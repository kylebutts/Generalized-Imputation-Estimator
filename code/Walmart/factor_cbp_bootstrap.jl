# Event Study Estimator for Walmart data 

# %% 
#| label: setup
using CSV, DataFrames
using FixedEffectModels, CategoricalArrays, GLM
using Optim, LineSearches
using ForwardDiff
using StatsBase, Statistics, Random, Distributions
using LinearAlgebra, Calculus
using SparseArrays, BlockDiagonals
using PrettyTables, Printf
using Plots
using Dates

# Set parameters
YEARS = collect(1977:1999)
T = length(YEARS)
T0 = 1985

rel_years = begin
  x = sample[sample.rel_year .!= -Inf, :rel_year]
  collect(minimum(x):maximum(x))
end

# %%
file = "data/County_Business_Patterns/sample_basker_YEARS_$(minimum(YEARS))_$(maximum(YEARS))_T0_$(T0).csv"
sample = DataFrame(CSV.File(file))

# %%
# Order by fips, then by year
sort!(sample, [:fips, :year])

# %% Median employment
median(sample[sample.year .== 1977, :retail_emp])
median(sample[sample.year .== 1977, :wholesale_emp])

# Estimate Factor Model Imputation TE function ---------------------------------
# %% New `m_thetas`
function m_thetas(
  parms::AbstractVector,
  p::Int64,
  ytilde::Vector{Float64},
  idx_control::Vector{Int64},
  W::Matrix{Float64},
  n_units::Int64,
  T::Int64,
)

  # Since I transposed W for quicker indexing
  n_instruments = size(W, 1)
  N_inf = length(idx_control)
  P_D_inf = (N_inf / n_units)

  # H(Θ)
  if p > 0
    theta = parms[1:((T - p) * p)]
    HΘ = [sparse(I, T - p, T - p); reshape(theta, T - p, p)']
  else
    HΘ = sparse(I, T, T)
  end

  # Initialize ms
  ms = zeros(Float64, n_units, (T - p) * n_instruments)

  # This code assumes balanced panel that is sorted by [:fips, :year] 
  # This gives a huge speedup
  for i in idx_control
    HΘytilde_i = (HΘ' * ytilde[(1 + (i - 1) * T):(i * T)])
    W_i = W[:, i]

    ms[i, :] = kron(HΘytilde_i, W_i) * P_D_inf
  end

  return ms
end

# %% New `m_theta_bar`
function m_theta_bar(
  parms::AbstractVector,
  p::Int64,
  ytilde::Vector{Float64},
  idx_control::Vector{Int64},
  W::Matrix{Float64},
  n_units::Int64,
  T::Int64,
)

  # Since I transposed W for quicker indexing
  n_instruments = size(W, 1)

  # P(D_\infty = 1)
  N_inf = length(idx_control)
  P_D_inf = (N_inf / n_units)

  # H(Θ)
  if p > 0
    theta = parms[1:((T - p) * p)]
    HΘ = [I(T - p); reshape(theta, T - p, p)']
  else
    HΘ = Matrix(1.0 * I(T))
  end

  # This code assumes balanced panel that is sorted by [:fips, :year] 
  # This gives a huge speedup
  # Do first loop manually
  i = idx_control[1]
  prod_i = HΘ' * view(ytilde, (1 + (i - 1) * T):(i * T))
  # prod_i = HΘ' * ytilde[(1 + (i-1)*T):(i*T)]
  m_bar = kron(prod_i, W[:, i]) * P_D_inf

  # Units 2-n_\infty
  for i in idx_control[2:end]
    # Hand do the += kron
    mul!(prod_i, HΘ', view(ytilde, (1 + (i - 1) * T):(i * T)))

    m = 1
    for j in axes(prod_i, 1)
      for k in 1:n_instruments
        m_bar[m] += prod_i[j] * W[k, i] * P_D_inf
        m += 1
      end
    end
  end

  # Divide by N for average
  m_bar ./= n_units

  return m_bar
end

# %% Objective function
function obj_theta(
  parms::AbstractVector,
  p::Int64,
  weight::Matrix{Float64},
  ytilde::Vector{Float64},
  idx_control::Vector{Int64},
  W::Matrix{Float64},
  n_units::Int64,
  T::Int64,
)
  m_bar = m_theta_bar(parms, p, ytilde, idx_control, W, n_units, T)
  return n_units * (m_bar' * weight * m_bar)
end

# %% Brown and Butts - Within transformation
function within_transform(df, outcome_var, g_var, t_var, T0)
  sample = copy(df)

  sample.y = sample[!, outcome_var]

  # Within Transformation ------------------------------------------------------
  # never-treated cross-sectional averages of y
  sample = (x -> transform(x, [:y, g_var] => ((y, g) -> mean(y[g .== Inf])) => :yinf_t))(
    (x -> groupby(x, :year))(sample)
  )

  # unit-averages pre-T0
  sample = (
    x -> transform(x, [:y, t_var] => ((y, year) -> mean(y[year .<= T0])) => :yi_pre)
  )(
    (x -> groupby(x, :fips))(sample)
  )

  # never-treated average pre-T0
  sample.yinf_pre .= mean(
    sample[(sample[!, g_var] .== Inf) .* (sample[!, t_var] .<= T0), :y]
  )

  # Create ytilde
  ytilde = sample.y - sample.yinf_t - sample.yi_pre + sample.yinf_pre

  return ytilde
end

# %% Function to estimate a treatment effect
function est_te(sample, p, outcome_var; export_y0=false, outcome="", quick_plot=false)
  sort!(sample, [:fips, :year])
  ytilde = within_transform(sample, outcome_var, :g, :year, T0)
  sample.ytilde = ytilde

  ## GMM estimate factor -------------------------------------------------------
  instruments = [
    :share_pop_ind_manuf,
    :share_pop_poverty_78_below,
    :share_pop_poverty_78_above,
    :share_pop_emp_private,
    :share_pop_emp_government,
    :share_school_col,
    :share_school_hs,
  ]
  ytilde = sample.ytilde
  fips = sample.fips
  year = sample.year
  rel_years = sort(unique(sample[(sample.g .!= Inf), :rel_year]))

  # Find untreated units
  n_units = length(unique(fips))
  idx_control = findall(sample.g[year .== 1977] .== Inf)
  N_inf = length(idx_control)

  # Number of instruments
  n_instruments = length(instruments)
  W = Matrix(sample[year .== 1977, instruments])
  W = Matrix(transpose(W))

  # Number of covariates
  k = 0

  # Estimate theta -------------------------------------------------------------
  # Step 1. Initial estimate
  parms = [ones((T - p) * p); zeros(length(rel_years))]
  weight = Matrix(1.0 * I((T - p) * n_instruments))

  if p > 0
    optres = optimize(
      x -> obj_theta(x, p, weight, ytilde, idx_control, W, n_units, T),
      parms,
      BFGS(; linesearch=BackTracking()),
      Optim.Options(; g_abstol=eps(Float64));
      autodiff=:forward,
    )

    parms_hat = Optim.minimizer(optres)
  else
    parms_hat = parms
  end

  # Step 2. Optimal weighting matrix 
  ms = m_thetas(parms_hat, p, ytilde, idx_control, W, n_units, T)
  Ωinv_theta = pinv(1 / n_units * (ms' * ms))

  # Ωinv_theta = zeros(147, 147)
  # for i in 1:n_units
  #   Ωinv_theta += ms[i, :] * ms[i, :]'
  # end
  # Ωinv_theta = 1 / n_units * Ωinv_theta

  if p > 0
    optres_opt = optimize(
      x -> obj_theta(x, p, Ωinv_theta, ytilde, idx_control, W, n_units, T),
      parms_hat,
      BFGS(; linesearch=BackTracking()),
      Optim.Options(; g_abstol=eps(Float64));
      autodiff=:forward,
    )

    parms_hat_opt = Optim.minimizer(optres_opt)
    theta_hat_opt = parms_hat_opt[1:((T - p) * p)]
    J = Optim.minimum(optres_opt)
  else
    J = obj_theta(parms_hat, p, Ωinv_theta, ytilde, idx_control, W, n_units, T)
  end

  # Testing for $p$ following Ahn, Lee, and Schmidt (2013)
  # Χ^2((T - p_0) (q - p_0) - k)
  # Where T is the number of periods, p_0 is the number of factors under the null, q is the number of instruments, and k is the number of covariates
  dist_hansen_sargent = Chisq((T - p) * (n_instruments - p) - k)
  p_value_hansen_sargent = 1 - cdf(dist_hansen_sargent, J)

  # Impute untreated potential outcomes
  sample.y0hat = zeros(size(sample, 1))
  if p > 0
    FΘ = [reshape(theta_hat_opt, T - p, p); (-1 * I(p))]
    FΘ = Matrix(FΘ)
    (FΘ * inv(FΘ' * FΘ) * FΘ')

    theta_hat_opt
    treated_ids = unique(sample[sample.g .< Inf, :fips])

    for id in treated_ids
      unit = sample[(sample.fips .== id), :]
      g_id = Int(unit.g[1])

      y_pre = unit[(unit.year .< g_id), :ytilde]
      FΘ_pre = FΘ[1:(g_id - minimum(YEARS)), :]

      sample[(sample.fips .== id), :y0hat] = FΘ * inv(FΘ_pre' * FΘ_pre) * FΘ_pre' * y_pre
    end
  end

  # Estimate τ^ℓ ---------------------------------------------------------------
  sample.tau_hat = sample.ytilde .- sample.y0hat
  if export_y0 == true
    outfile = "estimates/est_y0_outcome_$(outcome)_qld_p_$(p).csv"
    CSV.write(outfile, sample[:, [:fips, :year, :g, :rel_year, :ytilde, :y0hat, :tau_hat]])
  end

  # Factor Estimate
  sample.rel_year_cat = CategoricalArray(sample.rel_year)
  lm_tau = lm(@formula(tau_hat ~ 0 + rel_year_cat), sample[sample.g .< Inf, :])

  # Plot point estimates
  if quick_plot == true
    plot(sort(rel_years), coef(lm_tau); marker=:circle)
  end

  return (
    est=[coef(lm_tau), stderror(lm_tau)],
    sargan=(J, p_value_hansen_sargent),
    imputation=[sample.ytilde, sample.y0hat],
  )
end

# %% 
# p = 2
# outcome_var = :log_retail_emp

# %% Estimate
tau_hat_retail_p_0 = est_te(
  sample, 0, :log_retail_emp; export_y0=true, outcome="retail", quick_plot=true
)
tau_hat_retail_p_1 = est_te(
  sample, 1, :log_retail_emp; export_y0=true, outcome="retail", quick_plot=true
)
tau_hat_retail_p_2 = est_te(
  sample, 2, :log_retail_emp; export_y0=true, outcome="retail", quick_plot=true
)
tau_hat_retail_p_3 = est_te(
  sample, 3, :log_retail_emp; export_y0=true, outcome="retail", quick_plot=true
)
tau_hat_retail_p_4 = est_te(
  sample, 4, :log_retail_emp; export_y0=true, outcome="retail", quick_plot=true
)

# %% 
tau_hat_wholesale_p_0 = est_te(
  sample, 0, :log_wholesale_emp; export_y0=true, outcome="wholesale"
)
tau_hat_wholesale_p_1 = est_te(
  sample, 1, :log_wholesale_emp; export_y0=true, outcome="wholesale"
)
tau_hat_wholesale_p_2 = est_te(
  sample, 2, :log_wholesale_emp; export_y0=true, outcome="wholesale"
)
tau_hat_wholesale_p_3 = est_te(
  sample, 3, :log_wholesale_emp; export_y0=true, outcome="wholesale"
)
tau_hat_wholesale_p_4 = est_te(
  sample, 4, :log_wholesale_emp; export_y0=true, outcome="wholesale"
)

# %% 
# Ahn, Lee, and Schmidt procedure for selecting p
println("""
Checking for p for retail employment: 
p = 0; Hansen-Sargant Test p-value = $(tau_hat_retail_p_0.sargan[2])
p = 1; Hansen-Sargant Test p-value = $(tau_hat_retail_p_1.sargan[2])
p = 2; Hansen-Sargant Test p-value = $(tau_hat_retail_p_2.sargan[2])
p = 3; Hansen-Sargant Test p-value = $(tau_hat_retail_p_3.sargan[2])
p = 4; Hansen-Sargant Test p-value = $(tau_hat_retail_p_4.sargan[2])

Checking for p for wholesale employment: 
p = 0; Hansen-Sargant Test p-value = $(tau_hat_wholesale_p_0.sargan[2])
p = 1; Hansen-Sargant Test p-value = $(tau_hat_wholesale_p_1.sargan[2])
p = 2; Hansen-Sargant Test p-value = $(tau_hat_wholesale_p_2.sargan[2])
p = 3; Hansen-Sargant Test p-value = $(tau_hat_wholesale_p_3.sargan[2])
p = 4; Hansen-Sargant Test p-value = $(tau_hat_wholesale_p_4.sargan[2])
""")

# %%
plot(rel_years, tau_hat_wholesale_p_1.est[1]; marker=:circle)

# %% Chose p based on p-values
retail_p_values = (
  tau_hat_retail_p_0.sargan[2],
  tau_hat_retail_p_1.sargan[2],
  tau_hat_retail_p_2.sargan[2],
  tau_hat_retail_p_3.sargan[2],
  tau_hat_retail_p_4.sargan[2],
)
retail_p = findfirst(retail_p_values .> 0.1) - 1

retail_est = eval(Symbol("tau_hat_retail_p_$(retail_p)"))

wholesale_p_values = (
  tau_hat_wholesale_p_0.sargan[2],
  tau_hat_wholesale_p_1.sargan[2],
  tau_hat_wholesale_p_2.sargan[2],
  tau_hat_wholesale_p_3.sargan[2],
  tau_hat_wholesale_p_4.sargan[2],
)
wholesale_p = findfirst(wholesale_p_values .> 0.1) - 1

wholesale_est = eval(Symbol("tau_hat_wholesale_p_$(wholesale_p)"))

# %% 
plot(rel_years, retail_est.est[1]; marker=:circle)
plot(rel_years, wholesale_est.est[1]; marker=:circle)

# %% Export naive standard errors (ignoring estimation of factors)
est_retail_naive_se = DataFrame(;
  rel_year=collect(rel_years), estimate=retail_est.est[1], std_error=retail_est.est[2]
)
CSV.write("estimates/qld_p_$(retail_p)_est_retail_naive_se.csv", est_retail_naive_se)

est_wholesale_naive_se = DataFrame(;
  rel_year=collect(rel_years), estimate=wholesale_est.est[1], std_error=wholesale_est.est[2]
)
CSV.write(
  "estimates/qld_p_$(wholesale_p)_est_wholesale_naive_se.csv", est_wholesale_naive_se
)

# %% Retail bootstrap (takes about 10 minutes)
# display current time
start_time = Dates.now()
println("Current Time: $(Dates.format(start_time, "HH:MM:SS"))")
B = 5000
ids = unique(sample[:, :fips])
retail_ests = zeros(B, length(rel_years))
retail_ests[1, :] = retail_est.est[1]'

# In case estimation fails
b = 2
while b < B
  if (b % 50 == 0)
    println("On Bootstrap Iteration: $(b)")
  end

  resample_ids = StatsBase.sample(ids, length(ids))

  resample_df = vcat([sample[sample.fips .== id, :] for id in resample_ids]...)

  # Assuming balanced
  resample_df.fips = repeat(1:length(ids); inner=T)

  try
    boot_est = est_te(
      resample_df,
      retail_p,
      :log_retail_emp;
      export_y0=false,
      outcome="retail",
      quick_plot=false,
    )
    retail_ests[b, :] = boot_est.est[1]'
    b += 1
  catch
    println("Error on bootstrap iteration: $(b)")
  end
end

# Export 1000 bootstraps
CSV.write(
  "estimates/qld_p_$(retail_p)_est_retail_B_$(B).csv",
  DataFrame(retail_ests, ["tau$(i)" for i in rel_years]),
)
end_time = Dates.now()
print("Current Time: $(Dates.format(end_time, "HH:MM:SS"))")

# %% Wholesale bootstrap (takes about 10 minutes)

# display current time
start_time = Dates.now()
println("Current Time: $(Dates.format(start_time, "HH:MM:SS"))")
B = 5000
ids = unique(sample[:, :fips])
wholesale_ests = zeros(B, length(rel_years))
wholesale_ests[1, :] = wholesale_est.est[1]'

# In case estimation fails
b = 2
while b < B
  if (b % 50 == 0)
    println("On Bootstrap Iteration: $(b)")
  end

  resample_ids = StatsBase.sample(ids, length(ids))

  resample_df = vcat([sample[sample.fips .== id, :] for id in resample_ids]...)

  # Assuming balanced
  resample_df.fips = repeat(1:length(ids); inner=T)

  # est_te(
  #   resample_df, wholesale_p, :log_wholesale_emp,
  #   export_y0=false, outcome="wholesale", quick_plot=false
  # ).est[1]

  try
    boot_est = est_te(
      resample_df,
      wholesale_p,
      :log_wholesale_emp;
      export_y0=false,
      outcome="wholesale",
      quick_plot=false,
    )
    # If any is NaN
    # if !(any(isnan(boot_est.est[1])))
    wholesale_ests[b, :] = boot_est.est[1]'
    b += 1
    # end
  catch
    println("Error on bootstrap iteration: $(b)")
  end
end

# Export 1000 bootstraps
CSV.write(
  "estimates/qld_p_$(wholesale_p)_est_wholesale_B_$(B).csv",
  DataFrame(wholesale_ests, ["tau$(i)" for i in rel_years]),
)
end_time = Dates.now()
print("Current Time: $(Dates.format(end_time, "HH:MM:SS"))")
