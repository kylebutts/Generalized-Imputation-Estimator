# ------------------------------------------------------------------------------
# Event Study Estimator for Walmart data 


using Random, Distributions, Statistics, LinearAlgebra
using Plots, PrettyTables, GLM, CategoricalArrays
using CSV, DataFrames, Chain, FixedEffectModels, Printf
using Optim, LineSearches, RCall, Calculus, BlockDiagonals
import ForwardDiff, StatsBase

# Set parameters
YEARS = collect(1977:1999)
T0 = 1985
file = "data/sample_basker_YEARS_$(minimum(YEARS))_$(maximum(YEARS))_T0_$(T0).csv"
T = length(YEARS)
p = 1
outcome = "wholesale"
outcome_var = "log_$(outcome)_emp"


# Load data --------------------------------------------------------------------
sample = DataFrame(CSV.File(file))
sample.y = sample[:, outcome_var]
sample.share_pop_ind_manuf = sample.share_pop_ind_manuf_nondurable + sample.share_pop_ind_manuf_durable

# median retail employment
sample[sample.year .== 1977, :retail_emp] |> median
sample[sample.year .== 1977, :wholesale_emp] |> median


# Estimate Factor Model Imputation TE function ---------------------------------
function est_te(sample)
  # Within Transformation ------------------------------------------------------

  # never-treated cross-sectional averages of y
  sample = sample |>
  (x -> groupby(x, :year)) |>
  (x -> transform(x,
    [:y, :g] => ((y, g) -> mean(y[g .== Inf])) => :yinf_t
  ))

  # unit-averages pre-T0
  sample = sample |>
  x -> groupby(x, :fips) |>
  x -> transform(x,
    [:y, :year] => ((y, year) -> mean(y[year .<= T0])) => :yi_pre
  )

  # never-treated average pre-T0
  sample.yinf_pre .= mean(sample[(sample.g .== Inf) .* (sample.year .<= T0), :y])

  # Create ytilde
  sample.ytilde = sample.y - sample.yinf_t - sample.yi_pre + sample.yinf_pre


  ## GMM estimate factor -------------------------------------------------------
  instruments = [
    :share_pop_ind_manuf,
    :share_pop_poverty_78_below, :share_pop_poverty_78_above, 
    :share_pop_emp_private, :share_pop_emp_government, 
    :share_school_col, :share_school_hs
  ]

  n_instruments = length(instruments)
  n_units = sample.fips |> unique |> length
  rel_years = sample[(sample.g .!= Inf), :rel_year] |> unique |> sort

  # Estimate theta -------------------------------------------------------------

  function m_thetas(parms::AbstractVector, p::Int64)
    
    if p > 0
      theta = parms[1:((T - p) * p)] 
      HΘ = [I(T - p); reshape(theta, T-p, p)']
    else 
      HΘ = 1.0 * I(T) |> Matrix
    end

    m_theta = zeros(n_units, (T - p) * n_instruments)
    N_inf = 0
    i = 1
    for id in unique(sample.fips)
      W = sample[sample.fips .== id, instruments][1, :] |> Vector
      ytilde = sample[sample.fips .== id, :ytilde]

      D_inf = sample[sample.fips .== id, :g] |> 
        (x -> (x[1] == Inf))
      
      N_inf += D_inf

      m_theta[i, :] = D_inf * kron((HΘ' * ytilde), W)
      i += 1
    end

    # divide by 1 / P(D_i,inf) 
    m_theta = m_theta / (N_inf / n_units)

    return m_theta
  end

  function m_theta_bar(parms::AbstractVector, ytilde::Vector{Float64}, W::Matrix{Float64}, g::Vector{Float64}, p::Int64)

    if p > Int64(0)
      # theta = parms[1:((T - p) * p)] 
      HΘ = [
        I(T - p); 
        reshape(parms[1:((T - p) * p)], T-p, p)'
      ]
    else 
      HΘ = Matrix(1.0 * I(T))
    end

    m_theta = zeros(Float64, (T - p) * n_instruments)  
    N_inf::Int64 = 0
    
    for i in 1:n_units
      D_inf = (g[1 + (i-1) * T] == Inf)

      if D_inf
        N_inf += 1

        m_theta += kron(
          (HΘ' * ytilde[(1 + (i-1) * T):(i * T)]), 
          W[1 + (i-1)*T, :]
        )
      end
    end

    # 1 / P(D_i,inf) 
    # m_theta = m_theta / (N_inf / n_units)

    # return the avearage value of m_theta, m_theta_bar
    # m_theta = m_theta / n_units

    return m_theta / N_inf
  end

  function obj_theta(parms::AbstractVector, weight::Matrix{Float64}, ytilde::Vector{Float64}, W::Matrix{Float64}, g::Vector{Float64}, p::Int64)
    m_bar = m_theta_bar(parms, ytilde, W, g, p)
    return (m_bar' * weight * m_bar)
  end

  # Step 1. Initial estimate
  parms = [ones((T - p) * p); zeros(length(rel_years))]
  weight = 1.0 * I((T - p) * n_instruments) |> Matrix

  if p > 0
    optres = optimize(
      x -> obj_theta(x, weight, sample.ytilde, Matrix(sample[:, instruments]), sample.g, p), parms, 
      LBFGS(); autodiff = :forward
    )

    parms_hat = Optim.minimizer(optres)
  else 
    obj_theta(parms, weight, sample.ytilde, Matrix(sample[:, instruments]), sample.g, p)
  end

  # obj_value = Optim.minimum(optres)
  # print("p = $(p); objective function value = $(obj_value)")

  # Step 2. Optimal weighting matrix 
  Ωinv_theta = pinv( (m_thetas(parms_hat, p)' * m_thetas(parms_hat, p)) / n_units ) 

  if p > 0
    optres_opt = optimize(
      x -> obj_theta(x, Ωinv_theta, sample.ytilde, Matrix(sample[:, instruments]), sample.g, p), parms_hat, 
      BFGS(; linesearch = BackTracking()), Optim.Options(iterations=10000); 
      autodiff = :forward
    )

    parms_hat_opt = Optim.minimizer(optres_opt)
    theta_hat_opt = parms_hat_opt[1:((T - p) * p)] 
    J = Optim.minimum(optres_opt)
  else 
    J = obj_theta(parms, Ωinv_theta, sample.ytilde, Matrix(sample[:, instruments]), sample.g, p)
  end

  Optim.minimum(optres_opt)
  parms_hat_opt = Optim.minimizer(optres_opt)
  # parms_hat_opt = Optim.minimizer(optres)
  theta_hat_opt = parms_hat_opt[1:((T - p) * p)] 

  # obj_value = Optim.minimum(optres_opt)
  # println("p = $(p); objective function value = $(obj_value)")


  # Estimate τ^ℓ ---------------------------------------------------------------
  FΘ = [reshape(theta_hat_opt, T-p, p); (-1 * I(p))]
  FΘ = Matrix(FΘ)
  (FΘ * inv(FΘ' * FΘ) * FΘ')

  theta_hat_opt
  treated_ids = unique(sample[sample.g .< Inf, :fips])

  sample.y0hat = zeros(size(sample, 1))

  for id in treated_ids
    unit = sample[(sample.fips .== id), :]
    g_id = Int(unit.g[1])
      
    y_pre = unit[(unit.year .< g_id), :ytilde]  
    FΘ_pre = FΘ[1:(g_id - minimum(YEARS)), :]

    sample[(sample.fips .== id), :y0hat] = 
      FΘ * inv(FΘ_pre' * FΘ_pre) * FΘ_pre' * y_pre
  end

  sample.tau_hat = sample.ytilde .- sample.y0hat

  # ATT(1986, t)
  # sample[(sample.g .== 1986), :] |>
  #   (x -> groupby(x, :year)) |>
  #   (x -> combine(x, :tau_hat => mean)) |>
  #   (x -> x.tau_hat_mean)

  # Factor Estimate
  sample.rel_year_cat = CategoricalArray(sample.rel_year)
  lm_tau = lm(
    @formula(tau_hat ~ 0 + rel_year_cat), 
    sample[sample.g .< Inf, :]
  )


  return coef(lm_tau)
end


tau_hat = est_te(sample)

B = 1000
ids = sample[:, :fips] |> unique
ests = zeros(B, 36)
ests[1, :] = tau_hat

# Starting at 2 to not write over point estimate
for i in 2:B
  println("On Bootstrap Iteration: $(i)")

  resample_ids = StatsBase.sample(ids, length(ids))
  
  temp = sample[sample.fips .== resample_ids[1], :]
  temp.fips .= 1
  resample_df = temp
  for j in 2:length(resample_ids)
    temp = sample[sample.fips .== resample_ids[j], :]
    temp.fips .= j
    resample_df = [resample_df; temp]
  end

  ests[i, :] = est_te(resample_df)'
end

# Export 1000 bootstraps
CSV.write(
  "data/factor_est_$(outcome)_YEARS_$(minimum(YEARS))_$(maximum(YEARS))_T0_$(T0)_B_$(B)_p_$(p).csv", 
  DataFrame(ests, ["tau$(i)" for i in -22:13])
)

print("done")


# se_bootstrap = std(ests, dims = 1)
upper = zeros(size(ests)[2])
lower = zeros(size(ests)[2])
se_bootstrap = zeros(size(ests)[2])
for i in 1:size(ests)[2]
  println(i)
  col = ests[:, i]
  sort!(col)

  upper[i] = quantile(col, 0.975)
  lower[i] = quantile(col, 0.025)
  se_bootstrap[i] = std(col)
end





# Helpers ----------------------------------------------------------------------

est_did2s = DataFrame(CSV.File("data/did2s_est_$(outcome)_YEARS_1977_1999_T0_1985_B_1000.csv"))



# Plot TE Estimates ------------------------------------------------------------


# Extract point estimates 
rel_year = sample[sample.g .!= Inf, :rel_year] |> 
  unique |> 
  sort

did2s_est = est_did2s.estimate
did2s_errors = est_did2s.std_error
did2s_upper = est_did2s.upper
did2s_lower = est_did2s.lower
# factor_est = coef(lm_tau)
# factor_errors = coeftable(lm_tau).cols[2]
factor_est = tau_hat
factor_errors = se_bootstrap

# Linear extrapolation
pre_ests = DataFrame(
  rel_year   = rel_year[(idx_pre) .& (rel_year .>= -15)], 
  did2s_est  = did2s_est[(idx_pre) .& (rel_year .>= -15)],
  factor_est = factor_est[(idx_pre) .& (rel_year .>= -15)]
)
did2s_extrapolation = lm(
  @formula(did2s_est ~ 1 + rel_year), 
  pre_ests
) |> coef
factor_extrapolation = lm(
  @formula(factor_est ~ 1 + rel_year), 
  pre_ests
) |> coef


# idx_pre  = (rel_year .<  0) .& (rel_year .>= -15)
# idx_post = (rel_year .>= 0) .& (rel_year .<  9)
idx_pre  = (rel_year .<  0)
idx_post = (rel_year .>= 0)

pgfplotsx()
# gr(size = (600, 300))

did2s_plot = hline([0], linestyle = :dash, linecolor = :black)
# title!("TWFE Model Imputation")
xlabel!("Event Time")
ylabel!("Coefficient")
ylims!(-0.4, 0.1)
Plots.abline!(
  did2s_extrapolation[2], did2s_extrapolation[1],
  color = RGBA(154/255, 36/255, 21/255, 1)
)
scatter!(
  rel_year[idx_pre],
  did2s_est[idx_pre], 
  # yerror = 1.96 .* did2s_errors[idx_pre],
  yerror = (
    did2s_est[idx_pre] - did2s_lower[idx_pre], 
    did2s_upper[idx_pre] - did2s_est[idx_pre]
  ),
  legend = false,
  markersize = 2, 
  markerstrokecolor = RGBA(154/255, 36/255, 21/255, 1),
  markercolor = RGBA(154/255, 36/255, 21/255, 1)
)
scatter!(
  rel_year[idx_post],
  did2s_est[idx_post],
  # yerror = 1.96 .* did2s_errors[idx_post],
  yerror = (
    did2s_est[idx_post] - did2s_lower[idx_post], 
    did2s_upper[idx_post] - did2s_est[idx_post]
  ),
  legend = false,
  markersize = 2, 
  markerstrokecolor = RGBA(16/255, 120/255, 149/255, 1),
  markercolor = RGBA(16/255, 120/255, 149/255, 1)
)


factor_plot = hline([0], linestyle = :dash, linecolor = :black)
# title!("Factor Model Imputation")
xlabel!("Event Time")
ylabel!("Coefficient")
ylims!(-0.4, 0.1)
Plots.abline!(
  factor_extrapolation[2], factor_extrapolation[1],
  color = RGBA(154/255, 36/255, 21/255, 1)
)
scatter!(
  rel_year[idx_pre],
  factor_est[idx_pre],
  yerror = (
    factor_est[idx_pre] - lower[idx_pre], 
    upper[idx_pre] - factor_est[idx_pre]
  ),
  # yerror = 1.96 .* se_bootstrap[idx_pre],
  legend = false,
  markersize = 2, 
  markerstrokecolor = RGBA(154/255, 36/255, 21/255, 1),
  markercolor = RGBA(154/255, 36/255, 21/255, 1)
)
scatter!(
  rel_year[idx_post],
  factor_est[idx_post],
  yerror = (
    factor_est[idx_post] - lower[idx_post], 
    upper[idx_post] - factor_est[idx_post]
  ),
  # yerror = 1.96 .* se_bootstrap[idx_post],
  legend = false,
  markersize = 2, 
  markerstrokecolor = RGBA(16/255, 120/255, 149/255, 1),
  markercolor = RGBA(16/255, 120/255, 149/255, 1)
)

print("plotting YEARS: $(minimum(YEARS)):$(maximum(YEARS)); T0: $(T0); p: $(p)\n")


did2s_plot = plot(
  # did2s_plot, size = (600, 300) .* 0.75
  did2s_plot, size = (600, 300) .* 1.25
)
factor_plot = plot(
  # factor_plot, size = (600, 300) .* 0.75
  factor_plot, size = (600, 300) .* 1.25
)

combined_plot = plot(
  did2s_plot, factor_plot, layout = (1,2), legend = false, size = (1200, 400)
) 

# savefig(did2s_plot, "figures/did2s_$(outcome)_p_$(p)_bootstrap_$(B).tex")
# savefig(factor_plot, "figures/factor_$(outcome)_p_$(p)_bootstrap_$(B).tex")
savefig(did2s_plot, "figures/did2s_$(outcome)_p_$(p)_bootstrap_$(B).tex")
savefig(factor_plot, "figures/factor_$(outcome)_p_$(p)_bootstrap_$(B).tex")


