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
file = "data/County_Business_Patterns/sample_basker_YEARS_$(minimum(YEARS))_$(maximum(YEARS))_T0_$(T0).csv"
T = length(YEARS)

# Note p was selected by Ahn, Lee, and Schmidt test by calling est_te for each p and looking at the objective function value
outcome = "wholesale"
p = 1
# outcome = "retail"
# p = 3
outcome_var = "log_$(outcome)_emp"


# Load data --------------------------------------------------------------------
sample = DataFrame(CSV.File(file))
sample.y = sample[:, outcome_var]
sample.share_pop_ind_manuf = sample.share_pop_ind_manuf_nondurable + sample.share_pop_ind_manuf_durable

# median retail employment
sample[sample.year .== 1977, :retail_emp] |> median
sample[sample.year .== 1977, :wholesale_emp] |> median


# Estimate Factor Model Imputation TE function ---------------------------------
function est_te(sample; return_se = false, export_y0 = false)
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

  # Test for p following Ahn, Lee, and Schmidt (2013)
  # obj_value = Optim.minimum(optres)
  # print("p = $(p); objective function value = $(obj_value)")

  # Step 2. Optimal weighting matrix 
  Ωinv_theta = pinv( (m_thetas(parms_hat, p)' * m_thetas(parms_hat, p)) / n_units ) 

  if p > 0
    optres_opt = optimize(
      x -> obj_theta(x, Ωinv_theta, sample.ytilde, Matrix(sample[:, instruments]), sample.g, p), parms_hat, 
      SimulatedAnnealing(),
      Optim.Options(iterations=100000)
      # BFGS(; linesearch = BackTracking()), 
      # Optim.Options(iterations=10000); 
      # autodiff = :forward
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

  # Testing for $p$ following Ahn, Lee, and Schmidt (2013)
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

  if export_y0 == true
    outfile = "estimates/est_y0_outcome_$(outcome)_qld_p_$(p).csv"
    CSV.write(
      outfile, sample[:, [:fips, :year, :ytilde, :y0hat, :tau_hat]]
    )
  end

  # Factor Estimate
  sample.rel_year_cat = CategoricalArray(sample.rel_year)
  lm_tau = lm(
    @formula(tau_hat ~ 0 + rel_year_cat), 
    sample[sample.g .< Inf, :]
  )

  if return_se
    return coef(lm_tau), stderror(lm_tau)
  else
    return coef(lm_tau)
  end
end

tau_hat = est_te(sample, return_se = true, export_y0 = true)

est_naive_se = DataFrame(
  rel_year = collect(-22:13),
  estimate = tau_hat[1],
  std_error = tau_hat[2]
)
CSV.write(
  "data/qld_p_$(p)_est_$(outcome)_naive_se.csv", 
  est_naive_se
)

B = 1000
ids = sample[:, :fips] |> unique
ests = zeros(B, 36)
ests[1, :] = tau_hat[1]

# Bootstrap standard errors ----------------------------------------------------
# Starting at 2 to not write over point estimate
for i in 2:B
  println("On Bootstrap Iteration: $(i)")

  resample_ids = StatsBase.sample(ids, length(ids))
  
  temp = sample[sample.fips .== resample_ids[1], :]
  temp.fips .= 1
  resample_df = temp
  for j in 2:eachindex(resample_ids)
    temp = sample[sample.fips .== resample_ids[j], :]
    temp.fips .= j
    resample_df = [resample_df; temp]
  end

  ests[i, :] = est_te(resample_df)'
end

# Export 1000 bootstraps
CSV.write(
  "data/qld_p_$(p)_est_$(outcome)_B_$(B).csv", 
  DataFrame(ests, ["tau$(i)" for i in -22:13])
)

