# ------------------------------------------------------------------------------
# Event Study Estimator for Walmart data 


using Random, Distributions, Statistics, LinearAlgebra
using Plots, LaTeXStrings, PrettyTables, GLM, CategoricalArrays
using CSV, DataFrames, Chain, FixedEffectModels, Printf
using Optim, LineSearches, RCall, Calculus, BlockDiagonals
import ForwardDiff, StatsBase

# Set parameters
YEARS = collect(1977:1999)
T0 = 1985
file = "data/sample_basker_YEARS_$(minimum(YEARS))_$(maximum(YEARS))_T0_$(T0).csv"
T = length(YEARS)

ests = [(p = 2, outcome = "retail"), (p = 1, outcome = "wholesale")]


for est in ests 

  outcome = est.outcome
  p = est.p
  outcome_var = "log_$(outcome)_emp"


  # Load data ------------------------------------------------------------------
  sample = DataFrame(CSV.File(file))
  sample.y = sample[:, outcome_var]
  sample.share_pop_ind_manuf = sample.share_pop_ind_manuf_nondurable + sample.share_pop_ind_manuf_durable

  # median retail employment
  sample[sample.year .== 1977, :retail_emp] |> median


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


  ## Estimate τ^ℓ --------------------------------------------------------------
  FΘ = [reshape(theta_hat_opt, T-p, p); (-1 * I(p))]
  FΘ = Matrix(FΘ)
  (FΘ * inv(FΘ' * FΘ) * FΘ')

  theta_hat_opt
  treated_ids = unique(sample[sample.g .< Inf, :fips])

  sample.ytilde_0_hat = zeros(size(sample, 1))

  for id in treated_ids
    unit = sample[(sample.fips .== id), :]
    g_id = Int(unit.g[1])
      
    y_pre = unit[(unit.year .< g_id), :ytilde]  
    FΘ_pre = FΘ[1:(g_id - minimum(YEARS)), :]

    sample[(sample.fips .== id), :ytilde_0_hat] = 
      FΘ * inv(FΘ_pre' * FΘ_pre) * FΘ_pre' * y_pre
  end

  sample.tau_hat = sample.ytilde .- sample.ytilde_0_hat



  collapsed = sample |>
    (x -> groupby(x, [:g, :year])) |>
    (x -> combine(x, 
      :y => mean
    ))

  # plot(
  #   collapsed.year, collapsed.y_mean, group = collapsed.g
  # )
    
  ## Synthetic Control estimates
  treated = sample[sample.g .!= Inf, :]

  treated.y_0_hat = 
    treated.ytilde_0_hat .+ treated.yinf_t .+ treated.yi_pre .- treated.yinf_pre

  collapsed = treated |>
    (x -> groupby(x, [:rel_year])) |>
    (x -> combine(x, 
      :y => mean,
      :y_0_hat => mean,
      :ytilde => mean,
      :ytilde_0_hat => mean
    ))

  CSV.write("data/synthetic_$(outcome)_p_$(p).csv", collapsed)

end











