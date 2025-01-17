# %% 
using CSV
using DataFrames
using Chain
using Optim, LineSearches
using Statistics
using LinearAlgebra
using SparseArrays
using Distributions
using CategoricalArrays

# %% 
# Estimate Factor Model Imputation TE function ---------------------------------
function m_theta_bar(
  params::AbstractVector,
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
    HΘ = [I(T - p); reshape(params, T - p, p)']
  else
    HΘ = Matrix(1.0 * I(T))
  end

  # This code assumes balanced panel that is sorted by [:fips, :year] 
  # This gives a huge speedup
  # Do first loop manually
  i = idx_control[1]
  prod_i = HΘ' * view(ytilde, (1 + (i - 1) * T):(i * T))
  # prod_i = HΘ' * ytilde[(1 + (i-1)*T):(i*T)]
  m_bar = kron(prod_i, W[:, i])

  # Units 2-n_\infty
  for i in idx_control[2:end]
    # Hand do the += kron
    mul!(prod_i, HΘ', view(ytilde, (1 + (i - 1) * T):(i * T)))

    m = 1
    for j in axes(prod_i, 1)
      for k in 1:n_instruments
        m_bar[m] += prod_i[j] * W[k, i]
        m += 1
      end
    end
  end

  m_bar *= P_D_inf
  # Divide by N for average
  m_bar ./= n_units

  return m_bar
end

function m_thetas(
  params::AbstractVector,
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
    theta = params[1:((T - p) * p)]
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

function obj_theta(
  params::AbstractVector,
  p::Int64,
  weight::Matrix{Float64},
  ytilde::Vector{Float64},
  idx_control::Vector{Int64},
  W::Matrix{Float64},
  n_units::Int64,
  T::Int64,
)
  m_bar = m_theta_bar(params, p, ytilde, idx_control, W, n_units, T)
  return n_units * (m_bar' * weight * m_bar)
end

function est_theta_hat(
  p::Int64,
  ytilde::Vector{Float64},
  idx_control::Vector{Int64},
  W::Matrix{Float64},
  n_units::Int64,
  T::Int64,
)
  n_instruments = size(W, 1) # Since I transposed W for quicker indexing

  params = ones((T - p) * p)
  weight = Matrix(1.0 * I((T - p) * n_instruments))

  if p > 0
    optres = optimize(
      x -> obj_theta(x, p, weight, ytilde, idx_control, W, n_units, T),
      params,
      BFGS(; linesearch=BackTracking()),
      Optim.Options(; g_abstol=eps(Float64));
      autodiff=:forward,
    )
    params_hat = Optim.minimizer(optres)

  else
    params_hat = params
  end

  ms = m_thetas(params_hat, p, ytilde, idx_control, W, n_units, T)
  Ωinv_theta = pinv((ms' * ms) / n_units)

  if p > 0
    # Having problems with numeric stability when doing all the 1/N_inf correctly
    RESCALE = 1000
    optres_opt = optimize(
      x -> 1 / RESCALE * obj_theta(x, p, Ωinv_theta, ytilde, idx_control, W, n_units, T),
      params_hat,
      BFGS(; linesearch=BackTracking()),
      Optim.Options(; iterations=10000);
      autodiff=:forward,
    )

    params_hat_opt = Optim.minimizer(optres_opt)
    J = Optim.minimum(optres_opt) * RESCALE
  else
    params_hat_opt = params
    J = obj_theta(params_hat_opt, p, Ωinv_theta, ytilde, idx_control, W, n_units, T)
  end

  # Testing for $p$ following Ahn, Lee, and Schmidt (2013)
  # Χ^2((T - p_0) (q - p_0) - k)
  # Where T is the number of periods, p_0 is the number of factors under the null, q is the number of instruments, and k is the number of additional covariates = 0
  k = 0
  if p == n_instruments
    p_value_hansen_sargent = 10
  else
    dist_hansen_sargent = Chisq((T - p) * (n_instruments - p) - k)
    p_value_hansen_sargent = 1 - cdf(dist_hansen_sargent, J)
  end

  return (params_hat_opt=params_hat_opt, J=J, p_value_hansen_sargent=p_value_hansen_sargent)
end

# %% Function to estimate a treatment effect
function est_F_qld(df; T0, p=-1)

  ## GMM estimate factor -------------------------------------------------------
  y = df.y
  id = df.id
  t = df.t
  T = maximum(df.t)

  # Find untreated units
  n_units = length(unique(df.id))
  idx_control = findall(df.treated[df.t .== 1] .== false)

  # Number of instruments
  instruments = [:W1, :W2]
  n_instruments = length(instruments)
  W = Matrix(df[df.t .== 1, instruments])
  W = Matrix(transpose(W))

  # Estimate theta -------------------------------------------------------------
  # Step 1. Initial estimate

  res = missing
  if p == -1
    p_star = -1

    for curr_p in 1:n_instruments
      # print("$(curr_p)")
      res = est_theta_hat(curr_p, y, idx_control, W, n_units, T)
      # print(": $(res.p_value_hansen_sargent)\n")

      p_star = curr_p
      if (res.p_value_hansen_sargent > 0.10) || (p + 1 == n_instruments)
        break
      end
    end

    p = p_star
  else
    res = est_theta_hat(p, y, idx_control, W, n_units, T)
  end
  theta_hat_opt = res.params_hat_opt
  FΘ = [reshape(theta_hat_opt, T - p, p); (-1 * I(p))]
  return FΘ
end
