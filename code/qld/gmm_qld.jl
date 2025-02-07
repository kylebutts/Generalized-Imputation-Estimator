"""
    m_theta_bar_fast(parameters, p, ymat_W)

These are the quasi-long differencing moments
1/N \\sum_i w_i ⊗ [H(Θ)' * y_i]

# Details
This uses some linear algebra tricks since we only care about the average of the unit moments.
By multiplying the ``T x N`` matrix of ``Y_{it}`` by the ``N x L`` matrix ``W``, 
we get a ``T x L`` matrix with element ``(t, l)`` consisting of ``\\sum_i y_{it} w_{i,l}``.

Then, we can multiply this by ``H(Θ)`` to get moments (instead of doing ``H(Θ) y_i`` first). 
"""
function mbar_theta(
  theta::AbstractVector,
  p::Int64, # Number of factors
  N::Int64,
  ymat_W::AbstractMatrix,
)
  if p < 0 || !isa(p, Int)
    error("p must be a non-negative integer")
  end
  T = size(ymat_W, 1)

  if p > 0
    # \sum_i H(\Tilde) y_i ⊗ w_i
    qld_ymat_W = ymat_W[1:(T - p), :] + reshape(theta, T - p, p) * ymat_W[(T - p + 1):T, :]
  elseif p == 0
    qld_ymat_W = ymat_W
  end

  # m_bar = 1 / N * 1 / (N_inf / N) * vec(transpose(qld_ymat_W))
  # m_bar = 1 / N_inf * vec(transpose(qld_ymat_W))
  m_bar = 1 / N * vec(transpose(qld_ymat_W))
  return m_bar
end

function obj_theta(
  theta::AbstractVector,
  weight::AbstractMatrix,
  p::Int64, # Number of factors
  N::Int64,
  ymat_W::AbstractMatrix,
)
  mbar = mbar_theta(theta, p, N, ymat_W)
  # Multiply by N^2 for Sargan J-statistic
  # Doing 1/N because it helps with numeric stability
  return 1 / N * (mbar' * weight * mbar)
end

# Matrix of m(z_i, θ)
function ms_theta(
  theta::AbstractVector,
  p::Int64, # Number of factors
  ymat::AbstractMatrix,
  W::AbstractMatrix,
  idx_control::AbstractVector,
)
  if p < 0 || !isa(p, Int)
    error("p must be a non-negative integer")
  end

  # Since I transposed W for quicker indexing
  N, N_instruments = size(W)
  T = size(ymat, 1)
  # N_inf = length(idx_control)

  # H(Θ)
  if p > 0
    Hθ = [I(T - p) reshape(theta, T - p, p)]
  else
    Hθ = Matrix(I(T))
  end

  ms = zeros(N_instruments * (T - p), N)
  for i in idx_control
    ms[:, i] .= kron(Hθ * ymat[:, i], W[i, :])
  end
  # ms *= 1 / (N_inf / N)
  return ms
end


"""
   gmm_qld(p, ymat, W, idx_control)

Estimates the quasi-long differencing estimator 

`ymat` is a T x N matrix with element y_{it}
`W` is a N x K matrix of instruments for each unit
`p` is the number of factors to estimate. `p` must be <= K.
`idx_control` denotes which columns of `ymat` and rows of `W` to use in estimation of the QLD problem

"""
function gmm_qld_p_known(
  p::Int64, # Number of factors
  ymat::AbstractMatrix,
  W::AbstractMatrix,
  idx_control::AbstractVector,
)
  T, N = size(ymat)
  N_instruments = size(W, 2)
  N_inf = length(idx_control)
  ymat_W_inf = ymat[:, idx_control] * W[idx_control, :]

  ## First-stage initial consistent estimator ----
  theta_initial = ones((T - p) * p)
  W_initial = Matrix(1.0 * I((T - p) * N_instruments))

  gmm_methods = [
    BFGS(; linesearch=LineSearches.BackTracking()), BFGS(), NewtonTrustRegion()
  ]
  if p > 0
    optres = nothing
    for method in gmm_methods
      optres = Optim.optimize(
        x -> obj_theta(x, W_initial, p, N_inf, ymat_W_inf),
        theta_initial,
        method,
        Optim.Options(; allow_f_increases=true, g_abstol=eps(Float64));
        autodiff=:forward,
      )
      if Optim.converged(optres)
        break
      end
    end

    if !Optim.converged(optres)
      error("Failed to converge in $(Optim.iterations(optres)) iterations")
    end
    theta_hat = Optim.minimizer(optres)
  else
    theta_hat = theta_initial
  end

  ## Calculate optimal weighting matrix
  ms_hat = ms_theta(theta_hat, p, ymat, W, idx_control)
  W_opt = pinv(1 / N_inf * ms_hat * ms_hat')

  ## Second-stage asymptotically efficient estimator ----
  if p > 0
    optres = nothing
    for method in gmm_methods
      optres = Optim.optimize(
        x -> obj_theta(x, W_opt, p, N_inf, ymat_W_inf),
        theta_hat,
        method,
        Optim.Options(; allow_f_increases=true, g_abstol=eps(Float64));
        autodiff=:forward,
      )
      if Optim.converged(optres)
        break
      end
    end

    if !Optim.converged(optres)
      error("Failed to converge in $(Optim.iterations(optres)) iterations")
    end
    theta_hat_opt = Optim.minimizer(optres)
  else
    theta_hat_opt = theta_hat
  end

  J = N_inf^2 * obj_theta(theta_hat_opt, W_opt, p, N_inf, ymat_W_inf)

  # Testing for $p$ following Ahn, Lee, and Schmidt (2013)
  # Χ^2((T - p_0) (q - p_0) - k)
  # Where T is the number of periods, p_0 is the number of factors under the null, q is the number of instruments, and k is the number of covariates
  k = 0
  if p < N_instruments
    dist_hansen_sargent = Distributions.Chisq((T - p) * (N_instruments - p) - k)
    p_value_hansen_sargent = 1.0 - Distributions.cdf(dist_hansen_sargent, J)
  else
    p_value_hansen_sargent = 1.0
  end

  Mbar_theta = ForwardDiff.jacobian(x -> mbar_theta(x, p, N_inf, ymat_W_inf), theta_hat_opt)

  # # Switch from conditional to unconditional moments
  # Mbar_theta *= 1 / (1 / N_inf) * (1 / N) * (1 / (N_inf / N))
  # Mbar_theta *= 1 / (1 / N_inf) * (1 / N_inf)
  # Mbar_theta *= 1

  # # Switch from conditional to unconditional moments
  # W_opt *= N / N_inf

  return theta_hat_opt, W_opt, Mbar_theta, J, p_value_hansen_sargent
end
