# QLD Moment Functions ---------------------------------------------------------
"""
    m_theta_bar_fast(parameters, p, ymat_W)

These are the quasi-long differencing moments
1/N \\sum_i w_i ⊗ [H(Θ)' * y_i]

Note you must call `vec(transpose(m_theta_bar_fast(...)))` to get the correct moment vector.

# Details
This uses some linear algebra tricks since we only care about the average of the unit moments.
By multiplying the ``T x N`` matrix of ``Y_{it}`` by the ``N x L`` matrix ``W``, 
we get a ``T x L`` matrix with element ``(t, l)`` consisting of ``\\sum_i y_{it} w_{i,l}``.

Then, we can multiply this by ``H(Θ)`` to get moments (instead of doing ``H(Θ) y_i`` first). 
"""
function m_theta_bar_fast(
  theta::AbstractVector,
  p::Int64, # Number of factors
  ymat_W::AbstractMatrix, # never-treated only
)
  T = size(ymat_W, 1)
  # H(Θ)
  if p > 0
    # Note that the I just takes the first T-p rows, to simplify mat mult
    qld_ymat_W =
      ymat_W[1:(T-p), :] + reshape(theta, T - p, p) * ymat_W[(T-p+1):T, :]
  else
    qld_ymat_W = ymat_W
  end

  return qld_ymat_W
end

# This function is written more idiomatically and is mainly used to test the fast implementation
function m_theta_bar_idiomatic(
  theta::AbstractVector,
  p::Int64,
  ymat::AbstractMatrix,
  W::AbstractMatrix,
)
  # Since I transposed W for quicker indexing
  N_units, N_instruments = size(W)
  T = size(ymat, 1)

  # H(Θ)
  if p > 0
    HΘ = [I(T - p) reshape(theta, T - p, p)]
  else
    HΘ = Matrix(I(T))
  end

  m_bar_correct = zeros(N_instruments * (T - p))
  for i in axes(ymat, 2)
    m_bar_correct += kron(Hθ * ymat[:, i], W[i, :])
  end

  return (1 / N_units) * m_bar_correct
end

function obj_theta(
  theta::AbstractVector,
  weight::Matrix{Float64},
  p::Int64, # Number of factors
  ymat_W::AbstractMatrix, # Note this is already scaled
)
  m_bar = vec(transpose(m_theta_bar_fast(theta, p, ymat_W)))
  # Multiply by `N_units` since `m_bar` divides by `1/N_units`
  return (m_bar' * weight * m_bar)
end


# Matrix of m(z_i, θ)
function ms_theta(
  theta::AbstractVector,
  p::Int64,
  ymat::AbstractMatrix,
  W::AbstractMatrix,
)
  # Since I transposed W for quicker indexing
  N_units, N_instruments = size(W)
  T = size(ymat, 1)

  # H(Θ)
  if p > 0
    Hθ = [I(T - p) reshape(theta, T - p, p)]
  else
    Hθ = Matrix(I(T))
  end

  ms = zeros(N_instruments * (T - p), N_units)
  for i in axes(ymat, 2)
    ms[:, i] .= kron(Hθ * ymat[:, i], W[i, :])
  end
  return ms
end



# QLD GMM Procedure ------------------------------------------------------------
function gmm_qld(ymat, W, p)
  """
      gmm_qld(ymat, W, p)

    Estimates the quasi-long differencing estimator 

    `ymat` is a T x N matrix with element y_{it}
    `W` is a N x K matrix of instruments for each unit
    `p` is the number of factors to estimate. `p` must be <= K.

  """
  T, N_units = size(ymat)
  N_instruments = size(W, 2)
  k = 0

  # each column consists of \sum_i y_{it} ⊗ w_i for each t
  # multiplied by 1/N to keep from having to do this every time `m_theta_bar` is estimated
  ymat_W = (ymat * W) / N_units

  ## First stage ----
  theta_initial = ones((T - p) * p)
  weight = Matrix(1.0 * I((T - p) * N_instruments))

  # mbar = m_theta_bar_fast(theta_initial, p, ymat_W)
  # @show mbar[:, 1]
  # Mbar_initial = ForwardDiff.jacobian(
  #   x -> vec(transpose(m_theta_bar_fast(x, p, ymat_W))), theta_initial
  # )
  # @show Mbar_initial[:, 1]

  if p > 0
    optres = Optim.optimize(
      x -> obj_theta(x, weight, p, ymat_W),
      theta_initial,
      BFGS(; linesearch=LineSearches.BackTracking()),
      Optim.Options(; g_abstol=eps(Float64));
      autodiff=:forward,
    )
    theta_hat = Optim.minimizer(optres)
  else
    theta_hat = theta_initial
  end

  ## Estimate optimal weight matrix using first-stage estimates ----
  ms = ms_theta(theta_hat, p, ymat, W)
  Wopt = pinv(ms * ms' / N_units)

  ## Second stage ----
  if p > 0
    # Having problems with numeric stability when doing all the 1/N_inf correctly
    RESCALE = 15.266221542835233 * 83.84513366993183
    optres_opt = Optim.optimize(
      x -> 1 / RESCALE * obj_theta(x, Wopt, p, ymat_W),
      theta_hat,
      BFGS(; linesearch=LineSearches.BackTracking()),
      Optim.Options(; iterations=10000);
      autodiff=:forward,
    )

    theta_hat_opt = Optim.minimizer(optres_opt)
    J = N_units * Optim.minimum(optres_opt) * RESCALE
  else
    theta_hat_opt = theta_hat
    J = N_units * obj_theta(theta_hat_opt, Wopt, p, ymat_W)
  end

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

  Mbar_theta = ForwardDiff.jacobian(
    x -> vec(transpose(m_theta_bar_fast(x, p, ymat_W))), theta_hat_opt
  )

  return theta_hat_opt, Wopt, Mbar_theta, J, p_value_hansen_sargent
end


