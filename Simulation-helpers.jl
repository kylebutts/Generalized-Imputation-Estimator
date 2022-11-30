
# DGP for Simulations

function generate_data(rng, f, include_factor = true, parallel_trends = false, instrument_noise = 1)

  # f = 1.0:T
  # p = 1


  id = repeat(1:N, inner=T)
  t = repeat(1:T, outer=N)

  # factor loadings and unit FEs
  g = zeros(N * p, 1)
  W = zeros(N * p, 1)
  mu = zeros(N, 1)

  for i in 1:N
    idx = ((i-1)*p+1):(i*p)
    mu[idx, :] = 2 * randn(rng, 1)

    d = MvNormal(diag(mu[idx, :]), I(p))
    g[idx, :] = rand(rng, d, 1)'
    
    # Generate time-invariant covariate
    W[idx, :] = g[idx,:] + sqrt(instrument_noise) * randn(rng, 1)
  end

  # time FEs
  eta = zeros(T, 1)
  eta[1, :] = randn(rng, 1)
  for i in 2:T
    eta[i, :] = 0.75 * eta[i-1, :] + randn(rng, 1)
  end

  # Generate treatment, increasing probability in g2
  g_range = maximum(g[:, 1]) - minimum(g[:, 1])
  g_trans = (g ./ g_range)

  if(parallel_trends == false) 
    # Probability of treatment increasing in g2
    prob = 0.5 .+ 1 .* g_trans

    # Unconditional probability equal to 0.5
    prob = prob .* 0.5 ./ mean(prob)
  else
    prob = fill(0.5, N)
  end
  

  ever_treated = rand(rng, N) .>= prob  
  ever_treated = repeat(vec(ever_treated), inner = T)
  treat = (ever_treated .== 1) .& (t .> T0)


  
  # Correlation matrix for errors
  C = zeros(T, T)
  for t = 1:T
    for s = 1:T
      C[t, s] = 0.75^abs(t - s)
    end
  end
  d = MvNormal(zeros(T), C)
  u = rand(rng, d, N) |>
    x -> reduce(vcat, x)

  # Generating y(0)
  time_fe = kron(fill(1, N), eta)
  unit_fe = kron(mu, fill(1, T))

  if(include_factor == true)
    y0 = unit_fe .+ time_fe .+ kron(I(N), f) * g .+ u
  else 
    y0 = unit_fe .+ time_fe .+ u
  end
  
  y0 = vec(y0)

  # Generate y
  y = y0

  tau = treat .* repeat(true_te', N) .* (1 .+ unit_fe)

  idx_treat_6 = treat .* (t .== 6)
  idx_treat_7 = treat .* (t .== 7)
  idx_treat_8 = treat .* (t .== 8)
  
  tau[idx_treat_6] = true_te[6] .* (1 .+ standardize(ZScoreTransform, tau[idx_treat_6], dims = 1)[:,1])
  tau[idx_treat_7] = true_te[7] .* (1 .+ standardize(ZScoreTransform, tau[idx_treat_7], dims = 1)[:,1])
  tau[idx_treat_8] = true_te[8] .* (1 .+ standardize(ZScoreTransform, tau[idx_treat_8], dims = 1)[:,1])

  # mean(tau[8 .+ idx_0])

  y = y .+ tau


  W = kron(W, fill(1, T))

  data = DataFrame(
    id=id, t=t, treat=treat, ever_treated=ever_treated,
    y0=vec(y0), y=vec(y),
    W=vec(W)
  )

end

# Note: this is not general, it hard codes some things from the simulation
# CCE model, not currently used

function est_ccedid(data)
  
  sort!(data, [:id, :t])

  X_vars = [:X1, :X2]
  X = [data.X1 data.X2]
  y = data.y
  id = data.id
  t = data.t
  ever_treated = data.ever_treated


  # Fhat = cross-sectional averages of X for control group
  Fhat = data |>
    x -> x[(x.ever_treated.==0), :] |>
    x -> groupby(x, :t) |>
    x -> combine(x,
      X_vars .=> (x -> mean(x)) .=> X_vars
    ) |>
    x -> sort(x, :t) |>
    x -> x[:, X_vars] |>
    x -> Matrix(x)

  Fpre = Fhat[1:T0, :]
  Fpost = Fhat[(T0+1):T, :]

  M_Fpre = I(T0) - Fpre * inv(Fpre' * Fpre) * Fpre'

  # Estimate CCEP for beta_hat
  B = zeros(K, K)
  A = zeros(K, 1)
  for i = 1:N
    Xi = Matrix(data[
      (data.id .== i) .& (data.t .<= T0), X_vars
    ])
    yi = data[
      (data.id .== i) .& (data.t .<= T0), :y
    ]
    B = B + Xi' * M_Fpre * Xi
    A = A + Xi' * M_Fpre * yi
  end

  # CCEP estimator using pre-treatment X's for all groups
  bhat = inv(B) * A

  # # imputed factor loadings and covariates
  # # (NOTE: The dimension depends on how many factor proxies. For now it's K because I only use Xbar)
  treated_ids = data |>
    x -> x[x.ever_treated .== 1, :id] |>
    x -> unique(x)
  
  ghat = zeros(Float32, N * K, 1)
  Xhat = zeros(Float32, N * T, K)
  y0hat = zeros(Float32, N * T)
  ccedid = zeros(Float32, T - T0)
  for i in treated_ids
    idx_ghat = ((i-1)*K+1):(i*K)
    idx_i_pre = ((i-1)*T+1):((i-1)*T+T0)
    idx_i_post = ((i-1)*T+T0+1):(i*T)

    # impute factor loadings
    ghat[idx_ghat, :] =
      inv(Fpre' * Fpre) * Fpre' *
      (y[idx_i_pre, :] - X[idx_i_pre, :] * bhat)

    # impute X(0)
    Xhat[idx_i_post, :] =
      Fpost * inv(Fpre' * Fpre) * Fpre' *
      X[idx_i_pre, :]

    # impute y(0)
    y0hat[idx_i_post] =
      Xhat[idx_i_post, :] * bhat +
      Fpost * ghat[((i-1)*K+1):(i*K), :]

    # estimate TE
    ccedid = ccedid + 
      (y[idx_i_post] - y0hat[idx_i_post]) / N1
  end

  return (bhat=bhat', ccedid=ccedid)
end

# TWFE OLS
function est_twfe(data)
  data.rel_year = (data.t .- (T0 + 1)) .* (data.ever_treated .== 1)
  est = reg(
    data,
    @formula(y ~ rel_year + fe(id) + fe(t)),
    contrasts=Dict(:rel_year => DummyCoding(base=-1))
  )

  post = occursin.(r"rel\_year: [0-9]+", coefnames(est))

  return coef(est)[post]
end

# TWFE OLS with w_i β_t

function est_twfe_w_covariates(data)
  data.rel_year = 
    (data.t .- (T0 + 1)) .* (data.ever_treated .== 1) .+
    -1 * (data.ever_treated .== 0)
  est = reg(
    data,
    @formula(y ~ rel_year + fe(t) & W + fe(id) + fe(t)),
    contrasts=Dict(:rel_year => DummyCoding(base=-1)),
    save=:fe
  )

  post = occursin.(r"rel\_year: [0-9]+", coefnames(est))

  return coef(est)[post]
end

# {did2s} from R

function est_did2s_R(data)
  data.rel_year = 
    (data.t .- (T0 + 1)) .* (data.ever_treated .== 1) .+
    -1 * (data.ever_treated .== 0)

  rcopy(R"
    coef = did2s::did2s(
      data = $(data), 
      yname = 'y', 
      treatment = 'treat', 
      cluster_var = 'id', 
      first_stage = ~ 0 | id + t, 
      second_stage = ~ i(rel_year, ref = c(-1)),
      verbose = FALSE
    ) |> 
      coef()

    as.numeric(coef[
      stringr::str_detect(names(coef), 'rel_year::[0-9]+')
    ])
  ")
end

# {did2s} from R adding in w_i β_t

function est_did2s_covariates_R(data)
  data.rel_year = 
    (data.t .- (T0 + 1)) .* (data.ever_treated .== 1) .+
    -1 * (data.ever_treated .== 0)

  rcopy(R"
    coef = did2s::did2s(
      data = $(data), 
      yname = 'y', 
      treatment = 'treat', 
      cluster_var = 'id', 
      first_stage = ~ i(t, W) | id + t, 
      second_stage = ~ i(rel_year, ref = c(-1)),
      verbose = FALSE
    ) |> 
      coef()

    as.numeric(coef[
      stringr::str_detect(names(coef), 'rel_year::[0-9]+')
    ])
  ")
end


# Take matrix and make into latex table

function matrix_to_table(mat::Matrix{String}) 
  (r, c) = size(mat)
  tex = ""

  for i in 1:r
    for j in 1:c
      sep = (j == 1) ? "" : " & "
      tex *= "$(sep)$(mat[i, j])"
    end
    tex *= " \\\\ \n"
  end

  return tex
end


# Generate AR(1) process

sample_start(c::Real, ϕ::Real, σ::Real, rng::AbstractRNG) =
  c / (1 - ϕ) + randn(rng) * σ / sqrt(1 - ϕ^2)

function ar1!(v::AbstractVector{<:AbstractFloat}, c::Real, ϕ::Real, σ::Real;
              rng::AbstractRNG=Random.default_rng(),
              start::Real=sample_start(c, ϕ, σ, rng))
  σ < 0 && throw(ArgumentError("σ is not allowed to be negative"))
  abs(ϕ) < 1 || throw(ArgumentError("|ϕ| < 1 is required"))
  x = start
  n = length(v)
  if n > 0
      v[1] = x
      for i in 2:n
          x = c + x * ϕ + randn(rng) * σ
          v[i] = x
      end
  end
  return v
end

function ar1(n::Integer, c::Real, ϕ::Real, σ::Real;
             rng::AbstractRNG=Random.default_rng(),
             start::Real=sample_start(c, ϕ, σ, rng))
  return ar1!(Vector{Float64}(undef, n), c, ϕ, σ, rng=rng, start=start)
end
