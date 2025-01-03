# Estimate Factor Model Imputation TE function ---------------------------------
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

# QLD GMM Functions ------------------------------------------------------------
"""
    fast_m_theta_bar(parameters, p, ymat_W_inf)

These are the quasi-long differencing moments
1/N \\sum_i w_i ⊗ [H(Θ)' * y_i]

Note you must call `vec(transpose(fast_m_theta_bar(...)))` to get the correct moment vector.

# Details
This uses some linear algebra tricks since we only care about the average of the unit moments.
By multiplying the ``T x N_∞`` matrix of ``Y_{it}`` by the ``N_∞ x L`` matrix ``W``, 
we get a ``T x L`` matrix with element ``(t, l)`` consisting of ``\\sum_i y_{it} w_{i,l}``.

Then, we can multiply this by ``H(Θ)`` to get moments (instead of doing ``H(Θ) y_i`` first). 
"""
function fast_m_theta_bar(
  parameters::AbstractVector,
  p::Int64, # Number of factors
  ymat_W_inf::AbstractMatrix, # never-treated only
)
  T = size(ymat_W_inf, 1)
  # H(Θ)
  if p > 0
    # Note that the I just takes the first T-p rows, to simplify mat mult
    theta = parameters[1:((T - p) * p)]
    qld_ymat_W_inf =
      ymat_W_inf[1:(T - p), :] + reshape(theta, T - p, p) * ymat_W_inf[(T - p + 1):T, :]
  else
    qld_ymat_W_inf = ymat_W_inf
  end

  return qld_ymat_W_inf
end

# This function is written more idiomatically and is mainly used to test the fast implementation
function m_theta_bar_idiomatic(
  parameters::AbstractVector,
  p::Int64,
  N_units::Int64,
  ymat_inf::AbstractMatrix,
  W_inf::AbstractMatrix,
)
  # Since I transposed W for quicker indexing
  N_inf, N_instruments = size(W_inf)

  # H(Θ)
  if p > 0
    theta = parameters[1:((T - p) * p)]
    HΘ = [I(T - p) reshape(theta, T - p, p)]
  else
    HΘ = Matrix(I(T))
  end

  P_D_inf = (N_inf / N_units)
  m_bar_correct = zeros(N_instruments * (T - p))
  for i in axes(ymat_inf, 2)
    m_bar_correct += P_D_inf * kron(Hθ * ymat_inf[:, i], W_inf[i, :])
  end

  return (1 / N_units) * m_bar_correct
end

# # The idea of this version is to store `ymat_W_inf[1:(T - p), :]` once and 
# # then just do Θ * ymat_W_inf_end calculation. It is faster, but 
# # doesn't seem to matter for jacobian which dominates run time.
# # I think it is not worth the extra complexity
# 
# function fast_m_theta_bar_v2(
#   parameters::AbstractVector,
#   p::Int64, # Number of factors
#   T::Int64,
#   ymat_W_inf_end::AbstractMatrix, # never-treated only 
# )
#   if p == 0
#     error("This function only works for p > 0.")
#   end
# 
#   # H(Θ)
#   # Hθ = reshape(parameters[1:((T - p) * p)], T - p, p)
#   # The result is the moments, but needs to be reshaped by `vec(m_bar_mat')`
#   return reshape(parameters[1:((T - p) * p)], T - p, p) * ymat_W_inf_end
# end

"""
    cov_m_theta(parameters, p, N_units, ymat_inf, W_inf)

  This function will only be called once to estimate the weight matrix for GMM, 
  so does not need to be optimized

  It estimates ``1/N \\sum_i m_i * m_i'``. The weights matrix to GMM will be `pinv(cov_m_theta(...))`.
"""
function cov_m_theta(
  parameters::AbstractVector,
  p::Int64,
  N_units::Int64,
  ymat_inf::AbstractMatrix,
  W_inf::AbstractMatrix,
)
  # Since I transposed W for quicker indexing
  N_inf, N_instruments = size(W_inf)

  # H(Θ)
  if p > 0
    theta = parameters[1:((T - p) * p)]
    HΘ = [I(T - p) reshape(theta, T - p, p)]
  else
    HΘ = Matrix(I(T))
  end

  ms = zeros(N_inf, N_instruments * (T - p))
  for i in axes(ymat_inf, 2)
    ms[i, :] .= (N_inf / N_units) * kron(Hθ * ymat_inf[:, i], W_inf[i, :])
  end
  return (ms' * ms) / N_units
end

function ms_theta(
  parameters::AbstractVector,
  p::Int64,
  N_units::Int64,
  ymat_inf::AbstractMatrix,
  W_inf::AbstractMatrix,
)
  # Since I transposed W for quicker indexing
  N_inf, N_instruments = size(W_inf)

  # H(Θ)
  if p > 0
    theta = parameters[1:((T - p) * p)]
    Hθ = [I(T - p) reshape(theta, T - p, p)]
  else
    Hθ = Matrix(I(T))
  end

  ms = zeros(N_inf, N_instruments * (T - p))
  for i in axes(ymat_inf, 2)
    ms[i, :] .= kron(Hθ * ymat_inf[:, i], W_inf[i, :])
  end
  ms *= (N_inf / N_units) # P(D_∞ = 1)
  return ms
end

function obj_theta(
  parameters::AbstractVector,
  weight::Matrix{Float64},
  p::Int64, # Number of factors
  N_units::Int64,
  ymat_W_inf::AbstractMatrix, # Note this is already scaled
)
  m_bar = vec(transpose(fast_m_theta_bar(parameters, p, ymat_W_inf)))
  # Multiply by `N_units` since `m_bar` divides by `1/N_units`
  return N_units * (m_bar' * weight * m_bar)
end

# Factor Imputation Event-study Functions --------------------------------------
function m_bar_tau_fix_theta(
  parameters::AbstractVector,
  theta_opt::AbstractVector,
  ymat::AbstractMatrix,
  g::AbstractVector,
  rel_year_shift::AbstractVector,
)
  N_units = size(ymat, 2)
  theta = theta_opt
  tau = parameters[(((T - p) * p) + 1):end]
  Fhat = [reshape(theta, T - p, p); -1.0 * I(p)]

  # Impute y0_hat
  y_diff = zeros(typeof(parameters[1]), size(ymat))
  g_tau = zeros(typeof(parameters[1]), (N_units, length(tau)))
  for i in axes(ymat, 2)
    if g[i] == Inf
      # For now, no need to compute never-treated y(0)
      # y0_hat[:, i] .= Fhat * ((Fhat' * Fhat) \ (Fhat' * ymat[:, i]))
    else
      T0_i = convert(Int64, g[i])
      Fhat_pre = Fhat[1:(T0_i), :]
      y_diff[:, i] .=
        ymat[:, i] - Fhat * ((Fhat_pre' * Fhat_pre) \ (Fhat_pre' * ymat[1:(T0_i), i]))

      rel_year_shift_i = convert(Vector{Int64}, rel_year_shift[(1 + (i - 1) * T):(i * T)])
      g_tau[i, rel_year_shift_i] .= y_diff[:, i] - tau[rel_year_shift_i]
    end
  end

  return vec(mean(g_tau; dims=1))
end

function obj_tau(
  parameters::AbstractVector,
  theta_opt::AbstractVector,
  ymat::AbstractMatrix,
  g::AbstractVector,
  rel_year_shift::AbstractVector,
)
  N_units = size(ymat, 2)
  m_bar = m_bar_tau_fix_theta(parameters, theta_opt, ymat, g, rel_year_shift)
  return N_units * (m_bar' * m_bar)
end

function m_bar_tau(
  parameters::AbstractVector,
  ymat::AbstractMatrix,
  g::AbstractVector,
  rel_year_shift::AbstractVector,
)
  T, N_units = size(ymat)
  theta = parameters[1:((T - p) * p)]
  tau = parameters[(((T - p) * p) + 1):end]
  Fhat = [reshape(theta, T - p, p); -1.0 * I(p)]

  # Impute y0_hat
  y_diff = zeros(typeof(parameters[1]), size(ymat))
  g_tau = zeros(typeof(parameters[1]), (N_units, length(tau)))
  for i in axes(ymat, 2)
    if g[i] == Inf
      # For now, no need to compute never-treated y(0)
      # y0_hat[:, i] .= Fhat * ((Fhat' * Fhat) \ (Fhat' * ymat[:, i]))
    else
      T0_i = convert(Int64, g[i])
      Fhat_pre = Fhat[1:(T0_i), :]
      y_diff[:, i] .=
        ymat[:, i] - Fhat * ((Fhat_pre' * Fhat_pre) \ (Fhat_pre' * ymat[1:(T0_i), i]))

      rel_year_shift_i = convert(Vector{Int64}, rel_year_shift[(1 + (i - 1) * T):(i * T)])
      g_tau[i, rel_year_shift_i] .= y_diff[:, i] - tau[rel_year_shift_i]
    end
  end

  return vec(mean(g_tau; dims=1))
end

function ms_tau(
  parameters::AbstractVector,
  ymat::AbstractMatrix,
  g::AbstractVector,
  rel_year_shift::AbstractVector,
)
  N_units = size(ymat, 2)
  theta = parameters[1:((T - p) * p)]
  tau = parameters[(((T - p) * p) + 1):end]
  Fhat = [reshape(theta, T - p, p); -1.0 * I(p)]

  # Impute y0_hat
  y_diff = zeros(typeof(parameters[1]), size(ymat))
  g_tau = zeros(typeof(parameters[1]), (N_units, length(tau)))
  for i in axes(ymat, 2)
    if g[i] == Inf
      # For now, no need to compute never-treated y(0)
      # y0_hat[:, i] .= Fhat * ((Fhat' * Fhat) \ (Fhat' * ymat[:, i]))
    else
      T0_i = convert(Int64, g[i])
      Fhat_pre = Fhat[1:(T0_i), :]
      y_diff[:, i] .=
        ymat[:, i] - Fhat * ((Fhat_pre' * Fhat_pre) \ (Fhat_pre' * ymat[1:(T0_i), i]))

      rel_year_shift_i = convert(Vector{Int64}, rel_year_shift[(1 + (i - 1) * T):(i * T)])
      g_tau[i, rel_year_shift_i] .= y_diff[:, i] - tau[rel_year_shift_i]
    end
  end

  return g_tau
end








#
#
#
#
# %% Function to estimate a treatment effect
function est_te(df, p, outcome_var; export_y0=false, outcome="", quick_plot=false)
  # Sort just in case
  sort!(df, [:fips, :year])
  ytilde = within_transform(df, outcome_var, :g, :year, T0)
  df.ytilde = ytilde

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
  ytilde = df.ytilde
  fips = df.fips
  year = df.year
  rel_years = sort(unique(df[(df.g .!= Inf), :rel_year]))

  # Find untreated units
  N_units = length(unique(fips))
  idx_control = findall(df.g[year .== 1977] .== Inf)
  N_inf = length(idx_control)

  # Number of instruments
  N_instruments = length(instruments)
  W = Matrix(df[year .== 1977, instruments])
  W = Matrix(transpose(W))

  # Number of covariates
  k = 0

  # Estimate theta -------------------------------------------------------------
  # Step 1. Initial estimate
  parameters = [ones((T - p) * p); zeros(length(rel_years))]
  weight = Matrix(1.0 * I((T - p) * N_instruments))

  # EDIT ----
  # mean(m_thetas(parameters,  p, ytilde, idx_control, W, N_units, T), dims = 1)
  # m_theta_bar(parameters, p, ytilde, idx_control, W, N_units, T)
  # obj_theta(
  #   parameters, p, weight, ytilde, idx_control, W, N_units, T,
  # )
  # END EDIT ----

  if p > 0
    optres = optimize(
      x -> obj_theta(x, p, weight, ytilde, idx_control, W, N_units, T),
      parameters,
      BFGS(; linesearch=BackTracking()),
      Optim.Options(; g_abstol=eps(Float64));
      autodiff=:forward,
    )

    parameters_hat = Optim.minimizer(optres)
  else
    parameters_hat = parameters
  end

  # Step 2. Optimal weighting matrix 
  ms = m_thetas(parameters_hat, p, ytilde, idx_control, W, N_units, T)
  return (ms' * ms) / N_units
  Ωinv_theta = pinv((ms' * ms) / N_units)

  if p > 0
    # Having problems with numeric stability when doing all the 1/N_inf correctly
    RESCALE = 15.266221542835233 * 83.84513366993183
    optres_opt = optimize(
      x -> 1 / RESCALE * obj_theta(x, p, Ωinv_theta, ytilde, idx_control, W, N_units, T),
      parameters_hat,
      BFGS(; linesearch=BackTracking()),
      Optim.Options(; iterations=10000);
      autodiff=:forward,
    )

    parameters_hat_opt = Optim.minimizer(optres_opt)
    theta_hat_opt = parameters_hat_opt[1:((T - p) * p)]
    J = Optim.minimum(optres_opt) * RESCALE
  else
    J = obj_theta(parameters_hat, p, Ωinv_theta, ytilde, idx_control, W, N_units, T)
  end

  # Testing for $p$ following Ahn, Lee, and Schmidt (2013)
  # Χ^2((T - p_0) (q - p_0) - k)
  # Where T is the number of periods, p_0 is the number of factors under the null, q is the number of instruments, and k is the number of covariates
  dist_hansen_sargent = Chisq((T - p) * (N_instruments - p) - k)
  p_value_hansen_sargent = 1 - cdf(dist_hansen_sargent, J)

  # Impute untreated potential outcomes
  df.y0hat = zeros(size(df, 1))

  if p > 0
    FΘ = [reshape(theta_hat_opt, T - p, p); (-1 * I(p))]
    FΘ = Matrix(FΘ)

    treated_ids = unique(df[df.g .< Inf, :fips])
    for id in treated_ids
      unit = df[(df.fips .== id), :]
      g_id = Int(unit.g[1])

      y_pre = unit[(unit.year .< g_id), :ytilde]
      FΘ_pre = FΘ[1:(g_id - minimum(YEARS)), :]

      df[(df.fips .== id), :y0hat] = FΘ * ((FΘ_pre' * FΘ_pre) \ (FΘ_pre' * y_pre))
    end
  end

  # Estimate τ^ℓ ---------------------------------------------------------------
  df.y_diff = df.ytilde .- df.y0hat

  if export_y0 == true
    @assert outcome != "" "Outcome needs to be specified with `export_y0`"

    outfile = "estimates/est_y0_qld.csv"
    out = DataFrame(CSV.File(outfile))

    ytilde_name = Symbol(string(outcome_var) * "_tilde")
    ytilde_0hat_name = Symbol(string(outcome_var) * "_tilde_0hat_qld_p_$(p)")

    merge = select(df, :fips, :year, :ytilde => ytilde_name, :y0hat => ytilde_0hat_name)
    if hasproperty(out, ytilde_name)
      select!(merge, Not(ytilde_name))
    end
    if hasproperty(out, ytilde_0hat_name)
      select!(out, Not(ytilde_0hat_name))
    end

    out = innerjoin(out, merge; on=[:fips, :year])

    CSV.write(outfile, out)
  end

  df.rel_year_cat = CategoricalArray(df.rel_year)
  lm_tau = lm(@formula(y_diff ~ 0 + rel_year_cat), df[df.g .< Inf, :])

  # Plot point estimates
  if quick_plot == true
    plot(sort(rel_years), coef(lm_tau); marker=:circle)
  end

  return (
    est=[coef(lm_tau), stderror(lm_tau)],
    sargan=(J, p_value_hansen_sargent),
    imputation=[df.ytilde, df.y0hat],
  )
end
