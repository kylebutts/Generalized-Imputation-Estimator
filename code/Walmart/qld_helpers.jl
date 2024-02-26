# Estimate Factor Model Imputation TE function ---------------------------------
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
function est_te(df, p, outcome_var; export_y0=false, outcome="", quick_plot=false)
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
  n_units = length(unique(fips))
  idx_control = findall(df.g[year .== 1977] .== Inf)
  N_inf = length(idx_control)

  # Number of instruments
  n_instruments = length(instruments)
  W = Matrix(df[year .== 1977, instruments])
  W = Matrix(transpose(W))

  # Number of covariates
  k = 0

  # Estimate theta -------------------------------------------------------------
  # Step 1. Initial estimate
  parms = [ones((T - p) * p); zeros(length(rel_years))]
  weight = Matrix(1.0 * I((T - p) * n_instruments))

  # EDIT ----
  # mean(m_thetas(parms,  p, ytilde, idx_control, W, n_units, T), dims = 1)
  # m_theta_bar(parms, p, ytilde, idx_control, W, n_units, T)
  # obj_theta(
  #   parms, p, weight, ytilde, idx_control, W, n_units, T,
  # )
  # END EDIT ----

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
    # Having problems with numeric stability when doing all the 1/N_inf correctly
    optres_opt = optimize(
      x ->
        1 / 15.266221542835233 *
        obj_theta(x, p, Ωinv_theta / 83.84513366993183, ytilde, idx_control, W, n_units, T),
      parms_hat,
      BFGS(; linesearch=BackTracking()),
      Optim.Options(; iterations=10000);
      autodiff=:forward,
    )

    parms_hat_opt = Optim.minimizer(optres_opt)
    theta_hat_opt = parms_hat_opt[1:((T - p) * p)]
    J = Optim.minimum(optres_opt) * 15.266221542835233 * 83.84513366993183
  else
    J = obj_theta(parms_hat, p, Ωinv_theta, ytilde, idx_control, W, n_units, T)
  end

  # Testing for $p$ following Ahn, Lee, and Schmidt (2013)
  # Χ^2((T - p_0) (q - p_0) - k)
  # Where T is the number of periods, p_0 is the number of factors under the null, q is the number of instruments, and k is the number of covariates
  dist_hansen_sargent = Chisq((T - p) * (n_instruments - p) - k)
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
  df.tau_hat = df.ytilde .- df.y0hat
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

  # Factor Estimate
  df.rel_year_cat = CategoricalArray(df.rel_year)
  lm_tau = lm(@formula(tau_hat ~ 0 + rel_year_cat), df[df.g .< Inf, :])

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
