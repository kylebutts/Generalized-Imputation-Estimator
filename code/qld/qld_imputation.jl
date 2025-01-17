function qld_imputation(
  df;
  y::Union{String,Symbol},
  id::Union{String,Symbol},
  t::Union{String,Symbol},
  g::Union{String,Symbol},
  W::Union{Vector{String},Vector{Symbol}},
  do_within_transform::Bool,
  p::Union{Int64,Real},
  type::String="dynamic",
  return_y0::Bool,
  return_naive_se::Bool,
)
  # y = :log_retail_emp
  # id = :fips
  # t = :year
  # g = :g
  DataFrames.sort!(df, [g, id, t])

  # Check if the panel is balanced (approximately)
  panel_counts = DataFrames.combine(DataFrames.groupby(df, [id]), DataFrames.nrow => :count)
  @assert length(unique(panel_counts.count)) == 1 "Panel is not balanced"
  id_name = Symbol(id)
  t_name = Symbol(t)
  g_name = Symbol(g)
  y = df[!, y]
  id = df[!, id]
  t = df[!, t]
  g = df[!, g]

  # 
  uniq_t = unique(t)
  uniq_g = unique(g)
  T = length(uniq_t)
  T_0 = minimum(uniq_g) - 1
  N_pre = count(x -> x .<= T_0, uniq_t)

  min_t = minimum(t)
  g_shift = (g .- min_t)[t .== min_t]
  uniq_g_shift = unique(g_shift[g_shift .!== Inf])
  rel_year = t - g
  uniq_rel_years = unique(rel_year[rel_year .!== -Inf])
  uniq_rel_years = sort(uniq_rel_years)

  N_units = length(unique(id))
  idx_control = findall(g_shift .== Inf)

  # N x L matrix of instruments
  W = Matrix(df[t .== min_t, W])
  N_instruments = size(W, 2)

  # Number of strictly exogenous covariates
  # k = 0

  p = convert(Int64, p)
  @assert p >= -1 "`p` must be an integer >= 0. If you want to select p based on the data, use `p = -1`"
  @assert p <= N_pre - 1 "`p` must be smaller than the number of time periods before *any* unit is treated"
  @assert p <= N_instruments "The number of instruments must be >= `p`"
  @assert type ∈ ["gt", "dynamic", "overall"] "Three options: \"gt\", \"dynamic\" and \"overall\""
  @assert minimum(g_shift) > p

  # T x N matrix of $y_{it}$
  ymat = reshape(y, T, N_units)
  if do_within_transform == true
    ymat = within_transform(ymat, idx_control, N_pre)
  end

  # ----
  # Two-step GMM Estimation of QLD parameters
  # Using only `idx_control` to estimate the factors
  if p >= 0
    theta_hat_opt, Wopt, Mbar_theta, J, p_value_hansen_sargent = gmm_qld(
      ymat[:, idx_control], W[idx_control, :], p
    )
    # TODO:
    # Currently, this multiplies by N_inf a bunch of places instead of N_units, so need to be careful to multiply moments the correct way
    # In short, Need a bunch more mult of N and N_∞
  elseif p == -1
    p = 0
    # @info "Selecting p based on Hansen-Sargent statistic"
    while p <= N_pre - 1
      theta_hat_opt, Wopt, Mbar_theta, J, p_value_hansen_sargent = gmm_qld(
        ymat[:, idx_control], W[idx_control, :], p
      )

      # Note that if p == N_instruments, p_value will be returned as 1
      if p_value_hansen_sargent >= 0.10
        break
      end
      p += 1
    end
    # @info "Selected p=$(p) based on Hansen-Sargent statistic"
  end

  # ----
  # Estimate \tau(g,t) using imputation
  tau_gt_hat, N_tau_gt = estimate_tau_gt(theta_hat_opt, p, ymat, g_shift)

  # ----
  # Asymptotic variance estimation
  #
  # Following Newey and McFadden (1994), the influence function of a two-step 
  # estimator with a method of moments second-step is given in (6.6).
  # 
  # IF(m(z, γ)) = influence-function of first-step GMM parameters
  # g(z, τ, γ) = second-step moments (ms_tau)
  # G_γ(z, τ, γ) = ∇_γ gbar(z, τ, γ) (jacobian of mbar_tau)
  # G_τ(z, τ, γ) = ∇_τ gbar(z, τ, γ) (jacobian of mbar_tau)
  #
  # Two-step influence function is given by:
  # G_τ^{-1} [g(z) - G_γ IF(m(z, γ))]
  # 
  # `Mbar_theta` is returned by `gmm_qld`
  m_theta = ms_theta(theta_hat_opt, p, ymat, W)
  g_tau = ms_tau_gt(theta_hat_opt, p, tau_gt_hat, ymat, g_shift)

  # Gbar_tau = `-1.0 * I(length(tau_gt_hat))`, so I can safely ignore (crossprod will cancel out the negatives)
  # Gbar_theta = ForwardDiff.jacobian(
  #   x -> m_tau_gt_bar(x, p, tau_gt_hat, ymat, g_shift), theta_hat_opt
  # )
  Gbar_theta = ForwardDiff.jacobian(
    # TODO: use m_tau_gt_bar
    x -> mean(ms_tau_gt(x, p, tau_gt_hat, ymat, g_shift); dims=1),
    theta_hat_opt,
  )
  IF_theta =
    Gbar_theta * pinv(Mbar_theta' * Wopt * Mbar_theta) * Mbar_theta' * Wopt * m_theta

  # IF_tau = 1 / N_units * g_tau

  vcov_tau_gt = g_tau' * g_tau + (IF_theta * IF_theta')
  vcov_tau_gt /= N_units^2

  # Naive ses
  vcov_tau_gt_naive = g_tau' * g_tau
  vcov_tau_gt_naive /= N_units^2

  # vcov_tau_gt = (IF_tau' * IF_tau)

  # Aggregate effects if needed
  if type == "gt"
    return tau_gt_hat, vcov_tau_gt
  elseif type == "dynamic"
    # aggte to dynamic ATT (event-study)
    mat_agg_es = zeros(length(uniq_rel_years), length(tau_gt_hat))
    i = 1
    for curr_g in uniq_g_shift
      for t in 1:T
        es_idx = convert(Int, -1 * minimum(uniq_rel_years) + (t - curr_g))
        mat_agg_es[es_idx, i] = N_tau_gt[i]
        i += 1
      end
    end
    # Normalize each row by the row's sum
    mat_agg_es = mat_agg_es ./ sum(mat_agg_es; dims=2)
    tau_es_hat = mat_agg_es * tau_gt_hat
    vcov_tau_es = mat_agg_es * vcov_tau_gt * mat_agg_es'

    vcov_tau_es_naive = mat_agg_es * vcov_tau_gt_naive * mat_agg_es'

    if return_y0 == true

      impute_df = df[:, [id_name, t_name, g_name]]
      if do_within_transform == true
        impute_df.ytilde0_hat = vec(impute_y0(theta_hat_opt, p, ymat, g_shift))
        impute_df.ytilde = vec(ymat)
      else
        impute_df.y0_hat = vec(impute_y0(theta_hat_opt, p, ymat, g_shift))
        impute_df.y = vec(ymat)
      end

      if return_naive_se == true
        return uniq_rel_years, tau_es_hat, vcov_tau_es, impute_df, vcov_tau_es_naive
      else
        return uniq_rel_years, tau_es_hat, vcov_tau_es, impute_df
      end
    else
      if return_naive_se == true
        return uniq_rel_years, tau_es_hat, vcov_tau_es, vcov_tau_es_naive
      else
        return uniq_rel_years, tau_es_hat, vcov_tau_es
      end
    end
  elseif type == "overall"
    # aggte to overall ATT
    mat_agg_overall = zeros(1, length(tau_gt_hat))
    i = 1
    for curr_g in uniq_g_shift
      for t in 1:T
        if t >= curr_g
          mat_agg_overall[1, i] = N_tau_gt[i]
        end
        i += 1
      end
    end
    # Normalize each row by the row's sum
    mat_agg_overall = mat_agg_overall ./ sum(mat_agg_overall; dims=2)
    tau_overall_hat = mat_agg_overall * tau_gt_hat
    vcov_tau_overall = mat_agg_overall * vcov_tau_gt * mat_agg_overall'

    return tau_overall_hat, vcov_tau_overall
  end

end
