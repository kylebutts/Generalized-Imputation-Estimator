function qld_imputation(
  df;
  y::Union{String,Symbol},
  id::Union{String,Symbol},
  t::Union{String,Symbol},
  g::Union{String,Symbol},
  W::Union{String,Symbol,Vector{String},Vector{Symbol}},
  do_within_transform::Bool,
  p::Union{Int64,Real},
  type::String="dynamic",
  return_y0::Bool=false,
  return_naive_se::Bool=false,
)
  #
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

  N = length(unique(id))
  idx_control = findall(g_shift .== Inf)
  N_inf = length(idx_control)

  N_tau_gt = zeros(T * length(uniq_g_shift))
  for (l, curr_g) in enumerate(uniq_g_shift)
    curr_g = convert(Int64, curr_g)
    curr_idx = findall(g_shift .== curr_g) # All units with this g
    for i in curr_idx
      gt_idx = (1 + ((l - 1) * T)):(l * T)
      N_tau_gt[gt_idx] .+= 1
    end
  end

  # N x L matrix of instruments
  W = df[t .== min_t, W]
  if (ndims(W) == 1)
    W = reshape(W, length(W), 1)
  else
    W = Matrix(W)
  end
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
  ymat = reshape(y, T, N)
  if do_within_transform == true
    ymat = within_transform(ymat, idx_control, N_pre)
  end
  T, N = size(ymat)

  # Estimate Quasi-long differencing estimator
  # ----
  # Two-step GMM Estimation of QLD parameters
  # Using only `idx_control` to estimate the factors
  if p >= 0
    theta_hat_opt, W_opt, Mbar_theta, J, p_value_hansen_sargent = gmm_qld_p_known(
      p, # Number of factors
      ymat,
      W,
      idx_control,
    )
  elseif p == -1
    p = Int(0)
    # @info "Selecting p based on Hansen-Sargent statistic"
    while p <= N_pre - 1
      # @info "Trying p=$(p)"
      theta_hat_opt, W_opt, Mbar_theta, J, p_value_hansen_sargent = gmm_qld_p_known(
        p, # Number of factors
        ymat,
        W,
        idx_control,
      )

      # Note that if p == N_instruments, p_value will be returned as 1
      if p_value_hansen_sargent >= 0.10
        break
      end
      p += 1
    end
    # @info "Selected p=$(p) based on Hansen-Sargent statistic"
  end


  # Estimate τ(g,t) parameters
  # ----
  tau_gt_hat, N_tau_gt = estimate_tau_gt(theta_hat_opt, p, ymat, g_shift)

  # Estimate VCOV of τ(g,t) 
  # ----
  ms = ms_theta(theta_hat_opt, p, ymat, W, idx_control)
  ms *= 1 / (N_inf / N)
  gs = ms_tau_gt(theta_hat_opt, tau_gt_hat, p, ymat, g_shift)

  Gbar_theta = ForwardDiff.jacobian(
    x -> mean(ms_tau_gt(x, tau_gt_hat, p, ymat, g_shift); dims=1), theta_hat_opt
  )
  # Gbar_tau = -1 * I(length(tau_gt_hat))
  # Gbar_tau = ForwardDiff.jacobian(
  #   x -> mean(ms_tau_gt(theta_hat_opt, x, p, ymat, g_shift); dims=1), tau_gt_hat
  # )

  IF_tau = (1 / sqrt(N) * gs')
  IF_theta =
    Gbar_theta *
    pinv(Mbar_theta' * W_opt * Mbar_theta) *
    Mbar_theta' *
    W_opt *
    (1 / sqrt(N) * ms)

  vcov_tau_gt_naive = 1 / N * (IF_tau * IF_tau')
  vcov_tau_gt = 1 / N * (IF_theta + IF_tau) * (IF_theta + IF_tau)'


  # Aggregate effects if needed
  if type == "gt"
    ret = (estimate=tau_gt_hat, vcov=vcov_tau_gt, selected_p=p)
    if return_naive_se
      ret = merge(ret, (vcov_naive=vcov_tau_gt_naive,))
    end
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

    ret = (rel_year=uniq_rel_years, estimate=tau_es_hat, vcov=vcov_tau_es, selected_p=p)
    if return_naive_se
      ret = merge(ret, (vcov_naive=vcov_tau_es_naive,))
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
    vcov_tau_overall_naive = mat_agg_overall * vcov_tau_gt_naive * mat_agg_overall'

    ret = (estimate=tau_overall_hat, vcov=vcov_tau_overall, selected_p=p)
    if return_naive_se
      ret = merge(ret, (vcov_naive=vcov_tau_overall_naive,))
    end
  end

  if return_y0 == true
    impute_df = df[:, [id_name, t_name, g_name]]
    if do_within_transform == true
      impute_df.ytilde0_hat = vec(impute_y0(theta_hat_opt, p, ymat, g_shift))
      impute_df.ytilde = vec(ymat)
    else
      impute_df.y0_hat = vec(impute_y0(theta_hat_opt, p, ymat, g_shift))
      impute_df.y = vec(ymat)
    end
    ret = merge(ret, (impute_df=impute_df,))
  end

  return ret
end
