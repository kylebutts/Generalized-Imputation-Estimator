"""
    qld_imputation(df; y, id, t, g, W, do_within_transform, p, type="dynamic", return_y0=false, vcov_type="pointwise")

Estimate treatment effects using Quasi-Long Differencing (QLD) imputation method for panel data with staggered adoption.

# Arguments
- `df`: A DataFrame containing the panel data.
- `y::Union{String,Symbol}`: The outcome variable name.
- `id::Union{String,Symbol}`: The unit identifier variable name.
- `t::Union{String,Symbol}`: The time period variable name.
- `g::Union{String,Symbol}`: The treatment group variable name (timing of treatment). Units never treated should have `g` set to `Inf`.
- `W::Union{String,Symbol,Vector{String},Vector{Symbol}}`: Variable(s) to use as instruments.
- `do_within_transform::Bool`: Whether to apply within-unit transformation to the outcome variable.
- `p::Union{Int64,Real}`: Number of factors to use in the model. If `p = -1`, the number of factors is selected based on Hansen-Sargent statistic.
- `type::String="dynamic"`: The type of treatment effect to estimate:
  - `"gt"`: Group-time specific treatment effects.
  - `"dynamic"`: Event study effects relative to treatment timing.
  - `"overall"`: Overall average treatment effect.
- `return_y0::Bool=false`: Whether to return the imputed counterfactual outcomes.
- `vcov_type::String="pointwise"`: The type of variance-covariance matrix to compute:
  - `"pointwise"`: Standard pointwise asymptotic inference.
  - `"uniform"`: Multiplier bootstrap for sup-t uniform inference.
  - `"naive"`: Naive standard errors ignoring first-stage estimation.

# Returns
A Dictionary containing:
- `:estimate`: The estimated treatment effects.
- `:vcov_type`: The type of variance-covariance matrix computed.
- `:selected_p`: The number of factors used.

Additional returned elements depend on `vcov_type` and `type`:
- If `vcov_type = "uniform"`: `:se` (standard errors) and `:crit_val` (critical values).
- If `vcov_type = "pointwise"` or `"naive"`: `:vcov` (variance-covariance matrix).
- If `type = "gt"`: `:gt_index` (group-time indices) and `:N_tau_gt` (number of units for each group-time pair).
- If `type = "dynamic"`: `:rel_year` (relative years to treatment).
- If `return_y0 = true`: `:impute_df` (DataFrame with imputed counterfactual outcomes).

# Notes
- Requires a balanced panel.
- Treatment timing must be strictly after period `p` for all treated units.
- The number of instruments must be at least `p`.
"""
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
  vcov_type="pointwise",
)
  @assert vcov_type in ["pointwise", "uniform", "naive"] "vcov_type must be one of 'pointwise', 'uniform', or 'naive'"

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
  gt_index = zeros(T * length(uniq_g_shift), 2)
  for (l, curr_g) in enumerate(uniq_g_shift)
    curr_g = convert(Int64, curr_g)
    curr_idx = findall(g_shift .== curr_g) # All units with this g
    for i in curr_idx
      g_idx = (1 + ((l - 1) * T)):(l * T)
      N_tau_gt[g_idx] .+= 1
      gt_index[g_idx, 1] .= curr_g + min_t
      gt_index[g_idx, 2] .= uniq_t
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
  tau_gt_hat, _ = estimate_tau_gt(theta_hat_opt, p, ymat, g_shift)

  # Estimate VCOV of τ(g,t) 
  # ----
  # Multiply to make unconditional moments
  ms = 1 / (N_inf / N) * ms_theta(theta_hat_opt, p, ymat, W, idx_control)
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

  if vcov_type == "naive"
    IF = IF_tau
  else
    IF = IF_tau + IF_theta
  end

  # Aggregate effects if needed
  ret = Dict(:vcov_type => vcov_type, :selected_p => p)

  if type == "gt"
    ret[:gt_index] = gt_index
    ret[:N_tau_gt] = N_tau_gt
    ret[:estimate] = tau_gt_hat
    ret[:inf_func] = IF

    if vcov_type == "uniform"
      se_tau_gt, crit_val_tau_gt = mboot(1 / sqrt(N) * IF')
      ret[:se] = se_tau_gt
      ret[:crit_val] = crit_val_tau_gt
    else
      vcov_tau_gt = 1 / N * (IF * IF')
      ret[:vcov] = vcov_tau_gt
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

    IF_es = mat_agg_es * IF
    tau_es_hat = mat_agg_es * tau_gt_hat
    ret[:rel_year] = uniq_rel_years
    ret[:estimate] = tau_es_hat
    ret[:inf_func] = IF_es

    if vcov_type == "uniform"
      se_tau_es, crit_val_tau_es = mboot(1 / sqrt(N) * IF_es')
      ret[:se] = se_tau_es
      ret[:crit_val] = crit_val_tau_es
    else
      vcov_tau_es = 1 / N * (IF_es * IF_es')
      ret[:vcov] = vcov_tau_es
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

    IF_overall = mat_agg_overall * IF
    tau_overall_hat = mat_agg_overall * tau_gt_hat
    ret[:estimate] = tau_overall_hat
    ret[:inf_func] = IF_overall

    if vcov_type == "uniform"
      se_tau_overall, crit_val_tau_overall = mboot(1 / sqrt(N) * IF_overall')
      ret[:se] = se_tau_overall
      ret[:crit_val] = crit_val_tau_overall
    else
      vcov_tau_overall = 1 / N * (IF_overall * IF_overall')
      ret[:vcov] = vcov_tau_overall
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
    ret[:impute_df] = impute_df
  end

  return ret
end
