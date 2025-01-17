function estimate_tau_gt(
  theta::AbstractVector, p::Int, ymat::AbstractMatrix, g_shift::AbstractVector
)
  T = size(ymat, 1)
  Fhat = [reshape(theta, T - p, p); -1.0 * I(p)]
  uniq_g_shift = sort(unique(g_shift))
  uniq_g_shift = filter(x -> x .!== Inf, uniq_g_shift)

  # Impute y0_hat
  y_diff = zeros(size(ymat))
  tau_gt_hat = zeros(T * length(uniq_g_shift))
  N_tau_gt = zeros(T * length(uniq_g_shift))

  for (l, curr_g) in enumerate(uniq_g_shift)
    curr_idx = findall(g_shift .== curr_g) # All units with this g
    Fhat_pre = Fhat[1:convert(Int64, curr_g), :]
    PFhat = Fhat * ((Fhat_pre' * Fhat_pre) \ Fhat_pre')
    for i in curr_idx
      y_diff[:, i] = ymat[:, i] - PFhat * ymat[1:convert(Int64, curr_g), i]
      tau_gt_hat[(1 + ((l - 1) * T)):(l * T)] += y_diff[:, i]
      N_tau_gt[(1 + ((l - 1) * T)):(l * T)] .+= 1
    end
  end

  tau_gt_hat = tau_gt_hat ./ N_tau_gt
  return tau_gt_hat, N_tau_gt
end

function impute_y0(
  theta::AbstractVector, p::Int, ymat::AbstractMatrix, g_shift::AbstractVector
)
  T = size(ymat, 1)
  Fhat = [reshape(theta, T - p, p); -1.0 * I(p)]
  uniq_g_shift = sort(unique(g_shift))
  uniq_g_shift = filter(x -> x .!== Inf, uniq_g_shift)

  # Impute y0_hat
  y0_hat = zeros(size(ymat))

  for (l, curr_g) in enumerate(uniq_g_shift)
    curr_idx = findall(g_shift .== curr_g) # All units with this g
    Fhat_pre = Fhat[1:convert(Int64, curr_g), :]
    PFhat = Fhat * ((Fhat_pre' * Fhat_pre) \ Fhat_pre')
    for i in curr_idx
      y0_hat[:, i] = PFhat * ymat[1:convert(Int64, curr_g), i]
    end
  end

  return y0_hat
end

function ms_tau_gt(
  theta::AbstractVector,
  p::Int,
  tau_gt::AbstractVector,
  ymat::AbstractMatrix,
  g_shift::AbstractVector,
)
  T, N_units = size(ymat)
  Fhat = [reshape(theta, T - p, p); -1.0 * I(p)]
  uniq_g_shift = sort(unique(g_shift))
  uniq_g_shift = filter(x -> x .!== Inf, uniq_g_shift)

  # Impute y0_hat
  y_diff = zeros(eltype(theta), size(ymat, 1), size(ymat, 2))
  ms = zeros(eltype(theta), N_units, length(tau_gt))
  N_tau_gt = zeros(length(tau_gt))
  for (l, curr_g) in enumerate(uniq_g_shift)
    curr_g = convert(Int64, curr_g)
    curr_idx = findall(g_shift .== curr_g) # All units with this g
    Fhat_pre = Fhat[1:convert(Int64, curr_g), :]
    PFhat = Fhat * ((Fhat_pre' * Fhat_pre) \ Fhat_pre')
    for i in curr_idx
      y_diff[:, i] = ymat[:, i] - PFhat * ymat[1:curr_g, i]
      gt_idx = (1 + ((l - 1) * T)):(l * T)
      ms[i, gt_idx] = y_diff[:, i] - tau_gt[gt_idx]
      N_tau_gt[gt_idx] .+= 1
    end
  end

  ms ./= (N_tau_gt / N_units)'
  return ms
end

function m_tau_gt_bar(
  theta::AbstractVector,
  p::Int,
  tau_gt::AbstractVector,
  ymat::AbstractMatrix,
  g_shift::AbstractVector,
)
  T, N_units = size(ymat)
  Fhat = [reshape(theta, T - p, p); -1.0 * I(p)]
  uniq_g_shift = sort(unique(g_shift))
  uniq_g_shift = filter(x -> x .!== Inf, uniq_g_shift)

  # Impute y0_hat
  y_diff = zeros(eltype(theta), size(ymat, 1), size(ymat, 2))
  m_bar = zeros(eltype(theta), length(tau_gt))
  N_tau_gt = zeros(length(tau_gt))
  for (l, curr_g) in enumerate(uniq_g_shift)
    curr_g = convert(Int64, curr_g)
    curr_idx = findall(g_shift .== curr_g) # All units with this g
    Fhat_pre = Fhat[1:convert(Int64, curr_g), :]
    PFhat = Fhat * ((Fhat_pre' * Fhat_pre) \ Fhat_pre')
    for i in curr_idx
      y_diff[:, i] = ymat[:, i] - PFhat * ymat[1:curr_g, i]
      gt_idx = (1 + ((l - 1) * T)):(l * T)
      m_bar[gt_idx] .+= y_diff[:, i]
      N_tau_gt[gt_idx] .+= 1
    end
  end

  m_bar .-= tau_gt
  m_bar ./= (N_tau_gt / N_units)
  return m_bar
end

