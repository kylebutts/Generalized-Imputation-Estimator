# Brown and Butts - Within transformation --------------------------------------
function within_transform(ymat, idx_control, N_pre)
  # never-treated cross-sectional averages of y
  ymat .-= mean(ymat[:, idx_control]; dims=2)

  # unit-averages pre-T_0
  ymat .-= mean(ymat[1:N_pre, :]; dims=1)

  # never-treated average pre-T_0
  ymat .+= mean(ymat[1:N_pre, idx_control])
  return ymat
end
