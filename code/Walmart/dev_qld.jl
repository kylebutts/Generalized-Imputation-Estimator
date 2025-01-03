# Event Study Estimator for Walmart data 
# Really important: This assumes balanced panels
using BlockDiagonals
using Calculus
using CategoricalArrays
using CSV
using DataFrames
using Dates
using Distributions
using FixedEffectModels
using ForwardDiff
using GLM
using LinearAlgebra
using LineSearches
using Optim
using Plots
using ProjectRoot
using Random
using SparseArrays
using Statistics
using StatsBase
using Kronecker
using BenchmarkTools

here = @projectroot()
cd(here)
include("$here/code/Walmart/dev_qld_helpers.jl")

# Set parameters
YEARS = collect(1977:1999)
T = length(YEARS)
T0 = 1985
outcome_var = :log_retail_emp
p = 2

# %%
df = DataFrame(
  CSV.File(
    "data/County_Business_Patterns/sample_basker_YEARS_$(minimum(YEARS))_$(maximum(YEARS))_T0_$(T0).csv",
  ),
)

# ----
#
#
#
# Process Data 
sort!(df, [:fips, :year])
ytilde = within_transform(df, outcome_var, :g, :year, T0)
df.ytilde = ytilde
fips = df.fips
year = df.year

# _shift
g = df.g[year .== 1977] .- minimum(year) # indexed 1:T
keep_obs = df.g .!= Inf
uniq_rel_years = sort(unique(df[(df.g .!= Inf), :rel_year]))
# 1 ... length(uniq_rel_years)
rel_year_shift = 1 .+ abs(minimum(uniq_rel_years)) .+ df.rel_year

# # Second-stage regression (event-study indicators)
# X_ss = similar(df.ytilde, length(df.ytilde), length(uniq_rel_years))
# for (i, rel_yr) in enumerate(uniq_rel_years)
#   X_ss[:, i] = df.rel_year .== rel_yr
# end

# Find untreated units
N_units = length(unique(fips))
idx_control = findall(df.g[year .== 1977] .== Inf)
N_inf = length(idx_control)

# Number of instruments
instruments = [
  :share_pop_ind_manuf,
  :share_pop_poverty_78_below,
  :share_pop_poverty_78_above,
  :share_pop_emp_private,
  :share_pop_emp_government,
  :share_school_col,
  :share_school_hs,
]
N_instruments = length(instruments)
W = Matrix(df[year .== 1977, instruments])

# Number of covariates
k = 0

# Fills in by column:
# reshape([1 2 3 4 5 6 7 8 9], 3, 3)
ymat = reshape(ytilde, T, N_units)
ymat_inf = ymat[:, idx_control]
W_inf = W[idx_control, :]

# each column consists of \sum_i y_{it} ⊗ w_i for each t
# multiplied by P(g_i = \infty) * 1/N 
ymat_W_inf = ymat_inf * W_inf * (N_inf / N_units) / N_units

# # TEMP ----
# # Ensure correctness
# m1 = m_theta_bar_idiomatic(parameters, p, N_units, ymat_inf, W_inf)
# m2 = fast_m_theta_bar(parameters, p, ymat_W_inf)
# m1 ≈ vec(transpose(m2))
# 
# # benchmark
# # 192.869 μs
# @benchmark m_theta_bar_idiomatic(parameters, p, N_units, ymat_inf, W_inf)
# # 447.649 ns = 0.447649 μs
# @benchmark fast_m_theta_bar(parameters, p, ymat_W_inf)
# 
# @benchmark ForwardDiff.jacobian(
#   parameters -> m_theta_bar_idiomatic(parameters, p, N_units, ymat_inf, W_inf), parameters
# )
# @benchmark ForwardDiff.jacobian(
#   parameters -> vec(transpose(fast_m_theta_bar(parameters, p, ymat_W_inf))), parameters
# )
# # End TEMP ----

# ----
#
#
#
# Two-step GMM Estimation of QLD parameters
parameters = [ones((T - p) * p); zeros(length(uniq_rel_years))]
weight = Matrix(1.0 * I((T - p) * N_instruments))

if p > 0
  optres = optimize(
    x -> obj_theta(x, weight, p, N_units, ymat_W_inf),
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
ms = ms_theta(parameters_hat, p, N_units, ymat_inf, W_inf)
Wopt = pinv(ms' * ms / N_units)

if p > 0
  # Having problems with numeric stability when doing all the 1/N_inf correctly
  RESCALE = 15.266221542835233 * 83.84513366993183
  optres_opt = optimize(
    x -> 1 / RESCALE * obj_theta(x, Wopt, p, N_units, ymat_W_inf),
    parameters_hat,
    BFGS(; linesearch=BackTracking()),
    Optim.Options(; iterations=10000);
    autodiff=:forward,
  )

  parameters_hat_opt = Optim.minimizer(optres_opt)
  theta_hat_opt = parameters_hat_opt[1:((T - p) * p)]
  J = Optim.minimum(optres_opt) * RESCALE
else
  J = obj_theta(parameters_hat, p, Ωinv_theta, ytilde, idx_control, W, n_units, T)
end

# Testing for $p$ following Ahn, Lee, and Schmidt (2013)
# Χ^2((T - p_0) (q - p_0) - k)
# Where T is the number of periods, p_0 is the number of factors under the null, q is the number of instruments, and k is the number of covariates
dist_hansen_sargent = Chisq((T - p) * (N_instruments - p) - k)
p_value_hansen_sargent = 1 - cdf(dist_hansen_sargent, J)

# ----
#
#
#
# Estimate event-study

# TEMP ----
Hθ = [I(T - p) reshape(parameters_hat[1:((T - p) * p)], T - p, p)];
Fθ = [reshape(parameters_hat[1:((T - p) * p)], T - p, p); -1.0 * I(p)];
# Confirm I have constructed the QLD matrices correctly
Hθ * Fθ
# END TEMP ----

theta_opt = parameters_hat_opt[1:((T - p) * p)]
optres_tau = optimize(
  x -> obj_tau(x, theta_opt, ymat, g, rel_year_shift),
  parameters_hat_opt,
  BFGS(; linesearch=BackTracking()),
  Optim.Options(; iterations=10000);
  autodiff=:forward,
)

params_final = Optim.minimizer(optres_tau)
tau_hat = params_final[(((T - p) * p) + 1):end]

tau_hat[uniq_rel_years .< 0]
tau_hat[uniq_rel_years .>= 0]

# ----
#
#
#
# Asymptotic variance estimation
# Following Newey and McFadden (1994)
# The second-step is a method of moments estimator using a first-step GMM estimate. 
# This is derived in (6.6) of Newey and McFadden and expanding the GMM first-order conditions from the first-step.
# 
# IF(m(z, γ)) = influence-function of first-step GMM parameters
# g(z, τ, γ) = second-step moments (ms_tau)
# G_γ(z, τ, γ) = ∇_γ gbar(z, τ, γ) (jacobian of mbar_tau)
# G_τ(z, τ, γ) = ∇_τ gbar(z, τ, γ) (jacobian of mbar_tau)
#
# Two-step influence function is given by:
# G_τ^{-1} [g(z) - G_γ IF(m(z, γ))]
Mbar_theta = ForwardDiff.jacobian(
  x -> vec(transpose(fast_m_theta_bar(x, p, ymat_W_inf))), params_final
)
Mbar_theta = Mbar_theta[:, 1:length(theta_opt)] # only w.r.t theta
m_theta = zeros(N_units, size(Mbar_theta, 1))
m_theta[idx_control, :] = ms_theta(params_final, p, N_units, ymat_inf, W_inf)

IF_theta_solver = pinv(Mbar_theta' * Wopt * Mbar_theta) * Mbar_theta' * Wopt
# IF_theta_solver * m_theta[4, :]


m_tau = ms_tau(params_final, ymat, g, rel_year_shift)
Gbar = ForwardDiff.jacobian(x -> m_bar_tau(x, ymat, g, rel_year_shift), params_final)
Gbar_theta = Gbar[:, 1:length(theta_opt)]
Gbar_tau = Gbar[:, (1 + length(theta_opt)):end]

IF_tau_solver = pinv(Gbar_tau)
# IF_tau_solver * m_tau[1, :]

IF_tau = zeros(size(m_tau, 1), size(m_tau, 2))
IF_tau_ignore_fs = zeros(size(m_tau, 1), size(m_tau, 2))
for i in axes(m_theta, 1)
  IF_tau[i, :] =
    IF_tau_solver *
    (1 / N_units * m_tau[i, :] + Gbar_theta * IF_theta_solver * 1 / N_units * m_theta[i, :])
  IF_tau_ignore_fs[i, :] = IF_tau_solver * (1 / N_units * m_tau[i, :])
end
# Note: I have double-checked all my Ns and I have it all correct. 
# See (6.6) of Newey and McFadden

function get_ses(IF)
  return sqrt.(diag(IF' * IF))
end
get_ses(IF_tau)
get_ses(IF_tau_ignore_fs)
# get_ses(IF_tau) ./ get_ses(IF_tau_ignore_fs)


