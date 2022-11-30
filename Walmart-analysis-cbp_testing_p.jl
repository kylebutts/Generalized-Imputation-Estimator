# ------------------------------------------------------------------------------
# Event Study Estimator for Walmart data 

using Random, Distributions, Statistics, LinearAlgebra
using Plots, PrettyTables, GLM, CategoricalArrays
using CSV, DataFrames, Chain, FixedEffectModels, Printf
using Optim, LineSearches, RCall, Calculus, BlockDiagonals
using BenchmarkTools
import ForwardDiff

# Set parameters
YEARS = collect(1977:1999)
T0 = 1985
# YEARS = collect(1988:1999)
# T0 = 1995
file = "data/sample_basker_YEARS_$(minimum(YEARS))_$(maximum(YEARS))_T0_$(T0).csv"
T = length(YEARS)
p = 2
outcome = "retail"
outcome_var = "log_$(outcome)_emp"
lims = outcome == "retail" ? [-0.175, 0.3] : [-0.4, 0.2]


# Load data --------------------------------------------------------------------

sample = DataFrame(CSV.File(file))
sample.y = sample[:, outcome_var]
sample.share_pop_ind_manuf = sample.share_pop_ind_manuf_nondurable + sample.share_pop_ind_manuf_durable

sort!(sample, [:state_fips, :county_fips, :year])

# median retail employment
# sample[sample.year .== 1977, :retail_emp] |> median
# sample[sample.year .== 1977, :wholesale_emp] |> median

# Within Transformation --------------------------------------------------------

# never-treated cross-sectional averages of y
sample = sample |>
  (x -> groupby(x, :year)) |>
  (x -> transform(x,
    [:y, :g] => ((y, g) -> mean(y[g .== Inf])) => :yinf_t
  ))

# unit-averages pre-T0
sample = sample |>
  x -> groupby(x, :fips) |>
  x -> transform(x,
    [:y, :year] => ((y, year) -> mean(y[year .<= T0])) => :yi_pre
  )

# never-treated average pre-T0
sample.yinf_pre .= mean(sample[(sample.g .== Inf) .* (sample.year .<= T0), :y])

# Create ytilde
sample.ytilde = sample.y - sample.yinf_t - sample.yi_pre + sample.yinf_pre


# Export for Nick
# untreated = sample[sample.g .== Inf, [:ytilde, :share_pop_ind_manuf,
# :share_pop_poverty_78_below, :share_pop_poverty_78_above, 
# :share_pop_emp_private, :share_pop_emp_government, :share_school_col, :share_school_hs]]
# CSV.write("data/untreated_YEARS_$(minimum(YEARS))_$(maximum(YEARS))_T0_$(T0).csv", untreated)
# 
# treated = sample[sample.g .!= Inf, [:ytilde, :rel_year]]
# CSV.write("data/treated_YEARS_$(minimum(YEARS))_$(maximum(YEARS))_T0_$(T0).csv", treated)
# 
# g1986 = sample[sample.g .!= 1986, [:ytilde]]
# CSV.write("data/g1986_YEARS_$(minimum(YEARS))_$(maximum(YEARS))_T0_$(T0).csv", g1986)



## GMM estimate factor ---------------------------------------------------------

# untreated = sample[(sample.g .== Inf), :]  
# Θ = ones(T - p, p)

instruments = [
  :share_pop_ind_manuf,
  :share_pop_poverty_78_below, :share_pop_poverty_78_above, 
  :share_pop_emp_private, :share_pop_emp_government, 
  # :share_pop_income_79_bin_1, :share_pop_income_79_bin_2,
  # :share_pop_income_79_bin_3, :share_pop_income_79_bin_4,
  # :share_pop_income_79_bin_5, :share_pop_income_79_bin_6,
  # :share_pop_income_79_bin_7, :share_pop_income_79_bin_8,
  # :share_pop_income_79_bin_9, :share_pop_income_79_bin_10,
  # :share_pop_income_79_bin_11, :share_pop_income_79_bin_12,
  # :share_pop_income_79_bin_13, :share_pop_income_79_bin_14,
  # :share_pop_income_79_bin_15, :share_pop_income_79_bin_16,
  :share_school_col, :share_school_hs
]

n_instruments = length(instruments)
n_units = sample.fips |> unique |> length
rel_years = sample[(sample.g .!= Inf), :rel_year] |> unique |> sort

# Estimate theta ---------------------------------------------------------------

function m_thetas(parms::AbstractVector, p::Int64)
  
  if p > 0
    theta = parms[1:((T - p) * p)] 
    HΘ = [I(T - p); reshape(theta, T-p, p)']
  else 
    HΘ = 1.0 * I(T) |> Matrix
  end

  m_theta = zeros(n_units, (T - p) * n_instruments)
  N_inf = 0
  i = 1
  for id in unique(sample.fips)
    W = sample[sample.fips .== id, instruments][1, :] |> Vector
    ytilde = sample[sample.fips .== id, :ytilde]

    D_inf = sample[sample.fips .== id, :g] |> 
      (x -> (x[1] == Inf))
    
    N_inf += D_inf

    m_theta[i, :] = D_inf * kron((HΘ' * ytilde), W)
    i += 1
  end

  # divide by 1 / P(D_i,inf) 
  m_theta = m_theta / (N_inf / n_units)

  return m_theta
end

function m_theta_bar(parms::AbstractVector, ytilde::Vector{Float64}, W::Matrix{Float64}, g::Vector{Float64}, p::Int64)

  if p > 0
    theta = parms[1:((T - p) * p)] 
    HΘ = [I(T - p); reshape(theta, T-p, p)']
  else 
    HΘ = 1.0 * I(T) |> Matrix
  end

  m_theta = zeros((T - p) * n_instruments)  
  N_inf = 0
  
  for i in 1:n_units
    W_unit = W[1 + (i-1)*T, :]
    ytilde_unit = ytilde[(1 + (i-1)*T):(i*T)]

    D_inf = g[1 + (i-1)*T] == Inf
    N_inf += D_inf

    m_theta_idx = D_inf * kron((HΘ' * ytilde_unit), W_unit)

    m_theta = m_theta + m_theta_idx
  end

  # 1 / P(D_i,inf) 
  # m_theta = m_theta / (N_inf / n_units)

  # return the avearage value of m_theta, m_theta_bar
  # m_theta = m_theta / n_units

  m_theta = m_theta / N_inf

  return m_theta
end

function obj_theta(parms::AbstractVector, weight::Matrix{Float64}, ytilde::Vector{Float64}, W::Matrix{Float64}, g::Vector{Float64}, p::Int64)
  m_bar = m_theta_bar(parms, ytilde, W, g, p)
  return (m_bar' * weight * m_bar)
end

# overidentifying_restrictions = (T - p) * (n_instruments - p)
parms = [ones((T - p) * p); zeros(length(rel_years))]
weight = 1.0 * I((T - p) * n_instruments) |> Matrix

# mean(m_thetas(parms), dims = 1)
# m_theta_bar(parms)

# Step 1. Initial estimate
if p > 0
  optres = optimize(
    x -> obj_theta(x, weight, sample.ytilde, Matrix(sample[:, instruments]), sample.g, p), parms, 
    LBFGS(); autodiff = :forward
  )

  parms_hat = Optim.minimizer(optres)
  theta_hat = parms_hat[1:((T - p) * p)] 
else 
  obj_theta(parms, weight, sample.ytilde, Matrix(sample[:, instruments]), sample.g, p)
end

# obj_value = Optim.minimum(optres)
# print("p = $(p); objective function value = $(obj_value)")

# Step 2. Optimal weighting matrix 
Ωinv_theta = pinv( (m_thetas(parms_hat, p)' * m_thetas(parms_hat, p)) / n_units ) 

if p > 0
  optres_opt = optimize(
    x -> obj_theta(x, Ωinv_theta, sample.ytilde, Matrix(sample[:, instruments]), sample.g, p), parms_hat, 
    BFGS(; linesearch = BackTracking()), Optim.Options(iterations=10000); 
    autodiff = :forward
  )

  parms_hat_opt = Optim.minimizer(optres_opt)
  theta_hat_opt = parms_hat_opt[1:((T - p) * p)] 
  
  J = Optim.minimum(optres_opt)
else 
  J = obj_theta(parms, Ωinv_theta, sample.ytilde, Matrix(sample[:, instruments]), sample.g, p)
end

# From Equation (7) and Propoition 1 in Ahn, Lee, and Schmidt (2013)
p_value = 1 - cdf(Chisq((T - p) * (n_instruments - p)), J * n_units)
println("p = $(p); objective function value = $(J); p-value = $(p_value)")

# Retail - 11/23
# p = 0; objective function value = 0.19375520523380194; p-value = 1.565445430551371e-5
# p = 1; objective function value = 0.14650594060704658; p-value = 0.0012370648020560981
# p = 2; objective function value = 0.09517373124970177; p-value = 0.1326886245943757
# p = 3; objective function value = 0.0666109340899677; p-value = 0.33383328557602265
# p = 4; objective function value = 0.02629459726616495; p-value = 0.9944891251456437

# Wholesale
# p = 0; objective function value = 0.14890180055470742; p-value = 0.049099915394573035
# p = 1; objective function value = 0.10518532685717184; p-value = 0.4000556784516629
# p = 2; objective function value = 0.0660008619398643; p-value = 0.9238535531764781
# p = 3; objective function value = 0.03513302983950854; p-value = 0.999395353135183


# Estimate τ^ℓ ---------------------------------------------------------------
FΘ = [reshape(theta_hat_opt, T-p, p); (-1 * I(p))]
FΘ = Matrix(FΘ)
(FΘ * inv(FΘ' * FΘ) * FΘ')

theta_hat_opt
treated_ids = unique(sample[sample.g .< Inf, :fips])

sample.y0hat = zeros(size(sample, 1))

for id in treated_ids
  unit = sample[(sample.fips .== id), :]
  g_id = Int(unit.g[1])
    
  y_pre = unit[(unit.year .< g_id), :ytilde]  
  FΘ_pre = FΘ[1:(g_id - minimum(YEARS)), :]

  sample[(sample.fips .== id), :y0hat] = 
    FΘ * inv(FΘ_pre' * FΘ_pre) * FΘ_pre' * y_pre
end

sample.tau_hat = sample.ytilde .- sample.y0hat

# ATT(1986, t)
# sample[(sample.g .== 1986), :] |>
#   (x -> groupby(x, :year)) |>
#   (x -> combine(x, :tau_hat => mean)) |>
#   (x -> x.tau_hat_mean)

# Factor Estimate
sample.rel_year_cat = CategoricalArray(sample.rel_year)
lm_tau = lm(
  @formula(tau_hat ~ 0 + rel_year_cat), 
  sample[sample.g .< Inf, :]
)


# Extract point estimates 
rel_year = collect(-22:13)

idx_pre  = (rel_year .<  0)
idx_post = (rel_year .>= 0)
factor_estimate = coef(lm_tau)
factor_se = coeftable(lm_tau).cols[2]

factor_plot = plot(size = (600, 400) .* 0.75)
hline!([0], linestyle = :dash, linecolor = :black)
# title!("Factor Model Imputation")
xlabel!("Event Time")
ylabel!("Coefficient")
ylims!(lims[1], lims[2])
scatter!(
  rel_year[idx_pre],
  factor_estimate[idx_pre],
  yerror = 1.96 .* factor_se[idx_pre],
  # yerror = (
  #   factor_estimate[idx_pre] - factor_lower_95[idx_pre], 
  #   factor_upper_95[idx_pre] - factor_estimate[idx_pre]
  # ),
  legend = false,
  markersize = 2, 
  markerstrokecolor = RGBA(154/255, 36/255, 21/255, 1),
  markercolor = RGBA(154/255, 36/255, 21/255, 1)
)
scatter!(
  rel_year[idx_post],
  factor_estimate[idx_post],
  yerror = 1.96 .* factor_se[idx_post],
  # yerror = (
  #   factor_estimate[idx_post] - factor_lower_95[idx_post], 
  #   factor_upper_95[idx_post] - factor_estimate[idx_post]
  # ),
  legend = false,
  markersize = 2, 
  markerstrokecolor = RGBA(16/255, 120/255, 149/255, 1),
  markercolor = RGBA(16/255, 120/255, 149/255, 1)
)

savefig(factor_plot, "figures/factor_$(outcome)_p_$(p)_naive_se.tex")




