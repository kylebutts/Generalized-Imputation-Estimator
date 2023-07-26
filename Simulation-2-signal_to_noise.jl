# ------------------------------------------------------------------------------
# Simulation showing the effectiveness of including w_i β_t in TWFE imputation
# as a function of signal-to-noise ratio of w_i to the factor loading

using Random, Distributions, Statistics, StatsBase, LinearAlgebra
using Plots, PrettyTables
using DataFrames, Chain, FixedEffectModels, Printf
using Optim, ForwardDiff
using CSV

# NOTE: need {did2s} and {fixest} installed in R to run this code
using RCall 
@rlibrary did2s
@rlibrary fixest

include("Simulation-helpers.jl")
include("Simulation-factor_imputation.jl")

# Simulations ------------------------------------------------------------------

rng = MersenneTwister(12);

# Same factors
f = 1:8

sims = DataFrame(
  signal_to_noise = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1] 
)
sims.instrument_noise = 1 ./ sims.signal_to_noise .- 1

# parameters
N = 200;
T = 8;
T0 = 5;
b = [1.0, 1.0];
K = 2;
S = 250;
# S = 20;
s = 1;
p = 1;
true_te = [0 0 0 0 0 1 2 3];


nsims = length(sims.instrument_noise)
did2s_cov_bias = zeros(nsims, T - T0)
did2s_cov_bias_sd = zeros(nsims, T - T0)
generalized_bias = zeros(nsims, T - T0)
generalized_bias_sd = zeros(nsims, T - T0)

# for i in 1:1
for i in 1:nsims
  
  sim = sims[i, :]

  include_factor = true
  parallel_trends = false
  instrument_noise = sim.instrument_noise
  println("\n\nStarting on: 
  Include Factor = $(include_factor),
  Parallel Trends = $(parallel_trends), 
  Instrument Noise = $(instrument_noise)
  ")

  did2s_w_covariates = zeros(Float32, T-T0, S);
  did2s = zeros(Float32, T-T0, S);
  generalized = zeros(Float32, T-T0, S);

  s = 1
  while s <= S
    (s % 20 == 0) && print("\n.")

    # Generate data ----------------------------------------------------------
    data = generate_data(rng, f, include_factor, parallel_trends, instrument_noise)
    
    try
      # did2s w/ covariates Estimator ------------------------------------------
      did2s_w_covariates[:, s] = est_did2s_covariates_R(data)
      
      # Factor Imputation Estimator --------------------------------------------
      generalized[:, s] = est_factor_imputation(data)

      s = s+1
      print(".")
    catch
      print("x")
    end
  end

  for j in 1:(T-T0)
    did2s_cov_bias[i, j] = mean(did2s_w_covariates[j, :] .- true_te[T0 + j])
    did2s_cov_bias_sd[i, j] = std(did2s_w_covariates[j, :] .- true_te[T0 + j])
    generalized_bias[i, j] = mean(generalized[j, :] .- true_te[T0 + j])
    generalized_bias_sd[i, j] = std(generalized[j, :] .- true_te[T0 + j])
  end
end


sims.did2s_cov_bias_tau6 = did2s_cov_bias[:, 1]
sims.did2s_cov_bias_tau7 = did2s_cov_bias[:, 2]
sims.did2s_cov_bias_tau8 = did2s_cov_bias[:, 3]
sims.did2s_cov_bias_sd_tau6 = did2s_cov_bias_sd[:, 1]
sims.did2s_cov_bias_sd_tau7 = did2s_cov_bias_sd[:, 2]
sims.did2s_cov_bias_sd_tau8 = did2s_cov_bias_sd[:, 3]
sims.generalized_bias_tau6 = generalized_bias[:, 1]
sims.generalized_bias_tau7 = generalized_bias[:, 2]
sims.generalized_bias_tau8 = generalized_bias[:, 3]
sims.generalized_bias_sd_tau6 = generalized_bias_sd[:, 1]
sims.generalized_bias_sd_tau7 = generalized_bias_sd[:, 2]
sims.generalized_bias_sd_tau8 = generalized_bias_sd[:, 3]

sims

# Export results 
CSV.write("data/simulation-2.csv", sims)

# Plot -------------------------------------------------------------------------

using Plots, DataFrames, CSV, LaTeXStrings
# sims = DataFrame(CSV.File("simulation-2.csv"))

sims.did2s_cov_bias_tau8_upper = 
  sims.did2s_cov_bias_tau8 .+ 1.96 * sims.did2s_cov_bias_sd_tau8
sims.did2s_cov_bias_tau8_lower = 
  sims.did2s_cov_bias_tau8 .- 1.96 * sims.did2s_cov_bias_sd_tau8
sims.generalized_bias_tau8_upper = 
  sims.generalized_bias_tau8 .+ 1.96 * sims.generalized_bias_sd_tau8
sims.generalized_bias_tau8_lower = 
  sims.generalized_bias_tau8 .- 1.96 * sims.generalized_bias_sd_tau8


pgfplotsx()
# gr(size = (600, 300))

alice = RGBA(16/255, 120/255, 149/255, 1)
ruby = RGBA(154/255, 36/255, 21/255, 1)

bias_plot = plot(legend = :bottomright, size = (600, 300) .* 0.75)
xlabel!("Signal to Noise Ratio")
ylabel!("Average Bias")
plot!(
  sims.signal_to_noise,
  sims.did2s_cov_bias_tau8, 
  ribbon = 1.96 .* sims.did2s_cov_bias_sd_tau8,
  fill = alice,
  labels = ""
)
plot!(
  sims.signal_to_noise,
  sims.did2s_cov_bias_tau8, 
  linewidth = 1,
  linecolor = alice,
  markershape = :circle,
  markercolor = alice,
  markersize = 2.2,
  markerstrokewidth = 0,
  markerstrokealpha = 0,
  labels = L"TWFE Imputation with $w_i \beta_t$"
)
plot!(
  sims.signal_to_noise,
  sims.generalized_bias_tau8, 
  ribbon = 1.96 .* sims.generalized_bias_sd_tau8,
  fill = ruby,
  labels = ""
)
plot!(
  sims.signal_to_noise,
  sims.generalized_bias_tau8,
  linewidth = 1,
  linecolor = ruby,
  markershape = :square,
  markercolor = ruby,
  markersize = 2.2,
  markerstrokewidth = 0,
  markerstrokealpha = 0,
  labels = "Factor Model Imputation"
)

plot(bias_plot, size = (600, 400) .* 1)

# savefig(bias_plot, "figures/monte-bias_signal_to_noise.tex")

