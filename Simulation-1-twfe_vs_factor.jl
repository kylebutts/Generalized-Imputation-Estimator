# ------------------------------------------------------------------------------
# 10/18/2022: Preliminary sims for CCEDID paper to compare imputed CCE 
# estimator to TWFE. Currently there is no treatment effect

using Random, Distributions, Statistics, LinearAlgebra, StatsBase
using Plots, PrettyTables
using DataFrames, Chain, FixedEffectModels, Printf
using Optim, ForwardDiff

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

# parameters
N = 200
T = 8
T0 = 5
b = [1.0, 1.0]
K = 2
S = 150;
s = 1;
p = 1;
true_te = [0 0 0 0 0 1 2 3];


sims = DataFrame(
  include_factor = [false, true, true],
  parallel_trends = [true, true, false],
  instrument_noise = [1, 1, 1]
)

for sim in eachrow(sims)
 
  include_factor = sim.include_factor
  parallel_trends = sim.parallel_trends
  instrument_noise = sim.instrument_noise
  println("Starting on: 
  Include Factor = $(include_factor),
  Parallel Trends = $(parallel_trends), 
  Instrument Noise = $(instrument_noise)
  ")

  twfe = zeros(Float32, T-T0, S);
  did2s_est = zeros(Float32, T-T0, S);
  did2s_w_covariates = zeros(Float32, T-T0, S);
  generalized = zeros(Float32, T-T0, S);

  for s in 1:S
    # print(".")
    
    # Generate data ------------------------------------------------------------
    data = generate_data(rng, f, include_factor, parallel_trends, instrument_noise)
    # data[[1:3; 598:600], :]


    # TWFE Estimator -----------------------------------------------------------
    twfe[:, s] = est_twfe(data)

    # did2s Imputation Estimator -----------------------------------------------
    did2s_est[:, s] = est_did2s_R(data)

    # did2s w/ covariates Estimator --------------------------------------------
    did2s_w_covariates[:, s] = est_did2s_covariates_R(data)

    # Factor Imputation Estimator ----------------------------------------------
    generalized[:, s] = est_factor_imputation(data)

  end


  # Loop through each col
  results = fill("", 4, 7)
  results[:, 1] = [
    "TWFE"; "TWFE Imputation"; "TWFE Imputation with Covariates"; 
    "Factor Imputation" 
  ]

  for i in 1:(T-T0)
    results[1, i*2] = 
      @sprintf "%.2f" mean(twfe[i, :] .- true_te[T0 + i])
    results[1, i*2 + 1] = 
      @sprintf "%.2f" sum((twfe[i, :] .- true_te[T0 + i]).^2) / S

    # results[2, i*2] = 
    #   @sprintf "%.2f" mean(twfe_w_covariates[i, :] .- true_te[T0 + i])
    # results[2, i*2 + 1] = 
    #   @sprintf "%.2f" sum((twfe_w_covariates[i, :] .- true_te[T0 + i]).^2) / S

    results[2, i*2] = 
      @sprintf "%.3f" mean(did2s_est[i, :] .- true_te[T0 + i])
    results[2, i*2 + 1] = 
      @sprintf "%.3f" sum((did2s_est[i, :] .- true_te[T0 + i]).^2) / S

    results[3, i*2] = 
      @sprintf "%.3f" mean(did2s_w_covariates[i, :] .- true_te[T0 + i])
    results[3, i*2 + 1] = 
      @sprintf "%.3f" sum((did2s_w_covariates[i, :] .- true_te[T0 + i]).^2) / S

    results[4, i*2] = 
      @sprintf "%.3f" mean(generalized[i, :] .- true_te[T0 + i])
    results[4, i*2 + 1] = 
      @sprintf "%.3f" sum((generalized[i, :] .- true_te[T0 + i]).^2) / S
  end


  pretty_table(
    results, 
    header = ["Estimator", "Bias τ_1", "MSE τ_1", "Bias τ_2", "MSE τ_2", "Bias τ_3", "MSE τ_3"]
  )

  # Export results
  tex = matrix_to_table(results)

  dgp = include_factor == true ? "factor" : "twfe"

  filename = "tables/simulation-dgp_$(dgp)-pt_$(parallel_trends)-iv_noise_$(instrument_noise).tex"

  if !isfile(filename) 
    touch(filename)
  end

  open(filename, "w") do file
    write(file, tex)
  end
end






