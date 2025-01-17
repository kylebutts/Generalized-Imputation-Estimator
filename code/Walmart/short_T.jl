#' ---
#' title: 'Short T0'
#' --- 
#' This script uses just 1986-1992 treated units (trimming sample a bunch). 
#' This might improve efficiency with more “never treated” while highlighting 
#' the advantage of the method (short panels)

# %% 
#| label: setup
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

here = @projectroot()
cd(here)
include("$here/code/qld_helpers.jl")

# Set parameters
YEARS = collect(1977:1999)
T0 = 1985
file = "data/County_Business_Patterns/sample_basker_YEARS_$(minimum(YEARS))_$(maximum(YEARS))_T0_$(T0).csv"

main = DataFrame(CSV.File(file))

# Sample to years 1977:1992
df = main[main.year .<= 1992, :]
df[df.g .> 1992, :rel_year] .= -Inf
df[df.g .> 1992, :g] .= Inf

# %% 
YEARS = collect(1977:1992)
T = length(YEARS)
rel_years = begin
  x = df[df.rel_year .!= -Inf, :rel_year]
  collect(minimum(x):maximum(x))
end

# %%
# Order by fips, then by year
sort!(df, [:fips, :year])

# %% Estimate
tau_hat_retail_p_0 = est_te(
  df, 0, :log_retail_emp; export_y0=true, outcome="retail", quick_plot=true
)
tau_hat_retail_p_1 = est_te(
  df, 1, :log_retail_emp; export_y0=true, outcome="retail", quick_plot=true
)
tau_hat_retail_p_2 = est_te(
  df, 2, :log_retail_emp; export_y0=true, outcome="retail", quick_plot=true
)
tau_hat_retail_p_3 = est_te(
  df, 3, :log_retail_emp; export_y0=true, outcome="retail", quick_plot=true
)
tau_hat_retail_p_4 = est_te(
  df, 4, :log_retail_emp; export_y0=true, outcome="retail", quick_plot=true
)

# %% 
tau_hat_wholesale_p_0 = est_te(
  df, 0, :log_wholesale_emp; export_y0=true, outcome="wholesale"
)
tau_hat_wholesale_p_1 = est_te(
  df, 1, :log_wholesale_emp; export_y0=true, outcome="wholesale"
)
tau_hat_wholesale_p_2 = est_te(
  df, 2, :log_wholesale_emp; export_y0=true, outcome="wholesale"
)
tau_hat_wholesale_p_3 = est_te(
  df, 3, :log_wholesale_emp; export_y0=true, outcome="wholesale"
)
tau_hat_wholesale_p_4 = est_te(
  df, 4, :log_wholesale_emp; export_y0=true, outcome="wholesale"
)

# %% 
# Ahn, Lee, and Schmidt procedure for selecting p
println("""
Checking for p for retail employment: 
p = 0; Hansen-Sargant Test p-value = $(tau_hat_retail_p_0.sargan[2])
p = 1; Hansen-Sargant Test p-value = $(tau_hat_retail_p_1.sargan[2])
p = 2; Hansen-Sargant Test p-value = $(tau_hat_retail_p_2.sargan[2])
p = 3; Hansen-Sargant Test p-value = $(tau_hat_retail_p_3.sargan[2])
p = 4; Hansen-Sargant Test p-value = $(tau_hat_retail_p_4.sargan[2])

Checking for p for wholesale employment: 
p = 0; Hansen-Sargant Test p-value = $(tau_hat_wholesale_p_0.sargan[2])
p = 1; Hansen-Sargant Test p-value = $(tau_hat_wholesale_p_1.sargan[2])
p = 2; Hansen-Sargant Test p-value = $(tau_hat_wholesale_p_2.sargan[2])
p = 3; Hansen-Sargant Test p-value = $(tau_hat_wholesale_p_3.sargan[2])
p = 4; Hansen-Sargant Test p-value = $(tau_hat_wholesale_p_4.sargan[2])
""")

# %% Chose p based on p-values
retail_p_values = (
  tau_hat_retail_p_0.sargan[2],
  tau_hat_retail_p_1.sargan[2],
  tau_hat_retail_p_2.sargan[2],
  tau_hat_retail_p_3.sargan[2],
  tau_hat_retail_p_4.sargan[2],
)
retail_p = findfirst(retail_p_values .> 0.1) - 1
retail_est = eval(Symbol("tau_hat_retail_p_$(retail_p)"))

wholesale_p_values = (
  tau_hat_wholesale_p_0.sargan[2],
  tau_hat_wholesale_p_1.sargan[2],
  tau_hat_wholesale_p_2.sargan[2],
  tau_hat_wholesale_p_3.sargan[2],
  tau_hat_wholesale_p_4.sargan[2],
)
wholesale_p = findfirst(wholesale_p_values .> 0.1) - 1
wholesale_est = eval(Symbol("tau_hat_wholesale_p_$(wholesale_p)"))


# %% Retail bootstrap (takes about 10 minutes)
# display current time
start_time = Dates.now()
println("Current Time: $(Dates.format(start_time, "HH:MM:SS"))")
Random.seed!(123)

B = 1000
ids = unique(df[:, :fips])
retail_ests = zeros(B, length(rel_years))
retail_ests[1, :] = retail_est.est[1]'

b = 2
while b < B
  b % 50 == 0 && println("On Bootstrap Iteration: $(b)")

  resample_ids = StatsBase.sample(ids, length(ids))
  resample_df = vcat([df[df.fips .== id, :] for id in resample_ids]...)
  resample_df.fips = repeat(1:length(ids); inner=T)

  try
    boot_est = est_te(
      resample_df,
      retail_p,
      :log_retail_emp;
      export_y0=false,
      outcome="retail",
      quick_plot=false,
    )
    retail_ests[b, :] = boot_est.est[1]'
    b += 1
  catch
    println("Error on bootstrap iteration: $(b)")
  end
end

CSV.write(
  "estimates/short_t_qld_p_$(retail_p)_est_retail_B_$(B).csv",
  DataFrame(retail_ests, ["tau$(i)" for i in rel_years]),
)
end_time = Dates.now()
print("Current Time: $(Dates.format(end_time, "HH:MM:SS"))")

# %% Wholesale bootstrap (takes about 10 minutes)

# display current time
start_time = Dates.now()
println("Current Time: $(Dates.format(start_time, "HH:MM:SS"))")
B = 1000
ids = unique(df[:, :fips])
wholesale_ests = zeros(B, length(rel_years))
wholesale_ests[1, :] = wholesale_est.est[1]'

b = 2
while b < B
  b % 50 == 0 && println("On Bootstrap Iteration: $(b)")

  resample_ids = StatsBase.sample(ids, length(ids))
  resample_df = vcat([df[df.fips .== id, :] for id in resample_ids]...)
  resample_df.fips = repeat(1:length(ids); inner=T)

  try
    boot_est = est_te(
      resample_df,
      wholesale_p,
      :log_wholesale_emp;
      export_y0=false,
      outcome="wholesale",
      quick_plot=false,
    )
    wholesale_ests[b, :] = boot_est.est[1]'
    b += 1
  catch
    println("Error on bootstrap iteration: $(b)")
  end
end

CSV.write(
  "estimates/short_t_qld_p_$(wholesale_p)_est_wholesale_B_$(B).csv",
  DataFrame(wholesale_ests, ["tau$(i)" for i in rel_years]),
)
end_time = Dates.now()
print("Current Time: $(Dates.format(end_time, "HH:MM:SS"))")
