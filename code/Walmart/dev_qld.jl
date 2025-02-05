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
# include("$here/code/Walmart/dev_qld_helpers.jl")

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

include("../qld/QLD.jl")
include("../qld/within_transform.jl")
include("../qld/attgt.jl")
include("../qld/gmm_qld.jl")
include("../qld/qld_imputation.jl")


ret_retail = qld_imputation(
  df;
  y=:log_retail_emp,
  id=:fips,
  t=:year,
  g=:g,
  W=[
    :share_pop_ind_manuf,
    :share_pop_poverty_78_below,
    :share_pop_poverty_78_above,
    :share_pop_emp_private,
    :share_pop_emp_government,
    :share_school_col,
    :share_school_hs,
  ],
  do_within_transform=true,
  p=Int(-1),
  type="dynamic",
  return_y0=true,
  return_naive_se=true,
)

plot(ret_retail[1], ret_retail[2]; ribbon=1.96 .* sqrt.(diag(ret_retail[3])), fillalpha=0.3)



