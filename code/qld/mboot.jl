
"""
    mboot(IF::AbstractMatrix{<:Real}, B::Integer=1000; level::Real=0.95)

Performs a multiplier bootstrap using Rademacher weights.

Calculates bootstrap standard errors and the critical value for constructing
uniform confidence bands based on the maximum absolute t-statistic.

# Arguments
- `IF::AbstractMatrix{<:Real}`: N x K matrix of influence functions, where N is
  the number of observations and K is the number of parameters.
- `B::Integer=1000`: The number of bootstrap iterations.
- `level::Real=0.95`: The confidence level for the critical value calculation.

# Returns
- `NamedTuple{(:se, :crit_val)}`: A named tuple containing:
    - `se`: A vector of bootstrap standard errors (K x 1).
    - `crit_val`: The critical value for the uniform confidence band.

# Details
Follows the approach similar to the `mboot.R` function from the `did` R package,
using an IQR-based standard error calculation.

# Examples
```julia
# Example usage:
N = 100
K = 5
IF = randn(N, K) # Example influence function matrix
result = mboot(IF, 1000)
println("Bootstrap SEs: ", result.se)
println("Critical Value (95%): ", result.crit_val)

# Use BenchmarkTools for performance check
# using BenchmarkTools
# @benchmark mboot(IF, 1000)
"""
function mboot(IF::AbstractMatrix{<:Real}, B::Integer=1000; level::Real=0.95)
  N, K = size(IF)

  # Pre-allocate matrix for bootstrap estimates (B x K)
  bs_ests = Matrix{Float64}(undef, B, K)

  # Multiplier bootstrap loop
  for b in 1:B
    # Generate N Rademacher weights (+1 or -1)
    rad = rand([-1, 1], N)

    # Calculate bootstrap estimate for iteration b: (K x N) * (N x 1) -> K x 1
    bs_ests[b, :] = transpose(transpose(IF) * rad)
  end

  # Calculate Bootstrap standard errors using IQR method
  se_bs = Vector{Float64}(undef, K)
  qnorm_diff =
    quantile(Distributions.Normal(), 0.75) - quantile(Distributions.Normal(), 0.25)
  for k in 1:K
    # Calculate quantiles, handling potential NaNs (though unlikely here)
    q75 = quantile(bs_ests[:, k], 0.75)
    q25 = quantile(bs_ests[:, k], 0.25)
    se_bs[k] = (q75 - q25) / qnorm_diff
  end

  # Calculate absolute t-statistics for each draw
  max_abs_t_stats_per_draw = Vector{Float64}(undef, B)
  for b in 1:B
    max_abs_t_stats_per_draw[b] = maximum(abs.(bs_ests[b, :] ./ se_bs))
  end

  # Calculate the critical value (quantile of the max t-stats)
  crit_val = quantile(max_abs_t_stats_per_draw, level)

  return (se=se_bs, crit_val=crit_val)
end
