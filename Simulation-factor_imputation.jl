
# Note: this is not general, it hard codes some things from the simulation
# A more general implemenation for staggered timing is in Walmart application
function est_factor_imputation(data)
  p = 1
  
  # Remove unit and time fixed effects as detailed in the paper
  iota_Pre = [ones(T0, 1); zeros(T - T0, 1)]
  iota = ones(T, 1)
  Imp = iota * inv(iota_Pre' * iota_Pre) * iota_Pre'
  M = I(T) - Imp

  # untreated cross-sectional averages of y
  y0bar = data |>
    x -> x[(x.ever_treated .== 0), :] |>
    x -> groupby(x, :t) |>
    x -> combine(x,
      :y => (mean) => :y
    ) |>
    x -> sort(x, :t) |>
    x -> x[:, :y]

  # Create y1tilde
  data |>
    x -> groupby(x, :id) |>
    x -> transform!(x,
      :y => (y -> vec(M * (y - y0bar))) => :y1tilde
    )

  untreated = data |> 
    x -> x[(x.ever_treated .== 0), :]


  # GMM estimate factor
  # Θ = ones(T - p, p)
  function moment(Θ) 
    HΘ = [I(T - p); Θ']


    transformed = untreated |>
      x -> groupby(x, :id) |>
      x -> combine(x,
        [:y1tilde, :W] => 
          ((y, w) -> (HΘ' * y) * w[1]) => 
          :Hy1tildeW,
      ) 
    
    return transformed[:, :Hy1tildeW]' * transformed[:, :Hy1tildeW] 
  end

  optres = optimize(moment, ones(T - p, p), LBFGS(); autodiff = :forward)
  Θhat = optres.minimizer

  FΘ = [Θhat; (-1 * I(p))]
  FΘ = Matrix(FΘ)
  FΘ_pre = FΘ[1:T0, :]
  FΘ_post = FΘ[(T0 + 1):end, :]

  treated_ids = data |> 
    x -> x[(x.ever_treated .== 1), :id] |>
    unique 

  data.y0hat = zeros(size(data, 1))
  for id in treated_ids
    y_pre  = data[(data.id .== id) .& (data.t .<= T0), :y1tilde]  

    data[(data.id .== id), :y0hat] = 
      FΘ * inv(FΘ_pre' * FΘ_pre) * FΘ_pre' * y_pre
  end

  taus = data |> 
    x -> x[(x.ever_treated .== 1), :] |>
    x -> groupby(x, :t) |>
    x -> combine(x, 
      [:y1tilde, :y0hat] => 
        ((y, y0hat) -> mean(y - y0hat)) => 
        :tau
    ) |>
    x -> x[x.t .> T0, :tau]

  return taus
end
