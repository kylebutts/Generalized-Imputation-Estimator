library(tidyverse)

# %%
#' Draw eps_t = rho eps_{t-1} + v_t
#'
#' @param T Integer number of observations to draw
#' @param sd Standard deviation of AR(1) process
#' @param mean Mean of AR(1) process
#' @param rho Vector of autocorrelation coefficients
rAR <- function(T, mean = 0, sd = 1, rho = 0.75) {
  x <- arima.sim(
    list(ar = rho),
    n = T,
    rand.gen = function(n) rnorm(n, mean = 0, sd = sd)
  )
  # Make variance = var_eps
  x <- as.numeric(x) * sqrt(1 - rho^2)
  # Make mean = mean
  x <- x + mean
  return(x)
}

# rAR <- function(T, mean = 0, sd = 1, rho = 0.75) {
#   x <- arima.sim(
#     list(ar = rho),
#     n = T,
#     # Make variance = var_eps
#     rand.gen = function(n) rnorm(n, mean = 0, sd = sd * sqrt(1 - rho^2))
#   )
#   x <- as.numeric(x) + mean
#   return(x)
# }

# %%
gen_time_base <- function(T0 = 12, T = T0 + 1) {
  # True treatment effects
  true_te <- c(rep(0, T0), 1:(T - T0))
  post <- c(rep(FALSE, T0), rep(TRUE, T - T0))

  # Common factors
  F1 <- 1:T

  # twfe_F1 <- rep(1, T)
  # unit-fixed effects

  F2 <- rAR(T = T, mean = 0, sd = 0.5, rho = 0.75)

  time_effects <- tibble(
    t = 1:T, post = post, true_te = true_te,
    F1 = F1, F2 = F2
  )
  return(time_effects)
}

# %%
#' Generate treatment, increasing probability in g2 if `parallel_trends == FALSE`
gen_treat_probs <- function(g1, parallel_trends = TRUE) {
  if (parallel_trends == TRUE) {
    p <- as.numeric(
      rep(0.5, length(g1))
    )
  } else {
    p <- as.numeric(g1 > 1.0)
  }
  return(p)
}

# This assumes df is sorted by (id, t) before being called
gen_b <- function(N = 200, time_effects, twfe = FALSE, parallel_trends = FALSE, ar_error_term = FALSE, instrument_noise = 1) {
  
  unit_info <- tibble(id = 1:N) |>
    mutate(
      # Factor loadings
      mu = rnorm(n(), mean = 1, sd = 1),
      g1 = mu,
      # if \gamma_2 = 1, then these are time fixed-effects
      g2 = if (.env$twfe == TRUE) {
        1
      } else {
        rnorm(n(), mean = 1, sd = 1)
      },

      # Time-invariant covariate proxies for the factor loadings
      W1 = g1 + rnorm(n(), sd = sqrt(.env$instrument_noise)),
      W2 = g2 + rnorm(n(), sd = sqrt(.env$instrument_noise)),

      # Treatment
      prob_treat = gen_treat_probs(g1 = g1, parallel_trends = .env$parallel_trends),
      treated = runif(n()) <= prob_treat,
    )

  F1 = if (twfe == TRUE) {
    rep(1, nrow(time_effects))
  } else {
    time_effects$F1 
  }
  F2 = time_effects$F2
  F = cbind(F1, F2)
  
  # Expand out over T0
  df <- unit_info |>
    rowwise(everything()) |>
    reframe(
      t = time_effects$t,
      post = time_effects$post,
      true_te = time_effects$true_te,
      F1 = F[, 1],
      F2 = F[, 2],
      factor_prod = as.numeric(F %*% rbind(g1, g2))
    ) |>
    mutate(
      .by = id,
      g = if_else(
        treated[1] == TRUE, min(t[post == TRUE]), Inf
      )
    ) |>
    mutate(
      treat = treated * post,
      rel_year = ifelse(treated, t - g, -10)
    ) |>
    select(
      id, t, everything()
    )

  if (ar_error_term == FALSE) {
    df <- df |> mutate(eps = rnorm(n(), mean = 0, sd = 1))
  } else {
    df <- df |> mutate(.by = id, eps = rAR(n(), mean = 0, sd = 1, rho = 0.75))
  }
  df <- df |> mutate(y = treat * true_te + factor_prod + eps)

  return(df)
}
