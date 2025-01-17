library(tidyverse)

# %%
#' @param n Integer number of observations to draw
#' @param rho Vector of autocorrelation coefficients
rAR <- function(n, rho = 0.75) {
  as.numeric(arima.sim(list(ar = rho), n = n, rand.gen = rnorm))
}

# %%
gen_df_base <- function(N = 200, T0 = 4, T = 7, twfe = FALSE) {
  # True treatment effects
  true_te <- c(rep(0, T0), 1:(T - T0))
  post <- c(rep(FALSE, T0), rep(TRUE, T - T0))

  # Common factors
  if (twfe == FALSE) {
    F1 <- 1:T
  } else if (twfe == TRUE) {
    # unit-fixed effects
    F1 <- rep(1, T)
  }
  F2 <- rAR(n = T, rho = 0.75)
  F <- cbind(F1, F2)

  time_effects <- cbind(
    tibble(t = 1:T, post = post, true_te = true_te, eta = F2),
    F
  )

  # F and eta are non-random across simulations
  df <- expand_grid(id = 1:N, t = 1:T) |>
    left_join(time_effects, by = "t")

  df <- df |> arrange(id, t)

  return(df)
}

# %%
# This assumes df is sorted by (id, t) before being called
gen_b <- function(df, T0, F, twfe = FALSE, parallel_trends = FALSE, ar_error_term = FALSE, instrument_noise = 1) {
  start <- Sys.time()

  #' Generate treatment, increasing probability in g2 if `parallel_trends == FALSE`
  gen_treat_probs <- function(g1, parallel_trends = TRUE) {
    if (parallel_trends == TRUE) {
      p <- as.numeric(
        rep(0.5, length(g1))
      )
    } else {
      # p <- (g1 - min(g1)) / (max(g1) - min(g1))
      # p <- 0.5 / mean(p) * p
      p <- as.numeric(g1 > 0)
    }
    return(p)
  }

  ids <- unique(df$id)
  n_ids <- length(ids)
  unit_info <- tibble(
    id = ids,
    # unit fixed effects
    mu = rnorm(n_ids, sd = 1),
  )
  unit_info <- unit_info |>
    mutate(
      g1 = if (.env$twfe == TRUE) {
        rnorm(n(), mean = mu, sd = 1)
      } else {
        rnorm(n(), mean = mu, sd = 1)
      },
      # if \gamma_2 = 1, then these are time fixed-effects
      g2 = if (.env$twfe == TRUE) {
        1
      } else {
        rnorm(n(), mean = mu, sd = 1)
      },
      W1 = g1 + rnorm(n_ids, sd = sqrt(.env$instrument_noise)),
      W2 = g2 + rnorm(n_ids, sd = sqrt(.env$instrument_noise)),
      prob_treat = gen_treat_probs(g1 = g1, parallel_trends = .env$parallel_trends),
      treated = runif(n_ids) <= prob_treat,
      g = ifelse(treated == TRUE, (T0 + 1), Inf)
    )

  unit_info <- unit_info |>
    rowwise(everything()) |>
    reframe(
      t = 1:nrow(F),
      factor_prod = as.numeric(F %*% rbind(g1, g2))
    )

  df <- df |>
    left_join(unit_info, by = join_by(id, t)) |>
    mutate(
      treat = as.numeric(treated & post),
      rel_year = ifelse(treated, t - (T0 + 1), -10)
    )

  if (ar_error_term == FALSE) {
    df <- df |> mutate(eps = rnorm(n(), 0, 1))
  } else {
    df <- df |> mutate(eps = rAR(n(), 0.75), .by = id)
  }

  if (twfe == TRUE) {
    df <- df |>
      mutate(
        y = treated * true_te + mu + eta + eps,
        F1 = 1,
        F2 = eta
      )
  } else {
    df <- mutate(df, y = treated * true_te + factor_prod + eps)
  }

  return(df)
}
