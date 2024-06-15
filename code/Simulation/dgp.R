library(tidyverse)

# %%
#' @param n Integer number of observations to draw
#' @param rho Vector of autocorrelation coefficients
rAR <- function(n, rho = 0.75) {
  as.numeric(arima.sim(list(ar = rho), n = n, rand.gen = rnorm))
}

# %%
gen_df_base <- function(N = 200, T0 = 4, T = 7) {
  # True treatment effects
  true_te <- c(rep(0, T0), 1:(T - T0))
  post <- c(rep(FALSE, T0), rep(TRUE, T - T0))

  # time FEs
  eta <- rAR(T, c(0.75))

  # Common factors
  f1 <- 1:T
  f2 <- rAR(n = T, rho = 0.5)
  F <- cbind(f1, f2)

  time_effects <- tibble(t = 1:T, post = post, true_te = true_te, eta = eta)
  for (i in seq_len(ncol(F))) {
    time_effects[[paste0("F", i)]] <- F[, i]
  }

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
      # make between 0/1
      p <- (g1 - min(g1)) / (max(g1) - min(g1))
      p <- 0.5 / mean(p) * p
    }
    return(p)
  }

  ids <- unique(df$id)
  n_ids <- length(ids)
  unit_info <- tibble(
      id = ids,
      mu = rnorm(n_ids, sd = 1),
      g1 = rnorm(n_ids, mean = mu, sd = 1),
      g2 = rnorm(n_ids, mean = mu, sd = 1),
      W1 = g1 + rnorm(n_ids, sd = sqrt(.env$instrument_noise)),
      W2 = g2 + rnorm(n_ids, sd = sqrt(.env$instrument_noise)),
      prob_treat = gen_treat_probs(g1 = g1, parallel_trends = .env$parallel_trends),
      treated = runif(n_ids) <= prob_treat
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
    df <- mutate(df, y = treated * true_te + mu + eta + eps)
  } else {
    df <- mutate(df, y = treated * true_te + factor_prod + eps)
  }

  return(df)
}
