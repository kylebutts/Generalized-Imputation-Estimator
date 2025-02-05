# %%
each_row <- function(data) {
  unname(lapply(
    split(data, seq_len(nrow(data))),
    as.list
  ))
}

# %% 
run_simulation <- function(B, dgps, estimators, seed = NULL) {
  set.seed(seed)
  max_T0 <- 12
  time_effects <- gen_time_base(T0 = max_T0)

  ests <- imap(
    each_row(dgps),
    function(dgp, i_dgp) {
      curr_time_effects <- time_effects |>
        slice_tail(n = dgp$T0 + 1)

      true_te <- curr_time_effects |>
        filter(post == TRUE) |>
        slice(1, .by = t) |>
        arrange(t) |>
        pull(true_te)
      rel_year <- (1:length(true_te)) - 1

      cat(sprintf("Simulation %02i (`.` = 10): \n", i_dgp))
      res <- list_rbind(map(1:B, function(b) {
        if (b %% 10 == 0) {
          if (b %% 250 == 0) {
            cat(".\n")
          } else {
            cat(".")
          }
        }

        df_b <- gen_b(
          N = dgp$N,
          time_effects = curr_time_effects,
          twfe = dgp$twfe,
          parallel_trends = dgp$parallel_trends,
          ar_error_term = dgp$ar_error_term,
          instrument_noise = dgp$instrument_noise
        )

        estimates <- list_rbind(map(
          each_row(estimators),
          function(estimator) {
            est_function <- estimator$est_function[[1]]
            res <- purrr::possibly(
              est_function,
              otherwise = NULL, quiet = TRUE
            )(df_b)
            
            if (is.null(res)) {
              warning(sprintf("On DGP %s and iteration %s, %s failed", dgp$dgp_num, b, estimator$estimator_short))
              return(res)
            }

            res <- res |> 
              mutate(estimator = estimator$estimator, .before = 1) |>
              mutate(true_te = .env$true_te, .before = estimate)
            return(res)
          }
        ))
        
        # record DGP information
        estimates <- estimates |>
          mutate(
            b = .env$b,
            dgps[i_dgp, ],
            .before = 1
          )
      }))
      cat("\n")

      return(res)
    }
  )

  list_rbind(ests)
}

