# %%
extract_tt_latex_body <- function(x) {
  paste0(
    tinytable:::build_tt(x, "latex")@body,
    collapse = "\n"
  )
}

# %%
each_row <- function(data) {
  lapply(
    split(data, seq_len(nrow(data))),
    as.list
  )
}

run_simulation <- function(B, dgps, estimators, do_inference = FALSE, cores = 8, seed = NULL) {
  # plan(multicore, workers = cores)
  set.seed(seed)
  # ests <- future_map(
  ests <- map(
    each_row(dgps),
    function(dgp) {
      df_base <- gen_df_base(N = dgp$N, T0 = dgp$T0, T = dgp$T, twfe = dgp$twfe)
      
      F <- df_base |>
        filter(id == 1) |>
        arrange(t) |>
        select(F1, F2) |>
        as.matrix()

      true_te <- df_base |>
        filter(post == TRUE) |>
        slice(1, .by = t) |>
        arrange(t) |>
        pull(true_te)

      rel_year <- (1:length(true_te)) - 1

      cat("Simulation progress (`.` = 10): \n")
      res <- list_rbind(map(1:B, function(b) {
        if (b %% 10 == 0) {
          if (b %% 200 == 0) {
            cat(".\n")
          } else {
            cat(".")
          }
        }

        df_b <- gen_b(
          df = df_base, T0 = dgp$T0, F = F,
          twfe = dgp$twfe, parallel_trends = dgp$parallel_trends,
          ar_error_term = dgp$ar_error_term, instrument_noise = dgp$instrument_noise
        )

        estimates <- list_rbind(map(
          each_row(estimators), 
          function(estimator) {
            est_function <- estimator$est_function[[1]]

            if (estimator$estimator_short %in% c("synth", "augsynth")) { 
              res <- purrr::possibly(
                est_function, otherwise = NULL, quiet = TRUE
              )(df_b, do_inference)
            } else if (estimator$estimator_short == "qld") {
              res <- purrr::possibly(
                est_function, otherwise = NULL, quiet = TRUE
              )(df_b)

              # if (res$selected_p[1] == 0) {
              #   message("FOUND ONE")
              #   write_csv(df_b, here("data/Simulations/df_selected_p_0.csv"))
              # }
              
            } else {
              res <- purrr::possibly(
                est_function, otherwise = NULL, quiet = TRUE
              )(df_b)
            }
            if (is.null(res)) {
              warning(sprintf("On DGP %s and iteration %s, %s failed", dgp$dgp_num, b, estimator$estimator_short))
              return(res)
            }
            
            res <- res |> mutate(
              estimator = estimator$estimator, 
              true_te = .env$true_te, 
              .before = 1
            ) 
            return(res)
          }
        ))


        estimates <- estimates |>
          mutate(
            b = .env$b, 
            dgp_num = dgp$dgp_num, 
            T_0 = dgp$T0,
            .before = 1
          )
      }))
      cat("\n")

      return(res)
    }
    # },
    # .options = furrr_options(seed = seed)
  )

  list_rbind(ests)
}

# %% 
summarize_ests <- function(ests, dgps, dgp_num = 1, T_0 = 4) {
  summary <- ests |>
    filter(dgp_num == .env$dgp_num, T_0 == .env$T_0) |>
    summarize(
      .by = c(dgp_num, estimator, rel_year),
      mean_estimate = mean(estimate),
      mean_bias = mean(estimate - true_te),
      mse = mean((estimate - true_te)^2),
      rmse = sqrt(mean((estimate - true_te)^2)),
      coverage = mean(true_te >= ci_lower & true_te <= ci_upper)
    ) |>
    pivot_wider(
      id_cols = c(dgp_num, estimator),
      names_from = rel_year,
      names_glue = "{.value}_{rel_year}",
      values_from = c(mean_estimate, mean_bias, rmse, coverage)
    )

  rel_years <- ests |>
    filter(dgp_num == .env$dgp_num) |>
    with(unique(rel_year))
  cols <- c("estimator")
  new_col_names <- c("Estimator")
  groups <- list()
  for (i in seq_along(rel_years)) {
    rel_year <- rel_years[i]
    cols <- c(cols, paste0("mean_bias_", rel_year), paste0("rmse_", rel_year), paste0("coverage_", rel_year))

    groups[[sprintf("$\\\\tau^{%s}$", rel_year)]] <- 1 + 3*(i-1) + (1:3)

    new_col_names <- c(new_col_names, "Bias", "RMSE", "Coverage")
  }

  summary <- summary[, cols]
  colnames(summary) <- new_col_names

  caption <- sprintf("DGP %s", dgp_num)
  dgp_details <- dgps |> filter(dgp_num == .env$dgp_num)
  notes <- with(dgp_details, sprintf(
    "twfe = %s; parallel_trends = %s; ar_error_term = %s",
    twfe, parallel_trends, ar_error_term
  ))

  summary |>
    tt(caption = caption, notes = notes) |>
    group_tt(j = groups) |>
    format_tt(j = "(Bias|RMSE)", sprintf = "%0.2f") |>
    format_tt(j = "(Coverage)", sprintf = "%0.3f")
}

# %%
summarize_ests_wide <- function(ests, which_rel_year = 0, latex_body = FALSE) {
  summary <- ests |>
    summarize(
      .by = c(dgp_num, estimator, rel_year),
      mean_estimate = mean(estimate),
      bias = mean(estimate - true_te),
      mse = mean((estimate - true_te)^2),
      rmse = sqrt(mean((estimate - true_te)^2)),
      coverage = 100 * mean(true_te >= ci_lower & true_te <= ci_upper)
    ) |>
    filter(rel_year == .env$which_rel_year) |>
    pivot_wider(
      id_cols = c(estimator),
      names_from = dgp_num,
      names_glue = "{.value}_{dgp_num}",
      values_from = c(bias, rmse, coverage)
    )

  dgps <- ests$dgp_num |> unique()
  cols <- c("estimator")
  new_col_names <- c("Estimator")
  groups <- list()
  for (i in seq_along(dgps)) {
    dgp <- dgps[[i]]
    cols <- c(cols, paste0("bias_", dgp), paste0("rmse_", dgp), paste0("coverage_", dgp))
    new_col_names <- c(new_col_names, "Bias", "RMSE", "Coverage")
    groups[[paste0("DGP ", dgps[[i]])]] <-  1 + 3*(i-1) + (1:3)
  }
  summary <- summary[, cols]
  colnames(summary) <- new_col_names

  table <- summary |>
    tt() |>
    format_tt(j = ("(Bias|RMSE)"), sprintf = "%0.2f") |>
    format_tt(j = ("(Coverage)"), sprintf = "%0.1f\\%%") |>
    group_tt(j = groups)

  return(table)
}
