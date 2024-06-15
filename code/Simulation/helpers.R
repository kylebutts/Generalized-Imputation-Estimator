# %%
extract_tt_latex_body <- function(x) {
  paste0(
    tinytable:::build_tt(x, "latex")@body,
    collapse = "\n"
  )
}

# %%
each_row <- function(data) {
  split(data, seq_len(nrow(data)))
}

run_simulation <- function(B, dgps, estimators, cores = 8, seed = NULL) {
  plan(multicore, workers = cores)
  ests <- future_map(
    each_row(dgps),
    function(dgp) {
      df_base <- gen_df_base(N = dgp$N, T0 = dgp$T0, T = dgp$T)
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

      list_rbind(map(1:B, function(b) {
        df_b <- gen_b(
          df = df_base, T0 = dgp$T0, F = F,
          twfe = dgp$twfe, parallel_trends = dgp$parallel_trends,
          ar_error_term = dgp$ar_error_term, instrument_noise = dgp$instrument_noise
        )

        list_rbind(map(each_row(estimators), function(estimator) {
          est_function <- estimator$est_function[[1]]
          est <- quietly({
            est_function(df_b)
          })

          tibble(
            dgp_num = dgp$dgp_num,
            estimator = estimator$estimator,
            rel_year = rel_year,
            true_te = true_te,
            est = est
          )
        }))
      }))
    },
    .options = furrr_options(seed = seed)
  )
  plan(multicore, workers = cores)

  list_rbind(ests)
}

# %%
summarize_ests <- function(ests, dgps, dgp_num = 1) {
  summary <- ests |>
    filter(dgp_num == .env$dgp_num) |>
    summarize(
      .by = c(dgp_num, estimator, rel_year),
      mean_estimate = mean(est),
      mean_bias = mean(est - true_te),
      mse = mean((est - true_te)^2)
    ) |>
    pivot_wider(
      id_cols = c(dgp_num, estimator),
      names_from = rel_year,
      names_glue = "{.value}_{rel_year}",
      values_from = c(mean_estimate, mean_bias, mse)
    )

  rel_years <- ests |>
    filter(dgp_num == .env$dgp_num) |>
    with(unique(rel_year))
  cols <- c("estimator")
  new_col_names <- c("Estimator")
  groups <- list()
  for (i in seq_along(rel_years)) {
    rel_year <- rel_years[i]
    cols <- c(cols, paste0("mean_bias_", rel_year), paste0("mse_", rel_year))

    groups[[sprintf("$\\\\tau^{%s}$", rel_year)]] <- (2 * i):(2 * i + 1)

    new_col_names <- c(new_col_names, "Mean Bias", "MSE")
    # new_col_names = c(new_col_names,
    #   sprintf("Mean Bias $\\\\tau_{%s}$", rel_year),
    #   sprintf("MSE $\\\\tau_{%s}$", rel_year)
    # )
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
    format_tt(j = "(Mean|MSE)", sprintf = "%0.2f")
}

# %%
summarize_ests_wide <- function(ests, which_rel_year = 0, latex_body = FALSE) {
  summary <- ests |>
    summarize(
      .by = c(dgp_num, estimator, estimator_short,, rel_year),
      mean_estimate = mean(est),
      mean_bias = mean(est - true_te),
      mse = mean((est - true_te)^2)
    ) |>
    filter(rel_year == .env$which_rel_year) |>
    pivot_wider(
      id_cols = c(estimator),
      names_from = dgp_num,
      names_glue = "{.value}_{dgp_num}",
      values_from = c(mean_bias, mse)
    )

  dgps <- ests$dgp_num |> unique()
  cols <- c("estimator")
  new_col_names <- c("Estimator")
  groups <- list()
  for (i in seq_along(dgps)) {
    dgp <- dgps[[i]]
    cols <- c(cols, paste0("mean_bias_", dgp), paste0("mse_", dgp))
    new_col_names <- c(new_col_names, "Mean Bias", "MSE")
    groups[[paste0("DGP ", dgps[[i]])]] <- (2 * i):(2 * i + 1)
  }
  summary <- summary[, cols]
  colnames(summary) <- new_col_names

  table <- summary |>
    tt() |>
    format_tt(j = unlist(groups), sprintf = "%0.2f") |>
    group_tt(j = groups)

  return(table)
}
