# %%
library(tidyverse)
library(here)
library(fs)
library(kfbmisc)
library(tinytable)

fs::dir_create(here("data/Simulations"))
fs::dir_create(here("out/tables/simulation-1/"))
source(here("code/Simulation/dgp.R"))
source(here("code/Simulation/estimators.R"))
source(here("code/Simulation/run_simulation.R"))
source(here("code/helpers/extract_tt_latex_body.R"))

RUN_SIMULATION <- FALSE

estimators <- tibble(
  est_function = rlang::list2(
    function(df) {
      est_dim(df)
    },
    function(df) {
      est_twfe(df)
    },
    function(df) {
      est_twfe_covs(df)
    },
    function(df) {
      est_augsynth(df)
    },
    function(df) {
      est_gsynth(df, force = "none", p = 2L)
    },
    function(df) {
      est_gsynth(df, force = "none", p = c(0L, 3L))
    },
    # function(df) {
    #   est_gsynth(df, force = "none", p = 0L)
    # },
    # function(df) {
    #   est_gsynth(df, force = "none", p = 1L)
    # },
    function(df) {
      est_qld_F_known(df, do_within_transform = FALSE)
    },
    function(df) {
      est_qld(df, do_within_transform = FALSE, p = 2)
    },
    function(df) {
      est_qld(df, do_within_transform = FALSE, p = -1)
    }
  ),
  estimator = c(
    "Difference-in-means",
    "TWFE",
    "TWFE with $\\bm{w}_i' \\bm{\\beta}_t$",
    "Augmented Synthetic Control",
    "Generalized Synth ($p$ known)",
    "Generalized Synth ($p$ estimated)",
    # "Generalized Synth ($p = 0$)",
    # "Generalized Synth ($p = 1$)",
    "Factor Imputation ($\\bm{F}$ known)",
    "QLD ($p$ known)",
    "QLD ($p$ estimated)"
  ),
  estimator_short = c(
    "dim",
    "twfe",
    "twfe_covs",
    "augsynth",
    "gsynth_p_known",
    "gsynth",
    # "gsynth_p_0",
    # "gsynth_p_1",
    "qld_F_known",
    "qld_p_known",
    "qld"
  )
)

# fmt: skip
dgps <- tribble(
  ~dgp_num, ~N, ~T0, ~twfe, ~parallel_trends, ~ar_error_term, ~instrument_noise,
  01, 300L, 4L, TRUE, TRUE, FALSE, 1,
  02, 300L, 4L, FALSE, TRUE, FALSE, 1,
  03, 300L, 4L, FALSE, FALSE, FALSE, 1,
  04, 300L, 4L, FALSE, FALSE, TRUE, 1,
  01, 300L, 12L, TRUE, TRUE, FALSE, 1,
  02, 300L, 12L, FALSE, TRUE, FALSE, 1,
  03, 300L, 12L, FALSE, FALSE, FALSE, 1,
  04, 300L, 12L, FALSE, FALSE, TRUE, 1,
)

B <- 2000
# B <- 100

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run simulation ----
if (RUN_SIMULATION == TRUE) {
  tictoc::tic()
  ests <- run_simulation(B, dgps, estimators, seed = 20250121)
  tictoc::toc()

  write_csv(ests, here("data/Simulations/simulation_1_ests.csv"))
}


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Report on simulation ----
ests <- read_csv(
  here("data/Simulations/simulation_1_ests.csv"),
  show_col_types = FALSE
)

# %%
# quick detailed sum
ests |>
  summarize(
    .by = c(dgp_num, T0, estimator),
    n = n(),
    bias = mean(estimate - true_te, na.rm = TRUE),
    rmse = sqrt(mean((estimate - true_te)^2)),
    coverage = mean(true_te >= ci_lower & true_te <= ci_upper),
    avg_std_error = mean(std.error, na.rm = TRUE),
    avg_ci_length = mean(ci_upper - ci_lower, na.rm = TRUE),
    avg_est = mean(estimate, na.rm = TRUE),
    est_05th = quantile(estimate, 0.05),
    est_50th = quantile(estimate, 0.50),
    est_95th = quantile(estimate, 0.95),
    pct_p_0 = mean(selected_p == 0, na.rm = TRUE),
    pct_p_1 = mean(selected_p == 1, na.rm = TRUE),
    pct_p_2 = mean(selected_p == 2, na.rm = TRUE),
    pct_p_3 = mean(selected_p == 3, na.rm = TRUE),
    mean_p = mean(selected_p, na.rm = TRUE)
  ) |>
  print(n = 100)

ests |>
  summarize(
    .by = c(dgp_num, T0, estimator),
    n = n(),
    bias = mean(estimate - true_te, na.rm = TRUE),
    rmse = sqrt(mean((estimate - true_te)^2)),
    coverage = mean(true_te >= ci_lower & true_te <= ci_upper),
    avg_std_error = mean(std.error, na.rm = TRUE),
    avg_ci_length = mean(ci_upper - ci_lower, na.rm = TRUE),
    avg_est = mean(estimate, na.rm = TRUE),
    est_05th = quantile(estimate, 0.05),
    est_50th = quantile(estimate, 0.50),
    est_95th = quantile(estimate, 0.95),
    pct_p_0 = mean(selected_p == 0, na.rm = TRUE),
    pct_p_1 = mean(selected_p == 1, na.rm = TRUE),
    pct_p_2 = mean(selected_p == 2, na.rm = TRUE),
    pct_p_3 = mean(selected_p == 3, na.rm = TRUE),
    mean_p = mean(selected_p, na.rm = TRUE)
  ) |>
  mutate(estimator = str_replace(estimator, "Generalized Synth", "gsynth")) |>
  filter(str_detect(estimator, "gsynth")) |>
  left_join(dgps |> select(dgp_num, N), by = "dgp_num") |>
  # select(-dgp_num) |>
  select(dgp_num, T0, N, estimator, B = n, everything()) |>
  print(n = 100)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Export ----
summarize_ests_wide <- function(ests, which_rel_year = 0) {
  summary <- ests |>
    filter(rel_year == .env$which_rel_year) |>
    summarize(
      .by = c(dgp_num, estimator),
      mean_estimate = mean(estimate),
      bias = mean(estimate - true_te),
      mse = mean((estimate - true_te)^2),
      rmse = sqrt(mean((estimate - true_te)^2))
    ) |>
    pivot_wider(
      id_cols = c(estimator),
      names_from = dgp_num,
      names_glue = "{.value}_{dgp_num}",
      values_from = c(bias, rmse)
    )

  dgp_nums <- ests$dgp_num |> unique()
  cols <- c("estimator")
  new_col_names <- c("Estimator")
  groups <- list()
  for (i in seq_along(dgp_nums)) {
    dgp <- dgp_nums[[i]]
    cols <- c(
      cols,
      paste0("bias_", dgp),
      paste0("rmse_", dgp)
    )
    new_col_names <- c(new_col_names, "Bias", "RMSE")
    groups[[paste0("DGP ", dgp_nums[[i]])]] <- 1 + 2 * (i - 1) + (1:2)
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

summarize_ests_wide_ci <- function(ests, which_rel_year = 0) {
  summary <- ests |>
    filter(rel_year == .env$which_rel_year) |>
    summarize(
      .by = c(dgp_num, estimator),
      bias = mean(estimate - true_te),
      rmse = sqrt(mean((estimate - true_te)^2)),
      coverage = mean(true_te >= ci_lower & true_te <= ci_upper),
      average_ci_length = mean(ci_upper - ci_lower),
      sd_ci_length = sd(ci_upper - ci_lower)
    ) |>
    pivot_wider(
      id_cols = c(estimator),
      names_from = dgp_num,
      names_glue = "{.value}_{dgp_num}",
      values_from = c(bias, rmse, coverage, average_ci_length, sd_ci_length)
    )

  dgp_nums <- ests$dgp_num |> unique()
  cols <- c("estimator")
  new_col_names <- c("Estimator")
  groups <- list()
  for (i in seq_along(dgp_nums)) {
    dgp <- dgp_nums[[i]]
    cols <- c(
      cols,
      paste0("bias_", dgp),
      paste0("rmse_", dgp),
      paste0("coverage_", dgp),
      paste0("average_ci_length_", dgp),
      paste0("sd_ci_length_", dgp)
    )

    groups[[paste0("DGP ", dgp_nums[[i]])]] <- 1 + 5 * (i - 1) + (1:5)

    new_col_names <- c(
      new_col_names,
      "Bias",
      "RMSE",
      "Coverage",
      "Avg. CI Length",
      "SD CI Length"
    )
  }

  summary <- summary[, cols]
  colnames(summary) <- new_col_names

  caption <- sprintf("DGP %s; T_0 = %s", dgp_num, T0)
  dgp_details <- dgps |> dplyr::filter(dgp_num == .env$dgp_num, T0 == .env$T0)
  notes <- with(
    dgp_details,
    sprintf(
      "twfe = %s; parallel_trends = %s; ar_error_term = %s",
      twfe,
      parallel_trends,
      ar_error_term
    )
  )

  summary |>
    tt(caption = caption, notes = notes) |>
    group_tt(j = groups) |>
    format_tt(j = "(Bias|RMSE|Coverage|CI Length)", sprintf = "%0.2f") |>
    format_tt(j = "(Coverage)", sprintf = "%0.3f")
}

# %%
cat("\n\n\n")
ests |>
  filter(estimator != "Difference-in-means") |>
  filter(between(dgp_num, 1, 4), T0 == 4) |>
  summarize_ests_wide(which_rel_year = 0) |>
  print("markdown") |>
  extract_tt_latex_body() |>
  cat(file = here("out/tables/simulation-1/T0_4.tex"))

cat("\n\n\n")
ests |>
  filter(estimator != "Difference-in-means") |>
  filter(between(dgp_num, 1, 4), T0 == 12) |>
  summarize_ests_wide(which_rel_year = 0) |>
  print("markdown") |>
  extract_tt_latex_body() |>
  cat(file = here("out/tables/simulation-1/T0_12.tex"))


cat("\n\n\n")
ests |>
  filter(estimator != "Difference-in-means") |>
  filter(between(dgp_num, 1, 4), T0 == 4) |>
  summarize_ests_wide_ci(which_rel_year = 0) |>
  # print("markdown") |>
  extract_tt_latex_body() |>
  cat(file = here("out/tables/simulation-1/T0_4_ci.tex"))

cat("\n\n\n")
ests |>
  filter(estimator != "Difference-in-means") |>
  filter(between(dgp_num, 1, 4), T0 == 12) |>
  summarize_ests_wide_ci(which_rel_year = 0) |>
  # print("markdown") |>
  extract_tt_latex_body() |>
  cat(file = here("out/tables/simulation-1/T0_12_ci.tex"))
