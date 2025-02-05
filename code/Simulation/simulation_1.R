# %%
library(tidyverse)
library(fixest)
library(here)
library(fs)
library(furrr)
library(glue)

library(tinytable)
options(tinytable_print_output = "markdown")

fs::dir_create(here("data/Simulations"))
fs::dir_create(here("out/tables/simulation-1/"))
source(here("code/Simulation/dgp.R"))
source(here("code/Simulation/estimators.R"))
source(here("code/Simulation/run_simulation.R"))
source(here("code/Simulation/summarize.R"))

RUN_SIMULATION <- TRUE

# %%
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
    "TWFE with $\\bm{w}_i \\beta_t$",
    "Augmented Synthetic Control",
    "Generalized Synth ($p$ known)",
    "Generalized Synth ($p$ estimated)",
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
    "qld_F_known",
    "qld_p_known",
    "qld"
  )
)

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

# %%
if (RUN_SIMULATION == TRUE) {
  tictoc::tic()
  ests <- run_simulation(B, dgps, estimators, seed = 20250121)
  tictoc::toc()

  write_csv(ests, here("data/Simulations/simulation_1_ests.csv"))
}


# %%
## Report on simulation
ests <- read_csv(here("data/Simulations/simulation_1_ests.csv"), show_col_types = FALSE)

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

# %%
for (curr_dgp_num in dgps$dgp_num) {
  for (curr_T0 in c(4L, 12L)) {
    out <- here(glue("out/tables/simulation-1/dgp{curr_dgp_num}_T0_{curr_T0}.tex"))
    ests |> 
      filter(estimator != "Difference-in-means") |>
      summarize_ests(dgps, dgp_num = curr_dgp_num, T0 = curr_T0) |>
      # print() |>
      extract_tt_latex_body() |>
      cat(file = out)
    # cat()
  }
}

cat("\n\n\n")
out <- here(glue("out/tables/simulation-1/T0_4.tex"))
ests |>
  filter(estimator != "Difference-in-means") |>
  filter(between(dgp_num, 1, 4), T0 == 4) |>
  summarize_ests_wide(which_rel_year = 0) |>
  print() |>
  extract_tt_latex_body() |>
  cat(file = out)

cat("\n\n\n")
out <- here(glue("out/tables/simulation-1/T0_12.tex"))
ests |>
  filter(estimator != "Difference-in-means") |>
  filter(between(dgp_num, 1, 4), T0 == 12) |>
  summarize_ests_wide(which_rel_year = 0) |>
  print() |>
  extract_tt_latex_body() |>
  cat(file = out)
