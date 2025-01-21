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
source(here("code/Simulation/helpers.R"))

RUN_SIMULATION <- FALSE
DO_INFERENCE <- TRUE

# %%
estimators <- tibble(
  est_function = list(
    function(df) {
      est_twfe(df)
    },
    function(df) {
      est_twfe_covs(df)
    },
    function(df, do_inference = FALSE) {
      est_augsynth(df, do_inference)
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
  ~dgp_num, ~N, ~T0, ~T, ~twfe, ~parallel_trends, ~ar_error_term, ~instrument_noise,
  1, 200L, 4L, 7L, TRUE, TRUE, FALSE, 1,
  2, 200L, 4L, 7L, FALSE, TRUE, FALSE, 1,
  3, 200L, 4L, 7L, FALSE, FALSE, FALSE, 1,
  4, 200L, 4L, 7L, FALSE, FALSE, TRUE, 1,
  1, 200L, 12L, 15L, TRUE, TRUE, FALSE, 1,
  2, 200L, 12L, 15L, FALSE, TRUE, FALSE, 1,
  3, 200L, 12L, 15L, FALSE, FALSE, FALSE, 1,
  4, 200L, 12L, 15L, FALSE, FALSE, TRUE, 1,
)

B <- 1000

# %%
if (RUN_SIMULATION == TRUE) {
  tictoc::tic()
  ests <- run_simulation(B, dgps, estimators, do_inference = DO_INFERENCE, seed = 20250121)
  tictoc::toc()

  if (DO_INFERENCE == FALSE) {
    write_csv(ests, here("data/Simulations/simulation_1_ests_no_inference.csv"))
  } else {
    write_csv(ests, here("data/Simulations/simulation_1_ests.csv"))
  }
}

#' ## Report on simulation
# %%
ests <- read_csv(here("data/Simulations/simulation_1_ests.csv"), show_col_types = FALSE)
# ests <- read_csv(here("data/Simulations/simulation_1_ests_no_inference.csv"), show_col_types = FALSE)

for (curr_dgp_num in dgps$dgp_num) {
  for (curr_T_0 in c(4L, 12L)) {
    out <- here(glue("out/tables/simulation-1/dgp{curr_dgp_num}_T0_{curr_T_0}.tex"))

    summarize_ests(ests, dgps, dgp_num = curr_dgp_num, T_0 = curr_T_0) |>
      # print() |>
      extract_tt_latex_body() |>
      cat(file = out)
    # cat()
  }
}

cat("\n\n\n")
out <- here(glue("out/tables/simulation-1/T0_4.tex"))
ests |>
  filter(between(dgp_num, 1, 4), T_0 == 4) |>
  summarize_ests_wide(which_rel_year = 0) |>
  print() |>
  extract_tt_latex_body() |>
  cat(file = out)

cat("\n\n\n")
out <- here(glue("out/tables/simulation-1/T0_12.tex"))
ests |>
  filter(between(dgp_num, 1, 4), T_0 == 12) |>
  summarize_ests_wide(which_rel_year = 0) |>
  print() |>
  extract_tt_latex_body() |>
  cat(file = out)
