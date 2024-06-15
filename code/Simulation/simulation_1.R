# %%
library(tidyverse)
library(collapse)
library(fixest)
library(arrow)
library(here)
library(fs)
library(yaml)
library(glue)
library(stringmagic)
library(furrr)

library(tinytable)
options(tinytable_print_output = "markdown")

dir_create(here("data/Simulations"))
dir_create(here("out/tables/simulation-1/"))
source(here("code/Simulation/dgp.R"))
source(here("code/Simulation/estimators.R"))
source(here("code/Simulation/helpers.R"))

RUN_SIMULATION <- FALSE

# %%
estimators <- tibble(
  est_function = list(
    function(df) {
      est_twfe(df)
    },
    function(df) {
      est_twfe_covs(df)
    },
    function(df) {
      est_synth(df)
    },
    function(df) {
      est_augsynth(df)
    },
    function(df) {
      est_gsynth(df, p = 2L)
    },
    function(df) {
      est_gsynth(df, p = c(0L, 3L))
    },
    function(df) {
      est_qld_F_known(df, within_transform = FALSE)
    },
    function(df) {
      est_qld(df, within_transform = FALSE, p = 2)
    },
    function(df) {
      est_qld(df, within_transform = FALSE, p = -1)
    }
  ),
  estimator = c(
    "TWFE",
    "TWFE with $\\bm{W}_i \\beta_t$",
    "Synthetic Control",
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
    "synth",
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
  5, 200L, 12L, 15L, TRUE, TRUE, FALSE, 1,
  6, 200L, 12L, 15L, FALSE, TRUE, FALSE, 1,
  7, 200L, 12L, 15L, FALSE, FALSE, FALSE, 1,
  8, 200L, 12L, 15L, FALSE, FALSE, TRUE, 1,
)

B <- 2500

# %% 
if (RUN_SIMULATION == TRUE) {
  tictoc::tic()
  ests <- run_simulation(B, dgps, estimators, cores = 8, seed = 20240518)
  tictoc::toc()
  
  write_csv(ests, here("data/Simulations/simulation_1_ests.csv"))
}

#' ## Report on simulation
# %%
ests <- read_csv(here("data/Simulations/simulation_1_ests.csv"), show_col_types = FALSE)

ests <- ests |>
  left_join(
    estimators |> select(estimator, estimator_short)
  )

for (curr_dgp_num in dgps$dgp_num) {
  out = here(glue("out/tables/simulation-1/dgp{curr_dgp_num}.tex"))
  
  ests |> 
    summarize_ests(dgps, dgp_num = curr_dgp_num) |>
    # print() |>
    extract_tt_latex_body() |>
    cat(file = out)
}

cat("\n\n\n")
out = here(glue("out/tables/simulation-1/T0_4.tex"))
ests |> 
  filter(between(dgp_num, 1, 4)) |>
  summarize_ests_wide(which_rel_year = 0) |>
  print() |>
  extract_tt_latex_body() |>
  cat(file = out)

cat("\n\n\n")
out = here(glue("out/tables/simulation-1/T0_10.tex"))
ests |> 
  filter(between(dgp_num, 5, 8)) |>
  summarize_ests_wide(which_rel_year = 0) |>
  print() |>
  extract_tt_latex_body() |>
  cat(file = out)

