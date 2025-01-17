# %%
library(tidyverse)
library(fixest)
library(here)
library(fs)
library(furrr)
library(glue)

library(tinytable)
options(tinytable_print_output = "markdown")

library(tikzDevice)
tikzDevice::setTikzDefaults()
default_packages <- getOption("tikzLatexPackages")
options("tikzLatexPackages" = c(default_packages, "\\usepackage{bm}\n"))

dir_create(here("data/Simulations"))
dir_create(here("out/tables/simulation-2/"))
source(here("code/Simulation/dgp.R"))
source(here("code/Simulation/estimators.R"))
source(here("code/Simulation/helpers.R"))

# Flags
RUN_SIMULATION <- TRUE

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
      est_qld(df, within_transform = FALSE, p = -1)
    }
  ),
  estimator = c(
    "TWFE",
    "TWFE with $\\bm{w}_i \\beta_t$",
    "QLD ($p$ estimated)"
  ),
  estimator_short = c(
    "twfe",
    "twfe_covs",
    "qld"
  )
)

dgps <- tribble(
  ~dgp_num, ~N, ~T0, ~T, ~twfe, ~parallel_trends, ~ar_error_term, ~signal_to_noise,
  1, 200L, 4L, 7L, FALSE, FALSE, FALSE, 0.01,
  2, 200L, 4L, 7L, FALSE, FALSE, FALSE, 0.025,
  3, 200L, 4L, 7L, FALSE, FALSE, FALSE, 0.05,
  4, 200L, 4L, 7L, FALSE, FALSE, FALSE, 0.1,
  5, 200L, 4L, 7L, FALSE, FALSE, FALSE, 0.2,
  6, 200L, 4L, 7L, FALSE, FALSE, FALSE, 0.3,
  7, 200L, 4L, 7L, FALSE, FALSE, FALSE, 0.4,
  8, 200L, 4L, 7L, FALSE, FALSE, FALSE, 0.5,
  9, 200L, 4L, 7L, FALSE, FALSE, FALSE, 0.6,
  10, 200L, 4L, 7L, FALSE, FALSE, FALSE, 0.7,
  11, 200L, 4L, 7L, FALSE, FALSE, FALSE, 0.8,
  12, 200L, 4L, 7L, FALSE, FALSE, FALSE, 0.9,
  13, 200L, 4L, 7L, FALSE, FALSE, FALSE, 1.0,
)

# stnr = v_\gamma / (v_\gamma + v_\xi)
# v_\xi = v_\gamma / stnr - v_\gamma
var_gamma <- 1
dgps$instrument_noise <- (var_gamma / dgps$signal_to_noise) - var_gamma

# B <- 25
B <- 5000

# %%
if (RUN_SIMULATION == TRUE) {
  tictoc::tic()
  ests <- run_simulation(B, dgps, estimators, cores = 8, seed = 20240518)
  tictoc::toc()

  write_csv(ests, here("data/Simulations/simulation_2_ests.csv"))
}

#' ## Report on simulation
# %%
ests <- read_csv(here("data/Simulations/simulation_2_ests.csv"), show_col_types = FALSE)

ests <- ests |>
  left_join(
    estimators |> select(estimator, estimator_short)
  )

summary <- ests |>
  filter(rel_year == 0) |>
  filter(estimator_short != "twfe") |>
  mutate(bias = est - true_te) |>
  summarize(
    .by = c(dgp_num, estimator),
    mean_bias = mean(bias),
    std_error_bias = sd(bias),
    bias_empirical_ci_upper = quantile(bias, 0.975),
    bias_empirical_ci_lower = quantile(bias, 0.025)
  ) |>
  left_join(
    dgps |> select(dgp_num, signal_to_noise, instrument_noise),
    by = "dgp_num"
  )

avg_bias_twfe <- ests |>
  filter(rel_year == 0) |>
  filter(estimator_short == "twfe") |>
  with(mean(est - true_te))


(plot_signal_to_noise <- ggplot() +
  # Average bias of TWFE
  # geom_hline(yintercept = avg_bias_twfe, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(
      x = signal_to_noise, y = mean_bias,
      color = estimator, shape = estimator
    ),
    data = summary,
    size = 3
  ) +
  geom_line(
    aes(
      x = signal_to_noise, y = mean_bias,
      color = estimator, group = estimator
    ),
    data = summary,
    linewidth = 1.5
  ) +
  geom_ribbon(
    aes(
      x = signal_to_noise, y = mean_bias,
      color = estimator,
      fill = estimator, group = estimator,
      # ymin = mean_bias - 1.96 * std_error_bias,
      # ymax = mean_bias + 1.96 * std_error_bias
      ymin = bias_empirical_ci_lower,
      ymax = bias_empirical_ci_upper
    ),
    data = summary,
    alpha = 0.2, show.legend = FALSE
  ) +
  # Average bias of TWFE
  # annotate(
  #   "label",
  #   x = 1, hjust = 1, y = avg_bias_twfe + 0.25,
  #   label = string_magic("bias of TWFE"),
  #   label.r = unit(0, "pt"), label.size = 0,
  #   label.padding = unit(0, "pt"),
  #   size = 5, fill = "white"
  # ) +
  labs(
    x = "Signal to Noise Ratio",
    y = "Bias",
    color = NULL, shape = NULL,
    fill = NULL
  ) +
  scale_color_manual(
    values = c("black", "#107895"),
    guide = guide_legend(
      byrow = TRUE,
      override.aes = list(size = 0)
    )
  ) +
  scale_fill_manual(
    values = c("black", "#107895"),
  ) +
  kfbmisc::theme_kyle(base_size = 16) +
  theme(
    legend.text = element_text(size = rel(1 / 1.2)),
    legend.title = element_text(size = rel(1.2)),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.2),
    legend.key.spacing.x = unit(0, "pt"),
    legend.key.spacing.y = unit(4, "pt"),
    legend.background = element_rect(colour = "black", linewidth = 0.5),
    legend.margin = margin(t = 4, r = 0, b = 4, l = 4, unit = "pt")
  ))


# %%
(plot_signal_to_noise_errorbars <- ggplot() +
  # Average bias of TWFE
  # geom_hline(yintercept = avg_bias_twfe, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(
    aes(
      x = signal_to_noise + 0.002 * grepl("QLD", estimator) - 0.002 * grepl("TWFE", estimator),
      y = mean_bias,
      color = estimator, shape = estimator
    ),
    data = summary,
    size = 3
  ) +
  geom_errorbar(
    aes(
      x = signal_to_noise + 0.002 * grepl("QLD", estimator) - 0.002 * grepl("TWFE", estimator),
      y = mean_bias,
      color = estimator,
      group = estimator,
      ymin = bias_empirical_ci_lower,
      ymax = bias_empirical_ci_upper
    ),
    data = summary,
    linewidth = 1.25
  ) +
  # Average bias of TWFE
  # annotate(
  #   "label",
  #   x = 1, hjust = 1, y = avg_bias_twfe + 0.25,
  #   label = string_magic("bias of TWFE"),
  #   label.r = unit(0, "pt"), label.size = 0,
  #   label.padding = unit(0, "pt"),
  #   size = 5, fill = "white"
  # ) +
  labs(
    x = "Signal to Noise Ratio",
    y = "Bias",
    color = NULL, shape = NULL,
    fill = NULL
  ) +
  scale_color_manual(
    values = c("black", "#107895"),
    guide = guide_legend(
      byrow = TRUE,
      override.aes = list(size = 0)
    )
  ) +
  scale_fill_manual(
    values = c("black", "#107895"),
  ) +
  kfbmisc::theme_kyle(base_size = 16) +
  theme(
    legend.text = element_text(size = rel(1 / 1.2)),
    legend.title = element_text(size = rel(1.2)),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.2),
    legend.key.spacing.x = unit(0, "pt"),
    legend.key.spacing.y = unit(4, "pt"),
    legend.background = element_rect(colour = "black", linewidth = 0.5),
    legend.margin = margin(t = 4, r = 0, b = 4, l = 4, unit = "pt")
  ))

# %%
kfbmisc::tikzsave(
  here("out/figures/simulation-2/bias_signal_to_noise.pdf"),
  plot_signal_to_noise,
  width = 10, height = 5
)
kfbmisc::tikzsave(
  here("out/figures/simulation-2/bias_signal_to_noise_errorbars.pdf"),
  plot_signal_to_noise_errorbars,
  width = 10, height = 5
)
