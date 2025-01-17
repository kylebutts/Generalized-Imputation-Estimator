#' ---
#' title: "Create figures for Walmart application"
#' author: "Kyle Butts and Nicholas Brown"
#' ---

# %% setup
#| warning: false
library(tidyverse)
library(data.table)
library(ggplot2)
library(here)
library(fixest)
library(kfbmisc)


# Load Estimates ---------------------------------------------------------------
# %%
short_t_qld_retail <- fread(here("estimates/short_t_qld_p_3_est_retail_B_5000.csv"))
short_t_qld_retail <- as.matrix(short_t_qld_retail)
short_t_qld_retail <- data.table(
  rel_year = -15L:6L,
  estimate = short_t_qld_retail[1, ],
  std_error = drop(apply(short_t_qld_retail, 2, sd)),
  lower = drop(apply(short_t_qld_retail, 2, quantile, probs = 0.025)),
  upper = drop(apply(short_t_qld_retail, 2, quantile, probs = 0.975))
)
short_t_qld_retail$pre <- short_t_qld_retail$rel_year < 0
short_t_qld_retail$lower_normal <-
  short_t_qld_retail$estimate - 1.96 * short_t_qld_retail$std_error
short_t_qld_retail$upper_normal <-
  short_t_qld_retail$estimate + 1.96 * short_t_qld_retail$std_error

# %%
# Estimate pre-trend
pre_short_t_qld_retail <- coef(feols(
  estimate ~ rel_year, short_t_qld_retail[rel_year < 0 & rel_year >= -15]
))

(plot_short_t_qld_retail <- ggplot() +
  geom_abline(
    slope = 0, intercept = 0,
    linetype = "dashed", linewidth = 0.8
  ) +
  geom_abline(
    slope = pre_short_t_qld_retail[2], intercept = pre_short_t_qld_retail[1],
    color = "#9A2415",
    linewidth = 2
  ) +
  geom_point(
    aes(x = rel_year, y = estimate, color = pre),
    data = short_t_qld_retail,
    size = 3
  ) +
  # Empirical 95%
  # geom_errorbar(
  #   aes(x = rel_year, ymin = lower, ymax = upper, color = pre),
  #   data = short_t_qld_retail,
  #   width = 0.6, linewidth = 2
  # ) +
  # 1.96 * standard error
  geom_errorbar(
    aes(x = rel_year, ymin = lower_normal, ymax = upper_normal, color = pre),
    data = short_t_qld_retail,
    width = 0.6, linewidth = 2
  ) +
  labs(
    x = "Event Time", y = NULL
  ) +
  scale_y_continuous(limits = c(-0.2, 0.3), expand = c(0,0)) +
  scale_x_continuous(limits = c(-22, 13)) +
  scale_color_manual(
    values = c("#107895", "#9A2415"),
    guide = "none"
  ) +
  kfbmisc::theme_kyle(base_size = 16)
)

# %% 
kfbmisc::tikzsave(
  here("out/figures/short_t_qld_retail.pdf"),
  plot_short_t_qld_retail, width = 10, height = 5
)

# %% 
short_t_qld_wholesale <- fread(here("estimates/short_t_qld_p_2_est_wholesale_B_5000.csv"))
short_t_qld_wholesale <- as.matrix(short_t_qld_wholesale)
short_t_qld_wholesale <- short_t_qld_wholesale[apply(short_t_qld_wholesale, 1, \(x) !any(is.nan(x))), ]
short_t_qld_wholesale <- data.table(
  rel_year = -15L:6L,
  estimate = short_t_qld_wholesale[1, ],
  std_error = drop(apply(short_t_qld_wholesale, 2, sd)),
  lower = drop(apply(short_t_qld_wholesale, 2, quantile, probs = 0.025)),
  upper = drop(apply(short_t_qld_wholesale, 2, quantile, probs = 0.975))
)
short_t_qld_wholesale$pre <- short_t_qld_wholesale$rel_year < 0
short_t_qld_wholesale$lower_normal <-
  short_t_qld_wholesale$estimate - 1.96 * short_t_qld_wholesale$std_error
short_t_qld_wholesale$upper_normal <-
  short_t_qld_wholesale$estimate + 1.96 * short_t_qld_wholesale$std_error

# %%
# Estimate pre-trend
pre_short_t_qld_wholesale <- coef(feols(
  estimate ~ rel_year, short_t_qld_wholesale[rel_year < 0 & rel_year >= -15]
))


(plot_short_t_qld_wholesale <- ggplot() +
  geom_abline(
    slope = 0, intercept = 0,
    linetype = "dashed", linewidth = 0.8
  ) +
  geom_abline(
    slope = pre_short_t_qld_wholesale[2], intercept = pre_short_t_qld_wholesale[1],
    color = "#9A2415",
    linewidth = 2
  ) +
  geom_point(
    aes(x = rel_year, y = estimate, color = pre),
    data = short_t_qld_wholesale,
    size = 3
  ) +
  # Empirical 95%
  # geom_errorbar(
  #   aes(x = rel_year, ymin = lower, ymax = upper, color = pre),
  #   data = short_t_qld_wholesale,
  #   width = 0.6, linewidth = 2
  # ) +
  # 1.96 * standard error
  geom_errorbar(
    aes(x = rel_year, ymin = lower_normal, ymax = upper_normal, color = pre),
    data = short_t_qld_wholesale,
    width = 0.6, linewidth = 2
  ) +
  labs(
    x = "Event Time", y = NULL
  ) +
  scale_y_continuous(limits = c(-0.5, 0.2), expand = c(0,0)) +
  scale_x_continuous(limits = c(-22, 13)) +
  scale_color_manual(
    values = c("#107895", "#9A2415"),
    guide = "none"
  ) +
  kfbmisc::theme_kyle(base_size = 16)
)

# %% 
kfbmisc::tikzsave(
  here("out/figures/short_t_qld_wholesale.pdf"),
  plot_short_t_qld_wholesale, width = 10, height = 5
)
