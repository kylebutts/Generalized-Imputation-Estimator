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
options(readr.show_col_types = FALSE)


# Load Estimates ---------------------------------------------------------------
# %%
#| label: "Load QLD"
qld_retail <- read_csv(here("estimates/qld_p_2_est_retail_B_5000.csv"))
qld_retail <- tibble(
  rel_year = colnames(qld_retail) |> str_replace("tau", "") |> as.numeric(),
  pre = rel_year < 0,
  estimate = as.numeric(qld_retail[1, ]),
  std_error = drop(apply(qld_retail, 2, sd)),
  lower = drop(apply(qld_retail, 2, quantile, probs = 0.025)),
  upper = drop(apply(qld_retail, 2, quantile, probs = 0.975)),
  lower_normal = estimate - 1.96 * std_error,
  upper_normal = estimate + 1.96 * std_error
)

qld_wholesale <- read_csv(here("estimates/qld_p_1_est_wholesale_B_5000.csv"))
qld_wholesale <- tibble(
  rel_year = colnames(qld_wholesale) |> str_replace("tau", "") |> as.numeric(),
  pre = rel_year < 0,
  estimate = as.numeric(qld_wholesale[1, ]),
  std_error = drop(apply(qld_wholesale, 2, sd)),
  lower = drop(apply(qld_wholesale, 2, quantile, probs = 0.025)),
  upper = drop(apply(qld_wholesale, 2, quantile, probs = 0.975)),
  lower_normal = estimate - 1.96 * std_error,
  upper_normal = estimate + 1.96 * std_error
)

# %%
#| label: "Load did2s"
did2s_retail <- read_csv(here("estimates/did2s_est_retail_B_1000.csv"))
did2s_retail$pre <- did2s_retail$rel_year < 0

did2s_wholesale <- read_csv(here("estimates/did2s_est_wholesale_B_1000.csv"))
did2s_wholesale$pre <- did2s_wholesale$rel_year < 0

# %%
#| label: "Load did2s (with QLD instr. as covs)"
did2s_covs_retail <- read_csv(here("estimates/did2s_est_covs_retail.csv"))
did2s_covs_retail$pre <- did2s_covs_retail$rel_year < 0

did2s_covs_wholesale <- read_csv(here("estimates/did2s_est_covs_wholesale.csv"))
did2s_covs_wholesale$pre <- did2s_covs_wholesale$rel_year < 0

# %%
#| label: "Load PCA"
pca_retail <- read_csv(here("estimates/pca_est_retail.csv"))
pca_retail <- tibble(
  rel_year = -22L:13L, estimate = as.numeric(pca_retail[1, ])
)
pca_retail$pre <- pca_retail$rel_year < 0

pca_wholesale <- read_csv(here("estimates/pca_est_wholesale.csv"))
pca_wholesale <- tibble(
  rel_year = -22L:13L, estimate = as.numeric(pca_wholesale[1, ])
)
pca_wholesale$pre <- pca_wholesale$rel_year < 0

# %%
#| label: "Load CCE"

cce_retail <- read_csv(here("estimates/cce_est_retail.csv"))
cce_retail <- tibble(
  rel_year = -22L:13L, estimate = as.numeric(cce_retail[1, ])
)
cce_retail$pre <- cce_retail$rel_year < 0

cce_wholesale <- read_csv(here("estimates/cce_est_wholesale.csv"))
cce_wholesale <- tibble(
  rel_year = -22L:13L, estimate = as.numeric(cce_wholesale[1, ])
)
cce_wholesale$pre <- cce_wholesale$rel_year < 0

# %%
#| label: "Load naive SEs"
# Load Naive SEs
qld_retail_naive_se <- 
  read_csv(here("estimates/qld_p_2_est_retail_naive_se.csv")) |>
  mutate(
    pre = rel_year < 0,
    lower = estimate - 1.96 * std_error,
    upper = estimate + 1.96 * std_error
  )

qld_wholesale_naive_se <-   
  read_csv(here("estimates/qld_p_1_est_wholesale_naive_se.csv")) |>
  mutate(
    pre = rel_year < 0,
    lower = estimate - 1.96 * std_error,
    upper = estimate + 1.96 * std_error
  )

# %%
#| label: "Load Synthetic Control Estimates"
est_y0_qld = read_csv(here("estimates/est_y0_qld.csv"))
synth_retail <- 
  est_y0_qld |>
  filter(g != Inf) |>
  summarize(
    ytilde_mean = mean(log_retail_emp_tilde, na.rm = TRUE),
    ytilde_0hat_mean = mean(log_retail_emp_tilde_0hat_qld_p_2, na.rm = TRUE),
    .by = c(rel_year)
  )
synth_retail <- rbind(
  synth_retail |>
    select(rel_year, mean = ytilde_mean) |>
    mutate(group = r'(Average of $\tilde{y}_{it}$)'),
  synth_retail |>
    select(rel_year, mean = ytilde_0hat_mean) |>
    mutate(group = r'(Average of $\hat{\tilde{y}}_{it}(0)$)')
)

synth_wholesale <- 
  est_y0_qld |>
  filter(g != Inf) |>
  summarize(
    ytilde_mean = mean(log_wholesale_emp_tilde, na.rm = TRUE),
    ytilde_0hat_mean = mean(log_wholesale_emp_tilde_0hat_qld_p_1, na.rm = TRUE),
    .by = c(rel_year)
  )
synth_wholesale <- rbind(
  synth_wholesale |>
    select(rel_year, mean = ytilde_mean) |>
    mutate(group = r'(Average of $\tilde{y}_{it}$)'),
  synth_wholesale |>
    select(rel_year, mean = ytilde_0hat_mean) |>
    mutate(group = r'(Average of $\hat{\tilde{y}}_{it}(0)$)')
)

# Figures ----------------------------------------------------------------------
# %%
#| label: "did2s retail"
# Estimate pre-trend
pre_did2s_retail <- coef(feols(
  estimate ~ rel_year, did2s_retail |> filter(rel_year < 0 & rel_year >= -15)
))

(plot_did2s_retail <- ggplot() +
  geom_abline(
    slope = 0, intercept = 0,
    linetype = "dashed", linewidth = 0.8
  ) +
  geom_abline(
    slope = pre_did2s_retail[2], intercept = pre_did2s_retail[1],
    color = "#9A2415",
    linewidth = 2
  ) +
  geom_point(
    aes(x = rel_year, y = estimate, color = pre),
    data = did2s_retail,
    size = 3
  ) +
  # Empirical 95%
  geom_errorbar(
    aes(x = rel_year, ymin = lower, ymax = upper, color = pre),
    data = did2s_retail,
    width = 0.6, linewidth = 2
  ) +
  labs(
    x = "Event Time", y = NULL
  ) +
  scale_y_continuous(limits = c(-0.2, 0.3)) +
  scale_color_manual(
    values = c("#107895", "#9A2415"),
    guide = "none"
  ) +
  kfbmisc::theme_kyle(base_size = 16) +
  theme(legend.text = element_text(size = rel(1))))

# %%
#| label: "did2s wholesale"
pre_did2s_wholesale <- feols(
  estimate ~ rel_year, did2s_wholesale |> filter(rel_year < 0 & rel_year >= -15)
) |>
  coef()

(plot_did2s_wholesale <- ggplot() +
  geom_abline(
    slope = 0, intercept = 0,
    linetype = "dashed", linewidth = 0.8
  ) +
  geom_abline(
    slope = pre_did2s_wholesale[2], intercept = pre_did2s_wholesale[1],
    color = "#9A2415",
    linewidth = 2
  ) +
  geom_point(
    aes(x = rel_year, y = estimate, color = pre),
    data = did2s_wholesale,
    size = 3
  ) +
  # Empirical 95%
  geom_errorbar(
    aes(x = rel_year, ymin = lower, ymax = upper, color = pre),
    data = did2s_wholesale,
    width = 0.6, linewidth = 2
  ) +
  labs(
    x = "Event Time", y = NULL
  ) +
  scale_y_continuous(limits = c(-0.5, 0.2)) +
  scale_color_manual(
    values = c("#107895", "#9A2415"),
    guide = "none"
  ) +
  kfbmisc::theme_kyle(base_size = 16) +
  theme(legend.text = element_text(size = rel(1))))

# %%
#| label: "qld retail"
# Estimate pre-trend
pre_qld_retail <- coef(feols(
  estimate ~ rel_year, qld_retail |> filter(rel_year < 0 & rel_year >= -15)
))

(plot_qld_retail <- ggplot() +
  geom_abline(
    slope = 0, intercept = 0,
    linetype = "dashed", linewidth = 0.8
  ) +
  geom_abline(
    slope = pre_qld_retail[2], intercept = pre_qld_retail[1],
    color = "#9A2415",
    linewidth = 2
  ) +
  geom_point(
    aes(x = rel_year, y = estimate, color = pre),
    data = qld_retail,
    size = 3
  ) +
  # Empirical 95%
  geom_errorbar(
    aes(x = rel_year, ymin = lower, ymax = upper, color = pre),
    data = qld_retail,
    width = 0.6, linewidth = 2
  ) +
  labs(
    x = "Event Time", y = NULL
  ) +
  scale_y_continuous(limits = c(-0.2, 0.3)) +
  scale_color_manual(
    values = c("#107895", "#9A2415"),
    guide = "none"
  ) +
  kfbmisc::theme_kyle(base_size = 16) +
  theme(
    legend.text = element_text(size = rel(1))
  ))

# %%
#| label: "qld wholesale"
pre_qld_wholesale <- coef(feols(
  estimate ~ rel_year, qld_wholesale |> filter(rel_year < 0 & rel_year >= -15)
))

(plot_qld_wholesale <- ggplot() +
  geom_abline(
    slope = 0, intercept = 0,
    linetype = "dashed", linewidth = 0.8
  ) +
  geom_abline(
    slope = pre_qld_wholesale[2], intercept = pre_qld_wholesale[1],
    color = "#9A2415",
    linewidth = 2
  ) +
  geom_point(
    aes(x = rel_year, y = estimate, color = pre),
    data = qld_wholesale,
    size = 3
  ) +
  # Empirical 95%
  geom_errorbar(
    aes(x = rel_year, ymin = lower, ymax = upper, color = pre),
    data = qld_wholesale,
    width = 0.6, linewidth = 2
  ) +
  labs(
    x = "Event Time", y = NULL
  ) +
  scale_y_continuous(limits = c(-0.5, 0.2)) +
  scale_color_manual(
    values = c("#107895", "#9A2415"),
    guide = "none"
  ) +
  kfbmisc::theme_kyle(base_size = 16) +
  theme(
    legend.text = element_text(size = rel(1))
  ))

# %%
#| label: "Many Estimators"
cce_retail$group <- "Common Correlated Effects"
pca_retail$group <- "Principal Components"
qld_retail$group <- "Quasi-Long Differencing"
retail_estimators <- data.table::rbindlist(
  list(cce_retail, pca_retail, qld_retail),
  fill = TRUE
)

cce_wholesale$group <- "Common Correlated Effects"
pca_wholesale$group <- "Principal Components"
qld_wholesale$group <- "Quasi-Long Differencing"
wholesale_estimators <- data.table::rbindlist(
  list(cce_wholesale, pca_wholesale, qld_wholesale),
  fill = TRUE
)

# %%
(plot_retail_many_estimators <- ggplot() +
  geom_abline(
    slope = 0, intercept = 0,
    linetype = "dashed", linewidth = 0.8
  ) +
  geom_point(
    aes(x = rel_year, y = estimate, color = group),
    data = retail_estimators,
    size = 3
  ) +
  geom_line(
    aes(x = rel_year, y = estimate, color = group, group = group),
    data = retail_estimators,
    linewidth = 2
  ) +
  scale_y_continuous(limits = c(-0.2, 0.18)) +
  scale_color_manual(
    values = c(
      "Quasi-Long Differencing" = "grey80",
      "Common Correlated Effects" = "grey20",
      "Principal Components" = "grey50"
    ),
    guide = guide_legend(
      override.aes = list(size = 0, linewidth = 2.2),
      byrow = TRUE
    )
  ) +
  labs(
    x = "Event Time", y = NULL, color = NULL, linetype = NULL
  ) +
  kfbmisc::theme_kyle(base_size = 16) +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.2)),
    legend.position = "inside",
    legend.position.inside = c(0.25, 0.8),
    legend.spacing.y = unit(10, "pt"),
    legend.background = element_rect(colour = "black", linewidth = 0.5),
    legend.margin = margin(t = 2, r = 8, b = 8, l = 8, unit = "pt")
  ))

# %%
(plot_wholesale_many_estimators <- ggplot() +
  geom_abline(
    slope = 0, intercept = 0,
    linetype = "dashed", linewidth = 0.8
  ) +
  geom_point(
    aes(x = rel_year, y = estimate, color = group),
    data = wholesale_estimators,
    size = 3
  ) +
  geom_line(
    aes(x = rel_year, y = estimate, color = group, group = group),
    data = wholesale_estimators,
    linewidth = 2
  ) +
  scale_y_continuous(limits = c(-0.2, 0.18)) +
  scale_color_manual(
    values = c(
      "Quasi-Long Differencing" = "grey80",
      "Common Correlated Effects" = "grey20",
      "Principal Components" = "grey50"
    ),
    guide = guide_legend(
      override.aes = list(size = 0, linewidth = 2.2),
      byrow = TRUE
    )
  ) +
  labs(
    x = "Event Time", y = NULL, color = NULL, linetype = NULL
  ) +
  kfbmisc::theme_kyle(base_size = 16) +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.2)),
    legend.position = "inside",
    legend.position.inside = c(0.25, 0.8),
    legend.spacing.y = unit(10, "pt"),
    legend.background = element_rect(colour = "black", linewidth = 0.5),
    legend.margin = margin(t = 2, r = 8, b = 8, l = 8, unit = "pt")
  ))

# %%
#| label: "QLD naive estimates"
(plot_qld_retail_naive_se <- ggplot() +
  geom_abline(
    slope = 0, intercept = 0,
    linetype = "dashed", linewidth = 0.8
  ) +
  geom_point(
    aes(x = rel_year, y = estimate, color = pre),
    data = qld_retail_naive_se,
    size = 3
  ) +
  geom_errorbar(
    aes(
      x = rel_year, ymin = lower, ymax = upper, color = pre
    ),
    data = qld_retail_naive_se,
    width = 0.6, linewidth = 2
  ) +
  labs(
    x = "Event Time", y = NULL
  ) +
  scale_y_continuous(limits = c(-0.2, 0.3)) +
  scale_color_manual(
    values = c("#107895", "#9A2415"),
    guide = "none"
  ) +
  kfbmisc::theme_kyle(base_size = 16) +
  theme(legend.text = element_text(size = rel(1))))

(plot_qld_wholesale_naive_se <- ggplot() +
  geom_abline(
    slope = 0, intercept = 0,
    linetype = "dashed", linewidth = 0.8
  ) +
  geom_point(
    aes(x = rel_year, y = estimate, color = pre),
    data = qld_wholesale,
    size = 3
  ) +
  geom_errorbar(
    aes(
      x = rel_year,
      ymin = lower,
      ymax = upper,
      color = pre
    ),
    data = qld_wholesale_naive_se,
    width = 0.6, linewidth = 2
  ) +
  labs(
    x = "Event Time", y = NULL
  ) +
  scale_y_continuous(limits = c(-0.4, 0.2)) +
  scale_color_manual(
    values = c("#107895", "#9A2415"),
    guide = "none"
  ) +
  kfbmisc::theme_kyle(base_size = 16) +
  theme(legend.text = element_text(size = rel(1))))

# %%
#| label: "Synthetic Control QLD Plots"
(plot_synth_retail <- ggplot() +
  geom_line(
    aes(x = rel_year, y = mean, color = group),
    data = synth_retail,
    linewidth = 2
  ) +
  labs(
    x = "Event Time", y = NULL,
    color = NULL
  ) +
  scale_y_continuous(limits = c(-0.1, 0.3)) +
  scale_color_manual(
    values = c("#107895", "#9A2415"),
    guide = guide_legend(byrow = TRUE)
  ) +
  kfbmisc::theme_kyle(base_size = 16) +
  theme(legend.text = element_text(size = rel(1))) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.2, 0.8),
    legend.background = element_rect(colour = "black", linewidth = 0.5),
    legend.spacing.y = unit(8, "pt"),
    legend.margin = margin(t = 2, r = 8, b = 8, l = 8, unit = "pt")
  ))

(plot_synth_wholesale <- ggplot() +
  geom_line(
    aes(x = rel_year, y = mean, color = group),
    data = synth_wholesale,
    linewidth = 2
  ) +
  labs(
    x = "Event Time", y = NULL,
    color = NULL
  ) +
  scale_y_continuous(limits = c(-0.1, 0.3)) +
  scale_color_manual(
    values = c("#107895", "#9A2415"),
    guide = guide_legend(byrow = TRUE)
  ) +
  kfbmisc::theme_kyle(base_size = 16) +
  theme(legend.text = element_text(size = rel(1))) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.2, 0.8),
    legend.background = element_rect(
      colour = "black", linewidth = 0.5
    ),
    legend.spacing.y = unit(8, "pt"),
    legend.margin = margin(t = 2, r = 8, b = 8, l = 8, unit = "pt")
  ))

# %%
#| label: "Covs. vs. QLD Plot"
retail_covs_ests <- bind_rows(
  did2s_retail |> mutate(group = "TWFE Imputation"),
  did2s_covs_retail |> mutate(group = "TWFE Imputation (with $w_i \\beta_t$)"),
  qld_retail |> mutate(group = "Quasi-Long Differencing")
)
wholesale_covs_ests <- bind_rows(
  did2s_wholesale |> mutate(group = "TWFE Imputation"),
  did2s_covs_wholesale |> mutate(group = "TWFE Imputation (with $w_i \\beta_t$)"),
  qld_wholesale |> mutate(group = "Quasi-Long Differencing")
)

# %%
(plot_retail_covs <- ggplot() +
  geom_abline(
    slope = 0, intercept = 0,
    linetype = "dashed", linewidth = 0.8
  ) +
  geom_point(
    aes(x = rel_year, y = estimate, color = group),
    data = retail_covs_ests,
    size = 3
  ) +
  geom_line(
    aes(x = rel_year, y = estimate, color = group, group = group),
    data = retail_covs_ests,
    linewidth = 2
  ) +
  scale_y_continuous(limits = c(-0.2, 0.24)) +
  scale_color_manual(
    values = c(
      "TWFE Imputation" = "grey80",
      "TWFE Imputation (with $w_i \\beta_t$)" = "grey50",
      "Quasi-Long Differencing" = "grey20"
    ),
    guide = guide_legend(
      override.aes = list(size = 0, linewidth = 2.2),
      byrow = TRUE
    )
  ) +
  labs(
    x = "Event Time", y = NULL, color = NULL, linetype = NULL
  ) +
  kfbmisc::theme_kyle(base_size = 16) +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.2)),
    legend.position = "inside",
    legend.position.inside = c(0.25, 0.8),
    legend.spacing.y = unit(10, "pt"),
    legend.background = element_rect(colour = "black", linewidth = 0.5),
    legend.margin = margin(t = 2, r = 8, b = 8, l = 8, unit = "pt")
  ))

# %%
(plot_wholesale_covs <- ggplot() +
  geom_abline(
    slope = 0, intercept = 0,
    linetype = "dashed", linewidth = 0.8
  ) +
  geom_point(
    aes(x = rel_year, y = estimate, color = group),
    data = wholesale_covs_ests,
    size = 3
  ) +
  geom_line(
    aes(x = rel_year, y = estimate, color = group, group = group),
    data = wholesale_covs_ests,
    linewidth = 2
  ) +
  scale_y_continuous(limits = c(-0.2, 0.24)) +
  scale_color_manual(
    values = c(
      "TWFE Imputation" = "grey80",
      "TWFE Imputation (with $w_i \\beta_t$)" = "grey50",
      "Quasi-Long Differencing" = "grey20"
    ),
    guide = guide_legend(
      override.aes = list(size = 0, linewidth = 2.2),
      byrow = TRUE
    )
  ) +
  labs(
    x = "Event Time", y = NULL, color = NULL, linetype = NULL
  ) +
  kfbmisc::theme_kyle(base_size = 16) +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.2)),
    legend.position = "inside",
    legend.position.inside = c(0.25, 0.8),
    legend.spacing.y = unit(10, "pt"),
    legend.background = element_rect(colour = "black", linewidth = 0.5),
    legend.margin = margin(t = 2, r = 8, b = 8, l = 8, unit = "pt")
  ))


# %%
cce_retail$group <- "Common Correlated Effects"
pca_retail$group <- "Principal Components"
qld_retail$group <- "Quasi-Long Differencing"
retail_estimators <- data.table::rbindlist(
  list(cce_retail, pca_retail, qld_retail),
  fill = TRUE
)
cce_wholesale$group <- "Common Correlated Effects"
pca_wholesale$group <- "Principal Components"
qld_wholesale$group <- "Quasi-Long Differencing"
wholesale_estimators <- data.table::rbindlist(
  list(cce_wholesale, pca_wholesale, qld_wholesale),
  fill = TRUE
)

(plot_retail_many_estimators <- ggplot() +
  geom_abline(
    slope = 0, intercept = 0,
    linetype = "dashed", linewidth = 0.8
  ) +
  geom_point(
    aes(x = rel_year, y = estimate, color = group),
    data = retail_estimators,
    size = 3
  ) +
  geom_line(
    aes(x = rel_year, y = estimate, color = group, group = group),
    data = retail_estimators,
    linewidth = 2
  ) +
  scale_y_continuous(limits = c(-0.2, 0.18)) +
  scale_color_manual(
    values = c(
      "Quasi-Long Differencing" = "grey80",
      "Common Correlated Effects" = "grey20",
      "Principal Components" = "grey50"
    ),
    guide = guide_legend(
      override.aes = list(size = 0, linewidth = 2.2),
      byrow = TRUE
    )
  ) +
  labs(
    x = "Event Time", y = NULL, color = NULL, linetype = NULL
  ) +
  kfbmisc::theme_kyle(base_size = 16) +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.2)),
    legend.position = "inside",
    legend.position.inside = c(0.25, 0.8),
    legend.spacing.y = unit(10, "pt"),
    legend.background = element_rect(colour = "black", linewidth = 0.5),
    legend.margin = margin(t = 2, r = 8, b = 8, l = 8, unit = "pt")
  ))

(plot_wholesale_many_estimators <- ggplot() +
  geom_abline(
    slope = 0, intercept = 0,
    linetype = "dashed", linewidth = 0.8
  ) +
  geom_point(
    aes(x = rel_year, y = estimate, color = group),
    data = wholesale_estimators,
    size = 3
  ) +
  geom_line(
    aes(x = rel_year, y = estimate, color = group, group = group),
    data = wholesale_estimators,
    linewidth = 2
  ) +
  scale_y_continuous(limits = c(-0.2, 0.18)) +
  scale_color_manual(
    values = c(
      "Quasi-Long Differencing" = "grey80",
      "Common Correlated Effects" = "grey20",
      "Principal Components" = "grey50"
    ),
    guide = guide_legend(
      override.aes = list(size = 0, linewidth = 2.2),
      byrow = TRUE
    )
  ) +
  labs(
    x = "Event Time", y = NULL, color = NULL, linetype = NULL
  ) +
  kfbmisc::theme_kyle(base_size = 16) +
  theme(
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.2)),
    legend.position = "inside",
    legend.position.inside = c(0.25, 0.8),
    legend.spacing.y = unit(10, "pt"),
    legend.background = element_rect(colour = "black", linewidth = 0.5),
    legend.margin = margin(t = 2, r = 8, b = 8, l = 8, unit = "pt")
  ))


# Export plots -----------------------------------------------------------------
# %%
## did2s
kfbmisc::tikzsave(
  here("out/figures/did2s_retail.pdf"),
  plot_did2s_retail,
  width = 10, height = 5
)
kfbmisc::tikzsave(
  here("out/figures/did2s_wholesale.pdf"),
  plot_did2s_wholesale,
  width = 10, height = 5
)

## QLD
kfbmisc::tikzsave(
  here("out/figures/qld_retail.pdf"),
  plot_qld_retail,
  width = 10, height = 5
)
kfbmisc::tikzsave(
  here("out/figures/qld_wholesale.pdf"),
  plot_qld_wholesale,
  width = 10, height = 5
)

## Many Factor Estimators
kfbmisc::tikzsave(
  here("out/figures/retail_many_estimators.pdf"),
  plot_retail_many_estimators,
  width = 10, height = 5
)
kfbmisc::tikzsave(
  here("out/figures/wholesale_many_estimators.pdf"),
  plot_wholesale_many_estimators,
  width = 10, height = 5
)

## Covariates
kfbmisc::tikzsave(
  here("out/figures/retail_covs.pdf"),
  plot_retail_covs,
  width = 10, height = 5
)
kfbmisc::tikzsave(
  here("out/figures/wholesale_covs.pdf"),
  plot_wholesale_covs,
  width = 10, height = 5
)

## QLD Naive SEs
kfbmisc::tikzsave(
  here("out/figures/qld_retail_naive_se.pdf"),
  plot_qld_retail_naive_se,
  width = 10, height = 5
)
kfbmisc::tikzsave(
  here("out/figures/qld_wholesale_naive_se.pdf"),
  plot_qld_wholesale_naive_se,
  width = 10, height = 5
)

## Synthetic Control
kfbmisc::tikzsave(
  here("out/figures/synth_retail.pdf"),
  plot_synth_retail,
  width = 10, height = 5
)
kfbmisc::tikzsave(
  here("out/figures/synth_wholesale.pdf"),
  plot_synth_wholesale,
  width = 10, height = 5
)
