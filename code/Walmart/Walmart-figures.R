library(data.table)
library(ggplot2)
library(here)
library(fixest)
library(tikzDevice)

source(here("figures/convert_tex_to_pdf.R"))

# Load Estimates ---------------------------------------------------------------

## QLD -------------------------------------------------------------------------

qld_retail = fread(here("data/qld_p_2_est_retail_B_1000.csv"))
qld_retail = as.matrix(qld_retail)
qld_wholesale = fread(here("data/qld_p_1_est_wholesale_B_1000.csv"))
qld_wholesale = as.matrix(qld_wholesale)

qld_retail = data.table(
  rel_year = -22L:13L,
  estimate = qld_retail[1, ],
  std_error = drop(apply(qld_retail, 2, sd)),
  lower = drop(apply(qld_retail, 2, quantile, probs = 0.025)),
  upper = drop(apply(qld_retail, 2, quantile, probs = 0.975))
)
qld_wholesale = data.table(
  rel_year = -22L:13L,
  estimate = qld_wholesale[1, ],
  std_error = drop(apply(qld_wholesale, 2, sd)),
  lower = drop(apply(qld_wholesale, 2, quantile, probs = 0.025)),
  upper = drop(apply(qld_wholesale, 2, quantile, probs = 0.975))
)

qld_retail$pre = qld_retail$rel_year < 0
qld_wholesale$pre = qld_wholesale$rel_year < 0

## did2s -----------------------------------------------------------------------
did2s_retail = fread(here("estimates/did2s_est_retail_B_1000.csv"))
did2s_wholesale = fread(here("estimates/did2s_est_wholesale_B_1000.csv"))

did2s_retail$rel_year = -22L:13L
did2s_wholesale$rel_year = -22L:13L

did2s_retail$pre = did2s_retail$rel_year < 0
did2s_wholesale$pre = did2s_wholesale$rel_year < 0

## PCA -------------------------------------------------------------------------

pca_retail = fread(here("estimates/pca_est_retail.csv"))
pca_wholesale = fread(here("estimates/pca_est_wholesale.csv"))

pca_retail = data.table(
  rel_year = -22L:13L,
  estimate = as.numeric(pca_retail[1, ])
)
pca_wholesale = data.table(
  rel_year = -22L:13L,
  estimate = as.numeric(pca_wholesale[1, ])
)

pca_retail$pre = pca_retail$rel_year < 0
pca_wholesale$pre = pca_wholesale$rel_year < 0

## CCE -------------------------------------------------------------------------

cce_retail = fread(here("estimates/cce_est_retail.csv"))
cce_wholesale = fread(here("estimates/cce_est_wholesale.csv"))

cce_retail = data.table(
  rel_year = -22L:13L,
  estimate = as.numeric(cce_retail[1, ])
)
cce_wholesale = data.table(
  rel_year = -22L:13L,
  estimate = as.numeric(cce_wholesale[1, ])
)

cce_retail$pre = cce_retail$rel_year < 0
cce_wholesale$pre = cce_wholesale$rel_year < 0

# Figures ----------------------------------------------------------------------

## did2s estimates -------------------------------------------------------------

# Estimate pre-trend
pre_did2s_retail = coef(feols(
  estimate ~ rel_year, did2s_retail[rel_year < 0 & rel_year >= -15]
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

pre_did2s_wholesale = feols(
  estimate ~ rel_year, did2s_wholesale[rel_year < 0 & rel_year >= -15]
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
  geom_errorbar(
    aes(x = rel_year, ymin = lower, ymax = upper, color = pre),
    data = did2s_wholesale,
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

## qld estimates ---------------------------------------------------------------

# Estimate pre-trend
pre_qld_retail = coef(feols(
  estimate ~ rel_year, qld_retail[rel_year < 0 & rel_year >= -15]
))

plot_qld_retail <- ggplot() + 
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
  theme(legend.text = element_text(size = rel(1)))

pre_qld_wholesale = coef(feols(
  estimate ~ rel_year, qld_wholesale[rel_year < 0 & rel_year >= -15]
))

plot_qld_wholesale <- ggplot() + 
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
  geom_errorbar(
    aes(x = rel_year, ymin = lower, ymax = upper, color = pre),
    data = qld_wholesale,
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
  theme(legend.text = element_text(size = rel(1)))

## Many Estimators -------------------------------------------------------------

cce_retail$group = "Common Correlated Effects"
pca_retail$group = "Principal Components"
qld_retail$group = "Quasi-Long Differencing"
retail_estimators = rbindlist(
  list(cce_retail, pca_retail, qld_retail), 
  fill = TRUE
)
cce_wholesale$group = "Common Correlated Effects"
pca_wholesale$group = "Principal Components"
qld_wholesale$group = "Quasi-Long Differencing"
wholesale_estimators = rbindlist(
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
    legend.position = c(0.25, 0.8),
    legend.spacing.y = unit(10, "pt"),
    legend.background = element_rect(colour = "black", linewidth = 0.5),
    legend.margin = margin(t = 2, r = 8, b = 8, l = 8, unit = "pt"),
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
    legend.position = c(0.25, 0.8),
    legend.spacing.y = unit(10, "pt"),
    legend.background = element_rect(colour = "black", linewidth = 0.5),
    legend.margin = margin(t = 2, r = 8, b = 8, l = 8, unit = "pt"),
  ))

## QLD naive estimates ---------------------------------------------------------

# Load Naive SEs
qld_retail_naive_se = fread(here("data/factor_est_retail_p_2_naive_se.csv"))
qld_wholesale_naive_se = fread(here("data/factor_est_wholesale_p_1_naive_se.csv"))

qld_retail_naive_se$pre = qld_retail_naive_se$rel_year < 0
qld_wholesale_naive_se$pre = qld_wholesale_naive_se$rel_year < 0

qld_retail_naive_se$lower = qld_retail_naive_se$estimate - 1.96 * qld_retail_naive_se$std_error
qld_retail_naive_se$upper = qld_retail_naive_se$estimate + 1.96 * qld_retail_naive_se$std_error
qld_wholesale_naive_se$lower = qld_wholesale_naive_se$estimate - 1.96 * qld_wholesale_naive_se$std_error
qld_wholesale_naive_se$upper = qld_wholesale_naive_se$estimate + 1.96 * qld_wholesale_naive_se$std_error

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
      x = rel_year, ymin = lower, ymax = upper, color = pre),
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
## Synthetic Control QLD Plots -------------------------------------------------

# Load Synthetic Control Estimates
synth_retail = fread(here("data/synthetic_retail_p_2.csv"))
synth_retail = rbind(
  synth_retail[, .(rel_year, mean = ytilde_mean, group = r'(Average of $\tilde{y}_{it}$)')],
  synth_retail[, .(rel_year, mean = ytilde_0_hat_mean, group = r'(Average of $\hat{\tilde{y}}_{it}(0)$)')]
)
synth_wholesale = fread(here("data/synthetic_wholesale_p_1.csv"))
synth_wholesale = rbind(
  synth_wholesale[, .(rel_year, mean = ytilde_mean, group = r'(Average of $\tilde{y}_{it}$)')],
  synth_wholesale[, .(rel_year, mean = ytilde_0_hat_mean, group = r'(Average of $\hat{\tilde{y}}_{it}(0)$)')]
)

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
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(colour = "black", linewidth = 0.5),
    legend.spacing.y = unit(8, "pt"),
    legend.margin = margin(t = 2, r = 8, b = 8, l = 8, unit = "pt"),
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
    legend.position = c(0.2, 0.8),
    legend.background = element_rect(
      colour = "black", linewidth = 0.5
    ),
    legend.spacing.y = unit(8, "pt"),
    legend.margin = margin(t = 2, r = 8, b = 8, l = 8, unit = "pt"),
  ))
# Export plots -----------------------------------------------------------------

## did2s 
tikzDevice::tikz(
  here("figures/plot_did2s_retail.tex"),
  width = 10, height = 5, standAlone = FALSE
)
plot(plot_did2s_retail)
dev.off()
compile_tikz(here("figures/plot_did2s_retail.tex"))

tikzDevice::tikz(
  here("figures/plot_did2s_wholesale.tex"),
  width = 10, height = 5, standAlone = FALSE
)
plot(plot_did2s_wholesale)
dev.off()
compile_tikz(here("figures/plot_did2s_wholesale.tex"))

## QLD 
tikzDevice::tikz(
  here("figures/plot_qld_retail.tex"),
  width = 10, height = 5, standAlone = FALSE
)
plot(plot_qld_retail)
dev.off()
compile_tikz(here("figures/plot_qld_retail.tex"))

tikzDevice::tikz(
  here("figures/plot_qld_wholesale.tex"),
  width = 10, height = 5, standAlone = FALSE
)
plot(plot_qld_wholesale)
dev.off()
compile_tikz(here("figures/plot_qld_wholesale.tex"))

## Many Factor Estimators
tikzDevice::tikz(
  here("figures/plot_retail_many_estimators.tex"),
  width = 10, height = 5, standAlone = FALSE
)
plot(plot_retail_many_estimators)
dev.off()
compile_tikz(here("figures/plot_retail_many_estimators.tex"))

tikzDevice::tikz(
  here("figures/plot_wholesale_many_estimators.tex"),
  width = 10, height = 5, standAlone = FALSE
)
plot(plot_wholesale_many_estimators)
dev.off()
compile_tikz(here("figures/plot_wholesale_many_estimators.tex"))

## QLD Naive SEs
tikzDevice::tikz(
  here("figures/plot_qld_retail_naive_se.tex"),
  width = 10, height = 5, standAlone = FALSE
)
plot(plot_qld_retail_naive_se)
dev.off()
compile_tikz(here("figures/plot_qld_retail_naive_se.tex"))

tikzDevice::tikz(
  here("figures/plot_qld_wholesale_naive_se.tex"),
  width = 10, height = 5, standAlone = FALSE
)
plot(plot_qld_wholesale_naive_se)
dev.off()
compile_tikz(here("figures/plot_qld_wholesale_naive_se.tex"))

## Synthetic Control
tikzDevice::tikz(
  here("figures/plot_synth_retail.tex"),
  width = 10, height = 5, standAlone = FALSE
)
plot(plot_synth_retail)
dev.off()
compile_tikz(here("figures/plot_synth_retail.tex"))

tikzDevice::tikz(
  here("figures/plot_synth_wholesale.tex"),
  width = 10, height = 5, standAlone = FALSE
)
plot(plot_synth_wholesale)
dev.off()
compile_tikz(here("figures/plot_synth_wholesale.tex"))


