library(data.table)
library(ggplot2)
library(kfbmisc)
library(here)

# Load Estimates ---------------------------------------------------------------

est_qld_retail = fread(here("data/factor_est_retail_YEARS_1977_1999_T0_1985_B_1000_p_2.csv"))[1, ]
est_cce_retail = fread(here("data/cce_est_retail.csv"))
est_pca_retail = fread(here("data/pca_est_retail.csv"))
est_qld_wholesale = fread(here("data/factor_est_wholesale_YEARS_1977_1999_T0_1985_B_1000_p_1.csv"))[1, ]
est_cce_wholesale = fread(here("data/cce_est_wholesale.csv"))
est_pca_wholesale = fread(here("data/pca_est_wholesale.csv"))

format_data = function(est, estimator, outcome) {
  est = suppressWarnings(melt(est))
  est$estimator = estimator
  est$outcome = outcome
  est$term = as.numeric(gsub("tau", "", est$variable))

  est$x = est$term
  if (estimator == "Quasi-Long Differencing") {
    est$x = est$x + 0.2
  } else if (estimator == "Common Correlated Effects") {
    est$x = est$x - 0.2
  }

  return(est)
}

ests_retail = rbindlist(list(
  format_data(est_qld_retail, "Quasi-Long Differencing", "Retail"), 
  format_data(est_cce_retail, "Common Correlated Effects", "Retail"), 
  format_data(est_pca_retail, "Principal Components", "Retail")
))
ests_wholesale = rbindlist(list(
  format_data(est_qld_wholesale, "Quasi-Long Differencing", "Wholesale"), 
  format_data(est_cce_wholesale, "Common Correlated Effects", "Wholesale"), 
  format_data(est_pca_wholesale, "Principal Components", "Wholesale")
))


# Plot estimators --------------------------------------------------------------

(plot_retail = ggplot(ests_retail) +
  geom_point(
    aes(
      x = x, y = value, color = estimator
    ),
    size = 2.5
  ) + 
  geom_line(
    aes(
      x = x, y = value, color = estimator,
      group = estimator
    ),
    linewidth = 1.2
  ) +
  scale_y_continuous(limits = c(-0.05, 0.125)) +
  scale_color_manual(
    values = c("Quasi-Long Differencing" = "black", 
               "Common Correlated Effects" = "#9a2515", 
               "Principal Components" = "#5601A4"),
  ) + 
  guides(
    color = guide_legend(override.aes = list(size = 2.5, linewidth = 1.2))
  ) +
  labs(
    x = "Event Time", y = "Coefficient", color = "Factor Estimator"
  ) +
  kfbmisc::theme_kyle(base_size = 16) + 
  theme(
    legend.position = c(0.25, 0.75),
    legend.background = element_rect(fill = "white"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  ))

ggsave(
  here("figures/retail_many_estimators.pdf"), 
  plot_retail, width = 12, height = 8
)

(plot_wholesale = ggplot(ests_wholesale) +
  geom_point(
    aes(
      x = x, y = value, color = estimator
    ),
    size = 2.5
  ) + 
  geom_line(
    aes(
      x = x, y = value, color = estimator,
      group = estimator
    ),
    linewidth = 1.2
  ) +
  scale_y_continuous(limits = c(-0.2, 0.2)) +
  scale_color_manual(
    values = c("Quasi-Long Differencing" = "black", 
               "Common Correlated Effects" = "#9a2515", 
               "Principal Components" = "#5601A4"),
  ) + 
  guides(
    color = guide_legend(override.aes = list(size = 2.5, linewidth = 1.2))
  ) +
  labs(
    x = "Event Time", y = "Coefficient", color = "Factor Estimator"
  ) +
  kfbmisc::theme_kyle(base_size = 16) + 
  theme(
    legend.position = c(0.25, 0.75),
    legend.background = element_rect(fill = "white"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  ))

ggsave(
  here("figures/wholesale_many_estimators.pdf"), 
  plot_wholesale, width = 12, height = 8
)




