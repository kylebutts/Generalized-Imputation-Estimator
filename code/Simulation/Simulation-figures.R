library(data.table)
library(ggplot2)
library(here)
library(fixest)
library(tikzDevice)

source(here("figures/convert_tex_to_pdf.R"))

# Load Estimates ---------------------------------------------------------------

sim2 = fread(here("estimates/simulation-2.csv"))

sim2_points = rbind(
  data.table(
    model = r'(TWFE Imputation with $w_i \beta_t$)', 
    signal_to_noise = sim2$signal_to_noise,
    estimate = sim2$did2s_cov_bias_tau8, 
    std_error = sim2$did2s_cov_bias_sd_tau8
  ), 
  data.table(
    model = r'(Factor Model Imputation)', 
    signal_to_noise = sim2$signal_to_noise,
    estimate = sim2$generalized_bias_tau8, 
    std_error = sim2$generalized_bias_sd_tau8
  )
)

# Figures ----------------------------------------------------------------------

## Signal to Noise Ratio -------------------------------------------------------

(plot_signal_to_noise = ggplot(
    aes(
      x = signal_to_noise, y = estimate, 
      color = model, shape = model, 
      fill = model, group = model
    ),
    data = sim2_points
  ) +
  geom_point(size = 3) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(
    aes(
      ymin = estimate - 1.96 * std_error, 
      ymax = estimate + 1.96 * std_error
    ),
    alpha = 0.2, show.legend = FALSE
  ) +
  labs(
    x = "Signal to Noise Ratio",
    y = "Average Bias",
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
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1.2)),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.2),
    legend.spacing.y = unit(10, "pt"),
    legend.background = element_rect(colour = "black", linewidth = 0.5),
    legend.margin = margin(t = 2, r = 8, b = 8, l = 8, unit = "pt")
  ))


# Export -----------------------------------------------------------------------
tikzDevice::tikz(
  here("figures/simulation-bias_signal_to_noise.tex"),
  width = 10, height = 5, standAlone = FALSE
)
plot(plot_signal_to_noise)
dev.off()
compile_tikz(here("figures/simulation-bias_signal_to_noise.tex"))
