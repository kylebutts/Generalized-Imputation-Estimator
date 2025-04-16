# %%
library(tidyverse)
library(JuliaCall)
julia_setup()
julia_source("code/qld/QLD.jl")
julia_source("code/qld/gmm_qld.jl")
julia_source("code/qld/attgt.jl")
julia_source("code/qld/within_transform.jl")
julia_source("code/qld/qld_imputation.jl")

mboot <- function(IF, B = 1000) {
  # Following https://github.com/bcallaway11/did/blob/master/R/mboot.R
  # IF = N x K matrix of influence function
  N = nrow(IF)

  # Rademacher multiplier bootstrap
  bs_ests <- do.call(
    "rbind",
    purrr::map(1:B, function(i) {
      rad <- 1 - 2 * (runif(N) > 0.5)
      t(t(IF) %*% rad)
    })
  )

  # Bootstrap standard error
  # se_bs <- sqrt(diag(var(bs_ests)))
  se_bs <- apply(
    bs_ests,
    2,
    function(b)
      (quantile(b, .75, type = 1, na.rm = T) -
        quantile(b, .25, type = 1, na.rm = T)) /
        (qnorm(.75) - qnorm(.25))
  )

  # Each row divided by se
  abs_t_stats <- abs(t(t(bs_ests) / se_bs))

  # Max t-stat for each bootstrap
  max_abs_t_stat_per_draw <- apply(abs_t_stats, 1, max)

  crit_val <- quantile(max_abs_t_stat_per_draw, 0.95)

  return(list(se = se_bs, crit_val = crit_val))
}

# %%
df <- read_csv("data/Simulations/df_ex_dgp_3.csv", show_col_types = FALSE)
qld_est_analytic <- julia_call(
  "qld_imputation",
  df,
  y = "y",
  id = "id",
  t = "t",
  g = "g",
  W = c("W1", "W2"),
  do_within_transform = FALSE,
  p = -1L,
  type = "gt",
  vcov_type = "analytic"
)
qld_est_analytic$vcov |> diag() |> sqrt()

qld_est <- julia_call(
  "qld_imputation",
  df,
  y = "y",
  id = "id",
  t = "t",
  g = "g",
  W = c("W1", "W2"),
  do_within_transform = FALSE,
  p = -1L,
  type = "gt",
  vcov_type = "uniform"
)
IF <- t(qld_est$inf_func)
N <- nrow(IF)
mboot(1 / sqrt(N) * IF, B = 10000)

qld_est$se

# %%
# %%
cbp <- read_csv(
  "data/County_Business_Patterns/sample_basker_YEARS_1977_1999_T0_1985.csv",
  show_col_types = FALSE
)
qld_est_analytic <- julia_call(
  "qld_imputation",
  cbp,
  y = "log_retail_emp",
  id = "fips",
  t = "year",
  g = "g",
  W = c(
    "share_pop_ind_manuf",
    "share_pop_poverty_78_below",
    "share_pop_poverty_78_above",
    "share_pop_emp_private",
    "share_pop_emp_government",
    "share_school_col",
    "share_school_hs"
  ),
  do_within_transform = TRUE,
  p = -1L,
  type = "dynamic",
  vcov_type = "analytic"
)
qld_est_analytic$vcov |> diag() |> sqrt()

qld_est <- julia_call(
  "qld_imputation",
  cbp,
  y = "log_retail_emp",
  id = "fips",
  t = "year",
  g = "g",
  W = c(
    "share_pop_ind_manuf",
    "share_pop_poverty_78_below",
    "share_pop_poverty_78_above",
    "share_pop_emp_private",
    "share_pop_emp_government",
    "share_school_col",
    "share_school_hs"
  ),
  do_within_transform = TRUE,
  p = -1L,
  type = "dynamic",
  vcov_type = "uniform"
)
IF <- t(qld_est$inf_func)
N <- nrow(IF)
se_uniform <- mboot(1 / sqrt(N) * IF, B = 10000)$se
crit_val <- mboot(1 / sqrt(N) * IF, B = 10000)$crit_val
qld_est_analytic$vcov |> diag() |> sqrt()

ggplot() +
  geom_point(aes(
    x = qld_est$rel_year,
    y = 1.96 * qld_est_analytic$vcov |> diag() |> sqrt(),
    color = "Analytic SE * 1.96"
  )) +
  geom_point(aes(
    x = qld_est$rel_year,
    y = 1.96 * se_uniform,
    color = "Uniform SE * 1.96"
  )) +
  geom_point(aes(
    x = qld_est$rel_year,
    y = crit_val * se_uniform,
    color = "Uniform SE * crit_val"
  )) +
  kfbmisc::theme_kyle(base_size = 14, legend = "top")
