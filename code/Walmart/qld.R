# %% 
library(tidyverse)
library(JuliaCall)
library(here)
library(glue)
library(broom)

JuliaCall::julia_setup()
JuliaCall::julia_source(here("code/qld/qld_imputation.jl"))
JuliaCall::julia_source(here("code/qld/attgt.jl"))
JuliaCall::julia_source(here("code/qld/gmm_qld.jl"))
JuliaCall::julia_source(here("code/qld/within_transform.jl"))

df <- read_csv(
  "data/County_Business_Patterns/sample_basker_YEARS_1977_1999_T0_1985.csv",
  show_col_types = FALSE
)

#' Fact on median employment
med_retail_emp <- df |> filter(year == 1977) |> with(median(retail_emp))
med_wholesale_emp <- df |> filter(year == 1977) |> with(median(wholesale_emp))

# %% 
tidy.qld_dynamic <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  stopifnot(isTRUE(conf.int) || isFALSE(conf.int))
  stopifnot(conf.level > 0 && conf.level < 1)
  result <- data.frame(
    rel_year = x[[1]],
    estimate = x[[2]],
    std_error = sqrt(diag(x[[3]])),
    std_error_naive = sqrt(diag(x[[6]]))
  )
  if (isTRUE(conf.int)) {
    moe <- abs(qnorm((1 - conf.level) / 2))
    result$lower  <- result$estimate - moe * result$std_error
    result$upper <- result$estimate + moe * result$std_error
  }
  return(result)
}

# %% 
res_retail <- julia_call("qld_imputation",
  df,
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
  return_y0 = TRUE,
  return_naive_se = TRUE
)
class(res_retail) <- c("qld_dynamic", "qld", "list")
tidy_retail <- broom::tidy(res_retail, conf.int = TRUE, conf.level = 0.95)

res_wholesale <- julia_call("qld_imputation",
  df,
  y = "log_wholesale_emp",
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
  p = 1L,
  type = "dynamic",
  return_y0 = TRUE,
  return_naive_se = TRUE
)
class(res_wholesale) <- c("qld_dynamic", "qld", "list")
tidy_wholesale <- broom::tidy(res_wholesale, conf.int = TRUE, conf.level = 0.95)

# %% 
# Check thinks look right
ggplot() + 
  geom_hline(yintercept = 0) +
  geom_point(
    aes(x = rel_year, y = estimate),
    data = tidy_retail
  )+ 
  geom_errorbar(
    aes(x = rel_year, ymin = lower, ymax = upper),
    data = tidy_retail
  )

# %% 
write_csv(
  tidy_retail,
  here("estimates/qld_est_retail.csv")
)
write_csv(
  tidy_wholesale,
  here("estimates/qld_est_wholesale.csv")
)
write_csv(
  res_retail[[5]], 
  here("estimates/est_y0_outcome_retail_qld.csv")
)
write_csv(
  res_wholesale[[5]], 
  here("estimates/est_y0_outcome_wholesale_qld.csv")
)

