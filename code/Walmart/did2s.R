library(tidyverse)
library(did2s)
library(glue)
library(here)
library(broom)

df <- read_csv(
  here("data/County_Business_Patterns/sample_basker_YEARS_1977_1999_T0_1985.csv"), 
  show_col_types = FALSE
)

process_did2s <- function(res) {
  broom::tidy(res) |>
    filter(str_detect(term, "rel_year::")) |>
    mutate(
      rel_year = term |>
        str_replace("rel_year::", "") |>
        as.numeric(),
    ) |>
    select(rel_year, estimate, std_error = std.error) |>
    mutate(
      lower = estimate - 1.96 * std_error,
      upper = estimate + 1.96 * std_error
    )
}

res_retail <- did2s::did2s(
  data = df,
  yname = "log_retail_emp",
  treatment = "any_open",
  cluster_var = "fips",
  first_stage = ~ 0 | fips + year,
  second_stage = ~ i(rel_year, ref = c(-Inf)),
  verbose = FALSE
)
# coefplot(res_retail)

res_covs_retail <- did2s::did2s(
  data = df,
  yname = "log_retail_emp",
  treatment = "any_open",
  cluster_var = "fips",
  first_stage = ~ i(year, share_pop_ind_manuf) + i(year, share_pop_poverty_78_below) + i(year, share_pop_poverty_78_above) + i(year, share_pop_emp_private) + i(year, share_pop_emp_government) + i(year, share_school_col) + i(year, share_school_hs) | fips + year,
  second_stage = ~ i(rel_year, ref = c(-Inf)),
  verbose = FALSE
)
# coefplot(res_covs_retail)

res_wholesale <- did2s::did2s(
  data = df,
  yname = "log_wholesale_emp",
  treatment = "any_open",
  cluster_var = "fips",
  first_stage = ~ 0 | fips + year,
  second_stage = ~ i(rel_year, ref = c(-Inf)),
  verbose = FALSE
)
# coefplot(res_wholesale)

res_covs_wholesale <- did2s::did2s(
  data = df,
  yname = "log_wholesale_emp",
  treatment = "any_open",
  cluster_var = "fips",
  first_stage = ~ i(year, share_pop_ind_manuf) + i(year, share_pop_poverty_78_below) + i(year, share_pop_poverty_78_above) + i(year, share_pop_emp_private) + i(year, share_pop_emp_government) + i(year, share_school_col) + i(year, share_school_hs) | fips + year,
  second_stage = ~ i(rel_year, ref = c(-Inf)),
  verbose = FALSE
)
# coefplot(res_covs_wholesale)


## Export
write_csv(
  process_did2s(res_retail),
  here(glue("estimates/did2s_est_retail.csv"))
)
write_csv(
  process_did2s(res_covs_retail),
  here(glue("estimates/did2s_est_covs_retail.csv"))
)
write_csv(
  process_did2s(res_wholesale),
  here(glue("estimates/did2s_est_wholesale.csv"))
)
write_csv(
  process_did2s(res_covs_wholesale),
  here(glue("estimates/did2s_est_covs_wholesale.csv"))
)

