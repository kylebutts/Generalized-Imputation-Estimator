library(data.table)
library(did2s)
library(glue)
library(here)

# outcome <- "retail"
outcome <- "wholesale"
sample <- fread(here("data/County_Business_Patterns/sample_basker_YEARS_1977_1999_T0_1985.csv"))

sample$y <- sample[[glue("log_{outcome}_emp")]]
sample$share_pop_ind_manuf <- sample$share_pop_ind_manuf_durable + sample$share_pop_ind_manuf_nondurable

res <- did2s::did2s(
  data = as.data.frame(sample),
  yname = "y",
  treatment = "any_open",
  cluster_var = "fips",
  first_stage = ~ 0 | fips + year,
  second_stage = ~ i(rel_year, ref = c(-Inf)),
  verbose = FALSE
)
coefplot(res)

res_covs = did2s::did2s(
  data = as.data.frame(sample),
  yname = "y",
  treatment = "any_open",
  cluster_var = "fips",
  first_stage = ~ i(year, share_pop_ind_manuf) + i(year, share_pop_poverty_78_below) + i(year, share_pop_poverty_78_above) + i(year, share_pop_emp_private) + i(year, share_pop_emp_government) + i(year, share_school_col) + i(year, share_school_hs) | fips + year,
  second_stage = ~ i(rel_year, ref = c(-Inf)),
  verbose = FALSE
)
coefplot(res_covs)

bootstraps <- did2s(
  data = sample,
  yname = "y",
  treatment = "any_open",
  cluster_var = "fips",
  first_stage = ~ 0 | fips + year,
  second_stage = ~ i(rel_year, ref = c(-Inf)),
  bootstrap = TRUE,
  n_bootstraps = 1000,
  return_bootstrap = TRUE,
  verbose = FALSE
)

upper = apply(bootstraps, 2, \(x) quantile(x, 0.975))
lower = apply(bootstraps, 2, \(x) quantile(x, 0.025))
std_error = apply(bootstraps, 2, \(x) sd(x))
estimate = coef(res)

est_did2s = data.table(
  rel_year = -22L:13L, estimate = estimate, std_error = std_error, upper = upper, lower = lower
)

## Export
fwrite(
  est_did2s, 
  file = here(glue("estimates/did2s_est_{outcome}_B_1000.csv"))
)

## Covs export
est_covs = as.data.table(coeftable(res_covs))
est_covs = est_covs |>
  _[, .(
    rel_year = -22L:13L,
    estimate = Estimate,
    std_error = `Std. Error`
  )] |>
  _[, let(
    upper = estimate + 1.96 * std_error,
    lower = estimate - 1.96 * std_error
  )]

fwrite(
  est_covs, 
  file = here(glue("estimates/did2s_est_covs_{outcome}.csv"))
)
