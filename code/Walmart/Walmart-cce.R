library(data.table)
library(did2s)
library(glue)
library(here)
library(gsynth)


outcome <- "wholesale"
# outcome <- "retail"
sample <- fread(here("data/County_Business_Patterns/sample_basker_YEARS_1977_1999_T0_1985.csv"))

sample$y <- sample[[glue("log_{outcome}_emp")]]
sample$share_pop_ind_manuf <- sample$share_pop_ind_manuf_durable + sample$share_pop_ind_manuf_nondurable
sample$any_open = as.numeric(sample$any_open)

xvars = c("log_manufacturing_emp", "log_construction_emp", "log_agriculture_emp", "log_healthcare_emp")

## CCE pooled estimate of \beta hat --------------------------------------------

T0 <- 1985

# CSAs of y and Xs
Fhat <- sample |>
  DT(
    g == Inf,
    lapply(.SD, \(x) mean(x)),
    by = year,
    .SDcols = c("y", xvars)
  ) |>
  DT(,
    as.matrix(.SD),
    .SDcols = c("y", xvars)
  )

# Add constant to Fhat
Fhat = cbind(Fhat, rep(1, nrow(Fhat)))

# CCEP estimator for \hat{\beta} using pre-treatment X's for all groups
N_pre_T0 <- (T0 - min(sample$year) + 1)
Fpre <- Fhat[1:N_pre_T0, ]
tFpreFpreinv = MASS::ginv(crossprod(Fpre))
M_Fpre <- diag(N_pre_T0) - Fpre %*% tFpreFpreinv %*% t(Fpre)

B <- matrix(0, nrow = length(xvars), ncol = length(xvars))
A <- matrix(0, nrow = length(xvars), ncol = 1)

for (id in unique(sample[, fips])) {
  Xi_pre <- sample |>
    DT(
      fips == id & year <= T0,
      as.matrix(.SD),
      .SDcols = xvars
    )

  yi_pre <- sample[
    fips == id & year <= T0,
    y
  ]

  B <- B + t(Xi_pre) %*% M_Fpre %*% Xi_pre
  A <- A + t(Xi_pre) %*% M_Fpre %*% yi_pre
}

bhat <- MASS::ginv(B) %*% A

## Imputation of covariates and outcome ----------------------------------------
for (fips_id in unique(sample[g < Inf, fips])) {

  idx = sample$fips == fips_id

  t <- sample[idx, ]$year

  Xi <- sample[idx, ] |>
    _[
      ,
      as.matrix(.SD),
      .SDcols = xvars
    ]

  Xi_pre <- Xi[t <= T0, , drop = FALSE]
  yi_pre <- sample[idx, ][year <= T0, y]

  # impute y(0) using X
  gi_hat <- tFpreFpreinv %*% t(Fpre) %*% (yi_pre - Xi_pre %*% bhat)
  y0hat <- (Xi %*% bhat) + (Fhat %*% gi_hat)

  sample[idx, "y0hat"] = y0hat
}

# Difference variables: Z_{it} - \hat{Z}_{it}(0)
sample$tau_hat <- sample$y - sample$y0hat

est = feols(
  tau_hat ~ 0 + i(rel_year), 
  sample[g < Inf, ]
)  
coefplot(est)

## Export estimates ------------------------------------------------------------

fwrite(
  sample[, c("fips", "year", "y", "y0hat", "tau_hat")], 
  here(glue("estimates/est_y0_outcome_{outcome}_cce.csv"))
)

c(
  paste0(
    names(coef(est)) |> gsub("rel_year::", "tau", x = _), 
    collapse = ","
  ),
  "\n",
  paste0(coef(est), collapse = ",")
) |>
  # cat()
  cat(file = here(paste0("estimates/cce_est_", outcome, ".csv")))






