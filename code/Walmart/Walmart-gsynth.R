library(data.table)
library(did2s)
library(glue)
library(here)
library(gsynth)

source(here("code/Walmart/helpers.R"))
options(max.print = 1000)


outcome <- "wholesale"
# outcome <- "retail"
sample <- fread(here("data/County_Business_Patterns/sample_basker_YEARS_1977_1999_T0_1985.csv"))

# Number of factors
# r = 3
cat(paste0("Outcome variable: log_", outcome, "_emp"))

setorder(sample, fips, year)

T0 <- sample[, min(g)] - 1

sample$y <- sample[[glue("log_{outcome}_emp")]]
sample$any_open = as.numeric(sample$any_open)

sample$share_pop_ind_manuf = sample$share_pop_ind_manuf_nondurable + sample$share_pop_ind_manuf_durable

# Double-demean y
sample$y_tilde <- double_demean(sample[, "y"], sample$fips, sample$year, sample$any_open)
sample[, .(fips, year, any_open, y_tilde)]

# Projected principal components
for (t in unique(sample$year)) {
  est = feols(
    y_tilde ~ share_pop_ind_manuf + share_pop_poverty_78_below + share_pop_poverty_78_above + share_pop_emp_private + share_pop_emp_government + share_school_col + share_school_hs, 
    sample[year == t & g == Inf, ]
  )
  sample[year == t, "y_tilde_phat"] = predict(est, newdata = sample[year == t, ])
}

# out_gsynth_twoway <- gsynth(
#   y ~ any_open, data = sample, 
#   index = c("fips", "year"), force = "two-way", 
#   se = TRUE,
#   CV = TRUE, r = c(0, 5)
# )
# plot(out_gsynth_twoway)

# GSYNTH (MSE to select number of factors) -------------------------------------

out_gsynth <- gsynth(
  y_tilde ~ any_open, data = sample, 
  index = c("fips", "year"), force = "none", 
  se = TRUE, nboots = 50,
  CV = TRUE, r = c(0, 5)# , 
)
plot(out_gsynth)
# out_gsynth$att

# Number of factors
r = out_gsynth$r.cv

# Bai (2009) and Chan and Kwok (2022)
untreated_ids = sample[g == Inf, unique(fips)]
Nco <- length(untreated_ids)
TT <- length(unique(sample$year))
ymat = matrix(nrow = Nco, ncol = TT)

for (i in seq_along(untreated_ids)) { 
  id = untreated_ids[i]
  ymat[i, ] = sample[fips == id, y_tilde]
}

# gsynth:::panel_factor(ymat, r = r)
# demean by the grand mean 
ymat = ymat - mean(ymat)
tYY = tcrossprod(ymat) / ncol(ymat) * nrow(ymat)
eigen = eigen(tYY)
Fhat = crossprod(ymat, eigen$vectors[, 1:r, drop = FALSE]) / ncol(ymat)

# Data Check:
# Fhat_gsynth = gsynth:::panel_factor(ymat, r = r)$lambda
# Fhat_gsynth - Fhat %*% solve(crossprod(Fhat)) %*% crossprod(Fhat, Fhat_gsynth)
for (fips_id in unique(sample[g < Inf, fips])) {
  idx = sample$fips == fips_id
  last_untreated_year = sample[idx, ][any_open == 0, max(year)]
  
  N_pre <- (last_untreated_year - min(sample$year) + 1)
  Fpre <- Fhat[1:N_pre, ]
  tFpreFpreinv = MASS::ginv(crossprod(Fpre))

  yi_pre <- sample[idx, ][year <= last_untreated_year, y_tilde]

  # impute y(0) using X
  gi_hat <- tFpreFpreinv %*% t(Fpre) %*% (yi_pre)
  y0hat <- (Fhat %*% gi_hat)

  sample[idx, "y0hat"] = y0hat
}

# Difference variables: Z_{it} - \hat{Z}_{it}(0)
sample$tau_hat <- sample$y_tilde - sample$y0hat

est = feols(
  tau_hat ~ 0 + i(rel_year), 
  sample[g < Inf, ]
)  
coefplot(est)

## Export estimates ------------------------------------------------------------

fwrite(
  sample[, c("fips", "year", "y_tilde", "y0hat", "tau_hat")], 
  here(glue("estimates/est_y0_outcome_{outcome}_pca.csv"))
)

# Export to .csv
coef(est)
c(
  paste0(
    names(coef(est)) |> gsub("rel_year::", "tau", x = _), 
    collapse = ","
  ),
  "\n",
  paste0(coef(est), collapse = ",")
) |>
  # cat()
  cat(file = here(paste0("estimates/pca_est_", outcome, ".csv")))


