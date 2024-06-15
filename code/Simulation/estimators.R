# %%
quietly <- function(code) {
  warnings <- character()
  wHandler <- function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }

  messages <- character()
  mHandler <- function(m) {
    messages <<- c(messages, conditionMessage(m))
    invokeRestart("muffleMessage")
  }

  temp <- file()
  sink(temp)
  on.exit({
    sink()
    close(temp)
  })

  result <- withCallingHandlers(
    code,
    warning = wHandler,
    message = mHandler
  )
  
  return(result)
}


# %%
est_twfe <- function(df) {
  fs <- feols(
    y ~ 0 | id + t,
    df |> filter(treat == 0)
  )
  df$ytilde <- df$y - predict(fs, newdata = df)
  ss <- feols(
    ytilde ~ i(rel_year, ref = -10), df
  )

  est <- coef(ss, keep = "rel_year")
  est <- est[c("rel_year::0", "rel_year::1", "rel_year::2")]
  return(est)
}

# %%
est_twfe_covs <- function(df) {
  fs <- feols(
    y ~ 0 + i(t, W1) + i(t, W2) | id + t,
    df |> filter(treat == 0),
    warn = FALSE, notes = FALSE
  )
  df$ytilde <- df$y - predict(fs, newdata = df)
  ss <- feols(
    ytilde ~ i(rel_year, ref = -10), df
  )

  est <- coef(ss, keep = "rel_year")
  as.numeric(str_replace(names(est), "rel_year::", ""))

  est <- est[c("rel_year::0", "rel_year::1", "rel_year::2")]
  return(est)
}

# %%
library(augsynth)
est_synth <- function(df) {
  T0 = df |> filter(post == FALSE) |> with(max(t))

  est <- augsynth::augsynth(
    y ~ treat,
    unit = id, time = t, data = df,
    progfunc = "None", scm = T
  )

  estimates <- summary(est, inf = FALSE)$att
  estimates$Estimate[match((T0 + 1):(T0 + 3), estimates$Time)]
}

# %%
library(augsynth)
est_augsynth <- function(df) {
  T0 = df |> filter(post == FALSE) |> with(max(t))

  est <- augsynth::augsynth(
    y ~ treat,
    unit = id, time = t, data = df,
    progfunc = "Ridge"
  )

  estimates <- summary(est, inf = FALSE)$att
  estimates$Estimate[match((T0 + 1):(T0 + 3), estimates$Time)]
}

# %% 
est_matrix_completion <- function(df) {
  T0 = df |> filter(post == FALSE) |> with(max(t))

  est <- augsynth::augsynth(
    y ~ treat,
    unit = id, time = t, data = df,
    progfunc = "MCP"
  )

  estimates <- summary(est, inf = FALSE)$att
  estimates$Estimate[match((T0 + 1):(T0 + 3), estimates$Time)]
}

# %%
# library(gsynth)
source(here("code/Simulation/gsynth_default_patch.R"))
est_gsynth <- function(df, force = "none", p = c(0L, 3L)) {
  T0 = df |> filter(post == FALSE) |> with(max(t))

  est <- patched_gsynth(
    Y = "y", D = "treat", index = c("id", "year"),
    data = df |> rename(year = t),
    force = force, se = FALSE,
    r = p, CV = ifelse(length(p) == 1, FALSE, TRUE),
    min.T0 = 3
  )$att  
  est[match((T0 + 1):(T0 + 3), names(est))]
}


# %%
compute_within_transform <- function(df) {
  T0 = df |> filter(post == FALSE) |> with(max(t))

  df |>
    mutate(ybar_inf_t = mean(y[treated == 0]), .by = t) |>
    mutate(ybar_i_pre = mean(y[t <= .env$T0]), .by = id) |>
    mutate(ybar_inf_pre = mean(y[treated == 0 & t <= .env$T0])) |>
    mutate(ytilde = y - ybar_inf_t - ybar_i_pre + ybar_inf_pre) |>
    pull(ytilde)
}

# %%
est_qld_F_known <- function(df, within_transform = FALSE) {
  T0 = df |> filter(post == FALSE) |> with(max(t))

  if (within_transform == TRUE) {
    df$y <- compute_within_transform(df, T0 = T0)
  }

  F <- with(df |> filter(id == 1), cbind(F1, F2))
  Fpre <- F[1:T0, ]

  imputation_mat <- F %*% solve(crossprod(Fpre), t(Fpre))

  # Y_mat = sparse_mat_from_ijx(df$t, df$id, df$y)
  df <- df |>
    arrange(id, t) |>
    mutate(
      te_hat = y - as.numeric(.env$imputation_mat %*% y[t <= .env$T0]),
      .by = id
    )

  est <- feols(
    te_hat ~ i(rel_year, ref = -10),
    data = df
  )
  est |>
    coef(keep = "rel_year") |>
    _[c("rel_year::0", "rel_year::1", "rel_year::2")]
}

# %% 
# Estimate p with `p = NULL`, otherwise give int
library(JuliaCall)
julia_setup()
julia_source(here("code/Simulation/estimate_qld.jl"))

est_qld <- function(df, within_transform = FALSE, p = -1L) {
  T0 = df |> filter(post == FALSE) |> with(max(t))

  if (within_transform == TRUE) {
    df$y <- compute_within_transform(df, T0 = T0)
  }

  julia_assign("df", df)
  julia_assign("T0", T0)
  julia_assign("p", as.integer(p))
  F <- julia_eval("est_F_qld(df; T0 = T0, p = p)")
  Fpre <- F[1:T0, ]

  imputation_mat <- F %*% solve(crossprod(Fpre), t(Fpre))

  df <- df |>
    arrange(id, t) |>
    mutate(
      te_hat = y - as.numeric(.env$imputation_mat %*% y[t <= .env$T0]),
      .by = id
    )

  est <- feols(
    te_hat ~ i(rel_year, ref = -10),
    data = df
  )
  est |>
    coef(keep = "rel_year") |>
    _[c("rel_year::0", "rel_year::1", "rel_year::2")]
}

