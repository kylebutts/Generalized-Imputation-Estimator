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
library(did2s)
est_twfe <- function(df) {
  # if (do_inference == TRUE) { 
  #   
  # } else {
  #   fs <- feols(
  #     y ~ 0 | id + t,
  #     df |> dplyr::filter(treat == 0)
  #   )
  #   df$ytilde <- df$y - predict(fs, newdata = df)
  #   ss <- feols(
  #     ytilde ~ i(rel_year, ref = -10), df
  #   )
  #   est <- broom::tidy(ss) |> 
  #     filter(str_detect(term, "rel_year::")) |> 
  #     mutate(rel_year = as.numeric(str_replace(term, "rel_year::", ""))) |>
  #     filter(rel_year >= 0) 
  #     select(rel_year, estimate)
  # }
  ss <- did2s::did2s(
    data = df, 
    yname = "y",
    first_stage = ~ 0 | id + t,
    second_stage = ~ i(rel_year, ref = -10),
    treatment = "treat",
    cluster_var = "id",
    verbose = FALSE
  )
  est <- broom::tidy(ss) |> 
    filter(str_detect(term, "rel_year::")) |> 
    mutate(rel_year = as.numeric(str_replace(term, "rel_year::", ""))) |>
    filter(rel_year >= 0) |> 
    mutate(
      ci_lower = estimate - 1.96 * std.error,
      ci_upper = estimate + 1.96 * std.error
    ) |> 
    select(rel_year, estimate, std.error, ci_lower, ci_upper)
  
  return(est)
}

# %%
est_twfe_covs <- function(df) {
  # if (do_inference == TRUE) { 
  #   
  # } else {
  #   fs <- feols(
  #     y ~ 0 + i(t, W1) + i(t, W2) | id + t,
  #     df |> filter(treat == 0),
  #     warn = FALSE, notes = FALSE
  #   )
  #   df$ytilde <- df$y - predict(fs, newdata = df)
  #   ss <- feols(
  #     ytilde ~ i(rel_year, ref = -10), df
  #   )  
  #   est <- broom::tidy(ss) |> 
  #     filter(str_detect(term, "rel_year::")) |> 
  #     mutate(rel_year = as.numeric(str_replace(term, "rel_year::", ""))) |>
  #     filter(rel_year >= 0) 
  #     select(rel_year, estimate)
  # }
  ss <- did2s::did2s(
    data = df, 
    yname = "y",
    first_stage = ~ 0 + i(t, W1) + i(t, W2) | id + t,
    second_stage = ~ i(rel_year, ref = -10),
    treatment = "treat",
    cluster_var = "id",
    verbose = FALSE
  )
  est <- broom::tidy(ss) |> 
    filter(str_detect(term, "rel_year::")) |> 
    mutate(rel_year = as.numeric(str_replace(term, "rel_year::", ""))) |>
    filter(rel_year >= 0) |> 
    mutate(
      ci_lower = estimate - 1.96 * std.error,
      ci_upper = estimate + 1.96 * std.error
    ) |> 
    select(rel_year, estimate, std.error, ci_lower, ci_upper)
  
  return(est)
}

# %%
library(augsynth)
est_synth <- function(df, do_inference = FALSE) {
  T0 = df |> filter(post == FALSE) |> with(max(t))

  synth_est <- quietly({
    augsynth::augsynth(
      y ~ treat,
      unit = id, time = t, data = df,
      progfunc = "None", scm = T
    )
  })

  # summ <- summary(synth_est, inf = do_inference, inf_type = "conformal")$att
  summ <- summary(synth_est, inf = do_inference, inf_type = "jackknife+")$att
  if (do_inference == TRUE) {
    est <- summ |>
      as_tibble() |> 
      mutate(rel_year = Time - (T0 + 1)) |>
      filter(rel_year >= 0) |> 
      select(rel_year, estimate = Estimate, ci_lower = lower_bound, ci_upper = upper_bound) |> 
      mutate(std.error = (ci_upper - estimate) / 1.96, .after = "estimate")
  } else {
    est <- summ |>
      as_tibble() |> 
      mutate(rel_year = Time - (T0 + 1)) |>
      filter(rel_year >= 0) |> 
      select(rel_year, estimate = Estimate)
  }

  return(est)
  }

# %%
library(augsynth)
est_augsynth <- function(df, do_inference = FALSE) {
  T0 = df |> filter(post == FALSE) |> with(max(t))

  augsynth_est <- quietly({
    augsynth::augsynth(
      y ~ treat,
      unit = id, time = t, data = df,
      progfunc = "Ridge"
    )
  })

  # summ <- summary(augsynth_est, inf = do_inference, inf_type = "conformal")$att
  summ <- summary(augsynth_est, inf = do_inference, inf_type = "jackknife+")$att

  if (do_inference == TRUE) {
    est <- summ |>
      as_tibble() |> 
      mutate(rel_year = Time - (T0 + 1)) |>
      filter(rel_year >= 0) |> 
      select(rel_year, estimate = Estimate, ci_lower = lower_bound, ci_upper = upper_bound) |> 
      mutate(std.error = (ci_upper - estimate) / 1.96, .after = "estimate")
  } else {
    est <- summ |>
      as_tibble() |> 
      mutate(rel_year = Time - (T0 + 1)) |>
      filter(rel_year >= 0) |> 
      select(rel_year, estimate = Estimate)
  }

  return(est)
  }

# %% 
# WARNING: Very slow to conduct inference
library(augsynth)
est_matrix_completion <- function(df, do_inference = FALSE) {
  T0 = df |> filter(post == FALSE) |> with(max(t))

  MCP_est <- quietly({
    augsynth::augsynth(
      y ~ treat,
      unit = id, time = t, data = df,
      progfunc = "MCP"
    )
  })

  summ <- summary(MCP_est, inf = do_inference)$att
  if (do_inference == TRUE) {
    est <- summ |>
      as_tibble() |> 
      mutate(rel_year = Time - (T0 + 1)) |>
      filter(rel_year >= 0) |> 
      select(rel_year, estimate = Estimate, ci_lower = lower_bound, ci_upper = upper_bound) |> 
      mutate(std.error = (ci_upper - estimate) / 1.96, .after = "estimate")
  } else {
    est <- summ |>
      as_tibble() |> 
      mutate(rel_year = Time - (T0 + 1)) |>
      filter(rel_year >= 0) |> 
      select(rel_year, estimate = Estimate)
  }

  return(est)
}

# %%
library(gsynth)
library(parallel)
library(doParallel)
library(doRNG)
source(here::here("code/Simulation/gsynth_default_patch.R"))
est_gsynth <- function(df, force = "none", p = c(0L, 3L)) {
  T0 = df |> filter(post == FALSE) |> with(max(t))

  gsynth_est <- quietly({
    patched_gsynth(
      Y = "y", D = "treat", index = c("id", "year"),
      data = df |> rename(year = t),
      force = force,
      r = p, CV = ifelse(length(p) == 1, FALSE, TRUE),
      se = TRUE, parallel = FALSE,
      min.T0 = 3
    )
  })
  
  est <- gsynth_est$est.att |>
    as.data.frame() |>
    as_tibble(rownames = "rel_year") |> 
    mutate(rel_year = as.numeric(rel_year) - 1) |> # "1" is year of treatment
    filter(rel_year >= 0) |>
    mutate(selected_p = gsynth_est$r.cv) |>
    select(rel_year, estimate = ATT, std.error = `S.E.`, ci_lower = `CI.lower`, ci_upper = `CI.upper`, selected_p) 
  
  return(est)
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
est_qld_F_known <- function(df, do_within_transform = FALSE) {
  T0 = df |> filter(post == FALSE) |> with(max(t))

  if (do_within_transform == TRUE) {
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

  ss <- feols(
    te_hat ~ i(rel_year, ref = -10),
    data = df, cluster = ~id
  )
  
  est <- broom::tidy(ss) |> 
    filter(str_detect(term, "rel_year::")) |> 
    mutate(rel_year = as.numeric(str_replace(term, "rel_year::", ""))) |>
    filter(rel_year >= 0) |> 
    mutate(
      ci_lower = estimate - 1.96 * std.error,
      ci_upper = estimate + 1.96 * std.error
    ) |> 
    select(rel_year, estimate, std.error, ci_lower, ci_upper)
  
  return(est)
}

# %% 
# Estimate p with `p = NULL`, otherwise give int
library(JuliaCall)
julia_setup()
julia_source(here("code/qld/QLD.jl"))
julia_source(here("code/qld/qld_imputation.jl"))
julia_source(here("code/qld/attgt.jl"))
julia_source(here("code/qld/gmm_qld.jl"))
julia_source(here("code/qld/within_transform.jl"))

est_qld <- function(df, do_within_transform = FALSE, p = -1L) {
  p <- as.integer(p)
  qld_est <- julia_call("qld_imputation",
    df,
    y = "y",
    id = "id",
    t = "t",
    g = "g",
    W = c("W1", "W2"),
    do_within_transform = do_within_transform,
    p = p,
    type = "dynamic"
  )
  
  est <- tibble(
    rel_year = qld_est[[1]],
    estimate = qld_est[[2]],
    std.error = sqrt(diag(qld_est[[3]]))
  ) |> 
    filter(rel_year >= 0) |> 
    mutate(
      ci_lower = estimate - 1.96 * std.error,
      ci_upper = estimate + 1.96 * std.error
    ) |> 
    mutate(selected_p = qld_est[[4]]) |>
    select(rel_year, estimate, std.error, ci_lower, ci_upper, selected_p)

  return(est)
}

