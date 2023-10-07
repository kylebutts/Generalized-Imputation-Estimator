#' Double-demeaning procedure (while preserving factor structure)
#' 
#' @param X a N-by-K matrix of variables
#' @param i a N-vector of individual identifiers
#' @param t a N-vector of time identifiers
#' @param d a N-vector of 0/1s indicating when treatment is on
#' 
#' @return The matrix X double-demeaned
double_demean <- function(X, i, t, d) {

  if (!all(d %in% c(TRUE, FALSE, 0, 1))) {
    stop("d must be a vector of 0/1s or a vector of TRUE/FALSEs")
  } 
  
  dt = as.data.table(X)
  vars = colnames(dt)
  dt$i = i
  dt$t = as.numeric(t)
  dt$d = as.numeric(d)

  T0 = dt[d == 1, min(t)] - 1
  
  dt[, 
    g := ifelse(
      any(d == 1), 
      min(t[d == 1]), 
      Inf
    ), 
    by = i
  ]
  
  for (var in vars) {
    
    dt$y = dt[[var]]
    
    # never-treated cross-sectional averages of y
    dt[, 
      yinf_t := mean(y[g == Inf]), 
      by = t
    ]
    # unit-averages pre-T0
    dt[, 
      yi_pre := mean(y[t <= T0]), 
      by = i
    ]
    # never-treated average pre-T0
    dt[, 
      yinf_pre := mean(y[g == Inf & t <= T0])
    ]
    # Create ytilde
    dt[[paste0(var, "_tilde")]] = 
      dt$y - dt$yinf_t - dt$yi_pre + dt$yinf_pre
  }

  as.matrix(
    dt[, .SD, .SDcols = paste0(vars, "_tilde")]
  )
}

#' Takes T-by-r estimate of the estimated factors and imputes y
#' 
#' @param F a T-by-r matrix of estimated factors. Ordered by time
#' @param dt a data.table of the original data. Will be modified in place
#' @param y String. Variable name denoting the outcome variable. Can be 
#'  double-demeaned
#' @param i String. Variable name denoting the individual identifiers
#' @param t String. Variable name denoting the time identifiers
#' @param d String. Variable name denoting the 0/1s indicating when treatment is on
factor_to_imputation <- function(F, dt, y, i, t, d, balanced = TRUE) {

  if (!all(d %in% c(TRUE, FALSE, 0, 1))) {
    stop("d must be a vector of 0/1s or a vector of TRUE/FALSEs")
  } 
  
  i = "fips"
  t = "year"
  d = "any_open"
  y = "log_retail_emp"

  # Need to sort!
  dt = as.data.table(dt)
  setorderv(dt, c(i, t))

  T0 = min(dt[dt[[d]] == 1, ][[t]]) - 1

  g = dt[, 
    g := ifelse(
      any(d == 1), 
      min(t[d == 1]), 
      Inf
    ), 
    by = i,
    env = list(i = i, d = d, t = t)
  ][, V1]


  # Group by group to be efficient
  if (balanced == TRUE) {

    dt$y0_hat = 0
    for (g_id in setdiff(unique(g), Inf)) {
      N_periods_pre <- sum(unique(dt$t) < g_id)
      Fpre = F[1:N_periods_pre, , drop = FALSE]
      P_F_Fpre = F %*% MASS::ginv(Matrix::crossprod(Fpre)) %*% Matrix::t(Fpre)

      dt[,
        y0_hat := ifelse(
          g == g_id, 
          P_F_Fpre %*% y[t < g_id],
          y0_hat
        ),
        by = i
      ]
    }

  }

  return(dt)
}




