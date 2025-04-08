extract_tt_latex_body <- function(x) {
  stopifnot(inherits(x, "tinytable"))
  paste0(
    tinytable:::build_tt(x, "latex")@body,
    collapse = "\n"
  )
}
