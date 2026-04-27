#' anesrake compatibility wrapper
#'
#' Maps anesrake parameter names to harvest(). Unknown params produce warnings.
#' @param inputter Data frame (maps to \code{data}).
#' @param targets Named list of target proportions.
#' @param weightvec Starting weights or NULL.
#' @param caseid Ignored (row label not used by calibration).
#' @param pctlim Deprecated; maps to \code{convergence["pct"]}.
#' @param cap Upper weight cap; maps to \code{max_weight}.
#' @param choosemethod "rake" or "nrake"; maps to \code{method="lbfgsb"} with warning.
#' @param type Ignored.
#' @param nlim Max iterations; maps to \code{max_iterations}.
#' @param iterate Ignored (always iterates).
#' @param threads Ignored (single-threaded v1).
#' @param ... Unknown params produce a warning then are ignored.
#' @return Same as \code{harvest()}.
#' @export
anesrake <- function(inputter, targets, weightvec = NULL, caseid = NULL,
                     pctlim = 0.05, cap = 5, choosemethod = "rake",
                     type = "pctlim", nlim = 500L, iterate = TRUE,
                     threads = 1L, ...) {
  dots <- list(...)
  if (length(dots) > 0)
    warning("anesrake: ignoring unknown arguments: ",
            paste(names(dots), collapse = ", "))
  if (!is.null(caseid))
    message("anesrake: caseid is ignored (not used by leafblower calibration)")
  if (!is.null(pctlim))
    warning("anesrake: pctlim=", pctlim, " is ignored ",
            "(pct-based tolerance is not implemented in leafblower). ",
            "The default tol_abs=1e-6 is used. ",
            "Pass convergence=list(absolute=...) to harvest() to control tolerance.")
  harvest(
    data           = inputter,
    target         = targets,
    start_weights  = weightvec,
    max_weight     = cap,
    method         = choosemethod,
    max_iterations = as.integer(nlim)
  )
}
