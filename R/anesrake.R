#' anesrake compatibility wrapper
#'
#' Maps anesrake parameter names to harvest(). Unknown params produce warnings.
#' @param inputter Data frame (maps to \code{data}).
#' @param targets Named list of target proportions.
#' @param weightvec Starting weights or NULL.
#' @param caseid Ignored (row label not used by calibration).
#' @param pctlim Deprecated; maps to \code{convergence["pct"]}.
#' @param cap Upper weight cap; maps to \code{max_weight}.
#' @param choosemethod Algorithm; \code{"rake"} and \code{"nrake"} are silently
#'   remapped to \code{"ieppa"} (the recommended method).
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

  # F6: map pctlim → convergence[["pct"]] instead of silently dropping it
  convergence <- list()
  if (!missing(pctlim) && !is.null(pctlim)) {
    warning("anesrake: 'pctlim' is deprecated; use convergence = list(pct = pctlim)")
    convergence[["pct"]] <- pctlim
  }

  # F7: silently remap legacy "rake"/"nrake" to "ieppa" — avoids deprecation warning.
  # tolower() guards against caller-passed "Rake"/"NRAKE" etc.
  if (tolower(choosemethod) %in% c("rake", "nrake")) choosemethod <- "ieppa"

  harvest(
    data           = inputter,
    target         = targets,
    start_weights  = weightvec,
    max_weight     = cap,
    method         = choosemethod,
    max_iterations = as.integer(nlim),
    convergence    = convergence
  )
}
