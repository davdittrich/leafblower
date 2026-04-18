#' Kish design effect (1-argument) or Henry-Valliant (4-argument)
#'
#' When called with only \code{weights}, computes the Kish (1992) design effect:
#' \code{n * sum(w^2) / sum(w)^2}.
#'
#' When called with all four arguments, computes the Henry & Valliant (2015)
#' calibration design effect: \code{var_weighted(y) / var_unweighted(y)}.
#'
#' @param weights Numeric vector of calibrated weights.
#' @param outcome Numeric outcome vector (optional; 4-arg form only).
#' @param data Data frame used in calibration (optional; 4-arg form only).
#' @param target Named list of target proportions (optional; 4-arg form only).
#' @return Numeric scalar design effect.
#' @export
design_effect <- function(weights, outcome = NULL, data = NULL, target = NULL) {
  if (is.null(outcome)) {
    n <- length(weights)
    return(n * sum(weights^2) / sum(weights)^2)
  }
  n <- length(weights)
  if (length(outcome) != n)
    stop("design_effect: 'outcome' length must equal length(weights)")
  y_bar_w <- sum(weights * outcome) / sum(weights)
  var_w   <- sum(weights * (outcome - y_bar_w)^2) / sum(weights)
  var_u   <- var(outcome)
  if (var_u < 1e-20 || var_w < 1e-20) return(1.0)
  var_w / var_u
}

#' Effective sample size
#'
#' @param weights Numeric vector of calibrated weights.
#' @return Numeric scalar: n / design_effect(weights).
#' @export
effective_sample_size <- function(weights) {
  length(weights) / design_effect(weights)
}
