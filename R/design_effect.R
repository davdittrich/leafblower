#' Kish design effect (1-argument) or Cochran (1977) calibration design effect (4-argument)
#'
#' When called with only \code{weights}, computes the Kish (1992) design effect:
#' \code{n * sum(w^2) / sum(w)^2}.
#'
#' When called with all four arguments, computes the Henry & Valliant (2015)
#' calibration design effect ratio \code{var_weighted(y) / var_unweighted(y)},
#' using Cochran (1977) §4.5 unbiased weighted variance for the numerator
#' (denominator \code{(sum(w)^2 - sum(w^2)) / sum(w)}, reduces to Bessel
#' \code{(n-1)} under uniform weights; matches \code{survey::svydesign}).
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
  # Cochran (1977) §4.5: unbiased weighted variance of y under reliability weights.
  # Reduces to Bessel (n-1) denominator when weights are uniform; matches survey::svydesign.
  var_w_denom <- (sum(weights)^2 - sum(weights^2)) / sum(weights)
  if (var_w_denom <= 0) return(1.0)  # degenerate (single observation or all-zero weights)
  var_w   <- sum(weights * (outcome - y_bar_w)^2) / var_w_denom
  var_u   <- var(outcome)  # Bessel (n-1) denominator; Cochran var_w now matches this basis
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
