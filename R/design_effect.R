#' Kish (1965) design effect (1-argument) or Henry & Valliant (2015) calibration design effect (4-argument)
#'
#' When called with only \code{weights}, computes the Kish (1965) weighting design effect
#' \code{deff_K(w) = n * sum(w^2) / sum(w)^2}.
#'
#' When called with all four arguments, computes the Henry & Valliant (2015) calibration
#' design effect (Eq 3.5, zero-correlation approximation):
#' \deqn{deff_H \approx deff_K(w) \cdot \hat{\sigma}^2_u / \hat{\sigma}^2_y}
#' where \eqn{\hat{u}_i = y_i - x_i^T \hat{\beta}} are the GREG residuals from a survey-
#' weighted least-squares fit of \code{outcome} on the calibration covariate matrix
#' \eqn{X} (built from \code{data} + \code{target}), and \eqn{\hat{\sigma}^2_u, \hat{\sigma}^2_y}
#' are survey-weighted population unit variances per Eq 2.7. Falls back to \code{deff_K(w)}
#' when \code{target} is empty, when \eqn{\hat{\sigma}^2_y} is zero, or when calibration
#' margins are rank-deficient (with warning).
#'
#' All numerical computation is performed in C++17 via \code{.Call(C_rk_design_effect, ...)};
#' R + Python wrappers produce byte-identical output by construction.
#'
#' @param weights Numeric vector of calibrated weights. Must be finite, no NA, sum > 0.
#' @param outcome Numeric outcome vector (optional; 4-arg form only). Must match
#'   \code{length(weights)}, no NA, all finite.
#' @param data Data frame with the calibration covariate columns named in \code{target}.
#'   Must have \code{nrow(data) == length(weights)}. Each column named in \code{target}
#'   must contain only levels present in \code{names(target[[var]])} and no NAs.
#' @param target Named list of target proportions per categorical calibration margin
#'   (one named entry per margin variable in \code{data}). May be \code{list()}.
#' @return Numeric scalar design effect (\code{deff_K} for 1-arg, \code{deff_H} for 4-arg).
#' @references
#'   Kish, L. (1965). \emph{Survey Sampling}. Wiley.
#'
#'   Henry, K.A. and Valliant, R. (2015). A design effect measure for calibration weighting
#'   in single-stage samples. \emph{Survey Methodology}, 41(2), 315-331.
#'   Statistics Canada Catalogue No. 12-001-X. Equation 3.5.
#' @export
design_effect <- function(weights, outcome = NULL, data = NULL, target = NULL) {
  if (is.null(outcome)) {
    if (!is.numeric(weights))
      stop("design_effect: weights must be a numeric vector")
    if (anyNA(weights))
      stop("design_effect: weights must not contain NA values")
    if (!all(is.finite(weights)))
      stop("design_effect: weights must all be finite")
    if (sum(weights) <= 0)
      stop("design_effect: sum(weights) must be positive")
    res <- .Call(C_rk_design_effect, as.double(weights), NULL, NULL, NULL, 0L)
    return(res$deff_K)
  }
  if (is.null(data) || is.null(target))
    stop("design_effect: 4-argument form requires both 'data' and 'target'")
  n <- length(weights)
  if (length(outcome) != n) stop("design_effect: 'outcome' length must equal length(weights)")
  if (nrow(data) != n) stop("design_effect: nrow(data) must equal length(weights)")
  K <- length(target)
  if (K == 0L) {
    res <- .Call(C_rk_design_effect, as.double(weights), as.double(outcome), NULL, NULL, 0L)
    return(res$deff_H)
  }
  missing_cols <- setdiff(names(target), names(data))
  if (length(missing_cols) > 0L)
    stop(sprintf("design_effect: data is missing column(s) named in target: %s",
                 paste(missing_cols, collapse = ", ")))
  data_codes <- integer(n * K)
  cat_counts  <- integer(K)
  for (k in seq_len(K)) {
    var  <- names(target)[k]
    levs <- names(target[[k]])
    col  <- data[[var]]
    if (anyNA(col))
      stop(sprintf("design_effect: data[['%s']] contains NA", var))
    bad <- setdiff(unique(col), levs)
    if (length(bad) > 0L)
      stop(sprintf("design_effect: data[['%s']] has level(s) {%s} not in target",
                   var, paste(bad, collapse = ", ")))
    f <- factor(col, levels = levs)
    code_vec <- as.integer(f) - 1L   # 0-based for C
    data_codes[seq(k, n * K, by = K)] <- code_vec
    cat_counts[k] <- length(levs)
  }
  res <- .Call(C_rk_design_effect, as.double(weights), as.double(outcome),
               data_codes, cat_counts, as.integer(K))
  if (isTRUE(res$rank_def == 1L))
    warning("design_effect: calibration margins rank-deficient; deff_H = deff_K")
  res$deff_H
}

#' Effective sample size given calibrated weights
#'
#' Computes \code{length(weights) / design_effect(weights)} (Kish 1965 form).
#'
#' @param weights Numeric vector of calibrated weights. Must be finite, no NA, sum > 0.
#' @return Numeric scalar effective sample size.
#' @references
#'   Kish, L. (1965). \emph{Survey Sampling}. Wiley.
#' @export
effective_sample_size <- function(weights) {
  length(weights) / design_effect(weights)
}
