#' Weighted proportions
#'
#' @param x Factor or character vector.
#' @param weights Numeric weights, same length as x.
#' @return Named numeric vector of weighted proportions summing to 1.
#' @export
weighted_pct <- function(x, weights) {
  lvls    <- if (is.factor(x)) levels(x) else sort(unique(x[!is.na(x)]))
  total_w <- sum(weights[!is.na(x)])
  out <- vapply(lvls, function(lv) {
    mask <- !is.na(x) & (as.character(x) == lv)
    if (total_w > 0) sum(weights[mask]) / total_w else 0.0
  }, numeric(1))
  stats::setNames(out, lvls)
}
