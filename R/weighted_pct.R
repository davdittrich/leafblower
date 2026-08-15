#' Weighted proportions
#'
#' @param x Factor or character vector.
#' @param weights Numeric weights, same length as x.
#' @return Named numeric vector of weighted proportions summing to 1.
#' @export
weighted_pct <- function(x, weights) {
  lvls    <- if (is.factor(x)) levels(x) else sort(unique(x[!is.na(x)]))
  total_w <- sum(weights[!is.na(x)])
  out <- numeric(length(lvls))
  if (total_w > 0) {
    # xc1s.13(f): one split+sum pass instead of L masked scans. tapply groups obs by
    # level (NA obs -> NA index, excluded) and totals each group with base sum(), so the
    # per-level sums are bit-identical to the old sum(weights[!is.na(x) &
    # as.character(x) == lv]) -- same subset, same long-double accumulation (NOT rowsum,
    # whose C kernel accumulates in plain double). Empty levels return NA -> 0.
    s <- tapply(weights, factor(as.character(x), levels = lvls), sum)
    s[is.na(s)] <- 0
    out <- as.numeric(s) / total_w
  }
  stats::setNames(out, lvls)
}
