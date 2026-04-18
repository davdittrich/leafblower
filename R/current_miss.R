#' Current calibration miss (matching autumn's export name)
#'
#' @param data Data frame.
#' @param target Named list of target proportions.
#' @param weights Numeric weight vector.
#' @return Named numeric vector of max absolute errors per variable.
#' @export
get_current_miss <- function(data, target, weights) {
  vapply(names(target), function(v) {
    col <- data[[v]]
    tgt <- target[[v]]
    W   <- sum(weights[!is.na(col)])
    errs <- vapply(names(tgt), function(lv) {
      mask <- !is.na(col) & (as.character(col) == lv)
      prop <- if (W > 0) sum(weights[mask]) / W else 0.0
      abs(prop - tgt[[lv]])
    }, numeric(1))
    if (length(errs) == 0L) return(0.0)
    max(errs)
  }, numeric(1))
}
