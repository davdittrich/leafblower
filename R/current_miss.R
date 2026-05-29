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
    # Encode NA observations into the literal "NA" bin so an injected
    # add_na_proportion target is counted (as.character(NA) is NA_character_,
    # which never equals the string "NA"). Denominator spans all observations
    # so shares sum to 1.
    col_char <- as.character(col)
    col_char[is.na(col)] <- "NA"
    W   <- sum(weights)
    errs <- vapply(names(tgt), function(lv) {
      mask <- col_char == lv
      prop <- if (W > 0) sum(weights[mask]) / W else 0.0
      abs(prop - tgt[[lv]])
    }, numeric(1))
    if (length(errs) == 0L) return(0.0)
    max(errs)
  }, numeric(1))
}
