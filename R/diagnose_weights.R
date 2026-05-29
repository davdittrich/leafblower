#' Diagnose calibration quality
#'
#' Returns a data frame comparing original and weighted marginals to targets.
#'
#' @param data Data frame used in calibration.
#' @param target Named list of target proportions (same format as harvest()).
#' @param weights Numeric vector of calibrated weights, length nrow(data).
#' @return Data frame with columns: variable, level, prop_original, prop_weighted,
#'   target, error_original, error_weighted.
#' @export
diagnose_weights <- function(data, target, weights) {
  if (!is.list(target))
    stop("target must be a named list of named numeric vectors")
  if (length(weights) != nrow(data))
    stop("weights length must equal nrow(data)")

  total_rows <- sum(vapply(target, length, integer(1L)))
  rows <- vector("list", total_rows)
  row_idx <- 0L

  for (varname in names(target)) {
    col <- data[[varname]]
    if (is.null(col)) stop("Variable '", varname, "' not found in data")
    tgt      <- target[[varname]]
    is_na    <- is.na(col)
    has_na_bin <- "NA" %in% names(tgt)
    if (has_na_bin) {
      # add_na_proportion case: fill NA observations into the literal "NA" bin
      # so the injected target is counted (as.character(NA) is NA_character_,
      # which never matches the string "NA"). All-obs denominators so the
      # prop_original / prop_weighted shares sum to 1 across bins (incl. NA).
      col_char <- as.character(col)
      col_char[is_na] <- "NA"
      n_total <- length(col)
      w_total <- sum(weights)
    } else {
      # No "NA" bin (common case): exclude NA observations from BOTH the level
      # masks and the denominators. harvest drops NA/gid<0 obs from the
      # marginal constraints (raking.cpp: `if (g>=0)`), so named-level shares
      # must be measured over non-NA obs only — otherwise they sum to <1 and
      # produce a spurious error_weighted ~ -na_frac*target on every level for
      # well-calibrated data.
      col_char <- as.character(col)
      col_char[is_na] <- NA_character_
      n_total <- sum(!is_na)
      w_total <- sum(weights[!is_na])
    }

    for (lvl in names(tgt)) {
      mask      <- !is.na(col_char) & col_char == lvl
      prop_orig <- if (n_total > 0L) sum(mask) / n_total else 0.0
      prop_wtd  <- if (w_total > 0.0) sum(weights[mask]) / w_total else 0.0
      tgt_val   <- tgt[[lvl]]
      row_idx   <- row_idx + 1L
      rows[[row_idx]] <- data.frame(
        variable       = varname,
        level          = lvl,
        prop_original  = prop_orig,
        prop_weighted  = prop_wtd,
        target         = tgt_val,
        error_original = prop_orig - tgt_val,
        error_weighted = prop_wtd - tgt_val,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
