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
    is_na <- is.na(col)
    has_na_bin <- "NA" %in% names(tgt)
    if (has_na_bin) {
      # add_na_proportion case: the injected "NA" bin CONFLATES true-missing rows
      # with rows whose value is the literal string "NA" (is_na OR char=="NA"),
      # matching the solver's documented encoding (harvest.R:130-131: "real NA
      # values AND a literal ... 'NA' will collide — both mapped to the injected
      # NA bin"; solver fill harvest.R:475 `key[is.na(key)] <- "NA"`). Denominator
      # spans ALL observations, so shares (incl. NA) sum to 1 when every obs falls
      # in a named level or the NA bin (an out-of-vocabulary value lands in no bin
      # -> Σ<1). Supersedes 4ihf.4 (mask-only), which under-reported the NA bin
      # and broke Σshares==1 even with no OOV rows [4ihf.5].
      col_char <- as.character(col)
      col_char[is_na] <- NA_character_
      W <- sum(weights)
    } else {
      # No "NA" bin (common case): exclude NA observations from BOTH the level
      # masks and the denominator. harvest drops NA/gid<0 obs from the marginal
      # constraints (raking.cpp: `if (g>=0)`), so the named-level shares must be
      # measured over non-NA obs only — otherwise they sum to <1 and produce a
      # spurious miss ~ na_frac*target on every level for calibrated data.
      col_char <- as.character(col)
      col_char[is_na] <- NA_character_
      W <- sum(weights[!is_na])
    }
    errs <- vapply(names(tgt), function(lv) {
      # The injected NA bin (lv == "NA" with has_na_bin) CONFLATES true-missings
      # with literal-"NA" rows (is_na OR char=="NA"), matching the solver
      # encoding (harvest.R:130-131,475). Real levels match the NA-cleared char
      # vector; a literal-"NA" row falls into the NA bin, not its own level.
      mask <- if (has_na_bin && lv == "NA")
                is_na | (!is.na(col) & as.character(col) == "NA")
              else !is.na(col_char) & col_char == lv
      prop <- if (W > 0) sum(weights[mask]) / W else 0.0
      abs(prop - tgt[[lv]])
    }, numeric(1))
    if (length(errs) == 0L) return(0.0)
    max(errs)
  }, numeric(1))
}
