#' Current calibration miss (matching autumn's export name)
#'
#' @param data Data frame.
#' @param target Named list of target proportions.
#' @param weights Numeric weight vector.
#' @return Named numeric vector of max absolute errors per variable.
#' @details An \code{"NA"} entry in \code{target} is treated as the missing-data
#'   bin: a row counts toward it if it is \code{NA} \emph{or} its value is the
#'   literal string \code{"NA"} (conflation), matching how
#'   \code{harvest(add_na_proportion = TRUE)} encodes the injected NA bin. For a
#'   hand-built target that names a real category \code{"NA"} without
#'   \code{add_na_proportion}, genuinely-missing rows are also counted there (a
#'   documented collision; see \code{?harvest}).
#' @export
get_current_miss <- function(data, target, weights) {
  vapply(names(target), function(v) {
    col <- data[[v]]
    # CR-F5 (dtkn.5): fail loudly on an absent variable instead of silently
    # returning max(tgt) as a fabricated "miss" (NULL col -> logical(0) -> W=0).
    # Matches diagnose_weights().
    if (is.null(col)) stop("Variable '", v, "' not found in data")
    tgt <- target[[v]]
    is_na <- is.na(col)
    has_na_bin <- "NA" %in% names(tgt)
    if (has_na_bin) {
      # add_na_proportion case: the injected "NA" bin CONFLATES true-missing rows
      # with rows whose value is the literal string "NA" (is_na OR char=="NA"),
      # matching the solver's documented encoding (harvest.R:130-131: "real NA
      # values AND a literal ... 'NA' will collide -- both mapped to the injected
      # NA bin"; solver fill harvest.R:475 `key[is.na(key)] <- "NA"`). Denominator
      # spans ALL observations, so shares (incl. NA) sum to 1 when every obs falls
      # in a named level or the NA bin (an out-of-vocabulary value lands in no bin
      # -> sum < 1). Supersedes 4ihf.4 (mask-only), which under-reported the NA bin
      # and broke sum(shares)==1 even with no OOV rows [4ihf.5].
      col_char <- as.character(col)
      col_char[is_na] <- NA_character_
      W <- sum(weights)
    } else {
      # No "NA" bin (common case): exclude NA observations from BOTH the level
      # masks and the denominator. harvest drops NA/gid<0 obs from the marginal
      # constraints (raking.cpp: `if (g>=0)`), so the named-level shares must be
      # measured over non-NA obs only -- otherwise they sum to <1 and produce a
      # spurious miss ~ na_frac*target on every level for calibrated data.
      col_char <- as.character(col)
      col_char[is_na] <- NA_character_
      W <- sum(weights[!is_na])
    }
    # xc1s.13(g): one split+sum pass for per-level weighted sums (was an O(n*L) mask
    # per level). tapply groups by level (NA/out-of-vocab excluded) and totals each with
    # base sum(), so the sums are bit-identical to the old sum(weights[!is.na(col_char)
    # & col_char==lv]) -- long-double accumulation (NOT rowsum, which sums in plain
    # double). The injected "NA" bin (is_na OR literal char=="NA") is overridden with
    # its conflation mask (harvest.R:130-131,475). A literal-"NA" row lands in the NA
    # bin, not its own level.
    lvs   <- names(tgt)
    ws    <- tapply(weights, factor(col_char, levels = lvs), sum)
    ws[is.na(ws)] <- 0
    wsums <- as.numeric(ws)
    if (has_na_bin) {
      na_i <- match("NA", lvs)
      wsums[na_i] <- sum(weights[.encode_na_bin_mask(col, TRUE)])
    }
    errs <- vapply(seq_along(lvs), function(i) {
      prop <- if (W > 0) wsums[i] / W else 0.0
      abs(prop - tgt[[lvs[i]]])
    }, numeric(1))
    if (length(errs) == 0L) return(0.0)
    max(errs)
  }, numeric(1))
}
