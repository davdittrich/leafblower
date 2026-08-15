#' Diagnose calibration quality
#'
#' Returns a data frame comparing original and weighted marginals to targets.
#'
#' @param data Data frame used in calibration.
#' @param target Named list of target proportions (same format as harvest()).
#' @param weights Numeric vector of calibrated weights, length nrow(data).
#' @return Data frame with columns: variable, level, prop_original, prop_weighted,
#'   target, error_original, error_weighted.
#' @details An \code{"NA"} entry in \code{target} is treated as the missing-data
#'   bin: a row counts toward it if it is \code{NA} \emph{or} its value is the
#'   literal string \code{"NA"} (conflation), matching how
#'   \code{harvest(add_na_proportion = TRUE)} encodes the injected NA bin. For a
#'   hand-built target that names a real category \code{"NA"} without
#'   \code{add_na_proportion}, genuinely-missing rows are also counted there (a
#'   documented collision; see \code{?harvest}).
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
      # add_na_proportion case: the injected "NA" bin CONFLATES true-missing rows
      # with rows whose value is the literal string "NA" (is_na OR char=="NA"),
      # matching the solver's documented encoding (harvest.R:130-131: "real NA
      # values AND a literal ... 'NA' will collide -- both mapped to the injected
      # NA bin"; solver fill harvest.R:475). All-obs denominators, so prop_original
      # / prop_weighted sum to 1 across bins (incl. NA) when every obs falls in a
      # named level or the NA bin (out-of-vocabulary values land in no bin -> sum < 1).
      # Supersedes 4ihf.4 (mask-only), which under-reported the NA bin and broke
      # sum(shares)==1 even with no OOV rows [4ihf.5].
      col_char <- as.character(col)
      col_char[is_na] <- NA_character_
      n_total <- length(col)
      w_total <- sum(weights)
    } else {
      # No "NA" bin (common case): exclude NA observations from BOTH the level
      # masks and the denominators. harvest drops NA/gid<0 obs from the
      # marginal constraints (raking.cpp: `if (g>=0)`), so named-level shares
      # must be measured over non-NA obs only -- otherwise they sum to <1 and
      # produce a spurious error_weighted ~ -na_frac*target on every level for
      # well-calibrated data.
      col_char <- as.character(col)
      col_char[is_na] <- NA_character_
      n_total <- sum(!is_na)
      w_total <- sum(weights[!is_na])
    }

    # xc1s.13(g): precompute all per-level counts + weighted sums in one pass (was an
    # O(n*L) mask per level). tabulate ignores NA & out-of-vocab => integer counts equal
    # the old sum(mask). Weighted sums use tapply+base sum() (long-double accumulation,
    # bit-identical to the old sum(weights[mask]); NOT rowsum, whose C kernel sums in
    # plain double). The injected "NA" bin (is_na OR literal char=="NA") is then
    # overridden with its conflation mask.
    lvs    <- names(tgt)
    counts <- tabulate(match(col_char, lvs), nbins = length(lvs))
    ws     <- tapply(weights, factor(col_char, levels = lvs), sum)
    ws[is.na(ws)] <- 0
    wsums  <- as.numeric(ws)
    if (has_na_bin) {
      na_mask <- .encode_na_bin_mask(col, TRUE)
      na_i    <- match("NA", lvs)
      counts[na_i] <- sum(na_mask)
      wsums[na_i]  <- sum(weights[na_mask])
    }
    names(counts) <- lvs
    names(wsums)  <- lvs

    for (lvl in names(tgt)) {
      prop_orig <- if (n_total > 0L) counts[[lvl]] / n_total else 0.0
      prop_wtd  <- if (w_total > 0.0) wsums[[lvl]] / w_total else 0.0
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
