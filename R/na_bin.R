#' Encode the injected NA-bin mask for a diagnostic column
#'
#' Internal helper shared by \code{diagnose_weights()} and
#' \code{get_current_miss()}. Given a data column and whether the
#' corresponding target has an injected \code{"NA"} bin
#' (\code{target_has_na = "NA" \%in\% names(target[[varname]])}), returns the
#' boolean mask of observations that fall into that NA bin: rows that are
#' genuinely missing (\code{is.na(col)}) \emph{or} whose value is the literal
#' string \code{"NA"} (conflation), matching how
#' \code{harvest(add_na_proportion = TRUE)} encodes the injected NA bin
#' (see harvest.R:130-131, 475).
#'
#' When \code{target_has_na} is \code{FALSE} the injected-bin concept does
#' not apply; the returned mask is simply \code{is.na(col)} (genuinely-
#' missing rows), for callers that need to exclude NA observations from
#' level masks and denominators.
#'
#' @param col A data column (any atomic vector).
#' @param target_has_na Logical: does the corresponding target have an
#'   \code{"NA"} bin?
#' @return Logical vector, same length as \code{col}.
#' @keywords internal
.encode_na_bin_mask <- function(col, target_has_na) {
  is_na <- is.na(col)
  if (!target_has_na) return(is_na)
  is_na | (!is_na & as.character(col) == "NA")
}
