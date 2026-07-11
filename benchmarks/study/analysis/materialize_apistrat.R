# STUDY-BRANCH-ONLY-DO-NOT-MERGE
#
# Materializes the R-only survey::apistrat dataset to a parquet file so the
# Python metrics aggregator (which cannot load R packages) can consume it.
# Data export ONLY -- no metric computation happens here (avoids R/Py drift).
#
# Spec: benchmarks/study/spec/canonical_survey_apistrat.json
#   data_ref:       pkg:survey::apistrat
#   design_weights: column:pw
#   margins:        ["stype"]
#   targets.stype:  {E: 4421, H: 755, M: 1018}

suppressPackageStartupMessages({
  library(survey)
  library(arrow)
  library(jsonlite)
})

spec_path <- "benchmarks/study/spec/canonical_survey_apistrat.json"
out_dir <- "benchmarks/study/results/_problem_data"
out_path <- file.path(out_dir, "canonical_survey_apistrat.parquet")

spec <- fromJSON(spec_path)
stopifnot(identical(spec$data_ref, "pkg:survey::apistrat"))

margin_cols <- spec$margins
design_weight_col <- sub("^column:", "", spec$design_weights)

data(api, package = "survey")
df <- apistrat

select_cols <- unique(c(margin_cols, design_weight_col))
missing_cols <- setdiff(select_cols, names(df))
if (length(missing_cols) > 0) {
  stop("apistrat missing expected columns: ", paste(missing_cols, collapse = ", "))
}

out <- df[, select_cols, drop = FALSE]

# problem_io convention: margin columns are cast to character.
for (col in margin_cols) {
  out[[col]] <- as.character(out[[col]])
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
arrow::write_parquet(out, out_path)

# --- Verification against spec ---
if (!is.null(spec$n)) {
  stopifnot(nrow(df) == spec$n)
}

for (col in margin_cols) {
  target_levels <- names(spec$targets[[col]])
  observed_levels <- sort(unique(out[[col]]))
  stopifnot(identical(sort(target_levels), observed_levels))
}

cat(sprintf(
  "materialize_apistrat.R: wrote %d rows x %d cols to %s\n",
  nrow(out), ncol(out), out_path
))
cat(sprintf("  margin columns: %s\n", paste(margin_cols, collapse = ", ")))
cat(sprintf("  design weight column: %s\n", design_weight_col))
for (col in margin_cols) {
  cat(sprintf("  %s levels: %s\n", col, paste(sort(unique(out[[col]])), collapse = ", ")))
}
