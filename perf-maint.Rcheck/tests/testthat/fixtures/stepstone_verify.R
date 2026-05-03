#!/usr/bin/env Rscript
# Verification-only: report n + per-margin counts from the saved parquet.
# No solver run. User confirms data integrity before reference benchmark.
#
# Expected per ORIGINAL Rmd (Salary-Data/code/22-weighting-create-weights-2025-firmsize.Rmd):
#   n after imputation: 1,248,521 (line 484 comment)
#   K = 9 margins
#   margin order (lines 457-466):
#     rk_i_loc_time_age10_gender, rk_i_loc_wz, rk_age10_gender,
#     rk_i_loc_time_gender, rk_i_loc_fsize, rk_i_loc_Abschluss_gender,
#     rk_gender_time, rk_time, rk_gender
#   ESS: 760,219.5 (line 494 comment) → DEFF ≈ 1.64
#
# Known discrepancy from prior session: the saved parquet reports n=1,582,732
# (NOT 1,248,521). Likely causes: different wave cutoff, different imputation,
# or regenerated data. Verification report will make the divergence visible.

suppressPackageStartupMessages({
  library(arrow); library(jsonlite)
})

DATA_PATH    <- "benchmarks/stepstone_fulldata_bench_data.parquet"
TARGETS_PATH <- "benchmarks/stepstone_fulldata_bench_targets.json"

stopifnot(file.exists(DATA_PATH), file.exists(TARGETS_PATH))

df  <- read_parquet(DATA_PATH)
tgt <- fromJSON(TARGETS_PATH)
tgt <- lapply(tgt, function(t) unlist(t))

cat("=== Parquet + targets inventory ===\n")
cat(sprintf("data file:    %s (%.1f MB)\n", DATA_PATH, file.info(DATA_PATH)$size/1024^2))
cat(sprintf("targets file: %s (%.1f kB)\n", TARGETS_PATH, file.info(TARGETS_PATH)$size/1024))
cat(sprintf("rows in parquet: %s\n", format(nrow(df), big.mark=",")))
cat(sprintf("cols in parquet: %d\n", ncol(df)))
cat(sprintf("column names: %s\n", paste(names(df), collapse=", ")))

cat(sprintf("\n=== Rmd expected ===\n"))
cat("  n = 1,248,521  (Rmd line 484 comment)\n")
cat("  K = 9 margins\n")
cat("  ESS = 760,219.5  (Rmd line 494 comment)  →  DEFF ≈ 1.64\n")

cat(sprintf("\n=== Delta vs Rmd ===\n"))
cat(sprintf("  n: parquet=%s - Rmd=1,248,521 = %+d (%+.2f%%)\n",
    format(nrow(df), big.mark=","),
    nrow(df) - 1248521L,
    100 * (nrow(df) - 1248521L) / 1248521L))

# Margin order
expected_order <- c(
  "rk_i_loc_time_age10_gender", "rk_i_loc_wz", "rk_age10_gender",
  "rk_i_loc_time_gender", "rk_i_loc_fsize", "rk_i_loc_Abschluss_gender",
  "rk_gender_time", "rk_time", "rk_gender"
)
cat(sprintf("\n=== Margin order (Rmd line 456-466) ===\n"))
if (identical(names(tgt), expected_order)) {
  cat("  OK: all 9 margins in expected order.\n")
} else {
  cat("  MISMATCH:\n")
  cat(sprintf("    parquet: %s\n", paste(names(tgt), collapse=", ")))
  cat(sprintf("    Rmd:     %s\n", paste(expected_order, collapse=", ")))
}

# Per-margin counts: raw (as in JSON) + filtered (drop cells absent from data + zero-obs)
cat(sprintf("\n=== Per-margin cell counts (raw targets JSON) ===\n"))
raw_counts <- sapply(tgt, length)
for (nm in names(tgt)) cat(sprintf("  %-32s %d cats\n", nm, length(tgt[[nm]])))
cat(sprintf("  TOTAL (raw): %d\n", sum(raw_counts)))

cat(sprintf("\n=== Per-margin filtered (drop absent + zero-obs, renormalized) ===\n"))
filt_counts <- integer(length(tgt))
names(filt_counts) <- names(tgt)
for (i in seq_along(tgt)) {
  nm <- names(tgt)[i]
  sv <- table(df[[nm]])
  present <- names(sv)[sv > 0]
  t_ <- tgt[[nm]]
  keep <- names(t_) %in% present
  dropped_missing <- sum(!keep)
  t_ <- t_[keep]
  n_per_cat <- as.integer(sv[names(t_)])
  zero_obs <- which(n_per_cat == 0 & t_ > 0)
  dropped_zero <- length(zero_obs)
  if (dropped_zero > 0) t_ <- t_[-zero_obs]
  filt_counts[i] <- length(t_)
  cat(sprintf("  %-32s raw=%3d filt=%3d  (dropped: absent=%d, zero_obs=%d)\n",
      nm, length(tgt[[nm]]), filt_counts[i], dropped_missing, dropped_zero))
}
cat(sprintf("  TOTAL (filtered): %d\n", sum(filt_counts)))

cat(sprintf("\n=== Prior-session value for comparison ===\n"))
cat("  /tmp/stepstone_3algo_bench.R output (earlier this session):\n")
cat("    stepstone full: n=1,582,732 K=9 total_cats=836\n")
cat("    cat counts: 408, 204, 12, 68, 68, 68, 4, 2, 2\n")
cat(sprintf("  This run filtered counts: %s (total=%d)\n",
    paste(filt_counts, collapse=", "), sum(filt_counts)))

cat(sprintf("\n=== Data-type sanity ===\n"))
for (nm in names(tgt)) {
  col <- df[[nm]]
  cat(sprintf("  %-32s class=%s unique_in_data=%d\n",
      nm, paste(class(col), collapse=","), length(unique(col))))
}

cat("\nDone. Review then authorize benchmark run.\n")
