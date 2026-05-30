#!/usr/bin/env Rscript
# WI-0b: basin-overlap kill-switch for Epic-C ORIS warm-start.
#
# Purpose: falsify the warm-start hypothesis BEFORE writing any code.
# If running ORIS and Newton-KL on the same fixture (current master,
# no warm-start) lands at per-observation weights within 1e-3 (log-ratio),
# warm-start would be a functional no-op and the entire epic is wasted work.
#
# Saves benchmarks/results/basin_overlap_killswitch.csv and exits with
# status 0 on PROCEED, status 2 on ABORT.

Sys.setenv(OMP_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(leafblower)
  library(arrow)
  library(jsonlite)
})
set.seed(1L)

# Path-resolve fixtures (testthat-style: works from repo root or worktree).
data_path <- if (file.exists("benchmarks/stepstone_bench_data.parquet")) {
  "benchmarks/stepstone_bench_data.parquet"
} else {
  "../../benchmarks/stepstone_bench_data.parquet"
}
tgt_path <- if (file.exists("benchmarks/stepstone_bench_targets.json")) {
  "benchmarks/stepstone_bench_targets.json"
} else {
  "../../benchmarks/stepstone_bench_targets.json"
}

df <- read_parquet(data_path)
df$uuid <- NULL
tgt <- lapply(jsonlite::fromJSON(tgt_path), function(t) {
  t <- unlist(t)
  t / sum(t)
})
for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])

# ORIS + SRAA acceleration (current Epic-C baseline pipeline).
w_ie <- suppressWarnings(harvest(
  df, tgt,
  method = "oris", accelerate = TRUE,
  max_weight = 5, min_weight = 0,
  max_iterations = 200L,
  attach_weights = FALSE, verbose = 0L
))

# Newton-KL on current master (no warm-start exists yet).
w_nk <- suppressWarnings(harvest(
  df, tgt,
  method = "newton_kl",
  max_weight = 5, min_weight = 0,
  max_iterations = 200L,
  attach_weights = FALSE, verbose = 0L
))

# Element-wise log-ratio of weights. min_weight=0 means a solver may
# place an observation at the boundary; log(0/x)=-Inf, log(x/0)=+Inf.
# We report finite-only summary statistics plus a count of non-finite
# entries so the CSV stays numeric.
log_ratio <- log(w_ie / w_nk)
abs_lr <- abs(log_ratio)
finite_mask <- is.finite(abs_lr)
abs_lr_finite <- abs_lr[finite_mask]
stats <- list(
  max_abs_log_ratio    = if (length(abs_lr_finite)) max(abs_lr_finite) else NA_real_,
  median_abs_log_ratio = if (length(abs_lr_finite)) median(abs_lr_finite) else NA_real_,
  p95_abs_log_ratio    = if (length(abs_lr_finite)) quantile(abs_lr_finite, 0.95, names = FALSE) else NA_real_,
  n_close              = sum(finite_mask & abs_lr < 1e-3),
  n_nonfinite          = sum(!finite_mask),
  n_total              = length(log_ratio)
)

dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
write.csv(
  data.frame(metric = names(stats), value = unlist(stats)),
  "benchmarks/results/basin_overlap_killswitch.csv",
  row.names = FALSE
)

cat(sprintf(
  "[basin_overlap] max_abs_log_ratio=%.4e median=%.4e p95=%.4e\n",
  stats$max_abs_log_ratio, stats$median_abs_log_ratio, stats$p95_abs_log_ratio
))
cat(sprintf(
  "[basin_overlap] n_close (|log_ratio|<1e-3) / n_total = %d / %d ; n_nonfinite=%d\n",
  stats$n_close, stats$n_total, stats$n_nonfinite
))

# Decision: ABORT only if basins are tight everywhere (no boundary
# disagreement AND finite max-log-ratio below 1e-3). Any non-finite
# entry means the two solvers disagree at a boundary — basins differ.
if (stats$n_nonfinite == 0L &&
    !is.na(stats$max_abs_log_ratio) &&
    stats$max_abs_log_ratio < 1e-3) {
  cat(sprintf(
    "ABORT: max_log_ratio=%.4e < 1e-3 — warm-start is functional no-op, epic wasted work\n",
    stats$max_abs_log_ratio
  ))
  quit(status = 2)
} else {
  cat(sprintf(
    "PROCEED: max_log_ratio=%.4e (finite) ; n_nonfinite=%d\n",
    stats$max_abs_log_ratio, stats$n_nonfinite
  ))
  quit(status = 0)
}
