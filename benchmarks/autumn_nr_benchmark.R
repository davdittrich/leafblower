#!/usr/bin/env Rscript
# benchmarks/autumn_nr_benchmark.R
#
# Benchmark autumn::harvest method="calibrate" (logit distance, enforces box)
# vs method="rake" (IPF, post-hoc clip) on full Stepstone data (n=1,582,732).
#
# method="rake":      standard IPF — already measured: 2144s, NO strict box
# method="calibrate": logit distance with max_weight<Inf → enforces box in-loop
#                     (replaces the removed method="nr")

suppressPackageStartupMessages({
  library(arrow)
  library(autumn)
  library(jsonlite)
})

ROOT         <- normalizePath(".")
DATA_PATH    <- file.path(ROOT, "benchmarks", "stepstone_fulldata_bench_data.parquet")
TARGETS_PATH <- file.path(ROOT, "benchmarks", "stepstone_fulldata_bench_targets.json")

cat("=== autumn method='nr' benchmark on full Stepstone data ===\n")

# ── load data ─────────────────────────────────────────────────────────────────
cat("Loading data...\n")
df      <- arrow::read_parquet(DATA_PATH)
targets_raw <- jsonlite::read_json(TARGETS_PATH, simplifyVector = FALSE)
targets <- lapply(targets_raw, function(v) {
  vals <- unlist(v)
  vals / sum(vals)  # renormalise to fix float round-trip
})

n <- nrow(df)
cat(sprintf("n = %s | margins = %d | total categories = %d\n\n",
    format(n, big.mark = ","),
    length(targets),
    sum(sapply(targets, length))))

MAX_ITER    <- 3000L
MAX_WEIGHT  <- 5
TOL         <- 1e-3
AUTUMN_CONV <- c(pct = 1e-15, absolute = TOL)

time_one <- function(expr) {
  t0  <- proc.time()["elapsed"]
  res <- force(expr)
  list(result = res, ms = (proc.time()["elapsed"] - t0) * 1000)
}

# ── warmup (minimal iterations) ───────────────────────────────────────────────
cat("Warmup (max_iterations=10)...\n")
invisible(autumn::harvest(df, targets, method = "calibrate",
                          max_weight = MAX_WEIGHT, max_iterations = 10L,
                          convergence = AUTUMN_CONV))
cat("Warmup done.\n\n")

# ── timed run ─────────────────────────────────────────────────────────────────
cat("--- autumn::harvest method='calibrate' (logit distance, enforces box) ---\n")
cat(sprintf("Settings: max_iter=%d  tol=%.0e  max_weight=%.1f\n\n",
            MAX_ITER, TOL, MAX_WEIGHT))

r_nr <- time_one(autumn::harvest(df, targets, method = "calibrate",
                                 max_weight = MAX_WEIGHT,
                                 max_iterations = MAX_ITER,
                                 convergence = AUTUMN_CONV))

w_nr <- r_nr$result$weights
cat(sprintf("  time:     %.0f ms  (%.1f s)\n", r_nr$ms, r_nr$ms / 1000))
cat(sprintf("  weights:  min=%.3f  med=%.3f  max=%.3f\n",
            min(w_nr), median(w_nr), max(w_nr)))
cat(sprintf("  w > 5:    %d observations (%.2f%%)\n",
            sum(w_nr > MAX_WEIGHT), 100 * mean(w_nr > MAX_WEIGHT)))
cat(sprintf("  DEFF:     %.3f\n", autumn::design_effect(r_nr$result)))
cat(sprintf("  ESS:      %.0f / %d\n\n",
            autumn::effective_sample_size(r_nr$result), n))

# ── comparison summary ─────────────────────────────────────────────────────────
rake_ms <- 2144373  # from previous benchmark run
lb_r_ms <-  133296
lb_py_ms <- 127000  # Python leafblower (median of 3 runs)

cat("=== Summary (n=1,582,732, 9 margins, 836 categories, tol=1e-3, max_weight=5) ===\n")
cat(sprintf("  autumn rake:       %8.0f ms (%6.1f s)  | strict box: NO\n",
            rake_ms, rake_ms / 1000))
cat(sprintf("  autumn calibrate:  %8.0f ms (%6.1f s)  | strict box: YES (logit)\n",
            r_nr$ms, r_nr$ms / 1000))
cat(sprintf("  R leafblower:      %8.0f ms (%6.1f s)  | strict box: YES\n",
            lb_r_ms, lb_r_ms / 1000))
cat(sprintf("  Python leafblower: %8.0f ms (%6.1f s)  | strict box: YES\n",
            lb_py_ms, lb_py_ms / 1000))
cat(sprintf("  calibrate speedup vs rake:          %.1fx\n", rake_ms / r_nr$ms))
cat(sprintf("  leafblower (R) speedup vs calibrate: %.1fx\n", r_nr$ms / lb_r_ms))
