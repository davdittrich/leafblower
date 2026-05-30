#!/usr/bin/env Rscript
# benchmarks/stepstone_fulldata_homotopy.R
#
# WU-6: stepstone-fulldata merge gate.
# Compares 5 leafblower configs + autumn reference on 1.58M-row data.
# Gate: leafblower_AB errRp (internal convergence measure) <= 1.60e-3.
# Note: tapply max_err and errRp are DIFFERENT metrics. The plan gate uses errRp.

suppressPackageStartupMessages({
  library(arrow)
  library(leafblower)
  library(jsonlite)
})

BENCH_DIR <- "benchmarks"
DATA_F    <- file.path(BENCH_DIR, "stepstone_fulldata_bench_data.parquet")
TGT_F     <- file.path(BENCH_DIR, "stepstone_fulldata_bench_targets.json")
RPT_F     <- file.path(BENCH_DIR, "stepstone_fulldata_homotopy_report.rds")

cat("=== WU-6: stepstone-fulldata homotopy merge gate ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Load data
cat("Loading data...\n")
data <- arrow::read_parquet(DATA_F)
cat(sprintf("  n = %s rows\n", format(nrow(data), big.mark = ",")))

# Load and renormalize targets (jsonlite returns nested lists)
target_raw <- jsonlite::fromJSON(TGT_F)
target <- lapply(target_raw, function(x) {
  v <- unlist(x)
  v / sum(v)
})
cat(sprintf("  margins: %d  |  total categories: %d\n\n",
    length(target), sum(sapply(target, length))))

# tapply max_err (external calibration check)
max_err_of <- function(w, data, target) {
  errs <- numeric(length(target))
  for (k in seq_along(target)) {
    tab <- tapply(w, data[[names(target)[k]]], sum) / sum(w)
    # Name-aligned subtraction: tapply alphabetizes; JSON target preserves
    # original order. Positional subtraction inflates error for reordered
    # categories (e.g. German compound names in rk_gender_time).
    errs[k] <- max(abs(tab[names(target[[k]])] - target[[k]]), na.rm = TRUE)
  }
  max(errs)
}

# Autumn reference
run_autumn <- function() {
  cat("--- autumn::harvest (reference) ---\n")
  t0 <- proc.time()["elapsed"]
  suppressWarnings(
    res <- autumn::harvest(
      data, target = target,
      max_weight   = 5,
      convergence  = c(pct = 1e-15, absolute = 1e-10),
      max_iterations = 3000,
      accelerate   = TRUE,
      attach_weights = FALSE
    )
  )
  elapsed <- proc.time()["elapsed"] - t0
  w <- as.numeric(res)
  me <- max_err_of(w, data, target)
  cat(sprintf("  max_err (tapply): %.3e  |  time: %.1fs\n\n", me, elapsed))
  list(config = "autumn", max_err = me, errRp = NA_real_, elapsed = elapsed, weights = w)
}

# leafblower runner — captures both tapply max_err and internal errRp
run_lb <- function(label, homotopy_levels, homotopy_start_factor,
                   homotopy_end_factor, scheduler,
                   homotopy_budget_p = 0.5,
                   eta_schedule = "fixed", eta_start = 1.0, eta_end = 1.0,
                   eta_schedule_power = 0.5,
                   max_iterations = 3000L) {
  cat(sprintf("--- leafblower::%s ---\n", label))
  cat(sprintf("    levels=%d  start_factor=%g  end_factor=%g  scheduler=%s  budget_p=%g\n",
      homotopy_levels, homotopy_start_factor, homotopy_end_factor,
      scheduler, homotopy_budget_p))
  t0 <- proc.time()["elapsed"]
  suppressWarnings(
    w_raw <- leafblower::harvest(
      data, target,
      max_weight          = 5,
      method              = "oris",
      max_iterations      = max_iterations,
      convergence         = list(absolute = 1e-10),
      homotopy_levels     = homotopy_levels,
      homotopy_start_factor = homotopy_start_factor,
      homotopy_end_factor = homotopy_end_factor,
      homotopy_budget_p   = homotopy_budget_p,
      scheduler           = scheduler,
      eta_schedule        = eta_schedule,
      eta_start           = eta_start,
      eta_end             = eta_end,
      eta_schedule_power  = eta_schedule_power,
      attach_weights      = FALSE
    )
  )
  elapsed <- proc.time()["elapsed"] - t0
  # Extract internal errRp from result attribute
  res_attr  <- attr(w_raw, "result")
  errRp     <- if (!is.null(res_attr)) res_attr$max_error else NA_real_
  iterations <- if (!is.null(res_attr)) res_attr$iterations else NA_integer_
  w <- as.numeric(w_raw)
  me <- max_err_of(w, data, target)
  cat(sprintf("  max_err (tapply): %.3e  |  errRp (internal): %.3e  |  iters: %d  |  time: %.1fs\n\n",
      me, errRp, iterations, elapsed))
  list(config = label, max_err = me, errRp = errRp, elapsed = elapsed,
       iterations = iterations, weights = w)
}

results <- list()

# Autumn reference
results[["autumn"]] <- run_autumn()

# Config A: baseline (no homotopy, round-robin)
results[["leafblower_A"]] <- run_lb(
  "leafblower_A (baseline)",
  homotopy_levels = 1L, homotopy_start_factor = 1.0, homotopy_end_factor = 1.0,
  scheduler = "round_robin"
)

# Config B: greedy scheduler only
results[["leafblower_B"]] <- run_lb(
  "leafblower_B (greedy only)",
  homotopy_levels = 1L, homotopy_start_factor = 1.0, homotopy_end_factor = 1.0,
  scheduler = "greedy"
)

# Config AB: homotopy_levels=5, homotopy_start_factor=10, greedy (primary gate config)
results[["leafblower_AB"]] <- run_lb(
  "leafblower_AB (homotopy+greedy, gate config)",
  homotopy_levels = 5L, homotopy_start_factor = 10.0, homotopy_end_factor = 1.0,
  scheduler = "greedy", homotopy_budget_p = 0.5
)

# Config AB2: wider homotopy, more levels
results[["leafblower_AB2"]] <- run_lb(
  "leafblower_AB2 (levels=10, start=20)",
  homotopy_levels = 10L, homotopy_start_factor = 20.0, homotopy_end_factor = 1.0,
  scheduler = "greedy", homotopy_budget_p = 0.5
)

# Config ABE: AB + Tang dynamic-eta
results[["leafblower_ABE"]] <- run_lb(
  "leafblower_ABE (homotopy+greedy+tang-eta)",
  homotopy_levels = 5L, homotopy_start_factor = 10.0, homotopy_end_factor = 1.0,
  scheduler = "greedy", homotopy_budget_p = 0.5,
  eta_schedule = "tang_dynamic", eta_start = 20.0, eta_end = 1.0,
  eta_schedule_power = 0.5
)

# Compute Pearson r vs autumn reference
w_autumn <- results[["autumn"]]$weights
pearson <- c(
  AB  = cor(results[["leafblower_AB"]]$weights,  w_autumn),
  AB2 = cor(results[["leafblower_AB2"]]$weights, w_autumn),
  ABE = cor(results[["leafblower_ABE"]]$weights, w_autumn)
)

# Build report data frame
rpt <- data.frame(
  config     = names(results),
  max_err    = sapply(results, `[[`, "max_err"),
  errRp      = sapply(results, `[[`, "errRp"),
  elapsed_s  = sapply(results, `[[`, "elapsed"),
  stringsAsFactors = FALSE
)
attr(rpt, "pearson_r") <- pearson

cat("=== Summary ===\n")
print(rpt[, c("config", "max_err", "errRp", "elapsed_s")], row.names = FALSE)
cat("\nPearson r vs autumn:\n")
print(pearson)

# Gate check — errRp is the internal convergence metric (comparable to lb_max_error in reference)
ab_errRp <- results[["leafblower_AB"]]$errRp
ab_me    <- results[["leafblower_AB"]]$max_err
cat(sprintf("\nGATE (errRp): leafblower_AB errRp = %.3e  (threshold <= 1.60e-3)\n", ab_errRp))
cat(sprintf("INFO (tapply): leafblower_AB max_err = %.3e  (reference baseline: %.3e)\n",
    ab_me, results[["leafblower_A"]]$max_err))
ref_errRp <- 0.00222  # lb_max_error from tests/testthat/fixtures/stepstone_reference_summary.rds
cat(sprintf("Reference errRp at git 9a97cc8: %.3e\n", ref_errRp))

if (!is.na(ab_errRp) && ab_errRp <= 1.60e-3) {
  cat("GATE: PASS\n")
} else {
  cat("GATE: FAIL — best errRp by config:\n")
  for (nm in names(results)) {
    eRp <- results[[nm]]$errRp
    if (!is.na(eRp)) cat(sprintf("  %s: errRp=%.3e\n", nm, eRp))
  }
}

saveRDS(rpt, RPT_F)
cat(sprintf("\nReport saved to %s\n", RPT_F))
