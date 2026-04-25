#!/usr/bin/env Rscript
# Generates tests/testthat/fixtures/raking_obs_reference_stepstone.rds
# Run BEFORE Plan B migration to capture obs-level raking reference.
suppressPackageStartupMessages({
  library(leafblower); library(arrow); library(jsonlite)
})
data <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
tgt_raw <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
target <- lapply(tgt_raw, function(x) { v <- unlist(x); v / sum(v) })

cat("Running obs-level raking to record reference...\n")
t0 <- proc.time()["elapsed"]
suppressWarnings(
  w <- leafblower::harvest(data, target, method="raking",
                           max_weight=5, max_iterations=500, attach_weights=FALSE)
)
elapsed_obs <- proc.time()["elapsed"] - t0

r <- attr(w, "result")
ref <- list(
  max_error     = r$max_error,
  best_error    = r$best_error,
  iterations    = r$iterations,
  status        = r$status,
  elapsed_obs   = elapsed_obs,
  version       = as.character(packageVersion("leafblower")),
  date          = Sys.Date()
)
cat(sprintf("  max_error=%.4e  elapsed=%.1fs  iterations=%d\n",
            ref$max_error, ref$elapsed_obs, ref$iterations))
saveRDS(ref, "tests/testthat/fixtures/raking_obs_reference_stepstone.rds")
cat("Saved.\n")
