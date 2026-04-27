#!/usr/bin/env Rscript
# Generates tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds
# Run: Rscript data-raw/gen_ieppa_kl_ref.R
# Records iEPPA's KL divergence at best_iter on stepstone-fulldata.
# Used by A1 test: sinkhorn KL at convergence must be <= this value.
suppressPackageStartupMessages({
  library(leafblower)
  library(arrow)
  library(jsonlite)
})

data <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
tgt_raw <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
target <- lapply(tgt_raw, function(x) { v <- unlist(x); v / sum(v) })

cat("Running iEPPA on stepstone-fulldata to record KL at best_iter...\n")
suppressWarnings(
  w <- leafblower::harvest(data, target, method = "ieppa",
                           max_weight = 5, max_iterations = 3000,
                           attach_weights = FALSE)
)
r <- attr(w, "result")

ref <- list(
  kl_at_best_iter = r$convergence_used$solver_objective,
  best_iter       = r$best_iter,
  max_error       = r$max_error,
  ieppa_version   = as.character(packageVersion("leafblower")),
  date            = Sys.Date()
)
cat(sprintf("  KL at best_iter=%d: %.4e\n", ref$best_iter, ref$kl_at_best_iter))
cat(sprintf("  max_error: %.4e\n", ref$max_error))

saveRDS(ref, "tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds")
cat("Saved to tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds\n")
