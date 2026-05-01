#!/usr/bin/env Rscript
# WI-1 empirical K_warm=8 validation.
# Sweep IEPPA+SRAA on stepstone for K_warm in {1,2,4,8,16,32}.
# Verify K_warm=8 reaches max_err <= 1e-2 (basin-of-attraction proxy for
# warm-start utility — K=8 SRAA iterations should drive lf close enough to
# the fixed point that Newton's quadratic-region radius is plausibly entered).
Sys.setenv(OMP_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(leafblower)
  library(arrow)
  library(jsonlite)
})
set.seed(1L)

data_path <- {
  if (file.exists("benchmarks/stepstone_bench_data.parquet"))
    "benchmarks/stepstone_bench_data.parquet"
  else "../../benchmarks/stepstone_bench_data.parquet"
}
tgt_path <- {
  if (file.exists("benchmarks/stepstone_bench_targets.json"))
    "benchmarks/stepstone_bench_targets.json"
  else "../../benchmarks/stepstone_bench_targets.json"
}
df  <- read_parquet(data_path); df$uuid <- NULL
tgt <- lapply(jsonlite::fromJSON(tgt_path), function(t) { t <- unlist(t); t / sum(t) })
for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])

K_warms <- c(1L, 2L, 4L, 8L, 16L, 32L)
results <- data.frame(K_warm = integer(0), max_err = numeric(0), iters = integer(0))
for (K in K_warms) {
  r <- suppressWarnings(harvest(df, tgt, method = "ieppa", accelerate = TRUE,
    max_weight = 5, min_weight = 0, max_iterations = as.integer(K),
    attach_weights = FALSE, verbose = 0L))
  res <- attr(r, "result")
  results <- rbind(results, data.frame(K_warm = K,
                                       max_err = res$max_error,
                                       iters   = res$iterations))
}
dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
write.csv(results, "benchmarks/results/k_warm_sweep.csv", row.names = FALSE)
cat("K_warm sweep:\n"); print(results)

k8_max_err <- results$max_err[results$K_warm == 8L]
cat(sprintf("\nK_warm=8 max_err = %.3e (threshold for PROCEED: <= 1e-2)\n", k8_max_err))
if (is.na(k8_max_err) || k8_max_err > 1e-2) {
  cat(sprintf("ABORT: K_warm=8 max_err=%.3e > 1e-2 -- IEPPA contraction insufficient on stepstone\n",
              k8_max_err))
  quit(status = 2)
}
cat("PROCEED: K_warm=8 reaches the basin-of-attraction proxy threshold\n")
