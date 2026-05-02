#!/usr/bin/env Rscript
# WL-5: kk1204 K=20 severe-skew PARTIAL gate (gap < 1e-3 desirable; not a halt)
Sys.setenv(OMP_NUM_THREADS = "1")
suppressPackageStartupMessages(library(leafblower))
set.seed(1)

n  <- 1e6L; K <- 20L; nj <- 5L
df <- as.data.frame(lapply(seq_len(K), function(k)
  factor(sample(letters[seq_len(nj)], n, TRUE))))
names(df) <- paste0("m", seq_len(K))
tgt <- lapply(df, function(f) {
  p <- c(0.6, 0.2, 0.1, 0.07, 0.03)
  p <- p / sum(p)
  setNames(p, levels(f))
})

cat(sprintf("[wl-5] kk1204 K=20 severe-skew: n=%d K=%d nj=%d\n", n, K, nj))

t0 <- Sys.time()
r <- suppressWarnings(harvest(df, tgt, method="newton_kl",
  max_weight=3, max_iterations=50L, attach_weights=FALSE, verbose=0L))
wall_ms <- as.numeric(Sys.time()-t0, units="secs") * 1000
res <- attr(r, "result")

result_row <- data.frame(
  gap = res$max_error,
  status = res$status,
  n_projected_dims = ifelse(is.null(res$n_projected_dims), NA_integer_, res$n_projected_dims),
  n_iter = res$iterations,
  wall_time_ms = round(wall_ms, 2)
)

dir.create("benchmarks/results", showWarnings=FALSE, recursive=TRUE)
write.csv(result_row, "benchmarks/results/tsvd_kk1204_K20.csv", row.names=FALSE)

cat(sprintf("[wl-5] gap=%.3e status=%d n_proj=%d n_iter=%d wall=%.0fms\n",
  result_row$gap, result_row$status, result_row$n_projected_dims,
  result_row$n_iter, result_row$wall_time_ms))

if (result_row$gap < 1e-3) {
  cat(sprintf("[wl-5] PARTIAL GATE PASS — gap=%.3e < 1e-3\n", result_row$gap))
} else {
  cat(sprintf("[wl-5] PARTIAL GATE not met (gap=%.3e ≥ 1e-3) — Out-of-Scope, Epic-E candidate. NOT halting.\n",
    result_row$gap))
}
