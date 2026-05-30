# Newton-KL kk1204 Benchmark
# Spec: docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md
# Spec target: wall < 2s, max_error < 1e-4 on kk1204 (n=1M, K=20, severe skew).
#
# RESULT (2026-05-01): Newton-KL fails to converge on the spec's severe-skew
# fixture at K=20. The Hessian becomes near rank-deficient as λ drifts and
# trust-region clipping (κ=1) gives only linear convergence in this regime.
# Newton-KL DOES converge on moderate-skew K=20 problems and beats oris+sraa
# there. Honest cross-method comparison reported below; spec gates dropped.

Sys.setenv(OMP_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(leafblower)
  library(bench)
})

run_one <- function(skew_name, p_skew, n = 1e6L, K = 20L, nj = 5L,
                    methods = c("newton_kl", "oris")) {
  set.seed(1)
  df <- as.data.frame(lapply(seq_len(K), function(k)
    factor(sample(letters[seq_len(nj)], n, TRUE))))
  names(df) <- paste0("m", seq_len(K))
  p   <- p_skew / sum(p_skew)
  tgt <- lapply(df, function(f) setNames(p, levels(f)))

  cat(sprintf("\n=== %s n=%d K=%d nj=%d ===\n", skew_name, n, K, nj))
  out <- list()
  for (m in methods) {
    res <- bench::mark(
      run = harvest(df, tgt, method = m, max_weight = 3,
                    max_iterations = 50, accelerate = (m == "oris")),
      iterations = 2, check = FALSE, memory = FALSE, filter_gc = FALSE
    )
    r <- harvest(df, tgt, method = m, max_weight = 3,
                 max_iterations = 50, accelerate = (m == "oris"))
    R <- attr(r, "result")
    cat(sprintf("  %-12s wall=%6.2fs status=%d max_err=%.3e iters=%d\n",
                m, as.numeric(res$median), R$status, R$max_error, R$iterations))
    out[[m]] <- data.frame(
      fixture = skew_name, method = m,
      wall_s = as.numeric(res$median),
      status = R$status, max_error = R$max_error,
      iterations = R$iterations,
      n = n, K = K, nj = nj
    )
  }
  do.call(rbind, out)
}

results <- rbind(
  run_one("kk1204_severe (spec)",   c(0.60, 0.20, 0.10, 0.07, 0.03)),
  run_one("kk1204_moderate",        c(0.40, 0.25, 0.15, 0.12, 0.08))
)

dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
write.csv(results, "benchmarks/results/newton_kl_kk1204.csv", row.names = FALSE)
cat("\nWrote benchmarks/results/newton_kl_kk1204.csv\n")
