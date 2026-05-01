#!/usr/bin/env Rscript
# BLSE routing sweep: localize ieppa/raking crossover across compression ratio.
# Runs a 3x3 grid: mc_ratio (M_cell/n) x n.

suppressPackageStartupMessages({
  library(leafblower); library(bench)
})

# Grid: mc_ratio = M_cell/n (compression = n/M_cell = 1/mc_ratio)
# mc_ratio=0.5 → 2x compression (at threshold), 0.1 → 10x, 0.02 → 50x
mc_ratios <- c(0.5, 0.1, 0.02)  # n/M_cell = 2, 10, 50
ns        <- c(10000L, 50000L, 200000L)

results <- list()
for (n in ns) {
  for (mc in mc_ratios) {
    M_cell_approx <- max(1L, round(n * mc))
    # Build a problem with approx M_cell cells:
    # Use 2 margins with nj1*nj2 = M_cell_approx
    nj1 <- max(2L, round(sqrt(M_cell_approx)))
    nj2 <- max(2L, round(M_cell_approx / nj1))
    set.seed(42L)
    df <- data.frame(
      v1 = factor(sample(nj1, n, replace=TRUE)),
      v2 = factor(sample(nj2, n, replace=TRUE))
    )
    tgt <- list(
      v1 = setNames(rep(1/nj1, nj1), as.character(1:nj1)),
      v2 = setNames(rep(1/nj2, nj2), as.character(1:nj2))
    )
    compression <- n / (nj1 * nj2)
    label <- sprintf("n=%d mc=%.2f compression=%.1fx", n, mc, compression)
    cat(sprintf("  %-52s ...", label)); flush.console()
    t_ieppa <- system.time(suppressWarnings(
      harvest(df, tgt, method="ieppa",  max_iterations=200L, verbose=0L)
    ))["elapsed"]
    t_raking <- system.time(suppressWarnings(
      harvest(df, tgt, method="raking", max_iterations=200L, verbose=0L)
    ))["elapsed"]
    winner <- if (t_ieppa < t_raking) "ieppa" else "raking"
    cat(sprintf("  ieppa=%.2fs raking=%.2fs  winner=%s\n",
      t_ieppa, t_raking, winner))
    results[[length(results)+1]] <- list(
      n=n, mc_ratio=mc, compression=compression,
      t_ieppa=t_ieppa, t_raking=t_raking, winner=winner
    )
  }
}

cat("\n=== Routing crossover summary ===\n")
for (r in results) {
  cat(sprintf("  n=%6d  compr=%5.1fx  ieppa=%5.2fs  raking=%5.2fs  winner=%-6s\n",
    r$n, r$compression, r$t_ieppa, r$t_raking, r$winner))
}

saveRDS(results, "benchmarks/blse_routing_sweep.rds")
cat("\nSaved to benchmarks/blse_routing_sweep.rds\n")
