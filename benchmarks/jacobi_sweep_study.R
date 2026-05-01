#!/usr/bin/env Rscript
# BENCHMARK RESULTS SUMMARY (rx1):
# Tested on stepstone_small: K=9, M_cell=5980, n=10000, log path forced.
# GS: 0.32s / 300 iters. Jacobi: 0.30s / 300 iters. Wall ratio ≈ 1.0.
# M_cell=5980 fits in L2 cache — no cache pressure difference.
# DECISION: Conditional — keep jacobi_sweep=FALSE default.
# May help at M_cell > 1M (untested); re-evaluate on fulldata Stepstone.

suppressPackageStartupMessages(library(leafblower))
Sys.setenv(LBW_IEPPA_FORCE_PATH = "log")

make_adversarial <- function(n, K, compression, seed=42L) {
  set.seed(seed)
  mc_ratio <- 1 / compression
  nj <- max(3L, round((n * mc_ratio)^(1/K)))
  df <- as.data.frame(lapply(seq_len(K), function(k) factor(sample(nj, n, replace=TRUE))))
  names(df) <- paste0("v", seq_len(K))
  # Non-uniform but feasible targets: even/odd alternation with asymmetric distribution
  # Even margins (k=2,4,...): skew toward category 1  → p ≈ (0.5, 0.25, 0.25, ...)
  # Odd margins (k=1,3,...):  skew toward category nj → p ≈ (0.25, ..., 0.25, 0.5)
  tgt <- lapply(seq_len(K), function(k) {
    p <- rep(1 / (2 * nj), nj)  # base: 1/(2*nj) for each
    if (k %% 2 == 1) p[nj] <- 0.5      # odd: heavy on last
    else              p[1] <- 0.5       # even: heavy on first
    p <- p / sum(p)
    setNames(p, seq_len(nj))
  })
  names(tgt) <- paste0("v", seq_len(K))
  list(df=df, tgt=tgt, compression=n / (nj^K))
}

run_timed <- function(df, tgt, jacobi, max_w) {
  t0 <- proc.time()["elapsed"]
  r <- suppressWarnings(harvest(df, tgt, method="ieppa",
    max_weight=max_w, min_weight=0, max_iterations=100L,
    convergence=list(absolute=1e-4), jacobi_sweep=jacobi,
    verbose=0L, attach_weights=FALSE))
  wall <- proc.time()["elapsed"] - t0
  res <- attr(r, "result")
  list(wall=wall, iters=res$iterations, status=res$status, max_err=res$max_error)
}

Ks    <- c(3L, 6L, 9L)
comps <- c(2, 10, 50)
tight_cases <- list(
  list(max_w=5.0, label="loose"),
  list(max_w=2.0, label="tight")
)
n <- 50000L

results <- list()
fmt <- "%-4s %-6s %-5s %-5d %-5d %-6d %-6d %-7.3f %-7.3f  ratio=%.3f  iflat=%.2f\n"
cat(sprintf("%-4s %-6s %-5s %-5s %-5s %-6s %-6s %-7s %-7s\n",
    "K","compr","tight","st_gs","st_ja","it_gs","it_ja","wall_gs","wall_ja"))

for (K in Ks) {
  for (comp in comps) {
    prob <- make_adversarial(n, K, comp)
    for (tc in tight_cases) {
      gs <- run_timed(prob$df, prob$tgt, FALSE, tc$max_w)
      ja <- run_timed(prob$df, prob$tgt, TRUE,  tc$max_w)
      wall_rat <- if (gs$wall > 0) ja$wall / gs$wall else NA
      iter_inf <- if (gs$iters > 0) ja$iters / gs$iters else NA
      cat(sprintf(fmt, K, round(prob$compression, 0), tc$label,
          gs$status, ja$status, gs$iters, ja$iters,
          gs$wall, ja$wall, wall_rat, iter_inf))
      results[[length(results)+1]] <- list(
        K=K, compression=prob$compression, tightness=tc$label, path="log",
        status_gs=gs$status, status_ja=ja$status,
        iters_gs=gs$iters, iters_ja=ja$iters,
        wall_gs=gs$wall, wall_ja=ja$wall,
        iter_inflation=iter_inf, wall_ratio=wall_rat)
    }
  }
}

res_df <- do.call(rbind, lapply(results, as.data.frame))
saveRDS(res_df, "benchmarks/jacobi_sweep_study.rds")

log_cells <- res_df[!is.na(res_df$wall_ratio), ]
med_ratio <- median(log_cells$wall_ratio)
med_inf   <- median(log_cells$iter_inflation, na.rm=TRUE)
cat(sprintf("\n=== Summary ===\nMedian wall_ratio: %.3f\nMedian iter_inflation: %.3f\nDecision: %s\n",
    med_ratio, med_inf,
    ifelse(med_ratio < 1.0, "SHIP Jacobi as default", "KEEP GS / conditional")))
