# BLSE 3x3 Routing Sweep Implementation Plan

**Goal:** Localize the AUTO-router breakpoint between ieppa and raking by sweeping `(M_cell/n, n)` on a 3x3 grid and recording per-cell winner.
**Architecture:** AUTO dispatch currently routes on `M_cell/n` ratio (kLinearSpaceThreshold ≈ 2.0 per ticket context). A single-point probe at `M_cell/n=1.0` cannot localize the threshold; this benchmark records iters and per-iter cost across `M_cell/n ∈ {0.01, 0.1, 0.5}` × `n ∈ {1e4, 1e5, 1e6}` to derive numeric T1 (ratio) and T2 (n floor).
**Tech Stack:** R, `bench::mark()`, leafblower::harvest, no C++ changes.

**Mechanism:** R-only benchmark script writing rds + a small report.
**Forbidden:** changing AUTO dispatch logic in this ticket; using `system.time()` in place of `bench::mark()`; pivoting to fix the threshold here (separate ticket).

---

## Task T1: Author `benchmarks/blse_routing_sweep.R`

Steps:

1. Read `benchmarks/ieppa_vs_raking_bench.R` (existing pattern). Reuse its data-generation helper if present; otherwise:

```r
# benchmarks/blse_routing_sweep.R
# Sweep (M_cell/n, n) over 3x3 grid; record per-method iters, time, max_err.
# Output: benchmarks/blse_routing_sweep.rds + a markdown summary.
suppressPackageStartupMessages({
  library(leafblower); library(bench); library(dplyr)
})

set.seed(20260501)

make_problem <- function(n, mc_ratio) {
  M_cell_target <- max(1L, as.integer(round(mc_ratio * n)))
  # Pick K and cat_counts so prod(cat_counts) ~ M_cell_target.
  # Use K=4 with balanced category counts.
  K  <- 4L
  cc <- rep(max(2L, as.integer(round(M_cell_target ^ (1/K)))), K)
  group_ids <- replicate(K, sample.int(cc[1], n, replace = TRUE) - 1L,
                         simplify = FALSE)
  targets   <- lapply(cc, function(k) {
    p <- runif(k); p / sum(p)
  })
  list(weights = rep(1.0, n), group_ids = group_ids,
       cat_counts = cc, targets = targets,
       n = n, K = K, M_cell_target = M_cell_target,
       mc_ratio = mc_ratio)
}

run_one <- function(prob, method) {
  t <- bench::mark(
    res = harvest(prob$weights, prob$group_ids, prob$cat_counts, prob$targets,
                  method = method, max_iter = 500L, tol_pct = 1e-4,
                  verbose = 0L),
    iterations = 5L, check = FALSE)
  data.frame(
    n         = prob$n,
    mc_ratio  = prob$mc_ratio,
    method    = method,
    median_s  = as.numeric(t$median),
    mem_alloc = as.numeric(t$mem_alloc),
    iters     = res$iterations,
    max_err   = res$max_err,
    status    = res$status)
}

grid <- expand.grid(
  n        = c(1e4, 1e5, 1e6),
  mc_ratio = c(0.01, 0.1, 0.5))

results <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  prob <- make_problem(as.integer(grid$n[i]), grid$mc_ratio[i])
  rbind(
    run_one(prob, "ieppa"),
    run_one(prob, "raking"))
}))

saveRDS(results, "benchmarks/blse_routing_sweep.rds")

# Per-cell winner
winner <- results |>
  group_by(n, mc_ratio) |>
  summarise(winner = method[which.min(median_s)],
            speedup = max(median_s) / min(median_s),
            .groups = "drop")

cat("\n## BLSE routing winner per cell\n\n")
print(winner)
saveRDS(winner, "benchmarks/blse_routing_sweep_winner.rds")
```

Confidence: 80 — pattern matches existing `benchmarks/ieppa_vs_raking_bench.R`; need to verify exact harvest() arg names (max_iter vs maxit) at write time.

---

## Task T2: Run + record

```bash
Rscript benchmarks/blse_routing_sweep.R 2>&1 | tee benchmarks/blse_routing_sweep.log
```

Pass criteria:
- 9 cells × 2 methods = 18 rows in the rds.
- Each row has finite `median_s` and `max_err < 1e-2` (i.e. both methods solved).
- Markdown winner table written to log.

If any cell fails to converge for either method, increase `max_iter` to 2000 and re-run that cell. Do NOT silently lower tolerance.

---

## Task T3: Append summary to ticket

`bd update leafblower-2jw --comment "$(cat benchmarks/blse_routing_sweep.log)"`

The summary table is the deliverable. Setting numeric T1/T2 thresholds is **out of scope** — file a follow-up ticket if the sweep recommends a value different from current `kLinearSpaceThreshold`.
