#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(leafblower))

Sys.setenv(LBW_IEPPA_FORCE_PATH = "log")

# Compression control: M_cell ≈ n * mc_ratio
# K margins, nj categories each: M_cell ≈ nj^K → nj = (n * mc_ratio)^(1/K)
make_problem <- function(n, K, compression, seed = 42L) {
  set.seed(seed)
  mc_ratio <- 1 / compression
  nj <- max(2L, round((n * mc_ratio)^(1/K)))
  df <- as.data.frame(lapply(seq_len(K), function(k) factor(sample(nj, n, replace=TRUE))))
  names(df) <- paste0("v", seq_len(K))
  tgt <- setNames(lapply(seq_len(K), function(k) setNames(rep(1/nj, nj), seq_len(nj))),
                  paste0("v", seq_len(K)))
  list(df=df, tgt=tgt, M_cell_approx=nj^K, compression=n/(nj^K))
}

run_timed <- function(df, tgt, jacobi, max_w) {
  t0 <- proc.time()["elapsed"]
  r <- suppressWarnings(harvest(df, tgt, method="ieppa",
                                max_weight=max_w, min_weight=0,
                                max_iterations=1000L,
                                convergence=list(absolute=1e-4),
                                jacobi_sweep=jacobi, verbose=0L,
                                attach_weights=FALSE))
  wall <- proc.time()["elapsed"] - t0
  res <- attr(r, "result")
  list(wall=wall, iters=res$iterations, status=res$status, max_err=res$max_error)
}

Ks     <- c(3L, 6L, 9L)
comps  <- c(2, 10, 50)
tight  <- list(loose=list(max_w=5.0, label="loose"),
               tight=list(label="tight", max_w=1.5))
n      <- 30000L

results <- list()
cat(sprintf("%-6s %-6s %-6s %-5s %-5s %-6s %-6s %-8s %-8s\n",
    "K","compr","tight","st_gs","st_ja","it_gs","it_ja","wall_gs","wall_ja"))

for (K in Ks) {
  for (comp in comps) {
    prob <- make_problem(n, K, comp)
    for (t in tight) {
      gs <- run_timed(prob$df, prob$tgt, FALSE, t$max_w)
      ja <- run_timed(prob$df, prob$tgt, TRUE,  t$max_w)
      iter_inf <- if (gs$iters > 0) ja$iters / gs$iters else NA
      wall_rat <- if (gs$wall  > 0) ja$wall  / gs$wall  else NA
      cat(sprintf("%-6d %-6.0fx %-6s %-5d %-5d %-6d %-6d %-8.3f %-8.3f  ratio=%.3f\n",
          K, prob$compression, t$label,
          gs$status, ja$status, gs$iters, ja$iters,
          gs$wall, ja$wall, wall_rat))
      results[[length(results)+1]] <- list(
        K=K, compression=prob$compression, tightness=t$label,
        path="log",  # forced by env var
        status_gs=gs$status, status_ja=ja$status,
        iters_gs=gs$iters, iters_ja=ja$iters,
        wall_gs=gs$wall, wall_ja=ja$wall,
        iter_inflation=iter_inf, wall_ratio=wall_rat
      )
    }
  }
}

res_df <- do.call(rbind, lapply(results, as.data.frame))
saveRDS(res_df, "benchmarks/jacobi_sweep_study.rds")
cat("\n=== Decision summary ===\n")
log_cells <- res_df[res_df$path == "log" & !is.na(res_df$wall_ratio), ]
cat(sprintf("Log-path cells: %d\n", nrow(log_cells)))
cat(sprintf("Median wall_ratio: %.3f\n", median(log_cells$wall_ratio)))
cat(sprintf("Median iter_inflation: %.3f\n", median(log_cells$iter_inflation, na.rm=TRUE)))
cat(sprintf("Decision gate (wall_ratio < 1.0): %s\n",
    ifelse(median(log_cells$wall_ratio) < 1.0, "SHIP Jacobi as default", "KEEP GS (revert flag)")))
