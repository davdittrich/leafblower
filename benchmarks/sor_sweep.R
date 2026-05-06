#!/usr/bin/env Rscript
# benchmarks/sor_sweep.R
# SOR cost/benefit sweep: 3×3 design (bounds × K) × SOR on/off
# Purpose: confirm whether SOR default (auto=TRUE) is net positive across its design envelope.

suppressPackageStartupMessages(library(leafblower))
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ─── DGP ─────────────────────────────────────────────────────────────────────
make_oscillating_frame <- function(n, K, max_wt) {
  # Opposite sampling/target probs (40pp mismatch) forces large weight corrections
  # and pushes weights against the max_weight boundary
  set.seed(42 + 1000L * K + as.integer(max_wt * 10))  # seed = 42 + K*1000 + max_wt*10
  nms <- paste0("v", seq_len(K))
  tgt <- lapply(seq_len(K), function(k) c(A=0.7, B=0.3))
  names(tgt) <- nms
  df  <- as.data.frame(lapply(nms, function(v)
    factor(sample(c("A","B"), n, replace=TRUE, prob=c(0.3,0.7)), levels=names(tgt[[v]]))
  ))
  names(df) <- nms
  list(df=df, tgt=tgt)
}

# Pre-validate DGP: check weights hit boundary
validate_dgp <- function(df, tgt, max_wt) {
  res <- harvest(df, tgt, method="ieppa", max_weight=max_wt,
                 convergence=list(tol=1e-4), verbose=0L, attach_weights=FALSE)
  boundary_pct <- mean(res >= max_wt * 0.9)
  if (boundary_pct < 0.01) {
    message(sprintf("  Warning: only %.1f%% weights near boundary; DGP may not stress SOR", boundary_pct*100))
  } else {
    cat(sprintf("  DGP boundary check: %.1f%% weights > 90%% of max_wt\n", boundary_pct*100))
  }
  invisible(boundary_pct)
}

# ─── Grid ────────────────────────────────────────────────────────────────────
bounds_grid <- c(tight=2, medium=5, loose=10)
k_grid      <- c(low=3, medium=6, high=9)
NREPS       <- 10L
CONV        <- list(tol=1e-4)
N           <- 10000L

SOR_ON  <- list(auto=TRUE,  omega_min=0.3)
SOR_OFF <- NULL  # NULL disables SOR (parse_sor(NULL) → enabled=0)

cat(sprintf("SOR sweep: %d bounds × %d K × 2 variants × %d reps = %d calls\n",
            length(bounds_grid), length(k_grid), 2L, length(bounds_grid)*length(k_grid)*2L*NREPS))

results <- vector("list", length(bounds_grid) * length(k_grid))
idx <- 0L
for (bname in names(bounds_grid)) {
  max_wt <- bounds_grid[[bname]]
  for (kname in names(k_grid)) {
    K <- k_grid[[kname]]
    idx <- idx + 1L
    label <- sprintf("bounds=%s(%.0f) K=%s(%d)", bname, max_wt, kname, K)
    cat(sprintf("\n[%d/9] %s\n", idx, label))

    frm  <- make_oscillating_frame(N, K, max_wt)
    df   <- frm$df; tgt <- frm$tgt
    validate_dgp(df, tgt, max_wt)

    run_n <- function(sor_cfg) {
      times <- numeric(NREPS); iters <- integer(NREPS)
      for (i in seq_len(NREPS)) {
        t0 <- proc.time()["elapsed"]
        res <- harvest(df, tgt, method="ieppa", max_weight=max_wt,
                       convergence=CONV, sor=sor_cfg, verbose=0L, attach_weights=FALSE)
        times[i] <- (proc.time()["elapsed"] - t0) * 1000
        iters[i] <- attr(res, "iterations") %||% NA_integer_
      }
      list(med_ms=median(times), med_iters=median(iters))
    }

    on  <- run_n(SOR_ON)
    off <- run_n(SOR_OFF)

    delta_iters <- (on$med_iters - off$med_iters) / pmax(off$med_iters, 1)
    delta_time  <- on$med_ms / pmax(off$med_ms, 0.001)
    hurts <- (delta_iters >= 0.10) || (delta_time > 1.10)

    cat(sprintf("  SOR_on:  %5.0fms %4.0f iters\n", on$med_ms, on$med_iters))
    cat(sprintf("  SOR_off: %5.0fms %4.0f iters\n", off$med_ms, off$med_iters))
    cat(sprintf("  Δiters=%.2f  time_ratio=%.2f  hurts=%s\n",
                delta_iters, delta_time, if(hurts)"YES" else "no"))

    results[[idx]] <- data.frame(
      condition=label, bounds=bname, max_wt=max_wt, K_level=kname, K=K,
      sor_on_ms=round(on$med_ms,1), sor_off_ms=round(off$med_ms,1),
      sor_on_iters=on$med_iters, sor_off_iters=off$med_iters,
      delta_iters=round(delta_iters,3), time_ratio=round(delta_time,3),
      hurts=hurts, stringsAsFactors=FALSE
    )
  }
}

tab <- do.call(rbind, results)
n_hurts <- sum(tab$hurts)
cat(sprintf("\n=== SOR sweep complete: %d/9 conditions SOR hurts (threshold >=3) ===\n", n_hurts))
cat(sprintf("Verdict: %s\n",
    if(n_hurts >= 3) "SOR HURTS in >=3/9 conditions — file follow-up to reconsider default"
    else "SOR is net positive or neutral — default confirmed"))
print(tab[,c("condition","sor_on_iters","sor_off_iters","delta_iters","time_ratio","hurts")],
      row.names=FALSE)

outfile <- "benchmarks/results/sor_sweep.csv"
write.csv(tab, outfile, row.names=FALSE)
cat(sprintf("Saved: %s\n", outfile))
