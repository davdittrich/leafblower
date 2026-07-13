#!/usr/bin/env Rscript
# STUDY-BRANCH-ONLY-DO-NOT-MERGE
# test_replicate_weights.R -- tests for benchmarks/study/analysis/replicate_weights.R
# (WU leafblower-2ouc.14.1).
#
# Usage: Rscript benchmarks/study/analysis/test_replicate_weights.R

here <- function(...) file.path("benchmarks", "study", ...)
source(here("analysis", "replicate_weights.R"))

failures <- 0L
check <- function(desc, cond) {
  if (isTRUE(cond)) {
    cat(sprintf("  PASS: %s\n", desc))
  } else {
    cat(sprintf("  FAIL: %s\n", desc))
    failures <<- failures + 1L
  }
}

set.seed(42L)
n <- 2000L
design_weights <- runif(n, 0.5, 2)
B <- 20L
seed <- 12345L

cat("== output shape ==\n")
W <- rm_make_replicate_weights(design_weights, B, seed)
check("is a matrix", is.matrix(W))
check("nrow == n", nrow(W) == n)
check("ncol == B", ncol(W) == B)

cat("== finite and non-negative ==\n")
check("all entries finite", all(is.finite(W)))
check("all entries >= 0", all(W >= 0))

cat("== mean-preservation (bootstrap noise tolerance) ==\n")
target_mean <- mean(design_weights)
tol <- 0.05 * target_mean
col_means <- colMeans(W)
check("all column means within tol of design mean",
      all(abs(col_means - target_mean) <= tol))
cat(sprintf("  (target_mean=%.6f, range of col_means=[%.6f, %.6f], tol=%.6f)\n",
            target_mean, min(col_means), max(col_means), tol))

cat("== byte-determinism ==\n")
W1 <- rm_make_replicate_weights(design_weights, B, seed)
W2 <- rm_make_replicate_weights(design_weights, B, seed)
check("same seed => identical matrices", identical(W1, W2))
W3 <- rm_make_replicate_weights(design_weights, B, seed + 1L)
check("different seed => different matrix", !identical(W1, W3))

cat("== RNG non-perturbation ==\n")
set.seed(999L)
snapshot_before <- .Random.seed
invisible(rm_make_replicate_weights(design_weights, B, seed))
snapshot_after <- .Random.seed
check("caller's .Random.seed unchanged after call",
      identical(snapshot_before, snapshot_after))
# also confirm the *stream* is unperturbed: draw a value, reset, call fn, draw again
set.seed(999L)
draw_before <- runif(1)
set.seed(999L)
invisible(rm_make_replicate_weights(design_weights, B, seed))
draw_after <- runif(1)
check("caller's RNG stream unperturbed (same next draw)",
      identical(draw_before, draw_after))

cat("== arg validation ==\n")
bad_weights_zero <- design_weights
bad_weights_zero[1] <- 0
err1 <- tryCatch({
  rm_make_replicate_weights(bad_weights_zero, B, seed)
  NULL
}, error = function(e) conditionMessage(e))
check("non-positive weight raises error", is.character(err1))

bad_weights_neg <- design_weights
bad_weights_neg[1] <- -1
err2 <- tryCatch({
  rm_make_replicate_weights(bad_weights_neg, B, seed)
  NULL
}, error = function(e) conditionMessage(e))
check("negative weight raises error", is.character(err2))

bad_weights_nonfinite <- design_weights
bad_weights_nonfinite[1] <- NA_real_
err3 <- tryCatch({
  rm_make_replicate_weights(bad_weights_nonfinite, B, seed)
  NULL
}, error = function(e) conditionMessage(e))
check("non-finite weight raises error", is.character(err3))

err4 <- tryCatch({
  rm_make_replicate_weights(c(1.0), B, seed)
  NULL
}, error = function(e) conditionMessage(e))
check("n < 2 raises error", is.character(err4))

err5 <- tryCatch({
  rm_make_replicate_weights(design_weights, B, seed, type = "JK1")
  NULL
}, error = function(e) conditionMessage(e))
check("unsupported type raises error", is.character(err5) && grepl("unsupported type", err5))

cat(sprintf("\nRESULT: %d failure(s)\n", failures))
if (failures > 0L) quit(status = 1L)
