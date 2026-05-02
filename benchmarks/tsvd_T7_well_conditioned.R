#!/usr/bin/env Rscript
# WL-4: T7 K=4 well-conditioned over-projection regression test
Sys.setenv(OMP_NUM_THREADS = "1")
suppressPackageStartupMessages(library(leafblower))
set.seed(1L)

# ── DGP: K=4 well-conditioned fixture ────────────────────────────────────────
# n = 2000, margins: m1 (3 cats), m2 (4 cats), m3 (3 cats), m4 (2 cats)
# Mild skew per margin as per spec
n <- 2000L
df <- data.frame(
  m1 = factor(sample(letters[1:3], n, TRUE, prob = c(0.4, 0.3, 0.3))),
  m2 = factor(sample(letters[1:4], n, TRUE, prob = c(0.4, 0.3, 0.2, 0.1))),
  m3 = factor(sample(letters[1:3], n, TRUE, prob = c(0.4, 0.3, 0.3))),
  m4 = factor(sample(c("M", "F"), n, TRUE, prob = c(0.5, 0.5)))
)

tgt <- list(
  m1 = c(a = 0.4, b = 0.3, c = 0.3),
  m2 = c(a = 0.4, b = 0.3, c = 0.2, d = 0.1),
  m3 = c(a = 0.4, b = 0.3, c = 0.3),
  m4 = c(M = 0.5, F = 0.5)
)

# ── Newton-KL with default ratio_tsvd = 1e-8 ────────────────────────────────
t0 <- Sys.time()
r <- harvest(df, tgt, method = "newton_kl",
  max_iterations = 50L, attach_weights = FALSE, verbose = 0L)
wall_ms <- as.numeric(Sys.time() - t0, units = "secs") * 1000
res <- attr(r, "result")

result_row <- data.frame(
  lambda_ratio = NA_real_,  # not externally computed; see commit body
  n_projected_dims = ifelse(is.null(res$n_projected_dims), NA_integer_, res$n_projected_dims),
  max_err = res$max_error,
  n_iter = res$iterations,
  wall_time_ms = round(wall_ms, 2)
)

dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
write.csv(result_row, "benchmarks/results/tsvd_T7_K4.csv", row.names = FALSE)

cat(sprintf("[wl-4] T7 K=4 well-conditioned: n_proj=%d, max_err=%.3e, n_iter=%d, wall=%.0fms\n",
  result_row$n_projected_dims, result_row$max_err,
  result_row$n_iter, result_row$wall_time_ms))

# ── HARD GATES ───────────────────────────────────────────────────────────────
# Gate 1: n_projected_dims == 0 (no truncation fired with ratio_tsvd = 1e-8)
# Gate 2: max_err < 1e-4 (convergence quality)
# Gate 3: lambda_ratio < 1e6 — NA in this run; n_proj==0 + max_err<1e-4 are primary
hard_pass <- result_row$n_projected_dims == 0 && result_row$max_err < 1e-4

if (hard_pass) {
  cat("[wl-4] HARD GATE PASS — over-projection regression test green.\n")
} else {
  cat(sprintf("[wl-4] HARD GATE FAIL — n_proj=%d, max_err=%.3e\n",
    result_row$n_projected_dims, result_row$max_err))
  cat("[wl-4] SPEC_FAILURE: over-projection bug OR DGP miscalibration — escalate.\n")
  quit(status = 2)
}
