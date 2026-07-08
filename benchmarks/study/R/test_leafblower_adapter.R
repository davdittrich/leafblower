#!/usr/bin/env Rscript
# benchmarks/study/R/test_leafblower_adapter.R -- WU-7 (leafblower-2ouc.8)
# Contract-shape + golden tests for benchmarks/study/R/leafblower_adapter.R.
#
# Usage (repo root, single-thread BLAS per CLAUDE.md determinism rule):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#     Rscript benchmarks/study/R/test_leafblower_adapter.R

here <- function(...) file.path("benchmarks", "study", ...)
source(here("common", "problem_io.R"))
source(here("common", "metrics.R"))
source(here("R", "leafblower_adapter.R"))

fail_count <- 0L
check <- function(desc, cond) {
  ok <- isTRUE(cond)
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", desc))
  if (!ok) fail_count <<- fail_count + 1L
  invisible(ok)
}
near <- function(a, b, tol = 1e-4) all(abs(a - b) <= tol)

# Currently-installed build (see .wolf/... / final report for which flags are
# in effect). The portable<->native flag-swap-and-rebuild cycle is a separate
# manual DoD step (spec/contract.md Section 3); this suite validates adapter
# CORRECTNESS against whichever build is installed when it runs, tagging its
# own weights_ref filenames accordingly.
BUILD <- Sys.getenv("LBW_BUILD_TAG", "native")

toy <- load_problem_spec(here("spec", "toy_inline.json"))

KEYS <- c("weights_ref", "iterations", "status", "converged", "error_message",
          "wall_time_s", "peak_rss_bytes")
STATUS_ENUM <- c("converged", "no_conv", "infeasible", "bound_violation",
                  "bad_arg", "budget", "stall", "error")

## ---------------------------------------------------------------------------
## Contract-shape checks (all 9 methods)
## ---------------------------------------------------------------------------
cat("== contract-shape (all 9 methods) ==\n")
for (m in .LBW_METHODS) {
  res <- run_leafblower(toy, m, BUILD)
  check(sprintf("%s: return has exactly the 7 contract keys", m),
        identical(sort(names(res)), sort(KEYS)))
  check(sprintf("%s: status in harmonized enum", m), res$status %in% STATUS_ENUM)
  check(sprintf("%s: weights_ref file exists on disk", m),
        file.exists(file.path(.REPO_ROOT, res$weights_ref)))
  w <- arrow::read_parquet(file.path(.REPO_ROOT, res$weights_ref))$weight
  check(sprintf("%s: weights_ref has n=%d rows", m, nrow(toy$data)),
        length(w) == nrow(toy$data))
  check(sprintf("%s: wall_time_s is a positive finite number", m),
        is.numeric(res$wall_time_s) && is.finite(res$wall_time_s) && res$wall_time_s > 0)
  check(sprintf("%s: peak_rss_bytes is NA or a positive integer", m),
        is.na(res$peak_rss_bytes) || res$peak_rss_bytes > 0)
  check(sprintf("%s: converged is a logical scalar", m),
        is.logical(res$converged) && length(res$converged) == 1L && !is.na(res$converged))
}

## ---------------------------------------------------------------------------
## LEAFBLOWER_R_ADAPTERS registry: 9 entries, leafblower_<method>_r keys
## ---------------------------------------------------------------------------
cat("\n== registry ==\n")
check("registry has 9 entries", length(LEAFBLOWER_R_ADAPTERS) == 9L)
check("registry keys are leafblower_<method>_r",
      identical(sort(names(LEAFBLOWER_R_ADAPTERS)), sort(paste0("leafblower_", .LBW_METHODS, "_r"))))
reg_res <- LEAFBLOWER_R_ADAPTERS[["leafblower_raking_r"]](toy, BUILD)
check("registry entry is callable and contract-shaped",
      identical(sort(names(reg_res)), sort(KEYS)))

## ---------------------------------------------------------------------------
## Golden: toy_inline hand-derived closed form.
##
## toy_inline (spec/toy_inline.json): single margin "grp", 2 disjoint groups
## A={row1,row2}, B={row3,row4}, design_weights=(1,1,2,2), targets A=B=0.5.
## leafblower enforces Sigma(w) = n = 4 (row count, CLAUDE.md "Output
## weights" convention) at exit, regardless of sum(design_weights)=6.
##
## For a SINGLE margin split into disjoint groups, the constraint set is two
## independent linear equalities sum_{i in A} w_i = 0.5*4 = 2 and
## sum_{i in B} w_i = 0.5*4 = 2. Any KL/entropic-family projection from a
## starting point d onto a per-group-sum-only constraint is the UNIQUE
## multiplicative rescale within each group: w_i = d_i * (target_group_mass /
## d_group_sum). Here scale_A = 2 / (1+1) = 1, scale_B = 2 / (2+2) = 0.5, so
##   w = (1*1, 1*1, 2*0.5, 2*0.5) = (1, 1, 1, 1) EXACTLY.
## This closed form holds for every solver in KL_NATIVE_FAMILIES
## (metrics.R): oris, oris_soft, raking, sinkhorn, greenkhorn, newton_kl.
## logit/chebyshev/greg use a different divergence with no assumed common
## closed form here, so they get a looser (but still fully derived) check:
## convergence achieved, Sigma(w)=n, bounds respected.
## ---------------------------------------------------------------------------
cat("\n== golden: toy_inline hand-derived closed form ==\n")
kl_native_methods <- intersect(.LBW_METHODS, KL_NATIVE_FAMILIES)
other_methods <- setdiff(.LBW_METHODS, KL_NATIVE_FAMILIES)

for (m in kl_native_methods) {
  res <- run_leafblower(toy, m, BUILD)
  w <- arrow::read_parquet(file.path(.REPO_ROOT, res$weights_ref))$weight
  check(sprintf("%s: recomputed converged == TRUE", m), isTRUE(res$converged))
  check(sprintf("%s: weights == (1,1,1,1) (exact KL-family projection)", m),
        near(w, c(1, 1, 1, 1), tol = 1e-4))
}

for (m in other_methods) {
  res <- run_leafblower(toy, m, BUILD)
  w <- arrow::read_parquet(file.path(.REPO_ROOT, res$weights_ref))$weight
  check(sprintf("%s: recomputed converged == TRUE", m), isTRUE(res$converged))
  check(sprintf("%s: Sigma(w) == n == 4", m), near(sum(w), 4, tol = 1e-4))
  check(sprintf("%s: bounds respected", m),
        all(w >= toy$bounds$min - 1e-9) && all(w <= toy$bounds$max + 1e-9))
}

## ---------------------------------------------------------------------------
## Error classification: bad_arg (R-level pre-C-call validation stop(), no
## "leafblower:" prefix -- max_weight < min_weight is rejected before the
## solver ever runs).
## ---------------------------------------------------------------------------
cat("\n== error classification ==\n")
bad_problem <- toy
bad_problem$bounds <- list(min = 10, max = 0)
res_bad <- run_leafblower(bad_problem, "raking", BUILD)
check("max_weight<min_weight -> status bad_arg", identical(res_bad$status, "bad_arg"))
check("bad_arg run converged == FALSE", identical(res_bad$converged, FALSE))
check("bad_arg run still writes a weights_ref (never dangling)",
      file.exists(file.path(.REPO_ROOT, res_bad$weights_ref)))
w_bad <- arrow::read_parquet(file.path(.REPO_ROOT, res_bad$weights_ref))$weight
check("bad_arg sentinel is all-NaN, length n", length(w_bad) == nrow(toy$data) && all(is.na(w_bad)))
check("bad_arg error_message is non-empty", nzchar(res_bad$error_message %||% ""))

cat(sprintf("\n%s: %d assertion(s) failed.\n", if (fail_count == 0) "GOLDEN PASS" else "GOLDEN FAIL", fail_count))
quit(status = if (fail_count == 0) 0L else 1L, save = "no")
