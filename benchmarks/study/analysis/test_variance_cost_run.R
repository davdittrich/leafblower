#!/usr/bin/env Rscript
# STUDY-BRANCH-ONLY-DO-NOT-MERGE
# test_variance_cost_run.R -- WU-13 (leafblower-2ouc.14.2) contract + audit
# tests for benchmarks/study/analysis/variance_cost_run.R.
#
# Usage (repo root, single-thread BLAS per CLAUDE.md determinism rule):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#     Rscript benchmarks/study/analysis/test_variance_cost_run.R

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

here <- function(...) file.path("benchmarks", "study", ...)
source(here("analysis", "replicate_weights.R"))
source(here("common", "problem_io.R"))
source(here("common", "instance_family.R"))
install_gen_resolver()
source(here("analysis", "variance_cost_run.R"))

fail_count <- 0L
check <- function(desc, cond) {
  ok <- isTRUE(cond)
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", desc))
  if (!ok) fail_count <<- fail_count + 1L
  invisible(ok)
}

## ---------------------------------------------------------------------------
## Fixture: smallest n=1000 instance-family spec on disk. SKIP cleanly (not a
## failure) if none is present.
## ---------------------------------------------------------------------------
spec_dir <- here("spec", "instance_family")
n1000_specs <- sort(list.files(spec_dir, pattern = "^if_n1000_.*\\.json$", full.names = TRUE))

if (length(n1000_specs) == 0L) {
  cat("SKIP: no if_n1000_*.json instance-family spec found under ", spec_dir, " -- nothing to test.\n")
  quit(status = 0L, save = "no")
}

spec_path <- n1000_specs[[1]]
cat(sprintf("using fixture spec: %s\n", spec_path))
problem <- load_problem_spec(spec_path)
bounded <- is.finite(problem$bounds$min) && is.finite(problem$bounds$max)

## ---------------------------------------------------------------------------
## 1-2. Basic contract shape: 4 rows, expected columns, finite non-negative
## wall times, status/converged populated.
## ---------------------------------------------------------------------------
cat("== basic contract shape ==\n")
B <- 4L
cell <- rm_variance_cost_cell(problem, "leafblower_raking_r", B = B, seed = 123L)

EXPECTED_COLS <- c("solver", "problem", "rep", "wall_time_s", "status",
                    "converged", "iterations", "repwt_fingerprint")

check("rm_variance_cost_cell returns a data.frame", is.data.frame(cell))
check(sprintf("returns exactly B=%d rows", B), nrow(cell) == B)
check("has all expected columns", all(EXPECTED_COLS %in% names(cell)))
check("rep column is 1..B", identical(cell$rep, seq_len(B)))
check("solver column is the requested solver_id",
      all(cell$solver == "leafblower_raking_r"))
check("problem column equals problem$id", all(cell$problem == problem$id))
check("wall_time_s all finite and >= 0",
      all(is.finite(cell$wall_time_s)) && all(cell$wall_time_s >= 0))
check("status is a non-empty character for every row",
      is.character(cell$status) && all(nzchar(cell$status)))
check("converged is logical, non-NA for every row",
      is.logical(cell$converged) && !any(is.na(cell$converged)))
check("repwt_fingerprint is a non-empty character, identical across rows",
      is.character(cell$repwt_fingerprint) && all(nzchar(cell$repwt_fingerprint)) &&
        length(unique(cell$repwt_fingerprint)) == 1L)

## ---------------------------------------------------------------------------
## 3-4. Fair-comparison audit: a second solver on the SAME (problem, seed)
## must get the IDENTICAL repwt_fingerprint (same replicate weights).
## ---------------------------------------------------------------------------
cat("\n== fair-comparison audit (shared replicate weights across solvers) ==\n")
second_solver <- if (bounded) "leafblower_logit_r" else "survey_calibrate_raking"
cell2 <- rm_variance_cost_cell(problem, second_solver, B = B, seed = 123L)

check(sprintf("second solver ('%s') also returns B=%d rows", second_solver, B), nrow(cell2) == B)
check("repwt_fingerprint IDENTICAL across the two different solvers (same problem+seed)",
      identical(cell$repwt_fingerprint[[1]], cell2$repwt_fingerprint[[1]]))

## ---------------------------------------------------------------------------
## 5. Determinism: same (problem, solver, seed) -> identical fingerprint;
## different seed -> different fingerprint.
## ---------------------------------------------------------------------------
cat("\n== determinism ==\n")
cell_repeat <- rm_variance_cost_cell(problem, "leafblower_raking_r", B = B, seed = 123L)
cell_diff_seed <- rm_variance_cost_cell(problem, "leafblower_raking_r", B = B, seed = 456L)

check("same (problem, solver, seed) -> identical repwt_fingerprint",
      identical(cell$repwt_fingerprint[[1]], cell_repeat$repwt_fingerprint[[1]]))
check("different seed -> different repwt_fingerprint",
      !identical(cell$repwt_fingerprint[[1]], cell_diff_seed$repwt_fingerprint[[1]]))

## ---------------------------------------------------------------------------
## 6. COLD check: variance_cost_run.R must never pass start_weights= to
## harvest() -- each replicate re-solve starts from the reference measure,
## no warm start carried over from a previous replicate or solve.
## ---------------------------------------------------------------------------
cat("\n== cold-solve check ==\n")
src_path <- here("analysis", "variance_cost_run.R")
src_lines <- readLines(src_path)
# Strip comments before checking -- the source file's own header PROSE
# documents the cold-solve invariant by naming start_weights (e.g. "no
# start_weights= argument passed"), which must not itself trip this check;
# only an actual code reference (harvest(..., start_weights = ...)) should.
code_only <- sub("#.*$", "", src_lines)
check("variance_cost_run.R never references start_weights in actual code (cold-solve only)",
      !any(grepl("start_weights", code_only, fixed = TRUE)))

cat(sprintf("\n%s: %d assertion(s) failed.\n", if (fail_count == 0) "ALL PASS" else "FAIL", fail_count))
quit(status = if (fail_count == 0) 0L else 1L, save = "no")
