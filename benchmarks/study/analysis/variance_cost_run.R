#!/usr/bin/env Rscript
# STUDY-BRANCH-ONLY-DO-NOT-MERGE
# variance_cost_run.R -- WU-13 (leafblower-2ouc.14.2) COLD per-replicate
# re-solve cost of a solver on a problem: for each of B Rao-Wu-Yue rescaled
# bootstrap replicate design-weight vectors (leafblower-2ouc.14.1's
# rm_make_replicate_weights()), COLD-solve the same problem and record the
# per-replicate wall time and outcome. This is the "variance axis" cost --
# how much wall time producing B replicate estimates for a variance
# calculation actually costs, as opposed to the single-solve cost measured
# by benchmarks/study/R/*.R's contract-v2 harness.
#
# Timer discipline mirrors .comp_run (competitors.R) / run_leafblower
# (leafblower_adapter.R) EXACTLY: t0 <- Sys.time() immediately before the
# solve, wall <- difftime(Sys.time(), t0, "secs") immediately after -- no
# parquet/weights I/O inside the timed region.
#
# COLD means: each replicate re-solve starts from the reference measure
# (design_weights = that replicate's weight vector) with NO warm-start
# (no start_weights= argument passed to harvest()). This file intentionally
# never passes start_weights to leafblower::harvest().
#
# Fair comparison: the SAME n x B replicate-weight matrix (same problem,
# same seed) is generated once per rm_variance_cost_cell() call and reused
# for every replicate of that call -- calling this function for two
# different solver_ids with the same (problem, seed) yields byte-identical
# replicate weights (see `repwt_fingerprint` below), which is the audit key
# for cross-solver fairness.

# Single-thread BLAS for determinism (CLAUDE.md) -- MUST be set before any
# BLAS-linked package (leafblower, arrow, survey, sampling, ...) is loaded.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

suppressMessages(library(leafblower))

.vcr_here <- tryCatch(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                       error = function(e) ".")
if (is.na(.vcr_here) || !nzchar(.vcr_here)) .vcr_here <- "."
.vcr_common <- file.path(.vcr_here, "..", "common")
.vcr_R <- file.path(.vcr_here, "..", "R")
source(file.path(.vcr_common, "problem_io.R"))       # load_problem_spec
source(file.path(.vcr_common, "instance_family.R"))  # install_gen_resolver
source(file.path(.vcr_here, "replicate_weights.R"))  # rm_make_replicate_weights
install_gen_resolver()

#' Deterministic fingerprint of a replicate-weight matrix -- the fair-
#' comparison audit key: two calls with the SAME repwt matrix (same problem +
#' seed) must produce the SAME fingerprint, regardless of solver_id.
.vcr_repwt_fingerprint <- function(repwt) {
  if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(repwt)
  } else {
    paste(dim(repwt), format(sum(repwt), digits = 17), format(sum(repwt^2), digits = 17))
  }
}

# leafblower_{raking,greg,logit,oris}_r -- public harvest() entry point,
# design_weights = the replicate's weight vector (COLD: no start_weights=,
# so the dual starts fresh at the reference measure every replicate).
.vcr_lbw_solve <- function(method) {
  function(p) {
    res <- leafblower::harvest(
      data = p$data, target = p$targets, method = method,
      min_weight = p$bounds$min, max_weight = p$bounds$max,
      design_weights = p$design_weights,
      convergence = list(absolute = p$tol), attach_weights = FALSE
    )
    result <- attr(res, "result")
    status_code <- result$status
    status <- unname(.LBW_STATUS_MAP[as.character(status_code)])
    if (is.na(status)) status <- "error"
    iterations <- attr(res, "iterations")
    list(status = status, converged = identical(status, "converged"),
         iterations = if (is.null(iterations)) NA_integer_ else as.integer(iterations))
  }
}

# .LBW_STATUS_MAP (RK_* code -> harmonized enum) is defined in
# leafblower_adapter.R; source it here (read-only reuse) so .vcr_lbw_solve
# above can map harvest()'s attr(res,"result")$status the SAME way
# run_leafblower() does, without duplicating the map.
source(file.path(.vcr_R, "leafblower_adapter.R"))
source(file.path(.vcr_R, "competitors.R"))

.VCR_SOLVERS <- list(
  leafblower_raking_r = .vcr_lbw_solve("raking"),
  leafblower_greg_r    = .vcr_lbw_solve("greg"),
  leafblower_logit_r   = .vcr_lbw_solve("logit"),
  leafblower_oris_r    = .vcr_lbw_solve("oris"),
  survey_calibrate_raking = function(p) .vcr_from_comp(.comp_survey_solve(p, "raking")),
  survey_calibrate_linear = function(p) .vcr_from_comp(.comp_survey_solve(p, "linear")),
  survey_calibrate_logit  = function(p) .vcr_from_comp(.comp_survey_solve(p, "logit")),
  sampling_calib_linear = function(p) .vcr_from_comp(.comp_sampling_solve(p, "linear")),
  sampling_calib_logit  = function(p) .vcr_from_comp(.comp_sampling_solve(p, "logit")),
  anesrake = function(p) .vcr_from_comp(.comp_anesrake_solve(p))
)

# Adapt a competitors.R solve_fn's list(weights, iterations, status,
# error_message) return into the (status, converged, iterations) triple this
# file records -- converged is read from the harmonized status enum, same
# rule .comp_run / run_leafblower apply (status == "converged").
.vcr_from_comp <- function(r) {
  iterations <- r$iterations
  if (is.null(iterations) || length(iterations) == 0 || is.na(iterations)) {
    iterations <- NA_integer_
  } else {
    iterations <- as.integer(iterations)
  }
  list(status = r$status, converged = identical(r$status, "converged"), iterations = iterations)
}

#' Measure COLD per-replicate re-solve cost of `solver_id` on `problem` over
#' B Rao-Wu-Yue rescaled bootstrap replicates.
#'
#' @param problem   a problem object (common/problem_io.R's load_problem_spec())
#' @param solver_id one of names(.VCR_SOLVERS)
#' @param B         integer number of replicates, B >= 1
#' @param seed      integer RNG seed passed to rm_make_replicate_weights();
#'                  the SAME seed with the SAME problem produces the SAME
#'                  replicate-weight matrix for every solver_id (fair
#'                  comparison -- see repwt_fingerprint).
#' @return a data.frame with B rows: solver, problem, rep, wall_time_s,
#'   status, converged, iterations, repwt_fingerprint.
rm_variance_cost_cell <- function(problem, solver_id, B, seed) {
  thunk <- .VCR_SOLVERS[[solver_id]]
  if (is.null(thunk)) {
    stop(sprintf("rm_variance_cost_cell: unknown solver_id '%s' (known: %s)",
                 solver_id, paste(names(.VCR_SOLVERS), collapse = ", ")))
  }

  base_dw <- problem$design_weights
  n <- nrow(problem$data)
  repwt <- rm_make_replicate_weights(base_dw, B, seed)  # n x B, shared across solvers
  fingerprint <- .vcr_repwt_fingerprint(repwt)

  rows <- vector("list", B)
  for (b in seq_len(B)) {
    p_b <- problem
    p_b$design_weights <- repwt[, b]

    t0 <- Sys.time()
    out <- tryCatch(thunk(p_b), error = function(e) {
      list(status = "error", converged = FALSE, iterations = NA_integer_)
    })
    wall_time_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    rows[[b]] <- data.frame(
      solver = solver_id,
      problem = problem$id,
      rep = b,
      wall_time_s = wall_time_s,
      status = out$status,
      converged = isTRUE(out$converged),
      iterations = if (is.null(out$iterations) || is.na(out$iterations)) NA_integer_ else as.integer(out$iterations),
      repwt_fingerprint = fingerprint,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

## ---------------------------------------------------------------------------
## CLI main: only runs when this file is invoked directly via Rscript (not
## when source()'d by the test suite), matching capture_one.R's convention.
## ---------------------------------------------------------------------------
.vcr_is_main <- function() {
  file_args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_args) == 0L) return(FALSE)
  # Only the LAST --file= is the top-level script Rscript was invoked with;
  # source()-ing this file from another script (e.g. the test suite) does
  # NOT add a new --file= entry, so basename-matching the last one against
  # this file's own name distinguishes "run directly" from "sourced".
  basename(sub("^--file=", "", file_args[length(file_args)])) == "variance_cost_run.R"
}

if (.vcr_is_main()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 5L) {
    stop("usage: variance_cost_run.R <spec_path> <solver_id> <B> <seed> <out_shard_path>")
  }
  spec_path <- args[[1]]; solver_id <- args[[2]]
  B <- as.integer(args[[3]]); seed <- as.integer(args[[4]]); out_shard_path <- args[[5]]

  p <- load_problem_spec(spec_path)
  cell <- rm_variance_cost_cell(p, solver_id, B, seed)
  dir.create(dirname(out_shard_path), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(cell, out_shard_path)
}
