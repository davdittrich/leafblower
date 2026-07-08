#!/usr/bin/env Rscript
# benchmarks/study/R/leafblower_adapter.R -- WU-7 (leafblower-2ouc.8).
#
# Uniform-contract R adapters wrapping leafblower's OWN 9 solver methods
# (oris, oris_soft, raking, sinkhorn, greenkhorn, chebyshev, greg, logit,
# newton_kl), timed under TWO build variants (portable / native).
#
# STRICT SEPARATION (user constraint, 2026-07-08): this file calls leafblower
# ONLY through its PUBLIC `leafblower::harvest()` entry point. No leafblower
# source file (src/, r_bridge.cpp, R/, DESCRIPTION, NAMESPACE) is read or
# edited here, and there is no solve-only entry point -- leafblower is timed
# end-to-end (groupby/cell-build inside the timer) exactly like a whole-unit
# competitor (contract.md Section 1.1).
#
# Contract v2 (benchmarks/study/spec/contract.md): every `run_leafblower()`
# call below returns EXACTLY
#   {weights_ref, iterations, status, converged, error_message,
#    wall_time_s, peak_rss_bytes}
# `status` is leafblower's OWN native RK_* status code mapped to the 8-value
# harmonized enum (status_enum.json), with a harness-side bound_violation
# reclassification layered on top. `converged` is ALWAYS independently
# recomputed here from common/metrics.R's margin_stats() against
# problem$tol -- never the solver's self-reported status (contract.md
# Section 2.4).
#
# `build` (portable | native) is NOT part of the frozen 7-key return -- per
# the WU-7 ticket's own "Format" field, leafblower is the only adapter with
# a build axis, so `run_leafblower(problem, method, build)` takes it as an
# explicit third argument; it surfaces in output only via the weights_ref
# filename (runs_schema.json: `build` is a driver-added column, not part of
# an adapter's own return schema).

# Path convention matches benchmarks/study/common/test_problem_io.R's here():
# cwd == repo root for every command in this project (CLAUDE.md "Working
# dir: /home/dd/Gemini/leafblower"), so paths are resolved from getwd()
# rather than commandArgs() -- robust to this file being source()'d from a
# directory other than benchmarks/study/R (e.g. a future WU-8 driver script).
.REPO_ROOT <- normalizePath(getwd())
.STUDY_DIR <- file.path(.REPO_ROOT, "benchmarks", "study")

source(file.path(.STUDY_DIR, "common", "metrics.R"))

WEIGHTS_DIR <- file.path(.REPO_ROOT, "weights")

.STATUS_ENUM <- jsonlite::fromJSON(file.path(.STUDY_DIR, "spec", "status_enum.json"))[["$defs"]][["StatusEnum"]][["enum"]]

.LBW_METHODS <- c("oris", "oris_soft", "raking", "sinkhorn", "greenkhorn",
                   "chebyshev", "greg", "logit", "newton_kl")

# RK_* status code (leafblower.h) -> harmonized enum. This is a DIRECT 1:1
# map for the codes harvest() RETURNS normally (0,1,4,5). Codes 2 (infeasible)
# and 3 (bad_arg) are ONLY ever thrown via stop() (R/harvest.R:724-728) --
# they never reach this map; see .lbw_classify_error() below.
.LBW_STATUS_MAP <- c(`0` = "converged", `1` = "no_conv", `4` = "budget", `5` = "stall")

#' Best-effort in-process high-water-mark RSS, in bytes.
#' contract.md Section 2.7: VmHWM sampled from /proc/self/status (Linux).
.lbw_peak_rss_bytes <- function() {
  lines <- tryCatch(readLines("/proc/self/status"), error = function(e) character(0))
  hwm <- grep("^VmHWM:", lines, value = TRUE)
  if (length(hwm) == 0L) return(NA_integer_)
  kb <- as.numeric(sub("^VmHWM:\\s*([0-9]+)\\s*kB.*$", "\\1", hwm[1]))
  as.integer(kb * 1024)
}

#' Classify an uncaught harvest() stop() into the harmonized status enum.
#' harvest()'s ONLY two custom-prefixed stop()s (R/harvest.R:724-728) are:
#'   status==3 -> "leafblower: invalid arguments -- <msg>"
#'   status==2 -> "leafblower: <msg>" / "leafblower: infeasible problem"
#' Every other stop() is an R-level pre-C-call validation error (no
#' "leafblower:" prefix, e.g. bad max_iterations/min_weight/max_weight) --
#' classified bad_arg, matching status 3's semantics (invalid user-supplied
#' argument caught before the solver ever ran).
.lbw_classify_error <- function(msg) {
  if (grepl("invalid arguments", msg, fixed = TRUE)) return("bad_arg")
  if (grepl("^leafblower:", msg)) return("infeasible")
  "bad_arg"
}

.lbw_weights_path <- function(solver_id, problem_id, thread, build) {
  dir.create(WEIGHTS_DIR, recursive = TRUE, showWarnings = FALSE)
  file.path(WEIGHTS_DIR, sprintf("%s__%s__t%d__%s.parquet", solver_id, problem_id, thread, build))
}

.lbw_write_weights <- function(weights, solver_id, problem_id, thread, build) {
  path <- .lbw_weights_path(solver_id, problem_id, thread, build)
  arrow::write_parquet(data.frame(weight = as.numeric(weights)), path)
  sub(paste0("^", .REPO_ROOT, "/"), "", path)
}

#' contract.md Section 2.1: hard-failure runs still write a length-n
#' all-NaN sentinel vector so weights_ref is never dangling.
.lbw_nan_sentinel <- function(n, solver_id, problem_id, thread, build) {
  .lbw_write_weights(rep(NA_real_, n), solver_id, problem_id, thread, build)
}

.lbw_bound_violation <- function(weights, problem, atol = 1e-9) {
  lo <- problem$bounds$min
  hi <- problem$bounds$max
  any(weights < lo - atol) || any(weights > hi + atol)
}

#' contract.md Section 2.4: converged is ALWAYS harness-recomputed here via
#' common/metrics.R's margin_stats() uniform margin-L-infinity check against
#' problem$tol -- never the solver's own self-reported status.
.lbw_recompute_converged <- function(weights, problem) {
  groups <- stats::setNames(lapply(problem$margins, function(m) problem$data[[m]]), problem$margins)
  ms <- margin_stats(weights, groups, problem$targets)
  list(converged = isTRUE(ms$margin_linf <= problem$tol), margin_linf = ms$margin_linf)
}

.lbw_thread <- function() as.integer(Sys.getenv("OMP_NUM_THREADS", "1"))

#' Run one leafblower method on `problem` through the public harvest() API
#' and return the frozen 7-key contract v2 result list.
#'
#' @param problem  a problem object as returned by common/problem_io.R's loader
#' @param method   one of .LBW_METHODS
#' @param build    "portable" | "native" -- tags the weights_ref filename only
#' @param thread   thread-count tag for the weights_ref filename; defaults to
#'                 OMP_NUM_THREADS (single-thread-BLAS convention, CLAUDE.md)
run_leafblower <- function(problem, method, build, thread = NULL) {
  stopifnot(method %in% .LBW_METHODS)
  if (is.null(thread)) thread <- .lbw_thread()
  solver_id <- paste0("leafblower_", method)
  n  <- nrow(problem$data)
  # contract.md Section 2.6: high-res in-process timer. Sys.time() (gettimeofday,
  # microsecond resolution on Linux) rather than proc.time()[["elapsed"]] (clock
  # ticks, ~10ms resolution) -- the latter reads exactly 0 for a solve this fast
  # on a 4-row toy problem, which is not a meaningful measurement.
  t0 <- Sys.time()

  warnings_seen <- character(0)
  res <- withCallingHandlers(
    tryCatch(
      leafblower::harvest(
        data           = problem$data,
        target         = problem$targets,
        method         = method,
        min_weight     = problem$bounds$min,
        max_weight     = problem$bounds$max,
        design_weights = problem$design_weights,
        convergence    = list(absolute = problem$tol),
        attach_weights = FALSE
      ),
      error = function(e) e
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  # contract.md Section 2.6: single timer stop at weights-out (harvest returned,
  # success or error). The recompute/bound-check/parquet write below are post-run
  # harness diagnostics (Section 2.4) and MUST NOT be inside the timer.
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (inherits(res, "error")) {
    status <- .lbw_classify_error(conditionMessage(res))
    ref    <- .lbw_nan_sentinel(n, solver_id, problem$id, thread, build)
    return(list(
      weights_ref = ref, iterations = NA_integer_, status = status,
      converged = FALSE, error_message = conditionMessage(res),
      wall_time_s = as.numeric(wall), peak_rss_bytes = .lbw_peak_rss_bytes()
    ))
  }

  result      <- attr(res, "result")
  status_code <- result$status
  status      <- unname(.LBW_STATUS_MAP[as.character(status_code)])
  if (is.na(status)) status <- "error"  # unmapped/unexpected returned code -- defensive only

  weights <- as.numeric(res)
  if (.lbw_bound_violation(weights, problem)) status <- "bound_violation"

  cv  <- .lbw_recompute_converged(weights, problem)
  ref  <- .lbw_write_weights(weights, solver_id, problem$id, thread, build)
  iterations <- attr(res, "iterations")

  list(
    weights_ref = ref,
    iterations  = if (is.null(iterations)) NA_integer_ else as.integer(iterations),
    status      = status,
    converged   = cv$converged,
    error_message = if (length(warnings_seen) == 0L) NULL else paste(warnings_seen, collapse = "; "),
    wall_time_s = as.numeric(wall),
    peak_rss_bytes = .lbw_peak_rss_bytes()
  )
}

LEAFBLOWER_R_ADAPTERS <- stats::setNames(
  lapply(.LBW_METHODS, function(m) {
    force(m)
    function(problem, build) run_leafblower(problem, m, build)
  }),
  paste0("leafblower_", .LBW_METHODS, "_r")
)
