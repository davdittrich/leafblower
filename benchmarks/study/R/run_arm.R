#!/usr/bin/env Rscript
# benchmarks/study/R/run_arm.R -- WU-8 (leafblower-2ouc.9).
#
# Registry-driven run-matrix driver for the R arm (competitors.R's 19
# ranked competitors, and leafblower_adapter.R's 9 methods via public
# leafblower::harvest()). Mirrors python/run_arm.py 1:1 -- see that file's
# module docstring for the full measurement-protocol rationale (DESIGN.md
# Sec4/Sec5): fresh subprocess per cell, registry-gated applicability,
# single wall_time_s per contract v2, data-loaded RSS baseline, hardware-
# isolation capture (warn-not-fail), randomized cell order, ott-jax
# excluded from ranked timing, pre-run frozen-tag hook, build is a
# pass-through tag only (no CMakeLists/Makevars edits here).
#
# STRICT SEPARATION (user constraint, 2026-07-08): this file, and every
# file it sources, reaches leafblower ONLY through leafblower_adapter.R's
# public `leafblower::harvest()` wrapper. No leafblower source file is
# read or edited by this WU.
#
# Path convention matches leafblower_adapter.R / run_matrix.R: cwd == repo
# root for every command in this project.
#
# Usage (repo root):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#     Rscript benchmarks/study/R/run_arm.R --smoke
#   Rscript benchmarks/study/R/run_arm.R --sync-registry

.RUN_ARM_REPO_ROOT <- normalizePath(getwd())
source(file.path(.RUN_ARM_REPO_ROOT, "benchmarks", "study", "common", "run_matrix.R"))

RUN_ARM_PATH <- file.path(STUDY_DIR, "R", "run_arm.R")

CONTRACT_KEYS <- c("weights_ref", "iterations", "status", "converged",
                    "error_message", "wall_time_s", "peak_rss_bytes")
STATUS_ENUM <- c("converged", "no_conv", "infeasible", "bound_violation",
                  "bad_arg", "budget", "stall", "error", "dnf")
RUNS_ROW_KEYS <- c("solver", "problem", "thread", "build", "rep",
                    "weights_ref", "iterations", "status", "converged",
                    "error_message", "wall_time_s", "peak_rss_bytes", "trajectory_ref")

## ----------------------------------------------------------------------
## Small shared helpers
## ----------------------------------------------------------------------

.set_thread_env <- function(thread) {
  # MUST run before leafblower_adapter.R/competitors.R (and thus
  # leafblower::/arrow's BLAS-linked deps) are sourced in THIS process --
  # CLAUDE.md single-thread-BLAS determinism rule, extended to the
  # per-cell thread sweep.
  Sys.setenv(OMP_NUM_THREADS = as.character(thread),
             OPENBLAS_NUM_THREADS = as.character(thread),
             MKL_NUM_THREADS = as.character(thread))
}

.peak_rss_bytes <- function() {
  lines <- tryCatch(readLines("/proc/self/status"), error = function(e) character(0))
  hwm <- grep("^VmHWM:", lines, value = TRUE)
  if (length(hwm) == 0L) return(NA_real_)
  as.numeric(sub("^VmHWM:\\s*([0-9]+)\\s*kB.*$", "\\1", hwm[1])) * 1024
}

.get <- function(r, k, default = NULL) {
  v <- r[[k]]
  if (is.null(v)) default else v
}

`%+%` <- function(a, b) paste0(a, b)

## ----------------------------------------------------------------------
## Worker mode: ONE fresh subprocess, ONE (solver,problem,thread,build) cell.
## ----------------------------------------------------------------------

#' Lazily source the adapter file that owns `solver_id` (deferred exactly
#' like leafblower_adapter.py/competitors.py's per-call imports in
#' run_arm.py's _resolve_adapter -- NOT eagerly at file top, so worker mode
#' can set OMP_/OPENBLAS_/MKL_NUM_THREADS before any BLAS-linked package
#' (arrow, leafblower::) loads).
.resolve_adapter <- function(solver_id) {
  if (grepl("^leafblower_.*_r$", solver_id)) {
    if (!exists("LEAFBLOWER_R_ADAPTERS")) source(file.path(STUDY_DIR, "R", "leafblower_adapter.R"))
    fn <- LEAFBLOWER_R_ADAPTERS[[solver_id]]
    if (is.null(fn)) stop(sprintf("unknown leafblower R solver id %s", solver_id))
    return(function(problem, build) fn(problem, build))
  }
  fn_name <- paste0("run_", solver_id)
  if (!exists(fn_name, mode = "function")) source(file.path(STUDY_DIR, "R", "competitors.R"))
  if (!exists(fn_name, mode = "function")) {
    stop(sprintf("no run_%s() adapter found in competitors.R -- registry entry has no implementation",
                  solver_id))
  }
  fn <- get(fn_name, mode = "function")
  function(problem, build) fn(problem)  # competitor adapters hardcode build="na" (frozen)
}

.subsample <- function(problem, n_cap, seed) {
  n <- nrow(problem$data)
  if (is.null(n_cap) || n <= n_cap) return(problem)
  set.seed(seed)
  idx <- sample.int(n, n_cap)
  out <- problem
  out$data <- problem$data[idx, , drop = FALSE]
  if (!is.null(problem$design_weights)) out$design_weights <- problem$design_weights[idx]
  out
}

run_worker <- function(cell_path, result_path) {
  cell <- jsonlite::fromJSON(cell_path, simplifyVector = TRUE)
  .set_thread_env(cell$thread)  # MUST precede any source() below (BLAS-linked deps)

  source(file.path(STUDY_DIR, "common", "problem_io.R"))
  # Instance-family specs carry a gen:instance_family?... data_ref that
  # load_problem_spec() cannot resolve alone (raises NotImplemented). Source
  # instance_family.R AFTER problem_io.R and install_gen_resolver() -- it
  # monkey-patches the just-sourced .pio_resolve_data_ref, is idempotent, and
  # intercepts ONLY gen:instance_family refs (a no-op for the 4 static specs).
  source(file.path(STUDY_DIR, "common", "instance_family.R"))
  install_gen_resolver()
  fn <- .resolve_adapter(cell$solver)  # lazily sources the owning adapter file
  problem <- load_problem_spec(rm_resolve_spec_path(cell$problem))
  problem <- .subsample(problem, cell$n_cap, if (is.null(cell$seed)) 0L else cell$seed)

  # Universal structural-infeasibility short-circuit: if a target category has 0
  # sample observations, the cell is infeasible for EVERY solver -> record one
  # "infeasible" row, do NOT run the solver (no misleading crash/error/fake-conv).
  infeas <- .structural_infeasible_cats(problem)
  if (length(infeas)) {
    jsonlite::write_json(list(rows = list(.infeasible_row(cell, infeas))), result_path,
                         auto_unbox = TRUE, null = "null")
    return(invisible())
  }

  warmups <- as.integer(cell$warmups)
  min_reps <- as.integer(cell$reps)
  min_total_duration <- as.numeric(cell$min_total_duration)
  max_reps <- as.integer(cell$max_reps)

  timed <- list()
  calls <- 0L
  cumulative <- 0.0
  repeat {
    res <- fn(problem, cell$build)
    calls <- calls + 1L
    if (calls <= warmups) next
    timed[[length(timed) + 1L]] <- res
    cumulative <- cumulative + as.numeric(res$wall_time_s)
    if (length(timed) >= min_reps && cumulative >= min_total_duration) break
    if (length(timed) >= max_reps) break
  }

  weights_refs <- unique(vapply(timed, function(r) r$weights_ref, character(1)))
  if (length(weights_refs) > 1L) {
    message(sprintf("WARN: cell %s/%s/t%s/%s: weights_ref differs across reps (%s) -- contract.md "
                     %+% "Sec2.1 expects a single shared file per cell; recording as observed.",
                     cell$solver, cell$problem, cell$thread, cell$build,
                     paste(weights_refs, collapse = ", ")))
  }

  peaks <- vapply(timed, function(r) { v <- .get(r, "peak_rss_bytes", NA_real_); as.numeric(v) }, double(1))
  cell_peak_rss <- if (all(is.na(peaks))) .peak_rss_bytes() else max(peaks, na.rm = TRUE)

  rows <- vector("list", length(timed))
  for (i in seq_along(timed)) {
    r <- timed[[i]]
    if (!identical(sort(names(r)), sort(CONTRACT_KEYS))) {
      stop(sprintf("%s: adapter key drift %s", cell$solver, paste(sort(names(r)), collapse = ",")))
    }
    if (!(r$status %in% STATUS_ENUM)) stop(sprintf("%s: status %s not harmonized", cell$solver, r$status))
    rows[[i]] <- list(
      solver = cell$solver, problem = cell$problem, thread = as.integer(cell$thread),
      build = cell$build, rep = i - 1L,
      weights_ref = r$weights_ref,
      iterations = if (is.na(.get(r, "iterations", NA))) NULL else as.integer(r$iterations),
      status = r$status, converged = isTRUE(r$converged),
      error_message = .get(r, "error_message", NULL),
      wall_time_s = as.numeric(r$wall_time_s), peak_rss_bytes = cell_peak_rss,
      # Per-iteration margin-error trajectory (RQ3): no currently-shipped
      # adapter exposes it beyond the frozen 7-key contract (asserted
      # above) -- honest always-null hook, see run_arm.py's mirror comment.
      trajectory_ref = NULL
    )
  }

  jsonlite::write_json(list(rows = rows), result_path, auto_unbox = TRUE, null = "null")
}

run_baseline_worker <- function(solver_id, result_path) {
  .set_thread_env(1L)
  source(file.path(STUDY_DIR, "common", "problem_io.R"))
  t0 <- Sys.time()
  if (startsWith(solver_id, "leafblower_")) {
    source(file.path(STUDY_DIR, "R", "leafblower_adapter.R"))
  } else {
    source(file.path(STUDY_DIR, "R", "competitors.R"))
  }
  load_problem_spec(file.path(STUDY_DIR, "spec", "toy_inline.json"))
  load_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  jsonlite::write_json(
    list(solver = solver_id, load_time_s = load_s, peak_rss_bytes = .peak_rss_bytes()),
    result_path, auto_unbox = TRUE, null = "null"
  )
}

## ----------------------------------------------------------------------
## runs.parquet row assembly
## ----------------------------------------------------------------------

.rows_to_df <- function(rows) {
  n <- length(rows)
  solver <- character(n); problem <- character(n); thread <- integer(n)
  build <- character(n); rep_ <- integer(n); weights_ref <- character(n)
  iterations <- integer(n); status <- character(n); converged <- logical(n)
  error_message <- character(n); wall_time_s <- double(n)
  peak_rss_bytes <- double(n); trajectory_ref <- character(n)
  for (i in seq_len(n)) {
    r <- rows[[i]]
    solver[i] <- .get(r, "solver", NA_character_)
    problem[i] <- .get(r, "problem", NA_character_)
    thread[i] <- as.integer(.get(r, "thread", NA))
    build[i] <- .get(r, "build", NA_character_)
    rep_[i] <- as.integer(.get(r, "rep", NA))
    weights_ref[i] <- .get(r, "weights_ref", NA_character_)
    it <- .get(r, "iterations", NA); iterations[i] <- if (is.null(it) || is.na(it)) NA_integer_ else as.integer(it)
    status[i] <- .get(r, "status", NA_character_)
    converged[i] <- isTRUE(.get(r, "converged", NA))
    em <- .get(r, "error_message", NA_character_); error_message[i] <- if (is.null(em)) NA_character_ else as.character(em)
    wall_time_s[i] <- as.numeric(.get(r, "wall_time_s", NA))
    peak_rss_bytes[i] <- as.numeric(.get(r, "peak_rss_bytes", NA))
    tr <- .get(r, "trajectory_ref", NA_character_); trajectory_ref[i] <- if (is.null(tr)) NA_character_ else as.character(tr)
  }
  data.frame(solver = solver, problem = problem, thread = thread, build = build, rep = rep_,
             weights_ref = weights_ref, iterations = iterations, status = status,
             converged = converged, error_message = error_message, wall_time_s = wall_time_s,
             peak_rss_bytes = peak_rss_bytes, trajectory_ref = trajectory_ref,
             stringsAsFactors = FALSE)
}

## ----------------------------------------------------------------------
## Orchestrator mode
## ----------------------------------------------------------------------

#' Fresh-subprocess spawn with a wall-clock cap (a hung competitor must not
#' stall the whole matrix). R's system2() has no native timeout argument;
#' coreutils `timeout` is used when available (nzchar(Sys.which("timeout")))
#' and skipped (uncapped, with a one-time warning) otherwise.
#'
#' `cpus` (a taskset -c CPU list string, e.g. "3" or "4,5,6,7"): when set AND
#' `taskset` is available, the worker is PINNED to exactly those CPUs. This
#' CONFINES a solver whose backend ignores the single-thread env vars (JAX/XLA for
#' ott_jax_sinkhorn, Rust rayon for cvxpy CLARABEL) to its assigned cores -- a
#' 2800%-CPU runaway that would otherwise starve every co-scheduled worker and
#' POISON the timing is contained to 100% x |cpus| of its own cores. A thread=k
#' cell is given k physical cores so its k BLAS threads are not strangled on one.
#' taskset exec()s the child, so the timeout/exit-124 semantics are unchanged.
.spawn <- function(args, timeout_s = 120, cpus = NULL) {
  old <- getwd()
  on.exit(setwd(old))
  setwd(REPO_ROOT)
  have_timeout <- nzchar(Sys.which("timeout"))
  base <- if (have_timeout) {
    c("timeout", sprintf("%ss", timeout_s), "Rscript", shQuote(RUN_ARM_PATH), args)
  } else {
    c("Rscript", shQuote(RUN_ARM_PATH), args)
  }
  if (!is.null(cpus) && nzchar(cpus) && nzchar(Sys.which("taskset"))) {
    cmd <- "taskset"
    full_args <- c("-c", cpus, base)
  } else {
    cmd <- base[1]
    full_args <- base[-1]
  }
  out <- suppressWarnings(system2(cmd, full_args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  status <- if (is.null(status)) 0L else status
  if (have_timeout && identical(status, 124L)) {
    return(list(status = status, output = sprintf("cell timed out after %ss", timeout_s)))
  }
  list(status = status, output = paste(out, collapse = "\n"))
}

#' Right-censored DNF row for a cell killed at the wall-clock budget (coreutils
#' `timeout` exit 124). Recorded, NOT dropped, so a timeout is distinguishable
#' from a cell that never ran. Weights are UNDEFINED for an unfinished cell -> a
#' length-1 all-NaN sentinel (non-dangling; scoring skips weights for
#' status=='dnf'). Mirrors run_arm.py::_dnf_row.
.dnf_row <- function(cell, budget_s) {
  ref_rel <- file.path("benchmarks", "study", "results", "weights",
                       sprintf("%s__%s__t%d__%s.parquet",
                               cell$solver, cell$problem, as.integer(cell$thread), cell$build))
  dir.create(dirname(ref_rel), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(data.frame(weight = NA_real_), ref_rel)
  list(
    solver = cell$solver, problem = cell$problem, thread = as.integer(cell$thread),
    build = cell$build, rep = 0L, weights_ref = ref_rel,
    iterations = NULL, status = "dnf", converged = FALSE,
    error_message = sprintf("exceeded --cell-timeout=%ss (right-censored DNF)", budget_s),
    wall_time_s = as.numeric(budget_s), peak_rss_bytes = NA_real_, trajectory_ref = NULL
  )
}

#' Structural infeasibility: a margin category with target>0 but ZERO sample
#' observations (no row carries it) cannot be reached by ANY reweighting.
#' Extreme-skew/high-cardinality small-n instances realize only a subset of the
#' declared categories. Detected on the MATERIALIZED (post-subsample) problem the
#' solver would see, so the universal short-circuit below applies uniformly to
#' every solver -- competitors AND leafblower -- for a clean apples-to-apples
#' "infeasible" classification (the head-to-head runs on feasible cells). Mirrors
#' run_arm.py::structural_infeasible_cats.
.structural_infeasible_cats <- function(problem) {
  bad <- character(0)
  for (m in problem$margins) {
    tg <- problem$targets[[m]]
    present <- unique(as.character(problem$data[[m]]))
    miss <- names(tg)[tg > 0 & !(names(tg) %in% present)]
    if (length(miss)) bad <- c(bad, paste0(m, ".", miss))
  }
  bad
}

#' Row synthesized when a cell is structurally infeasible: no solver runs; weights
#' are undefined -> length-1 all-NaN sentinel (scoring skips), status="infeasible".
.infeasible_row <- function(cell, cats) {
  ref_rel <- file.path("benchmarks", "study", "results", "weights",
                       sprintf("%s__%s__t%d__%s.parquet",
                               cell$solver, cell$problem, as.integer(cell$thread), cell$build))
  dir.create(dirname(ref_rel), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(data.frame(weight = NA_real_), ref_rel)
  list(
    solver = cell$solver, problem = cell$problem, thread = as.integer(cell$thread),
    build = cell$build, rep = 0L, weights_ref = ref_rel,
    iterations = NULL, status = "infeasible", converged = FALSE,
    error_message = sprintf("structurally infeasible: target>0 but 0 observations for %s%s",
                            paste(head(cats, 3), collapse = ", "),
                            if (length(cats) > 3) sprintf(" (+%d more)", length(cats) - 3) else ""),
    wall_time_s = 0.0, peak_rss_bytes = NA_real_, trajectory_ref = NULL
  )
}

#' Row synthesized when a worker DIES without writing a result (non-timeout, e.g.
#' OOM-kill / allocation failure / hard crash): recorded (not dropped to
#' cell_failures) so a crashed cell is distinguishable from one that never ran --
#' the WU-9 no-selective-reporting guarantee (extends the .dnf_row treatment from
#' timeouts to crashes). status="error"; weights undefined -> NaN sentinel; timing
#' unknown -> NA. exit_status 137 == SIGKILL (typically OOM). Mirrors run_arm.py.
.crash_row <- function(cell, exit_status, reason) {
  ref_rel <- file.path("benchmarks", "study", "results", "weights",
                       sprintf("%s__%s__t%d__%s.parquet",
                               cell$solver, cell$problem, as.integer(cell$thread), cell$build))
  dir.create(dirname(ref_rel), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(data.frame(weight = NA_real_), ref_rel)
  oom <- if (identical(as.integer(exit_status), 137L)) " (exit 137 = SIGKILL, likely OOM)" else ""
  list(
    solver = cell$solver, problem = cell$problem, thread = as.integer(cell$thread),
    build = cell$build, rep = 0L, weights_ref = ref_rel,
    iterations = NULL, status = "error", converged = FALSE,
    error_message = sprintf("worker died without result: exit=%s%s; captured: %s",
                            as.character(exit_status), oom, substr(gsub("\\s+", " ", reason), 1, 160)),
    wall_time_s = NA_real_, peak_rss_bytes = NA_real_, trajectory_ref = NULL
  )
}

#' Solver ids resolving to a runnable R adapter this arm: the run_<id> competitor
#' functions (competitors.R) + the 9 leafblower_*_r methods. A registry entry
#' without one (e.g. `cvxr_reference`, the convex anchor produced separately via
#' produce_ref_rows.R -- NOT a matrix competitor) is gated OUT, else
#' .resolve_adapter stop()s on it.
.r_available_adapters <- function(registry) {
  source(file.path(STUDY_DIR, "R", "competitors.R"))
  if (!exists("LEAFBLOWER_R_ADAPTERS")) source(file.path(STUDY_DIR, "R", "leafblower_adapter.R"))
  ids <- names(registry$solvers)
  keep <- Filter(function(id) {
    if (grepl("^leafblower_", id)) return(id %in% names(LEAFBLOWER_R_ADAPTERS))
    exists(paste0("run_", id), mode = "function")
  }, ids)
  unlist(keep)
}

orchestrate <- function(opts) {
  registry <- rm_load_registry()
  # Full matrix problem set: 4 static specs UNION the WU-3 instance family (32
  # instances). The instance family is the only arena declaring ot/newton_kl
  # objective_families -- what gives the OT / Newton-KL solvers their applicable
  # problems. Also materializes the spec/instance_family/<id>.json the worker loads.
  # regenerate=FALSE: use the COMMITTED (frozen) instance specs -- regenerating
  # would overwrite them and dirty the frozen runnable tree (WU-9). The committed
  # spec/instance_family/*.json ARE the frozen truth; regen only on a missing dir.
  problem_specs <- rm_load_all_problem_specs(regenerate = FALSE)
  hyperparams <- rm_load_hyperparams()
  available <- .r_available_adapters(registry)

  if (isTRUE(opts$sync_registry)) {
    registry <- rm_compute_applicable_problems(registry, problem_specs,
                                               available = available, arm = "R")
    jsonlite::write_json(registry, file.path(STUDY_DIR, "registry.json"),
                          auto_unbox = TRUE, pretty = TRUE, null = "null")
    # --sync-registry is SYNC-ONLY: refresh applicable_problems and exit, never
    # launching a run. A bounded --smoke may be combined (sync then smoke); a
    # bare --sync-registry must NOT fall through to the full production matrix
    # (WU-11, post-freeze).
    if (!isTRUE(opts$smoke)) {
      return(list(synced = TRUE, n_solvers = length(registry$solvers),
                  registry_path = file.path(STUDY_DIR, "registry.json")))
    }
  }

  tag_status <- rm_assert_frozen_tag()
  # WU-9: a scored/timed run must execute against the frozen runnable tree. When
  # --assert-runnable-tag is passed (WU-11 scored launcher), hard-stop on a
  # dirty/drifted benchmarks/study tree. The rehearsal omits it (unfrozen by design).
  runnable_status <- NULL
  if (!is.null(opts$assert_runnable_tag)) {
    runnable_status <- rm_assert_runnable_frozen(opts$assert_runnable_tag)
    message(sprintf("WU-9 runnable-tree gate OK: tree matches signed tag '%s'", opts$assert_runnable_tag))
  }
  env_info <- rm_capture_environment(pin_core = opts$pin_core)

  # --smoke bounds the instance family to n <= opts$n: a huge-n synthetic
  # instance would burn minutes in the pure-R generator only to be subsampled
  # to opts$n rows. Static specs (rm_instance_family_n is NULL) always kept.
  matrix_specs <- problem_specs
  if (isTRUE(opts$smoke)) {
    keep <- Filter(function(pid) {
      nn <- rm_instance_family_n(problem_specs[[pid]])
      is.null(nn) || nn <= opts$n
    }, names(problem_specs))
    matrix_specs <- problem_specs[keep]
  }

  threads <- as.integer(strsplit(opts$threads, ",")[[1]])
  cells <- rm_build_matrix(registry, matrix_specs, arm = "R", threads = threads, rng_seed = opts$seed, available = available)
  if (!is.null(opts$solvers)) cells <- Filter(function(c) grepl(opts$solvers, c$solver), cells)
  if (!is.null(opts$problems)) cells <- Filter(function(c) grepl(opts$problems, c$problem), cells)
  if (!is.null(opts$max_cells)) cells <- cells[seq_len(min(length(cells), opts$max_cells))]

  out_dir <- opts$out
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  work_dir <- file.path(out_dir, "_cells")
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  warmups <- if (isTRUE(opts$smoke)) 0L else opts$warmups
  reps <- opts$reps
  min_total_duration <- if (isTRUE(opts$smoke)) 0.0 else opts$min_total_duration
  n_cap <- if (isTRUE(opts$smoke)) opts$n else NULL

  # Parallel cell execution (--jobs N; memory-gated EXTERNALLY by run_all.sh:
  # nonheavy jobs<=run_config default_concurrency, heavy jobs=1 so only one
  # >=0.5*RAM worker is ever alive). Each cell is a fresh --worker subprocess
  # under `timeout`; DNF/crash rows are SYNTHESIZED and PERSISTED to
  # result_<key>.json so a killed run RESUMES (skips completed cells, incl.
  # hour-long right-censored DNFs). Files keyed by cell IDENTITY
  # (solver__problem__t<thread>__build), NOT loop index, so resume is stable
  # across a re-run whose cell order/count differs.
  jobs <- if (is.null(opts$jobs)) 1L else max(1L, as.integer(opts$jobs))
  # CPU-affinity SLOT POOL (fixes the idx%%n_cores collision + HT contention).
  # `pool` = pinnable logical CPUs, ONE per physical core (the launcher passes
  # --pin-cpus so two workers never land on HT siblings of one physical core;
  # fallback = 0..detectCores-1 with the documented NA return guarded). Each slot
  # owns `width` CPUs (width = widest thread count this run) so a thread=4 STRETCH
  # cell gets 4 real cores, not 4 threads strangled on 1. Cells are partitioned
  # round-robin into `jobs` slots; slot s runs its cells SERIALLY pinned to its
  # fixed CPUs => two concurrently-timed workers NEVER share a CPU (collision-free
  # by construction; independent of mclapply's internal chunking).
  pool <- if (!is.null(opts$pin_cpus) && nzchar(opts$pin_cpus)) {
    as.integer(strsplit(opts$pin_cpus, ",", fixed = TRUE)[[1]])
  } else {
    nc <- parallel::detectCores(); nc <- if (is.na(nc)) 1L else max(1L, as.integer(nc))
    seq.int(0L, nc - 1L)
  }
  width <- max(1L, max(threads))
  n_slots <- max(1L, length(pool) %/% width)
  jobs <- min(jobs, n_slots)
  slot_cpus <- function(s) {  # s in 0..jobs-1 -> "c1,c2,..." of `width` CPUs
    lo <- s * width + 1L
    paste(pool[lo:(lo + width - 1L)], collapse = ",")
  }
  cell_key <- function(cell) sprintf("%s__%s__t%d__%s", cell$solver, cell$problem,
                                     as.integer(cell$thread), cell$build)
  run_one <- function(cell, cpus) {
    key <- cell_key(cell)
    cell_path <- file.path(work_dir, sprintf("cell_%s.json", key))
    result_path <- file.path(work_dir, sprintf("result_%s.json", key))
    if (file.exists(result_path)) {  # resume: reuse a prior (incl. DNF/crash) result
      rr <- tryCatch(jsonlite::fromJSON(result_path, simplifyVector = FALSE),
                     error = function(e) NULL)
      if (!is.null(rr) && !is.null(rr$rows)) return(rr$rows)
    }
    cell_full <- c(cell, list(warmups = warmups, reps = reps, max_reps = opts$max_reps,
                              min_total_duration = min_total_duration, n_cap = n_cap,
                              seed = opts$seed))
    jsonlite::write_json(cell_full, cell_path, auto_unbox = TRUE, null = "null")
    r <- .spawn(c("--worker", "--cell", shQuote(cell_path), "--result-out", shQuote(result_path)),
                timeout_s = opts$cell_timeout, cpus = cpus)
    if (!identical(r$status, 0L) || !file.exists(result_path)) {
      row <- if (identical(r$status, 124L)) {
        # Right-censored DNF (coreutils `timeout` exit 124): a real row, never dropped.
        message(sprintf("DNF: cell %s exceeded %ss budget (right-censored)",
                        key, opts$cell_timeout))
        .dnf_row(cell, opts$cell_timeout)
      } else {
        # Worker died without a result (OOM-kill / hard crash): a real error row,
        # never dropped (WU-9 no-selective-reporting).
        message(sprintf("CRASH: cell %s died (exit=%s) -> recorded as error row",
                        key, r$status))
        .crash_row(cell, r$status, r$output)
      }
      # Persist the synthesized row so a resume skips this cell (esp. 1h DNFs).
      jsonlite::write_json(list(rows = list(row)), result_path, auto_unbox = TRUE, null = "null")
      return(list(row))
    }
    jsonlite::fromJSON(result_path, simplifyVector = FALSE)$rows
  }
  # Slot s owns cells {s+1, s+1+jobs, ...} (round-robin) and runs them serially,
  # pinned to slot_cpus(s). One fork per slot => at most `jobs` live workers on
  # `jobs` disjoint CPU sets.
  run_slot <- function(s) {
    cpus <- slot_cpus(s)
    idxs <- which((seq_along(cells) - 1L) %% jobs == s)
    do.call(c, lapply(idxs, function(i) run_one(cells[[i]], cpus)))
  }
  res_list <- if (jobs > 1L) {
    parallel::mclapply(0:(jobs - 1L), run_slot, mc.cores = jobs, mc.preschedule = FALSE)
  } else {
    list(run_slot(0L))
  }
  all_rows <- list()
  for (rl in res_list) {
    if (inherits(rl, "try-error") || is.null(rl)) next  # fork died: cell re-runs on resume
    all_rows <- c(all_rows, rl)
  }
  # cell_failures == crash (status='error') rows; dnf is a separate right-censor.
  cell_failures <- Filter(function(r) identical(r$status, "error"), all_rows)

  baseline_reps <- list()
  for (cell in cells) {
    pkg <- rm_solver_package(cell$solver, hyperparams)
    if (is.null(baseline_reps[[pkg]])) baseline_reps[[pkg]] <- cell$solver
  }
  baselines <- list()
  for (pkg in names(baseline_reps)) {
    solver_id <- baseline_reps[[pkg]]
    result_path <- file.path(work_dir, sprintf("baseline_%s.json", pkg))
    r <- .spawn(c("--baseline-rss", "--solver", shQuote(solver_id), "--result-out", shQuote(result_path)),
                timeout_s = opts$cell_timeout)
    if (identical(r$status, 0L) && file.exists(result_path)) {
      b <- jsonlite::fromJSON(result_path, simplifyVector = TRUE)
      b$package <- pkg
      baselines[[length(baselines) + 1L]] <- b
    } else {
      message(sprintf("WARN: baseline for package %s (via %s) failed: %s", pkg, solver_id, r$output))
    }
  }

  runs_df <- .rows_to_df(all_rows)
  if (nrow(runs_df) > 0L) {
    runs_df$build <- factor(runs_df$build)
    runs_df$status <- factor(runs_df$status, levels = sort(STATUS_ENUM))
  }
  runs_path <- file.path(out_dir, "runs.parquet")
  arrow::write_parquet(runs_df, runs_path)

  jsonlite::write_json(
    c(env_info, list(frozen_tag = tag_status, baselines = baselines,
                      n_cells = length(cells), n_cell_failures = length(cell_failures))),
    file.path(out_dir, "environment.json"), auto_unbox = TRUE, pretty = TRUE, null = "null"
  )

  list(n_cells = length(cells), n_rows = nrow(runs_df), n_failures = length(cell_failures),
       runs_path = runs_path, cell_failures = cell_failures)
}

## ----------------------------------------------------------------------
## CLI parsing
## ----------------------------------------------------------------------

.parse_args <- function(argv) {
  opts <- list(worker = FALSE, baseline_rss = FALSE, cell = NULL, solver = NULL, result_out = NULL,
               smoke = FALSE, n = 5000L, reps = 2L, warmups = 2L, max_reps = 200L,
               min_total_duration = 0.5, threads = "1,4", pin_core = NULL, pin_cpus = NULL, seed = NULL,
               out = file.path(STUDY_DIR, "results"), sync_registry = FALSE, solvers = NULL, problems = NULL, max_cells = NULL,
               cell_timeout = 120, assert_runnable_tag = NULL, jobs = 1L)
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[i]
    val_next <- function() { i <<- i + 1L; argv[i] }
    if (a == "--worker") opts$worker <- TRUE
    else if (a == "--baseline-rss") opts$baseline_rss <- TRUE
    else if (a == "--cell") opts$cell <- val_next()
    else if (a == "--solver") opts$solver <- val_next()
    else if (a == "--result-out") opts$result_out <- val_next()
    else if (a == "--smoke") opts$smoke <- TRUE
    else if (a == "--n") opts$n <- as.integer(val_next())
    else if (a == "--reps") opts$reps <- as.integer(val_next())
    else if (a == "--warmups") opts$warmups <- as.integer(val_next())
    else if (a == "--max-reps") opts$max_reps <- as.integer(val_next())
    else if (a == "--min-total-duration") opts$min_total_duration <- as.numeric(val_next())
    else if (a == "--threads") opts$threads <- val_next()
    else if (a == "--pin-core") opts$pin_core <- as.integer(val_next())
    else if (a == "--pin-cpus") opts$pin_cpus <- val_next()
    else if (a == "--seed") opts$seed <- as.integer(val_next())
    else if (a == "--out") opts$out <- val_next()
    else if (a == "--sync-registry") opts$sync_registry <- TRUE
    else if (a == "--solvers") opts$solvers <- val_next()
    else if (a == "--problems") opts$problems <- val_next()
    else if (a == "--assert-runnable-tag") opts$assert_runnable_tag <- val_next()
    else if (a == "--max-cells") opts$max_cells <- as.integer(val_next())
    else if (a == "--cell-timeout") opts$cell_timeout <- as.numeric(val_next())
    else if (a == "--jobs") opts$jobs <- as.integer(val_next())
    else stop(sprintf("run_arm.R: unrecognized argument %s", a))
    i <- i + 1L
  }
  opts
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  opts <- .parse_args(argv)
  if (isTRUE(opts$worker)) {
    run_worker(opts$cell, opts$result_out)
    return(invisible(NULL))
  }
  if (isTRUE(opts$baseline_rss)) {
    run_baseline_worker(opts$solver, opts$result_out)
    return(invisible(NULL))
  }
  result <- orchestrate(opts)
  if (isTRUE(result$synced)) {
    # --sync-registry sync-only path: no run, print the sync summary.
    cat(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE), "\n")
    return(invisible(result))
  }
  cat(jsonlite::toJSON(list(n_cells = result$n_cells, n_rows = result$n_rows,
                             n_failures = result$n_failures, runs_path = result$runs_path),
                        auto_unbox = TRUE, pretty = TRUE), "\n")
  invisible(result)
}

.is_main <- function() {
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa) == 0L) return(FALSE)
  invoked <- normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE)
  identical(invoked, normalizePath(RUN_ARM_PATH, mustWork = FALSE))
}

if (.is_main()) {
  invisible(main())
}
