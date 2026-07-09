#!/usr/bin/env Rscript
# benchmarks/study/common/run_matrix.R -- WU-8 (leafblower-2ouc.9).
#
# Pure-logic helpers shared by R/run_arm.R (mirrors common/run_matrix.py's
# Python implementation 1:1) for the registry-driven run-matrix driver:
# applicability gating, run-cell enumeration + randomization, hardware-
# isolation capture, the pre-run frozen-tag SHA-check hook, and timing
# statistics (median/percentile/bootstrap CI).
#
# This file intentionally does NOT source any adapter, arrow, or leafblower
# -- it depends only on jsonlite (already a transitive dep of the shipped
# adapters) so both the driver's orchestrator and its unit tests can
# exercise gating/environment/timing logic without paying subprocess or
# package-load cost.
#
# STRICT SEPARATION (user constraint 2026-07-08): sources nothing from
# leafblower's own src/, r_bridge.cpp, R/, DESCRIPTION, NAMESPACE. No
# registry "applicability pin" field (arm/families/bounds/K_max/builds) is
# ever mutated here -- only `applicable_problems`, which contract.md/
# registry_schema.json document as the WU-10-installed-version-gated,
# cross-checkable derived field this WU is explicitly permitted to
# extend/validate.

.RM_THIS_DIR <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                               mustWork = FALSE)
if (is.na(.RM_THIS_DIR) || !nzchar(.RM_THIS_DIR)) {
  .RM_THIS_DIR <- normalizePath(file.path(getwd(), "benchmarks", "study", "common"), mustWork = FALSE)
}
STUDY_DIR <- normalizePath(file.path(.RM_THIS_DIR, ".."), mustWork = FALSE)
REPO_ROOT <- normalizePath(file.path(STUDY_DIR, "..", ".."), mustWork = FALSE)

FROZEN_TAG <- "benchmark-config-freeze-v1"

# Paths (relative to REPO_ROOT) whose git tree must match the frozen tag's
# tree before any timed cell runs. registry.json is DELIBERATELY excluded --
# contract.md/this WU's own ticket authorize extending/validating its
# `applicable_problems` field, so it is not part of the frozen perimeter.
FROZEN_PATHS <- c(
  "benchmarks/study/spec",
  "benchmarks/study/common/problem_io.R",
  "benchmarks/study/common/problem_io.py",
  "benchmarks/study/common/metrics.R",
  "benchmarks/study/common/metrics.py",
  "benchmarks/study/common/instance_family.R",
  "benchmarks/study/common/instance_family.py",
  "benchmarks/study/R/competitors.R",
  "benchmarks/study/R/leafblower_adapter.R",
  "benchmarks/study/python/competitors.py",
  "benchmarks/study/python/leafblower_adapter.py"
)

THREAD_SWEEP <- c(1L, 4L)

# ott-jax is excluded from ranked {1,4}-thread timing (DESIGN.md Sec2/Sec5,
# WU-OTT leafblower-2ouc.18): XLA ignores OMP_/OPENBLAS_/MKL_NUM_THREADS, so
# a thread-sweep row for it would misrepresent the isolation contract.
RANKED_TIMING_EXCLUDED <- c("ott_jax_sinkhorn")

## ----------------------------------------------------------------------
## Registry / spec loading
## ----------------------------------------------------------------------

rm_load_registry <- function(path = file.path(STUDY_DIR, "registry.json")) {
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

rm_load_problem_specs <- function(spec_dir = file.path(STUDY_DIR, "spec")) {
  skip <- c("schema.json", "runs_schema.json", "registry_schema.json", "status_enum.json",
            "hyperparams.json", "tol_mapping.json",
            "instance_family.json")  # WU-3 manifest (not a spec; per-instance specs live in instance_family/)
  files <- sort(list.files(spec_dir, pattern = "\\.json$", full.names = TRUE))
  files <- files[!basename(files) %in% skip]
  out <- list()
  for (f in files) {
    spec <- jsonlite::fromJSON(f, simplifyVector = FALSE)
    out[[spec$id]] <- spec
  }
  out
}

rm_load_instance_family_specs <- function(spec_dir = file.path(STUDY_DIR, "spec"),
                                          regenerate = TRUE) {
  # Materialize (WU-3 generate_instance_family_specs) + load the 32-instance
  # synthetic family from spec/instance_family/<id>.json. instance_family.R is
  # sourced lazily here (not at run_matrix.R top) so the static-only unit tests
  # keep run_matrix.R's dependency surface at jsonlite. Generation is cheap for
  # any n -- gen: data is resolved lazily at solve time, never here.
  source(file.path(STUDY_DIR, "common", "instance_family.R"))
  inst_dir <- file.path(spec_dir, "instance_family")
  if (isTRUE(regenerate) || !dir.exists(inst_dir)) {
    generate_instance_family_specs(spec_dir)
  }
  files <- sort(list.files(inst_dir, pattern = "\\.json$", full.names = TRUE))
  out <- list()
  for (f in files) {
    spec <- jsonlite::fromJSON(f, simplifyVector = FALSE)
    out[[spec$id]] <- spec
  }
  out
}

rm_load_all_problem_specs <- function(spec_dir = file.path(STUDY_DIR, "spec"),
                                      include_instance_family = TRUE,
                                      regenerate = TRUE) {
  # MATRIX problem set: 4 canonical static specs UNION the WU-3 instance family.
  # rm_load_problem_specs() stays static-only for callers/tests that must see
  # ONLY the canonical 4; the instance family is the sole arena declaring
  # ot/newton_kl objective_families (giving the OT / Newton-KL solvers their
  # applicable problems), so folding it into the static loader would change what
  # static-only consumers observe.
  specs <- rm_load_problem_specs(spec_dir)
  if (isTRUE(include_instance_family)) {
    fam <- rm_load_instance_family_specs(spec_dir, regenerate = regenerate)
    for (id in names(fam)) specs[[id]] <- fam[[id]]
  }
  specs
}

rm_instance_family_n <- function(spec) {
  # Row count n parsed from a gen:instance_family?n=... data_ref; NULL for any
  # non-instance-family spec. Lets --smoke bound the family to small instances.
  data_ref <- spec$data_ref
  if (is.null(data_ref) || !startsWith(data_ref, "gen:instance_family")) return(NULL)
  query <- sub("^[^?]*\\?", "", data_ref)
  for (part in strsplit(query, "&", fixed = TRUE)[[1]]) {
    kv <- strsplit(part, "=", fixed = TRUE)[[1]]
    if (identical(kv[1], "n")) {
      v <- suppressWarnings(as.integer(kv[2]))
      return(if (is.na(v)) NULL else v)
    }
  }
  NULL
}

rm_resolve_spec_path <- function(problem_id, spec_dir = file.path(STUDY_DIR, "spec")) {
  # Canonical spec/<id>.json first, then generated spec/instance_family/<id>.json.
  direct <- file.path(spec_dir, paste0(problem_id, ".json"))
  if (file.exists(direct)) return(direct)
  nested <- file.path(spec_dir, "instance_family", paste0(problem_id, ".json"))
  if (file.exists(nested)) return(nested)
  stop(sprintf("no spec JSON for problem id '%s' under %s or %s/instance_family",
               problem_id, spec_dir, spec_dir))
}

rm_load_hyperparams <- function(path = file.path(STUDY_DIR, "spec", "hyperparams.json")) {
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

rm_solver_package <- function(solver_id, hyperparams) {
  if (startsWith(solver_id, "leafblower_")) return("leafblower")
  for (tbl in c("R_competitors", "Python_competitors")) {
    entry <- hyperparams[[tbl]][[solver_id]]
    if (!is.null(entry) && !is.null(entry$package)) return(as.character(entry$package))
  }
  solver_id
}

## ----------------------------------------------------------------------
## Applicability gating (contract.md Sec4 / registry_schema.json)
## ----------------------------------------------------------------------

rm_problem_bounds_kind <- function(problem_bounds) {
  mx <- problem_bounds$max
  if (is.null(mx) || is.na(mx) || is.infinite(mx)) return("unbounded")
  "bounded"
}

rm_is_applicable <- function(entry, problem_spec) {
  families <- unlist(entry$families)
  obj_families <- unlist(problem_spec$objective_families)
  if (length(intersect(families, obj_families)) == 0L) return(FALSE)
  kind <- rm_problem_bounds_kind(problem_spec$bounds)
  if (!(entry$bounds == "both" || entry$bounds == kind)) return(FALSE)
  k_max <- entry$K_max
  if (!is.null(k_max) && as.integer(problem_spec$K) > as.integer(k_max)) return(FALSE)
  data_ref <- problem_spec$data_ref
  if (!is.null(data_ref) && identical(entry$arm, "python") && grepl("^pkg:", data_ref)) {
    # problem_io.py's _PKG_LOADERS is empty by design (DESIGN.md Sec3):
    # pkg: data_refs are R-package datasets (survey/sampling/anesrake/
    # ebal/optweight/balance), resolvable generically by problem_io.R's
    # requireNamespace()-based loader but NOT by Python (no pkg loader
    # registered). Gate python-arm out rather than schedule a cell that
    # is guaranteed to fail with NotImplementedError at solve time.
    return(FALSE)
  }
  TRUE
}

rm_compute_applicable_problems <- function(registry, problem_specs,
                                           available = NULL, arm = NULL) {
  # `available`+`arm` gate out THIS-ARM registry solvers with no runnable adapter
  # (empty applicable_problems), so the matrix never schedules a cell that would
  # error at adapter-resolution time. Solvers of OTHER arms are left UNTOUCHED
  # (their own driver syncs them) -- a single-arm sync never clobbers the other.
  out <- registry
  ids <- sort(names(problem_specs))
  for (solver_id in names(out$solvers)) {
    entry <- out$solvers[[solver_id]]
    if (!is.null(arm) && !identical(entry$arm, arm)) next  # other arm: untouched
    if (!is.null(available) && !(solver_id %in% available)) {
      out$solvers[[solver_id]]$applicable_problems <- list()
      next
    }
    applicable <- Filter(function(pid) rm_is_applicable(entry, problem_specs[[pid]]), ids)
    out$solvers[[solver_id]]$applicable_problems <- as.list(applicable)
  }
  out
}

## ----------------------------------------------------------------------
## Run-cell enumeration
## ----------------------------------------------------------------------

rm_build_matrix <- function(registry, problem_specs, arm, threads = THREAD_SWEEP,
                             rng_seed = NULL, exclude_ranked_timing = TRUE,
                             available = NULL) {
  cells <- list()
  for (solver_id in names(registry$solvers)) {
    entry <- registry$solvers[[solver_id]]
    if (!identical(entry$arm, arm)) next
    if (!is.null(available) && !(solver_id %in% available)) next  # no runnable adapter this arm
    if (exclude_ranked_timing && solver_id %in% RANKED_TIMING_EXCLUDED) next
    applicable_ids <- unlist(entry$applicable_problems)
    for (pid in applicable_ids) {
      spec <- problem_specs[[pid]]
      if (is.null(spec) || !rm_is_applicable(entry, spec)) next  # defensive re-check
      for (build in unlist(entry$builds)) {
        for (thread in threads) {
          cells[[length(cells) + 1L]] <- list(solver = solver_id, problem = pid,
                                               thread = as.integer(thread), build = build)
        }
      }
    }
  }
  if (length(cells) > 1L) {
    if (!is.null(rng_seed)) set.seed(rng_seed)
    cells <- cells[sample.int(length(cells))]
  }
  cells
}

rm_weights_store_filename <- function(solver, problem, thread, build) {
  sprintf("%s__%s__t%d__%s.parquet", solver, problem, as.integer(thread), build)
}

## ----------------------------------------------------------------------
## Pre-run frozen-tag SHA-check hook
## ----------------------------------------------------------------------

.rm_git <- function(args, repo_root) {
  system2("git", args, stdout = TRUE, stderr = TRUE, wd = NULL)
}

.rm_git_in <- function(args, repo_root) {
  old <- getwd()
  on.exit(setwd(old))
  setwd(repo_root)
  out <- suppressWarnings(system2("git", args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(status = if (is.null(status)) 0L else status, stdout = out)
}

rm_assert_frozen_tag <- function(tag = FROZEN_TAG, paths = FROZEN_PATHS, repo_root = REPO_ROOT) {
  exists_res <- .rm_git_in(c("rev-parse", "-q", "--verify", paste0("refs/tags/", tag)), repo_root)
  if (!identical(exists_res$status, 0L)) {
    msg <- sprintf("WARN: frozen tag '%s' not found -- freeze not yet cut (WU-9 wires enforcement); proceeding unchecked.", tag)
    message(msg)
    return(list(tag_exists = FALSE, ok = TRUE, mismatches = list(), message = msg))
  }
  mismatches <- character(0)
  for (p in paths) {
    tag_obj <- .rm_git_in(c("rev-parse", paste0(tag, ":", p)), repo_root)
    head_obj <- .rm_git_in(c("rev-parse", paste0("HEAD:", p)), repo_root)
    if (!identical(tag_obj$status, 0L) || !identical(head_obj$status, 0L)) {
      mismatches <- c(mismatches, p)
    } else if (!identical(tag_obj$stdout, head_obj$stdout)) {
      mismatches <- c(mismatches, p)
    }
  }
  if (length(mismatches) > 0L) {
    stop(sprintf(
      "SPEC_FAILURE: frozen tag '%s' exists but HEAD has drifted from it under: %s -- halting per CLAUDE.md Anti-Pivot rule (no workaround; re-tag+re-rehearse per contract.md Sec6).",
      tag, paste(mismatches, collapse = ", ")), call. = FALSE)
  }
  list(tag_exists = TRUE, ok = TRUE, mismatches = list())
}

## ----------------------------------------------------------------------
## Hardware isolation capture (DESIGN.md Sec5 Blocker F)
## ----------------------------------------------------------------------

.rm_read_text <- function(path) {
  if (!file.exists(path)) return(NULL)
  v <- tryCatch(readLines(path, warn = FALSE), error = function(e) NULL)
  if (is.null(v) || length(v) == 0L) return(NULL)
  trimws(v[1])
}

.rm_cpu_governors <- function() {
  files <- Sys.glob("/sys/devices/system/cpu/cpu[0-9]*/cpufreq/scaling_governor")
  vals <- vapply(files, .rm_read_text, character(1))
  vals[!vapply(vals, is.null, logical(1))]
}

.rm_turbo_boost_enabled <- function() {
  boost <- .rm_read_text("/sys/devices/system/cpu/cpufreq/boost")
  if (!is.null(boost)) return(identical(boost, "1"))
  no_turbo <- .rm_read_text("/sys/devices/system/cpu/intel_pstate/no_turbo")
  if (!is.null(no_turbo)) return(identical(no_turbo, "0"))
  NA
}

rm_capture_environment <- function(pin_core = NULL) {
  governors <- .rm_cpu_governors()
  governor <- if (length(governors) == 0L) NULL else if (length(unique(governors)) == 1L) governors[1] else "mixed"
  turbo <- .rm_turbo_boost_enabled()
  taskset <- nzchar(Sys.which("taskset"))
  numactl <- nzchar(Sys.which("numactl"))

  governor_ok <- isTRUE(identical(governor, "performance"))
  turbo_ok <- isTRUE(!is.na(turbo) && !turbo)
  pinning_active <- !is.null(pin_core)

  state <- list(
    cpu_governor_all_cores = governor,
    cpu_governor_required_for_timed_runs = "performance",
    cpu_governor_conformant = governor_ok,
    turbo_boost_enabled = if (is.na(turbo)) NULL else turbo,
    turbo_boost_required_for_timed_runs = "disabled",
    turbo_boost_conformant = turbo_ok,
    core_pinning_active = pinning_active,
    core_pinning_pin_core = pin_core,
    core_pinning_mechanism_available = (taskset || numactl),
    core_pinning_mechanism = if (taskset) "taskset" else if (numactl) "numactl" else NULL
  )
  isolated <- governor_ok && turbo_ok && pinning_active
  if (!isolated) {
    gaps <- c(if (!governor_ok) "governor", if (!turbo_ok) "turbo", if (!pinning_active) "pinning")
    message(sprintf(
      "WARN: hardware not fully isolated for timed runs (gaps: %s); see environment.json 'hardware_state_for_timed_runs' -- proceeding (warn-not-fail per WU-8 scope; root-level fix is WU-11's job).",
      paste(gaps, collapse = ", ")))
  }

  list(
    captured_at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
    captured_by = "WU-8 run_arm.R (leafblower-2ouc.9)",
    arm = "R",
    hardware_state_for_timed_runs = state,
    isolated = isolated
  )
}

## ----------------------------------------------------------------------
## Timing statistics: median + 5/95 percentile + percentile-bootstrap CI
## ----------------------------------------------------------------------

.rm_percentile <- function(sorted_vals, q) {
  n <- length(sorted_vals)
  if (n == 1L) return(sorted_vals[1])
  pos <- q * (n - 1) + 1  # 1-indexed
  lo <- floor(pos); hi <- ceiling(pos)
  if (lo == hi) return(sorted_vals[lo])
  frac <- pos - lo
  sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac
}

rm_summarize_timing <- function(reps, n_boot = 2000L, ci = 0.90, rng_seed = NULL) {
  if (length(reps) == 0L) stop("rm_summarize_timing: reps must be non-empty")
  reps_sorted <- sort(reps)
  med <- stats::median(reps)
  if (length(reps) == 1L) {
    return(list(median = med, p05 = reps[1], p95 = reps[1],
                boot_ci_lo = reps[1], boot_ci_hi = reps[1], n_reps = 1L))
  }
  p05 <- .rm_percentile(reps_sorted, 0.05)
  p95 <- .rm_percentile(reps_sorted, 0.95)
  if (!is.null(rng_seed)) set.seed(rng_seed)
  boot_medians <- sort(replicate(n_boot, stats::median(sample(reps, length(reps), replace = TRUE))))
  alpha <- (1 - ci) / 2
  list(
    median = med, p05 = p05, p95 = p95,
    boot_ci_lo = .rm_percentile(boot_medians, alpha),
    boot_ci_hi = .rm_percentile(boot_medians, 1 - alpha),
    n_reps = length(reps)
  )
}
