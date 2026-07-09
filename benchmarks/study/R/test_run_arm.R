#!/usr/bin/env Rscript
# benchmarks/study/R/test_run_arm.R -- WU-8 (leafblower-2ouc.9).
#
# Manual check()/cat() TDD harness (matches test_leafblower_adapter.R's
# established convention -- no testthat framework used in this dir). Mirrors
# python/test_run_arm.py 1:1: applicability gating against the real
# registry.json + real problem specs, cell-matrix enumeration, the
# frozen-tag hook against the real repo state, hardware-isolation capture
# shape, adapter resolution (including the disclosed cvxr_reference/
# samplics registry entries with no adapter implementation), and two REAL
# end-to-end integration tests: a single worker cell (toy_inline +
# leafblower_oris_r) and a bounded --smoke run, schema-checked against
# runs_schema.json.

.REPO_ROOT <- normalizePath(getwd())
.STUDY_DIR <- file.path(.REPO_ROOT, "benchmarks", "study")
source(file.path(.STUDY_DIR, "R", "run_arm.R"))

.n_pass <- 0L
.n_fail <- 0L
check <- function(desc, cond) {
  ok <- isTRUE(cond)
  if (ok) .n_pass <<- .n_pass + 1L else .n_fail <<- .n_fail + 1L
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", desc))
  if (!ok) cat(sprintf("       cond was: %s\n", paste(deparse(substitute(cond)), collapse = " ")))
}
check_error <- function(desc, expr) {
  ok <- inherits(tryCatch({ force(expr); NULL }, error = function(e) e), "error")
  if (ok) .n_pass <<- .n_pass + 1L else .n_fail <<- .n_fail + 1L
  cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", desc))
}

registry_raw <- rm_load_registry()
problem_specs <- rm_load_problem_specs()
hyperparams <- rm_load_hyperparams()
registry <- rm_compute_applicable_problems(registry_raw, problem_specs)

## ----------------------------------------------------------------------
## Applicability gating
## ----------------------------------------------------------------------

{
  entry <- registry$solvers[["ipfr"]]
  check("ipfr (bounded-only) applicable to stepstone_bounded",
        rm_is_applicable(entry, problem_specs[["stepstone_bounded"]]))
  check("ipfr (bounded-only) NOT applicable to stepstone_unbounded",
        !rm_is_applicable(entry, problem_specs[["stepstone_unbounded"]]))
}

{
  # 4-family taxonomy realignment (2026-07-09): stepstone_bounded/unbounded
  # now declare 'minimax' in objective_families -- leafblower_chebyshev_r
  # (families=['minimax'], bounds='both') is applicable to both stepstones,
  # and still gates out on the kl-only toy_inline/canonical_survey_apistrat.
  entry <- registry$solvers[["leafblower_chebyshev_r"]]
  check("leafblower_chebyshev_r (minimax family) applicable to stepstone_bounded (declares minimax)",
        rm_is_applicable(entry, problem_specs[["stepstone_bounded"]]))
  check("leafblower_chebyshev_r (minimax family) applicable to stepstone_unbounded (declares minimax)",
        rm_is_applicable(entry, problem_specs[["stepstone_unbounded"]]))
  check("leafblower_chebyshev_r (minimax family) gates out on toy_inline (kl-only, no minimax)",
        !rm_is_applicable(entry, problem_specs[["toy_inline"]]))
  check("leafblower_chebyshev_r (minimax family) gates out on canonical_survey_apistrat (kl-only, no minimax)",
        !rm_is_applicable(entry, problem_specs[["canonical_survey_apistrat"]]))
}

{
  # ott_jax_sinkhorn is a KL-family algorithm (method_class='sinkhorn'); 'ot'
  # is retired as an objective_families value (4-family taxonomy realignment,
  # 2026-07-09) -- synthetic probe specs must use 'kl', matching the registry
  # pin, to isolate the K_max gating branch from family matching.
  entry <- registry$solvers[["ott_jax_sinkhorn"]]
  check("ott_jax_sinkhorn families pin == kl (KL-family unification; 'ot' retired as an objective family)",
        identical(unlist(entry$families), "kl"))
  check("ott_jax_sinkhorn K_max pin == 2", identical(entry$K_max, 2L) || identical(entry$K_max, 2))
  synthetic_low_k <- list(objective_families = list("kl"), bounds = list(min = 0, max = NULL), K = 2)
  synthetic_high_k <- list(objective_families = list("kl"), bounds = list(min = 0, max = NULL), K = 9)
  check("ott_jax_sinkhorn applicable to synthetic same-family K<=K_max spec",
        rm_is_applicable(entry, synthetic_low_k))
  check("ott_jax_sinkhorn NOT applicable to synthetic same-family K>K_max spec",
        !rm_is_applicable(entry, synthetic_high_k))
}

{
  # canonical_survey_apistrat has data_ref='pkg:survey::apistrat'. problem_io.py's
  # _PKG_LOADERS is empty (no Python loader registered) -- scheduling this spec
  # on any python-arm solver is guaranteed to fail with NotImplementedError at
  # solve time, so rm_is_applicable() must gate it out regardless of
  # family/bounds/K match. Regression test for a real full-scale --smoke run
  # failure (weightipy/canonical_survey_apistrat). problem_io.R resolves pkg:
  # generically via requireNamespace() -- no such restriction on the R arm.
  spec <- problem_specs[["canonical_survey_apistrat"]]
  check("canonical_survey_apistrat data_ref is pkg:survey::apistrat",
        identical(spec$data_ref, "pkg:survey::apistrat"))
  weightipy <- registry$solvers[["weightipy"]]
  check("weightipy is python arm", identical(weightipy$arm, "python"))
  check("python-arm solver (weightipy) gated OUT of pkg: data_ref spec (canonical_survey_apistrat)",
        !rm_is_applicable(weightipy, spec))
  ebal <- registry$solvers[["ebal"]]
  check("ebal is R arm", identical(ebal$arm, "R"))
  check("R-arm solver (ebal, same family/bounds shape as weightipy) STAYS applicable to pkg: data_ref spec",
        rm_is_applicable(ebal, spec))
}

{
  # 4-family taxonomy realignment (2026-07-09): available specs now declare
  # exactly the authoritative {kl,chi2,logit,minimax} union -- stepstone_
  # bounded/unbounded added 'minimax'; 'ot'/'newton_kl' are retired as
  # objective_families values (they are method_class/algorithm-class values).
  seen_families <- unique(unlist(lapply(problem_specs, function(s) unlist(s$objective_families))))
  check("available specs declare exactly the 4-family taxonomy union {kl,chi2,logit,minimax}",
        setequal(seen_families, c("kl", "chi2", "logit", "minimax")))
  check("no available spec declares retired 'ot' family",
        !("ot" %in% seen_families))
  check("no available spec declares retired 'newton_kl' family",
        !("newton_kl" %in% seen_families))
}

{
  ok <- TRUE
  for (sid in names(registry_raw$solvers)) {
    e0 <- registry_raw$solvers[[sid]]
    e1 <- registry$solvers[[sid]]
    for (k in c("arm", "families", "bounds", "K_max", "builds")) {
      if (!identical(e0[[k]], e1[[k]])) ok <- FALSE
    }
  }
  check("rm_compute_applicable_problems preserves all frozen pin fields", ok)
}

## ----------------------------------------------------------------------
## Cell-matrix enumeration
## ----------------------------------------------------------------------

cells_r <- rm_build_matrix(registry, problem_specs, arm = "R", threads = c(1L, 4L), rng_seed = 1L)
check("build_matrix R arm produces >0 cells", length(cells_r) > 0L)
check("build_matrix R arm excludes ott_jax_sinkhorn (python-only anyway, defensive)",
      !("ott_jax_sinkhorn" %in% vapply(cells_r, function(c) c$solver, character(1))))
check("build_matrix R arm covers thread sweep {1,4}",
      setequal(unique(vapply(cells_r, function(c) c$thread, integer(1))), c(1L, 4L)))
check("build_matrix R arm excludes python-arm solvers",
      length(intersect(vapply(cells_r, function(c) c$solver, character(1)),
                        names(Filter(function(e) identical(e$arm, "python"), registry$solvers)))) == 0L)

lbw_cells_r <- Filter(function(c) startsWith(c$solver, "leafblower_"), cells_r)
check("build_matrix leafblower R rows carry build=='native'",
      length(lbw_cells_r) > 0L && all(vapply(lbw_cells_r, function(c) identical(c$build, "native"), logical(1))))

cells_r_seed2 <- rm_build_matrix(registry, problem_specs, arm = "R", threads = c(1L, 4L), rng_seed = 2L)
check("build_matrix randomizes order across seeds",
      length(cells_r) == length(cells_r_seed2) &&
        !identical(vapply(cells_r, function(c) paste(c$solver, c$problem, c$thread, c$build), character(1)),
                   vapply(cells_r_seed2, function(c) paste(c$solver, c$problem, c$thread, c$build), character(1))))

## ----------------------------------------------------------------------
## Frozen-tag hook + hardware isolation
## ----------------------------------------------------------------------

tag_status <- rm_assert_frozen_tag()  # raises on real drift -- reaching here means ok
check("assert_frozen_tag against real repo: ok", isTRUE(tag_status$ok))

env <- rm_capture_environment()
state <- env$hardware_state_for_timed_runs
check("capture_environment: governor field present", "cpu_governor_all_cores" %in% names(state))
check("capture_environment: turbo field present", "turbo_boost_enabled" %in% names(state))
check("capture_environment: pinning field present", "core_pinning_active" %in% names(state))
check("capture_environment: isolated is logical", is.logical(env$isolated))

single_stat <- rm_summarize_timing(c(0.1))
check("summarize_timing: single rep n_reps==1", identical(single_stat$n_reps, 1L))
multi_stat <- rm_summarize_timing(c(0.10, 0.12, 0.11, 0.13, 0.09), rng_seed = 7L)
check("summarize_timing: p05 <= median <= p95",
      multi_stat$p05 <= multi_stat$median && multi_stat$median <= multi_stat$p95)
check("summarize_timing: boot_ci_lo <= boot_ci_hi", multi_stat$boot_ci_lo <= multi_stat$boot_ci_hi)

## ----------------------------------------------------------------------
## Adapter resolution (including disclosed registry/adapter gaps)
## ----------------------------------------------------------------------

check("resolve_adapter: leafblower_oris_r resolves to a function",
      is.function(.resolve_adapter("leafblower_oris_r")))
check("resolve_adapter: ipfr (R competitor) resolves to a function",
      is.function(.resolve_adapter("ipfr")))
check("cvxr_reference has no run_cvxr_reference() in competitors.R (disclosed gap; not fixed here)",
      !exists("run_cvxr_reference", mode = "function"))
check_error("resolve_adapter: cvxr_reference raises (no adapter implementation)",
            .resolve_adapter("cvxr_reference"))

## ----------------------------------------------------------------------
## Real end-to-end integration: one worker cell
## ----------------------------------------------------------------------

.tmp_dir <- tempfile("run_arm_test_")
dir.create(.tmp_dir)

cell_path <- file.path(.tmp_dir, "cell.json")
result_path <- file.path(.tmp_dir, "result.json")
jsonlite::write_json(
  list(solver = "leafblower_oris_r", problem = "toy_inline", thread = 1L, build = "native",
       warmups = 1L, reps = 2L, max_reps = 5L, min_total_duration = 0.0, n_cap = NULL, seed = 0L),
  cell_path, auto_unbox = TRUE, null = "null"
)

old_wd <- getwd()
setwd(REPO_ROOT)
worker_res <- suppressWarnings(system2(
  "Rscript", c(shQuote(RUN_ARM_PATH), "--worker", "--cell", shQuote(cell_path),
               "--result-out", shQuote(result_path)),
  stdout = TRUE, stderr = TRUE
))
worker_status <- attr(worker_res, "status")
setwd(old_wd)

check(sprintf("worker subprocess (leafblower_oris_r/toy_inline) exits 0 [output: %s]",
              paste(utils::head(worker_res, 5), collapse = " | ")),
      is.null(worker_status) || identical(worker_status, 0L))

if (file.exists(result_path)) {
  out <- jsonlite::fromJSON(result_path, simplifyVector = FALSE)
  rows <- out$rows
  check("worker cell: reps=2 -> 2 timed rows (warmups discarded)", length(rows) == 2L)
  if (length(rows) == 2L) {
    row1 <- rows[[1]]
    check("worker row: exact RUNS_ROW_KEYS shape", identical(sort(names(row1)), sort(RUNS_ROW_KEYS)))
    check("worker row: solver id correct", identical(row1$solver, "leafblower_oris_r"))
    check("worker row: thread==1", identical(as.integer(row1$thread), 1L))
    check("worker row: build=='native'", identical(row1$build, "native"))
    check("worker row: status in harmonized enum", row1$status %in% STATUS_ENUM)
    check("worker row: wall_time_s > 0", as.numeric(row1$wall_time_s) > 0)
    check("worker row: weights_ref file exists on disk",
          file.exists(file.path(REPO_ROOT, row1$weights_ref)))
    check("worker row: weights_ref filename uses bare method (frozen adapter convention: 'leafblower_oris', not 'leafblower_oris_r')",
          identical(basename(row1$weights_ref), "leafblower_oris__toy_inline__t1__native.parquet"))
    refs <- unique(vapply(rows, function(r) r$weights_ref, character(1)))
    check("worker cell: weights_ref identical across reps (deterministic solve)", length(refs) == 1L)
  }
} else {
  check("worker cell produced a result.json", FALSE)
}

## ----------------------------------------------------------------------
## Real end-to-end integration: bounded --smoke run, schema-checked
## ----------------------------------------------------------------------

smoke_out_dir <- file.path(.tmp_dir, "results")
old_wd <- getwd()
setwd(REPO_ROOT)
smoke_res <- suppressWarnings(system2(
  "Rscript",
  c(shQuote(RUN_ARM_PATH), "--smoke", "--solvers", shQuote("^(leafblower_oris_r|ipfr)$"),
    "--threads", "1", "--out", shQuote(smoke_out_dir), "--seed", "3"),
  stdout = TRUE, stderr = TRUE
))
smoke_status <- attr(smoke_res, "status")
setwd(old_wd)

check(sprintf("smoke run exits 0 [output: %s]", paste(utils::tail(smoke_res, 8), collapse = " | ")),
      is.null(smoke_status) || identical(smoke_status, 0L))

runs_path <- file.path(smoke_out_dir, "runs.parquet")
check("smoke run: runs.parquet exists", file.exists(runs_path))

if (file.exists(runs_path)) {
  df <- arrow::read_parquet(runs_path)
  check("smoke run: runs.parquet has >0 rows", nrow(df) > 0L)

  schema <- jsonlite::fromJSON(file.path(.STUDY_DIR, "spec", "runs_schema.json"), simplifyVector = FALSE)
  required <- unlist(schema$items$required)
  check("smoke run: runs.parquet has all required columns", all(required %in% names(df)))
  check("smoke run: status values within harmonized enum",
        all(as.character(df$status) %in% STATUS_ENUM))
  check("smoke run: thread values within {1,4}", all(df$thread %in% c(1L, 4L)))
  check("smoke run: build values within {portable,native,na}",
        all(as.character(df$build) %in% c("portable", "native", "na")))
  check("smoke run: wall_time_s all positive", all(df$wall_time_s > 0))
  check("smoke run: rep starts at 0", min(df$rep) == 0L)
}

env_path <- file.path(smoke_out_dir, "environment.json")
check("smoke run: environment.json exists", file.exists(env_path))
if (file.exists(env_path)) {
  env_out <- jsonlite::fromJSON(env_path, simplifyVector = FALSE)
  check("smoke run environment.json: hardware_state_for_timed_runs present",
        !is.null(env_out$hardware_state_for_timed_runs))
  check("smoke run environment.json: frozen_tag present", !is.null(env_out$frozen_tag))
  check("smoke run environment.json: baselines present", !is.null(env_out$baselines))
}

## ----------------------------------------------------------------------

cat(sprintf("\n%d passed, %d failed\n", .n_pass, .n_fail))
if (.n_fail > 0L) quit(status = 1L)
