# Testing Patterns

**Analysis Date:** 2026-08-14

## Test Framework

**R:**
- `testthat` (>= 3.0.0), declared in `DESCRIPTION:19` Suggests.
- Entry point: `tests/testthat.R` — `library(testthat); library(leafblower); test_check("leafblower")`.
- 94 test files under `tests/testthat/test-*.R`.
- **Caveat:** `DESCRIPTION` has NO `Config/testthat/edition: 3` line. Despite the "testthat v3" claim in `CLAUDE.md`, the suite runs under **edition 2 semantics** (edition 3 is opt-in per package). Do not rely on 3e-only behaviour (`expect_snapshot`, stricter `expect_equal` via waldo, deprecation of `expect_equivalent`) without adding that field first.

**Python:**
- `pytest` (`python/pyproject.toml:26`, `[project.optional-dependencies] test`).
- 19 test files inside the package dir `python/leafblower/test_*.py` (so the editable install discovers them). They are excluded from the built wheel via `wheel.exclude = ["leafblower/test_*.py"]` (`python/pyproject.toml:36`).
- Plus `tests/test_parity_weights.py` at repo root.

**C++:** `tests/cpp/` exists but is **empty** — no C++ unit-test framework is wired. Core behaviour is covered only through the R and Python bindings.

## Run Commands (exact)

```bash
# 1. R build gate — this, NOT devtools::install
R CMD INSTALL --preclean .

# 2. R tests
Rscript -e "devtools::test()"
# or, matching the enforcement command:
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", stop_on_failure=TRUE)'

# 3. Python install — uv-managed venv, NO bare pip
cd python && uv pip install -e . --reinstall-package leafblower

# 4. Python + parity tests
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -m pytest -q

# 5. Optional heavy stepstone no-regression gate (local only)
LBW_BENCH_GATE=1 Rscript -e 'testthat::test_file("tests/testthat/test-bench-gate.R")'
```

Use `.venv/bin/python -m pytest`, never bare `python`/`pytest` — a stale shadow `.so` in `~/.local` gets imported instead of the freshly built extension.

## Single-Thread BLAS Requirement

`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1` and `MKL_NUM_THREADS=1` must ALL be exported TOGETHER, and in Python **before `import numpy`**.

**Why:** multi-threaded BLAS/OpenMP reductions sum in a nondeterministic order, so floating-point results differ run to run at the last few ulps. The R↔Python parity assertions are `np.allclose(..., rtol=1e-6, atol=0.0)` against the *same* C++ core — any thread-order drift makes them flaky, and the stepstone benchmark comparison drifts likewise.

The repo enforces this defensively: `python/conftest.py:4-7` sets all three via `os.environ.setdefault` *before* the numpy import, and `.coverage-thresholds.json` bakes them into the blocking enforcement command.

## Definition of Done / Quality Gate

`.coverage-thresholds.json` declares `coverage_model: "behavioral"` with all line/branch/function/statement thresholds `null`. There is **no `covr` and no `pytest-cov`** — do not add a `--cov-fail-under` gate; it would only cover the thin Python layer.

The blocking gate (`enforcement.command`, `blocking: true`) is:
`R CMD INSTALL --preclean .` → R testthat 0 FAIL → `uv pip install -e .` → Python pytest 0 FAIL, single-thread BLAS throughout. Stepstone no-regression is opt-in via `LBW_BENCH_GATE=1`.

## Test File Organization

**R** — flat `tests/testthat/`, one file per behaviour, three naming families:
- Feature: `test-harvest.R`, `test-oris.R`, `test-raking.R`, `test-newton-kl.R`, `test-design.R`.
- Ticket-scoped regression: `test-cr-d16-nbounds.R`, `test-cr-e14-start-weights.R`, `test-cr-f11-margin-kl-zero-obs.R`, `test-xc1s13-vectorization.R` — file name encodes the review finding it locks.
- Invariant/contract: `test-clamp-contract.R`, `test-returned-weights-invariant.R`, `test-unit-bounds-status-consistency.R`, `test-calibration-result-names.R`.

No `helper-*.R` files exist — every test is self-sufficient and constructs its own data inline.

**Python** — mirrors the same families inside `python/leafblower/`: `test_cr_d16_nbounds.py`, `test_cr_e14_start_weights.py`, plus the parity family `test_solver_parity.py`, `test_harvest_na_parity.py`, `test_diagnose_na_parity.py`, `test_design_effect_parity.py`. Every R-side ticket test with a Python-visible surface has a same-named Python counterpart.

**R-side driver scripts for parity:** `tests/parity/run_parity_r.R`, `run_oris_soft_r.R`, `run_chebyshev_r.R`.

## Test Structure

R tests use inline synthetic data with a pinned seed and an explicit convergence spec:

```r
# tests/testthat/test-harvest.R:15-23
test_that("default routing selects ORIS for large complexity", {
  set.seed(1)
  n   <- 200000L
  df  <- data.frame(x = factor(sample(c("a","b","c"), n, replace=TRUE)))
  tgt <- list(x = c(a=0.33, b=0.34, c=0.33))
  result <- harvest(df, tgt, convergence = list(absolute = 1e-6))
  expect_identical(attr(result, "algorithm"), "oris")
})
```

Conventions visible throughout:
- `set.seed()` first, always.
- Pass an **explicit `convergence = list(...)`** so the stopping point is fixed by the spec, not by a per-method default that may change.
- Assert on the documented public surface: returned weights, `attr(result, "algorithm")`, `attr(result, "result")$<field>`.
- Comment the before/after behaviour a regression test locks (`test-harvest.R:44-46`: "Before fix: tol_abs ignored … After fix: tol_abs=0.3 forwarded").

**Expectation mix** (counts across `tests/testthat/`): `expect_equal` 276, `expect_true` 214, `expect_error` 104, `expect_lt` 91, `expect_lte` 41, `expect_gt` 39, `expect_false` 34, `expect_gte` 24, `expect_no_error` 20, `expect_identical` 16, `expect_warning` 14, `expect_no_warning` 4. Numerical results are asserted with inequality bounds (`expect_lt(max_err, tol)`), never bare equality.

**No snapshot tests** — `expect_snapshot` appears in zero files. Reference comparisons go through `.rds` fixtures instead.

## Parity Tests (the distinguishing pattern)

`python/leafblower/test_solver_parity.py` documents the four-step protocol in its module docstring (`:7-11`):

1. Build the Python result with an explicit convergence spec.
2. Build the R result by shelling out — `subprocess.run(["Rscript", "-e", r_script], ...)` (`:149-150`), asserting `proc.returncode == 0` with `proc.stderr` in the message (`:153`).
3. **Precheck both sides converged** (`max_error < _CONV_TOL`, `_CONV_TOL = 0.01` at `:107`) before comparing — so a non-converged run fails loudly instead of silently skipping the assertion (`:90`, `:173-177`).
4. `assert np.allclose(w_r, w_py, rtol=1e-6, atol=0.0)` (`:184`).

Fixture data is **embedded as Python list constants** (`_A`, `_B` at `:43-70`), generated once from R `set.seed(42)`, so both sides run identical data with no filesystem side-effects.

One case deliberately passes NO explicit rule (`convergence={}` / `list()`) to lock the *per-method default* resolution across bindings: `test_logit_default_rule_parity` (`:213`).

Documented non-obvious cases are annotated rather than special-cased: `greenkhorn`/`sinkhorn` are entropic so `sum(w) != n` by construction; `greg` is one-shot so its `max_error ~0.005` is a chi2-scaled residual, not non-convergence.

## Fixtures

`tests/testthat/fixtures/` holds the reference artifacts, loaded with `readRDS()` (`test-calibration-solvers.R:48`, `test-best-iterate.R:66`):
- Stepstone benchmark references: `stepstone_reference.rds`, `stepstone_reference_summary.rds`, `stepstone_best_error_ref.rds`, `stepstone_reference_autumn_only.rds`, `stepstone_small.parquet`, `stepstone_small_targets.rds`.
- Solver references: `oris_kl_reference_stepstone.rds`, `oris_pre_alm_ref.rds`, `raking_obs_reference_stepstone.rds`, `raking_squarem_baseline.rds`.
- Ship-gate JSON baselines: `oris_shipgate_reference.json`, `oris_fixed_omega_baseline.json`.
- Regeneration scripts live beside the data: `stepstone_reference_run.R`, `oris_shipgate_fixture.R`, `oris_baseline_snapshot.R`, `stepstone_verify.R`.

Top-level `tests/testthat/task1_ref.rds` and `task2_oris_ref.rds` are additional pinned references.

Fixture-dependent tests guard on existence and skip rather than error:

```r
# tests/testthat/test-best-iterate.R:49-54
skip_on_cran()
skip_if(!file.exists(ref_path))
skip_if(!file.exists(fx) || !file.exists(tg))
```

## Skips

- `skip_on_cran()` on anything slow or fixture-backed.
- `skip_if(Sys.getenv("LBW_BENCH_GATE") == "")` gates the whole benchmark suite (`test-bench-gate.R:4, 21`).
- `skip_if(Sys.getenv("CI") != "")` excludes timing-sensitive checks from CI (`test-bench-gate.R:30`).
- `skip_if_not_installed("leafblower")` guards binding-level checks (`test-calibration-result-names.R:2`).

## Import-Path Guards

Both conftests exist solely to stop pytest importing the *source* tree (which lacks the compiled `_leafblower` extension) instead of the installed wheel:
- `conftest.py` (repo root) removes `<root>/python` from `sys.path`.
- `python/conftest.py` removes its own directory from `sys.path`, after setting the three BLAS env vars.

## Adding Tests for a New Solver

Steps 6–8 of the 8-step new-solver checklist are test artifacts and all three are required:
6. R test fixture (`.rds`) under `tests/testthat/fixtures/`.
7. Python parity test in `python/leafblower/test_*_parity.py` (pattern: `test_solver_parity.py`).
8. Benchmark fixture for the stepstone regression gate.

## Known Gaps

- `tests/cpp/` is empty — no direct unit tests of `src/` internals; all coverage is via bindings.
- `tests/testthat/_problems/` holds five quarantined files (`test-convergence-criteria-203.R`, `test-ieppa-nonuniform-d-28.R`, `test-ieppa-nonuniform-d-29.R`, `test-raking-89.R`, `test-raking-92.R`) that are not collected by `test_check()`.
- No coverage instrumentation of any kind; coverage claims cannot be measured today.
- No `Config/testthat/edition: 3` (see above).

---

*Testing analysis: 2026-08-14*
