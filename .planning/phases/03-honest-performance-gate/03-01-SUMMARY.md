---
phase: 03-honest-performance-gate
plan: 01
subsystem: testing
tags: [r, benchmarking, oris_soft, survey-calibrate, convergence, honest-performance-gate]

# Dependency graph
requires:
  - phase: 02-one-engine-not-two
    provides: oris_soft as the single authoritative solver for the headline claim (D-04)
provides:
  - "benchmarks/oris_soft_vs_competitors.R: fresh, narrow measurement script for oris_soft vs survey::calibrate on the medium_100k_5margins fixture, each arm run under its own canonical convergence config"
  - "benchmarks/results/oris_soft_vs_competitors.csv: frozen 17-column schema, one row per arm, measured fresh in this run"
  - "benchmarks/results/oris_soft_vs_competitors_env.txt: machine/BLAS/version provenance for the published figures"
  - "benchmarks/run_honest_gate.sh: D-12 single entry point, regenerates the artefacts and runs the opt-in stepstone gate"
  - "tests/testthat/test-bench-gate.R honest-gate test_that() block: opt-in (LBW_BENCH_GATE=1) pass/fail regression on bound compliance and accuracy (D-10/D-11)"
affects: [03-02, 03-03, 03-04]

actuals:
  tokens: 9200
  tasks: 3
  commits: 5

tech-stack:
  added: []
  patterns:
    - "Convergence config must match each solver's own canonical/default stopping rule, not a config borrowed from the competitor being measured against — matching rules artificially rather than matching intent (comparable accuracy target) reintroduces the 'stopped early' confound this phase exists to eliminate."

key-files:
  created:
    - benchmarks/oris_soft_vs_competitors.R
    - benchmarks/run_honest_gate.sh
  modified:
    - benchmarks/results/oris_soft_vs_competitors.csv
    - benchmarks/results/oris_soft_vs_competitors_env.txt
    - tests/testthat/test-bench-gate.R

key-decisions:
  - "oris_soft's convergence config uses its canonical default (metric='marginal_kl', rule='improvement', tol=0.001), passed explicitly rather than via convergence=list(), so the call-site intent is visible. survey::calibrate keeps its own canonical epsilon=1e-3. The two arms are each measured under their own real-world-representative stopping behavior, not an artificially matched rule."
  - "run_honest_gate.sh exports NOT_CRAN=true alongside the thread variables and LBW_BENCH_GATE: without it every skip_on_cran()-gated bench-gate assertion silently no-ops under testthat::test_dir(), independent of LBW_BENCH_GATE."
  - "The new honest-gate test anchors its CSV path on testthat::test_path() rather than a bare relative string, because a filtered test_dir(filter='bench-gate') run never executes test-algo-selection.R's cwd-resetting setwd() and a bare path silently false-skips."

patterns-established:
  - "Confound control extends to convergence RULE choice, not just weight self-reporting: forcing a solver into a competitor's native stopping rule (rather than giving it a comparable target under its own rule) can silently bias a speed/accuracy comparison."
  - "Opt-in testthat gates behind skip_on_cran() require NOT_CRAN=true to actually run under a direct testthat::test_dir() call (devtools::test() sets it internally; a bare test_dir() invocation, as used by wrapper scripts and CI-style commands, does not)."
  - "File paths read inside a filtered testthat::test_dir(filter=...) run should anchor on testthat::test_path(), not a bare repo-relative string — the filter can exclude whichever file happens to reset the working directory for the full suite, leaving cwd at tests/testthat/ for the filtered subset."

requirements-completed: []  # US-003/KPI-04 close across the full plan set (03-01..03-04), not this plan alone; 03-01 supplies the measured fact and both gate mechanisms (D-12 wrapper, D-10/D-11 regression assertion) the later plans build on.

coverage:
  - id: D1
    description: "medium_100k_5margins fixture measured fresh: leafblower_oris_soft (canonical marginal_kl/improvement convergence) vs survey_calibrate (canonical epsilon=1e-3), identical per-observation bounds, identically-computed max_error/max_w/n_eff."
    requirement: "US-003"
    verification:
      - kind: unit
        ref: "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript benchmarks/oris_soft_vs_competitors.R"
        status: pass
      - kind: unit
        ref: "Rscript -e stopifnot(...) CSV schema/row-count/bound checks (task 1 <verify> block)"
        status: pass
    human_judgment: true
    rationale: "Tracer feedback gate (execute-plan.md): the measured numbers are a substantive claim about the package's headline performance and require human sign-off before Task 2/3 (wrapper script, regression gate) build on top of them."
  - id: D2
    description: "One command (benchmarks/run_honest_gate.sh) runs the pre-existing stepstone regression gate and regenerates the new comparison artefacts under the single-thread BLAS envelope (D-12)."
    requirement: "KPI-04"
    verification:
      - kind: unit
        ref: "test -x benchmarks/run_honest_gate.sh && bash -n benchmarks/run_honest_gate.sh && CI=1 bash benchmarks/run_honest_gate.sh && test -f benchmarks/results/oris_soft_vs_competitors.csv"
        status: pass
    human_judgment: false
  - id: D3
    description: "With LBW_BENCH_GATE=1 a testthat assertion reads the CSV and fails when the oris_soft row breaches its bound or accuracy floor; unset, it skips (D-10/D-11)."
    requirement: "US-003"
    verification:
      - kind: unit
        ref: "CI=1 LBW_BENCH_GATE=1 NOT_CRAN=true Rscript -e testthat::test_dir(filter='bench-gate', stop_on_failure=TRUE) -- honest gate: line printed, PASS 5"
        status: pass
      - kind: unit
        ref: "Negative control: perturbed max_w -> FAIL 1 with expected message; restored -> PASS 5"
        status: pass
      - kind: unit
        ref: "devtools::test() default run (no LBW_BENCH_GATE): FAIL 0 | WARN 141 | SKIP 12 | PASS 1839; new assertion SKIP"
        status: pass
    human_judgment: false

duration: ~55min (includes checkpoint fix cycle)
completed: 2026-08-15
status: complete
---

# Phase 3 Plan 1: Honest Performance Gate — Tracer Slice (measurement, wrapper, regression gate) Summary

**oris_soft measured 3.35e-5 max margin error vs survey::calibrate's 1.49e-4 on a 100K-row/5-margin fixture, each solver run under its own canonical convergence rule; one wrapper command regenerates the CSV and an opt-in testthat assertion regresses against it.**

## Performance

- **Duration:** ~55 min (original measurement + checkpoint fix cycle + tasks 2-3)
- **Tasks:** 3 of 3 complete
- **Files modified:** 5 (2 created: `run_honest_gate.sh`, `oris_soft_vs_competitors.R`; 2 regenerated: the CSV, the env.txt; 1 appended-to: `test-bench-gate.R`)

## Accomplishments
- `benchmarks/oris_soft_vs_competitors.R` created: single-thread-BLAS-enforced, fresh (not aggregated from `benchmarks/study/`) measurement of `oris_soft` vs `survey::calibrate` on the `medium_100k_5margins` fixture (100,000 rows, K=5 margins, 4 categories/margin, `bounds_mode="unit"`/`max_weight=3` matched to `survey`'s `bounds=c(0,3)` on `g`).
- Shared `margin_max_error()` helper computes accuracy identically for both arms from their returned weight vectors — never from a package's self-reported convergence number.
- Determinism guard: script refuses to measure (non-zero exit, names the missing variable) unless `OMP_NUM_THREADS`/`OPENBLAS_NUM_THREADS`/`MKL_NUM_THREADS` are all `"1"`.
- **Methodological fix (this checkpoint cycle):** the `leafblower_oris_soft` arm originally passed `convergence = list(absolute = 1e-3)`, which per `R/harvest.R`'s roxygen contract sets `metric="max_err"`/`rule="threshold"` — `survey::calibrate`'s own native IPF stopping behavior, not `oris_soft`'s canonical `metric="marginal_kl"`/`rule="improvement"` default. This forced `oris_soft` to stop the instant it first crossed the 1e-3 threshold instead of continuing to refine under its own plateau-detection rule. Corrected to pass the canonical config explicitly: `convergence = list(metric = "marginal_kl", rule = "improvement", tol = 0.001)`. `survey::calibrate`'s `epsilon = 1e-3` was left unchanged — it is that arm's own canonical stopping parameter.

## Measured Numbers

### Original (flawed) measurement — commit `700b5bc`

Convergence requested for `leafblower_oris_soft`: `convergence=list(absolute=1e-3)` (forces `metric="max_err"`, `rule="threshold"` — foreign to `oris_soft`, native to `survey`).

| arm | wall_s | max_error | max_w | min_w | deff | n_eff | iterations | status | ok |
|---|---|---|---|---|---|---|---|---|---|
| leafblower_oris_soft | 0.0407 | 2.580e-04 | 3.000 | 0.1303 | 1.4801 | 67563.7 | 10 | 0 | TRUE |
| survey_calibrate | 0.5215 | 1.487e-04 | 2.998 | 0.1296 | 1.4829 | 67433.8 | NA | — | TRUE |

Under this config `oris_soft`'s `max_error` (2.58e-4) was LOOSER than `survey`'s (1.49e-4) — the symptom that surfaced the bug: `oris_soft` had stopped as soon as it crossed the 1e-3 target, not because it could not do better.

### Corrected measurement — commit `42f0441`

Convergence requested for `leafblower_oris_soft`: `convergence=list(metric='marginal_kl',rule='improvement',tol=0.001)` (oris_soft's own canonical default). `survey_calibrate` unchanged: `epsilon=1e-3`.

| arm | wall_s | max_error | max_w | min_w | deff | n_eff | iterations | status | ok |
|---|---|---|---|---|---|---|---|---|---|
| leafblower_oris_soft | 0.0418 | 3.347e-05 | 3.000 | 0.1299 | 1.4817 | 67489.4 | 30 | 5 | TRUE |
| survey_calibrate | 0.4902 | 1.487e-04 | 2.998 | 0.1296 | 1.4829 | 67433.8 | NA | — | TRUE |

Under its own canonical convergence rule, `oris_soft` needed more iterations (30 vs 10) and slightly more wall time (0.0418s vs 0.0407s, still ~11.7x faster than `survey`'s 0.4902s) but achieved a max margin error **4.4x tighter** than `survey::calibrate` (3.347e-05 vs 1.487e-04) — consistent with `oris_soft` being the documented higher-precision solver for this package's headline claim (D-04), and no longer an artifact of an unnaturally early stop. `status=5` (constrained-optimum plateau — 2,525 weights are legitimately clamped at the `max_weight=3` bound) is a usable result per the task's own `ok_lb` gate (`status %in% c(0L, 5L)`).

Fixture and machine metadata (both measurements): `m_cell=1024`, `m_cell_over_n=0.01024`, `n=100000`, `n_margins=5`, `n_categories=4`. Machine: AMD Ryzen 9 9950X3D 16-Core Processor, R 4.6.1 (2026-06-24), BLAS `/usr/lib/libblas.so.3.12.0`, LAPACK `/usr/lib/liblapack.so.3.12.0`, `leafblower` 0.1.0, `survey` 4.5, all three thread env vars `=1`.

## Task Commits

1. **Task 1: End-to-end measurement — medium-scale class, oris_soft vs survey::calibrate** - `700b5bc` (feat) — original measurement, later found to use a methodologically flawed convergence config for the `leafblower_oris_soft` arm.
2. **Task 1 fix: canonical convergence** - `42f0441` (fix) — corrected the convergence config to `oris_soft`'s own default and re-measured; superseded the flawed numbers with the table above.
3. **Task 1 summary (tracer checkpoint)** - `2591d31` (docs) — recorded the corrected measurement; human approved at the tracer feedback gate.
4. **Task 2: One-command wrapper (D-12)** - `362d841` (feat) — `benchmarks/run_honest_gate.sh`, mirroring `run_allmethod.sh`'s shape: runs the pre-existing opt-in stepstone gate then regenerates `oris_soft_vs_competitors.csv`/`_env.txt`, all under the single-thread BLAS envelope. Verified end-to-end (`CI=1 bash benchmarks/run_honest_gate.sh`, exit 0, fresh CSV with matching accuracy).
5. **Task 3: Opt-in regression assertion (D-10/D-11)** - `c11762d` (test) — new `test_that()` block in `tests/testthat/test-bench-gate.R` gating on `LBW_BENCH_GATE=1`, asserting bound compliance and the 1e-3 accuracy floor on the `leafblower_oris_soft`/`medium_100k_5margins` row. Negative control verified (perturbed `max_w` → FAIL; restored → PASS). Full `devtools::test()` default run: `[ FAIL 0 | WARN 141 | SKIP 12 | PASS 1839 ]`, new assertion correctly SKIP (opt-in, D-11).

All three tasks complete; plan done.

## Files Created/Modified
- `benchmarks/oris_soft_vs_competitors.R` - measurement script; `leafblower_oris_soft` arm now requests its canonical `marginal_kl`/`improvement` convergence explicitly instead of the `absolute` shorthand that forced `survey`'s native `max_err`/`threshold` rule onto it.
- `benchmarks/results/oris_soft_vs_competitors.csv` - regenerated (task 1 fix, and again via task 2's wrapper — accuracy unchanged run to run, wall-time varies by <2% noise).
- `benchmarks/results/oris_soft_vs_competitors_env.txt` - regenerated (byte-identical content: same machine/versions).
- `benchmarks/run_honest_gate.sh` (task 2, new) - D-12 single entry point: stepstone regression gate then the new measurement, single-thread BLAS envelope, `NOT_CRAN=true` exported (see deviation below).
- `tests/testthat/test-bench-gate.R` (task 3) - new opt-in `test_that()` block; kk1204 block at the end of the file left byte-unchanged (diff is additions-only).

## Task 2/3 Measured Numbers (wrapper-regenerated, post-Task-1 approval)

| arm | wall_s | max_error | max_w | min_w | n_eff |
|---|---|---|---|---|---|
| leafblower_oris_soft | 0.0411 | 3.347e-05 | 3.000 | 0.1299 | 67489.4 |
| survey_calibrate | 0.4830 | 1.487e-04 | 2.998 | 0.1296 | 67433.8 |

Matches task 1's corrected numbers within run-to-run wall-time noise; accuracy figures are bit-identical (deterministic given the fixed seed and single-thread BLAS).

## Negative-Control Evidence (Task 3 acceptance criterion)

1. **FAIL:** perturbed `benchmarks/results/oris_soft_vs_competitors.csv`'s `leafblower_oris_soft` row `max_w` from `3` to `4` (above `max_weight=3`), re-ran `CI=1 LBW_BENCH_GATE=1 NOT_CRAN=true Rscript -e "testthat::test_dir('tests/testthat', filter='bench-gate', stop_on_failure=TRUE)"`: `FAILURE: 'test-bench-gate.R:50:3' — Expected max_w: per-observation bound honoured <= r$max_weight + 1e-10. Actual comparison: 4.0 > 3.0`, `[ FAIL 1 | WARN 0 | SKIP 3 | PASS 4 ]`.
2. **PASS:** restored the CSV from backup, re-ran the identical command: `[ FAIL 0 | WARN 0 | SKIP 3 | PASS 5 ]`, `honest gate: wall_s=0.0411 max_error=3.347e-05 max_w=3.0000 n_eff=67489.4` printed.

## Decisions Made
- Do not force a result: the fix's instruction was to re-measure honestly and report the actual numbers even if `oris_soft` turned out not to be tighter under its canonical config. It did turn out tighter (3.35e-5 vs 1.49e-4), so no further adjustment was made or needed.
- Kept `survey::calibrate`'s `epsilon=1e-3` unchanged — each arm is measured under its own real-world-representative convergence behavior, not a rule matched across packages.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug, caught at tracer checkpoint] Wrong convergence config forced oris_soft into a foreign stopping rule**
- **Found during:** Task 1 tracer feedback gate (human review of `700b5bc`)
- **Issue:** `convergence = list(absolute = 1e-3)` sets `metric="max_err"`/`rule="threshold"` per `R/harvest.R`'s roxygen contract — `survey::calibrate`'s own native IPF stopping rule, not `oris_soft`'s canonical `marginal_kl`/`improvement` default. This made `oris_soft` stop the instant it crossed the 1e-3 target rather than continuing to refine, producing a measured `max_error` (2.58e-4) looser than `survey`'s (1.49e-4) — the exact "stopped early" confound the phase's determinism protocol exists to eliminate.
- **Fix:** Changed the `leafblower_oris_soft` arm's `convergence` argument to `list(metric = "marginal_kl", rule = "improvement", tol = 0.001)` (oris_soft's canonical default, made explicit at the call site). Updated the in-script comment and `note_lb` string to record the corrected config and cite this fix. Re-ran the script under the mandated single-thread BLAS envelope and re-verified the task 1 `<verify>` block, the negative-threads-unset guard, and the two negative greps (`margin_max_error` call count, absence of the unreleased-package name).
- **Files modified:** `benchmarks/oris_soft_vs_competitors.R`, `benchmarks/results/oris_soft_vs_competitors.csv`
- **Verification:** Re-ran `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript benchmarks/oris_soft_vs_competitors.R` (exit 0, both artefacts written) and the task's own `stopifnot()` verify block (schema, one row per arm, `max_w <= max_weight + 1e-10`, `m_cell_over_n < 0.5`) — all pass. Confirmed the determinism guard still refuses to run with `OPENBLAS_NUM_THREADS` unset.
- **Committed in:** `42f0441` (own atomic commit, not squashed into `700b5bc` per project convention: local-only repo, prefer a new commit over amend unless the user explicitly asks for amend).

---

**2. [Rule 3 - blocking issue] `testthat::test_dir(filter='bench-gate')` never actually asserts anything without `NOT_CRAN=true`**
- **Found during:** Task 2, verifying `benchmarks/run_honest_gate.sh` actually exercises the opt-in stepstone gate.
- **Issue:** Every `bench-gate` test opens with `skip_on_cran()`. `devtools::test()` sets `NOT_CRAN=true` internally via pkgload's `test_check()` machinery, but a direct `testthat::test_dir()` call (what the wrapper and the plan's own `<verify>` blocks use) does not. Without it, `skip_on_cran()` fires before the `LBW_BENCH_GATE` check is ever reached — the entire opt-in gate silently reports SKIP regardless of `LBW_BENCH_GATE=1`, which is a false-positive "gate passed" that never actually gated anything. Confirmed empirically: identical invocation without `NOT_CRAN=true` produces `Reason: On CRAN` for all three `bench-gate` tests; with it, the file-existence/env-var skip reasons and the `honest gate:` PASS/FAIL line appear.
- **Fix:** Added `export NOT_CRAN=true` to `benchmarks/run_honest_gate.sh` alongside the three thread variables and `LBW_BENCH_GATE`. Used `NOT_CRAN=true` in the ad hoc verification of Task 3's `<verify>` command (not itself a file edit — the plan's literal verify text omits it, which is the artifact of the same gap).
- **Files modified:** `benchmarks/run_honest_gate.sh`.
- **Verification:** `CI=1 bash benchmarks/run_honest_gate.sh` (exit 0, both stepstone gates and the new measurement step run) and the Task 3 negative-control FAIL/PASS pair above, both of which depend on this fix to reach the assertions at all.
- **Committed in:** `362d841`.

**3. [Rule 1 - bug] Bare-relative `csv_path` silently false-skips under a filtered `test_dir()` call**
- **Found during:** Task 3, following the plan's literal instruction to use a bare `"benchmarks/results/oris_soft_vs_competitors.csv"` string matching the sibling stepstone tests' idiom.
- **Issue:** `testthat::test_dir()` sets the working directory to the directory it is called on (here `tests/testthat/`) for the duration of the run. The full test suite only sees the repo root as cwd because `test-algo-selection.R` performs a top-level `setwd(rprojroot::find_root(...))` that runs once per session and persists for every file executed after it alphabetically — but `filter='bench-gate'` (specified in the plan's own Task 2/3 `<verify>` blocks) excludes that file, so cwd stays at `tests/testthat/` for the whole filtered run and any bare `"benchmarks/..."` path resolves to a nonexistent location, silently false-skipping (this is also why the two *pre-existing* stepstone tests always skip under the filtered invocation — confirmed, but out of scope to fix, since that predates this plan and those files are untouched here).
- **Fix:** Anchored `csv_path` on `testthat::test_path()` (cwd-agnostic by design — resolves to `.` under the filtered run and to the real relative offset under `devtools::test()`), the same pattern `test-sinkhorn-invariants.R` already uses in this directory: `pkg_root <- normalizePath(file.path(testthat::test_path(), "../.."))`. Kept the skip idiom (`skip_on_cran()`, then `skip_if(LBW_BENCH_GATE=="")`, then `skip_if(!file.exists(csv_path))`) verbatim per D-11 — only the path *resolution* changed, not the gate convention.
- **Files modified:** `tests/testthat/test-bench-gate.R`.
- **Verification:** Filtered `test_dir(filter='bench-gate')` run with `LBW_BENCH_GATE=1`/`NOT_CRAN=true` finds the CSV and prints `honest gate: ...` (would otherwise SKIP with `!file.exists(csv_path)`, exactly the failure mode this fix avoids); confirmed under `devtools::test(filter="bench-gate")` too (no explicit `NOT_CRAN`/gate — correctly SKIPs for the *gate* reason, not a path-resolution reason).
- **Committed in:** `c11762d`.

---

**Total deviations:** 3 auto-fixed (1 Rule 1 bug caught by human review at the tracer checkpoint; 1 Rule 3 blocking-issue fix making the opt-in gate actually assert instead of silently no-op; 1 Rule 1 bug fix making the new gate's file-existence check actually find the CSV under the plan's own specified filtered invocation).
**Impact on plan:** Corrects the plan's single most important number before it becomes load-bearing for the rest of the phase (03-02's large-scale figure and 03-03's published threshold both build on this task's methodology), and closes two environment-dependent gaps that would have made the "opt-in regressable gate" (D-10) a false-positive no-op rather than an actual assertion. No scope creep — same two arms, same fixture, same schema, same skip-idiom shape; only the invocation environment (`NOT_CRAN`) and path-resolution mechanism (`test_path()` vs a bare string) changed, both load-bearing for D-10's "hard pass/fail" requirement to be true rather than nominal.

## Issues Encountered
None beyond the deviation above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Plan 03-01 complete: `benchmarks/oris_soft_vs_competitors.R` measures fresh, `benchmarks/run_honest_gate.sh` is the one-command entry point (D-12), and `tests/testthat/test-bench-gate.R` carries an opt-in, verified-to-bite regression assertion (D-10/D-11).
- 03-02 depends on this task's `margin_max_error()` helper, `run_input_class()` fixture/measurement split, and the frozen 17-column CSV schema staying stable — all unchanged by tasks 2/3.
- 03-03 depends on the `honest gate:` assertion site existing with its wall-time-ceiling comment marker (added here, threshold left unset per D-06/D-10) — 03-03 adds the threshold after a human sees the measured values.
- Carry-forward for later plans in this phase: `testthat::test_dir(filter=...)` runs (used by this plan's own wrapper/verify commands and likely reused in 03-02/03-03/03-04's verify blocks) require `NOT_CRAN=true` to actually execute `skip_on_cran()`-gated assertions, and any new file-path reads inside such filtered test files should anchor on `testthat::test_path()` rather than a bare `"benchmarks/..."` string — both gaps are now fixed at the wrapper/assertion level added by this plan, but a *new* opt-in test added in a later plan that doesn't reuse `run_honest_gate.sh` would need the same two fixes independently.

---
*Phase: 03-honest-performance-gate*
*Completed: 2026-08-15*

## Self-Check: PASSED

All 5 task commits (`700b5bc`, `42f0441`, `2591d31`, `362d841`, `c11762d`) confirmed present via `git log --oneline --all`. All 6 claimed files confirmed present on disk: `benchmarks/oris_soft_vs_competitors.R`, `benchmarks/run_honest_gate.sh`, `benchmarks/results/oris_soft_vs_competitors.csv`, `benchmarks/results/oris_soft_vs_competitors_env.txt`, `tests/testthat/test-bench-gate.R`, `.planning/phases/03-honest-performance-gate/03-01-SUMMARY.md`.
