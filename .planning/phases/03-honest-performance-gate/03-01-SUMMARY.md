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
affects: [03-02, 03-03, 03-04]

actuals:
  tokens: 3808
  tasks: 1
  commits: 2

tech-stack:
  added: []
  patterns:
    - "Convergence config must match each solver's own canonical/default stopping rule, not a config borrowed from the competitor being measured against — matching rules artificially rather than matching intent (comparable accuracy target) reintroduces the 'stopped early' confound this phase exists to eliminate."

key-files:
  created:
    - benchmarks/oris_soft_vs_competitors.R
  modified:
    - benchmarks/results/oris_soft_vs_competitors.csv
    - benchmarks/results/oris_soft_vs_competitors_env.txt

key-decisions:
  - "oris_soft's convergence config uses its canonical default (metric='marginal_kl', rule='improvement', tol=0.001), passed explicitly rather than via convergence=list(), so the call-site intent is visible. survey::calibrate keeps its own canonical epsilon=1e-3. The two arms are each measured under their own real-world-representative stopping behavior, not an artificially matched rule."

patterns-established:
  - "Confound control extends to convergence RULE choice, not just weight self-reporting: forcing a solver into a competitor's native stopping rule (rather than giving it a comparable target under its own rule) can silently bias a speed/accuracy comparison."

requirements-completed: []  # Task 1 only; US-003/KPI-04 close across the full plan set, not this task alone.

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

duration: ~25min (includes checkpoint fix cycle)
completed: 2026-08-15
status: halted
---

# Phase 3 Plan 1: Honest Performance Gate — Task 1 (medium-scale measurement) Summary

**oris_soft measured 3.35e-5 max margin error vs survey::calibrate's 1.49e-4 on a 100K-row/5-margin fixture, each solver run under its own canonical convergence rule — corrected after a tracer-checkpoint review caught the first measurement forcing oris_soft into survey's native stopping rule.**

## Performance

- **Duration:** ~25 min (original measurement + checkpoint fix cycle)
- **Tasks:** 1 of 3 (Task 1 complete and re-verified after fix; Tasks 2-3 not started — halted at the tracer feedback gate pending renewed human approval)
- **Files modified:** 3 (1 created, 2 regenerated)

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

Tasks 2 and 3 (wrapper script, opt-in regression gate) have not started — this plan is halted at the tracer feedback gate pending renewed human approval of the corrected numbers.

## Files Created/Modified
- `benchmarks/oris_soft_vs_competitors.R` - measurement script; `leafblower_oris_soft` arm now requests its canonical `marginal_kl`/`improvement` convergence explicitly instead of the `absolute` shorthand that forced `survey`'s native `max_err`/`threshold` rule onto it.
- `benchmarks/results/oris_soft_vs_competitors.csv` - regenerated with the corrected measurement (see table above).
- `benchmarks/results/oris_soft_vs_competitors_env.txt` - regenerated (byte-identical content: same machine/versions).

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

**Total deviations:** 1 auto-fixed (Rule 1 - bug, caught by human review at the tracer checkpoint before it could propagate into Tasks 2/3 or the published headline figure).
**Impact on plan:** Corrects the plan's single most important number before it becomes load-bearing for the rest of the phase (03-02's large-scale figure and 03-03's published threshold both build on this task's methodology). No scope creep — same two arms, same fixture, same schema.

## Issues Encountered
None beyond the deviation above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Task 1's measurement is now methodologically sound and re-verified; the corrected numbers are ready for renewed tracer-checkpoint human approval.
- Tasks 2 (`benchmarks/run_honest_gate.sh` wrapper) and 3 (opt-in `testthat` regression gate) are NOT started and must not proceed without a fresh approval on the corrected numbers above (checkpoint re-fires because the underlying measurement changed).
- 03-02 and 03-03 depend on this task's `margin_max_error()` helper and CSV schema staying stable; both are unchanged by this fix.

---
*Phase: 03-honest-performance-gate*
*Completed: 2026-08-15 (Task 1 only; plan halted at tracer checkpoint)*
