---
phase: 03-honest-performance-gate
plan: 03
subsystem: testing
tags: [r, benchmarking, oris_soft, testthat, honest-performance-gate, documentation]

# Dependency graph
requires:
  - phase: 03-honest-performance-gate
    provides: "03-01: tracer measurement + benchmarks/run_honest_gate.sh + honest-gate test_that() block (bound/accuracy assertions, wall-time ceiling deferred). 03-02: full measured table across medium/large/known-limit classes and the full D-07 competitor set."
provides:
  - "tests/testthat/test-bench-gate.R: honest-gate block extended with a regressable wall_s<=0.5s ceiling and n_eff>=60000 floor (D-06/D-10 paired framing), each with an in-comment derivation and independently negative-control-verified."
  - "docs/performance.md: the single linked methodology page every published figure points to — headline claim, results tables, input classes, machine/version provenance, methodology, reproduction command, known limit, competitors."
affects: [03-04]

actuals:
  tokens: 3555
  tasks: 3
  commits: 2

tech-stack:
  added: []
  patterns:
    - "A wall-time gate ceiling should be derived from an explicit argument (e.g. half the slowest competitor's measured time) rather than a fixed multiplier of the measurement itself, and the derivation recorded in-comment with source file, machine line and headroom factor — a threshold with no recorded derivation is the stale-promise pattern this phase retires."
    - "Restoring a perturbed benchmark CSV during a negative-control check must restore from git (git checkout -- <path>), not from a manually-copied backup file — a stale backup left over from an earlier session silently reintroduced wrong data (caught before commit; see Issues Encountered)."

key-files:
  created:
    - docs/performance.md
  modified:
    - tests/testthat/test-bench-gate.R

key-decisions:
  - "Task 1 checkpoint resolved by the developer (verbatim, recorded below): paired framing selected — wall_s<=0.5s AND n_eff>=60000 both asserted alongside the pre-existing bound/accuracy checks."
  - "Wall-time ceiling 0.5s derived as half of the slowest doc-named competitor's measured time (ReGenesees_e_calibrate, 0.9051s) — ~11.7x headroom over the measured 0.0427s, staying meaningfully faster than every competitor under adverse host conditions."
  - "n_eff floor 60000 derived as ~11% headroom below the measured 67489.4, matching the pre-existing max_error<=1e-3 floor's own margin below its measured 3.35e-05."
  - "docs/performance.md's known-limit section presents 03-02's own known_limit_k20_uniform measurement (milder-skew fixture) as the primary figures, and quotes the leafblower-ylsy closure's conclusion (severe-skew fixture) as corroborating qualitative context only — explicitly not conflating the two fixtures' numbers, per 03-02-SUMMARY's carried-forward caveat."

patterns-established:
  - "A paired headline claim (speed + accuracy + bound + n_eff together) is asserted as multiple independent expect_*() calls in one test_that() block, each with its own label= naming the ceiling and class, so a single failing assertion identifies exactly which promise broke without opening the source."

requirements-completed: []  # US-003/KPI-04 close across the full plan set (03-01..03-04), not this plan alone.

coverage:
  - id: D1
    description: "Task 1 checkpoint: developer selected the 'paired' framing (speed never asserted without accuracy/bound/n_eff), with the numeric ceiling (wall_s<=0.5s) and floor (n_eff>=60000) and their derivations, on the record before any code was touched."
    requirement: "US-003"
    verification:
      - kind: manual_procedural
        ref: "Developer's response recorded verbatim in this SUMMARY's 'Task 1: Resolved Decision' section"
        status: pass
    human_judgment: true
    rationale: "A one-way published threshold decision is the human judgment itself — not something automation verifies after the fact."
  - id: D2
    description: "tests/testthat/test-bench-gate.R's honest-gate block extended with wall_s<=0.5s and n_eff>=60000 assertions, each with an in-comment derivation (source CSV/env.txt, headroom argument), alongside the pre-existing bound/accuracy/finiteness checks. Gate stays opt-in; kk1204 block untouched."
    requirement: "US-003"
    verification:
      - kind: unit
        ref: "CI=1 LBW_BENCH_GATE=1 NOT_CRAN=true Rscript -e testthat::test_dir(filter='bench-gate', stop_on_failure=TRUE) -> [ FAIL 0 | WARN 0 | SKIP 3 | PASS 7 ], honest gate: wall_s=0.0427 max_error=3.347e-05 max_w=3.0000 n_eff=67489.4 printed"
        status: pass
      - kind: unit
        ref: "Negative control 1 (wall_s): perturbed CSV wall_s to 0.6 -> FAIL 'wall_s: <=0.5s ceiling...' 0.60 > 0.50; restored via git checkout -- -> PASS again"
        status: pass
      - kind: unit
        ref: "Negative control 2 (n_eff): perturbed CSV n_eff to 50000 -> FAIL 'n_eff: >=60000 floor...' 50000.0 < 60000.0; restored via git checkout -- -> PASS again"
        status: pass
      - kind: unit
        ref: "Opt-in check: CI=1 Rscript -e testthat::test_dir(filter='bench-gate') (LBW_BENCH_GATE unset) | grep -c 'honest gate:' -> 0"
        status: pass
    human_judgment: false
  - id: D3
    description: "docs/performance.md published: headline claim (paired framing), per-class results tables, input-class shapes, verbatim machine/version provenance, methodology (bounds_mode='unit' rationale, independent re-grading of returned weights), one-command reproduction, known-limit section cross-linked to both investigation docs without restating them, competitors section grounded in docs/methods/oris.md."
    requirement: "KPI-04"
    verification:
      - kind: unit
        ref: "Task 3's automated verify block (Rscript stopifnot chain over docs/performance.md vs the CSV's input_class values, run_honest_gate.sh/LBW_BENCH_GATE/methods-oris.md/bounds_mode presence) -> 'OK docs page'"
        status: pass
      - kind: unit
        ref: "Acceptance-criteria greps: BLAS|LAPACK>=1 (6), CPU model string matches env.txt exactly (1), autumn=0, uuid|/home/dd/stepstone=0, bibliography keys=0 with methods/oris.md link present (1)"
        status: pass
    human_judgment: true
    rationale: "SC1/SC4's 'every figure traceable, no detached speed number' requirement needs a human read of the rendered page (done below in the number-to-CSV-cell mapping), not just grep coverage."
  - id: D4
    description: "Project Definition of Done unbroken after both tasks: R CMD INSTALL --preclean . succeeds, devtools::test() reports 0 FAIL at the same counts as the 03-02 baseline, benchmarks/run_honest_gate.sh still exits 0 and regenerates every figure the page publishes."
    requirement: "US-003"
    verification:
      - kind: unit
        ref: "R CMD INSTALL --preclean . (exit 0) && OMP/OPENBLAS/MKL_NUM_THREADS=1 Rscript -e devtools::test() -> [ FAIL 0 | WARN 141 | SKIP 12 | PASS 1839 ] (unchanged from 03-02's baseline)"
        status: pass
      - kind: unit
        ref: "CI=1 bash benchmarks/run_honest_gate.sh -> exit 0; regenerated CSV/env.txt then restored via git checkout -- to keep this plan's diff scoped to tests/testthat/test-bench-gate.R and docs/performance.md only"
        status: pass
    human_judgment: false

duration: ~30min
completed: 2026-08-15
status: complete
---

# Phase 3 Plan 3: Honest Performance Gate — Headline Metric, Hard Ceiling, Methodology Page Summary

**Developer selected the "paired" headline framing (speed never asserted alone); `tests/testthat/test-bench-gate.R` now regresses `wall_s<=0.5s` and `n_eff>=60000` on `leafblower_oris_soft`/`medium_100k_5margins` alongside the existing bound/accuracy checks, and `docs/performance.md` publishes every measured figure with its solver, class, machine, and one-command reproduction path.**

## Task 1: Resolved Decision (recorded verbatim, per output spec)

> **Decision: "paired" framing.**
> - Gate asserts wall_s + max_error + max_w/min_w + n_eff together (4 assertions) on `leafblower_oris_soft` × `medium_100k_5margins`, extending the existing `honest gate:` block in `tests/testthat/test-bench-gate.R:28-56` per Task 2's action text.
> - Wall-time ceiling: **0.5s** — derived as "half of the slowest doc-named competitor's measured time (ReGenesees 0.9051s), giving ~11.7x headroom over the measured 0.0427s while staying meaningfully faster than every competitor even under adverse host conditions." Record this derivation in the code comment per Task 2's action text.
> - n_eff floor: **60000** — measured 67489.4, ~11% headroom below measurement, matching the existing `max_error <= 1e-3` floor's own margin below its measured 3.35e-05 (a floor clearly cleared today but that would catch a real accuracy/deff regression, not host noise).
> - Keep the three existing assertions (`max_w`, `min_w`, `max_error <= 1e-3`) unchanged; add `wall_s <= 0.5` and `n_eff >= 60000` alongside them (four/five total assertions).
> - Input class confirmed: `medium_100k_5margins` (the class the existing gate block cleanly extends into).

## Performance

- **Duration:** ~30 min (continuation from a resolved checkpoint)
- **Tasks:** 3 of 3 complete (task 1: checkpoint resolved by developer, no code; tasks 2-3: implemented)
- **Files modified:** 2 (`tests/testthat/test-bench-gate.R`, `docs/performance.md` created)

## Accomplishments

- **Task 2 (hard gate):** `tests/testthat/test-bench-gate.R`'s `honest gate:` block extended with `expect_lte(r$wall_s, 0.5, ...)` and `expect_gte(r$n_eff, 60000, ...)`, each with an in-comment derivation (measured value, competitor comparison or existing-floor analogy, machine line, headroom factor). The placeholder comment deferring the ceiling to this plan was removed. Both assertions verified independently to bite (see Negative Controls below) and independently to pass on the true measurement. Gate remains opt-in (`LBW_BENCH_GATE` unset -> 0 `honest gate:` lines). The pre-existing kk1204 block is byte-unchanged (`git diff` shows no hunk touching it).
- **Task 3 (docs/performance.md):** Created as the single linked methodology page D-13 calls for. Contains, in the required order: headline claim (paired), per-class results tables (every wall-time figure carries `max_error`/`max_w`/`n_eff` in the same row), input-class shapes, verbatim machine/version transcription, methodology (what was held constant, why `bounds_mode="unit"`), one-command reproduction, known-limit section (cross-linked to both investigation docs, not restating them, explicitly not conflating 03-02's milder-skew measurement with the ylsy closure's severe-skew numbers), and competitors section grounded in `docs/methods/oris.md`'s practitioner table.

## Negative Controls (Task 2 acceptance criterion)

Both new assertions verified independently to bite and to clear, restoring the CSV from `git checkout --` between each perturbation (not a manual backup — see Issues Encountered):

1. **wall_s ceiling — FAIL:** perturbed `leafblower_oris_soft`/`medium_100k_5margins`'s `wall_s` to `0.6` -> `FAILURE: 'test-bench-gate.R:65:3' — Expected wall_s: <=0.5s ceiling on leafblower_oris_soft/medium_100k_5margins <= 0.5. Actual comparison: 0.60 > 0.50`, `[ FAIL 1 | WARN 0 | SKIP 3 | PASS 6 ]`.
2. **wall_s ceiling — PASS:** restored via `git checkout -- benchmarks/results/oris_soft_vs_competitors.csv`, re-ran: `honest gate: wall_s=0.0427 ...`, `[ FAIL 0 | WARN 0 | SKIP 3 | PASS 7 ]`.
3. **n_eff floor — FAIL:** perturbed the same row's `n_eff` to `50000` -> `FAILURE: 'test-bench-gate.R:72:3' — Expected n_eff: >=60000 floor on leafblower_oris_soft/medium_100k_5margins >= 60000. Actual comparison: 50000.0 < 60000.0`, `[ FAIL 1 | WARN 0 | SKIP 3 | PASS 6 ]`.
4. **n_eff floor — PASS:** restored via `git checkout --`, re-ran: `honest gate: wall_s=0.0427 max_error=3.347e-05 max_w=3.0000 n_eff=67489.4`, `[ FAIL 0 | WARN 0 | SKIP 3 | PASS 7 ]`.

## docs/performance.md — Number-to-CSV-Cell Mapping (required by output spec)

Every published figure, mapped to the CSV cell (or env.txt line) it was transcribed from:

| Published number | Source |
|---|---|
| Headline: wall_s=0.0427, max_error=3.35e-05, max_w=3.0000, min_w=0.1299, n_eff=67489.4 | CSV row: `medium_100k_5margins` / `leafblower_oris_soft` |
| Results table, `medium_100k_5margins`, `survey::calibrate`: 0.4932 / 1.487e-04 / 2.9984 / 0.1296 / 67433.8 | CSV row: `medium_100k_5margins` / `survey_calibrate` |
| Results table, `medium_100k_5margins`, `icarus::calibration`: 0.4200 / 6.196e-09 / 2.8729 / 0.0551 / 69075.3 | CSV row: `medium_100k_5margins` / `icarus_calibration` |
| Results table, `medium_100k_5margins`, `ReGenesees::e.calibrate`: 0.9051 / 1.097e-08 / 3.0000 / 0.1298 / 67478.2 | CSV row: `medium_100k_5margins` / `ReGenesees_e_calibrate` |
| Results table, `large_stepstone_fulldata`, `leafblower_oris_soft`: 3.5459 / 9.422e-03 / 3.0000 / ~0 / 904107.9 | CSV row: `large_stepstone_fulldata` / `leafblower_oris_soft` |
| `large_stepstone_fulldata` competitor rows: all `—` | CSV rows: `large_stepstone_fulldata` / {survey,icarus,ReGenesees}, all `ok=FALSE` |
| Results table, `known_limit_k20_uniform`, `leafblower_oris_soft`: 7.3934 / 5.229e-03 / 3.0000 / 1.487e-04 / 205120.3 | CSV row: `known_limit_k20_uniform` / `leafblower_oris_soft` |
| Results table, `known_limit_k20_uniform`, `leafblower_raking_accelerated`: 3.3888 / 1.516e-03 / 3.0000 / 0.0000 / 182761.1 | CSV row: `known_limit_k20_uniform` / `leafblower_raking_accelerated` |
| Input classes table: n/margins/categories/m_cell/m_cell_over_n for all three classes | CSV columns `n`, `n_margins`, `n_categories`, `m_cell`, `m_cell_over_n` (one row per class, constant across arms) |
| Machine/versions table (R version, platform, BLAS, LAPACK, CPU model, thread vars, package versions) | `benchmarks/results/oris_soft_vs_competitors_env.txt`, verbatim |
| "9.9 GB" dense-matrix projection | Computed from CSV's `large_stepstone_fulldata` row: `n=1582732 x n_categories=836 x 8 bytes`; note field of the `ok=FALSE` competitor rows states the same figure |
| Known-limit ylsy closure quote (DEFF 8000-14000, n_eff 71-118) | `.planning/phases/03-honest-performance-gate/03-CONTEXT.md` D-02, itself the recorded verbatim close reason for beads ticket `leafblower-ylsy` |
| kk1204 commit reference `3effd3a` | `docs/investigations/2026-04-23-kk1204-convergence.md` header |

## Task Commits

1. **Task 2: Complete the hard gate with the selected ceiling** - `34a95d1` (test)
2. **Task 3: docs/performance.md — the methodology page every figure links to** - `1afef90` (docs)

_Task 1 (checkpoint:decision) recorded no commit — the developer's selection above is the task's only output, consumed by tasks 2-3._

## Files Created/Modified

- `tests/testthat/test-bench-gate.R` - `honest gate:` block gains `wall_s<=0.5s` and `n_eff>=60000` assertions with in-comment derivations; placeholder deferral comment removed; kk1204 block byte-unchanged.
- `docs/performance.md` (new) - the linked methodology page: headline claim, results tables, input classes, machine/versions, methodology, reproduction, known limit, competitors.

## Decisions Made

See "Task 1: Resolved Decision" above (verbatim) and `key-decisions` in frontmatter.

## Deviations from Plan

None - plan executed exactly as written for tasks 2 and 3. Task 1's decision (recorded above) was made by the developer per the plan's own blocking-checkpoint design, not a deviation.

## Issues Encountered

**Stale backup file during negative-control verification (self-caught, no plan/code impact).** While running the first negative control for task 2, an initial attempt to snapshot `benchmarks/results/oris_soft_vs_competitors.csv` via `cp` inside a shell script that also contained a blocked `python3` heredoc command silently did not execute (the harness blocked the whole script before the `cp` line ran). A subsequent restore from that path picked up a stale file left over from an earlier, unrelated session (matching 03-01's older 2-row measurement, not 03-02's current 10-row CSV) instead of the actual original. This was caught immediately by diffing the restored file against `git diff` (which showed an unexpected 10-line delta instead of "no diff"), and corrected by restoring with `git checkout -- benchmarks/results/oris_soft_vs_competitors.csv` — the sanctioned single-file restore pattern — before any commit was made. No incorrect data reached a commit; the plan's own negative-control verification step is what caught this. All subsequent restores in this plan (including after running `benchmarks/run_honest_gate.sh` for the DoD check, which regenerates the CSV with fresh run-to-run wall-time noise) used `git checkout --` exclusively.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- Plan 03-03 complete: the headline metric and gate ceiling are on the record (developer-selected), the claim is regressable (`wall_s<=0.5s` and `n_eff>=60000` both verified to bite and to clear independently), and `docs/performance.md` is the single page every figure links to.
- 03-04 depends on this plan's published ceiling and `docs/performance.md` existing before it writes the one-line README claim and points to this page, rewrites the US-003/KPI-04 REQUIREMENTS.md rows, adds the `tasks/prd-leafblower-core.md` supersede markers, and re-gates the kk1204 block (untouched by this plan, confirmed via `git diff`).
- Carry-forward: docs/performance.md's known-limit section deliberately keeps 03-02's own `known_limit_k20_uniform` measurement (milder skew) separate from the `leafblower-ylsy` closure's quoted numbers (severe skew) — any future plan citing "the kk1204 known limit" should preserve this distinction rather than treating the two fixtures' figures as interchangeable.

---
*Phase: 03-honest-performance-gate*
*Completed: 2026-08-15*

## Self-Check: PASSED

All 2 task commits (`34a95d1`, `1afef90`) confirmed present via `git log --oneline --all`. All 3 claimed files confirmed present on disk: `tests/testthat/test-bench-gate.R`, `docs/performance.md`, `.planning/phases/03-honest-performance-gate/03-03-SUMMARY.md`.
