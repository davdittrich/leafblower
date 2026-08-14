---
phase: 01-verification-coverage-closed
plan: 03
subsystem: testing
tags: [testthat, r, property-testing, bound-invariant, kpi-02]

# Dependency graph
requires: []
provides:
  - "tests/testthat/test-bound-property.R: one-dataset tracer (Task 1) plus a
    50-dataset stratified property sweep (Task 2) asserting the KPI-02 weight-
    bound invariant unconditionally at 1e-10 across the eight bounds-enforcing
    solvers"
affects: [phase-2-solver-work, kpi-02-tracking]

# Actuals (#2632) — chars/4 over the realized diff, this session's work only (Task 2).
actuals:
  tokens: 9500
  tasks: 1
  commits: 1

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Fixed-seed literal per dataset (set.seed(N) inline before each .bp_gen()
      call) rather than a seed vector consumed in a loop -- keeps every dataset
      independently replayable and grep-auditable (D-02)."
    - "Non-vacuity witness: track whether the clamp counter (n_bounds_clamped)
      fired on at least one swept solver per dataset, assert the count against
      a floor with real headroom below the measured value."

key-files:
  created: []
  modified:
    - "tests/testthat/test-bound-property.R"

key-decisions:
  - "Stratified 50 datasets 17/16/17 across skewed-weights-only, sparse-cells-
    only, and both, per D-04 -- counts chosen to keep sparse-cells-only (the
    weaker clamp driver on its own) below half the sweep while still fully
    represented."
  - "Non-vacuity witness measured per-dataset (clamp fired on >=1 of 8 solvers),
    not per (dataset, solver) pair -- matches the plan's '40 of the 50' wording
    literally."
  - "newton_kl excluded from the Task 2 sweep, per plan -- its shipped
    reporting-not-clamping contract is adjudicated in Task 3 (blocking
    checkpoint, not yet resolved this run)."

patterns-established:
  - "Design-weight mixture generator (.bp_gen): lognormal bulk (rlnorm) plus a
    Cauchy heavy-tailed contaminant on a small row fraction (abs(rcauchy(...))),
    never Gaussian -- reusable for any future bound-stress fixture (D-03)."

requirements-completed: []  # KPI-02/SC4 NOT yet closed -- Task 3 (blocking checkpoint) and Task 4 remain.

coverage:
  - id: D1
    description: "50 fixed, stratified, heavy-tailed datasets assert min(w) >= min_weight and max(w) <= max_weight within 1e-10, unconditionally, across the 8 bounds-enforcing solvers"
    requirement: "KPI-02"
    verification:
      - kind: unit
        ref: "tests/testthat/test-bound-property.R#bound invariant holds across 50 fixed stratified datasets and 8 solvers (KPI-02)"
        status: pass
    human_judgment: false
  - id: D2
    description: "newton_kl's treatment under KPI-02's unconditional wording (pin the shipped report-only contract vs. hold the requirement as written and accept a red test)"
    verification: []
    human_judgment: true
    rationale: "Task 3 is an explicit blocking checkpoint:decision (gate=\"blocking\") requiring a developer choice between two documented options with measured evidence; not resolvable by the executor. Not yet decided this run."

# Metrics
duration: ~25min (Task 2 only, this session)
completed: 2026-08-15
status: halted
---

# Phase 1 Plan 03: KPI-02 Bound Property Test (Tasks 1-2) Summary

**50 fixed, stratified, heavy-tailed R testthat datasets prove the weight-bound invariant unconditionally at 1e-10 across all eight bounds-enforcing solvers; the ninth solver's contract conflict awaits a blocking developer decision (Task 3).**

## Performance

- **Duration:** ~25 min (this session, Task 2 only; Task 1 was executed and human-verified in a prior session)
- **Tasks:** 2 of 4 complete (Task 1 tracer, Task 2 stratified sweep); Task 3 (checkpoint:decision) reached and returned to the developer; Task 4 not started
- **Files modified:** 1 (`tests/testthat/test-bound-property.R`)

## Accomplishments
- Generalized the Task 1 tracer into a fixed-seed, 50-dataset generator (`.bp_gen`) stratified 17/16/17 into skewed-weights-only, sparse-cells-only, and both.
- Swept the eight bounds-enforcing solvers (`oris`, `raking`, `sinkhorn`, `chebyshev`, `greg`, `oris_soft`, `greenkhorn`, `logit`) against all 50 datasets, asserting `min(w) >= 0.2 - 1e-10` and `max(w) <= 5 + 1e-10` unconditionally (no convergence precheck), per D-05.
- Added a non-vacuity witness: the clamp must fire on at least 40 of the 50 datasets (via `n_bounds_clamped` on at least one swept solver). Measured at **50/50** with the shipped parameters; verified externally that neutralizing the Cauchy contaminant (a constant `x1.0` multiplier in place of the heavy-tailed draw) drops the count to **34/50**, confirming the floor is load-bearing, not decorative.
- Full R suite (94 files) still reports `FAILED=0 ERROR=0` after the expansion; new file alone runs in ~1.3s (804 passing expectations), well inside the "fast" requirement.
- Reached Task 3, the `checkpoint:decision` on `newton_kl`'s treatment under KPI-02, and returned control to the developer per this run's explicit instruction rather than auto-selecting an option.

## Task Commits

Each task was committed atomically:

1. **Task 1: One dataset, one solver — bound invariant end to end (tracer)** - `7f2ee42` (test) — completed and human-verified in a prior session.
2. **Task 2: Expand to 50 fixed stratified datasets across the clamping solvers** - `2c9c2b1` (test)

**Plan metadata:** not yet committed — plan is not complete (Task 3 checkpoint pending, Task 4 not started). No STATE.md/ROADMAP.md update in this run per orchestrator instruction (wave owns those writes after all plans in the wave complete).

## Files Created/Modified
- `tests/testthat/test-bound-property.R` — Task 1's one-dataset tracer (unchanged) plus Task 2's `.bp_gen` generator, 50 literal `set.seed(...)` dataset-construction lines (17 skewed-weights-only / 16 sparse-cells-only / 17 both), and the stratified sweep `test_that` block with its non-vacuity witness.

## Decisions Made
- **Stratum parameters (Claude's Discretion per D-04/CONTEXT.md):**
  - `skewed-weights-only` (17, seeds 2001-2017): n=1000, category marginals `A=0.40/B=0.35/C=0.25` (near-balanced, no sparse cell), `y=P0.55/Q0.45`, `sdlog=1.3`, 6% Cauchy contaminant (`location=20, scale=15`).
  - `sparse-cells-only` (16, seeds 3001-3016): n=1000, category marginals `A=0.90/B=0.07/C=0.03` (a 3% minority cell), `y=P0.5/Q0.5`, `sdlog=0.4` (mild), 2% Cauchy contaminant (`location=10, scale=8`) — deliberately milder weight skew so the failure axis under stress is cell sparsity, not weight distribution.
  - `both` (17, seeds 4001-4017): n=1200, category marginals `A=0.80/B=0.15/C=0.05`, `y=P0.55/Q0.45`, `sdlog=1.2`, 5% Cauchy contaminant (`location=20, scale=15`) — identical shape to the Task 1 tracer fixture, replicated across 17 seeds.
  - All strata use bounds `min_weight=0.2, max_weight=5` (same as Task 1).
- **Non-vacuity witness granularity:** engagement tracked per-dataset (fired on >=1 of 8 solvers), not per (dataset, solver) pair — matches the plan's literal "engaged on at least 40 of the 50" wording. Measured 50/50 engaged; floor set at 40 for headroom while remaining a real regression signal (verified: contaminant-neutralized run drops to 34/50, below the floor).
- **Kept Task 1's original tracer test_that block intact** rather than folding it into the generalized generator — it is a cheap (<0.1s), independently-readable sanity check with its own hand-picked fixture, and nothing in Task 2's acceptance criteria required removing it.
- **Did not resolve Task 3.** Per this run's explicit instruction, the checkpoint:decision is presented to the developer rather than auto-selected, even though the plan's own `gate="blocking"` (not `blocking-human`) would otherwise permit auto-selection under an active auto-mode config.

## Deviations from Plan

None — Task 2 executed exactly as specified. No Rule 1-4 auto-fixes were needed; the acceptance criteria were satisfied on the first parameter set validated externally (prototype sweep, not committed to the repo).

## Issues Encountered
None.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness

**Blocked on Task 3.** This plan cannot be marked complete, and KPI-02/SC4 cannot be marked satisfied, until the developer selects a `newton_kl` treatment (option-a: pin the shipped report-only contract and ticket the requirement-scope exception; option-b: hold KPI-02 as written for all nine solvers and accept a red test until the solver is fixed; or a described alternative). See the checkpoint returned alongside this summary for full context, measured evidence, and options.

Once Task 3 resolves, Task 4 implements the selected treatment (test-only, no `src/`/`R/`/`python/` changes), appends the decision outcome to `leafblower-og7d.5` via `bd comment`, and this SUMMARY should be regenerated with `status: complete`, `requirements-completed: [KPI-02, SC4]`, and the plan-metadata commit.

---
*Phase: 01-verification-coverage-closed*
*Completed: partial — halted at Task 3 checkpoint, 2026-08-15*
