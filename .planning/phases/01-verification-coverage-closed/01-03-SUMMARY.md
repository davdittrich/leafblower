---
phase: 01-verification-coverage-closed
plan: 03
subsystem: testing
tags: [testthat, r, property-testing, bound-invariant, kpi-02, newton_kl]

# Dependency graph
requires: []
provides:
  - "tests/testthat/test-bound-property.R: one-dataset tracer (Task 1), a
    50-dataset stratified property sweep across the 8 bounds-enforcing
    solvers (Task 2), and a dedicated report-not-clamp contract assertion
    pinning newton_kl's documented KPI-02 exception (Task 3/4,
    leafblower-og7d.5)"
affects: [phase-2-solver-work, kpi-02-tracking]

# Actuals (#2632) — chars/4 over the realized diff, this session's work only (Task 3/4).
actuals:
  tokens: 813
  tasks: 2
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
    - "Documented-exception assertion for a solver whose shipped contract
      conflicts with a requirement's literal wording: pin the contract with a
      conditional assertion (IF violating THEN non-OK status AND non-zero
      violation count) rather than silently excluding the solver from
      coverage."

key-files:
  created: []
  modified:
    - "tests/testthat/test-bound-property.R"

key-decisions:
  - "Stratified 50 datasets 17/16/17 across skewed-weights-only, sparse-cells-
    only, and both, per D-04 -- counts chosen to keep sparse-cells-only (the
    weaker clamp driver on its own) below half the sweep while still fully
    represented."
  - "Non-vacuity witness granularity: measured per-dataset (clamp fired on
    >=1 of 8 solvers), not per (dataset, solver) pair -- matches the plan's
    '40 of the 50' wording literally."
  - "Task 3 checkpoint:decision resolved as option-a (developer decision):
    pin newton_kl's shipped report-not-clamp contract as a documented KPI-02
    exception rather than holding the requirement's unconditional wording
    literally across all nine solvers and accepting a permanently red test.
    Recorded on leafblower-og7d.5 via bd comment."
  - "newton_kl's Task 4 assertion reuses the Task 1 tracer fixture
    (.bound_fixture_1, seed 7) rather than a new generator -- it is the exact
    fixture shape measured at plan time to drive the solver out of bounds,
    and reusing it keeps the decision test's non-vacuity guarantee tied to
    already-audited code."

patterns-established:
  - "Design-weight mixture generator (.bp_gen): lognormal bulk (rlnorm) plus a
    Cauchy heavy-tailed contaminant on a small row fraction (abs(rcauchy(...))),
    never Gaussian -- reusable for any future bound-stress fixture (D-03)."

requirements-completed: [KPI-02, SC4]

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
    description: "newton_kl's treatment under KPI-02's unconditional wording, adjudicated by the developer (option-a: pin the shipped report-not-clamp contract as a documented exception) and pinned by an assertion that fails on silent violation"
    requirement: "KPI-02"
    verification:
      - kind: unit
        ref: "tests/testthat/test-bound-property.R#newton_kl reports (not clamps) bound violations -- documented KPI-02 exception (leafblower-og7d.5)"
        status: pass
    human_judgment: false

# Metrics
duration: ~25min (Tasks 1-2, prior session) + ~15min (Tasks 3-4, this session) = ~40min total
completed: 2026-08-15
status: complete
---

# Phase 1 Plan 03: KPI-02 Bound Property Test Summary

**50 fixed, stratified, heavy-tailed R testthat datasets prove the weight-bound invariant unconditionally at 1e-10 across eight bounds-enforcing solvers; the ninth solver, newton_kl, is pinned to its documented report-not-clamp contract as an explicit, developer-adjudicated KPI-02 exception (leafblower-og7d.5).**

## Performance

- **Duration:** ~25 min (Tasks 1-2, prior session) + ~15 min (Tasks 3-4, this session)
- **Tasks:** 4 of 4 complete
- **Files modified:** 1 (`tests/testthat/test-bound-property.R`)

## Accomplishments
- Generalized the Task 1 tracer into a fixed-seed, 50-dataset generator (`.bp_gen`) stratified 17/16/17 into skewed-weights-only, sparse-cells-only, and both.
- Swept the eight bounds-enforcing solvers (`oris`, `raking`, `sinkhorn`, `chebyshev`, `greg`, `oris_soft`, `greenkhorn`, `logit`) against all 50 datasets, asserting `min(w) >= 0.2 - 1e-10` and `max(w) <= 5 + 1e-10` unconditionally (no convergence precheck), per D-05.
- Added a non-vacuity witness: the clamp must fire on at least 40 of the 50 datasets. Measured at **50/50** with the shipped parameters; verified externally that neutralizing the Cauchy contaminant drops the count to **34/50**, confirming the floor is load-bearing.
- Resolved the Task 3 `checkpoint:decision` (developer selected option-a) and implemented it in Task 4: `newton_kl` is excluded from the unconditional 8-solver sweep and instead gets its own assertion pinning its shipped write-guard contract — on the Task 1 tracer fixture (seed 7), which reproducibly drives it out of bounds (`min(w)≈0.005`, `max(w)≈130`), the returned result MUST carry `status != RK_OK` and `n_bounds_violated > 0`. Verified load-bearing by temporarily forcing an OK-status expectation (fails as designed) and restoring.
- Appended the Task 3 decision outcome to the pre-filed ticket `leafblower-og7d.5` via `bd comment` (no duplicate ticket filed).
- Full R suite still reports `FAILED=0 ERROR=0` after both expansions (808 passing expectations in the file alone).

## Task Commits

Each task was committed atomically:

1. **Task 1: One dataset, one solver — bound invariant end to end (tracer)** - `7f2ee42` (test) — completed and human-verified in a prior session.
2. **Task 2: Expand to 50 fixed stratified datasets across the clamping solvers** - `2c9c2b1` (test)
3. **Task 3: Decide how newton_kl is treated under KPI-02** - `checkpoint:decision`, resolved by developer as option-a (no code commit; decision recorded in this plan and on the ticket)
4. **Task 3/4 implementation: pin newton_kl's report-not-clamp contract** - `a11ae8c` (test) — implements the option-a decision and includes the ticket-ID comment per Task 4's requirement

**Plan metadata:** this SUMMARY.md commit (docs) — final metadata commit for this plan. STATE.md and ROADMAP.md are intentionally NOT touched by this run; the orchestrator updates those after all plans in this wave complete.

## Files Created/Modified
- `tests/testthat/test-bound-property.R` — Task 1's one-dataset tracer, Task 2's `.bp_gen` generator with 50 stratified datasets and the sweep `test_that` block, plus Task 3/4's dedicated `newton_kl` report-not-clamp contract assertion.

## Decisions Made

- **Stratum parameters (Claude's Discretion per D-04/CONTEXT.md):** unchanged from the prior session — see Task 2 details below.
  - `skewed-weights-only` (17, seeds 2001-2017): n=1000, category marginals `A=0.40/B=0.35/C=0.25`, `y=P0.55/Q0.45`, `sdlog=1.3`, 6% Cauchy contaminant (`location=20, scale=15`).
  - `sparse-cells-only` (16, seeds 3001-3016): n=1000, category marginals `A=0.90/B=0.07/C=0.03`, `y=P0.5/Q0.5`, `sdlog=0.4`, 2% Cauchy contaminant (`location=10, scale=8`).
  - `both` (17, seeds 4001-4017): n=1200, category marginals `A=0.80/B=0.15/C=0.05`, `y=P0.55/Q0.45`, `sdlog=1.2`, 5% Cauchy contaminant (`location=20, scale=15`).
  - All strata use bounds `min_weight=0.2, max_weight=5`.
- **Task 3 decision — option-a (developer-selected, this session):** pin `newton_kl`'s shipped report-not-clamp contract (`RK_ERR_NOCONV` status + `n_bounds_violated > 0`, documented in `test-newton-bounds-write-guard.R` / leafblower-73d7) as a recorded exception to KPI-02's literal unconditional wording, rather than holding the requirement as written for all nine solvers and accepting a permanently red suite (option-b). Rejected option-b because: (1) it would violate this project's local-only "complete = committed + gates green" definition, (2) the fix belongs in the C++ core, which this test-layer-only phase may not touch, and (3) the suspected mechanism overlaps the dual-dispatch path Phase 2 is scheduled to unify — fixing it here risks being thrown away.
- **Task 4 fixture choice:** reused `.bound_fixture_1()` (the Task 1 tracer, `set.seed(7L)`) for the newton_kl decision assertion rather than a new generator. This is the exact fixture shape measured at plan time (and re-measured this session: `min(w)=0.004960, max(w)=129.532699, status=1, n_bounds_violated=547, n_bounds_clamped=0`) to reproducibly drive the solver out of bounds — reusing it keeps the assertion's non-vacuity tied to an already-read, already-audited fixture instead of introducing a fourth generator into the file.
- **Assertion shape:** conditional, not solver-specific-magic-number: `IF the returned vector violates the bounds THEN status must be non-OK AND n_bounds_violated must be > 0`, plus `n_bounds_clamped == 0` to positively confirm the report-not-clamp mechanism (not a coincidental non-violation). A `expect_true(violates, ...)` guard makes the whole test fail loudly (not silently pass) if a future change to the fixture or the solver stops triggering the violation this assertion exists to pin.
- Kept Task 1's original tracer `test_that` block intact rather than folding it into the generalized generator — unchanged rationale from the prior session.

## Deviations from Plan

None — Tasks 3 and 4 executed exactly per the developer's option-a selection and the plan's Task 4 instructions. No Rule 1-4 auto-fixes were needed.

## Issues Encountered

None. The measured violation magnitudes this session (`min≈0.005, max≈130`) differ numerically from the ticket's original plan-time measurement (`min≈0.005, max≈297`, same seed) — both are far outside `[0.2, 5]` and both trigger the identical reporting contract (`status=1`, `n_bounds_violated>0`, `n_bounds_clamped=0`); the magnitude difference does not affect the assertion, which checks the contract, not a specific numeric value.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness

**Plan complete.** KPI-02 and ROADMAP SC4 are satisfied: the bound invariant is proven unconditionally across the 8 bounds-enforcing solvers on 50 fixed stratified datasets, and the ninth solver's documented contract conflict is adjudicated, pinned by a regression-catching assertion, and ticketed (`leafblower-og7d.5`, scheduled no earlier than Phase 2). No `src/`, `R/`, or `python/` files were touched — the phase's test-layer-only boundary held throughout all four tasks.

---
*Phase: 01-verification-coverage-closed*
*Completed: 2026-08-15*

## Self-Check: PASSED

- FOUND: `tests/testthat/test-bound-property.R`
- FOUND: `.planning/phases/01-verification-coverage-closed/01-03-SUMMARY.md`
- FOUND commit: `7f2ee42` (Task 1)
- FOUND commit: `2c9c2b1` (Task 2)
- FOUND commit: `a11ae8c` (Task 3/4)
