---
phase: 01-verification-coverage-closed
plan: 02
subsystem: testing
tags: [pytest, parity, raking, sinkhorn, convergence-rule, leafblower-6uhm]

# Dependency graph
requires: []
provides:
  - "raking and sinkhorn covered by the explicit-spec weight-parity checks in python/leafblower/test_solver_parity.py (was logit/greg/newton_kl/chebyshev/greenkhorn only)"
  - "raking and sinkhorn covered by the per-method default-rule-resolution lock, mirroring test_logit_default_rule_parity"
  - "_CONV_TOL rationale block extended with measured raking/sinkhorn max_error"
affects: [02-dispatch-unification]

# Actuals (#2632)
actuals:
  tokens: 1450
  tasks: 2
  commits: 2

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Default-rule lock pattern applied to a second and third solver (raking, sinkhorn), confirming test_logit_default_rule_parity's shape generalizes: assert the resolved rule string on BOTH bindings separately, then delegate to _assert_parity for weights — because a fixture that converges to ~1e-16 cannot distinguish stopping rules by weights alone"

key-files:
  created: []
  modified:
    - python/leafblower/test_solver_parity.py

key-decisions:
  - "Hardcoded the observed 'improvement' literal on each binding's rule assertion separately (mirroring test_logit_default_rule_parity exactly), rather than comparing py_rule == r_rule directly — verified by acceptance-criteria's load-bearing check (temporarily breaking one side's expected literal fails exactly that test, not both)"
  - "Measured raking (5.551e-17) and sinkhorn (1.110e-16) max_error before touching _CONV_TOL — both sit far inside the 0.01 threshold, so greg (4.956e-03) remains the binding case and the threshold needs no change, per D-06/D-07's measure-before-touching discipline carried over from 01-01"

patterns-established: []

requirements-completed: [SC2]

coverage:
  - id: D1
    description: "test_raking_parity and test_sinkhorn_parity: explicit-spec (rule=improvement, tol=0.001) weight parity, one-line wrappers delegating to the existing _assert_parity helper"
    requirement: "SC2"
    verification:
      - kind: unit
        ref: "python/leafblower/test_solver_parity.py::test_raking_parity, ::test_sinkhorn_parity (10 passed in file, was 6)"
        status: pass
    human_judgment: false
  - id: D2
    description: "test_raking_default_rule_parity and test_sinkhorn_default_rule_parity: resolved default rule ('improvement') asserted on both R and Python bindings separately, then weight parity, mirroring test_logit_default_rule_parity"
    requirement: "SC2"
    verification:
      - kind: unit
        ref: "python/leafblower/test_solver_parity.py::test_raking_default_rule_parity, ::test_sinkhorn_default_rule_parity"
        status: pass
      - kind: unit
        ref: "load-bearing check: temporarily set expected raking rule literal to a wrong value, confirmed exactly test_raking_default_rule_parity failed, restored, re-confirmed 10/10 green"
        status: pass
    human_judgment: false

duration: 20min
completed: 2026-08-15
status: complete
---

# Phase 1 Plan 2: Raking/Sinkhorn Convergence-Rule Parity Coverage Summary

**Extended `test_solver_parity.py`'s weight-parity and default-rule-lock coverage from five solvers to seven by adding raking and sinkhorn, closing ROADMAP SC2 with four new tests that reuse the file's existing four-step protocol unchanged.**

## Performance

- **Duration:** ~20 min
- **Started:** 2026-08-15 (session start)
- **Completed:** 2026-08-15
- **Tasks:** 2/2
- **Files modified:** 1 (`python/leafblower/test_solver_parity.py`)

## Accomplishments

- `test_raking_parity` and `test_sinkhorn_parity` added as one-line wrappers around the existing `_assert_parity` helper, exactly the shape of `test_chebyshev_parity`/`test_greenkhorn_parity`. No new subprocess helper, comparison function, or epsilon-diff helper was introduced.
- `_CONV_TOL` rationale comment block extended with raking (`5.551e-17`) and sinkhorn (`1.110e-16`) measured `max_error` — the two tightest of every solver in the file. greg (`4.956e-03`) remains the binding case; the `0.01` threshold is unchanged.
- `test_raking_default_rule_parity` and `test_sinkhorn_default_rule_parity` added, mirroring `test_logit_default_rule_parity`'s structure exactly: assert the resolved rule string (`"improvement"`) on the Python binding via `convergence_used`, assert it separately on the R binding, and only then delegate to `_assert_parity` for the weight comparison. Both methods resolve to `rule=improvement` on both bindings on this fixture.
- Module docstring updated (title line + default-rule paragraph) to name raking and sinkhorn alongside logit as methods whose default-rule resolution is locked.

## Task Commits

Each task was committed atomically:

1. **Task 1: Add explicit-spec weight parity for raking and sinkhorn (SC2)** - `1d11b3f` (feat)
2. **Task 2: Lock the default-rule resolution for raking and sinkhorn (SC2)** - `0789008` (feat)

**Plan metadata:** (this commit, following this SUMMARY)

## Files Created/Modified

- `python/leafblower/test_solver_parity.py` - module docstring extended to name raking/sinkhorn; `_CONV_TOL` rationale block extended with two measured values; four new test functions (`test_raking_parity`, `test_sinkhorn_parity`, `test_raking_default_rule_parity`, `test_sinkhorn_default_rule_parity`) added, all delegating to the pre-existing `_assert_parity`/`_run_py`/`_run_r` triplet.

## Measured Values (this fixture, 2026-08-15)

| Method | max\|error\| (explicit spec) | Resolved default rule (both bindings) |
|---|---|---|
| raking | 5.551e-17 | improvement |
| sinkhorn | 1.110e-16 | improvement |

Explicit-spec weight-diff (`max\|w_py − w_r\|`): raking `2.22e-16`, sinkhorn `4.44e-16`. Default-rule-path weight-diff: raking `2.22e-16`, sinkhorn `4.44e-16`. Both five orders of magnitude inside the `rtol=1e-6` parity bound.

## Decisions Made

- **Hardcode the observed literal, not a cross-binding comparison:** initially wrote the default-rule tests comparing `py_rule == r_rule` directly, then rewrote to hardcode `"improvement"` on each side separately (mirroring `test_logit_default_rule_parity`'s exact shape) after re-reading the plan's acceptance criteria: "temporarily changing one expected rule string to a wrong value makes exactly that test fail." A direct `py_rule == r_rule` comparison would pass even if both bindings drifted to the same wrong value together — the hardcoded-literal form is strictly more load-bearing and was verified by the temporary-break-and-restore check (see coverage D2).
- **Measure before touching `_CONV_TOL`:** ran raking/sinkhorn through the existing protocol first (`5.551e-17`, `1.110e-16`), confirmed both sit far inside the `0.01` threshold with greg remaining the binding case, then added the values to the rationale comment — no threshold change was needed or made.
- **Split into two atomic commits per task** despite both tasks touching the same file: assembled all changes first, then isolated Task 1's diff (parity wrappers + `_CONV_TOL` extension) by temporarily removing Task 2's default-rule tests and re-verifying 8/10 tests before the first commit, then restored and re-verified 10/10 before the second commit — keeping the task-level commit granularity the plan's task_commit_protocol requires.

## Deviations from Plan

None — plan executed exactly as written. Both tasks matched their `<action>` and `<acceptance_criteria>` blocks. The only refinement was the hardcoded-literal vs. cross-binding-comparison choice documented above, which is a direct application of the plan's own acceptance criteria (load-bearing assertion requirement), not a deviation from it.

## Issues Encountered

None new. The full Python suite (`cd python && ... pytest -q`) reports the same single pre-existing, out-of-scope failure already documented in `01-01-SUMMARY.md` and `.planning/phases/01-verification-coverage-closed/deferred-items.md`: `leafblower/test_trajectory_csv_smoke.py::test_trajectory_csv_smoke` fails because the CSV header assertion (`iter,errRp`) doesn't match the actual header (`iter,errRp,marginal_kl`). Confirmed unrelated: `test_trajectory_csv_smoke.py` and `src/oris_trajectory.cpp` were untouched by either of this plan's commits. Not fixed here — out of scope per this plan's single-file scope and the phase's test-layer-only boundary. Total suite count moved from 152 to 156 (1 failed, 155 passed), consistent with +4 new tests.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- ROADMAP SC2 satisfied: raking and sinkhorn are covered by both the explicit-spec weight-parity check and the default-rule-resolution lock, so a per-method default-rule divergence between bindings now fails the suite for raking and sinkhorn, not only for logit.
- Combined with 01-01, `test_solver_parity.py` and `test_parity_weights.py` together give Phase 2 (dispatch unification) a wider parity net across more solvers before the dual-dispatch rewire lands.
- The pre-existing `test_trajectory_csv_smoke.py` failure (see Issues Encountered) still needs its own beads ticket before the overall Python suite is fully green — unchanged from 01-01's finding, not introduced or worsened by this plan.

---
*Phase: 01-verification-coverage-closed*
*Completed: 2026-08-15*

## Self-Check: PASSED

- FOUND: `python/leafblower/test_solver_parity.py`
- FOUND: `.planning/phases/01-verification-coverage-closed/01-02-SUMMARY.md`
- FOUND: commit `1d11b3f`
- FOUND: commit `0789008`
