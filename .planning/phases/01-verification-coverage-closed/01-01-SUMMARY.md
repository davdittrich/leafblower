---
phase: 01-verification-coverage-closed
plan: 01
subsystem: testing
tags: [pytest, parity, rk_algorithm_t, leafblower-x7n8]

# Dependency graph
requires: []
provides:
  - "python/leafblower/test_parity_weights.py collected by the blocking Python gate (152 tests, was 141)"
  - "R-vs-Python weight-vector parity across all nine non-AUTO rk_algorithm_t solvers"
  - "Uniform 1e-10 parity tolerance with no unexplained per-method relaxation"
affects: [02-dispatch-unification]

# Actuals (#2632)
actuals:
  tokens: 1164
  tasks: 3
  commits: 3

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Import-time assertion on R-helper path existence converts a wrong REPO_ROOT into a loud pytest collection error instead of a silent pytest.skip"

key-files:
  created: []
  modified:
    - python/leafblower/test_parity_weights.py

key-decisions:
  - "Relocated via git mv (rename-detectable) rather than copy+delete, per D-08"
  - "Measured logit's actual R-vs-Python divergence (5.33e-15) before touching its tolerance, per D-06 — no mechanism justified the prior 1e-6 relaxation, so it was removed rather than kept and documented"

patterns-established:
  - "Module-level `assert <path>.exists()` immediately after computing a REPO_ROOT-derived path, so wrong relocation arithmetic fails collection loudly rather than degrading every dependent test to skip"

requirements-completed: [leafblower-x7n8, SC1, SC3]

coverage:
  - id: D1
    description: "8 weight-vector parity tests relocated from tests/ into python/leafblower/, now collected by the blocking Definition-of-Done gate (149 tests, up from 141)"
    requirement: "leafblower-x7n8"
    verification:
      - kind: unit
        ref: "cd python && pytest --collect-only -q (149 tests collected)"
        status: pass
      - kind: unit
        ref: "python/leafblower/test_parity_weights.py (8 passed, 0 skipped after relocation)"
        status: pass
    human_judgment: false
  - id: D2
    description: "test_weight_parity parametrized over all nine non-AUTO rk_algorithm_t solvers (was six)"
    requirement: "SC1"
    verification:
      - kind: unit
        ref: "python/leafblower/test_parity_weights.py::test_weight_parity[chebyshev|greg|oris_soft] (152 tests collected, 11 passed in file)"
        status: pass
    human_judgment: false
  - id: D3
    description: "logit's unexplained per-method tolerance (1e-6 special case) removed; every method uses a uniform 1e-10 bound with an adjacent measured-evidence comment"
    requirement: "SC3"
    verification:
      - kind: unit
        ref: "python/leafblower/test_parity_weights.py::test_weight_parity[logit] (max diff 5.33e-15 < 1e-10)"
        status: pass
    human_judgment: false

duration: 25min
completed: 2026-08-15
status: complete
---

# Phase 1 Plan 1: Weight-Vector Parity Coverage Closure Summary

**Relocated the orphaned 8-test weight-parity file into the blocking pytest gate, extended it to all nine shipped solvers, and replaced logit's unexplained 1e-6 tolerance special-case with a uniform 1e-10 bound backed by a measured 5.33e-15 divergence.**

## Performance

- **Duration:** ~25 min
- **Started:** 2026-08-15 (session start)
- **Completed:** 2026-08-15
- **Tasks:** 3/3
- **Files modified:** 1 (relocated: `tests/test_parity_weights.py` → `python/leafblower/test_parity_weights.py`)

## Accomplishments
- `tests/test_parity_weights.py`'s 8 R-vs-Python weight-vector parity tests, previously outside the pytest rootdir the Definition-of-Done gate walks, now collect and run inside it (`cd python && pytest --collect-only -q` reads 149, then 152 after Task 2).
- `test_weight_parity` extended from 6 to all 9 non-AUTO `rk_algorithm_t` solvers (added `chebyshev`, `greg`, `oris_soft`), each measured R-vs-Python on the shared four-margin synthetic fixture.
- Every parity tolerance in the file is now uniform (`1e-10`) — the prior `1e-6 if method == "logit" else 1e-10` ternary is gone, replaced by a comment citing the measured value, the date, the range across all nine solvers, and the two candidate mechanisms considered and ruled out.

## Task Commits

Each task was committed atomically:

1. **Task 1: Relocate the parity file into the gate's collection root (leafblower-x7n8, D-08, D-09)** - `4d362af` (feat)
2. **Task 2: Extend the parity matrix to all nine non-AUTO solvers (SC1)** - `1bd3d14` (feat)
3. **Task 3: Retire the unexplained logit parity tolerance (SC3, D-06, D-07)** - `23e7605` (fix)

**Plan metadata:** (this commit, following this SUMMARY)

## Files Created/Modified
- `python/leafblower/test_parity_weights.py` - relocated from `tests/test_parity_weights.py` (git-detected rename); `REPO_ROOT` recomputed for the new depth; import-fallback removed in favour of the bare `leafblower` import; import-time existence assertions added on all three R-helper paths; parametrize list extended to 9 solvers; logit tolerance special-case removed in favour of a uniform, measured-and-commented `1e-10` bound.

## Measured Per-Method Parity (baseline for future regression, 2026-08-15, this fixture)

| Method | max\|w_py − w_r\| |
|---|---|
| greenkhorn | 4.66e-15 |
| logit | 5.33e-15 |
| raking | 4.88e-15 |
| oris | 4.88e-15 |
| sinkhorn | 4.88e-15 |
| newton_kl | 9.77e-15 |
| chebyshev | 4.88e-15 |
| greg | 5.11e-15 |
| oris_soft | 4.88e-15 |

All nine sit 4.66e-15..9.77e-15 — five orders of magnitude inside the uniform `1e-10` bound. Any future run landing materially above this range on the same fixture is a genuine regression, not noise.

## Decisions Made
- **D-08 relocation via `git mv`:** kept the move git-detectable as a rename (`git diff -M` confirms `tests/test_parity_weights.py => python/leafblower/test_parity_weights.py`), not an add+delete pair.
- **D-06/D-07 measure-before-touching:** ran the logit case under the pre-existing 1e-6 tolerance first, recorded 5.33e-15 (matching the plan's pre-recorded 5.329e-15 planning measurement), and only then removed the special case — no ticket was needed since the measured value did not exceed 1e-10.
- **Import-time loud-failure verification:** per the acceptance criteria, deliberately broke `parents[2]` to `parents[1]`, confirmed collection ERRORs (not skips) with a clear `AssertionError` pointing at the wrong resolved path, then restored the correct arithmetic and re-verified 149 tests collected.

## Deviations from Plan

None — plan executed exactly as written. All three tasks matched their `<action>` and `<acceptance_criteria>` blocks without requiring auto-fixes, architectural changes, or scope expansion.

## Issues Encountered

None during execution of this plan's own tasks. One **pre-existing, out-of-scope** failure was discovered while running the plan's own `<verification>` block (`cd python && pytest -q`): `leafblower/test_trajectory_csv_smoke.py::test_trajectory_csv_smoke` fails because the actual CSV header (`iter,errRp,marginal_kl`) no longer matches its hardcoded assertion (`iter,errRp`). Confirmed pre-existing and unrelated — that file and `src/oris_trajectory.cpp` were untouched by any commit in this plan, and the header change predates this plan's HEAD (last touched in `f4b56bd`). Not fixed here per the phase's test-layer-only boundary and this plan's file scope (`tests/test_parity_weights.py` / `python/leafblower/test_parity_weights.py` only). Logged to `.planning/phases/01-verification-coverage-closed/deferred-items.md` for a follow-up beads ticket.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- `leafblower-x7n8` is closable: the gate-collection gap is closed, all 8 relocated tests run inside the blocking Python command.
- ROADMAP SC1 and SC3 satisfied for the weight-vector parity file.
- Phase 2 (dispatch unification) now has a wider net: all nine solvers are R-vs-Python weight-checked inside the blocking gate before the dual-dispatch rewire lands.
- One pre-existing, unrelated test failure (`test_trajectory_csv_smoke.py`) needs its own beads ticket before the overall Python suite is fully green — see Issues Encountered above.

---
*Phase: 01-verification-coverage-closed*
*Completed: 2026-08-15*

## Self-Check: PASSED

- FOUND: `python/leafblower/test_parity_weights.py`
- CONFIRMED-GONE: `tests/test_parity_weights.py`
- FOUND: `.planning/phases/01-verification-coverage-closed/01-01-SUMMARY.md`
- FOUND: commit `4d362af`
- FOUND: commit `1bd3d14`
- FOUND: commit `23e7605`
