# Phase 2: One Engine, Not Two - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-08-15
**Phase:** 2-One Engine, Not Two
**Areas discussed:** Dispatch migration strategy, SC3 water-fill duplication scope, SC4 build-list divergence gate

---

## Dispatch migration strategy (leafblower-rywn, P0)

| Option | Description | Selected |
|--------|-------------|----------|
| Incremental, solver-by-solver | One solver migrated + verified per task/commit before the next; a regression bisects to one solver | ✓ |
| Atomic single cutover | All 9 solvers rewired through `rk_calibrate()` in one commit | |

**User's choice:** Incremental, solver-by-solver.
**Notes:** Given P0 risk, no CI, and RAII/ABI tripwires, per-solver verification (R↔Python
parity + full DoD gate) before advancing was preferred over a single large commit. A
separate rollback/bake-in safety-net option (feature-flag fallback to old dispatch) was
offered as its own discussion area but not selected — git-revert at solver granularity is
the accepted rollback mechanism (recorded as D-02).

---

## SC3 water-fill duplication scope

| Option | Description | Selected |
|--------|-------------|----------|
| Verify-and-close | Treat SC3 as already satisfied by existing shared-helper architecture; add a verification task instead of a code change | ✓ |
| I have a specific duplication in mind | User names the actual remaining duplication | |

**User's choice:** Verify-and-close.
**Notes:** Presented with direct-read evidence that every solver (oris, raking, chebyshev,
greenkhorn, greg, logit_calib, sinkhorn) already routes through the shared
`finalize_weights`/`finalize_weights_buf` in `calib_dispatch.hpp`, including an explicit
"Delegate to the single source" comment in `oris_finalize.cpp` documenting a prior
consolidation. Raking's `water_fill_cat` is different math (inline IPF box projection),
not the SC3 target. D-04 records that a genuine remaining duplication found during deeper
research/planning supersedes this and should be flagged, not silently dropped.

---

## SC4 build-list divergence gate

| Option | Description | Selected |
|--------|-------------|----------|
| Wired into the DoD gate | Test assertion running as part of `.coverage-thresholds.json`'s `enforcement.command`, fires automatically | ✓ |
| Separate manual script/Makefile target | Standalone check run by hand, not part of the automatic gate | |

**User's choice:** Wired into the DoD gate.
**Notes:** No CI exists in this repo — the DoD gate command IS the enforcement mechanism,
so a manual-only check would be easy to forget. Left to researcher/planner discretion
whether the assertion lives in the R or Python test suite.

---

## Claude's Discretion

- Whether the SC4 sync-check assertion lives in the R or Python test suite.
- Exact task boundaries and ordering for the solver-by-solver migration (which solver
  first, how AUTO-routing consolidation is sequenced relative to the 9 individual solver
  migrations).

## Deferred Ideas

None. The rollback/bake-in safety-net area was presented as a discussion option but not
selected by the user — not a scope-creep deferral, just an unselected option (see D-02 in
CONTEXT.md for the accepted alternative).
