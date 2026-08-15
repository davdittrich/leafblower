---
phase: 02-one-engine-not-two
plan: 06
subsystem: api
tags: [cpp17, dispatch, r-bridge, newton-kl, calib_dispatch]

requires:
  - phase: 02-one-engine-not-two
    plan: 05
    provides: "lbw::DispatchResult / lbw::dispatch_solver shared table extended through oris, oris_soft — this plan adds the ninth and final named-method arm"
provides:
  - "RK_ALG_NEWTON_KL case arm in lbw::dispatch_solver (calib_dispatch.hpp), including the two newton-only diagnostics (n_projected_dims, lm_mu_final) and the empty-best-weights zero-fill sentinel"
  - "newton_kl reached from both c_api.cpp::rk_calibrate() and r_bridge.cpp::C_rk_calibrate() through the same dispatch call"
  - "Confirmation that all 9 non-AUTO dispatch_solver arms assign stall_kind, and that it reaches harvest.R's production convergence_stall_kind consumer unbroken"
affects: [phase-02-plan-07-auto-consolidation, phase-02-plan-08-single-dispatch-site-test]

actuals:
  tokens: 2452
  tasks: 2
  commits: 1

tech-stack:
  added: []
  patterns:
    - "Reused the existing pack_dispatch_oris_extras_c (added in plan 02-05) for c_api.cpp's newton_kl branch instead of writing a near-duplicate 'newton extras' packer: DispatchResult's default-constructed ORIS-only-field values are byte-identical to what the pre-migration pack_newton_result_c hardcoded as the non-ORIS reset, so calling the existing helper reproduces the reset exactly with zero new code."

key-files:
  created: []
  modified:
    - src/calib_dispatch.hpp
    - src/c_api.cpp
    - src/r_bridge.cpp

key-decisions:
  - "c_api.cpp's explicit RK_ALG_NEWTON_KL branch: verified field-by-field that DispatchResult's default constructor values for every ORIS-only field (min_alpha_seen=1.0, final_alpha=1.0, homotopy_final_factor=1.0, all 9 sor_* / 4 alm_* / remaining homotopy fields at 0/0.0) exactly match what the pre-migration pack_newton_result_c hardcoded as the post-fallback reset. Since the new RK_ALG_NEWTON_KL dispatch arm never touches those DispatchResult members, calling the existing pack_dispatch_oris_extras_c (added in plan 02-05, previously only used by the two oris branches) after pack_dispatch_result_c reproduces byte-identical output to the old explicit reset block, with no new function. pack_newton_result_c itself is left fully intact and still used by the AUTO-fallback branch (c_api.cpp:606), which this plan does not touch."
  - "r_bridge.cpp's explicit newton_kl branch now sets res_n_bounds_clamped explicitly (dres.n_bounds_clamped), which pre-migration was set implicitly via pack_solver_result's has_n_bounds trait on NewtonCalibResult — same value (0, since newton_calib.cpp never sets n_bounds_clamped), just now assigned via the dispatch-result copy pattern shared by all 9 migrated solvers instead of a template trait. No behavior change."
  - "Both res_n_projected_dims and res_lm_mu_final (the two superset-only R diagnostics DispatchResult already declared as scaffolding in plan 02-01) are now populated inside the RK_ALG_NEWTON_KL dispatch arm and copied into the res_* locals at the r_bridge.cpp call site, exactly matching pre-migration values."
  - "Did not add a shared finalize_weights call to newton_calib.cpp, per the plan's explicit prohibition — the unit-mode redistribution gap remains a separately-tracked pre-existing issue (leafblower-og7d.5), unrelated to this dispatch refactor."

requirements-completed: [US-004]

coverage:
  - id: D1
    description: "R and Python both reach lbw::newton_calibrate through lbw::dispatch_solver for method newton_kl"
    requirement: US-004
    verification:
      - kind: unit
        ref: "tests/testthat filter=\"newton\" (5 files): 0 FAIL, 25 PASS, 2 WARN (pre-existing, unrelated to migration); each of the 3 tripwire files (test-newton-kl.R, test-newton-tsvd-projection.R, test-newton-kl-tsvd-ratio.R) also run individually: 0 FAIL / 0 ERROR each"
        status: pass
      - kind: integration
        ref: "python -m pytest -k newton: 8 passed, 0 failed"
        status: pass
    human_judgment: false
  - id: D2
    description: "R's newton_kl result still contains n_projected_dims and a finite lm_mu_final on the converging fixture, still surfaces n_bounds_violated, and the empty-best-weight zero-fill sentinel still fires when the >5% violation guard trips"
    requirement: US-004
    verification:
      - kind: unit
        ref: "test-newton-kl.R T1/T2 (n_projected_dims/lm_mu_final finite), T4 (>5% violation -> RK_ERR_NOCONV + zero-filled sentinel), test-newton-tsvd-projection.R (n_projected_dims pinned against pre-migration value) — all pass"
        status: pass
    human_judgment: false
  - id: D3
    description: "Every explicitly-named method now reaches its solver through lbw::dispatch_solver; only the AUTO branch retains its own routing copy"
    verification:
      - kind: other
        ref: "grep -c 'case RK_ALG_' src/calib_dispatch.hpp: 9 (SINKHORN, GREG, GREENKHORN, LOGIT, CHEBYSHEV, RAKING, ORIS, ORIS_SOFT, NEWTON_KL) — all 9 non-AUTO enum values now have a dispatch arm; c_api.cpp/r_bridge.cpp's method='auto' branches still call lbw::newton_calibrate/oris_solve/raking_solve directly, untouched (plan 07's scope)"
        status: pass
    human_judgment: false
  - id: D4
    description: "Every arm of lbw::dispatch_solver assigns stall_kind, and it reaches harvest.R's production convergence_stall_kind consumer (accelerate heuristic) unbroken"
    verification:
      - kind: other
        ref: "grep -c 'out.stall_kind' src/calib_dispatch.hpp: 9, one per case arm, each = res.base.stall_kind; r_bridge.cpp copies dres.stall_kind into res_stall_kind at all 9 dispatch call sites (plus the AUTO path via pack_solver_result), which SET_VECTOR_ELT writes to SEXP element 48 as convergence_stall_kind — the exact field R/harvest.R:702 reads"
        status: pass
    human_judgment: false
  - id: D5
    description: "Full DoD gate green after newton_kl migration, including the stepstone benchmark gate"
    verification:
      - kind: unit
        ref: "R CMD INSTALL --preclean . succeeds; full testthat suite 0 FAIL / 1833 PASS (identical to 02-05 baseline)"
        status: pass
      - kind: integration
        ref: "python -m pytest -q (single-thread BLAS): 159 passed, 0 failed (identical to 02-05 baseline)"
        status: pass
      - kind: other
        ref: "LBW_BENCH_GATE=1 NOT_CRAN=true testthat filter=\"bench-gate\": kk1204 gate status=0 iters=10 best_error=-7.376e-14 time=1.5s; 3 PASS / 0 FAIL / 2 SKIP (2 skips are pre-existing missing local report/rds fixtures, unrelated to this migration)"
        status: pass
    human_judgment: false
  - id: D6
    description: "leafblower.h unchanged; no new src/*.cpp translation unit; exactly one commit for the solver"
    verification:
      - kind: other
        ref: "git diff --stat src/leafblower.h empty; ls src/*.cpp | wc -l unchanged at 18; git log shows one feat(02-06) commit touching exactly the three planned files"
        status: pass
    human_judgment: false

duration: ~20min
completed: 2026-08-15
status: complete
---

# Phase 2 Plan 6: Migrate newton_kl onto the shared dispatch table Summary

**Added the `RK_ALG_NEWTON_KL` case arm to `lbw::dispatch_solver` — the last named method still on the per-bridge method-string chain — and confirmed all 9 non-AUTO solver arms populate `stall_kind` for harvest.R's production accelerate heuristic; full DoD gate plus stepstone benchmark green.**

## Performance

- **Duration:** ~20min
- **Completed:** 2026-08-15
- **Tasks:** 2/2 completed
- **Files modified:** 3 (`src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`)

## Accomplishments

- `lbw::dispatch_solver` now covers all 9 non-AUTO algorithm slots: sinkhorn, greg, greenkhorn, logit, chebyshev, raking, oris, oris_soft (plans 02-01/02-03/02-04/02-05), plus **newton_kl** (this plan). Only the `RK_ALG_AUTO` routing logic (plan 07) still bypasses the shared table.
- **Task 1 (migrate newton_kl):** the new `RK_ALG_NEWTON_KL` arm calls `lbw::newton_calibrate`, copies the standard base-field set (mirroring the other 8 arms), and additionally populates the two newton-only diagnostics `n_projected_dims`/`lm_mu_final` (fields `DispatchResult` already declared as scaffolding since plan 02-01's inventory). The arm centralizes the empty-best-weights zero-fill sentinel (length `st.n`) that fires when newton's >5% bounds-violation guard leaves `res.base.best_weights` empty — previously duplicated per-callsite on both bridges.
- **c_api.cpp:** replaced the explicit `RK_ALG_NEWTON_KL` branch's `lbw::newton_calibrate` + `pack_newton_result_c` call with `lbw::dispatch_solver` + `pack_dispatch_result_c` + a reuse of the existing `pack_dispatch_oris_extras_c` (added in plan 02-05). This reuse is exact, not approximate: `DispatchResult`'s default-constructed ORIS-only-field values are byte-identical to what `pack_newton_result_c` hardcoded as the post-newton reset, and the new dispatch arm never touches those fields — so no new "newton extras" packer was needed. `pack_newton_result_c` itself is untouched and still serves the AUTO-fallback branch (out of scope, plan 07).
- **r_bridge.cpp:** replaced the explicit `newton_kl` method-string branch's `lbw::newton_calibrate` + `pack_solver_result` call with `lbw::dispatch_solver` + a full field copy into the existing `res_*` locals — same pattern as the other 8 migrated solvers — plus the two R-only diagnostics `res_n_projected_dims`/`res_lm_mu_final`.
- **Task 2 (full DoD gate + stall_kind production-path check):** verified by reading that all 9 case arms in `lbw::dispatch_solver` assign `out.stall_kind = res.base.stall_kind` (grep confirms 9 arms, 9 assignments — 1:1), and that `r_bridge.cpp` copies `dres.stall_kind` into `res_stall_kind` at every dispatch call site, which is packed into SEXP element 48 as `convergence_stall_kind` — the exact field `harvest.R:702` reads to drive its accelerate/stall-classification heuristic. No silent regression risk on this production-critical diagnostic.
- Full DoD gate green: R testthat 0 FAIL / 1833 PASS (141 WARN, 13 SKIP — identical to 02-05 baseline); Python pytest 159 passed / 0 failed; stepstone `bench-gate` filter 3 PASS / 0 FAIL (kk1204 gate: status=0, 10 iters, best_error=-7.376e-14, 1.5s — no regression).
- `leafblower.h` untouched; `ls src/*.cpp | wc -l` unchanged at 18 (no new translation unit). No shared `finalize_weights` call added to `newton_calib.cpp` (explicitly prohibited by the plan — the unit-mode redistribution gap remains tracked separately at `leafblower-og7d.5`).

## Per-Solver Migration Status (all 9 non-AUTO `rk_algorithm_t` values)

| Solver | Enum | Dispatch arm added | Plan |
|---|---|---|---|
| sinkhorn | `RK_ALG_SINKHORN` | Yes | 02-01 |
| greg | `RK_ALG_GREG` | Yes | 02-03 |
| greenkhorn | `RK_ALG_GREENKHORN` | Yes | 02-03 |
| logit | `RK_ALG_LOGIT` | Yes | 02-03 |
| chebyshev | `RK_ALG_CHEBYSHEV` | Yes | 02-04 |
| raking | `RK_ALG_RAKING` | Yes | 02-04 |
| oris | `RK_ALG_ORIS` | Yes | 02-05 |
| oris_soft | `RK_ALG_ORIS_SOFT` | Yes | 02-05 |
| newton_kl | `RK_ALG_NEWTON_KL` | Yes | **02-06 (this plan)** |
| *(auto)* | `RK_ALG_AUTO` | No — pure routing logic, deferred by design | 07 |

All 9 named methods are migrated. Plan 07 starts from this fully-migrated-except-AUTO state.

## Task Commits

1. **Task 1: Migrate newton_kl** — `097b57b` (feat) — `src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`.
2. **Task 2: Full DoD gate + stall-kind production-path check** — no source edits (gate run + read-based verification only); no commit.

## Files Created/Modified

- `src/calib_dispatch.hpp` — added `#include "newton_calib.hpp"`; added the `RK_ALG_NEWTON_KL` case arm to `dispatch_solver` (copies the standard base-field set, the two newton-only diagnostics `n_projected_dims`/`lm_mu_final`, and applies the empty-best-weights zero-fill sentinel); updated the function's doc comment to state all 9 non-AUTO slots are now covered.
- `src/c_api.cpp` — replaced the `RK_ALG_NEWTON_KL` branch body with `lbw::dispatch_solver` + `pack_dispatch_result_c` + `pack_dispatch_oris_extras_c` (reused, not duplicated). `pack_newton_result_c` left intact, still called by the untouched AUTO-fallback branch.
- `src/r_bridge.cpp` — replaced the `"newton_kl"` method-string branch's inline `lbw::newton_calibrate` + `pack_solver_result` call with `lbw::dispatch_solver(RK_ALG_NEWTON_KL, st, dres)` + a full field copy into the existing `res_*` locals, matching the other 8 migrated branches' pattern.

## Decisions Made

See `key-decisions` in frontmatter — summarized: (1) reused `pack_dispatch_oris_extras_c` for c_api.cpp's newton_kl branch rather than writing a near-duplicate helper, verified field-by-field to be an exact reproduction of the pre-migration reset; (2) `res_n_bounds_clamped` is now set explicitly at the r_bridge.cpp call site rather than implicitly via a template trait — same value, no behavior change; (3) both newton-only R diagnostics populated inside the new dispatch arm; (4) no `finalize_weights` change to `newton_calib.cpp`, per the plan's explicit prohibition.

## Deviations from Plan

None — both tasks' `<action>`/`<verify>`/`<acceptance_criteria>` were followed as written. The plan's `read_first` line numbers for `r_bridge.cpp`'s newton_kl branch and `c_api.cpp`'s `pack_newton_result_c`/`RK_ALG_NEWTON_KL` arm were stale (shifted by plan 02-05's edits); located the actual branches via `grep -n "newton"` instead of the plan's literal line ranges — a Rule 3 navigational correction, not a scope change.

## Issues Encountered

None. No architectural decision needed — the migration recipe from plan 02-03's objective applied directly; the one non-obvious step (preserving c_api.cpp's ORIS-only-field reset for the explicit newton_kl branch) resolved cleanly by reusing plan 02-05's existing helper after confirming field-by-field equivalence.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- All 9 named methods (sinkhorn, greg, greenkhorn, logit, chebyshev, raking, oris, oris_soft, newton_kl) now route through `lbw::dispatch_solver`. Only `RK_ALG_AUTO`'s routing logic remains outside the shared table — both bridges' `method="auto"` branches still call `lbw::oris_solve`/`lbw::raking_solve`/`lbw::newton_calibrate` directly (untouched by this or any prior plan in this phase), exactly as plan 06's objective anticipated.
- Plan 07 (AUTO consolidation) can now route every branch of its severe-skew/moderate-skew/zero-compression decision tree through `lbw::dispatch_solver` calls, since every solver it might select is already migrated.
- Plan 08's `test_single_dispatch_site.py` (single-dispatch-site regression guard) has a fully-migrated table to assert against for every non-AUTO enum value.
- No blockers for plan 02-07.

## Self-Check: PASSED

All claimed files exist (`src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`, this SUMMARY.md) and commit `097b57b` is present in `git log`.
