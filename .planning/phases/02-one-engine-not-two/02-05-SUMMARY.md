---
phase: 02-one-engine-not-two
plan: 05
subsystem: api
tags: [cpp17, dispatch, r-bridge, oris, oris-soft, alm, sor, calib_dispatch]

requires:
  - phase: 02-one-engine-not-two
    plan: 04
    provides: "lbw::DispatchResult / lbw::dispatch_solver shared table extended to chebyshev, raking — this plan extends it further"
provides:
  - "RK_ALG_ORIS, RK_ALG_ORIS_SOFT case arms in lbw::dispatch_solver (calib_dispatch.hpp), including the full oris-family diagnostic surface: SOR (9 fields), ALM (4 fields), homotopy/SRAA-alpha/aa_accepted_count/sraa_demoted"
  - "oris, oris_soft reached from both c_api.cpp::rk_calibrate() and r_bridge.cpp::C_rk_calibrate() through the same dispatch call"
affects: [phase-02-plans-06-through-08]

actuals:
  tokens: 6766
  tasks: 2
  commits: 2

tech-stack:
  added: []
  patterns:
    - "Extended plan 02-01's shared-dispatch-table pattern to the largest superset-only diagnostic surface in the package (SOR/ALM/homotopy/SRAA-alpha oris-family fields)."
    - "New c_api.cpp helper pack_dispatch_oris_extras_c: a second, narrower pack function called ONLY for RK_ALG_ORIS/RK_ALG_ORIS_SOFT, immediately after the shared pack_dispatch_result_c — keeps the shared pack function's field set unchanged for every other already-migrated solver (sinkhorn/greg/greenkhorn/logit/raking/chebyshev), which must keep leaving rk_result_init's memset/1.0 defaults untouched for these ORIS-only ABI fields."

key-files:
  created: []
  modified:
    - src/calib_dispatch.hpp
    - src/c_api.cpp
    - src/r_bridge.cpp

key-decisions:
  - "The plan's Task 1 action text enumerated only a subset of the oris-family DispatchResult members to add (alm_*, aa_accepted_count, and 7 of the 9 sor_* fields, omitting sor_min_omega/sor_n_damped and all of n_xcur_writes_per_iter_last/min_alpha_seen/final_alpha/homotopy_levels_used/homotopy_final_factor/greedy_sweeps_taken/eta_final). Grepping the R SEXP-list construction (res_list indices 5-13, 20-21) confirmed EVERY one of these omitted fields is exported into R's 49-element result list today via pack_oris_result. Adding only the enumerated subset would have left the omitted fields at DispatchResult's post-migration defaults for method='oris'/'oris_soft' — a silent correctness regression, directly contradicting the plan's own must_haves truth ('the SOR diagnostics ... with unchanged values', unqualified = all 9, not 7). Added the full set (Rule 2: auto-add missing critical functionality) rather than the plan's literal enumerated list."
  - "rk_result_t (leafblower.h) already declares every oris-family field EXCEPT aa_accepted_count (confirmed via grep — sor_min_omega, sor_n_damped, alm_capacity_mu_final, n_xcur_writes_per_iter_last, min_alpha_seen, final_alpha, homotopy_levels_used, homotopy_final_factor, greedy_sweeps_taken, eta_final all exist in the frozen ABI). pack_oris_result_c (the pre-migration C-side function) wrote all of them directly. Extending the SHARED pack_dispatch_result_c to also copy these fields would have overwritten rk_result_init's defaults (0.0/0, not the oris-only 1.0/0 pattern) for sinkhorn/greg/greenkhorn/logit/raking/chebyshev — a regression for six already-migrated solvers, out of this plan's scope. Instead added a second, oris-only helper (pack_dispatch_oris_extras_c) called only from the RK_ALG_ORIS / RK_ALG_ORIS_SOFT branches, preserving byte-identical rk_result_t output for every other solver."
  - "aa_accepted_count is genuinely superset-only: absent from leafblower.h, and the pre-migration pack_oris_result_c never copied it either (confirmed by reading its full body before deletion). Not a migration-introduced gap — matches the pre-existing 6th superset-only field flagged in plan 02-01 (bd comment on leafblower-rywn). Left unfixed, as directed (out of scope for this plan)."
  - "pack_oris_result vs pack_oris_result_c field diff (Task 1's required acceptance-criterion comparison): no VALUE divergence found. The only differences are fields the frozen rk_result_t ABI structurally lacks — stall_kind (absent from leafblower.h entirely, confirmed by grep), aa_accepted_count, sraa_demoted, and message/solver_message (ORISResult carries no message field; pack_oris_result_c never synthesized one, matching the R lambda's has_message-trait-false '\\0' behavior). Every field BOTH functions carry (SOR ×9, ALM ×4, homotopy ×4, n_bounds ×2, all convergence/quality-metric fields) is copied with an identical source expression on both sides. Adopted the R behaviour per the task's instruction, which was moot here since there was nothing to reconcile."
  - "capacity_penalty / estimate_M_cell comparison (Task 2's required deliverable for plan 07): both bridges resolve capacity_penalty identically — same formula (M_cell_est/n via lbw::estimate_M_cell, falling back to 1.0 when n<=0) on the same inputs, so an explicit oris_soft request produces the SAME auto-resolved capacity_mu on R and C/Python. Neither bridge double-computes estimate_M_cell for an explicit oris_soft request: c_api.cpp's RK_ALG_ORIS_SOFT branch has exactly one call site (the AUTO-routing block never routes to ORIS_SOFT, so its own estimate_M_cell call is never reached on the same invocation); r_bridge.cpp computes it in one unconditional pre-dispatch block (m_cell_est_cache) that already runs for EVERY method call (not just oris_soft — a pre-existing, out-of-scope waste for non-capacity-consuming methods, per RESEARCH.md Pitfall 3), and the AUTO-branch's second site reuses the cached value rather than recomputing. RESEARCH.md Pitfall 3's genuine double-compute risk only exists on r_bridge.cpp's method='auto' path, which this plan's must_haves truth (single evaluation per solve) is trivially satisfied by for oris_soft, and which this plan does not touch — the finding is recorded here for plan 07's consolidation, not fixed."
  - "Removed c_api.cpp's pack_oris_result_c entirely once both of its call sites (RK_ALG_ORIS in Task 1, RK_ALG_ORIS_SOFT in Task 2) were migrated within this same plan — it is a fully dead static function I orphaned by removing its last callers, not a pre-existing orphan, so cleanup falls within this plan's scope per the project's 'clean up only own newly created orphans' rule."

requirements-completed: [US-004]

coverage:
  - id: D1
    description: "R and Python both reach lbw::oris_solve for method oris and method oris_soft through lbw::dispatch_solver"
    requirement: US-004
    verification:
      - kind: unit
        ref: "tests/testthat filter=\"oris\" (11 oris-family test files): 0 FAIL, 138 PASS, 12 WARN (pre-existing); full suite run after both task commits"
        status: pass
      - kind: integration
        ref: "python/leafblower/ pytest -k oris (9 passed) and full pytest suite after Task 2 (159 passed, 0 failed)"
        status: pass
    human_judgment: false
  - id: D2
    description: "R's oris and oris_soft results still report sraa_demoted, aa_accepted_count, all 9 SOR diagnostics, and the four ALM diagnostics with unchanged values; oris_auto_selected/use_admm_capacity assignments stay caller-side"
    requirement: US-004
    verification:
      - kind: unit
        ref: "full R testthat suite green (0 FAIL, 1833 PASS, 141 WARN, 13 SKIP — identical counts to the 02-04 baseline) after both migrations; res_* locals feeding the 49-element R result list re-read unchanged at the R-list-packing tail (res_list indices 5-13, 20-28, 37-41, 47)"
        status: pass
      - kind: other
        ref: "dispatch_solver's RK_ALG_ORIS / RK_ALG_ORIS_SOFT arms copy all 9 sor_* fields, all 4 alm_* fields, aa_accepted_count, sraa_demoted, and the 7 homotopy/alpha/n_xcur fields the plan's Task 1 action text under-enumerated (see key-decisions) — grepped against pack_oris_result's full field list to confirm 1:1 coverage"
        status: pass
    human_judgment: false
  - id: D3
    description: "The M_cell estimate is computed at most once per solve on either bridge for an explicit oris_soft request; Rf_error unwinding still releases every heap-backed local on both the error and success path"
    verification:
      - kind: other
        ref: "c_api.cpp's RK_ALG_ORIS_SOFT branch has one estimate_M_cell call site per invocation (AUTO never routes to ORIS_SOFT); r_bridge.cpp's m_cell_est_cache guard (unchanged, not touched by this plan) already memoizes across its two call sites"
        status: pass
      - kind: other
        ref: "No new heap-backed member added to DispatchResult (all 20 new fields are scalar int/double/bool); both existing Rf_error swap-release blocks (r_bridge.cpp) already cover dres.best_weights unconditionally and require no new swap line"
        status: pass
    human_judgment: false
  - id: D4
    description: "rk_result_t stays 536 bytes; leafblower.h unchanged; no new src/*.cpp translation unit; oris.cpp/oris_finalize.cpp/oris_trajectory.cpp not modified"
    verification:
      - kind: other
        ref: "git diff --stat src/leafblower.h empty across both commits; ls src/*.cpp | wc -l unchanged at 18; git diff --stat src/oris.cpp src/oris_finalize.cpp src/oris_trajectory.cpp empty"
        status: pass
    human_judgment: false
  - id: D5
    description: "Full DoD gate green after oris_soft migration, including the stepstone benchmark gate"
    verification:
      - kind: unit
        ref: "R CMD INSTALL --preclean . succeeds; full testthat suite 0 FAIL / 1833 PASS (identical to 02-04 baseline)"
        status: pass
      - kind: integration
        ref: "python -m pytest -q (single-thread BLAS): 159 passed, 0 failed (identical to 02-04 baseline)"
        status: pass
      - kind: other
        ref: "LBW_BENCH_GATE=1 NOT_CRAN=true testthat filter=\"bench-gate\": kk1204 gate status=0 iters=10 best_error=-7.4e-14 time=1.5s; 3 PASS / 0 FAIL / 2 SKIP (2 skips are pre-existing missing local report/rds fixtures, unrelated to this migration)"
        status: pass
    human_judgment: false

duration: ~35min
completed: 2026-08-15
status: complete
---

# Phase 2 Plan 5: Migrate oris and oris_soft onto the shared dispatch table Summary

**Added `RK_ALG_ORIS` and `RK_ALG_ORIS_SOFT` case arms to `lbw::dispatch_solver` — the package's largest superset-only diagnostic surface (9 SOR fields, 4 ALM fields, 4 homotopy/alpha fields, `aa_accepted_count`, `sraa_demoted`) — one commit per solver, full DoD gate plus stepstone benchmark green after both.**

## Performance

- **Duration:** ~35min
- **Completed:** 2026-08-15
- **Tasks:** 2/2 completed
- **Files modified:** 3 (`src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`)

## Accomplishments

- `lbw::dispatch_solver` now covers 8 of 12 algorithm slots (sinkhorn, greg, greenkhorn, logit, chebyshev, raking from plans 02-01/02-03/02-04, plus oris and oris_soft from this plan). Only `newton_kl` and the AUTO branch still bypass the shared table.
- **oris (Task 1)**: `RK_ALG_ORIS` case arm absorbs what r_bridge.cpp's `pack_oris_result` lambda and c_api.cpp's `pack_oris_result_c` did today. Diffed the two field-by-field per the task's acceptance criterion — no value divergence found (see `key-decisions`); the only differences are fields the frozen `rk_result_t` ABI structurally lacks (`stall_kind`, `aa_accepted_count`, `sraa_demoted`, message).
- **oris_soft (Task 2)**: `RK_ALG_ORIS_SOFT` case arm is identical to `RK_ALG_ORIS` plus the four ALM diagnostics, which only this path populates. `st.use_admm_capacity`, `st.oris_auto_selected`, and the `capacity_penalty` auto-resolution (`estimate_M_cell`) all stay on the caller side, as directed.
- The plan's Task 1 action text enumerated only a partial DispatchResult field list (7 of 9 SOR fields, omitting `sor_min_omega`/`sor_n_damped` and all 7 of `n_xcur_writes_per_iter_last`/`min_alpha_seen`/`final_alpha`/`homotopy_levels_used`/`homotopy_final_factor`/`greedy_sweeps_taken`/`eta_final`). Grepped the R SEXP-list construction to confirm every omitted field is genuinely exported into R's 49-element result today; added the FULL set (Rule 2 — auto-add missing critical functionality) rather than the plan's literal enumeration, since the omission would have silently regressed those fields to stale defaults for `method="oris"`/`"oris_soft"`.
- `rk_result_t` (leafblower.h) already carries every oris-family field except `aa_accepted_count`. Rather than extending the SHARED `pack_dispatch_result_c` (which would have overwritten `rk_result_init`'s defaults for the six already-migrated non-oris solvers — a regression outside this plan's scope), added a second, oris-only helper `pack_dispatch_oris_extras_c`, called only from the two oris branches immediately after `pack_dispatch_result_c`.
- `aa_accepted_count` remains genuinely superset-only (absent from `leafblower.h`; the pre-migration `pack_oris_result_c` never surfaced it either) — a pre-existing gap flagged in plan 02-01, not fixed here (out of scope).
- Removed c_api.cpp's `pack_oris_result_c` entirely: both of its call sites were migrated across this plan's two tasks, leaving it fully dead — cleaned up as an orphan this plan itself created, not a pre-existing one.
- `capacity_penalty`/`estimate_M_cell` comparison for plan 07 (both bridges resolve to the same value via the same formula on the same inputs; neither double-computes for an explicit `oris_soft` request; RESEARCH.md Pitfall 3's genuine double-compute risk lives entirely in r_bridge.cpp's untouched `method="auto"` path) — full finding recorded in `key-decisions`.
- `src/oris.cpp`, `src/oris_finalize.cpp`, `src/oris_trajectory.cpp` untouched — no per-iteration hot-loop code moved across a TU boundary (no LTO in this build).
- `leafblower.h` untouched; `ls src/*.cpp | wc -l` unchanged at 18 (no new translation unit).
- Full DoD gate green after Task 2: R testthat 0 FAIL / 1833 PASS (141 WARN, 13 SKIP — identical to the 02-04 baseline); Python pytest 159 passed / 0 failed; stepstone `bench-gate` filter 3 PASS / 0 FAIL (kk1204 gate: status=0, 10 iters, best_error=-7.4e-14, 1.5s — no regression).

## Task Commits

1. **Task 1: Migrate oris (the default method)** — `8e4d23a` (feat) — `src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`.
2. **Task 2: Migrate oris_soft, then run the full DoD gate** — `75b7236` (feat) — `src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`.

## Files Created/Modified

- `src/calib_dispatch.hpp` — extended `DispatchResult` with the full oris-family diagnostic set (20 new scalar members: `n_xcur_writes_per_iter_last`, `min_alpha_seen`, `final_alpha`, `homotopy_levels_used`, `homotopy_final_factor`, `greedy_sweeps_taken`, `eta_final`, 9× `sor_*`, `aa_accepted_count`, 4× `alm_*`); added `RK_ALG_ORIS` and `RK_ALG_ORIS_SOFT` case arms (each calls `lbw::oris_solve`, copies the ~40-field block, `RK_ALG_ORIS_SOFT` additionally copies the 4 ALM fields); updated the function's doc comment to list current coverage (8 of 12 slots).
- `src/c_api.cpp` — removed the now-dead `pack_oris_result_c`; added `pack_dispatch_oris_extras_c` (copies the oris-only ABI fields `pack_dispatch_result_c` deliberately skips); rewired the default/`RK_ALG_ORIS` branch and the `RK_ALG_ORIS_SOFT` branch to call `lbw::dispatch_solver` + `pack_dispatch_result_c` + `pack_dispatch_oris_extras_c`, preserving the fall-through into the AUTO-fallback tail (unchanged control-flow shape — oris/oris_soft still participate in the `p->algorithm == RK_ALG_AUTO` fallback gate check, a no-op for explicit calls). Capacity-penalty resolution and `use_admm_capacity` assignment for `oris_soft` unchanged.
- `src/r_bridge.cpp` — replaced the default/oris branch's inline `lbw::oris_solve` + `pack_oris_result` call with `lbw::dispatch_solver(RK_ALG_ORIS, st, dres)` + a full field copy into the existing `res_*` locals; replaced the `"oris_soft"` branch's inline call with `lbw::dispatch_solver(RK_ALG_ORIS_SOFT, st, dres)` + the same field copy plus the 4 ALM fields. `st.oris_auto_selected` / `st.use_admm_capacity` assignments and the `m_cell_est_cache` capacity-resolution block are all unchanged, still running before dispatch.

## Decisions Made

See `key-decisions` in frontmatter — summarized: (1) completed the plan's under-enumerated DispatchResult field list to the full set actually exported to R (Rule 2); (2) added a dedicated `pack_dispatch_oris_extras_c` rather than extending the shared pack function, to avoid regressing the six other already-migrated solvers; (3) `aa_accepted_count` stays superset-only, matching pre-migration behavior, not fixed; (4) no `pack_oris_result` vs `pack_oris_result_c` value divergence found — differences are purely ABI-structural (fields the frozen `rk_result_t` never carried); (5) `capacity_penalty`/`estimate_M_cell` resolves identically on both bridges and is evaluated at most once per explicit `oris_soft` solve; (6) removed the now-fully-dead `pack_oris_result_c`.

## Deviations from Plan

**[Rule 2 — auto-add missing critical functionality] Completed the DispatchResult field list beyond the plan's literal enumeration.** Task 1's action text listed only `alm_capacity_mu_final`, `alm_n_growth_events`, `alm_max_dual_norm`, `alm_sum_drift`, `aa_accepted_count`, and 7 of the 9 `sor_*` fields (`sor_omega_mean`, `sor_any_latched`, `sor_n_pinned_fb`, `sor_n_warmup_fb`, `sor_n_conv_fb`, `sor_n_resid_grew`, `sor_n_monotone_cd`) — omitting `sor_min_omega`, `sor_n_damped`, and all 7 of `n_xcur_writes_per_iter_last`/`min_alpha_seen`/`final_alpha`/`homotopy_levels_used`/`homotopy_final_factor`/`greedy_sweeps_taken`/`eta_final`. Grepping `r_bridge.cpp`'s SEXP-list construction (`res_list` indices 5-13, 20-21) confirmed every one of these omitted fields is exported into R's 49-element result list today. Implementing only the enumerated subset would have left these fields at stale defaults for `method="oris"`/`"oris_soft"` — a silent correctness regression directly contradicting the plan's own must_haves truth ("the SOR diagnostics ... with unchanged values", which is unqualified and therefore means all 9, not 7). Added the complete set. Files: `src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`. Commits: `8e4d23a`, `75b7236`.

No other deviations — both tasks' `<action>`/`<verify>`/`<acceptance_criteria>` were otherwise followed exactly.

## Issues Encountered

None beyond the deviation above. No architectural decision needed — `ORISResult` fit within `DispatchResult`'s existing shape once extended; no `oris.cpp`/`oris_finalize.cpp`/`oris_trajectory.cpp` change was needed or made.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- 8 of 12 algorithm slots now route through `lbw::dispatch_solver` (sinkhorn, greg, greenkhorn, logit, chebyshev, raking, oris, oris_soft). Plan 06 migrates `newton_kl`, the last field-bearing solver; the AUTO branch (r_bridge.cpp's `strcmp(method_str, "auto")` block and c_api.cpp's `RK_ALG_AUTO` routing) remains unmigrated by design — it is pure routing logic layered on top of the now-shared per-solver arms, deferred to a later plan per this plan's own objective statement.
- The `capacity_penalty`/`estimate_M_cell` finding (both bridges agree, at-most-once evaluation for an explicit `oris_soft` request; the double-compute risk lives entirely in r_bridge.cpp's untouched `method="auto"` pre-dispatch block) is the input plan 07 needs to consolidate the estimate into a single evaluation without introducing a double-compute or an eager compute on paths that never needed it.
- No blockers for plan 02-06.

## Self-Check: PASSED

All claimed files exist (`src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`, this SUMMARY.md) and commits `8e4d23a`, `75b7236` are present in `git log`.
