---
phase: 02-one-engine-not-two
plan: 07
subsystem: api
tags: [cpp17, dispatch, r-bridge, c-api, auto-routing, calib_dispatch]

requires:
  - phase: 02-one-engine-not-two
    plan: 06
    provides: "lbw::dispatch_solver covering all 9 non-AUTO rk_algorithm_t values — this plan adds the final RK_ALG_AUTO routing decision and collapses the R bridge's per-method branching that called dispatch_solver 9 separate times into one call site"
provides:
  - "lbw::route_auto() (calib_dispatch.hpp) — the single AUTO routing DECISION (compression ratio, K, target_skew thresholds), shared by both bridges; never solves, only decides"
  - "lbw::resolve_m_cell_est() — the single lazily-cached lbw::estimate_M_cell wrapper both bridges' AUTO routing and oris_soft capacity_penalty auto-resolution share"
  - "lbw::kAlgNames / lbw::kAlgNamesLen — the single 12-entry enum-to-name table (reserved-slot holes at 2/7, two static_asserts), replacing both bridges' independently-drifted copies"
  - "C_rk_calibrate (r_bridge.cpp) now dispatches purely off the already-resolved/validated rk_algorithm_t enum through ONE lbw::dispatch_solver() call for every explicit method, plus route_auto()'s up-to-two calls for AUTO — no strcmp(method_str, ...) branching left"
  - "st.use_admm_capacity moved into dispatch_solver's RK_ALG_ORIS_SOFT arm (both bridges set it identically, no caller context needed)"
affects: [phase-02-plan-08-single-dispatch-site-test]

actuals:
  tokens: 17596
  tasks: 3
  commits: 3

tech-stack:
  added: []
  patterns:
    - "route_auto() takes raw (n, K, group_ids, cat_counts, targets) rather than a CalibState, because c_api.cpp's algorithm-resolution switch runs BEFORE its CalibState is built (r_bridge.cpp's CalibState already exists by the time it calls route_auto and just passes st.n/st.K/st.group_ids/st.cat_counts/st.targets through) — the two bridges' call sites differ in build order, not in what they pass."
    - "A comprehensive dres.* -> res_* field copy (all ~40 DispatchResult fields, unconditional) after a single lbw::dispatch_solver() call is behaviorally identical to per-algorithm-selective copying, because DispatchResult's default-constructed values for fields an algorithm doesn't touch are BYTE-IDENTICAL to the res_* locals' own top-of-function declared defaults — confirmed field-by-field for every one of the 9 non-AUTO algorithms before relying on it (Task 3)."
    - "State that depends only on the resolved enum value (st.use_admm_capacity for RK_ALG_ORIS_SOFT) moves into dispatch_solver's case arm; state that depends on WHY that enum was chosen (st.oris_auto_selected — true only when AUTO itself picked ORIS, per types.hpp's documented contract) stays on the caller side, because dispatch_solver's (alg, st, out) signature has no way to see the caller's routing history."

key-files:
  created: []
  modified:
    - src/calib_dispatch.hpp
    - src/c_api.cpp
    - src/r_bridge.cpp

key-decisions:
  - "AUTO-implementation diff table (Task 1, required even if identical): two genuine, both behaviorally-inert, divergences were found while tabulating r_bridge.cpp's vs c_api.cpp's pre-migration inline AUTO routing before writing route_auto(). (1) st.oris_auto_selected: c_api.cpp set it true unconditionally for ANY AUTO pick (including newton_kl/raking); r_bridge.cpp only set it true when AUTO actually selected ORIS, matching types.hpp's own documented contract ('true iff AUTO routing selected ORIS'). Only oris_solve()'s verbose '[AUTO->ORIS] ' log prefix reads this field — no result-field impact. (2) accelerate restore on auto-fallback: r_bridge.cpp's severe-skew branch forces st.accelerate=true for ORIS+SRAA then restores the user's original value (CR-D5) before a newton_kl fallback solve; c_api.cpp's fallback never restored it. newton_calib.cpp never reads st.accelerate, confirmed inert (st goes out of scope with no further reads). Both fixed by adopting R's fixture-pinned behavior in both bridges, per the plan's explicit instruction; recorded on leafblower-uqyf (filed and closed same session, since the fix already landed in the Task 1 commit)."
  - "kAlgNames observability finding (Task 2, required even if identical): grepped R/Python test suites and R/harvest.R for the two tables' differing spellings. r_bridge.cpp's lowercase-canonical names (matching kAlgMap's keys) feed a STRUCTURED, tested SEXP result field (res_list element 3, R's algorithm_used string) — these cannot change and are what the unified table uses. c_api.cpp's different spellings ('ORIS'/'ORIS-soft'/'(reserved)'/'(gap)') only ever fed rk_result_t's default message free text; nothing greps for that capitalization, so c_api.cpp converges onto the canonical names — a message-text-only change with zero test impact."
  - "Literal grep-count acceptance criteria (both files) return non-zero because both bridges correctly REFERENCE the single lbw::kAlgNames/lbw::dispatch_solver symbols from calib_dispatch.hpp — the criteria's intent (no LOCAL duplicate definition, no redundant per-branch call) is met; see Coverage D2/D4 below for the exact counts and why they are non-zero."
  - "Both Rf_error swap-release blocks (r_bridge.cpp:831-839 error path, :875-881 success path) re-read after Task 3's collapse: no new heap-backed function-scope local was introduced. weights_backup/accel_backup (AUTO branch) are try-block-scoped, not function-scope — normal C++ stack unwinding destroys them on a thrown exception before the deferred Rf_error() call, and they are simply out of scope (no leak risk) on the success path. The pre-existing covered set (pre_error, solver_error, group_ids, cat_counts, tgt_storage, targets, weights, res_best_weights, dres.best_weights) is unchanged and complete on both paths."

requirements-completed: [US-004]

coverage:
  - id: D1
    description: "lbw::route_auto() is the single AUTO routing decision (compression/skew thresholds preserved exactly); both bridges call it instead of maintaining independent inline copies"
    requirement: US-004
    verification:
      - kind: unit
        ref: "tests/testthat filter=\"auto\": 0 FAIL, 30 PASS; test-cr-d5-auto-fallback-fields.R: 0 FAIL, 14 PASS"
        status: pass
      - kind: integration
        ref: "python -m pytest -q: 159 passed, 0 failed (AUTO path exercised via c_api.cpp's rk_calibrate through the Python bindings)"
        status: pass
    human_judgment: false
  - id: D2
    description: "route_auto exists once in calib_dispatch.hpp; both bridges reference it (grep -c 'route_auto' >= 1 in all three files) — literal grep counts are 3/3/3 (1 usage line + comments mentioning it by name), not the raw definition count, since the plan's own acceptance criterion phrasing counts any occurrence"
    verification:
      - kind: other
        ref: "grep -c 'route_auto' src/calib_dispatch.hpp = 3 (1 definition + 2 comment mentions); src/c_api.cpp = 3 (1 call site line 281 + 2 comment mentions); src/r_bridge.cpp = 3 (1 call site line 621 + 2 comment mentions) — all >= 1, satisfying the acceptance criterion's literal '>= 1' wording"
        status: pass
    human_judgment: false
  - id: D3
    description: "One 12-entry enum-to-name table (lbw::kAlgNames, calib_dispatch.hpp) with reserved-slot holes at 2/7 and both static_asserts; both bridges' own hand-maintained copies deleted"
    requirement: US-004
    verification:
      - kind: other
        ref: "grep -A15 'inline constexpr const char\\* kAlgNames' src/calib_dispatch.hpp shows 12 entries, slots 2 and 7 empty, both static_asserts present; grep for a LOCAL 'static const char* kAlgNames' or 'static constexpr const char* kAlgNames' definition in src/r_bridge.cpp and src/c_api.cpp returns zero matches (both files' own copies were deleted, not just aliased)"
        status: pass
      - kind: unit
        ref: "full R testthat suite: 0 FAIL, 0 ERROR (141 WARN, unchanged from 02-06 baseline)"
        status: pass
    human_judgment: false
  - id: D4
    description: "R bridge's method-string chain collapsed to one lbw::dispatch_solver() call per explicit method (plus route_auto's up to two for AUTO); zero strcmp(method_str, ...) branches remain"
    requirement: US-004
    verification:
      - kind: other
        ref: "grep -v '^ *//' src/r_bridge.cpp | grep -c 'strcmp(method_str' = 0; actual lbw::dispatch_solver(...) CALL SITES in r_bridge.cpp = 3 (line 621 AUTO primary, line 678 AUTO fallback, line 739 the one unified explicit-method call) — literal 'grep -c dispatch_solver' returns 13 because 10 of those lines are prose comments explaining the design, not call sites"
        status: pass
      - kind: unit
        ref: "full R testthat suite: 0 FAIL, 0 ERROR; tests/testthat/test-bridge-length-checks.R (CR-D10 unrecognized-method-string test): 0 FAIL, 21 PASS — the exact error message and behavior survived the collapse unchanged"
        status: pass
      - kind: integration
        ref: "python -m pytest -q: 159 passed, 0 failed"
        status: pass
      - kind: other
        ref: "LBW_BENCH_GATE=1 NOT_CRAN=true testthat filter=\"bench-gate\": kk1204 gate status=0, iters=10, best_error=-7.376e-14, time=1.5s — byte-identical to the 02-06 baseline, confirming zero regression from the collapse"
        status: pass
    human_judgment: false
  - id: D5
    description: "leafblower.h unchanged; no new src/*.cpp translation unit"
    verification:
      - kind: other
        ref: "git diff --stat 2718671 HEAD -- src/leafblower.h: empty (0 lines changed); ls src/*.cpp | wc -l: 18 (unchanged from 02-06 baseline)"
        status: pass
    human_judgment: false
---

# Phase 2 Plan 7: Consolidate AUTO routing and unify the enum-to-name table Summary

**Added `lbw::route_auto()` and one `lbw::kAlgNames` table to `calib_dispatch.hpp`, then collapsed R's 9-branch `strcmp(method_str, ...)` chain into a single `lbw::dispatch_solver()` call — SC1's "adding a routing rule requires editing ONE dispatch site" now holds for AUTO routing too, closing the last gap `leafblower-rywn` opened.**

## Performance

- **Duration:** ~50min
- **Completed:** 2026-08-15
- **Tasks:** 3/3 completed
- **Files modified:** 3 (`src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`)

## Accomplishments

- **Task 1 — `lbw::route_auto()`:** the single AUTO routing DECISION (K/compression-ratio/target_skew thresholds, preserved verbatim — exact-integer `>=` 0.9 comparison, `K>=5`, `target_skew>5`) now lives once in `calib_dispatch.hpp`, returning the chosen `rk_algorithm_t` plus a `force_accelerate` flag; the caller then dispatches via `lbw::dispatch_solver`, which is what makes the auto-fallback path a *second* `dispatch_solver` call instead of nested routing logic. `lbw::resolve_m_cell_est()` factors out the lazily-cached `estimate_M_cell` wrapper both `route_auto` and both bridges' `oris_soft` `capacity_penalty` auto-resolution now share (eb79.15), so the O(n·K) estimate runs at most once per call regardless of how many callers touch it. Both bridges' inline AUTO copies (r_bridge.cpp's 135-line branch, c_api.cpp's 60-line switch case) are deleted; two now-orphaned R-bridge lambdas (`pack_solver_result`/`pack_oris_result`) and their `has_message`/`has_n_bounds` trait detectors were removed as a consequence, and c_api.cpp's `pack_newton_result_c` was retired (its fallback now uses a fresh `DispatchResult` whose default-constructed ORIS-extra fields reproduce the old hardcoded reset exactly, per plan 02-06's established pattern).
- **Task 2 — one enum-to-name table:** `lbw::kAlgNames`/`lbw::kAlgNamesLen`, 12 positional entries with the documented holes at slot 2 (removed LBFGSB) and slot 7 (withdrawn GRAKE), guarded by the two `static_assert`s. Both bridges' independently-drifted copies (r_bridge.cpp's lowercase canonical spellings vs c_api.cpp's capitalized display spellings) are deleted. Grepped R/Python test suites and `R/harvest.R` before choosing which spelling wins: r_bridge.cpp's names feed a *structured* SEXP result field (`algorithm_used`), c_api.cpp's only ever fed a free-text default message nothing tests — so the shared table uses r_bridge.cpp's canonical names, and c_api.cpp's default message text changes cosmetically (e.g. `"ORIS: 5 iters..."` → `"oris: 5 iters..."`) with zero test impact.
- **Task 3 — collapse the R bridge's dispatch chain:** discovered that `p.algorithm` is *already* resolved and validated against `kAlgMap` upfront in `C_rk_calibrate` (CR-D10), before the `strcmp` chain even runs — so the chain was pure redundant re-derivation. Replaced it with `const rk_algorithm_t alg = p.algorithm;` and one unified `if (alg == RK_ALG_AUTO) { ... } else { lbw::dispatch_solver(alg, st, dres); <one comprehensive field copy> }`, collapsing 8 near-identical ~35-line branches (sinkhorn/greg/greenkhorn/logit/chebyshev/newton_kl/oris_soft/oris) into one. Verified field-by-field before trusting the "copy everything unconditionally" design: `DispatchResult`'s default-constructed values for fields an algorithm doesn't touch are byte-identical to the `res_*` locals' own top-of-function defaults, so an unconditional copy reproduces every algorithm's exact pre-migration output. `st.use_admm_capacity` (identical on both bridges, no caller context needed) moved into `dispatch_solver`'s `RK_ALG_ORIS_SOFT` arm; `st.oris_auto_selected` stays caller-side since it depends on *why* `RK_ALG_ORIS` was chosen (AUTO vs. explicit), information `dispatch_solver`'s signature can't see. All 3 remaining "mirrors r_bridge.cpp" comments in c_api.cpp retired.
- Full DoD gate green throughout: R testthat 0 FAIL / 0 ERROR (141 WARN, unchanged from 02-06 baseline); Python pytest 159 passed / 0 failed; `LBW_BENCH_GATE=1` kk1204 gate `status=0`, `best_error=-7.376e-14`, `1.5s` — byte-identical to the 02-06 baseline, confirming zero performance regression from the collapse.
- `leafblower.h` untouched; `ls src/*.cpp | wc -l` unchanged at 18 (no new translation unit).

## AUTO-Implementation Diff Table (Task 1, required per acceptance criteria)

| Aspect | r_bridge.cpp (pre-migration) | c_api.cpp (pre-migration) | Divergence? | Resolution |
|---|---|---|---|---|
| 0.9 compression threshold | `M_cell_est*10 >= n*9` (exact int) | `M_cell_est*10 >= n*9` (exact int) | No | Preserved verbatim in `route_auto` |
| K>=5 gate | `K >= 5` | `K >= 5` | No | Preserved verbatim |
| target_skew formula | `max_target / max(min_target, 1e-12)` | identical | No | Preserved verbatim |
| Severe-skew threshold | `target_skew > 5.0` | `target_skew > 5.0` | No | Preserved verbatim |
| `M_cell_est` caching | lazy, shared with `oris_soft` capacity site via `m_cell_est_cache` | recomputed independently, no cache | Computation-redundancy only (Pitfall 3), not a value divergence | `route_auto`'s in/out `m_cell_est_cache` param unifies both onto the lazy-cache policy; c_api.cpp's `oris_soft` capacity site also switched to `resolve_m_cell_est` for symmetry (inert there — AUTO never selects `oris_soft`) |
| `st.oris_auto_selected` on AUTO | `true` only when AUTO selected ORIS (severe-skew or compressed branch) | `true` unconditionally for ANY AUTO pick (including newton_kl/raking) | **Yes — genuine, behaviorally inert** | Adopted R's precise contract (matches `types.hpp`'s own doc comment) in both bridges; only affects `oris_solve()`'s verbose log prefix, no result field |
| `accelerate` restore on fallback | `st.accelerate = accel_backup` before newton_kl fallback (CR-D5) | never restored | **Yes — genuine, behaviorally inert** | Adopted R's behavior; `newton_calib.cpp` never reads `st.accelerate`, confirmed inert |
| Fallback trigger | `NOCONV \|\| BUDGET` | `NOCONV \|\| BUDGET` | No | Preserved verbatim |
| Fallback algorithm | always `newton_kl` | always `newton_kl` | No | Preserved verbatim |

Both genuine divergences recorded and closed on `leafblower-uqyf` per the plan's instruction ("file a bd ticket recording the... change so the divergence is a tracked fix, not a silent one").

## Task Commits

1. **Task 1: One AUTO routing decision** — `cac0cfb` (feat) — `src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`.
2. **Task 2: One enum-to-name table** — `8843806` (feat) — `src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`.
3. **Task 3: Collapse the R bridge's method-string chain** — `5938d08` (feat) — `src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`.

## Files Created/Modified

- `src/calib_dispatch.hpp` — added `lbw::resolve_m_cell_est()`, `lbw::AutoRouteResult`/`lbw::route_auto()`, `lbw::kAlgNames`/`lbw::kAlgNamesLen`; moved `st.use_admm_capacity = true;` into the `RK_ALG_ORIS_SOFT` case arm.
- `src/c_api.cpp` — algorithm-resolution switch's `RK_ALG_AUTO` case now calls `lbw::route_auto`; the auto-fallback block now dispatches through `lbw::dispatch_solver` into a fresh `DispatchResult` instead of `lbw::newton_calibrate` + the now-deleted `pack_newton_result_c`; `st.oris_auto_selected` fixed to the precise per-algorithm contract; `accel_backup` now captured and restored on fallback; own `kAlgNames` table deleted (reads `lbw::kAlgNames`); `oris_soft`'s capacity resolution reads `lbw::resolve_m_cell_est`; all 3 "mirrors r_bridge.cpp" comments retired.
- `src/r_bridge.cpp` — `"auto"` branch now calls `lbw::route_auto` + `lbw::dispatch_solver` (primary + fallback) instead of the 135-line inline AUTO copy; the entire remaining 8-branch `strcmp(method_str, ...)` chain collapsed into `if (alg == RK_ALG_AUTO) {...} else { lbw::dispatch_solver(alg, st, dres); <one field copy> }`; `pack_solver_result`/`pack_oris_result` lambdas and `has_message`/`has_n_bounds` traits removed (orphaned by Task 1); own `kAlgNames` table deleted (reads `lbw::kAlgNames`); site-A `capacity_penalty` resolution reads `lbw::resolve_m_cell_est`.

## Decisions Made

See `key-decisions` in frontmatter — summarized: (1) `route_auto` takes raw scalars/pointers rather than a `CalibState`, since c_api.cpp resolves the algorithm before building its `CalibState`; (2) the "copy every `DispatchResult` field unconditionally" design for Task 3's unified branch is safe because non-applicable fields are byte-identical defaults on both sides — verified field-by-field, not assumed; (3) `st.use_admm_capacity` moves into the shared dispatch arm (no caller context needed), `st.oris_auto_selected` does not (depends on caller-side routing history the shared arm's signature can't see).

## Deviations from Plan

None beyond what the plan itself anticipated and scoped (the AUTO-implementation diff table, the observability grep, and the caller-vs-shared-arm split for `use_admm_capacity`/`oris_auto_selected` were all explicitly called for in the plan's task actions). Two Rule-1 bug fixes were made per the plan's own explicit instruction to "adopt the R behaviour... and file a bd ticket" for genuine AUTO-implementation divergences — see the diff table above and `leafblower-uqyf` (filed and closed in this session, since the fix landed in the same commit that found it).

## Issues Encountered

None. The largest risk (Task 3's field-copy unification silently dropping a per-algorithm value) was retired by reading each `dispatch_solver` case arm and each removed branch's field list side-by-side before writing the unified copy, not by assuming it would work — confirmed correct by the full R suite staying at exactly the same PASS/WARN counts as the 02-06 baseline (no new WARN, no FAIL) and the bench-gate producing byte-identical numbers.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- SC1 ("adding a solver, a result field, or a routing rule requires editing ONE dispatch site") is now fully satisfied: all 9 named methods (plan 01-06) plus AUTO routing (this plan) route through `lbw::dispatch_solver`/`lbw::route_auto`; the R bridge has zero per-method branching left.
- Plan 08's `test_single_dispatch_site.py` (single-dispatch-site regression guard) has a fully-consolidated table + routing function to assert against, including the AUTO path.
- No blockers for plan 02-08.

## Self-Check: PASSED

All claimed files exist (`src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`, this SUMMARY.md) and commits `cac0cfb`, `8843806`, `5938d08` are present in `git log`.
