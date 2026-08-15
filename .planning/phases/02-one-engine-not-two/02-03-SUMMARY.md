---
phase: 02-one-engine-not-two
plan: 03
subsystem: api
tags: [cpp17, dispatch, r-bridge, greg, greenkhorn, logit, calib_dispatch]

requires:
  - phase: 02-one-engine-not-two
    plan: 01
    provides: "lbw::DispatchResult / lbw::dispatch_solver shared table (sinkhorn tracer) — this plan extends it"
provides:
  - "RK_ALG_GREG, RK_ALG_GREENKHORN, RK_ALG_LOGIT case arms in lbw::dispatch_solver (calib_dispatch.hpp)"
  - "greg, greenkhorn, logit reached from both c_api.cpp::rk_calibrate() and r_bridge.cpp::C_rk_calibrate() through the same dispatch call"
affects: [phase-02-plans-04-through-08]

actuals:
  tokens: 4775
  tasks: 3
  commits: 3

tech-stack:
  added: []
  patterns:
    - "Extended plan 02-01's shared-dispatch-table pattern to three more solvers, one case arm each, one commit each (D-01 solver-by-solver granularity)."
    - "Best-iterate solvers (greenkhorn, logit) keep their own existing copy-into-output step at each FFI boundary unchanged — dispatch_solver only moves best_weights out of the solver result; it does not unify the two callers' copy logic (explicit non-goal per the plan's Task 2 action text)."

key-files:
  created: []
  modified:
    - src/calib_dispatch.hpp
    - src/c_api.cpp
    - src/r_bridge.cpp

key-decisions:
  - "Ran the full R testthat suite (not testthat's filter=\"greg\"/\"greenkhorn\" narrowing) for Task 1 and Task 2 verification. testthat's test_dir(filter=...) matches only file basenames (stripped of test- and .R), and no file in tests/testthat/ has \"greg\" as a filename substring (test-calibration-solvers.R etc. carry greg tests but don't match) — running the plan's literal filter=\"greg\" command reproducibly errors with \"No test files found.\" (verified). Filed as a deviation (Rule 3 — the plan's own verify command doesn't run); the full suite already gates on every FAIL, so this is verification-equivalent, not weaker."
  - "The two stale 'mirrors r_bridge.cpp:806-810' line-number comments in c_api.cpp's greenkhorn/logit arms (already pointing at sinkhorn's block, not the real centralized copy-back, before this plan started) were replaced with a description of the actual target — r_bridge.cpp's centralized post-dispatch-chain GREENKHORN/LOGIT weights copy-back block — instead of a fresh numeric range that would go stale the next time either file's line count shifts."

requirements-completed: [US-004]

coverage:
  - id: D1
    description: "R and Python both reach lbw::greg_solve, lbw::greenkhorn_solve and the logit solver through lbw::dispatch_solver"
    requirement: US-004
    verification:
      - kind: unit
        ref: "tests/testthat full suite (0 FAIL, 1833 PASS) run after each of the 3 task commits"
        status: pass
      - kind: integration
        ref: "python/leafblower/test_solver_parity.py -k 'greg or greenkhorn or logit' (4 passed, both sides converged)"
        status: pass
    human_judgment: false
  - id: D2
    description: "R's harvest() result for greg, greenkhorn and logit has the same element names and values as before the migration (rtol=1e-6 R-vs-Python)"
    requirement: US-004
    verification:
      - kind: unit
        ref: "full R testthat suite green after logit migration (0 FAIL, 1833 PASS, 141 WARN, 13 SKIP — same counts pre/post migration)"
        status: pass
      - kind: integration
        ref: "full Python pytest suite (159 passed, 0 failed)"
        status: pass
    human_judgment: false
  - id: D3
    description: "greenkhorn's and logit's best-iterate weights still land in the returned weight vector on both bridges"
    requirement: US-004
    verification:
      - kind: unit
        ref: "greenkhorn/logit both remain best-iterate solvers post-migration: dispatch_solver moves best_weights out; c_api.cpp's own copy-into-weights[] step and r_bridge.cpp's centralized post-dispatch GREENKHORN/LOGIT copy-back block are both unchanged and still fire (test_greenkhorn_parity, test_logit_parity, test_logit_default_rule_parity all pass with real, non-zero weights)"
        status: pass
    human_judgment: false
  - id: D4
    description: "Rf_error unwinding still releases every heap-backed local on both the error and the success path after each branch is rewired"
    verification:
      - kind: other
        ref: "Both r_bridge.cpp swap-release blocks (error path, success path) already covered dres.best_weights as a function-scope local before this plan (plan 02-01); greg/greenkhorn/logit reuse the same function-scope dres, no new heap-backed member added, both blocks re-read unchanged"
        status: pass
    human_judgment: false
  - id: D5
    description: "Input validation still runs before dispatch and is unchanged"
    verification:
      - kind: other
        ref: "greenkhorn's min_weight<max_weight R-side guard (throw before dispatch_solver call) and validate_calibrate_inputs (c_api.cpp, runs before alg dispatch) both untouched — grep confirms no diff to either check"
        status: pass
    human_judgment: false

duration: ~30min
completed: 2026-08-15
status: complete
---

# Phase 2 Plan 3: Migrate greg, greenkhorn, logit onto the shared dispatch table Summary

**Added `RK_ALG_GREG`, `RK_ALG_GREENKHORN` and `RK_ALG_LOGIT` case arms to `lbw::dispatch_solver`, routing all three solvers through the shared table from both `c_api.cpp` and `r_bridge.cpp` — one commit per solver, full DoD gate green after all three.**

## Performance

- **Duration:** ~30min
- **Completed:** 2026-08-15
- **Tasks:** 3/3 completed
- **Files modified:** 3 (`src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`)

## Accomplishments

- `lbw::dispatch_solver` now covers 4 of 12 algorithm slots (sinkhorn from plan 02-01, plus greg/greenkhorn/logit from this plan) — the shared `{enum -> solver -> result}` table grows by one case arm per solver, matching D-01's solver-by-solver granularity exactly.
- **greg** (no superset-only field, writes `st.weights` in-place): migrated arm returns `dres.status` directly in `c_api.cpp` and mirrors sinkhorn's shape in `r_bridge.cpp`, preserving the pre-existing asymmetry where greg does NOT fall through to the generic message-synthesis tail (unlike greenkhorn/logit) — this was read and preserved, not "fixed," since AUTO never routes to greg/greenkhorn/logit so the skipped auto-fallback check was always dead code for them.
- **greenkhorn** (best-iterate solver): `dispatch_solver` moves `best_weights` out of the solver result; both FFI boundaries keep their own pre-existing copy-into-output step (`c_api.cpp`'s copy into the flat `weights[]` buffer; `r_bridge.cpp`'s centralized post-dispatch-chain copy-back for `GREENKHORN`/`LOGIT`) completely unchanged, per the plan's explicit instruction not to unify them in this task.
- **logit** (also best-iterate): identical shape to greenkhorn's migration; `grake_norm` convergence metric confirmed still flowing through `DispatchResult.grake_norm` (already a pre-existing member from the sinkhorn tracer, needed no new field).
- The `strcmp` chain in `C_rk_calibrate` is three branches shorter in substance (still three `strcmp` arms, but each is now a ~28-line copy-from-`dres` block instead of a direct solver call + `pack_solver_result`) — matches the plan's `<verification>` claim.
- `leafblower.h` untouched; `ls src/*.cpp | wc -l` unchanged at 18 (no new translation unit).
- Full DoD gate green after Task 3: R testthat 0 FAIL / 1833 PASS (141 WARN, 13 SKIP — unchanged counts from before this plan); Python pytest 159 passed / 0 failed; `test_solver_parity.py -k 'greg or greenkhorn or logit'` 4/4 passed, both sides converged.

## Task Commits

1. **Task 1: Migrate greg** — `39b99d3` (feat) — `src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`.
2. **Task 2: Migrate greenkhorn** — `c055f8d` (feat) — same 3 files.
3. **Task 3: Migrate logit, then run the full DoD gate** — `6ce0ee9` (feat) — same 3 files.

## Files Created/Modified

- `src/calib_dispatch.hpp` — added `#include "greg.hpp"`, `#include "greenkhorn.hpp"`, `#include "logit_calib.hpp"`; added `RK_ALG_GREG`, `RK_ALG_GREENKHORN`, `RK_ALG_LOGIT` case arms to `lbw::dispatch_solver` (same 24-field copy shape as the existing `RK_ALG_SINKHORN` arm — no new `DispatchResult` member needed, since none of the three solvers carries a field beyond what sinkhorn already required); updated the function's doc comment to list current coverage.
- `src/c_api.cpp` — `RK_ALG_GREG` arm now calls `lbw::dispatch_solver` + `pack_dispatch_result_c` and returns directly (matching the original early-return shape); `RK_ALG_GREENKHORN`/`RK_ALG_LOGIT` arms now call `lbw::dispatch_solver` into a block-scoped `DispatchResult` (`dres_gk`/`dres_lg`), keep their pre-existing best_weights-into-`weights[]` copy, and still fall through to the shared post-switch tail (auto-fallback check, generic message synthesis) exactly as before. The two stale "mirrors r_bridge.cpp:806-810" comments were rewritten to describe the actual target block instead of a numeric range.
- `src/r_bridge.cpp` — `"greg"`, `"greenkhorn"`, `"logit"` `strcmp` branches each replaced with a `lbw::dispatch_solver(RK_ALG_*, st, dres)` call + a copy into the existing `res_*` locals, mirroring sinkhorn's already-migrated branch shape exactly. greenkhorn's `min_weight < max_weight` R-side validation guard stays exactly where it was (before the dispatch call, unchanged). The function-scope `dres` local and both `Rf_error` swap-release blocks (error path, success path) were re-read and required no change — they already cover `dres.best_weights` from plan 02-01.

## Decisions Made

- **Full-suite verification instead of the plan's literal `filter="greg"`/`filter="greenkhorn"` testthat commands.** Direct execution of `test_dir("tests/testthat", filter="greg", stop_on_failure=TRUE)` fails with `Error: No test files found.` — testthat's `filter` matches only file basenames (after stripping `test-`/`.R`), and no file in `tests/testthat/` has `"greg"` as a filename substring (greg-covering tests live in `test-calibration-solvers.R`, `test-clamp-contract.R`, etc., none named with "greg"). Verified this reproduces on the unmodified plan-provided command before treating it as a Rule 3 blocking issue. Ran the full `test_dir("tests/testthat", stop_on_failure=TRUE)` suite for Task 1 and Task 2 verification instead — it gates on every FAIL including greg/greenkhorn-specific tests, so this is verification-equivalent to (in fact strictly broader than) what the plan's filter was intended to check, not a weakening. `filter="greenkhorn"` DOES match (test-greenkhorn-best-metric.R, test-greenkhorn-exit.R, test-greenkhorn-l1-weight.R, test-stall-kind-greenkhorn.R) and was not the blocking case, but the full suite was run for both tasks for consistency and because Task 3's own verify command already runs the full suite unconditionally.
- **Rewrote rather than renumbered the stale "mirrors r_bridge.cpp:806-810" comments.** Both instances (c_api.cpp's greenkhorn and logit arms) already pointed at the wrong location before this plan started — line 806-810 in the current r_bridge.cpp sits inside the (already-migrated, plan 02-01) sinkhorn branch, not anywhere near greenkhorn/logit's real copy-back logic (which is a single centralized block after the entire dispatch chain, not a per-branch block at any fixed line range that would map 1:1 to a per-branch comment in c_api.cpp). Since this plan's own edits shift every subsequent line number in both files, writing a fresh numeric range would go stale immediately upon plan 04's edits to c_api.cpp/r_bridge.cpp. Replaced with a plain description of the real target block instead — this doesn't reintroduce the staleness the migration recipe's step 5 exists to prevent.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking issue] Plan's testthat filter verify commands for Task 1/Task 2 match zero files**
- **Found during:** Task 1, running the plan's literal `<verify>` command for greg.
- **Issue:** `test_dir("tests/testthat", filter="greg", stop_on_failure=TRUE)` errors with `No test files found.` because testthat's `filter` argument matches file basenames, not file contents or test descriptions, and no file is named with a "greg" substring.
- **Fix:** Ran the full `test_dir("tests/testthat", stop_on_failure=TRUE)` suite instead for both Task 1 and Task 2, confirming 0 FAIL / 1833 PASS each time (same counts as the pre-migration baseline in 02-02-SUMMARY.md).
- **Files modified:** none (verification-only; no source change from this finding).
- **Verification:** Full suite green after both task commits; Task 3's own `<verify>` command (which already runs the full suite unconditionally) confirms this again after logit's migration.
- **Committed in:** No separate commit — reflected in the verification evidence for `39b99d3` and `c055f8d`.

---

**Total deviations:** 1 auto-fixed (Rule 3). No architectural changes, no scope creep — the plan's `must_haves.artifacts` and `<verification>` claims are all met as specified; only the literal test-filter *command text* in two `<verify>` blocks turned out not to exercise any file, which the full-suite run supersedes without weakening coverage.

## Issues Encountered

None beyond the testthat filter finding above. No new field drift, no new architectural decision needed — greg/greenkhorn/logit's result structs (`GregResult`, `GreenkornResult`, `LogitCalibResult`) all fit within the `DispatchResult` shape plan 02-01 already established for sinkhorn; no `M_cell` field surfacing was needed since neither bridge exposed `GregResult::M_cell`/`LogitCalibResult::M_cell` before this plan (confirmed via grep — out of scope, unchanged).

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- 4 of 12 algorithm slots now route through `lbw::dispatch_solver` (sinkhorn, greg, greenkhorn, logit). Plans 04-06 migrate the field-bearing solvers (oris/oris_soft/newton_kl/chebyshev/raking) against this now-twice-proven pattern.
- The `RK_ALG_GREENKHORN`/`RK_ALG_LOGIT` copy-into-output duplication the plan's Task 2 action explicitly deferred ("Do NOT unify or delete either copy in this task even though they now read from the same struct member; that is a follow-on simplification") remains open for whichever future plan wants to pursue it — not blocking, not regressed.
- No blockers for plan 02-04.

---
*Phase: 02-one-engine-not-two*
*Completed: 2026-08-15*

## Self-Check: PASSED

All claimed files exist (`src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`, this SUMMARY.md) and commits `39b99d3`, `c055f8d`, `6ce0ee9` are present in `git log`.
