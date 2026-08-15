---
phase: 02-one-engine-not-two
plan: 01
subsystem: api
tags: [cpp17, dispatch, r-bridge, pybind11, calib_dispatch, sinkhorn]

requires:
  - phase: 01-verification-coverage
    provides: regression net (testthat + pytest suites) the migration is verified against
provides:
  - "lbw::DispatchResult — neutral, ABI-unconstrained internal result struct (calib_dispatch.hpp)"
  - "lbw::dispatch_solver() — shared enum-dispatch table, RK_ALG_SINKHORN case only"
  - "sinkhorn reached from both c_api.cpp::rk_calibrate() and r_bridge.cpp::C_rk_calibrate() through the same dispatch call"
  - "field inventory (rk_result_t vs R's 49-element SEXP result) posted as a bd comment on leafblower-rywn, including one previously-undocumented drift field (aa_accepted_count)"
affects: [phase-02-plans-02-through-07]

actuals:
  tokens: 3533
  tasks: 3
  commits: 1

tech-stack:
  added: []
  patterns:
    - "Shared dispatch table + neutral result struct in calib_dispatch.hpp (mirrors the existing finalize_weights/finalize_weights_buf shared-helper pattern) — the target shape every remaining solver migration (plans 02-07) converges on."
    - "Narrowing copy at each FFI boundary: c_api.cpp marshals DispatchResult -> rk_result_t (ABI-frozen subset); r_bridge.cpp marshals DispatchResult -> the unchanged 49-element SEXP list (full superset)."

key-files:
  created: []
  modified:
    - src/calib_dispatch.hpp
    - src/c_api.cpp
    - src/r_bridge.cpp

key-decisions:
  - "Dropped the `const CellTable&` parameter from dispatch_solver's plan-specified signature — every solver already builds its own CellTable internally from CalibState (verified: sinkhorn.cpp:84, oris.cpp:126, raking.cpp:62); no caller has one available and none is needed for sinkhorn. Documented as a deviation (Rule 3) rather than a silent respec."
  - "Committed Task 2's code change and Task 3's gate verification as ONE commit, per the plan's own Task 3 instruction ('commit the migration as one commit so a later git revert of this solver is clean') — not per-task, for this plan only."
  - "solver_message on DispatchResult is a fixed char[256] (matching r_bridge's res_solver_message and rk_result_t's message field), not std::string — keeps the RAII swap-release surface to the single std::vector<double> best_weights member the plan explicitly flagged."

requirements-completed: [US-004]

coverage:
  - id: D1
    description: "R's harvest(method=\"sinkhorn\") and Python's harvest(method=\"sinkhorn\") reach lbw::sinkhorn_solve through the same lbw::dispatch_solver() call in calib_dispatch.hpp"
    requirement: US-004
    verification:
      - kind: unit
        ref: "tests/testthat/test-sinkhorn-invariants.R (filter=sinkhorn, 1 passed)"
        status: pass
      - kind: integration
        ref: "python -m pytest -k sinkhorn (4 passed: test_solver_parity.py, test_parity_weights.py)"
        status: pass
    human_judgment: false
  - id: D2
    description: "Every field R's harvest() returned for sinkhorn before this plan is still returned, identical values, rtol=1e-6 R-vs-Python"
    requirement: US-004
    verification:
      - kind: unit
        ref: "full R testthat suite (0 FAIL, 1833 PASS)"
        status: pass
      - kind: integration
        ref: "full Python pytest suite (156 passed, 0 failed)"
        status: pass
    human_judgment: false
  - id: D3
    description: "sizeof(rk_params_t)==264 and sizeof(rk_result_t)==536 ABI static_asserts still hold"
    requirement: US-004
    verification:
      - kind: unit
        ref: "R CMD INSTALL --preclean . (static_asserts compile as part of the build; build succeeded)"
        status: pass
    human_judgment: false
  - id: D4
    description: "Task 1 field inventory posted as a bd comment on leafblower-rywn, evidence-driven (quoted from source, none invented)"
    verification:
      - kind: other
        ref: "bd show leafblower-rywn | grep -c n_projected_dims -> 2"
        status: pass
    human_judgment: false

duration: ~25min
completed: 2026-08-15
status: complete
---

# Phase 2 Plan 1: Shared Dispatch Table (Sinkhorn Tracer) Summary

**Built `lbw::DispatchResult` + `lbw::dispatch_solver()` in `calib_dispatch.hpp` and routed sinkhorn through it from both `c_api.cpp` and `r_bridge.cpp`, proving the shared-table pattern (SC1, leafblower-rywn) on one low-risk solver before the remaining eight migrate.**

## Performance

- **Duration:** ~25min (not explicitly timestamped at plan start; estimated from session length)
- **Completed:** 2026-08-15
- **Tasks:** 3/3 completed
- **Files modified:** 3 (`src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`)

## Accomplishments

- One shared `{enum -> solver -> result}` dispatch site now exists in `calib_dispatch.hpp` — no new `.cpp` translation unit, header-only per the no-LTO / `CORE_SOURCES`-sync constraint.
- Both R's `C_rk_calibrate()` and the C API's `rk_calibrate()` reach `lbw::sinkhorn_solve` through the identical `lbw::dispatch_solver()` call instead of two independently hand-maintained code paths.
- Full field-drift inventory (rk_result_t's 43 fields vs. R's 49-element SEXP result) posted as a `bd comment` on `leafblower-rywn`, confirming the 5 previously-known superset-only fields (`n_projected_dims`, `lm_mu_final`, `sraa_demoted`, `convergence_stall_kind`, `best_weights`) **and discovering a 6th undocumented one**: `aa_accepted_count` (ORISResult-only, present in R's SEXP result at element 41, absent from both `rk_result_t` and Python's `result_dict`). Filed as a finding, not fixed in this plan (out of scope — sinkhorn has no `aa_accepted_count`).
- Full Definition-of-Done gate green: R testthat 0 FAIL / 1833 PASS; Python pytest 156 passed / 0 failed; sinkhorn-specific parity (`test_solver_parity.py`, `test_parity_weights.py`) 4/4 passed.

## Task Commits

1. **Task 1: Field inventory** — no commit (no source file modified; inventory posted as a `bd comment` on `leafblower-rywn` per the task's own scope).
2. **Task 2 + Task 3: Shared dispatch table wired for sinkhorn + full DoD gate** — `a8b2e57` (feat) — combined into a single commit per Task 3's explicit instruction ("commit the migration as one commit so a later `git revert` of this solver is clean").

_Note: Tasks 2 and 3 land in one commit by plan design — Task 3 has no source edits of its own, only gate verification gating whether Task 2's change is committed._

## Files Created/Modified

- `src/calib_dispatch.hpp` — added `lbw::DispatchResult` (neutral, ABI-unconstrained result struct) and `lbw::dispatch_solver()` (RK_ALG_SINKHORN case only; every other enum value is currently a no-op). Added `#include "sinkhorn.hpp"`.
- `src/c_api.cpp` — added `#include "calib_dispatch.hpp"`; added `pack_dispatch_result_c()` (narrows `DispatchResult` into the ABI-frozen `rk_result_t`); `RK_ALG_SINKHORN` arm now calls `lbw::dispatch_solver` + the narrowing pack instead of `lbw::sinkhorn_solve` + `pack_solver_result` directly.
- `src/r_bridge.cpp` — added `#include "calib_dispatch.hpp"`; added a function-scope `lbw::DispatchResult dres` local (declared alongside the other `res_*` locals so its heap-backed `best_weights` member is covered by both `Rf_error` swap-release blocks); the `sinkhorn` `strcmp` branch now calls `lbw::dispatch_solver` and copies `dres`'s fields into the existing `res_*` locals — the 49-element SEXP-packing block is untouched.

## Decisions Made

- **Dropped `const CellTable&` from `dispatch_solver`'s signature.** The plan's frontmatter and Task 2 action text both specify `dispatch_solver(rk_algorithm_t, CalibState&, const CellTable&, DispatchResult&)`. Direct read of `sinkhorn.cpp:84`, `oris.cpp:126`, and `raking.cpp:62` confirms every solver already builds its own `CellTable` internally from `CalibState` — no call site in `c_api.cpp` or `r_bridge.cpp` has a `CellTable` available at the dispatch call, and sinkhorn's dispatch doesn't need one externally. Implemented `dispatch_solver(rk_algorithm_t alg, CalibState& st, DispatchResult& out)` without it. Add the parameter in whichever future plan (03-06) first has a genuine consumer for it at the dispatch level.
- **Single commit for Tasks 2+3.** Task 3's own `<action>` text says: "If the R suite and the Python parity suite are both green, commit the migration as one commit so a later `git revert` of this solver is clean (D-01/D-02: git history IS the rollback mechanism)." Held Task 2's code change uncommitted until Task 3's full DoD gate passed, then made one commit — this is a plan-specified exception to the default per-task commit cadence, not a deviation.
- **`solver_message` as `char[256]`, not `std::string`.** Matches `r_bridge.cpp`'s existing `res_solver_message` local and `rk_result_t`'s `message` field exactly, and keeps the RAII swap-release surface to exactly the one heap-backed member (`best_weights`) the plan's Pitfall-4 discussion calls out — no second heap allocation to track through the `Rf_error` unwind blocks.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking issue] `dispatch_solver`'s literal signature is unsatisfiable without dead code**
- **Found during:** Task 2
- **Issue:** The plan's frontmatter artifacts table and Task 2's `<action>` both specify `dispatch_solver(rk_algorithm_t, CalibState&, const CellTable&, DispatchResult&)`. `CellTable` has no default constructor path cheap enough to fabricate at a call site (`build_cell_table` requires real `group_ids`/`cat_counts`/`weights`), and neither `c_api.cpp` nor `r_bridge.cpp` has a `CellTable` available in their `rk_calibrate`/`C_rk_calibrate` bodies today — every solver (sinkhorn included, `sinkhorn.cpp:84`) builds its own internally. Passing the literal signature would force constructing a real `CellTable` at both call sites purely to satisfy an unused reference parameter.
- **Fix:** Implemented `dispatch_solver(rk_algorithm_t alg, CalibState& st, DispatchResult& out)` — dropped the `CellTable` parameter. Documented in the commit message and here.
- **Files modified:** `src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`
- **Verification:** R build + full testthat/pytest suites green with the 3-parameter signature; no caller needed a `CellTable`.
- **Committed in:** `a8b2e57` (part of the task commit)

---

**Total deviations:** 1 auto-fixed (Rule 3).
**Impact on plan:** Signature-only; no behavioral change, no scope creep. The dispatch table's actual purpose — one shared enum-dispatch-and-extraction site both bridges call — is fully realized as specified.

## Issues Encountered

- **New drift finding beyond the plan's named 5 superset-only fields:** direct comparison of every `res_*` local in `r_bridge.cpp` against `rk_result_t`'s field list (`leafblower.h:112-161`) turned up a 6th field absent from the C ABI struct: `aa_accepted_count` (`ORISResult.aa_accepted_count`, `oris.hpp:42`), exposed to R at SEXP result element 41 but never packed by `c_api.cpp`'s `pack_oris_result_c` nor Python's `_bindings.cpp` `result_dict`. This is irrelevant to sinkhorn (which has no `aa_accepted_count`) and out of scope for this plan; recorded on the `leafblower-rywn` `bd comment` per that ticket's own DoD ("Any pre-existing field drift found in step 1 is filed as its own ticket, not fixed here") so whichever future plan migrates `oris`/`oris_soft` inherits the finding instead of rediscovering it.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- The shared-table pattern is proven end-to-end (both bridges, real solver, full DoD gate green) — plans 02-07 can each add one `case` arm to `dispatch_solver` and extend `DispatchResult` with whatever solver-specific fields that solver needs (`aa_accepted_count`, `sor_*`, `alm_*`, `homotopy_*` for oris/oris_soft; nothing extra for raking/greg/greenkhorn/logit/chebyshev beyond what already exists).
- `leafblower-rywn`'s Step 1 field inventory is now posted and complete (including the newly-found `aa_accepted_count` drift) — later migration plans can read it directly instead of re-deriving.
- No blockers for the next solver migration plan.

---
*Phase: 02-one-engine-not-two*
*Completed: 2026-08-15*

## Self-Check: PASSED

All claimed files exist (`src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`, this SUMMARY.md) and commit `a8b2e57` is present in `git log`.
