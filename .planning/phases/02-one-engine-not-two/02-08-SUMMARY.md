---
phase: 02-one-engine-not-two
plan: 08
subsystem: build-and-test-infra
tags: [pytest, r-testthat, regression-guard, sc1, sc5, phase-gate, dispatch]

requires:
  - phase: 02-one-engine-not-two
    plan: 07
    provides: "the fully-consolidated single dispatch site (lbw::dispatch_solver / lbw::route_auto in calib_dispatch.hpp, R bridge's strcmp chain collapsed) this plan's guard test asserts against"
  - phase: 02-one-engine-not-two
    plan: 02
    provides: "the SC2/SC3/SC4 regression guards (test_core_sources_sync.py, test_finalize_weights_sync.py) whose text-scan pattern this plan's SC1 guard copies"
provides:
  - "python/leafblower/test_single_dispatch_site.py — SC1 now enforced by a pytest text-scan assertion, not only satisfied by the current source (RED/GREEN verified with a temporary probe)"
  - "Full Definition-of-Done gate proven green with recorded numbers (SC5): R testthat 0 FAIL/1833 PASS, Python pytest 160 passed/0 failed, single-thread BLAS"
  - "Stepstone benchmark (kk1204, LBW_BENCH_GATE=1) proven unregressed against the 02-07 baseline with measured numbers, not asserted"
  - "leafblower-rywn (P0 dispatch-unification epic) closed with DoD evidence"
affects: []

actuals:
  tokens: 700
  tasks: 2
  commits: 1

tech-stack:
  added: []
  patterns:
    - "Third sibling text-scan pytest guard (after test_core_sources_sync.py/SC4 and test_finalize_weights_sync.py/SC3): no compiled-module import, comment-stripped before assertion, RED/GREEN verified via a temporary probe line inserted into and removed from the scanned source."

key-files:
  created:
    - python/leafblower/test_single_dispatch_site.py
  modified: []

key-decisions:
  - "Task 1's <behavior> text specified 'at most two dispatch_solver call sites (the primary and the AUTO fallback)'. Direct read of the post-plan-07 src/r_bridge.cpp shows 3 non-comment call sites: the AUTO primary (line 621), the AUTO newton_kl fallback (line 678), and the one unified explicit-method call (line 739) -- the plan's own Task 2 acceptance criteria (Rule 1 elsewhere) and 02-07-SUMMARY.md's own D4 coverage entry both independently confirm this same count of 3. The plan's stated bound omitted the explicit-method branch's call. Implemented and committed the test asserting <= 3, the count that actually matches the consolidated architecture, and documented the correction in the test's own docstring and in the task commit message -- asserting the plan's literal '<=2' would make Task 1's own acceptance criterion ('the new test passes on the post-plan-07 tree') impossible to satisfy against real code."

requirements-completed: [US-004]

coverage:
  - id: D1
    description: "SC1 defended by a pytest guard: zero strcmp(method_str, ...) occurrences (comment-excluded) and a bounded number of lbw::dispatch_solver call sites in src/r_bridge.cpp"
    requirement: US-004
    verification:
      - kind: unit
        ref: "python/leafblower/test_single_dispatch_site.py::test_r_bridge_has_no_method_dispatch_chain -- RED (temporary strcmp(method_str, \"raking\") probe line: FAILED as expected) and GREEN (unmodified tree: PASSED) both verified; comment-exclusion also verified (same probe moved into a // comment: PASSED)"
        status: pass
    human_judgment: false
  - id: D2
    description: "Full Definition-of-Done gate green with single-thread BLAS: R build, R testthat, Python reinstall, Python pytest"
    requirement: US-004
    verification:
      - kind: unit
        ref: "R CMD INSTALL --preclean . succeeded; OMP/OPENBLAS/MKL_NUM_THREADS=1 Rscript devtools test_dir: 0 FAIL, 1833 PASS, 141 WARN, 13 SKIP"
        status: pass
      - kind: integration
        ref: "cd python && uv pip install -e . --reinstall-package leafblower succeeded; OMP/OPENBLAS/MKL_NUM_THREADS=1 .venv/bin/python -m pytest -q: 160 passed, 0 failed (159 pre-plan baseline + this plan's 1 new guard)"
        status: pass
    human_judgment: false
  - id: D3
    description: "Stepstone benchmark (kk1204, LBW_BENCH_GATE=1) shows no regression against the pre-plan (02-07) baseline"
    requirement: US-004
    verification:
      - kind: other
        ref: "LBW_BENCH_GATE=1 NOT_CRAN=true testthat filter=\"bench-gate\": 0 FAIL, 3 PASS, 2 SKIP (pre-existing missing local report/rds fixtures, unrelated); kk1204 gate: status=0 iters=10 best_error=-7.376e-14 time=1.5s -- byte-identical to 02-07-SUMMARY.md's recorded baseline (status=0 iters=10 best_error=-7.376e-14 time=1.5s), expected since this plan changed no C++ source"
        status: pass
    human_judgment: false
  - id: D4
    description: "Human confirms the phase's user-visible surface (all 9 solvers + AUTO, R-only fields, Python result dict) is unchanged"
    verification:
      - kind: manual_procedural
        ref: "harvest() run for all 9 explicit methods + auto in R (identical attr(result) field set incl. all 4 R-only fields + aa_accepted_count) and the same 9 methods in Python (auto correctly rejected as R-only per _harvest.py:381-383, unchanged pre-existing contract; identical 36-key result dicts, zero R-only-field leakage); stepstone delta (byte-identical to 02-07 baseline) judged acceptable"
        status: pass
    human_judgment: true
    rationale: "Requires running harvest() interactively in R and Python and comparing result shapes/field presence by eye against the pre-phase baseline, plus a subjective judgment on benchmark noise -- the plan's Task 3 explicitly gates this on human sign-off (checkpoint:human-verify, gate=\"blocking\"). User approved 2026-08-15 after independent corroboration (test-newton-kl.R + test-cr-d5-auto-fallback-fields.R: 0 FAIL/29 PASS; test_solver_parity.py + test_parity_weights.py: 21 passed/0 failed)."

duration: ~30min
completed: 2026-08-15
status: complete
---

# Phase 2 Plan 8: Phase gate — SC1 test enforcement + SC5 proof Summary

**SC1 is now defended by a pytest regression guard (not just satisfied by the current source); SC5 is proven with recorded numbers — full DoD gate green (R 0 FAIL/1833 PASS, Python 160 passed/0 failed) and the stepstone benchmark byte-identical to the pre-plan baseline. Human sign-off on the unchanged user-visible surface obtained — phase 2 (One Engine, Not Two) is complete.**

## Performance

- **Duration:** ~30min
- **Completed:** 2026-08-15
- **Tasks:** 3/3 completed
- **Files modified:** 1 (`python/leafblower/test_single_dispatch_site.py`, new)

## Accomplishments

- **Task 1 — SC1 guard:** `python/leafblower/test_single_dispatch_site.py::test_r_bridge_has_no_method_dispatch_chain` reads `src/r_bridge.cpp`, strips `//`-comment lines, and asserts (a) zero occurrences of `strcmp(method_str` in the remainder and (b) at most 3 `lbw::dispatch_solver(` call sites. RED/GREEN verified with a temporary probe: inserting `strcmp(method_str, "raking") == 0` as live code made the test FAIL with the exact drift-explaining message; moving the identical text into a `//` comment made it pass again (comment-exclusion confirmed); the probe was then fully removed and `git status --porcelain src/` confirmed empty before commit. The test imports no compiled module — pure text scan, same shape as its SC3/SC4 siblings from plan 02.
- **Task 2 — SC5 proof:** Ran the exact `.coverage-thresholds.json` `enforcement.command` sequence with single-thread BLAS exported throughout: `R CMD INSTALL --preclean .` succeeded; R testthat — 0 FAIL, 1833 PASS, 141 WARN, 13 SKIP (matches the established 02-04/02-05 baseline exactly); Python reinstall succeeded; Python pytest — 160 passed, 0 failed (159 baseline + this plan's 1 new guard test). Then ran the stepstone benchmark gate (`LBW_BENCH_GATE=1 NOT_CRAN=true`, `filter="bench-gate"`): kk1204 gate `status=0`, `iters=10`, `best_error=-7.376e-14`, `time=1.5s` — byte-identical to 02-07-SUMMARY.md's recorded post-AUTO-consolidation baseline, confirming zero regression (expected, since this plan added a Python-only test file and touched no C++ source). No gate step was skipped.
- `leafblower-rywn` (the P0 dispatch-unification epic tracking SC1 since plan 01) closed with a comment enumerating evidence against its own 9-item Definition of Done — all satisfied except one literal sub-clause (see Deviations below).
- **Task 3 — human sign-off:** the checkpoint was surfaced (per this plan's `autonomous: false`) and NOT auto-approved on an agent's say-so — a first message claiming approval explicitly disclosed it came from the orchestrator running the checks itself, not the user, and was rejected on that basis alone (no agent message substitutes for the user's own approval on a `gate="blocking"` checkpoint). Independently re-ran the field-presence evidence (`test-newton-kl.R` + `test-cr-d5-auto-fallback-fields.R`: 0 FAIL/29 PASS; `test_solver_parity.py` + `test_parity_weights.py`: 21 passed/0 failed) before presenting it. The user then approved directly. R's `harvest()` across all 9 explicit methods + `auto` returns an identical `attr(result)` field set including all 4 R-only fields (`n_projected_dims`, `lm_mu_final`, `sraa_demoted`, `convergence_stall_kind`) plus `aa_accepted_count`; Python's `harvest()` across the same 9 methods (`auto` correctly rejected — R-only, per `_harvest.py:381-383`, a pre-existing contract this phase did not touch) returns identical 36-key result dicts with zero leakage of the R-only fields, confirming the C ABI struct stayed frozen. Stepstone delta (byte-identical to the 02-07 baseline) judged acceptable.

## Task Commits

1. **Task 1: SC1 guard — a second dispatch chain cannot come back unnoticed** — `3f6d4f1` (test) — `python/leafblower/test_single_dispatch_site.py`.
2. **Task 2: SC5 — full DoD gate and stepstone no-regression** — no source edits (gate run only, per the plan's own `<files>` spec); evidence recorded above and on `leafblower-rywn`'s closing comment.
3. **Task 3: Human confirmation that the user-visible surface is unchanged** — `6d93d8d` (docs, Tasks 1-2 metadata) + this plan-metadata commit (checkpoint resolution) — no source edits, verification-only.

**Plan metadata:** `6d93d8d` (Tasks 1-2 interim) + final commit below (Task 3 resolution + phase close).

## Files Created/Modified

- `python/leafblower/test_single_dispatch_site.py` (new) — SC1 single-dispatch-site regression guard.

## Decisions Made

See `key-decisions` in frontmatter: the plan's Task 1 `<behavior>` text stated an incorrect bound ("at most two" `dispatch_solver` call sites); the actual, architecturally-correct count is 3 (AUTO primary + AUTO newton_kl fallback + the one unified explicit-method call), independently confirmed by 02-07-SUMMARY.md's own D4 coverage entry. Implemented the test against the real, correct invariant rather than the plan's miscounted one, since asserting the plan's literal number would make the guard itself always fail on the very tree it exists to protect.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Corrected the plan's stated `dispatch_solver` call-site bound from 2 to 3**
- **Found during:** Task 1
- **Issue:** The plan's `<behavior>` said "at most two `dispatch_solver` call sites (the primary and the AUTO fallback)". Grepping the post-plan-07 `src/r_bridge.cpp` (`grep -nF "lbw::dispatch_solver("`) shows 3 non-comment call sites, not 2 — the plan's count omitted the non-AUTO explicit-method branch's own single call (line 739), counting only the two calls inside the `if (alg == RK_ALG_AUTO)` branch (lines 621, 678).
- **Fix:** Wrote the assertion against the real, verified count (`<= 3`), with the rationale documented in the test's docstring and inline assert message, citing 02-07-SUMMARY.md's D4 coverage entry as the independent cross-check.
- **Files modified:** `python/leafblower/test_single_dispatch_site.py`.
- **Verification:** Test passes GREEN on the unmodified post-plan-07 tree (required by Task 1's own acceptance criteria); would have been an immediately-failing, un-committable test had the plan's literal "2" been used.
- **Committed in:** `3f6d4f1` (part of the Task 1 commit).

**2. [Rule 3 - out of literal scope, low-risk] `leafblower-rywn`'s own DoD item "Any pre-existing field drift found in step 1 is filed as its OWN ticket" is not literally satisfied**
- **Found during:** closing `leafblower-rywn` after Task 2's DoD proof.
- **Issue:** Plan 02-01 found a 6th superset-only field (`aa_accepted_count`) beyond the 5 already named in `leafblower-rywn`/RESEARCH.md, and — per 02-01-SUMMARY.md's own recorded decision — documented it on `leafblower-rywn`'s own comment thread rather than filing a new, separate `bd` ticket. This satisfies the *intent* (traceable, not silently dropped) but not the DoD's literal wording ("its OWN ticket").
- **Fix:** None taken — this was plan 02-01's decision, out of this plan's task scope, and is explicitly not a blocker to SC1 (dispatch unification, which is what `leafblower-rywn` tracks) — `aa_accepted_count`'s C-API/Python omission is a separate, pre-existing field-parity gap unrelated to whether R and Python share one dispatch site. Recorded honestly in `leafblower-rywn`'s closing comment rather than silently claimed as fully compliant.
- **Files modified:** none (bd ticket only).
- **Committed in:** not applicable (bd state, not a git commit — `.beads/issues.jsonl`'s update is committed alongside this SUMMARY).

---

**Total deviations:** 2 (1 Rule 1 bug-fix in a test assertion; 1 Rule 3 scope note, no code change).
**Impact on plan:** Both preserve the plan's actual intent (a correctly-enforced SC1 guard; an honestly-closed epic). No scope creep — no solver code touched, no new architecture.

## Issues Encountered

None beyond the plan's own Task 1 miscount (see Deviations). The DoD gate, R build, and stepstone benchmark all ran clean on the first attempt with no debugging required.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- SC1 and SC5 both have durable evidence in the repository: SC1 via `test_single_dispatch_site.py` (rides the existing DoD gate automatically), SC5 via the recorded DoD-gate and stepstone numbers above.
- **Task 3 resolved — user approved 2026-08-15** after a first "approved" claim was rejected for coming from the orchestrator's own verification run rather than the user (see Accomplishments). Phase 2 (`One Engine, Not Two`) is complete: all 5 success criteria (SC1-SC5) have durable artifacts or recorded measurements — no phase success criterion rests on an unverified assertion.
- Phase 3 (`Honest Performance Gate`, per `ROADMAP.md`) is next. `leafblower-2ouc`'s `benchmarks/` infrastructure is the reuse target already noted in STATE.md's carried-forward decisions.

## Self-Check: PASSED

All claimed files exist:
- `python/leafblower/test_single_dispatch_site.py` — FOUND
- this SUMMARY.md — FOUND

All claimed commits present in `git log`: `3f6d4f1`, `6d93d8d`.
