---
phase: 02-one-engine-not-two
plan: 02
subsystem: build-and-test-infra
tags: [pytest, cmake, regression-guard, documentation, sc2, sc3, sc4]

requires:
  - phase: 02-one-engine-not-two
    plan: 01
    provides: shared dispatch table pattern this plan's guards protect (calib_dispatch.hpp)
provides:
  - "python/leafblower/test_core_sources_sync.py — src/*.cpp vs python/CMakeLists.txt CORE_SOURCES drift now a pytest assertion (SC4)"
  - "python/leafblower/test_finalize_weights_sync.py — 7-solver shared-finalize delegation + single-definition-site guard (SC3)"
  - "R↔Python -O optimization-level asymmetry documented in CLAUDE.md, python/CMakeLists.txt, and test_solver_parity.py header (SC2, leafblower-qzto closed)"
affects: [phase-02-plans-03-through-08]

actuals:
  tokens: 9200
  tasks: 3
  commits: 4

tech-stack:
  added: []
  patterns:
    - "Text-scan pytest modules (no compiled artifact, no BLAS dependence) for build-list and delegation-sync regression guards — rides the existing DoD gate per D-05, no separate script/Makefile target."
    - "RED/GREEN verified via temporary probe files (src/zzz_*.cpp, src/zzz_*.hpp) created, asserted-to-fail, then deleted before commit — same TDD discipline as a compiled-code RED/GREEN cycle, applied to text-scan assertions."

key-files:
  created:
    - python/leafblower/test_core_sources_sync.py
    - python/leafblower/test_finalize_weights_sync.py
  modified:
    - CLAUDE.md
    - python/CMakeLists.txt
    - python/leafblower/test_solver_parity.py

key-decisions:
  - "Did not file a duplicate bd ticket for the newton_kl unit-mode bounds gap. leafblower-og7d.5 already fully documents the exact same gap (same fixture, same measurement, same 'no fix before Phase 2' scheduling) — added a cross-referencing bd comment instead and used og7d.5's own ID in the test source, satisfying the acceptance criterion ('a bd ticket exists... its ID appears in the test source') without ticket sprawl."
  - "CLAUDE.md commit isolated from a pre-existing, unrelated, already-dirty working-tree deletion (a 'bd remember' bullet under Beads Rules, absent since before this plan started). Built the commit content from HEAD + only this plan's one-line edit via a temp-file swap, so the commit carries exactly the -O asymmetry documentation change and nothing else; the pre-existing unrelated deletion remains uncommitted in the working tree, untouched, per the executor's explicit instruction not to stage/restore/touch it."

requirements-completed: []

coverage:
  - id: D1
    description: "A build-list divergence (src/*.cpp vs python/CMakeLists.txt CORE_SOURCES) fails the DoD gate as a pytest assertion naming the drifted file, not as an undefined-symbol link error"
    requirement: US-004
    verification:
      - kind: unit
        ref: "python/leafblower/test_core_sources_sync.py::test_core_sources_matches_src_glob — RED (temp src/zzz_probe.cpp) and GREEN (unmodified tree) both verified"
        status: pass
    human_judgment: false
  - id: D2
    description: "Every 7 bounded solvers still call the shared calib_dispatch.hpp finalize_weights/finalize_weights_buf helper; the helper is defined exactly once; newton_calib.cpp's exclusion is explicit and reasoned in the test source"
    requirement: US-004
    verification:
      - kind: unit
        ref: "python/leafblower/test_finalize_weights_sync.py — 2 tests, both RED (temp probe files) and GREEN (unmodified tree) verified"
        status: pass
    human_judgment: false
  - id: D3
    description: "The R vs Python -O optimization-level asymmetry is documented as a deliberate, bounded decision at CLAUDE.md, python/CMakeLists.txt, and the parity assertion's own file header — no compiler flag changed on either side"
    requirement: US-004
    verification:
      - kind: other
        ref: "grep -c check_make_vars CLAUDE.md python/CMakeLists.txt -> 1 each; git diff on src/Makevars.in and configure empty; python/CMakeLists.txt diff comment-only"
        status: pass
    human_judgment: false
  - id: D4
    description: "Full DoD gate green after all three tasks: R build + testthat, Python reinstall + pytest, single-thread BLAS"
    verification:
      - kind: unit
        ref: "R CMD INSTALL --preclean . succeeded; Rscript -e 'devtools::test()' -> 0 FAIL, 1839 PASS, 11 SKIP"
        status: pass
      - kind: integration
        ref: "cd python && uv pip install -e . --reinstall-package leafblower succeeded; OMP/OPENBLAS/MKL_NUM_THREADS=1 .venv/bin/python -m pytest -> 159 passed, 0 failed"
        status: pass
    human_judgment: false

duration: ~35min
completed: 2026-08-15
status: complete
---

# Phase 2 Plan 2: Guards — SC4/SC3/SC2 Summary

**Closed the three non-dispatch success criteria (SC4 build-list sync, SC3 shared-finalize delegation, SC2 optimization-level asymmetry) as regression-prevention pytest guards and documentation, verifying each with a RED/GREEN probe cycle before committing — zero solver code touched.**

## Performance

- **Duration:** ~35min
- **Completed:** 2026-08-15
- **Tasks:** 3/3 completed
- **Files modified:** 5 (2 new pytest modules, CLAUDE.md, python/CMakeLists.txt, test_solver_parity.py)

## Accomplishments

- **SC4 (build-list sync):** `python/leafblower/test_core_sources_sync.py` asserts `src/*.cpp` (minus `r_bridge.cpp`) equals `python/CMakeLists.txt`'s `CORE_SOURCES` set. Verified RED with a temporary `src/zzz_probe.cpp` (failed, named in the assertion message) and GREEN on the unmodified tree (17 files each side, zero drift).
- **SC3 (shared-finalize delegation):** `python/leafblower/test_finalize_weights_sync.py` asserts (a) all 7 bounded solvers (`oris_finalize`, `raking`, `chebyshev`, `greenkhorn`, `greg`, `logit_calib`, `sinkhorn`) still call `finalize_weights`, and (b) `finalize_weights_buf` is defined exactly once, in `calib_dispatch.hpp`. `newton_calib.cpp`'s exclusion is stated in the test source, not implied by omission, with the reason (smooth-dual Newton, no box-constrained inner step, report-not-clamp T4 contract) and cross-referenced to the pre-existing `leafblower-og7d.5`. Verified RED with two independent temporary probes (a non-delegating solver entry; a second `finalize_weights_buf` definition site) and GREEN on the unmodified tree.
- **SC2 (optimization-level asymmetry):** Documented at all three places a reader can hit it — `CLAUDE.md`'s bullet now scopes the "no `-O`" claim to R and states Python's `-O3` explicitly; a new comment at `python/CMakeLists.txt`'s `-O3` line names the R-side counterpart and the CRAN rule; `test_solver_parity.py`'s file header records the asymmetry and states the tolerances are the bound on how much it may move a result. No compiler flag changed on either side. `leafblower-qzto` closed with the resolution comment.
- Full DoD gate green after all three tasks: R (`R CMD INSTALL --preclean .` succeeded; `devtools::test()` — 0 FAIL, 1839 PASS, 11 SKIP) and Python (`uv pip install -e . --reinstall-package leafblower` succeeded; pytest — 159 passed, 0 failed, single-thread BLAS).

## Task Commits

1. **Task 1: SC4 build-list divergence gate** — `4a44340` (test) — `python/leafblower/test_core_sources_sync.py`.
2. **Task 2: SC3 shared-finalize delegation gate** — `d1db12e` (test) — `python/leafblower/test_finalize_weights_sync.py`.
3. **Task 3: SC2 documentation** — split into two commits by necessity (see Decisions Made): `3225799` (docs, `CLAUDE.md` only, isolated from an unrelated pre-existing working-tree deletion) and `33c6326` (docs, `python/CMakeLists.txt` + `python/leafblower/test_solver_parity.py`).

## Files Created/Modified

- `python/leafblower/test_core_sources_sync.py` (new) — SC4 build-list sync assertion.
- `python/leafblower/test_finalize_weights_sync.py` (new) — SC3 shared-finalize delegation + single-definition-site assertions.
- `CLAUDE.md` — optimization bullet rescoped to name both build sites and both `-O` levels.
- `python/CMakeLists.txt` — comment added at the `-O3` line naming the R-side counterpart and the CRAN constraint; no flag changed.
- `python/leafblower/test_solver_parity.py` — file-header note on the asymmetry and what the parity tolerances bound; no test logic changed.

## Decisions Made

- **No duplicate bd ticket for the `newton_kl` gap.** The plan's Task 2 action instructs filing a new ticket referencing `leafblower-og7d.5`. Direct read of that ticket showed it already fully records the identical gap — same fixture (seed 7/11/23), same measured out-of-bounds magnitudes, same conclusion (T4 report-not-clamp contract is firing as designed; scheduled "no earlier than Phase 2," i.e., now), still open. Filing a second ticket for the same gap would be pure duplication. Added a cross-referencing `bd comment` to `leafblower-og7d.5` instead, linking it to this plan's SC3 test, and used `leafblower-og7d.5`'s own ID in the test source's comment. This satisfies the acceptance criterion's literal wording ("A `bd` ticket exists for the `newton_kl` unit-mode bounds gap and its ID appears in the test source") without ticket sprawl.
- **CLAUDE.md commit isolated from a pre-existing unrelated deletion.** The working tree's `CLAUDE.md` had a "Use `bd remember`..." bullet already missing (uncommitted, present at `git show HEAD:CLAUDE.md`, absent from the working file since before this plan started — confirmed by comparing against the CLAUDE.md content shown to this agent at conversation start, which also lacked it). Committing the file as-is would have swept that unrelated deletion into this plan's commit, violating the executor's instruction not to stage/restore/touch pre-existing dirty-tree state. Built the commit content by applying only this plan's one-line edit on top of `git show HEAD:CLAUDE.md`, staged and committed that, then restored the working-tree file to its prior state (this plan's edit + the pre-existing unrelated deletion, still uncommitted) — so the commit is surgical and the working tree's other dirt is left exactly as found.

## Deviations from Plan

None beyond the two decisions documented above, both scope-preserving (Rule 3: avoided a redundant duplicate ticket / avoided committing unrelated dirty-tree state; neither changes what the plan's success criteria require).

## Issues Encountered

None. `newton_kl`'s existing exclusion was already fully anticipated by CONTEXT.md D-03/D-04 and RESEARCH.md Pitfall 2/Assumption A1 — no new duplication or architectural finding surfaced during Task 2's read (D-04's trigger condition did not fire).

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- SC2, SC3, SC4 are closed for Phase 2. Remaining phase work is SC1 (dispatch unification) across plans 02-03 through 02-08.
- The two new pytest modules ride the existing DoD gate automatically (`.coverage-thresholds.json`'s `enforcement.command`) — no separate script or CI step to remember.
- `leafblower-qzto` is closed. `leafblower-og7d.5` remains open (by design — no fix in Phase 2's scope for the `newton_kl` gap itself; it is tracked, not resolved).
- No blockers for plan 02-03 (migrate `greg`, `greenkhorn`, `logit`).

---
*Phase: 02-one-engine-not-two*
*Completed: 2026-08-15*

## Self-Check: PASSED

All claimed files exist:
- `python/leafblower/test_core_sources_sync.py` — FOUND
- `python/leafblower/test_finalize_weights_sync.py` — FOUND
- `CLAUDE.md`, `python/CMakeLists.txt`, `python/leafblower/test_solver_parity.py` — modified, present

All claimed commits present in `git log`: `4a44340`, `d1db12e`, `3225799`, `33c6326`.
