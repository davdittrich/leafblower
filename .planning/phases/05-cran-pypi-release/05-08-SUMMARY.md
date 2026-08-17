---
phase: 05-cran-pypi-release
plan: 08
subsystem: testing
tags: [design-effect, statistical-correctness, calibration, greg, henry-valliant, pratools, tdd]

# Dependency graph
requires:
  - phase: 05-cran-pypi-release
    provides: "05-05's r-universe publish, which surfaced the R CMD check ERROR (build 32069801748) this plan root-causes"
provides:
  - "src/design_effect.cpp's 4-argument path now returns a mathematically valid Henry & Valliant (2015) Eq-3.5 statistic (constant column in the calibration design matrix)"
  - "A permanent deff_H <= deff_K regression invariant, independent of any third-party package"
  - "An actually-executing (no skip) PracTools parity test with a var(w)=0-safe oracle"
affects: [05-09, 07-cran-submission]

# Actuals (#2632)
actuals:
  tokens: 3963
  tasks: 3
  commits: 3

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Nested-model invariant testing: deff_H <= deff_K is derived from the constant-only model being nested in the calibration model (weighted SSE <= weighted SST) and pinned as a permanent regression guard, not just satisfied by current source"
    - "Third-party oracle built from primitives defined at the fixture's degenerate point (var(w)=0), with the direct upstream call kept as a conditional (is.finite-gated) tripwire rather than dropped"

key-files:
  created: []
  modified:
    - src/design_effect.cpp
    - tests/testthat/test-design-pratools-parity.R
    - tests/testthat/test-design.R
    - python/leafblower/test_design_effect.py
    - R/design_effect.R
    - man/design_effect.Rd
    - NEWS.md

key-decisions:
  - "Task 1 measured the exact RED numbers locally with PracTools 1.7.5 actually installed (not simulated): d_lbw=8.293020900270720 vs oracle=0.197671531130964 (~42x divergence)"
  - "PracTools::deffH is NaN under uniform weights (upstream 0/0 removable singularity in its Eq-3.4 correction term, R/deffH.R:22-25, version-independent) -- kept as a conditional is.finite() tripwire rather than the primary assertion, which now routes through PracTools::deffK + a glm()-fit WLS Eq-3.5 ratio (a QR path independent of the C++ core's Cholesky path)"
  - "Task 3's measured relative difference between the corrected core and the re-derived hand-rolled oracle is 1.974e-12 (was asserted at 1e-6 pre-fix, never actually measured at that tightness); tolerance tightened to 1e-11 -- the next power of ten above the measurement, per the plan's explicit anti-inflation instruction"
  - "roxygen2::roxygenise() on this machine is v8.0.0 vs. the repo's pinned RoxygenNote 7.3.3, so it rewrote DESCRIPTION's Config/roxygen2/version field on regeneration; reverted that unrelated toolchain-version drift with git checkout, keeping only the doc-text delta in man/design_effect.Rd (git diff --stat confirms roxygen-shaped changes only)"

patterns-established:
  - "GREG residual design matrices must include the constant column: dropping one reference level per calibration margin without an intercept removes a dimension from the model space, making the reference cell's residual equal to its raw outcome and inflating sigma^2_u -- documented in-comment at the fix site, citing Henry & Valliant (2015) Eq 3.5 and leafblower-xfz4"

requirements-completed: [US-010, US-008, KPI-05, KPI-06]

coverage:
  - id: D1
    description: "PracTools parity test executes (0 skips) and fails against the current buggy build with a quantified, non-NaN divergence (RED)"
    requirement: US-010
    verification:
      - kind: unit
        ref: "tests/testthat/test-design-pratools-parity.R::design_effect 4-arg matches PracTools::deffK/glm Eq-3.5 oracle"
        status: pass
    human_judgment: false
  - id: D2
    description: "design_effect() 4-argument path builds a calibration design matrix that spans the constant vector, so its GREG residuals are true GREG residuals and the corrected value (0.1976715311) matches an independent PracTools/glm oracle"
    requirement: US-008
    verification:
      - kind: unit
        ref: "src/design_effect.cpp (p_acc seeded at 1, column 0 filled with 1.0)"
        status: pass
      - kind: unit
        ref: "tests/testthat/test-design-pratools-parity.R (both blocks, 0 skipped, 0 failed)"
        status: pass
    human_judgment: false
  - id: D3
    description: "Permanent deff_H <= deff_K invariant guard pinned in the suite (would have caught leafblower-xfz4 with no third-party package)"
    requirement: US-008
    verification:
      - kind: unit
        ref: "tests/testthat/test-design.R::H&V Eq 3.5 invariant: deff_H <= deff_K (leafblower-xfz4)"
        status: pass
    human_judgment: false
  - id: D4
    description: "R and Python still return byte-identical deff_H for the same input after the C++ core change; empty-target degenerate path still returns exactly deff_K"
    requirement: KPI-06
    verification:
      - kind: unit
        ref: "python/leafblower/test_design_effect_parity.py::test_r_python_parity_4arg_same_c_entry"
        status: pass
      - kind: unit
        ref: "tests/testthat/test-design.R::design_effect 4-arg empty target -> deff_K (degenerate)"
        status: pass
    human_judgment: false
  - id: D5
    description: "Full Definition-of-Done gate green with the fix applied: R build + full testthat suite, and Python parity suite under single-thread BLAS"
    requirement: KPI-05
    verification:
      - kind: integration
        ref: "R CMD INSTALL --preclean . && Rscript -e 'devtools::test()' -> 0 FAIL/1840 PASS/141 WARN/12 SKIP"
        status: pass
      - kind: integration
        ref: "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest -> 161 passed/0 failed"
        status: pass
    human_judgment: false

duration: ~25min
completed: 2026-08-18
status: complete
---

# Phase 05 Plan 08: Design Effect Constant-Column Fix Summary

**Fixed a real statistical-correctness bug in the published `design_effect()`: its calibration design matrix omitted the constant column, so `deff_H` could (and did, by ~42x) exceed `deff_K` — mathematically impossible once the intercept is nested in the model; corrected seed-2024 fixture value moves from 8.2930209003 to 0.1976715311, matching an independent PracTools/glm oracle, with a permanent `deff_H <= deff_K` invariant now pinned in the suite.**

## Performance

- **Duration:** ~25 min
- **Tasks:** 3/3 completed
- **Files modified:** 7

## Accomplishments

- Reproduced the divergence locally with PracTools 1.7.5 actually installed (was previously a silent `skip_if_not_installed` no-op): current build returns 8.293020900270720, an independent `PracTools::deffK * glm-WLS` Eq-3.5 oracle returns 0.197671531130964 — a ~42x divergence, and `PracTools::deffH` itself is NaN under uniform weights (an upstream, version-independent 0/0 removable singularity, kept as a conditional tripwire rather than the primary assertion).
- Root-caused and fixed `src/design_effect.cpp`: the calibration design matrix `X` was built from reference-coded margin dummies with **no constant column**, so `u_i = y_i - x_i^T β̂` was not a true GREG residual (the reference cell's residual equaled its raw outcome). Seeded `p_acc` at 1 and filled column 0 with the constant vector; margin dummies shift right by one column. Corrected value: 0.1976715311, matching the oracle to within the measured tolerance.
- Added a permanent `deff_H <= deff_K` invariant guard (`tests/testthat/test-design.R`), derived from model nesting (constant-only model nested in the calibration model ⇒ weighted SSE ≤ weighted SST), on both the parity fixture and the pre-existing RVAL.4 n=20/K=2 fixture — this is the check that would have caught leafblower-xfz4 with no third-party package at all.
- Re-derived every downstream fixture, oracle, and doc for the corrected column space: the hand-rolled Eq-3.5 oracle in `test-design-pratools-parity.R` (measured relative difference 1.974e-12, tolerance tightened 1e-6 → 1e-11), the 3-level perfect-fit fixture in `test-design.R` and its Python mirror (beta_hat now `[10, 10, 20]`, all residuals ≈0, `deff_H = 0` was `1/3`), roxygen + regenerated `man/design_effect.Rd`, and a `NEWS.md` entry recording the superseded value.
- Full Definition-of-Done gate confirmed green after the fix: R `devtools::test()` 0 FAIL/1840 PASS/141 WARN/12 SKIP (PracTools parity file: 0 skipped); Python `pytest` under single-thread BLAS 161 passed/0 failed, including `test_design_effect_parity.py` (R and Python still agree byte-for-byte on the shared C++ core).

## Task Commits

Each task was committed atomically:

1. **Task 1: RED — reproduce the divergence end-to-end with PracTools installed** - `60f2a57` (test)
2. **Task 2: GREEN — add the constant column so residuals are true GREG residuals** - `30842f4` (fix)
3. **Task 3: Propagate — re-derive fixtures, oracle, roxygen, NEWS** - `aac45d6` (docs)

_No plan-metadata commit: STATE.md/ROADMAP.md updates are owned by the orchestrator per this plan's execution instructions._

## Files Created/Modified

- `src/design_effect.cpp` - Calibration design matrix now spans the constant vector (`p_acc` seeded at 1, column 0 filled with 1.0, margin dummies shift right by one column)
- `tests/testthat/test-design-pratools-parity.R` - Both blocks rewritten: oracle 1 now `PracTools::deffK * glm-WLS` (var(w)=0-safe) with a conditional `deffH` tripwire; oracle 2's hand-rolled `X2` gains the constant column, tolerance tightened to the measured value
- `tests/testthat/test-design.R` - Perfect-fit fixture re-derived (`deff_H = 0`, absolute bound); new permanent `deff_H <= deff_K` invariant guard
- `python/leafblower/test_design_effect.py` - Perfect-fit fixture mirror re-derived to match R
- `R/design_effect.R` - Roxygen states X is a constant column + reference-coded dummies and why that bounds `deff_H <= deff_K`
- `man/design_effect.Rd` - Regenerated via roxygen2 (doc-text delta only; reverted an unrelated `RoxygenNote`/`Config/roxygen2/version` toolchain drift from the local roxygen2 8.0.0 vs. pinned 7.3.3)
- `NEWS.md` - Entry under 0.1.1 recording the corrected numeric behavior and the superseded 8.2930209003 value

## Decisions Made

- Kept `PracTools::deffH` as a conditional (`is.finite`-gated) tripwire rather than dropping it, so the test self-heals if upstream ever removes its NaN-producing singularity, without asserting on a NaN.
- Set Task 3's hand-rolled-oracle tolerance to 1e-11 (next power of ten above the measured 1.974e-12 relative difference), not the plan's ceiling of 1e-6, per the plan's explicit instruction not to raise tolerance above what the measurement justifies.
- Reverted `roxygen2::roxygenise()`'s incidental `DESCRIPTION` `Config/roxygen2/version` bump (local roxygen2 8.0.0 vs. repo-pinned `RoxygenNote: 7.3.3`) — out of this plan's scope, not touched.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking, plan-sequencing gap] Task 2's stated acceptance criteria conflicted with its own action text**

- **Found during:** Task 2
- **Issue:** Task 2's acceptance criteria required `tests/testthat/test-design-pratools-parity.R` to report 0 failures after the C++ fix. But that file's *second* (hand-rolled) block still built its oracle matrix `X2` without a constant column — unavoidably diverging from the now-corrected core — and the plan's own "Artifacts this plan produces" table and Task 3 action item 1 explicitly assign fixing that block to Task 3, while Task 2's action text says "do not touch [test-design.R's perfect-fit fixture]... in this task" (silent on this file's second block, but the same logic applies: fixing it early would violate "do not touch the C++ again to make them pass" scope and duplicate Task 3's stated responsibility).
- **Fix:** Verified Task 2's actual `<done>` criterion directly (corrected numeric value matches the Task-1 oracle to 1e-8; the first parity block passes; the new invariant guard passes) instead of the literal 0-failures-on-that-file bound, and proceeded to Task 3, which fixes the second block as designed.
- **Files modified:** None beyond Task 2's planned scope (`src/design_effect.cpp`, `tests/testthat/test-design.R`).
- **Verification:** `tests/testthat/test-design-pratools-parity.R` reports 0 failures/0 skipped once Task 3 completes (confirmed).
- **Committed in:** `30842f4` (documented inline in the commit message)

**2. [Rule 3 - Blocking] roxygen2 version drift**

- **Found during:** Task 3
- **Issue:** `roxygen2::roxygenise()` on this machine (v8.0.0) rewrote `DESCRIPTION`'s `RoxygenNote: 7.3.3` field to `Config/roxygen2/version: 8.0.0`, an unrelated toolchain-version change out of this plan's scope.
- **Fix:** `git checkout -- DESCRIPTION` to revert the version-field drift, keeping only the doc-text delta in `man/design_effect.Rd` (confirmed via `git diff --stat` to be roxygen-shaped, not hand-edited).
- **Files modified:** `DESCRIPTION` (reverted, not committed).
- **Verification:** `git diff man/design_effect.Rd` shows only the intended `\eqn{}` doc-text addition.
- **Committed in:** N/A (revert, no commit)

---

**Total deviations:** 2 auto-fixed (2 Rule 3). **Impact on plan:** No scope creep — both were plan-sequencing/toolchain frictions resolved by following the plan's own stated Task-3 responsibilities and by reverting an out-of-scope incidental file change; the design_effect() fix itself, the invariant guard, and all fixture/doc/NEWS propagation were executed exactly as specified.

## Issues Encountered

None beyond the two deviations above.

## User Setup Required

None - PracTools installed from CRAN into the local R library during Task 1 (already a declared `Suggests` entry, covered by 05-RESEARCH.md's Package Legitimacy Audit per this plan's threat model T-05-08-03; no new package added to any manifest).

## Next Phase Readiness

The r-universe `R CMD check` ERROR from build 32069801748 (commit 8920cf22) cannot recur: the assertion that produced it now passes against a mathematically correct statistic, verified locally with PracTools installed (not skipped). SC1/KPI-05's remaining blocker moves from "a real numeric defect" to "CI must actually install the Suggests" — 05-09's scope, per this plan's stated success criteria. Any consumer who recorded a 4-argument `deff_H` value from `design_effect()` under 0.1.0 or an earlier 0.1.1 build should treat it as superseded (NEWS.md entry filed).
