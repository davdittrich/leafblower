---
phase: 05-cran-pypi-release
plan: 09
subsystem: testing
tags: [ci, r-cmd-check, suggests, dependency-resolution, cran-comments, github-actions]

# Dependency graph
requires:
  - phase: 05-cran-pypi-release
    provides: "05-08's design_effect() constant-column fix, which this plan's CI change now actually exercises via PracTools instead of skipping it"
provides:
  - "DESCRIPTION Suggests field with no unresolvable entry (autumn removed -- GitHub-only, unsatisfiable from any mainstream repo)"
  - "r-check.yml with no dependency-resolution masking: an unavailable Suggests package now fails the check instead of silently widening the skip set"
  - "A real, durable CRAN-representative R CMD check result (0 errors, 0 warnings) with PracTools/survey/arrow/DiceKriging actually installed and their guarded tests actually executing, proven on a real GitHub Actions run (32080650203), not just locally"
  - "cran-comments.md restated against that Suggests-complete result, with an explicit disclosure of the earlier masked-CI configuration and what corrected it"
affects: [07-cran-submission]

# Actuals (#2632)
actuals:
  tokens: 2500
  tasks: 3
  commits: 3

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "CI dependency-resolution masking is structurally worse than a per-package skip_if_not_installed guard: a hand-maintained 'hard' dependency-class override or a FORCE_SUGGESTS=false env var downgrades EVERY unavailable Suggests package to a skip, not just the one it was added to route around -- the fix is deleting the override entirely and letting R CMD check's own default missing-Suggests error do the guarding, so nothing bespoke has to be written or maintained"

key-files:
  created: []
  modified:
    - DESCRIPTION
    - .github/workflows/r-check.yml
    - cran-comments.md

key-decisions:
  - "Both new local NOTEs found during Task 1's Suggests-complete check (CRAN incoming feasibility 'New submission', HTML manual 'V8 unavailable') are NOT caused by the now-installed Suggests -- verified autumn:: is called nowhere outside .Rbuildignore'd benchmarks/, and both NOTEs are inherent to --as-cran (new-submission boilerplate) or this local machine's R-CMD-check-tooling state (V8 R package absent, unrelated to any package Suggests), not new defects introduced by this plan. Documented, not root-cause-fixed (V8 install falls under Rule 3's package-manager-install exclusion and isn't part of DESCRIPTION anyway)."
  - "CI's R CMD check result stayed at 2 NOTEs (compilation flags, HTML manual) before and after the fix -- no new NOTE appeared on GitHub Actions, confirming the extra local NOTEs are local-environment artifacts, not CI-visible defects from the Suggests change"
  - "Also corrected cran-comments.md's stale tarball filename (leafblower_0.1.0.tar.gz -> leafblower_0.1.1.tar.gz, matching the actual DESCRIPTION Version) while restating the check-evidence section, since Task 3 already required rewriting that exact line for accuracy"

patterns-established: []

requirements-completed: [US-010, KPI-05]

coverage:
  - id: D1
    description: "autumn (>= 0.2.0) removed from DESCRIPTION Suggests; every remaining Suggests entry is CRAN-resolvable and installed locally, and a Suggests-complete R CMD check --as-cran (default missing-Suggests guard, no skip override) reports 0 errors, 0 warnings on the built tarball"
    requirement: US-010
    verification:
      - kind: integration
        ref: "R CMD build . && R CMD check --as-cran leafblower_0.1.1.tar.gz -> Status: 3 NOTEs, 0 errors, 0 warnings, 0 'is not installed' lines in 00check.log"
        status: pass
    human_judgment: false
  - id: D2
    description: "r-check.yml no longer masks Suggests installation (_R_CHECK_FORCE_SUGGESTS_ override and dependencies: 'hard' key both deleted); a real GitHub Actions run under the corrected workflow installs PracTools/survey/arrow/DiceKriging and their guarded tests execute rather than skip"
    requirement: KPI-05
    verification:
      - kind: e2e
        ref: "GitHub Actions run 32080650203 (davdittrich/leafblower, workflow r-check.yml) -- conclusion success; log shows 'Installed DiceKriging 1.6.1', 'Installed PracTools 1.7.5', 'Installed arrow 25.0.0', 'Installed survey 4.5'; testthat summary [ FAIL 0 | WARN 130 | SKIP 31 | PASS 1794 ], 0 skips attributable to a missing package; R CMD check Status: 2 NOTEs (0 errors, 0 warnings)"
        status: pass
    human_judgment: false
  - id: D3
    description: "cran-comments.md restated against the Suggests-complete result: cites GitHub Actions run 32080650203 verbatim, states the Suggests set was installed, explains every NOTE reported by the local Suggests-complete check, discloses the earlier masked-CI configuration and its correction, and makes no claim about a Windows check result"
    requirement: US-010
    verification:
      - kind: unit
        ref: "grep -c 32080650203 cran-comments.md -> 4; grep -ci suggests cran-comments.md -> 19; grep -in windows cran-comments.md -> no match; all 3 local NOTE headings (CRAN incoming feasibility, compilation flags used, HTML version of manual) have matching explanation paragraphs"
        status: pass
    human_judgment: false

duration: ~22min
completed: 2026-08-18
status: complete
---

# Phase 05 Plan 09: Un-mask the Suggests Skip Summary

**Removed the unresolvable `autumn` Suggests entry and the two CI settings that were silently downgrading every unavailable Suggests package (not just `autumn`) to a skip -- a real GitHub Actions run (32080650203) now shows PracTools, survey, arrow and DiceKriging actually installing and their ~25 guarded tests actually executing, with `R CMD check --as-cran` still 0 errors/0 warnings.**

## Performance

- **Duration:** ~22 min
- **Tasks:** 3/3 completed
- **Files modified:** 3

## Accomplishments

- Deleted `autumn (>= 0.2.0)` from `DESCRIPTION`'s Suggests field -- confirmed by repo-wide search that no shipped `R/`, `src/`, or `man/` file calls `autumn::` (only roxygen prose in `R/harvest.R`/`man/harvest.Rd`, comments in `src/leafblower.h`, and `.Rbuildignore`d `benchmarks/` reference it). All 12 remaining Suggests entries verified CRAN-resolvable (`available.packages()`) and locally installed. A Suggests-complete `R CMD check --as-cran` on the built tarball, with the default missing-Suggests guard restored (no skip override), reports 0 errors, 0 warnings, 3 NOTEs, and 0 "is not installed" lines in `00check.log`.
- Deleted the job-level `_R_CHECK_FORCE_SUGGESTS_: false` env override and the `setup-r-dependencies` `dependencies: '"hard"'` key from `.github/workflows/r-check.yml` -- both existed only to route around `autumn`'s unresolvability, but together downgraded every unavailable Suggests package to a silent skip. Pushed and watched GitHub Actions run **32080650203** to completion (`success`): log confirms `✔ Installed DiceKriging 1.6.1`, `✔ Installed PracTools 1.7.5`, `✔ Installed arrow 25.0.0`, `✔ Installed survey 4.5`; testthat summary `[ FAIL 0 | WARN 130 | SKIP 31 | PASS 1794 ]` with all 31 skips traced to `On CRAN` gates or local-only benchmark-data guards -- none reference a missing package (0 "is not installed" lines); `R CMD check` Status stayed `2 NOTEs` (0 errors, 0 warnings), unchanged from before the fix.
- Restated `cran-comments.md`: cites run 32080650203, states the full Suggests set installed in both local and CI environments, adds a "Suggests availability correction" section disclosing the earlier masked-CI configuration and the commits (`8e28df1`, `790341e`) that fixed it, and updates the NOTE inventory to match what the Suggests-complete runs actually reported (the `_R_CHECK_FORCE_SUGGESTS_=false`/`skip_if_not_installed`-guarded framing from the old text no longer applies -- those tests execute for real now). Makes no Windows claim (deferred to 05-10).

## Task Commits

Each task was committed atomically:

1. **Task 1: Make every declared Suggests resolvable, then prove a Suggests-complete check locally** - `8e28df1` (fix)
2. **Task 2: Delete the two lines that masked the skip, and prove on real CI that the Suggests tests execute** - `790341e` (fix)
3. **Task 3: Restate cran-comments.md against the Suggests-complete result** - `cbecdab` (docs)

_No plan-metadata commit: STATE.md/ROADMAP.md updates are owned by the orchestrator per this plan's execution instructions._

## Files Created/Modified

- `DESCRIPTION` - `autumn (>= 0.2.0)` removed from Suggests; other 12 entries untouched
- `.github/workflows/r-check.yml` - Job-level `_R_CHECK_FORCE_SUGGESTS_: false` and `setup-r-dependencies`'s `dependencies: '"hard"'` both deleted; the rest of the `with:` block (`extra-packages`, `needs`, `install-pandoc`) untouched
- `cran-comments.md` - Restated check evidence against the Suggests-complete run (local and CI 32080650203), new "Suggests availability correction" disclosure section, corrected tarball filename to `leafblower_0.1.1.tar.gz`, updated NOTE explanations

## Decisions Made

- Investigated but did not "fix" two new local NOTEs (CRAN incoming feasibility "New submission", HTML manual "V8 unavailable") found during Task 1's Suggests-complete run: confirmed neither is caused by the Suggests change (autumn's removal only simplifies the incoming-feasibility NOTE's text, doesn't create it; the V8 NOTE is R CMD check's own tooling state, unrelated to any package Suggests). Documented both in `cran-comments.md` per Task 3 rather than installing V8 (package-manager installs are Rule 3-excluded, and V8 isn't a declared Suggests entry anyway).
- Corrected `cran-comments.md`'s stale `leafblower_0.1.0.tar.gz` tarball filename to `leafblower_0.1.1.tar.gz` while rewriting the check-evidence section for accuracy (matches DESCRIPTION's actual `Version: 0.1.1`).
- Treated Task 1 (`type="tracer"`) as auto-verified after its automated `<verify>` passed (0 errors/0 warnings/0 "not installed" lines) and proceeded directly to Task 2 rather than pausing for a separate human-verify checkpoint -- this is an infrastructure/CI-config plan with fully automated, non-visual verification criteria, consistent with the session's active auto-mode framing.

## Deviations from Plan

None - plan executed exactly as written. The two extra local NOTEs discovered during Task 1 were investigated as the plan's own action text required ("if a NEW note appears that is caused by the now-installed Suggests, record it verbatim...") and confirmed NOT to be caused by the Suggests change, so no fix was needed -- this was anticipated plan behavior, not a deviation.

## Issues Encountered

None.

## User Setup Required

None - `gh auth status` was already authenticated and this clone could already push to `origin` (the plan's Task 2 precondition), set up in an earlier session (05-05).

## Next Phase Readiness

SC1/KPI-05's "0 errors, 0 warnings" claim is now durable: it comes from a check run under CRAN-representative conditions (every declared Suggests installed, missing-Suggests treated as an error) rather than a check that was quietly skipping ~25+ tests. The configuration that produced the misleading green result is gone, and R CMD check's own built-in guard now prevents its return -- no bespoke masking-detection code was needed or added. The remaining Phase 5 gap is SC6's Windows compile failure, out of this plan's scope and assigned to 05-10.

## Self-Check: PASSED

All 3 modified files confirmed present on disk (`DESCRIPTION`, `.github/workflows/r-check.yml`, `cran-comments.md`); all 3 commits (`8e28df1`, `790341e`, `cbecdab`) confirmed present in `git log --oneline`.

---
*Phase: 05-cran-pypi-release*
*Completed: 2026-08-18*
