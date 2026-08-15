---
phase: 04-truthful-surface
plan: 01
subsystem: api
tags: [r, harvest, input-validation, breaking-change]

requires: []
provides:
  - "harvest(..., weights=w) now raises a hard stop() naming design_weights= instead of silently falling into ... and returning a plausible-but-wrong unweighted result"
  - "RVAL.4 regression test locking the guard in place"
affects: [04-02]

actuals:
  tokens: 1113
  tasks: 2
  commits: 1

tech-stack:
  added: []
  patterns:
    - "RVAL.2 dots-check pattern: a named-typo guard (stop()) inserted ahead of the generic unknown-argument warning, so a plausible misnamed argument fails hard instead of being silently absorbed"

key-files:
  created: []
  modified:
    - R/harvest.R
    - tests/testthat/test-harvest-rval.R
    - tests/testthat/test-logit.R
    - NEWS.md

key-decisions:
  - "D-04 resolved as option-a (hard stop()) at the Task 1 checkpoint: harvest(..., weights=w) now errors naming design_weights=, rather than only sharpening the existing non-fatal warning — a warning remains silently suppressible (the one real in-repo caller already wrapped its call in suppressWarnings()), which would leave SC3's defect reachable."
  - "eb79.18's suppressWarnings() was removed (not just the weights= rename) after an empirical run under design_weights= confirmed zero warnings fire — per the plan's explicit instruction not to strip it by inspection alone, only after confirming."

requirements-completed: []

coverage:
  - id: D1
    description: "harvest(df, tgt, weights=runif(n)) raises an error (not a warning) whose message names design_weights="
    verification:
      - kind: unit
        ref: "tests/testthat/test-harvest-rval.R#RVAL.4: bare weights= arg errors naming design_weights"
        status: pass
    human_judgment: false
  - id: D2
    description: "eb79.18 (test-logit.R) still converges under the design_weights= rename, and existing RVAL.1-3/META.2/dtkn.6 tests are unaffected by the new guard"
    verification:
      - kind: unit
        ref: "tests/testthat/test-logit.R#eb79.18: consistent collinear with HETEROGENEOUS design weights still reaches feasibility"
        status: pass
      - kind: unit
        ref: "tests/testthat/test-harvest-rval.R (RVAL.1-3, META.2, dtkn.6 blocks)"
        status: pass
    human_judgment: false

duration: 25min
completed: 2026-08-15
status: complete
---

# Phase 04 Plan 01: harvest() weights= hard-stop() Summary

**`harvest(..., weights = w)` now raises a hard `stop()` naming `design_weights=` instead of silently falling into `...` and returning a plausible-but-wrong unweighted result (SC3, D-04 option-a).**

## Performance

- **Duration:** 25 min
- **Started:** 2026-08-15T15:26:00Z (approx, per checkpoint resume)
- **Completed:** 2026-08-15T15:51:34Z
- **Tasks:** 2 (Task 1 checkpoint:decision resolved externally as option-a; Task 2 implementation)
- **Files modified:** 4

## Accomplishments
- `R/harvest.R`'s RVAL.2 dots-check gains a `stop()` guard that fires on a bare `weights=` name in `...`, ahead of the existing generic unknown-argument warning, naming `design_weights=` in the error message.
- New `RVAL.4` regression test in `tests/testthat/test-harvest-rval.R` locks the hard-error behavior in place.
- The one real in-repo caller of the old silent fallback (`eb79.18` in `tests/testthat/test-logit.R`) is renamed from `weights = base_w` to `design_weights = base_w` in the same commit — DoD gate never went red.
- `NEWS.md`'s `## Breaking changes` block gains a 4th bullet documenting the change.

## Task Commits

1. **Task 1: Confirm the weights= breaking-change decision (D-04)** — `checkpoint:decision`, resolved externally as option-a (hard `stop()`); no code/files of its own, no commit.
2. **Task 2: Implement the weights= guard, RVAL.4 test, and eb79.18 rename** - `5839b16` (feat)

**Plan metadata:** (this commit)

## Files Created/Modified
- `R/harvest.R` - RVAL.2 block gains an `if ("weights" %in% names(dots)) stop(...)` guard naming `design_weights=`
- `tests/testthat/test-harvest-rval.R` - new `RVAL.4` test appended after the existing RVAL/META/dtkn blocks
- `tests/testthat/test-logit.R` - `eb79.18`'s `weights = base_w` renamed to `design_weights = base_w`; `suppressWarnings()` dropped after confirming empirically zero warnings fire under the rename
- `NEWS.md` - one new bullet in the `## Breaking changes` block (the block matching the plan's "3 prior bullets" description, under the second `# leafblower (development version)` header)

## Decisions Made
- D-04 resolved option-a (hard `stop()`) per the Task 1 checkpoint, matching SC3's literal wording ("raises an informative error").
- `suppressWarnings()` removed from `eb79.18` (beyond the minimal `weights=` → `design_weights=` rename) after an empirical `withCallingHandlers` run confirmed zero warnings fire under the renamed call — the plan's action text permitted removal only on empirical confirmation, not by inspection.

## Deviations from Plan

None - plan executed exactly as written (Task 2's `<action>` text was authored specifically for the option-a decision already resolved at the checkpoint).

## Issues Encountered
None.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- SC3 closed. `04-02-PLAN.md`'s lower-risk doc/comment fixes are unblocked (this plan was the phase's tracer task, proven first per TRACER_MODE).
- Full `devtools::test()` run: `[ FAIL 0 | WARN 141 | SKIP 13 | PASS 1837]` — all 141 warnings and 13 skips pre-exist this plan's changes (unrelated fixtures/benchmarks), verified via targeted `filter='harvest-rval'` (11 pass, 1 pre-existing unrelated warning) and `filter='logit'` (54 pass, 0 warnings) runs.

---
*Phase: 04-truthful-surface*
*Completed: 2026-08-15*
