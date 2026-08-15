---
phase: 04-truthful-surface
reviewed: 2026-08-15T16:01:20Z
depth: standard
files_reviewed: 7
files_reviewed_list:
  - R/harvest.R
  - tests/testthat/test-harvest-rval.R
  - tests/testthat/test-logit.R
  - NEWS.md
  - docs/raking.md
  - src/leafblower.h
  - CLAUDE.md
findings:
  critical: 0
  warning: 1
  info: 1
  total: 2
status: issues_found
---

# Phase 4: Code Review Report

**Reviewed:** 2026-08-15T16:01:20Z
**Depth:** standard
**Files Reviewed:** 7
**Status:** issues_found

## Summary

Two plans landed: (1) `R/harvest.R` gained a hard `stop()` guard rejecting a bare
`weights=` argument, naming `design_weights=`; (2) doc-truthfulness fixes to
`docs/raking.md` (deleted §8.2/§12 ORIS/L-BFGS-B misattribution) and
`src/leafblower.h`/`CLAUDE.md` (annotated `rk_algorithm_t` slot 7 as removed
`RK_ALG_GRAKE`).

**Guard ordering (the review's primary focus): correct.** The `weights=` `stop()`
guard (`R/harvest.R:309-312`) executes before the generic RVAL.2 unknown-arg
`warning()` (`R/harvest.R:313-315`), so a bare `weights=` typo cannot be silently
swallowed into the non-fatal warning path. Confirmed by direct trace and by the
new `RVAL.4` test (`tests/testthat/test-harvest-rval.R:150-162`), which asserts
`expect_error(..., regexp = "design_weights")`. The error message itself is
specific and actionable ("did you mean 'design_weights'? harvest() takes
per-observation design weights via design_weights=, not weights="). The
`test-logit.R` rename (`weights=base_w` → `design_weights=base_w`, line 195) is
the one in-repo caller that would otherwise now hard-error; the rename is
correct and matches the formal parameter name in `harvest()`'s signature
(`R/harvest.R:298`).

**Enum comment change: verified comment-only, no ABI/value-shift risk.**
`rk_algorithm_t` in `src/leafblower.h` uses explicit `= N` assignments for every
member (`RK_ALG_RAKING = 3`, `RK_ALG_ORIS_SOFT = 8`, etc.), so the new
`/* 7 = removed (was RK_ALG_GRAKE) */` comment at line 49 does not shift any
enclosing enumerator's value — `RK_ALG_ORIS_SOFT` stays `8` before and after.
The `EXPECTED_RK_PARAMS_BYTES`/`EXPECTED_RK_RESULT_BYTES` static-assert tripwires
are also untouched by this commit (confirmed via `git show 2a06f1f` — only a
comment line and one `CLAUDE.md` prose line were touched in `src/leafblower.h`).

**`docs/raking.md` misattribution deletion: correctly scoped, exit criteria met.**
`git show 3143594` confirms both the §8.2 "ORIS Solver (renamed from iEPPA)" and
§12 "Synthesis" passages (which attributed an unimplemented entropic-proximal-point
loop to ORIS, and recommended the removed L-BFGS-B algorithm as "the definitive
state-of-the-art") are gone, with no replacement text, matching the commit
message's stated intent. The remaining 9 `L-BFGS-B` mentions in §6.2/§11 are a
documented, deliberate out-of-scope decision (per `.planning/phases/04-truthful-surface/04-CONTEXT.md`
D-01/D-02 and `04-02-SUMMARY.md`) — general-numerical-technique literature
discussion, not a leafblower capability claim — so this review does not
re-flag them as a defect.

One finding below is new and not accounted for in the phase's planning
artifacts.

## Warnings

### WR-01: Undocumented, unrelated deletion in the CLAUDE.md commit

**File:** `CLAUDE.md:27` (pre-commit `2a06f1f`; line now absent)
**Issue:** Commit `2a06f1f` ("docs(leafblower.h,CLAUDE.md): annotate rk_algorithm_t
slot 7 as removed GRAKE") also silently deletes an unrelated bullet:
```
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files
```
from the "Rules" section under "Beads Issue Tracker". This line is not
mentioned anywhere in the commit message (which describes only the slot-7
enum annotation and CLAUDE.md's slot-2 sentence extension), is not referenced
in any of the phase-04 planning artifacts (`04-CONTEXT.md`, `04-02-PLAN.md`,
`04-02-SUMMARY.md`, `04-RESEARCH.md`, `04-DISCUSSION-LOG.md` — grepped, zero
hits for "bd remember" or "MEMORY.md"), and has no replacement text anywhere
else in the file (`grep -n "bd remember\|MEMORY.md" CLAUDE.md` → 0 matches
post-commit). This is a real, standing project convention ("use `bd remember`
for persistent knowledge, not MEMORY.md files") that has now been silently
dropped from the agent-facing instructions with no stated rationale, violating
the repo's own "Surgical Changes: touch only what strictly required" rule
(`CLAUDE.md` §2, still present) and creating risk that future agent sessions
default back to writing MEMORY.md files or otherwise lose this governance rule
with no audit trail explaining why.
**Fix:** Restore the line (or, if its removal was intentional, add a commit
note / beads ticket explaining why the `bd remember`/MEMORY.md convention was
dropped, and confirm no other doc still asserts it should be followed):
```diff
 - Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
 - Run `bd prime` for detailed command reference and session close protocol
+- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files
```

## Info

### IN-01: `weights=` guard is exact-match only; does not cover common near-typos

**File:** `R/harvest.R:309-312`
**Issue:** The guard checks `"weights" %in% names(dots)` exactly. A caller who
writes `harvest(df, tgt, weight = w)` (singular, matching neither the
`weight_column` formal nor the guarded name) or `wt = w` still falls through
to the generic non-fatal `warning()` at line 313-315, reproducing the same
silently-wrong-unweighted-result failure mode SC3/D-04 was written to close —
just under a different misspelling. This is explicitly the one typo form
identified by the ticket (D-04: "the one in-repo caller of the old fallback"
was `weights=`), so it is not a gap in what was scoped, only a residual
partial coverage of the underlying failure class.
**Fix:** None required for this phase; if broader typo coverage is wanted,
file a follow-up ticket to match on a small alias set (e.g.
`c("weight", "wt", "weights")`) rather than silently expanding this phase's
diff.

---

_Reviewed: 2026-08-15T16:01:20Z_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
