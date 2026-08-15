---
phase: 01-verification-coverage-closed
plan: 04
subsystem: testing
tags: [testthat, r, edition-3, waldo, sc5, deprecation]

# Dependency graph
requires: ["01-03"]
provides:
  - "DESCRIPTION Config/testthat/edition: 3 -- the package now runs under the
    testthat edition CLAUDE.md documents, closing ROADMAP SC5 and
    leafblower-og7d."
  - "11 test files with the deprecated per-file context() declaration
    removed (edition 3 derives grouping from filename)."
  - "2 test files (test-autocollapse.R, test-harvest-rval.R) with
    expect_warning() scoped via an outer suppressWarnings() so an
    incidental, unrelated diagnostic warning from harvest() no longer
    leaks past the assertion under 3e's non-swallowing semantics."
affects: [ci-gate, future-test-authoring]

# Actuals (#2632) -- chars/4 over the realized diff, this plan's work only.
actuals:
  tokens: 2249
  tasks: 2
  commits: 3

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "expect_warning() under testthat 3e only muffles the ONE warning
      matching its regexp/class; any other warning raised by the same
      expression propagates instead of being silently swallowed as it was
      under 2e. Fix: wrap the whole expect_warning() call in an OUTER
      suppressWarnings() (never the reverse -- an inner suppressWarnings()
      around the harvest() call would consume the target warning too,
      before expect_warning()'s own handler ever sees it). Verified by
      direct probe: expect_warning() still fails correctly when its target
      warning never fires, with or without the outer suppressWarnings()."
    - "Fix-first-then-flip (D-12), generalized: a fix discovered only by
      running the NEW edition (here, the expect_warning scoping bug) still
      belongs in a commit that lands and stays green under the OLD edition
      first, if it is edition-2-safe (verified: identical FAILED/ERROR/
      PASSED/WARNING counts under 2e with or without the fix). That keeps
      the final flip commit isolated to the metadata field alone, which is
      what a bisect -- and this plan's own acceptance criteria -- expect."

key-files:
  created: []
  modified:
    - "DESCRIPTION"
    - "tests/testthat/test-alm-config-grouping.R"
    - "tests/testthat/test-cell-table.R"
    - "tests/testthat/test-compare.R"
    - "tests/testthat/test-oris-b12-fallback-best-reset.R"
    - "tests/testthat/test-oris-b13-best-error-honesty.R"
    - "tests/testthat/test-oris-dispatch.R"
    - "tests/testthat/test-oris-faithful.R"
    - "tests/testthat/test-oris-nonuniform-d.R"
    - "tests/testthat/test-oris.R"
    - "tests/testthat/test-oris-sraa-log-path.R"
    - "tests/testthat/test-oris-sraa.R"
    - "tests/testthat/test-autocollapse.R"
    - "tests/testthat/test-harvest-rval.R"

key-decisions:
  - "Plan's Task 2 diagnosis (\"the +5 warning delta is from context()
    deprecation notices\") was disproved by direct A/B measurement: after
    Task 1 removed all 11 context() calls, edition 3 still reported
    WARNING=146, identical to the pre-cleanup number. The real +5 came from
    testthat 3e's expect_warning() no longer swallowing non-matching
    warnings raised by the same expression -- a genuine D-11 finding,
    investigated individually rather than triaged into a blanket '3e
    artifact' bucket."
  - "Landed the expect_warning scoping fix as its OWN commit, before the
    edition-3 flip, rather than folding it into the flip commit -- keeps
    Task 2's acceptance criterion (\"git diff --name-only HEAD~1 lists
    DESCRIPTION and nothing else\") literally true, and keeps D-12's
    fix-first-then-flip bisectability intact for this newly-discovered fix
    too. Verified the fix is edition-2-safe (no observable behavior change
    under 2e) before committing it ahead of the flip."
  - "Fix direction: outer suppressWarnings() around expect_warning(), not
    the reverse. Verified empirically with a two-warning probe function
    that an inner suppressWarnings() would consume the target warning
    before expect_warning() ever saw it (breaking the assertion), while
    the outer-wrap correctly lets expect_warning() catch+verify its target
    first and only silences what escapes past it."

patterns-established:
  - "Baseline-vs-observed reconciliation before touching anything: Task 1's
    action explicitly required matching the plan's edition-2 baseline row
    before any edit. PASSED reported 1833, not the plan's 1025 -- reconciled
    as 1025 + Plan 03's ~808 new bound-property expectations (01-03-SUMMARY
    confirms), not a discrepancy requiring a halt."

requirements-completed: [SC5, leafblower-og7d]

coverage:
  - id: D1
    description: "Rscript -e 'testthat::test_dir(...)' runs under testthat edition 3, matching CLAUDE.md's documented claim"
    requirement: "SC5"
    verification:
      - kind: build
        ref: "DESCRIPTION Config/testthat/edition: 3 + testthat::edition_get() == 3 at test-run time"
        status: pass
    human_judgment: false
  - id: D2
    description: "The deprecation fixes and the edition flip are separate, independently-green commits so a bisect can distinguish mechanical churn from real fallout"
    requirement: "D-12"
    verification:
      - kind: manual
        ref: "git log --oneline: da7bff9 (context() removal) -> d6c2c5c (expect_warning scoping fix) -> dc57317 (DESCRIPTION-only flip); each commit's own suite run recorded green"
        status: pass
    human_judgment: false
  - id: D3
    description: "Every newly-surfaced warning under the edition flip investigated as a finding, not blanket-triaged"
    requirement: "D-11"
    verification:
      - kind: manual
        ref: "A/B diff of sorted warning message text (edition 2 vs edition 3, both post-Task-1) isolated the exact 4 call sites and the exact mechanism (expect_warning non-swallowing); fixed by scoping, not by loosening or ignoring"
        status: pass
    human_judgment: false

# Metrics
duration: ~70min
completed: 2026-08-15
status: complete
---

# Phase 1 Plan 04: testthat Edition 3 Flip Summary

**`DESCRIPTION` now carries `Config/testthat/edition: 3`, making CLAUDE.md's documented "testthat v3" claim true; getting there required not just the planned deprecation cleanup but discovering and fixing a real D-11 finding — four `expect_warning()` call sites that 3e's stricter semantics exposed as accidentally swallowing an unrelated diagnostic warning under 2e.**

## Performance

- **Duration:** ~70 min
- **Tasks:** 2 of 2 complete (plus one D-11-driven fix commit between them)
- **Files modified:** 14 (11 test files' `context()` removed, 2 test files' `expect_warning()` scoped, `DESCRIPTION`)

## Accomplishments

- Enumerated the deprecated-`context()` files at execution time (`grep -l '^context(' tests/testthat/*.R`) rather than trusting the plan's frontmatter snapshot — the 11-file list matched exactly.
- Established the edition-2 baseline first per Task 1's own instruction: `FAILED=0 ERROR=0 PASSED=1833 WARNING=141`. `PASSED=1833` (not the plan's measured `1025`) reconciled immediately as `1025 + ~808` — Plan 03's bound-property expectations, landed after this plan's baseline was measured, per 01-03-SUMMARY.
- Removed the `context(...)` line (and its stranded blank line) from all 11 files, verified still edition 2, still green, committed as its own change (`da7bff9`).
- Added `Config/testthat/edition: 3` and re-ran the suite: `WARNING=146`, exceeding the plan's `<=141` acceptance ceiling — the same 146 as the plan's own pre-cleanup "edition 3, no other change" baseline row. This proved the plan's causal diagnosis wrong: Task 1's `context()` removal had **zero** effect on the edition-3 warning count.
- Investigated per D-11 rather than accepting or blanket-triaging: sorted and diffed the full warning-message list between an edition-2 run and the edition-3 run (both post-Task-1). The +5 delta was 4 call sites (`test-autocollapse.R:26,48,57`, `test-harvest-rval.R:72`) where `expect_warning(harvest(...), "<target regexp>")` also incidentally triggers `harvest()`'s sparse-category or fixed-point-convergence diagnostic warning. Edition 2's `expect_warning()` swallows every warning raised in its expression; edition 3's only swallows the one matching its regexp/class and lets the rest propagate — a documented, non-buggy 2e-vs-3e behavior change, confirmed with a minimal two-warning probe function before touching any test file.
- Fixed by wrapping each `expect_warning()` call in an *outer* `suppressWarnings()` — verified by probe that this ordering (never the reverse) preserves assertion strength: `expect_warning()`'s own handler still fires first and still fails the test if the target warning never occurs; the outer wrap only catches what escapes past it.
- Verified the fix is edition-2-safe (identical `FAILED/ERROR/PASSED/WARNING` to the untouched baseline under edition 2) and committed it as its own commit *before* the flip (`d6c2c5c`), preserving both D-12's fix-first-then-flip bisectability and Task 2's own acceptance criterion that the flip commit's diff contain `DESCRIPTION` alone.
- Re-added the field, re-verified `EDITION=3 FAILED=0 ERROR=0 PASSED=1833 WARNING=141` — back at the edition-2 floor — and `R CMD INSTALL --preclean .` completing `DONE (leafblower)`, then committed the DESCRIPTION-only flip (`dc57317`).
- Ran the complete `.coverage-thresholds.json` `enforcement.command` end to end after all three commits: R half `FAIL 0 | WARN 141 | SKIP 13 | PASS 1833`; Python half `156 passed` under single-thread BLAS. Both green.

## Task Commits

Each unit of work was committed atomically, in fix-first-then-flip order (D-12):

1. **Task 1: Remove the deprecated per-file context declarations while still on edition 2 (D-12)** — `da7bff9` (test)
2. **D-11-driven fix (discovered during Task 2's investigation, landed ahead of the flip): scope `expect_warning()` to its target warning** — `d6c2c5c` (test)
3. **Task 2: Opt into testthat edition 3 (SC5, D-10, D-11, D-12)** — `dc57317` (chore) — DESCRIPTION only

**Plan metadata:** this SUMMARY.md commit (docs) — final metadata commit for this plan. STATE.md and ROADMAP.md are intentionally NOT touched by this run; the orchestrator updates those after all plans in this wave complete.

## Files Created/Modified

- `DESCRIPTION` — new `Config/testthat/edition: 3` DCF field.
- `tests/testthat/test-alm-config-grouping.R`, `test-cell-table.R`, `test-compare.R`, `test-oris-b12-fallback-best-reset.R`, `test-oris-b13-best-error-honesty.R`, `test-oris-dispatch.R`, `test-oris-faithful.R`, `test-oris-nonuniform-d.R`, `test-oris.R`, `test-oris-sraa-log-path.R`, `test-oris-sraa.R` — deprecated `context(...)` line (and its stranded blank line) deleted; no other change.
- `tests/testthat/test-autocollapse.R` — 3 `expect_warning()` sites wrapped in an outer `suppressWarnings()`.
- `tests/testthat/test-harvest-rval.R` — 1 `expect_warning()` site wrapped in an outer `suppressWarnings()`.

## Decisions Made

- **D-11 finding, root cause, not symptom:** the plan's own measured baseline table correctly recorded the numbers (146 vs 141) but misattributed the mechanism to `context()` deprecation notices. Re-running edition 3 with `context()` already removed and still seeing 146 falsified that attribution immediately — investigated with a message-level diff rather than assuming the plan's stated cause was correct. See `tech-stack.patterns` above for the actual mechanism and the fix.
- **Fix ordering:** `suppressWarnings(expect_warning(...))`, never `expect_warning(suppressWarnings(...))`. R's calling-handler stack invokes the most-recently-established (innermost) handler first; nesting `suppressWarnings()` inside `expect_warning()`'s argument would let it consume the target warning via its own `muffleWarning` restart before `expect_warning()`'s outer-established handler ever observes it, silently breaking the assertion. Verified both orderings with a two-warning probe function before editing any real test.
- **Commit sequencing:** the scoping fix could have been folded into the flip commit (Task 2's action text anticipates "if any assertion newly fails... fix... commit" as part of that task). Chose instead to verify it was edition-2-safe and land it as an independent commit ahead of the flip, because (1) Task 2's own acceptance criteria literally require the flip commit's diff to be `DESCRIPTION` and nothing else, and (2) D-12's whole rationale — isolating mechanical/hygiene churn from true edition semantics — applies exactly as well to a fix that happens to only be *observable* under the new edition but is not *caused* by it.
- **No changes to `harvest()` or any `src/`/`R/` core code.** The sparse-category and fixed-point-convergence diagnostic warnings are legitimate, documented `harvest()` output (the messages themselves point at `result$diagnostics$sparseness` and suggest `accelerate=TRUE`/`method='newton_kl'`); nothing was hidden or wrong in the solver layer. The test-layer-only boundary this phase operates under held throughout.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] `expect_warning()` leaked incidental warnings under testthat 3e, exceeding the plan's WARNING<=141 acceptance ceiling**
- **Found during:** Task 2, after adding the `Config/testthat/edition: 3` field and re-running the suite (`WARNING=146`, ceiling `141`).
- **Issue:** 4 `expect_warning()` call sites across 2 test files pass because `harvest()` also incidentally raises an unrelated diagnostic warning (sparse-category detection or fixed-point convergence) in the same call. Edition 2's `expect_warning()` silently swallows every warning in its wrapped expression, hiding the incidental one; edition 3 only swallows the one matching the regexp, letting the rest propagate and inflate the WARNING count.
- **Fix:** Wrapped each `expect_warning()` call in an outer `suppressWarnings()`, verified by direct probe not to weaken the assertion (still fails if the target warning is absent) and verified edition-2-safe (identical counts to the untouched baseline) before landing.
- **Files modified:** `tests/testthat/test-autocollapse.R`, `tests/testthat/test-harvest-rval.R`
- **Commit:** `d6c2c5c`

## Issues Encountered

None beyond the D-11 finding documented above, which was investigated to a definitive mechanism (not left as an unexplained artifact) and fixed at the test layer without touching `src/`, `R/`, or `python/`.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

**Plan complete.** ROADMAP SC5 is satisfied: the suite runs under testthat edition 3, `DESCRIPTION` and `CLAUDE.md` agree, and `leafblower-og7d` is closable. The complete `.coverage-thresholds.json` `enforcement.command` (R build + testthat + Python parity, single-thread BLAS) runs green end to end: R `FAIL 0 | WARN 141 | SKIP 13 | PASS 1833`; Python `156 passed`. Every commit in this plan's chain (`da7bff9` -> `d6c2c5c` -> `dc57317`) is independently green, honoring the project's local-only "complete = committed + gates green" definition and D-12's bisectability requirement.

---
*Phase: 01-verification-coverage-closed*
*Completed: 2026-08-15*

## Self-Check: PASSED

- FOUND: `DESCRIPTION`
- FOUND: `tests/testthat/test-alm-config-grouping.R` (and the other 10 `context()`-cleanup files)
- FOUND: `tests/testthat/test-autocollapse.R`
- FOUND: `tests/testthat/test-harvest-rval.R`
- FOUND: `.planning/phases/01-verification-coverage-closed/01-04-SUMMARY.md`
- FOUND commit: `da7bff9` (Task 1)
- FOUND commit: `d6c2c5c` (D-11 fix)
- FOUND commit: `dc57317` (Task 2)
