---
phase: 05-cran-pypi-release
plan: 01
subsystem: infra
tags: [r-package, cran, Rbuildignore, roxygen2, ascii-portability, news-md]

# Dependency graph
requires:
  - phase: 04-docs-audit-and-correction
    provides: corrected docs/ (no ORIS/L-BFGS-B misattribution), clean AGENTS.md audit trail
provides:
  - hygiene-clean git tree (no tracked dev artifacts, no leafblower.Rcheck/ tree)
  - .Rbuildignore covering every tracked-and-untracked non-package path R CMD build was
    sweeping into the tarball (58 patterns total, up from 23)
  - ASCII-only R/*.R source and man/harvest.Rd (LaTeX PDF-manual ERROR + non-ASCII WARNING
    both eliminated at the root cause)
  - parseable NEWS.md (three "(development)" headers retitled to "leafblower 0.1.0")
  - cran-comments.md documenting a real, non-placeholder R CMD check --as-cran result
affects: [05-02, 05-03, 05-04, 05-05]

# Actuals (#2632)
actuals:
  tokens: 91906
  tasks: 2
  commits: 2

# Tech tracking
tech-stack:
  added: []
  patterns:
    - ".Rbuildignore must exclude by working-tree presence, not git-tracked status — R CMD
      build scans the filesystem, not the git index, so untracked dev-tool cruft (agent
      caches, JS tooling stubs, benchmark output dirs) ships in the tarball unless
      Rbuildignore'd too"
    - "R_MAKEVARS_USER=/dev/null isolates a check run from a developer's personal
      ~/.R/Makevars (e.g. -march=native), producing a result representative of what a
      stranger's machine (or CRAN's own build farm) would actually see"

key-files:
  created:
    - cran-comments.md
  modified:
    - .Rbuildignore
    - DESCRIPTION
    - NEWS.md
    - R/harvest.R
    - R/anesrake.R
    - R/current_miss.R
    - R/diagnose_weights.R
    - R/weighted_pct.R
    - man/harvest.Rd

key-decisions:
  - "Extended .Rbuildignore 21 patterns beyond the plan's declared 14 after the real R CMD
    build run surfaced untracked dev-tool cruft (graphify-out/ cache, .planning/, .tldr/,
    .gemini*/, hidden agent-tool dirs, stray node_modules/package.json, data-raw/, a
    benchmark weights/ output dir) being swept into the tarball — none of this was git-rm'd
    since none of it (except data-raw/, standard R convention) was git-tracked"
  - "Fixed the LaTeX PDF-manual ERROR and non-ASCII WARNING at the root: replaced every
    non-ASCII math/typography character (Sigma, Delta, em/en-dash, >=, ~=, ->, x, lambda,
    ||...||, ^2) in R/harvest.R's roxygen comments and 4 other R/*.R files with ASCII
    equivalents, then hand-mirrored the same substitution into the already-generated
    man/harvest.Rd rather than re-running roxygen2 (installed 8.0.0 vs. DESCRIPTION's pinned
    RoxygenNote 7.3.3 — regenerating risked an unrelated formatting diff)"
  - "Retitled NEWS.md's three '(development)' headers to 'leafblower 0.1.0' (verified via
    tools:::.build_news_db_from_package_NEWS_md — R's version-number regex doesn't match
    free text) — content of every bullet left untouched; this is the first release, so
    collapsing accumulated dev history under one version header is standard practice"
  - "Added arrow/callr/jsonlite/withr to DESCRIPTION Suggests (used by tests/, undeclared)"
  - "Did NOT install checkbashisms, HTML Tidy, or the R V8 package (each would eliminate one
    remaining WARNING/NOTE) — installing check-environment tooling is a package-manager
    install, explicitly excluded from Rule 3 auto-fix, and this session's sandbox blocks
    pacman regardless; documented as a local-environment gap in cran-comments.md that
    05-03/05-04's CI matrix is expected to close"
  - "Did use an already-present, unrelated pandoc binary (bundled with a local Quarto
    install at /opt/quarto/bin/tools/x86_64/pandoc) via a one-off PATH prepend for the
    check run — this closed half the top-level-files WARNING (the pandoc-unavailable half);
    not a new install, just pointing at a binary already on disk"
  - "Ran the check with _R_CHECK_FORCE_SUGGESTS_=false and R_MAKEVARS_USER=/dev/null
    prepended to the plan's literal verify command — PracTools isn't installed locally,
    this machine's autumn is pinned at 0.1 (Suggests needs >=0.2.0), and ~/.R/Makevars
    injects -march=native for every package on this machine; none of the three are
    leafblower defects, and R's own check-failure message names FORCE_SUGGESTS as the
    sanctioned workaround"

requirements-completed: [US-010, KPI-05]

coverage:
  - id: D1
    description: "Hygiene-clean git tree: the 11 named dev-artifact strays and the tracked leafblower.Rcheck/ tree (120 files) removed from git tracking, with matching .Rbuildignore patterns as defense-in-depth"
    requirement: US-010
    verification:
      - kind: other
        ref: "git ls-files | grep -Ec '^(cell_table_92c4f45\\.(cpp|hpp)|ieppa_92c4f45\\.cpp|patch_raking\\.py|patch_wolfe\\.py|test_output\\.log|leafblower_0\\.1\\.0\\.tar\\.gz|REVIEW_FINDINGS\\.md|code-review-findings\\.md|baseline_bench\\.R|security-review-2026-04-27-ieppa\\.json|leafblower\\.Rcheck/)' -> 0"
        status: pass
    human_judgment: false
  - id: D2
    description: "R CMD check --as-cran on the real, cleaned tarball: 0 ERRORs, down from 1 ERROR + 4 WARNINGs + 6 NOTEs on the first real run to 1 WARNING + 3 NOTEs (checkbashisms/tidy/V8 tooling gaps + the expected new-submission NOTE + a local-Makeconf compilation-flags NOTE) — does not literally reach the plan's stated 'at most the new-submission NOTE' bar"
    requirement: KPI-05
    verification:
      - kind: other
        ref: "R CMD build . && _R_CHECK_FORCE_SUGGESTS_=false R_MAKEVARS_USER=/dev/null R CMD check --as-cran leafblower_0.1.0.tar.gz -> Status: 1 WARNING, 3 NOTEs (0 ERRORs)"
        status: fail
    human_judgment: true
    rationale: "The plan's must_have truth ('0 errors and 0 warnings, with at most the new-submission NOTE') is not literally met. The 0-ERROR / 0-non-portability-WARNING bar IS met after fixing a real LaTeX-manual ERROR, a real non-ASCII WARNING, and a real unstated-test-dependency WARNING at the root cause. The one remaining WARNING (checkbashisms absent) and two of three NOTEs (tidy/V8 absent; this machine's own R Makeconf/-march flags) are local check-environment tooling gaps this session is procedurally barred from closing (Rule 3 excludes package-manager installs; the sandbox also blocks pacman directly) — not package defects. A human should confirm this reasoning is acceptable, or run the check on a machine/CI runner where checkbashisms/tidy/V8 are installed to get the literal 0-WARNING/1-NOTE result 05-03/05-04's CI matrix is expected to produce."
  - id: D3
    description: "cran-comments.md documents the real observed result and quotes the no-O3/LAPACK build rationale verbatim from source"
    requirement: US-010
    verification:
      - kind: other
        ref: "grep -c '^## Submission type$|^## Test environments$|^## R CMD check results$|^## Notes on build configuration$|find_package(LAPACK REQUIRED)' cran-comments.md -> all 1"
        status: pass
    human_judgment: false

duration: 19min
completed: 2026-08-15
status: complete
---

# Phase 05 Plan 01: CRAN Tracer — Hygiene Cleanup and Real R CMD Check Summary

**Removed all git-tracked dev artifacts (11 strays + 120-file leafblower.Rcheck/ tree),
closed a LaTeX-manual-breaking non-ASCII-Unicode ERROR and a non-ASCII WARNING at the root
in R/harvest.R (and four sibling R files), made NEWS.md parseable, and ran the real
`R CMD build . && R CMD check --as-cran` cycle end-to-end for the first time in this
project's history — down from 1 ERROR/4 WARNINGs/6 NOTEs to 0 ERRORs/1 WARNING/3 NOTEs,
with the remaining findings traced to this machine's missing check-tooling (checkbashisms,
tidy, V8) rather than package defects.**

## Performance

- **Duration:** 19 min
- **Started:** 2026-08-15T17:14:03Z
- **Completed:** 2026-08-15T17:33:01Z
- **Tasks:** 2
- **Files modified:** 149 (140 in Task 1's commit, 1 new file in Task 2's commit, 8 more
  than the plan's declared `[".Rbuildignore", "cran-comments.md"]` due to Rule 1/2 fixes
  found by actually running the check — see Deviations)

## Accomplishments

- `git rm`'d the 11 named root-level dev-artifact strays and the entire tracked
  `leafblower.Rcheck/` output tree (120 files); `.Rbuildignore` grew from 23 to 58 lines,
  covering both the removed tracked strays and 21 more untracked working-tree paths (agent
  tool caches, JS tooling stubs, `.planning/`, `.tldr/`, `graphify-out/`, a benchmark
  `weights/` output dir, `data-raw/`) that `R CMD build` — which scans the filesystem, not
  the git index — was sweeping into the tarball regardless of git-tracked status.
- Ran the real `R CMD build . && R CMD check --as-cran` cycle three times as fixes landed,
  going from `1 ERROR, 4 WARNINGs, 6 NOTEs` (first real run against the cleaned tree) to
  `0 ERRORs, 1 WARNING, 3 NOTEs` (final run). Fixed at the root: a LaTeX PDF-manual ERROR
  and a non-ASCII-characters WARNING (both caused by Unicode math symbols — Sigma, Delta,
  em-dash, `>=`, `~=`, `->`, `x`, `lambda`, `||...||`, `^2` — in roxygen comments), an
  unstated-test-dependencies WARNING (`arrow`/`callr`/`jsonlite`/`withr` missing from
  DESCRIPTION Suggests), and a NEWS.md "no news entries found" NOTE.
- Wrote `cran-comments.md` with the real, non-placeholder check result, explaining which
  remaining findings are local check-environment tooling gaps (vs. package defects) and
  quoting the no-`-O3`/hard-LAPACK-dependency build-flag rationale verbatim from
  `src/Makevars.in` and `python/CMakeLists.txt`.

## Task Commits

1. **Task 1: Remove tracked dev artifacts, extend .Rbuildignore, prove R CMD check --as-cran
   clean** - `4798d1b` (chore)
2. **Task 2: Write cran-comments.md against the real Task 1 result** - `5bddbfa` (docs)

_Note: Task 1's single commit bundles the hygiene `git rm` + `.Rbuildignore` extension with
the ASCII-portability and NEWS.md/DESCRIPTION fixes discovered while actually running the
check — these are interdependent (the fixes only exist because Task 1's own action ran the
real build/check cycle) and were verified together before any commit, per the tracer task's
"prove it end-to-end, once, for real" objective._

## Files Created/Modified

- `.Rbuildignore` - 23 -> 58 lines; 14 patterns from the plan (11 stray files + 3 kept-tracked
  tooling config) plus 21 more found by the real build run
- `cran-comments.md` - new; real R CMD check result + build-flag rationale
- `DESCRIPTION` - added arrow/callr/jsonlite/withr to Suggests
- `NEWS.md` - 3 headers retitled `# leafblower 0.1.0` (was `# leafblower (development[...])`)
- `R/harvest.R`, `R/anesrake.R`, `R/current_miss.R`, `R/diagnose_weights.R`,
  `R/weighted_pct.R` - non-ASCII Unicode characters in roxygen/inline comments replaced with
  ASCII equivalents
- `man/harvest.Rd` - hand-mirrored the same non-ASCII substitution (not regenerated via
  roxygen2, to avoid an installed-version-mismatch diff — see Decisions)
- 11 root-level files removed from git tracking (see plan frontmatter for the full list)
- `leafblower.Rcheck/` (120 tracked files) removed from git tracking

## Decisions Made

See `key-decisions` in frontmatter for the full list with rationale. Summary: extended
`.Rbuildignore` well past the plan's declared 14 patterns once the real build run showed
what was actually landing in the tarball; fixed the LaTeX/non-ASCII/NEWS/Suggests issues at
the root cause rather than suppressing them; explicitly declined to install missing
check-environment tools (checkbashisms/tidy/V8) as out of Rule-3 scope, documenting the gap
instead; used env-var overrides (`_R_CHECK_FORCE_SUGGESTS_`, `R_MAKEVARS_USER`) to get a
check result representative of a stranger's machine rather than this one developer box.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 2 - Missing Critical] .Rbuildignore didn't cover untracked dev-tool cruft**
- **Found during:** Task 1, first real `R CMD build` run
- **Issue:** `R CMD build` scans the working tree, not git's index — 14MB of
  `graphify-out/cache/`, `.planning/` (1.8M, 134 tracked files never in the plan's stray
  list), `.tldr/`, `.gemini/`, `.gemini-bridge/`, `__pycache__/`, and a handful of
  untracked hidden agent-tool dirs/JS-tooling stubs (`.beads.gate.lock`, `.crush`, `.gsd`,
  `.lean-ctx`, `.pytest_cache`, `.serena`, `node_modules`, `package.json`,
  `package-lock.json`, `Cholesky`, `sor_corun_aa`, `weights`) were all landing in the
  tarball, some triggering "storing paths of more than 100 bytes is not portable" build
  warnings.
- **Fix:** Added 21 more `.Rbuildignore` patterns beyond the plan's 14 (`data-raw/` also
  added, matching standard R packaging convention for generator scripts).
- **Files modified:** `.Rbuildignore`
- **Verification:** Re-ran `R CMD build .`; zero "not portable" warnings; tarball shrank
  from 3.9MB to 344KB.
- **Committed in:** `4798d1b`

**2. [Rule 1 - Bug] Non-ASCII Unicode math symbols broke LaTeX PDF-manual generation**
- **Found during:** Task 1, second real `R CMD check --as-cran` run
- **Issue:** `checking PDF version of manual ... ERROR` — LaTeX errors on Sigma/approx/gte/
  Delta characters that roxygen carried from `R/harvest.R`'s doc comments into
  `man/harvest.Rd`; `checking code files for non-ASCII characters ... WARNING` on the same
  root cause, plus the same character classes found (via `tools::showNonASCIIfile`) in 4
  more R/*.R files.
- **Fix:** Replaced every non-ASCII character (em/en-dash, Sigma, Delta, `->`, `>=`, `||`,
  `~=`, `x`, `lambda`, `^2`, `.`, `=>`) with an ASCII equivalent across all 5 files;
  hand-mirrored the identical substitution into the already-generated `man/harvest.Rd`
  rather than re-running `roxygen2::roxygenize()` (installed 8.0.0 vs. DESCRIPTION's pinned
  `RoxygenNote: 7.3.3` — regenerating risked collateral formatting changes unrelated to this
  fix).
- **Files modified:** `R/harvest.R`, `R/anesrake.R`, `R/current_miss.R`,
  `R/diagnose_weights.R`, `R/weighted_pct.R`, `man/harvest.Rd`
- **Verification:** `tools::showNonASCIIfile()` returns empty for all 5 files;
  `checking PDF version of manual ... OK` and `checking code files for non-ASCII
  characters ... OK` on the next check run.
- **Committed in:** `4798d1b`

**3. [Rule 1 - Bug] NEWS.md unparseable — "No news entries found" NOTE**
- **Found during:** Task 1, second real check run
- **Issue:** `checking package subdirectories ... NOTE` — R's news parser
  (`tools:::.build_news_db_from_package_NEWS_md`) requires top-level headings to contain a
  valid version number; NEWS.md's three headings read `# leafblower (development)` /
  `(development version)`, matching no version, so the whole file parsed to zero entries.
- **Fix:** Retitled all three headings to `# leafblower 0.1.0` (this is the first release;
  content of every bullet under each heading left untouched).
- **Files modified:** `NEWS.md`
- **Verification:** `tools:::.build_news_db_from_package_NEWS_md("NEWS.md")` returns 8 rows
  (was NULL); `checking package subdirectories ... OK` on the next check run.
- **Committed in:** `4798d1b`

**4. [Rule 2 - Missing Critical] Undeclared test dependencies**
- **Found during:** Task 1, second real check run
- **Issue:** `checking for unstated dependencies in 'tests' ... WARNING` — `arrow`, `callr`,
  `jsonlite`, `withr` used across `tests/parity/` and several `tests/testthat/*.R` files but
  absent from `DESCRIPTION`'s `Suggests` field.
- **Fix:** Added all four to `Suggests`.
- **Files modified:** `DESCRIPTION`
- **Verification:** `checking for unstated dependencies in 'tests' ... OK` on the next check
  run.
- **Committed in:** `4798d1b`

**5. [Rule 3 - Blocking, environmental] -march=native NOTE traced to developer's own
~/.R/Makevars, not the package**
- **Found during:** Task 1, first real check run
- **Issue:** `checking compilation flags used ... NOTE` listed `-march=native -mtune=native`
  among other flags — neither `src/Makevars.in` nor `configure` sets these; grep confirmed
  they come from this developer's personal `~/.R/Makevars` (applies to every R package built
  on this machine).
- **Fix:** Re-ran the check with `R_MAKEVARS_USER=/dev/null`, the R-sanctioned mechanism to
  bypass a personal Makevars file for one invocation (does not modify `~/.R/Makevars`
  itself). This produced a check representative of a stranger's/CRAN's build machine rather
  than this one.
- **Files modified:** none (environment variable, not a file change)
- **Verification:** Re-run confirmed `-march=native`/`-mtune=native` gone from the NOTE;
  the remaining flags in that NOTE (`-march=x86-64`, `-Wformat`, etc.) were then traced via
  `R CMD config CXXFLAGS` to this Arch Linux R installation's own system `Makeconf` — also
  not package-controlled. `-mavx2` is the one package-set, load-bearing, already-documented
  flag in the list.
- **Committed in:** n/a (documented in `cran-comments.md`, `5bddbfa`)

**6. [Rule 3 - Blocking, environmental, NOT auto-fixed per Rule-3 exclusion] checkbashisms /
tidy / V8 absent from this local machine**
- **Found during:** Task 1, all three real check runs
- **Issue:** `checking top-level files ... WARNING` (`checkbashisms` script absent;
  `pandoc` initially also absent) and `checking HTML version of manual ... NOTE` (`tidy`
  and R's `V8` package absent) — these are optional tools `R CMD check` itself uses to
  validate/render, not package dependencies.
- **Fix:** Did NOT install `checkbashisms` (devscripts) or HTML Tidy — both would require a
  system package-manager install, explicitly excluded from Rule 3 auto-fix, and this
  session's sandbox independently blocks `pacman` invocation. Did point the check at an
  already-present `pandoc` binary bundled with a local Quarto install
  (`/opt/quarto/bin/tools/x86_64/pandoc`, added to `PATH` for one check invocation only —
  not a new install), which resolved the pandoc half of the top-level-files WARNING.
- **Files modified:** none
- **Verification:** Confirmed via `which`/`find` that `checkbashisms`/`tidy` genuinely don't
  exist anywhere reachable on this machine; documented as a known local-environment gap in
  `cran-comments.md`, expected to be closed by the 05-03/05-04 CI matrix (standard runner
  image, these tools present).
- **Committed in:** n/a (documented in `cran-comments.md`, `5bddbfa`)

---

**Total deviations:** 6 (4 auto-fixed via Rule 1/2, 2 documented-not-fixed via Rule 3's
package-install exclusion).
**Impact on plan:** All 4 auto-fixes are real correctness/portability fixes required for an
honest CRAN-representative check, not scope creep. The 2 documented-not-fixed items are the
reason Task 1's literal "0 warnings, at most 1 NOTE" acceptance criterion is not fully met
in this local sandboxed environment — see `## Known Stubs` below and coverage item D2's
`human_judgment: true` for the explicit gap.

## Known Stubs

None — no placeholder data or unwired UI. The one open item is a verification gap, not a
stub: Task 1's `<verify>` combined-command bar ("0 errors, 0 warnings, at most the
new-submission NOTE") is not literally met — the real, final local result is
`0 ERRORs, 1 WARNING, 3 NOTEs`, with 1 WARNING + 2 of 3 NOTEs traced to this machine's
absent check-tooling (checkbashisms, tidy, V8) rather than a package defect (see coverage
item D2 and deviation 6 above). `05-03`/`05-04`'s CI matrix is expected to produce the
literal result on a standard runner image.

## Issues Encountered

The plan's own `<verify>` automated command, run literally as written (no env var
overrides), returns `1 ERROR` (missing/outdated Suggests packages `PracTools`/`autumn`) —
this is a pre-existing local-machine package-library state, not something Task 1's action
list anticipated. `_R_CHECK_FORCE_SUGGESTS_=false` is R's own documented workaround
(printed directly in the ERROR message) and does not touch package code; used it to get a
representative check. See deviation 5 for the parallel `R_MAKEVARS_USER` finding.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

- The hygiene-clean tree and extended `.Rbuildignore` are the foundation 05-03 (CI) and
  05-04 (release automation) build on, per the plan's own stated purpose.
- `cran-comments.md` exists with real content; 05-05 (or whichever plan handles the actual
  submission) should re-run the check on the CI matrix's standard runner and update the
  "Test environments"/"R CMD check results" sections with that result once available, since
  it is expected to close the checkbashisms/tidy/V8 gaps this plan could not.
- A future plan should decide whether to also non-ASCII-clean `tests/testthat/*.R` (many
  files still carry Unicode math symbols in comments) — `R CMD check` does not scan `tests/`
  for this, so it is not currently blocking, but it is the same class of issue fixed in
  `R/*.R` here.

---
*Phase: 05-cran-pypi-release*
*Completed: 2026-08-15*
