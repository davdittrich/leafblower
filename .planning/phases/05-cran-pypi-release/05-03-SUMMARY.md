---
phase: 05-cran-pypi-release
plan: 03
subsystem: infra
tags: [ci, github-actions, r-cmd-check, release-hygiene, regression-guard]

# Dependency graph
requires:
  - phase: 05-cran-pypi-release
    plan: 01
    provides: hygiene-clean git tree, extended .Rbuildignore, real local R CMD check --as-cran result to reproduce in CI
provides:
  - .github/workflows/r-check.yml — CI job reproducing 05-01's local R CMD build/check --as-cran
    cycle plus an ongoing tarball-hygiene regression guard
affects: [05-05]

# Actuals (#2632)
actuals:
  tokens: 357
  tasks: 2
  commits: 1

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Third-party GitHub Action supply chain pinned to full commit SHA (not floating tag),
      per RESEARCH.md's Known Threat Patterns table — same discipline as a flagged PyPI/npm
      package, confirmed via a checkpoint:human-verify before being written into the workflow"
    - "Ongoing regression guard implemented as a CI step that re-runs the exact local grep
      from 05-01's verify command, not a new/divergent check"

key-files:
  created:
    - .github/workflows/r-check.yml
  modified: []

key-decisions:
  - "Pinned all three r-lib/actions references (setup-r, setup-r-dependencies,
    check-r-package) to commit d3c5be51b12e724e68f33216ca3c148b66d5f0b6 (tag v2) after
    explicit human verification via the checkpoint:human-verify Task 0 — the human approved
    via the interactive AskUserQuestion UI, confirming the org/repo and SHA before any
    pinning happened."
  - "Hygiene guard step placed BEFORE the check step and copies 05-01's exact grep pattern
    verbatim (11 stray-file patterns + leafblower.Rcheck/) rather than reinventing an
    equivalent check, so CI and the proven local verification never drift apart."
  - "args: 'c(\"--as-cran\")' with no --no-manual, matching 05-01's full manual/vignette
    build exactly; error-on: '\"warning\"' so any WARNING fails the job (check-r-package's
    default of error-on error would let a WARNING through silently)."
  - "No Python setup step in this job — R and Python CI stay in separate jobs/files per
    RESEARCH.md Pattern 2 (05-04 owns the Python side)."

requirements-completed: [US-010, KPI-05, KPI-06]

coverage:
  - id: D1
    description: ".github/workflows/r-check.yml exists, is valid YAML, defines exactly one
      job r-check, and references all three r-lib/actions steps pinned to
      d3c5be51b12e724e68f33216ca3c148b66d5f0b6"
    requirement: US-010
    verification:
      - kind: automated
        ref: "Rscript -e 'wf <- yaml::read_yaml(\".github/workflows/r-check.yml\"); stopifnot(is.list(wf$jobs$\"r-check\")); cat(\"valid YAML, job present\\n\")' -> valid YAML, job present"
        status: pass
      - kind: automated
        ref: "grep -c d3c5be51b12e724e68f33216ca3c148b66d5f0b6 r-check.yml -> 3 references"
        status: pass
    human_judgment: false
  - id: D2
    description: "Check step reproduces 05-01's exact flags: --as-cran present, --no-manual
      absent, error-on set to warning (not the check-r-package default of error)"
    requirement: KPI-05
    verification:
      - kind: automated
        ref: "grep -q -- '--as-cran' r-check.yml && ! grep -q -- '--no-manual' r-check.yml && grep -q 'error-on:.*warning' r-check.yml -> structure OK"
        status: pass
    human_judgment: false
  - id: D3
    description: "Tarball hygiene guard step re-runs 05-01's stray-file grep as an ongoing
      regression check, running before the check step"
    requirement: KPI-06
    verification:
      - kind: automated
        ref: "grep -c 'leafblower\\\\.Rcheck' r-check.yml -> 2 references (guard step + trigger comment context)"
        status: pass
    human_judgment: false
---

# Phase 05 Plan 03: r-check CI Workflow Summary

**Added `.github/workflows/r-check.yml`, a GitHub Actions job that reproduces 05-01's
proven local `R CMD build . && R CMD check --as-cran` cycle on pinned r-lib/actions SHAs,
with an ongoing tarball-hygiene regression guard step running before the check.**

## Performance

- **Duration:** ~5 min (continuation from a checkpoint:human-verify approval)
- **Completed:** 2026-08-15T17:53:31Z
- **Tasks:** 2 (Task 0 checkpoint approved in prior agent turn; Task 1 authored + committed
  this turn; Task 2 structural verification this turn, no commit — read-verify only)
- **Files modified:** 1 (new file, matches plan's declared `files_modified`)

## Accomplishments

- Task 0 (checkpoint:human-verify, gate=blocking-human): human approved pinning
  r-lib/actions at commit `d3c5be51b12e724e68f33216ca3c148b66d5f0b6` (tag `v2`) via the
  interactive AskUserQuestion UI before any file was written.
- Task 1: Authored `.github/workflows/r-check.yml` — triggers on `push`/`pull_request`, one
  job `r-check` on `ubuntu-latest`: `actions/checkout@v4`, then
  `r-lib/actions/setup-r@d3c5be51...`, `r-lib/actions/setup-r-dependencies@d3c5be51...`
  (`extra-packages: any::rcmdcheck`, `needs: check`), a "Tarball hygiene guard" step
  reproducing 05-01's exact grep against the 11 named strays + `leafblower.Rcheck/`, then
  `r-lib/actions/check-r-package@d3c5be51...` with `args: 'c("--as-cran")'` (no
  `--no-manual`) and `error-on: '"warning"'`. Verified valid YAML with the `r-check` job
  present.
- Task 2: Structural verification pass (no source changes) confirmed `--as-cran` present,
  `--no-manual` absent, `error-on:` set to the literal `warning`, and both `push` and
  `pull_request` triggers present.

## Task Commits

1. **Task 1: Author .github/workflows/r-check.yml** - `74c26e5` (feat)
2. **Task 2: Structural verification of the r-check.yml workflow** - no commit (read-verify
   only, no source-file changes per the plan's own task action)

## Files Created/Modified

- `.github/workflows/r-check.yml` — new; 35 lines; CI job reproducing 05-01's local
  `R CMD check --as-cran` sequence with pinned r-lib/actions SHAs and a hygiene regression
  guard.

## Decisions Made

See `key-decisions` in frontmatter. Summary: full-SHA pinning of all three r-lib/actions
references only after explicit human verification (Task 0, gate=blocking-human, approved via
AskUserQuestion UI — not orchestrator self-approval); hygiene guard step copies 05-01's exact
grep pattern verbatim rather than a divergent equivalent; `--as-cran` with full manual build
and `error-on: warning` matching 05-01's local run exactly; R and Python CI intentionally kept
in separate jobs (05-04 owns Python).

## Deviations from Plan

None — plan executed exactly as written, including the human-verify checkpoint gate.

## Known Stubs

None.

## Issues Encountered

None.

## User Setup Required

None — no external service configuration required. Note: this repository has no git remote
(local-only per project convention), so the workflow's first real GitHub Actions execution
will happen once a remote exists (tracked in 05-05's phase-gate summary, per this plan's
`<verification>` section). Local structural/YAML verification is the full extent of what is
checkable in this environment.

## Next Phase Readiness

- `.github/workflows/r-check.yml` is authored and structurally verified; ready for 05-04
  (Python CI job) and 05-05 (phase-gate summary) to build on.
- No blockers for 05-04/05-05.

## Self-Check: PASSED

- FOUND: .github/workflows/r-check.yml
- FOUND: commit 74c26e5

---
*Phase: 05-cran-pypi-release*
*Completed: 2026-08-15*
