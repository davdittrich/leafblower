---
phase: 05-cran-pypi-release
plan: 02
subsystem: infra
tags: [python, pytest, version-sync, release-hygiene, regression-guard]

# Dependency graph
requires:
  - phase: 05-cran-pypi-release
    plan: 01
    provides: hygiene-clean git tree, extended .Rbuildignore, real R CMD check --as-cran result
provides:
  - python/leafblower/test_version_sync.py — regression test guarding DESCRIPTION vs.
    python/pyproject.toml version drift, runs in every local DoD gate and CI
affects: [05-03, 05-04, 05-05]

# Actuals (#2632)
actuals:
  tokens: 317
  tasks: 1
  commits: 1

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Build/config-drift guards live as plain pytest assertions (regex-extract two files,
      assert equality) rather than TOML-parsing dependencies or CI-only scripts — same
      pattern as test_core_sources_sync.py (SC4, phase 02-02)"

key-files:
  created:
    - python/leafblower/test_version_sync.py
  modified: []

key-decisions:
  - "Followed test_core_sources_sync.py's structure exactly (REPO_ROOT = parents[2], one
    test function, docstring shape citing ticket/SC/date, plain Path.read_text(), no
    mocking/fixtures) rather than inventing a new pattern for the same class of guard test"
  - "Plain regex extraction for both files' single-line version formats, no TOML-parsing
    dependency — avoids the tomllib (3.11+) vs tomli (<3.11) split per RESEARCH.md's Don't
    Hand-Roll table"
  - "No RED/GREEN split exercised: this is a regression guard against an already-in-sync
    state (both files read 0.1.0 today), not new production behavior — same shape as its
    analog test_core_sources_sync.py, which also ships already-passing. The plan's own
    <behavior> block frames the drift-detection path as 'reasoned, not separately
    exercised.' Single test(...) commit, no accompanying feat(...) commit."

requirements-completed: [US-010]

coverage:
  - id: D1
    description: "python/leafblower/test_version_sync.py exists, passes against the current in-sync 0.1.0/0.1.0 state, and is included in standard pytest collection with no special flag"
    requirement: US-010
    verification:
      - kind: automated
        ref: "cd python && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest leafblower/test_version_sync.py -v -> 1 passed"
        status: pass
    human_judgment: false

duration: 4min
completed: 2026-08-15
status: complete
---

# Phase 05 Plan 02: Version-Sync Regression Guard Summary

**Added `python/leafblower/test_version_sync.py`, a pytest assertion that regex-extracts
DESCRIPTION's `Version:` field and `python/pyproject.toml`'s `version = "..."` field and
fails loudly if they diverge — closing the phase's SC5 (version drift caught by a check,
not by a reader).**

## Performance

- **Duration:** 4 min
- **Started:** 2026-08-15T17:34:00Z (approx.)
- **Completed:** 2026-08-15T17:38:09Z
- **Tasks:** 1
- **Files modified:** 1 (new file, matches plan's declared `files_modified`)

## Accomplishments

- Wrote `python/leafblower/test_version_sync.py` following `test_core_sources_sync.py`'s
  exact structure: `REPO_ROOT = Path(__file__).resolve().parents[2]`, one test function
  `test_description_pyproject_version_match()`, docstring citing the ticket/SC/date and the
  verified baseline, plain `Path.read_text()` with no mocking or fixtures.
- Verified the regexes against both files' actual current single-line formats: `DESCRIPTION`
  line 3 (`Version: 0.1.0`), `python/pyproject.toml` line 10 (`version = "0.1.0"`).
- Ran the test: 1 collected, 1 passed.

## Task Commits

1. **Task 1: Add DESCRIPTION/pyproject.toml version-sync test** - `c4752e9` (test)

## Files Created/Modified

- `python/leafblower/test_version_sync.py` — new; 30 lines; regex-extracts and compares
  DESCRIPTION `Version:` vs. pyproject.toml `version = "..."`.

## Decisions Made

See `key-decisions` in frontmatter. Summary: mirrored the existing
`test_core_sources_sync.py` build/config-drift-guard pattern exactly rather than inventing a
new shape; used plain regex extraction (no TOML dependency); shipped as a single
`test(...)` commit since there is no production code to implement for this pure
regression-guard test — the plan's own `<behavior>` block frames the drift-detection path as
reasoned, not separately exercised.

## Deviations from Plan

None — plan executed exactly as written.

## Known Stubs

None.

## Issues Encountered

None.

## User Setup Required

None — no external service configuration required.

## Next Phase Readiness

- The version-sync guard is now part of the standard pytest collection, so 05-03 (CI matrix)
  will exercise it automatically on every run without additional wiring.
- No blockers for 05-03/05-04/05-05.

---
*Phase: 05-cran-pypi-release*
*Completed: 2026-08-15*
