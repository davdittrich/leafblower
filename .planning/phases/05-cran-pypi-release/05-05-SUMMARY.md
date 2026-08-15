---
phase: 05-cran-pypi-release
plan: 05
subsystem: infra
tags: [cran, pypi, github-actions, ci, release-hygiene, phase-gate, git-remote]

# Dependency graph
requires:
  - phase: 05-cran-pypi-release
    plan: 01
    provides: hygiene-clean git tree, extended .Rbuildignore, real local R CMD check --as-cran result, cran-comments.md
  - phase: 05-cran-pypi-release
    plan: 02
    provides: python/leafblower/test_version_sync.py (DESCRIPTION vs pyproject.toml drift guard)
  - phase: 05-cran-pypi-release
    plan: 03
    provides: .github/workflows/r-check.yml (authored, structurally verified, never executed)
  - phase: 05-cran-pypi-release
    plan: 04
    provides: .github/workflows/python-wheels.yml + [tool.cibuildwheel] config (authored, structurally verified, never executed) plus the sdist force-include fix
provides:
  - "A real GitHub remote (https://github.com/davdittrich/leafblower, public) — the phase's
    prior CI-only gap closed by actually pushing, not by a local-only sign-off"
  - "Both .github/workflows/r-check.yml and .github/workflows/python-wheels.yml executed for
    real on GitHub Actions and are green: r-check 0 errors/0 warnings/2 NOTEs, 1744 tests
    pass; python-wheels builds+twine-checks+imports cleanly on ubuntu-latest + macos-14
    across Python 3.9-3.13"
  - "A real test-code bug fix (unguarded library(DiceKriging) in test-algo-selection.R) found
    only by real CI execution, not any prior local run"
  - "cran-comments.md updated with the real, final CI outcome (not a projected/expected one)"
affects: []

# Actuals (#2632)
actuals:
  tokens: 2263
  tasks: 3
  commits: 10

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "Real CI execution surfaces defects no amount of local structural verification can:
      an uninstallable Suggests package, a pandoc auto-detection quirk, an implicit
      testthat extra-package assumption, an unguarded library() call in a benchmark-sourcing
      test file, and a runner that silently never schedules (macos-13) — none of these were
      visible from 05-01 through 05-04's local-only proofs"
    - "A checkpoint's two offered options are not exhaustive: the human chose a third path
      (set up the git remote now) that produced strictly better evidence than either
      'accept the CI-only gap' or 'defer to a follow-up ticket'"

key-files:
  created: []
  modified:
    - .gitignore
    - .github/workflows/r-check.yml
    - .github/workflows/python-wheels.yml
    - tests/testthat/test-algo-selection.R
    - cran-comments.md

key-decisions:
  - "User overrode the plan's Task 3 checkpoint (which offered only 'approve local-only
    closure' or 'name what should happen next') by choosing to set up the GitHub remote
    immediately, converting SC3/SC4's stated CI-only limitation into real, green CI evidence
    within the same session rather than deferring it"
  - "Per project CLAUDE.md's pre-push-to-new-remote hygiene audit, untracked .wolf/, .beads/,
    .claude/, .metaswarm/, .tldr/ before the first push (194535d) rather than pushing them to
    a now-public repository"
  - "Dropped macos-13 from the python-wheels CI matrix (bc73cb0) after the runner failed to
    schedule across 25+ minutes on every CI run on this account — consistent with GitHub's
    documented phase-out of Intel macOS runners, not a leafblower-side defect. User-approved
    scope reduction: macos-14 (arm64) stays, x86_64-macOS wheel coverage is no longer proven
    by CI"
  - "Fixed the r-check.yml autumn/testthat/pandoc CI failures (fc1fd3c, 82bb159, c1986c2)
    without weakening what the check enforces — each fix targets a real
    CI-environment-vs-local-environment gap (autumn unresolvable from any repo available to
    the runner; check-r-package's pandoc auto-detection vs. explicit r-lib/actions/setup-
    pandoc; testthat needing to stay an explicit extra-package despite already being a hard
    Imports dependency)"
  - "Installed TinyTeX via r-lib/actions/setup-tinytex (12c803f) so CI's R CMD check builds
    the real PDF manual rather than skip it — user explicitly chose this over --no-manual,
    matching 05-01's local proof exactly (no --no-manual flag anywhere in the check
    invocation)"
  - "Marked WINDOWS.md ledger entry #1 (05-01's 1 WARNING/3 NOTE local result, open,
    'expected to close on CI matrix') as fixed — real CI now returns 0 errors/0 warnings/
    2 NOTEs, both explained in cran-comments.md, satisfying SC1's literal clause ('at most
    the new-submission note... cran-comments.md explains any remaining note') in spirit even
    though the 2 remaining NOTEs (-mavx2 flags, HTML-tidy-absent) are not literally the
    new-submission NOTE"

requirements-completed: [US-010, US-008, KPI-05, KPI-06]

coverage:
  - id: D1
    description: "Full local DoD gate green: R CMD INSTALL --preclean . + devtools::test() (0 FAIL), uv pip install -e . + single-thread-BLAS pytest (0 failed), test_version_sync.py present and passing in the pytest collection"
    requirement: US-010
    verification:
      - kind: automated
        ref: "R CMD INSTALL --preclean . && Rscript -e 'devtools::test()' -> 0 FAIL; cd python && uv pip install -e . --reinstall-package leafblower && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest -v -> 0 failed, test_version_sync.py::test_description_pyproject_version_match passed"
        status: pass
    human_judgment: false
  - id: D2
    description: "Combined final tree re-confirmed hygiene-clean and CRAN-check-clean after 05-02/05-03/05-04's edits: hygiene grep 0, both workflow YAMLs valid, cran-comments.md present"
    requirement: KPI-05
    verification:
      - kind: automated
        ref: "git ls-files hygiene grep -> 0; Rscript -e 'yaml::read_yaml(\".github/workflows/r-check.yml\"); yaml::read_yaml(\".github/workflows/python-wheels.yml\")' -> both valid YAML; test -f cran-comments.md -> present"
        status: pass
    human_judgment: false
  - id: D3
    description: "r-check.yml executed for real on GitHub Actions (not merely authored/structurally verified) and is green: 0 errors, 0 warnings, 2 NOTEs (-mavx2 compilation flags — expected/documented; HTML-manual tidy missing — environment gap, harmless), all 1744 tests pass"
    requirement: KPI-05
    verification:
      - kind: e2e
        ref: "https://github.com/davdittrich/leafblower/actions/runs/31908234869"
        status: pass
    human_judgment: false
  - id: D4
    description: "python-wheels.yml executed for real on GitHub Actions across ubuntu-latest + macos-14 (arm64), Python 3.9-3.13: wheels build, pass twine check, import + calibrate cleanly on all 5 versions on both platforms. macos-13 (x86_64) dropped from the matrix — runner never scheduled across 25+ min on every run, a GitHub-side Intel-macOS-runner phase-out, not a package defect; x86_64-macOS wheel coverage is genuinely unproven"
    requirement: KPI-06
    verification:
      - kind: e2e
        ref: "https://github.com/davdittrich/leafblower/actions/runs/31908234870"
        status: pass
    human_judgment: false
  - id: D5
    description: "Human sign-off on phase closure — resolved not by approving the checkpoint's stated CI-only gap, but by directing the orchestrator to set up a real git remote and run real CI to close the gap outright, which then happened and is green"
    verification: []
    human_judgment: true
    rationale: "A checkpoint gate=blocking requires the user's own decision; the user's decision here was a third option not offered by the checkpoint's two stated choices, and produced strictly stronger evidence (real green CI) than either offered option would have."

duration: ~65min (Task 1/2 verification + checkpoint resolution + post-checkpoint remote setup, hygiene audit, and 6 CI-fix commits)
completed: 2026-08-15
status: complete
---

# Phase 5 Plan 5: CRAN + PyPI Release Phase Gate Summary

**Phase gate closed with real green CI, not a documented local-only limitation: pushed to a
new GitHub remote, untracked dev-tooling artifacts first per CLAUDE.md's pre-push audit, then
iteratively fixed six real CI-only defects (an unresolvable Suggests dependency, a pandoc
auto-detection quirk, an implicit testthat extra-package assumption, an unguarded
`library(DiceKriging)` in a benchmark-sourcing test file, a missing PDF-manual toolchain, and
a chronically-unscheduled macos-13 runner) until both `.github/workflows/r-check.yml` (0
errors/0 warnings/2 NOTEs, 1744 tests) and `.github/workflows/python-wheels.yml` (wheels
build+twine-check+import clean on ubuntu-latest + macos-14, Python 3.9-3.13) ran for real and
turned green.**

## Performance

- **Duration:** ~65 min total (Task 1 DoD gate + Task 2 hygiene re-confirmation, both
  verification-only; Task 3 checkpoint raised and resolved; then ~45 min of orchestrator-
  driven remote setup, hygiene audit, and 6 iterative CI-fix commits after the human's
  "set up a git remote now" decision)
- **Tasks:** 3 (Task 1 and 2 automated verification, Task 3 checkpoint resolved by the user
  choosing a path not offered by either of the checkpoint's two stated options)
- **Files modified (post-checkpoint work only — Tasks 1/2 modified nothing, per plan):** 5
  (`.gitignore`, `.github/workflows/r-check.yml`, `.github/workflows/python-wheels.yml`,
  `tests/testthat/test-algo-selection.R`, `cran-comments.md`)

## Accomplishments

- **Task 1 (Full local DoD gate):** Ran the project's exact Definition-of-Done gate from
  `CLAUDE.md` § Build & Test. R suite: 0 FAIL. Python pytest (single-thread BLAS): 0 failed,
  with `test_version_sync.py::test_description_pyproject_version_match` present and passing
  in the collection, confirming 05-02's guard is live in the standard test run rather than
  orphaned.
- **Task 2 (Re-confirm hygiene and CRAN-check-clean):** Re-ran 05-01's combined hygiene grep
  and CRAN-check verification against the final tree after 05-02/05-03/05-04's edits — 0
  tracked strays, both `.github/workflows/*.yml` files parse as valid YAML, `cran-comments.md`
  present with real (non-placeholder) content.
- **Task 3 (checkpoint, resolved):** The plan's checkpoint offered exactly two paths — approve
  local-only closure with the CI-only gap disclosed, or name a follow-up. The user chose a
  third: set up the git remote immediately. This closed the gap for real rather than deferring
  it:
  1. Audited tracked dev-tooling artifacts per `CLAUDE.md`'s pre-push-to-new-remote rule;
     untracked `.wolf/`, `.beads/`, `.claude/`, `.metaswarm/`, `.tldr/` and added them to
     `.gitignore` (`194535d`) before the repository became public.
  2. Created `https://github.com/davdittrich/leafblower` (public) via `gh repo create`,
     pushed `master` and all tags.
  3. Both workflows actually ran on GitHub Actions. Iteratively fixed six real,
     CI-environment-only defects, each committed separately:
     - `fc1fd3c`, `82bb159`, `c1986c2`: `autumn` (an uninstallable-on-any-repo Suggests
       package) resolution failures, a pandoc auto-detection quirk, and restoring `testthat`
       as an explicit `extra-packages` entry for `r-lib/actions/setup-r-dependencies` — none
       weaken what the check enforces; each targets a real CI-environment-vs-local gap.
     - `61fa6e0`: a genuine test-code bug — `tests/testthat/test-algo-selection.R` sourced a
       benchmark file with an unguarded `library(DiceKriging)`; added
       `skip_if_not_installed()` guards matching the project's established pattern, visible
       only because a runner without `DiceKriging` installed actually executed this file.
     - `12c803f`: installed TinyTeX (`r-lib/actions/setup-tinytex`) so `R CMD check` builds
       the real PDF manual on `ubuntu-latest` (no `pdflatex` by default) — user chose this
       over `--no-manual`, keeping CI's check identical in scope to 05-01's local proof.
     - `bc73cb0`: dropped `macos-13` from the wheel matrix — the runner never scheduled after
       25+ minutes across every CI run on this account, consistent with GitHub phasing out
       Intel macOS runners; kept `ubuntu-latest` + `macos-14` (arm64). User-approved scope
       reduction, honestly reflected in `cran-comments.md` and in this summary's SC3/SC4
       coverage below.
     - `2ac7870`: updated `cran-comments.md` with the real, final CI outcome, replacing the
       "expected to close" language from 05-01 with the actual measured result.

## Task Commits

1. **Task 1: Full local DoD gate** — verification-only, no commit (`22c3cde` from the prior
   session's continuation was the last code fix ahead of this run)
2. **Task 2: Re-confirm hygiene and CRAN-check-clean** — verification-only, no commit
3. **Task 3 follow-on (post-checkpoint):**
   - `194535d` (chore) — untrack local dev-tooling artifacts before first remote push
   - `fc1fd3c` (fix) — skip uninstallable Suggests (`autumn`) in r-check workflow
   - `82bb159` (fix) — skip pandoc-need auto-detection in r-check workflow
   - `c1986c2` (fix) — restore testthat as explicit extra-package in r-check
   - `61fa6e0` (fix) — `skip_if_not_installed()` guard for algo-selection bench source
   - `12c803f` (fix) — install TinyTeX so R CMD check can build the PDF manual
   - `bc73cb0` (fix) — drop macos-13 from wheel matrix (runner never schedules)
   - `2ac7870` (docs) — record real CI outcome in cran-comments.md
   - `26eb97c` (chore) — trigger Actions workflow indexing on first push

## Files Created/Modified

- `.gitignore` — added `.wolf/`, `.beads/`, `.claude/`, `.metaswarm/`, `.tldr/` before the
  first push to a public remote
- `.github/workflows/r-check.yml` — 7-line diff: skip unresolvable `autumn` Suggests, drop
  pandoc auto-detection reliance, restore explicit `testthat` extra-package, add TinyTeX setup
- `.github/workflows/python-wheels.yml` — 1-line diff: `macos-13` removed from the build
  matrix
- `tests/testthat/test-algo-selection.R` — `skip_if_not_installed("DiceKriging")` guard added
  around the benchmark-file source that pulled in an unguarded `library(DiceKriging)`
- `cran-comments.md` — Test environments and R CMD check results sections rewritten with the
  real CI run URLs and the actual 0-error/0-warning/2-NOTE outcome, replacing 05-01's
  "expected to close" projection

## Decisions Made

See `key-decisions` in frontmatter. Summary: the user's choice to set up the git remote
immediately, rather than accept either checkpoint option, converted a documented limitation
into closed, verified scope within the same session. All six CI-fix commits target real
environment gaps or a real test bug, none weaken what either workflow enforces. `macos-13` was
dropped after empirical evidence (25+ minutes of non-scheduling across every run) rather than
guessed at.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Unguarded `library(DiceKriging)` in test-algo-selection.R**
- **Found during:** Post-checkpoint CI run of `r-check.yml` (real GitHub Actions execution;
  never surfaced by any of this project's prior local runs, since `DiceKriging` happened to
  be installed on the developer's own machine)
- **Issue:** `tests/testthat/test-algo-selection.R` sourced a benchmark file containing an
  unconditional `library(DiceKriging)` call; the CI runner does not have `DiceKriging`
  installed (it is an optional `Suggests`-only comparison package), so the test file errored
  at source-time rather than skipping gracefully.
- **Fix:** Added `skip_if_not_installed("DiceKriging")` guards matching the project's
  established pattern used elsewhere for the same class of optional dependency.
- **Files modified:** `tests/testthat/test-algo-selection.R`
- **Commit:** `61fa6e0`

**2. [Rule 3 - Blocking, environmental] `autumn` Suggests package unresolvable on any repo
reachable by the CI runner**
- **Found during:** Post-checkpoint CI run of `r-check.yml`
- **Issue:** `r-lib/actions/setup-r-dependencies`'s default dependency-resolution step tried
  to install `autumn` (an optional `Suggests` comparison package, itself not on CRAN) and
  failed — no repository configured on the runner can resolve it.
- **Fix:** Adjusted the workflow's dependency step to skip `autumn` explicitly rather than
  let the whole dependency-install step fail, while keeping `testthat` (a hard runtime
  dependency of the test suite) as an explicit extra-package so it is never silently dropped
  along with the optional Suggests it was bundled with.
- **Files modified:** `.github/workflows/r-check.yml`
- **Commits:** `fc1fd3c`, `82bb159`, `c1986c2` (three iterations to isolate the correct fix
  without weakening the check)

**3. [Rule 2 - Missing Critical] R CMD check's PDF-manual build has no LaTeX toolchain on
`ubuntu-latest` by default**
- **Found during:** Post-checkpoint CI run of `r-check.yml`
- **Issue:** `ubuntu-latest` ships no `pdflatex`; `R CMD check --as-cran` (no `--no-manual`,
  matching 05-01's local proof) therefore cannot build the manual it is asked to build.
- **Fix:** Added `r-lib/actions/setup-tinytex` before the check step, per the user's explicit
  choice of "install the toolchain" over "weaken the check with `--no-manual`."
- **Files modified:** `.github/workflows/r-check.yml`
- **Commit:** `12c803f`

**4. [Rule 4 - Architectural, user-approved] macos-13 runner never schedules**
- **Found during:** Post-checkpoint CI run of `python-wheels.yml`, repeated across every run
  on this account
- **Issue:** `macos-13` (Intel) never entered the queue in 25+ minutes across multiple full
  CI runs, while `ubuntu-latest` and `macos-14` both ran immediately every time — consistent
  with GitHub's documented phase-out of Intel macOS hosted runners on this account tier, not
  a defect in this project's workflow configuration.
- **Change:** Dropped `macos-13` from the wheel-build matrix; `ubuntu-latest` + `macos-14`
  (arm64) remain. This is a genuine scope reduction — x86_64 macOS wheel coverage is no
  longer proven by CI — presented to and approved by the user rather than silently dropped.
- **Files modified:** `.github/workflows/python-wheels.yml`
- **Commit:** `bc73cb0`

---

**Total deviations:** 4 (1 Rule 1 bug fix, 1 Rule 2 missing-critical fix, 1 Rule 3
environmental blocking fix requiring 3 iterations, 1 Rule 4 architectural/scope change,
explicitly presented to and approved by the user).
**Impact on plan:** All four are real defects or real scope decisions surfaced only by
actual CI execution — none is scope creep, and none weakens what either workflow checks.

## Known Stubs

None. One honestly-recorded coverage gap, not a stub: **x86_64 macOS (`macos-13`) wheel
coverage is unproven.** The CI matrix builds and verifies `ubuntu-latest` (manylinux) and
`macos-14` (arm64) across Python 3.9-3.13; Intel-macOS wheel-building mechanics were proven
once, locally, without Docker, in 05-04 (auditwheel/MKL vendoring), but never on CI and never
on the specific `macos-13` runner image, since it never scheduled. This affects US-010's and
KPI-06's exact wording — see Requirements Impact below.

## Issues Encountered

None beyond the six items already documented as deviations above — all were resolved within
this session.

## User Setup Required

None further. The one item this phase previously flagged as requiring a user decision — "set
up a git remote" — was resolved by the user during this session's checkpoint.

## Threat Flags

None new. The plan's own `T-05-09` (repudiation risk of silently claiming CI-matrix coverage
that was never executed) is now closed by the opposite outcome: CI genuinely ran, and this
summary states the one remaining coverage gap (`macos-13`/x86_64-macOS) explicitly rather
than folding it into a blanket "CI is green" claim.

## Requirements Impact

- **US-010** (CRAN + PyPI distributable): moves from **Partial** toward **Implemented** in
  substance — `cran-comments.md` documents a real CI run (0 errors/0 warnings/2 NOTEs, both
  explained), the hygiene-clean tarball is CI-reproduced, and the version-sync guard runs in
  every CI invocation. **REQUIREMENTS.md's literal text should be updated** to reflect the
  real CI outcome (see `.planning/REQUIREMENTS.md` traceability note below) — this summary
  records the evidence; REQUIREMENTS.md itself is updated in the same commit as this file.
- **US-008** (Python `harvest()` via compiled core): the previously-stated gap — "no wheel
  artefact has ever been built... no CI exists" — is now factually superseded: wheels build,
  pass `twine check`, and import+calibrate cleanly on `ubuntu-latest` + `macos-14` across
  Python 3.9-3.13, for real, on CI. **Residual, honestly stated:** `macos-13` (x86_64 macOS)
  is not part of the matrix that actually ran — US-008's underlying PRD text does not name a
  specific macOS architecture, so this does not block calling US-008's CI-matrix clause
  closed, but it is the one platform combination this phase cannot claim to have proven.
- **KPI-05** (CRAN check, 0 errors/0 warnings via `R CMD check --as-cran`): the CI run
  (`https://github.com/davdittrich/leafblower/actions/runs/31908234869`) reports **0 errors,
  0 warnings, 2 NOTEs** — literally satisfies the "0 errors and 0 warnings" clause. The 2
  NOTEs (`-mavx2` compilation flags, HTML-manual `tidy` absent) are not literally "the
  new-submission note" SC1 anticipated, but SC1's own text allows this: "`cran-comments.md`
  explains any remaining note" — both are explained there. **KPI-05 moves from Open to
  Implemented.**
- **KPI-06** (Python wheel installs on Linux/macOS, Python 3.9-3.13, via a CI matrix): the CI
  run (`https://github.com/davdittrich/leafblower/actions/runs/31908234870`) proves this for
  Linux (`ubuntu-latest`, manylinux-tagged) and macOS arm64 (`macos-14`) across all 5 targeted
  Python versions. **Residual, stated honestly:** macOS x86_64 (`macos-13`) coverage is
  unproven by CI — the runner never scheduled. KPI-06's literal text says "Linux/macOS"
  without naming an architecture, so the Linux+arm64-macOS evidence satisfies it as written;
  this summary and `cran-comments.md` both flag the x86_64 gap so it is not silently
  generalized into "macOS" without qualification. **KPI-06 moves from Open to Implemented**,
  with this residual gap recorded in REQUIREMENTS.md rather than omitted.

## Next Phase Readiness

- Phase 5 (CRAN + PyPI Release) is complete — this was the final phase in the v1.0 roadmap
  (`.planning/ROADMAP.md` lists Phases 1-5, all now `[x]`).
- No blocking follow-up work identified by this plan. The one recorded residual gap
  (`macos-13`/x86_64-macOS CI coverage) is documented here and in `cran-comments.md`, not
  filed as a new ticket — it reflects a GitHub-side runner-availability constraint outside
  this project's control, not an open defect to fix.
- `bd` (beads) is currently schema-mismatched on this machine (database v65 vs. this binary's
  v53) and could not be queried during this session to close/verify Phase 5's tracked
  tickets (`leafblower-l6h0`, `leafblower-dns3`, `leafblower-kk1.24`, `leafblower-kk1.24.3`).
  This is a pre-existing environment issue, not something introduced by this plan — flagged
  here so a future session runs the documented recovery (`BD_IGNORE_SCHEMA_SKEW=1` or the
  v1.2.2 recovery guide) before assuming those tickets are still open.

## Self-Check: PASSED

- FOUND: `.github/workflows/r-check.yml`
- FOUND: `.github/workflows/python-wheels.yml`
- FOUND: `cran-comments.md`
- FOUND: commit `194535d`
- FOUND: commit `fc1fd3c`
- FOUND: commit `82bb159`
- FOUND: commit `c1986c2`
- FOUND: commit `61fa6e0`
- FOUND: commit `12c803f`
- FOUND: commit `bc73cb0`
- FOUND: commit `2ac7870`

---
*Phase: 05-cran-pypi-release*
*Completed: 2026-08-15*
