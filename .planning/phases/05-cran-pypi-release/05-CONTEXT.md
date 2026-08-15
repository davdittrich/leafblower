# Phase 5: CRAN + PyPI Release - Context

**Gathered:** 2026-08-15
**Status:** Ready for planning

<domain>
## Phase Boundary

A survey analyst who has never seen this repository installs leafblower and gets a working,
self-contained package. This phase makes the package **installable and CI-verified** —
tarball hygiene, wheel build, version-sync check, and a CI matrix. It does NOT submit to
CRAN (deferred — see D-06) and does NOT write vignettes/user-facing documentation (deferred
— see Deferred Ideas).

</domain>

<decisions>
## Implementation Decisions

### Root strays (tracked dev artifacts)
- **D-01:** `cell_table_92c4f45.{cpp,hpp}`, `ieppa_92c4f45.cpp`, `patch_raking.py`,
  `patch_wolfe.py`, `test_output.log`, `leafblower_0.1.0.tar.gz`, `REVIEW_FINDINGS.md`,
  `code-review-findings.md` — confirmed **git-tracked** (not just untracked cruft;
  `leafblower.Rcheck/00_pkg_src/leafblower/patch_wolfe.py` is a stray build artifact and
  should not be committed either). Fix both ends: `git rm` them from tracking AND add
  `.Rbuildignore` patterns as defense-in-depth. — **Reversibility:** reversible — files
  remain in git history, `git rm` doesn't destroy anything, and `.Rbuildignore` patterns are
  trivially editable.

### CI matrix (KPI-06)
- **D-02:** No CI exists in this repo. Build a new GitHub Actions pipeline
  (`.github/workflows/`) rather than a documented manual matrix. The pipeline must build
  real, portable wheels via **cibuildwheel** across manylinux/macOS × Python 3.9–3.13 (not
  just a source-build smoke test) — this is what "pip install leafblower installs a
  self-contained wheel" in US-010/US-008 requires; source-build-only would verify
  compatibility but produce no distributable artifact. Also run `R CMD check --as-cran` in
  CI (still required — see D-06). — **Reversibility:** costly — a CI pipeline is durable
  infrastructure; ripping it out later means re-deciding the whole matrix strategy.

### Version sync (SC5)
- **D-03:** `DESCRIPTION` and `python/pyproject.toml` version drift is caught by a
  **test-suite assertion** (an R or Python test reads both files and fails on mismatch) —
  not a separate CI-only step, not a pre-commit hook. Runs on every DoD gate invocation
  (local and CI both), no new tooling. — **Reversibility:** reversible — a single test file.

### cran-comments.md content
- **D-04:** Write `cran-comments.md` even though CRAN submission itself is deferred (see
  D-06) — r-universe and any eventual CRAN submission both need it. Must explain, beyond
  new-submission boilerplate:
  - The deliberate **no `-O` flag** in `PKG_CXXFLAGS` (R supplies the optimization level via
    user/site Makevars; `tools:::.check_make_vars` rejects an `-O` flag in the package's own
    flags — preempt a reviewer flagging this as an oversight).
  - The **LAPACK hard dependency** (`find_package(... REQUIRED)` on the Python side; note
    the R-side `SystemRequirements`/`configure` behavior if relevant).
  — **Reversibility:** reversible — a markdown file, freely editable.

### Publication channel (scope-defining decision)
- **D-05:** Do **NOT** submit to CRAN this phase. Publish via **GitHub + r-universe** as an
  intermediate distribution channel instead. Vignettes and other documentation CRAN/general
  users would expect are not yet written — that's separate, substantial content work.
  — **Reversibility:** reversible — r-universe publication doesn't preclude a later CRAN
  submission; nothing about this choice is one-way.

- **D-06:** `R CMD check --as-cran` (0 errors/0 warnings, roadmap SC1) still applies this
  phase, unchanged. r-universe builds and checks packages the same way CRAN does, and it's
  the honest quality bar regardless of where the package is actually published. Only the
  CRAN web-form submission step itself is skipped. — **Reversibility:** reversible — running
  a stricter check than currently required costs nothing to later relax.

### Claude's Discretion
- Exact cibuildwheel config (which manylinux image, macOS runner versions, caching strategy)
  is an implementation detail for the researcher/planner to work out against current
  cibuildwheel docs — not discussed here.
- Whether the version-sync test lives in the R suite or Python suite is left to whichever
  the implementer finds cleaner to wire.

</decisions>

<canonical_refs>
## Canonical References

**Downstream agents MUST read these before planning or implementing.**

### Roadmap & requirements
- `.planning/ROADMAP.md` § Phase 5: CRAN + PyPI Release — success criteria, dependencies,
  bead IDs, and planning notes (LAPACK hard-required, `PKG_CXXFLAGS` no `-O`, no CI today,
  `.Rbuildignore` already excludes `cran-comments.md`, `tests/testthat/fixtures` is
  `.Rbuildignore`d so fixture-backed tests need `skip_if(!file.exists(...))` guards).
- `.planning/REQUIREMENTS.md` — US-010 (Distribution, Partial, "largest genuinely-unfinished
  requirement"), US-008 (Python wheel + 3.9-3.13 matrix, Partial), KPI-05 (CRAN check,
  Open), KPI-06 (Python wheel CI matrix, Open).
- `.planning/PROJECT.md` § Session Continuity, § Known Defects — current project state.

### Build configuration
- `.Rbuildignore` — current exclusion list; confirmed it does NOT yet exclude the 8 tracked
  root strays (D-01) and already excludes `cran-comments.md` (correct — CRAN wants that file
  outside the tarball, so create it, don't unexclude it).
- `python/CMakeLists.txt` — hard `find_package(LAPACK REQUIRED)`, `-O3` unconditional,
  explicit `CORE_SOURCES` list (does NOT glob `src/*.cpp` — new solver `.cpp` files must be
  added here manually).
- `DESCRIPTION` (Version: 0.1.0) / `python/pyproject.toml` (version = "0.1.0") — the two
  version strings D-03's test must compare.
- `CLAUDE.md` § Build & Test, § Architecture — `R CMD INSTALL --preclean .` is the R build
  gate; single-thread BLAS env vars required for deterministic parity/tests.

### Beads tickets
- `leafblower-l6h0` (P1) — `.Rbuildignore` gaps ship dev artifacts (D-01).
- `leafblower-dns3` (P2) — generated build artifacts tracked in git (D-01).
- `leafblower-kk1.24` (epic) / `leafblower-kk1.24.3` (T19 final gate) — CRAN submission
  package + PyPI wheel artifacts.

</canonical_refs>

<code_context>
## Existing Code Insights

### Reusable Assets
- `configure` — already has C++17→C++14 fallback logic; no changes needed for this phase.
- `pyproject.toml` — already on scikit-build-core; cibuildwheel builds on top of this
  without restructuring the Python build.

### Established Patterns
- Definition of Done gate (`CLAUDE.md`): `R CMD INSTALL --preclean .` + R testthat 0 FAIL +
  Python pytest 0 FAIL under single-thread BLAS. Any new CI workflow should invoke this same
  gate, not reinvent a parallel check.
- No `--cov-fail-under` / coverage gate is wired anywhere — don't add one as part of CI setup.

### Integration Points
- New GitHub Actions workflow(s) under `.github/workflows/` (currently absent — confirmed
  via `ls`).
- `.Rbuildignore` and a new/edited test file (version-sync assertion) are the two touch
  points besides the CI workflow itself.

</code_context>

<specifics>
## Specific Ideas

No specific UI/UX-style ideas — this is packaging/infra work. The concrete asks were: git-rm
the tracked strays, stand up cibuildwheel-based CI across Python 3.9-3.13, add a version-sync
test, write `cran-comments.md` explaining the no-`-O` and LAPACK points, and target
GitHub + r-universe rather than an actual CRAN submission this phase.

</specifics>

<deferred>
## Deferred Ideas

- **Vignettes and user-facing documentation** — needed before any real public release (CRAN
  or otherwise), but substantial content work. Deferred to its own future phase, not folded
  into Phase 5's packaging/distribution scope.
- **Actual CRAN submission** — `R CMD check --as-cran` still runs and must pass (D-06), but
  clicking submit on the CRAN web form is explicitly out of scope for this phase (D-05).
  File as a follow-up phase/ticket once docs/vignettes exist.

### Reviewed Todos (not folded)
None — `todo.match-phase 5` returned zero matches.

</deferred>

---

*Phase: 5-CRAN + PyPI Release*
*Context gathered: 2026-08-15*
