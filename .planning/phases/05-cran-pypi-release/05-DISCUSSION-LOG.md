# Phase 5: CRAN + PyPI Release - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-08-15
**Phase:** 5-CRAN + PyPI Release
**Areas discussed:** Root strays, CI matrix, Version sync, cran-comments content, Publication channel, CRAN check scope, Docs scope

---

## Root strays (tracked dev artifacts)

| Option | Description | Selected |
|--------|-------------|----------|
| git rm + .Rbuildignore both | Remove from git tracking AND add .Rbuildignore patterns as defense-in-depth | ✓ |
| .Rbuildignore only | Leave tracked in git, only exclude from tarball | |
| git rm only | Delete from git, no .Rbuildignore pattern | |

**User's choice:** git rm + .Rbuildignore both (Recommended option).
**Notes:** Confirmed via `git ls-files` that all 8 named strays are git-tracked, not just
untracked cruft — plus discovered a 9th, `leafblower.Rcheck/00_pkg_src/leafblower/patch_wolfe.py`.

---

## CI matrix (KPI-06)

| Option | Description | Selected |
|--------|-------------|----------|
| New GitHub Actions CI pipeline | .github/workflows/ running wheel build + import test across Python 3.9-3.13 | ✓ |
| Documented manual matrix | No CI infra, markdown doc records manual testing | |

**User's choice:** New GitHub Actions CI pipeline (Recommended option).
**Notes:** Confirmed via `ls .github/workflows/` that no CI exists at all. Later refined
(see "Wheel scope" below) to require cibuildwheel, not just a source-build smoke test.

---

## Version sync (SC5)

| Option | Description | Selected |
|--------|-------------|----------|
| Test-suite assertion | R or Python test reads both files, fails on mismatch | ✓ |
| CI-only check | Separate CI step compares version strings | |
| Pre-commit hook | Git hook blocks commit with only one version changed | |

**User's choice:** Test-suite assertion (Recommended option).
**Notes:** Confirmed both files currently read `0.1.0`.

---

## cran-comments content

| Option | Description | Selected |
|--------|-------------|----------|
| No -O flag by design | Note PKG_CXXFLAGS deliberately carries no -O level | ✓ |
| LAPACK hard dependency | Note SystemRequirements/configure LAPACK behavior | ✓ |
| Just new-submission boilerplate | Standard boilerplate only | |
| (free text) | Do not yet submit to CRAN — prepare only. Docs/vignettes needed first. Intermediate step: publish on GitHub + r-universe. | ✓ (drove D-05/D-06) |

**User's choice:** Multi-select — both "No -O flag" and "LAPACK hard dependency" content
requirements, PLUS a free-text answer that reframed the phase's publication scope entirely.
**Notes:** This answer surfaced two new gray areas (CRAN check scope, docs scope), asked as
a follow-up round below.

---

## Wheel scope

| Option | Description | Selected |
|--------|-------------|----------|
| cibuildwheel manylinux/macOS wheels | Real portable wheels via cibuildwheel across the matrix | ✓ |
| Source-build-only CI matrix | pytest against scikit-build-core install per Python version, no wheel artifact | |

**User's choice:** cibuildwheel manylinux/macOS wheels (Recommended option).
**Notes:** This is what "self-contained wheel" in US-010/US-008 literally requires.

---

## CRAN check scope

| Option | Description | Selected |
|--------|-------------|----------|
| Yes, still required | R CMD check --as-cran must still pass 0/0; only the web-form submission is skipped | ✓ |
| Relaxed check | Only plain R CMD check needed this phase | |

**User's choice:** Yes, still required (Recommended option).
**Notes:** Follow-up question triggered by the cran-comments free-text answer's scope
change (no CRAN submission this phase, r-universe instead).

---

## Docs scope

| Option | Description | Selected |
|--------|-------------|----------|
| Separate phase — file for later | Phase 5 stays focused on packaging/distribution mechanics | ✓ |
| In scope for Phase 5 | Vignettes/docs included in this phase's success criteria | |

**User's choice:** Separate phase — file for later (Recommended option).
**Notes:** Vignette/doc-writing is substantial content work — kept out to avoid scope creep
per the scope guardrail.

---

## Claude's Discretion

- Exact cibuildwheel configuration (manylinux image choice, macOS runner versions, caching).
- Whether the version-sync test lives in the R suite or the Python suite.

## Deferred Ideas

- Vignettes and user-facing documentation — needed before any real public release; own
  future phase.
- Actual CRAN web-form submission — deferred until docs/vignettes exist; R CMD check
  --as-cran still runs and must pass this phase, but nothing gets submitted.
