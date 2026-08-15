# Phase 5: CRAN + PyPI Release - Research

**Researched:** 2026-08-15
**Domain:** R package (CRAN-style) hygiene/check + Python wheel packaging (cibuildwheel) + CI (GitHub Actions) + version-drift prevention
**Confidence:** HIGH (tarball/tracked-file facts, DESCRIPTION/CMakeLists/Makevars contents — all read directly this session) / MEDIUM (cibuildwheel/CRAN process best practice — WebSearch against official docs, current as of 2026-08-15)

## Summary

This phase has two independently-verifiable deliverables (R tarball hygiene + `R CMD check
--as-cran` cleanliness; a self-contained multi-version Python wheel built in CI) plus one
small cross-cutting guard (version-sync test). None of it is exploratory — CONTEXT.md has
already locked the mechanism for four of five decisions (D-01 through D-04); this research
verifies the *current, exact state* those decisions must operate against, and fills in the
one area CONTEXT.md left to research/planner discretion: the concrete cibuildwheel
configuration.

The tarball-hygiene picture is worse than CONTEXT.md's D-01 already knew. Beyond the 8 named
strays plus the discovered `leafblower.Rcheck/00_pkg_src/leafblower/patch_wolfe.py`, this
session found the **entire `leafblower.Rcheck/` directory tree is git-tracked** — roughly 90
files including compiled `.so`/`.o` binaries, `.rdb`/`.rdx` compiled help databases, and a
full second copy of `tests/testthat/` frozen from a much earlier state of the package
(`test-lbfgsb.R`, `test-ieppa.R` — names that no longer exist in the live `tests/testthat/`
tree, confirming this is stale build output, not source). Separately, `git ls-files` at repo
root surfaces five more tracked files with no `.Rbuildignore` pattern covering them:
`baseline_bench.R`, `security-review-2026-04-27-ieppa.json`, `GEMINI.md`, `.mcp.json`,
`.coverage-thresholds.json`. R's build system does **not** blanket-exclude dotfiles or
arbitrary root files by default — only VCS/editor artifacts are excluded automatically — so
every one of these ships in the source tarball today. One root file, `cleanup`, that looks
stray at first glance is in fact a standard, recognized R package special file (a POSIX
shell script R runs after build/check to remove transient files) and must NOT be removed or
`.Rbuildignore`d.

The Python wheel side has a hard, unconditional `find_package(LAPACK REQUIRED)`
(`python/CMakeLists.txt:93`) with no fallback path — "self-contained" can only mean "LAPACK
vendored into the wheel by the platform-appropriate repair tool," not "no LAPACK dependency."
cibuildwheel's built-in repair step (auditwheel on Linux, delocate on macOS) does exactly
this, but the *source* of LAPACK differs sharply by platform: Linux manylinux containers need
`lapack-devel`/`openblas-devel` installed via `CIBW_BEFORE_ALL_LINUX` before the build even
starts (no LAPACK ships in the base image); macOS gets LAPACK for free from the OS-provided
Accelerate framework, but CMake's `FindLAPACK` module has a long history of not
auto-detecting Accelerate without an explicit `BLA_VENDOR=Apple` hint. Windows was not
selected in CONTEXT.md's CI-matrix decision (D-02 says "manylinux/macOS") and this is
correct to leave out of scope: `find_package(LAPACK REQUIRED)` has no free LAPACK source on
Windows (would need vcpkg/conda), so a Windows wheel is materially more work than the other
two platforms.

**Primary recommendation:** Use cibuildwheel (PyPA's own tool, GitHub-Actions-native, drives
manylinux/macOS wheel builds + auditwheel/delocate repair automatically) building on the
already-scikit-build-core `python/pyproject.toml` with no restructuring; gate Linux LAPACK
via `CIBW_BEFORE_ALL_LINUX = "yum install -y lapack-devel openblas-devel"` (manylinux2014) and
macOS LAPACK via `CIBW_ENVIRONMENT_MACOS` hints to force Accelerate; run `R CMD check
--as-cran` as a second, separate CI job via `r-lib/actions/setup-r` + `rcmdcheck`, not folded
into the same job as the Python build; implement the D-03 version-sync test as a regex read
of both `DESCRIPTION`'s `Version:` line and `pyproject.toml`'s `version = "..."` line — no
new dependency, no `tomllib`/`tomli` Python-version gating concern since the project supports
3.9 (pre-`tomllib`).

## Architectural Responsibility Map

This phase has no browser/API/DB tiers — it is packaging/build infrastructure. The
equivalent ownership split:

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Source tarball hygiene (strip dev artifacts) | git tracking (`git rm`) | `.Rbuildignore` (defense-in-depth) | D-01: fix both ends — git rm is the actual fix, `.Rbuildignore` guards against future re-adds |
| `R CMD check --as-cran` cleanliness | R package build (`DESCRIPTION`/`Makevars.in`/`configure`) | CI (r-lib/actions) | Check runs against the *built* tarball; CI just automates re-running it on every push |
| Python wheel self-containment | `python/CMakeLists.txt` (LAPACK link) | cibuildwheel repair step (auditwheel/delocate) | CMake finds+links LAPACK at build time; the repair step is what makes the *artifact* self-contained by vendoring the shared lib in |
| Multi-version (3.9-3.13) matrix | CI (`.github/workflows/`, cibuildwheel `CIBW_BUILD`) | — | No manual matrix chosen (D-02) — CI is the sole owner |
| Version-drift prevention | Test suite (R or Python, per D-03) | — | Deliberately NOT a CI-only or pre-commit mechanism — runs on every local DoD gate too |
| `cran-comments.md` content | Documentation (markdown, hand-written) | — | Not generated; explains build-flag/LAPACK decisions to a human CRAN/r-universe reviewer |

## User Constraints

<user_constraints>
### Locked Decisions (from CONTEXT.md — copied verbatim)

- **D-01:** `cell_table_92c4f45.{cpp,hpp}`, `ieppa_92c4f45.cpp`, `patch_raking.py`,
  `patch_wolfe.py`, `test_output.log`, `leafblower_0.1.0.tar.gz`, `REVIEW_FINDINGS.md`,
  `code-review-findings.md` — confirmed **git-tracked** (not just untracked cruft;
  `leafblower.Rcheck/00_pkg_src/leafblower/patch_wolfe.py` is a stray build artifact and
  should not be committed either). Fix both ends: `git rm` them from tracking AND add
  `.Rbuildignore` patterns as defense-in-depth.

- **D-02:** No CI exists in this repo. Build a new GitHub Actions pipeline
  (`.github/workflows/`) rather than a documented manual matrix. The pipeline must build
  real, portable wheels via **cibuildwheel** across manylinux/macOS × Python 3.9–3.13 (not
  just a source-build smoke test). Also run `R CMD check --as-cran` in CI (still required —
  see D-06).

- **D-03:** `DESCRIPTION` and `python/pyproject.toml` version drift is caught by a
  **test-suite assertion** (an R or Python test reads both files and fails on mismatch) —
  not a separate CI-only step, not a pre-commit hook. Runs on every DoD gate invocation
  (local and CI both), no new tooling.

- **D-04:** Write `cran-comments.md` even though CRAN submission itself is deferred (see
  D-06) — r-universe and any eventual CRAN submission both need it. Must explain, beyond
  new-submission boilerplate: the deliberate **no `-O` flag** in `PKG_CXXFLAGS`; the
  **LAPACK hard dependency** (`find_package(... REQUIRED)` on the Python side; note the
  R-side `SystemRequirements`/`configure` behavior if relevant).

- **D-05:** Do **NOT** submit to CRAN this phase. Publish via **GitHub + r-universe** as an
  intermediate distribution channel instead. Vignettes and other documentation CRAN/general
  users would expect are not yet written.

- **D-06:** `R CMD check --as-cran` (0 errors/0 warnings, roadmap SC1) still applies this
  phase, unchanged. r-universe builds and checks packages the same way CRAN does. Only the
  CRAN web-form submission step itself is skipped.

### Claude's Discretion

- Exact cibuildwheel config (which manylinux image, macOS runner versions, caching
  strategy) is an implementation detail for the researcher/planner to work out against
  current cibuildwheel docs — not discussed here.
- Whether the version-sync test lives in the R suite or Python suite is left to whichever
  the implementer finds cleaner to wire.

### Deferred Ideas (OUT OF SCOPE)

- Vignettes and user-facing documentation — needed before any real public release; own
  future phase.
- Actual CRAN web-form submission — deferred until docs/vignettes exist; `R CMD check
  --as-cran` still runs and must pass this phase, but nothing gets submitted.
</user_constraints>

## Phase Requirements

<phase_requirements>
| ID | Description | Research Support |
|----|-------------|------------------|
| US-010 | leafblower distributable on CRAN and PyPI (`install.packages()`/`pip install`) | Tarball-hygiene findings below (§ Runtime State Inventory), Standard Stack (cibuildwheel), Package Legitimacy Audit |
| US-008 (residual) | `pip install leafblower` installs a self-contained wheel, Python 3.9–3.13 test matrix | Don't Hand-Roll (cibuildwheel vs. manual manylinux), Common Pitfalls (LAPACK vendoring per-platform) |
| KPI-05 | `R CMD check --as-cran`: 0 errors, 0 warnings | Common Pitfalls (CRAN check hierarchy), Code Examples (`cran-comments.md` skeleton) |
| KPI-06 | Python wheel installs on Linux/macOS, Python 3.9–3.13, via a CI matrix | Architecture Patterns (GitHub Actions structure), State of the Art (cibuildwheel version range) |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| `cibuildwheel` | 3.x (latest stable at execution time — pin exact SHA/tag, see Pitfall below) [CITED: cibuildwheel.pypa.io] | Builds manylinux/macOS wheels across a Python version matrix, auto-repairs with auditwheel/delocate | The PyPA-maintained tool purpose-built for exactly this; used by NumPy/SciPy/pandas and the entire scientific-Python C-extension ecosystem [CITED: cibuildwheel.pypa.io] |
| `twine` | latest [CITED: pypi.org/project/twine] | Validates (`twine check`) and uploads wheel/sdist metadata | Standard PyPA upload tool; `twine check` catches metadata/README-render failures before any upload attempt [CITED: pydevtools.com/handbook/reference/twine] |
| `build` (PyPA `build`) | latest [CITED: pypi.org/project/build] | Builds sdist + wheel via PEP 517 from `pyproject.toml`, no setup.py invocation | Standard front-end for scikit-build-core-based projects; already implied by the existing `pyproject.toml` build-system table |
| `r-lib/actions/setup-r` + `rcmdcheck` | v2 (tag) [CITED: github.com/r-lib/actions] | Runs `R CMD check --as-cran` inside GitHub Actions | The tidyverse/r-lib-maintained standard for R CI; superior to hand-rolling an R install + check script |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| `auditwheel` | bundled inside cibuildwheel's manylinux repair step | Vendors shared libs (incl. LAPACK/OpenBLAS) into the Linux wheel, relabels manylinux tag | Automatic — cibuildwheel invokes it; no separate install/config unless overriding `CIBW_REPAIR_WHEEL_COMMAND_LINUX` |
| `delocate` | bundled inside cibuildwheel's macOS repair step | Same vendoring for macOS `.dylib` deps (LAPACK via Accelerate) | Automatic on macOS runners |
| `pypa/cibuildwheel` GitHub Action | `@v3` (major tag) [CITED: github.com/pypa/cibuildwheel] | Runs cibuildwheel directly as a CI step (`uses: pypa/cibuildwheel@v3.x`) instead of a manual `pip install cibuildwheel && cibuildwheel` invocation | Preferred over pip-installing cibuildwheel in a run step — matches the officially documented pattern |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| cibuildwheel | Hand-rolled manylinux Docker + auditwheel invocation per Python version | More control, dramatically more CI YAML to maintain and debug; cibuildwheel already encodes the manylinux-image-per-Python-version knowledge (D-02 explicitly rejected this — "documented manual matrix" was the non-chosen option) |
| cibuildwheel | `scikit-build-core` + plain `pip wheel .` in a matrix job, no repair step | Produces a wheel that is NOT self-contained on Linux — the LAPACK `.so` stays an external system dependency, silently violating "self-contained" (US-010/US-008's literal ask per CONTEXT.md's Wheel-scope decision) |
| `build`+`twine` two-step | `uv build` + `uv publish` | `uv` bundles both steps and is faster; but this repo's Python tooling directive (CLAUDE.md: "venv is uv-managed — NO pip") already mandates `uv` for the *dev* environment — `uv build`/`uv publish` is consistent with that and is a defensible substitution the planner may make; documented here as a legitimate lighter alternative, not a rejection of `twine check` (twine check itself has no uv equivalent and should still run as a distinct verification step regardless of which tool produces/uploads the artifact) |
| Version-sync test-suite assertion (D-03, locked) | `scikit-build-core` dynamic metadata (read `DESCRIPTION`'s `Version:` at build time via a custom Python metadata provider, drop the hardcoded `version = "0.1.0"` from `pyproject.toml` entirely) | Structurally eliminates drift instead of just detecting it — a genuinely stronger mechanism. Not chosen: D-03 already locked "test-suite assertion, no new tooling" specifically to avoid a metadata-provider script; documented here per CLAUDE.md's alternatives-considered mandate, not as an unresolved question |
| Version-sync test-suite assertion (D-03, locked) | `setuptools_scm`/`hatch-vcs` git-tag-derived versioning | Standard for many PyPI packages, but requires the git tag to be the source of truth instead of `DESCRIPTION` — inverts this project's existing convention where `DESCRIPTION` already IS the canonical version (R packages have no git-tag-version convention); would also require R-side tooling changes DESCRIPTION doesn't currently have. Rejected for the same reason as above — D-03 locked the simpler test assertion |

**Installation:**
```bash
# Python (build/release tooling — NOT a runtime dependency, dev-only)
uv pip install cibuildwheel twine build   # per CLAUDE.md's uv-only rule for python/.venv

# R (CI only — not a package dependency)
# r-lib/actions/setup-r + r-lib/actions/setup-r-dependencies handle this in CI;
# locally: Rscript -e 'install.packages(c("rcmdcheck"))' if running --as-cran manually
```

**Version verification:** `pip`/`npm view` were unavailable in this research session's
sandbox (no `pip`/`pip3` binary on PATH, no `npm`). Package identities and current-ness were
cross-checked via `mcp` package-legitimacy gate (below) and WebSearch against
`cibuildwheel.pypa.io`, `github.com/pypa/cibuildwheel`, `github.com/pypa/twine`,
`github.com/r-lib/actions` — all three PyPA tools confirmed to exist with recent releases;
exact pinned patch version must be verified at execution time with
`uv pip index versions cibuildwheel` (or equivalent) since this session had no working pip.

## Package Legitimacy Audit

| Package | Registry | Published (latest) | Downloads | Source Repo | Verdict | Disposition |
|---------|----------|---------------------|-----------|--------------|---------|-------------|
| `cibuildwheel` | PyPI | 2026-08-05 [VERIFIED: package-legitimacy gate] | unknown (lookup returned null, not zero) | `github.com/pypa/cibuildwheel` [VERIFIED: package-legitimacy gate] | SUS | Flagged by the automated gate on "too-new" + "unknown-downloads" signals — both are known false-positive triggers for a high-release-cadence, officially PyPA-maintained tool (org confirmed `pypa`). **Keep.** Planner must add a `checkpoint:human-verify` task before pinning the exact CI version, per protocol, despite this being a well-known first-party tool. |
| `twine` | PyPI | 2026-07-27 [VERIFIED: package-legitimacy gate] | unknown (lookup returned null) | `github.com/pypa/twine` [VERIFIED: package-legitimacy gate] | SUS | Same false-positive pattern (recent release + null download count). **Keep.** Planner must add `checkpoint:human-verify` before pinning version. |
| `build` | PyPI | 2026-04-30 [VERIFIED: package-legitimacy gate] | unknown | none returned by gate (actual repo is `github.com/pypa/build` per training knowledge — `[ASSUMED]`, gate's repo lookup failed to resolve it) | SUS | **Keep** — same PyPA-org pattern as the two above by cross-reference, but the repo-URL claim itself is `[ASSUMED]` (gate returned `null`, not a confirming URL). Planner must add `checkpoint:human-verify`. |

**Packages removed due to [SLOP] verdict:** none.
**Packages flagged as suspicious [SUS]:** `cibuildwheel`, `twine`, `build` — all three flagged
purely on "too-new"/"unknown-downloads" heuristics that misfire for high-cadence official
PyPA tooling; a human should still confirm the exact pinned version against
`pypi.org/project/<name>` at execution time before the planner locks a version string into
the CI workflow, since this research session had no working `pip`/`npm` to independently
query download counts.

## Architecture Patterns

### System Architecture Diagram

```
git push / PR
      |
      v
.github/workflows/  (new, D-02)
      |
      +--> job: r-check
      |      r-lib/actions/setup-r --> install deps --> R CMD build --> R CMD check --as-cran
      |      (fails CI on any error/warning; new-submission NOTE is the only tolerated NOTE)
      |
      +--> job: python-wheels (matrix: ubuntu-latest, macos-13, macos-14 x cp39..cp313)
      |      pypa/cibuildwheel@v3
      |        CIBW_BEFORE_ALL_LINUX: yum install lapack-devel openblas-devel   (manylinux only)
      |        CIBW_ENVIRONMENT_MACOS: BLA_VENDOR hint for Accelerate            (macOS only)
      |        --> cmake configure (find_package(LAPACK REQUIRED)) --> build --> auditwheel/delocate repair
      |      --> upload-artifact: wheelhouse/*.whl
      |
      +--> job: wheel-check (needs: python-wheels)
             download wheelhouse artifacts --> twine check wheelhouse/*
             --> install into clean venv (matrix 3.9-3.13) --> python -c "import leafblower; leafblower.harvest"

Local DoD gate (unchanged, both CI and dev machine):
  R CMD INSTALL --preclean . + testthat suite (incl. version-sync test, D-03)
  uv pip install -e . --reinstall-package leafblower + pytest (incl. version-sync test, D-03)
```

### Recommended Project Structure
```
.github/
└── workflows/
    ├── r-check.yml          # R CMD check --as-cran job
    └── python-wheels.yml    # cibuildwheel matrix + twine check + import smoke test
cran-comments.md              # new (D-04) — already .Rbuildignore'd, correctly so
tests/testthat/test-version-sync.R   # OR python/leafblower/test_version_sync.py (D-03, either location)
```

### Pattern 1: cibuildwheel platform-specific LAPACK provisioning
**What:** Different `CIBW_*` overrides per OS to satisfy the hard `find_package(LAPACK
REQUIRED)` at `python/CMakeLists.txt:93`.
**When to use:** Any CI matrix job building this wheel on Linux or macOS.
**Example:**
```toml
# pyproject.toml (or a separate cibuildwheel config block)
[tool.cibuildwheel]
build = "cp39-* cp310-* cp311-* cp312-* cp313-*"
skip = "*-musllinux*"          # not requested by CONTEXT.md (manylinux/macOS only)
test-command = "python -c \"import leafblower; leafblower.harvest\""

[tool.cibuildwheel.linux]
before-all = "yum install -y lapack-devel openblas-devel"   # manylinux2014 (yum-based)
# NOTE: verify the CIBW default manylinux image for each cp3x target before locking this —
# newer cibuildwheel defaults may use manylinux_2_28 (dnf-based) for some Python versions;
# confirm via `cibuildwheel --print-build-identifiers` at execution time [CITED: cibuildwheel.pypa.io/options]

[tool.cibuildwheel.macos]
environment = { CMAKE_ARGS = "-DBLA_VENDOR=Apple" }   # forces CMake's FindLAPACK to Accelerate
```
Source: [CITED: cibuildwheel.pypa.io/options], [CITED: cmake.org/cmake/help/latest/module/FindLAPACK.html]

### Pattern 2: two independent CI jobs, not one combined job
**What:** R check and Python wheel build run as separate GitHub Actions jobs (not sequential
steps in one job).
**When to use:** Always for a dual-language package like this one.
**Rationale:** The R and Python toolchains have no shared setup steps; a single failing R
check should not block wheel-matrix jobs from running (and vice versa), and the matrix
strategy (5 Python versions x 2-3 OS runners) only applies to the Python job.

### Anti-Patterns to Avoid
- **Source-build-only "wheel" job:** Running `pip wheel python/` inside a bare Ubuntu runner
  without cibuildwheel produces a wheel tagged for that exact glibc/Python ABI only, with
  LAPACK still an unvendored `.so` dependency — this fails "self-contained" and fails on any
  machine without a system LAPACK installed. CONTEXT.md's Wheel-scope decision already
  rejected this ("Source-build-only CI matrix" was the non-selected option).
- **Setting `DYLD_LIBRARY_PATH` directly in the cibuildwheel macOS job:** macOS SIP strips
  `DYLD_LIBRARY_PATH` before it reaches the `delocate-wheel` child process; a documented
  cibuildwheel workaround is to store the path under a different env var name and reference
  it explicitly inside `CIBW_REPAIR_WHEEL_COMMAND_MACOS` [CITED: github.com/pypa/cibuildwheel
  FAQ]. Only relevant if a custom repair command is needed — the default repair step does not
  hit this, but a custom one debugging LAPACK vendoring might.
- **Renormalizing `.Rbuildignore` before `git rm`-ing the tracked strays:** `.Rbuildignore`
  only controls what `R CMD build` puts in the tarball — it does NOT remove files from git
  history or the working tree. D-01 requires both; doing only the `.Rbuildignore` half still
  leaves `leafblower.Rcheck/` (with its stale `test-lbfgsb.R` etc.) bloating the git
  repository and available for accidental re-inclusion.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| manylinux-compatible Linux wheel with vendored LAPACK | Custom Docker image + manual `auditwheel repair` invocation | `pypa/cibuildwheel@v3` GitHub Action | cibuildwheel already encodes which manylinux image supports which CPython version, drives the repair step, and handles the artifact upload path — hundreds of scientific-Python projects rely on exactly this |
| macOS universal2/arm64+x86_64 wheel repair | Manual `lipo`/codesign scripting | cibuildwheel's built-in `delocate` repair step | Same tool the entire PyPA ecosystem uses; hand-rolling reintroduces the exact arm64-detection bugs the pypa/cibuildwheel issue tracker documents fixes for |
| Version-string parsing from `DESCRIPTION`/`pyproject.toml` | A TOML parser dependency (`tomli`) gated behind a Python-version check for the pre-3.11 matrix | Regex extraction (`Version:\s*(\S+)` / `version\s*=\s*"([^"]+)"`) — both files have a fixed, simple single-line format | Both target strings are simple `KEY: value`/`key = "value"` lines with no nested TOML structure to parse; a regex avoids adding a new dependency and avoids the `tomllib` (3.11+) vs. `tomli` (<3.11) split entirely, consistent with `ponytail`'s lazy-first ladder — rung 3 (stdlib) is unneeded when rung 6 (one line) already works |

**Key insight:** Both the wheel-repair problem and the version-drift-detection problem look
like they need a dedicated library/framework, but the actual mechanism in each case is
already either fully solved upstream (cibuildwheel) or trivially small (two regex reads) —
resist the urge to add TOML-parsing machinery for a one-line value.

## Runtime State Inventory

> Rename/refactor-adjacent: this phase moves tracked files out of git and into
> `.Rbuildignore`. Not a rename, but the same "what does grep miss" discipline applies —
> git tracking state is not visible to a plain filesystem `ls`, and `.Rbuildignore` alone
> does not touch it.

| Category | Items Found | Action Required |
|----------|-------------|------------------|
| **Git-tracked dev artifacts (root)** | `cell_table_92c4f45.cpp`, `cell_table_92c4f45.hpp`, `ieppa_92c4f45.cpp`, `patch_raking.py`, `patch_wolfe.py`, `test_output.log`, `leafblower_0.1.0.tar.gz`, `REVIEW_FINDINGS.md`, `code-review-findings.md` — all confirmed git-tracked via `git ls-files` [VERIFIED: `git ls-files` output, this session] | `git rm` each + add `.Rbuildignore` patterns (D-01) |
| **Git-tracked full R-check output tree** | `leafblower.Rcheck/` — confirmed **~90 tracked files** including compiled binaries (`src/*.o`, `src/leafblower.so`, `libs/leafblower.so`), compiled help (`R/leafblower.rdb`/`.rdx`, `help/*.rdb`/`.rdx`), and a **frozen historical copy** of `tests/testthat/` containing `test-lbfgsb.R` and `test-ieppa.R` — filenames that do not exist in the live `tests/testthat/` tree, confirming this is stale generated output, not source [VERIFIED: `git ls-files \| grep leafblower.Rcheck`, this session, full file list obtained] | `git rm -r leafblower.Rcheck/` — this is a strictly larger fix than D-01's single named `patch_wolfe.py` file inside it; the whole directory must go, not just that one file. `.Rbuildignore` already has `^leafblower\.Rcheck$` (line present) so the tarball itself is unaffected — this is a pure git-hygiene fix, not a CRAN-check fix |
| **Additional untracked-by-`.Rbuildignore` root files (not in D-01's list)** | `baseline_bench.R`, `security-review-2026-04-27-ieppa.json`, `GEMINI.md`, `.mcp.json`, `.coverage-thresholds.json` — all confirmed git-tracked at repo root with no matching `.Rbuildignore` pattern; R does not exclude dotfiles by default [VERIFIED: `git ls-files` root listing + `.Rbuildignore` content, this session; CITED: r-hub blog / R Packages (2e) confirming no default dotfile exclusion] | Decide per-file: `baseline_bench.R` and `security-review-2026-04-27-ieppa.json` look like the same class of dev artifact as D-01's list (add to `.Rbuildignore`, consider `git rm`); `GEMINI.md`/`.mcp.json`/`.coverage-thresholds.json` are tooling config that may be legitimate to keep tracked in git but should still be `.Rbuildignore`d (CRAN doesn't need them any more than it needs `CLAUDE.md`, which is already `.Rbuildignore`d) — planner should confirm disposition with the user if not obviously in the same bucket as D-01 |
| **Non-stray root file initially suspected** | `cleanup` — confirmed to be a standard R package special file (POSIX shell script; R's build/check machinery runs it to remove transient files post-build) [VERIFIED: `file cleanup` + content read, this session — `#!/bin/sh` / `rm -f src/Makevars`] | **No action** — do not `.Rbuildignore` or remove; it is expected packaging infrastructure, same tier as `configure` |
| **Fixture-backed tests and their guards** | 5 test files reference `tests/testthat/fixtures` (`.Rbuildignore`d): `test-calibration-solvers.R`, `test-sraa-global.R`, `test-eta-schedule.R`, `test-bench-gate.R`, `test-convergence-trajectory.R` — all 5 confirmed to already have `skip_if(!file.exists(...))` or `skip_if_not_installed(...)` guards on every fixture reference [VERIFIED: grep + read of all 5 files, this session] | **No action needed** — the ROADMAP's stated risk ("fixture-backed tests must keep their `skip_if` guards or the CRAN run errors") is already satisfied by the current code; nothing to fix here, just don't regress it |
| **Version strings** | `DESCRIPTION:3` reads `Version: 0.1.0`; `python/pyproject.toml:9` reads `version = "0.1.0"` — currently in sync [VERIFIED: `DESCRIPTION` and `pyproject.toml`, read this session] | Add the D-03 test now, while they're in sync, so any future edit that touches only one file fails immediately |
| **No CI / no `.github` directory** | Confirmed absent — `ls .github` exits non-zero [VERIFIED: `ls .github` / `ls .github/workflows`, this session] | Create `.github/workflows/` fresh (D-02) — no existing workflow to extend or conflict with |
| **No `Makevars.win`, no `configure.win`** | Confirmed absent from repo root and `src/` [VERIFIED: `ls` this session] | Not this phase's concern (Windows wheel is out of D-02's scope) but worth flagging: a future Windows R build (CRAN's win-builder) will fall back to `Makevars.in`'s `@...@` placeholders unsubstituted unless CRAN's own Windows toolchain handles the substitution — this is CRAN's win-builder's problem for a plain `configure`-based package (standard R Windows build pipeline: `configure.win` is only needed if the Unix `configure` script itself can't run on Windows, and CRAN's win-builder provides its own R toolchain that handles `Makevars.in` `@VAR@` substitution via `tools:::.win_config` conventions for packages without `configure.win`) — no action needed this phase, since D-05 has already deferred actual CRAN submission |

## Common Pitfalls

### Pitfall 1: Treating a CRAN "NOTE" as automatically acceptable
**What goes wrong:** Teams assume any NOTE (vs. ERROR/WARNING) is fine to ship.
**Why it happens:** The three-tier severity hierarchy looks like NOTEs are optional feedback.
**How to avoid:** CRAN's own guidance treats significant NOTEs the same as WARNINGs for
practical purposes — eliminate every NOTE that can be eliminated; only the "New submission"
NOTE is universally tolerated on a first release, and even that must be explicitly
acknowledged in `cran-comments.md` [CITED: r-pkgs.org/release.html, CRAN submission
checklist]. D-06 already locks "0 errors, 0 warnings ... at most the new-submission note" —
this pitfall is about not silently accumulating a second, third NOTE alongside it.
**Warning signs:** `R CMD check --as-cran` output listing more than one NOTE line.

### Pitfall 2: LAPACK vendoring succeeding on Linux, silently failing "self-contained" on macOS
**What goes wrong:** `find_package(LAPACK REQUIRED)` succeeds on a macOS CI runner (finds
Accelerate) but delocate has nothing to vendor because Accelerate is a system framework, not
a relocatable `.dylib` — this is actually fine (Accelerate ships with every macOS install),
but if CMake instead finds a Homebrew-installed OpenBLAS on the runner (common if `brew
install openblas` ran for some other reason), delocate WILL try to vendor that non-system
`.dylib`, and the result depends entirely on which LAPACK provider CMake happened to pick.
**Why it happens:** CMake's `FindLAPACK` searches multiple providers and picks the first
found; a CI runner's exact installed-package state affects the choice non-deterministically
run-to-run.
**How to avoid:** Explicitly pin the vendor with `-DBLA_VENDOR=Apple` (forces Accelerate,
skips the search) rather than relying on `find_package(LAPACK REQUIRED)`'s default search
order to consistently find the same provider [CITED: cmake.org FindLAPACK docs].
**Warning signs:** A macOS wheel that imports fine on the CI runner but fails
`import leafblower` on a clean machine with `Library not loaded: libopenblas.dylib` —
indicates a non-system LAPACK got linked and delocate either failed to vendor it or vendored
an incompatible arch slice.

### Pitfall 3: manylinux image mismatch across the Python version matrix
**What goes wrong:** cibuildwheel's default manylinux image selection can differ per CPython
version (e.g., older CPython targets defaulting to `manylinux2014`, newer ones to
`manylinux_2_28`), which changes which package manager (`yum` vs. `dnf`) is available inside
the build container — a single `before-all = "yum install ..."` line silently fails (or is
silently skipped) on a `dnf`-based image for a subset of the matrix.
**Why it happens:** cibuildwheel abstracts the image per-Python-version but the
`CIBW_BEFORE_ALL_LINUX` command is not abstracted — it runs verbatim inside whatever
container that Python version selected.
**How to avoid:** Either pin every target explicitly to the same manylinux image via
`manylinux-x86_64-image` overrides in `[tool.cibuildwheel.linux]`, or write the before-all
command to detect the package manager (`command -v yum && yum install ... || dnf install
...`) [CITED: cibuildwheel.pypa.io/options — overrides example].
**Warning signs:** Some Python versions in the matrix succeed, others fail at the
`find_package(LAPACK REQUIRED)` CMake step specifically (not at compile or link) — the
telltale sign the dev package never got installed on that particular image.

### Pitfall 4: `-O` flag creeping back into `PKG_CXXFLAGS` via a future edit
**What goes wrong:** A future contributor adds `-O2`/`-O3` to `src/Makevars.in`'s
`PKG_CXXFLAGS` line thinking it's a performance win, silently reintroducing the exact defect
`leafblower-qzto` already closed and D-04 requires `cran-comments.md` to proactively explain.
**Why it happens:** The reasoning ("R supplies `-O` via user/site `Makevars`, CRAN's
`tools:::.check_make_vars` rejects package-supplied `-O` flags") is documented in a comment
in `src/Makevars.in` but not enforced by any check.
**How to avoid:** `R CMD check --as-cran` itself catches this (the `tools:::.check_make_vars`
sub-check runs as part of `--as-cran`) — no separate guard needed, but `cran-comments.md`
should state this is intentional so a human reviewer doesn't "fix" it back in before the
check ever flags it.
**Warning signs:** A NOTE from `--as-cran` about non-portable flags in `PKG_CXXFLAGS`.

### Pitfall 5: pinning `pypa/cibuildwheel@v3` (major tag) vs. a specific SHA
**What goes wrong:** Using a moving major-version tag (`@v3`) in a CI workflow means the
exact cibuildwheel version — and therefore exact manylinux image defaults — can change
between CI runs without any change to this repo, causing "worked yesterday, fails today."
**Why it happens:** GitHub Actions convention favors major-tag pinning for convenience;
cibuildwheel's own examples use `@v3.4`-style tags.
**How to avoid:** For a release-critical pipeline, pin either an exact minor/patch tag
(`@v3.x.y`) or a commit SHA, and bump deliberately. cibuildwheel's own GitHub Action has
moved toward full-SHA pinning for its *own* dependencies as a security practice
[CITED: cibuildwheel changelog] — apply the same discipline to how this repo references it.
**Warning signs:** A wheel-matrix job failure with no corresponding change in this repo's
diff.

## Code Examples

### Version-sync test (R side, D-03 — one of two acceptable locations)
```r
# tests/testthat/test-version-sync.R
test_that("DESCRIPTION and pyproject.toml versions match", {
  desc_version <- read.dcf("../../DESCRIPTION")[1, "Version"]
  pyproject_lines <- readLines("../../python/pyproject.toml")
  py_version_line <- grep('^version\\s*=', pyproject_lines, value = TRUE)[1]
  py_version <- sub('^version\\s*=\\s*"([^"]+)".*$', '\\1', py_version_line)
  expect_equal(desc_version, py_version)
})
```
Path note: `test_path()`/relative paths inside `tests/testthat/` resolve against the package
root differently under `devtools::test()` vs. `testthat::test_dir()` — verify the correct
relative path at implementation time (existing tests in this repo already navigate this;
follow their established pattern rather than the illustrative `../../` above). [ASSUMED —
exact relative-path convention not verified against a concrete existing example this
session; planner/implementer should confirm against an existing testthat file that reads a
repo-root file, if one exists, or use `rprojroot`/`testthat::test_path()`'s documented
root-finding instead of a literal `../../`.]

### Version-sync test (Python side, D-03 — alternative location)
```python
# python/leafblower/test_version_sync.py
import re
from pathlib import Path

def test_description_pyproject_version_match():
    root = Path(__file__).resolve().parents[2]  # repo root
    desc_text = (root / "DESCRIPTION").read_text()
    desc_version = re.search(r"^Version:\s*(\S+)", desc_text, re.M).group(1)
    pyproject_text = (root / "python" / "pyproject.toml").read_text()
    py_version = re.search(r'^version\s*=\s*"([^"]+)"', pyproject_text, re.M).group(1)
    assert desc_version == py_version
```
[ASSUMED — exact repo-root path depth from `python/leafblower/` not independently verified
this session; confirm `parents[N]` depth against the actual installed/test-discovered file
location before use.]

### cran-comments.md skeleton (D-04)
```markdown
## Submission type

This is an intermediate GitHub + r-universe release, not a CRAN submission (see project
notes). `R CMD check --as-cran` is run and passes cleanly regardless, as the honest quality
bar for this stage of the project.

## Test environments

* local: R 4.x, <platform>
* r-universe / GitHub Actions: <platforms in CI matrix>

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission NOTE — expected for a first release.

## Notes on build configuration

* `PKG_CXXFLAGS` intentionally carries no `-O` flag. R supplies the optimization level via
  the user/site `Makevars` `$(CXXFLAGS)`; a package-supplied `-O` flag is rejected by
  `tools:::.check_make_vars` under `--as-cran`. This is deliberate, not an oversight.
* The package links `$(LAPACK_LIBS) $(BLAS_LIBS)` via R's own build-supplied macros
  (`src/Makevars.in`) — this is the standard R mechanism and requires no additional
  `SystemRequirements` declaration, since R itself is built against LAPACK. This is distinct
  from the Python bindings (`python/`, not part of the R package build), which have a hard
  `find_package(LAPACK REQUIRED)` CMake dependency with no fallback — irrelevant to the R
  CRAN check but documented here for reviewer context.
```
[ASSUMED structure — synthesized from CITED cran-comments.md conventions (r-pkgs.org,
observed example diffs), not a verbatim copy of any single source; content requirements
(no-`-O`, LAPACK) are the locked D-04 decision, not assumed.]

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Manual `auditwheel`/`delocate` invocation per platform in hand-written CI scripts | `cibuildwheel` as a single cross-platform driver | Established practice for several years, still current in 2026 [CITED: cibuildwheel.pypa.io] | Drastically less CI YAML; correctness of manylinux-image-to-Python-version mapping owned upstream |
| `python setup.py bdist_wheel` | PEP 517 `build` front-end (`python -m build`) reading `pyproject.toml` | setuptools-era pattern fully superseded; this repo never used the old pattern (already on scikit-build-core) | Not a migration for this repo — noted only because a wheel-build task description should not accidentally reintroduce `setup.py` |
| Username/password PyPI upload | Trusted publishing (OIDC) or `__token__`/API-token auth | PyPI removed password auth April 2024 [CITED: WebSearch summary of pypi.org changelog] | Any upload step (even if manual/deferred past this phase) must use a token or trusted publishing, never a password |

**Deprecated/outdated:**
- Plain `pip wheel` without a repair step for any package with a compiled extension linking
  a system library: produces a wheel that is not portable off the build machine. Not
  "deprecated" in the formal sense, just never sufficient for this project's LAPACK
  dependency.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | Exact relative path (`../../`) for the R version-sync test's `DESCRIPTION` read | Code Examples | Test fails to find the file under `testthat::test_dir()` vs `devtools::test()` — caught immediately by TDD RED-phase-first execution per CLAUDE.md, low risk |
| A2 | Exact `parents[N]` depth for the Python version-sync test's repo-root resolution | Code Examples | Same as A1 — caught immediately on first test run |
| A3 | cibuildwheel's exact current stable version number/patch (WebSearch could not resolve a precise number; classify-confidence gate marked all three PyPA tools SUS on heuristic signals, not on actual illegitimacy) | Standard Stack, Package Legitimacy Audit | Low — these are unambiguously the correct, well-known official tools by name and org; only the *exact pinned version string* needs execution-time confirmation, not the choice of tool |
| A4 | cran-comments.md's exact conventional section headers/structure | Code Examples | Low — CRAN accepts free-form markdown; the content requirements (D-04) are locked, not assumed, only the formatting scaffold is synthesized |
| A5 | CRAN win-builder's handling of `Makevars.in` `@VAR@` substitution without `configure.win` present | Runtime State Inventory | Not this phase's concern (D-05 defers actual CRAN submission) — flagged only as a forward-looking note, zero risk to this phase's success criteria |

**If this table is empty:** N/A — table is populated; five low-risk items, none blocking.

## Open Questions

1. **Disposition of the 5 additional untracked-by-`.Rbuildignore` root files found beyond D-01's list**
   - What we know: `baseline_bench.R`, `security-review-2026-04-27-ieppa.json`, `GEMINI.md`,
     `.mcp.json`, `.coverage-thresholds.json` are all git-tracked, none `.Rbuildignore`d,
     confirmed this session.
   - What's unclear: CONTEXT.md's D-01 explicitly named only 8 files + the discovered 9th
     (`patch_wolfe.py` inside `leafblower.Rcheck/`). It did not discuss these 5, and did not
     discuss the fact that ALL of `leafblower.Rcheck/` is tracked (not just that one file
     inside it).
   - Recommendation: Planner should treat "no development artifact" (phase Success
     Criterion 2) as the binding requirement and extend D-01's mechanical fix (`git rm` +
     `.Rbuildignore` pattern) to cover all 5 files plus the full `leafblower.Rcheck/`
     directory, since SC2 is broader than D-01's enumerated list and this is squarely within
     the phase's existing scope (not a scope expansion, just a more complete application of
     the same locked mechanism) — flag to user only if any of the 5 look intentionally
     tracked for a reason not visible from filename alone (e.g., `.mcp.json` might be
     legitimately needed by contributors and only needs `.Rbuildignore`, not `git rm`).

2. **Exact manylinux image per Python 3.9-3.13 target**
   - What we know: cibuildwheel selects a default manylinux image per Python version, and
     recent versions may split older/newer CPython targets across `manylinux2014` (yum) and
     `manylinux_2_28`/`manylinux_2_34` (dnf).
   - What's unclear: The exact current (2026-08-15) default mapping was not pinned down by
     WebSearch to a specific per-version table.
   - Recommendation: Planner/implementer should run `cibuildwheel --print-build-identifiers`
     (or check `cibuildwheel.pypa.io/options#linux-image` at execution time) before writing
     the final `before-all` command, and prefer the package-manager-detecting shell snippet
     from Pitfall 3 over a single hardcoded `yum install` line if the matrix spans both image
     families.

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|------------|-----------|---------|----------|
| `R` / `Rscript` | R CMD check, R tests | Yes | 4.6.1 [VERIFIED: `Rscript --version`, this session] | — |
| `pip`/`pip3` | Version verification of Python packaging tools | No — not on PATH in this sandbox [VERIFIED: `command -v pip3` failed, this session] | — | `uv pip index versions <pkg>` (per this project's uv-only convention) at implementation time |
| `.github/workflows/` (existing CI) | D-02's CI pipeline | No — directory does not exist | — | None needed — this phase creates it fresh |
| `Makevars.win` / `configure.win` | Windows R CRAN build | No — absent | — | Not required this phase (Windows out of D-02's scope; CRAN win-builder has its own fallback per Open Question / Runtime State Inventory note) |
| Docker (for local manylinux testing) | Optional local cibuildwheel dry-run before pushing to CI | Not probed this session | — | Not required — cibuildwheel can run entirely inside GitHub Actions without a local Docker requirement for this phase's scope |

**Missing dependencies with no fallback:** none — every gap above has a documented fallback
or is genuinely out of this phase's locked scope.

**Missing dependencies with fallback:** `pip`/`pip3` (use `uv` per CLAUDE.md), local
CI/Docker dry-run tooling (not required, CI itself is the verification surface).

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | testthat 3 (`Config/testthat/edition: 3` in `DESCRIPTION`, confirmed [VERIFIED: `DESCRIPTION`, this session]) / pytest (Python) |
| Config file | `DESCRIPTION` (R, testthat config lives there) / `python/pyproject.toml` `[project.optional-dependencies] test = ["pytest"]` [VERIFIED: `pyproject.toml`, this session] |
| Quick run command | `Rscript -e "testthat::test_dir('tests/testthat', filter='version-sync')"` / `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest python/leafblower/test_version_sync.py` |
| Full suite command | `R CMD INSTALL --preclean .` then `Rscript -e "devtools::test()"` / `uv pip install -e . --reinstall-package leafblower` then the pinned-thread-count pytest invocation (both per `CLAUDE.md` Build & Test) |

### Phase Requirements → Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| US-010 (SC1) | `R CMD check --as-cran` 0 errors/0 warnings | CI/manual | `R CMD build . && R CMD check --as-cran leafblower_*.tar.gz` | ❌ Wave 0 (new CI job) |
| US-010 (SC2) | Tarball has no dev artifact; `git ls-files` shows no tracked build output | manual/scripted | `git ls-files \| grep -E '<stray patterns>'` (should be empty) | ❌ Wave 0 (no existing check script for this) |
| US-008 (SC3) | `pip wheel python/` passes `twine check`; clean-env import works | CI | `cibuildwheel` build + `twine check wheelhouse/*` + `python -c "import leafblower; leafblower.harvest"` in a fresh venv | ❌ Wave 0 (new CI job) |
| US-008 (SC4) | Wheel imports/calibrates on Python 3.9-3.13 | CI matrix | Same wheel-check job, matrixed across `cp39`-`cp313` | ❌ Wave 0 |
| US-010 (SC5) | `DESCRIPTION`/`pyproject.toml` version match | unit | `testthat::test_dir(filter='version-sync')` or `pytest python/leafblower/test_version_sync.py` | ❌ Wave 0 (new test file, D-03) |

### Sampling Rate
- **Per task commit:** Quick run of whichever new test file was just added (version-sync
  test) plus the existing DoD quick-check habits already established in this repo.
- **Per wave merge:** Full DoD gate (`R CMD INSTALL --preclean .` + full testthat + full
  pytest, single-thread BLAS) — unchanged from existing project convention.
- **Phase gate:** `R CMD check --as-cran` clean (0/0 + only new-submission NOTE) AND a full
  cibuildwheel matrix run green AND `twine check` clean AND the version-sync test passing —
  all four are the phase's actual success criteria, not proxies for them.

### Wave 0 Gaps
- [ ] `.github/workflows/r-check.yml` — new, runs `R CMD check --as-cran` (covers SC1)
- [ ] `.github/workflows/python-wheels.yml` — new, runs cibuildwheel matrix + twine check +
      import smoke test (covers SC3, SC4)
- [ ] `tests/testthat/test-version-sync.R` OR `python/leafblower/test_version_sync.py` — new,
      one file only per D-03 (covers SC5)
- [ ] A tarball-hygiene verification step (either a one-off manual `git ls-files` audit
      documented in the plan, or a small script) — covers SC2; no existing test infrastructure
      for this in the repo

*(No existing test infrastructure covers any of this phase's 5 success criteria — all five
are genuinely new verification surface, consistent with this being the "largest
genuinely-unfinished requirement" per REQUIREMENTS.md's US-010 note.)*

## Security Domain

### Applicable ASVS Categories

| ASVS Category | Applies | Standard Control |
|---------------|---------|-----------------|
| V2 Authentication | No | This phase does not add any auth surface |
| V3 Session Management | No | N/A |
| V4 Access Control | No | N/A |
| V5 Input Validation | No | No new user-facing input surface — packaging/build only |
| V6 Cryptography | No | N/A |

### Known Threat Patterns for this stack

| Pattern | STRIDE | Standard Mitigation |
|---------|--------|---------------------|
| Malicious/compromised third-party GitHub Action supply-chain injection | Tampering | Pin Actions (including `pypa/cibuildwheel`, `r-lib/actions/*`) to a specific tag or commit SHA rather than a floating major tag, per Pitfall 5 above; this is the standard mitigation for any new CI pipeline, not specific to this project |
| CI-published PyPI credentials leaking via workflow logs or a compromised Action | Elevation of Privilege / Information Disclosure | If/when an actual `twine upload`/`pypi publish` step is added (not required this phase per D-05's CRAN-deferral, but the PyPI wheel-publish equivalent may still be in scope for US-008/US-010's PyPI half — clarify with planner), use PyPI Trusted Publishing (OIDC, no stored token) over a stored `PYPI_API_TOKEN` secret [CITED: State of the Art table, PyPI auth changes] |
| Auto-vendoring an unintended/incompatible LAPACK provider into a published wheel (see Pitfall 2) | Tampering (of the build artifact, not malicious but unintentional) | Pin `BLA_VENDOR` / manylinux image explicitly rather than relying on `find_package`'s default search order |

## Sources

### Primary (HIGH confidence — files read directly this session)
- `/home/dd/Gemini/leafblower/.Rbuildignore` — full exclusion list
- `/home/dd/Gemini/leafblower/DESCRIPTION` — version, testthat edition, SystemRequirements
- `/home/dd/Gemini/leafblower/python/pyproject.toml` — version, build-system, cibuildwheel-relevant scikit-build config
- `/home/dd/Gemini/leafblower/python/CMakeLists.txt` — LAPACK REQUIRED (line 93), `-O3` (line 99), CORE_SOURCES (lines 60-78), AVX2 gating (lines 37-52, 110-112)
- `/home/dd/Gemini/leafblower/src/Makevars.in` — PKG_CXXFLAGS/PKG_LIBS, no-`-O` comment block
- `/home/dd/Gemini/leafblower/configure` — C++17/C++14 fallback (lines 2-8)
- `/home/dd/Gemini/leafblower/cleanup` — confirmed standard R packaging special file
- `git ls-files` (this session) — full tracked-file audit at repo root and inside `leafblower.Rcheck/`
- `.planning/phases/05-cran-pypi-release/05-CONTEXT.md`, `.planning/REQUIREMENTS.md`, `.planning/STATE.md`

### Secondary (MEDIUM confidence — WebSearch verified against official docs)
- cibuildwheel.pypa.io/options, github.com/pypa/cibuildwheel — CIBW_BEFORE_ALL, overrides, macOS SIP/DYLD_LIBRARY_PATH workaround
- r-pkgs.org/release.html, CRAN submission checklist (cran.r-project.org/web/packages/submission_checklist.html) — CRAN check severity hierarchy, cran-comments.md conventions
- github.com/r-lib/actions — check-r-package / setup-r workflow patterns, v2.12.0 (2026-04-30)
- cmake.org/cmake/help/latest/module/FindLAPACK.html — BLA_VENDOR mechanism

### Tertiary (LOW confidence — WebSearch only, not cross-checked against a second source)
- Exact current cibuildwheel/twine/build pinned version numbers (package-legitimacy gate
  flagged all three SUS on heuristic signals; treat as needing execution-time reconfirmation,
  see Assumptions Log A3)
- Exact manylinux-image-per-Python-version default mapping (Open Question 2)

## Metadata

**Confidence breakdown:**
- Standard stack: MEDIUM — cibuildwheel/twine/build tool choices are HIGH confidence (unambiguous, well-known official tooling); exact version pins are LOW confidence pending execution-time `uv pip index versions` check
- Architecture: HIGH — CI job split and LAPACK-provisioning pattern verified against current official cibuildwheel docs this session
- Pitfalls: HIGH for tarball-hygiene/`-O`-flag pitfalls (all verified against this repo's actual tracked files and source); MEDIUM for manylinux-image and macOS-Accelerate pitfalls (WebSearch-sourced, not reproduced in a live CI run this session)

**Research date:** 2026-08-15
**Valid until:** 2026-09-15 (30 days — cibuildwheel/manylinux-image defaults and CRAN policy specifics move fast enough to warrant re-verification for a fast-moving packaging domain; the repo-state findings in Runtime State Inventory are valid until the next commit touching any of the listed files)
