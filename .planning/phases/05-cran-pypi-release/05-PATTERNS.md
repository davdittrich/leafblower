# Phase 5: CRAN + PyPI Release - Pattern Map

**Mapped:** 2026-08-15
**Files analyzed:** 6 (new/modified)
**Analogs found:** 4 / 6

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `python/leafblower/test_version_sync.py` (or `tests/testthat/test-version-sync.R`, D-03 either location) | test | file-I/O (read-and-compare) | `python/leafblower/test_core_sources_sync.py` | exact |
| `.Rbuildignore` (edit) | config | transform (exclusion-list append) | `.Rbuildignore` (itself — extend existing pattern list) | exact |
| `cran-comments.md` | config/docs | transform (static markdown) | none in-repo — RESEARCH.md skeleton is the source | no analog |
| `.github/workflows/r-check.yml` | config (CI) | event-driven (push/PR trigger) | none in-repo (no `.github/` exists) | no analog |
| `.github/workflows/python-wheels.yml` | config (CI) | event-driven (push/PR trigger) | none in-repo (no `.github/` exists) | no analog |
| Tracked stray removal (`git rm` of 9+ files, D-01) | — (no new file) | — | `.Rbuildignore` existing entries (pattern for what NOT to ship) | role-match |

## Pattern Assignments

### `python/leafblower/test_version_sync.py` (test, file-I/O)

**Analog:** `python/leafblower/test_core_sources_sync.py` (full file, 36 lines, read verbatim above)

This is not just the closest analog — it is the SAME class of test (build/config-drift
guard, read two files, regex-extract, assert set/string equality) already established at
`python/leafblower/test_core_sources_sync.py:1-36`. Copy its structure exactly:

**Repo-root resolution pattern** (line 21):
```python
REPO_ROOT = Path(__file__).resolve().parents[2]
```
This resolves `python/leafblower/test_core_sources_sync.py` → repo root via `parents[2]`
(file → `leafblower/` → `python/` → repo root). Confirms RESEARCH.md's Assumption A2 —
`parents[2]` is the correct, already-proven depth from this exact directory. Use identically
in `test_version_sync.py`.

**Regex-extraction pattern** (lines 28-30, adapted):
```python
def test_description_matches_pyproject_version():
    desc_text = (REPO_ROOT / "DESCRIPTION").read_text()
    desc_version = re.search(r"^Version:\s*(\S+)", desc_text, re.M).group(1)
    pyproject_text = (REPO_ROOT / "python" / "pyproject.toml").read_text()
    py_version = re.search(r'^version\s*=\s*"([^"]+)"', pyproject_text, re.M).group(1)
    assert desc_version == py_version, (
        f"DESCRIPTION Version={desc_version!r} != pyproject.toml version={py_version!r}"
    )
```

**Docstring pattern** (lines 1-16 of analog): opens with a `SC<N> (leafblower-<id> / phase
<N>)` ticket header, explains WHY the drift is dangerous and WHERE it currently surfaces
without the test, states the verified-baseline date. Follow the same structure — one
paragraph, ticket ID first line, "Verified baseline (date): X == Y — zero drift today" as
the closing line (mirrors analog lines 13-15).

**No mocking/fixtures needed** — same as analog: plain `Path.read_text()`, no pytest
fixtures, no parametrize (single assertion, matches analog's minimalism).

If D-03 is instead placed in the R suite (`tests/testthat/test-version-sync.R`), the
**path-resolution analog** is `tests/testthat/test-bench-gate.R` (excerpt already read):
```r
pkg_root <- normalizePath(file.path(testthat::test_path(), "../.."))
```
— `testthat::test_path()` is the established, cwd-agnostic root-finder in this repo (used
because `test-bench-gate.R`'s own comment explains bare relative paths break under a
filtered `testthat::test_dir()` run). Use this exact idiom, NOT the illustrative `../../`
literal from RESEARCH.md's Code Examples section (RESEARCH.md itself flags this as
unverified/A1 — this repo's actual codebase already has the answer).

---

### `.Rbuildignore` (config, transform)

**Analog:** itself — append to the existing pattern list, matching its exact style
(anchored regex, one pattern per line, no comments needed for simple filename excludes):
```
^\.beads$
^\.claude$
^cran-comments\.md$
^\.github$
```
New entries for D-01 strays follow the same `^literal\.ext$` anchoring, e.g.:
```
^cell_table_92c4f45\.cpp$
^cell_table_92c4f45\.hpp$
^ieppa_92c4f45\.cpp$
^patch_raking\.py$
^patch_wolfe\.py$
^test_output\.log$
^leafblower_0\.1\.0\.tar\.gz$
^REVIEW_FINDINGS\.md$
^code-review-findings\.md$
^baseline_bench\.R$
^security-review-2026-04-27-ieppa\.json$
```
Note: `^\.github$` is ALREADY present (line present in current file) — no new entry needed
for the CI workflow directory itself.

---

### `cran-comments.md` (docs, static content)

**No in-repo analog** — this file type does not exist anywhere in the repo yet (confirmed
already `.Rbuildignore`d at line `^cran-comments\.md$`, i.e., expected-to-exist-but-absent).
Use RESEARCH.md's `Code Examples` skeleton verbatim as the structural template (already
cites D-04's exact required content: no-`-O` explanation, LAPACK hard-dependency
explanation). Cross-reference the exact source-of-truth comments already in the codebase
rather than paraphrasing:
- No-`-O` rationale: `src/Makevars.in` comment block (read above, lines "PKG_CXXFLAGS =
  ... -O3 is intentionally NOT set: R supplies the user/site -O level via $(CXXFLAGS)...").
- LAPACK: `python/CMakeLists.txt` lines ~93-95 (`find_package(LAPACK REQUIRED)`) plus the
  R-side `src/Makevars.in` line `PKG_LIBS = ... $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` (R's
  own build-supplied macros — no separate `SystemRequirements` needed on the R side).

---

### `.github/workflows/r-check.yml` and `.github/workflows/python-wheels.yml` (CI config, event-driven)

**No in-repo analog** — `.github/` does not exist (confirmed via `ls` — directory absent).
No prior CI pipeline to pattern-match against in this codebase. Use RESEARCH.md's
`Architecture Patterns` § "System Architecture Diagram" and § "Pattern 1/Pattern 2" as the
structural template — two independent jobs (r-check, python-wheels), not one combined job
(RESEARCH.md Pattern 2 rationale: no shared toolchain setup, independent failure isolation).

The one hard constraint pulled from THIS repo (not external docs): both workflows must
invoke the SAME Definition-of-Done gate already established in `CLAUDE.md` § Build & Test —
`R CMD INSTALL --preclean .` (not `devtools::install`) and the pinned single-thread BLAS env
vars (`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`) for the Python job's test
step, to avoid parity/determinism drift between local and CI runs.

---

## Shared Patterns

### Repo-root path resolution (Python side)
**Source:** `python/leafblower/test_core_sources_sync.py:21` — `Path(__file__).resolve().parents[2]`
**Apply to:** `test_version_sync.py` and any other new Python test file that needs to read
repo-root files (`DESCRIPTION`, `pyproject.toml`).

### Repo-root path resolution (R side)
**Source:** `tests/testthat/test-bench-gate.R` — `normalizePath(file.path(testthat::test_path(), "../.."))`
**Apply to:** `test-version-sync.R` if the R-side location is chosen instead (D-03 discretion).
**Rule:** never use a bare relative path literal (`"../../DESCRIPTION"`) — breaks under a
filtered `testthat::test_dir(filter=...)` run per the analog's own inline comment.

### `.Rbuildignore` anchoring style
**Source:** `.Rbuildignore` (existing file, all 23 current lines)
**Apply to:** every new D-01 stray-exclusion entry — `^literal\.ext$` anchored regex, one
per line, no wildcard globs used anywhere in the existing file.

### No-`-O` / LAPACK build-flag rationale (single source of truth)
**Source:** `src/Makevars.in` (no-`-O` comment block) + `python/CMakeLists.txt` (`-O3`
comment block referencing "phase-02 SC2, leafblower-qzto") + `python/CMakeLists.txt:93`
(`find_package(LAPACK REQUIRED)`)
**Apply to:** `cran-comments.md` — do not re-derive the rationale, quote/paraphrase these
exact existing comments so the CI-facing doc and the source-code doc never drift apart.

## No Analog Found

| File | Role | Data Flow | Reason |
|------|------|-----------|--------|
| `cran-comments.md` | docs | static | No prior cran-comments file exists (project has never done a CRAN-style release); use RESEARCH.md's skeleton + this repo's own Makevars.in/CMakeLists.txt comments as content source |
| `.github/workflows/r-check.yml` | config (CI) | event-driven | No `.github/` directory exists at all — greenfield CI; use RESEARCH.md Architecture Patterns as template |
| `.github/workflows/python-wheels.yml` | config (CI) | event-driven | Same — greenfield CI |

## Metadata

**Analog search scope:** repo root, `tests/testthat/`, `python/leafblower/`, `python/CMakeLists.txt`, `src/Makevars.in`, `.Rbuildignore`, `.github/` (confirmed absent)
**Files scanned:** `.Rbuildignore`, `DESCRIPTION`, `python/pyproject.toml` (head), `python/CMakeLists.txt` (LAPACK/-O3 section), `src/Makevars.in` (full), `tests/testthat/test-bench-gate.R` (partial), `python/leafblower/test_core_sources_sync.py` (full), `python/leafblower/test_solver_parity.py` (docstring)
**Pattern extraction date:** 2026-08-15
