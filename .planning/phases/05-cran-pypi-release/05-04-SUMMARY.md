---
phase: 05-cran-pypi-release
plan: 04
subsystem: infra
tags: [cibuildwheel, pypi, wheel, manylinux, macos, lapack, scikit-build-core, sdist]

# Dependency graph
requires:
  - phase: 05-cran-pypi-release
    provides: "05-01's hygiene-clean tree and cran-comments.md baseline that this plan's Python-track wheel matrix parallels"
provides:
  - "[tool.cibuildwheel] config in python/pyproject.toml (cp39-cp313, manylinux/macOS, LAPACK-provisioned)"
  - ".github/workflows/python-wheels.yml (python-wheels build matrix + wheel-check twine/import gate)"
  - "a working sdist for the Python package -- python/pyproject.toml's [tool.scikit-build.sdist] force-include/exclude config and python/CMakeLists.txt's LBW_SRC_DIR, fixing a pre-existing 'uv build'/'python -m build' failure that would have also broken cibuildwheel's own sdist-based build path and any future real PyPI publish"
affects: [05-05]

# Actuals (#2632)
actuals:
  tokens: 1515
  tasks: 2
  commits: 2

# Tech tracking
tech-stack:
  added: [cibuildwheel 4.2.0, twine 7.0.0, build 1.5.0 (per-checkpoint approved), auditwheel]
  patterns:
    - "scikit-build-core's sdist walk can never reach outside the project root (python/); a
      CORE_SOURCES tree living one level up (../src, shared with the R build) must be
      force-included into the sdist under an internal src/ path, with sdist.exclude
      dropping the gitignored build byproducts (*.o, *.so, Makevars) that live alongside
      the real sources in that directory"
    - "CMakeLists.txt resolves a build-layout-detecting LBW_SRC_DIR (../src for a dev
      checkout, src for an extracted sdist) via EXISTS(...c_api.cpp) rather than hardcoding
      one relative path -- the same CMakeLists.txt now serves both layouts unmodified"

key-files:
  created:
    - .github/workflows/python-wheels.yml
  modified:
    - python/pyproject.toml
    - python/CMakeLists.txt

key-decisions:
  - "Fixed the sdist packaging bug at the root cause (force-include + CMake layout
    detection) rather than only running 'uv build --wheel' (wheel-only, skips the sdist
    round-trip) for Task 2's local proof -- the plan's own acceptance criteria requires
    'uv build' (both sdist+wheel) to exit 0, and a broken sdist is a real defect
    independent of this plan's cibuildwheel scope (PyPI publish uploads sdist+wheels
    together; D-05 defers publish but doesn't defer sdist correctness)"
  - "auditwheel repair vendored Intel MKL (this sandbox's find_package(LAPACK) provider,
    found via 'Looking for cheev_ - found' against /opt/intel/oneapi/mkl) rather than
    OpenBLAS -- proves the vendoring MECHANISM (rpath patch + leafblower.libs/ vendoring
    into a manylinux_2_39_x86_64-tagged wheel) works for whichever LAPACK provider CMake's
    find_package resolves; CI's manylinux container will resolve openblas/lapack-devel
    per Task 1's before-all instead, exercising the same mechanism against a different
    provider"
  - "Did not fix twine check's two pre-existing WARNINGs (missing long_description/
    long_description_content_type) -- unrelated to this plan's cibuildwheel/sdist scope,
    pyproject.toml carries no README reference at all; documented as a known gap, not a
    blocking finding (twine check still reports PASSED, not FAILED)"

requirements-completed: [US-010, US-008, KPI-06]

coverage:
  - id: D1
    description: "[tool.cibuildwheel] config authored: cp39-cp313 build matrix, musllinux
      skipped, linux before-all installs lapack-devel/openblas-devel (yum/dnf fallback),
      macos environment pins CMAKE_ARGS to -DBLA_VENDOR=Apple"
    requirement: US-010
    verification:
      - kind: other
        ref: "grep -c '\\[tool.cibuildwheel\\]' python/pyproject.toml -> 1; grep -c 'cp39-\\* cp310-\\* cp311-\\* cp312-\\* cp313-\\*' python/pyproject.toml -> 1"
        status: pass
    human_judgment: false
  - id: D2
    description: ".github/workflows/python-wheels.yml authored: python-wheels job
      (manylinux/macos-13/macos-14 matrix, pinned pypa/cibuildwheel@v4.2.0) + wheel-check
      job (needs: python-wheels; twine check + fresh-venv import across py3.9-3.13)"
    requirement: US-010
    verification:
      - kind: other
        ref: "Rscript -e 'wf <- yaml::read_yaml(\".github/workflows/python-wheels.yml\"); stopifnot(is.list(wf$jobs$\"python-wheels\"), is.list(wf$jobs$\"wheel-check\"))' -> valid YAML, both jobs present"
        status: pass
    human_judgment: false
  - id: D3
    description: "Build/repair/import mechanics proven locally without Docker: uv build
      (sdist+wheel) exits 0, twine check PASSED, auditwheel repair vendors LAPACK into a
      self-contained manylinux-tagged wheel"
    requirement: KPI-06
    verification:
      - kind: other
        ref: "cd python && uv build && .venv/bin/python -m twine check dist/* && .venv/bin/python -m auditwheel repair dist/*.whl -w wheelhouse/ -> sdist+wheel built, twine PASSED with warnings (missing long_description, pre-existing/out-of-scope), repaired wheel manylinux_2_39_x86_64 with leafblower.libs/{libmkl_core,libmkl_intel_lp64,libmkl_sequential} vendored"
        status: pass
    human_judgment: false
  - id: D4
    description: "C++ extension builds and fresh-venv-imports cleanly on all 5 targeted
      CPython versions (3.9-3.13) via uv-managed interpreters, proving the Python-version
      half of SC4 without a manylinux container"
    requirement: KPI-06
    verification:
      - kind: other
        ref: "for v in 3.9 3.10 3.11 3.12 3.13: uv build --python $v --wheel && uv pip install --python /tmp/wheel-check-$v/bin/python <wheel> && /tmp/wheel-check-$v/bin/python import_check.py -> 5/5 'import OK', exit 0"
        status: pass
    human_judgment: false

duration: ~35min
completed: 2026-08-15
status: complete
---

# Phase 05 Plan 04: Python Wheel CI Matrix (cibuildwheel) and Local Build Proof Summary

**Authored the cibuildwheel-based manylinux/macOS x Python 3.9-3.13 CI matrix
(`.github/workflows/python-wheels.yml`, `[tool.cibuildwheel]`) and, while proving its
underlying mechanics locally without Docker, found and fixed a pre-existing sdist
packaging bug (`../src` unreachable from scikit-build-core's project-root-scoped sdist
walk) that would have broken `uv build`/`python -m build` and any future real PyPI
publish, not just this plan's local verification.**

## Performance

- **Duration:** ~35 min
- **Tasks:** 2 (plus Task 0, a `checkpoint:human-verify` approved in a prior session --
  see Checkpoint below)
- **Files modified:** 3 (1 created, 2 modified)

## Accomplishments

- `python/pyproject.toml` carries a `[tool.cibuildwheel]` table (`build = "cp39-* cp310-*
  cp311-* cp312-* cp313-*"`, `skip = "*-musllinux*"`, `test-command` importing `leafblower`
  and touching `leafblower.harvest`), plus `[tool.cibuildwheel.linux]` (`before-all`
  detects `yum` vs. `dnf` before installing `lapack-devel openblas-devel`, since
  cibuildwheel's default manylinux image can differ per CPython target) and
  `[tool.cibuildwheel.macos]` (`CMAKE_ARGS = "-DBLA_VENDOR=Apple"`, pinning CMake's
  `FindLAPACK` to Accelerate instead of its non-deterministic default provider search).
- `.github/workflows/python-wheels.yml` (new): `python-wheels` job builds wheels across
  `ubuntu-latest`/`macos-13`/`macos-14` via `pypa/cibuildwheel@v4.2.0` (checkpoint-approved
  exact tag, not a floating major-version ref) and uploads `wheelhouse/*.whl`;
  `wheel-check` (needs: `python-wheels`) downloads them, runs `twine check`, and installs
  the matching wheel into a fresh venv across Python 3.9-3.13 to confirm
  `import leafblower; leafblower.harvest` works.
- Ran the real local build/repair/import chain, no Docker: `uv build` (sdist+wheel) exits
  0, `twine check` PASSED, `auditwheel show`/`repair` vendored the sandbox's LAPACK
  provider (Intel MKL, 3 shared libs) into `leafblower.libs/`, producing a self-contained
  `manylinux_2_39_x86_64`-tagged wheel. Then proved the same chain across all 5 targeted
  CPython versions via `uv`-managed interpreters (3.9, 3.10, 3.11, 3.12, 3.13): each
  built, installed into a fresh venv, and imported clean -- 5/5.
- Found the sdist bug on the FIRST `uv build` run (not `uv build --wheel`): scikit-build-
  core's sdist packaging walks only the project directory (`python/`) and can never
  reach `../src` (CORE_SOURCES lives one level up, shared with the R build). `python -m
  build`/`uv build`'s default "build wheel from sdist" path therefore failed with `Cannot
  find source file: ../src/c_api.cpp`. Fixed via `[tool.scikit-build.sdist.force-include]`
  (vendors `../src`'s 40 git-tracked files under sdist-internal `src/`) +
  `sdist.exclude` (drops the *.o/*.so/Makevars build byproducts sitting alongside them)
  + a `LBW_SRC_DIR` CMake variable that detects which of the two layouts (`../src` for a
  dev checkout, `src` for an extracted sdist) is present.

## Task Commits

1. **Task 0: Confirm cibuildwheel/twine/build legitimacy before pinning** - checkpoint
   only, no commit (approved in prior session via `AskUserQuestion` -- see Checkpoint)
2. **Task 1: Author [tool.cibuildwheel] config and python-wheels.yml** - `40ad7e9` (feat)
3. **Task 2: Prove the build/repair/import mechanics locally, across Python 3.9-3.13** -
   `4954139` (fix) -- committed as `fix` rather than `test`/`chore` since the sdist bug
   fix is the substantive change; the proof itself left no artifacts (build outputs
   cleaned up after verification, per the task's own "no committed artifact" reversibility
   note)

## Files Created/Modified

- `.github/workflows/python-wheels.yml` - new; `python-wheels` (cibuildwheel matrix) +
  `wheel-check` (twine check + clean-venv import) jobs
- `python/pyproject.toml` - `[tool.cibuildwheel]` / `[tool.cibuildwheel.linux]` /
  `[tool.cibuildwheel.macos]` tables added; `[tool.scikit-build.sdist]` +
  `[tool.scikit-build.sdist.force-include]` added to fix the sdist bug
- `python/CMakeLists.txt` - `LBW_SRC_DIR` variable added, `CORE_SOURCES` and
  `target_include_directories` switched from hardcoded `../src` to `${LBW_SRC_DIR}`

## Decisions Made

See `key-decisions` in frontmatter. Summary: fixed the sdist bug at the root (force-include
+ CMake layout detection) instead of narrowing the local proof to `uv build --wheel` only;
kept `twine check`'s two pre-existing metadata WARNINGs out of scope (unrelated to this
plan, still `PASSED` not `FAILED`); accepted MKL as the vendored LAPACK provider for the
local proof since the mechanism, not the specific library, is what this task verifies.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] `uv build` (sdist+wheel) failed: sdist cannot reach `../src`**
- **Found during:** Task 2, first `uv build` run (before any per-version loop)
- **Issue:** `uv build` (default: sdist, then wheel-from-sdist, matching `python -m
  build`'s standard behavior) failed with a CMake configure error: `Cannot find source
  file: ../src/c_api.cpp`. `uv build --wheel` (which builds directly from the live source
  tree, skipping the sdist round-trip) succeeded, isolating the cause to scikit-build-
  core's sdist packaging: its file walk starts at `Path()` (the project root, `python/`)
  and can never enumerate files outside it, so the generated sdist silently omitted all
  17 `CORE_SOURCES` files (and their headers) living in the sibling `../src` directory.
  This is not scoped to this plan's local verification -- cibuildwheel's own build
  frontend and any future real `twine upload`/PyPI-publish step would hit the identical
  failure, since both consume the sdist.
- **Fix:** Added `[tool.scikit-build.sdist.force-include]` (`"../src" = "src"`, vendoring
  the 40 git-tracked files under `../src` into the sdist at `src/`) and
  `[tool.scikit-build.sdist].exclude` (`src/*.o`, `src/*.so`, `src/Makevars`,
  `src/Makevars.bak`, `src/lbw_config.h` -- the gitignored build byproducts that
  physically sit alongside the real sources and would otherwise also get force-included).
  Added a `LBW_SRC_DIR` CMake variable to `python/CMakeLists.txt`
  (`if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../src/c_api.cpp)` picks `../src` for a dev
  checkout, else falls back to `src` for an extracted sdist) and switched `CORE_SOURCES`
  and `target_include_directories` to use it instead of the hardcoded `../src`.
- **Files modified:** `python/pyproject.toml`, `python/CMakeLists.txt`
- **Verification:** `tar -tzf dist/*.tar.gz | grep '/src/'` shows exactly the 40
  git-tracked `src/` files, zero `.o`/`.so`/`Makevars` artifacts; `uv build` (full
  sdist+wheel) exits 0 and produces both `dist/leafblower-0.1.0.tar.gz` and a working
  wheel; a separate `cmake -S python -B <dir>` configure-only run against the real repo
  checkout (unaffected by the sdist changes) confirmed `LBW_SRC_DIR` still resolves to
  `../src` there, with no "cannot find source" regression for the existing dev workflow
  (`uv pip install -e . --reinstall-package leafblower`, per CLAUDE.md).
- **Committed in:** `4954139`

---

**Total deviations:** 1 (Rule 3, auto-fixed).
**Impact on plan:** The fix is a genuine correctness requirement for Task 2's literal
acceptance criteria (`uv build` exits 0) and for the wider phase goal (a real, working
Python packaging pipeline) -- not scope creep. It also closes a latent defect that would
have surfaced downstream in cibuildwheel's own build-frontend or a real PyPI publish,
neither of which this plan's own verify command would otherwise have caught before CI/
publish time.

## Issues Encountered

`uv run twine check ...` triggered `uv`'s project-sync machinery, which recreated
`python/.venv` from scratch (dropping the `twine`/`auditwheel` just installed via `uv pip
install`) because `twine` is not a declared project dependency. Worked around by invoking
`.venv/bin/python -m twine`/`.venv/bin/python -m auditwheel` directly against the venv
`uv pip install` populated, rather than `uv run`. `uv venv` also does not install `pip` by
default; used `uv pip install --python <venv>/bin/python <wheel>` instead of
`<venv>/bin/pip install` for the per-version fresh-venv installs.

## User Setup Required

None - no external service configuration required.

## Checkpoint: Task 0 (human-verify, gate=blocking-human)

Approved in a prior session via the interactive `AskUserQuestion` UI, relayed faithfully
by the orchestrator (not an orchestrator self-approval): pinning `cibuildwheel` 4.2.0
(`pypa/cibuildwheel@v4.2.0`), `twine` 7.0.0, and `build` 1.5.0 -- all three PyPA-maintained,
flagged `[SUS]` by the automated legitimacy gate purely on "too-new release"/"unknown
downloads" heuristics for high-cadence official tooling, with `build`'s
`github.com/pypa/build` repo link specifically `[ASSUMED]` (unverified by the automated
gate's repo lookup). User selected "Approved (Recommended)".

## Known Stubs

None -- no placeholder data or unwired functionality. One documented gap: `twine check`
reports `PASSED with warnings` (missing `long_description`/`long_description_content_type`
in `python/pyproject.toml`'s `[project]` table) -- pre-existing, unrelated to this plan's
scope, does not block `twine check`'s PASS verdict.

## Next Phase Readiness

- `[tool.cibuildwheel]` and `.github/workflows/python-wheels.yml` are ready for CI to
  exercise once this repository gains a git remote (this project is local-only per
  CLAUDE.md's Session Completion section) -- the manylinux/macOS OS-and-glibc portability
  layer stays genuinely untested until then, as stated honestly in the plan's own
  `<verification>` section.
- The sdist fix (`force-include` + `LBW_SRC_DIR`) is durable infrastructure any future
  real `twine upload`/PyPI-publish step (deferred by D-05) will depend on -- without it,
  publishing a source distribution alongside the wheels would have failed at build time.
- 05-05 (the phase-gate summary) should note the remaining CI-only gap (manylinux/macOS
  OS-and-glibc portability) alongside 05-01's parallel checkbashisms/tidy/V8 gap, both
  expected to close once the CI matrix actually runs on a standard runner image.

---
*Phase: 05-cran-pypi-release*
*Completed: 2026-08-15*

## Self-Check: PASSED

All claimed files (`python/pyproject.toml`, `.github/workflows/python-wheels.yml`,
`python/CMakeLists.txt`) and commits (`40ad7e9`, `4954139`) verified present.
