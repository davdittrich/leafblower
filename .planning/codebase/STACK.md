# Technology Stack

**Analysis Date:** 2026-08-14

## Languages

**Primary:**
- C++17 — the entire numerical core in `src/` (~11k LOC of `.cpp`/`.hpp`). Standard is probed by `configure` (`-std=c++17`, falling back to `-std=c++14` with a warning) and fixed by `set(CMAKE_CXX_STANDARD 17)` in `python/CMakeLists.txt`.
- R (>= no explicit floor in `DESCRIPTION`) — user-facing API in `R/`: `R/harvest.R`, `R/anesrake.R`, `R/design_effect.R`, `R/diagnose_weights.R`, `R/current_miss.R`, `R/na_bin.R`, `R/weighted_pct.R`, `R/zzz.R`.
- Python >= 3.9 (`python/pyproject.toml` `requires-python`) — thin wrapper layer `python/leafblower/_harvest.py`, `python/leafblower/_design_effect.py`, `python/leafblower/__init__.py`.

**Secondary:**
- Shell — `configure` (POSIX `/bin/sh` feature-detection script).
- CMake >= 3.18 — `python/CMakeLists.txt`.
- Make — `src/Makevars.in` (configure-substituted into `src/Makevars`).

## Runtime

**Environment:**
- R with a compiled shared object (`useDynLib(leafblower, .registration = TRUE)` in `NAMESPACE`), `NeedsCompilation: yes`.
- CPython >= 3.9 loading the pybind11 extension module `_leafblower` (`python/leafblower/_bindings.cpp`).

**Package Manager:**
- R: base `R CMD INSTALL` (no renv/packrat lockfile present).
- Python: `uv` — the project venv is uv-managed; `uv pip install -e . --reinstall-package leafblower` per `CLAUDE.md`. No `uv.lock` or `requirements.txt` at repo root.
- Lockfile: missing for both languages. A root `package-lock.json` exists but `package.json` is literally `{}` — Node tooling is vestigial, not part of the build.

## Frameworks

**Core:**
- pybind11 >= 2.11 — Python/C++ binding layer (`find_package(pybind11 2.11 REQUIRED)`, `pybind11_add_module(_leafblower ...)`).
- scikit-build-core >= 0.8 — PEP 517 build backend (`[build-system]` in `python/pyproject.toml`, `build-backend = "scikit_build_core.build"`).
- R `.Call`/`Rinternals` bridge — `src/r_bridge.cpp` (excluded from the Python build; it is the only `src/*.cpp` NOT in `CORE_SOURCES`).

**Testing:**
- testthat (>= 3.0.0, edition 3) — R suite in `tests/testthat/` (~100 `test-*.R` files), driver `tests/testthat.R`.
- pytest — Python suite lives *inside* the package dir (`python/leafblower/test_*.py`, 20 files) so an editable install discovers it; excluded from wheels via `wheel.exclude = ["leafblower/test_*.py"]`.
- Cross-language parity harness: `tests/parity/run_parity_r.R`, `run_oris_soft_r.R`, `run_chebyshev_r.R`, `tests/test_parity_weights.py`, root `conftest.py`.
- C++-level tests: `tests/cpp/`.

**Build/Dev:**
- `configure` + `src/Makevars.in` — autoconf-style hand-rolled feature probes substituting `@CXXFLAGS_STD@ @OMP_FLAGS@ @SIMD_FLAGS@ @MAVX2_FLAG@ @MVEC_LIBS@`.
- roxygen2 7.3.3 (`RoxygenNote`) — generates `man/` and `NAMESPACE`.

## Key Dependencies

**Critical:**
- LAPACK / BLAS — linked in both build sites: `PKG_LIBS = @MVEC_LIBS@ $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` (`src/Makevars.in`) and `find_package(LAPACK REQUIRED)` (`python/CMakeLists.txt`). Used by `src/calib_linalg.cpp`.
- numpy >= 1.21, pandas >= 1.3 — the only Python runtime dependencies (`python/pyproject.toml`).
- R `stats` — sole R import (`Imports: stats`; `importFrom("stats", "setNames")`).
- glibc `libmvec` (`-lmvec`, symbol `_ZGVdN4v_exp`) — optional AVX2 4-wide vector `exp`, gates the `_mm256_*` intrinsics in `src/lbw_math.hpp` behind `LBW_HAS_GLIBC_MVEC`.
- OpenMP — optional (`SystemRequirements: C++17 compiler, OpenMP (optional, for SIMD exp vectorisation)`). Probed in two tiers: `-fopenmp-simd` hints only, or full `SHLIB_OPENMP_CXXFLAGS` threading.

**Suggested (R, test/benchmark only):**
- `autumn (>= 0.2.0)` — the package this is a drop-in replacement for (`harvest()` API compatibility reference).
- `testthat`, `bench`, `lhs`, `DiceKriging`, `ggplot2`, `rprojroot`, `survey`, `PracTools (>= 1.4.0)` — reference implementations for parity tests and benchmark tooling.

## Configuration

**Environment:**
- No `.env`, `.envrc`, or secrets file exists in the repo.
- Determinism env vars are load-bearing for tests/benchmarks: `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` must all be `1` together (in Python, exported before `import numpy`), else R↔Python parity drifts.
- `LBW_BENCH_GATE=1` opts into the heavy stepstone no-regression benchmark gate.

**Build:**
- `configure` → writes `src/Makevars` and the R-side `src/lbw_config.h` (gitignored, machine-specific).
- `python/CMakeLists.txt` + `python/lbw_config.h.in` → generates a *hermetic* `lbw_config.h` in the CMake binary dir, force-included via `-include` so a stale R-generated `src/lbw_config.h` cannot shadow it.
- Optimization policy (verified in both files): the package sets **no `-O` level** for the R build — R supplies it via `$(CXXFLAGS)` in `$(ALL_CXXFLAGS)`; `tools:::.check_make_vars` rejects `-O*` in `PKG_CXXFLAGS`. The **Python** build does set `-O3` explicitly (`target_compile_options(_leafblower PRIVATE -O3)`), since CRAN rules do not apply there. `-mavx2` is added on `x86_64|AMD64` only, in both build sites.
- **Two build sites for `src/*.cpp`:** R auto-globs `src/*.cpp` (the `PKG_SOURCES` list in `Makevars.in` is explicitly commented "DECORATIVE ONLY"). The Python build uses an explicit `CORE_SOURCES` list (17 files, all of `src/*.cpp` except `r_bridge.cpp`). A new `src/*.cpp` must be added to `CORE_SOURCES` or the pybind11 link fails with undefined symbols.
- No LTO: `-flto` appears in neither `configure` nor `Makevars.in`.
- Python build defines `LBW_NO_R` to compile the core without R headers.

**Build & test commands (from `CLAUDE.md`, matching the files):**
```bash
R CMD INSTALL --preclean .                                        # R build gate (NOT devtools::install)
Rscript -e "devtools::test()"                                     # R tests
cd python && uv pip install -e . --reinstall-package leafblower    # Python build (uv, never pip)
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest
```

## Platform Requirements

**Development:**
- C++17 compiler (GCC/Clang), R with dev toolchain, Python >= 3.9 with `uv`, CMake >= 3.18, LAPACK/BLAS headers and libs.
- Full SIMD/`libmvec` fast path requires glibc >= 2.22 on x86_64 with AVX2. ARM and non-AVX2 builds compile cleanly with the intrinsics `#ifdef`-guarded out.

**Production:**
- Distributed as a source R package (`leafblower_0.1.0.tar.gz`) and a Python wheel built by scikit-build-core. No container, no server, no deployment target — this is a library.

## Tooling Directories (not part of the stack)

`.wolf/`, `.beads/`, `.metaswarm/`, `graphify-out/`, `tasks/`, `benchmarks/study/`, `node_modules/` are agent/issue-tracking/analysis tooling, not build inputs. Several `*.Rcheck/` dirs and `src/*.o`/`src/leafblower.so` are build artifacts checked into the working tree.

---

*Stack analysis: 2026-08-14*
