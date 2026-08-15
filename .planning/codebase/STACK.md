# Technology Stack

**Analysis Date:** 2026-08-15

## Languages

**Primary:**
- C++17 - Core calibration solvers and algorithms (`src/*.cpp`, `src/*.hpp`)
- R - User-facing package and wrapper functions (`R/*.R`)
- Python - Parallel bindings via pybind11 (`python/leafblower/*.py`, `python/leafblower/_bindings.cpp`)

**Secondary:**
- C - C API surface for cross-language compatibility (`src/leafblower.h`, `src/c_api.cpp`)
- Shell - Build configuration (`configure` script)

## Runtime

**Environment:**
- R: 3.6.0 or later (implied by NeedsCompilation: yes in DESCRIPTION)
- Python: 3.9, 3.10, 3.11, 3.12, 3.13 (per `python/pyproject.toml` classifiers)
- C++17 compiler required (auto-detected by `configure`, fallback to C++14 available)

**Package Manager:**
- R: Built-in (via `R CMD INSTALL`; `install.packages()` on CRAN)
- Python: uv (managed via `~/.venv/`, project uses uv for isolated venv)
- No lockfiles committed; dependencies specified declaratively only

## Frameworks

**Core:**
- No standalone framework; C++17 standard library only
- pybind11 2.11+ - Bindings generator (`python/CMakeLists.txt:7`)

**Build/Dev:**
- R CMD (R package build system) - Primary R build gate
- CMake 3.18+ - Python extension build via scikit-build-core
- scikit-build-core 0.8+ - Build backend for Python wheel (`python/pyproject.toml:2-3`)

**Testing:**
- R: testthat 3.0.0+ (`DESCRIPTION` Suggests)
- Python: pytest (via `python/pyproject.toml` optional-dependencies)
- No continuous integration configured (local-only repo; no remote)

## Key Dependencies

**Critical (Core Computation):**
- LAPACK (system library) - Linear algebra routines for Cholesky solve, eigenvector decomposition (`src/calib_linalg.cpp`, `src/newton_calib.cpp`)
  - `dpotrf`, `dpotrs`, `dsyevd` called directly
- BLAS (system library) - Dense matrix operations
  - Referenced via `$(BLAS_LIBS)` in `src/Makevars.in`

**Python Layer:**
- numpy 1.21+ - Numerical arrays and vectorization (`python/leafblower/_harvest.py`)
- pandas 1.3+ - Data frame operations and group-by aggregation (`python/leafblower/_harvest.py`)

**R Layer:**
- stats (base R) - Imported in `DESCRIPTION`
- autumn (>= 0.2.0) - Suggested; optional dependency for harvest() compatibility

**Supplementary (Optional, Testing/Benchmarking):**
- bench - R benchmarking framework
- lhs - Latin hypercube sampling
- DiceKriging - Gaussian process kriging
- ggplot2 - Visualization
- rprojroot - Project root detection
- survey - Survey design functions
- PracTools 1.4.0+ - Practical tools for complex surveys

**Optional SIMD/Vectorization:**
- glibc libmvec (system library, x86_64 only) - Vectorized `_ZGVdN4v_exp` (AVX2 4-wide double exp)
  - Auto-detected by `configure`; present on modern Linux x86_64 systems
  - Fallback: scalar `exp()` from libm

## Configuration

**Environment:**
- Environment variables (if set before build):
  - `CXX` - C++ compiler command (used by `configure` for feature detection)
  - `CXXFLAGS` - User/site compilation flags (user-supplied, not set by package)
  - `SHLIB_OPENMP_CXXFLAGS` - R's OpenMP flag template (auto-detected by R)
  - `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` - Threadpool control (REQUIRED for deterministic R↔Python parity; set to 1 for testing)

**Build:**
- `src/Makevars.in` - R-side compiler flags (substituted by `configure`)
  - Decorated placeholders: `@CXXFLAGS_STD@`, `@OMP_FLAGS@`, `@SIMD_FLAGS@`, `@MAVX2_FLAG@`, `@MVEC_LIBS@`
  - No `-O` optimization flags (user-supplied via `~/.R/Makevars`)
  - `configure` generates `src/lbw_config.h` with feature detection (`LBW_HAS_OMP_SIMD`, `LBW_HAS_OMP`, `LBW_HAS_GLIBC_MVEC`)

- `python/CMakeLists.txt` - Python extension build
  - Hermetic feature detection (re-probes to ensure Python module independence from R build artifacts)
  - Forces `-O3 -mavx2` (x86_64 only; ARM/no-AVX2 omits `-mavx2`)
  - Generates `python/build/_/generated/lbw_config.h` (CMake-generated, isolated)

- `python/pyproject.toml` - Python packaging metadata
  - `scikit-build-core` backend with `cmake.build-type = "Release"`
  - Excludes test files from wheel: `wheel.exclude = ["leafblower/test_*.py"]`

## Platform Requirements

**Development:**
- C++17 capable compiler (g++, clang, Apple Clang)
- CMake 3.18+ (for Python builds)
- R development headers (linux: `r-base-dev`, macOS: XCode Command Line Tools)
- LAPACK/BLAS development libraries (linux: `libblas-dev`, `liblapack-dev`; macOS: Accelerate framework)
- Python 3.9+ with development headers (`python3-dev` on linux)
- glibc 2.22+ with libmvec (optional; present on most modern Linux x86_64)
- OpenMP runtime (optional; auto-detected; improves SIMD vectorization hints)

**Production:**
- Deployment: R package via CRAN binary or source build
  - R 3.6+, compiled against system LAPACK/BLAS
- Python package via PyPI wheel or source
  - Python 3.9-3.13, bundled C++ core
- No server deployment; local compute library only
- Single-machine, single-threaded by default (thread control via environment variables)

## Special Notes

**No Hard-Coded Optimization Flags:**
- Package follows R Packaging Guidelines (WRE §1.2.1): does NOT set `-O` levels
- User/site `-O` level comes via R's `$(CXXFLAGS)` environment
- Python build sets `-O3 -mavx2` to match R's expected performance tier
- A user wanting `-O3` in R sets it in `~/.R/Makevars`

**Feature Detection at Build Time:**
- `./configure` probes for C++17, OpenMP, `-fopenmp-simd`, glibc libmvec
- Output written to `src/lbw_config.h` (used by R build)
- Python CMake re-probes independently (`CMakeLists.txt` hermetic feature detection) to avoid dependency on gitignored R artifacts

**Parity Requirement (R↔Python):**
- Deterministic weights require single-thread BLAS: export `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` before any R/Python import
- Both builds link same LAPACK/BLAS system libraries but may differ on OpenMP availability
- Test suite runs under thread constraint to guarantee bit-parity

**No LTO (Link-Time Optimization):**
- Solvers and helper functions intentionally kept in single or paired translation units
- Cross-TU inlining (required for hot-path performance in iterative solvers) is disabled by design
- `oris.cpp` keeps main solve loop co-located; cold finalize logic split to `oris_finalize.cpp`

---

*Stack analysis: 2026-08-15*
