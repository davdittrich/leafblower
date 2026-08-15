# External Integrations

**Analysis Date:** 2026-08-15

## Overview

**leafblower** is a standalone computational library with **no external API integrations, no cloud dependencies, and no network requirements.** All computation is local; all I/O is file-based or in-memory.

## APIs & External Services

**None.** The package does not call any external APIs, cloud services, or web endpoints.

**No authentication providers, no webhooks, no external monitoring.**

## Data Storage

**Local Filesystem Only:**
- Input: data frames / arrays passed in-memory (R data.frame, Python pandas.DataFrame, or C arrays via C API)
- Output: calibrated weights returned as vectors in-memory; no automatic file persistence
- Benchmarks: reference fixtures stored as serialized files (`.rds` for R, `.json` for benchmark metadata)
  - `benchmarks/stepstone_bench_targets.json` - Target definitions
  - `benchmarks/stepstone_fulldata_bench_targets.json` - Full-dataset benchmarks
  - Test fixtures in `tests/testthat/fixtures/` - Pre-computed .rds files for regression testing

**No databases, no persistent backend, no cache server.**

**Optional trajectory logging:**
- `oris.cpp` writes trajectory CSV to user-specified file path (if requested via C API)
  - Uses C stdio (`fopen`/`fprintf`/`fclose`), not C++ streams
  - Deterministic output: one row per iteration, tab-separated values

## Authentication & Identity

**Not applicable.** No user accounts, no API keys, no session management.

## Monitoring & Observability

**Error Tracking:**
- None; errors reported via return codes and status messages in `rk_result_t` struct (`leafblower.h`)
- Solver messages set in `rk_result_t.message[]` (256-byte string)

**Logs:**
- Optional callback logging via `rk_params_t.log_fn` (C API)
  - Function pointer: `void (*log_fn)(const char* msg, void* ctx)`
  - R wrapper prints to stderr on `verbose >= 1`
  - Python wrapper logs via Python's logging module (or stderr on verbose)
- No persistent log files; all output transient

**Metrics:**
- Solver diagnostics returned in `rk_result_t`: iterations, convergence metric, error progression, SOR adaptation stats
- Accessed by user code; not sent anywhere

## CI/CD & Deployment

**Hosting:**
- Not hosted. This is a library; users install via CRAN (R) or PyPI (Python).
- Source repository is local-only (no remote per `CLAUDE.md`: "No git remote (local-only)")

**CI Pipeline:**
- None configured in repo
- Local quality gates (defined in `CLAUDE.md`):
  - R: `R CMD INSTALL --preclean .` + `Rscript -e "devtools::test()"`
  - Python: `pytest` with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`
  - Optional stepstone no-regression check: `LBW_BENCH_GATE=1` (heavy, local-only)

**Distribution:**
- CRAN (R package): binary or source builds on submission
- PyPI (Python package): wheels via scikit-build-core + GitHub Actions (user-managed; not in this repo)
- No continuous deployment

## Environment Configuration

**Required environment variables (for deterministic testing):**
- `OMP_NUM_THREADS=1` - Disable OpenMP threading
- `OPENBLAS_NUM_THREADS=1` - Disable OpenBLAS threading
- `MKL_NUM_THREADS=1` - Disable Intel MKL threading
- (Set all three TOGETHER, BEFORE importing Python/R, to guarantee bit-parity R↔Python)

**Optional:**
- `LBW_BENCH_GATE=1` - Enable expensive stepstone benchmark regression gate (local use only)
- `CXX` - Override C++ compiler (for `./configure` feature detection)
- `CXXFLAGS` - Override compiler flags (R uses these; package does NOT add `-O`)

**Secrets location:**
- No secrets. No credentials stored, no API keys, no .env files required.

## Webhooks & Callbacks

**Incoming:**
- None

**Outgoing:**
- Optional solver progress callback (C API only):
  - `rk_params_t.log_fn` user-supplied function pointer
  - Called during solver iterations for logging / monitoring purposes
  - No external network usage

## R Package Dependencies

**Hard Imports:**
- `stats` (base R) - Standard distributions, variance/covariance functions

**Suggested (optional):**
- `autumn (>= 0.2.0)` - For harvest() API compatibility comparison
- `testthat (>= 3.0.0)` - Unit test framework
- `bench` - Benchmarking tools
- `lhs` - Latin hypercube sampling (benchmarks)
- `DiceKriging` - Gaussian process fitting (benchmarks)
- `ggplot2` - Visualization (benchmarks/diagnostics)
- `rprojroot` - Project directory detection
- `survey` - Survey design / variance functions (tests)
- `PracTools (>= 1.4.0)` - Practical survey tools (tests)

## Python Package Dependencies

**Core:**
- `numpy >= 1.21` - Numerical arrays
- `pandas >= 1.3` - Data manipulation

**Build:**
- `scikit-build-core >= 0.8` - Build backend
- `pybind11 >= 2.11` - C++ bindings

**Optional:**
- `pytest` - Test framework (installed via `test` extra: `pip install -e ".[test]"`)

## System Library Requirements

**Required (system-provided):**
- LAPACK (liblapack.so / Accelerate on macOS)
- BLAS (libblas.so / Accelerate on macOS)
- libm (math library, standard on all POSIX systems)

**Optional (system-provided, auto-detected):**
- libmvec (glibc >= 2.22, x86_64) - Vectorized exponential for AVX2 SIMD
- OpenMP runtime (libgomp for GCC, libomp for Clang)
- `-fopenmp-simd` compiler support (SIMD vectorization hints)

## No External Service Dependencies

**Summary:**
- No HTTP/REST calls
- No database connections
- No message queues
- No cloud storage
- No SaaS integrations
- No third-party authentication
- No rate-limited APIs
- No subscription services

The library runs entirely in-process, on the user's machine, with no network access or external dependencies beyond standard system libraries (LAPACK/BLAS).

---

*Integration audit: 2026-08-15*
