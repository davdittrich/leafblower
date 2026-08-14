# External Integrations

**Analysis Date:** 2026-08-14

## APIs & External Services

**None.** This is an offline numerical library. No HTTP client, no SDK, no network call exists in `R/`, `python/leafblower/`, or `src/`. Python imports are limited to `numpy`, `pandas`, `json`, `math`, `os`, `pathlib`, `subprocess`, `sys`, `warnings`, `typing` plus internal modules.

**Cross-language "integration" (the real coupling surface):**
- C ABI: `src/leafblower.h` — the stable C API surface (`rk_algorithm_t` enum, `CalibResult`) consumed by both bridges. Algorithm slot 2 is reserved (LBFGSB removed); do not reuse.
- R → C++: `src/r_bridge.cpp` via `.Call` with `useDynLib(leafblower, .registration = TRUE)` (`NAMESPACE`).
- Python → C++: `python/leafblower/_bindings.cpp` via pybind11, module `_leafblower`, wrapped by `python/leafblower/_harvest.py` and `python/leafblower/_design_effect.py`.
- R ↔ Python parity: `tests/parity/run_parity_r.R`, `run_oris_soft_r.R`, `run_chebyshev_r.R` are invoked from Python (`subprocess`) by `tests/test_parity_weights.py` and the `python/leafblower/test_*_parity.py` suite (rtol=1e-6).

## Data Storage

**Databases:**
- None. No DB driver, connection string, or ORM anywhere in the source tree.

**File Storage:**
- Local filesystem only. Test fixtures are R serialized objects: `tests/testthat/task1_ref.rds`, `tests/testthat/task2_oris_ref.rds`, `tests/testthat/fixtures/`. Raw data generation lives in `data-raw/`.
- Optional CSV trajectory output (`python/leafblower/test_trajectory_csv_smoke.py`, `src/oris_trajectory.cpp`).

**Caching:**
- None.

## Authentication & Identity

- Not applicable. A calibration library with no user model, no sessions, no auth.

## Monitoring & Observability

**Error Tracking:**
- None. Errors surface as R conditions / Python exceptions raised from the bridges; validation lives in `src/calib_validate.cpp` and `src/validation.cpp`.

**Logs:**
- No logging framework. Diagnostics are returned as result fields (`R/diagnose_weights.R`, `src/design_effect.cpp`) rather than emitted.

## CI/CD & Deployment

**Hosting:**
- None. Library only.

**CI Pipeline:**
- None. No `.github/workflows/`, no `.gitlab-ci.yml`, no `Jenkinsfile`, no `Dockerfile`.
- Quality gating is local and manual per `CLAUDE.md`: `R CMD INSTALL --preclean .` + testthat (0 FAIL) + pytest (0 FAIL) under single-thread BLAS; stepstone benchmark no-regression opt-in via `LBW_BENCH_GATE=1`.

**Git remote:**
- None — the repository is local-only. Work is complete when committed locally; there is nothing to push to.

## Environment Configuration

**No secrets, no credentials.** No `.env*`, `credentials.*`, key, or token file exists in the tree.

**Behavior-affecting env vars:**
- `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` — must be set to `1` together for deterministic R↔Python parity and benchmarks (in Python, before `import numpy`).
- `LBW_BENCH_GATE=1` — enables the heavy stepstone benchmark regression gate.
- `CXX`, `SHLIB_OPENMP_CXXFLAGS` — consumed at configure time by `configure` for compiler/OpenMP feature probing.
- `~/.R/Makevars` — where a user sets `-O3`; the package deliberately supplies no `-O` level for the R build.

## Webhooks & Callbacks

**Incoming:** None.

**Outgoing:** None.

## Optional Third-Party R Packages (test-time only)

Declared in `DESCRIPTION` `Suggests`, used as reference implementations rather than runtime integrations:
- `autumn (>= 0.2.0)` — the `harvest()` API this package is a drop-in replacement for.
- `survey`, `PracTools (>= 1.4.0)` — independent calibration/design-effect references for parity tests.
- `bench`, `lhs`, `DiceKriging`, `ggplot2`, `rprojroot` — benchmarking and study tooling (`benchmarks/`, `tools/check_research_isolation.R`).

---

*Integration audit: 2026-08-14*
