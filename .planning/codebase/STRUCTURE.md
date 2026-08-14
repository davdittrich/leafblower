# Codebase Structure

**Analysis Date:** 2026-08-15

## Directory Layout

```
/home/dd/Gemini/leafblower/
├── src/                          # C++17 core library (compiled by both R and Python builds)
│   ├── leafblower.h              # C99-clean public API header (ABI contract)
│   ├── c_api.cpp                 # Main dispatcher: algorithm selection, result packing
│   ├── types.hpp                 # Shared structs (CalibState, CalibResult, enums)
│   ├── validation.hpp/.cpp       # Input validation (dimensions, bounds, K_MAX)
│   │
│   ├── cell_table.hpp/.cpp       # Cross-tabulation: cell mapping, M_cell estimation
│   ├── calib_dispatch.hpp        # Convergence rules & metrics (shared inline helpers)
│   ├── calib_validate.hpp        # Bounds checking utils (used by solvers)
│   ├── calib_linalg.hpp/.cpp     # Matrix ops: Cholesky, TSVD, eigendecomposition
│   ├── lbw_math.hpp              # Numeric kernels: logsumexp, stable arithmetic
│   ├── lbw_config.h              # Build-time configuration (OpenMP, compiler flags)
│   │
│   ├── oris.hpp/.cpp             # ORIS (block-coordinate descent, box-constrained Sinkhorn)
│   ├── oris_internal.hpp         # ORIS internals (gates, SOR, step control)
│   ├── oris_finalize.cpp         # ORIS weight finalization + bounds clamping
│   ├── oris_trajectory.cpp       # ORIS trajectory logging (CSV output)
│   │
│   ├── raking.hpp/.cpp           # Raking (IPF with acceleration)
│   ├── sinkhorn.hpp/.cpp         # Sinkhorn (entropy-minimizing IPF)
│   ├── greg.hpp/.cpp             # GREG (Newton QP, chi-square minimization)
│   ├── chebyshev.hpp/.cpp        # Chebyshev (interior-point LP, minimax)
│   ├── greenkhorn.hpp/.cpp       # Greenkhorn (greedy coordinate descent IPF)
│   ├── logit_calib.hpp/.cpp      # Logit Newton calibration (Deville-Sarndal)
│   ├── newton_calib.hpp/.cpp     # Newton-KL (smooth dual, TSVD damping)
│   ├── design_effect.hpp/.cpp    # Kish & Henry-Valliant design effect
│   │
│   └── r_bridge.cpp              # R entry point (.Call interface, roxygen2 wrapping)
│
├── R/                            # R user-facing code
│   ├── harvest.R                 # Main wrapper function (S3, result wrapping, diagnostics)
│   ├── design_effect.R           # Design effect wrapper + CI computation
│   ├── diagnose_weights.R        # Weight diagnostics (bounds violations, summary stats)
│   ├── na_bin.R                  # NA binning utility
│   ├── weighted_pct.R            # Weighted percentage computation
│   ├── anesrake.R                # anesrake-compatible backward-compat wrapper
│   ├── current_miss.R            # Missing data imputation (experimental)
│   └── zzz.R                     # Package initialization (.onLoad, roxygen2 hook)
│
├── python/                       # Python build system & bindings
│   ├── CMakeLists.txt            # scikit-build2 configuration (compiles src/*.cpp)
│   ├── pyproject.toml            # PEP 517 build config (setuptools, scikit-build2)
│   ├── leafblower/               # Python package
│   │   ├── __init__.py           # Module entry (imports _leafblower.so, exposes calibrate)
│   │   ├── _bindings.cpp         # pybind11 module (calibrate function, callbacks)
│   │   ├── _harvest.py           # High-level Python API (calibrate wrapper, result unpacking)
│   │   ├── _design_effect.py     # Design effect wrapper
│   │   ├── test_*.py             # Unit tests (30+ files: regression tests, parity checks, CR fixes)
│   │   └── __pycache__/          # Bytecode cache (ignored by .gitignore)
│   ├── parity/                   # R↔Python parity fixtures (reference .rds files)
│   └── .venv/                    # uv-managed Python venv (not committed)
│
├── tests/                        # R & Python test suite
│   ├── testthat/                 # R testthat v3 tests
│   │   ├── fixtures/             # .rds reference files (cached solver outputs)
│   │   ├── test_*.R              # 50+ test files (algorithm correctness, regression, CR fixes)
│   │   ├── test-calibration-solvers.R        # Main solver parity suite
│   │   ├── test-convergence-criteria.R       # Convergence rule/metric tests
│   │   ├── test-calib-linalg.R               # Linear algebra correctness
│   │   ├── test-design-effect.R              # Design effect computation
│   │   └── ...
│   ├── testthat.R               # testthat runner (single line: library(testthat); ...)
│   ├── cpp/                      # C++ unit tests (legacy; mostly dormant)
│   ├── parity/                   # Python parity runner
│   └── test_parity_weights.py    # Python weight parity checks
│
├── benchmarks/                   # Performance/regression gate
│   ├── stepstone.R               # Benchmark harness (LBW_BENCH_GATE=1)
│   ├── cases/                    # Benchmark case definitions (n, K, data complexity)
│   └── results/                  # Historical benchmark results (csv)
│
├── docs/                         # Documentation & specs
│   ├── superpowers/              # Phase design docs, architecture notes
│   ├── guides/                   # User guides (WIP)
│   └── ...
│
├── .planning/                    # GSD planning artifacts
│   ├── codebase/                 # This directory (ARCHITECTURE.md, STRUCTURE.md, etc.)
│   └── phases/                   # Phase plans, specs
│
├── .beads/                       # Beads issue tracker (local-only)
│   ├── issues.jsonl              # Issue database
│   ├── plans/                    # Linked work plans
│   └── ...
│
├── .claude/                      # Claude Code local config
│   ├── settings.json             # Tool permissions, model routing
│   ├── rules/                    # Project-specific rules
│   └── ...
│
├── .metaswarm/                   # Orchestrated execution config
│   ├── project-profile.json      # Metaswarm profile (phases, gating)
│   └── external-tools.yaml       # Delegated executor config (gemini, codex)
│
├── DESCRIPTION                   # R package metadata (CRAN-style)
├── NAMESPACE                     # R export list (roxygen2 artifact)
├── LICENSE                       # MIT license
├── CLAUDE.md                     # Project instructions (this file, checked in)
├── .gitignore                    # VCS exclusions (build artifacts, .venv, .Rcheck)
├── .coverage-thresholds.json     # Quality gate config (0 FAIL enforcement)
├── configure                     # GNU autoconf script (optional; R uses Makevars)
└── README.md                     # User-facing intro (if present)
```

## Directory Purposes

**`src/`:**
- Purpose: C++17 core library (single compiled translation unit per solver, linked by both R and Python)
- Contains: Algorithm implementations, calibration logic, shared infrastructure
- Key files: `leafblower.h` (public API), `c_api.cpp` (dispatcher), all solver `*.cpp`
- R build: auto-globs all `src/*.cpp` (see `Makevars.in`)
- Python build: explicit `CORE_SOURCES` list in `CMakeLists.txt` (must update when adding new solver)
- Output: Static library (`src/.libs/leafblower.a`), linked into both R package and Python wheel

**`R/`:**
- Purpose: R-facing wrappers, roxygen2 documentation, user-convenient functions
- Contains: `harvest()` S3 function, design effect wrapper, weight diagnostics
- Generated: `NAMESPACE`, `man/*.Rd` files (roxygen2 artifact)
- Build: `R CMD INSTALL --preclean .` compiles src/ then installs R/ functions
- Pattern: Each `.R` file is a separate module; roxygen2 docstrings (`#' @export`) drive NAMESPACE

**`python/`:**
- Purpose: Python build system, pybind11 bindings, high-level Python API
- Contains: `CMakeLists.txt` (scikit-build2), `_bindings.cpp` (pybind11 module), `leafblower/` package
- Build: `pip install -e .` or `uv pip install -e .` invokes scikit-build2, compiles src/, creates `_leafblower.so`, installs package
- Output: `_leafblower.so` (compiled module), `leafblower/` (Python code)
- Pattern: `__init__.py` imports `_leafblower` and re-exports `calibrate()`; `_harvest.py` wraps for high-level API

**`tests/testthat/`:**
- Purpose: R test suite (testthat v3)
- Contains: 50+ test files covering solver correctness, convergence rules, regression tests, CR fixes
- Fixtures: `.rds` files with reference outputs (cached solver results)
- Run: `devtools::test()` or `Rscript -e "devtools::test()"`
- Output: Pass/fail per test; test report in console

**`python/leafblower/test_*.py`:**
- Purpose: Python test suite (pytest-based)
- Contains: 30+ test files mirroring R testthat suite; parity checks (R vs. Python); CR-specific regression tests
- Run: `pytest` or `OMP_NUM_THREADS=1 pytest` (single-thread BLAS for determinism)
- Output: Pass/fail; can measure coverage via `pytest --cov`

**`benchmarks/`:**
- Purpose: Performance regression gate (LBW_BENCH_GATE=1 opt-in)
- Contains: `stepstone.R` harness, benchmark case definitions, historical results
- Run: `source("benchmarks/stepstone.R")` or integrated into CI
- Threshold: Time/iteration parity within 5% of baseline (project-specific)

**`.planning/codebase/`:**
- Purpose: GSD codebase maps (this directory)
- Contains: ARCHITECTURE.md, STRUCTURE.md, CONVENTIONS.md, TESTING.md, CONCERNS.md
- Generated by: `/gsd-map-codebase` skill
- Consumed by: `/gsd-plan-phase`, `/gsd-execute-phase` skills

**`.beads/`:**
- Purpose: Local issue tracker (no remote)
- Contains: `issues.jsonl` (work items), plans, session state
- Commands: `bd ready`, `bd show <id>`, `bd update <id> --claim`, `bd close <id>`
- Manual: `bd prime` for full workflow

## Key File Locations

**Entry Points:**
- `src/leafblower.h`: C API (rk_calibrate, rk_params_t, rk_result_t, enums)
- `src/c_api.cpp`: Algorithm dispatcher (line ~190 for main rk_calibrate)
- `src/r_bridge.cpp`: R .Call entry point (C_rk_calibrate native routine, ~400 lines)
- `python/leafblower/_bindings.cpp`: Python pybind11 module (calibrate function, ~300 lines)

**Configuration:**
- `DESCRIPTION`: R package metadata (version, deps, build flags)
- `python/pyproject.toml`: Python package metadata (version, build config)
- `python/CMakeLists.txt`: Compilation config (CORE_SOURCES list, include paths)
- `src/lbw_config.h`: Build-time flags (OpenMP, SIMD)
- `configure`: GNU autoconf (optional; R prefers Makevars)

**Core Logic:**
- `src/oris.cpp`: ORIS solver (main iterative algorithm, ~800 lines)
- `src/oris_internal.hpp`: ORIS step control, gate logic
- `src/raking.cpp`: Raking solver (~600 lines)
- `src/calib_dispatch.hpp`: Convergence helpers (used by all solvers)
- `src/cell_table.cpp`: Cross-tabulation builder (~400 lines)
- `src/validation.cpp`: Input validation (~300 lines)

**Testing:**
- `tests/testthat/test-calibration-solvers.R`: Main R solver parity suite (~1500 lines)
- `python/leafblower/test_solver_parity.py`: Main Python solver parity suite
- `python/leafblower/test_python.py`: Python-specific integration tests (~800 lines)
- `benchmarks/stepstone.R`: Performance regression harness

## Naming Conventions

**Files:**
- Solvers: `<algorithm>.hpp`, `<algorithm>.cpp` (e.g., `oris.hpp`, `raking.cpp`)
- Helpers: `<area>.hpp`, `<area>.cpp` (e.g., `calib_dispatch.hpp`, `cell_table.cpp`)
- R functions: `<domain>.R` (e.g., `harvest.R`, `design_effect.R`)
- Tests: `test-<area>.R` (R) or `test_<area>.py` (Python)
- C API types: `rk_*` prefix (e.g., `rk_calibrate`, `rk_params_t`, `rk_algorithm_t`)
- C++ namespace: `lbw` (leafblower internal)

**Functions:**
- Public C API: lowercase with underscores (e.g., `rk_calibrate`, `rk_design_effect`)
- Solver entry points: `<name>_solve(CalibState&)` → `<Name>Result` (e.g., `oris_solve`, `raking_solve`)
- Helpers: lowercase (e.g., `build_cell_table`, `validate_inputs`, `select_metric`)
- R functions: lowercase (e.g., `harvest`, `design_effect`)

**Types:**
- Public C structs: `rk_*` or `rk_*_t` (e.g., `rk_params_t`, `rk_result_t`)
- C++ classes/structs: PascalCase (e.g., `CalibState`, `ORISResult`, `CellTable`)
- Enums: PascalCase with `_t` suffix (e.g., `rk_algorithm_t`, `rk_bounds_mode_t`)
- C++ enums (ns lbw): PascalCase (e.g., `CalibMetric`, `CalibRule`)

**Variables:**
- camelCase for local variables (e.g., `cellTable`, `maxError`)
- UPPER_CASE for constants (e.g., `K_MAX`, `kErrCheckInterval`)
- Prefixes: `n_*` for counts, `is_*` for bools, `*_of` for mappings

## Where to Add New Code

**New Solver Algorithm:**
1. Primary implementation: `src/<algorithm>.hpp` (result struct + solver interface)
2. Core loop: `src/<algorithm>.cpp` (solver implementation, ~800 lines)
3. Header entry point: `src/leafblower.h` — add enum value to `rk_algorithm_t` (ensure slot 2 remains reserved for removed LBFGSB)
4. Dispatcher: `src/c_api.cpp` — add branch to rk_calibrate (algorithm selection switch, result packing)
5. Algorithm name: Add to `kAlgNames[]` array in `src/c_api.cpp` (must match enum order)
6. R wrapper: `R/harvest.R` — add algorithm string to switch statement
7. Python wrapper: `python/leafblower/_bindings.cpp` (params dict → rk_params_t)
8. R tests: `tests/testthat/test-calibration-solvers.R` — add fixture + parity check
9. Python tests: `python/leafblower/test_solver_parity.py` — add parity test
10. Build system: `python/CMakeLists.txt` — add `src/<algorithm>.cpp` to `CORE_SOURCES` list (R auto-globs; Python does NOT)

**New Convergence Rule:**
1. Enum: Add to `CalibRule` in `src/types.hpp`
2. Logic: `src/calib_dispatch.hpp::apply_rule()` — add case to switch
3. R param: `R/harvest.R` — map string name to enum
4. Python param: `python/leafblower/_bindings.cpp` — dict key handling
5. Tests: `tests/testthat/test-convergence-criteria.R` + Python parity

**New Calibration Metric:**
1. Enum: Add to `CalibMetric` in `src/types.hpp`
2. Selection: `src/calib_dispatch.hpp::select_metric()` — add parameter and case
3. Computation: Solver must compute all 7 metrics (max_err, mean_err, kl, chi2, grake_norm, l1_weight, marginal_kl) at each check; metric selection happens at check time (no solver-specific filtering)
4. Result field: Already exists in `CalibResult`; no struct change needed

**New Configuration Parameter:**
1. Struct field: Add to `rk_params_t` in `src/leafblower.h`
2. ABI tripwire: Update `EXPECTED_RK_PARAMS_BYTES` (track padding)
3. Initialization: Add default to `rk_params_init()` in `c_api.cpp`
4. Marshaling: R bridge (`r_bridge.cpp`) and Python bindings (`_bindings.cpp`) must extract from their input formats
5. Solver usage: Solvers read from `CalibState` (which mirrors rk_params_t); ensure field is transferred during state marshaling
6. Documentation: Update comment in header and CLAUDE.md architecture notes

**New Utility/Helper (e.g., stabilized computation):**
- Shared across solvers: `src/lbw_math.hpp` (inline) or `src/calib_linalg.cpp` (function)
- Solver-specific: `src/<algorithm>_internal.hpp` or static fn in `src/<algorithm>.cpp`
- Matrix ops: `src/calib_linalg.hpp` and `src/calib_linalg.cpp`

**New Test (R):**
- Location: `tests/testthat/test-<domain>.R`
- Pattern: `test_that("<feature>", { expect_*(...) })`
- Fixtures: Save reference outputs to `tests/testthat/fixtures/<name>.rds` via `saveRDS(result, "tests/testthat/fixtures/<name>.rds")`
- Run: `devtools::test()` or `testthat::test_dir("tests/testthat")`

**New Test (Python):**
- Location: `python/leafblower/test_<domain>.py`
- Pattern: `def test_<feature>(): assert ... or pytest.approx(...)`
- Run: `pytest python/leafblower/test_<domain>.py` or full suite via `pytest`
- Parity: Compare R output to Python via `test_solver_parity.py` pattern

## Special Directories

**`tests/testthat/fixtures/`:**
- Purpose: Cached reference outputs (.rds files)
- Generated: Solver runs with known inputs produce reference outputs
- Usage: Regression tests compare new solver runs to fixture outputs
- Committed: Yes (ensures reproducibility across time)
- Pattern: `fixture_<case>_<algorithm>_<config>.rds`

**`benchmarks/results/`:**
- Purpose: Historical benchmark data (CSV)
- Generated: `stepstone.R` harness appends new run results
- Usage: Performance regression detection (time/iter drift > 5%)
- Committed: Selective (not all runs; major versions only)

**`.planning/codebase/`:**
- Purpose: GSD documentation (ARCHITECTURE.md, STRUCTURE.md, etc.)
- Generated: `/gsd-map-codebase` skill
- Committed: Yes (input to `/gsd-plan-phase` and `/gsd-execute-phase`)
- Versioning: Date stamps in headers (refreshed on each remap)

**`.beads/`:**
- Purpose: Issue tracking state
- Generated: `bd` commands update `issues.jsonl`
- Committed: Yes (work history)
- Usage: `bd ready`, `bd show <id>`, `bd close <id>` for workflow

**`.lean-ctx/`:**
- Purpose: MCP lean-ctx cache (context compression, symbol indices)
- Generated: Auto-built on first lean-ctx use
- Committed: No (.gitignore excludes)
- Cleared: Safe to delete; will regenerate on next tool use

**`leafblower.Rcheck/`:**
- Purpose: R CMD check output (validation, NOT a build artifact)
- Generated: `R CMD INSTALL --preclean .` → writes to .Rcheck
- Committed: No (.gitignore excludes)
- Cleaned: `rm -rf *.Rcheck` between builds

**`python/.venv/`:**
- Purpose: Python virtual environment (uv-managed)
- Generated: `uv sync` from `pyproject.toml`
- Committed: No (.gitignore excludes)
- Regenerate: `cd python && uv sync` (not installed via pip; uses lockfile)

---

*Structure analysis: 2026-08-15*
