# Codebase Structure

**Analysis Date:** 2026-08-14

## Directory Layout

```
leafblower/
├── src/                 # C++17 numerical core + C ABI + R bridge (the whole engine)
├── R/                   # R user layer (8 files, marshalling + validation only)
├── man/                 # roxygen2-generated .Rd docs
├── python/              # Python package: pybind11 bindings + scikit-build
│   ├── CMakeLists.txt   # explicit CORE_SOURCES list (does NOT glob src/)
│   ├── pyproject.toml   # scikit-build-core; version must be synced with DESCRIPTION
│   └── leafblower/      # _bindings.cpp, _harvest.py, _design_effect.py, test_*.py
├── tests/
│   ├── testthat/        # 117 R test files + fixtures/ + .rds references
│   ├── parity/          # R driver scripts for R↔Python parity runs
│   ├── cpp/             # (empty)
│   └── testthat.R       # testthat entry
├── benchmarks/          # stepstone / algo-selection benchmark scripts + result .rds/.csv
├── docs/                # methods notes, investigations, iEPPA specs, superpowers specs
├── data-raw/            # raw data prep scripts
├── tools/               # maintenance R scripts (e.g. check_research_isolation.R)
├── tasks/               # PRD + lessons
├── .beads/              # beads issue tracker (issues.jsonl, plans, hooks)
├── configure            # autoconf-style feature probes → src/Makevars, src/lbw_config.h
├── DESCRIPTION          # R package metadata + version
├── NAMESPACE            # useDynLib(.registration=TRUE) + 7 exports
└── CLAUDE.md            # authoritative project invariants
```

## Directory Purposes

**`src/`:**
- Purpose: the entire numerical implementation. No numerics live anywhere else.
- Contains: 18 `.cpp`, 17 `.hpp`, the public C header, build inputs.
- Key files: `leafblower.h` (C ABI), `c_api.cpp` (ABI impl + AUTO routing),
  `r_bridge.cpp` (SEXP layer), `calib_dispatch.hpp` (shared helpers),
  `cell_table.hpp/.cpp`, `types.hpp`, `Makevars.in`.
- Note: build artifacts (`*.o`, `leafblower.so`, generated `Makevars`, `lbw_config.h`) live
  here untracked after a build.

**`R/`:**
- Purpose: R-facing API, argument validation, factor encoding, result shaping.
- Key files: `harvest.R` (1269 lines — the main entry plus 10 private parsers),
  `design_effect.R`, `diagnose_weights.R`, `anesrake.R` (compat shim),
  `current_miss.R`, `weighted_pct.R`, `na_bin.R`, `zzz.R`.

**`python/leafblower/`:**
- Purpose: Python-facing API + the pybind11 module source + the test suite (tests are
  co-located with the package, not in a separate `tests/` tree).
- Key files: `_bindings.cpp` (322 lines), `_harvest.py` (874), `_design_effect.py` (116),
  `__init__.py` (re-exports 4 symbols).

**`tests/testthat/`:**
- Purpose: the R gate. One file per behaviour area, plus `.rds` reference fixtures
  (`task1_ref.rds`, `task2_oris_ref.rds`, `fixtures/`).

**`benchmarks/`:**
- Purpose: stepstone regression gate and algorithm-selection studies. Mixed R and Python
  drivers plus checked-in `.rds`/`.csv`/`.pdf` results. Subdirs are per-ticket
  (`2apm/`, `yh0l/`, `study/`, `results/`).

## Key File Locations

**Entry Points:**
- `R/harvest.R:266`: `harvest()` — R calibration entry.
- `python/leafblower/_harvest.py:237`: `harvest()` — Python calibration entry.
- `src/c_api.cpp:256`: `rk_calibrate()` — C ABI entry.
- `src/r_bridge.cpp:234`: `C_rk_calibrate` — the 39-arg `.Call` entry.
- `python/leafblower/_bindings.cpp:27`: `PYBIND11_MODULE(_leafblower, m)`.

**Configuration:**
- `configure`: probes glibc libmvec/AVX2 and SIMD, substitutes `@MAVX2_FLAG@` etc.
- `src/Makevars.in`: R build flags; `PKG_SOURCES` is decorative (R globs `src/*.cpp`).
- `python/CMakeLists.txt`: hermetic re-probe of the same macros + explicit `CORE_SOURCES`.
- `DESCRIPTION` and `python/pyproject.toml`: versions, bumped manually in lockstep.
- `.coverage-thresholds.json`: the behavioural quality gate command.

**Core Logic:**
- `src/leafblower.h`: types, enums, return codes, ABI size tripwires.
- `src/calib_dispatch.hpp`: metrics, convergence rules, setup, `finalize_weights`.
- `src/cell_table.hpp` / `src/cell_table.cpp`: CellTable.
- `src/types.hpp`: `CalibState`, `CalibResult`, overlay config structs.
- Solvers: `oris.cpp` (+`oris_finalize.cpp`, `oris_trajectory.cpp`, `oris_internal.hpp`),
  `raking.cpp`, `sinkhorn.cpp`, `greenkhorn.cpp`, `chebyshev.cpp`, `greg.cpp`,
  `logit_calib.cpp` (+`logit.cpp`/`logit.hpp` link-function helpers), `newton_calib.cpp`.
- Support: `calib_linalg.cpp` (LAPACK wrappers), `lbw_math.hpp` (bulk exp/log, AVX2),
  `sraa.hpp` (Anderson acceleration), `validation.cpp`, `calib_validate.cpp`,
  `design_effect.cpp`.

**Testing:**
- `tests/testthat/test-*.R`: R suite; `tests/testthat.R` is the runner hook.
- `python/leafblower/test_*.py`: pytest suite (parity + Python-layer).
- `tests/parity/run_*.R`: R-side drivers producing references for parity comparison.
- `tests/test_parity_weights.py`, `conftest.py` (root and `python/`).

## Naming Conventions

**Files:**
- Solvers: `<algorithm>.cpp` + `<algorithm>.hpp`, lowercase snake (`newton_calib.cpp`).
- Cold TU splits of one solver: `<solver>_<role>.cpp` + private `<solver>_internal.hpp`
  (`oris_finalize.cpp`, `oris_trajectory.cpp`, `oris_internal.hpp`).
- Shared headers: `calib_*.hpp` for the calibration-wide layer, `lbw_*.h[pp]` for
  math/config primitives.
- R: one file per exported function, named after it (`harvest.R`, `design_effect.R`).
- Python: private modules are underscore-prefixed (`_harvest.py`, `_bindings.cpp`).
- R tests: `test-<kebab-case-area>.R`. Ticket-scoped regressions embed the ticket id
  (`test-cr-d16-nbounds.R`).
- Python tests: `test_<snake_case>.py`, ticket-scoped as `test_cr_<id>_<area>.py`.

**Symbols:**
- C ABI: `rk_` prefix for functions and types, `RK_` for enum values and macros.
- C++ core: everything in `namespace lbw`; types `PascalCase` (`CalibState`, `ORISResult`),
  functions `snake_case` (`oris_solve`, `finalize_weights_buf`), constants `kPascalCase`
  (`kMinSafeTotalWeight`, `kUnboundedSentinel`, `K_MAX`).
- R `.Call` symbols: `C_<c_function_name>` (`C_rk_calibrate`), registered in
  `src/r_bridge.cpp:182-188`.
- R private helpers: bare snake_case, not exported (`parse_target`, `map_method`); dot-prefixed
  for internal-only (`.encode_na_bin_mask`).

## Where to Add New Code

**New solver (8 required steps — see `CLAUDE.md`):**
- Implementation: `src/<name>.cpp` + `src/<name>.hpp` in `namespace lbw`, result struct
  embedding `CalibResult base`.
- REQUIRED: add the `.cpp` to `CORE_SOURCES` in `python/CMakeLists.txt` (R globs, Python does not).
- Enum value: `rk_algorithm_t` in `src/leafblower.h` (never slot 2 or 7) + `kAlgNames` in
  `src/c_api.cpp:30`.
- Dispatch: a branch in BOTH `src/c_api.cpp:414+` (enum) and `src/r_bridge.cpp:620+` (string),
  plus `map_method()` in `R/harvest.R:975` and the Python method map.
- Shared helpers: `src/calib_dispatch.hpp`, never in the solver file.
- Tests: `tests/testthat/test-<name>.R` with an `.rds` fixture, and
  `python/leafblower/test_solver_parity.py` coverage.
- Benchmark fixture for the stepstone regression gate under `benchmarks/`.

**New shared C++ helper:**
- General calibration logic → `src/calib_dispatch.hpp` (header-only, `inline`).
- CellTable-specific → `src/cell_table.hpp`.
- Linear algebra → `src/calib_linalg.hpp`; elementwise math kernels → `src/lbw_math.hpp`.

**New R user-facing function:**
- `R/<function_name>.R`, roxygen2 block, add to `NAMESPACE` export list, regenerate `man/`.

**New Python user-facing function:**
- `python/leafblower/_<name>.py`, re-export in `python/leafblower/__init__.py`.

**New C ABI field:**
- Add to `rk_params_t` / `rk_result_t` in `src/leafblower.h`, update the
  `EXPECTED_RK_*_BYTES` tripwire and its layout comment, then thread it through
  `rk_params_init`, `CalibState`, both dispatch tables, both bridges, and both language layers.

**Splitting a hot TU:**
- Only COLD (once-per-solve) code may move out — there is no LTO. Follow the `oris.cpp`
  pattern: hot loop stays, finalization/trajectory move, private declarations in
  `<solver>_internal.hpp`.

## Special Directories

**`man/`:** roxygen2-generated. Committed. Do not hand-edit; a stray regenerated
`man/dot-*.Rd` is a known artifact to clean before committing.

**`.beads/`:** issue tracker state (`issues.jsonl`, `plans/`, `hooks/`). Committed, but a hook
re-stages `issues.jsonl` — always commit with an explicit pathspec, never `git add -A`.

**`benchmarks/`:** results (`.rds`, `.csv`, `.pdf`) are committed alongside the drivers.

**`*.Rcheck/`, `node_modules/`, `python/.venv/`, `__pycache__/`, `graphify-out/`,
`src/*.o`, `src/leafblower.so`, `src/Makevars`, `src/lbw_config.h`:** generated, not committed.

**Untracked strays at repo root:** `cell_table_92c4f45.cpp/.hpp`, `ieppa_92c4f45.cpp`,
`patch_raking.py`, `patch_wolfe.py`, `test_output.log`, `leafblower_0.1.0.tar.gz`,
`code-review-findings.md`, `REVIEW_FINDINGS.md` — snapshots and scratch, not part of the build.

---

*Structure analysis: 2026-08-14*
