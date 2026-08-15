<!-- refreshed: 2026-08-15 -->
# Architecture

**Analysis Date:** 2026-08-15

## System Overview

```text
┌──────────────────────────────────────────────────────────────────────────┐
│                        Public C API Layer                                 │
│                    `src/leafblower.h` (C99-clean)                        │
│  rk_calibrate() → rk_params_t + rk_result_t (C-compatible structs)       │
└─────────────────────────────┬──────────────────────────────────────────┘
                              │
         ┌────────────────────┼─────────────────────┐
         │                    │                     │
         ▼                    ▼                     ▼
┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐
│   R Bridge      │  │  Python Bridge  │  │  C++ Internal   │
│  `r_bridge.cpp` │  │ `_bindings.cpp` │  │   Solvers       │
│  (Rcpp/S3)      │  │  (pybind11)     │  │ `*.cpp`/`*.hpp` │
└────────┬────────┘  └────────┬────────┘  └────────┬────────┘
         │                    │                    │
         └────────────────────┼────────────────────┘
                              │
                              ▼
        ┌─────────────────────────────────────────┐
        │  Algorithm Dispatch Layer               │
        │  `c_api.cpp` — route algorithm enum     │
        │  to concrete solver implementations     │
        └────────────┬────────────────────────────┘
                     │
        ┌────────────┴──────────────────────────────────────────┐
        │                                                       │
        ▼                                                       ▼
┌──────────────────────────────────┐        ┌──────────────────────────────┐
│   Iterative Solvers (BCD/IPF)    │        │   Batch/Direct Solvers       │
│  • ORIS `oris.cpp`               │        │  • GREG `greg.cpp`           │
│  • ORIS-Soft `oris.cpp` (ADMM)   │        │  • Chebyshev `chebyshev.cpp` │
│  • Raking `raking.cpp`           │        │  • Newton-KL                 │
│  • Sinkhorn `sinkhorn.cpp`       │        │    `newton_calib.cpp`        │
│  • Greenkhorn `greenkhorn.cpp`   │        │  • Logit `logit_calib.cpp`   │
│                                  │        │                              │
└────────────┬─────────────────────┘        └──────────────┬───────────────┘
             │                                             │
             └────────────────┬────────────────────────────┘
                              │
                              ▼
        ┌─────────────────────────────────────────┐
        │   Shared Infrastructure                 │
        │  • `types.hpp` — structs, enums         │
        │  • `calib_dispatch.hpp` — convergence   │
        │  • `calib_validate.hpp` — bounds check  │
        │  • `calib_linalg.cpp` — matrix ops     │
        │  • `cell_table.cpp` — cross-tab build   │
        │  • `design_effect.cpp` — Kish/H&V      │
        │  • `validation.cpp` — input validation  │
        │  • `lbw_math.hpp` — math kernels       │
        │  • `lbw_config.h` — build config       │
        └─────────────────────────────────────────┘
```

## Component Responsibilities

| Component | Responsibility | File |
|-----------|----------------|------|
| **C API** | Public entry point, ABI contract, parameter/result marshaling | `src/c_api.cpp`, `src/leafblower.h` |
| **C_rk_calibrate** | Main R entry point, S3 wrapping, list-to-struct conversion | `src/r_bridge.cpp` |
| **Python calibrate()** | Main Python entry point, NumPy array handling, GIL management | `python/leafblower/_bindings.cpp` |
| **Algorithm selector** | Route algorithm enum (AUTO, ORIS, RAKING, etc.) to solver | `src/c_api.cpp` (rk_calibrate) |
| **ORIS solver** | Over-relaxed iterative scaling (box-constrained Sinkhorn BCD) | `src/oris.cpp`, `src/oris.hpp`, `src/oris_internal.hpp`, `src/oris_finalize.cpp`, `src/oris_trajectory.cpp` |
| **Raking solver** | Raking (IPF, proportional margins adjustment) with acceleration | `src/raking.cpp`, `src/raking.hpp` |
| **Sinkhorn solver** | Entropy-minimizing IPF (KL Bregman Dykstra) | `src/sinkhorn.cpp`, `src/sinkhorn.hpp` |
| **GREG solver** | Generalized regression estimator (Newton QP, chi-square) | `src/greg.cpp`, `src/greg.hpp` |
| **Chebyshev solver** | Interior-point LP-based solver (minimax norm) | `src/chebyshev.cpp`, `src/chebyshev.hpp` |
| **Greenkhorn solver** | Greedy coordinate-descent IPF (autumn::harvest style) | `src/greenkhorn.cpp`, `src/greenkhorn.hpp` |
| **Logit solver** | Deville-Sarndal (1992) logit Newton calibration | `src/logit_calib.cpp`, `src/logit_calib.hpp` |
| **Newton-KL solver** | Newton-KL smooth dual (zero-compression regime, TSVD damping) | `src/newton_calib.cpp`, `src/newton_calib.hpp` |
| **Calibration dispatch** | Convergence metrics/rules, error checking, metric selection | `src/calib_dispatch.hpp` |
| **Cell table builder** | Cross-tabulation (cell = unique margin combination) | `src/cell_table.cpp`, `src/cell_table.hpp` |
| **Validation** | Input range checks (group IDs, category counts, bounds) | `src/validation.cpp`, `src/validation.hpp` |
| **Linear algebra** | Matrix solve, Cholesky, TSVD, eigendecomposition | `src/calib_linalg.cpp`, `src/calib_linalg.hpp` |
| **Design effect** | Kish (1965) and Henry-Valliant (2015) design effect | `src/design_effect.cpp`, `src/design_effect.hpp` |
| **Math kernels** | Numeric stability helpers, log-sum-exp, safe arithmetic | `src/lbw_math.hpp` |

## Pattern Overview

**Overall:** Pluggable solver architecture with shared convergence and validation infrastructure.

**Key Characteristics:**
- **Single solver entry** — `rk_calibrate(n, K, weights, group_ids, cat_counts, targets, params, result)` dispatches to 9 distinct algorithms via `rk_algorithm_t` enum
- **Shared result packing** — `CalibResult` base struct (status, iterations, metrics) embedded in every solver-specific result type; `c_api.cpp` and `r_bridge.cpp` unpack solver results into C-compatible `rk_result_t`
- **Algorithm-agnostic convergence** — `calib_dispatch.hpp::apply_rule()` and `select_metric()` reduce repetitive convergence switching across ORIS, raking, etc.
- **Cross-language support** — same C core wrapped by R (Rcpp + roxygen2) and Python (pybind11 + scikit-build), parity validated via co-located test suites
- **Bounds enforcement** — two modes: `RK_BOUNDS_CELL` (cell-aggregate, default) and `RK_BOUNDS_UNIT` (per-observation water-filling via ORIS)

## Layers

**C API Layer:**
- Purpose: Language-agnostic FFI contract, C99-clean public header
- Location: `src/leafblower.h`, `src/c_api.cpp`
- Contains: Enum definitions (`rk_algorithm_t`, `rk_bounds_mode_t`), struct definitions (`rk_params_t`, `rk_result_t`), parameter initialization, result initialization, main `rk_calibrate()` entry point, ABI tripwires (static_assert on struct sizes)
- Depends on: Nothing — header-only C dependencies
- Used by: R bridge, Python bindings, any external C caller

**R Bridge Layer:**
- Purpose: Translate R S3 objects → C structs, manage RAII cleanup (deferred throw on bad scalar args), wrap results as R list
- Location: `src/r_bridge.cpp`
- Contains: `C_rk_calibrate()` native routine (exported via `NAMESPACE`), scalar extraction helpers, struct marshaling, roxygen2 wrapper `harvest()`
- Depends on: C API, R internals (`<R.h>`, `<Rinternals.h>`), `types.hpp` (enum names for dispatch)
- Used by: R code via `.Call()`, roxygen2 wrapping

**Python Binding Layer:**
- Purpose: NumPy array marshaling, GIL release during solve, callback trampoline for logging
- Location: `python/leafblower/_bindings.cpp`
- Contains: `calibrate()` pybind11 module function, array size validation, GIL re-acquisition for py_log_trampoline
- Depends on: C API, pybind11 headers
- Used by: `python/leafblower/__init__.py` (imports the compiled `_leafblower` module)

**Algorithm Dispatch Layer:**
- Purpose: Route `params.algorithm` enum to concrete solver, handle AUTO routing heuristics, pack solver results into rk_result_t
- Location: `src/c_api.cpp` (main rk_calibrate function, lines ~170-550)
- Contains: Algorithm name lookup table (kAlgNames), AUTO heuristic (M_cell estimation → ORIS if M_cell > threshold, else raking), solver-specific result packing (pack_solver_result, pack_oris_result_c, pack_newton_result_c)
- Depends on: All solver headers, calib_dispatch, cell_table, validation
- Used by: C API entry point

**Solver Layer (9 algorithms):**
- Purpose: Distinct mathematical formulation, each with own convergence logic and result struct
- Locations: `src/oris.cpp`, `src/raking.cpp`, `src/sinkhorn.cpp`, `src/greg.cpp`, `src/chebyshev.cpp`, `src/greenkhorn.cpp`, `src/logit_calib.cpp`, `src/newton_calib.cpp`, plus ORIS internals (`oris_finalize.cpp`, `oris_trajectory.cpp`, `oris_internal.hpp`)
- Pattern: Each solver `<name>.cpp` defines struct `<Name>Result : public CalibResult { ... solver-specific fields ... }` and function `<Name>Result <name>_solve(CalibState&)`
- Shared utilities: `calib_dispatch::apply_rule()`, `calib_dispatch::select_metric()`, `calib_validate` bounds checking, `calib_linalg` for matrix ops
- Depends on: types.hpp, calib_dispatch, calib_validate, calib_linalg, lbw_math.hpp, cell_table.hpp
- Used by: Algorithm dispatch layer

**Convergence Helpers:**
- Purpose: Eliminate triplicated convergence switches across solvers
- Location: `src/calib_dispatch.hpp` (inline functions)
- Contains: `select_metric(CalibMetric, max_err, mean_err, ...)` — pick active metric; `apply_rule(CalibRule, curr, prev, tol)` — test convergence criterion (THRESHOLD/IMPROVEMENT/PLATEAU)
- Invariant: all rules require finite curr and tol > 0; IMPROVEMENT skips first check (prev=inf); PLATEAU skips when prev non-finite
- Used by: All iterative solvers (ORIS, raking, sinkhorn, greenkhorn, logit, chebyshev)

**Cell Table Builder:**
- Purpose: Build cross-tabulation (every unique combination of margin categories = one cell); enables matrix-free IPF and compression diagnostics
- Location: `src/cell_table.cpp`, `src/cell_table.hpp`
- Contains: `build_cell_table()` — returns CellTable with M_cell count, cell_of mapping (obs → cell), n_per_cell counts, g_per_cell (margin category per cell); `estimate_M_cell()` — O(n) fast path for AUTO routing; `build_cells_per_cat()` — inverted index (cell → obs per margin category)
- Depends on: types.hpp, leafblower.h (K_MAX=64 limit)
- Used by: c_api (AUTO routing, cell count diagnostic), all solvers for matrix-free computation

**Validation Layer:**
- Purpose: Input sanitization before solver dispatch
- Location: `src/validation.cpp`, `src/validation.hpp`
- Contains: `validate_inputs()` — checks group_ids range ([-1, cat_counts[k])), cat_counts > 0, weights > 0, dimension consistency, K ≤ K_MAX
- Depends on: leafblower.h (return codes)
- Used by: c_api, validation tests

**Math/Numerics:**
- Purpose: Cancellation-free arithmetic, stable log-sum-exp, adaptive finite-difference, condition number estimation
- Location: `src/lbw_math.hpp` (inline), `src/calib_linalg.cpp`/.hpp
- Contains: `logsumexp()`, `norm_weighted()`, TSVD truncation, Cholesky solve, regularized solve with ridge penalty, spectral condition number
- Depends on: Eigen (externally linked), LAPACK/BLAS (system)
- Used by: All solvers for matrix ops and stable computation

## Data Flow

### Primary Solver Path (e.g., ORIS)

1. **Entry** `rk_calibrate(n, K, weights, group_ids, cat_counts, targets, params, result)` (`src/c_api.cpp:~190`)
2. **Validation** → `validation::validate_inputs()` checks array bounds, K ≤ 64 (`src/validation.cpp`)
3. **Algorithm selection** → if `params.algorithm == RK_ALG_AUTO`, estimate M_cell and pick ORIS or raking (`src/c_api.cpp:~220`)
4. **Cell table build** → `build_cell_table(n, K, group_ids, cat_counts, weights, cell_table)` maps each obs to a cell (`src/cell_table.cpp`)
5. **State marshaling** → pack `group_ids`, `targets`, convergence config into `CalibState` struct (`src/c_api.cpp:~250`)
6. **Solver invoke** → `oris_solve(state)` runs block-coordinate descent on Sinkhorn dual (`src/oris.cpp`)
   - Inner loop: cycle through each margin, update IPF scaling factors
   - Convergence check: `calib_dispatch::apply_rule(params.rule, curr_metric, prev_metric, params.pct_tol)` every `kErrCheckInterval` iters
   - On convergence: cap iterations, store best_weights snapshot
7. **Result packing** → `pack_oris_result_c(result, oris_result)` copies ORISResult.base → rk_result_t and ORIS diagnostics (`src/c_api.cpp:~95`)
8. **Return** → status code (RK_OK / RK_ERR_NOCONV / ...) to caller

### Algorithm Selection (AUTO Routing)

1. `estimate_M_cell()` counts unique cells, capped at n
2. If M_cell > threshold (heuristic: ~n/10 or algorithm-specific): route to ORIS (sparse problem, compression saves iterations)
3. Else: route to raking (dense problem, IPF is competitive)
4. Diagnostic: `result.homotopy_levels_used` = 1 if AUTO selected a single-pass algorithm

### Convergence Check Flow

```
Every kErrCheckInterval (default 10 iters):
  1. Compute calibration metrics: max_err, mean_err, kl, chi2, grake_norm, l1_weight_change
  2. select_metric(params.metric) → extract the active metric
  3. apply_rule(params.rule, curr, prev_metric, tol) → test convergence:
     - THRESHOLD: curr < tol
     - IMPROVEMENT: |curr - prev| / prev < tol (skip if prev ≤ 1e-15 or non-finite)
     - PLATEAU: curr >= prev * (1 - tol) [metric did NOT improve by ≥tol fraction]
  4. If converged (via selected rule) and params.stop_when=ANY: exit
  5. If multiple rules active (stop_when=ALL): test all, exit when all fire
  6. If not converged: update prev_metric, continue looping
```

### Python/R Bridge Data Flow

**R:**
- `harvest(..., algorithm="oris")` (R wrapper, `R/harvest.R`)
- ↓ roxygen2 → `.Call("C_rk_calibrate", ...)`
- ↓ `src/r_bridge.cpp::C_rk_calibrate()` extracts scalars, builds rk_params_t
- ↓ `rk_calibrate()` (C API, solves)
- ↓ `src/r_bridge.cpp` packs rk_result_t → R list
- ← Returns S3 `harvest_result` list

**Python:**
- `leafblower.calibrate(n, K, weights_np, group_ids_list, ...)` (Python wrapper)
- ↓ `python/leafblower/_bindings.cpp::calibrate()` (pybind11)
- ↓ Validates NumPy array sizes, builds std::vector<const int32_t*> pointers
- ↓ Releases GIL (py::gil_scoped_release) during rk_calibrate
- ↓ `rk_calibrate()` (C API, solves); py_log_trampoline re-acquires GIL on callback
- ↓ Packs rk_result_t → Python dict
- ← Returns (weights_np, result_dict) tuple

**State Management:**
- Weights: updated in-place during solve (C array modified, reflected back to caller)
- Result: struct populated by solver, copied to caller's output struct (no aliasing)
- Cell table: local (destroyed on return); not exposed to R/Python

## Key Abstractions

**CalibResult (base struct):**
- Purpose: Shared result fields across all solvers (status, iterations, metrics, convergence metadata)
- Examples: `ORISResult { CalibResult base; ...; }`, `RakingResult { CalibResult base; ...; }`
- Pattern: Each solver's result struct inherits (embeds) CalibResult and extends it with solver-specific diagnostics (e.g., ORIS adds sor_min_omega, alm_capacity_mu_final)
- File: `src/types.hpp`

**CalibState (solver configuration):**
- Purpose: Unified input struct for all solvers (weights, group_ids, targets, convergence config, etc.)
- Contains: n, K, pointers to arrays, bounds, tolerances, convergence criterion (metric/rule), homotopy/scheduler/eta config, ALM config, ridge penalty, etc.
- Pattern: Built once by c_api, passed by reference to solver; solver reads (does not modify except weights)
- File: `src/types.hpp`

**CellTable (compression mapping):**
- Purpose: Fast cross-tabulation; enables IPF without dense matrix allocation
- Contains: M_cell (unique cell count), cell_of (n-elem obs→cell map), n_per_cell (M_cell-elem counts), g_per_cell (K×M_cell category grid)
- Pattern: Built by `build_cell_table()`; passed to solver; solver loops over cells, not all n obs
- File: `src/cell_table.hpp`

**rk_params_t (C-exposed config):**
- Purpose: Language-agnostic parameter struct; ABI-frozen (264 bytes, checked via static_assert)
- Contains: min_weight, max_weight, tolerances, algorithm enum, bounds_mode, overlay knobs (homotopy, scheduler, eta, SOR, ALM), convergence (metric/rule/thresholds)
- Pattern: Caller fills via defaults (rk_params_init) and selective overrides; solver reads (no modification)
- File: `src/leafblower.h`

**rk_result_t (C-exposed output):**
- Purpose: Language-agnostic result struct; ABI-frozen (536 bytes)
- Contains: status, iterations, max_error, mean_error, kl, chi2, algorithm_used, convergence metadata, ORIS-specific diagnostics (sor_*, alm_*, homotopy_*), design-effect diagnostics
- Pattern: Zero-initialized by caller; solver writes via pack functions; caller reads (read-only after return)
- File: `src/leafblower.h`

## Entry Points

**C Entry Point:**
- Location: `src/c_api.cpp::rk_calibrate()`
- Triggers: Called by R bridge (via .Call) or Python bindings
- Responsibilities: Validate inputs, route algorithm, build cell table, marshal CalibState, invoke solver, pack result

**R Entry Point:**
- Location: `src/r_bridge.cpp::C_rk_calibrate()`
- Triggers: Invoked via `.Call("C_rk_calibrate", ...)` from R roxygen2 wrapper `harvest()`
- Responsibilities: Extract scalars from SEXP, build rk_params_t, call rk_calibrate(), pack rk_result_t → R list

**Python Entry Point:**
- Location: `python/leafblower/_bindings.cpp::calibrate()` (pybind11 module)
- Triggers: Called from `python/leafblower/__init__.py` via `_leafblower.calibrate(...)`
- Responsibilities: Validate NumPy arrays, build std::vector pointers, release GIL, call rk_calibrate(), pack rk_result_t → Python dict

## Architectural Constraints

- **Threading:** Single-threaded per-solve; inner loops use OpenMP `#pragma omp simd` for vectorization (optional, not required); GIL is released during solve in Python
- **Global state:** None — all state is passed via CalibState struct; no module-level singletons
- **Circular imports:** None detected; dependency DAG is acyclic (solvers → calib_dispatch → types)
- **ABI stability:** rk_params_t and rk_result_t sizes frozen (264B and 536B); padding fields reserved for future extensions (e.g., sor_omega_mode_id repurposed _pad slot)
- **Cell limit:** K ≤ 64 (prevents unbounded memory; enforced by build_cell_table and estimate_M_cell)
- **Weight bounds enforcement:** Two distinct modes (cell vs. unit); cell mode does NOT clamp per-observation weights (diagnostic only); unit mode uses intra-cell water-filling to enforce per-obs bounds
- **No cancellation policy:** All variance/covariance computed without subtraction of like-magnitude terms (e.g., `p*(1-p)` not `p - p*p`)

## Anti-Patterns

### Algorithm Repeated in Multiple Solvers

**What happens:** Convergence checking logic (THRESHOLD/IMPROVEMENT/PLATEAU decision trees) was copy-pasted across ORIS, raking, sinkhorn, etc., leading to maintenance burden and divergence bugs.

**Why it's wrong:** Any change to convergence semantics must be applied to 6+ solvers; bugs in one solver may not trigger in others; future solver additions will copy stale logic.

**Do this instead:** Use `calib_dispatch.hpp::apply_rule()` and `select_metric()` inline helpers. These are generic and reduce the switch blocks to a single call site per solver.

### Solver Results With No Shared Base

**What happens:** Early solver result structs (ORISResult, RakingResult) had different memory layouts; packing them into rk_result_t required separate branches per solver.

**Why it's wrong:** Result marshaling scales O(# solvers) and is error-prone when adding new fields; branch coverage is hard to test; new solvers risk omitting fields.

**Do this instead:** Every solver result embeds `CalibResult base` as the first field; `pack_solver_result()` template unpacks the base uniformly; solver-specific extras (ORIS sor_*, ALM fields) are handled via SFINAE (if constexpr has_n_bounds_c<R>::value).

### Hard-Coded Algorithm Limits in Multiple Places

**What happens:** Algorithm AUTO heuristic (M_cell threshold to choose ORIS vs. raking) was hand-tuned in c_api.cpp and duplicated in R wrapper logic.

**Why it's wrong:** Tuning the threshold requires changes in two places; desyncs silently lead to divergent routing decisions.

**Do this instead:** Define the threshold in one place (cell_table.cpp::estimate_M_cell or c_api.cpp); expose it via a parameter or comment; document the heuristic rationale (e.g., "ORIS preferred for M_cell > n/10 to leverage cell compression").

## Error Handling

**Strategy:** Return codes (RK_OK, RK_ERR_NOCONV, RK_ERR_INFEAS, RK_ERR_BADARG, RK_ERR_BUDGET, RK_ERR_STALL) propagated from solver → c_api → bridge → caller. No exceptions thrown across language boundaries.

**Patterns:**
- **Validation errors** (RK_ERR_BADARG): Input range check fails → solver returns immediately; message set in rk_result_t.message
- **Convergence failure** (RK_ERR_NOCONV): Solver hits max_iterations without convergence; best_weights stored; caller can inspect best_iter and best_error
- **Infeasibility** (RK_ERR_INFEAS): Cell has positive target but zero input weight; no solution exists
- **Budget exhausted** (RK_ERR_BUDGET): Loss decreasing at max_iterations; suggest increasing max_iterations
- **Stall** (RK_ERR_STALL): Plateau detected (metric ≥ prev * (1 - tol)); at constrained optimum; weights are valid
- **Message field:** 256-char null-terminated string; populated on error (c_api, r_bridge only); ORISResult has no message field (synthesized at packing stage if empty)

## Cross-Cutting Concerns

**Logging:** Optional callback (`log_fn`/`log_ctx` in CalibState); R calls REprintf(), Python uses py_log_trampoline to re-acquire GIL and invoke Python callable. Default: silent.

**Validation:** Centralized in `validation.cpp::validate_inputs()` (called once per rk_calibrate); checks dimensions, category codes, K ≤ 64, bounds feasibility. Solver-specific checks (e.g., TSVD ratio > 0) done at solver entry.

**Authentication/Authorization:** Not applicable — C library with no I/O or network.

**Design tradeoffs:** 
- Cell table built eagerly (not lazy) to fail fast on validation and enable AUTO routing heuristics
- Solver results are not modified after return (caller owns best_weights snapshot); enables safe sharing with async threads (immutable after read)
- GIL released during solve in Python to allow concurrent calls from other threads; r_bridge has no threading constraint (R is single-threaded at C level)

---

*Architecture analysis: 2026-08-15*
