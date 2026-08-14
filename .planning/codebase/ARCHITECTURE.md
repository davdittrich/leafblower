<!-- refreshed: 2026-08-14 -->
# Architecture

**Analysis Date:** 2026-08-14

## System Overview

```text
┌─────────────────────────────────────────────────────────────┐
│                    User-facing language layers               │
├──────────────────────────────┬──────────────────────────────┤
│   R package (S3-free lists)  │   Python package (pandas)    │
│   `R/harvest.R`              │   `python/leafblower/`       │
│   `R/design_effect.R`        │     `_harvest.py`            │
│   `R/diagnose_weights.R`     │     `_design_effect.py`      │
└──────────────┬───────────────┴───────────────┬──────────────┘
               │ .Call(C_rk_calibrate, …)      │ _leafblower.calibrate(…)
               ▼                                ▼
┌──────────────────────────────┬──────────────────────────────┐
│  R bridge (SEXP marshalling) │  pybind11 bridge             │
│  `src/r_bridge.cpp`          │  `python/leafblower/`        │
│  39-arg .Call, own dispatch  │    `_bindings.cpp`           │
│  on method STRING            │  calls the C API             │
└──────────────┬───────────────┴───────────────┬──────────────┘
               │                                │
               │            ┌───────────────────┘
               ▼            ▼
┌─────────────────────────────────────────────────────────────┐
│  C ABI  `src/leafblower.h` / `src/c_api.cpp`                │
│  rk_calibrate() · rk_design_effect() · rk_params_init()     │
│  AUTO routing · rk_params_t → lbw::CalibState · result pack │
└─────────────────────────────┬───────────────────────────────┘
                              ▼
┌─────────────────────────────────────────────────────────────┐
│  Shared solver-support layer (header-only, `lbw` namespace) │
│  `src/calib_dispatch.hpp`  metrics, convergence rules,      │
│                            setup, finalize_weights          │
│  `src/cell_table.hpp`      CellTable build / M_cell estimate│
│  `src/types.hpp`           CalibState, CalibResult, configs │
│  `src/calib_linalg.hpp` `src/lbw_math.hpp` `src/sraa.hpp`   │
└─────────────────────────────┬───────────────────────────────┘
                              ▼
┌─────────────────────────────────────────────────────────────┐
│  Solver translation units (one per algorithm)               │
│  oris · raking · sinkhorn · greenkhorn · chebyshev · greg   │
│  logit_calib · newton_calib                                 │
└─────────────────────────────┬───────────────────────────────┘
                              ▼
┌─────────────────────────────────────────────────────────────┐
│  Weight finalization — the single exit funnel               │
│  `lbw::finalize_weights[_buf]` (`src/calib_dispatch.hpp`)   │
│  normalize Σw=n → bounds_mode dispatch                      │
└─────────────────────────────────────────────────────────────┘
```

## Component Responsibilities

| Component | Responsibility | File |
|-----------|----------------|------|
| `harvest()` (R) | User entry: validate, encode margins to 0-based int codes, marshal 39 args, post-process result | `R/harvest.R:266` |
| `harvest()` (Python) | Same contract over pandas; parses convergence/SOR dicts | `python/leafblower/_harvest.py:237` |
| R bridge | SEXP ⇄ C++ marshalling, its OWN string-keyed solver dispatch, RAII-safe `Rf_error` unwinding | `src/r_bridge.cpp:234` |
| pybind11 bridge | Python ⇄ C API marshalling; exports `calibrate`, `_design_effect` | `python/leafblower/_bindings.cpp:27` |
| C ABI | Stable `extern "C"` surface + AUTO routing + result packing | `src/leafblower.h`, `src/c_api.cpp:256` |
| Cell table | Compress n observations into M_cell distinct margin-crossing cells | `src/cell_table.hpp:20`, `src/cell_table.cpp` |
| Dispatch helpers | `select_metric`, `apply_rule`, `check_convergence`, `solver_setup_ct`, `finalize_weights` | `src/calib_dispatch.hpp` |
| ORIS solver | Over-Relaxed Iterative Scaling, default/AUTO primary; hot loop | `src/oris.cpp` (+ `oris_finalize.cpp`, `oris_trajectory.cpp`) |
| Other solvers | raking, sinkhorn, greenkhorn, chebyshev, greg, logit, newton_kl | `src/raking.cpp` … `src/newton_calib.cpp` |
| Design effect | Kish + Henry&Valliant deff, independent C entry | `src/design_effect.cpp`, `src/leafblower.h:191` |

## Pattern Overview

**Overall:** Single numerical core, thin polyglot bridges. A C ABI is the waist of the
hourglass; both language layers are marshalling code only — no numerics live in `R/` or
`python/leafblower/*.py`.

**Key Characteristics:**
- Cell-table compression: all solvers work in cell space (`M_cell` ≤ n), then expand to
  observations once at exit (`expand_obs`, `src/calib_dispatch.hpp:294`).
- Solver polymorphism by convention, not virtual dispatch: each solver is a free function
  `lbw::<name>_solve(CalibState&) -> <Name>Result`, and every result struct embeds
  `CalibResult base` (`src/types.hpp:87`). Packing is done by C++ templates
  (`pack_solver_result`, `src/c_api.cpp:57`) plus SFINAE traits for optional fields.
- Header-only shared layer: `calib_dispatch.hpp` is `inline`-everything so hot helpers
  inline into each solver TU (there is no LTO — see Constraints).
- ABI is frozen with `static_assert` tripwires on struct sizes
  (`src/leafblower.h:245`, `:269`).

## Layers

**Language layer (R):**
- Purpose: argument validation, factor/character → integer margin codes, NA-bin encoding,
  auto-collapse of rare categories, result naming and diagnostics.
- Location: `R/` (`harvest.R`, `anesrake.R`, `design_effect.R`, `diagnose_weights.R`,
  `current_miss.R`, `weighted_pct.R`, `na_bin.R`, `zzz.R`).
- Depends on: registered C symbols via `useDynLib(leafblower, .registration = TRUE)` (`NAMESPACE:1`).
- Used by: package users; exported set is 7 functions (`NAMESPACE`).

**Language layer (Python):**
- Purpose: same contract over pandas DataFrames.
- Location: `python/leafblower/_harvest.py`, `_design_effect.py`, `__init__.py`.
- Depends on: `from ._leafblower import calibrate` (`python/leafblower/_harvest.py:8`).

**Bridge layer:**
- `src/r_bridge.cpp` — 1274 lines, registers 6 `.Call` entries (`src/r_bridge.cpp:182-188`);
  `C_rk_calibrate` takes 39 SEXP args. It is EXCLUDED from the Python build.
- `python/leafblower/_bindings.cpp` — pybind11 module `_leafblower`.

**C ABI layer:**
- `src/leafblower.h` — `rk_params_t`, `rk_result_t`, `rk_algorithm_t`, return codes.
- `src/c_api.cpp` — `rk_calibrate` (`:256`), `rk_design_effect` (`:590`), defaults (`:177`).

**Shared solver-support layer (`lbw` namespace, header-only):**
- `src/calib_dispatch.hpp` — canonical home for anything shared by 2+ solvers.
- `src/cell_table.hpp` / `.cpp` — CellTable-specific helpers only.
- `src/types.hpp` — `CalibState` (solver input) and `CalibResult` (solver output base).
- `src/calib_linalg.hpp`, `src/lbw_math.hpp` (bulk `exp`/`log`, AVX2 libmvec path),
  `src/sraa.hpp` (Anderson acceleration), `src/calib_validate.hpp`, `src/validation.hpp`.

**Solver layer:**
- One `.cpp` + `.hpp` pair per algorithm, all in `namespace lbw`. ORIS is split across three
  TUs (`oris.cpp`, `oris_finalize.cpp`, `oris_trajectory.cpp`, shared `oris_internal.hpp`).

## Data Flow

### Primary calibration path (R)

1. `harvest(data, target, …)` validates and encodes each margin column to a 0-based
   `integer` vector (`-1` = NA/OOV) (`R/harvest.R:266`, encoding block `R/harvest.R:554`).
2. `map_method()` normalizes the method string; convergence/SOR dicts parsed
   (`R/harvest.R:975`, `:989`, `:1072`).
3. `.Call(C_rk_calibrate, …)` with 39 arguments (`R/harvest.R:600`).
4. `C_rk_calibrate` builds `lbw::CalibState st`, then dispatches on the method STRING
   (`src/r_bridge.cpp:620` onward) directly into `lbw::<solver>_solve(st)` — it does NOT go
   through `rk_calibrate`.
5. Solver: `solver_setup_ct(st, ct, …)` builds the CellTable, `X_init`, cell bounds and runs
   pre-entry validation (`src/calib_dispatch.hpp:511`); the iteration loop calls
   `compute_cell_metrics` (`:247`) and `check_convergence` (`:204`) each check interval, with
   `BestIterTracker` (`:136`) holding the best iterate.
6. Exit: `expand_obs` (cell → obs, no clamp) then `lbw::finalize_weights` /
   `finalize_weights_buf` (`src/calib_dispatch.hpp:359`).
7. Bridge packs fields via the `pack_solver_result` lambda (`src/r_bridge.cpp:571`) and
   returns a named list; R re-nests convergence fields and maps status → message
   (`R/harvest.R:700-770`).

### Primary calibration path (Python / direct C)

1. `leafblower.harvest(df, targets, …)` (`python/leafblower/_harvest.py:237`).
2. `_leafblower.calibrate(…)` (`python/leafblower/_bindings.cpp:30`).
3. `rk_calibrate` (`src/c_api.cpp:256`): resolves algorithm, validates
   (`validate_inputs`, `:336`), rejects non-finite homotopy factors (`:343`), builds
   `CalibState` (`:352`), dispatches on the `rk_algorithm_t` enum (`:414-548`).
4. Same solver → `finalize_weights` path as above.
5. AUTO-only fallback: on `RK_ERR_NOCONV` / `RK_ERR_BUDGET` the original weights are
   restored and `newton_calibrate` is re-run (`src/c_api.cpp:554-570`).

### AUTO algorithm routing (`src/c_api.cpp:292-330`)

1. `estimate_M_cell(n, K, group_ids, cat_counts)` — O(n) cell-count estimate.
2. `M_cell*10 < n*9` (compression ratio < 0.9) → `RK_ALG_ORIS`.
3. Else `K < 5` → `RK_ALG_RAKING`.
4. Else compute `target_skew = max_T / max(min_T, 1e-12)`; `> 5.0` → `RK_ALG_ORIS` with
   `st.accelerate = true` (SRAA), otherwise `RK_ALG_NEWTON_KL`.

**State Management:**
- `CalibState` (`src/types.hpp:116`) is the single mutable solver input; `st.weights` is a
  caller-owned `double*` mutated in place. All other pointers are borrowed, never owned.
- Results are values (`ORISResult`, `RakingResult`, …) returned by value, moved into the
  bridge.

## Key Abstractions

**`CalibState`:**
- Purpose: everything a solver needs — data pointers, bounds, tolerances, every overlay
  config (homotopy, scheduler, eta schedule, convergence, SOR, ALM).
- File: `src/types.hpp:116`.

**`CalibResult` (`res.base`):**
- Purpose: the field block every solver result carries. VERIFIED: all solver structs embed
  it as a member named `base` (`src/oris.hpp:8`), and every consumer reads `res.base.status`
  etc. (`src/c_api.cpp:60-79`, `src/r_bridge.cpp:571`). Direct `res.status` does not compile
  / is not present — the CLAUDE.md invariant holds.

**`CellTable`:**
- Purpose: compress observations to distinct margin-crossing cells; `cell_of[i]`,
  `n_per_cell[c]`, `g_per_cell[k][c]`.
- File: `src/cell_table.hpp:7`. `K_MAX = 64` (`:37`).

**`BestIterTracker`:**
- Purpose: uniform best-iterate bookkeeping, replacing per-solver ad-hoc variables.
- File: `src/calib_dispatch.hpp:136`.

## Entry Points

**`harvest()` (R):** `R/harvest.R:266` — main calibration API.
**`harvest()` (Python):** `python/leafblower/_harvest.py:237`.
**`rk_calibrate()` (C ABI):** `src/c_api.cpp:256`, declared `src/leafblower.h:226`.
**`C_rk_calibrate` (.Call):** `src/r_bridge.cpp:234`, registered `src/r_bridge.cpp:185`.
**`_leafblower.calibrate` (pybind11):** `python/leafblower/_bindings.cpp:30`.
**`rk_design_effect()`:** `src/c_api.cpp:590`; R side `R/design_effect.R:35`.
**`R_init_leafblower` symbol table:** `src/r_bridge.cpp:182-188` (6 entries).

## Architectural Constraints

- **No LTO.** VERIFIED: `-flto` appears in neither `configure` nor `src/Makevars.in`.
  Consequence: cross-TU calls do not inline. When splitting a TU, move only COLD
  (once-per-solve) code out; the hot per-iteration loop must stay with its callers. `oris` is
  split exactly this way — `oris_solve` stays in `oris.cpp`, finalization and trajectory dump
  moved to `oris_finalize.cpp` / `oris_trajectory.cpp`.
- **No `-O` level set by the package.** VERIFIED: `src/Makevars.in` states it explicitly and
  `PKG_CXXFLAGS = @CXXFLAGS_STD@ @OMP_FLAGS@ @SIMD_FLAGS@ @MAVX2_FLAG@ -I. -DSTRICT_R_HEADERS`
  contains no `-O` (CRAN's `tools:::.check_make_vars` "pflags" check would reject it). The
  Python build DOES set `-O3` (`python/CMakeLists.txt`), so the two builds differ in
  optimization unless the user sets `-O3` in `~/.R/Makevars`.
- **Two independent source lists.** R auto-globs `src/*.cpp` and IGNORES `PKG_SOURCES`
  (`src/Makevars.in` says so verbatim). `python/CMakeLists.txt` uses an explicit
  `CORE_SOURCES` list (17 files; `r_bridge.cpp` excluded). A new `src/*.cpp` MUST be added to
  `CORE_SOURCES` or the pybind11 link fails.
- **Algorithm slot 2 reserved.** VERIFIED: `src/leafblower.h:44` `/* 2 = removed (was
  RK_ALG_LBFGSB) */`, and `src/c_api.cpp:33` maps index 2 to `"(reserved)"`. Slot 7 is a
  further gap. Do not reuse either.
- **ABI frozen by `static_assert`.** `sizeof(rk_params_t) == 264`, `sizeof(rk_result_t) == 536`
  (`src/leafblower.h:245`, `:269`). Adding a field requires updating the tripwire AND every
  consumer.
- **`K_MAX = 64`** margins (`src/cell_table.hpp:37`); `build_cell_table` returns -1 above it.
- **Threading:** single-threaded. `LBW_HAS_OMP` is fixed to 0 in the Python build
  (`python/CMakeLists.txt`) to preserve R↔Python parity; determinism additionally requires
  `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`.
- **Global state:** none. No module-level mutable singletons; `R/zzz.R:5` `.onLoad` is a no-op.
- **`LBW_NO_R`:** the Python build defines it; `types.hpp:152` and `calib_dispatch.hpp:185`
  switch between `REprintf` and `fprintf(stderr)` on it. CRAN forbids R packages writing
  directly to stderr, so any new diagnostic print must be `#ifndef LBW_NO_R`-guarded.

## Anti-Patterns

### Renormalizing weights after the bounds pass

**What happens:** scaling weights to Σw=n after `finalize_weights` has already water-filled.
**Why it's wrong:** it re-pushes clamped observations back above `max_weight`, silently
voiding the `bounds_mode="unit"` guarantee. The R wrapper explicitly documents that it does
NOT normalize for this reason (`R/harvest.R:~745`), and the C contract repeats it
(`src/c_api.cpp:248-255`).
**Do this instead:** the sanctioned order is normalize → bounds, both inside
`finalize_weights_buf` (`src/calib_dispatch.hpp:363-367` scale to n, `:369` bounds_mode
dispatch). VERIFIED: the invariant holds in code. Degenerate `total_w <= kMinSafeTotalWeight`
(1e-100, `:321`) is left unscaled to avoid subnormal overflow.

### Clamping per observation during cell → obs expansion

**What happens:** applying `[min_weight, max_weight]` inside `expand_obs`.
**Why it's wrong:** the cell contract is on the aggregate `X[c] <= U_cell`; per-obs clamping
there distorts marginals and breaks Σw=n (measured 13pp margin drift, Σw 3392 vs n 4000 —
`src/calib_dispatch.hpp:284-293`).
**Do this instead:** expand unclamped, then let `finalize_weights` apply the bounds_mode
contract.

### Adding shared logic to an individual solver file

**What happens:** a helper needed by two solvers grows a private copy in each.
**Why it's wrong:** this repo already paid for it — `oris_finalize.cpp:152-154` records that
the normalize/clamp/water-fill core was duplicated and had to be collapsed back into
`finalize_weights_buf`.
**Do this instead:** put it in `src/calib_dispatch.hpp` (general) or `src/cell_table.hpp`
(CellTable-specific).

### Two divergent dispatch tables

**What happens:** `src/r_bridge.cpp:620+` dispatches on a method STRING into the solvers
directly, while `src/c_api.cpp:414` dispatches on the `rk_algorithm_t` enum. Both then
hand-pack `rk_result_t`/SEXP fields.
**Why it's wrong:** every new solver, every new result field and every routing tweak must be
made twice, and the comment trail (`src/c_api.cpp:458`, `:473` "mirrors r_bridge.cpp:806-810";
`:505` "mirrors r_bridge.cpp:628-657") shows the two have already drifted and been re-synced
repeatedly.
**Do this instead:** when touching either dispatch path, grep the mirrored line ranges in the
other and change both in the same commit. Adding a solver is an 8-step checklist for this
reason (`CLAUDE.md`).

## Error Handling

**Strategy:** integer status codes end-to-end; no exceptions cross the ABI.

**Patterns:**
- Codes: `RK_OK`(0), `RK_ERR_NOCONV`(1), `RK_ERR_INFEAS`(2), `RK_ERR_BADARG`(3),
  `RK_ERR_BUDGET`(4), `RK_ERR_STALL`(5) — `src/leafblower.h:32-38`.
- `c_api.cpp` propagates the solver status verbatim (`src/leafblower.h:218`).
- `RK_ERR_STALL` means "valid weights at a constrained optimum"; `base.stall_kind`
  (1=weight-change, 2=KL, `src/types.hpp:112`) tells R which message to emit
  (`R/harvest.R:~715`).
- Unit-mode re-gating: `regate_unit_status` (`src/calib_dispatch.hpp:460`) demotes RK_OK →
  RK_ERR_STALL when post-finalize margin error exceeds `st.tol_abs`.
- R side: status 2/3 → `stop()`, status 4 → budget advice with stall heuristic
  (`R/harvest.R:~735-765`).
- C++ exceptions inside `C_rk_calibrate` are caught and converted to a string, so all RAII
  objects are destroyed before `Rf_error` longjmps (`src/r_bridge.cpp:~597`, R-exts §5.5).

## Cross-Cutting Concerns

**Logging:** `CalibState::log()` (`src/types.hpp:147`) — gated on `verbose > 0`, routed to a
caller-supplied `log_fn`, else `REprintf` (R) / `fprintf(stderr)` (`LBW_NO_R`).
**Validation:** three tiers — language layer (`R/harvest.R:443-540`,
`python/leafblower/_harvest.py:13`), ABI layer (`lbw::validate_calibrate_inputs`,
`src/validation.cpp`), solver pre-entry (`calib_validate_preentry`, `src/calib_validate.cpp`,
invoked from `solver_setup_ct`). Direct C-ABI callers get tiers 2 and 3 only — several branches
carry explicit "direct C API callers bypass R-layer validation" notes (`src/c_api.cpp:454`, `:469`).
**Numerical stability:** cancellation-free formulations are a project rule; `src/lbw_math.hpp`
provides `bulk_log` / `bulk_scaled_exp` with an AVX2 glibc-libmvec path gated on
`LBW_HAS_GLIBC_MVEC` (probed by `configure` for R and re-probed by `python/CMakeLists.txt` for
Python — the parity-critical macro).

---

*Architecture analysis: 2026-08-14*
