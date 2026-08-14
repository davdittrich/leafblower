# Coding Conventions

**Analysis Date:** 2026-08-14

Three layers, three convention sets, one shared C++17 core:

| Layer | Location | Language |
|-------|----------|----------|
| Core | `src/` | C++17, `lbw` namespace |
| R API | `R/`, `src/r_bridge.cpp` | R + roxygen2 |
| Python API | `python/leafblower/`, `python/leafblower/_bindings.cpp` | Python 3.9+ + pybind11 |

## Naming Patterns

**C++ files (`src/`):**
- One solver per file pair: `<solver>.cpp` + `<solver>.hpp` (`oris.cpp`/`oris.hpp`, `raking.cpp`, `chebyshev.cpp`, `greenkhorn.cpp`, `greg.cpp`, `logit.cpp`, `newton_calib.cpp`, `sinkhorn.cpp`).
- Cold-path splits get a suffix and a private header: `oris_finalize.cpp`, `oris_trajectory.cpp`, `oris_internal.hpp`.
- Shared headers use a domain prefix: `calib_dispatch.hpp`, `calib_linalg.hpp`, `calib_validate.hpp`, `lbw_math.hpp`.
- Public C ABI header: `src/leafblower.h`.

**C++ symbols:**
- Types: `PascalCase` with an `Lbw` suffix on config structs — `HomotopyConfigLbw`, `SchedulerConfigLbw`, `EtaScheduleConfigLbw` (`src/types.hpp:15-36`).
- Scoped enums: `enum class SchedulerMode : int { ROUND_ROBIN = 0, GREEDY = 1 }` (`src/types.hpp:22`) — always `: int`, always explicit values (values are ABI).
- C-ABI enums/typedefs: `rk_`-prefixed, `RK_`-prefixed members — `rk_algorithm_t`, `rk_scheduler_t`, `RK_SCHED_ROUND_ROBIN` (`src/leafblower.h:18-19`).
- Compile-time constants: `k`-prefixed `inline constexpr double` — `kUnboundedSentinel` (`src/calib_dispatch.hpp:312`), `kMinSafeTotalWeight` (`:321`), `kMetricEps` (`:252`).
- Functions: `snake_case` (`finalize_weights`, `select_metric`, `oris_solve`).

**R (`R/`):**
- One exported concept per file: `harvest.R`, `design_effect.R`, `diagnose_weights.R`, `anesrake.R`, `weighted_pct.R`, `current_miss.R`, `na_bin.R`, plus `zzz.R`.
- Functions and arguments: `snake_case` (`design_effect`, `effective_sample_size`, `get_current_miss`, `min_weight`, `add_na_proportion`).
- `.Call` entry points: `C_rk_*` registered symbols (`C_rk_calibrate`, `C_rk_design_effect`), never string-name `.Call`.

**Python (`python/leafblower/`):**
- Private implementation modules are underscore-prefixed: `_harvest.py`, `_design_effect.py`, `_bindings.cpp` → `_leafblower` extension.
- `__init__.py` is a 5-line re-export façade with an explicit `__all__`.
- Module-private constants `_UPPER_SNAKE`: `_TARGET_SUM_TOL`, `_METRIC_MAP`, `_RULE_MAP`, `_KNOWN_CONVERGENCE_KEYS` (`python/leafblower/_harvest.py:10-52`).
- Private helpers `_snake_case`: `_validate_pos_scalar`, `_compute_sparseness_diag`, `_parse_convergence`.

## Code Style

**Formatting:**
- No formatter or linter is configured — no `.lintr`, no `clang-format`, no `ruff`/`black` config. Match the surrounding file exactly.
- R: 2-space indent, `<-` assignment, aligned `=` in multi-line calls, single-statement `if` bodies on the next line without braces (`R/design_effect.R:39-48`).
- C++: 4-space indent, aligned struct-field initializers in a column block (`src/oris.hpp:9-45`), section banners drawn with box characters (`// ── ORIS diagnostics ──`).
- Python: `from __future__ import annotations` at top, type hints on public and helper signatures, keyword-only args after `*` (`python/leafblower/_harvest.py:1-14`).

**Optimization flags:** the package sets no `-O` level of its own — `configure` and `src/Makevars.in` leave it to R's `$(CXXFLAGS)`. Never add `-O*` to `PKG_CXXFLAGS` (CRAN `tools:::.check_make_vars` rejects it).

## Numerical Style — No Cancellations

Project rule: compute variances, covariances and residuals cancellation-free. Never write a difference of near-equal quantities where an algebraically equal product exists.

Canonical example, `src/newton_calib.cpp:417-429`:

```cpp
// the DIAGONAL variance p_a − p_a² is a catastrophic cancellation as p_a→1 …
// Compute it cancellation-free as p_a·(1−p_a) — one rounding (~0.5 ulp)
H[row_off + a] = G[a] * (1.0 - G[a]);        // cancellation-free variance
for (int b = a + 1; b < n_lam; ++b)
    H[row_off + b] -= G[a] * G[b];            // off-diagonal covariance
```

The same file documents why the off-diagonals stay a plain subtraction (independent quantities, not self-cancellation) and applies the stable quadratic root instead of `(−b+√disc)/(2a)` (`src/newton_calib.cpp:237`). Log-sum-exp shifts subtract `u_max` and note that the factor cancels in all ratios (`:356`, `:857`). Rows with zero design weight are excluded from LSE sums to avoid `0*inf`.

**Corollary:** every non-obvious numerical rewrite carries an inline comment naming the failure mode it avoids and the ticket ID that introduced it (`CR-C18 (kxna.18)`).

## Comments

- Comments explain *why*, and are dense in solver math. A comment that pins a verified-correct formula is a guard — `CLAUDE.md` lists two formulas (chebyshev Mehrotra corrector, oris ALM Newton step) that must NOT be "fixed".
- Ticket IDs are cited inline (`CR-F7 (dtkn.7)`, `CR-D9 (j7x8.9)`, `y2ks.8`) so a reader can recover the rationale from beads.
- Struct fields carry trailing comments stating units and when the field is inert (`src/oris.hpp:14-21`).
- R: roxygen2 blocks are exhaustive — `R/harvest.R` opens with ~250 lines of `#'` covering every argument, per-method defaults and references. `@references` cite the primary literature (`R/design_effect.R:29-34`).
- Python docstrings mirror the R rule they implement, with R file:line pointers (`python/leafblower/_harvest.py:15-27`).

## Error Handling

**R — base R, not `cli`.** There is no `cli` usage anywhere in `R/` (verified: zero `cli::` matches across all eight files). Use:
- `stop("<fn>: <what> ...; got: ", value)` for contract violations, with `call. = FALSE` where the caller frame is unhelpful (`R/harvest.R:508`).
- `warning(sprintf(...))` for degraded-but-proceeding paths (unknown args ignored, deprecated method names, budget exhaustion) (`R/harvest.R:306, 357, 409`).
- Messages are prefixed with the function or package name and always echo the offending value.
- Validate at the R boundary before `.Call` — both entry paths, not just one (`R/design_effect.R:37-48`).

**C++ / R bridge — RAII-safe longjmp.** `Rf_error` longjmps past C++ destructors (R-exts §5.5). The convention in `src/r_bridge.cpp` is: set a `pre_error` string, let the enclosing scope destroy every RAII object at its `}`, then call `Rf_error` (`src/r_bridge.cpp:270, 642-648`). Direct `Rf_error` is permitted only where no heap-backed local is live (`:46, :52`).

**Python — `ValueError` mirroring R text.** `_harvest.py` raises `ValueError` (31 sites) and `_design_effect.py` (10 sites) with messages formatted to match the R `stop()` string, and `warnings.warn(..., UserWarning, stacklevel=2)` where R warns (`python/leafblower/_harvest.py:28-32`). Parity of *behaviour on bad input* is part of the contract, not just parity of weights.

**Sentinels over exceptions in the core:** the solver core returns status codes through `CalibResult`; degenerate cases are handled numerically (`total_w < kMinSafeTotalWeight` is left unscaled to avoid subnormal overflow) rather than by throwing.

## Result-Struct Access

`CalibResult` fields live at `res.base.*`, not `res.*` (since ztid.4). Direct field access on the outer struct compiles and breaks silently. Solver-specific extras sit alongside `base` in the wrapper struct (`src/oris.hpp:7-9`).

## Module Design

**C++ placement rules (hard constraints, no LTO in this build):**
- Shared solver helpers go in `calib_dispatch.hpp` — the canonical home. Do NOT add shared logic to an individual solver file.
- CellTable-specific helpers go in `cell_table.hpp`.
- Only COLD (once-per-solve) code may move to a new TU. HOT per-iteration loops and `static inline` kernels must stay co-located with their caller — cross-TU calls do not inline without LTO, and `-flto` is absent from `configure`/`Makevars`.

**Two build sites for `src/*.cpp`:** R auto-globs the directory (the `PKG_SOURCES` list in `Makevars.in` is decorative). The Python build does not glob — a new `.cpp` MUST be added to `CORE_SOURCES` in `python/CMakeLists.txt` or the pybind11 link fails with undefined symbols.

**Exports:**
- R: explicit `export()` lines in `NAMESPACE` (8 exported functions) plus `useDynLib(leafblower, .registration = TRUE)`.
- Python: explicit `__all__` in `__init__.py`; everything else is underscore-private.

## Lambda Capture Footgun

In `raking.cpp`, declare bool guards BEFORE the `auto F_eval = [&]` definition. `[&]` captures only variables in scope at the definition site, not at the call site.

## Weight-Finalization Order

`Σw = n` is enforced at exit via `lbw::finalize_weights[_buf]` in `calib_dispatch.hpp`: a single pre-bounds scale to `n`, THEN `bounds_mode` dispatch (`cell` = count-only, `unit` = per-cell water-fill). The sanctioned order is normalize → bounds. Renormalizing AFTER water-fill is FORBIDDEN — it silently breaks the `bounds_mode="unit"` clamps.

## Argument-Name Footguns

- `design_weights=` is the `harvest()` argument for per-observation design weights (`d_i` in `Z = Σ d_i exp(u_i)`). There is no `weights=` argument; it lands in `...` and is silently ignored.
- Version numbers are duplicated in `DESCRIPTION` and `python/pyproject.toml` with no automation — bump both by hand.
- Algorithm slot 2 in `rk_algorithm_t` is reserved (LBFGSB removed); do not reuse.

---

*Convention analysis: 2026-08-14*
