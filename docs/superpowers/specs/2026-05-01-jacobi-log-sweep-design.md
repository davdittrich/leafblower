# iEPPA Log-Path Jacobi Sweep — Design Spec

## Problem

The iEPPA log-path margin sweep (`apply_single_margin_log`) uses **Gauss-Seidel** semantics: after updating `lf[k][j]` for margin k, `cell_lf[c]` is patched in-place for every cell c in bucket (k, j). These patches are scattered writes (bucket-order), causing poor cache utilization for large M_cell.

**Jacobi** alternative: freeze `cell_lf` at outer-iteration start, sweep all K margins using the frozen snapshot (no mid-sweep patches), then rebuild `cell_lf` once at iter end with a sequential O(K·M_cell) pass. Cache-friendly, but typically requires ~2× more outer iterations to converge. Net wall-time change is empirically uncertain — this design builds and benchmarks both variants to decide.

## Scope

- **In scope**: log path only. The linear path already has O(1) leave-one-out via `X_cur[c] / f_lin[k][j]`; Jacobi adds no value there.
- **Out of scope**: linear path, SRAA, SOR, homotopy, ieppa_soft ALM.

## Changes

### 1. `src/types.hpp` — add flag to `CalibState`

```cpp
bool jacobi_log = false;  // log-path: freeze cell_lf at iter start (Jacobi) vs patch inline (GS)
```

### 2. `R/harvest.R` — new parameter

```r
jacobi_sweep = FALSE   # bool; only affects log path (high-compression surveys n/M_cell >= 2)
```

Wired to C bridge as `as.integer(jacobi_sweep)`. Validated: must be TRUE/FALSE.

### 3. `src/r_bridge.cpp` — wire to CalibState

```cpp
st.jacobi_log = (INTEGER(jacobi_sweep_sexp)[0] != 0);
```

### 4. `src/ieppa.cpp` — Jacobi implementation in log-path outer loop

**New workspace vector** (declared at homotopy-level scope, outside iter loop):
```cpp
std::vector<double> cell_lf_frozen;  // size ct.M_cell; only allocated when jacobi_log
if (st.jacobi_log) cell_lf_frozen.assign(ct.M_cell, 0.0);
```

**At top of each outer iteration** (before margin sweeps, log path only):
```cpp
if (st.jacobi_log && !use_linear) {
    std::copy(cell_lf.begin(), cell_lf.end(), cell_lf_frozen.begin());
}
```

**Inside `apply_single_margin_log`** — use `cell_lf_frozen` instead of `cell_lf` when flag is set:
```cpp
// compute S_kj using frozen or live cell_lf depending on flag
const double* clf = st.jacobi_log ? cell_lf_frozen.data() : cell_lf.data();
// ... use clf[c] in bucket sum ...
```

**Skip mid-sweep cell_lf patches** when Jacobi:
```cpp
if (!st.jacobi_log) {
    // existing GS patch: update cell_lf[c] for cells in bucket
    cell_lf[c] += lf_delta;
}
```

**End-of-sweep cell_lf rebuild** when Jacobi (replaces the K scattered patch operations):
```cpp
if (st.jacobi_log && !use_linear) {
    // Sequential O(K*M_cell) rebuild
    std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
    for (int k = 0; k < st.K; k++) {
        const int off = cat_offset[k];
        const int* gk = ct.g_per_cell[k].data();
        for (int c = 0; c < ct.M_cell; c++) {
            int g = gk[c];
            if (g >= 0) cell_lf[c] += lf[off + g];
        }
    }
}
```

### 5. `benchmarks/jacobi_sweep_study.R` — benchmark script

3×3×2 = 18 cells:
- K ∈ {3, 6, 9}
- compression ∈ {2×, 10×, 50×} (M_cell/n ∈ {0.5, 0.1, 0.02})
- tightness ∈ {loose: max_weight=5, tight: max_weight=1.5}

For each cell, run both `jacobi_sweep=FALSE` (GS, baseline) and `jacobi_sweep=TRUE` (Jacobi) with `convergence=list(absolute=1e-4)`. Record `iters`, `wall_s`. Derive:
- `iter_inflation = iters_jacobi / iters_GS`
- `wall_ratio = wall_jacobi / wall_GS`

Save results to `benchmarks/jacobi_sweep_study.rds`.

## Decision Gate

After benchmark:
- **Ship Jacobi as default** (`jacobi_log = true`): `median(wall_ratio) < 1.0` across grid cells where log path was active (compression ≥ 2×).
- **Remove flag, close as not-worth-it**: otherwise.
- **Flag retained but off-by-default**: if wall_ratio shows improvement only at specific (K, compression) — expose as `jacobi_sweep = FALSE` with note in docs.

## Files Changed

| File | Change |
|---|---|
| `src/types.hpp` | `bool jacobi_log = false` on `CalibState` |
| `src/ieppa.cpp` | `cell_lf_frozen` workspace, conditional freeze/rebuild in log path |
| `src/r_bridge.cpp` | wire `jacobi_sweep_sexp` → `st.jacobi_log` |
| `R/harvest.R` | `jacobi_sweep = FALSE` parameter |
| `benchmarks/jacobi_sweep_study.R` | benchmark script (new) |

## Correctness Guarantee

Jacobi is mathematically equivalent to a different fixed-point iteration (all-margins-old vs one-at-a-time-new). It is NOT identical to GS and may converge to a slightly different fixed-point. Acceptance criterion: `max_error` at convergence must be within `1e-4` of GS result on the same fixture (verified in T3 test).

## Testing

- T2: unit test: run same problem GS vs Jacobi to tol=1e-4; assert `max_error` within 1e-4 of each other and both status==0.
- T3: benchmark 18-cell grid (output: decision CSV).
- T4: if shipping: existing test suite must remain FAIL 0.
