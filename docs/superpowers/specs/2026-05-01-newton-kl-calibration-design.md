# Newton-KL Calibration Solver — Design Spec

## Problem

The K=20, zero-compression regime (M_cell/n ≈ 1) exhibits pure Sinkhorn slow-rate convergence (errRp ~ iter^{-0.5}) in iEPPA. At 500 iterations, errRp stalls at ~1.15e-3. Target: converge to 1e-4 in < 3s.

The KL-calibration dual is smooth and K×max_cats-dimensional (~80 variables for K=20, 5 cats/margin). Full Newton on this dual converges in 5-10 steps (quadratically). Per-step cost is O(n × K²/2) — one pass over n observations accumulating the Fisher information matrix. For K=20, n=1M: ~190M ops → ~0.1s/step → ~0.5-1s total.

## Scope

- New standalone solver: `method="newton"` in `harvest()`
- Designed to replace lbfgsb for the zero-compression regime
- AUTO routing: K ≥ 5 AND M_cell/n ≥ 0.9
- Linear path only (log path not relevant at M_cell/n ≥ 0.9 since n/M_cell ≤ 1.1 < 2.0)

## Math

**Dual objective** (maximize over λ):

```
g(λ) = log Z(λ) - Σ_{k,j} T_kj λ_kj
Z(λ) = Σ_i d_i exp(u_i(λ))
u_i(λ) = Σ_k λ_{k, j_k(i)}    [j_k(i) = category of obs i for margin k]
```

**Reference elimination:** fix λ_{k,0} = 0 (one category per margin = reference). Free variables: `n_λ = Σ_k (cat_counts[k]-1)`.

**Gradient (marginal error):**
```
∇g_{kj} = S_kj(λ)/Z - T_kj    [j ≥ 1; reference contributions implicit]
```

**Hessian (Fisher information):**
```
H_{(k1,j1),(k2,j2)} = S_{k1j1,k2j2}(λ)/Z - [S_k1j1(λ)/Z] × [S_k2j2(λ)/Z]
```
where `S_{k1j1,k2j2}(λ) = Σ_{i: j_k1(i)=j1, j_k2(i)=j2} f_i(λ)`.

**Newton step:** `δ = H^{-1} ∇g`, `λ ← λ + α δ` (Armijo line search for α).

## Single-Pass Per-Step Algorithm

```
f_i = d_i × exp(Σ_k λ_{k, j_k(i)})       // calibrated weights (unnormalized)
Z   = Σ_i f_i
G   = Σ_i f_i × A_i / Z - T               // gradient; A_i is free-variable indicator vector
H   = Σ_i f_i × A_i A_i^T / Z - G G^T    // Hessian (Fisher information)
δ   = LDLT(H)^{-1} G                      // Newton direction (n_λ × n_λ solve)
α   = Armijo(g, λ, δ)                      // line search
λ  += α δ
```

All of G, H accumulated in **one pass** over n observations.

## Cost Analysis

| Operation | Cost | K=20, n=1M |
|---|---|---|
| exp(u_i) per obs | O(n × K) = 20M | ~0.02s |
| Gradient accumulation | O(n × K) = 20M | ~0.02s |
| Hessian accumulation | O(n × K²/2) = 190M | ~0.05s |
| LDLT solve (n_λ×n_λ) | O(n_λ³/3) = 170K | negligible |
| Line search (3 evals) | O(n × K) × 3 = 60M | ~0.06s |
| **Total per step** | | **~0.15s** |
| **10 steps** | | **~1.5s** |

Convergence: quadratic local, ~5-10 steps to 1e-4 from cold start.

## Files

| File | Change |
|---|---|
| `src/newton_calib.hpp` | `NewtonCalibResult` struct (embeds `CalibResult base`), `newton_calibrate()` declaration |
| `src/newton_calib.cpp` | Full implementation: Hessian accumulation, LDLT solve, Armijo line search, convergence check |
| `src/leafblower.h` | `RK_ALG_NEWTON = 11` in `rk_algorithm_t` enum |
| `src/c_api.cpp` | Dispatch arm for `RK_ALG_NEWTON` |
| `src/r_bridge.cpp` | Dispatch arm, expose `NewtonCalibResult` fields |
| `R/harvest.R` | `method="newton"` in `match.arg`, AUTO routing logic |

**Reuse from existing infrastructure:**
- `lbw::ldlt_factor_inplace` + `lbw::ldlt_solve` from `calib_linalg.hpp`
- `lbw::solver_setup_ct` for cell-table + validation preamble
- `lbw::mark_converged` for result struct population
- `CalibState` unchanged (no new fields needed)

## Result Struct

```cpp
struct NewtonCalibResult {
    CalibResult base;           // status, iterations, max_error, convergence_*
    int    n_lambda     = 0;    // dual dimension (= n_λ free variables)
    double dual_gap     = 0.0;  // primal-dual gap at exit
    double step_norm    = 0.0;  // ||δ|| of last Newton step
    double line_alpha   = 1.0;  // final Armijo step size
};
```

## Convergence Criterion

Per-iteration (each Newton step counts as one iteration):
1. **Primary:** `||∇g||_∞ < tol_abs` (gradient norm below user tolerance)
2. **Secondary:** `max_error < pct_tol` (marginal error in user's pct_tol range)
3. **Tertiary:** step norm `||δ|| < 1e-12` (Newton step collapsed)

Max iterations: `min(max_iterations, 50)` — Newton should never need more than 50 steps.

## AUTO Routing

In `c_api.cpp` AUTO selection, add before existing rules:

```cpp
// Newton: zero-compression regime (M_cell/n ≥ 0.9, K ≥ 5)
// Full Newton on n_λ-dim dual outperforms iEPPA/raking/lbfgsb here.
if (static_cast<double>(M_cell_est) / n >= 0.9 && K >= 5) {
    alg = RK_ALG_NEWTON;
}
```

## Bounds and Constraints

- `max_weight` / `min_weight`: enforced in the primal (calibrated weights `w_i = f_i/Z`) by adding a capacity constraint to the dual. The standard KL-calibration dual assumes no weight bounds. With bounds, the problem is a constrained LP; Newton no longer applies directly.
- **Design decision:** For the zero-compression regime (kk1204), bounds are typically loose (`max_weight=3`). Newton applies to the unconstrained dual and the primal solution is checked post-hoc. If any w_i violates bounds: fall back to lbfgsb.
- Fallback condition: `fraction_violated > 0.01` after Newton convergence → report `RK_ERR_NOCONV` with status note.

## Testing

1. `test-newton-kl.R`: K=3 trivial problem, status==0, max_error < 1e-6
2. `test-newton-kl.R`: K=9 stepstone_small, status==0, competitive with greenkhorn
3. Manual benchmark: kk1204 (K=20, n=1M) → target < 3s, max_error < 1e-4
4. Full regression: FAIL 0 after adding newton dispatch

## Non-Goals

- Not a drop-in replacement for lbfgsb on compressed problems (M_cell << n)
- Not compatible with homotopy, SOR, SRAA overlays (Newton is its own convergence mechanism)
- No greedy scheduler (Newton has no sweep structure)
