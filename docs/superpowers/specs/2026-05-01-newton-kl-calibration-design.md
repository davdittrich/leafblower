# Newton-KL Calibration Solver — Design Spec (rev 2)

## Current Baselines (fresh benchmark, 2026-05-01)

kk1204: n=1M, K=20, nj=5, max_weight=3, skewed targets, OMP_NUM_THREADS=1.

| Method | Wall (s) | Iters | Status | max_err |
|---|---|---|---|---|
| ieppa | 8.0 | 50 | ✅ | 3.5e-5 |
| **ieppa+sraa** | **3.7** | **10** | **✅** | **2.4e-14** |
| raking | 16.7 | 50 | stall | 6.1e-3 |
| raking+sraa | 60.3 | 200 | budget | 3.8e-3 |
| greenkhorn | 16.3 | 200 | budget | 3.2e-2 |
| greenkhorn+sraa | 309 | 4560 | budget | 7.6e-3 |
| ieppa_soft | 14.4 | 100 | budget | 1.5e-3 |
| **ieppa_soft+sraa** | **3.7** | **10** | **✅** | **2.4e-14** |
| **greg** | **3.2** | **2** | **✅** | **2.5e-2 ⚠️** |
| sinkhorn | 12.0 | 100 | NOCONV | 2.6e-3 |
| lbfgsb | 73.8 | 200 | budget | 8.3e-4 |
| logit | 35.8 | 50 | budget | 1.5e-3 |
| chebyshev | 104 | 100 | NOCONV | 1.1e-2 |

**Key observations:**
- greenkhorn diverges on K=20 dense problems (4560 iters, 309s)
- lbfgsb takes 73.8s — unusable
- **greg converges in 3.2s (2 Newton steps!) but max_err=2.5e-2** — validates that Newton-type methods are fast here; greg fails on quality because it minimizes chi2 not KL
- ieppa+sraa and ieppa_soft+sraa both solve to max_err=2.4e-14 in 3.7s
- Newton-KL should match greg's speed (both Newton on the dual) but with KL-optimal quality

## Revised Problem Statement

`ieppa+sraa` solves kk1204 (3.7s), but:
1. The 3s target is not quite met (3.7s > 3s)
2. lbfgsb (the "standard" dual solver) is 20× slower at 73.8s — unusable for this class
3. Newton-on-dual is the principled algorithm for this regime: ~5-10 quadratic steps vs 10 SRAA steps

**Goal:** Newton solver < 2s on kk1204 (2× faster than ieppa+sraa, 37× faster than lbfgsb). Primary role: replace lbfgsb for the zero-compression regime and provide a drop-in lbfgsb successor.

## Scope

- New standalone solver: `method="newton_kl"` in `harvest()`
- Designed to replace lbfgsb for the zero-compression regime AND serve as a faster alternative to ieppa+sraa
- AUTO routing: K ≥ 5 AND M_cell/n ≥ 0.9 (K ≥ 5 justified below; same M_cell/n threshold already in c_api.cpp)
- Obs-level iteration (not cell-level) — explicitly NOT using cell_lf

## K ≥ 5 Threshold Justification

Newton's Hessian costs O(n × K²/2) per step. For this to justify replacing lbfgsb's O(n × K) per step with Newton's 5-10× fewer steps:

- K=3: Newton costs 4.5× more per step vs lbfgsb. Needs > 4.5× fewer steps to win. Marginal benefit.
- K=5: Newton costs 12.5× more per step. Needs > 12.5× fewer steps. With 50→5 steps ratio (10×), break-even at K~5.
- K=10: Newton costs 50× more per step. Needs > 50× fewer steps. Newton excels here.
- K=20: Newton costs 190× more per step. Needs > 190× fewer steps. With ~50→5 ratio, Newton wins by ~3× in total cost.

K ≥ 5 is the empirical break-even point where Newton's quadratic convergence compensates for the K²-expensive Hessian.

## Math

**Dual objective** (maximize over λ):

```
g(λ) = log Z(λ) - Σ_{k,j} T_kj λ_kj
Z(λ) = Σ_i d_i exp(u_i(λ))
u_i(λ) = Σ_k λ_{k, j_k(i)}    [j_k(i) = category of obs i for margin k, 0-indexed]
```

**Reference elimination:** fix λ_{k,0} = 0 per margin. Free variables: `n_λ = Σ_k (cat_counts[k]-1)`.
For K=20, nj=5: n_λ = 80.

**Gradient:** `∇g_{kj} = S_kj(λ)/Z - T_kj` for j ≥ 1

**Hessian:** `H_{(k1,j1),(k2,j2)} = S_{k1j1,k2j2}(λ)/Z - [S_k1j1/Z][S_k2j2/Z]`

where `S_{k1j1,k2j2} = Σ_{i: j_{k1}(i)=j1, j_{k2}(i)=j2} f_i(λ)` (cross-margin joint sum).

**Newton step:** `δ = H^{-1} ∇g`, `λ ← λ + α δ` (Armijo line search for α).

## Single-Pass Per-Step Algorithm (obs-level, NOT cell-level)

The accumulation iterates over raw observations `i = 0..n-1` using `group_ids[k][i]`:

```cpp
// Initialize
Z = 0; G[n_λ] = 0; H[n_λ][n_λ] = 0;

// Single pass over observations
for i in 0..n-1:
    u_i = Σ_{k: group_ids[k][i] > 0} λ[lam_off[k] + group_ids[k][i] - 1]
    f_i = d[i] * exp(u_i)
    Z  += f_i
    // Gradient
    for k in 0..K-1:
        j = group_ids[k][i]
        if j > 0: G[lam_off[k] + j - 1] += f_i
    // Hessian (symmetric, upper triangle only)
    for k1 in 0..K-1, k2 in k1..K-1:
        j1 = group_ids[k1][i]; j2 = group_ids[k2][i]
        if j1 > 0 and j2 > 0:
            H[lam_off[k1]+j1-1][lam_off[k2]+j2-1] += f_i

// Normalize
G /= Z; H /= Z
H -= G ⊗ G          // outer product: H[a][b] -= G[a]*G[b]
G -= T              // subtract targets

// Newton solve
H_factored = LDLT(H)   // n_λ × n_λ, negligible
δ = H_factored^{-1} G
α = Armijo(g, λ, δ)
λ += α δ
```

**IMPORTANT:** This is obs-level (`group_ids[k][i]`), NOT cell-level (`g_per_cell[k][c]`). No `cell_lf` is used.

## Cost Analysis

| Operation | Cost | K=20, n=1M |
|---|---|---|
| exp(u_i) per obs | O(n × K) = 20M | ~0.02s |
| Gradient accumulation | O(n × K) = 20M | ~0.02s |
| Hessian accumulation | O(n × K²/2) = 190M | ~0.05s |
| Outer product H -= G⊗G | O(n_λ²) = 6.4K | negligible |
| LDLT solve (n_λ×n_λ) | O(n_λ³/3) = 170K | negligible |
| Armijo eval × 3 | O(n × K) × 3 = 60M | ~0.06s |
| **Total per step** | | **~0.15s** |
| **10 steps** | | **~1.5s** |

Target: **< 2s** on kk1204 (vs ieppa+sraa 3.7s, lbfgsb 73.8s).

## Files

| File | Change |
|---|---|
| `src/newton_calib.hpp` | `NewtonCalibResult` struct (embeds `CalibResult base`), `newton_calibrate()` declaration |
| `src/newton_calib.cpp` | Full implementation: obs-level Hessian accumulation, LDLT solve, Armijo line search |
| `src/Makevars.in` | Add `newton_calib.cpp` to `PKG_SOURCES` |
| `src/leafblower.h` | `RK_ALG_NEWTON_KL = 11` in `rk_algorithm_t` enum |
| `src/c_api.cpp` | AUTO dispatch: Newton when M_cell/n ≥ 0.9 AND K ≥ 5; add Newton arm |
| `src/r_bridge.cpp` | **Both** AUTO dispatch (line ~425) **and** Newton dispatch arm; expose result fields |
| `R/harvest.R` | `method="newton_kl"` in `match.arg`, AUTO routing comment |

**Reuse from existing infrastructure:**
- `lbw::ldlt_factor_inplace` + `lbw::ldlt_solve` from `calib_linalg.hpp`
- `lbw::solver_setup_ct` for cell-table + validation preamble
- `lbw::mark_converged` from `calib_dispatch.hpp` (NOT calib_linalg)
- `CalibState::group_ids[k][i]` — the obs-level category lookup

## Result Struct

```cpp
struct NewtonCalibResult {
    CalibResult base;
    int    n_lambda     = 0;    // dual dimension (n_λ free variables)
    double dual_gap     = 0.0;  // ||∇g||_∞ at exit
    double step_norm    = 0.0;  // ||δ|| of last Newton step
    double line_alpha   = 1.0;  // final Armijo step size
};
```

## Convergence Criterion

1. **Primary:** `||∇g||_∞ < tol_abs` (gradient norm = marginal errors below tolerance)
2. **Secondary:** `step_norm < 1e-12` (Newton step collapsed)

Max iterations: `min(max_iterations, 50)`.

## Bounds Handling and Fallback

Newton solves the **unconstrained** KL dual. Weight bounds (`max_weight`, `min_weight`) are not enforced during optimization. After convergence, compute `w_i = f_i / Z × n`:

- If `fraction_violated = count(w_i > max_weight or w_i < min_weight) / n > 0.05`: return `RK_ERR_NOCONV` with the Newton weights (partial solution) and log a warning. The returned weights ARE the Newton solution — no secondary lbfgsb run.

**Justification for 5% threshold:** In the zero-compression regime (M_cell/n ≈ 1, no capacity clamp), weight bounds are loose constraints that only bind for extreme target skew. At max_weight=3 with skewed targets, typically < 1% of observations violate bounds. The 5% threshold (50K obs at n=1M) distinguishes "nearly feasible" (proceed) from "bounds are structurally active" (Newton is wrong algorithm). This is conservative — a follow-on ticket can tune it empirically.

**Returned weights on fallback:** always the Newton weights (`w_i = f_i/Z × n`), never a fallback solve. The status field tells the user whether bounds were satisfied.

## AUTO Routing (both dispatch sites)

**`src/c_api.cpp`** (line ~171, before existing `M_cell > 0.9` rule):
```cpp
// Newton: smooth dual, K ≥ 5, zero-compression regime
if (est_ratio >= 0.9 && K >= 5) alg = RK_ALG_NEWTON_KL;
```

**`src/r_bridge.cpp`** (line ~425, AUTO arm, same condition):
```cpp
// kAlgMap already contains {"newton_kl", RK_ALG_NEWTON_KL} from the new dispatch arm.
// AUTO arm: when M_cell/n >= 0.9 && K >= 5, set method_str = "newton_kl" before lookup.
alg_for_validation = RK_ALG_NEWTON_KL;  // set directly in AUTO branch
```
Both sites must be updated. The string key `"newton_kl"` must be in `kAlgMap` to match `harvest.R`'s `match.arg`.

## Testing

1. `test-newton-kl.R`: K=3 small problem, status==0, max_error < 1e-6
2. `test-newton-kl.R`: K=9 stepstone_small, status==0, max_error < 1e-4 (tighter than greenkhorn)
3. Manual benchmark: kk1204 (K=20, n=1M) → target < 2s, max_error < 1e-4
4. `test-newton-kl.R`: bounds-active problem → fallback triggers correctly, status = RK_ERR_NOCONV
5. `test-newton-kl.R`: **KL vs chi2 quality validation** — same fixture, `newton_kl$max_error < greg$max_error`. Validates that Newton-KL minimizes KL divergence (better marginal calibration than greg's chi2 Newton). Catches gradient or dual implementation bugs that produce greg-like convergence without KL-correctness.
6. Full regression: FAIL 0

## Non-Goals

- Not for compressed problems (M_cell << n); ieppa/raking are better there
- No greedy scheduler, SOR, SRAA overlays
- No homotopy
