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

**Newton step (LM-damped):** `H_damped · δ = ∇g` where `H_damped[a,a] = max(H[a,a]·(1+lm_mu), lm_mu·d_floor)`, `d_floor = mean(diag(H))`; `λ -= α δ` (Armijo line search; gain-ratio adaptive lm_mu — see "Levenberg-Marquardt damping" subsection below).

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

### Levenberg-Marquardt damping

Sample H rank-collapses dynamically as λ drifts (severe-skew K=20 case): only a small subset of obs has meaningful `exp(u_i − u_max)`, so the empirical Hessian effectively builds from a shrinking subset. Plain Newton overshoots, gradient stalls. Solution: scale-invariant damping with additive floor.

**Damped Hessian:**
- `H_damped[a,a] = max(H[a,a]·(1+lm_mu), lm_mu·d_floor)` where `d_floor = mean(diag(H))`.
- `(1+lm_mu)` preserves per-coordinate scale where H is well-conditioned.
- `max(..., lm_mu·d_floor)` provides a non-zero damping floor where `H[a,a]→0` (rank-collapsed direction).

**Adaptive `lm_mu` via Marquardt gain ratio:**
```
ρ = (g_curr − g_trial) / (α·G·δ − ½·α²·δᵀHδ)
```
- `ρ > 0.75` AND Armijo full accept (α=1) ⇒ `lm_mu ← max(lm_mu/3, 1e-12)` (recover Newton near optimum).
- `ρ < 0.25` OR Armijo line search exhausted ⇒ `lm_mu ← min(lm_mu·10, 1e12)` (poor model — back off toward gradient descent).
- Otherwise: keep `lm_mu` unchanged.

**Init:** `lm_mu = 1.0`. Bounds: `[1e-12, 1e12]`.

**Failed line search:** retry the LDLT solve with `lm_mu *= 10` in the same iter (max 3 retries). Do NOT advance λ with a tiny-step fallback. If three consecutive iterations fail line search, return `RK_ERR_NOCONV` with `lm_mu_final` set to the saturated value (used by AUTO-routing safety patch on Epic-D ESCAPE verdict).

See plan rev2 (`docs/superpowers/plans/2026-05-01-newton-kl-lm-plan.md`) for the full rationale and the failure modes ruled out (multiplicative-only damping, μ·I spherical damping, trust-region clip, tiny-step fallback).

### Target homotopy

LM damping stabilizes Newton inside its convergence basin but cannot deform the basin itself. Severe-skew K=20 problems and stepstone-style overlapping-margin fixtures have rank-deficient Hessians whose dual landscape has a finite-gap plateau — Newton converges to that plateau but no further. **Mechanism:** target homotopy provides a warm-start path: each intermediate problem's `λ*` is a starting point for the next, keeping Newton inside its quadratic basin throughout the descent to `λ*(T_0)`. Rank deficiency in `H` persists at every ε, but Newton's basin is wide enough when `λ_init` is close to `λ*`.

**Schedule:** ε ∈ {0.5, 0.1, 0.03, 0.01, 0.003, 0.0}.

**Per-ε targets:** `T_eps[k][j] = (1-ε)·T[k][j] + ε·(1/cat_counts[k])` (per-margin uniform; NOT joint uniform — joint would explode cell count).

**Feasibility guard:** before the first ε iteration, verify every (k,j) with `T_eps[k][j] > 0` has at least one observation `cat_k(i) == j` with `d_i > 0`. If any (k,j) fails, abort with `RK_ERR_BADARG` (smoothed target unreachable from sample).

**Inner solve:** existing LM-damped Newton iter loop, max 20 iters per ε. λ carries across ε boundaries (warm start).

**`lm_mu` schedule:** at first ε, `lm_mu = 1e-3`. At each subsequent ε, `lm_mu = max(lm_mu_final_prev / 3, 1e-6)` — keeps recovered damping while avoiding saturation lock-in.

**Final ε=0:** runs the same inner with original targets; expected to converge in 3-5 Newton iters from a near-optimal warm start.

**Recovery:** weights computed from final ε=0's converged λ. No homotopy bias in returned weights.

**New diagnostic:** `n_homotopy_levels_used` field on `NewtonCalibResult` records the count of ε levels that ran ≥1 inner Newton iter.

**AUTO routing impact:** none in the success/PARTIAL paths (homotopy is internal to `newton_calibrate`; signature unchanged). On ESCAPE verdict, AUTO routing in `c_api.cpp` adds a target-skew gate to route severely-skewed K≥5 problems to ieppa+sraa (homotopy code stays in `newton_calib.cpp` for moderate-skew users).

See plan rev 2: `docs/superpowers/plans/2026-05-01-newton-kl-homotopy-plan.md`.

### IEPPA warm-start

LM damping + best-iterate fallback stabilize Newton inside its convergence basin, but rank-deficient `H_pre` (overlapping-margin samples like stepstone) leaves a basin floor that Newton alone cannot cross from a `λ=0` cold start. **IEPPA (iterative proportional fitting / coordinate descent in log-factor space)** cannot overshoot across margin boundaries — its updates always land `λ` inside the correct basin. Running IEPPA for `K_warm` sweeps before LM-Newton produces a `λ_init` close enough to `λ*(T_0)` that Newton's quadratic basin captures the rest.

**Mechanism:** `K_warm = 8` IEPPA+SRAA sweeps with the original target `T` (NOT a smoothed `T_eps`). Capture the lf vector at IEPPA's best-iterate (`lf_best`, mirroring `W_best`). Convert to Newton's `λ` via `λ[lam_off_newton[k] + j - 1] = lf_best[cat_offset_ieppa[k] + j] - lf_best[cat_offset_ieppa[k] + 0]` for `j ≥ 1`. The constant offset `C = sum_k lf_best[cat_offset_ieppa[k] + 0]` shifts away under Newton's LSE; weights are scale-equivalent.

**Three-tier fallback:** SRAA-fail (non-finite lf or status ≥ 1) → plain IEPPA retry (`accelerate=FALSE`) → cold start (`λ = 0`). Strictly additive — warm-start can only improve, never regress. Diagnostic field `n_warmstart_iters_used` records which tier fired (8 = SRAA, 8 = plain IEPPA, 0 = cold).

**Recovery:** weights computed from final ε=0 Newton's `λ` (no homotopy smoothing). After warm-start, Newton runs the existing LM-damped inner with `lm_mu = max(lm_mu_after_first_iter / 3, 1e-6)` (see "Levenberg-Marquardt damping" subsection for the rest of the inner-loop semantics).

**Diagnostic field:** `n_warmstart_iters_used` (integer) on `NewtonCalibResult`, surfaced via `r_bridge.cpp pack_solver_result`. Replaces the now-stale `n_homotopy_levels_used` (Epic-B target homotopy was BLOCKED; see `docs/investigations/2026-05-01-newton-kl-homotopy-result.md`).

**AUTO routing impact:** none. Warm-start is internal to `newton_calibrate`; `c_api.cpp` AUTO selection logic is independent of inner-solver mechanics.

See plan rev 3: `docs/superpowers/plans/2026-05-01-newton-kl-ieppa-warmstart-plan.md`.

> **ERRATUM (2026-05-02):** This subsection was added in commit `32fcee6` documenting
> Epic-C's intended IEPPA warm-start mechanism. Empirical falsification at Epic-C
> closure (verdict BLOCKED): warm-start regresses T2 stepstone K=9 from 2.79e-4 (cold)
> to 4.39e-4. IEPPA's basin floor sits OUTSIDE Newton's quadratic basin on stepstone;
> warm-start lands Newton in a worse position than cold start. See
> `docs/investigations/2026-05-02-newton-kl-ieppa-warmstart-result.md` for full
> empirical analysis. Section retained for historical reference and potential future
> revival; NOT current behavior.

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
