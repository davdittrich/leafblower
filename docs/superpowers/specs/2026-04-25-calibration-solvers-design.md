# Calibration Solvers Redesign — Spec

**Date:** 2026-04-25
**Status:** Draft v1

## Problem

leafblower has three calibration solvers (iEPPA, raking, lbfgsb) but:
1. iEPPA uses Euclidean water-fill for capacity enforcement, breaking KL monotone descent — it never reaches its own theoretical minimum.
2. No solver targets the true minimum of max_err, mean_err, chi2, or grake_norm specifically.
3. Raking operates obs-level O(n) when cell-table O(M_cell) would give 100–300× speedup.

## Goals

1. Add `method="sinkhorn"` — true KL minimum subject to hard capacity constraints.
2. Add `method="chebyshev"` — true L∞ marginal error minimum (LP).
3. Add `method="greg"` — true chi2 minimum (Newton QP, GREG calibration).
4. Add `method="grake"` — true grake_norm minimum (normalized Chebyshev, LP).
5. Migrate `method="raking"` to cell-table — same algorithm, 100–300× speedup.
6. Change `method="ieppa"` default convergence from `max_err+improvement` to `kl+improvement` (consistent with what iEPPA actually minimizes).
7. All methods: hard capacity bounds satisfied at solver exit. No post-hoc clamping.
8. All speed-critical code in C++. No external LP/QP library dependencies.
9. Python parity via existing pybind11 path.

## Non-goals

- Replacing lbfgsb (it has its own use case: non-multiplicative calibration).
- Supporting non-cell-table obs-level paths for new methods (cell-table always).

---

## Architecture

### §1 Shared cell-table exit expansion

All cell-table solvers compute **cell mass multipliers** `M[c]` (ratio of final to initial cell mass). Obs-level weights are recovered at exit:

```
w_i = d_i × M[cell_of[i]]
normalize: w_i *= n / Σ w_i
```

where `d_i = start_weights[i]` (initial per-obs weight, may be non-uniform). Within-cell variation of d_i is preserved exactly: all obs in cell c get the same multiplier but their pre-existing d_i differences are maintained. This is the iEPPA pattern already in use (ieppa.cpp:1039).

**Capacity constraints** are box constraints on `X[c] = M[c] × X_init[c]`:
```
L_c = min_weight × n_per_cell[c]
U_c = max_weight × n_per_cell[c]
```
These are hard inequality constraints in each solver's formulation — not post-hoc. All methods satisfy `L_c ≤ X_final[c] ≤ U_c` at the fixed point.

### §2 Shared normal equations kernel (`src/calib_linalg.hpp`)

The **marginal incidence matrix A** (rows = margin-category pairs (k,j), cols = cells c) is binary: A[(k,j), c] = 1 iff cell c belongs to bucket (k,j). This matrix underlies all LP/QP/Newton steps.

The critical shared computation:
```
N = A × D × Aᵀ    (shape: n_cats_total × n_cats_total)
```
where D is diagonal on cell masses. For stepstone: N is 836×836. Computing N costs O(M_cell × K) using `cells_by_margin_cat`. Cholesky factor of N costs O(836³/3) ≈ negligible.

```cpp
// src/calib_linalg.hpp
namespace lbw {
void compute_normal_equations(const CellTable& ct, const double* D, double* N,
                               const int* cat_offset, int n_cats_total);
void cholesky_factor_inplace(double* N, int n);   // in-place LDLT
void cholesky_solve(const double* L, double* rhs, int n); // forward/back sub
} // namespace lbw
```

### §3 New method constants

```c
/* leafblower.h additions */
#define RK_ALG_SINKHORN   4   /* true KL min, Bregman Dykstra */
#define RK_ALG_CHEBYSHEV  5   /* true L∞ marginal error min, LP */
#define RK_ALG_GREG       6   /* true chi2 min, Newton QP */
#define RK_ALG_GRAKE      7   /* true grake_norm min, LP */
```

All methods return the standard `rk_result_t`. Two new result fields:
```c
double objective;        /* value of minimized metric at convergence */
int    objective_metric; /* CalibMetric enum: which metric was minimized */
```

---

## Method Designs

### §4 method="sinkhorn" — KL Bregman Dykstra

**Minimizes:** KL(X || X_init) subject to margin constraints + hard capacity bounds.

**Why iEPPA fails:** iEPPA's water-fill is a Euclidean projection onto the capacity box. In KL geometry, the correct projection is multiplicative (log-domain). Euclidean Dykstra on a KL problem breaks monotone descent, causing the two-accumulation-point oscillation (Gietl-Fröhlich 2013).

**Algorithm:** Alternating KL Bregman projections (Benamou 2015, Chizat 2016):

```
Initialize: lf[k][j] = 0  ∀k,j    (log Sinkhorn factors)
            a[c] = 0      ∀c      (log Dykstra corrections for capacity)
            X[c] = X_init[c]

Each outer iteration:
  ① K Sinkhorn sweeps (KL projection onto each margin constraint):
     For k = 1..K, j = 1..cat_counts[k]:
       S_kj = Σ_{c∈bucket(k,j)} X[c]
       ratio = T_kj × n / S_kj          (ratio > 0 for feasible margins)
       lf[k][j] += log(ratio)
       X[c] *= ratio  ∀c ∈ bucket(k,j)

  ② KL projection onto capacity box (replaces Euclidean water-fill):
     Bisection on μ ∈ ℝ s.t. Σ_c clamp(X[c]·exp(a[c]+μ), L_c, U_c) = n
     X_proj[c] = clamp(X[c]·exp(a[c]+μ), L_c, U_c)
     a[c] += log(X[c]) - log(X_proj[c])   ← log-domain Dykstra correction
     X[c] = X_proj[c]

  Check convergence (CalibRule on active CalibMetric).
  KL guaranteed monotone decreasing each outer iteration.
```

**Bisection bounds for μ:** μ ∈ [log(L_min/X_max), log(U_max/X_min)]. Bisection to tol=1e-12 needs ~40 iterations over M_cell evaluations: O(M_cell × 40) per capacity step.

**Convergence:** Geometric O(κ^t) to unique KL minimum satisfying all constraints. κ < 1 depends on capacity tightness.

**Convergence metric:** default `kl+improvement+0.001`. Supports all 6 metrics via improvement/plateau rules (same framework as other methods). `kl+improvement` is consistent with what the algorithm actually minimizes.

**Storage:** lf[K×max_cats], a[M_cell] — same order as iEPPA.

**Complexity per outer iter:** O(K × M_cell) Sinkhorn + O(M_cell × 40) bisection = O(K × M_cell). Asymptotically same as iEPPA.

**References:** Benamou et al. SIAM J. Sci. Comput. 2015; Chizat et al. Math. Comp. 2016; implemented in POT/GeomLoss for general OT with capacity.

---

### §5 method="chebyshev" — LP for L∞ calibration (custom IPM)

**Minimizes:** δ = max_{k,j} |S_kj/n − T_kj| subject to hard capacity bounds.

**LP formulation:**
```
min  δ
s.t. Σ_{c∈(k,j)} X[c] − δ·n ≤  T_kj·n     ∀k,j  (2×Σcat_counts constraints)
    −Σ_{c∈(k,j)} X[c] − δ·n ≤ −T_kj·n     ∀k,j
     L_c ≤ X[c] ≤ U_c                       ∀c    (M_cell box constraints)
     Σ X[c] = n                                    (normalization)
Variables: X[0..M_cell-1], δ  →  M_cell+1 variables
```

**Custom primal-dual IPM:**

The IPM normal equations for this LP have the structure:
```
N_eff × Δλ = rhs
where N_eff = A × D_X × Aᵀ + δ_component
```
D_X is a diagonal matrix derived from interior-point barrier weights on X[c]. The N_eff matrix is n_cats_total × n_cats_total — computed using `compute_normal_equations` from §2.

IPM iterations:
1. Compute D_X from current (X, δ, μ) primal-dual variables: O(M_cell)
2. Compute N_eff = A × D_X × Aᵀ: O(M_cell × K)
3. Cholesky factor N_eff (836×836): O(836³/3)
4. Solve for (ΔX, Δδ, Δμ) via back-substitution: O(836²)
5. Line search (Mehrotra predictor-corrector): O(M_cell)

**Total per iteration:** O(M_cell × K). ~25–40 IPM iterations. Total: O(40 × M_cell × K).

For stepstone: 40 × 6000 × 9 = 2.16M ops → sub-millisecond.

**Convergence:** to true L∞ minimum. `res.objective = δ*` (achievable L∞ error). `res.objective_metric = MAX_ERR`.

---

### §6 method="greg" — Newton QP for chi2 (GREG calibration)

**Minimizes:** Σ_c (X[c] − X_init[c])² / X_init[c] = Σ_i (w_i/d_i − 1)² d_i subject to hard capacity.

**Unconstrained optimum (no active bounds):**

KKT conditions → closed form for λ ∈ ℝ^{n_cats_total}:
```
N × λ = b   where  N = A × diag(X_init) × Aᵀ,  b = T·n − A·X_init
X*[c] = X_init[c] × (1 + Aᵀ·λ evaluated at cell c)
       = X_init[c] × (1 + Σ_k λ_{(k, group_of_c_in_k)})
```

One linear solve in N gives exact solution when no bounds are active. Same N as §2 with D = X_init.

**With active bounds (projected Newton):**
Active set method: identify cells at bounds, solve reduced Newton system. 3–5 Newton iterations empirically (survey::grake does this at obs-level). Cell-table version: O(K² × M_cell) vs O(K² × n).

**Speedup vs obs-level raking:** M_cell/n ratio. stepstone: 263×.

**Convergence:** `res.objective = chi2*`. `res.objective_metric = CHI2`.

**Reference:** Deville-Särnäl 1992; survey::grake Newton implementation.

---

### §7 method="grake" — LP for normalized Chebyshev

**Minimizes:** δ = max_{k,j} |S_kj − T_kj·n| / (1 + T_kj·n) subject to hard capacity.

Identical LP structure to §5 but with per-(k,j) scaling of the δ constraint:
```
Σ_{c∈(k,j)} X[c] − δ·(1 + T_kj·n) ≤  T_kj·n    ∀k,j
−Σ_{c∈(k,j)} X[c] − δ·(1 + T_kj·n) ≤ −T_kj·n    ∀k,j
```

Same IPM solver as §5 with rescaled right-hand sides. `res.objective_metric = GRAKE_NORM`.

Matches survey::grake's convergence criterion exactly at the fixed point.

---

### §8 method="raking" migration to cell-table

**Algorithm:** unchanged cyclic IPF multiplies cell masses. Same as current obs-level raking but indexed by cell.

**Per-obs sub-weights (non-uniform d_i):** Within a cell, d_i vary. Cell multiplier M[c] is computed at cell level; obs weights recovered at exit as `w_i = d_i × M[c]`. Exact: no within-cell information lost.

**Cell-level IPF iteration:**
```
For k = 1..K, j = 1..cat_counts[k]:
  X_kj = Σ_{c∈bucket(k,j)} X[c]          ← O(cells_in_bucket)
  scale = T_kj × n / X_kj
  X[c] *= scale  ∀c ∈ bucket(k,j)         ← uniform multiplier per bucket
```

At exit: sum constraint and margin constraints satisfied to the same level as obs-level raking. Capacity bounds enforced via Dykstra on X[c] (same as current raking but per-cell).

**Complexity:** O(K × M_cell) per iteration. stepstone speedup: 263× → expected ~0.2s vs current ~50s.

**Convergence:** same metric framework as other methods. Default `kl+improvement` (raking IS a KL minimizer for feasible problems without bounds; with bounds it's approximate KL).

---

### §9 method="ieppa" — default metric change

Change default from `max_err+improvement+0.001` to `kl+improvement+0.001`.

Rationale: iEPPA is a Sinkhorn-type algorithm — it minimizes KL divergence from the initial distribution. Reporting convergence via `kl+improvement` is semantically consistent with the algorithm's actual objective. `max_err+improvement` was chosen as a proxy for calibration quality, but KL is more directly aligned with what iEPPA computes.

**Backward compat:** `convergence=list(absolute=1e-6)` still works as before.

---

## §10 File Structure

| File | Role | New/Modified |
|---|---|---|
| `src/calib_linalg.hpp` | Normal equations, Cholesky: shared by greg/chebyshev/grake | NEW |
| `src/sinkhorn.cpp` + `.hpp` | KL Bregman Dykstra solver | NEW |
| `src/chebyshev.cpp` + `.hpp` | LP IPM for L∞ | NEW |
| `src/greg.cpp` + `.hpp` | Newton QP for chi2 | NEW |
| `src/grake.cpp` + `.hpp` | LP IPM for grake_norm | NEW |
| `src/raking.cpp` | Migrate to cell-table | MODIFIED |
| `src/ieppa.cpp` | Change default metric to kl | MODIFIED |
| `src/c_api.cpp` | Wire new methods; new result fields | MODIFIED |
| `src/r_bridge.cpp` | New method SEXPs | MODIFIED |
| `src/leafblower.h` | New method constants; objective/objective_metric | MODIFIED |
| `R/harvest.R` | method= enum, default update | MODIFIED |
| `python/leafblower/_harvest.py` | Mirror | MODIFIED |

---

## §11 Acceptance Criteria

**A1 (sinkhorn KL minimum):** On stepstone-fulldata, `method="sinkhorn"` KL at convergence ≤ KL from `method="ieppa"` at best_iter. `best_iter == last_iter` (monotone, no overshoot).

**A2 (chebyshev L∞ minimum):** On stepstone-fulldata, `method="chebyshev"` max_err ≤ autumn's max_err (autumn minimizes approximately; chebyshev minimizes exactly subject to hard bounds).

**A3 (greg chi2 minimum):** On smooth synthetic, `method="greg"` chi2 ≤ all other methods' chi2 at convergence.

**A4 (grake exact):** `method="grake"` grake_norm at convergence matches `survey::calibrate(epsilon=1e-10)` within 1%.

**A5 (raking speedup):** `method="raking"` on stepstone-fulldata < 1s (vs current ~50s). Same max_err as current obs-level raking (same fixed point).

**A6 (hard bounds all methods):** All methods: `all(weights >= min_weight - 1e-10) && all(weights <= max_weight + 1e-10)` at exit.

**A7 (best_iter == last_iter):** For sinkhorn, chebyshev, greg, grake: `result$best_iter == result$iterations` (monotone algorithms don't overshoot their optimum).

**A8 (CRAN gate):** devtools::test() FAIL 0; R CMD check --as-cran 0 ERROR 0 WARNING.

---

## §12 Resolved Design Questions

| Question | Resolution |
|---|---|
| Cell-table for raking non-uniform weights | `d_i × M[c]` exit expansion — exact |
| LP solver dependency | Custom IPM in C++, no external deps |
| iEPPA default metric | Change to kl+improvement |
| sinkhorn convergence criteria | All 6 metrics via improvement/plateau; kl default |
| Water-fill replacement | Log-domain Dykstra (multiplicative corrections) |
| Normal equations kernel | Shared calib_linalg.hpp |
| L1_weight metric for raking/sinkhorn | Applicable; same exit expansion gives l1_weight_change |
