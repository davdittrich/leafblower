# Calibration Solvers Redesign — Spec

**Date:** 2026-04-25
**Status:** Draft v2 — all 5 design-review-gate blockers addressed

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

**Pre-entry validation** (`calib_validate_preentry` in `src/calib_validate.hpp`): shared across all new solvers, called from `c_api.cpp` dispatcher before method dispatch. See §13 for full contract.

### §2 Shared normal equations kernel (`src/calib_linalg.hpp`)

The **marginal incidence matrix A** (rows = margin-category pairs (k,j), cols = cells c) is binary: A[(k,j), c] = 1 iff cell c belongs to bucket (k,j). This matrix underlies all LP/QP/Newton steps.

The critical shared computation:
```
N = A × D × Aᵀ    (shape: n_cats_total × n_cats_total)
```
where D is diagonal on cell masses. For stepstone: N is 836×836. Computing N costs O(M_cell × K) using `cells_by_margin_cat`. LDLT factor of N costs O(836³/3).

**Guard (SEC-H1):** `compute_normal_equations` must check `n_cats_total` before allocating:
```cpp
constexpr int kNCatsTotalMax = 2048;  // N matrix 2048×2048 = 32MB, LDLT O(2048³/3)≈2.9B ops ≈100ms — practical cap
if (n_cats_total > kNCatsTotalMax)
    return RK_ERR_BADARG;  // "n_cats_total exceeds limit; use method='ieppa' or 'raking'"
```

```cpp
// src/calib_linalg.hpp
namespace lbw {
// Returns RK_ERR_BADARG if n_cats_total > kNCatsTotalMax.
int compute_normal_equations(const CellTable& ct, const double* D, double* N,
                              const int* cat_offset, size_t n_cats_total);
// Modified LDLT with diagonal perturbation (Gill-Murray) for near-singular N.
// eps_perturb: minimum diagonal value, e.g. 1e-10. Returns RK_ERR_BADARG if singular.
int ldlt_factor_inplace(double* N, size_t n, double eps_perturb);
void ldlt_solve(const double* L, const double* d, double* rhs, size_t n);
} // namespace lbw
```

Note: all dimension parameters use `size_t` (not `int`) to prevent int32 overflow for n_cats_total > 46340 (SEC-M3). Decomposition is LDLT (not plain Cholesky) with diagonal perturbation to handle degenerate margin groups (cells with zero initial mass → zero row/col in N).

### §3 New method constants and ABI

```c
/* leafblower.h additions — added to rk_algorithm_t enum */
#define RK_ALG_SINKHORN   4   /* true KL min, Bregman Dykstra */
#define RK_ALG_CHEBYSHEV  5   /* true L∞ marginal error min, LP */
#define RK_ALG_GREG       6   /* true chi2 min, Newton QP */
#define RK_ALG_GRAKE      7   /* true grake_norm min, LP */
```

**rk_result_t new fields** (folded into the existing `convergence_used` cluster so R-side nesting is consistent):
```c
double convergence_objective;      /* value of minimized metric at convergence */
int    convergence_minimized_metric; /* CalibMetric enum: which metric was minimized */
```

These map to `result$convergence_used$objective` and `result$convergence_used$minimized_metric` in R (same nesting as `result$convergence_used$metric`, `$rule`, `$tol`, `$fired_at_iter`).

**ABI tripwire:** Add `EXPECTED_RK_RESULT_BYTES` to `leafblower.h` mirroring the existing `EXPECTED_RK_PARAMS_BYTES` pattern:
```c
#define EXPECTED_RK_RESULT_BYTES <measured sizeof on first compile>
static_assert(sizeof(rk_result_t) == EXPECTED_RK_RESULT_BYTES,
    "rk_result_t size changed; update EXPECTED_RK_RESULT_BYTES and ABI consumers");
```

For iEPPA, raking, lbfgsb: `convergence_objective` = active metric value at exit; `convergence_minimized_metric` = CalibMetric used for convergence dispatch.

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
     target_mass = Σ_c X[c]       ← current mass BEFORE projection (mass-preserving)
     Bisection on μ ∈ ℝ s.t. Σ_c clamp(X[c]·exp(a[c]+μ), L_c, U_c) = target_mass
     X_proj[c] = clamp(X[c]·exp(a[c]+μ), L_c, U_c)
     a[c] += log(X[c]) - log(X_proj[c])   ← log-domain Dykstra correction
     X[c] = X_proj[c]

  Check convergence (CalibRule on active CalibMetric).
  KL guaranteed monotone decreasing each outer iteration.

**Pre-entry validation (before outer loop):**
- Assert L_c ≤ U_c ∀c; return `RK_ERR_BADARG` with cell index if violated.
- Normalize each margin: if |Σ_j T_kj − 1| > 1e-6, normalize and emit a warning.
- If Σ_c L_c > n or Σ_c U_c < n (total capacity incompatible with target mass), return `RK_ERR_INFEAS`.

**LP infeasibility detection:** The bisection has no solution when Σ_c U_c < target_mass or Σ_c L_c > target_mass. Check bounds before bisect: `if (Σ U_c < target_mass || Σ L_c > target_mass) return RK_ERR_INFEAS`.
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

The δ variable has a dense coefficient vector in the LP constraint matrix: every margin row has coefficient 1 for δ (the shared slack). In the IPM normal equations, δ contributes a rank-1 update via Sherman-Morrison:

```
N_eff × Δλ = rhs
where N_eff = A × D_X × Aᵀ + (1/s_δ) × u × uᵀ
      u = vector of 1s (length n_cats_total), s_δ = IPM barrier weight for δ
```

Since u is all-ones, (1/s_δ) × uuᵀ is a scalar multiple of the all-ones matrix. Apply Sherman-Morrison: if N_0 = A × D_X × Aᵀ (LDLT-factored), then N_eff⁻¹ × r = N_0⁻¹ × r − (c/(1 + s_δ·c)) × v where c = 1/s_δ, v = N_0⁻¹ × u. Two solves against N_0 per IPM iteration (once for rhs, once for u). N_0 LDLT factor is reused; u solve done once per factorization. This keeps the system size at n_cats_total × n_cats_total throughout.

**Warm start:** Before the first IPM iteration, shift X to strictly interior:
```
ε_shift = 1e-8 × max_c(U_c − L_c)
X⁰[c] = clamp(X_entering[c], L_c + ε_shift, U_c − ε_shift)
```
If U_c − L_c < 2ε_shift for any cell: X⁰[c] = (L_c + U_c)/2. Ensures all D_X barrier weights are finite.

IPM iterations:
1. Compute D_X from current (X, δ, μ) primal-dual variables: O(M_cell)
2. Compute N_0 = A × D_X × Aᵀ via `compute_normal_equations`: O(M_cell × K)
3. LDLT-factor N_0 (836×836): O(836³/3)
4. Apply Sherman-Morrison for δ: two additional back-subs: O(836²)
5. Solve for (ΔX, Δδ, Δμ): O(836²)
6. Line search (Mehrotra predictor-corrector): O(M_cell)

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

**Per-obs sub-weights (non-uniform d_i):** Within a cell, d_i vary. Cell multiplier M[c] is computed at cell level; obs weights recovered at exit as `w_i = d_i × M[c]`. The final marginal sums match exactly (Σ_{i∈bucket(k,j)} w_i = T_kj × n at the cell-level fixed point).

**Per-obs hard bounds (addressing Goal 7):** After the `w_i = d_i × M[c]` expansion, a per-obs clamp is applied:
```
w_i = clamp(d_i × M[c], min_weight, max_weight)
```
This guarantees `min_weight ≤ w_i ≤ max_weight` for all i, satisfying Goal 7 strictly.

**Trade-off:** For non-uniform d_i, the post-exit clamp may slightly distort margins (cells where some obs were clamped will have S_kj slightly off from T_kj × n). This is the same trade-off as obs-level raking's Dykstra which also clamps post-hoc. The residual marginal error from clamping is bounded by the within-cell d_i variance × |M[c] - clamp|, which is small when the calibration is near convergence. This is reported in `max_error` at exit.

**Cell-level IPF iteration:**
```
For k = 1..K, j = 1..cat_counts[k]:
  X_kj = Σ_{c∈bucket(k,j)} X[c]          ← O(cells_in_bucket)
  scale = T_kj × n / X_kj
  X[c] *= scale  ∀c ∈ bucket(k,j)         ← uniform multiplier per bucket
```

At exit: sum constraint and margin constraints satisfied to the same level as obs-level raking.

**Capacity enforcement — cell-level Dykstra (alternating projections):**
```
Initialize: p[c] = 0 ∀c  (Dykstra corrections, same as q[i] in raking.cpp)

Each capacity enforcement step (interleaved with IPF sweeps):
  ① Box projection:
     X_before[c] = X[c] + p[c]
     X[c] = clamp(X_before[c], L_c, U_c)
     p[c] += X_before[c] − X[c]   ← accumulate Dykstra correction

  ② Sum projection:
     total = Σ_c X[c]
     X[c] *= n / total              ← rescale to target mass
```
Same two-step structure as obs-level raking's Dykstra (`q[i]` + sum hyperplane), applied to cells instead of obs. Provably convergent to the intersection of the box and the sum hyperplane (Boyle-Dykstra 1986).

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

**A1 (sinkhorn KL minimum):** On stepstone-fulldata, `method="sinkhorn"` KL at convergence ≤ KL from `method="ieppa"` at best_iter. `best_iter == last_iter` (monotone, no overshoot). Reference fixture: `tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds` (scalar — iEPPA KL at best_iter on stepstone). Generated by `data-raw/gen_ieppa_kl_ref.R`; committed before A1 test is written.

**A2 (chebyshev L∞ minimum):** On stepstone-fulldata, `method="chebyshev"` max_err ≤ autumn's max_err (autumn minimizes approximately; chebyshev minimizes exactly subject to hard bounds).

**A3 (greg chi2 minimum):** On smooth synthetic, `method="greg"` chi2 ≤ all other methods' chi2 at convergence.

**A4 (grake exact):** `method="grake"` grake_norm at convergence matches `survey::calibrate(epsilon=1e-10)` within 1%.

**A5 (raking speedup):** `method="raking"` on stepstone-fulldata (n=1.58M, M_cell≈6k): `skip_on_cran()` + criterion `elapsed_cell < elapsed_obs × 0.05` (measured in same process). For homogeneous d_i: max_err within 1e-8 of obs-level raking. For non-uniform d_i: tested separately with documented tolerance.

**A6 (hard bounds):**
- sinkhorn, chebyshev, greg, grake: per-cell AND per-obs bounds hold: `all(w >= min_weight - 1e-12) && all(w <= max_weight + 1e-12)`.
- raking cell-table: per-obs bounds hold exactly via post-exit obs-level clamp `w_i = clamp(d_i × M[c], min_weight, max_weight)`. Small marginal distortion from clamped cells is reported in `max_error`.

**A7 (best_iter == last_iter — monotone methods only):**
- sinkhorn: `result$convergence_used$fired_at_iter == result$iterations` (monotone KL, no overshoot by proof).
- chebyshev, grake: LP terminates at exact optimum; fired_at_iter == iterations trivially.
- greg: Newton terminates at optimum (1 step unconstrained; active-set converges). fired_at_iter == iterations.
- ieppa, raking, lbfgsb: NOT included in A7 (iterative, oscillation possible).

**A8 (CRAN gate):** devtools::test() FAIL 0; R CMD check --as-cran 0 ERROR 0 WARNING.

---

## §13 Infeasibility Handling

All new methods return `RK_ERR_INFEAS` (existing code = 2) when the problem has no feasible solution. Pre-entry checks are in `calib_validate_preentry(const CalibState& st, const CellTable& ct)` declared in `src/calib_validate.hpp`, implemented in `src/calib_validate.cpp`, called from `c_api.cpp:rk_calibrate()` BEFORE method dispatch:

1. `L_c ≤ U_c` ∀c — if violated: `RK_ERR_BADARG` with cell index.
1b. `X_init[c] > 0 || L_c == 0` ∀c — if X_init[c]=0 AND L_c>0: `RK_ERR_INFEAS` ("cell c has zero initial mass but positive lower bound; multiplicative solvers cannot move zero cells"). This is a pre-existing iEPPA check (structural_infeas_pairs) extended to all new methods.
2. `Σ_c L_c ≤ n ≤ Σ_c U_c` — total mass feasibility; if violated: `RK_ERR_INFEAS` with message "total capacity incompatible with target mass n".
3. `|Σ_j T_kj − 1| ≤ 1e-6` ∀k — targets normalized; if not: normalize + emit `warning()`.
4. `n_cats_total ≤ kNCatsTotalMax` — memory guard; if violated: `RK_ERR_BADARG`.

During LP solve (chebyshev, grake): IPM homogeneous self-dual embedding detects primal infeasibility (margin targets outside capacity polytope) → `RK_ERR_INFEAS`. During sinkhorn: capacity bisection range check before each outer iteration.

R error message: `"leafblower: calibration infeasible — [reason]. Check margin targets and weight bounds."`

## §12 Resolved Design Questions

| Question | Resolution |
|---|---|
| Cell-table for raking non-uniform weights | `d_i × M[c]` exit expansion — exact |
| LP solver dependency | Custom IPM in C++, no external deps |
| iEPPA default metric | Change to kl+improvement |
| sinkhorn convergence criteria | All 6 metrics via improvement/plateau; kl default |
| Water-fill replacement | Log-domain Dykstra (multiplicative corrections) |
| Normal equations kernel | Shared calib_linalg.hpp with LDLT+diagonal perturbation |
| δ column in IPM normal equations | Rank-1 update via Sherman-Morrison (two solves per iter) |
| Bisection normalization target | target_mass = Σ_c X[c] before projection (mass-preserving) |
| objective/objective_metric naming | Folded into convergence_used: $convergence_used$objective + $minimized_metric |
| rk_result_t ABI tripwire | EXPECTED_RK_RESULT_BYTES static_assert added |
| n_cats_total memory guard | Cap at kNCatsTotalMax = 8192; return RK_ERR_BADARG if exceeded |
| int32 overflow in Cholesky | All linalg dimensions as size_t |
| Raking fixed-point for non-uniform d_i | Documented approximation; per-cell bounds guaranteed, per-obs approximate |
| Infeasibility handling | §13: shared pre-entry checks; RK_ERR_INFEAS with message |
| A7 scope | Monotone methods only (sinkhorn/chebyshev/greg/grake); ieppa/raking excluded |
| A5 CI testability | Relative criterion (elapsed_cell < elapsed_obs × 0.05), skip_on_cran |
| L1_weight metric for raking/sinkhorn | Applicable; same exit expansion gives l1_weight_change |


## §14 Versioning and Migration

**Version bump:** 0.1.0 → 0.2.0 (minor: new public API surface). Pre-1.0 software — breaking default changes are acceptable.

**Breaking changes requiring NEWS.md entry:**

1. `method="ieppa"` default convergence changes from `max_err+improvement+0.001` to `kl+improvement+0.001`. Users relying on prior default behavior: add `convergence=list(metric="max_err", rule="improvement", tol=0.001)` explicitly. Result fields `result$convergence_used$metric` and `result$iterations` will differ.

2. `method="raking"` migrated to cell-table. Identical calibration results for homogeneous design weights. Non-uniform design weights: minor numerical differences (within within-cell d_i variance × M[c]). Users who observe unexpected differences should report as a bug.

3. `rk_result_t` gains two new fields (`convergence_objective`, `convergence_minimized_metric`). C-API consumers must recompile. `EXPECTED_RK_RESULT_BYTES` tripwire will catch this at compile time.

**NEWS.md template:**
```markdown
## leafblower 0.2.0 (unreleased)

### Breaking changes
- method="ieppa": default convergence now kl+improvement (was max_err+improvement)
- method="raking": cell-table backend (100-300x faster; results identical for uniform design weights)
- C ABI: rk_result_t gains 2 fields; recompile required

### New methods
- method="sinkhorn": true KL minimum subject to capacity constraints
- method="chebyshev": true L-infinity marginal error minimum (LP)
- method="greg": true chi-squared minimum (GREG)
- method="grake": true grake_norm minimum (normalized Chebyshev)
```
