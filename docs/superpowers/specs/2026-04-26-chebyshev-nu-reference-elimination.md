# Chebyshev ν Fix: Reference Category Elimination

**Date**: 2026-04-26
**Status**: Approved for implementation
**Tickets**: leafblower-8i6f (root cause), leafblower-2f1p (original nu-fix filed during Plan E)
**File**: `src/chebyshev.cpp` only

---

## Problem

The chebyshev/grake LP IPM minimizes δ = max_m |S[m]/W − T[m]| via interior-point. The normalization constraint Σ X[c] = n_d is NOT in the Newton system because the normalization dual ν is **algebraically degenerate**: for any margin k with Σ_j T[k][j] = 1 (standard case), the normalization equation is a linear combination of the margin equations, so `schur_nu = D_ν − e^T · N_eff^{-1} · e ≈ 0`.

Without ν, the Newton step is not sum-preserving. Σ X[c] drifts from n_d across iterations. On K≥4 complex margin problems (e.g. stepstone K=9), W → 0 after ~500 iterations — all-zero metrics, misleading NOCONV.

The renorm workaround (cheb-renorm branch) prevents collapse but does not produce meaningful convergence because the Newton direction is still wrong.

---

## Root Cause

For margin k with J_k categories (each cell belongs to exactly one category per k), the margin constraint rows satisfy:

```
Σ_{j=0}^{J_k-1} A[cat_offset[k]+j, c] = 1  for all cells c
```

This means: Σ_{j in margin k} e[m] = Σ_c D_eff[c] = D_ν for each k.

Equivalently: the vector `e` (used in the ν second Sherman-Morrison) is in the span of the rows of N_0 = A · D_eff · A^T. Therefore `N_eff^{-1} · e` ≈ N_0^{-1} · e, and `e^T · N_0^{-1} · e = D_ν`, giving `schur_nu ≈ 0`.

---

## Fix: Reference Category Elimination

For each margin k, drop the **last category** (index J_k − 1) from the Newton system. This removes the linearly dependent row, making the reduced system full-rank. The ν second Sherman-Morrison then works correctly.

### Key Properties

- All K · Σ J_k LP slack variables (s_up/s_dn) and duals (y_up/y_dn) remain at full size — IPM tracking unaffected.
- Reference category calibration quality is enforced by difference: `S[k, j_ref] = W − Σ_{j≠j_ref} S[k,j]`, automatically correct when ν forces Σ X[c] = n_d.
- Reference multiplier `dlambda_ref[k] = 0` per Newton step (implicitly). This is the standard identification in regression/ANOVA.

---

## Edge Case: Single-Category Margins

If `cat_counts[k] == 1` for any margin k, that margin has only one category. Dropping it as reference removes all Newton correction for that margin, leaving a 0-contribution row in the reduced system (or a 0×0 system for K=1). 

**Rule**: When `cat_counts[k] < 2`, do NOT eliminate a reference for that margin. Treat all categories as non-reference. The margin is already trivially satisfied (only one category, T[k][0]=1, so X must sum to W) — ν alone handles the normalization.

This adjusts `nct_red` dynamically:
```
nct_red = nct - count(k : cat_counts[k] >= 2)
```

---

## Data Structures

```cpp
// Computed once per IPM call, after cat_offset is built:
// nct_red = nct - number of margins with >=2 categories (reference elimination only there)
int nct_red_count = 0;
for (int k = 0; k < st.K; k++)
    if (st.cat_counts[k] >= 2) nct_red_count++;
const int nct_red = nct - nct_red_count;

std::vector<bool> is_ref(nct, false);          // true for last category of each margin k (only when cat_counts[k]>=2)
std::vector<int>  full_to_red(nct, -1);        // -1 for reference margins
std::vector<int>  red_to_full(nct_red);        // reduced → full index

int nr = 0;
for (int k = 0; k < st.K; k++) {
    for (int j = 0; j < st.cat_counts[k]; j++) {
        int m = cat_offset[k] + j;
        if (st.cat_counts[k] >= 2 && j == st.cat_counts[k] - 1) {
            is_ref[m] = true;  // reference: last category (only for multi-category margins)
        } else {
            full_to_red[m] = nr;
            red_to_full[nr++] = m;
        }
    }
}
// Invariant: nr == nct_red
```

---

## Newton System Changes

### Normal equations (reduced)

Replace current `N0` (nct × nct) with `N_red` (nct_red × nct_red):

```cpp
std::vector<double> N_red((size_t)nct_red * (size_t)nct_red);
```

`compute_normal_equations` is called with the reduced index mapping. Two options:
- **Option A (preferred)**: Add a `compute_normal_equations_reduced` variant that accepts `full_to_red` and skips is_ref margins. Avoids extra allocation.
- **Option B**: Build N_red from N0 by selecting non-reference rows/cols post-hoc.

Use Option A: O(K · M_cell) pass, same complexity as current.

### Reduced vectors

```cpp
// Replaces current w_sol[nct] — reduced form only
std::vector<double> w_sol_red(nct_red);   // N_red^{-1} · rhs_red
std::vector<double> u_red(nct_red);       // w_kj for non-reference m
std::vector<double> v_red(nct_red);       // N_eff_red^{-1} · u_red  (first SM)
std::vector<double> e_red(nct_red);       // Σ_{c∈m} D_eff[c]  (ν SM vector)
std::vector<double> w_e_red(nct_red);     // N_eff_red^{-1} · e_red  (second SM)
```

### RHS (reduced)

```cpp
for (int nr_idx = 0; nr_idx < nct_red; nr_idx++) {
    int m = red_to_full[nr_idx];
    rhs_red[nr_idx] = -(S[m] - T_flat[m]*W)
                      + D_marg[m] * (rmu_up[m]/s_up[m] - rmu_dn[m]/s_dn[m]);
}
```

### Three back-solves (same as nu-fix plan, but on reduced system)

```
1. LDLT(N_red)  (size nct_red × nct_red)
2. w_sol_red = N_red^{-1} · rhs_red
3. v_red = N_red^{-1} · u_red     ← first SM (δ)
4. e_red[nr] = Σ_{c∈red_to_full[nr]} D_eff[c]
5. w_e_red = N_red^{-1} · e_red   ← third back-solve (ν)
6. Apply first SM correction to w_sol_red and w_e_red
```

### ν computation (non-degenerate)

```cpp
double eTw_e = 0.0, eTdlambda = 0.0;
for (int nr = 0; nr < nct_red; nr++) {
    eTw_e     += e_red[nr] * w_e_red[nr];
    eTdlambda += e_red[nr] * dlambda_red[nr];  // dlambda_red = SM-corrected w_sol_red
}
const double D_nu       = /* Σ_c D_eff[c] */;
const double schur_nu   = D_nu - eTw_e;
// schur_nu > 0 guaranteed by reference elimination (key invariant).
// Verbose diagnostic: if (st.verbose >= 2 && iter == 1) log schur_nu value.
// Tightened guard: 1e-8 (schur_nu is O(D_ν) = O(n_d) >> 1e-8 for any realistic problem).
const double r_nu       = W - n_d;             // normalization residual
const double dnu        = (schur_nu > 1e-8)
                        ? (-r_nu - eTdlambda) / schur_nu : 0.0;

// Correct dlambda_red
for (int nr = 0; nr < nct_red; nr++) dlambda_red[nr] -= dnu * w_e_red[nr];
```

### dX reconstruction

```cpp
for (int c = 0; c < ct.M_cell; c++) {
    double sum_dlam = 0.0;
    for (int k = 0; k < st.K; k++) {
        int g = ct.g_per_cell[k][c];
        if (g < 0 || g >= st.cat_counts[k]) continue;
        int m  = cat_offset[k] + g;
        int nr = full_to_red[m];
        if (nr >= 0) sum_dlam += dlambda_red[nr];  // reference: nr == -1, skip
    }
    dX[c] = D_eff[c] * (sum_dlam + dnu);  // Δν applies to ALL cells
}
```

### dS update (reduced)

```cpp
std::fill(delta_S.begin(), delta_S.end(), 0.0);  // delta_S still full size nct
for (int c = 0; c < ct.M_cell; c++) {
    for (int k = 0; k < st.K; k++) {
        int g = ct.g_per_cell[k][c];
        if (g < 0 || g >= st.cat_counts[k]) continue;
        int m = cat_offset[k] + g;
        delta_S[m] += dX[c];   // full nct — reference margins also get dS
    }
}
```

The reference margin `delta_S[m_ref]` is computed but not used in the Newton solve (it contributes to slack updates). This is correct — reference margin slacks s_up/s_dn are still updated each iteration.

### w_dot_dlambda (recompute after ν correction)

```cpp
double w_dot_dlambda = 0.0;
for (int nr = 0; nr < nct_red; nr++)
    w_dot_dlambda += w_kj[red_to_full[nr]] * dlambda_red[nr];
// d_delta = -(w_dot_dlambda + ...) / ...  [existing formula]
```

---

## dlambda storage

`dlambda_red` replaces `dlambda` in the Newton step. The existing `dlambda` (size nct) can be repurposed as a full-size vector used only for slack/dual updates:

```cpp
// After ν correction, expand to full for dual updates:
std::fill(dlambda_full.begin(), dlambda_full.end(), 0.0);
for (int nr = 0; nr < nct_red; nr++) dlambda_full[red_to_full[nr]] = dlambda_red[nr];
// dlambda_full[is_ref[m]] remains 0.0 (reference: no Newton contribution)
```

Then the existing dY_up/dY_dn computation uses `dlambda_full`.

---

## Slack/dual updates — unchanged

All nct s_up/s_dn, y_up/y_dn updated with full `delta_S` and `dlambda_full`. No change to complementarity tracking, line search, mu computation.

---

## compute_normal_equations_reduced

New variant (or extend existing) in `src/calib_linalg.hpp`:

```cpp
// Compute N = A_red × diag(D) × A_red^T where A_red excludes reference margins.
// full_to_red: nct-length array, -1 for reference margins.
// nct_red: number of non-reference constraints.
int compute_normal_equations_reduced(
    const CellTable& ct,
    const double* D,
    double* N,               // nct_red × nct_red output
    const int* cat_offset,   // full cat_offset, size K
    int K,
    size_t nct_red,
    const int* full_to_red   // full → reduced index (-1 = reference)
) noexcept;
```

Or: pass `is_ref` bool array. Either works.

---

## kNCatsTotalMax

Still applies to full `nct`. Reference elimination reduces the Newton system but doesn't change the LP constraint count. Guard remains at `nct > kNCatsTotalMax → RK_ERR_BADARG`.

---

## Verification

### In-scope (required for this ticket)

1. **E1/E2 correctness tests** (existing, `test-calibration-solvers.R`) — must stay GREEN after fix. E1: chebyshev max_err ≤ raking max_err + 1e-6. E2: grake grake_norm ≤ raking grake_norm + 1e-6.
2. **schur_nu diagnostic**: On first IPM iteration (`iter == 0`), log `schur_nu` value when `verbose >= 2`. Add `expect_gt(schur_nu_logged, 1e-6)` in a test fixture that captures log output. This verifies the non-degeneracy claim. File: `tests/testthat/test-calib-linalg.R` (new test).
3. **Single-category margin test**: Add test with `a = factor(c("x","x","x"))` (1 category). Verify no crash and sensible output (chevyshev with K=1, J=1 should trivially converge or return INFEAS, not crash).

### Deferred (follow-on tickets)

4. **Stepstone K=9 convergence**: chebyshev/grake should converge (status=0). Deferred to `leafblower-yned` (spike validation ticket).
5. **kk1204 K=20 convergence**: run after fix. Deferred to P1 gate `leafblower-kk1.20.4`.

---

## Files Changed

- `src/chebyshev.cpp` — main changes (Newton system, reduced bookkeeping, ν fix)
- `src/calib_linalg.hpp` + `src/calib_linalg.cpp` — add `compute_normal_equations_reduced`

---

## Out of Scope

- Renorm workaround in `cheb-renorm` branch — superseded by this fix; discard that branch after this ships.
- Python/R API changes — none required.
- Other solvers — raking/sinkhorn/greg unaffected.
