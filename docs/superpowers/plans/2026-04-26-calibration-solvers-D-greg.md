# Calibration Solvers — Plan D: calib_linalg + method="greg" (Newton QP for chi2)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the shared linear algebra kernel (`calib_linalg.cpp`) and `method="greg"` — the true chi2 minimum (GREG calibration, Deville-Särnäl 1992). Removes the `RK_ALG_GREG` stub and makes tests GREEN.

**Architecture:** `calib_linalg.cpp` implements `compute_normal_equations` + `ldlt_factor_inplace` + `ldlt_solve`. `greg.cpp` + `greg.hpp` implement one-shot Newton QP with active-set for bounds. chebyshev/grake deferred to Plan E (LP-IPM).

**Spec:** `docs/superpowers/specs/2026-04-25-calibration-solvers-design.md` §2, §6

**Baseline:** FAIL 0 | PASS 361 | SKIP 5

---

## File Structure

| File | Action |
|---|---|
| `src/calib_linalg.cpp` | New — compute_normal_equations + ldlt_factor_inplace + ldlt_solve |
| `src/calib_linalg.hpp` | Update NOTE comments (remove "Plan D" forward refs) |
| `src/greg.hpp` | New — GregResult struct + greg_solve declaration |
| `src/greg.cpp` | New — Newton QP implementation |
| `src/c_api.cpp` | Replace RK_ALG_GREG stub with greg_solve dispatch |
| `src/r_bridge.cpp` | Remove "greg" from stub block; add greg dispatch |
| `src/Makevars` + `Makevars.in` | Add calib_linalg.cpp + greg.cpp |
| `tests/testthat/test-calibration-solvers.R` | Add D1 test (greg chi2 ≤ other methods) |

---

## Key algorithms

### compute_normal_equations

N is `n_cats_total × n_cats_total`, stored row-major. `cat_offset[k]` = starting row index for margin k.

For each cell c (with g_per_cell[k][c] = j_k for each margin k):
- Cell c contributes D[c] to N[row1, row2] for every pair (k1,j1), (k2,j2) where c ∈ bucket(k1,j1) AND c ∈ bucket(k2,j2)
- Equivalently: for each cell c, for each pair of margins (k1,k2), add D[c] to N[cat_offset[k1]+g_k1[c], cat_offset[k2]+g_k2[c]]

O(M_cell × K²) — for stepstone (K=9, M_cell=6k): 486k ops.

```cpp
int compute_normal_equations(const CellTable& ct, const double* D, double* N,
                              const int* cat_offset, size_t n_cats_total)
{
    if (n_cats_total > (size_t)kNCatsTotalMax) return RK_ERR_BADARG;
    std::fill(N, N + n_cats_total * n_cats_total, 0.0);
    for (int c = 0; c < ct.M_cell; c++) {
        if (D[c] <= 0.0) continue;
        for (int k1 = 0; k1 < ct.K_stored; k1++) {  // NOTE: check actual field name
            int j1 = ct.g_per_cell[k1][c];
            if (j1 < 0) continue;
            size_t row = (size_t)cat_offset[k1] + (size_t)j1;
            for (int k2 = 0; k2 < ct.K_stored; k2++) {
                int j2 = ct.g_per_cell[k2][c];
                if (j2 < 0) continue;
                size_t col = (size_t)cat_offset[k2] + (size_t)j2;
                N[row * n_cats_total + col] += D[c];
            }
        }
    }
    return RK_OK;
}
```

**IMPORTANT before implementing:** check exact field name for K in CellTable:
```bash
grep -n "K_stored\|int K\b\|\.K\b" src/cell_table.hpp src/cell_table.cpp | head -10
```
CellTable may not store K — the outer `K` comes from CalibState. Use `st.K` passed to greg_solve, not a CellTable field.

### ldlt_factor_inplace

Modified LDLT with Gill-Murray diagonal perturbation. Overwrites A with L (lower triangular, unit diagonal) and d (diagonal, in-place on A's diagonal):

```cpp
int ldlt_factor_inplace(double* A, size_t n, double eps_perturb)
{
    if (n > (size_t)kNCatsTotalMax) return RK_ERR_BADARG;
    for (size_t j = 0; j < n; j++) {
        // d[j] = A[j][j] - sum_{k<j} d[k] * L[j][k]^2
        double d_j = A[j * n + j];
        for (size_t k = 0; k < j; k++) {
            double l_jk = A[j * n + k];
            double d_k  = A[k * n + k];
            d_j -= d_k * l_jk * l_jk;
        }
        d_j = std::max(d_j, eps_perturb);  // Gill-Murray perturbation
        A[j * n + j] = d_j;
        // L[i][j] = (A[i][j] - sum_{k<j} d[k] * L[i][k] * L[j][k]) / d[j]  for i > j
        for (size_t i = j + 1; i < n; i++) {
            double l_ij = A[i * n + j];
            for (size_t k = 0; k < j; k++) {
                double d_k  = A[k * n + k];
                l_ij -= d_k * A[i * n + k] * A[j * n + k];
            }
            A[i * n + j] = l_ij / d_j;
        }
    }
    return RK_OK;
}
```

Storage convention: `A[i*n+j]` (row-major). After factorization:
- Diagonal `A[j*n+j]` holds d[j]
- Lower triangle `A[i*n+j]` (i>j) holds L[i][j]
- Upper triangle unused

### ldlt_solve

Forward-backward substitution: solves LDLᵀx = b, overwrites b with x.

```cpp
void ldlt_solve(const double* A, size_t n, double* b)
{
    // Forward: solve L y = b  (L unit lower triangular)
    for (size_t i = 1; i < n; i++)
        for (size_t j = 0; j < i; j++)
            b[i] -= A[i * n + j] * b[j];
    // Scale: solve D z = y
    for (size_t i = 0; i < n; i++)
        b[i] /= A[i * n + i];
    // Backward: solve Lᵀ x = z  (Lᵀ unit upper triangular)
    for (size_t i = n - 1; i < n; i--) {  // NOTE: i is size_t, loop carefully
        for (size_t j = i + 1; j < n; j++)
            b[i] -= A[j * n + i] * b[j];
    }
}
```

**NOTE:** the backward loop `for (size_t i = n-1; i < n; i--)` works because size_t wraps on underflow — when i=0 and decremented, it becomes SIZE_MAX > n, exiting the loop. This is correct but looks odd. Alternative: use `int i` for the backward loop to avoid confusion.

---

## Task 1 — calib_linalg.cpp

**Ticket:**
```bash
bd create --title="feat(linalg): implement compute_normal_equations + LDLT" --type=feature --priority=1 2>&1 | tail -2
```

### Step 1.1: Create src/calib_linalg.cpp

**Before writing:** read `src/cell_table.hpp` to confirm the CellTable struct field for K. The function receives K as a separate parameter (from CalibState), NOT from CellTable. Confirm: `grep -n "int K\b\|K_stored\|K =" src/cell_table.hpp | head -5`

```cpp
#include "calib_linalg.hpp"
#include "leafblower.h"
#include <algorithm>
#include <cstring>
#include <cmath>

namespace lbw {

int compute_normal_equations(const CellTable& ct,
                              const double* D,
                              double* N,
                              const int* cat_offset,
                              int K,                   // number of margins
                              size_t n_cats_total)
{
    if (n_cats_total > (size_t)kNCatsTotalMax) return RK_ERR_BADARG;
    std::fill(N, N + n_cats_total * n_cats_total, 0.0);
    for (int c = 0; c < ct.M_cell; c++) {
        if (D[c] <= 0.0) continue;
        for (int k1 = 0; k1 < K; k1++) {
            int j1 = ct.g_per_cell[k1][c];
            if (j1 < 0) continue;
            size_t row = (size_t)cat_offset[k1] + (size_t)j1;
            for (int k2 = 0; k2 < K; k2++) {
                int j2 = ct.g_per_cell[k2][c];
                if (j2 < 0) continue;
                size_t col = (size_t)cat_offset[k2] + (size_t)j2;
                N[row * n_cats_total + col] += D[c];
            }
        }
    }
    return RK_OK;
}

int ldlt_factor_inplace(double* A, size_t n, double eps_perturb)
{
    if (n > (size_t)kNCatsTotalMax) return RK_ERR_BADARG;
    for (size_t j = 0; j < n; j++) {
        double d_j = A[j * n + j];
        for (size_t k = 0; k < j; k++) {
            double l_jk = A[j * n + k];
            double d_k  = A[k * n + k];
            d_j -= d_k * l_jk * l_jk;
        }
        d_j = std::max(d_j, eps_perturb);
        A[j * n + j] = d_j;
        for (size_t i = j + 1; i < n; i++) {
            double s = A[i * n + j];
            for (size_t k = 0; k < j; k++)
                s -= A[k * n + k] * A[i * n + k] * A[j * n + k];
            A[i * n + j] = s / d_j;
        }
    }
    return RK_OK;
}

void ldlt_solve(const double* A, size_t n, double* b)
{
    // Forward: L y = b
    for (size_t i = 1; i < n; i++)
        for (size_t j = 0; j < i; j++)
            b[i] -= A[i * n + j] * b[j];
    // Scale: D z = y
    for (size_t i = 0; i < n; i++)
        b[i] /= A[i * n + i];
    // Backward: Lᵀ x = z (careful with size_t underflow)
    for (int i = (int)n - 1; i >= 0; i--)
        for (int j = i + 1; j < (int)n; j++)
            b[i] -= A[(size_t)j * n + (size_t)i] * b[j];
}

} // namespace lbw
```

**NOTE on `ldlt_solve` signature:** the existing `calib_linalg.hpp` declares:
```cpp
void ldlt_solve(const double* L, const double* d_diag, double* b, size_t n);
```
This uses separate `L` and `d_diag` arrays. Since `ldlt_factor_inplace` stores both in-place in the same matrix `A`, the implementation should either:
- Change the `.hpp` signature to `void ldlt_solve(const double* A, size_t n, double* b)` (combined storage), OR
- Keep the separate L/d_diag signature and unpack in greg_solve

**Read the existing .hpp declaration first** and match the signature exactly. If the .hpp says `(L, d_diag, b, n)`, implement with those parameters.

### Step 1.2: Update calib_linalg.hpp — MANDATORY edits

CellTable does NOT store K (verified: `src/cell_table.hpp` exposes M_cell, cell_of, n_per_cell, g_per_cell — no K). Therefore:

**Edit 1 — add `int K` to `compute_normal_equations` declaration:**
```cpp
// OLD:
int compute_normal_equations(const CellTable& ct,
                              const double* D,
                              double* N,
                              const int* cat_offset,
                              size_t n_cats_total);

// NEW:
int compute_normal_equations(const CellTable& ct,
                              const double* D,
                              double* N,
                              const int* cat_offset,
                              int K,
                              size_t n_cats_total);
```

**Edit 2 — rewrite `ldlt_solve` to combined in-place storage (matches `ldlt_factor_inplace`):**
```cpp
// OLD:
void ldlt_solve(const double* L, const double* d_diag, double* b, size_t n);

// NEW — A is the in-place factored matrix from ldlt_factor_inplace:
void ldlt_solve(const double* A, size_t n, double* b);
```

**Edit 3 — remove "NOTE: implementation in Plan D" comments** from all three function docstrings.

### Step 1.3: Add to Makevars AND Makevars.in
```bash
cat src/Makevars
```
Add `calib_linalg.cpp` to `PKG_SOURCES` in BOTH files.

### Step 1.4: Build gate
```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```
Expected: `* DONE (leafblower)`.

### Step 1.5: Commit
```bash
git add src/calib_linalg.cpp src/calib_linalg.hpp src/Makevars src/Makevars.in
git commit -m "feat(linalg): compute_normal_equations + LDLT factorization/solve

compute_normal_equations: O(M_cell*K^2) cell-table normal equations.
ldlt_factor_inplace: Gill-Murray modified LDLT (diagonal ≥ eps_perturb).
ldlt_solve: forward-back substitution on in-place factored matrix.
Used by greg (Plan D) and chebyshev/grake (Plan E)."
bd close <ticket-id>
```

---

## Task 2 — greg.hpp + greg.cpp

**Ticket:**
```bash
bd create --title="feat(greg): Newton QP for chi2 calibration (GREG)" --type=feature --priority=1 2>&1 | tail -2
```

### Step 2.1: Create src/greg.hpp

```cpp
#pragma once
#include "types.hpp"
#include "cell_table.hpp"
#include <vector>
#include <limits>

namespace lbw {

struct GregResult {
    int    status        = RK_ERR_NOCONV;
    int    iterations    = 0;  // Newton iterations (typically 1 unconstrained, ≤5 with bounds)
    double max_error     = 1.0;
    double mean_error    = 0.0;
    double kl            = 0.0;
    double chi2          = 0.0;
    double grake_norm    = 0.0;
    double l1_weight_change = 0.0;
    int    convergence_metric = static_cast<int>(CalibMetric::CHI2);
    int    convergence_rule   = 0;
    double convergence_tol    = 0.0;
    int    convergence_iter   = 1;
    double convergence_objective        = std::numeric_limits<double>::infinity();
    int    convergence_minimized_metric = static_cast<int>(CalibMetric::CHI2);
    double best_error   = std::numeric_limits<double>::infinity();
    int    best_iter    = 1;
    std::vector<double> best_weights;
    int    M_cell       = 0;
    char   message[256] = {};
};

// Newton QP for chi2 calibration (GREG — Deville-Sarnal 1992).
// Minimizes Σ_c (X[c] - X_init[c])^2 / X_init[c] subject to margin constraints + capacity.
// One Newton iteration (exact solution when no bounds active).
// Active-set for bounds: ≤5 reduced Newton iterations.
GregResult greg_solve(CalibState& st);

} // namespace lbw
```

### Step 2.2: Create src/greg.cpp

```cpp
#include "greg.hpp"
#include "calib_linalg.hpp"
#include "calib_dispatch.hpp"
#include "leafblower.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <limits>

namespace lbw {

GregResult greg_solve(CalibState& st) {
    GregResult res;
    res.status = RK_ERR_NOCONV;

    CellTable ct;
    if (build_cell_table(st.n, st.K, st.group_ids, st.cat_counts, st.weights, ct) != 0) {
        res.status = RK_ERR_BADARG;
        std::strncpy(res.message, "build_cell_table failed", sizeof(res.message) - 1);
        return res;
    }
    res.M_cell = ct.M_cell;

    // Initial cell masses and capacity bounds
    std::vector<double> X_init(ct.M_cell, 0.0);
    for (int i = 0; i < st.n; i++) X_init[ct.cell_of[i]] += st.weights[i];

    const double lo = st.min_weight;
    const double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    std::vector<double> L_cell(ct.M_cell), U_cell(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        L_cell[c] = lo * static_cast<double>(ct.n_per_cell[c]);
        U_cell[c] = hi * static_cast<double>(ct.n_per_cell[c]);
    }

    // cat_offset[k] = starting index for margin k in the n_cats_total vector
    std::vector<int> cat_offset(st.K);
    int n_cats_total = 0;
    for (int k = 0; k < st.K; k++) { cat_offset[k] = n_cats_total; n_cats_total += st.cat_counts[k]; }

    if (n_cats_total > kNCatsTotalMax) {
        res.status = RK_ERR_BADARG;
        std::snprintf(res.message, sizeof(res.message),
                      "n_cats_total=%d exceeds limit %d; use method='ieppa'",
                      n_cats_total, kNCatsTotalMax);
        return res;
    }

    // Active-set: cells at lower bound (fixed_lo[c]=true) or upper (fixed_hi[c]=true)
    std::vector<bool> fixed_lo(ct.M_cell, false), fixed_hi(ct.M_cell, false);

    // X = current iterate (starts at X_init)
    std::vector<double> X(X_init);

    const int max_cats = *std::max_element(st.cat_counts, st.cat_counts + st.K);
    std::vector<double> bucket_b(max_cats);  // hoisted: reused across Newton iterations

    static constexpr int kMaxNewtonIters = 10;
    static constexpr double kEps = 1e-10;

    for (int newton_iter = 0; newton_iter < kMaxNewtonIters; newton_iter++) {
        res.iterations = newton_iter + 1;

        // Effective D: X_init[c] for free cells, 0 for fixed cells (excluded from N)
        std::vector<double> D_eff(ct.M_cell, 0.0);
        for (int c = 0; c < ct.M_cell; c++) {
            if (!fixed_lo[c] && !fixed_hi[c] && X_init[c] > kEps)
                D_eff[c] = X_init[c];
        }

        // N = A × diag(D_eff) × Aᵀ  (only free cells contribute)
        std::vector<double> N((size_t)n_cats_total * (size_t)n_cats_total, 0.0);
        if (compute_normal_equations(ct, D_eff.data(), N.data(),
                                     cat_offset.data(), st.K,
                                     (size_t)n_cats_total) != RK_OK) {
            res.status = RK_ERR_BADARG; return res;
        }

        // b[k][j] = T_kj * n - Σ_{c∈bucket(k,j)} X[c]  (marginal defect relative to fixed n)
        // GREG uses fixed proportions T_kj; target COUNT = T_kj * n (original sum(X_init)).
        // Use st.n (original target mass), not current W, so normalization is consistent
        // across active-set iterations regardless of which cells are clamped.
        // O(K × M_cell) scatter-accumulate — one pass over cells per margin.
        const double n_total = static_cast<double>(st.n);
        std::vector<double> b((size_t)n_cats_total, 0.0);
        // First: scatter X[c] into each margin's buckets in O(M_cell) per margin
        std::fill(bucket_b.begin(), bucket_b.begin() + max_cats, 0.0);  // reuse hoisted vector
        for (int k = 0; k < st.K; k++) {
            std::fill(bucket_b.begin(), bucket_b.begin() + st.cat_counts[k], 0.0);
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) bucket_b[g] += X[c];
            }
            for (int j = 0; j < st.cat_counts[k]; j++)
                b[(size_t)cat_offset[k] + (size_t)j] = st.targets[k][j] * n_total - bucket_b[j];
        }

        // LDLT factor N and solve N*lambda = b
        if (ldlt_factor_inplace(N.data(), (size_t)n_cats_total, 1e-10) != RK_OK) {
            res.status = RK_ERR_BADARG; return res;
        }
        ldlt_solve(N.data(), (size_t)n_cats_total, b.data());
        const std::vector<double>& lambda = b;

        // Newton update: X_new[c] = X_init[c] * (1 + Σ_k lambda[cat_offset[k] + g_k[c]])
        // For free cells only. Fixed cells stay at L_c or U_c.
        bool any_clamped = false;
        for (int c = 0; c < ct.M_cell; c++) {
            if (fixed_lo[c]) { X[c] = L_cell[c]; continue; }
            if (fixed_hi[c]) { X[c] = U_cell[c]; continue; }
            double sum_lambda = 0.0;
            for (int k = 0; k < st.K; k++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k])
                    sum_lambda += lambda[(size_t)cat_offset[k] + (size_t)g];
            }
            double X_new = X_init[c] > kEps ? X_init[c] * (1.0 + sum_lambda) : 0.0;
            if (X_new < L_cell[c] - 1e-10) {
                X[c] = L_cell[c];
                fixed_lo[c] = true;
                any_clamped = true;
            } else if (X_new > U_cell[c] + 1e-10) {
                X[c] = U_cell[c];
                fixed_hi[c] = true;
                any_clamped = true;
            } else {
                X[c] = std::clamp(X_new, L_cell[c], U_cell[c]);
            }
        }

        if (!any_clamped) {
            // No new bound violations: exact solution found
            res.status = RK_OK;
            res.convergence_iter = newton_iter + 1;
            break;
        }
        // Else: iterate with active-set updated
    }

    // If loop exhausted without `any_clamped == false` break, active-set did not converge cleanly.
    // Report NOCONV rather than silently claiming OK. The final X is still a valid feasible
    // point (all cells clamped within bounds), so metrics are computed and returned.
    // Callers can still use the weights at NOCONV — they satisfy bounds, just not exact margins.

    // Compute all 6 metrics at exit (same as raking/sinkhorn)
    double W = 0.0;
    for (int c = 0; c < ct.M_cell; c++) W += X[c];
    std::vector<double> bucket(*std::max_element(st.cat_counts, st.cat_counts + st.K));
    constexpr double kMetricEps = 1e-10, kChi2Floor = 1.0;
    double errRp = 0.0, mean_err_sum = 0.0, kl_max = 0.0, chi2_total = 0.0, grake_norm = 0.0;
    for (int k = 0; k < st.K; k++) {
        const int nj = st.cat_counts[k];
        std::fill(bucket.begin(), bucket.begin() + nj, 0.0);
        for (int c = 0; c < ct.M_cell; c++) {
            int g = ct.g_per_cell[k][c];
            if (g >= 0 && g < nj) bucket[g] += X[c];
        }
        double max_k = 0.0, kl_k = 0.0;
        for (int j = 0; j < nj; j++) {
            double S_p = bucket[j] / W, T = st.targets[k][j];
            double err = std::fabs(S_p - T);
            if (err > max_k) max_k = err;
            if (err > errRp) errRp = err;
            if (T > 0.0) kl_k += T * std::log((T + kMetricEps) / (S_p + kMetricEps));
            double obs = bucket[j], pop_kj = T * W;
            chi2_total += (obs - pop_kj) * (obs - pop_kj) / (pop_kj + kChi2Floor);
            double nm = std::fabs(obs - pop_kj) / (1.0 + std::fabs(pop_kj));
            if (nm > grake_norm) grake_norm = nm;
        }
        mean_err_sum += max_k;
        if (kl_k > kl_max) kl_max = kl_k;
    }
    res.max_error   = errRp;
    res.kl          = kl_max;
    res.chi2        = chi2_total;
    res.mean_error  = mean_err_sum / (st.K > 0 ? st.K : 1);
    res.grake_norm  = grake_norm;
    res.convergence_objective = chi2_total;
    res.best_error = chi2_total;

    // Exit: obs expansion + clamp (same as raking/sinkhorn, no renorm)
    const double hi_obs = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        double mult = (X_init[c] > kEps) ? X[c] / X_init[c] : 1.0;
        st.weights[i] = std::clamp(st.weights[i] * mult, lo, hi_obs);
    }

    res.best_weights.resize(st.n);
    std::copy(st.weights, st.weights + st.n, res.best_weights.begin());

    return res;
}

} // namespace lbw
```

**CRITICAL before writing:** verify `ldlt_solve` signature in `calib_linalg.hpp` and match it exactly. The current stub says `void ldlt_solve(const double* L, const double* d_diag, double* b, size_t n)` — if we store L and d in-place in the same matrix, we need to update the signature. Read the actual .hpp and decide: either update the .hpp or use separate arrays.

### Step 2.3: Add to Makevars
Add `greg.cpp` to `PKG_SOURCES` in both Makevars files.

### Step 2.4: Wire into c_api.cpp and r_bridge.cpp

**c_api.cpp:** add `#include "greg.hpp"`. In the dispatch, replace the `RK_ALG_GREG` stub case:
```cpp
case RK_ALG_GREG: {
    auto gres = lbw::greg_solve(st);
    if (result) {
        result->status                       = gres.status;
        result->iterations                   = gres.iterations;
        result->max_error                    = gres.max_error;
        result->convergence_metric           = gres.convergence_metric;
        result->convergence_rule             = gres.convergence_rule;
        result->convergence_tol              = gres.convergence_tol;
        result->convergence_iter             = gres.convergence_iter;
        result->convergence_objective        = gres.convergence_objective;
        result->convergence_minimized_metric = gres.convergence_minimized_metric;
        result->best_error                   = gres.best_error;
        result->best_iter                    = gres.best_iter;
        result->mean_error                   = gres.mean_error;
        result->kl                           = gres.kl;
        result->chi2                         = gres.chi2;
        result->grake_norm                   = gres.grake_norm;
        result->l1_weight_change             = gres.l1_weight_change;
        result->algorithm_used               = static_cast<int>(alg);
        std::strncpy(result->message, gres.message, sizeof(result->message) - 1);
    }
    return gres.status;
}
```

**r_bridge.cpp:** remove "greg" from the stub block. Add dispatch after sinkhorn branch:
```cpp
} else if (strcmp(method_str, "greg") == 0) {
    auto res = lbw::greg_solve(st);
    pack_solver_result(res);
    res_status     = res.status;
    res_iterations = res.iterations;
    res_max_error  = res.max_error;
    res_alg_used   = static_cast<int>(RK_ALG_GREG);
    res_mean_error = res.mean_error;
    res_kl         = res.kl;
    res_chi2       = res.chi2;
    res_best_error = res.best_error;
    res_best_iter  = res.best_iter;
    if (!res.best_weights.empty())
        res_best_weights = std::move(res.best_weights);
    else
        res_best_weights.assign(st.n, 0.0);
```

Also update `alg_name` ternary to include "greg".

### Step 2.5: Update T1a — remove "greg" from stub loop

In `tests/testthat/test-calibration-solvers.R`: change `for (m in c("chebyshev", "greg", "grake"))` to `for (m in c("chebyshev", "grake"))`.

### Step 2.6: Add D1 test

```r
test_that("D1: greg achieves chi2 <= other methods on synthetic", {
  set.seed(5)
  n <- 400
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace=TRUE)),
    b = factor(sample(c("1","2"), n, replace=TRUE))
  )
  target <- list(a=c("1"=0.4,"2"=0.4,"3"=0.2), b=c("1"=0.6,"2"=0.4))
  w_greg <- leafblower::harvest(data, target, min_weight=0.2, max_weight=5, method="greg",
                                attach_weights=FALSE)
  w_rake <- leafblower::harvest(data, target, max_weight=5, method="raking",
                                max_iterations=500, attach_weights=FALSE)
  r_greg <- attr(w_greg, "result")
  r_rake <- attr(w_rake, "result")
  expect_equal(r_greg$status, 0L, info="greg must converge")
  expect_lte(r_greg$chi2, r_rake$chi2 + 1e-6,
             label="greg chi2 <= raking chi2 (greg is chi2 minimizer)")
  expect_true(all(w_greg >= 1/5 - 1e-10 & w_greg <= 5 + 1e-10),
              info="greg bounds must hold")
})
```

### Step 2.7: Build gate + smoke test + regression
```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
Rscript -e '
  set.seed(5); n <- 200L
  data <- data.frame(a=factor(sample(c("1","2","3"),n,T)))
  target <- list(a=c("1"=0.5,"2"=0.3,"3"=0.2))
  w <- leafblower::harvest(data, target, max_weight=5, method="greg", attach_weights=FALSE)
  r <- attr(w,"result")
  cat(sprintf("status=%d chi2=%.4e max_error=%.4e iter=%d\n",
              r$status, r$chi2, r$max_error, r$iterations))
  stopifnot(r$status == 0, r$max_error < 0.1)
' 2>&1 | grep -v Welcome
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 363.

### Step 2.8: Commit
```bash
git add src/greg.hpp src/greg.cpp src/c_api.cpp src/r_bridge.cpp \
        src/Makevars src/Makevars.in tests/testthat/test-calibration-solvers.R
git commit -m "$(cat <<'EOF'
feat(greg): Newton QP for chi2 calibration (GREG, Deville-Sarnal 1992)

greg_solve: one Newton iteration (exact when bounds inactive).
Active-set for bounds: ≤10 reduced iterations. Uses LDLT from
calib_linalg.cpp. D1 test: greg chi2 <= raking chi2.
T1a updated: greg removed from stub loop.
EOF
)"
bd close <ticket-id>
```

---

## Final Verification

- [ ] FAIL 0, PASS ≥ 363
- [ ] `grep "case RK_ALG_GREG" src/c_api.cpp` → dispatches to greg_solve
- [ ] `grep "greg" src/r_bridge.cpp | grep "Rf_error"` → 0
- [ ] D1 test: greg chi2 ≤ raking chi2
- [ ] `method="chebyshev"` and `method="grake"` still stub (Plan E)
