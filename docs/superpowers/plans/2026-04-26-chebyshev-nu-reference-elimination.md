# Chebyshev ν Fix: Reference Category Elimination — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix chebyshev/grake LP IPM so that Σ X[c] = n_d is enforced in the Newton system, enabling convergence on K≥4 complex margin problems where the current code produces W→0 collapse and NaN metrics.

**Architecture:** Reference category elimination removes one algebraically-redundant LP constraint row per multi-category margin from the Newton system, making schur_nu > 0 (non-degenerate). A second Sherman-Morrison step then adds the normalization dual ν correctly, forcing the Newton step to be sum-preserving.

**Tech Stack:** C++17, R (devtools::test), leafblower LDLT solver (calib_linalg.hpp).

---

## Baseline

```bash
Rscript -e 'devtools::test()' 2>&1 | tail -4
# Expected: FAIL 0 | WARN 31 | SKIP 5 | PASS 378
```

Build gate: `R CMD INSTALL --preclean . 2>&1 | tail -3` → `* DONE`

---

## File Map

| File | Change |
|------|--------|
| `src/chebyshev.cpp` | Major: index bookkeeping, reduced Newton system, ν second SM |
| `src/calib_linalg.hpp` | Add `compute_normal_equations_reduced` declaration |
| `src/calib_linalg.cpp` | Add `compute_normal_equations_reduced` implementation |
| `tests/testthat/test-calib-linalg.R` | Add schur_nu diagnostic test + single-cat margin test |
| `tests/testthat/test-calibration-solvers.R` | E1/E2 must stay GREEN (existing) |

---

## Task 0: Fix W=0 guard for best_errRp (leafblower-ng5f)

**Files:** `src/chebyshev.cpp`

This fix prevents NaN post-loop metrics when W collapses. It's a prerequisite guard that must be in place before the ν fix, since both deal with W→0 scenarios.

- [ ] **Read** `src/chebyshev.cpp` lines 174-190

- [ ] **Find and replace** the metrics + best_errRp block (around lines 178-184). Current code:
```cpp
        CellMetrics cm;
        if (W > 1e-300) cm = lbw::compute_cell_metrics(st, ct, X, W, bucket_tmp);
        double errRp = cm.errRp, mean_err = cm.mean_err;
        double kl_max = cm.kl, chi2_total = cm.chi2, grake_norm = cm.grake_norm;
        double l1_weight = 0.0;  // not tracked per IPM step (no prev weights)

        // Track X with the best actual calibration error seen so far.
        ...
        if (errRp < best_errRp) { best_errRp = errRp; res.best_iter = iter+1; X_best = X; }
```

Replace — move best_errRp update inside the W guard:
```cpp
        CellMetrics cm;
        if (W > 1e-300) {
            cm = lbw::compute_cell_metrics(st, ct, X, W, bucket_tmp);
            if (cm.errRp < best_errRp) { best_errRp = cm.errRp; res.best_iter = iter+1; X_best = X; }
        }
        double errRp = cm.errRp, mean_err = cm.mean_err;
        double kl_max = cm.kl, chi2_total = cm.chi2, grake_norm = cm.grake_norm;
        double l1_weight = 0.0;  // not tracked per IPM step (no prev weights)
```
(Remove the old standalone `if (errRp < best_errRp)` line.)

- [ ] **Find post-loop block** (around line 430). Currently:
```cpp
    double W_final = 0.0;
    for (int c = 0; c < ct.M_cell; c++) W_final += X_out[c];
    {
        auto m2 = lbw::compute_cell_metrics(st, ct, X_out, W_final, bucket_tmp);
        res.max_error  = m2.errRp;
        res.kl         = m2.kl;
        res.chi2       = m2.chi2;
        res.mean_error = m2.mean_err;
        res.grake_norm = m2.grake_norm;
    }
```

Add W_final guard:
```cpp
    double W_final = 0.0;
    for (int c = 0; c < ct.M_cell; c++) W_final += X_out[c];
    if (W_final > 1e-300) {
        auto m2 = lbw::compute_cell_metrics(st, ct, X_out, W_final, bucket_tmp);
        res.max_error  = m2.errRp;
        res.kl         = m2.kl;
        res.chi2       = m2.chi2;
        res.mean_error = m2.mean_err;
        res.grake_norm = m2.grake_norm;
    }
    // else: X_best also collapsed — leave metrics at default 0.0, not NaN
```

- [ ] **Compile gate:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
# Expected: * DONE (leafblower)
```

- [ ] **Run tests:**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -4
# Expected: FAIL 0 | PASS 378
```

- [ ] **Commit:**
```bash
git add src/chebyshev.cpp
git commit -m "fix(chebyshev): guard best_errRp update and post-loop metrics against W=0

Closes: leafblower-ng5f"
```

---

## Task 1: TDD RED — Write failing tests

**Files:** `tests/testthat/test-calib-linalg.R`

Write tests that will fail until the ν fix is implemented.

- [ ] **Read** `tests/testthat/test-calib-linalg.R` (existing file)

- [ ] **Add** at the end of `test-calib-linalg.R`:

```r
# ──────────────────────────────────────────────────────────────────────────────
# chebyshev ν fix: reference elimination makes schur_nu > 0
# ──────────────────────────────────────────────────────────────────────────────

test_that("chebyshev schur_nu > 0 with verbose=2 (non-degeneracy diagnostic)", {
  set.seed(42); n <- 200L
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, TRUE)),
    b = factor(sample(c("1","2"),     n, TRUE))
  )
  target <- list(a = c("1"=0.4,"2"=0.4,"3"=0.2), b = c("1"=0.6,"2"=0.4))
  log_lines <- character(0)
  r <- leafblower::harvest(data, target, method = "chebyshev",
                           min_weight = 0.2, max_weight = 5,
                           max_iterations = 5, attach_weights = FALSE,
                           verbose = 2)
  # schur_nu is logged to stderr (R message channel) at verbose >= 2
  log_lines <- capture.output(
    leafblower::harvest(data, target, method = "chebyshev",
                        min_weight = 0.2, max_weight = 5,
                        max_iterations = 5, attach_weights = FALSE,
                        verbose = 2),
    type = "message"
  )
  schur_lines <- grep("schur_nu=", log_lines, value = TRUE)
  expect_gt(length(schur_lines), 0L, label = "schur_nu should be logged at verbose=2")
  if (length(schur_lines) > 0) {
    val <- as.numeric(regmatches(schur_lines[1],
                                 regexpr("[0-9.eE+-]+", schur_lines[1])))
    expect_gt(val, 1e-6, label = "schur_nu must be positive (non-degenerate)")
  }
})

test_that("chebyshev: single-category margin does not crash", {
  # K=2, margin b has only 1 category (trivial constraint)
  set.seed(99); n <- 100L
  data <- data.frame(
    a = factor(sample(c("1","2"), n, TRUE)),
    b = factor(rep("x", n))   # single category
  )
  target <- list(a = c("1"=0.5,"2"=0.5), b = c("x"=1.0))
  r <- leafblower::harvest(data, target, method = "chebyshev",
                           min_weight = 0.2, max_weight = 5,
                           max_iterations = 20, attach_weights = FALSE,
                           verbose = 0)
  # Must not crash; weights must be finite and non-negative
  expect_true(all(is.finite(r)), label = "chebyshev with 1-cat margin: weights finite")
  expect_true(all(r >= 0), label = "chebyshev with 1-cat margin: weights non-negative")
})
```

- [ ] **Run tests to verify RED:**
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calib-linalg.R")' 2>&1 | tail -5
```
Expected: schur_nu test FAILS (no log line). Single-cat test may pass or crash — either is acceptable RED state.

- [ ] **Commit RED tests:**
```bash
git add tests/testthat/test-calib-linalg.R
git commit -m "test(chebyshev): RED tests for schur_nu diagnostic + single-cat margin"
```

---

## Task 2: Add `compute_normal_equations_reduced` to calib_linalg

**Files:** `src/calib_linalg.hpp`, `src/calib_linalg.cpp`

- [ ] **Read** `src/calib_linalg.hpp` and `src/calib_linalg.cpp`

- [ ] **Add declaration** to `src/calib_linalg.hpp` after the existing `compute_normal_equations` declaration:

```cpp
// Compute N = A_red × diag(D) × A_red^T where A_red excludes reference margins.
// full_to_red[m] == -1 for reference margins (skipped); >= 0 for reduced index.
// nct_red: count of non-reference constraints (rows of A_red).
// N output: nct_red × nct_red, row-major.
int compute_normal_equations_reduced(
    const CellTable& ct,
    const double* D,
    double* N,
    const int* cat_offset,
    int K,
    size_t nct_red,
    const int* full_to_red
) noexcept;
```

- [ ] **Add implementation** to `src/calib_linalg.cpp` after the existing `compute_normal_equations` function:

```cpp
int compute_normal_equations_reduced(const CellTable& ct,
                                      const double* D,
                                      double* N,
                                      const int* cat_offset,
                                      int K,
                                      size_t nct_red,
                                      const int* full_to_red)
{
    std::fill(N, N + nct_red * nct_red, 0.0);
    for (int c = 0; c < ct.M_cell; c++) {
        if (D[c] <= 0.0) continue;
        for (int k1 = 0; k1 < K; k1++) {
            int j1 = ct.g_per_cell[k1][c];
            if (j1 < 0) continue;
            int m1 = cat_offset[k1] + j1;
            int r1 = full_to_red[m1];
            if (r1 < 0) continue;  // reference margin — skip
            for (int k2 = 0; k2 < K; k2++) {
                int j2 = ct.g_per_cell[k2][c];
                if (j2 < 0) continue;
                int m2 = cat_offset[k2] + j2;
                int r2 = full_to_red[m2];
                if (r2 < 0) continue;  // reference margin — skip
                N[static_cast<size_t>(r1) * nct_red + static_cast<size_t>(r2)] += D[c];
            }
        }
    }
    return RK_OK;
}
```

- [ ] **Compile gate:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
# Expected: * DONE (leafblower)
```

- [ ] **Commit:**
```bash
git add src/calib_linalg.hpp src/calib_linalg.cpp
git commit -m "feat(calib_linalg): add compute_normal_equations_reduced for nu-fix"
```

---

## Task 3: Add reference index bookkeeping to chebyshev.cpp

**Files:** `src/chebyshev.cpp`

- [ ] **Read** `src/chebyshev.cpp` lines 44-50 (cat_offset build) and lines 140-160 (hoisted work vectors, pre-loop)

- [ ] **After the cat_offset/nct build** (after the line `if (nct > kNCatsTotalMax) { ... }` block, before the w_kj/Tgt/T_flat initialization), add:

```cpp
    // Reference category elimination: drop last category per multi-cat margin.
    // Makes schur_nu = D_nu - e^T·N_eff^{-1}·e > 0 (non-degenerate).
    // Single-category margins (cat_counts[k]==1) are never eliminated.
    int nct_red_count = 0;
    for (int k = 0; k < st.K; k++)
        if (st.cat_counts[k] >= 2) nct_red_count++;
    const int nct_red = nct - nct_red_count;

    std::vector<bool> is_ref(nct, false);
    std::vector<int>  full_to_red(nct, -1);
    std::vector<int>  red_to_full(nct_red);
    {
        int nr = 0;
        for (int k = 0; k < st.K; k++) {
            for (int j = 0; j < st.cat_counts[k]; j++) {
                int m = cat_offset[k] + j;
                if (st.cat_counts[k] >= 2 && j == st.cat_counts[k] - 1) {
                    is_ref[m] = true;
                } else {
                    full_to_red[m] = nr;
                    red_to_full[nr++] = m;
                }
            }
        }
        // Invariant: nr == nct_red
    }
```

- [ ] **Replace the hoisted work vectors** (lines 142-148). Find this block:
```cpp
    std::vector<double> D_eff(ct.M_cell), D_marg(nct);
    std::vector<double> N0((size_t)nct*(size_t)nct);
    std::vector<double> u_vec(nct), v_vec(nct), rhs_v(nct), w_sol(nct);
    std::vector<double> dX(ct.M_cell);
    std::vector<double> dS_up(nct), dS_dn(nct);
    std::vector<double> dY_lo(ct.M_cell), dY_hi(ct.M_cell), dY_up(nct), dY_dn(nct);
    std::vector<double> dlambda(nct);   // hoisted: Sherman-Morrison result
    std::vector<double> delta_S(nct);   // hoisted: ΔS[m] = Σ_c A_mc*ΔX[c]
```

Replace with:
```cpp
    std::vector<double> D_eff(ct.M_cell), D_marg(nct);
    std::vector<double> N_red((size_t)nct_red * (size_t)nct_red);
    std::vector<double> u_red(nct_red), v_red(nct_red), rhs_red(nct_red), w_sol_red(nct_red);
    std::vector<double> e_red(nct_red), w_e_red(nct_red);
    std::vector<double> dlambda_red(nct_red);  // reduced Newton solution (post ν correction)
    std::vector<double> dX(ct.M_cell);
    std::vector<double> dS_up(nct), dS_dn(nct);
    std::vector<double> dY_lo(ct.M_cell), dY_hi(ct.M_cell), dY_up(nct), dY_dn(nct);
    std::vector<double> delta_S(nct);          // full size — reference margins included
```

- [ ] **Compile gate:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```
This will FAIL (N0, u_vec, etc. still referenced in loop body). That is expected — proceed to Task 4.

Actually — compile now to confirm _only_ the new symbols fail, not a syntax error:
```bash
R CMD INSTALL --preclean . 2>&1 | grep -E "error:|N0|u_vec|w_sol|dlambda\[" | head -15
```
Expected: errors mentioning `N0`, `u_vec`, `w_sol`, `dlambda` — not `full_to_red` or `N_red`.

---

## Task 4: Replace Newton solve with reduced system + first SM

**Files:** `src/chebyshev.cpp`

- [ ] **Find the N_0 compute + LDLT + SM block** in the IPM loop (around lines 234-288):

```cpp
        // N_0 = A * D_eff * A^T
        if (compute_normal_equations(ct, D_eff.data(), N0.data(), ...
        ...
        // Two back-solves (Sherman-Morrison)
        std::copy(rhs_v.begin(), rhs_v.end(), w_sol.begin());
        ldlt_solve(N0.data(), static_cast<size_t>(nct), w_sol.data());
        std::copy(u_vec.begin(), u_vec.end(), v_vec.begin());
        ldlt_solve(N0.data(), static_cast<size_t>(nct), v_vec.data());
        ...
        for (int m = 0; m < nct; m++) dlambda[m] = w_sol[m] - sm_coeff*v_vec[m];
```

- [ ] **Replace the RHS computation** (lines ~247-252). Find:
```cpp
        for (int m = 0; m < nct; m++) {
            double rmu_up = sigma_mu - s_up[m]*y_up[m];
            double rmu_dn = sigma_mu - s_dn[m]*y_dn[m];
            rhs_v[m] = -(S[m] - T_flat[m]*W)
                       + D_marg[m] * (rmu_up/s_up[m] - rmu_dn/s_dn[m]);
        }
```

Replace with reduced RHS:
```cpp
        for (int nr = 0; nr < nct_red; nr++) {
            int m = red_to_full[nr];
            double rmu_up = sigma_mu - s_up[m]*y_up[m];
            double rmu_dn = sigma_mu - s_dn[m]*y_dn[m];
            rhs_red[nr] = -(S[m] - T_flat[m]*W)
                          + D_marg[m] * (rmu_up/s_up[m] - rmu_dn/s_dn[m]);
        }
```

- [ ] **Keep** the r_delta_stat / margin_delta_center / Theta loops unchanged (they use all nct margins via y_up/y_dn/w_kj — no change needed).

- [ ] **Replace the SM setup + LDLT + two back-solves** (from `double alpha_sm = ...` through `dlambda[m] = w_sol[m] - sm_coeff*v_vec[m]`):
```cpp
        double alpha_sm = (Theta > 1e-300) ? 1.0/Theta : 0.0;
        for (int nr = 0; nr < nct_red; nr++) u_red[nr] = w_kj[red_to_full[nr]];

        // LDLT factor N_red (reduced, nct_red × nct_red)
        if (lbw::compute_normal_equations_reduced(ct, D_eff.data(), N_red.data(),
                                                  cat_offset.data(), st.K,
                                                  static_cast<size_t>(nct_red),
                                                  full_to_red.data()) != RK_OK) {
            res.status = RK_ERR_BADARG; return res;
        }
        if (ldlt_factor_inplace(N_red.data(), static_cast<size_t>(nct_red), kEpsLdlt) != RK_OK) {
            res.status = RK_ERR_BADARG; return res;
        }

        // Back-solve 1: w_sol_red = N_red^{-1} · rhs_red
        std::copy(rhs_red.begin(), rhs_red.end(), w_sol_red.begin());
        ldlt_solve(N_red.data(), static_cast<size_t>(nct_red), w_sol_red.data());

        // Back-solve 2: v_red = N_red^{-1} · u_red  (for first SM: δ correction)
        std::copy(u_red.begin(), u_red.end(), v_red.begin());
        ldlt_solve(N_red.data(), static_cast<size_t>(nct_red), v_red.data());

        double utv = 0.0, utw = 0.0;
        for (int nr = 0; nr < nct_red; nr++) { utv += u_red[nr]*v_red[nr]; utw += u_red[nr]*w_sol_red[nr]; }
        double sm_denom = 1.0 + alpha_sm*utv;
        double sm_coeff = (std::fabs(sm_denom) > 1e-300) ? (alpha_sm*utw/sm_denom) : 0.0;
        for (int nr = 0; nr < nct_red; nr++) dlambda_red[nr] = w_sol_red[nr] - sm_coeff*v_red[nr];
```

- [ ] **Compile gate:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```
Still expects FAIL (dX loop and w_dot_dlambda still use old `dlambda` / `N0`). Check errors are only on those remaining references.

---

## Task 5: Add ν second Sherman-Morrison

**Files:** `src/chebyshev.cpp`

After the first SM dlambda_red computation (end of Task 4), add the ν block:

- [ ] **Add after** `for (int nr = 0; nr < nct_red; nr++) dlambda_red[nr] = w_sol_red[nr] - sm_coeff*v_red[nr];`:

```cpp
        // Back-solve 3 (ν): e_red[nr] = Σ_{c in margin red_to_full[nr]} D_eff[c]
        std::fill(e_red.begin(), e_red.end(), 0.0);
        for (int c = 0; c < ct.M_cell; c++) {
            for (int k = 0; k < st.K; k++) {
                int g = ct.g_per_cell[k][c];
                if (g < 0 || g >= st.cat_counts[k]) continue;
                int m = cat_offset[k] + g;
                int nr = full_to_red[m];
                if (nr >= 0) e_red[nr] += D_eff[c];
            }
        }
        std::copy(e_red.begin(), e_red.end(), w_e_red.begin());
        ldlt_solve(N_red.data(), static_cast<size_t>(nct_red), w_e_red.data());
        // Apply first SM correction to w_e_red
        double ute = 0.0;
        for (int nr = 0; nr < nct_red; nr++) ute += u_red[nr] * w_e_red[nr];
        double sm_coeff_e = (std::fabs(sm_denom) > 1e-300) ? (alpha_sm*ute/sm_denom) : 0.0;
        for (int nr = 0; nr < nct_red; nr++) w_e_red[nr] -= sm_coeff_e * v_red[nr];

        // Compute ν: schur_nu = D_nu - e^T·w_e > 0 (guaranteed by reference elimination)
        double D_nu = 0.0;
        for (int c = 0; c < ct.M_cell; c++) D_nu += D_eff[c];
        double eTw_e = 0.0, eTdlambda = 0.0;
        for (int nr = 0; nr < nct_red; nr++) {
            eTw_e     += e_red[nr] * w_e_red[nr];
            eTdlambda += e_red[nr] * dlambda_red[nr];
        }
        const double schur_nu = D_nu - eTw_e;
        if (st.verbose >= 2 && iter == 0) {
            char msg[128];
            std::snprintf(msg, sizeof(msg), "chebyshev: schur_nu=%.4e (iter 0)", schur_nu);
            st.log(msg);
        }
        const double r_nu = W - n_d;
        const double dnu  = (schur_nu > 1e-8) ? (-r_nu - eTdlambda) / schur_nu : 0.0;
        // Correct dlambda_red with ν contribution
        for (int nr = 0; nr < nct_red; nr++) dlambda_red[nr] -= dnu * w_e_red[nr];
```

- [ ] **Compile gate:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```
Still expects FAIL on dX loop and w_dot_dlambda references to old `dlambda`.

---

## Task 6: Update dX, delta_S, w_dot_dlambda

**Files:** `src/chebyshev.cpp`

- [ ] **Replace dX computation** (find lines ~291-298):
```cpp
        // ΔX[c] = D_eff[c] * Σ_k Δλ[m_k]
        std::fill(dX.begin(), dX.end(), 0.0);
        for (int k = 0; k < st.K; k++)
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k])
                    dX[c] += dlambda[cat_offset[k]+g];
            }
        for (int c = 0; c < ct.M_cell; c++) dX[c] *= D_eff[c];
```

Replace with:
```cpp
        // ΔX[c] = D_eff[c] * (Σ_k Δλ_red[nr_k(c)] + Δν)
        // Reference margins (full_to_red[m]==-1) contribute 0 to Δλ sum.
        // Δν applies to ALL cells regardless of reference status.
        for (int c = 0; c < ct.M_cell; c++) {
            double sum_dlam = 0.0;
            for (int k = 0; k < st.K; k++) {
                int g = ct.g_per_cell[k][c];
                if (g < 0 || g >= st.cat_counts[k]) continue;
                int m = cat_offset[k] + g;
                int nr = full_to_red[m];
                if (nr >= 0) sum_dlam += dlambda_red[nr];
            }
            dX[c] = D_eff[c] * (sum_dlam + dnu);
        }
```

- [ ] **Replace w_dot_dlambda** (find lines ~301-302):
```cpp
        double w_dot_dlambda = 0.0;
        for (int m = 0; m < nct; m++) w_dot_dlambda += w_kj[m]*dlambda[m];
```

Replace with:
```cpp
        double w_dot_dlambda = 0.0;
        for (int nr = 0; nr < nct_red; nr++)
            w_dot_dlambda += w_kj[red_to_full[nr]] * dlambda_red[nr];
```

Note: Reference margins have dlambda_ref=0 implicitly. Their w_kj contribution to d_delta is 0 per Newton step. This is correct — ν handles the normalization instead.

- [ ] **Compile gate:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
# Expected: * DONE (leafblower)
```

- [ ] **If still failing**, grep for any remaining `dlambda[` or `N0` references:
```bash
grep -n "dlambda\[m\]\|dlambda\[cat\|N0\.\|u_vec\[m\]\|v_vec\[m\]\|w_sol\[m\]" src/chebyshev.cpp
```
Fix any remaining references by replacing with the reduced equivalents.

---

## Task 7: Run tests and commit

**Files:** `src/chebyshev.cpp`, `tests/testthat/test-calib-linalg.R`

- [ ] **Run full test suite:**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -6
```
Expected: FAIL 0 | PASS ≥ 379 (378 existing + at minimum the single-cat test green). E1 and E2 must be GREEN.

If tests fail:
- **E1/E2 fail**: The reduced Newton system is incorrect. Check dX reconstruction — verify Δν is included for all cells.
- **schur_nu test fails**: Verify the verbose log is written. Check `st.log(msg)` is called with `verbose >= 2`.
- **Single-cat test crashes**: Check the `is_ref` guard — `st.cat_counts[k] >= 2` must prevent dropping the only category.

- [ ] **Run quick benchmark** to verify no regression on simple K=2 problem:
```bash
Rscript -e '
suppressPackageStartupMessages(library(leafblower))
set.seed(11); n <- 200L
data <- data.frame(a=factor(sample(c("1","2","3"),n,TRUE)),
                   b=factor(sample(c("1","2"),n,TRUE)))
target <- list(a=c("1"=0.4,"2"=0.4,"3"=0.2), b=c("1"=0.6,"2"=0.4))
wc <- leafblower::harvest(data, target, min_weight=0.2, max_weight=5,
                          method="chebyshev", max_iterations=60, attach_weights=FALSE)
rc <- attr(wc,"result")
wr <- leafblower::harvest(data, target, min_weight=0.2, max_weight=5,
                          method="raking", max_iterations=500, attach_weights=FALSE)
rr <- attr(wr,"result")
cat("cheb max_err:", rc$max_error, " raking:", rr$max_error, "\n")
cat("cheb <= raking:", rc$max_error <= rr$max_error + 1e-6, "\n")
' 2>&1 | grep -v "Welcome\|Warning"
```
Expected: `cheb <= raking: TRUE`

- [ ] **Commit:**
```bash
git add src/chebyshev.cpp tests/testthat/test-calib-linalg.R
git commit -m "$(cat <<'EOF'
fix(chebyshev): reference elimination + nu second SM for sum-preserving Newton step

Normalization dual nu is algebraically degenerate for sum-to-1 targets
(normalization row = linear combination of margin rows, schur_nu=0).

Fix: drop last category per multi-cat margin from Newton system (reduced
nct_red = nct - K). N_eff_red is full-rank; schur_nu > 0. Third LDLT
back-solve gives N_eff_red^{-1}*e, then nu = (-r_nu - e^T*dlambda) /
schur_nu. Corrects dlambda_red and dX. Newton step is now sum-preserving.
Single-cat margins (cat_counts[k]<2) not eliminated.

Closes: leafblower-8i6f
EOF
)"
```

---

## Task 8: Close worktree and tickets

- [ ] **Close both tickets:**
```bash
bd close leafblower-ng5f leafblower-8i6f 2>&1
```

- [ ] **Delete cheb-renorm worktree** (superseded by this fix):
```bash
git worktree remove .worktrees/cheb-renorm --force
git branch -d cheb-renorm
```

- [ ] **Run final test suite:**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -4
# Expected: FAIL 0 | PASS ≥ 379
```

---

## Self-Review

**Spec coverage check:**

| Spec section | Plan task |
|---|---|
| Edge case: n_cats[k]=1 | Task 3 — `st.cat_counts[k] >= 2` guard in bookkeeping |
| compute_normal_equations_reduced | Task 2 |
| Reference bookkeeping | Task 3 |
| Reduced RHS | Task 4 |
| Three back-solves | Task 4 + Task 5 |
| ν computation with schur_nu > 1e-8 guard | Task 5 |
| schur_nu verbose log at iter==0 | Task 5 |
| dX reconstruction with Δν | Task 6 |
| w_dot_dlambda from reduced | Task 6 |
| W=0 guard for best_errRp | Task 0 |
| E1/E2 tests GREEN | Task 7 |
| schur_nu > 0 test | Task 1 + Task 7 |
| Single-cat margin test | Task 1 + Task 7 |
| Delete cheb-renorm | Task 8 |

**Gap check**: The `delta_S` computation (full nct, unchanged) — Task 6 leaves `delta_S` loop using `dX[c]` (which now includes Δν). This is correct — dS updates use the true primal step, not just the margin contribution.

**Type consistency**: `dlambda_red` (nct_red), `dlambda` removed, `dnu` (scalar) — consistent through Tasks 4-6.

**Placeholder scan**: None found.
