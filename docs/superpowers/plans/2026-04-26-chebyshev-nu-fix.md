# Fix chebyshev IPM: Add Normalization Dual ν (second Sherman-Morrison)

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use checkbox syntax.

**Goal:** Fix chebyshev/grake LP IPM correctness. `Σ X[c] = n` is an LP equality constraint but its dual `ν` is missing from the Newton system. Without it, ΔX is not sum-preserving and `W = Σ X[c]` drifts away from `n`, breaking the correspondence between the LP objective `δ` and the calibration metric `max_error`.

**Fix:** Add `ν` to the augmented normal equations via a second Sherman-Morrison rank-1 update (one additional LDLT back-solve per iteration).

**Tickets:** leafblower-ekxu, leafblower-nthn

**File:** `src/chebyshev.cpp` only.

**Baseline:** FAIL 0 | PASS 374 | SKIP 5

---

## Algorithm: Second Sherman-Morrison for ν

After computing `N_eff = N_0 + α_δ·u·u^T` (the δ SM update), the full augmented system is:

```
N_eff · Δλ + e · Δν = rhs_margin
e^T · Δλ + D_ν · Δν = -r_nu
```

where:
- `e[m] = Σ_{c∈m} D_eff[c]` — effective D column sum into margin space
- `D_ν = Σ_c D_eff[c]` — total effective weight
- `r_nu = Σ X[c] - n` — normalization primal residual (drift from target)

**Block elimination:**

1. Compute `w_e = N_eff^{-1} · e` (third LDLT back-solve + SM correction)
2. `schur_nu = D_ν - e^T · w_e` (Schur complement denominator for ν)
3. `eTdlambda = e^T · dlambda` (e dot the uncorrected dlambda = N_eff^{-1}·rhs)
4. `Δν = (-r_nu - eTdlambda) / schur_nu`
5. Correct: `dlambda[m] -= Δν · w_e[m]` for all m
6. Correct: `dX[c] += D_eff[c] · Δν` for all c (after existing dX computation)
7. Recompute `w_dot_dlambda` from corrected dlambda, then update `d_delta`

---

## Task 1 — TDD RED test

**Create ticket:** `bd create --title="test: chebyshev max_err <= raking (E1 correctness)" --type=task --priority=1 2>&1 | tail -2`

Restore E1/E2 to correctness assertions (currently wiring-only). Find the existing E1/E2 tests:

```bash
grep -n "E1\|E2.*grake\|chebyshev.*max_err\|grake.*grake_norm" tests/testthat/test-calibration-solvers.R | head -10
```

Replace the wiring assertions with correctness assertions:
```r
# E1: chebyshev max_err <= raking max_err
expect_equal(r_cheb$status, 0L)
expect_lte(r_cheb$max_error, r_rake$max_error + 1e-6,
           label="chebyshev max_err <= raking max_err")
expect_true(all(w_cheb >= 0.2 - 1e-10 & w_cheb <= 5 + 1e-10))

# E2: grake grake_norm <= raking grake_norm
expect_equal(r_grake$status, 0L)
expect_lte(r_grake$grake_norm, r_rake$grake_norm + 1e-6,
           label="grake grake_norm <= raking grake_norm")
expect_true(all(w_grake >= 0.2 - 1e-10 & w_grake <= 5 + 1e-10))
```

Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "E1|E2|FAIL" | head -5`

Expected: E1 and E2 FAIL (chebyshev max_err > raking = RED state).

Commit RED test.

---

## Task 2 — Fix src/chebyshev.cpp

**Ticket:** leafblower-ekxu (claim it)

Read the current IPM loop to find the Sherman-Morrison section:
```bash
grep -n "utv\|utw\|sm_coeff\|dlambda\|w_dot_dlambda\|d_delta" src/chebyshev.cpp | head -20
```

### Step 2.1: Add e_vec and D_nu to hoisted work vectors

Find the hoisted work vectors block (before the for loop):
```bash
grep -n "dX\|dS_up\|w_sol\|v_vec\|u_vec" src/chebyshev.cpp | head -8
```

Add after the existing vectors:
```cpp
    std::vector<double> e_vec(nct);    // hoisted: e[m] = Σ_{c∈m} D_eff[c]
    std::vector<double> w_e(nct);     // hoisted: N_eff^{-1} · e_vec
```

### Step 2.2: Compute e_vec and D_nu after D_eff is ready

Find where `D_eff[c]` is fully computed (after the inner loop). Add immediately after:
```cpp
        // Normalization dual ν correction: e[m] = Σ_{c∈m} D_eff[c]
        std::fill(e_vec.begin(), e_vec.end(), 0.0);
        double D_nu = 0.0;
        for (int c = 0; c < ct.M_cell; c++) {
            D_nu += D_eff[c];
            for (int k = 0; k < st.K; k++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k])
                    e_vec[cat_offset[k]+g] += D_eff[c];
            }
        }
        const double r_nu = W - n_d;   // normalization primal residual
```

(Note: `W` is already computed in the convergence metrics block above. Verify it's in scope here. If not, compute `W` here too: `double W_now = 0.0; for (int c...)`)

### Step 2.3: Add third back-solve for N_eff^{-1}·e

Find the two existing back-solves (for w_sol and v_vec). After them, add:
```cpp
        // Third back-solve: w_e = N_eff^{-1} · e_vec (normalization SM correction)
        std::copy(e_vec.begin(), e_vec.end(), w_e.begin());
        ldlt_solve(N0.data(), static_cast<size_t>(nct), w_e.data());
        double ute = 0.0;
        for (int m = 0; m < nct; m++) ute += u_vec[m] * w_e[m];
        double sm_coeff_e = (std::fabs(sm_denom) > 1e-300)
                          ? (alpha_sm * ute / sm_denom) : 0.0;
        for (int m = 0; m < nct; m++) w_e[m] -= sm_coeff_e * v_vec[m];
```

(Note: `sm_denom = 1 + alpha_sm * utv` is already computed for the δ SM. Reuse it.)

### Step 2.4: Compute Δν and correct dlambda

Find where `dlambda[m]` is computed from the SM formula. After it, add:
```cpp
        // Compute Δν via block elimination: ν corrects for normalization drift
        double eTw_e = 0.0, eTdlambda = 0.0;
        for (int m = 0; m < nct; m++) {
            eTw_e     += e_vec[m] * w_e[m];
            eTdlambda += e_vec[m] * dlambda[m];
        }
        const double schur_nu = D_nu - eTw_e;
        const double dnu = (std::fabs(schur_nu) > 1e-14)
                         ? (-r_nu - eTdlambda) / schur_nu : 0.0;
        for (int m = 0; m < nct; m++) dlambda[m] -= dnu * w_e[m];
```

### Step 2.5: Add Δν to ΔX and recompute d_delta

Find where `dX[c] *= D_eff[c]` is applied. After it, add:
```cpp
        // Apply ν correction to ΔX: ensures Σ ΔX = 0 (sum-preserving)
        for (int c = 0; c < ct.M_cell; c++) dX[c] += D_eff[c] * dnu;
```

Find where `w_dot_dlambda` is computed. Change it to compute AFTER dlambda correction:
```cpp
        // w_dot_dlambda uses corrected dlambda (post-ν correction)
        double w_dot_dlambda = 0.0;
        for (int m = 0; m < nct; m++) w_dot_dlambda += w_kj[m] * dlambda[m];
```

If `w_dot_dlambda` is currently computed BEFORE the dlambda correction, move it after step 2.4. If it's already after, no change needed.

### Step 2.6: Build gate
```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```
Expected: `* DONE`. If compile errors, fix before proceeding.

### Step 2.7: Test correctness
```bash
Rscript -e '
  set.seed(11); n <- 200L
  data <- data.frame(
    a=factor(sample(c("1","2","3"),n,TRUE)),
    b=factor(sample(c("1","2"),n,TRUE))
  )
  target <- list(a=c("1"=0.4,"2"=0.4,"3"=0.2), b=c("1"=0.6,"2"=0.4))
  wc <- leafblower::harvest(data, target, min_weight=0.2, max_weight=5,
                            method="chebyshev", max_iterations=60, attach_weights=FALSE)
  rc <- attr(wc,"result")
  wr <- leafblower::harvest(data, target, min_weight=0.2, max_weight=5,
                            method="raking", max_iterations=500, attach_weights=FALSE)
  rr <- attr(wr,"result")
  cat("cheb: max_err=", rc$max_error, "iter=", rc$iterations, "\n")
  cat("raking: max_err=", rr$max_error, "\n")
  cat("cheb<=raking:", rc$max_error <= rr$max_error + 1e-6, "\n")
' 2>&1 | grep -v Welcome
```
Expected: `cheb<=raking: TRUE`.

### Step 2.8: Full regression
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: E1 and E2 GREEN. FAIL 0, PASS ≥ 374.

### Step 2.9: Commit
```bash
git add src/chebyshev.cpp tests/testthat/test-calibration-solvers.R
git commit -m "$(cat <<'EOF'
fix(chebyshev): add normalization dual ν via second Sherman-Morrison

LP has Σ X[c]=n equality constraint; its dual ν was missing from the
Newton system. Without ν, ΔX is not sum-preserving and W drifts away
from n, causing errRp ≠ LP objective δ. Fix: compute e[m]=Σ D_eff[c]
in bucket m, third LDLT back-solve gives N_eff⁻¹·e, then Δν from
Schur complement. Corrects dlambda and dX. E1/E2 now GREEN.
EOF
)"
bd close leafblower-ekxu 2>/dev/null || true
bd close leafblower-nthn 2>/dev/null || true
```

---

## Final Verification

- [ ] `Rscript -e 'devtools::test()' 2>&1 | tail -3` → FAIL 0, PASS ≥ 374, E1+E2 GREEN
- [ ] chebyshev max_err ≤ raking max_err on synthetic
- [ ] grake grake_norm ≤ raking grake_norm on synthetic
- [ ] `grep "r_nu\|D_nu\|e_vec\|w_e\|dnu\|schur_nu" src/chebyshev.cpp | wc -l` → > 10 (fix present)
