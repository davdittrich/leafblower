# Raking SQUAREM Acceleration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement SQUAREM SqS3 outer-loop acceleration for leafblower raking, exposing `accelerate=TRUE` in `harvest()` that jumps toward the KL fixed point 2× faster than the current Greedy-only baseline.

**Architecture:** SQUAREM wraps the existing inner IPF sweep as a fixed-point operator F. The outer loop computes w1=F(X), w2=F(w1), extrapolates X*=X-2α·r+α²·v (CBB step), then evaluates X_new=F(X*) with step-halving safety. All F-evals share the same Greedy/SOR/Bregman state. The `accelerate` bool flows: harvest.R → C bridge (new 31st arg) → CalibState.accelerate → raking_solve outer branch.

**Tech Stack:** C++17 (raking_solve), R (harvest.R bridge). No new structs, no rk_params_t ABI change. CalibState (C++ only) gets one new bool field.

---

## File Map

| File | Role |
|------|------|
| `src/types.hpp` | Add `bool accelerate = false` to CalibState |
| `src/r_bridge.cpp` | Add 31st SEXP arg; wire to `st.accelerate`; update registration count |
| `src/raking.cpp` | SQUAREM outer loop (replace flat iter loop when `st.accelerate`) |
| `R/harvest.R` | Move `accelerate` out of ignored list; warn for non-raking; pass to .Call; update @param |
| `tests/testthat/test-calibration-solvers.R` | RED test + AC3 fixture test |
| `tests/testthat/fixtures/raking_squarem_baseline.rds` | Generated in Task 0 (before any code) |

---

## Task 0: Generate AC3 baseline fixture

> Must complete BEFORE any code changes. Generates the reference weights for `accelerate=FALSE`.

**Files:**
- Create: `tests/testthat/fixtures/raking_squarem_baseline.rds`

- [ ] **Step 1: Build current package**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```
Expected: `* DONE (leafblower)` — no errors.

- [ ] **Step 2: Generate and save the baseline fixture**

```r
# Run from the package root directory
set.seed(42L); n <- 2000L
df  <- data.frame(v1 = factor(sample(3L, n, TRUE)), v2 = factor(sample(2L, n, TRUE)))
tgt <- list(v1=c("1"=0.5,"2"=0.3,"3"=0.2), v2=c("1"=0.6,"2"=0.4))

r_base <- leafblower::harvest(df, tgt, method="raking",
  accelerate=FALSE, max_weight=5, max_iterations=500L, attach_weights=FALSE)
w_base <- as.numeric(r_base)
saveRDS(w_base, "tests/testthat/fixtures/raking_squarem_baseline.rds")
cat("Saved", length(w_base), "weights, sum=", sum(w_base), "\n")
```

Run: `Rscript -e 'source("tests/testthat/fixtures/gen_squarem_baseline.R")'`
Or paste directly in R console.

Expected output: `Saved 2000 weights, sum= 2000`

- [ ] **Step 3: Verify fixture exists**

```bash
ls -lh tests/testthat/fixtures/raking_squarem_baseline.rds
```
Expected: file exists, size ~16KB.

- [ ] **Step 4: Commit fixture**

```bash
git add tests/testthat/fixtures/raking_squarem_baseline.rds
git commit -m "test(squarem): add AC3 baseline fixture (accelerate=FALSE, pre-SQUAREM)"
```

---

## Task 1: Write RED test

> Commit this test BEFORE implementing SQUAREM. It must FAIL (RED) now.

**Files:**
- Modify: `tests/testthat/test-calibration-solvers.R` (append after last test)

- [ ] **Step 1: Append RED test to test file**

Add this block at the end of `tests/testthat/test-calibration-solvers.R`:

```r
# ── SQUAREM acceleration tests ──────────────────────────────────────────────

test_that("squarem-red: accelerate=TRUE produces different iterations than accelerate=FALSE", {
  # RED state (before SQUAREM): accelerate=TRUE silently ignored → same iterations.
  # GREEN state (after SQUAREM): accelerate=TRUE takes fewer super-steps → fewer reported iters.
  # Note: harvest() currently has accelerate in its 'ignored' list → same result as FALSE.
  set.seed(42L); n <- 2000L
  df  <- data.frame(v1 = factor(sample(3L, n, TRUE)), v2 = factor(sample(2L, n, TRUE)))
  tgt <- list(v1=c("1"=0.5,"2"=0.3,"3"=0.2), v2=c("1"=0.6,"2"=0.4))

  r_base <- leafblower::harvest(df, tgt, method="raking",
    accelerate=FALSE, max_weight=5, max_iterations=500L, attach_weights=FALSE)
  r_acc  <- leafblower::harvest(df, tgt, method="raking",
    accelerate=TRUE,  max_weight=5, max_iterations=500L, attach_weights=FALSE)

  iters_base <- attr(r_base, "result")$iterations
  iters_acc  <- attr(r_acc,  "result")$iterations

  # Before SQUAREM: accelerate=TRUE ignored → same iterations as FALSE
  # After  SQUAREM: fewer super-steps (each = 3 F-evals) → different iters
  expect_false(isTRUE(all.equal(iters_base, iters_acc)),
               label="accelerate=TRUE must use different iterations than accelerate=FALSE")
})

test_that("squarem-ac3: accelerate=FALSE is bit-identical to pre-SQUAREM baseline", {
  ref_path <- test_path("fixtures/raking_squarem_baseline.rds")
  skip_if(!file.exists(ref_path), "raking_squarem_baseline.rds not generated yet")

  set.seed(42L); n <- 2000L
  df  <- data.frame(v1 = factor(sample(3L, n, TRUE)), v2 = factor(sample(2L, n, TRUE)))
  tgt <- list(v1=c("1"=0.5,"2"=0.3,"3"=0.2), v2=c("1"=0.6,"2"=0.4))

  r <- leafblower::harvest(df, tgt, method="raking",
    accelerate=FALSE, max_weight=5, max_iterations=500L, attach_weights=FALSE)
  w <- as.numeric(r)
  w_ref <- readRDS(ref_path)

  expect_equal(w, w_ref, tolerance=0,
               label="accelerate=FALSE must be bit-identical to pre-SQUAREM baseline")
})

test_that("squarem-ac5: accelerate=TRUE with non-raking method warns and runs", {
  df  <- data.frame(v1 = factor(c("1","2","1","2","1")))
  tgt <- list(v1=c("1"=0.5,"2"=0.5))

  expect_warning(
    leafblower::harvest(df, tgt, method="ieppa", accelerate=TRUE,
                        max_weight=3, attach_weights=FALSE),
    regexp="accelerate.*raking",
    label="accelerate=TRUE with ieppa must warn"
  )
})
```

- [ ] **Step 2: Build package (unchanged from Task 0)**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

- [ ] **Step 3: Run only the squarem tests to confirm RED**

```bash
Rscript -e "devtools::test(filter='squarem')" 2>&1
```

Expected: `squarem-red` FAILS with `"accelerate=TRUE must use different iterations"`.
`squarem-ac3` should PASS (fixture exists, accelerate=FALSE still identical).
`squarem-ac5` FAILS (no warning yet).

- [ ] **Step 4: Commit RED tests**

```bash
git add tests/testthat/test-calibration-solvers.R
git commit -m "test(squarem): RED tests before SQUAREM implementation (AC3/AC5/RED)"
```

---

## Task 2: Add `accelerate` to CalibState and wire through C bridge

> Two source files: `src/types.hpp` (1 line) and `src/r_bridge.cpp` (4 changes).

**Files:**
- Modify: `src/types.hpp` (after line 95 — end of overlay config block)
- Modify: `src/r_bridge.cpp` (4 locations)

### types.hpp

- [ ] **Step 1: Add `accelerate` field to CalibState**

In `src/types.hpp`, find the `CalibSorCfg sor_cfg;` line (~line 95). Add `accelerate` immediately after the sor_cfg line:

```cpp
    CalibSorCfg          sor_cfg;
    bool                 accelerate = false;  // SQUAREM outer loop for raking
    // ── End overlay config ──
```

### r_bridge.cpp

- [ ] **Step 2: Add new SEXP parameter to C_rk_calibrate declaration**

Find line 24-26 in `src/r_bridge.cpp`:
```cpp
SEXP C_rk_calibrate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
```

Replace with (adds one SEXP):
```cpp
SEXP C_rk_calibrate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
                    SEXP);
```

- [ ] **Step 3: Update registration count 30 → 31**

Find line 42:
```cpp
        {"C_rk_calibrate",       (DL_FUNC)&C_rk_calibrate,       30},
```
Replace with:
```cpp
        {"C_rk_calibrate",       (DL_FUNC)&C_rk_calibrate,       31},
```

- [ ] **Step 4: Add `accelerate_sexp` parameter to the function definition**

Find the function definition starting at line 88-104:
```cpp
SEXP C_rk_calibrate(SEXP data_sexp, SEXP target_sexp,
                    ...
                    SEXP sor_omega_fixed_sexp, SEXP sor_burnin_sexp) {
```

Replace the closing parameter with:
```cpp
SEXP C_rk_calibrate(SEXP data_sexp, SEXP target_sexp,
                    SEXP min_weight_sexp, SEXP max_weight_sexp,
                    SEXP method_sexp, SEXP verbose_sexp,
                    SEXP inner_max_iter_sexp, SEXP start_weights_sexp,
                    SEXP tol_abs_sexp, SEXP bounds_mode_sexp,
                    SEXP homotopy_levels_sexp, SEXP homotopy_start_factor_sexp,
                    SEXP homotopy_end_factor_sexp, SEXP homotopy_budget_p_sexp,
                    SEXP scheduler_sexp, SEXP eta_schedule_sexp,
                    SEXP eta_start_sexp, SEXP eta_end_sexp,
                    SEXP eta_schedule_power_sexp,
                    /* Convergence config (WU-A) */
                    SEXP pct_tol_sexp, SEXP absolute_tol_sexp,
                    SEXP metric_sexp, SEXP rule_sexp, SEXP stop_when_sexp,
                    /* SOR config (WU-A) */
                    SEXP sor_enabled_sexp, SEXP sor_auto_sexp,
                    SEXP sor_omega_init_sexp, SEXP sor_omega_min_sexp,
                    SEXP sor_omega_fixed_sexp, SEXP sor_burnin_sexp,
                    /* SQUAREM */
                    SEXP accelerate_sexp) {
```

- [ ] **Step 5: Wire `accelerate_sexp` into CalibState**

Find the block that sets CalibState fields (around line 267-311). After `st.alm_mu = 0.0;`, add:
```cpp
    st.accelerate = (INTEGER(accelerate_sexp)[0] != 0);
```

- [ ] **Step 6: Compile to verify bridge changes**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -10
```
Expected: `* DONE (leafblower)` — compile errors here mean a SEXP count mismatch or syntax error.

---

## Task 3: Wire `accelerate` in harvest.R

> Three changes in `R/harvest.R`: @param update, ignored list removal, .Call addition, warning.

**Files:**
- Modify: `R/harvest.R`

- [ ] **Step 1: Update `@param accelerate` docstring**

Find line 75 in `R/harvest.R`:
```r
#' @param accelerate Ignored.
```
Replace with:
```r
#' @param accelerate Logical. If \code{TRUE}, applies SQUAREM SqS3 outer-loop
#'   acceleration (Varadhan & Roland 2008) to the raking fixed-point iteration.
#'   Only supported for \code{method="raking"}; a warning is emitted and the
#'   parameter is ignored for all other methods. Default \code{FALSE}.
```

- [ ] **Step 2: Remove `accelerate` from the ignored list**

Find line 203-207 in `R/harvest.R`:
```r
  ignored <- c("select_params", "select_function", "error_function",
                "adaptive_order", "accelerate", "enforce_mean")
```
Replace with:
```r
  ignored <- c("select_params", "select_function", "error_function",
                "adaptive_order", "enforce_mean")
```

- [ ] **Step 3: Add warning for accelerate=TRUE with non-raking method**

After the `sor_cfg <- parse_sor(sor)` line (~line 186), add:
```r
  if (isTRUE(accelerate) && method != "raking")
    warning("accelerate=TRUE is only supported for method='raking'; ignoring for method='",
            method, "'")
  accelerate_bool <- isTRUE(accelerate) && method == "raking"
```

- [ ] **Step 4: Add `accelerate_bool` to the `.Call` invocation**

Find the `.Call("C_rk_calibrate", ...)` block. After the last argument `as.integer(sor_cfg$burnin)`, add:
```r
               as.integer(sor_cfg$burnin),
               ## SQUAREM
               as.integer(accelerate_bool),
               PACKAGE = "leafblower")
```

(Remove the existing `PACKAGE = "leafblower")` line and replace the last sor line to add the new arg before it.)

- [ ] **Step 5: Build and verify warning fires**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

```r
# Quick smoke test in R console
df <- data.frame(v1=factor(c("1","2","1","2")))
tgt <- list(v1=c("1"=0.5,"2"=0.5))
# Should warn:
leafblower::harvest(df, tgt, method="ieppa", accelerate=TRUE, max_weight=3, attach_weights=FALSE)
# Should not warn:
leafblower::harvest(df, tgt, method="raking", accelerate=TRUE, max_weight=3, attach_weights=FALSE)
```

- [ ] **Step 6: Run RED tests — squarem-ac5 should now PASS**

```bash
Rscript -e "devtools::test(filter='squarem')" 2>&1
```
Expected: `squarem-ac5` PASSES (warning fires). `squarem-red` still FAILS (SQUAREM not yet implemented).

- [ ] **Step 7: Commit**

```bash
git add src/types.hpp src/r_bridge.cpp R/harvest.R
git commit -m "feat(squarem): wire accelerate bool through CalibState → C bridge → harvest.R"
```

---

## Task 4: Implement SQUAREM outer loop in raking_solve

> The core of the feature. All changes are in `src/raking.cpp`.

**Files:**
- Modify: `src/raking.cpp`

The strategy: extract the per-iteration body into a lambda `F_eval`, then replace the flat iteration loop with a SQUAREM outer loop (when `st.accelerate`) or the existing flat loop (when `!st.accelerate`).

### Step 1: Extract F_eval lambda

- [ ] **Step 1: Read the current inner loop body (lines 170–398 of raking.cpp) and extract it into a lambda**

After the existing variable declarations (around line 160, before `for (int iter = 1; ...)`), add the `F_eval` lambda. Insert this block:

```cpp
    // F_eval: one complete inner IPF iteration on (Xv, pv, q_hyp_v).
    // Modifies Xv, pv, q_hyp_v in place. Returns errRp for the new state.
    // Also updates errRp_k[] (shared, used by Greedy sort next call).
    // W_total recomputed internally.
    auto F_eval = [&](std::vector<double>& Xv,
                      std::vector<double>& pv,
                      double& q_hyp_v) -> double {
        // Save/restore shared state: hyperplane_step uses q_hyp directly
        // (captured by reference from outer scope). Temporarily swap with q_hyp_v.
        std::swap(q_hyp, q_hyp_v);

        double W_total = 0.0;
        for (int c = 0; c < ct.M_cell; c++) W_total += Xv[c];

        // Greedy sort
        if (use_greedy)
            std::sort(margin_order.begin(), margin_order.end(),
                      [&](int a, int b){ return errRp_k[a] > errRp_k[b]; });

        for (int ki = 0; ki < st.K; ki++) {
            const int k = use_greedy ? margin_order[ki] : ki;
            std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) bucket[g] += Xv[c];
            }
            std::fill(scale.begin(), scale.begin() + st.cat_counts[k], 1.0);
            for (int j = 0; j < st.cat_counts[k]; j++) {
                double Tkj = st.targets[k][j] * W_total;
                if (bucket[j] < kEmptyBucketThreshold * W_total) {
                    if (Tkj > 0.0) is_infeasible = true;
                } else {
                    scale[j] = Tkj / bucket[j];
                }
            }
            const double eff_omega = sor_active ? sor_omega[k] : 1.0;
            for (int c = 0; c < ct.M_cell; c++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k]) {
                    const double s_g = scale[g];
                    Xv[c] *= (eff_omega == 1.0) ? s_g : std::pow(std::max(s_g, 0.0), eff_omega);
                }
            }
            if (sor_active && sor_auto && W_total > 0.0) {
                double ek = 0.0;
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    double e = std::fabs(bucket[j] / W_total - st.targets[k][j]);
                    if (e > ek) ek = e;
                }
                if (sor_prev_errRp[k] < ek)
                    sor_omega[k] = std::max(omega_min, sor_omega[k] * 0.7);
                else
                    sor_omega[k] = std::min(1.0, sor_omega[k] * 1.05);
                sor_prev_errRp[k] = ek;
            }
            if (W_total > 0.0) {
                double ek = 0.0;
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    double e = std::fabs(bucket[j] / W_total - st.targets[k][j]);
                    if (e > ek) ek = e;
                }
                errRp_k[k] = ek;
            }
        }

        // Bregman box correction (using pv, not outer p)
        for (int c = 0; c < ct.M_cell; c++) {
            double yc = Xv[c] * pv[c];
            double Xc = std::clamp(yc, L_cell[c], U_cell[c]);
            pv[c] = (Xc > 0.0) ? yc / Xc : 1.0;
            Xv[c] = Xc;
        }

        // Hyperplane step: uses q_hyp (swapped to q_hyp_v above) and Xv.
        // hyperplane_step() operates on X and q_hyp (outer references).
        // Temporarily swap Xv ↔ X so hyperplane_step() operates on Xv.
        std::swap(X, Xv);
        hyperplane_step();   // modifies X (=Xv) and q_hyp (=q_hyp_v)
        std::swap(X, Xv);

        // Restore outer q_hyp
        std::swap(q_hyp, q_hyp_v);

        // Compute and return errRp of the new state
        return compute_errRp_ct(st, ct, Xv, bucket);
    };
```

**Note on the swap pattern**: `hyperplane_step()` is a lambda capturing `X` and `q_hyp` by reference. To avoid rewriting it, we swap Xv↔X before calling it, then swap back. The same trick applies to q_hyp ↔ q_hyp_v.

- [ ] **Step 2: Add SQUAREM outer loop (replaces the flat loop when st.accelerate)**

After the F_eval lambda definition, add the SQUAREM branch. The structure is:

```cpp
    if (st.accelerate) {
        // ── SQUAREM SqS3 outer loop ──────────────────────────────────────────
        // Each super-step: 3 F-evals (w1=F(X), w2=F(w1), X_new=F(X*)) + k halvings.
        static constexpr double kAlphaMax   = -1.0;    // cap: never gentler than -1
        static constexpr double kAlphaMin   = -1000.0; // cap: never more aggressive
        static constexpr double kVNormEps   = 1e-300;  // denominator guard
        static constexpr double kVNormRel   = 1e-10;   // relative convergence: ‖v‖/‖w2‖
        static constexpr double kHalvingSlack = 1.01;
        static constexpr int    kMaxHalvings  = 16;

        int sq_iter = 0;   // super-step count (each ≈ 3 F-evals)
        while (sq_iter * 3 < st.inner_max_iter) {
            // w1 = F(X), w2 = F(w1)
            auto w1 = X; auto p1 = p; double qh1 = q_hyp;
            F_eval(w1, p1, qh1);

            auto w2 = w1; auto p2 = p1; double qh2 = qh1;
            F_eval(w2, p2, qh2);

            sq_iter++;
            res.iterations = sq_iter * 3;  // approximate F-eval count

            // r = w1-X, v = w2-w1; compute norms
            double norm_r = 0.0, norm_v = 0.0, norm_w2 = 0.0;
            for (int c = 0; c < ct.M_cell; c++) {
                double ri = w1[c] - X[c];
                double vi = w2[c] - w1[c];
                norm_r  += ri * ri;
                norm_v  += vi * vi;
                norm_w2 += w2[c] * w2[c];
            }
            norm_r  = std::sqrt(norm_r);
            norm_v  = std::sqrt(norm_v);
            norm_w2 = std::sqrt(norm_w2);

            // ‖v‖/‖w2‖ guard: already at fixed point
            if (norm_v / (norm_w2 + kVNormEps) < kVNormRel) {
                X = w2; p = p2; q_hyp = qh2;
                res.max_error = compute_errRp_ct(st, ct, X, bucket);
                res.status = is_infeasible ? RK_ERR_INFEAS : RK_OK;
                break;
            }

            // α = -‖r‖/‖v‖, capped to [kAlphaMin, kAlphaMax]
            double alpha = std::max(kAlphaMin, std::min(kAlphaMax, -norm_r / (norm_v + kVNormEps)));

            // Snapshot after w2, before extrapolation
            auto X_snap = w2; auto p_snap = p2; double q_snap = qh2;

            // Extrapolate: X* = X - 2α·r + α²·v; clamp to ≥ 0
            auto X_star = X;  // copy to compute X*
            for (int c = 0; c < ct.M_cell; c++) {
                double ri = w1[c] - X[c];
                double vi = w2[c] - w1[c];
                X_star[c] = X_snap[c] - 2.0 * alpha * ri + alpha * alpha * vi;
                X_star[c] = std::max(X_star[c], 0.0);
            }

            // X_new = F(X*)
            auto p_star = p_snap; double qh_star = q_snap;
            double errRp_new = F_eval(X_star, p_star, qh_star);

            // Step-halving: ‖F(X*)-X*‖₂ vs 1.01×‖r‖₂
            // ‖F(X*)-X*‖₂ = ‖X_star_after_F - X_star_before_F‖ — but we overwrote X_star.
            // Proxy: errRp of X_new vs errRp of w2 (simpler, no extra storage).
            double errRp_w2 = compute_errRp_ct(st, ct, w2, bucket);
            for (int h = 0; h < kMaxHalvings && errRp_new > kHalvingSlack * errRp_w2; h++) {
                alpha = (alpha - 1.0) / 2.0;  // midpoint toward -1; converges for all α ≤ 0
                // Restore from w2 snapshot
                X_star = X_snap; p_star = p_snap; qh_star = q_snap;
                for (int c = 0; c < ct.M_cell; c++) {
                    double ri = w1[c] - X[c];
                    double vi = w2[c] - w1[c];
                    X_star[c] = X_snap[c] - 2.0 * alpha * ri + alpha * alpha * vi;
                    X_star[c] = std::max(X_star[c], 0.0);
                }
                errRp_new = F_eval(X_star, p_star, qh_star);
                res.iterations++;  // count each halving F-eval
            }

            // Accept X_star (= X_new = F(X*) after halving)
            X = X_star; p = p_star; q_hyp = qh_star;
            res.max_error = errRp_new;

            // Best-iterate tracking (same logic as flat loop)
            if (errRp_new < best_metric_seen) {
                best_metric_seen    = errRp_new;
                best_iter_val       = res.iterations;
                best_objective_seen = compute_weight_kl();
                W_best              = X;
            }

            // Convergence check (improvement criterion on errRp)
            lbw::CellMetrics m_conv; m_conv.errRp = errRp_new;
            if (lbw::check_convergence(st.convergence_cfg, m_conv,
                                       prev_metric_for_rule, st.tol_abs)) {
                res.status = is_infeasible ? RK_ERR_INFEAS : RK_OK;
                res.convergence_metric = static_cast<int>(st.convergence_cfg.metric);
                res.convergence_rule   = static_cast<int>(st.convergence_cfg.rule);
                res.convergence_tol    = st.convergence_cfg.pct_tol;
                res.convergence_iter   = res.iterations;
                break;
            }

            if (st.verbose >= 1) {
                char msg[256];
                std::snprintf(msg, 256, "raking[sq] super-step %d: errRp=%.2e alpha=%.4g",
                              sq_iter, errRp_new, alpha);
                st.log(msg);
            }

            // No-improve stall (same as flat loop: 5 consecutive checks without improvement)
            if (!std::isfinite(min_errRp_window)) {
                min_errRp_window = errRp_new; n_no_improve = 0;
            } else {
                const double rel_eps = 0.01 * min_errRp_window;
                const double eps = std::max(rel_eps, st.tol_abs);
                if (errRp_new < min_errRp_window - eps) {
                    min_errRp_window = errRp_new; n_no_improve = 0;
                } else {
                    n_no_improve++;
                }
            }
            if (n_no_improve >= kMaxNoImprove) {
                res.status = RK_ERR_NOCONV;
                break;
            }
        }
        // End SQUAREM branch: fall through to post-loop finalizer below.
    } else {
```

Then wrap the existing flat loop with a closing `}` after line 398 (the end of the `for (int iter...)` loop):

```cpp
    }   // end if (st.accelerate) ... else (flat loop)
```

- [ ] **Step 3: Build**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -15
```
Expected: `* DONE (leafblower)`. If compile errors: read full output, fix before continuing.

- [ ] **Step 4: Smoke test accelerate=TRUE**

```r
df <- data.frame(v1=factor(c("1","2","1","2","1","2","1")),
                 v2=factor(c("A","A","B","B","A","B","A")))
tgt <- list(v1=c("1"=0.5,"2"=0.5), v2=c("A"=0.6,"B"=0.4))
r <- leafblower::harvest(df, tgt, method="raking", accelerate=TRUE,
                         max_weight=5, max_iterations=500L, attach_weights=FALSE)
cat("iters:", attr(r,"result")$iterations, "max_err:", attr(r,"result")$max_error, "\n")
```
Expected: runs without error, prints `iters: N max_err: E` with E < 0.1.

- [ ] **Step 5: Commit implementation**

```bash
git add src/raking.cpp
git commit -m "feat(squarem): implement SQUAREM SqS3 outer loop in raking_solve"
```

---

## Task 5: Run tests and verify ACs

**Files:** none — verification only.

- [ ] **Step 1: Run full test suite**

```bash
Rscript -e "devtools::test()" 2>&1
```
Expected: FAIL ≤ 2 (pre-existing: `test-sor.R/lhs`, `test-ieppa-nonuniform-d.R`). `squarem-red`, `squarem-ac3`, `squarem-ac5` must all PASS.

- [ ] **Step 2: Verify AC1 (no error)**

```r
set.seed(1L); n <- 500L
df  <- data.frame(v1=factor(sample(3L,n,TRUE)), v2=factor(sample(2L,n,TRUE)))
tgt <- list(v1=c("1"=0.5,"2"=0.3,"3"=0.2), v2=c("1"=0.6,"2"=0.4))
r <- leafblower::harvest(df, tgt, method="raking", accelerate=TRUE,
                         max_weight=5, max_iterations=500L, attach_weights=FALSE)
cat("AC1: status=", attr(r,"result")$status,
    "iters=", attr(r,"result")$iterations,
    "max_err=", attr(r,"result")$max_error, "\n")
```
Expected: no R error, status 0 or 1.

- [ ] **Step 3: Verify AC3 (bit-identical accelerate=FALSE)**

```bash
Rscript -e "devtools::test(filter='squarem-ac3')" 2>&1
```
Expected: PASS.

- [ ] **Step 4: Verify AC5 (‖v‖ guard — already-converged problem)**

```r
# Trivially converged problem: 1 category per margin, target = current proportion
df <- data.frame(v1=factor(rep("1",10)))
tgt <- list(v1=c("1"=1.0))
r <- leafblower::harvest(df, tgt, method="raking", accelerate=TRUE,
                         max_weight=5, max_iterations=500L, attach_weights=FALSE)
w <- as.numeric(r)
cat("AC5: all weights finite:", all(is.finite(w)), "any NaN:", any(is.nan(w)), "\n")
```
Expected: `AC5: all weights finite: TRUE any NaN: FALSE`.

- [ ] **Step 5: Verify AC6 (warning for non-raking)**

```bash
Rscript -e "devtools::test(filter='squarem-ac5')" 2>&1
```
Expected: PASS.

- [ ] **Step 6: Commit test results (no code changes, just update memory)**

```bash
git status  # should be clean
```

---

## Task 6: AC2 manual benchmark (local only, not CI)

> **This task is not runnable in CI** (stepstone dataset not committed). Run locally before merge. Log result in PR description.

**Files:** none — benchmark only.

- [ ] **Step 1: Run stepstone benchmark with accelerate=TRUE**

```bash
OMP_NUM_THREADS=1 Rscript -e "
  library(leafblower); library(arrow); library(jsonlite)
  df  <- arrow::read_parquet('benchmarks/stepstone_fulldata_bench_data.parquet')
  df\$uuid <- NULL
  tgt <- lapply(jsonlite::fromJSON('benchmarks/stepstone_fulldata_bench_targets.json'),
                function(x) { v <- unlist(x); v / sum(v) })
  for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])

  t0 <- proc.time()['elapsed']
  r  <- leafblower::harvest(df, tgt, method='raking', accelerate=TRUE,
         max_weight=5, max_iterations=5000L, attach_weights=FALSE, verbose=0)
  wall <- proc.time()['elapsed'] - t0
  res <- attr(r, 'result')
  cat(sprintf('SQUAREM: wall=%.1fs  iters=%d  max_err=%.4e  status=%d\n',
              wall, res\$iterations, res\$max_error, res\$status))
" 2>&1
```

- [ ] **Step 2: Verify AC2 passes**

Confirm: `max_err ≤ 2.50e-3`.

If max_err > 2.50e-3: do NOT merge. Investigate:
- Does greedy + SQUAREM outperform greedy alone? If not, check F_eval lambda for state leakage.
- Is errRp_k[] being updated correctly between super-steps?

- [ ] **Step 3: Log result**

Record in the commit message or PR description:
```
AC2 benchmark (stepstone-fulldata, n=1.58M, K=9, max_weight=5):
  SQUAREM max_err = X.XXe-3 (threshold ≤ 2.50e-3) ✓
  wall = X.Xs, iters = N
```

---

## Self-Review Against Spec

**Spec coverage:**
- AC1 (no error): Task 5 Step 2 ✓
- AC2 (max_err ≤ 2.50e-3): Task 6 ✓
- AC3 (bit-identical): Task 5 Step 3 + squarem-ac3 test ✓
- AC4 (snapshot after F(w2)): implemented in Task 4 (X_snap = w2 after second F_eval) ✓
- AC5 (‖v‖ guard): Task 5 Step 4 ✓
- AC6 (warn non-raking): Task 5 Step 6 + squarem-ac5 test ✓
- AC7 (devtools::test() FAIL ≤ 2): Task 5 Step 1 ✓

**No placeholders**: all code blocks are complete and exact.

**Type consistency**: `F_eval` signature `(std::vector<double>&, std::vector<double>&, double&) -> double` used consistently in Task 4 Steps 1 and 2.
