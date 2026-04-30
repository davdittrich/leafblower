# Critical Code Review Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix all issues identified in the 2026-04-28 critical code review: two blocking SQUAREM correctness bugs, one blocking docstring mismatch, and seven required/suggested fixes across `src/raking.cpp`, `R/harvest.R`, `src/r_bridge.cpp`, and `src/types.hpp`.

**Architecture:** Fixes are grouped by compile dependency: raking.cpp changes first (three tasks, three compile gates), then pure-R harvest.R changes (two tasks), then r_bridge.cpp changes (one compile gate), then types.hpp comment. Each task is independently releasable.

**Tech Stack:** C++17 (`src/raking.cpp`, `src/r_bridge.cpp`, `src/types.hpp`), R (`R/harvest.R`), testthat.

---

## Issue Index (for traceability)

| ID | File | Lines | Description |
|----|------|-------|-------------|
| C1 | `src/raking.cpp` | 367, 439 | SQUAREM convergence exits return INFEAS when converging |
| C2 | `src/raking.cpp` | 436 | SQUAREM `m_conv` only populates `errRp` — kl/chi2/etc = 0 |
| C3 | `src/raking.cpp` | 54–58 | Docstring contradicts itself — says Dykstra IS and IS NOT used |
| R1 | `R/harvest.R` | 9 | `@param method` says "IPF + Dykstra box projection" (stale) |
| R2 | `R/harvest.R` | 343–345 | SQUAREM stall warning says "errRp plateau" (wrong) |
| R3 | `R/harvest.R` | 306 | `convergence_reason = "stall_errRp"` for SQUAREM (wrong) |
| R4 | `src/r_bridge.cpp` | 235–241 | `p.algorithm` missing sinkhorn → defaults to RK_ALG_IEPPA |
| R5 | `R/harvest.R` | 155 | `design_weights` param undocumented |
| R6 | `R/harvest.R` | 382–386 | `"algorithm"` attr missing when `attach_weights=FALSE` |
| S1 | `R/harvest.R` | 352–365 | PCT stall detection fires on intentional l1_weight convergence |
| S2 | `src/r_bridge.cpp` | 456–516 | `ieppa_soft` buried inside chebyshev/grake else block |
| S3 | `src/r_bridge.cpp` | 649–661 | Unprotected SEXPs in test function |
| S4 | `src/raking.cpp` | 337 | Misleading `(void)errRp_w1` comment |
| S5 | `src/types.hpp` | 68 | `burnin` dead for raking, not documented |
| S6 | `src/raking.cpp` | 60, 275 | `kEmptyBucketThreshold` used as both absolute and relative |

---

## File Map

| File | Tasks |
|------|-------|
| `tests/testthat/test-calibration-solvers.R` | Task 0: RED tests |
| `src/raking.cpp` | Task 1 (C1), Task 2 (C2), Task 3 (C3+S4+S6) |
| `R/harvest.R` | Task 4 (R2+R3), Task 5 (R1+R5+R6+S1) |
| `src/r_bridge.cpp` | Task 6 (R4+S2+S3) |
| `src/types.hpp` | Task 7 (S5) |

---

## Task 0: Write RED Tests (Before Any Code)

**Files:**
- Modify: `tests/testthat/test-calibration-solvers.R` (append)

- [ ] **Step 1: Append RED test at end of test file**

Only squarem-c2 is a genuine RED test (fails deterministically pre-fix). squarem-c1 is a regression guard added in Task 1 after the fix.

```r
# ── Critical review fix tests ─────────────────────────────────────────────────

test_that("squarem-c2: kl metric with accelerate=TRUE runs correct convergence", {
  # C2 RED: m_conv only sets errRp; kl=0.0 (default) → check_convergence fires on
  # iter 1-3 via IMPROVEMENT rule (`curr <= 1e-15` → trivially converged).
  # C2 GREEN: compute_cell_metrics populates kl correctly → proper kl convergence.
  set.seed(42L)
  df <- data.frame(v1 = factor(c(rep("A", 80L), rep("B", 20L))))
  tgt <- list(v1 = c("A" = 0.5, "B" = 0.5))  # imbalanced input → kl > 0 at start

  w <- suppressWarnings(
    leafblower::harvest(df, tgt, method = "raking", accelerate = TRUE,
                        max_weight = 5, max_iterations = 200L,
                        convergence = list(metric = "kl", rule = "improvement", tol = 0.001),
                        attach_weights = FALSE))
  r <- attr(w, "result")

  # RED: r$iterations <= 3 (kl=0 default causes instant convergence at first super-step)
  # GREEN: r$iterations > 3 (real kl convergence runs multiple super-steps)
  expect_gt(r$iterations, 3L,
            label = "C2: SQUAREM kl-metric must run >3 F-evals (not fire on kl=0 default)")
})
```

- [ ] **Step 2: Verify squarem-c2 is RED**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-calibration-solvers.R')" 2>&1 | grep -E "squarem-c2|PASS|FAIL" | tail -5
```

Expected: `squarem-c2` FAIL — `r$iterations <= 3` (kl=0 causes instant convergence pre-fix).

- [ ] **Step 3: Run full suite — verify existing tests still pass (syntax gate)**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: FAIL ≤ 4 (pre-existing 3 + squarem-c2 RED). Confirms the appended test is syntactically valid R.

- [ ] **Step 4: Commit RED test**

```bash
git add tests/testthat/test-calibration-solvers.R
git commit -m "test(squarem-c2): RED test — kl metric with accelerate=TRUE must not fire on kl=0 default"
```

---

## Task 1: Fix C1 — SQUAREM Convergence Exits Return RK_OK Unconditionally

**Files:**
- Modify: `src/raking.cpp` (lines 367, 439)

**Context:** The flat loop at line 575 correctly uses `res.status = RK_OK` unconditionally after the same refactor. Two SQUAREM convergence exits were missed.

- [ ] **Step 1: Fix fixed-point guard convergence exit (line 367)**

Find (inside `if (st.accelerate)` → `while` loop → fixed-point guard block):

```cpp
                if (norm_v / (norm_w2 + kVNormEps) < kVNormRel) {
                    X = w2;
                    res.max_error        = errRp_w2;
                    res.status           = is_infeasible ? RK_ERR_INFEAS : RK_OK;
                    res.convergence_iter = f_eval_count;
                    break;
                }
```

Replace with:

```cpp
                if (norm_v / (norm_w2 + kVNormEps) < kVNormRel) {
                    X = w2;
                    res.max_error        = errRp_w2;
                    res.status           = RK_OK;
                    res.convergence_iter = f_eval_count;
                    break;
                }
```

- [ ] **Step 2: Fix check_convergence break convergence exit (line 439)**

Find (inside same SQUAREM while loop, after `// Convergence criterion`):

```cpp
                if (lbw::check_convergence(st.convergence_cfg, m_conv,
                                           prev_metric_for_rule, st.tol_abs)) {
                    res.status             = is_infeasible ? RK_ERR_INFEAS : RK_OK;
                    res.convergence_metric = static_cast<int>(st.convergence_cfg.metric);
```

Replace `res.status` line only:

```cpp
                if (lbw::check_convergence(st.convergence_cfg, m_conv,
                                           prev_metric_for_rule, st.tol_abs)) {
                    res.status             = RK_OK;
                    res.convergence_metric = static_cast<int>(st.convergence_cfg.metric);
```

- [ ] **Step 3: Build**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 4: Add squarem-c1 regression guard to test file**

squarem-c1 is a regression guard (not a TDD RED test — C1 is hard to trigger deterministically). Add it now, after the fix, to lock in correct behavior:

```r
test_that("squarem-c1: feasible tight-bounds SQUAREM must not return status=INFEAS", {
  # Regression guard for C1 fix: SQUAREM convergence exits use unconditional RK_OK.
  # Water-fill may transiently set is_infeasible=true on extrapolated iterates;
  # this must never override convergence status on a genuinely feasible problem.
  set.seed(99L)
  df <- data.frame(
    v1 = factor(c(rep("A", 80L), rep("B", 20L))),
    v2 = factor(sample(c("x", "y"), 100L, TRUE))
  )
  tgt <- list(v1 = c("A" = 0.3, "B" = 0.7), v2 = c("x" = 0.5, "y" = 0.5))

  w_flat <- suppressWarnings(
    leafblower::harvest(df, tgt, method = "raking", accelerate = FALSE,
                        max_weight = 4, max_iterations = 500L, attach_weights = FALSE))
  skip_if(attr(w_flat, "result")$status == 2L,
          "flat raking reports INFEAS — problem is genuinely infeasible; choose different seed")

  expect_no_error(
    w_sq <- suppressWarnings(
      leafblower::harvest(df, tgt, method = "raking", accelerate = TRUE,
                          max_weight = 4, max_iterations = 500L, attach_weights = FALSE)),
    message = "C1: SQUAREM must not crash on a feasible tight-bounds problem"
  )
  expect_false(attr(w_sq, "result")$status == 2L,
               label = "C1: SQUAREM status must not be INFEAS=2 for a feasible problem")
})
```

Run to confirm PASS post-fix:

```bash
Rscript -e "testthat::test_file('tests/testthat/test-calibration-solvers.R')" 2>&1 | grep -E "squarem-c[12]|PASS|FAIL" | tail -5
```

Expected: `squarem-c1` PASS (or SKIP), `squarem-c2` PASS.

- [ ] **Step 5: Run full test suite — verify no regressions**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: FAIL ≤ 3 (pre-existing).

- [ ] **Step 6: Commit**

```bash
git add src/raking.cpp tests/testthat/test-calibration-solvers.R
git commit -m "fix(squarem): convergence exits use unconditional RK_OK — transient is_infeasible must not override converged status"
```

---

## Task 2: Fix C2 — SQUAREM `m_conv` Full Metric Propagation

**Files:**
- Modify: `src/raking.cpp` (around line 436)

**PREREQUISITE: Task 1 must be applied first.** The "Before" block below shows the post-Task-1 state (`res.status = RK_OK`). If Task 1 has not run, the actual file contains `res.status = is_infeasible ? RK_ERR_INFEAS : RK_OK` at that line. Do not execute this task in isolation — run Task 1 and its compile gate before starting here.

**Context:** The flat loop at lines 566–569 populates all CellMetrics fields. SQUAREM only set `errRp`. `compute_cell_metrics()` is already available from `calib_dispatch.hpp` (included). The weight-change `wchange` (obs-level L1) is moved to before the convergence check so it can populate `m_conv.l1`.

**Design note — why `compute_cell_metrics()` instead of a targeted single-metric compute:**
`compute_cell_metrics` is a single O(K×M_cell) pass that computes all 5 metrics. A targeted approach (switch on `cfg.metric`, compute only the needed field) would require ~30 lines of duplicated bucket-aggregation code per metric. The overhead is negligible: each SQUAREM super-step already performs 3+ full F_evals (each O(K×M_cell)), making the extra metrics pass <5% of super-step cost. The flat loop already calls the same comprehensive metrics block unconditionally (lines 511–551). Consistency with the flat loop and avoiding code duplication outweigh the marginal cost. No simpler targeted approach produces a meaningfully faster result at SQUAREM's granularity.

- [ ] **Step 1: Move wchange computation to before check_convergence and wire into m_conv**

Find the block after `res.iterations = f_eval_count;` in the SQUAREM while loop (currently around line 435–473). Replace the entire `// Convergence criterion` through `if (!fell_back) X_prev_sq = X;` block:

**Before (find this exact block):**

```cpp
                // Convergence criterion
                lbw::CellMetrics m_conv; m_conv.errRp = errRp_new;
                if (lbw::check_convergence(st.convergence_cfg, m_conv,
                                           prev_metric_for_rule, st.tol_abs)) {
                    res.status             = RK_OK;
                    res.convergence_metric = static_cast<int>(st.convergence_cfg.metric);
                    res.convergence_rule   = static_cast<int>(st.convergence_cfg.rule);
                    res.convergence_tol    = st.convergence_cfg.pct_tol;
                    res.convergence_iter   = f_eval_count;
                    break;
                }

                if (st.verbose >= 1) {
                    char msg[256];
                    std::snprintf(msg, 256, "raking[sq] f_eval=%d errRp=%.2e alpha=%.4g",
                                  f_eval_count, errRp_new, alpha);
                    st.log(msg);
                }

                // Weight-change stall: obs-level L1 Δw goes to zero at the fixed point.
                // Tried KL stall after geometry fix — gives identical result (81 F-evals,
                // same max_err), confirming accepted iterate KL is now approximately monotone.
                // Weight-change kept: equivalent result, more robust (no log(0) risk).
                // Skip snapshot update on fell_back: X=w2 → X_prev_sq=w2 → wchange=0
                // next iter → spurious stall after 5 consecutive fell_back super-steps.
                double wchange = 0.0;
                for (int c = 0; c < ct.M_cell; c++)
                    wchange += std::fabs(X[c] - X_prev_sq[c]) / static_cast<double>(ct.n_per_cell[c]);
                wchange /= static_cast<double>(st.n);

                if (!std::isfinite(min_loss_window)) {
                    min_loss_window = wchange; n_no_improve = 0;
                } else if (wchange < min_loss_window * (1.0 - st.convergence_cfg.pct_tol)) {
                    min_loss_window = wchange; n_no_improve = 0;
                } else {
                    n_no_improve++;
                }
                if (n_no_improve >= kMaxNoImprove) { res.status = RK_ERR_STALL; break; }

                if (!fell_back) X_prev_sq = X;
```

**After (replace with):**

```cpp
                // Weight-change (obs-level L1): computed here so m_conv.l1 can use it
                // for metric="l1_weight" convergence. Also used for stall detection below.
                // Skip snapshot update on fell_back (prevents wchange=0 spurious stall).
                double wchange = 0.0;
                for (int c = 0; c < ct.M_cell; c++)
                    wchange += std::fabs(X[c] - X_prev_sq[c]) / static_cast<double>(ct.n_per_cell[c]);
                wchange /= static_cast<double>(st.n);

                // Convergence criterion: compute all metrics so non-MAX_ERR metrics work correctly.
                // compute_cell_metrics fills errRp/mean_err/kl/chi2/grake_norm; override errRp
                // with F_eval's result (computed post-hyperplane-normalization, slightly more
                // accurate than re-computing from X). l1 = wchange (obs-level L1 Δw this step).
                {
                    double W_sq = 0.0;
                    for (int c = 0; c < ct.M_cell; c++) W_sq += X[c];
                    lbw::CellMetrics m_conv = lbw::compute_cell_metrics(st, ct, X, W_sq, bucket);
                    m_conv.errRp = errRp_new;
                    m_conv.l1    = wchange;
                    if (lbw::check_convergence(st.convergence_cfg, m_conv,
                                               prev_metric_for_rule, st.tol_abs)) {
                        res.status             = RK_OK;
                        res.convergence_metric = static_cast<int>(st.convergence_cfg.metric);
                        res.convergence_rule   = static_cast<int>(st.convergence_cfg.rule);
                        res.convergence_tol    = st.convergence_cfg.pct_tol;
                        res.convergence_iter   = f_eval_count;
                        break;
                    }
                }

                if (st.verbose >= 1) {
                    char msg[256];
                    std::snprintf(msg, 256, "raking[sq] f_eval=%d errRp=%.2e alpha=%.4g",
                                  f_eval_count, errRp_new, alpha);
                    st.log(msg);
                }

                // Weight-change stall: wchange computed above; reuse here.
                if (!std::isfinite(min_loss_window)) {
                    min_loss_window = wchange; n_no_improve = 0;
                } else if (wchange < min_loss_window * (1.0 - st.convergence_cfg.pct_tol)) {
                    min_loss_window = wchange; n_no_improve = 0;
                } else {
                    n_no_improve++;
                }
                if (n_no_improve >= kMaxNoImprove) { res.status = RK_ERR_STALL; break; }

                if (!fell_back) X_prev_sq = X;
```

- [ ] **Step 2: Build**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 3: Run C2 test**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-calibration-solvers.R')" 2>&1 | grep -E "squarem-c[12]|PASS|FAIL" | tail -5
```

Expected: both `squarem-c1` and `squarem-c2` PASS.

- [ ] **Step 4: Run full test suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: FAIL ≤ 3 (pre-existing).

- [ ] **Step 5: Commit**

```bash
git add src/raking.cpp
git commit -m "fix(squarem): populate all CellMetrics in m_conv before check_convergence — kl/chi2/grake_norm/l1 metrics now work correctly with accelerate=TRUE"
```

---

## Task 3: raking.cpp Cleanup — C3 Docstring + S4 Comment + S6 Constants

**Files:**
- Modify: `src/raking.cpp` (lines 54–58, 60, 275, 337)

No behavior change. Documentation and constant naming only.

- [ ] **Step 1: Fix C3 — rewrite stale function docstring (lines 54–58)**

Find:

```cpp
// Constrained raking solver: cyclic IPF for marginal projections + Dykstra box correction.
// Marginal step: pure IPF (Bregman/multiplicative projection — Euclidean Dykstra corrections
// diverge on multiplicative projections and are not used here).
// Box step: Dykstra additive correction q[i] prevents cycling at the [lo,hi]^n boundary.
// inner_max_iter is the single iteration budget; outer_max_iter is unused.
```

Replace with:

```cpp
// Constrained raking solver: cyclic IPF with per-category water-filling box projections.
// Each margin step applies water_fill_cat() — KL projection onto
// {Σ_c Xv[c]=T_kj, L_c≤Xv[c]≤U_c} (Csiszar-Tusnady 1984; autumn single_adjust).
// No Dykstra correction vectors; stateless F enables SQUAREM L2 step-halving.
// inner_max_iter is the single iteration budget; outer_max_iter is unused.
```

- [ ] **Step 2: Fix S6 — split `kEmptyBucketThreshold` into two named constants (line 60)**

Find (at top of `raking_solve`, with the other static constexpr declarations):

```cpp
    static constexpr double kEmptyBucketThreshold = 1e-15;
```

Replace with:

```cpp
    static constexpr double kAbsoluteZeroThreshold = 1e-15;  // bucket_j / free_sum is genuinely zero
    static constexpr double kRelativeZeroFraction  = 1e-15;  // bucket[j] < fraction * W_total
```

- [ ] **Step 3: Update kEmptyBucketThreshold usages in water_fill_cat**

Find and replace all three occurrences inside `water_fill_cat` lambda:

```cpp
// Occurrence 1 (line ~175):
        if (bucket_j < kEmptyBucketThreshold) {
// Replace:
        if (bucket_j < kAbsoluteZeroThreshold) {

// Occurrence 2 (line ~189 in pass loop):
            if (free_sum < kEmptyBucketThreshold) { is_infeasible = true; break; }
// Replace:
            if (free_sum < kAbsoluteZeroThreshold) { is_infeasible = true; break; }

// Occurrence 3 (line ~225 in post-loop fallback):
        const double m_final = (free_sum > kEmptyBucketThreshold && T_final > 0.0)
// Replace:
        const double m_final = (free_sum > kAbsoluteZeroThreshold && T_final > 0.0)
```

- [ ] **Step 4: Update kEmptyBucketThreshold usage in F_eval (line ~275)**

Find inside `F_eval` lambda, inside the margin sweep:

```cpp
                if (bucket[j] < kEmptyBucketThreshold * W_total) {
```

Replace:

```cpp
                if (bucket[j] < kRelativeZeroFraction * W_total) {
```

- [ ] **Step 5: Fix S4 — misleading `(void)errRp_w1` comment (line ~337)**

Find:

```cpp
                (void)errRp_w1;  // advances IPF side effects (errRp_k); value unused
```

Replace:

```cpp
                (void)errRp_w1;  // errRp_k updated inside F_eval but not consumed (use_greedy=false when accelerate=true)
```

- [ ] **Step 6: Build**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 7: Run full test suite — no regressions**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: FAIL ≤ 3 (pre-existing).

- [ ] **Step 8: Commit**

```bash
git add src/raking.cpp
git commit -m "docs(raking): fix stale Dykstra docstring; split kEmptyBucketThreshold into absolute/relative; clarify errRp_w1 comment"
```

---

## Task 4: harvest.R Semantic Fixes — R2 Warning + R3 convergence_reason

**Files:**
- Modify: `R/harvest.R` (lines 110–112, 306, 343–345)

Pure R changes. No recompile.

- [ ] **Step 1: Fix R3 — update `convergence_reason` from "stall_errRp" to "stall_wchange" (line 306)**

Find inside `calib_result$convergence_used <- list(...)` block:

```r
      else if (s == 5L && isTRUE(accelerate_bool))  "stall_errRp"
```

Replace:

```r
      else if (s == 5L && isTRUE(accelerate_bool))  "stall_wchange"
```

- [ ] **Step 2: Fix R3 — update roxygen @return doc for convergence_reason (lines 110–112)**

Find:

```r
#'           \code{"stall_kl"} (weight KL plateau — at constrained KL minimum),
#'           \code{"stall_errRp"} (SQUAREM errRp plateau — empirically near optimum),
```

Replace:

```r
#'           \code{"stall_kl"} (weight KL plateau — at constrained KL minimum),
#'           \code{"stall_wchange"} (SQUAREM weight-change plateau — at constrained optimum),
```

- [ ] **Step 3: Grep tests for any assertion on the old "stall_errRp" string**

```bash
grep -rn "stall_errRp\|stall_kl.*accelerate\|squarem-geo-ac1.*convergence_reason" \
     tests/testthat/ 2>&1
```

Expected: NO match for `stall_errRp` (confirming no test asserts the old string). The
`status-stall` test checks `"stall_kl"` but uses `accelerate=FALSE` — unaffected by R3.
`squarem-geo-ac1` checks only `max_error` (not `convergence_reason`) — also unaffected.

If ANY test asserts `"stall_errRp"`: update that assertion to `"stall_wchange"` before
proceeding, and commit with: `git commit -m "test: update stall_errRp assertion to stall_wchange (R3 rename)"`

- [ ] **Step 4: Fix R2 — rewrite incorrect SQUAREM stall warning (lines 343–345)**

Find:

```r
  if (calib_result$status == 5L && isTRUE(accelerate_bool))
    warning("leafblower: SQUAREM errRp plateau — weights are valid; ",
            "try accelerate=FALSE for KL-stall (reaches constrained KL minimum)")
```

Replace:

```r
  if (calib_result$status == 5L && isTRUE(accelerate_bool))
    warning("leafblower: SQUAREM weight-change plateau — at constrained optimum; ",
            "weights are valid; no further improvement is achievable")
```

- [ ] **Step 5: Verify convergence_reason for accelerate=TRUE**

```bash
Rscript -e "
  set.seed(42L)
  df  <- data.frame(v1=factor(c(rep('A',100L),rep('B',5L))),
                    v2=factor(sample(2L,105L,TRUE)))
  tgt <- list(v1=c('A'=0.9,'B'=0.1),v2=c('1'=0.5,'2'=0.5))
  w <- suppressWarnings(leafblower::harvest(df,tgt,method='raking',accelerate=TRUE,
         max_weight=5,max_iterations=50L,attach_weights=FALSE))
  r <- attr(w,'result')
  cat('reason:', r\$convergence_used\$convergence_reason, '\n')
  cat('status:', r\$status, '\n')
" 2>&1 | grep -E "^reason:|^status:"
```

Expected: `reason: stall_wchange` (or `criterion` if it converges without stalling).

- [ ] **Step 6: Run full test suite — verify squarem tests still pass**

```bash
Rscript -e "devtools::test()" 2>&1 | grep -E "squarem|FAIL|PASS" | tail -10
```

Expected: squarem-geo-smoke PASS; FAIL ≤ 3 (pre-existing).

- [ ] **Step 7: Commit**

```bash
git add R/harvest.R
git commit -m "fix(harvest): SQUAREM convergence_reason stall_errRp->stall_wchange; fix stall warning text (not errRp, not worse than flat)"
```

---

## Task 5: harvest.R Doc/API Fixes — R1 + R5 + R6 + S1

**Files:**
- Modify: `R/harvest.R` (lines 9, ~84, 352–365, 382–386)

Pure R changes. No recompile. No behavior change except R6 (adds missing attribute).

- [ ] **Step 1: Fix R1 — update `@param method` doc (line 9)**

Find:

```r
#'   \code{"raking"} (IPF + Dykstra box projection), \code{"lbfgsb"}
```

Replace:

```r
#'   \code{"raking"} (IPF + water-filling box projection, KL projection per Csiszar-Tusnady 1984), \code{"lbfgsb"}
```

- [ ] **Step 2: Fix R5 — add `@param design_weights` documentation**

Find (the last `@param` before `@return`):

```r
#' @param ... Additional arguments ignored.
#' @return data frame with weights column if \code{attach_weights=TRUE},
```

Replace:

```r
#' @param design_weights Optional design weights vector. When non-NULL and
#'   \code{start_weights} is NULL, used as starting weights (normalized to
#'   mean=1). Length must equal \code{nrow(data)}.
#' @param ... Additional arguments ignored.
#' @return data frame with weights column if \code{attach_weights=TRUE},
```

- [ ] **Step 3: Fix R6 — add `"algorithm"` attribute when `attach_weights=FALSE` (lines 382–386)**

Find:

```r
  if (!attach_weights) {
    attr(weights, "result") <- calib_result
    attr(weights, "iterations") <- calib_result$iterations
    return(weights)
  }
```

Replace:

```r
  if (!attach_weights) {
    attr(weights, "result")     <- calib_result
    attr(weights, "algorithm")  <- alg_used
    attr(weights, "iterations") <- calib_result$iterations
    return(weights)
  }
```

- [ ] **Step 4: Fix S1 — guard PCT stall detection to max_err / mean_err metrics only (lines 352–360)**

Find:

```r
  if (calib_result$status == 0L &&
      !is.null(conv$pct_tol) && conv$pct_tol > 0 &&
      !is.null(calib_result$max_error) &&
      calib_result$max_error > 10 * conv$pct_tol) {
```

Replace:

```r
  if (calib_result$status == 0L &&
      conv$metric %in% c("max_err", "mean_err") &&
      !is.null(conv$pct_tol) && conv$pct_tol > 0 &&
      !is.null(calib_result$max_error) &&
      calib_result$max_error > 10 * conv$pct_tol) {
```

- [ ] **Step 5: Add test_that regression test for R6 (algorithm attr symmetry)**

Append to `tests/testthat/test-calibration-solvers.R`:

```r
test_that("r6: algorithm attribute present for both attach_weights=TRUE and FALSE", {
  # R6: before fix, attr(r,"algorithm") was NULL when attach_weights=FALSE.
  # GREEN: both modes return a non-NULL, non-empty algorithm string.
  df  <- data.frame(v1 = factor(c("A", "B", "A")))
  tgt <- list(v1 = c("A" = 0.5, "B" = 0.5))

  w_detach <- leafblower::harvest(df, tgt, attach_weights = FALSE)
  w_attach  <- leafblower::harvest(df, tgt, attach_weights = TRUE)

  expect_false(is.null(attr(w_detach, "algorithm")),
               label = "R6: algorithm attr must be non-NULL when attach_weights=FALSE")
  expect_false(is.null(attr(w_attach,  "algorithm")),
               label = "algorithm attr must be non-NULL when attach_weights=TRUE")
  expect_equal(attr(w_detach, "algorithm"), attr(w_attach, "algorithm"),
               label = "R6: algorithm attr must be identical regardless of attach_weights")
})
```

Then run:

```bash
Rscript -e "testthat::test_file('tests/testthat/test-calibration-solvers.R')" 2>&1 | grep -E "r6|PASS|FAIL" | tail -5
```

Expected: `r6` PASS.

- [ ] **Step 7: Run full test suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: FAIL ≤ 3 (pre-existing).

- [ ] **Step 8: Commit**

```bash
git add R/harvest.R tests/testthat/test-calibration-solvers.R
git commit -m "fix(harvest): fix @param method doc (Dykstra->water-filling); add design_weights doc; add algorithm attr for attach_weights=FALSE; guard PCT stall warning to max_err/mean_err metrics"
```

---

## Task 6: r_bridge.cpp Fixes — R4 + S2 + S3

**Files:**
- Modify: `src/r_bridge.cpp` (lines 235–241, 456–516, 649–661)

Requires one compile gate.

- [ ] **Step 1: Fix R4 — add sinkhorn to `p.algorithm` switch (lines 235–241)**

Find:

```cpp
    else if (strcmp(method_str, "grake")     == 0) p.algorithm = RK_ALG_GRAKE;
    else if (strcmp(method_str, "auto")      == 0) p.algorithm = RK_ALG_AUTO;
    else                                            p.algorithm = RK_ALG_IEPPA;
```

Replace:

```cpp
    else if (strcmp(method_str, "grake")     == 0) p.algorithm = RK_ALG_GRAKE;
    else if (strcmp(method_str, "sinkhorn")  == 0) p.algorithm = RK_ALG_SINKHORN;
    else if (strcmp(method_str, "auto")      == 0) p.algorithm = RK_ALG_AUTO;
    else                                            p.algorithm = RK_ALG_IEPPA;
```

- [ ] **Step 2: Fix S2 — fix comment + indentation for ieppa_soft dispatch block (lines ~456–516)**

Find the comment at the start of the outer else block:

```cpp
    } else {
        // Shared dispatch for both chebyshev and grake (same solver, different variant).
        auto dispatch_cheb = [&](lbw::LpVariant variant, int alg_code) {
```

Replace comment only:

```cpp
    } else {
        // Dispatch for chebyshev, grake (shared solver), ieppa_soft, and default ieppa.
        auto dispatch_cheb = [&](lbw::LpVariant variant, int alg_code) {
```

Then fix the indentation of the `ieppa_soft` and default ieppa `else` blocks. Find (around line 474):

```cpp
        } else if (strcmp(method_str, "ieppa_soft") == 0) {
        st.ieppa_auto_selected = false;
        st.use_admm_capacity   = true;
```

Replace (add 8-space indent to match the chebyshev/grake blocks):

```cpp
        } else if (strcmp(method_str, "ieppa_soft") == 0) {
            st.ieppa_auto_selected = false;
            st.use_admm_capacity   = true;
```

And similarly fix the closing `} else {` block for default ieppa (around line 495):

```cpp
        } else {
        // Default / ieppa
        st.ieppa_auto_selected = (strcmp(method_str, "ieppa") != 0);
```

Replace:

```cpp
        } else {
            // Default / ieppa
            st.ieppa_auto_selected = (strcmp(method_str, "ieppa") != 0);
```

Apply consistent 12-space indentation to all statements inside both blocks to match the surrounding chebyshev/grake dispatch style. (Apply `clang-format` or manual indent to the affected region if easier.)

- [ ] **Step 3: Fix S3 — PROTECT SEXPs in C_leafblower_cell_table_probe (lines 649–661)**

Find inside `C_leafblower_cell_table_probe` (after `SEXP ret = PROTECT(...)`):

```cpp
    SEXP cell_of_sexp = Rf_allocVector(INTSXP, n);
    std::memcpy(INTEGER(cell_of_sexp), ct.cell_of.data(), n * sizeof(int));
    SET_VECTOR_ELT(ret, 1, cell_of_sexp);
    SEXP npc = Rf_allocVector(INTSXP, ct.M_cell);
    std::memcpy(INTEGER(npc), ct.n_per_cell.data(), ct.M_cell * sizeof(int));
    SET_VECTOR_ELT(ret, 2, npc);
    SEXP names = Rf_allocVector(STRSXP, 3);
    SET_STRING_ELT(names, 0, Rf_mkChar("M_cell"));
    SET_STRING_ELT(names, 1, Rf_mkChar("cell_of"));
    SET_STRING_ELT(names, 2, Rf_mkChar("n_per_cell"));
    Rf_setAttrib(ret, R_NamesSymbol, names);
    UNPROTECT(1);
    return ret;
```

Replace:

```cpp
    SEXP cell_of_sexp = PROTECT(Rf_allocVector(INTSXP, n));
    std::memcpy(INTEGER(cell_of_sexp), ct.cell_of.data(), n * sizeof(int));
    SET_VECTOR_ELT(ret, 1, cell_of_sexp);
    SEXP npc = PROTECT(Rf_allocVector(INTSXP, ct.M_cell));
    std::memcpy(INTEGER(npc), ct.n_per_cell.data(), ct.M_cell * sizeof(int));
    SET_VECTOR_ELT(ret, 2, npc);
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, Rf_mkChar("M_cell"));
    SET_STRING_ELT(names, 1, Rf_mkChar("cell_of"));
    SET_STRING_ELT(names, 2, Rf_mkChar("n_per_cell"));
    Rf_setAttrib(ret, R_NamesSymbol, names);
    UNPROTECT(4);  // ret + cell_of_sexp + npc + names
    return ret;
```

- [ ] **Step 4: Build**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 5: Smoke test sinkhorn p.algorithm**

```bash
Rscript -e "
  df  <- data.frame(v1=factor(c('A','B','A','B','A')))
  tgt <- list(v1=c('A'=0.5,'B'=0.5))
  suppressWarnings(w <- leafblower::harvest(df,tgt,method='sinkhorn',attach_weights=FALSE))
  r <- attr(w,'result')
  cat('algorithm_used:', r\$algorithm_used, '(expected 4 = RK_ALG_SINKHORN)\n')
" 2>&1 | grep algorithm_used
```

Expected: `algorithm_used: 4 (expected 4 = RK_ALG_SINKHORN)`

- [ ] **Step 6: Run full test suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: FAIL ≤ 3 (pre-existing).

- [ ] **Step 7: Commit**

```bash
git add src/r_bridge.cpp
git commit -m "fix(r_bridge): add sinkhorn to p.algorithm switch; fix ieppa_soft dispatch indentation; PROTECT child SEXPs in cell_table_probe"
```

---

## Task 7: types.hpp Documentation — S5

**Files:**
- Modify: `src/types.hpp` (line 68)

One-line comment addition. No recompile needed (comment-only change in a header only included by compiled units that are already built).

- [ ] **Step 1: Document `burnin` as iEPPA-only (line 68 in CalibSorCfg)**

Find:

```cpp
    int    burnin        = 20;
```

Replace:

```cpp
    int    burnin        = 20;  // iterations before SOR adaptation starts; iEPPA only (raking ignores)
```

- [ ] **Step 2: Verify file compiles (verify no stale .o from header-only change)**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 3: Run full test suite — final gate**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: FAIL ≤ 3 (pre-existing).

- [ ] **Step 4: Commit**

```bash
git add src/types.hpp
git commit -m "docs(types): note CalibSorCfg.burnin is iEPPA-only; raking ignores it"
```

---

## Self-Review Against Issue Index

| ID | Task | Covered? |
|----|------|----------|
| C1 | Task 1 | ✓ two exits, both fixed to RK_OK |
| C2 | Task 2 | ✓ compute_cell_metrics + wchange→l1 |
| C3 | Task 3 | ✓ docstring rewritten |
| R1 | Task 5 | ✓ @param method updated |
| R2 | Task 4 | ✓ warning text corrected |
| R3 | Task 4 | ✓ "stall_errRp"→"stall_wchange" + roxygen |
| R4 | Task 6 | ✓ sinkhorn added to p.algorithm switch |
| R5 | Task 5 | ✓ @param design_weights added |
| R6 | Task 5 | ✓ algorithm attr on attach_weights=FALSE |
| S1 | Task 5 | ✓ metric guard on PCT stall detection |
| S2 | Task 6 | ✓ comment + indentation fixed |
| S3 | Task 6 | ✓ PROTECT added, UNPROTECT(4) |
| S4 | Task 3 | ✓ errRp_w1 comment corrected |
| S5 | Task 7 | ✓ burnin comment added |
| S6 | Task 3 | ✓ split into two named constants |

**No placeholders**: all code blocks are complete. ✓  
**Type consistency**: `lbw::compute_cell_metrics` signature matches calib_dispatch.hpp. `UNPROTECT(4)` matches 4 PROTECTs in S3. ✓  
**Compile gates**: Tasks 1, 2, 3 (raking.cpp), Task 6 (r_bridge.cpp), Task 7 (types.hpp all have explicit compile+test steps. ✓
