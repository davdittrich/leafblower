# Performance + Correctness Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix 15 issues identified in the 2026-04-28 critical code review: 3 blocking correctness bugs in ieppa.cpp, 5 performance/required fixes split between ieppa.cpp and raking.cpp, 1 harvest.R API gap, and 6 simplification cleanups.

**Architecture:** Changes are ordered by compile dependency. C++ tasks get their own compile gate before advancing. ieppa.cpp issues (Tasks 1–3, 6) are sequenced first; raking.cpp performance (Tasks 4–5) second; R-only (Task 7) and cleanup (Task 8) last. No API changes.

**Tech Stack:** C++17 (`src/ieppa.cpp`, `src/raking.cpp`, `src/sinkhorn.cpp`, `src/calib_dispatch.hpp`), R (`R/harvest.R`), testthat.

---

## Issue Index

| ID | File | Lines | Description |
|----|------|-------|-------------|
| C1 | `src/ieppa.cpp` | 1192–1205 | PCT stall warning fires on RK_OK and for non-max_err metrics |
| C2 | `src/ieppa.cpp` | 1186 | RK_ERR_BUDGET=4 / RK_ERR_STALL=5 defined but never emitted |
| C3 | `src/ieppa.cpp` | 688–718 | Overflow fallback doesn't reset best-iterate tracking |
| R1 | `src/raking.cpp` | 381, 410 | Wasted `X_star = X_snap` copies immediately overwritten |
| R2 | `src/raking.cpp` | 336–392 | 5 × O(M_cell) heap allocs inside SQUAREM while loop |
| R3 | `src/raking.cpp` | 317, 451 | `compute_cell_metrics` duplicates O(K×M_cell) already done in F_eval |
| R4 | `src/ieppa.cpp` | 1211 | `convergence_tol` stores `pct_tol` even when `absolute_tol` fired |
| R5 | `src/ieppa.cpp` | 1001–1007 | `L1_WEIGHT` missing from `need_extra_metrics` gate |
| R6 | `R/harvest.R` | ~83 | `design_weights` silently ignored when `start_weights` also set |
| S1 | `src/raking.cpp` | 353, 441 | `1/n_per_cell[c]` divisions in hot loops; precompute |
| S2 | `src/ieppa.cpp` | 944–977 | Duplicate O(K×M_cell) S_lin pass in SOR adaptation block |
| S3 | `src/ieppa.cpp` | 1120–1121 | Magic number `2.302585` for ln(10) |
| S4 | `src/sinkhorn.cpp` | 204 | `!W_best.empty()` always true — dead condition |
| S5 | `src/calib_dispatch.hpp` | 99 | `apply_rule` side-effect (prev always updated) undocumented |
| S6 | `src/raking.cpp` | 405 | `(alpha + (-1.0)) / 2.0` — write as `(alpha - 1.0) / 2.0` |

---

## File Map

| File | Tasks |
|------|-------|
| `tests/testthat/test-calibration-solvers.R` | Task 0 (RED tests + regression guards) |
| `src/ieppa.cpp` | Tasks 1 (C1), 2 (C2), 3 (C3+regression), 6 (R4+R5+S3); S2 dropped |
| `src/raking.cpp` | Tasks 4 (R1+R2+S6), 5 (R3+S1) |
| `R/harvest.R` | Task 7 (R6) |
| `src/sinkhorn.cpp` | Task 8 (S4) |
| `src/calib_dispatch.hpp` | Task 8 (S5) |

---

## Task 0: Write RED Tests (before any code)

**Files:**
- Modify: `tests/testthat/test-calibration-solvers.R` (append)

- [ ] **Step 1: Append RED tests**

```r
# ── Performance + Correctness fix tests ─────────────────────────────────────

test_that("ieppa-c1-red: PCT stall Rprintf must not appear for metric=kl (not max_err/mean_err)", {
  # C1: The stall message goes via Rprintf (not R warning). capture.output(type="output") catches it.
  # C1 RED: pre-fix guard checks `metric != L1_WEIGHT`. metric=kl passes → message fires.
  # C1 GREEN: post-fix guard checks `metric in {max_err, mean_err}`. metric=kl blocked.
  # Reliability: use max_iterations=3L to GUARANTEE budget exhaustion AND max_error >> 0.01.
  # With max_iterations=3: max_error is still large (>> 10 * pct_tol = 0.01), so the
  # pre-fix guard fires deterministically.
  set.seed(42L)
  df <- data.frame(
    v1 = factor(c(rep("A", 80L), rep("B", 20L))),
    v2 = factor(c(rep("X", 40L), rep("Y", 60L)))
  )
  tgt <- list(v1 = c("A" = 0.5, "B" = 0.5), v2 = c("X" = 0.5, "Y" = 0.5))

  out <- capture.output(
    suppressWarnings(leafblower::harvest(
      df, tgt, method = "ieppa", verbose = 1,
      convergence = list(metric = "kl", rule = "improvement", tol = 0.001),
      max_iterations = 3L,  # force budget exhaustion → max_error guaranteed >> 0.01
      attach_weights = FALSE)),
    type = "output"
  )

  pct_stall_in_output <- any(grepl("PCT convergence stall", out, fixed = TRUE))

  # RED: TRUE (pre-fix guard: metric=kl != L1_WEIGHT, max_error >> 0.01 → message fires)
  # GREEN: FALSE (post-fix guard: metric=kl not in {max_err, mean_err} → message blocked)
  expect_false(pct_stall_in_output,
               label = "C1: PCT stall Rprintf must not appear for metric=kl (not max_err/mean_err)")
})

test_that("ieppa-c2-red: non-convergent iEPPA must return status 4 or 5 not 1", {
  # C2 RED: ieppa.cpp exits with RK_ERR_NOCONV=1 on budget exhaustion.
  # harvest.R maps status=4->"budget" and status=5->"stall_kl" but never receives them.
  # C2 GREEN: status=4 (budget exhausted) or status=5 (stall) returned.
  set.seed(99L)
  df <- data.frame(
    v1 = factor(c(rep("A", 80L), rep("B", 20L))),
    v2 = factor(c(rep("X", 60L), rep("Y", 40L)))
  )
  tgt <- list(v1 = c("A" = 0.2, "B" = 0.8), v2 = c("X" = 0.3, "Y" = 0.7))
  w <- suppressWarnings(
    leafblower::harvest(df, tgt, method = "ieppa", max_iterations = 2L,
                        attach_weights = FALSE))
  r <- attr(w, "result")
  # RED: r$status == 1L (RK_ERR_NOCONV)
  # GREEN: r$status %in% c(4L, 5L)
  expect_true(r$status %in% c(4L, 5L),
              label = "C2: iEPPA non-convergence must return budget(4) or stall(5), not legacy(1)")
})
```

- [ ] **Step 2: Verify RED state**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-calibration-solvers.R')" 2>&1 | \
  grep -E "ieppa-c[12]-red|PASS|FAIL" | tail -5
```

Expected: `ieppa-c1-red` FAIL (warning emitted), `ieppa-c2-red` FAIL (status==1).

- [ ] **Step 3: Run syntax gate**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```

Expected: FAIL ≤ 5 (3 pre-existing + 2 new RED). Confirms test file valid R.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-calibration-solvers.R
git commit -m "test(ieppa): RED tests for C1 (PCT stall false positive) and C2 (NOCONV not classified)"
```

---

## Task 1: Fix C1 — PCT Stall Warning Guard

**Files:**
- Modify: `src/ieppa.cpp` (lines 1192–1205)

**Context:** The warning fires when `cfg.metric != L1_WEIGHT && max_error >> pct_tol`. Bug: (a) fires even when `res.status == RK_OK` (successful convergence with marginal_kl metric); (b) excludes only L1_WEIGHT but not KL/CHI2/GRAKE_NORM/MARGINAL_KL. Fix: require `status != RK_OK` AND restrict to max_err/mean_err only.

- [ ] **Step 1: Fix the guard condition**

Find (lines 1192–1197):
```cpp
    {
        const auto& cfg = st.convergence_cfg;
        if (cfg.pct_tol > 0.0 &&
            cfg.metric != CalibMetric::L1_WEIGHT &&
            res.max_error > 10.0 * cfg.pct_tol &&
            st.log_fn != nullptr) {
```

Replace:
```cpp
    {
        const auto& cfg = st.convergence_cfg;
        if (res.status != RK_OK &&
            (cfg.metric == CalibMetric::MAX_ERR || cfg.metric == CalibMetric::MEAN_ERR) &&
            cfg.pct_tol > 0.0 &&
            res.max_error > 10.0 * cfg.pct_tol &&
            st.log_fn != nullptr) {
```

- [ ] **Step 2: Build**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 3: Run C1 RED test — must now be GREEN**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-calibration-solvers.R')" 2>&1 | \
  grep -E "ieppa-c1|PASS|FAIL" | tail -5
```

Expected: `ieppa-c1-red` PASS.

- [ ] **Step 4: Full suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```

Expected: FAIL ≤ 4 (3 pre-existing + ieppa-c2-red still RED).

- [ ] **Step 5: Commit**

```bash
git add src/ieppa.cpp
git commit -m "fix(ieppa): PCT stall warning requires status!=RK_OK and metric in {max_err,mean_err} — prevents false positives on successful marginal_kl convergence"
```

---

## Task 2: Fix C2 — ieppa.cpp Emits RK_ERR_BUDGET / RK_ERR_STALL

**Files:**
- Modify: `src/ieppa.cpp` (after line 1186)

**Context:** After the homotopy loop, `res.status == RK_ERR_NOCONV` when budget exhausted. `RK_ERR_BUDGET=4` (still descending) and `RK_ERR_STALL=5` (metric never improved) are defined but never emitted. Heuristic: if `best_metric_seen` is finite (metric improved at some point during the solve), it's BUDGET; if still `∞` (metric never decreased from its initial value), it's STALL.

- [ ] **Step 1: Classify NOCONV into BUDGET or STALL**

Find the line immediately after `}  // end homotopy level loop` at line ~1186. Insert after it:

```cpp
    // Classify RK_ERR_NOCONV into BUDGET or STALL.
    // RK_ERR_BUDGET (4): metric was improving; increase max_iterations.
    // RK_ERR_STALL  (5): metric never improved from initial → constrained optimum.
    if (res.status == RK_ERR_NOCONV) {
        res.status = std::isfinite(best_metric_seen) ? RK_ERR_BUDGET : RK_ERR_STALL;
    }
```

- [ ] **Step 2: Build**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 3: Run C2 RED test — must now be GREEN**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-calibration-solvers.R')" 2>&1 | \
  grep -E "ieppa-c2|PASS|FAIL" | tail -5
```

Expected: `ieppa-c2-red` PASS (status is 4 or 5, not 1).

- [ ] **Step 4: Verify convergence_reason in harvest.R**

```bash
Rscript -e "
  set.seed(99L)
  df <- data.frame(v1=factor(c(rep('A',80L),rep('B',20L))),
                   v2=factor(c(rep('X',60L),rep('Y',40L))))
  tgt <- list(v1=c('A'=0.2,'B'=0.8),v2=c('X'=0.3,'Y'=0.7))
  w <- suppressWarnings(leafblower::harvest(df,tgt,method='ieppa',max_iterations=2L,
         attach_weights=FALSE))
  r <- attr(w,'result')
  cat('status:', r\$status, '\n')
  cat('reason:', r\$convergence_used\$convergence_reason, '\n')
" 2>&1 | grep -E "^status:|^reason:"
```

Expected: `status: 4` and `reason: budget` (or `status: 5` and `reason: stall_kl`).

- [ ] **Step 5: Full suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```

Expected: FAIL ≤ 3 (pre-existing only).

- [ ] **Step 6: Commit**

```bash
git add src/ieppa.cpp
git commit -m "fix(ieppa): classify RK_ERR_NOCONV into BUDGET(4)/STALL(5) — finite best_metric_seen=BUDGET, never-improved=STALL"
```

---

## Task 3: Fix C3 — Reset Best-Iterate on Overflow Fallback

**Files:**
- Modify: `src/ieppa.cpp` (lines 688–718)
- Modify: `tests/testthat/test-calibration-solvers.R` (append regression guard)

**Context:** When linear-space overflow triggers the fallback to log-space (lines 688–718), the solver resets all state but leaves `best_metric_seen`, `W_best`, `best_iter_val` from the pre-fallback degenerate state. If any callers use `best_weights`, they receive corrupted weights silently.

- [ ] **Step 1: Reset best-iterate tracking in the fallback block**

Find (lines ~710–717, inside `if (overflow_trip && !linear_fallback_used)`):
```cpp
                // WU-B Fix 2: reset X_prev after fallback — X semantics changed (log-path).
                for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];
                if (st.verbose >= 1) {
                    st.log("iEPPA: linear-space overflow trip; fallback to log-space.");
                }
                continue;
```

Insert after the X_prev reset line and before the verbose log:
```cpp
                for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];
                // Reset best-iterate: pre-fallback W_best came from numerically
                // degenerate linear-space iterates; discard rather than return as "best".
                best_metric_seen = std::numeric_limits<double>::infinity();
                W_best.assign(ct.M_cell, 0.0);
                best_iter_val    = 0;
                if (st.verbose >= 1) {
                    st.log("iEPPA: linear-space overflow trip; fallback to log-space.");
                }
                continue;
```

- [ ] **Step 2: Add regression guard test**

Append to `tests/testthat/test-calibration-solvers.R`:
```r
test_that("ieppa-c3: best_weights has finite max after overflow fallback path", {
  # C3 regression guard: if overflow fallback doesn't reset best-iterate tracking,
  # best_weights can contain Inf/NaN from degenerate linear-space state.
  # This test verifies best_weights is all-finite regardless of fallback path.
  # (Cannot deterministically trigger overflow in unit test; verifies invariant holds.)
  set.seed(42L)
  df  <- data.frame(v1 = factor(c(rep("A", 100L), rep("B", 100L))))
  tgt <- list(v1 = c("A" = 0.5, "B" = 0.5))
  w <- suppressWarnings(
    leafblower::harvest(df, tgt, method = "ieppa", max_iterations = 200L,
                        attach_weights = FALSE))
  r <- attr(w, "result")
  if (is.finite(r$best_error)) {
    expect_true(all(is.finite(r$best_weights)),
                label = "C3: best_weights must be all-finite (no Inf/NaN from overflow)")
    expect_true(all(r$best_weights >= 0),
                label = "C3: best_weights must be non-negative")
  }
})
```

- [ ] **Step 3: Build**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 4: Run regression guard**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-calibration-solvers.R')" 2>&1 | \
  grep -E "ieppa-c3|PASS|FAIL" | tail -5
```

Expected: `ieppa-c3` PASS.

- [ ] **Step 5: Full suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```

Expected: FAIL ≤ 3 (pre-existing).

- [ ] **Step 6: Commit**

```bash
git add src/ieppa.cpp tests/testthat/test-calibration-solvers.R
git commit -m "fix(ieppa): reset best_metric_seen/W_best/best_iter_val on overflow fallback — pre-fallback snapshot from degenerate linear-space was returned as best_weights"
```

---

## Task 4: R1+R2+S6 — raking.cpp SQUAREM Buffer Pre-allocation + Wasted Copy Removal

**Files:**
- Modify: `src/raking.cpp` (SQUAREM while loop, lines ~327–480)

**Context:**
- R1: `auto X_star = X_snap` (line 381) and `X_star = X_snap` (line 410 in halving loop) are O(M_cell) copies immediately overwritten by the fill loop. For stepstone M_cell≈2M, worst case with kMaxHalvings=16: 17 × 16MB = 272MB wasted bandwidth per super-step.
- R2: Five O(M_cell) `auto` declarations inside the while loop trigger malloc/free per super-step. Pre-allocate before the loop.
- S6: `(alpha + (-1.0)) / 2.0` → `(alpha - 1.0) / 2.0`.

- [ ] **Step 1: Pre-allocate SQUAREM scratch buffers before the while loop**

Find (after `auto X_prev_sq = X;` before the while loop, around line 329):
```cpp
            int f_eval_count = 0;
            auto X_prev_sq = X;  // snapshot for weight-change stall; updated each accepted super-step
            while (f_eval_count + 3 <= st.inner_max_iter) {
```

Replace:
```cpp
            int f_eval_count = 0;
            auto X_prev_sq = X;  // snapshot for weight-change stall; updated each accepted super-step
            // Pre-allocate SQUAREM scratch — avoids 5 × O(M_cell) malloc/free per super-step.
            std::vector<double> sq_w1(ct.M_cell), sq_w2(ct.M_cell);
            std::vector<double> sq_X_snap(ct.M_cell), sq_X_star(ct.M_cell), sq_X_star_pre(ct.M_cell);
            while (f_eval_count + 3 <= st.inner_max_iter) {
```

- [ ] **Step 2: Replace `auto w1 = X` and `auto w2 = w1` with pre-allocated buffers**

Find (lines 336–343 inside the while loop):
```cpp
                auto w1 = X;
                double errRp_w1 = F_eval(w1);  ++f_eval_count;
                (void)errRp_w1;  // errRp_k updated inside F_eval but not consumed (use_greedy=false when accelerate=true)
                is_infeasible = infeas_before;

                auto w2 = w1;
                double errRp_w2 = F_eval(w2);  ++f_eval_count;
                is_infeasible = infeas_before;
```

Replace:
```cpp
                sq_w1 = X;
                double errRp_w1 = F_eval(sq_w1);  ++f_eval_count;
                (void)errRp_w1;  // errRp_k updated inside F_eval but not consumed (use_greedy=false when accelerate=true)
                is_infeasible = infeas_before;

                sq_w2 = sq_w1;
                double errRp_w2 = F_eval(sq_w2);  ++f_eval_count;
                is_infeasible = infeas_before;
```

- [ ] **Step 3: Update full obs-level norm loop (all w1/w2 references including w2_sq_obs)**

Find the COMPLETE obs-level norm loop (lines ~352–359). Must match the full block including the `w2_sq_obs` accumulator line:
```cpp
                for (int c = 0; c < ct.M_cell; c++) {
                    const double inv_nc = 1.0 / static_cast<double>(ct.n_per_cell[c]);
                    const double ri = w1[c] - X[c],  vi = w2[c] - w1[c];
                    r_sq_obs  += ri * ri * inv_nc;
                    v_sq_obs  += vi * vi * inv_nc;
                    w2_sq_obs += w2[c] * w2[c] * inv_nc;
                    v_sq_cell += vi * vi;
                }
```

Replace with (rename w1→sq_w1, w2→sq_w2; keep `1.0 / n_per_cell[c]` — Task 5 Step 3 upgrades it later):
```cpp
                for (int c = 0; c < ct.M_cell; c++) {
                    const double inv_nc = 1.0 / static_cast<double>(ct.n_per_cell[c]);
                    const double ri = sq_w1[c] - X[c],  vi = sq_w2[c] - sq_w1[c];
                    r_sq_obs  += ri * ri * inv_nc;
                    v_sq_obs  += vi * vi * inv_nc;
                    w2_sq_obs += sq_w2[c] * sq_w2[c] * inv_nc;
                    v_sq_cell += vi * vi;
                }
```

**Note:** `norm_w2` is `std::sqrt(w2_sq_obs)` declared after the loop — it does not need renaming since it's derived from the (now correctly updated) `w2_sq_obs` accumulator.

- [ ] **Step 4: Fix fixed-point guard to use sq_w2**

Find (line ~364–370):
```cpp
                if (norm_v / (norm_w2 + kVNormEps) < kVNormRel) {
                    X = w2;
                    res.max_error        = errRp_w2;
```

Replace:
```cpp
                if (norm_v / (norm_w2 + kVNormEps) < kVNormRel) {
                    X = sq_w2;
                    res.max_error        = errRp_w2;
```

- [ ] **Step 5: Replace `auto X_snap = w2` with pre-allocated buffer (no wasted copy)**

Find (line ~377–386):
```cpp
                // Snapshot: state at w2, before extrapolation
                auto X_snap = w2;

                // Extrapolate X* = X_snap - 2α·r + α²·v; clamp to ≥ 0
                auto X_star = X_snap;
                for (int c = 0; c < ct.M_cell; c++) {
                    double ri = w1[c] - X[c],  vi = w2[c] - w1[c];
                    X_star[c] = X_snap[c] - 2.0 * alpha * ri + alpha * alpha * vi;
                    if (X_star[c] < 0.0) X_star[c] = 0.0;
                }
```

Replace (no copy for X_snap — just point sq_X_snap at sq_w2 content; no wasted copy for X_star — fill loop overwrites all):
```cpp
                // Snapshot: state at w2, before extrapolation
                sq_X_snap = sq_w2;

                // Extrapolate X* = X_snap - 2α·r + α²·v; clamp to ≥ 0
                // No initial copy — the loop below fills every sq_X_star[c] unconditionally.
                for (int c = 0; c < ct.M_cell; c++) {
                    double ri = sq_w1[c] - X[c],  vi = sq_w2[c] - sq_w1[c];
                    sq_X_star[c] = sq_X_snap[c] - 2.0 * alpha * ri + alpha * alpha * vi;
                    if (sq_X_star[c] < 0.0) sq_X_star[c] = 0.0;
                }
```

- [ ] **Step 6: Replace X_star_pre and first F_eval(X_star)**

Find (lines ~392–398):
```cpp
                const double plain_resid = v_sq_cell;  // cell-level ‖v‖² — matches cand_resid geometry
                auto X_star_pre = X_star;
                double errRp_new = F_eval(X_star);  ++f_eval_count;
                double cand_resid = 0.0;
                for (int c = 0; c < ct.M_cell; c++) {
                    double d = X_star[c] - X_star_pre[c];
                    cand_resid += d * d;
                }
```

Replace:
```cpp
                const double plain_resid = v_sq_cell;  // cell-level ‖v‖² — matches cand_resid geometry
                sq_X_star_pre = sq_X_star;
                double errRp_new = F_eval(sq_X_star);  ++f_eval_count;
                double cand_resid = 0.0;
                for (int c = 0; c < ct.M_cell; c++) {
                    double d = sq_X_star[c] - sq_X_star_pre[c];
                    cand_resid += d * d;
                }
```

- [ ] **Step 7: Fix halving loop — remove wasted copy, update all references**

Find the halving loop (lines ~402–422):
```cpp
                bool fell_back = false;
                for (int h = 0; h < kMaxHalvings && cand_resid > kHalvingSlack * plain_resid; h++) {
                    is_infeasible = infeas_before;  // discard probe's infeasibility
                    alpha = (alpha + (-1.0)) / 2.0;  // midpoint toward -1 (autumn formula)
                    if (std::fabs(alpha - (-1.0)) < 1e-3) {
                        // Fell back to plain step
                        X = w2; errRp_new = errRp_w2; fell_back = true; break;
                    }
                    X_star = X_snap;
                    for (int c = 0; c < ct.M_cell; c++) {
                        double ri = w1[c] - X[c],  vi = w2[c] - w1[c];
                        X_star[c] = X_snap[c] - 2.0 * alpha * ri + alpha * alpha * vi;
                        if (X_star[c] < 0.0) X_star[c] = 0.0;
                    }
                    X_star_pre = X_star;
                    errRp_new = F_eval(X_star);  ++f_eval_count;
                    cand_resid = 0.0;
                    for (int c = 0; c < ct.M_cell; c++) {
                        double d = X_star[c] - X_star_pre[c];
                        cand_resid += d * d;
                    }
                }
                if (!fell_back) X = X_star;
```

Replace (S6 also: `(alpha - 1.0) / 2.0`; remove wasted `X_star = X_snap` copy):
```cpp
                bool fell_back = false;
                for (int h = 0; h < kMaxHalvings && cand_resid > kHalvingSlack * plain_resid; h++) {
                    is_infeasible = infeas_before;  // discard probe's infeasibility
                    alpha = (alpha - 1.0) / 2.0;  // midpoint toward -1 (autumn formula)
                    if (std::fabs(alpha - (-1.0)) < 1e-3) {
                        // Fell back to plain step
                        X = sq_w2; errRp_new = errRp_w2; fell_back = true; break;
                    }
                    // No copy of sq_X_snap needed — loop fills every sq_X_star[c] unconditionally.
                    for (int c = 0; c < ct.M_cell; c++) {
                        double ri = sq_w1[c] - X[c],  vi = sq_w2[c] - sq_w1[c];
                        sq_X_star[c] = sq_X_snap[c] - 2.0 * alpha * ri + alpha * alpha * vi;
                        if (sq_X_star[c] < 0.0) sq_X_star[c] = 0.0;
                    }
                    sq_X_star_pre = sq_X_star;
                    errRp_new = F_eval(sq_X_star);  ++f_eval_count;
                    cand_resid = 0.0;
                    for (int c = 0; c < ct.M_cell; c++) {
                        double d = sq_X_star[c] - sq_X_star_pre[c];
                        cand_resid += d * d;
                    }
                }
                if (!fell_back) X = sq_X_star;
```

- [ ] **Step 8: Update best-iterate tracking reference from X_star to sq_X_star**

Find (line ~427–432 — best-iterate block):
```cpp
                // Best-iterate tracking
                if (errRp_new < best_metric_seen) {
                    best_metric_seen    = errRp_new;
                    best_iter_val       = f_eval_count;
                    best_objective_seen = compute_weight_kl();
                    W_best              = X;
```

Verify `W_best = X` is correct here (X was just set to `sq_X_star` or `sq_w2` above). No change needed — `X` holds the accepted iterate.

- [ ] **Step 9: Build**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`.

If compile error referencing undeclared `w1`/`w2`/`X_star`/`X_snap`/`X_star_pre`:
```bash
grep -n "w1\[c\]\|w2\[c\]\|X_star\[c\]\|X_snap\[c\]\|X_star_pre" src/raking.cpp | head -20
```
Expected: **0 matches** (all usages renamed to `sq_*` prefix by steps 2–8). Any match = a missed rename: apply the `sq_` prefix and rebuild.

- [ ] **Step 10: Full suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```

Expected: FAIL ≤ 3 (pre-existing). If squarem tests fail: compare behavior vs pre-fix baseline by running the stepstone benchmark if available.

- [ ] **Step 11: Commit**

```bash
git add src/raking.cpp
git commit -m "perf(squarem): pre-allocate 5 scratch buffers; remove 17 wasted O(M_cell) copies per super-step — up to 272MB/step bandwidth saved on stepstone-scale problems"
```

---

## Task 5: R3+S1 — F_eval Populates CellMetrics; Precompute inv_n_per_cell

**Files:**
- Modify: `src/raking.cpp` (F_eval lambda, SQUAREM convergence block, pre-loop section)

**Context:**
- R3: `F_eval` ends with `compute_errRp_ct(st, ct, Xv, bucket)` — a full O(K×M_cell) pass. After accepting a SQUAREM step, Task 2's fix calls `compute_cell_metrics(...)` — another O(K×M_cell) pass. Replace `compute_errRp_ct` inside F_eval with `compute_cell_metrics` storing into `last_F_metrics`, then use `last_F_metrics` in the convergence block (one pass, not two).
- S1: `1.0 / n_per_cell[c]` computed per cell per super-step in two hot loops. Precompute as `inv_n_per_cell` vector once before the while loop.

- [ ] **Step 1: Declare `last_F_metrics` and `inv_n_per_cell` before the `if (st.accelerate)` block**

Find the line `bool is_infeasible = false;` (around line 116). After it, add:

```cpp
    // Pre-computed reciprocals for obs-level norm and wchange computations.
    // n_per_cell is constant across iterations; avoids one division per cell per super-step.
    std::vector<double> inv_n_per_cell(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++)
        inv_n_per_cell[c] = 1.0 / static_cast<double>(ct.n_per_cell[c]);
```

And inside `if (st.accelerate) { if (st.inner_max_iter >= 3) {` (before the while loop), add:

```cpp
            lbw::CellMetrics last_F_metrics;  // populated by last F_eval call; used for convergence check
```

- [ ] **Step 2: Replace `compute_errRp_ct` at the end of F_eval with `compute_cell_metrics`**

Find (last line of the F_eval lambda, around line 317):
```cpp
        return compute_errRp_ct(st, ct, Xv, bucket);
    };
```

Replace:
```cpp
        // compute_cell_metrics is a strict superset of compute_errRp_ct; uses W=n (post-hyperplane).
        // Stores all metrics in last_F_metrics — SQUAREM convergence check reuses them directly,
        // eliminating a second O(K×M_cell) pass after each accepted step.
        last_F_metrics = lbw::compute_cell_metrics(st, ct, Xv, static_cast<double>(st.n), bucket);
        return last_F_metrics.errRp;
    };
```

- [ ] **Step 3: Replace obs-level norm loop divisions with `inv_n_per_cell`**

Find (in the SQUAREM while loop, the norm computation block):
```cpp
                for (int c = 0; c < ct.M_cell; c++) {
                    const double inv_nc = 1.0 / static_cast<double>(ct.n_per_cell[c]);
                    const double ri = sq_w1[c] - X[c],  vi = sq_w2[c] - sq_w1[c];
                    r_sq_obs  += ri * ri * inv_nc;
                    v_sq_obs  += vi * vi * inv_nc;
                    w2_sq_obs += sq_w2[c] * sq_w2[c] * inv_nc;
                    v_sq_cell += vi * vi;
                }
```

Replace:
```cpp
                for (int c = 0; c < ct.M_cell; c++) {
                    const double inv_nc = inv_n_per_cell[c];
                    const double ri = sq_w1[c] - X[c],  vi = sq_w2[c] - sq_w1[c];
                    r_sq_obs  += ri * ri * inv_nc;
                    v_sq_obs  += vi * vi * inv_nc;
                    w2_sq_obs += sq_w2[c] * sq_w2[c] * inv_nc;
                    v_sq_cell += vi * vi;
                }
```

- [ ] **Step 4: Replace wchange loop division with `inv_n_per_cell`**

Find (the wchange computation block):
```cpp
                double wchange = 0.0;
                for (int c = 0; c < ct.M_cell; c++)
                    wchange += std::fabs(X[c] - X_prev_sq[c]) / static_cast<double>(ct.n_per_cell[c]);
                wchange /= static_cast<double>(st.n);
```

Replace:
```cpp
                double wchange = 0.0;
                for (int c = 0; c < ct.M_cell; c++)
                    wchange += std::fabs(X[c] - X_prev_sq[c]) * inv_n_per_cell[c];
                wchange /= static_cast<double>(st.n);
```

- [ ] **Step 5: Replace separate `compute_cell_metrics` call in SQUAREM convergence block with `last_F_metrics`**

Find (the convergence criterion block inside SQUAREM):
```cpp
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
```

Replace:
```cpp
                // Convergence criterion: reuse last_F_metrics from the accepted F_eval call.
                // F_eval now calls compute_cell_metrics internally (Task 5), so last_F_metrics
                // holds all metrics from the accepted step — no second O(K×M_cell) pass needed.
                {
                    lbw::CellMetrics m_conv = last_F_metrics;  // populated by accepted F_eval
                    m_conv.errRp = errRp_new;  // F_eval's errRp is authoritative (post-hyperplane)
                    m_conv.l1    = wchange;
                    if (lbw::check_convergence(st.convergence_cfg, m_conv,
                                               prev_metric_for_rule, st.tol_abs)) {
```

- [ ] **Step 6: Build**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`. If `last_F_metrics` is out of scope for F_eval lambda: move its declaration to before the F_eval lambda definition (before line ~165).

- [ ] **Step 7: Full suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```

Expected: FAIL ≤ 3 (pre-existing). squarem-c2 must remain GREEN.

- [ ] **Step 8: Commit**

```bash
git add src/raking.cpp
git commit -m "perf(squarem): F_eval populates last_F_metrics via compute_cell_metrics — eliminates second O(K×M_cell) pass per super-step; precompute inv_n_per_cell"
```

---

## Task 6: R4+R5+S3 — ieppa.cpp Diagnostic Correctness + ln10 Constant

**Files:**
- Modify: `src/ieppa.cpp` (lines 1001–1007, line 356 area, 1211, 1120–1121)

**Note on S2 (SOR S_lin dedup — dropped from scope):** The SOR adaptation block executes at line ~945, BEFORE the `need_extra_metrics` block at line ~1001. Passing values from `need_extra_metrics` to SOR would introduce a one-iteration lag. Fixing this requires block reordering which carries regression risk; tracked as a separate follow-up.

Three changes: R5 (L1_WEIGHT gate), R4 (convergence_tol flag), S3 (magic number). No behavior change on primary convergence path.

- [ ] **Step 1: R5 — Add `L1_WEIGHT` to `need_extra_metrics` gate (lines 1001–1007)**

Find:
```cpp
            const bool need_extra_metrics =
                (metric == lbw::CalibMetric::MEAN_ERR    ||
                 metric == lbw::CalibMetric::KL          ||
                 metric == lbw::CalibMetric::CHI2        ||
                 metric == lbw::CalibMetric::GRAKE_NORM  ||
                 metric == lbw::CalibMetric::MARGINAL_KL ||
                 iter_in_lvl == budget_lvl);
```

Replace:
```cpp
            const bool need_extra_metrics =
                (metric == lbw::CalibMetric::MEAN_ERR    ||
                 metric == lbw::CalibMetric::KL          ||
                 metric == lbw::CalibMetric::CHI2        ||
                 metric == lbw::CalibMetric::GRAKE_NORM  ||
                 metric == lbw::CalibMetric::MARGINAL_KL ||
                 metric == lbw::CalibMetric::L1_WEIGHT   ||
                 iter_in_lvl == budget_lvl);
```

- [ ] **Step 2: R4 — Fix `convergence_tol` to reflect which threshold fired (line 1211)**

Find (immediately after the `// WU-C: populate convergence diagnostics` comment, around line 1208–1212):
```cpp
    res.convergence_tol                = st.convergence_cfg.pct_tol;
```

Three sub-steps for R4:

**Step 2a** — Add `bool absolute_tol_fired = false;` before the outer homotopy `for` loop. Find the loop start:
```cpp
    for (int lvl = 0; lvl < N_levels && !homotopy_break; lvl++) {
```
(This is line 356 in ieppa.cpp — confirmed by grep.) Add immediately before it:
```cpp
    bool absolute_tol_fired = false;  // set when absolute_tol triggers convergence
```

**Step 2b** — Inside the convergence block, track when `converged_abs` fires. Find:
```cpp
                bool converged_abs = (cfg.absolute_tol > 0.0) && (curr_metric < cfg.absolute_tol);
```
Replace:
```cpp
                bool converged_abs = (cfg.absolute_tol > 0.0) && (curr_metric < cfg.absolute_tol);
                if (converged && converged_abs) absolute_tol_fired = true;
```

**Step 2c** — At line 1211, update the assignment:
```cpp
    res.convergence_tol = absolute_tol_fired
        ? st.convergence_cfg.absolute_tol : st.convergence_cfg.pct_tol;
```

- [ ] **Step 3: S3 — Replace magic number 2.302585 with named constant (lines 1120–1121)**

Find (inside the verbose=2 logging block):
```cpp
                    // log10(exp(v)) = v / ln(10)
                    std::snprintf(msg, sizeof(msg),
                                  "  margin=%d: log10(f) range [%.2f, %.2f]",
                                  k + 1,
                                  lf_min / 2.302585,
                                  lf_max / 2.302585);
```

Replace:
```cpp
                    static constexpr double kLn10 = 2.302585092994046;  // std::log(10.0)
                    std::snprintf(msg, sizeof(msg),
                                  "  margin=%d: log10(f) range [%.2f, %.2f]",
                                  k + 1,
                                  lf_min / kLn10,
                                  lf_max / kLn10);
```

- [ ] **Step 5: Build**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`.

- [ ] **Step 6: Full suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```

Expected: FAIL ≤ 3 (pre-existing).

- [ ] **Step 7: Commit**

```bash
git add src/ieppa.cpp
git commit -m "fix(ieppa): add L1_WEIGHT to need_extra_metrics gate; fix convergence_tol to reflect which threshold fired; name ln10 constant"
```

---

## Task 7: R6 — harvest.R design_weights Warning

**Files:**
- Modify: `R/harvest.R` (~line 83)

Pure R, no recompile.

- [ ] **Step 1: Add warning when design_weights silently ignored**

Find:
```r
  if (!is.null(design_weights) && is.null(start_weights)) {
    start_weights <- design_weights
  }
```

Replace:
```r
  if (!is.null(design_weights)) {
    if (!is.null(start_weights))
      warning("leafblower: both design_weights and start_weights supplied; design_weights ignored")
    else
      start_weights <- design_weights
  }
```

- [ ] **Step 2: Verify both paths**

```bash
Rscript -e "
  df  <- data.frame(v1=factor(c('A','B','A')))
  tgt <- list(v1=c('A'=0.5,'B'=0.5))
  # Path 1: design_weights only — should use it (no warning)
  w1 <- leafblower::harvest(df, tgt, design_weights=c(1,2,1), attach_weights=FALSE)
  cat('path1 weights:', round(as.numeric(w1), 3), '\n')
  # Path 2: both supplied — should warn
  expect_warn_str <- tryCatch({
    w2 <- leafblower::harvest(df, tgt, design_weights=c(1,2,1), start_weights=c(1,1,1),
                               attach_weights=FALSE)
    'no warning'
  }, warning=function(w) conditionMessage(w))
  cat('path2:', expect_warn_str, '\n')
" 2>&1 | grep -E "^path[12]:"
```

Expected:
- `path1 weights:` non-uniform values (design_weights used)
- `path2: leafblower: both design_weights and start_weights supplied; design_weights ignored`

- [ ] **Step 3: Full suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```

Expected: FAIL ≤ 3 (pre-existing).

- [ ] **Step 4: Commit**

```bash
git add R/harvest.R
git commit -m "fix(harvest): warn when design_weights silently ignored because start_weights also provided"
```

---

## Task 8: S4+S5 — sinkhorn.cpp Dead Branch + calib_dispatch.hpp Comment

**Files:**
- Modify: `src/sinkhorn.cpp` (line 204)
- Modify: `src/calib_dispatch.hpp` (line 99)

Pure simplification. No behavior change.

- [ ] **Step 1: S4 — Remove always-true `!W_best.empty()` condition (sinkhorn.cpp:204)**

Find (in sinkhorn.cpp):
```cpp
    if (std::isfinite(best_metric_seen) && !W_best.empty()) {
```

Replace:
```cpp
    if (std::isfinite(best_metric_seen)) {
```

`W_best` is declared as `std::vector<double> W_best(ct.M_cell, 0.0)` — always M_cell elements, never empty (M_cell ≥ 1 post-validation). The `!W_best.empty()` guard was defensive noise.

- [ ] **Step 2: S5 — Document `apply_rule` side effect (calib_dispatch.hpp:99)**

Find (at the end of `apply_rule`, before `return converged`):
```cpp
    prev = curr;   // always update so next call has a valid baseline
    return converged;
```

Replace:
```cpp
    prev = curr;   // always update — side effect: sliding-window baseline for next call
    return converged;
```

- [ ] **Step 3: Build**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 4: Full suite — final gate**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```

Expected: FAIL ≤ 3 (pre-existing).

- [ ] **Step 5: Commit**

```bash
git add src/sinkhorn.cpp src/calib_dispatch.hpp
git commit -m "simplify: remove always-true !W_best.empty() guard in sinkhorn; document apply_rule prev update side effect"
```

---

## Self-Review Against Issue Index

| ID | Task | Covered? |
|----|------|----------|
| C1 | Task 1 | ✓ `status != RK_OK` guard + restrict to `{MAX_ERR, MEAN_ERR}` |
| C2 | Task 2 | ✓ `BUDGET` when `best_metric_seen` finite; `STALL` otherwise |
| C3 | Task 3 | ✓ reset `best_metric_seen`, `W_best`, `best_iter_val` after overflow |
| R1 | Task 4 | ✓ wasted `X_star = X_snap` copies eliminated (initial + up to 16 halvings) |
| R2 | Task 4 | ✓ 5 pre-allocated buffers replace per-super-step malloc |
| R3 | Task 5 | ✓ F_eval calls `compute_cell_metrics`; SQUAREM reuses `last_F_metrics` |
| R4 | Task 6 | ✓ `absolute_tol_fired` flag; `convergence_tol = absolute_tol` when it fires |
| R5 | Task 6 | ✓ `L1_WEIGHT` added to `need_extra_metrics` |
| R6 | Task 7 | ✓ `warning()` when both `design_weights` and `start_weights` provided |
| S1 | Task 5 | ✓ `inv_n_per_cell` precomputed; divisions replaced by multiplications |
| S2 | Dropped | SOR block executes before need_extra_metrics — can't share values without one-iter lag. Tracked as follow-up. |
| S3 | Task 6 | ✓ `kLn10 = 2.302585092994046` named constant |
| S4 | Task 8 | ✓ `!W_best.empty()` removed from sinkhorn.cpp:204 |
| S5 | Task 8 | ✓ comment added to `apply_rule` prev update |
| S6 | Task 4 | ✓ `(alpha - 1.0) / 2.0` |

**No placeholders.** All code blocks complete. ✓

**Compile gates**: Tasks 1, 2, 3, 4, 5, 6, 8 each have `R CMD INSTALL --preclean` + `devtools::test()`. ✓

**Type consistency**: `lbw::CellMetrics last_F_metrics` declared before F_eval lambda, populated inside it (captured by reference), consumed in convergence block. Consistent across Tasks 5. ✓
