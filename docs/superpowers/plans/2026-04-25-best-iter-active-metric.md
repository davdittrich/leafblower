# Fix best_iter to Track Active Metric

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use checkbox syntax.

**Goal:** `best_iter` / `best_error` / `best_weights` currently always tracks the iterate with minimum `errRp` (max_err) regardless of the active convergence criterion. When `metric=kl`, best iterate should be minimum KL — not minimum errRp. Fix all three solvers.

**Architecture:** One WU. Replace `best_errRp_seen` with `best_metric_seen` using `select_metric(cfg.metric, ...)` in ieppa.cpp and raking.cpp. lbfgsb is batch (single iterate) — no change needed there.

**Tech Stack:** C++17 (ieppa.cpp, raking.cpp), R (no change), testthat.

**Ticket:** leafblower-g4oj

**Baseline:** FAIL 0 | PASS 332

---

## Key Facts

### ieppa.cpp current code (lines ~283–284, ~807–812)
```cpp
double best_errRp_seen = std::numeric_limits<double>::infinity();
int    best_iter_val   = 0;
// ...
if (errRp < best_errRp_seen) {
    best_errRp_seen = errRp;
    best_iter_val   = iter;
    for (int c ...) W_best[c] = ...;
}
```

### raking.cpp current code (lines ~131–132, ~207–209)
```cpp
double best_errRp_seen = std::numeric_limits<double>::infinity();
int    best_iter_val   = 0;
// ...
if (errRp < best_errRp_seen) {
    best_errRp_seen = errRp;
    best_iter_val   = iter;
```

### Ordering constraint
The best-iterate check at lines 807-812 in ieppa.cpp fires AFTER the metrics block (mean_err/kl/chi2/grake_norm/l1_weight are computed when need_extra_metrics=true). Need to ensure all 6 metric values are available when computing `best_metric`. Currently the check happens where `errRp`, `mean_err`, `kl_max`, `chi2_total`, `grake_norm`, `l1_weight` are all in scope.

Verify with: `grep -n "best_errRp_seen\|mean_err\|kl_max\|chi2_total\|grake_norm\|l1_weight" src/ieppa.cpp | grep "8[0-9][0-9]:"` — all should appear near line 800-812.

---

## Work Unit: Fix best-iterate in ieppa.cpp + raking.cpp

- [ ] **Step 1: Claim ticket**
```bash
bd update leafblower-g4oj --claim
```

- [ ] **Step 2: Write failing test (TDD)**

Append to `tests/testthat/test-best-iterate.R`:
```r
test_that("g4oj: best_iter tracks active metric (kl vs max_err differ)", {
  set.seed(42)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.4,"2"=0.4,"3"=0.2), b = c("1"=0.6,"2"=0.4))
  # Run with kl criterion — best_error is min KL, not min errRp
  w_kl <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                              max_iterations = 300,
                              convergence = list(metric="kl", rule="improvement", tol=0.001),
                              attach_weights = FALSE)
  r_kl <- attr(w_kl, "result")
  # Run with max_err criterion — best_error is min errRp
  w_me <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                              max_iterations = 300,
                              convergence = list(metric="max_err", rule="improvement", tol=0.001),
                              attach_weights = FALSE)
  r_me <- attr(w_me, "result")
  # KL run: best_error is a KL divergence value (positive, finite)
  expect_true(is.finite(r_kl$best_error))
  expect_gt(r_kl$best_error, 0)
  # max_err run: best_error <= max_error (standard best-iterate guarantee holds)
  expect_lte(r_me$best_error, r_me$max_error)
  # The two best_error values are different quantities (KL vs errRp) and should differ
  expect_false(isTRUE(all.equal(r_kl$best_error, r_me$best_error)),
               info = "KL best_error != max_err best_error (different metrics tracked)")
})
```

Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-best-iterate.R")' 2>&1 | tail -5`

This test will PASS before the fix (because `best_error <= max_error` is trivially true when both track errRp). The real semantic fix is verified by the benchmark comparison.

- [ ] **Step 3: Fix `src/ieppa.cpp`**

Find the `best_errRp_seen` declaration (~line 283) and rename:
```cpp
// OLD:
double best_errRp_seen = std::numeric_limits<double>::infinity();

// NEW:
double best_metric_seen = std::numeric_limits<double>::infinity();
```

Find the best-iterate update block (~lines 807–813) and change:
```cpp
// OLD:
if (errRp < best_errRp_seen) {
    best_errRp_seen = errRp;
    best_iter_val   = iter;
    for (int c = 0; c < ct.M_cell; c++) {
        W_best[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
    }
}

// NEW — track minimum of the active metric:
{
    const double curr_best_metric = lbw::select_metric(
        st.convergence_cfg.metric,
        errRp, mean_err, kl_max, chi2_total, grake_norm, l1_weight);
    if (std::isfinite(curr_best_metric) && curr_best_metric < best_metric_seen) {
        best_metric_seen = curr_best_metric;
        best_iter_val    = iter;
        for (int c = 0; c < ct.M_cell; c++) {
            W_best[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
        }
    }
}
```

Find the exit block (~lines 1103–1106) and update:
```cpp
// OLD:
res.best_error = best_errRp_seen;

// NEW:
res.best_error = best_metric_seen;
```

**IMPORTANT — exact placement (two-block approach, non-negotiable):**

`mean_err`, `kl_max`, `chi2_total`, `grake_norm`, `l1_weight` are 0.0 when gated. Using them outside the gate produces incorrect minima. The fix uses TWO separate update blocks:

**Block 1 — MAX_ERR only (stays at current line ~807, outside need_extra_metrics):**
`errRp` is always computed at kErrCheckInterval. MAX_ERR best-iterate fires at every check:
```cpp
// BLOCK 1: MAX_ERR best-iterate (errRp always valid here, outside gate).
// Remove old `if (errRp < best_errRp_seen)` and replace with:
if (st.convergence_cfg.metric == lbw::CalibMetric::MAX_ERR) {
    if (errRp < best_metric_seen) {
        best_metric_seen = errRp;
        best_iter_val    = iter;
        for (int c = 0; c < ct.M_cell; c++)
            W_best[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
    }
}
```

**Block 2 — all other metrics (INSIDE `if (need_extra_metrics && W_total > 0.0)`, at the end):**
After all 6 metrics are computed, update best when active metric is NOT MAX_ERR:
```cpp
// BLOCK 2: best-iterate for non-MAX_ERR metrics (all values valid here).
if (st.convergence_cfg.metric != lbw::CalibMetric::MAX_ERR) {
    const double curr_best = lbw::select_metric(
        st.convergence_cfg.metric,
        errRp, mean_err, kl_max, chi2_total, grake_norm, l1_weight);
    if (std::isfinite(curr_best) && curr_best < best_metric_seen) {
        best_metric_seen = curr_best;
        best_iter_val    = iter;
        for (int c = 0; c < ct.M_cell; c++)
            W_best[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
    }
}
```

This ensures:
- MAX_ERR: best tracked at every kErrCheckInterval (same frequency as before)
- All other metrics: best tracked when need_extra_metrics fires (convergence approach, final iter, about_to_converge, or explicit non-MAX_ERR criterion)

- [ ] **Step 4: Fix `src/raking.cpp`**

Find and rename `best_errRp_seen` → `best_metric_seen` (lines ~131–132, ~207–209, ~426).

Find the current update at lines ~207–209:
```cpp
// OLD:
if (errRp < best_errRp_seen) {
    best_errRp_seen = errRp;
    best_iter_val   = iter;
    for (int i ...) W_best_r[i] = w[i];  // obs-level copy
}
```

Replace with the same two-block approach as ieppa.cpp. In raking.cpp the `need_extra_metrics` gate is at line ~251 (`if (need_extra_metrics && W_total > 0.0) {`). The current update at ~207 is BEFORE this gate (errRp is always available there).

**Block 1 — MAX_ERR, at the old line ~207 location (before need_extra_metrics gate):**
```cpp
if (st.convergence_cfg.metric == lbw::CalibMetric::MAX_ERR) {
    if (errRp < best_metric_seen) {
        best_metric_seen = errRp;
        best_iter_val    = iter;
        for (int i = 0; i < st.n; i++) W_best_r[i] = w[i];
    }
}
```

**Block 2 — non-MAX_ERR, INSIDE `if (need_extra_metrics && W_total > 0.0)` at the end of that block:**
```cpp
if (st.convergence_cfg.metric != lbw::CalibMetric::MAX_ERR) {
    const double curr_best = lbw::select_metric(
        st.convergence_cfg.metric,
        errRp, mean_err, kl_max, chi2_total, grake_norm, l1_weight);
    if (std::isfinite(curr_best) && curr_best < best_metric_seen) {
        best_metric_seen = curr_best;
        best_iter_val    = iter;
        for (int i = 0; i < st.n; i++) W_best_r[i] = w[i];
    }
}
```

At exit (~line 426): `res.best_error = best_metric_seen;` (replace `res.best_error = best_errRp_seen;`)

- [ ] **Step 5: Build gate**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

- [ ] **Step 6: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 332.

Note: existing A5 test `expect_lte(result$best_error, result$max_error)` — this must still hold for MAX_ERR criterion. For other metrics, `best_error` is now the minimum of that metric (which could be > or < max_error). Verify existing A5 tests use MAX_ERR or update expectations.

- [ ] **Step 7: Commit**
```bash
git add src/ieppa.cpp src/raking.cpp tests/testthat/test-best-iterate.R
git commit -m "$(cat <<'EOF'
fix(best-iter): track active metric not always errRp

best_errRp_seen renamed best_metric_seen; updated via select_metric()
using the active CalibMetric. MAX_ERR updates at every kErrCheckInterval;
other metrics update only when need_extra_metrics=true (values valid).
best_error now reflects min of active criterion; best_weights is the
iterate that minimizes the user's chosen objective.
EOF
)"
```

- [ ] **Step 8: Close ticket**
```bash
bd close leafblower-g4oj
```

---

## Final Verification

- [ ] For MAX_ERR criterion: `best_error <= max_error` (unchanged semantics)
- [ ] For KL criterion: `best_error` is KL value, not errRp
- [ ] FAIL 0, PASS ≥ 332
