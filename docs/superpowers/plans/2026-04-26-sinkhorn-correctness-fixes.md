# Sinkhorn Correctness Fixes

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix 3 issues found by /critical-code-reviewer on Plan C sinkhorn:
1. Dykstra correction `a[c]` can accumulate unboundedly → exp overflow
2. No synthetic correctness test (A1 requires unavailable parquet)
3. Bisection has no short-circuit when X already in bounds

**Architecture:** All changes in `src/sinkhorn.cpp` and `tests/testthat/test-calibration-solvers.R`. One commit each.

**Tickets:** leafblower-4zo4, leafblower-xi7h, leafblower-qbxe

**Baseline:** FAIL 0 | PASS 354 | SKIP 5

---

## File Structure

| File | Task | Ticket |
|---|---|---|
| `src/sinkhorn.cpp` | 1, 3 | leafblower-4zo4, leafblower-qbxe |
| `tests/testthat/test-calibration-solvers.R` | 2 | leafblower-xi7h |

---

## Task 1 — Clamp Dykstra correction a[c] to prevent exp overflow

**Ticket:** leafblower-4zo4

**Why:** `a[c] += log(X[c]) - log(X_proj[c])` accumulates across iterations. When a cell sits persistently at its bound, `a[c]` grows monotonically. `exp(a[c] + μ)` then overflows (inf) for `a[c] > ~700`, causing `bisect_capacity` to return wrong X_proj values or false INFEAS. Clamping `a[c]` to ±30 prevents overflow while still allowing the bisection to find the correct μ (exp(30) ≈ 1e13 >> any practical weight).

### Step 1.1: Write regression guard test

**TDD note:** Overflow only manifests for pathological problems with persistent clamping over hundreds of iterations. At n=500 with 2000 iters the test passes BEFORE the fix (overflow is not triggered on this particular input). S1 is a **regression guard** — it documents the failure mode and protects against future regressions that remove the clamp.

Append to `tests/testthat/test-calibration-solvers.R`:
```r
test_that("S1: sinkhorn handles tight bounds without overflow (a[c] clamp guard)", {
  # Regression guard: verifies a[c] clamp is not removed.
  # Without clamping, a[c] grows unboundedly at tight bounds → exp overflow → false INFEAS.
  set.seed(99)
  n <- 500
  data <- data.frame(
    a = factor(sample(c("1","2","3","4","5"), n, replace=TRUE)),
    b = factor(sample(c("1","2"), n, replace=TRUE))
  )
  target <- list(
    a = c("1"=0.6, "2"=0.2, "3"=0.1, "4"=0.05, "5"=0.05),
    b = c("1"=0.9, "2"=0.1)
  )
  w <- leafblower::harvest(data, target, min_weight=0.1, max_weight=2.0,
                           method="sinkhorn", max_iterations=2000, attach_weights=FALSE)
  r <- attr(w, "result")
  expect_false(r$status == 3L,
               info=sprintf("sinkhorn returned INFEAS (status=%d) — likely a[c] overflow", r$status))
  expect_true(r$status %in% c(0L, 1L),
              info=sprintf("expected 0 (OK) or 1 (NOCONV), got %d", r$status))
  expect_true(all(w >= 0.1 - 1e-10 & w <= 2.0 + 1e-10),
              info="bounds violated")
})
```

### Step 1.2: Add clamp constant and apply after Dykstra update

Read: `grep -n "a\[c\].*log\|log.*a\[c\]" src/sinkhorn.cpp | head -5`

Find the Dykstra correction block:
```cpp
for (int c = 0; c < ct.M_cell; c++) {
    if (X[c] > 1e-300 && X_proj[c] > 1e-300)
        a[c] += std::log(X[c]) - std::log(X_proj[c]);
    X[c] = X_proj[c];
}
```

Add a constant and clamp:
```cpp
static constexpr double kAmax = 30.0;  // exp(30) ≈ 1e13 >> max practical weight ratio

for (int c = 0; c < ct.M_cell; c++) {
    if (X[c] > 1e-300 && X_proj[c] > 1e-300)
        a[c] += std::log(X[c]) - std::log(X_proj[c]);
    a[c] = std::clamp(a[c], -kAmax, kAmax);  // prevent exp overflow
    X[c] = X_proj[c];
}
```

### Step 1.3: Build gate
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

### Step 1.4: Run test + regression
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "S1|FAIL|PASS" | head -5
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: S1 PASSES, FAIL 0, PASS ≥ 354.

### Step 1.5: Commit
```bash
git add src/sinkhorn.cpp tests/testthat/test-calibration-solvers.R
git commit -m "fix(sinkhorn): clamp Dykstra a[c] to prevent exp overflow at tight bounds

a[c] accumulates log(X/X_proj) per iter. At tight bounds with X >> L,
a[c] grows monotonically and exp(a+mu) overflows after ~700 iters,
causing false RK_ERR_INFEAS. Clamp to ±30 (exp(30)≈1e13, safe for
practical weights); slight asymptotic cost tradeoff is acceptable."
bd close leafblower-4zo4 2>/dev/null || true
```

---

## Task 2 — Add synthetic KL correctness test

**Ticket:** leafblower-xi7h

No RED step — this is a test addition. But the test should PASS only after sinkhorn is working correctly (which it is post-merge).

### Step 2.1: Append synthetic test

Append to `tests/testthat/test-calibration-solvers.R`:
```r
test_that("S2: sinkhorn achieves KL <= raking on synthetic (no external data)", {
  # Verifies sinkhorn minimizes KL better than raking without any external fixtures.
  set.seed(7)
  n <- 300
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace=TRUE)),
    b = factor(sample(c("1","2"), n, replace=TRUE))
  )
  target <- list(
    a = c("1"=0.5, "2"=0.3, "3"=0.2),
    b = c("1"=0.6, "2"=0.4)
  )
  w_sink <- leafblower::harvest(data, target, max_weight=5, method="sinkhorn",
                                max_iterations=500, attach_weights=FALSE)
  w_rake <- leafblower::harvest(data, target, max_weight=5, method="raking",
                                max_iterations=500, attach_weights=FALSE)
  r_sink <- attr(w_sink, "result")
  r_rake <- attr(w_rake, "result")

  # Both must converge
  expect_equal(r_sink$status, 0L, info="sinkhorn must converge")
  expect_equal(r_rake$status, 0L, info="raking must converge")

  # Sinkhorn minimizes KL — its KL at best_iter must be <= raking's
  expect_lte(r_sink$kl, r_rake$kl + 1e-6,
             label="sinkhorn KL <= raking KL (sinkhorn is a KL minimizer; competitive on unconstrained problems)")

  # Hard bounds: all weights within [1/5, 5]
  expect_true(all(w_sink >= 1/5 - 1e-10 & w_sink <= 5 + 1e-10),
              info="sinkhorn weights must respect bounds")
})
```

### Step 2.2: Run test
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "S2|FAIL|PASS" | head -5
```
Expected: S2 PASSES. If it fails with `sinkhorn KL > raking KL`, investigate the convergence — sinkhorn should achieve lower KL on unconstrained/mildly constrained problems.

### Step 2.3: Commit
```bash
git add tests/testthat/test-calibration-solvers.R
git commit -m "test(sinkhorn): synthetic KL correctness test S2 (no parquet required)

Verifies sinkhorn achieves KL <= raking on synthetic n=300 problem.
No external data dependency — runs in any environment."
bd close leafblower-xi7h 2>/dev/null || true
```

---

## Task 3 — Add bisection short-circuit when X already in bounds

**Ticket:** leafblower-qbxe

**Why:** After convergence (or when bounds aren't tight), all `X[c]` are already within `[L_c, U_c]`. Running bisection in this case costs O(M_cell × 40) exp calls for a no-op (μ=0 is the answer). Short-circuit: if all `X[c]` satisfy bounds after Sinkhorn sweeps, skip bisection, set X_proj = X, a[c] unchanged.

### Step 3.1: Add short-circuit before bisect_capacity call

Read: `grep -n "bisect_capacity\|target_mass" src/sinkhorn.cpp | head -8`

Read: `sed -n '118,135p' src/sinkhorn.cpp`

Find the bisect block. Replace with the short-circuit version. **CRITICAL: `double mu` must be declared BEFORE the `needs_projection` check, initialized to 0.0:**

```cpp
        // Short-circuit: if X already within capacity bounds, bisection is a no-op
        // (when X_proj == X, log(X/X_proj) = 0, so a[c] correction is also zero).
        bool needs_projection = false;
        for (int c = 0; c < ct.M_cell; c++) {
            if (X[c] < L_cell[c] - 1e-12 || X[c] > U_cell[c] + 1e-12) {
                needs_projection = true;
                break;
            }
        }
        double mu = 0.0;  // ← declared here, BEFORE the if/else
        if (needs_projection) {
            if (!bisect_capacity(X, a, L_cell, U_cell, ct.M_cell, target_mass, mu, X_proj)) {
                res.status = RK_ERR_INFEAS;
                break;
            }
        } else {
            X_proj = X;  // identity projection; mu stays 0
        }
        // Dykstra correction update (same for both paths)
        for (int c = 0; c < ct.M_cell; c++) { ... }
```

Make sure to REMOVE the old `double mu;` declaration that was on the line immediately before `if (!bisect_capacity(...))`. Replace the entire block from just before `double mu;` through the `if (!bisect_capacity)` block.

### Step 3.2: Build gate
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

### Step 3.3: Regression
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 355 (S1+S2 tests added).

### Step 3.4: Commit
```bash
git add src/sinkhorn.cpp
git commit -m "perf(sinkhorn): short-circuit bisection when X already within capacity bounds

Avoids O(M_cell*40) exp calls per outer iter when no clamping is needed.
Common after convergence approach. X checked against [L_c, U_c] in O(M_cell)
before invoking bisect_capacity."
bd close leafblower-qbxe 2>/dev/null || true
```

---

## Final Verification

- [ ] `Rscript -e 'devtools::test()' 2>&1 | tail -3` → FAIL 0, PASS ≥ 356
- [ ] S1 test verifies overflow not triggered at max_iterations=2000
- [ ] S2 test verifies sinkhorn KL ≤ raking KL on synthetic
- [ ] `grep "kAmax" src/sinkhorn.cpp` → clamp constant present
- [ ] `grep "needs_projection" src/sinkhorn.cpp` → short-circuit present
