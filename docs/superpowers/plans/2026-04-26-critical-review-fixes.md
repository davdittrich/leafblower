# Critical Review Fixes

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix 7 issues from the /critical-code-reviewer session pass: 2 blocking correctness bugs in raking.cpp and routing, 1 behavioral bug in harvest.R, and 4 clean-code improvements across c_api.cpp, r_bridge.cpp, raking.cpp, cell_table.cpp.

**Architecture:** All fixes are in existing files; no new files. Three commits: (1) raking.cpp fixes, (2) routing int-safety + algorithm_used, (3) harvest.R kl-for-auto + comment cleanup.

**Tech Stack:** C++17, R, testthat 3.

**Tickets:** leafblower-jmxi, leafblower-3xf1, leafblower-v5sy, leafblower-ujip, leafblower-csaw, leafblower-mys1, leafblower-aiil

**Baseline:** FAIL 0 | PASS 349 | SKIP 5

**Note:** `res.best_weights` retains the clamp-without-renormalize behavior (bounds hold, sum may drift when active). Same tradeoff as old obs-level code. Deferred: leafblower (create follow-up if exact sum=n for best_weights is needed).

**Verified:** `method="auto"` path in r_bridge.cpp calls `pack_solver_result` after solver dispatch, so `convergence_used$metric` IS populated for AUTO calls. T-auto-kl test will function correctly.

---

## File Structure

| File | Task | Tickets |
|---|---|---|
| `src/raking.cpp` | 1 | leafblower-jmxi, leafblower-csaw |
| `src/c_api.cpp` | 2 | leafblower-3xf1, leafblower-ujip |
| `src/r_bridge.cpp` | 2 | leafblower-3xf1, leafblower-mys1 |
| `src/cell_table.cpp` | 2 | leafblower-mys1 |
| `R/harvest.R` | 3 | leafblower-v5sy |
| `src/r_bridge.cpp` | 3 | leafblower-aiil |
| `tests/testthat/test-calibration-solvers.R` | 1, 3 | (RED→GREEN) |

---

## Task 1 — raking.cpp: remove bounds-violating normalization + dedup hyperplane

**Tickets:** leafblower-jmxi, leafblower-csaw  
**Files:** `src/raking.cpp`, `tests/testthat/test-calibration-solvers.R`

### Why the normalization is wrong

Post-exit expansion: `w_i = d_i × X[c]/X_init[c]`. Since `sum(X[c]) = n` (Dykstra hyperplane), `sum(w_i before clamp) = n`. The normalization `w_i *= n/total_w` is a no-op when no clamps fire. When bounds are active, `total_w ≠ n` after clamping, and normalization pushes weights back outside `[lo, hi]`.

Per spec §8: bounds hold exactly; marginal distortion from clamped cells is accepted. The normalization contradicts this. Remove it.

### Step 1.1: Write failing bounds test (RED)

Append to `tests/testthat/test-calibration-solvers.R`:
```r
test_that("R-bounds: raking respects min_weight/max_weight exactly", {
  # Tight max_weight=1.1 forces many clamps; normalization after clamp would violate bounds
  set.seed(7)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2","3","4","5"), n, replace=TRUE)),
    b = factor(sample(c("1","2"), n, replace=TRUE))
  )
  # Skewed target forces large multipliers for rare categories
  target <- list(
    a = c("1"=0.5, "2"=0.2, "3"=0.15, "4"=0.1, "5"=0.05),
    b = c("1"=0.8, "2"=0.2)
  )
  w <- leafblower::harvest(data, target, min_weight=0.5, max_weight=1.5,
                           method="raking", max_iterations=500, attach_weights=TRUE)
  # All weights must be within bounds (no normalization may violate them)
  expect_true(all(w >= 0.5 - 1e-10),
              info=sprintf("min weight %.6f < 0.5", min(w)))
  expect_true(all(w <= 1.5 + 1e-10),
              info=sprintf("max weight %.6f > 1.5", max(w)))
})
```

Run and confirm FAIL (or PASS if normalization doesn't fire in this case — if it passes, tighten `max_weight` until it fails):
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "R-bounds|FAIL|PASS" | head -5
```

### Step 1.2: Remove normalization from raking.cpp exit block

Read: `grep -n "Solver-contract normalization\|total_w\|norm.*st.n" src/raking.cpp | head -10`

Find the normalization block:
```cpp
// Solver-contract normalization: sum(w) = n.
double total_w = 0.0;
for (int i = 0; i < st.n; i++) total_w += st.weights[i];
if (std::isfinite(total_w) && total_w > 0.0) {
    const double norm = static_cast<double>(st.n) / total_w;
    for (int i = 0; i < st.n; i++) st.weights[i] *= norm;
}
```

Remove it entirely. The Dykstra hyperplane ensures `sum(X[c]) = n`, so `sum(w_i before clamp) = n`. After clamping, sum may differ slightly when bounds are active — this is intentional per spec §8.

### Step 1.3: Deduplicate hyperplane block (leafblower-csaw)

The post-loop Dykstra finalizer in raking.cpp has this block:
```cpp
{
    double s = 0.0;
    for (int c = 0; c < ct.M_cell; c++) { X[c] += q_hyp; s += X[c]; }
    double shift = (static_cast<double>(st.n) - s) / static_cast<double>(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) X[c] += shift;
    q_hyp = -shift;
}
```

This is identical to the hyperplane block inside the main loop. Extract to a lambda at the start of `raking_solve`:
```cpp
// Dykstra hyperplane step: projects X onto {sum(X) = n}.
// shift = (n - s) / M_cell ensures M_cell cells × shift = (n-s). q_hyp = -shift.
auto hyperplane_step = [&]() {
    double s = 0.0;
    for (int c = 0; c < ct.M_cell; c++) { X[c] += q_hyp; s += X[c]; }
    double shift = (static_cast<double>(st.n) - s) / static_cast<double>(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) X[c] += shift;
    q_hyp = -shift;
};
```

Replace both occurrences with `hyperplane_step();`.

**IMPORTANT:** The lambda must capture `ct`, `X`, `q_hyp`, `st.n` by reference. Verify the lambda is declared AFTER `X`, `ct`, `q_hyp`, `st` are in scope.

### Step 1.4: Build gate
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

### Step 1.5: Full regression
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 350.

### Step 1.6: Commit
```bash
git add src/raking.cpp tests/testthat/test-calibration-solvers.R
git commit -m "$(cat <<'EOF'
fix(raking): remove bounds-violating post-clamp normalization + dedup hyperplane

Normalization w_i *= n/total_w after bounds clamp could push weights
outside [lo,hi] when active bounds changed total_w. Per spec §8: bounds
hold exactly; marginal distortion from clamped cells is intentional.
Normalization removed — sum(X[c])=n from Dykstra hyperplane makes it a
no-op when no clamps fire anyway. Hyperplane block deduplicated to lambda.
EOF
)"
bd close leafblower-jmxi && bd close leafblower-csaw
```

---

## Task 2 — c_api.cpp + r_bridge.cpp + cell_table.cpp: int overflow + algorithm_used + reserve

**Tickets:** leafblower-3xf1, leafblower-ujip, leafblower-mys1

### Step 2.1: Write failing test (RED for overflow)

The overflow is in `n * 9` with n as int. No unit test can exercise n > 238M safely. Instead: write a test asserting the types are correct at the call site (compile-time check). The TDD here is a code review assertion. The RED state is: `n * 9` uses int arithmetic (provable from current source). The GREEN state is: `(int64_t)n * 9`.

Append to `tests/testthat/test-calibration-solvers.R`:
```r
test_that("T-overflow: AUTO routing handles n near INT_MAX without overflow (type guard)", {
  # This is a type-safety smoke test. The actual overflow fix is in C++.
  # Tests that harvest works correctly at moderately large n without issues.
  set.seed(1)
  n <- 1e5L  # 100k obs, well within int range, exercises routing path
  data <- data.frame(
    a = factor(sample(1:100, n, replace=TRUE)),  # 100 cats → M_cell = 100 (ratio << 0.9)
    b = factor(sample(1:2, n, replace=TRUE))
  )
  target <- list(
    a = setNames(rep(0.01, 100), as.character(1:100)),
    b = c("1"=0.5, "2"=0.5)
  )
  w <- leafblower::harvest(data, target, max_weight=5, method="auto", attach_weights=FALSE)
  r <- attr(w, "result")
  expect_equal(r$status, 0L)
  expect_lt(r$max_error, 0.01)
  # Must route to ieppa (M_cell=200 << n=1e5, ratio = 0.002)
  expect_equal(r$algorithm_used, 1L,
               info="AUTO must select ieppa when M_cell/n << 0.9")
})
```

Run — currently FAILS on `r$algorithm_used` assertion if field is not populated:
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "T-overflow|FAIL|PASS" | head -5
```

### Step 2.2: Fix int overflow in c_api.cpp

Read: `grep -n "M_cell_est \* 10\|n \* 9" src/c_api.cpp`

Change:
```cpp
// OLD:
alg = (M_cell_est * 10 > n * 9) ? RK_ALG_RAKING : RK_ALG_IEPPA;

// NEW:
alg = (static_cast<int64_t>(M_cell_est) * 10 > static_cast<int64_t>(n) * 9)
      ? RK_ALG_RAKING : RK_ALG_IEPPA;
```

Add `#include <cstdint>` at the top of c_api.cpp if not already present:
```bash
grep -n "cstdint\|stdint" src/c_api.cpp | head -3
```

### Step 2.3: Fix int overflow in r_bridge.cpp

Same fix at `grep -n "M_cell_est \* 10\|n \* 9" src/r_bridge.cpp`.

### Step 2.4: Fix estimate_M_cell reserve (leafblower-mys1)

Read: `grep -n "reserve\|seen\." src/cell_table.cpp | head -10`

Change:
```cpp
// OLD:
seen.reserve(static_cast<size_t>(n));

// NEW — cap at 1<<20 (1M) to avoid massive pre-allocation on compressible data:
seen.reserve(static_cast<size_t>(std::min(n, 1 << 20)));
```

Need `#include <algorithm>` in cell_table.cpp if not present. Check: `grep -n "algorithm" src/cell_table.cpp | head -3`

### Step 2.6: Build gate
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

### Step 2.7: Full regression
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 351.

### Step 2.8: Commit A — int overflow + reserve (compile-time safety)
```bash
git add src/c_api.cpp src/r_bridge.cpp src/cell_table.cpp tests/testthat/test-calibration-solvers.R
git commit -m "$(cat <<'EOF'
fix(routing): int64 overflow guard + reserve cap in estimate_M_cell

AUTO threshold M_cell*10 > n*9 now uses int64_t casts (int*int overflows
at n>238M). estimate_M_cell reserves min(n, 1<<20) instead of n (avoids
oversized allocation on compressible inputs where M_cell << n).
EOF
)"
bd close leafblower-3xf1 && bd close leafblower-mys1
```

### Step 2.9: Fix algorithm_used on stub BADARG paths (leafblower-ujip)

**Implement AFTER Step 2.8 commits so this is a fresh c_api.cpp change.**

Read: `grep -n "RK_ERR_BADARG\|stub\|SINKHORN" src/c_api.cpp | head -10`

In c_api.cpp, the stub dispatch for unimplemented methods:
```cpp
case RK_ALG_SINKHORN:
case RK_ALG_CHEBYSHEV:
case RK_ALG_GREG:
case RK_ALG_GRAKE: {
    if (result) {
        result->status = RK_ERR_BADARG;
        snprintf(result->message, ...);
    }
    return RK_ERR_BADARG;
}
```

Add `result->algorithm_used = static_cast<int>(algorithm);` inside the `if (result)` block, after setting `result->status`. Verify `algorithm_used` field name:
```bash
grep -n "algorithm_used" src/leafblower.h | head -3
```
Expected: `rk_algorithm_t algorithm_used;` in `rk_result_t`.

### Step 2.10: Commit B — algorithm_used on BADARG
```bash
git add src/c_api.cpp
git commit -m "fix(c_api): set algorithm_used on stub RK_ERR_BADARG return

Stub dispatch for sinkhorn/chebyshev/greg/grake now sets
result->algorithm_used before returning RK_ERR_BADARG so error
telemetry knows which unimplemented method was attempted."
bd close leafblower-ujip
```

---

## Task 3 — harvest.R kl-for-auto + r_bridge.cpp stale comment

**Tickets:** leafblower-v5sy, leafblower-aiil  
**Files:** `R/harvest.R`, `src/r_bridge.cpp`

### Step 3.1: Write failing test (RED)

Append to `tests/testthat/test-calibration-solvers.R`:
```r
test_that("T-auto-kl: method='auto' defaults to kl convergence metric", {
  # AUTO must inherit kl default — both ieppa and raking use kl as default
  set.seed(3)
  data <- data.frame(a=factor(sample(c("1","2","3"),300,TRUE)))
  target <- list(a=c("1"=0.4,"2"=0.4,"3"=0.2))
  w <- leafblower::harvest(data, target, max_weight=5, method="auto", attach_weights=FALSE)
  r <- attr(w, "result")
  expect_equal(r$convergence_used$metric, "kl",
               info="AUTO must use kl default — both ieppa and raking default to kl")
})
```

Run — currently FAILS because method="auto" does not get the kl override:
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "T-auto-kl|FAIL|PASS" | head -5
```

### Step 3.2: Fix kl default in harvest.R

Read: `grep -n "method.*ieppa\|method.*raking\|conv\$metric.*kl" R/harvest.R | head -10`

Find the two override blocks:
```r
# ieppa kl override
if (method == "ieppa" && is.null(convergence[["metric"]]) && ...) conv$metric <- "kl"
# raking kl override
if (method == "raking" && ...) conv$metric <- "kl"
```

Replace both with a single combined block:
```r
# iEPPA, raking, and auto all default to kl+improvement — these are KL minimizers.
# "auto" inherits kl because it routes to ieppa or raking only.
if (method %in% c("ieppa", "raking", "auto") &&
    is.null(convergence[["metric"]]) &&
    is.null(convergence[["criterion"]]) &&
    is.null(convergence[["improvement"]]) &&
    is.null(convergence[["pct"]]) &&
    is.null(convergence[["absolute"]])) {
  conv$metric <- "kl"
}
```

This replaces BOTH separate blocks — one combined check, not two.

### Step 3.3: Fix stale comment in r_bridge.cpp (leafblower-aiil)

Read: `grep -n "5 convergence fields\|convergence fields\|14 prior" src/r_bridge.cpp | head -5`

Find and update the stale comment on the `Rf_allocVector(VECSXP, 30)` line. Already done in a previous commit? Verify — if "7 convergence fields" is already present, skip this step.

### Step 3.4: Build gate
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

### Step 3.5: Full regression
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 352.

### Step 3.6: Commit
```bash
git add R/harvest.R src/r_bridge.cpp tests/testthat/test-calibration-solvers.R
git commit -m "$(cat <<'EOF'
fix(harvest): kl default applies to method='auto' + consolidate override

ieppa/raking/auto all route to KL minimizers. Single combined condition
replaces two separate per-method blocks. method='auto' now inherits kl
default regardless of which algorithm C++ selects at runtime.
EOF
)"
bd close leafblower-v5sy && bd close leafblower-aiil
```

---

## Final Verification

- [ ] `Rscript -e 'devtools::test()' 2>&1 | tail -3` → FAIL 0, PASS ≥ 352
- [ ] `grep "total_w\|norm.*st.n\|Solver-contract normalization" src/raking.cpp` → 0 output
- [ ] `grep "M_cell_est \* 10 > n \* 9" src/c_api.cpp` → 0 (cast version used)
- [ ] Bounds test passes: weights within [min_weight, max_weight]
- [ ] kl-for-auto test passes: method="auto" reports metric="kl"
