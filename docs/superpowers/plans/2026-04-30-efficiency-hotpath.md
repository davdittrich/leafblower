# Hot-Path Efficiency Fixes — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans.

**Goal:** Eliminate per-iteration heap allocations and redundant O(K*M_cell) margin passes in the hot paths of ieppa, raking, lbfgsb, and sinkhorn solvers.

**Architecture:** Each task is an isolated mechanical fix (hoist alloc / guard computation / fuse loops). Algorithm logic does not change — only data layout and evaluation frequency. Ordered by estimated impact: lbfgsb/raking/ieppa HIGH items first, MEDIUM items last.

**Tech Stack:** C++17, R package build, testthat, benchmarks/stepstone_benchmark.R for perf verification

---

## Mechanism / Forbidden / Audit

- **Mechanism:** Hoist-out-of-loop allocs, single-pass metric fusion, guarded evaluation frequency
- **Forbidden:** Any change to convergence logic, loop bounds, solver math, or result struct fields. No new abstractions. No algorithmic substitutions.
- **Audit:** Output-identity tests: capture `harvest()` weights BEFORE and AFTER each change; assert `max(abs(before - after)) < 1e-10` on the stepstone_small fixture.

---

## Shared Test Fixture (used by ALL tasks)

### Step 0: Capture reference weights BEFORE any code change

The fixture is `tests/testthat/fixtures/stepstone_small.parquet` (parquet, not rds).
The companion targets file is `tests/testthat/fixtures/stepstone_small_targets.rds`.

**Before touching any C++ source**, run this capture script once and commit the `.rds` output:

```r
# capture_efficiency_ref.R — run ONCE before any code change, then git add + commit
library(leafblower)
library(arrow)

fx <- testthat::test_path("fixtures", "stepstone_small.parquet")
tg <- testthat::test_path("fixtures", "stepstone_small_targets.rds")
data   <- arrow::read_parquet(fx)
target <- readRDS(tg)

ref <- harvest(data = data, target = target, method = "lbfgsb")$weights
saveRDS(ref, testthat::test_path("fixtures", "lbfgsb_efficiency_ref.rds"))
cat("Reference captured:", length(ref), "weights\n")
```

```bash
# After running the script:
git add tests/testthat/fixtures/lbfgsb_efficiency_ref.rds
git commit -m "test(efficiency-hotpath): capture lbfgsb reference weights before fixes"
```

> **NOTE:** The reference `.rds` file must be committed **before ANY C++ source changes**. Tests skip automatically if the file is absent.

### Step 0 (also before touching any C++ source): Measure OLD binary performance

Run this immediately after capturing the reference weights, while the unmodified binary is still installed:

```bash
# Step 0: Measure OLD binary performance baseline (BEFORE any code change)
Rscript -e '
  library(bench); library(leafblower); library(arrow)
  fx <- "tests/testthat/fixtures/stepstone_small.parquet"
  tg <- "tests/testthat/fixtures/stepstone_small_targets.rds"
  data   <- arrow::read_parquet(fx)
  target <- readRDS(tg)
  cat("OLD lbfgsb:", format(bench::mark(
    harvest(data=data, target=target, method="lbfgsb"),
    iterations=20)$median), "\n")
  cat("OLD raking:", format(bench::mark(
    harvest(data=data, target=target, method="raking"),
    iterations=20)$median), "\n")
  cat("OLD ieppa:", format(bench::mark(
    harvest(data=data, target=target, method="ieppa",
            convergence=list(metric="MEAN_ERR", absolute=1e-6)),
    iterations=20)$median), "\n")
'
# Record OLD times. THEN apply code changes and rebuild.
```

Record the printed medians. After all tasks are applied and rebuilt, re-run the same script with the new binary. See the BENCHMARKING section for the full before/after protocol.

### Test file: `tests/testthat/test-efficiency-hotpath.R`

Add this test. It skips automatically until the reference is present, and will fail
only after code changes have broken numerical identity.

```r
test_that("lbfgsb output identity after efficiency fixes", {
  skip_if_not_installed("arrow")
  ref_path <- testthat::test_path("fixtures", "lbfgsb_efficiency_ref.rds")
  skip_if_not(file.exists(ref_path),
    "reference not yet captured — run capture_efficiency_ref.R before any code change")

  fx <- testthat::test_path("fixtures", "stepstone_small.parquet")
  tg <- testthat::test_path("fixtures", "stepstone_small_targets.rds")
  skip_if(!file.exists(fx) || !file.exists(tg))

  data   <- arrow::read_parquet(fx)
  target <- readRDS(tg)
  ref    <- readRDS(ref_path)

  got <- harvest(data = data, target = target, method = "lbfgsb")$weights
  expect_equal(got, ref, tolerance = 1e-10)
})
```

#### Raking metrics identity test (for Task 4)

Capture script — run once before any C++ change, then commit:

```r
# capture_raking_metrics_ref.R — run ONCE before Task 4, then git add + commit
library(leafblower); library(arrow)
fx <- testthat::test_path("fixtures", "stepstone_small.parquet")
tg <- testthat::test_path("fixtures", "stepstone_small_targets.rds")
data   <- arrow::read_parquet(fx)
target <- readRDS(tg)
r <- harvest(data = data, target = target, method = "raking")
saveRDS(list(max_error = r$max_error, iterations = r$iterations,
             mean_error = r$mean_error),
        testthat::test_path("fixtures", "raking_metrics_ref.rds"))
cat("Raking metrics captured\n")
```

```bash
git add tests/testthat/fixtures/raking_metrics_ref.rds
git commit -m "test(efficiency-hotpath): capture raking metrics reference (pre-Task4 baseline)"
```

Add to `tests/testthat/test-efficiency-hotpath.R`:

```r
test_that("raking convergence metrics unchanged after F_eval guard", {
  skip_if_not_installed("arrow")
  ref_path <- testthat::test_path("fixtures", "raking_metrics_ref.rds")
  skip_if_not(file.exists(ref_path),
    "raking metrics reference not captured — run capture_raking_metrics_ref.R before Task 4")
  fx <- testthat::test_path("fixtures", "stepstone_small.parquet")
  tg <- testthat::test_path("fixtures", "stepstone_small_targets.rds")
  skip_if(!file.exists(fx) || !file.exists(tg))
  data   <- arrow::read_parquet(fx)
  target <- readRDS(tg)
  ref    <- readRDS(ref_path)
  r      <- harvest(data = data, target = target, method = "raking")
  expect_equal(r$max_error,   ref$max_error,   tolerance = 1e-10)
  expect_equal(r$iterations,  ref$iterations)
  expect_equal(r$mean_error,  ref$mean_error,  tolerance = 1e-10)
})
```

> The reference `ref` is captured from the pre-change binary and committed. Each task's identity test runs the same `harvest()` call on the new binary and diffs against the committed reference. The test is not self-referential because `ref` comes from a different binary (pre-change .rds file).

---

## Task 1 (773f.1): Hoist S/S2 bucket allocs outside K-loops in `compute_final_weights_and_error`

**File:** `src/lbfgsb_solver.cpp`  
**Impact:** HIGH — two `vector<double>` allocs per K-iteration in a function called once per solve, but with K potentially large.

### Root cause

`compute_final_weights_and_error` (line 164) runs two separate K-loops. The first loop (line 184) allocates `std::vector<double> S(st.cat_counts[k], 0.0)` per iteration. The second loop (line 235) allocates `std::vector<double> S2(st.cat_counts[k], 0.0)` per iteration. Both are scratch buffers of varying size; both can be replaced with a single pre-allocated buffer sized to `max(cat_counts)`.

### Before (line 184–192)

```cpp
for (int k = 0; k < st.K; k++) {
    std::vector<double> S(st.cat_counts[k], 0.0);
    for (int i = 0; i < st.n; i++) {
        int g = st.group_ids[k][i];
        if (g >= 0) S[g] += st.weights[i];
    }
    for (int j = 0; j < st.cat_counts[k]; j++) {
        max_err = std::max(max_err, std::fabs(S[j] / Wn - st.targets[k][j]));
    }
}
```

### Before (line 234–258)

```cpp
for (int k = 0; k < st.K; k++) {
    std::vector<double> S2(st.cat_counts[k], 0.0);
    for (int i = 0; i < st.n; i++) {
        int g = st.group_ids[k][i];
        if (g >= 0) S2[g] += st.weights[i];
    }
    double max_k = 0.0;
    double kl_k  = 0.0;
    for (int j = 0; j < st.cat_counts[k]; j++) {
        double S_p = S2[j] / Wn;
        ...
        double obs    = S2[j];
        ...
    }
    ...
}
```

### After

Add before the `double Wn = 0.0;` line (line 181):

```cpp
// 773f.1: hoist bucket alloc outside K-loops — eliminates 2*K heap allocs per solve.
int max_cats = 0;
for (int k = 0; k < st.K; k++)
    if (st.cat_counts[k] > max_cats) max_cats = st.cat_counts[k];
std::vector<double> bucket(max_cats, 0.0);
```

First K-loop (replace `std::vector<double> S(...)` and all `S[g]`/`S[j]` refs):

```cpp
for (int k = 0; k < st.K; k++) {
    std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
    for (int i = 0; i < st.n; i++) {
        int g = st.group_ids[k][i];
        if (g >= 0) bucket[g] += st.weights[i];
    }
    for (int j = 0; j < st.cat_counts[k]; j++) {
        max_err = std::max(max_err, std::fabs(bucket[j] / Wn - st.targets[k][j]));
    }
}
```

Second K-loop (replace `std::vector<double> S2(...)` and all `S2[g]`/`S2[j]` refs):

```cpp
for (int k = 0; k < st.K; k++) {
    std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
    for (int i = 0; i < st.n; i++) {
        int g = st.group_ids[k][i];
        if (g >= 0) bucket[g] += st.weights[i];
    }
    double max_k = 0.0;
    double kl_k  = 0.0;
    for (int j = 0; j < st.cat_counts[k]; j++) {
        double S_p = bucket[j] / Wn;
        double T   = st.targets[k][j];
        double err = std::fabs(S_p - T);
        if (err > max_k) max_k = err;
        if (T > 0.0) {
            kl_k += T * std::log((T + kMetricEps) / (S_p + kMetricEps));
        }
        double obs    = bucket[j];
        double pop_kj = T * Wn;
        chi2_total += (obs - pop_kj) * (obs - pop_kj) / (pop_kj + kChi2Floor);
        double nm = std::fabs(obs - pop_kj) / (1.0 + std::fabs(pop_kj));
        if (nm > grake_norm) grake_norm = nm;
    }
    mean_err_sum += max_k;
    if (kl_k > kl_max) kl_max = kl_k;
}
```

### Unchanged components

`phi_and_grad`, `lbfgsb_solve_inner`, `wolfe_line_search`, `wolfe_zoom` — all untouched.

### Compile & test

```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean .
Rscript -e 'testthat::test_dir("tests/testthat")'
```

### Regressions prevented

No algorithm change. Both K-loops produce identical numerical output because `bucket` is `fill`ed to 0 before each `k`. The pre-allocated path cannot alias between loops since the second loop also fill-resets before use.

### Git commit message

```
perf(lbfgsb): hoist S/S2 bucket allocs outside K-loops in compute_final_weights_and_error

Eliminates 2*K heap allocations per solve by pre-sizing one scratch
vector to max(cat_counts) and std::fill-resetting per iteration.
No algorithmic change; output identity verified to 1e-10.
```

---

## Task 2 (773f.2): Hoist `log_denom_scratch` outside Wolfe bisection

**File:** `src/lbfgsb_solver.cpp`  
**Impact:** HIGH — `log_denom_scratch(fn.exponential ? 0 : st.n)` is allocated inside both `wolfe_zoom` (line 352) and `wolfe_line_search` (line 481), each of which is called up to 20 times per outer L-BFGS step and internally loops 20 more times. For logit link and large n, that is up to 400 heap allocs per outer iteration.

### Root cause

`wolfe_zoom` declares (line 352):
```cpp
std::vector<double> log_denom_scratch(fn.exponential ? 0 : st.n);
```

`wolfe_line_search` declares (line 481):
```cpp
std::vector<double> log_denom_scratch(fn.exponential ? 0 : st.n);
```

Both functions are called exclusively from `lbfgsb_solve_inner`. The buffer can be declared once there and passed by reference.

### Pre-implementation audit: confirm all wolfe_zoom call sites

Before implementation, run:

```bash
grep -n 'wolfe_zoom(' src/lbfgsb_solver.cpp
```

Expected output: **exactly 3 lines** — line 335 (the definition) plus lines 566 and 578 (both inside `wolfe_line_search`). No other call sites exist. If any additional call sites appear, they must also receive the `log_denom_scratch` parameter before proceeding.

### Call sites in `lbfgsb_solve_inner` (line 648, 566, 578)

```cpp
// Line 648 — primary call
wolfe_line_search(st, fn, off, T, d, lam, phi_curr, slope_0,
                  u, du, u_work, e_vec, dir, lam_new, grad_new, phi_new);

// Lines 566–569 (inside wolfe_line_search → calls wolfe_zoom)
return wolfe_zoom(st, fn, off, T, d, phi_0, slope_0,
                  alpha_prev, phi_prev, alpha,
                  u_base, du, lam, dir, u_work, e_vec,
                  lam_new, grad_new, phi_new);

// Lines 578–581 (second wolfe_zoom call site inside wolfe_line_search)
return wolfe_zoom(st, fn, off, T, d, phi_0, slope_0,
                  alpha, phi_trial, alpha_prev,
                  u_base, du, lam, dir, u_work, e_vec,
                  lam_new, grad_new, phi_new);
```

### After — function signature changes

`wolfe_zoom` — add `std::vector<double>& log_denom_scratch` as last parameter before `double& phi_new`:

```cpp
static double wolfe_zoom(
        const CalibState& st, const LinkFn& fn,
        const std::vector<int>& off, const std::vector<double>& T,
        const std::vector<double>& d, double phi_0, double slope_0,
        double alpha_lo, double phi_lo, double alpha_hi,
        const std::vector<double>& u_base, const std::vector<double>& du,
        const std::vector<double>& lam,
        const std::vector<double>& dir,
        std::vector<double>& u_work,
        std::vector<double>& e_vec,
        std::vector<double>& lam_new, std::vector<double>& grad_new,
        std::vector<double>& log_denom_scratch,   // 773f.2: hoisted
        double& phi_new) {
```

Remove the local declaration at line 352:
```cpp
// REMOVE:
std::vector<double> log_denom_scratch(fn.exponential ? 0 : st.n);
```

`wolfe_line_search` — add same parameter before `double& phi_new`:

```cpp
static double wolfe_line_search(
        const CalibState& st, const LinkFn& fn,
        const std::vector<int>& off, const std::vector<double>& T,
        const std::vector<double>& d,
        const std::vector<double>& lam, double phi_0, double slope_0,
        const std::vector<double>& u_base, const std::vector<double>& du,
        std::vector<double>& u_work,
        std::vector<double>& e_vec,
        const std::vector<double>& dir,
        std::vector<double>& lam_new, std::vector<double>& grad_new,
        std::vector<double>& log_denom_scratch,   // 773f.2: hoisted
        double& phi_new) {
```

Remove the local declaration at line 481:
```cpp
// REMOVE:
std::vector<double> log_denom_scratch(fn.exponential ? 0 : st.n);
```

Update the two `wolfe_zoom` call sites inside `wolfe_line_search` to pass `log_denom_scratch`:

```cpp
return wolfe_zoom(st, fn, off, T, d, phi_0, slope_0,
                  alpha_prev, phi_prev, alpha,
                  u_base, du, lam, dir, u_work, e_vec,
                  lam_new, grad_new, log_denom_scratch, phi_new);
// and:
return wolfe_zoom(st, fn, off, T, d, phi_0, slope_0,
                  alpha, phi_trial, alpha_prev,
                  u_base, du, lam, dir, u_work, e_vec,
                  lam_new, grad_new, log_denom_scratch, phi_new);
```

In `lbfgsb_solve_inner`, declare once after the existing scratch vectors (after line 618, before the iteration loop at line 629):

```cpp
// 773f.2: single allocation for logit log-denom scratch; passed into both Wolfe helpers.
std::vector<double> log_denom_scratch(fn.exponential ? 0 : st.n);
```

Update the `wolfe_line_search` call at line 648:

```cpp
wolfe_line_search(st, fn, off, T, d, lam, phi_curr, slope_0,
                  u, du, u_work, e_vec, dir, lam_new, grad_new,
                  log_denom_scratch, phi_new);
```

### Unchanged components

`phi_and_grad`, `phi_from_u`, `compute_du`, `compute_u`, `compute_final_weights_and_error` — untouched.

### Compile & test

```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean .
Rscript -e 'testthat::test_dir("tests/testthat")'
```

### Regressions prevented

`log_denom_scratch` is fully overwritten by `bulk_log` before any read in both functions. No state leaks between Wolfe trials — same guarantee as the local version. Exponential link path: `log_denom_scratch` has size 0 and is never accessed (existing guard `if (!fn.exponential)` before all use sites unchanged).

### Git commit message

```
perf(lbfgsb): hoist log_denom_scratch outside Wolfe bisection loops

Eliminates up to 400 heap allocs per outer L-BFGS step on logit link
by declaring log_denom_scratch once in lbfgsb_solve_inner and passing
by reference into wolfe_line_search and wolfe_zoom. No algorithmic
change; output identity verified to 1e-10.
```

---

## Task 3 (773f.4): Fuse ieppa extra-metrics sweep with errRp via `compute_cell_metrics`

**File:** `src/ieppa.cpp`  
**Impact:** HIGH — when `need_extra_metrics` is true, the convergence block (line 1122–1150) runs a second O(K×M_cell) accumulation pass immediately after the errRp pass (lines 978–1013) that already walks the same data. The two passes can be fused into one `compute_cell_metrics` call at the cost of refactoring the `errRp` accumulation out of the first pass.

### Root cause

The errRp pass (lines 975–1013) accumulates `errRp` and `marg_kl`. The extra-metrics pass (lines 1122–1150) accumulates `mean_err_sum`, `kl_max`, `chi2_total`, `grake_norm` — and replicates the same `S_lin` accumulation loop with identical O(K×M_cell) cost. `lbw::compute_cell_metrics` (calib_dispatch.hpp line 149) computes all six metrics in a single O(K×M_cell) pass and returns them in a `CellMetrics` struct.

### Verified `compute_cell_metrics` field mapping (calib_dispatch.hpp lines 149–189)

Read from source — not inferred:

| ieppa local variable | `CellMetrics` field | Equivalence |
|---|---|---|
| `errRp` (global max of per-margin max errors) | `cm.errRp` | **Bit-identical**: same formula `max over k,j of |bucket[j]/W - targets[k][j]|` with same W=W_total (passed as parameter, matches ieppa's `W_total`) |
| `mean_err_sum / st.K` (arithmetic mean of per-margin max errors) | `cm.mean_err` | `cm.mean_err = mean_sum / K` where `mean_sum += max_k` per margin — **identical formula**. Therefore `mean_err_sum = cm.mean_err * st.K` recovers `mean_sum` exactly. |
| `kl_max` (max over k of per-margin KL) | `cm.kl` | Same: `if (kl_k > m.kl) m.kl = kl_k` |
| `chi2_total` | `cm.chi2` | Same accumulation |
| `grake_norm` | `cm.grake_norm` | Same accumulation |
| `marg_kl` (sum over k of per-margin KL) | NOT in CellMetrics | Must retain separate accumulation — see "retained unchanged" note |

`CellMetrics` fields: `errRp`, `mean_err`, `kl`, `chi2`, `grake_norm` (`l1` is not populated by `compute_cell_metrics`; it remains computed inline).

### Before (schematic, lines 975–1150)

```cpp
// Pass 1: errRp + marg_kl (always runs)
double errRp = 0.0, marg_kl = 0.0;
if (use_linear) {
    for (int k = 0; k < st.K; k++) { ... S_lin accum ... errRp max ... marg_kl accum ... }
} else { ... }

// Pass 2: extra metrics (guarded by need_extra_metrics)
if (need_extra_metrics && W_total > 0.0) {
    for (int k = 0; k < st.K; k++) { ... S_lin accum again ... }
}
double mean_err = mean_err_sum / st.K;
```

### After

When `need_extra_metrics` is true AND `use_linear` is true AND `W_total > 0.0`, replace both passes with a single `compute_cell_metrics` call. The `marg_kl` value must still be computed — it is produced in the errRp pass and is NOT returned by `compute_cell_metrics`, so it must be retained as a separate single-pass accumulation for the `marg_kl` field only.

**Concrete change:** In the `if (need_extra_metrics && W_total > 0.0)` block (line 1122), replace the inner K-loop with:

```cpp
if (need_extra_metrics && W_total > 0.0) {
    // 773f.4: fuse extra-metrics sweep with errRp via compute_cell_metrics.
    // S_lin scratch is sized to max_cat (invariant from outer scope) — bucket
    // alias is valid since compute_cell_metrics uses its own bucket parameter.
    // W_total already computed above.
    const lbw::CellMetrics cm = lbw::compute_cell_metrics(st, ct, X, W_total, S_lin);
    // errRp already set above; cm.errRp is identical — use it to update res.
    errRp        = cm.errRp;     // keep consistent with cm
    mean_err_sum = cm.mean_err * static_cast<double>(st.K);  // back out for downstream mean_err calc
    kl_max       = cm.kl;
    chi2_total   = cm.chi2;
    grake_norm   = cm.grake_norm;
    // BEST-ITER UPDATE block 2 below uses mean_err_sum / st.K = cm.mean_err
}
```

> **IMPORTANT:** `compute_cell_metrics` requires a `bucket` scratch of size `>= max(cat_counts)`. In ieppa, `S_lin` is already sized to `max_cat` (verified by reading the outer declaration). Pass `S_lin` as the `bucket` argument. After the call, `S_lin` contents are overwritten but it is a scratch buffer — no downstream reader before next fill.

Replace the `mean_err` local with:
```cpp
double mean_err = (st.K > 0) ? (mean_err_sum / static_cast<double>(st.K)) : 0.0;
// After 773f.4: mean_err_sum = cm.mean_err * st.K when need_extra_metrics, so this is cm.mean_err.
```

The `marg_kl` computation in the first pass (lines 991–997, 1008–1011) is **retained unchanged** since `compute_cell_metrics` does not compute marginal KL.

### Unchanged components

The entire errRp sweep (lines 975–1013), the BEST-ITER update blocks, the `marg_kl` accumulation, the SOR omega block (lines 1047–1078), the `log-space` path. All untouched.

### Compile & test

```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean .
Rscript -e 'testthat::test_dir("tests/testthat")'
```

### Numerical identity — errRp bit-identity proof

ieppa errRp pass (linear path, lines 978–998):
```
for k in 0..K-1:
  fill S_lin[0..nj-1] = 0
  for c in 0..M_cell-1: S_lin[g_per_cell[k][c]] += X[c]
  for j in 0..nj-1: errRp = max(errRp, |S_lin[j]/W_total - targets[k][j]|)
```

`compute_cell_metrics` (calib_dispatch.hpp lines 155–175):
```
for k in 0..K-1:
  fill bucket[0..nj-1] = 0
  for c in 0..M_cell-1: bucket[g_per_cell[k][c]] += X[c]
  for j in 0..nj-1: m.errRp = max(m.errRp, |bucket[j]/W - targets[k][j]|)
```

The two loops are structurally identical with `W = W_total` (ieppa passes `W_total` as the
`W` parameter). Output is bit-identical for the linear path. The `errRp = cm.errRp`
assignment in Task 3 is therefore a no-op on the value — it just keeps the variable
consistent with `cm`.

### Git commit message

```
perf(ieppa): fuse extra-metrics sweep into single compute_cell_metrics call

Eliminates one redundant O(K*M_cell) accumulation pass per convergence
check when need_extra_metrics is true by delegating to
compute_cell_metrics. marg_kl retained as standalone accumulation.
No algorithmic change; output identity verified to 1e-10.
```

---

## Task 4 (773f.6): Guard raking `F_eval` full metrics by check-interval

**File:** `src/raking.cpp`  
**Impact:** HIGH — `F_eval` (line 348) calls `lbw::compute_cell_metrics` unconditionally on every iteration. The non-SRAA path (line 412) already gates on `kErrCheckInterval` for result reporting but calls `F_eval` every iteration. `F_eval` always computes all 6 metrics — full O(K×M_cell) sweep — even when the result is discarded.

### Root cause

`F_eval` (line 348):
```cpp
last_F_metrics = lbw::compute_cell_metrics(st, ct, Xv, static_cast<double>(st.n), bucket);
return last_F_metrics.errRp;
```

In the non-SRAA loop (line 412):
```cpp
for (int iter = 1; iter <= st.inner_max_iter; iter++) {
    double errRp = F_eval(X);
    if (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
        // use metrics
    }
}
```

On non-check iterations, `F_eval` performs the full `compute_cell_metrics` sweep and stores to `last_F_metrics`, but none of those results are read.

### Fast-path: errRp-only on non-check iterations

The `F_eval` lambda is defined at **line 270** and captures outer scope by `[&]`. Therefore
`f_eval_full_metrics` **must be declared at line 269** (the line immediately before
`auto F_eval = [&]`), not before the non-SRAA loop at line 411. The lambda captures it by
reference through `[&]`, so setting it in the non-SRAA loop still controls its value at
call time — but the declaration must precede the capture site at line 270.

**Step 1:** Declare `bool f_eval_full_metrics = true;` at **line 269** (just before `auto F_eval = [&]`).
Insert immediately before the `// last_F_metrics:` comment block (~line 261):

```cpp
// 773f.6: controls whether F_eval computes full metrics or errRp-only fast path.
// Declared here (before lambda) so [&] capture includes it.
// SRAA path: always true (sraa_step needs last_F_metrics).
// Non-SRAA path: set per-iteration before F_eval call.
bool f_eval_full_metrics = true;  // default true; non-SRAA loop overrides per-iter
```

**Step 2:** Convert the last two lines of `F_eval` (the `compute_cell_metrics` call and
return, currently near line 348) from:

```cpp
last_F_metrics = lbw::compute_cell_metrics(st, ct, Xv, static_cast<double>(st.n), bucket);
return last_F_metrics.errRp;
```

to:

```cpp
// 773f.6: on non-check iters, compute only errRp (fast path).
if (f_eval_full_metrics) {
    last_F_metrics = lbw::compute_cell_metrics(st, ct, Xv, static_cast<double>(st.n), bucket);
    return last_F_metrics.errRp;
}
// Fast path: compute errRp only — same O(K*M_cell) sweep, no chi2/kl/grake.
double errRp_fast = 0.0;
{
    double W_fast = static_cast<double>(st.n);
    for (int k = 0; k < st.K; k++) {
        const int nj = st.cat_counts[k];
        std::fill(bucket.begin(), bucket.begin() + nj, 0.0);
        for (int c = 0; c < ct.M_cell; c++) {
            int g = ct.g_per_cell[k][c];
            if (g >= 0 && g < nj) bucket[g] += Xv[c];
        }
        for (int j = 0; j < nj; j++) {
            double e = std::fabs(bucket[j] / W_fast - st.targets[k][j]);
            if (e > errRp_fast) errRp_fast = e;
        }
    }
}
return errRp_fast;
```

**Step 3:** In the non-SRAA loop (line 411), set `f_eval_full_metrics` per-iteration
before calling `F_eval`:

```cpp
for (int iter = 1; iter <= st.inner_max_iter; iter++) {
    res.iterations = iter;
    f_eval_full_metrics = (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter);
    double errRp = F_eval(X);
    if (f_eval_full_metrics) {
        res.max_error = errRp;
        // ... rest of check block unchanged ...
    }
}
```

**Step 4 (SRAA path):** `f_eval_full_metrics` is initialized to `true` above, so the SRAA
path (which calls `F_eval` via `lbw::sraa_step`) always runs the full path. No change needed
for SRAA correctness.

### Unchanged components

`F_eval` core logic (water-fill, hyperplane normalization), SRAA path behavior, convergence checks. All metric assignments in the check block remain identical.

### Compile & test

```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean .
Rscript -e 'testthat::test_dir("tests/testthat")'
```

### Identity note

On check iterations, `f_eval_full_metrics=true` → full `compute_cell_metrics` path → identical output. On non-check iterations, the returned `errRp_fast` value is assigned to `errRp` local but never written to `res` (the `if (iter == 1 || ...)` gate is false) — so result struct is unaffected.

### Git commit message

```
perf(raking): guard F_eval full metrics computation by check-interval

On non-check iterations, F_eval now computes only errRp instead of
all 6 metrics via compute_cell_metrics, reducing O(K*M_cell) metric
overhead by (kErrCheckInterval-1)/kErrCheckInterval per solve.
SRAA path unchanged. Output identity verified to 1e-10.
```

---

## Task 5 (773f.3): Hoist `per_margin_err` outside inner iter loop in `ieppa.cpp`

**File:** `src/ieppa.cpp`  
**Impact:** MEDIUM — `per_margin_err` is allocated twice per iteration when `use_greedy` is true: once on the linear path (line 667) and once on the log path (line 776). Both are inside the inner loop at line 442.

### Root cause

Linear-path greedy block (line 666–689):
```cpp
if (use_greedy) {
    std::vector<double> per_margin_err(st.K);
    ...
}
```

Log-path greedy block (line 776):
```cpp
std::vector<double> per_margin_err(st.K);
```

### After

`use_greedy` is declared at line 657 (`const bool use_greedy = ...`), which is INSIDE the
inner loop body (the loop starts at ~line 442). Therefore `use_greedy` is not in scope at
line 442. Use an unconditional declaration instead — `st.K` is small and the alloc cost of
an unconditionally-sized vector is negligible:

Add `std::vector<double> per_margin_err(st.K);` on the line immediately before
`for (int iter_in_lvl = 1; iter_in_lvl <= budget_lvl; iter_in_lvl++) {` (line ~442):

```cpp
// 773f.3: hoist per_margin_err outside inner iter loop — eliminates K-sized
// alloc on every iteration when use_greedy is true.
// Declared unconditionally (use_greedy not in scope here; st.K is small); written only when use_greedy.
std::vector<double> per_margin_err(st.K);
for (int iter_in_lvl = 1; iter_in_lvl <= budget_lvl; iter_in_lvl++) {
```

Remove the `std::vector<double> per_margin_err(st.K)` declarations at lines 667 and 776.
The vector is already declared and sized correctly; the `for` loops that follow are unchanged.

### Unchanged components

Everything. The greedy selection logic and loop body are identical.

### Compile & test

```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean .
Rscript -e 'testthat::test_dir("tests/testthat")'
```

### Git commit message

```
perf(ieppa): hoist per_margin_err outside inner iter loop

Eliminates one K-sized heap alloc per iteration on greedy scheduler
path (both linear and log branches) by declaring the vector once
before the loop. No algorithmic change; output identity verified to 1e-10.
```

---

## Task 6 (773f.5): Eliminate sinkhorn double O(M_cell) copy on non-projection path

**File:** `src/sinkhorn.cpp`  
**Impact:** MEDIUM — the `needs_projection=false` branch (line 200–205) executes `X_proj = X` (a full M_cell-element copy), and then line 214 executes `for (int c = 0; c < ct.M_cell; c++) X[c] = X_proj[c]` (another M_cell copy). On the non-projection path, `X_proj` is never read between these two copies, so the first copy is pure overhead.

### Root cause

```cpp
// Line 200–205
} else {
    X_proj = X;   // ← copy 1: O(M_cell), unnecessary
    // B10: ...
}
if (needs_projection) {
    // ... Dykstra update ...
}
for (int c = 0; c < ct.M_cell; c++) X[c] = X_proj[c];  // ← copy 2
```

When `needs_projection=false`, `X_proj = X` and then `X[c] = X_proj[c]` is a round-trip copy with no net effect (X → X_proj → X). Both copies are O(M_cell).

### Verification

`X_proj` is not read between lines 201 and 214 when `needs_projection=false`. The Dykstra `a[]` update block (lines 206–213) is fully guarded by `if (needs_projection)`. Confirmed by reading lines 185–214.

### After

Move the `X[c] = X_proj[c]` loop inside the `if (needs_projection)` block; remove the `X_proj = X` assignment from the else branch:

```cpp
bool needs_projection = false;
for (int c = 0; c < ct.M_cell; c++) {
    if (X[c] < L_cell[c] - 1e-12 || X[c] > U_cell[c] + 1e-12) {
        needs_projection = true;
        break;
    }
}
double mu = 0.0;
if (needs_projection) {
    if (!bisect_capacity_fast(X, exp_a.data(), L_cell, U_cell, ct.M_cell, target_mass, mu, X_proj)) {
        res.status = RK_ERR_INFEAS;
        break;
    }
    // Dykstra correction
    for (int c = 0; c < ct.M_cell; c++) {
        if (X[c] > 1e-300 && X_proj[c] > 1e-300)
            a[c] += std::log(X[c]) - std::log(X_proj[c]);
        a[c] = std::clamp(a[c], -kAmax, kAmax);
    }
    lbw::bulk_scaled_exp(1.0, a.data(), exp_a.data(), ct.M_cell);
    // 773f.5: apply projection result to X only on the projection path.
    for (int c = 0; c < ct.M_cell; c++) X[c] = X_proj[c];
}
// else: X unchanged — no copy needed.
// B10: do NOT zero a[] on the non-projection path.
```

Remove the original lines:
- `} else { X_proj = X; ... }` (the else branch with the unnecessary copy)
- `for (int c = 0; c < ct.M_cell; c++) X[c] = X_proj[c];` (the copy after the if/else)

### Unchanged components

`bisect_capacity_fast`, `exp_a`, the Dykstra `a[]` accumulation, convergence check, all metric reporting. Only the data-movement pattern changes.

### Compile & test

```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean .
Rscript -e 'testthat::test_dir("tests/testthat")'
```

### Git commit message

```
perf(sinkhorn): eliminate redundant O(M_cell) copy on non-projection path

When needs_projection=false, X_proj=X followed by X[c]=X_proj[c] was
a no-op round-trip copy. Move the X[c]=X_proj[c] loop inside the
needs_projection branch, eliminating 2*M_cell element copies per
non-binding iteration. No algorithmic change; output identity to 1e-10.
```

---

## Task 7 (773f.7): Fuse SOR omega adaptation `errRp` sweep with convergence sweep in `ieppa.cpp`

**File:** `src/ieppa.cpp`  
**Impact:** MEDIUM — when `sor_active && sor_auto && iter >= sor_burnin_v`, the SOR omega block (lines 1047–1078) runs a per-margin errRp sweep (O(K×M_cell)) to compute `errRp_k`. This sweep happens inside the same convergence-check block that already computed `errRp` by iterating over all K margins in the errRp pass (lines 975–1013).

### Root cause

The errRp pass (lines 975–1013) loops over `k=0..K-1` and accumulates `S_lin[j]` per margin, then takes max over j. This already computes the per-margin max error implicitly — but only stores the global max `errRp`. The SOR block (line 1048–1078) re-accumulates `S_lin` per margin to get `errRp_k`.

### After

During the convergence errRp sweep (lines 975–1013, linear path only), capture per-margin errRp into the existing `sor_prev_errRp` working array (or a local `per_k_errRp` vector pre-declared before the check block). Then the SOR omega block skips its O(K×M_cell) re-accumulation entirely.

**Step 1:** Declare `per_k_errRp_cache` and `per_k_errRp_valid` at the **top of the inner
iteration loop body** — after `res.iterations = iter;` at line ~444, and **outside** the
convergence-check `if` gate at line ~969
(`if (iter == 1 || iter % kErrCheckInterval == 0 || iter_in_lvl == budget_lvl)`).
Placing them inside that `if` block would make them invisible to the SOR block when the
convergence check does not fire.

```cpp
// 773f.7: per-margin errRp captured during the errRp sweep; reused by SOR block.
// Declared at outer iter loop body level — OUTSIDE the convergence-check if gate.
std::vector<double> per_k_errRp_cache(st.K, 0.0);
bool per_k_errRp_valid = false;
```

Reset `per_k_errRp_valid = false;` at the start of each iteration (immediately after the
declaration above, before any path-dependent code). This ensures stale cache from the
previous iteration is never used if the convergence check did not run.

**Step 2:** In the linear-path errRp loop (lines 978–998), capture per-margin max:

```cpp
for (int k = 0; k < st.K; k++) {
    const int nj = st.cat_counts[k];
    const int off = cat_offset[k];
    std::fill(S_lin.begin(), S_lin.begin() + nj, 0.0);
    const int* gk = ct.g_per_cell[k].data();
    for (int c = 0; c < ct.M_cell; c++) {
        int j = gk[c];
        if (j >= 0 && j < nj) S_lin[j] += X[c];
    }
    double errRp_k = 0.0;
    for (int j = 0; j < nj; j++) {
        double e = std::fabs(S_lin[j] / W_total - st.targets[k][j]);
        if (e > errRp) errRp = e;
        if (e > errRp_k) errRp_k = e;     // 773f.7: capture per-margin
        // marg_kl unchanged
        double tkj = st.targets[k][j];
        double skj = (W_total > 0.0) ? S_lin[j] / W_total : 0.0;
        if (tkj > 1e-300 && skj > 1e-300)
            marg_kl += tkj * std::log(tkj / skj);
    }
    per_k_errRp_cache[k] = errRp_k;   // 773f.7
}
per_k_errRp_valid = (use_linear && W_total > 0.0);
```

**Step 3:** In the SOR block (line 1047–1078), replace the inner K-loop with:

```cpp
if (sor_active && sor_auto && iter >= sor_burnin_v) {
    if (W_total > 0.0) {
        for (int k = 0; k < st.K; k++) {
            double errRp_k;
            if (per_k_errRp_valid) {
                // 773f.7: reuse cached per-margin errRp from the convergence sweep above.
                errRp_k = per_k_errRp_cache[k];
            } else {
                // Log path: re-accumulate (unchanged behavior).
                const int nj_k = st.cat_counts[k];
                std::fill(S_lin.begin(), S_lin.begin() + nj_k, 0.0);
                const int* gk_s = ct.g_per_cell[k].data();
                for (int c = 0; c < ct.M_cell; c++) {
                    int j = gk_s[c];
                    if (j >= 0 && j < nj_k) S_lin[j] += X[c];
                }
                errRp_k = 0.0;
                for (int j = 0; j < nj_k; j++) {
                    double e = std::fabs(S_lin[j] / W_total - st.targets[k][j]);
                    if (e > errRp_k) errRp_k = e;
                }
            }
            bool decreasing = (errRp_k < sor_prev_errRp[k]);
            bool sign_flip  = !decreasing && sor_prev_decreasing[k];
            if (sign_flip) {
                sor_omega[k] = std::max(omega_min_v, sor_omega[k] * 0.7);
                sor_n_damped++;
            } else if (decreasing) {
                sor_omega[k] = std::min(1.0, sor_omega[k] * 1.05);
            }
            if (sor_omega[k] < sor_min_omega) sor_min_omega = sor_omega[k];
            sor_prev_decreasing[k] = decreasing;
            sor_prev_errRp[k]      = errRp_k;
        }
    }
}
```

Remove the original inner K-loop body from lines 1049–1077 (replaced above).

### Unchanged components

Log-path SOR behavior (falls through to re-accumulation), omega update math, `sor_prev_errRp`, `sor_prev_decreasing`. All convergence reporting unchanged.

### Compile & test

```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean .
Rscript -e 'testthat::test_dir("tests/testthat")'
```

### Git commit message

```
perf(ieppa): fuse SOR omega errRp sweep with convergence sweep

On the linear path, per-margin errRp_k captured during the convergence
check pass is reused by the SOR omega adaptation block, eliminating
one redundant O(K*M_cell) accumulation per check iteration.
Log path unchanged. Output identity verified to 1e-10.
```

---

## BENCHMARKING

Run after completing Tasks 1, 2, 3, and 4 to measure actual speedup.

### Setup

Ensure the package is installed with the change applied:

```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean .
```

### Benchmark command

```bash
Rscript benchmarks/stepstone_benchmark.R
```

### What to measure (per task)

| Task | Solver/method flag | Expected saving |
|------|--------------------|-----------------|
| 1 (773f.1) | `method = "lbfgsb"` | 2K allocs per solve eliminated; negligible wall-time but visible in allocation profiler |
| 2 (773f.2) | `method = "lbfgsb"` | Up to 400n-allocs per outer iter on logit link; visible wall-time on large n |
| 3 (773f.4) | `method = "ieppa"` with non-MAX_ERR metric | One full O(K×M_cell) pass removed per check iter |
| 4 (773f.6) | `method = "raking"` | `(kErrCheckInterval-1)/kErrCheckInterval` fraction of `compute_cell_metrics` calls removed |

### Before/after benchmark (mandatory per CLAUDE.md §2)

`bench::mark(before=..., after=...)` with identical calls measures the SAME installed binary
and cannot detect changes. The correct procedure is sequential: measure old binary, install
new binary, measure new binary.

**Step 1: With OLD binary installed** (before any source change):

```bash
Rscript -e '
  library(bench); library(leafblower); library(arrow)
  fx <- "tests/testthat/fixtures/stepstone_small.parquet"
  tg <- "tests/testthat/fixtures/stepstone_small_targets.rds"
  data   <- arrow::read_parquet(fx)
  target <- readRDS(tg)

  cat("lbfgsb OLD:", format(bench::mark(
    harvest(data=data, target=target, method="lbfgsb"),
    iterations=20)$median), "\n")
  cat("ieppa OLD:", format(bench::mark(
    harvest(data=data, target=target, method="ieppa",
            convergence=list(metric="MEAN_ERR", absolute=1e-6)),
    iterations=20)$median), "\n")
  cat("raking OLD:", format(bench::mark(
    harvest(data=data, target=target, method="raking"),
    iterations=20)$median), "\n")
'
```

**Step 2: Apply source changes and rebuild:**

```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean .
```

**Step 3: With NEW binary installed** (same script, re-run):

```bash
Rscript -e '
  library(bench); library(leafblower); library(arrow)
  fx <- "tests/testthat/fixtures/stepstone_small.parquet"
  tg <- "tests/testthat/fixtures/stepstone_small_targets.rds"
  data   <- arrow::read_parquet(fx)
  target <- readRDS(tg)

  cat("lbfgsb NEW:", format(bench::mark(
    harvest(data=data, target=target, method="lbfgsb"),
    iterations=20)$median), "\n")
  cat("ieppa NEW:", format(bench::mark(
    harvest(data=data, target=target, method="ieppa",
            convergence=list(metric="MEAN_ERR", absolute=1e-6)),
    iterations=20)$median), "\n")
  cat("raking NEW:", format(bench::mark(
    harvest(data=data, target=target, method="raking"),
    iterations=20)$median), "\n")
'
```

**Step 4:** Report OLD vs NEW median for each method and compute ratio.

> Rationale: `bench::mark(before=X, after=X)` with identical expressions runs both against
> the currently-installed binary — it cannot measure the effect of a source change. The only
> valid comparison is OLD binary median vs NEW binary median measured separately.

### Acceptance criterion

No regression on any method. Tasks 2, 3, 4 expected to show ≥5% wall-time improvement on the stepstone_small fixture with `K ≥ 5` and `M_cell ≥ 10000`. Tasks 1, 5, 6, 7 may show negligible wall-time improvement but zero regression.
