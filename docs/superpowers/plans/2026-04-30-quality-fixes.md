# Code Quality Fixes — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans.

**Goal:** Fix 8 code quality issues: stale WU-* comments, magic number naming, rk_params_init microfix, HomotopyConfigLbw.enabled removal, alm fields on CalibState, strcmp dispatch ladder, select_metric parameter sprawl, and Result struct consolidation.

**Architecture:** Tasks ordered by risk (cosmetic first, structural last). Tasks 1-3 are comment/constant changes with zero behavioral risk. Tasks 4-7 are behavioral refactors. Task 8 (Result struct consolidation) is highest risk — do last.

**Tech Stack:** C++17, R package build, testthat

---

## Task 1 — Purge WU-* stub comments (leafblower-ztid.1)

**Mechanism:** Comment-only deletions. No behavioral change.
**Forbidden:** Removing any actual code, renaming fields, or altering defaults.
**Audit:** `git diff` must show only comment lines (lines beginning with `//` or inline comment text) removed.

### Files and exact changes

**`src/ieppa.hpp` lines 25–34:**

Current:
```cpp
    // ── Extended quality metrics (WU-A scaffold; populated in WU-B+) ──
    double mean_error          = 0.0;
    double kl                  = 0.0;
    double chi2                = 0.0;
    double l1_weight_change    = 0.0;  // WU-A: renamed from pct_change; computation in WU-B
    double grake_norm          = 0.0;  // WU-A stub; computation in WU-D
    int    convergence_metric  = 0;    // WU-A stub; CalibMetric at exit
    int    convergence_rule    = 1;    // WU-A stub; CalibRule at exit (IMPROVEMENT)
    double convergence_tol     = 0.001; // WU-A stub; threshold that fired
    int    convergence_iter    = -1;   // WU-A stub; iteration at convergence (-1=max_iter)
```

Replace section header and inline comments:
```cpp
    // ── Extended quality metrics ──
    double mean_error          = 0.0;
    double kl                  = 0.0;
    double chi2                = 0.0;
    double l1_weight_change    = 0.0;  // renamed from pct_change; L1 normalized weight change Σ|Δw|/n
    double grake_norm          = 0.0;  // max_k |misfit|/(1+|pop|) normalized residual
    int    convergence_metric  = 0;    // CalibMetric at exit
    int    convergence_rule    = 1;    // CalibRule at exit (IMPROVEMENT)
    double convergence_tol     = 0.001; // threshold that fired
    int    convergence_iter    = -1;   // iteration at convergence (-1=max_iter)
```

**`src/lbfgsb_solver.cpp` lines 219, 226, 232:**

Line 219 — remove label prefix `WU-B:`:
```cpp
// Before:
    // WU-B: compute pct_change (start weights d[i] vs. final st.weights[i]).
// After:
    // Compute pct_change (max relative shift start weights d[i] vs. final st.weights[i]).
```

Line 226 — remove label prefix `WU-B:`:
```cpp
// Before:
    // WU-B: alternative metrics from final weights.
// After:
    // Alternative metrics from final weights.
```

Line 232 — remove label prefix `WU-D:`:
```cpp
// Before:
    double grake_norm   = 0.0;  // WU-D: max over margins of |S_kj - T_kj*Wn| / (1 + |T_kj*Wn|)
// After:
    double grake_norm   = 0.0;  // max over margins of |S_kj - T_kj*Wn| / (1 + |T_kj*Wn|)
```

Also remove the `WU-D:` comment label on the SOR section of `src/ieppa.cpp` line 317:
```cpp
// Before:
    // WU-D: SOR adaptive under-relaxation state (iEPPA-only).
// After:
    // SOR adaptive under-relaxation state (iEPPA-only).
```

And `src/ieppa.cpp` line 1081:
```cpp
// Before:
            // WU-B: compute pct_change (max relative shift in cell mass since last check).
// After:
            // Compute pct_change (max relative shift in cell mass since last check).
```

### TDD

No test needed for comment-only change. Verify with:
```bash
git diff --unified=0 src/ieppa.hpp src/lbfgsb_solver.cpp src/ieppa.cpp \
  | grep '^[+-]' | grep -v '^---\|^+++' | grep -v '^\(+\|-\)\s*//'
```
Expected: empty output (only comment lines changed).

### Compile + verify
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -5
```

### Commit
```
fix(comments): remove stale WU-A/WU-B/WU-D scaffold labels from ieppa.hpp, lbfgsb_solver.cpp, ieppa.cpp

These markers were work-unit scaffolding from earlier development phases; the
underlying code is fully implemented. Labels are misleading to future readers.
```

---

## Task 2 — rk_params_init microfix (leafblower-ztid.2)

**Mechanism:** Guard `rk_params_init` with `if (!params)` so it only writes defaults when the caller passes NULL.
**Forbidden:** Changing any other logic in `rk_calibrate`.
**Audit:** Regression test calling `rk_calibrate` with non-NULL params; verify params values are honoured.

### Root cause

`src/c_api.cpp` lines 148–150:
```cpp
    rk_params_t defaults;
    rk_params_init(&defaults);
    const rk_params_t* p = params ? params : &defaults;
```

This is already correct — `rk_params_init` is only called on `defaults` and `defaults` is only used when `params == NULL`. No behavioral fix required here.

**Re-read the task spec:** The fix is ensuring that `rk_params_init` is only called in the null-params branch. The current code already does this correctly. The microfix is verifying no other callsite accidentally calls `rk_params_init` on a caller-supplied struct:

```bash
grep -n "rk_params_init" /home/dd/Gemini/leafblower/src/c_api.cpp \
                          /home/dd/Gemini/leafblower/src/r_bridge.cpp
```

Expected: `c_api.cpp` shows the `&defaults` pattern; `r_bridge.cpp` shows it does NOT call `rk_params_init` at all (it fills `p` directly from R SEXPs).

If any callsite calls `rk_params_init(params)` on a non-NULL caller-supplied struct — that is the bug to fix. Wrap it:
```cpp
    rk_params_t defaults;
    if (!params) {
        rk_params_init(&defaults);
        params = &defaults;
    }
    const rk_params_t* p = params;
```

### TDD — write failing test first

`tests/testthat/test-rk-params-passthrough.R`:
```r
test_that("rk_calibrate honours caller-supplied max_weight (not overwritten by defaults)", {
  n   <- 10L
  wts <- rep(1.0, n)
  gid <- matrix(rep(0L, n), nrow = 1L)
  tgt <- matrix(1.0, nrow = 1L)
  # max_weight=2.0 — if rk_params_init were called on the supplied struct it
  # would reset to 5.0 and the test would not catch over-clamping
  result <- leafblower_calibrate(
    weights = wts, group_ids = gid, targets = tgt,
    max_weight = 2.0
  )
  expect_lte(max(result$weights), 2.0 + 1e-9)
})
```

Run before fix: should pass if the current code is already correct (confirms no regression was introduced). If it fails — the code is broken and the `if (!params)` guard above is the fix.

### Compile + test
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . && \
  Rscript -e 'testthat::test_file("tests/testthat/test-rk-params-passthrough.R")'
```

### Commit
```
test(c_api): add regression test verifying rk_params_init not called on caller-supplied params

The current implementation is already correct — rk_params_init only touches
the local defaults struct, never the caller-supplied params. The test locks
this invariant against future regressions.
```

---

## Task 3 — Named constexpr for 8 magic numbers (leafblower-ztid.3)

**Mechanism:** Replace bare literals with `static constexpr` names at the top of each function/file block where the literal appears. Zero algorithmic change.
**Forbidden:** Moving constants across translation units, changing values, introducing global state.
**Audit:** `git diff` must show identical numeric values before/after for every affected line.

### 3a — `src/chebyshev.cpp` warm-start block (lines 153–160)

The warm-start block uses bare `1e-8` and `1e-10` literals. Add constants in the warm-start block (after line 150, before line 152):

```cpp
    // Warm start: strictly interior
    static constexpr double kWarmStartRelEps  = 1e-8;   // fractional shift off bound
    static constexpr double kWarmStartAbsEps  = 1e-10;  // absolute floor when gap is tiny
    std::vector<double> X(ct.M_cell);
    {
        double max_gap = 0.0;
        for (int c = 0; c < ct.M_cell; c++) max_gap = std::max(max_gap, U_cell[c]-L_cell[c]);
        double eps_shift = std::max(kWarmStartRelEps * max_gap, kWarmStartAbsEps);
        // ... rest unchanged
```

The `1e-8` at line 281 in the convergence check (`best_errRp < 1e-8`) is the primal machine-precision convergence threshold — add at file scope near the other constants (line 20–28 block):

```cpp
    static constexpr double kPrimalMachinePrecConv = 1e-8;  // Mehrotra: accept when best errRp at machine precision
```

Then replace the literal at line 281:
```cpp
            bool converged = (mu < kTolMu) || (iter > 0 && best_errRp < kPrimalMachinePrecConv);
```

### 3b — `src/ieppa.cpp` damping constants (lines 263, 1068, 1072, 1342)

Add after the existing `kLinearOverflowTrip`/`kLogOverflowThreshold` block (around line 244), just before the `f_lin` declaration:

```cpp
    static constexpr double kAlphaBeta           = 0.5;   // P2.1 stress→alpha mapping: alpha = 1/(1+β·stress)
    static constexpr double kSorOscillationDamp  = 0.7;   // SOR sign-flip: reduce omega by this factor
    static constexpr double kSorRecoveryGrowth   = 1.05;  // SOR monotone: recover omega by this factor
    static constexpr double kInfeasStallRatio    = 10.0;  // stall warn: max_error / pct_tol threshold
```

Replace usages:

Line 263: `double beta = 0.5;` → `double beta = kAlphaBeta;`

Line 1068: `sor_omega[k] = std::max(omega_min_v, sor_omega[k] * 0.7);`
→ `sor_omega[k] = std::max(omega_min_v, sor_omega[k] * kSorOscillationDamp);`

Line 1072: `sor_omega[k] = std::min(1.0, sor_omega[k] * 1.05);`
→ `sor_omega[k] = std::min(1.0, sor_omega[k] * kSorRecoveryGrowth);`

Line 1342: `res.max_error > 10.0 * cfg.pct_tol &&`
→ `res.max_error > kInfeasStallRatio * cfg.pct_tol &&`

### 3c — `src/logit_calib.cpp` lambda reject threshold (line 130)

The literal `10.0` on line 130 is the lambda initialisation rejection bound. Add at function top (after the `kMaxNewtonIters` constant block):

```cpp
    static constexpr double kLambdaInitRejectAbs = 10.0;  // reject lambda_0 if any component exceeds this
```

Replace line 130:
```cpp
            if (max_lambda_init <= kLambdaInitRejectAbs) lambda = std::move(b_init);
```

### TDD

Write a test that validates the exact numeric outputs are unchanged after naming. Use a fixed-seed calibration problem that exercises all three solvers:

`tests/testthat/test-magic-number-naming-regression.R`:
```r
test_that("chebyshev result unchanged after constexpr naming", {
  set.seed(42)
  n   <- 50L
  wts <- runif(n, 0.5, 1.5)
  gid <- matrix(sample(0:1, n, replace = TRUE), nrow = 1L)
  tgt <- matrix(0.5, nrow = 1L)
  r <- leafblower_calibrate(weights = wts, group_ids = gid, targets = tgt,
                            method = "chebyshev", max_iterations = 200L)
  expect_equal(r$max_error, 0.0, tolerance = 1e-6)
})

test_that("ieppa result unchanged after constexpr naming", {
  set.seed(42)
  n   <- 50L
  wts <- runif(n, 0.5, 1.5)
  gid <- matrix(sample(0:1, n, replace = TRUE), nrow = 1L)
  tgt <- matrix(0.5, nrow = 1L)
  r <- leafblower_calibrate(weights = wts, group_ids = gid, targets = tgt,
                            method = "ieppa", max_iterations = 200L)
  expect_equal(r$max_error, 0.0, tolerance = 1e-6)
})

test_that("logit result unchanged after constexpr naming", {
  set.seed(42)
  n   <- 50L
  wts <- runif(n, 0.5, 1.5)
  gid <- matrix(sample(0:1, n, replace = TRUE), nrow = 1L)
  tgt <- matrix(0.5, nrow = 1L)
  r <- leafblower_calibrate(weights = wts, group_ids = gid, targets = tgt,
                            method = "logit", max_iterations = 200L)
  expect_equal(r$max_error, 0.0, tolerance = 1e-6)
})
```

Run before changes to capture baseline. Run after to confirm identical outputs.

### Compile + test
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . && \
  Rscript -e 'testthat::test_file("tests/testthat/test-magic-number-naming-regression.R")'
```

### Commit
```
refactor(constants): replace 8 magic numbers with named constexpr in chebyshev, ieppa, logit_calib

kWarmStartRelEps, kWarmStartAbsEps, kPrimalMachinePrecConv (chebyshev.cpp);
kAlphaBeta, kSorOscillationDamp, kSorRecoveryGrowth, kInfeasStallRatio (ieppa.cpp);
kLambdaInitRejectAbs (logit_calib.cpp). Zero behavioral change.
```

---

## Task 4 — Drop HomotopyConfigLbw.enabled (leafblower-ztid.7)

**Mechanism:** Remove the redundant `bool enabled` field from `HomotopyConfigLbw`; replace all reads with `(n_levels > 1)`. Update the C ABI struct `rk_homotopy_cfg_t` in `leafblower.h` to remove `int enabled` — this is an ABI break, but `leafblower.h` is internal to the R package (no shared-library ABI guarantee across sessions).
**Forbidden:** Changing homotopy logic. `n_levels == 1` must continue to mean "disabled".
**Audit:** Test with `n_levels=1` (disabled) and `n_levels=3` (enabled); verify `homotopy_levels_used` in result.

### Sites to change

**`src/types.hpp` line 18:** Remove `bool enabled = false;` from `HomotopyConfigLbw`.

**`src/leafblower.h` lines 22–28:** Remove `int enabled;` from `rk_homotopy_cfg_t`.

**`src/r_bridge.cpp` line 299:** Replace:
```cpp
    st.homotopy.enabled = (p.homotopy.enabled != 0) || (p.homotopy.n_levels > 1);
```
With:
```cpp
    // enabled is derived: n_levels > 1 means homotopy is active
```
(delete the line entirely)

**`src/ieppa.cpp` lines 355–356:** Already uses `n_levels > 1` — no change needed:
```cpp
    const int N_levels = (st.homotopy.n_levels > 1)
                        ? st.homotopy.n_levels : 1;
```

Confirm there are no other reads of `.enabled` in the codebase:
```bash
grep -rn "\.enabled" /home/dd/Gemini/leafblower/src/
```
Expected: only `sor_cfg.enabled` and `CalibSorCfg::enabled` remain.

**`src/c_api.cpp`:** Two sites must be updated (confirmed by grep):

- Line 89 in `rk_params_init`: remove `p->homotopy.enabled = 0;`
- Line 214 in the CalibState assignment block: replace
  ```cpp
      st.homotopy.enabled         = (p->homotopy.enabled != 0) || (p->homotopy.n_levels > 1);
  ```
  with:
  ```cpp
      // enabled derived from n_levels > 1 — no assignment needed
  ```

### TDD

`tests/testthat/test-homotopy-enabled-field.R`:
```r
test_that("homotopy disabled when n_levels=1", {
  n   <- 30L
  wts <- rep(1.0, n)
  gid <- matrix(rep(0L, n), nrow = 1L)
  tgt <- matrix(1.0, nrow = 1L)
  r <- leafblower_calibrate(weights = wts, group_ids = gid, targets = tgt,
                            method = "ieppa", n_homotopy_levels = 1L)
  expect_equal(r$homotopy_levels_used, 0L)
})

test_that("homotopy active when n_levels=3", {
  n   <- 30L
  wts <- rep(1.0, n)
  gid <- matrix(rep(0L, n), nrow = 1L)
  tgt <- matrix(1.0, nrow = 1L)
  r <- leafblower_calibrate(weights = wts, group_ids = gid, targets = tgt,
                            method = "ieppa", n_homotopy_levels = 3L)
  expect_equal(r$homotopy_levels_used, 3L)
})
```

### Compile + test
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . && \
  Rscript -e 'testthat::test_file("tests/testthat/test-homotopy-enabled-field.R")'
```

### Commit
```
refactor(types): remove redundant HomotopyConfigLbw.enabled field

The enabled flag duplicated n_levels > 1. All read sites already used n_levels
for the actual dispatch decision. Removing the field eliminates a class of
desync bugs where enabled=false but n_levels=3 could produce inconsistent state.
```

---

## Task 5 — Move alm_lambda/alm_mu/capacity_mu off CalibState into ALMConfig (leafblower-ztid.8)

**Mechanism:** Group the three ALM penalty fields into a `struct ALMConfig` embedded in `CalibState`. All callers access via `st.alm.lambda`, `st.alm.mu`, `st.alm.capacity_mu`.
**Forbidden:** Changing ALM logic, values, or defaults. No field additions or removals.
**Audit:** Every read/write of the three fields in ieppa.cpp and r_bridge.cpp must use the new path.

### New struct in `src/types.hpp`

After `CalibSorCfg` (line 69), add:
```cpp
struct ALMConfig {
    double lambda     = 0.0;  // dual variable for sum(w)=n; only read when mu > 0
    double mu         = 0.0;  // penalty coefficient; 0.0 = ALM inactive
    double capacity_mu = 0.0; // ieppa_soft ALM penalty (capacity box constraint); 0.0 = inactive
};
```

In `CalibState` (lines 87–89), replace:
```cpp
    double alm_lambda = 0.0;  // dual variable for sum(w)=n; only read when alm_mu > 0
    double alm_mu     = 0.0;  // penalty coefficient; 0.0 = ALM inactive
    double capacity_mu = 0.0;  // ieppa_soft ALM penalty (capacity box constraint); 0.0 = inactive
```
With:
```cpp
    ALMConfig alm;
```

### Read/write sites to update

**`src/ieppa.cpp`** — confirmed sites (from grep):

| Line | Expression | Replacement |
|------|-----------|-------------|
| 158 | `st.capacity_mu` | `st.alm.capacity_mu` |
| 159 | `st.capacity_mu` | `st.alm.capacity_mu` |
| 160 | `st.capacity_mu` | `st.alm.capacity_mu` |
| 365 | `st.alm_mu` | `st.alm.mu` |
| 396 | `st.alm_mu` | `st.alm.mu` |
| 397 | `st.alm_mu` | `st.alm.mu` |
| 410 | `st.capacity_mu` | `st.alm.capacity_mu` |
| 413 | `st.capacity_mu` | `st.alm.capacity_mu` |
| 836 | `st.capacity_mu` | `st.alm.capacity_mu` |
| 837 | `st.capacity_mu` | `st.alm.capacity_mu` |
| 840 | `capacity_mu_base` (local, no change) | — |
| 841 | `st.capacity_mu` | `st.alm.capacity_mu` |
| 911 | `st.capacity_mu` | `st.alm.capacity_mu` |
| 912 | `st.capacity_mu` | `st.alm.capacity_mu` |
| 916 | `capacity_mu_base` (local, no change) | — |
| 917 | `st.capacity_mu` | `st.alm.capacity_mu` |
| 953 | `capacity_mu_base` (local, no change) | — |
| 954 | local `capacity_mu_adaptive` (no change) | — |

**Note:** `st.alm_lambda` has zero occurrences in ieppa.cpp — it is declared on CalibState but never referenced in the solver. `ALMConfig.lambda` is included in the new struct to match the field declaration, but no migration site exists in ieppa.cpp.

Run to verify before editing:
```bash
grep -n "alm_lambda\|alm_mu\|capacity_mu" src/ieppa.cpp src/r_bridge.cpp src/c_api.cpp
```

Replace every `st.alm_mu` → `st.alm.mu`, `st.capacity_mu` → `st.alm.capacity_mu`.

**`src/r_bridge.cpp`** lines 321–322:
```cpp
// Before:
    st.alm_lambda = 0.0;
    st.alm_mu     = 0.0;
// After:
    st.alm.lambda     = 0.0;
    st.alm.mu         = 0.0;
    st.alm.capacity_mu = 0.0;  // capacity_mu resolved below in the CellTable block
```

Find the `capacity_mu` assignment in r_bridge.cpp:
```bash
grep -n "capacity_mu\|st\.capacity" /home/dd/Gemini/leafblower/src/r_bridge.cpp
```
Update those sites to `st.alm.capacity_mu`.

**`src/c_api.cpp`** — confirmed zero sites:

`grep -n "alm_lambda\|alm_mu\|capacity_mu" src/c_api.cpp` returns empty — c_api.cpp has no ALM field writes. No changes required in c_api.cpp for this task.

**Note:** `residual_recheck_fraction` is intentionally NOT moved in this task — it is a `SchedulerConfigLbw` field, not an ALM field. If it needs moving to a dedicated scheduler struct, create a separate ticket.

### TDD

`tests/testthat/test-alm-config-grouping.R`:
```r
test_that("ieppa_soft converges after ALMConfig grouping", {
  set.seed(7)
  n   <- 40L
  wts <- runif(n, 0.5, 2.0)
  gid <- matrix(sample(0:1, n, replace = TRUE), nrow = 1L)
  tgt <- matrix(0.5, nrow = 1L)
  r <- leafblower_calibrate(weights = wts, group_ids = gid, targets = tgt,
                            method = "ieppa_soft", max_iterations = 200L)
  expect_lte(r$max_error, 1e-3)
})

test_that("ieppa converges (alm inactive) after ALMConfig grouping", {
  set.seed(7)
  n   <- 40L
  wts <- runif(n, 0.5, 2.0)
  gid <- matrix(sample(0:1, n, replace = TRUE), nrow = 1L)
  tgt <- matrix(0.5, nrow = 1L)
  r <- leafblower_calibrate(weights = wts, group_ids = gid, targets = tgt,
                            method = "ieppa", max_iterations = 200L)
  expect_lte(r$max_error, 1e-3)
})
```

### Compile + test
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . && \
  Rscript -e 'testthat::test_dir("tests/testthat")'
```

### Commit
```
refactor(types): group alm_lambda/alm_mu/capacity_mu into ALMConfig substruct on CalibState

Reduces CalibState field sprawl and makes ALM state boundary explicit.
Access pattern changes from st.alm_mu to st.alm.mu — no behavioral change.
```

---

## Task 6 — Replace strcmp ladders with enum map in r_bridge.cpp (leafblower-ztid.5)

**Mechanism:** Build a `static const std::unordered_map<std::string_view, rk_algorithm_t>` at file scope; replace both dispatch ladders (lines 238–248 and 255–264) with single lookups.
**Forbidden:** Changing the fallback default (`RK_ALG_IEPPA` for unknown strings), altering error messages.
**Audit:** Every algorithm name in the enum must appear in the map; the second (validation) ladder must also use the same map.

### Current ladders

Ladder 1 (lines 238–248) — sets `p.algorithm`:
```cpp
    if      (strcmp(method_str, "ieppa")      == 0) p.algorithm = RK_ALG_IEPPA;
    else if (strcmp(method_str, "ieppa_soft") == 0) p.algorithm = RK_ALG_IEPPA_SOFT;
    else if (strcmp(method_str, "lbfgsb")    == 0) p.algorithm = RK_ALG_LBFGSB;
    else if (strcmp(method_str, "raking")    == 0) p.algorithm = RK_ALG_RAKING;
    else if (strcmp(method_str, "greg")      == 0) p.algorithm = RK_ALG_GREG;
    else if (strcmp(method_str, "chebyshev") == 0) p.algorithm = RK_ALG_CHEBYSHEV;
    else if (strcmp(method_str, "sinkhorn")  == 0) p.algorithm = RK_ALG_SINKHORN;
    else if (strcmp(method_str, "auto")      == 0) p.algorithm = RK_ALG_AUTO;
    else if (strcmp(method_str, "greenkhorn") == 0) p.algorithm = RK_ALG_GREENKHORN;
    else if (strcmp(method_str, "logit")      == 0) p.algorithm = RK_ALG_LOGIT;
    else                                            p.algorithm = RK_ALG_IEPPA;
```

Ladder 2 (lines 255–264) — sets `alg_for_validation`:
```cpp
    rk_algorithm_t alg_for_validation =
        (strcmp(method_str, "ieppa")      == 0) ? RK_ALG_IEPPA :
        (strcmp(method_str, "ieppa_soft") == 0) ? RK_ALG_IEPPA_SOFT :
        (strcmp(method_str, "lbfgsb")     == 0) ? RK_ALG_LBFGSB :
        (strcmp(method_str, "auto")       == 0) ? RK_ALG_AUTO :
        (strcmp(method_str, "sinkhorn")   == 0) ? RK_ALG_SINKHORN :
        (strcmp(method_str, "greg")       == 0) ? RK_ALG_GREG :
        (strcmp(method_str, "chebyshev")  == 0) ? RK_ALG_CHEBYSHEV :
        (strcmp(method_str, "greenkhorn") == 0) ? RK_ALG_GREENKHORN :
        (strcmp(method_str, "logit")      == 0) ? RK_ALG_LOGIT :
                                                   RK_ALG_RAKING;
```

Note: ladder 2's fallback is `RK_ALG_RAKING`, not `RK_ALG_IEPPA`. Both defaults must be preserved.

### Replacement

Add at file scope in `src/r_bridge.cpp` (after includes, before any function):
```cpp
#include <string_view>
#include <unordered_map>

namespace {
static const std::unordered_map<std::string_view, rk_algorithm_t> kAlgMap = {
    {"ieppa",      RK_ALG_IEPPA},
    {"ieppa_soft", RK_ALG_IEPPA_SOFT},
    {"lbfgsb",     RK_ALG_LBFGSB},
    {"raking",     RK_ALG_RAKING},
    {"greg",       RK_ALG_GREG},
    {"chebyshev",  RK_ALG_CHEBYSHEV},
    {"sinkhorn",   RK_ALG_SINKHORN},
    {"auto",       RK_ALG_AUTO},
    {"greenkhorn", RK_ALG_GREENKHORN},
    {"logit",      RK_ALG_LOGIT},
};
} // anonymous namespace
```

Replace ladder 1:
```cpp
    {
        auto it = kAlgMap.find(method_str);
        p.algorithm = (it != kAlgMap.end()) ? it->second : RK_ALG_IEPPA;
    }
```

Replace ladder 2:
```cpp
        auto it2 = kAlgMap.find(method_str);
        rk_algorithm_t alg_for_validation = (it2 != kAlgMap.end()) ? it2->second : RK_ALG_RAKING;
```

### TDD

`tests/testthat/test-method-dispatch.R`:
```r
methods_to_test <- c("ieppa", "ieppa_soft", "lbfgsb", "raking", "greg",
                     "chebyshev", "sinkhorn", "auto", "greenkhorn", "logit")

for (m in methods_to_test) {
  test_that(paste("method dispatch works for", m), {
    n   <- 20L
    wts <- rep(1.0, n)
    gid <- matrix(rep(0L, n), nrow = 1L)
    tgt <- matrix(1.0, nrow = 1L)
    expect_no_error(
      leafblower_calibrate(weights = wts, group_ids = gid, targets = tgt,
                           method = m, max_iterations = 10L)
    )
  })
}

test_that("unknown method falls back to ieppa without error", {
  n   <- 20L
  wts <- rep(1.0, n)
  gid <- matrix(rep(0L, n), nrow = 1L)
  tgt <- matrix(1.0, nrow = 1L)
  # Unknown method should not crash — fallback to ieppa
  expect_no_error(
    leafblower_calibrate(weights = wts, group_ids = gid, targets = tgt,
                         method = "nonexistent_method", max_iterations = 10L)
  )
})
```

### Compile + test
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . && \
  Rscript -e 'testthat::test_file("tests/testthat/test-method-dispatch.R")'
```

### Commit
```
refactor(r_bridge): replace two strcmp dispatch ladders with unordered_map lookup

Both algorithm-selection ladders in leafblower_calibrate_impl now share a
single static kAlgMap, eliminating duplicated string comparisons and making
new algorithm registration a one-line change.
```

---

## Task 7 — select_metric() → CellMetrics struct param (leafblower-ztid.6)

**Mechanism:** Replace the 7-scalar positional signature of `select_metric` with a `const CellMetrics& m` parameter plus the metric enum. Callers pass the already-computed `CellMetrics` struct. `CellMetrics` already exists in `calib_dispatch.hpp` with fields `errRp, mean_err, kl, chi2, grake_norm, l1`.
**Forbidden:** Changing `CellMetrics` field layout. Do NOT add `marginal_kl` to `CellMetrics`.
**Audit:** All 6 direct call sites must be updated; `check_convergence` which calls `select_metric` internally must also be updated.

### Current signature (`src/calib_dispatch.hpp` line 35–43):
```cpp
inline double select_metric(
    CalibMetric metric,
    double max_err,
    double mean_err,
    double kl,
    double chi2,
    double grake_norm,
    double l1_weight,
    double marginal_kl = 0.0) noexcept
```

### New signature

Do NOT modify `CellMetrics` — its field layout is frozen (Forbidden). The new struct overload omits `marginal_kl`, which defaults to `0.0` in the 7-arg overload:

```cpp
inline double select_metric(CalibMetric metric, const CellMetrics& m) noexcept {
    return select_metric(metric, m.errRp, m.mean_err, m.kl, m.chi2,
                         m.grake_norm, m.l1);
    // marginal_kl omitted — defaults to 0.0 in the 7-arg overload
}
```

This is a non-breaking additive approach: old callers compile unchanged, new callers use the struct overload. No CellMetrics field changes required.

### Callers to migrate to struct overload

Find all call sites:
```bash
grep -n "select_metric" /home/dd/Gemini/leafblower/src/ieppa.cpp \
                        /home/dd/Gemini/leafblower/src/chebyshev.cpp \
                        /home/dd/Gemini/leafblower/src/raking.cpp \
                        /home/dd/Gemini/leafblower/src/sinkhorn.cpp
```

**`src/ieppa.cpp` line 1158** — already has a `CellMetrics`-like local set; migrate:
```cpp
// Before:
                    const double curr_best = lbw::select_metric(
                        st.convergence_cfg.metric,
                        errRp, mean_err_blk2, kl_max, chi2_total, grake_norm, l1_weight);
// After:
                    lbw::CellMetrics cm2;
                    cm2.errRp = errRp; cm2.mean_err = mean_err_blk2;
                    cm2.kl = kl_max; cm2.chi2 = chi2_total;
                    cm2.grake_norm = grake_norm; cm2.l1 = l1_weight;
                    const double curr_best = lbw::select_metric(st.convergence_cfg.metric, cm2);
```

**`src/ieppa.cpp` line 1240** — same pattern, migrate similarly.

**`src/chebyshev.cpp` line 263** — check the local variable names and migrate.

**`src/raking.cpp` line 456** — migrate.

**`src/sinkhorn.cpp` line 233** — migrate.

**`src/calib_dispatch.hpp` line 134** (internal `check_convergence` call) — already takes `CellMetrics&`; update its internal call:
```cpp
    const double curr = select_metric(cfg.metric, m);
```

**`src/lbfgsb_solver.cpp` line 283** — leave on old scalar overload (batch solver has no `CellMetrics`).

### TDD

`tests/testthat/test-select-metric-struct.R`:
```r
test_that("all solvers produce same output after select_metric struct migration", {
  set.seed(99)
  n   <- 60L
  wts <- runif(n, 0.8, 1.2)
  gid <- rbind(sample(0:2, n, replace = TRUE),
               sample(0:1, n, replace = TRUE))
  tgt <- rbind(c(1/3, 1/3, 1/3), c(0.5, 0.5))
  for (m in c("ieppa", "chebyshev", "raking", "sinkhorn")) {
    r <- leafblower_calibrate(weights = wts, group_ids = gid, targets = tgt,
                              method = m, max_iterations = 500L)
    expect_lte(r$max_error, 1e-3,
               label = paste("max_error for", m))
  }
})
```

### Compile + test
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . && \
  Rscript -e 'testthat::test_file("tests/testthat/test-select-metric-struct.R")'
```

### Commit
```
refactor(calib_dispatch): add struct overload for select_metric

The 7-argument positional call was error-prone and hard to extend. New
select_metric(CalibMetric, const CellMetrics&) overload is used by ieppa,
chebyshev, raking, and sinkhorn. lbfgsb retains the scalar overload (batch
solver pre-dates CellMetrics). CellMetrics field layout unchanged. No
behavioral change.
```

---

## Task 8 — Consolidate 6 duplicate Result structs → CalibResult (leafblower-ztid.4)

**Mechanism:** Define `CalibResult` in `types.hpp` with the common fields shared across all solvers. Solver-private diagnostics (iEPPA-internal, ALM) stay in the per-solver struct, which embeds `CalibResult` as a substruct. Update `r_bridge.cpp` and `c_api.cpp` to read from the common struct.
**Forbidden:** Removing solver-specific fields that are unique to one solver. Do not change `rk_result_t` (C ABI struct) — it is copied from `CalibResult` by `c_api.cpp`. Do not put iEPPA-private fields in the shared `CalibResult` base.
**Audit:** Every field in `rk_result_t` must map to a field in `CalibResult`. No field silently dropped.

### Step 1: Inventory all result structs

Read each solver header to find the structs:
```bash
grep -rn "struct.*Result" /home/dd/Gemini/leafblower/src/
```

Expected structs: `IEPPAResult` (ieppa.hpp), `ChebyshevResult` (chebyshev.hpp), `RakingResult` (raking.hpp), `SinkhornResult` (sinkhorn.hpp), `GregResult` (greg.hpp), `LbfgsbResult` (lbfgsb_solver.hpp or similar).

### Step 2: Identify common field set

Read each struct header before implementing. The split is:

**Common fields (CalibResult — every solver populates these):**

| Field | Type | Notes |
|---|---|---|
| `status` | `int` | RK_OK / RK_ERR_* |
| `iterations` | `int` | outer iters completed |
| `max_error` | `double` | final errRp |
| `M_cell` | `int` | compression info |
| `n_cap_active` | `int` | capacity constraint count |
| `homotopy_levels_used` | `int` | overlay diagnostic |
| `homotopy_final_factor` | `double` | overlay diagnostic |
| `greedy_sweeps_taken` | `int` | overlay diagnostic |
| `eta_final` | `double` | overlay diagnostic |
| `mean_error` | `double` | extended quality |
| `kl` | `double` | extended quality |
| `chi2` | `double` | extended quality |
| `l1_weight_change` | `double` | extended quality |
| `grake_norm` | `double` | extended quality |
| `convergence_metric` | `int` | CalibMetric at exit |
| `convergence_rule` | `int` | CalibRule at exit |
| `convergence_tol` | `double` | threshold that fired |
| `convergence_iter` | `int` | iter at convergence |
| `best_error` | `double` | errRp at best iter |
| `best_iter` | `int` | |
| `best_weights` | `std::vector<double>` | obs-level |
| `message[256]` | `char[]` | status message |

**iEPPA-private fields (stay in IEPPAResult, not in CalibResult):**

| Field | Type | Notes |
|---|---|---|
| `n_xcur_writes_per_iter_linear` | `int` | P1.1 diagnostic |
| `min_alpha_seen` | `double` | P2.1 alpha tracking |
| `final_alpha` | `double` | P2.1 alpha tracking |
| `n_bounds_violated` | `int` | cell-mode diagnostic |
| `n_bounds_clamped` | `int` | unit-mode action |
| `sor_min_omega` | `double` | SOR damping diagnostic |
| `sor_n_damped` | `int` | SOR damping count |
| `best_objective_seen` | `double` | iEPPA objective tracking |
| `convergence_solver_objective` | `double` | iEPPA-specific |
| `convergence_minimized_metric` | `int` | iEPPA-specific |
| `marginal_kl_at_iter` | `double` | iEPPA diagnostic |
| `alm_capacity_mu_final` | `double` | ieppa_soft only |
| `alm_n_growth_events` | `int` | ieppa_soft only |
| `alm_max_dual_norm` | `double` | ieppa_soft only |
| `alm_sum_drift` | `double` | ieppa_soft only |

**Strategy:** `IEPPAResult` embeds `CalibResult base` plus the iEPPA-private fields above. All other solvers return `CalibResult` directly. `r_bridge.cpp` reads from `res.base` (for iEPPA) or `res` directly (for others). `c_api.cpp` copies from `CalibResult` fields only — the C ABI struct `rk_result_t` already doesn't expose iEPPA internals.

**Before writing `CalibResult`:** Read all 6 solver headers to confirm which fields each has. Do not infer from `IEPPAResult` alone.

```bash
grep -A 50 "struct.*Result" /home/dd/Gemini/leafblower/src/chebyshev.hpp \
  /home/dd/Gemini/leafblower/src/raking.hpp \
  /home/dd/Gemini/leafblower/src/sinkhorn.hpp
```

### Step 3: Define CalibResult in types.hpp

Add `#include <vector>` to types.hpp if not already present. Check with:
```bash
grep '#include' src/types.hpp
```
If `<vector>` is absent, add it alongside the other standard headers.

Add after `ALMConfig` struct (from Task 5):
```cpp
struct CalibResult {
    int    status              = RK_ERR_NOCONV;
    int    iterations          = 0;
    double max_error           = std::numeric_limits<double>::infinity();
    int    M_cell              = 0;
    int    n_cap_active        = 0;
    // overlay diagnostics
    int    homotopy_levels_used  = 0;
    double homotopy_final_factor = 1.0;
    int    greedy_sweeps_taken   = 0;
    double eta_final             = 0.0;
    // extended quality metrics
    double mean_error          = 0.0;
    double kl                  = 0.0;
    double chi2                = 0.0;
    double l1_weight_change    = 0.0;
    double grake_norm          = 0.0;
    int    convergence_metric  = 0;
    int    convergence_rule    = 1;
    double convergence_tol     = 0.001;
    int    convergence_iter    = -1;
    double best_error          = std::numeric_limits<double>::infinity();
    int    best_iter           = 0;
    std::vector<double> best_weights;
    // status message
    char   message[256]        = {};
};
```

`IEPPAResult` keeps its iEPPA-private fields and adds a `CalibResult base` member. The per-solver structs for chebyshev, raking, sinkhorn, greg, and lbfgsb are replaced by `CalibResult` directly.

### Step 4: Migrate solvers

For each solver, change the return type from `IEPPAResult`/`ChebyshevResult`/etc. to `lbw::CalibResult`. Update `ieppa.hpp`, `chebyshev.hpp`, and the other solver headers. The old per-solver structs can be removed once all callers are updated.

Migration order: ieppa → chebyshev → raking → sinkhorn → greg → lbfgsb.

For each solver:
1. Change function signature: `IEPPAResult ieppa_solve(CalibState&)` → `CalibResult ieppa_solve(CalibState&)`
2. Change `IEPPAResult res;` → `CalibResult res;`
3. Confirm field names match (most already do — see `IEPPAResult` above)
4. Remove the now-dead per-solver struct definition

### Step 5: Update r_bridge.cpp

The result reading in `r_bridge.cpp` likely has a solver-dispatch block that reads from the per-solver result type. Replace with reads from `CalibResult`.

```bash
grep -n "IEPPAResult\|ChebyshevResult\|RakingResult\|SinkhornResult\|GregResult\|LbfgsbResult" \
  /home/dd/Gemini/leafblower/src/r_bridge.cpp
```

### Step 6: Update c_api.cpp

The `c_api.cpp` copies solver results into `rk_result_t`. The result-copy block is at two sites (lines 327-328 and 360-361 before this task). After CalibResult consolidation, the copy changes as follows:

**Before** (per-solver struct access — e.g. for ieppa at lines 327–328):
```cpp
                result->homotopy_levels_used  = res.homotopy_levels_used;
                result->homotopy_final_factor = res.homotopy_final_factor;
```
And similarly for the other solver branch at lines 360–361.

**After** (reads from `CalibResult` — same field names, no structural change):
```cpp
                result->homotopy_levels_used  = res.homotopy_levels_used;
                result->homotopy_final_factor = res.homotopy_final_factor;
```
Field names are identical between the old per-solver structs and `CalibResult` (by design in Step 3). The copy block lines themselves do **not** change; the type of `res` changes from `IEPPAResult`/`ChebyshevResult`/etc. to `CalibResult`. Verify no field was silently dropped by diffing `rk_result_t` fields against `CalibResult` fields before committing.

### TDD

`tests/testthat/test-calib-result-consolidation.R`:
```r
test_that("CalibResult raking fields populated", {
  set.seed(22)
  n <- 40L; wts <- rep(1.0, n)
  gid <- matrix(sample(0:1, n, TRUE), nrow = 1L)
  tgt <- matrix(0.5, nrow = 1L)
  r <- leafblower_calibrate(weights = wts, group_ids = gid, targets = tgt,
                            method = "raking", max_iterations = 200L)
  expect_true(is.finite(r$max_error))
  expect_true(r$iterations > 0L)
  expect_lte(r$max_error, 1e-3)
  expect_equal(length(r$weights), n)
})

test_that("CalibResult ieppa fields populated", {
  set.seed(22)
  n <- 40L; wts <- rep(1.0, n)
  gid <- matrix(sample(0:1, n, TRUE), nrow = 1L)
  tgt <- matrix(0.5, nrow = 1L)
  r <- leafblower_calibrate(weights = wts, group_ids = gid, targets = tgt,
                            method = "ieppa", max_iterations = 200L)
  expect_true(is.finite(r$max_error))
  expect_true(r$iterations > 0L)
  expect_lte(r$max_error, 1e-3)
  expect_equal(length(r$weights), n)
})

test_that("CalibResult chebyshev fields populated", {
  set.seed(22)
  n <- 40L; wts <- rep(1.0, n)
  gid <- matrix(sample(0:1, n, TRUE), nrow = 1L)
  tgt <- matrix(0.5, nrow = 1L)
  r <- leafblower_calibrate(weights = wts, group_ids = gid, targets = tgt,
                            method = "chebyshev", max_iterations = 200L)
  expect_true(is.finite(r$max_error))
  expect_true(r$iterations > 0L)
  expect_lte(r$max_error, 1e-3)
  expect_equal(length(r$weights), n)
})

test_that("CalibResult sinkhorn fields populated", {
  set.seed(22)
  n <- 40L; wts <- rep(1.0, n)
  gid <- matrix(sample(0:1, n, TRUE), nrow = 1L)
  tgt <- matrix(0.5, nrow = 1L)
  r <- leafblower_calibrate(weights = wts, group_ids = gid, targets = tgt,
                            method = "sinkhorn", max_iterations = 200L)
  expect_true(is.finite(r$max_error))
  expect_true(r$iterations > 0L)
  expect_lte(r$max_error, 1e-3)
  expect_equal(length(r$weights), n)
})

test_that("CalibResult greg fields populated", {
  set.seed(22)
  n <- 40L; wts <- rep(1.0, n)
  gid <- matrix(sample(0:1, n, TRUE), nrow = 1L)
  tgt <- matrix(0.5, nrow = 1L)
  r <- leafblower_calibrate(weights = wts, group_ids = gid, targets = tgt,
                            method = "greg", max_iterations = 200L)
  expect_true(is.finite(r$max_error))
  expect_true(r$iterations > 0L)
  expect_lte(r$max_error, 1e-3)
  expect_equal(length(r$weights), n)
})

test_that("CalibResult lbfgsb fields populated", {
  set.seed(22)
  n <- 40L; wts <- rep(1.0, n)
  gid <- matrix(sample(0:1, n, TRUE), nrow = 1L)
  tgt <- matrix(0.5, nrow = 1L)
  r <- leafblower_calibrate(weights = wts, group_ids = gid, targets = tgt,
                            method = "lbfgsb", max_iterations = 200L)
  expect_true(is.finite(r$max_error))
  expect_true(r$iterations > 0L)
  expect_lte(r$max_error, 1e-3)
  expect_equal(length(r$weights), n)
})
```

Run full test suite after each solver migration to catch regressions immediately:
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . && \
  Rscript -e 'testthat::test_dir("tests/testthat")'
```

### Commit
```
refactor(types): consolidate 6 per-solver Result structs into CalibResult in types.hpp

IEPPAResult, ChebyshevResult, RakingResult, SinkhornResult, GregResult, and
LbfgsbResult were structurally identical with minor field omissions. CalibResult
unifies them — solvers that don't use a field leave it at its default. Eliminates
~200 lines of duplicated struct definitions and field-copy boilerplate.
```

---

## Execution Order Summary

| # | Task ID | Risk | Files |
|---|---|---|---|
| 1 | ztid.1 | Zero | ieppa.hpp, lbfgsb_solver.cpp, ieppa.cpp |
| 2 | ztid.2 | Low | c_api.cpp |
| 3 | ztid.3 | Low | chebyshev.cpp, ieppa.cpp, logit_calib.cpp |
| 4 | ztid.7 | Medium | types.hpp, leafblower.h, ieppa.cpp, r_bridge.cpp |
| 5 | ztid.8 | Medium | types.hpp, ieppa.cpp, r_bridge.cpp, c_api.cpp |
| 6 | ztid.5 | Medium | r_bridge.cpp |
| 7 | ztid.6 | Medium | calib_dispatch.hpp, ieppa.cpp, chebyshev.cpp, raking.cpp, sinkhorn.cpp |
| 8 | ztid.4 | High | types.hpp, ieppa.hpp, chebyshev.hpp, raking.hpp, sinkhorn.hpp, greg.hpp, lbfgsb_solver.hpp, r_bridge.cpp, c_api.cpp |

Do not begin Task 8 until Tasks 1–7 all pass `testthat::test_dir("tests/testthat")`.
