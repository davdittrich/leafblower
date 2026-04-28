# Convergence Status Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `RK_ERR_BUDGET=4` and `RK_ERR_STALL=5` status codes; replace errRp-based stall detector in `raking_solve` with weight KL stall; add `convergence_reason` to harvest.R result.

**Architecture:** Single source of truth: `src/leafblower.h` defines codes → `src/raking.cpp` emits them → `R/harvest.R` maps to warnings + character field. No rk_result_t struct change; no ABI change beyond new `#define` values.

**Tech Stack:** C++ (`leafblower.h`, `raking.cpp`), R (`harvest.R`, 4 test files). No new structs, no bridge changes.

**Beads epic:** `leafblower-afbd`

**Spec:** `docs/superpowers/specs/2026-04-28-convergence-status-design.md`

---

## File Map

| File | Change |
|------|--------|
| `src/leafblower.h` | Add `RK_ERR_BUDGET=4`, `RK_ERR_STALL=5` |
| `src/raking.cpp` | Replace `min_errRp_window`/errRp stall with wkl stall; BUDGET at post-loop; STALL at stall-break |
| `R/harvest.R` | Status 4/5 warnings; `convergence_reason` field; `@return` doc update |
| `tests/testthat/test-calibration-solvers.R` | Line 204: `c(0L,1L)` → `c(0L,4L,5L)` |
| `tests/testthat/test-calib-linalg.R` | Line 47: same pattern |
| `tests/testthat/test-convergence-criteria.R` | Line 236: `status==1L` → `status %in% c(4L,5L)` |
| `tests/testthat/test-ieppa-persistent-infeas.R` | Comment update only |

---

## Task 0: Write RED tests

> Commit these BEFORE any implementation. They must fail (RED) against current code.

**Files:**
- Modify: `tests/testthat/test-calibration-solvers.R` (append)

- [ ] **Step 1: Append RED tests**

Add at end of `tests/testthat/test-calibration-solvers.R`:

```r
# ── Convergence status v2 tests ──────────────────────────────────────────────

test_that("status-budget: budget exhausted emits status=4, not status=1", {
  # Far-from-converged problem, tiny budget → budget exhausted (status=4).
  # RED before implementation: raking_solve still emits status=1 (RK_ERR_NOCONV).
  set.seed(99L)
  df <- data.frame(
    v1 = factor(sample(5L, 2000L, TRUE)),
    v2 = factor(sample(4L, 2000L, TRUE)),
    v3 = factor(sample(3L, 2000L, TRUE))
  )
  tgt <- list(
    v1 = setNames(rep(0.2, 5), as.character(1:5)),
    v2 = setNames(c(0.4, 0.3, 0.2, 0.1), as.character(1:4)),
    v3 = setNames(c(0.5, 0.3, 0.2), as.character(1:3))
  )
  r <- leafblower::harvest(df, tgt, method = "raking", accelerate = FALSE,
    max_weight = 5, max_iterations = 5L, attach_weights = FALSE)
  expect_equal(attr(r, "result")$status, 4L,
               label = "budget exhausted must return status=4 (RK_ERR_BUDGET)")
})

test_that("status-stall: wkl plateau emits status=5 with convergence_reason='stall_kl'", {
  # Constrained problem — all cells at U_cell in some categories.
  # KL plateau fires before budget (max_iterations=1000) → status=5.
  # RED before implementation: raking_solve still emits status=1.
  set.seed(7L)
  df <- data.frame(
    v1 = factor(sample(4L, 300L, TRUE)),
    v2 = factor(sample(3L, 300L, TRUE))
  )
  tgt <- list(
    v1 = c("1"=0.4,"2"=0.3,"3"=0.2,"4"=0.1),
    v2 = c("1"=0.5,"2"=0.3,"3"=0.2)
  )
  r <- leafblower::harvest(df, tgt, method = "raking", accelerate = FALSE,
    max_weight = 2, max_iterations = 1000L, attach_weights = FALSE)
  res <- attr(r, "result")
  expect_equal(res$status, 5L,
               label = "KL plateau must return status=5 (RK_ERR_STALL)")
  expect_equal(res$convergence_used$convergence_reason, "stall_kl",
               label = "convergence_reason must be 'stall_kl' for flat loop KL stall")
})

test_that("status-perfect: perfect calibration exits status=0 not status=5", {
  # 1-category 1-margin problem already calibrated: wkl=0 → RK_OK, not stall.
  # RED before implementation: might stall or error.
  df  <- data.frame(v1 = factor(rep("1", 20L)))
  tgt <- list(v1 = c("1" = 1.0))
  r <- leafblower::harvest(df, tgt, method = "raking", accelerate = FALSE,
    max_weight = 5, max_iterations = 500L, attach_weights = FALSE)
  expect_equal(attr(r, "result")$status, 0L,
               label = "perfect calibration must return status=0, not status=5")
})
```

- [ ] **Step 2: Run tests to confirm RED**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-calibration-solvers.R', filter='status')" 2>&1
```

Expected: `status-budget` FAILS (gets status=1), `status-stall` FAILS (gets status=1), `status-perfect` PASSES (already converges with status=0).

- [ ] **Step 3: Commit RED tests**

```bash
git add tests/testthat/test-calibration-solvers.R
git commit -m "test(status): RED tests for BUDGET/STALL codes before implementation"
```

---

## Task 1: Add status codes to leafblower.h

**Files:**
- Modify: `src/leafblower.h`

- [ ] **Step 1: Append new codes after existing definitions**

Find in `src/leafblower.h` (lines 31-35):
```c
/* ── Return codes ── */
#define RK_OK         0  /* Success */
#define RK_ERR_NOCONV 1  /* Did not converge within outer_max_iter */
#define RK_ERR_INFEAS 2  /* Infeasible: empty cell with positive target */
#define RK_ERR_BADARG 3  /* Invalid argument */
```

Replace with:
```c
/* ── Return codes ── */
#define RK_OK           0  /* Converged: improvement criterion satisfied */
#define RK_ERR_NOCONV   1  /* Legacy alias — no longer emitted by new solvers */
#define RK_ERR_INFEAS   2  /* Infeasible: empty cell with positive target */
#define RK_ERR_BADARG   3  /* Invalid argument */
#define RK_ERR_BUDGET   4  /* Budget exhausted while loss still decreasing; increase max_iterations */
#define RK_ERR_STALL    5  /* Loss function plateau — at constrained optimum; weights are valid */
```

- [ ] **Step 2: Compile gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```
Expected: `* DONE (leafblower)`.

- [ ] **Step 3: Commit**

```bash
git add src/leafblower.h
git commit -m "feat(status): add RK_ERR_BUDGET=4, RK_ERR_STALL=5 status codes"
```

---

## Task 2: Replace errRp stall with weight KL stall in raking_solve (flat loop)

**Files:**
- Modify: `src/raking.cpp`

Three changes in the flat loop (the `else { for (int iter ...) }` branch).

- [ ] **Step 1: Rename min_errRp_window → min_loss_window; initialize wkl tracking**

Find (around line 121):
```cpp
    // Descent monitor
    double min_errRp_window = std::numeric_limits<double>::infinity();
    int n_no_improve = 0;
```
Replace with:
```cpp
    // Descent monitor — tracks solver loss function (weight KL for flat loop).
    double min_loss_window = std::numeric_limits<double>::infinity();
    int n_no_improve = 0;
```

- [ ] **Step 2: Replace flat loop stall block with weight KL stall**

Find the flat loop stall block (around lines 544-553):
```cpp
                if (!std::isfinite(min_errRp_window)) {
                    min_errRp_window = errRp; n_no_improve = 0;
                } else {
                    const double rel_eps = 0.01 * min_errRp_window;
                    const double eps = std::max(rel_eps, st.tol_abs);
                    if (errRp < min_errRp_window - eps) {
                        min_errRp_window = errRp; n_no_improve = 0;
                    } else {
                        n_no_improve++;
                    }
                }
```
Replace with:
```cpp
                // Weight KL stall: monotone for water-filling IPF (Csiszar-Tusnady).
                // KL plateau ↔ constrained KL minimum — correct stall signal.
                // Guard: wkl ≤ tol_abs means effectively at optimum → converged (not stalled).
                const double wkl_flat = compute_weight_kl();
                if (wkl_flat <= st.tol_abs) {
                    res.status = RK_OK; res.convergence_iter = iter; break;
                }
                if (!std::isfinite(min_loss_window)) {
                    min_loss_window = wkl_flat; n_no_improve = 0;
                } else if (wkl_flat < min_loss_window * (1.0 - st.convergence_cfg.pct_tol)) {
                    min_loss_window = wkl_flat; n_no_improve = 0;
                } else {
                    n_no_improve++;
                }
```

- [ ] **Step 3: Change flat loop stall exit code**

Find (around line 579):
```cpp
                if (n_no_improve >= kMaxNoImprove) {
                    res.status = RK_ERR_NOCONV;
```
Replace with:
```cpp
                if (n_no_improve >= kMaxNoImprove) {
                    res.status = RK_ERR_STALL;
```

- [ ] **Step 4: Change SQUAREM stall exit code (errRp-based, but emit STALL not NOCONV)**

Find the SQUAREM stall (around line 457):
```cpp
                if (n_no_improve >= kMaxNoImprove) { res.status = RK_ERR_NOCONV; break; }
```
Replace with:
```cpp
                if (n_no_improve >= kMaxNoImprove) { res.status = RK_ERR_STALL; break; }
```

Note: SQUAREM keeps the errRp-based stall criterion (weight KL not monotone for extrapolation steps). Only the exit code changes.

Also update the SQUAREM stall tracking variable references from `min_errRp_window` → `min_loss_window` (5 occurrences in lines ~447-456):

Find in the SQUAREM stall block:
```cpp
                if (!std::isfinite(min_errRp_window)) {
                    min_errRp_window = errRp_new; n_no_improve = 0;
                } else {
                    const double eps = std::max(0.01 * min_errRp_window, st.tol_abs);
                    if (errRp_new < min_errRp_window - eps) {
                        min_errRp_window = errRp_new; n_no_improve = 0;
```
Replace with:
```cpp
                if (!std::isfinite(min_loss_window)) {
                    min_loss_window = errRp_new; n_no_improve = 0;
                } else {
                    const double eps = std::max(0.01 * min_loss_window, st.tol_abs);
                    if (errRp_new < min_loss_window - eps) {
                        min_loss_window = errRp_new; n_no_improve = 0;
```
(5 occurrences — all 5 must change; use replace_all on this block)

- [ ] **Step 5: Change budget exit code (post-loop)**

Find (around line 65):
```cpp
    res.status     = RK_ERR_NOCONV;
```
Replace with:
```cpp
    res.status     = RK_ERR_BUDGET;  // initial; overwritten by criterion/stall; remains if budget exhausted
```

- [ ] **Step 6: Compile gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -10
```
Expected: `* DONE (leafblower)`.

- [ ] **Step 7: Smoke test — verify basic raking still works**

```bash
Rscript -e "
  set.seed(1L); n <- 500L
  df  <- data.frame(v1=factor(sample(3L,n,TRUE)), v2=factor(sample(2L,n,TRUE)))
  tgt <- list(v1=c('1'=0.5,'2'=0.3,'3'=0.2), v2=c('1'=0.6,'2'=0.4))
  r <- leafblower::harvest(df, tgt, method='raking', max_weight=5,
       max_iterations=500L, attach_weights=FALSE)
  cat('status:', attr(r,'result')\$status,
      'max_err:', attr(r,'result')\$max_error, '\n')
  stopifnot(attr(r,'result')\$status %in% c(0L, 4L, 5L))
  cat('PASS\n')
"
```

- [ ] **Step 8: Run status RED tests — should now PASS**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-calibration-solvers.R', filter='status')" 2>&1
```

Expected: all 3 PASS (`status-budget`=4, `status-stall`=5, `status-perfect`=0).

- [ ] **Step 9: Commit**

```bash
git add src/raking.cpp
git commit -m "feat(status): replace errRp stall with weight KL stall; emit BUDGET/STALL codes"
```

---

## Task 3: Update harvest.R

**Files:**
- Modify: `R/harvest.R`

Three changes: (1) handle new status codes, (2) add `convergence_reason` to result, (3) update `@return` docs.

- [ ] **Step 1: Replace the existing status=1 warning AND add status 4/5 warnings**

Find (around line 323) the EXISTING status=1 block:
```r
  if (calib_result$status == 1L)
    warning("leafblower: calibration did not converge (max_error=",
            signif(calib_result$max_error, 3), "). Weights reflect last iterate.")
```

Replace it with (handles all soft-exit statuses in one place, no duplicates):
```r
  if (calib_result$status == 4L)
    warning("leafblower: budget exhausted — weights reflect best iterate; ",
            "increase max_iterations if further improvement is needed")
  if (calib_result$status == 5L && isTRUE(accelerate_bool))
    warning("leafblower: SQUAREM errRp plateau — weights are valid; ",
            "try accelerate=FALSE for guaranteed KL-minimum stall detection")
  if (calib_result$status == 5L && !isTRUE(accelerate_bool))
    warning("leafblower: loss function plateau — at constrained optimum given bounds; ",
            "weights are valid; no further improvement is achievable")
  if (calib_result$status == 1L)
    warning("leafblower: did not converge (legacy status code from solver not yet updated)")
```

- [ ] **Step 2: Add convergence_reason to convergence_used list**

Find the `calib_result$convergence_used <- list(...)` block (around line 288). Add `convergence_reason` as the last field:

```r
  calib_result$convergence_used <- list(
    metric           = .safe_lookup(.metric_names, calib_result$convergence_metric),
    rule             = .safe_lookup(.rule_names,   calib_result$convergence_rule),
    tol              = calib_result$convergence_tol,
    fired_at_iter    = calib_result$convergence_iter,
    solver_objective = calib_result$solver_objective,
    minimized_metric = .safe_lookup(.metric_names, calib_result$convergence_minimized_metric),
    convergence_reason = {
      s <- calib_result$status
      if      (is.null(s) || is.na(s))  NA_character_
      else if (s == 0L)                  "criterion"
      else if (s == 4L)                  "budget"
      else if (s == 5L && isTRUE(accelerate_bool)) "stall_errRp"
      else if (s == 5L)                  "stall_kl"
      else if (s == 2L)                  "infeasible"
      else if (s == 3L)                  "error"
      else                               "legacy"   # status=1 or unknown
    }
  )
```

- [ ] **Step 3: Update @return documentation**

Find the `@return` section in the roxygen block (around line 81 — the `result` list description). Add after the `best_weights` bullet:
```r
#'         \item \code{convergence_used$convergence_reason}: Character string.
#'           Why the solver exited: \code{"criterion"} (improvement criterion satisfied),
#'           \code{"budget"} (max_iterations exhausted, still improving — increase budget),
#'           \code{"stall_kl"} (weight KL plateau — at constrained KL minimum),
#'           \code{"stall_errRp"} (SQUAREM errRp plateau — empirically near optimum),
#'           \code{"infeasible"}, \code{"error"}, or \code{"legacy"}.
```

- [ ] **Step 4: Build and verify**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

```bash
Rscript -e "
  set.seed(1L); n <- 500L
  df  <- data.frame(v1=factor(sample(3L,n,TRUE)), v2=factor(sample(2L,n,TRUE)))
  tgt <- list(v1=c('1'=0.5,'2'=0.3,'3'=0.2), v2=c('1'=0.6,'2'=0.4))
  r <- leafblower::harvest(df, tgt, method='raking', max_weight=5,
       max_iterations=500L, attach_weights=FALSE)
  res <- attr(r,'result')
  cat('status:', res\$status, 'reason:', res\$convergence_used\$convergence_reason, '\n')
  stopifnot(!is.na(res\$convergence_used\$convergence_reason))
  cat('PASS: convergence_reason is non-NA\n')
"
```

- [ ] **Step 5: Commit**

```bash
git add R/harvest.R
git commit -m "feat(status): harvest.R — status 4/5 warnings, convergence_reason field, @return docs"
```

---

## Task 4: Update 4 test files

**Files:**
- Modify: 4 test files (mechanical updates to status assertions)

- [ ] **Step 1: test-calibration-solvers.R line 204**

Find:
```r
  expect_true(r$status %in% c(0L, 1L),
              info=sprintf("expected 0 (OK) or 1 (NOCONV), got %d", r$status))
```
Replace with:
```r
  expect_true(r$status %in% c(0L, 4L, 5L),
              info=sprintf("expected 0 (OK) or 4 (BUDGET) or 5 (STALL), got %d", r$status))
```

- [ ] **Step 2: test-calib-linalg.R line 47**

Find:
```r
  expect_true(r_greg$status %in% c(0L, 1L),
```
Replace with:
```r
  expect_true(r_greg$status %in% c(0L, 4L, 5L),
```

- [ ] **Step 3: test-convergence-criteria.R line 236**

Find:
```r
  if (result$status == 1L) {
    expect_lt(result$best_error, 0.9 * result$max_error)
  }
```
Replace with:
```r
  if (result$status %in% c(4L, 5L)) {
    expect_lt(result$best_error, 0.9 * result$max_error)
  }
```

- [ ] **Step 4: test-ieppa-persistent-infeas.R — comment update**

Find (line 63):
```r
  # hits max_iter -> RK_ERR_NOCONV with high errRp. This test guards that
```
Replace with:
```r
  # hits max_iter -> RK_ERR_BUDGET (status=4) with high errRp. This test guards that
```

- [ ] **Step 5: Run full test suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```
Expected: FAIL ≤ 3 (pre-existing).

- [ ] **Step 6: Commit**

```bash
git add tests/testthat/test-calibration-solvers.R \
        tests/testthat/test-calib-linalg.R \
        tests/testthat/test-convergence-criteria.R \
        tests/testthat/test-ieppa-persistent-infeas.R
git commit -m "test(status): update 4 test files — status %in% c(0L,4L,5L) instead of c(0L,1L)"
```

---

## Task 5: Verify ACs and local benchmark

**Files:** none — verification only.

- [ ] **Step 1: Run all status tests**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-calibration-solvers.R', filter='status')" 2>&1
```
Expected: `status-budget`, `status-stall`, `status-perfect` all PASS.

- [ ] **Step 2: AC7 — convergence_reason non-NA for all paths**

```bash
Rscript -e "
  # status=0 (criterion): normal converging problem
  df0 <- data.frame(v1=factor(c('1','2','1','2')))
  t0  <- list(v1=c('1'=0.5,'2'=0.5))
  r0  <- suppressWarnings(leafblower::harvest(df0, t0, method='raking',
           max_weight=5, max_iterations=500L, attach_weights=FALSE))
  cat('status=0 path: status=', attr(r0,'result')\$status,
      'reason=', attr(r0,'result')\$convergence_used\$convergence_reason, '\n')
  stopifnot(!is.na(attr(r0,'result')\$convergence_used\$convergence_reason))

  # status=4 (budget): far-from-converged, tiny budget
  set.seed(99L)
  df4 <- data.frame(v1=factor(sample(5L,2000L,TRUE)), v2=factor(sample(4L,2000L,TRUE)),
                    v3=factor(sample(3L,2000L,TRUE)))
  t4  <- list(v1=setNames(rep(0.2,5),as.character(1:5)),
              v2=setNames(c(0.4,0.3,0.2,0.1),as.character(1:4)),
              v3=setNames(c(0.5,0.3,0.2),as.character(1:3)))
  r4  <- suppressWarnings(leafblower::harvest(df4, t4, method='raking',
           max_weight=5, max_iterations=5L, attach_weights=FALSE))
  cat('status=4 path: status=', attr(r4,'result')\$status,
      'reason=', attr(r4,'result')\$convergence_used\$convergence_reason, '\n')
  stopifnot(!is.na(attr(r4,'result')\$convergence_used\$convergence_reason))

  # status=5 (stall): constrained problem, wkl plateau before budget
  set.seed(7L)
  df5 <- data.frame(v1=factor(sample(4L,300L,TRUE)), v2=factor(sample(3L,300L,TRUE)))
  t5  <- list(v1=c('1'=0.4,'2'=0.3,'3'=0.2,'4'=0.1), v2=c('1'=0.5,'2'=0.3,'3'=0.2))
  r5  <- suppressWarnings(leafblower::harvest(df5, t5, method='raking',
           max_weight=2, max_iterations=1000L, attach_weights=FALSE))
  cat('status=5 path: status=', attr(r5,'result')\$status,
      'reason=', attr(r5,'result')\$convergence_used\$convergence_reason, '\n')
  stopifnot(!is.na(attr(r5,'result')\$convergence_used\$convergence_reason))

  cat('AC7: PASS\n')
"
```

- [ ] **Step 3: AC8 — stepstone benchmark (local only)**

```bash
OMP_NUM_THREADS=1 Rscript -e "
  suppressPackageStartupMessages({library(arrow);library(jsonlite);library(leafblower)})
  df  <- arrow::read_parquet('benchmarks/stepstone_fulldata_bench_data.parquet')
  df\$uuid <- NULL
  tgt <- lapply(jsonlite::fromJSON('benchmarks/stepstone_fulldata_bench_targets.json'),
                function(t){t<-unlist(t);t/sum(t)})
  for(nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
  r <- suppressWarnings(leafblower::harvest(df, tgt, method='raking', accelerate=FALSE,
        max_weight=5, max_iterations=5000L, attach_weights=FALSE, verbose=0))
  res <- attr(r,'result')
  cat('status:', res\$status, 'reason:', res\$convergence_used\$convergence_reason,
      'max_err:', res\$max_error, '\n')
  stopifnot(res\$status == 5L)
  stopifnot(res\$max_error <= 2.97e-3)
  cat('AC8: PASS\n')
" 2>&1 | grep -v "^Welcome\|^Working"
```

Expected: `status: 5 reason: stall_kl max_err: X.XXe-3` where max_err ≤ 2.97e-3.

- [ ] **Step 4: Run full test suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```
Expected: FAIL ≤ 3 (pre-existing).

---

## Self-Review Against Spec

**AC coverage**:
- AC1 (criterion → 0): unchanged; Task 2 doesn't touch criterion path ✓
- AC2 (budget → 4): Task 2 Step 5 changes initial `res.status = RK_ERR_BUDGET` ✓
- AC3 (stall → 5): Task 2 Step 3 emits `RK_ERR_STALL` ✓
- AC4 (wkl=0 → 0): Task 2 Step 2 adds `if (wkl_flat <= st.tol_abs) → RK_OK` ✓
- AC5 (4/5 return weights, not stop): Task 3 Step 1 uses `warning()` not `stop()` ✓
- AC6 (2/3 still stop): Task 3 doesn't touch status 2/3 ✓
- AC7 (convergence_reason non-NA): Task 3 Step 2 maps all codes including `"legacy"` for status=1 ✓
- AC8 (stepstone stall_kl, max_err ≤ 2.97e-3): Task 5 Step 3 ✓ (local)
- AC9 (FAIL ≤ 3): Task 5 Step 4 ✓

**Type consistency**: `RK_ERR_BUDGET` and `RK_ERR_STALL` used consistently across leafblower.h and raking.cpp. `convergence_reason` is always character (never NA). ✓

**No placeholders**: all code blocks complete with exact text. ✓
