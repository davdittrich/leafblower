# SQUAREM Geometry Fix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix SQUAREM's CBB step geometry from cell-level to obs-level weighted norms; replace errRp stall with weight-change stall. Expected: SQUAREM reaches max_err=1.60e-3 (matching flat loop and autumn) instead of current 1.71e-3.

**Architecture:** Two changes inside the SQUAREM `while` loop only (`if (st.accelerate) { ... }`). Flat loop untouched. Single source file.

**Tech Stack:** C++17 (`src/raking.cpp`), R test file.

**Spec:** `docs/superpowers/specs/2026-04-28-squarem-geometry-fix.md`

---

## File Map

| File | Change |
|------|--------|
| `src/raking.cpp` | SQUAREM while loop: obs-level norms for α; cell-level norm for halving; weight-change stall; X_prev_sq snapshot |
| `tests/testthat/test-calibration-solvers.R` | Append RED test + AC1 skip_if test |

---

## Task 0: Write RED test (before any code)

**Files:**
- Modify: `tests/testthat/test-calibration-solvers.R` (append)

- [ ] **Step 1: Append RED test + AC1 test at end of test file**

```r
# ── SQUAREM geometry fix tests ───────────────────────────────────────────────

test_that("squarem-geo-red: heterogeneous cells → convergence_reason='stall_kl' after geo fix", {
  # RED before fix: cell-level α over-extrapolates → stall_errRp (errRp stall fires early).
  # GREEN after fix: obs-level α correct → stall_kl (weight-change stall, reaches KL min).
  #
  # Construct: 1 large category (500 obs) + 1 small (5 obs); target forces heavy IPF step.
  # Large cell dominates cell-level ‖r_cell‖² but not obs-level ‖r_obs‖².
  # α_cell >> α_obs → cell-level over-extrapolates, halving brings α → -1 → no acceleration.
  set.seed(42L)
  df  <- data.frame(v1 = factor(c(rep("A", 500L), rep("B", 5L))))
  tgt <- list(v1 = c("A" = 0.4, "B" = 0.6))  # B (5 obs) forced to 60% → huge IPF step

  r   <- suppressWarnings(leafblower::harvest(df, tgt, method = "raking", accelerate = TRUE,
           max_weight = 5, max_iterations = 200L, attach_weights = FALSE))
  res <- attr(r, "result")

  # RED: convergence_reason == "stall_errRp" (cell-level α fires errRp stall early)
  # GREEN: convergence_reason == "stall_kl" (obs-level α reaches true KL minimum)
  expect_equal(res$convergence_used$convergence_reason, "stall_kl",
               label = "obs-level geometry must reach KL stall, not errRp stall")
})

test_that("squarem-geo-ac1: stepstone SQUAREM reaches flat-loop quality", {
  skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
          "stepstone dataset not available (local-only benchmark)")
  skip_if(!requireNamespace("arrow", quietly = TRUE), "arrow not installed")
  skip_if(!requireNamespace("jsonlite", quietly = TRUE), "jsonlite not installed")

  df  <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
  df$uuid <- NULL
  tgt <- lapply(jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json"),
                function(t) { v <- unlist(t); v / sum(v) })
  for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])

  r   <- suppressWarnings(leafblower::harvest(df, tgt, method = "raking", accelerate = TRUE,
           max_weight = 5, min_weight = 0, max_iterations = 5000L,
           attach_weights = FALSE, verbose = 0L))
  res <- attr(r, "result")

  expect_lte(res$max_error, 1.60e-3,
             label = "SQUAREM+geo fix must reach flat-loop quality (max_err ≤ 1.60e-3)")
  expect_equal(res$convergence_used$convergence_reason, "stall_kl",
               label = "must stall on KL (weight-change), not errRp")
})
```

- [ ] **Step 2: Run tests to confirm RED state**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-calibration-solvers.R')" 2>&1 | grep -E "squarem-geo|PASS|FAIL" | tail -10
```

Expected:
- `squarem-geo-red`: FAIL (`convergence_reason` is `"stall_errRp"` not `"stall_kl"`)
- `squarem-geo-ac1`: SKIP (stepstone available) or FAIL if dataset present

- [ ] **Step 3: Commit RED tests**

```bash
git add tests/testthat/test-calibration-solvers.R
git commit -m "test(squarem-geo): RED tests before geometry fix — obs-level α + weight-change stall"
```

---

## Task 1: Fix norm computation + add X_prev_sq + weight-change stall

**Files:**
- Modify: `src/raking.cpp`

All changes inside the `while (f_eval_count + 3 <= st.inner_max_iter)` loop body.

### Sub-task A: Replace norm block with dual obs/cell norms

- [ ] **Step 1: Find and replace the norm computation block (around line 345)**

Find:
```cpp
                double norm_r = 0.0, norm_v = 0.0, norm_w2 = 0.0;
                for (int c = 0; c < ct.M_cell; c++) {
                    double ri = w1[c] - X[c],  vi = w2[c] - w1[c];
                    norm_r  += ri * ri;
                    norm_v  += vi * vi;
                    norm_w2 += w2[c] * w2[c];
                }
                norm_r = std::sqrt(norm_r);
                norm_v = std::sqrt(norm_v);
                norm_w2 = std::sqrt(norm_w2);
```

Replace with:
```cpp
                // Obs-level norms (1/n_c weighted) for α computation.
                // n_per_cell[c] >= 1 guaranteed by build_cell_table.
                // Cell-level norms kept separately for halving criterion (both sides must match).
                double r_sq_obs = 0.0, v_sq_obs = 0.0, w2_sq_obs = 0.0;
                double v_sq_cell = 0.0;  // cell-level ‖v‖² for halving (cand_resid is also cell-level)
                for (int c = 0; c < ct.M_cell; c++) {
                    const double inv_nc = 1.0 / static_cast<double>(ct.n_per_cell[c]);
                    const double ri = w1[c] - X[c],  vi = w2[c] - w1[c];
                    r_sq_obs  += ri * ri * inv_nc;
                    v_sq_obs  += vi * vi * inv_nc;
                    w2_sq_obs += w2[c] * w2[c] * inv_nc;
                    v_sq_cell += vi * vi;
                }
                const double norm_r  = std::sqrt(r_sq_obs);   // obs-level ‖r‖
                const double norm_v  = std::sqrt(v_sq_obs);   // obs-level ‖v‖ (for α + fixed-pt guard)
                const double norm_w2 = std::sqrt(w2_sq_obs);  // obs-level ‖w2‖ (for fixed-pt guard)
```

### Sub-task B: Fix plain_resid to use cell-level v_sq

- [ ] **Step 2: Fix the step-halving reference norm**

Find (around line 383):
```cpp
                const double plain_resid = norm_v * norm_v;  // ‖v‖²
```

Replace with:
```cpp
                const double plain_resid = v_sq_cell;  // cell-level ‖v‖² — matches cand_resid geometry
```

(Both `cand_resid` and `plain_resid` are now cell-level — dimensionally consistent.)

### Sub-task C: Add X_prev_sq snapshot + weight-change stall

- [ ] **Step 3: Add X_prev_sq declaration before the while loop**

Find the line `int f_eval_count = 0;` (inside `if (st.inner_max_iter >= 3)`). Add immediately after:
```cpp
            int f_eval_count = 0;
            auto X_prev_sq = X;  // snapshot for weight-change stall (obs-level L1 Δw)
```

- [ ] **Step 4: Replace errRp stall block with weight-change stall**

Find the errRp stall block (around lines 447-457):
```cpp
                // errRp stall for SQUAREM: KL is non-monotone for CBB extrapolation steps —
                // accepted iterate KL can increase even with step-halving, making KL stall
                // fire spuriously. errRp of accepted iterates is approximately non-increasing.
                if (!std::isfinite(min_loss_window)) {
                    min_loss_window = errRp_new; n_no_improve = 0;
                } else {
                    const double eps = std::max(0.01 * min_loss_window, st.tol_abs);
                    if (errRp_new < min_loss_window - eps) {
                        min_loss_window = errRp_new; n_no_improve = 0;
                    } else {
                        n_no_improve++;
                    }
                }
                if (n_no_improve >= kMaxNoImprove) { res.status = RK_ERR_STALL; break; }
```

Replace with:
```cpp
                // Weight-change stall for SQUAREM: obs-level L1 Δw goes to zero at the fixed
                // point regardless of errRp oscillation. Same sliding-window relative improvement
                // pattern as the flat loop KL stall. X_prev_sq is the previous accepted iterate.
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

                X_prev_sq = X;  // update snapshot for next super-step
```

### Sub-task D: Compile gate

- [ ] **Step 5: Build**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -10
```
Expected: `* DONE (leafblower)`. If errors: read full output. Fix before continuing.

- [ ] **Step 6: Smoke test — verify obs-level α differs from cell-level**

```bash
Rscript -e "
  set.seed(42L)
  df  <- data.frame(v1=factor(c(rep('A',500L), rep('B',5L))))
  tgt <- list(v1=c('A'=0.4,'B'=0.6))
  r   <- suppressWarnings(leafblower::harvest(df, tgt, method='raking', accelerate=TRUE,
           max_weight=5, max_iterations=200L, attach_weights=FALSE, verbose=1))
  res <- attr(r,'result')
  cat('reason:', res\$convergence_used\$convergence_reason, '\n')
  cat('status:', res\$status, 'iters:', res\$iterations, '\n')
" 2>&1 | grep -v "^Welcome\|^Working"
```

Expected: `reason: stall_kl` (GREEN — weight-change stall fires instead of errRp stall).

- [ ] **Step 7: Run squarem geometry RED tests**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-calibration-solvers.R')" 2>&1 | grep -E "squarem-geo|PASS|FAIL" | tail -5
```

Expected: `squarem-geo-red` PASSES (GREEN), `squarem-geo-ac1` PASS or SKIP.

- [ ] **Step 8: Run full test suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```

Expected: FAIL ≤ 3 (pre-existing).

- [ ] **Step 9: Commit**

```bash
git add src/raking.cpp
git commit -m "feat(squarem-geo): obs-level weighted norms for CBB α; cell-level for halving; weight-change stall"
```

---

## Task 2: Stepstone benchmark and verify AC1

> **Local only** — requires stepstone parquet. Run before declaring done.

- [ ] **Step 1: Run stepstone comparison**

```bash
OMP_NUM_THREADS=1 Rscript -e "
  suppressPackageStartupMessages({library(arrow);library(jsonlite);library(leafblower);library(autumn)})
  df  <- arrow::read_parquet('benchmarks/stepstone_fulldata_bench_data.parquet')
  df\$uuid <- NULL
  tgt <- lapply(jsonlite::fromJSON('benchmarks/stepstone_fulldata_bench_targets.json'),
                function(t){t<-unlist(t);t/sum(t)})
  for(nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
  fit<-function(w,df,tgt){W<-sum(w);max_err<-0
    for(nm in names(tgt)){lv<-names(tgt[[nm]]);S<-tapply(w,df[[nm]],sum,default=0)[lv];S[is.na(S)]<-0
      max_err<-max(max_err,max(abs(S/W-tgt[[nm]])))}; max_err}

  run<-function(label,...){cat(sprintf('%-32s',label));flush.console()
    t0<-proc.time()['elapsed']
    r<-suppressWarnings(if(label=='autumn')autumn::harvest(df,tgt,max_weight=5,min_weight=0,accelerate=TRUE)
       else leafblower::harvest(df,tgt,max_weight=5,min_weight=0,max_iterations=5000L,attach_weights=FALSE,verbose=0,...))
    wall<-proc.time()['elapsed']-t0
    w<-tryCatch(as.numeric(r),error=function(e)if(!is.null(r\$weights))as.numeric(r\$weights) else stop(e))
    res<-attr(r,'result')
    cat(sprintf(' %5.1fs  iters=%4s  reason=%-10s  max_err=%.4e  DEFF=%.4f\n',
      wall,ifelse(is.null(res\$iterations)|is.na(res\$iterations),'  —',res\$iterations),
      ifelse(is.null(res\$convergence_used\$convergence_reason),'—',res\$convergence_used\$convergence_reason),
      fit(w,df,tgt), length(w)*sum(w^2)/sum(w)^2))}

  run('raking (flat, wf)',       method='raking', accelerate=FALSE)
  run('raking+SQUAREM (geo fix)',method='raking', accelerate=TRUE)
  run('autumn',                  label='autumn')
" 2>&1 | grep -v "^Welcome\|^Working"
```

- [ ] **Step 2: Verify AC1 (max_err ≤ 1.60e-3)**

Confirm: SQUAREM reason=stall_kl AND max_err ≤ 1.60e-3.

If max_err > 1.60e-3: investigate. Check if wchange stall is firing too early (increase kMaxNoImprove) or if halving is still collapsing α (log verbose=1 α values).

---

## Self-Review Against Spec

**Spec coverage:**
- Fix 1 (obs-level norms for α): Task 1 Sub-task A ✓
- Step-halving stays cell-level: Task 1 Sub-task B (`plain_resid = v_sq_cell`) ✓
- Fixed-point guard uses obs-level: norm_v is obs-level after Sub-task A ✓
- Fix 2 (weight-change stall): Task 1 Sub-task C ✓
- X_prev_sq initialization before while: Task 1 Sub-task C Step 3 ✓
- RED test committed before implementation: Task 0 ✓
- AC1 skip_if-guarded test: Task 0 ✓
- AC4 flat-loop unchanged: changes in `if (st.accelerate)` branch only ✓

**No placeholders**: all code blocks complete. ✓

**Type consistency**: `v_sq_cell` (double) used for `plain_resid`; `cand_resid` computed from cell-level diffs — both cell-level, consistent. ✓
