# i0am-C: SRAA-m Correct Acceleration for All K Scales

**Date**: 2026-04-29
**Status**: Approved (rev4 — combined design)
**Ticket**: leafblower-i0am (continued)
**Files**: `src/sraa.hpp` (constants only), `src/greenkhorn.cpp`, `src/raking.cpp`,
           `tests/testthat/test-sraa-global.R`

---

## Problem

SRAA-m gives worse quality than plain on large-K overlapping-margin problems:

| Method | max_err | vs plain |
|---|---|---|
| greenkhorn plain | 1.57e-3 | baseline |
| greenkhorn+SRAA | 2.12e-3 | **+35% worse** |
| raking plain | 1.60e-3 | baseline |
| raking+SRAA | 1.75e-3 | +9% worse |

**Root cause:** `f_eval_sraa` sorts K margins ONCE at round entry (fixed sort). Plain
greenkhorn uses ADAPTIVE sort (re-sorts errRp after each step). These operators have
DIFFERENT fixed points — AA extrapolation finds the KL-optimal basin (2.12e-3), not the
max_err-optimal basin (1.57e-3). No per-step safeguard can distinguish the two basins in
early iterations.

---

## Design: Combined Fix (Tracks 2 + 1 simultaneously)

Two changes applied together — not sequential, not conditional:

### Change A: Adaptive-sort F_eval (eliminates basin multiplicity)

Change `f_eval_sraa` in `greenkhorn.cpp` from fixed-sort (sort once at entry) to
adaptive-sort (re-sort errRp after each greedy step, identical to plain greenkhorn).
This gives the AA operator a unique max_err-optimal fixed point.

```cpp
// REMOVE (sort once at entry):
std::iota(order_sraa.begin(), order_sraa.end(), 0);
std::stable_sort(order_sraa.begin(), order_sraa.end(),
    [&](int a, int b){ return errRp[a] > errRp[b]; });
for (int ki = 0; ki < K; ki++) greenkhorn_step(order_sraa[ki]);

// REPLACE WITH (re-sort after each step = plain greenkhorn):
for (int ki = 0; ki < K; ki++) {
    int k_star = (int)(std::max_element(errRp.begin(), errRp.end()) - errRp.begin());
    greenkhorn_step(k_star);
}
```

Also remove the now-dead `std::vector<int> order_sraa` declaration and `order_sraa.assign(K,0)` init.

`src/raking.cpp`: no change — raking's `F_eval` (line ~269) already uses adaptive order.

### Change B: Outer revert (recovery from transient AA overshoots)

After each SRAA outer iteration, compare outer quality (from adaptive-sort `errRp`) against
`best_errRp` (already tracked). If quality regresses > 10% above best for 5 consecutive
outer iterations: revert `X = X_best` and restart AA history.

**Why plateau gating is NOT needed (unlike the previous design):**
With Change A, F has a unique fixed point. Transient AA overshoots (from sparse early
history) produce quality detectably worse than the correct-basin `best_errRp`. The revert
catches them without needing to delay AA activation. No `aa_unlocked` field needed.

**New constants in `src/sraa.hpp`** (after `kSRAARestartGamma`):
```cpp
static constexpr double kSRAAOuterSlack       = 0.10; // 10% above best_errRp → stall
static constexpr int    kSRAAOuterStallWindow = 5;    // 5 stall iters → revert+restart
```

**`greenkhorn.cpp` outer loop addition** (inside `if (st.accelerate && K > 0)`, after errRp recompute):
```cpp
// Local at outer loop entry:
int outer_stall_count = 0;  // declared BEFORE the outer for-loop

// After existing X_best update block (which already computes curr_max):
// NOTE: Do NOT duplicate curr_max or best_errRp update — MERGE with existing block
if (curr_max > best_errRp * (1.0 + kSRAAOuterSlack)) {
    if (++outer_stall_count >= kSRAAOuterStallWindow) {
        X = X_best;               // revert to outer-quality best (copy, not swap)
        grk_sraa.clear();         // restart AA history
        outer_stall_count = 0;
        // Rebuild W, S_flat, errRp from reverted X
        W = 0.0; std::fill(S_flat.begin(), S_flat.end(), 0.0);
        for (int c = 0; c < M; c++) {
            W += X[c];
            for (int k = 0; k < K; k++) {
                int g = ct.g_per_cell[k][c];
                if (g >= 0 && g < st.cat_counts[k])
                    S_flat[k * S_stride + g] += X[c];
            }
        }
        for (int k = 0; k < K; k++) errRp[k] = compute_errRp_k(k);
        best_errRp = *std::max_element(errRp.begin(), errRp.end());
    }
} else { outer_stall_count = 0; }
```

**`raking.cpp` outer loop addition** (inside SRAA path, after existing `if (r.err_result < best_metric_seen)` block):
```cpp
// Local at SRAA block entry:
int rk_outer_stall_count = 0;

// After existing best_metric_seen update:
{
    double curr_quality_rk = r.err_result;
    if (curr_quality_rk > best_metric_seen * (1.0 + kSRAAOuterSlack)) {
        if (++rk_outer_stall_count >= kSRAAOuterStallWindow) {
            X = W_best;               // revert (no S_flat rebuild — F_eval handles it)
            rk_sraa.clear();
            rk_outer_stall_count = 0;
        }
    } else { rk_outer_stall_count = 0; }
}
```

---

## Files changed

| File | Change |
|---|---|
| `src/sraa.hpp` | Add 2 constants (`kSRAAOuterSlack`, `kSRAAOuterStallWindow`) |
| `src/greenkhorn.cpp` | Adaptive sort + remove `order_sraa` + outer revert |
| `src/raking.cpp` | Outer revert (scalar `r.err_result`; no S_flat rebuild) |
| `tests/testthat/test-sraa-global.R` | Add K=9 stepstone tests + outer revert test |

---

## Pre-implementation baseline (REQUIRED)

Before touching any code:
```bash
OMP_NUM_THREADS=1 Rscript benchmarks/stepstone_all_methods.R 2>&1 | grep -E "greenkhorn|raking"
```
Record `max_err` for greenkhorn+squarem and raking+squarem. Required in commit message.

---

## TDD Requirements

### T_sraa_adaptive_K9
```r
test_that("T_sraa_adaptive_K9: greenkhorn+AA K=9 max_err <= plain (stepstone)", {
  skip_if_not_installed("arrow"); skip_if_not_installed("jsonlite")
  skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
          "stepstone benchmark data not available")
  df  <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
  df$uuid <- NULL
  tgt <- lapply(jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json"),
                function(t) { t <- unlist(t); t / sum(t) })
  for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
  r_aa    <- suppressWarnings(harvest(df, tgt, method="greenkhorn", accelerate=TRUE,
                                       max_weight=5, min_weight=0, max_iterations=5000L,
                                       attach_weights=FALSE, verbose=0))
  r_plain <- suppressWarnings(harvest(df, tgt, method="greenkhorn", accelerate=FALSE,
                                       max_weight=5, min_weight=0, max_iterations=5000L,
                                       attach_weights=FALSE, verbose=0))
  me_aa    <- attr(r_aa,    "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  expect_lte(me_aa, me_plain * 1.001 + 1e-10,
    label=sprintf("K=9 SRAA (%.4e) must not exceed plain (%.4e)", me_aa, me_plain))
})
```

### T_sraa_raking_K9
```r
test_that("T_sraa_raking_K9: raking+AA K=9 max_err <= plain (stepstone)", {
  skip_if_not_installed("arrow"); skip_if_not_installed("jsonlite")
  skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
          "stepstone benchmark data not available")
  df  <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
  df$uuid <- NULL
  tgt <- lapply(jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json"),
                function(t) { t <- unlist(t); t / sum(t) })
  for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
  r_aa    <- suppressWarnings(harvest(df, tgt, method="raking", accelerate=TRUE,
                                       max_weight=5, min_weight=0, max_iterations=5000L,
                                       attach_weights=FALSE, verbose=0))
  r_plain <- suppressWarnings(harvest(df, tgt, method="raking", accelerate=FALSE,
                                       max_weight=5, min_weight=0, max_iterations=5000L,
                                       attach_weights=FALSE, verbose=0))
  me_aa    <- attr(r_aa,    "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  expect_lte(me_aa, me_plain * 1.001 + 1e-10,
    label=sprintf("K=9 raking+AA (%.4e) must not exceed plain (%.4e)", me_aa, me_plain))
})
```

### T_sraa_outer_revert (K=6 cross-margin — exercises the revert path)
```r
test_that("T_sraa_outer_revert: greenkhorn+AA K=6 quality <= plain with combined fix", {
  set.seed(42); n <- 8000L
  df <- data.frame(gender=factor(sample(c("M","F"),n,TRUE)),
                   time=factor(sample(1:3,n,TRUE)), age=factor(sample(1:4,n,TRUE)))
  df$gt <- factor(paste0(df$gender,df$time))
  df$ga <- factor(paste0(df$gender,df$age))
  df$ta <- factor(paste0(df$time,df$age))
  gt_t <- table(df$gt)/n; ga_t <- table(df$ga)/n; ta_t <- table(df$ta)/n
  tgt <- list(
    gender=c(M=0.48,F=0.52), time=setNames(c(0.4,0.35,0.25),1:3),
    age=setNames(c(0.3,0.25,0.25,0.2),1:4),
    gt={t<-setNames(as.numeric(gt_t)*c(0.95,1.02,0.98,1.03,0.97,1.05),names(gt_t));t/sum(t)},
    ga={t<-setNames(as.numeric(ga_t),names(ga_t));t/sum(t)},
    ta={t<-setNames(as.numeric(ta_t),names(ta_t));t/sum(t)})
  tgt <- lapply(tgt, function(t) t/sum(t))
  r_aa    <- suppressWarnings(harvest(df,tgt,method="greenkhorn",accelerate=TRUE,
                                       max_iterations=200L,attach_weights=FALSE))
  r_plain <- suppressWarnings(harvest(df,tgt,method="greenkhorn",accelerate=FALSE,
                                       max_iterations=200L,attach_weights=FALSE))
  expect_lte(attr(r_aa,"result")$max_error, attr(r_plain,"result")$max_error * 1.001 + 1e-10,
    label=sprintf("K=6 combined fix: AA (%.2e) must not exceed plain (%.2e)",
                  attr(r_aa,"result")$max_error, attr(r_plain,"result")$max_error))
})
```

### Existing tests must remain GREEN
T_sraa_grk (K=2 quality ≤ plain, faster), T_sraa_rk (K=2 raking), T_sraa_global (K=4),
T_sraa_ldlt_fallback, T_sraa_restart, T5-T8, T_logit_*

---

## Acceptance Criteria

| # | Criterion | Verify |
|---|-----------|--------|
| AC1 | T_sraa_grk GREEN (K=2, quality + speed preserved) | `devtools::test(filter="calibration-solvers")` |
| AC2 | T_sraa_rk GREEN (K=2 raking) | Same |
| AC3 | T_sraa_adaptive_K9 GREEN (K=9 greenkhorn) | Same (skip if parquet absent) |
| AC4 | T_sraa_raking_K9 GREEN (K=9 raking) | Same |
| AC5 | T_sraa_outer_revert GREEN (K=6 combined) | Same |
| AC6 | FAIL count = 3 (test-ieppa-nonuniform-d.R:28/:29, test-sor.R:18) | `devtools::test()` |
| AC7 (benchmark) | greenkhorn+AA iters < plain iters on K=9 | `Rscript benchmarks/stepstone_all_methods.R` |

---

## Calibration of constants

| Constant | Value | Derivation |
|---|---|---|
| `kSRAAOuterSlack` | 0.10 | K=9 basin escape = 35% regression. 10% catches genuine escape without false positives from normal transients (≤5% per step for adaptive operator). |
| `kSRAAOuterStallWindow` | 5 | 5 consecutive outer iterations above threshold; prevents false triggers from single transient; only ~0.5% of K=9 stepstone budget wasted before revert. |

---

## Out of Scope

- Plateau gating (`aa_unlocked` field, `kSRAAPlateauEps`, `kSRAAPlateauWindow`) — not needed with adaptive sort
- logit, ieppa, sinkhorn — do not use sraa_step, not affected
- Exposing constants to R user
