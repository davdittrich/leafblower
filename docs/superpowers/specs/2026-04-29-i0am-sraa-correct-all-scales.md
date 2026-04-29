# i0am-C: SRAA-m Correct Acceleration for All K Scales

**Date**: 2026-04-29
**Status**: Pending design review
**Ticket**: leafblower-i0am (continued)
**Files**: `src/sraa.hpp`, `src/greenkhorn.cpp`, `src/raking.cpp`,
           `tests/testthat/test-sraa-global.R`

---

## Problem

SRAA-m (Anderson Acceleration) gives worse quality than plain on large-K overlapping-margin
problems:

| Method | max_err | vs plain |
|---|---|---|
| greenkhorn plain | 1.57e-3 | baseline |
| greenkhorn+SRAA | 2.12e-3 | **+35% worse** |
| raking plain | 1.60e-3 | baseline |
| raking+SRAA | 1.75e-3 | +9% worse |

**Root cause (proven through investigation):** The fixed-point operator F used in f_eval_sraa
sorts K margins ONCE at round entry (stationary, fixed sort order). Plain greenkhorn uses
ADAPTIVE sort (re-sorts errRp after each step). These two operators have DIFFERENT fixed points:

- Adaptive-sort F → max_err-optimal fixed point (1.57e-3)
- Fixed-sort F → multiple fixed points including KL-optimal basin (2.12e-3)

AA extrapolation finds the KL-optimal basin of fixed-sort F. No per-step safeguard can
distinguish the two basins in early iterations because the wrong basin quality (2.12e-3)
is better than the unconverged correct-basin quality in early warmup.

**Failed approaches documented in leafblower-i0am notes:**
1. Global quality floor (best_err_seen from f_eval): fails — f_eval quality ≠ outer quality
2. Varadhan L2 residual: fails — both basins have similar fixed-point residuals
3. Strict L2 descent: fails — prevents AA entirely on K=9 (degrades to plain)
4. Outer revert without plateau gating: fails — wrong basin quality better than warmup quality

---

## Design: Two Complementary Tracks

### Track 2 (implement first, test empirically)
**Adaptive-sort F_eval** — change f_eval_sraa to re-sort errRp after EACH greedy step
(identical to plain greenkhorn). Gives AA a unique max_err-optimal fixed point.

### Track 1 (safety net, implement if Track 2 is insufficient)
**Plateau-Gated AA + Outer Revert** — AA disabled until outer quality improvement saturates
(correct basin established). Outer revert fires if AA escapes. Guarantees correctness at
cost of reduced large-K acceleration.

---

## Track 2: Adaptive-Sort F_eval

### What changes

**`src/greenkhorn.cpp`** — `f_eval_sraa` lambda, 5-line replacement:

```cpp
// BEFORE (fixed sort: sort once at round entry):
std::iota(order_sraa.begin(), order_sraa.end(), 0);
std::stable_sort(order_sraa.begin(), order_sraa.end(),
    [&](int a, int b){ return errRp[a] > errRp[b]; });
for (int ki = 0; ki < K; ki++) greenkhorn_step(order_sraa[ki]);

// AFTER (adaptive sort: re-sort errRp after each step, same as plain greenkhorn):
for (int ki = 0; ki < K; ki++) {
    int k_star = (int)(std::max_element(errRp.begin(), errRp.end()) - errRp.begin());
    greenkhorn_step(k_star);
}
```

**`src/raking.cpp`** — no change. Raking's `F_eval` lambda (line ~269) already performs
one full water-filling sweep in correct adaptive order (equivalent to K adaptive-sort steps).

**`src/sraa.hpp`** — no change.

### Why this works (or might work)

With adaptive sort, F is the same operator as plain greenkhorn's single-step update. This
operator has a UNIQUE fixed point (the max_err-optimal solution) because it always moves
toward the worst-error margin. AA extrapolation along this operator's trajectory stays
in the max_err-optimal basin.

### Risk

Breaks the "stationary operator" property that justifies AA convergence theory. AA becomes
a heuristic. Empirical validation required — AC3/AC4 tests confirm or deny.

### Outcome if Track 2 works

AC3 (K=9 stepstone greenkhorn+AA ≤ plain) and AC4 (raking+AA ≤ plain) pass. Track 1 is
not needed. Both methods get correct acceleration at all scales.

### Outcome if Track 2 fails

Track 1 is also implemented (belt-and-suspenders).

---

## Track 1: Plateau-Gated AA + Outer Revert

### Architecture

**New `SRAAState` field:**
```cpp
bool allow_aa = false;  // set externally by outer solver loop; prevents AA when false
```

**New constants in `sraa.hpp`:**
```cpp
static constexpr double kSRAAPlateauEps       = 1e-3;  // 0.1%/iter improvement → plateau
static constexpr int    kSRAAPlateauWindow     = 4;    // 4 consecutive plateau iters → AA enabled
static constexpr double kSRAAOuterSlack        = 0.10; // 10% above best_errRp → regression
static constexpr int    kSRAAOuterStallWindow  = 3;    // 3 regressed outer iters → revert
```

**`sraa_step` change** (one additional condition in step 5):
```cpp
// Not enough history OR outer loop has not enabled AA → plain step
if (state.count < kSRAAMinCount || !state.allow_aa) {
    std::swap(X, state.F_cur);
    return {false, 1, err_plain};
}
```
History still accumulates during plateau-detection phase (steps 3-4 run). AA fires only
when the outer loop sets `allow_aa = true`.

**Outer loop additions in `greenkhorn.cpp`** (after errRp recompute, within accelerate path):

```cpp
// Local variables at outer loop entry:
double prev_outer_quality = std::numeric_limits<double>::infinity();
int    plateau_count = 0, outer_stall_count = 0;

// ... (inside main outer loop, after errRp recompute) ...

if (st.accelerate && K > 0) {
    // Plateau detection: enable AA when outer quality improvement saturates
    if (std::isfinite(prev_outer_quality)) {
        double impr = (prev_outer_quality - curr_max) / std::max(prev_outer_quality, 1e-15);
        plateau_count = (impr < kSRAAPlateauEps) ? plateau_count + 1 : 0;
        if (plateau_count >= kSRAAPlateauWindow) grk_sraa.allow_aa = true;
    }
    prev_outer_quality = curr_max;

    // Outer revert: when AA escapes correct basin
    if (grk_sraa.allow_aa && curr_max > best_errRp * (1.0 + kSRAAOuterSlack)) {
        if (++outer_stall_count >= kSRAAOuterStallWindow) {
            X = X_best;                   // revert to outer-quality best
            grk_sraa.clear();             // restart AA history
            grk_sraa.allow_aa = false;    // re-require plateau before AA fires
            plateau_count = 0; outer_stall_count = 0;
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
            prev_outer_quality = best_errRp;
        }
    } else { outer_stall_count = 0; }
}
```

**Applies identically to `raking.cpp`** (same block, same variables, complements existing
`best_metric_seen`/`W_best` infrastructure).

### Why plateau-gating guarantees correctness

After `kSRAAPlateauWindow=4` outer iterations with improvement < 0.1%, the outer quality
is near the correct fixed point (max_err-optimal). At that point:

- `best_errRp` ≈ correct fixed point quality (e.g., 1.57e-3 for K=9 stepstone)
- If AA escapes to wrong basin (2.12e-3): curr_max = 2.12e-3 > 1.57e-3 × 1.10 = 1.73e-3
  → outer_stall_count increments → after 3 steps → revert to X_best

### Expected performance

| K | Plateau at (outer iters) | AA fires | Acceleration potential |
|---|---|---|---|
| 2 | ~7 outer iters | Early | High (3-5× faster) |
| 3-5 | ~40 outer iters | 80% through budget | Moderate |
| 9 (stepstone) | ~900 outer iters | 87% through budget | Low but CORRECT |

For K=9: correctness guaranteed at cost of minimal acceleration. AA fires near convergence
as a "polishing" step — limited benefit but quality ≤ plain. This is a correct tradeoff.

---

## Calibration of constants

| Constant | Value | Derivation |
|---|---|---|
| `kSRAAPlateauEps` | 1e-3 | For K=9 stepstone (1030 outer iters to converge 1.57e-3), improvement near convergence: Δq/q ≈ (1.65-1.57)/1.57/50 ≈ 1e-3/iter. Plateau when improvement < 0.1%. |
| `kSRAAPlateauWindow` | 4 | 4 consecutive iters below threshold filters transient quality oscillations (overlapping margins cause ≤2-5% per-step fluctuations). |
| `kSRAAOuterSlack` | 0.10 | K=9 AA escape = 35% regression. 10% slack catches escape without false positives from normal transients. |
| `kSRAAOuterStallWindow` | 3 | 3 consecutive regressed outer iters = minimal wasted budget before revert. |

---

## TDD Requirements

### T_sraa_adaptive_K9 (RED with current code, GREEN after Track 2 or Track 1)

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

### T_sraa_raking_K9 (same for raking)

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

### Existing tests must remain GREEN
T_sraa_grk (K=2 greenkhorn+AA faster and at least as good as plain)
T_sraa_rk (K=2 raking+AA)
T_sraa_global (K=4 regression guard)
T_sraa_ldlt_fallback, T_sraa_restart, T5-T8, T_logit_*

---

## Acceptance Criteria

| # | Criterion | Verify |
|---|-----------|--------|
| AC1 | T_sraa_grk GREEN (K=2, quality + speed) | `devtools::test(filter="calibration-solvers")` |
| AC2 | T_sraa_rk GREEN (K=2 raking) | Same |
| AC3 | T_sraa_adaptive_K9 GREEN (K=9 greenkhorn) | Same (skip if parquet absent) |
| AC4 | T_sraa_raking_K9 GREEN (K=9 raking) | Same |
| AC5 | FAIL count = 3 | `devtools::test()` |
| AC6 (benchmark) | greenkhorn+AA iters < plain iters on K=9 | `Rscript benchmarks/stepstone_all_methods.R` |

---

## Implementation Order

1. **Track 2 only first**: Change 5 lines in greenkhorn.cpp. Test AC3/AC4.
   - If both pass: stop here. Track 1 not needed.
   - If AC3 passes but AC4 fails (raking+AA): Track 1 needed for raking only.
   - If both fail: implement Track 1 for both.

2. **Track 1 as needed**: Add 30 lines to each failing solver + 1 field + 4 constants to sraa.hpp.

---

## Out of Scope

- Other methods (logit, ieppa, sinkhorn): do not use sraa_step, not affected
- Exposing plateau/revert parameters to R user (internal tuning)
- Adaptive m window size
