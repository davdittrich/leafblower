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

### Pre-implementation baseline (REQUIRED before touching code)

Run and record baseline values:
```bash
OMP_NUM_THREADS=1 Rscript benchmarks/stepstone_all_methods.R 2>&1 | grep -E "greenkhorn|raking"
```
Record `max_err` for greenkhorn+squarem and raking+squarem. These are the pre-fix baselines
that AC3/AC4 must beat. Without this baseline, a passing AC3 could be coincidental.

### What changes

**`src/greenkhorn.cpp`** — `f_eval_sraa` lambda, replace fixed-sort with adaptive-sort
AND remove the now-dead `order_sraa` vector:

```cpp
// REMOVE these declarations (order_sraa no longer needed):
// std::vector<int> order_sraa;
// if (st.accelerate && K > 0) { ... order_sraa.assign(K, 0); }

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

**New `SRAAState` field + methods:**
```cpp
bool aa_unlocked = false;  // lifecycle gate: set by enable_aa(), cleared by disable_aa()
                           // clear() does NOT reset this — outer loop owns the lifecycle

void enable_aa()  { aa_unlocked = true;  }
void disable_aa() { aa_unlocked = false; }
```

The field and methods are in SRAAState so that (a) naming communicates lifecycle intent,
not a feature flag, and (b) clear() semantics are explicit: history resets, gate does not.

**New constants in `sraa.hpp`:**
```cpp
static constexpr double kSRAAPlateauEps       = 1e-3;  // 0.1%/iter improvement → plateau
static constexpr int    kSRAAPlateauWindow     = 4;    // 4 consecutive plateau iters → AA enabled
static constexpr double kSRAAOuterSlack        = 0.10; // 10% above best_errRp → regression
static constexpr int    kSRAAOuterStallWindow  = 5;    // 5 consecutive regressed outer iters → revert
```

Note: `kSRAAOuterStallWindow = 5` (not 3). Near convergence (where plateau fires), a 3-step
window is too aggressive — normal AA warm-up transients can briefly elevate quality by 5-10%.
5 steps require sustained regression before triggering revert, reducing false positives.

**`sraa_step` change** (one additional condition in step 5):
```cpp
// Not enough history OR outer loop has not unlocked AA → plain step
if (state.count < kSRAAMinCount || !state.aa_unlocked) {
    std::swap(X, state.F_cur);
    return {false, 1, err_plain};
}
```
History still accumulates during plateau-detection phase (steps 3-4 run). AA fires only
when the outer loop calls `grk_sraa.enable_aa()`.

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
        if (plateau_count >= kSRAAPlateauWindow) grk_sraa.enable_aa();
    }
    prev_outer_quality = curr_max;

    // Outer revert: when AA escapes correct basin
    if (grk_sraa.aa_unlocked && curr_max > best_errRp * (1.0 + kSRAAOuterSlack)) {
        if (++outer_stall_count >= kSRAAOuterStallWindow) {
            X = X_best;                   // revert to outer-quality best
            grk_sraa.clear();             // restart AA history; clear() does NOT reset aa_unlocked
            grk_sraa.disable_aa();        // re-require plateau before AA fires
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

**Applies to `raking.cpp`** with one naming difference: raking uses `best_metric_seen`
(not `best_errRp`) and `W_best` (not `X_best`). The revert block substitutes these:
`X = W_best` → `X = grk_sraa_equiv.X_best` does not apply; for raking, `W_best` is the
iterate, so: `W = W_best; rk_sraa.clear(); rk_sraa.disable_aa()`. Do NOT introduce
a parallel X_best in raking — reuse W_best which already tracks the same quantity.

The plateau detection and outer_stall_count variables are local to the SRAA accelerate
block (not stored in SRAAState), so no naming conflict with existing raking infrastructure.

### Why plateau-gating guarantees correctness

After `kSRAAPlateauWindow=4` outer iterations with improvement < 0.1%, the outer quality
is near the correct fixed point (max_err-optimal). At that point:

- `best_errRp` ≈ correct fixed point quality (e.g., 1.57e-3 for K=9 stepstone)
- If AA escapes to wrong basin (2.12e-3): curr_max = 2.12e-3 > 1.57e-3 × 1.10 = 1.73e-3
  → outer_stall_count increments → after kSRAAOuterStallWindow=5 steps → revert to X_best

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
| `kSRAAOuterStallWindow` | 5 | 5 consecutive regressed outer iters before revert. AA fires near convergence where normal transients can briefly elevate quality 5-10%; 5 steps require sustained (not transient) regression. |

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

### T_sraa_plateau_gate (Track 1 only — verify AA is disabled before plateau)

```r
test_that("T_sraa_plateau_gate: greenkhorn+AA early iters use plain-only steps before plateau", {
  # With kSRAAPlateauWindow=4, AA must not fire in first 4 outer iterations.
  # Proxy: with max_iterations=4 outer iters × K, iters_aa must equal iters_plain
  # (all plain steps: K*1 each, no AA which gives K*2).
  set.seed(99); n <- 2000L
  df  <- data.frame(x=factor(sample(letters[1:3],n,TRUE)),
                    y=factor(sample(c("M","F"),n,TRUE)))
  tgt <- list(x=c(a=0.3,b=0.4,c=0.3), y=c(M=0.5,F=0.5))
  K_exp <- 2L
  r_aa    <- suppressWarnings(harvest(df,tgt,method="greenkhorn",accelerate=TRUE,
                                       max_iterations=4L,attach_weights=FALSE))
  r_plain <- suppressWarnings(harvest(df,tgt,method="greenkhorn",accelerate=FALSE,
                                       max_iterations=4L,attach_weights=FALSE))
  # Before plateau: AA disabled → same iter count as plain (K*1 per outer iter, not K*2)
  expect_equal(attr(r_aa,"result")$iterations, attr(r_plain,"result")$iterations,
    label="Before plateau, AA must not fire — iterations must match plain")
})
```

### T_sraa_outer_revert (Track 1 only — verify revert recovers quality)

```r
test_that("T_sraa_outer_revert: greenkhorn+AA K=4 quality recovers to <= plain after potential escape", {
  # Same K=4 problem from T_sraa_global. With plateau-gating + outer revert,
  # SRAA must not regress below plain even if AA escapes to a wrong basin.
  set.seed(42); n <- 8000L
  df <- data.frame(gender=factor(sample(c("M","F"),n,TRUE)),
                   time=factor(sample(1:3,n,TRUE)),
                   age=factor(sample(1:4,n,TRUE)))
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
  me_aa    <- attr(r_aa,   "result")$max_error
  me_plain <- attr(r_plain,"result")$max_error
  expect_lte(me_aa, me_plain * 1.001 + 1e-10,
    label=sprintf("K=6 revert: AA (%.2e) must not exceed plain (%.2e)", me_aa, me_plain))
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
| AC1 | T_sraa_grk GREEN (K=2, quality ≤ plain AND iters_aa < iters_plain) | `devtools::test(filter="calibration-solvers")` |
| AC2 | T_sraa_rk GREEN (K=2 raking quality ≤ plain) | Same |
| AC3 | T_sraa_adaptive_K9 GREEN (K=9 greenkhorn+AA ≤ plain) | Same (skip if parquet absent in CI) |
| AC4 | T_sraa_raking_K9 GREEN (K=9 raking+AA ≤ plain) | Same (skip if parquet absent) |
| AC5 | FAIL count = 3 (pre-existing: test-ieppa-nonuniform-d.R:28, :29, test-sor.R:18) | `devtools::test()` |
| AC6 (Track 1 only) | T_sraa_plateau_gate GREEN | `devtools::test(filter="calibration-solvers")` |
| AC7 (Track 1 only) | T_sraa_outer_revert GREEN | Same |
| AC8 (benchmark) | greenkhorn+AA iters < plain iters on K=9 | `Rscript benchmarks/stepstone_all_methods.R` |

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
