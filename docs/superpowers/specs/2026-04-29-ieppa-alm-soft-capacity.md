# iEPPA Soft: ALM Capacity Enforcement — Design Spec

**Date**: 2026-04-29
**Status**: Pending design review
**Tickets**: leafblower-scii (T3 bug), leafblower-bgx6 (alm_mu decision)
**Files**: `src/ieppa.cpp`, `src/ieppa.hpp`, `src/types.hpp`, `src/cell_table.cpp`,
          `src/r_bridge.cpp`, `src/c_api.cpp`, `R/harvest.R`,
          `tests/testthat/test-calibration-solvers.R`

---

## Problem

`method="ieppa_soft"` was specced (2026-04-27) as ADMM with z-update soft capacity enforcement. T3 acceptance test requires: `ieppa_soft max_err < ieppa max_err` on tight-bounds problems.

**Current state**: T3 fails. ieppa_soft produces bit-identical weights to ieppa hard-clamp:

```
ieppa     max_err: 0.0273525721455
ieppa_soft max_err: 0.0273525721455
weights identical: TRUE
```

**Root cause**: The existing ADMM z-update at ieppa.cpp:807 still hard-clamps internally (`z = clamp(X_tilde + u[c], L, U)` — z is always in [L, U]). The dual u[c] only shifts the clamp window. It cannot let X[c] temporarily violate bounds to escape local minima, so on initial iteration with u[c]=0 the result is identical to hard clamp. By the time u[c] accumulates meaningful signal (10+ iterations), the solver has already converged on the bad fixed point.

**User intent**: Augmented Lagrangian Method (ALM) — soft capacity penalty during optimization, exact bounds in final output. The hard cap should be **removed** during the iterative phase and **re-applied** as a single projection at the end.

---

## Design

### Decision Summary (from brainstorming 2026-04-29)

| # | Question | Decision |
|---|----------|----------|
| 1 | Scope | (1a) ALM only in `ieppa_soft`; `ieppa` keeps hard-clamp (backward compat) |
| 2 | Default `alm_mu` | (2c) Auto-compute as `M_cell/n`; user override via R param |
| 3 | Dynamic schedule | (3c) Tang-η across homotopy levels + adaptive μ growth on persistent violation |
| 4 | Final projection | (4a) Always hard-clamp after last iteration; renormalize once |
| 5 | u[c] reset policy | (5b) Reset only on homotopy level transitions; preserve across overflow fallback |
| 6 | Convergence criterion | (6c) Auto-increase μ when primal infeasibility persists ≥ K_persist iterations |
| 7 | R-facing API name | `capacity_penalty` (NULL = auto, positive scalar = manual) |

### Algorithm

**Per-cell ALM update** (replaces hard-clamp at ieppa.cpp:804-812):

```cpp
// ieppa_soft ALM capacity update — replaces P1.1 hard clamp
double z = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);   // proj target
double rho = st.alm_mu * X_tilde_c;                        // cell-normalized strength
X[c]  = (X_tilde_c + rho * z + u[c]) / (1.0 + rho);        // soft blend
u[c] += st.alm_mu * (X[c] - z);                            // dual accumulation
W[c]  = X[c] / X_tilde_c;
X_cur[c] = X[c];
```

**Mathematical derivation** (justifies blend formula): Solve

```
argmin_{X[c]>0}  KL(X[c] | X_tilde[c]) + λ[c]·(X[c] - z[c]) + (μ/2)·(X[c] - z[c])²
```

Setting gradient to zero with second-order Taylor at `X[c] ≈ X_tilde[c]` gives:

```
X[c] ≈ (X_tilde[c] + ρ·z[c] - λ[c]·X_tilde[c]/μ) / (1 + ρ),  ρ = μ·X_tilde[c]
```

With scaled dual `u[c] = -λ[c]·X_tilde[c]/μ` we get the closed form above. When `ρ → 0`: `X[c] = X_tilde[c]` (no constraint). When `ρ → ∞`: `X[c] = z[c]` (hard clamp recovered). When `ρ ≈ 1`: 50/50 blend, X[c] can lie outside [L, U] briefly.

**Dual update**: standard ALM ascent on Lagrangian, scaled. As `X[c] - z[c] → 0`: `u[c]` stops growing → at convergence X[c] is feasible.

**Final projection** (after last solver iteration, before return):

```cpp
double total = 0.0;
for (int c = 0; c < ct.M_cell; c++) {
    X[c] = std::clamp(X[c], L_cell[c], U_cell[c]);
    total += X[c];
}
const double rescale = static_cast<double>(st.n) / total;
for (int c = 0; c < ct.M_cell; c++) {
    X[c] *= rescale;
    // Re-clamp tiny renormalization drift; max drift bounded by |1 - rescale|·U_cell.
    X[c] = std::clamp(X[c], L_cell[c], U_cell[c]);
}
```

This guarantees `L_cell ≤ X[c] ≤ U_cell` exactly (modulo IEEE 754 rounding). Note the second clamp may break sum-to-n by O(eps·n); accept this drift as the cost of strict bounds adherence.

### Dynamic μ Schedule (combined Tang-η + adaptive growth)

**Macro (Tang-η, existing at ieppa.cpp:381-391)** — unchanged:

```cpp
if (st.eta_schedule.mode == EtaScheduleMode::TANG_DYNAMIC && N_levels > 1) {
    const double scaled_frac = std::pow(frac, st.eta_schedule.schedule_power);
    const double eta_i = st.eta_schedule.eta_start *
        std::pow(st.eta_schedule.eta_end / st.eta_schedule.eta_start, scaled_frac);
    if (alm_mu_base > 0.0) {
        st.alm_mu = eta_i * alm_mu_base;   // <-- already wired
    }
}
```

**Adaptive (per-iteration violation tracking)** — new:

```cpp
// Track primal infeasibility across iterations
constexpr int    kAlmPersistenceThreshold = 5;
constexpr double kAlmGrowthFactor         = 2.0;   // doubles when triggered
constexpr double kAlmMaxScale             = 1000.0; // cap at 1000x base
constexpr double kAlmTolPrimalScale       = 0.01;   // 1% of max_weight slack

// At each capacity block iteration:
double max_violation = 0.0;
for (int c = 0; c < ct.M_cell; c++) {
    double over  = X[c] - U_cell[c];
    double under = L_cell[c] - X[c];
    max_violation = std::max(max_violation, std::max(over, under));
}
const double tol_primal = kAlmTolPrimalScale * st.max_weight * mean_n_per_cell;

if (max_violation > tol_primal) {
    alm_violation_streak++;
    if (alm_violation_streak >= kAlmPersistenceThreshold &&
        st.alm_mu < alm_mu_base * kAlmMaxScale) {
        st.alm_mu *= kAlmGrowthFactor;
        alm_violation_streak = 0;  // reset after growth
    }
} else {
    alm_violation_streak = 0;
}
```

### Auto-Compute `alm_mu`

`cell_table.cpp build_cell_table` populates `out.alm_mu_auto = M_cell / n`. Read at solver entry:

```cpp
// In r_bridge.cpp / c_api.cpp before solver dispatch:
if (capacity_penalty_param == nullptr || *capacity_penalty_param <= 0.0) {
    p.alm_mu = ct.alm_mu_auto;       // computed by cell_table
} else {
    p.alm_mu = *capacity_penalty_param;
}
```

### `u[c]` Reset Policy

**Reset on**: homotopy level transitions ONLY.
**Preserve through**: overflow → log-path fallback (lines 692, 708).

```cpp
// At start of each homotopy level (NEW location, after level setup):
if (st.use_admm_capacity) std::fill(u.begin(), u.end(), 0.0);

// Remove existing resets at lines 692 and 708 (overflow fallback paths).
```

**Caveat**: u[c] was accumulated against linear-path X_tilde. After fallback to log-path, X_tilde scale may differ. Document this as a known limitation; rescaling u[c] across paths is out of scope.

### R API

```r
# In R/harvest.R:
harvest(df, target,
        method = "ieppa_soft",
        max_weight = 5,
        capacity_penalty = NULL,    # NULL = auto (M_cell/n); positive scalar = manual
        ...)
```

Propagation: harvest.R → `.Call("C_rk_calibrate", ...)` → `r_bridge.cpp` → `rk_params_t.alm_mu` (or new field `capacity_penalty`) → `CalibState.alm_mu` → ieppa.cpp.

**Validation in R layer**:

```r
if (!is.null(capacity_penalty)) {
  if (!is.numeric(capacity_penalty) || length(capacity_penalty) != 1L ||
      !is.finite(capacity_penalty) || capacity_penalty <= 0) {
    stop("capacity_penalty must be NULL or a positive finite scalar")
  }
}
```

---

## Files to Modify

| File | Change |
|------|--------|
| `src/cell_table.cpp` | `build_cell_table`: compute `out.alm_mu_auto = M_cell / n` |
| `src/cell_table.hpp` | Add `double alm_mu_auto` field to `CellTable` struct |
| `src/types.hpp` | (no change — `alm_mu` already on CalibState) |
| `src/ieppa.cpp` | Replace P1.1 capacity block (lines 800-816) with ALM update; add adaptive μ growth; remove u[] resets at 692, 708; add final projection |
| `src/ieppa.hpp` | (no change) |
| `src/r_bridge.cpp` | Read `capacity_penalty` SEXP; route to `p.alm_mu`; auto-fill from `ct.alm_mu_auto` if NULL |
| `src/c_api.cpp` | Same routing for direct C callers |
| `R/harvest.R` | Add `capacity_penalty = NULL` param; validate; pass through `.Call` |
| `tests/testthat/test-calibration-solvers.R` | T3 strengthened; T4–T7 new (see TDD) |

---

## TDD (RED tests before implementation)

### T3 (existing, currently failing — must turn GREEN):
ieppa_soft achieves strictly lower max_err than ieppa on tight-bounds. Already in spec 2026-04-27.

### T4 — capacity_penalty parameter accepted

```r
test_that("ieppa_soft: capacity_penalty=NULL routes to auto-computed value", {
  set.seed(4); n <- 1000L
  df <- data.frame(v=factor(sample(letters[1:3], n, TRUE)))
  tgt <- list(v=c(a=0.4, b=0.4, c=0.2))
  r1 <- harvest(df, tgt, method="ieppa_soft", capacity_penalty=NULL,
                attach_weights=FALSE)
  r2 <- harvest(df, tgt, method="ieppa_soft",
                attach_weights=FALSE)
  # NULL and unspecified should produce identical results
  expect_equal(as.numeric(r1), as.numeric(r2), tolerance=1e-12)
})

test_that("ieppa_soft: capacity_penalty rejects invalid input", {
  df <- data.frame(v=factor(c("a","a","b")))
  tgt <- list(v=c(a=0.5, b=0.5))
  expect_error(harvest(df, tgt, method="ieppa_soft", capacity_penalty=-1),
               "positive finite scalar")
  expect_error(harvest(df, tgt, method="ieppa_soft", capacity_penalty=0),
               "positive finite scalar")
  expect_error(harvest(df, tgt, method="ieppa_soft", capacity_penalty=c(1,2)),
               "positive finite scalar")
})
```

### T5 — Final bounds adherence (exact)

```r
test_that("ieppa_soft: final weights respect bounds exactly", {
  set.seed(5); n <- 5000L
  df <- data.frame(v=factor(sample(5, n, TRUE)))
  tgt <- list(v=setNames(c(0.4,0.3,0.15,0.1,0.05), as.character(1:5)))
  r <- harvest(df, tgt, method="ieppa_soft",
               max_weight=1.8, min_weight=0.1,
               max_iterations=300, attach_weights=FALSE)
  w <- as.numeric(r)
  # Strict inequality with no floating-point slack — final projection guarantees this
  expect_true(max(w) <= 1.8)
  expect_true(min(w) >= 0.1)
})
```

### T6 — ALM converges with auto μ on stepstone-like problem

```r
test_that("ieppa_soft: auto capacity_penalty achieves max_err <= ieppa max_err on stepstone", {
  skip_on_cran()
  if (!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet")) skip("no stepstone fixture")
  df <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet"); df$uuid <- NULL
  tgt <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
  tgt <- lapply(tgt, function(t) { t <- unlist(t); t / sum(t) })
  for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])

  r_hard <- harvest(df, tgt, method="ieppa",
                    max_weight=5, max_iterations=500, attach_weights=FALSE)
  r_soft <- harvest(df, tgt, method="ieppa_soft",
                    max_weight=5, max_iterations=500, attach_weights=FALSE)
  me_hard <- attr(r_hard, "result")$max_error
  me_soft <- attr(r_soft, "result")$max_error
  expect_lte(me_soft, me_hard)
})
```

### T7 — Adaptive μ growth fires on persistent infeasibility

```r
test_that("ieppa_soft: adaptive growth scales alm_mu when primal violation persists", {
  set.seed(7); n <- 2000L
  df <- data.frame(v=factor(sample(c("X","Y"), n, TRUE, prob=c(.3,.7))))
  tgt <- list(v=c(X=0.99, Y=0.01))   # adversarial: target 99% in minority cat
  r <- harvest(df, tgt, method="ieppa_soft",
               capacity_penalty=1e-6,    # start with tiny penalty
               max_weight=10, max_iterations=300, attach_weights=FALSE,
               verbose=2L)
  # Solver must converge: status in {OK, BUDGET, STALL} but NOT INFEAS
  s <- attr(r, "result")$status
  expect_true(s %in% c(0L, 4L, 5L))
  # Final weights respect bounds (T5-style)
  w <- as.numeric(r)
  expect_true(max(w) <= 10 + 1e-9)
})
```

---

## Acceptance Criteria

1. **T3 GREEN**: `me_soft < me_hard` on tight-bounds problem.
2. **T4 GREEN**: `capacity_penalty=NULL` routes to auto, validation rejects bad input.
3. **T5 GREEN**: Final weights satisfy bounds exactly (no IEEE drift outside [L, U]).
4. **T6 GREEN**: On stepstone, `ieppa_soft max_err ≤ ieppa max_err`.
5. **T7 GREEN**: Adaptive growth saves convergence on adversarial small-μ start.
6. **Backward compat**: `method="ieppa"` produces identical results before/after this change.
7. **No new failures**: `devtools::test()` regression count ≤ 3 (current baseline).
8. **Stepstone benchmark**: `ieppa_soft` matches or beats `ieppa` on max_err and weight KL.

---

## Risk / Limitations

1. **u[c] scale mismatch on log-path fallback**: dual is in linear-X scale; log-path uses log-X. Documented; rescaling out of scope. Worst case: the dual signal is wrong magnitude for one or two iterations until log-path stabilizes.

2. **Final clamp may break sum-to-n by O(eps·n)**: typical drift < 1e-10 · n. Acceptable for survey calibration where bounds adherence is the harder requirement.

3. **Adaptive μ growth caps at 1000× base**: pathological inputs (e.g., target=0.99 with 5% of obs available) may not converge. Solver returns BUDGET/STALL with best-iter weights. T7 documents this.

4. **Fixed `kAlmGrowthFactor=2.0`**: not exposed to user. If 2.0 is too aggressive (oscillation) or too conservative (slow growth), tune based on empirical results.

5. **Closed-form blend assumes ρ=μ·X_tilde balances KL Hessian**: derivation uses second-order Taylor, valid when X_tilde[c] is not near zero. For cells with X_tilde[c]→0 (structurally near-empty), the blend reduces to X[c]≈X_tilde[c] which is correct (no penalty applied to a zero cell).

---

## Out of Scope

- ALM for raking solver (separate spec, different geometry)
- Per-cell adaptive μ (single global μ for now)
- Removing `alm_lambda` field from CalibState (cleanup ticket if needed)
- Replacing ADMM in greg/sinkhorn (those use different convergence theory)
