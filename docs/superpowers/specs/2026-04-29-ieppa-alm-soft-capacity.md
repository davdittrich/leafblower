# iEPPA Soft: ALM Capacity Enforcement — Design Spec (rev 2)

**Date**: 2026-04-29
**Status**: Pending design review
**Tickets**: leafblower-scii (T3 bug), leafblower-bgx6 (alm_mu decision)
**Files**: `src/ieppa.cpp`, `src/ieppa.hpp`, `src/types.hpp`, `src/cell_table.cpp`,
          `src/cell_table.hpp`, `src/r_bridge.cpp`, `src/c_api.cpp`,
          `R/harvest.R`, `man/harvest.Rd` (auto-regenerated),
          `tests/testthat/test-calibration-solvers.R`

---

## Problem

`method="ieppa_soft"` was specced (2026-04-27) as ADMM with z-update soft capacity enforcement. T3 acceptance test requires `ieppa_soft max_err < ieppa max_err` on tight-bounds. T3 is failing — ieppa_soft produces bit-identical weights to ieppa hard-clamp:

```
ieppa     max_err: 0.0273525721455
ieppa_soft max_err: 0.0273525721455
weights identical: TRUE
```

**Root cause**: existing ADMM z-update at ieppa.cpp:807 still hard-clamps internally (`z = clamp(X̃ + u[c], L, U)`). The dual u[c] only shifts the clamp window. On first iteration with u[c]=0, result is identical to hard clamp. By the time u[c] accumulates signal, the solver has converged on the bad fixed point.

**User intent**: Augmented Lagrangian (ALM) — soft penalty during optimization (X[c] can briefly violate bounds), exact bounds in final output via single projection.

**When to use ieppa_soft over ieppa**: tight-bounds problems where `ieppa` exits with non-trivial max_err and the user wants the calibration to escape the local minimum at the constraint boundary. For loose bounds (max_weight ≥ 5× typical weight), `ieppa` is faster and gives equivalent results. Documented in harvest.Rd `@details`.

---

## Design

### Decision Summary

| # | Question | Decision |
|---|----------|----------|
| 1 | Scope | (1a) ALM only in `ieppa_soft`; `ieppa` keeps hard-clamp |
| 2 | Default `capacity_penalty` | (2c) Auto = `M_cell/n`; user override via R param |
| 3 | Dynamic schedule | (3c) Tang-η × persistent-adaptive growth (composed correctly) |
| 4 | Final projection | (4a) Always hard-clamp + bounded-iteration rescale |
| 5 | u[c] reset policy | **REVISED** (5b'): reset on level transitions AND on numerical fallback (linear→log) |
| 6 | Convergence criterion | (6c) Auto-grow μ on persistent infeasibility |
| 7 | R-facing API name | `capacity_penalty` (NULL = auto, positive scalar = manual) |
| 8 | Method name | Keep `ieppa_soft` (backward compat with 2026-04-27 spec) |

**Decision 5 changed in rev 2**: prior decision preserved u[c] across overflow fallback. Architect + Security review showed the linear-X scale of u[c] is invalid in log-path X̃, producing silent NaN or wrong-magnitude penalty. Reset on fallback accepts loss of convergence progress (the fallback is itself a reset event).

### Field naming (avoid lbfgsb collision)

`alm_lambda` and `alm_mu` in `CalibState` are **actively used by lbfgsb_solver.cpp** (lines 74-87, 378-392, 510-524) for the sum-to-n ALM. They must NOT be reused for ieppa's capacity ALM.

**New CalibState fields** (in `src/types.hpp`):
```cpp
double capacity_mu       = 0.0;  // ieppa_soft ALM penalty coefficient (0 = inactive)
// (per-cell dual `lambda_cell[c]` lives in solver-local std::vector, not CalibState)
```

The existing `use_admm_capacity` field is kept and still gates the ALM block (renamed semantically — see Implementation §2).

### Algorithm

**Notation**: For cell `c`, let `X̃[c]` = `X_tilde[c]` = unconstrained iterate from Sinkhorn sweeps. Let `λ[c]` = per-cell ALM dual (raw, no scaling). Let `μ` = `st.capacity_mu`.

**Per-cell ALM update** (replaces hard-clamp at ieppa.cpp:804-812):

```cpp
// ieppa_soft ALM capacity update
double z = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);   // proj target
double rho = st.capacity_mu * X_tilde_c;                  // dimensionless ratio (KL_Hess vs penalty_Hess)

// Primal (closed-form Newton step from X = X_tilde):
double X_alm = X_tilde_c * (1.0 - lambda_cell[c] + st.capacity_mu * z) / (1.0 + rho);

// NaN/Inf safety (decision: revert to hard clamp on numerical failure):
if (!std::isfinite(X_alm) || X_alm <= 0.0) {
    X_alm = z;  // fall back to hard-clamp value
}

X[c] = X_alm;
W[c] = X_alm / X_tilde_c;
X_cur[c] = X_alm;

// Dual ascent (raw, standard ALM):
lambda_cell[c] += st.capacity_mu * (X_alm - z);
// Bound dual to prevent runaway (10× max weight is plenty of dual range):
const double lambda_cap = 10.0 * st.max_weight;
lambda_cell[c] = std::clamp(lambda_cell[c], -lambda_cap, lambda_cap);
```

**Mathematical derivation** (replaces fuzzy reference in rev 1):

Per-cell objective: `f(X) = KL(X | X̃) + λ·(X − z) + (μ/2)·(X − z)²`

Substituting `s = X/X̃` (so `s=1` at `X = X̃`):
- `KL(X | X̃) = X̃·(s·log(s) − s + 1)` → `d²/ds² = X̃/s = X̃` at s=1
- Penalty `(μ/2)·(s·X̃ − z)²` → `d²/ds² = μ·X̃²` at s=1
- Total Hessian at s=1: `X̃·(1 + μ·X̃) = X̃·(1 + ρ)` where ρ = μ·X̃ (dimensionless)
- Gradient at s=1: `X̃·log(1) + X̃·(λ + μ·(X̃ − z)) = X̃·(λ + μ·(X̃ − z))`
- Newton step: `Δs = −(λ + μ·(X̃ − z)) / (1 + ρ)`
- `s_new = (1 + μ·X̃ − λ − μ·X̃ + μ·z) / (1 + ρ) = (1 − λ + μ·z) / (1 + ρ)`
- `X_new = X̃ · s_new = X̃·(1 − λ + μ·z) / (1 + ρ)`

**Limits**:
- `ρ → 0` (μ → 0): `X = X̃·(1 − λ)` — purely dual-driven, no constraint enforcement
- `ρ → ∞` (μ → ∞): `X → z` — recovers hard clamp (Newton dominated by penalty Hessian)
- `ρ ≈ 1`: balanced — KL Hessian ≈ penalty Hessian; step blends X̃ and z roughly equally
- `λ = 0`: `X = X̃·(1 + μ·z)/(1 + μ·X̃)` — pure penalty, no dual signal

**Standard ALM dual ascent**: `λ_new = λ + μ·(X_new − z)`. As `X − z → 0` at convergence, `λ` stops growing. λ has units mass⁻¹ (matches X̃⁻¹).

**Why `auto μ = M_cell/n`**: The KL Hessian per cell is `1/X[c] ≈ M_cell/n` (for X[c] near mean cell mass `n/M_cell`). To make the penalty Hessian comparable: `μ = M_cell/n` gives `ρ ≈ 1` — balanced 50/50 blend, neither hard-clamp-degenerate (ρ→∞) nor unconstrained (ρ→0). Justified by Boyd et al. (2011) §3.4 "Optimal step size selection for ADMM": for separable objectives with Hessian H, optimal μ ≈ √(λ_min·λ_max) of H. For our problem, λ_min ≈ 1/U[c] and λ_max ≈ 1/L[c] (or 1/X̃ for L=0); the geometric mean is bounded above by 1/X̃ ≈ M_cell/n.

### Final Projection (exact bounds, bounded sum drift)

```cpp
// Apply after solver loop exits, before returning result:
constexpr int kMaxRescaleIters = 3;
constexpr double kRescaleTol = 1e-12;

for (int c = 0; c < ct.M_cell; c++) {
    X[c] = std::clamp(X[c], L_cell[c], U_cell[c]);
}

// Iterate clamp+rescale until both stabilize OR max iters reached.
// Each iteration is monotone: rescale shifts weights, second clamp pulls back.
// Stable fixed point exists iff sum(L_cell) ≤ n ≤ sum(U_cell) (always true if upstream validation passed).
double prev_total = 0.0;
for (int iter = 0; iter < kMaxRescaleIters; iter++) {
    double total = 0.0;
    for (int c = 0; c < ct.M_cell; c++) total += X[c];
    if (std::abs(total - prev_total) < kRescaleTol * st.n) break;
    prev_total = total;
    if (total <= 0.0) break;  // pathological — accept output as-is
    const double rescale = static_cast<double>(st.n) / total;
    for (int c = 0; c < ct.M_cell; c++) {
        X[c] *= rescale;
        X[c] = std::clamp(X[c], L_cell[c], U_cell[c]);
    }
}
// Post-condition: max(L) ≤ X[c] ≤ max(U) exactly. Sum may differ from n by up to
// kRescaleTol * n + final_clamp_drift. Worst-case clamp drift is bounded by
// kMaxRescaleIters * sum_clamped_cells * U_cell. Logged at verbose >= 1 if drift > 1e-6 * n.
{
    double final_total = 0.0;
    for (int c = 0; c < ct.M_cell; c++) final_total += X[c];
    const double sum_drift = std::abs(final_total - static_cast<double>(st.n));
    if (sum_drift > 1e-6 * st.n && st.verbose >= 1) {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
            "[ieppa_soft] final projection sum drift = %.2e (target n = %d); "
            "bounds enforced exactly, sum-to-n approximate",
            sum_drift, st.n);
        st.log(msg);
    }
}
```

**Drift bound**: After clamp+rescale convergence, weights satisfy `L ≤ X[c] ≤ U` exactly (last operation is clamp). `sum(X)` may differ from `n` by up to the final-iteration rescale drift (bounded by `kRescaleTol * n` per iteration, plus single-iteration clamp pull-back of at most `sum_clamped_cells * eps`). Worst case ≈ `n * kRescaleTol * kMaxRescaleIters = 3e-12 * n` — well within survey calibration tolerance. Pathological inputs (sum(L) > n, infeasible) emit warning and accept as-is.

### Dynamic μ Schedule (composed correctly)

Persistent state across homotopy levels:
```cpp
double capacity_mu_base = st.capacity_mu;       // user-supplied (or auto from cell_table)
double capacity_mu_adaptive = capacity_mu_base; // adaptive scale, persists across levels
int    alm_violation_streak = 0;                // resets on level transition AND on decay
```

At each homotopy level entry (line 386-391, replace existing alm_mu logic):
```cpp
if (st.eta_schedule.mode == EtaScheduleMode::TANG_DYNAMIC && N_levels > 1) {
    const double scaled_frac = std::pow(frac, st.eta_schedule.schedule_power);
    const double eta_i = st.eta_schedule.eta_start *
        std::pow(st.eta_schedule.eta_end / st.eta_schedule.eta_start, scaled_frac);
    res.eta_final = eta_i;
    st.capacity_mu = eta_i * capacity_mu_adaptive;  // Tang-η × adaptive (compose, don't overwrite)
} else {
    st.capacity_mu = capacity_mu_adaptive;          // No Tang-η: pure adaptive
}
alm_violation_streak = 0;  // Reset streak on level transition
// Reset per-cell dual on level transition (was decision 5b, unchanged):
if (st.use_admm_capacity) std::fill(lambda_cell.begin(), lambda_cell.end(), 0.0);
```

Per-iteration adaptive growth (after capacity block):
```cpp
constexpr int    kAlmPersistenceThreshold = 5;
constexpr double kAlmGrowthFactor         = 2.0;
constexpr double kAlmMaxScale             = 1000.0;
constexpr double kAlmTolPrimalRel         = 0.01;  // 1% of max_weight slack

double max_violation = 0.0;
for (int c = 0; c < ct.M_cell; c++) {
    double over  = X[c] - U_cell[c];
    double under = L_cell[c] - X[c];
    double v = std::max(over, under);
    if (std::isfinite(v)) max_violation = std::max(max_violation, v);
}
const double tol_primal = kAlmTolPrimalRel * st.max_weight *
                          (static_cast<double>(st.n) / std::max(1, ct.M_cell));

if (max_violation > tol_primal) {
    alm_violation_streak++;
    if (alm_violation_streak >= kAlmPersistenceThreshold &&
        capacity_mu_adaptive < capacity_mu_base * kAlmMaxScale) {
        capacity_mu_adaptive *= kAlmGrowthFactor;
        st.capacity_mu = (Tang_active ? eta_i_current : 1.0) * capacity_mu_adaptive;
        alm_violation_streak = 0;
        if (st.verbose >= 2) {
            char msg[128];
            std::snprintf(msg, sizeof(msg),
                "[ieppa_soft] alm growth: capacity_mu_adaptive *= %.1f → %.4e",
                kAlmGrowthFactor, capacity_mu_adaptive);
            st.log(msg);
        }
    }
} else {
    alm_violation_streak = 0;
}
```

### u[c] Reset Policy (rev 2)

**Reset lambda_cell[c] on**: 
1. Homotopy level transitions (every level entry)
2. Linear→log path fallback (existing resets at lines 692, 708 — KEPT)
3. Solver entry (initialization)

**Preserve through**: nothing — accept loss of dual signal at fallback as the cost of correctness.

Rationale: linear-path X̃ ≠ log-path X̃ in scale; carrying linear-scale dual into log-path produces NaN or wrong-magnitude penalty (Architect #6, Security CS-2).

### Auto-Compute `capacity_mu`

In `src/cell_table.hpp`:
```cpp
struct CellTable {
    // ... existing fields ...
    double capacity_mu_auto = 0.0;  // auto-computed default for ieppa_soft ALM penalty
};
```

In `src/cell_table.cpp build_cell_table`, after M_cell is known and n is known:
```cpp
out.capacity_mu_auto = (n > 0 && out.M_cell > 0)
                       ? static_cast<double>(out.M_cell) / static_cast<double>(n)
                       : 1.0;  // fallback to safe default if degenerate
```

In `src/r_bridge.cpp` and `src/c_api.cpp`, after cell_table built and before solver dispatch:
```cpp
// Set capacity_mu from R param OR auto-compute:
if (Rf_isNull(capacity_penalty_sexp) || LENGTH(capacity_penalty_sexp) == 0) {
    p.capacity_mu = ct.capacity_mu_auto;
} else {
    if (LENGTH(capacity_penalty_sexp) != 1) {
        Rf_error("leafblower: capacity_penalty must be a single positive finite scalar");
    }
    const double v = REAL(capacity_penalty_sexp)[0];
    if (!R_FINITE(v) || v <= 0.0) {
        Rf_error("leafblower: capacity_penalty must be a positive finite scalar (NULL = auto = M_cell/n)");
    }
    if (v > 1e15) {
        Rf_error("leafblower: capacity_penalty=%.2e exceeds safe range (max 1e15)", v);
    }
    if (v < 1e-15) {
        Rf_warning("leafblower: capacity_penalty=%.2e is below recommended range; constraint enforcement may be ineffective", v);
    }
    p.capacity_mu = v;
}
```

### R API

```r
#' @param capacity_penalty Numeric, controls how strongly capacity bounds are
#'   enforced during ALM optimization in \code{method="ieppa_soft"}. Use
#'   \code{NULL} (default) for auto-computed value (M_cell / n) which gives a
#'   balanced 50/50 blend between unconstrained KL minimization and hard-clamp
#'   projection. Larger values force tighter constraint adherence at the cost
#'   of slower margin convergence; smaller values allow more constraint
#'   violation during optimization but are pulled back to the feasible set by
#'   the final projection. Ignored for methods other than \code{"ieppa_soft"}.
harvest <- function(data, target,
                    method = "auto",
                    max_weight = 5,
                    capacity_penalty = NULL,
                    ...)
```

Validation in R layer (defense-in-depth on top of C-layer validation above):
```r
if (!is.null(capacity_penalty)) {
  if (!is.numeric(capacity_penalty) || length(capacity_penalty) != 1L ||
      !is.finite(capacity_penalty) || capacity_penalty <= 0) {
    stop("capacity_penalty must be NULL (auto) or a positive finite scalar; ",
         "got: ", deparse(capacity_penalty))
  }
}
```

`.Call` argument insertion: `capacity_penalty` added as new positional argument
after `start_weights` (slot 8 → all subsequent slots shift by 1). The C bridge
unpacks via:
```c
SEXP capacity_penalty_sexp = CADR(...);  // exact slot determined by .Call signature
```
**Implementation note**: this requires updating ALL `.Call("C_rk_calibrate", ...)`
call sites — currently at `R/harvest.R:241-276` and any test or benchmark that
calls `.Call` directly. Grep verifies: `grep -rn 'C_rk_calibrate' R/ tests/ benchmarks/`.

### ALM Diagnostics Exposed

Add to `rk_result_t`:
```c
double alm_capacity_mu_final;   /* final capacity_mu after adaptive scaling; 0 if not ieppa_soft */
int    alm_n_growth_events;     /* count of adaptive growth fires; 0 if not ieppa_soft */
double alm_max_dual_norm;       /* max |lambda_cell[c]| at solver exit; for monitoring */
```

Exposed in `attr(r, "result")` so users can monitor:
- `alm_capacity_mu_final / capacity_penalty_input ≈ 1.0` → adaptive growth not needed
- `alm_capacity_mu_final / capacity_penalty_input ≥ 1000` → hit cap, increase capacity_penalty manually
- `alm_max_dual_norm` near `lambda_cap` → constraint very binding, may indicate infeasibility

### Files to Modify (with implementation order)

| Order | File | Change |
|-------|------|--------|
| 1 | `src/types.hpp` | Add `double capacity_mu = 0.0;` to CalibState |
| 2 | `src/cell_table.hpp` | Add `double capacity_mu_auto = 0.0;` to CellTable struct |
| 3 | `src/cell_table.cpp` | `build_cell_table`: compute `capacity_mu_auto = M_cell/n` (with n>0 guard) |
| 4 | `src/leafblower.h` | Add 3 `alm_*` fields to `rk_result_t`; bump `EXPECTED_RK_RESULT_BYTES`; update tripwire |
| 5 | `src/r_bridge.cpp` | Read `capacity_penalty` SEXP with proper Rf_isNull/LENGTH guards; expose 3 ALM diagnostics in result |
| 6 | `src/c_api.cpp` | Same routing; document validation contract for direct C callers (R-layer guards do not apply) |
| 7 | `src/ieppa.cpp` | Replace P1.1 capacity block with ALM update; add lambda_cell vector; add adaptive growth; update u[] reset policy; add final projection |
| 8 | `src/ieppa.hpp` | (no changes needed — IEPPAResult fields covered by rk_result_t additions) |
| 9 | `R/harvest.R` | Add `capacity_penalty = NULL` param + validation; pass through `.Call`; update `@param` Roxygen |
| 10 | `man/harvest.Rd` | Auto-regenerated by `devtools::document()` from harvest.R |
| 11 | `tests/testthat/test-calibration-solvers.R` | T3 strengthened; T4–T8 new |

**Compile order matters**: 1→2→3→4 must compile clean before 5/6 (downstream consumers). Then 7 before 8/9. Then 11 last.

---

## TDD (RED tests before implementation)

### T3 — Strengthened (currently weakened to `<=`, must restore strict `<` with epsilon)

```r
test_that("ieppa_soft: strictly better max_err than ieppa on tight-bounds problem", {
  set.seed(3); n <- 5000L
  df <- data.frame(v1=factor(sample(5, n, TRUE)))
  tgt <- list(v1=setNames(c(0.4,0.3,0.15,0.1,0.05), as.character(1:5)))
  r_hard <- harvest(df, tgt, method="ieppa",
    max_weight=1.8, min_weight=0, max_iterations=500, attach_weights=FALSE)
  r_soft <- harvest(df, tgt, method="ieppa_soft",
    max_weight=1.8, min_weight=0, max_iterations=500, attach_weights=FALSE)
  me_hard <- attr(r_hard, "result")$max_error
  me_soft <- attr(r_soft, "result")$max_error
  # Require ≥1e-6 absolute improvement to confirm ALM escaped the bad fixed point.
  # Tolerance: bit-identical weights would give difference < 1e-15. A tiny
  # floating-point delta is not a real algorithmic improvement.
  expect_lt(me_soft, me_hard - 1e-6,
    label = sprintf("ieppa_soft must beat ieppa by >=1e-6: hard=%.6e, soft=%.6e",
                    me_hard, me_soft))
})
```

### T4 — capacity_penalty parameter accepted

```r
test_that("ieppa_soft: capacity_penalty=NULL routes to auto-computed value", {
  set.seed(4); n <- 1000L
  df <- data.frame(v=factor(sample(letters[1:3], n, TRUE)))
  tgt <- list(v=c(a=0.4, b=0.4, c=0.2))
  r1 <- harvest(df, tgt, method="ieppa_soft", capacity_penalty=NULL,
                attach_weights=FALSE)
  r2 <- harvest(df, tgt, method="ieppa_soft", attach_weights=FALSE)
  expect_equal(as.numeric(r1), as.numeric(r2), tolerance=1e-12)
  # Auto value should be exposed:
  cm <- attr(r1, "result")$alm_capacity_mu_final
  expect_true(is.finite(cm) && cm > 0)
})

test_that("ieppa_soft: capacity_penalty rejects invalid input", {
  df <- data.frame(v=factor(c("a","a","b")))
  tgt <- list(v=c(a=0.5, b=0.5))
  expect_error(harvest(df, tgt, method="ieppa_soft", capacity_penalty=-1),
               "positive finite scalar")
  expect_error(harvest(df, tgt, method="ieppa_soft", capacity_penalty=0),
               "positive finite scalar")
  expect_error(harvest(df, tgt, method="ieppa_soft", capacity_penalty=Inf),
               "positive finite scalar")
  expect_error(harvest(df, tgt, method="ieppa_soft", capacity_penalty=NaN),
               "positive finite scalar")
  expect_error(harvest(df, tgt, method="ieppa_soft", capacity_penalty=c(1,2)),
               "positive finite scalar")
  expect_error(harvest(df, tgt, method="ieppa_soft", capacity_penalty="0.5"),
               "positive finite scalar")
})
```

### T5 — Final bounds adherence (exact, no IEEE slack)

```r
test_that("ieppa_soft: final weights respect bounds exactly via hard projection", {
  set.seed(5); n <- 5000L
  df <- data.frame(v=factor(sample(5, n, TRUE)))
  tgt <- list(v=setNames(c(0.4,0.3,0.15,0.1,0.05), as.character(1:5)))
  r <- harvest(df, tgt, method="ieppa_soft",
               max_weight=1.8, min_weight=0.1,
               max_iterations=300, attach_weights=FALSE)
  w <- as.numeric(r)
  # Final clamp guarantees IEEE-exact bounds:
  expect_true(max(w) <= 1.8)
  expect_true(min(w) >= 0.1)
  # Discriminating: weights must NOT be identical to ieppa hard-clamp result
  # (would mean ALM didn't run, just final clamp did the work):
  r_hard <- harvest(df, tgt, method="ieppa",
                    max_weight=1.8, min_weight=0.1,
                    max_iterations=300, attach_weights=FALSE)
  expect_false(isTRUE(all.equal(as.numeric(r), as.numeric(r_hard), tolerance=1e-10)),
    label="ieppa_soft weights must differ from ieppa hard-clamp (T3 escaped local min)")
})

test_that("ieppa_soft: degenerate asymmetric bounds — final projection bounded sum drift", {
  # Adversarial: tight upper bound, loose lower; many cells will need clamping.
  set.seed(15); n <- 2000L
  df <- data.frame(v=factor(sample(c("a","b"), n, TRUE, prob=c(.95,.05))))
  tgt <- list(v=c(a=0.3, b=0.7))   # need to upweight rare 'b' heavily
  r <- harvest(df, tgt, method="ieppa_soft",
               max_weight=8.0, min_weight=0.01,
               max_iterations=200, attach_weights=FALSE)
  w <- as.numeric(r)
  expect_true(max(w) <= 8.0)
  expect_true(min(w) >= 0.01)
  # Sum drift must be < 1e-6 of n (bounded by kMaxRescaleIters * kRescaleTol * n):
  expect_lt(abs(sum(w) - n), 1e-6 * n)
})
```

### T6 — ALM converges with auto μ on stepstone-like problem

```r
test_that("ieppa_soft: auto capacity_penalty achieves max_err <= ieppa on stepstone", {
  skip_on_cran()
  if (!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet")) skip("no fixture")
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
  # Allow tiny epsilon for non-determinism on a 1.58M-row problem:
  expect_lte(me_soft, me_hard + 1e-9)
})
```

### T7 — Adaptive μ growth fires AND is observable

```r
test_that("ieppa_soft: adaptive growth scales capacity_mu when primal violation persists", {
  set.seed(7); n <- 2000L
  df <- data.frame(v=factor(sample(c("X","Y"), n, TRUE, prob=c(.3,.7))))
  tgt <- list(v=c(X=0.99, Y=0.01))   # adversarial: target 99% in minority cat
  r <- harvest(df, tgt, method="ieppa_soft",
               capacity_penalty=1e-6,    # tiny start, must grow
               max_weight=10, max_iterations=300, attach_weights=FALSE)
  res <- attr(r, "result")
  expect_true(res$status %in% c(0L, 4L, 5L))   # not INFEAS
  # Adaptive growth must actually have fired:
  expect_gt(res$alm_capacity_mu_final, 1e-6 * 1.5,
    label="capacity_mu must have grown beyond initial 1e-6")
  expect_gt(res$alm_n_growth_events, 0L,
    label="adaptive growth must have fired ≥1 time")
  # Final weights respect bounds:
  w <- as.numeric(r)
  expect_true(max(w) <= 10)
})
```

### T8 — Wall-time budget (regression guard)

```r
test_that("ieppa_soft wall-time within 1.3x of ieppa on stepstone", {
  skip_on_cran()
  if (!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet")) skip("no fixture")
  df <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet"); df$uuid <- NULL
  tgt <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
  tgt <- lapply(tgt, function(t) { t <- unlist(t); t / sum(t) })
  for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])

  t_hard <- system.time(
    harvest(df, tgt, method="ieppa", max_iterations=300, attach_weights=FALSE)
  )["elapsed"]
  t_soft <- system.time(
    harvest(df, tgt, method="ieppa_soft", max_iterations=300, attach_weights=FALSE)
  )["elapsed"]
  expect_lt(t_soft, 1.3 * t_hard,
    label = sprintf("ieppa_soft (%.2fs) must be within 1.3x of ieppa (%.2fs)", t_soft, t_hard))
})
```

---

## Acceptance Criteria

1. **T3 GREEN**: `me_soft < me_hard - 1e-6` on tight-bounds problem (strict improvement, not noise).
2. **T4 GREEN**: `capacity_penalty=NULL` routes to auto; validation rejects all invalid inputs with descriptive message.
3. **T5 GREEN**: Final weights satisfy bounds exactly; sum drift < 1e-6 of n; ieppa_soft weights provably differ from ieppa.
4. **T6 GREEN**: On stepstone, `me_soft ≤ me_hard + 1e-9`.
5. **T7 GREEN**: Adaptive growth fires and is observable via `alm_n_growth_events > 0` and `alm_capacity_mu_final > capacity_penalty * 1.5`.
6. **T8 GREEN**: Wall time within 1.3× ieppa on stepstone.
7. **Backward compat**: `method="ieppa"` produces bit-identical results before/after this change. Test:
   ```r
   r_pre  <- readRDS("tests/testthat/fixtures/ieppa_pre_alm_ref.rds")  # captured pre-merge
   r_post <- harvest(df_ref, tgt_ref, method="ieppa", ...)
   expect_equal(as.numeric(r_post), r_pre$weights, tolerance=1e-12)
   ```
8. **No regressions**: `devtools::test()` FAIL count ≤ 3 (current baseline).
9. **harvest.Rd contains `capacity_penalty`**: `grep -c "capacity_penalty" man/harvest.Rd` returns ≥ 2 (param block + details).

---

## Risk / Limitations

1. **u[c] reset on log-path fallback**: accepts loss of accumulated dual (prior version preserved this; rev 2 inverts due to scale-mismatch risk). Convergence after fallback restarts from λ=0 with current best capacity_mu_adaptive. Acceptable: fallback is rare and itself a numerical reset event.

2. **Final clamp+rescale convergence**: bounded to `kMaxRescaleIters=3` iterations. For pathological inputs (sum(L) > n), inner loop terminates without convergence and emits warning. Bounds always exact; sum-to-n approximate.

3. **Adaptive growth ceiling at 1000×**: pathological inputs may hit cap and return BUDGET/STALL. T7 documents this case. User can manually pass higher `capacity_penalty` to compensate.

4. **ALM diagnostics add 3 fields to rk_result_t**: ABI break for any external C consumer. `EXPECTED_RK_RESULT_BYTES` static_assert catches at compile time. Expected new size: 448 + 24 = 472 bytes.

5. **Non-determinism risk**: ALM blend formula uses floating-point divisions with `1+ρ` denominator. Order of evaluation in compiler optimizations could introduce tiny variations. T6/T8 use epsilon tolerances to accommodate.

6. **lambda_cell vector lifetime**: allocated once at solver entry (size M_cell), zeroed on level transitions and overflow fallback. Memory cost: 8 * M_cell bytes — negligible (e.g., 8 MB for M_cell=1M).

---

## Out of Scope

- ALM for raking solver (separate spec, different geometry — water-fill-cat per bucket)
- Per-cell adaptive μ (single global μ for now; per-cell is a future enhancement)
- Removing dead `alm_lambda` field from CalibState — it's NOT dead (used by lbfgsb)
- ALM for greg/sinkhorn (those use different convergence theory)
- Runtime tuning of `kAlmGrowthFactor` / `kAlmPersistenceThreshold` / `kAlmMaxScale` — fixed values, will be empirically validated post-implementation; if problematic, separate ticket exposes them
