# iEPPA Soft: ALM Capacity Enforcement — Design Spec (rev 3)

**Date**: 2026-04-29
**Status**: Pending design review (rev 3 addresses all 5-agent gate findings from rev 2)
**Tickets**: leafblower-scii (T3 bug), leafblower-bgx6 (alm_mu decision)
**Files**: `src/ieppa.cpp`, `src/ieppa.hpp`, `src/types.hpp`, `src/leafblower.h`,
          `src/cell_table.cpp`, `src/cell_table.hpp`, `src/r_bridge.cpp`,
          `src/c_api.cpp`, `R/harvest.R`, `man/harvest.Rd` (auto-regenerated),
          `tests/testthat/test-calibration-solvers.R`,
          `tests/testthat/fixtures/ieppa_pre_alm_ref.rds` (Step 0 prerequisite)

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

**New CalibState field** (in `src/types.hpp`):
```cpp
double capacity_mu       = 0.0;  // ieppa_soft ALM penalty coefficient (0 = inactive)
// (per-cell dual `lambda_cell[c]` lives in solver-local std::vector, not CalibState)
```

**New rk_params_t field** (in `src/leafblower.h`, ABI break):
```cpp
double capacity_penalty;  /* ieppa_soft ALM penalty; 0.0 or negative → use cell_table auto */
```

The existing `use_admm_capacity` field on CalibState is kept and gates the ALM block — kept as `use_admm_capacity` for stability (no rename despite changed semantic; rename is out of scope to avoid churn).

**Plumbing chain** (data flow trace — every link must be implemented):

```
R user → harvest(capacity_penalty=NULL|scalar)
       → harvest.R validates, passes to .Call("C_rk_calibrate", ..., capacity_penalty_sexp, ...)
       → r_bridge.cpp C_rk_calibrate reads SEXP slot
       → r_bridge.cpp populates rk_params_t.capacity_penalty
       → r_bridge.cpp calls build_cell_table → CellTable.capacity_mu_auto = M_cell/n
       → r_bridge.cpp resolves p.capacity_penalty:
           if (capacity_penalty <= 0)  st.capacity_mu = ct.capacity_mu_auto
           else                         st.capacity_mu = capacity_penalty
       → r_bridge.cpp dispatch sets st.use_admm_capacity = true (for IEPPA_SOFT)
       → ieppa.cpp ieppa_solve reads st.capacity_mu and st.use_admm_capacity
       → ieppa.cpp capacity block uses st.capacity_mu in ALM formula
       → ieppa.cpp populates res fields (alm_capacity_mu_final, alm_n_growth_events, alm_max_dual_norm)
       → r_bridge.cpp packs result fields into rk_result_t
       → r_bridge.cpp exposes via attr(r, "result")$alm_*
       → R user reads attr(r, "result")$alm_capacity_mu_final
```

**Same chain for c_api.cpp direct C callers** — except they bypass R-layer validation (documented constraint).

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
// Bound dual to prevent runaway. λ has units mass⁻¹ (matches X̃⁻¹ = M_cell/n scale).
// Cap = 10·capacity_mu_base·max_weight matches the auto scale: when |λ|≈cap and
// |X − z|≈max_weight, the dual update Δλ = μ·(X − z) ≈ μ·max_weight cannot exceed cap
// in a single step. Prevents one-iteration runaway while leaving multi-iteration headroom.
const double lambda_cap = 10.0 * capacity_mu_base * st.max_weight;
lambda_cell[c] = std::clamp(lambda_cell[c], -lambda_cap, lambda_cap);
```

**Mathematical derivation** (linearized Newton step from s=1, NOT exact minimizer):

The exact ALM minimizer requires solving a transcendental equation (KL gradient is `log(s)`, penalty gradient is linear in s). The closed-form below is a single Newton step from s=1, which is exact when X[c] = X̃[c] and a good approximation when |X − X̃| / X̃ is small. Multiple Newton iterations per outer step are out of scope — single step has worked empirically for similar Bregman proximal solvers (Boyd 2011 §3.4).


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

Persistent state across homotopy levels (declared at solver entry, before homotopy loop):
```cpp
const double capacity_mu_base = st.capacity_mu;       // user-supplied (or auto from cell_table)
double capacity_mu_adaptive   = capacity_mu_base;     // adaptive scale, persists across levels
int    alm_violation_streak   = 0;                    // resets on level transition
double eta_i_current          = 1.0;                  // current Tang-η factor; persists for adaptive use
const bool tang_active = (st.eta_schedule.mode == EtaScheduleMode::TANG_DYNAMIC && N_levels > 1);
```

At each homotopy level entry (replace existing alm_mu logic at ieppa.cpp:381-391):
```cpp
if (tang_active) {
    const double scaled_frac = std::pow(frac, st.eta_schedule.schedule_power);
    eta_i_current = st.eta_schedule.eta_start *
        std::pow(st.eta_schedule.eta_end / st.eta_schedule.eta_start, scaled_frac);
    res.eta_final = eta_i_current;
    st.capacity_mu = eta_i_current * capacity_mu_adaptive;  // compose: Tang-η × adaptive
} else {
    eta_i_current = 1.0;
    st.capacity_mu = capacity_mu_adaptive;                   // no Tang-η: pure adaptive
}
alm_violation_streak = 0;  // reset streak on level transition
// Reset per-cell dual on level transition AND on linear→log fallback (rev 2 decision 5b'):
if (st.use_admm_capacity) std::fill(lambda_cell.begin(), lambda_cell.end(), 0.0);
```

Per-iteration adaptive growth (after capacity block, eta_i_current already declared above):
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
        st.capacity_mu = eta_i_current * capacity_mu_adaptive;  // immediate effect, eta_i tracked
        res.alm_n_growth_events++;
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

**`@param method` update** — must include `ieppa_soft` in the enumeration:
```r
#' @param method Calibration algorithm to use. One of: \code{"auto"} (default,
#'   selects based on data shape), \code{"ieppa"} (iterative entropic proximity
#'   projection with hard capacity clamp), \code{"ieppa_soft"} (iEPPA with
#'   augmented Lagrangian soft capacity enforcement; better than \code{"ieppa"}
#'   on tight-bounds problems where many cells hit \code{max_weight}),
#'   \code{"raking"} (cyclic IPF with water-filling), \code{"lbfgsb"}, etc.
```

**`@param capacity_penalty`** — new param with concrete user guidance:
```r
#' @param capacity_penalty Numeric, controls how strongly capacity bounds are
#'   enforced during ALM optimization in \code{method="ieppa_soft"}. Use
#'   \code{NULL} (default) for auto-computed value (M_cell / n) which gives a
#'   balanced 50/50 blend between unconstrained KL minimization and hard-clamp
#'   projection. Larger values force tighter constraint adherence at the cost
#'   of slower margin convergence; smaller values allow more constraint
#'   violation during optimization but are pulled back to the feasible set by
#'   the final projection.
#'
#'   \strong{Tuning guidance}: After running \code{ieppa_soft}, inspect
#'   \code{attr(result, "result")$alm_capacity_mu_final}. If
#'   \code{alm_capacity_mu_final / capacity_penalty_input >= 1000}, the adaptive
#'   schedule hit its growth ceiling — supply a manual \code{capacity_penalty}
#'   that is \code{10x} larger and re-run. If \code{alm_n_growth_events == 0}
#'   and \code{max_error} is acceptable, the auto value worked.
#'
#'   Ignored for methods other than \code{"ieppa_soft"}; passing it with another
#'   method emits a warning.
```

**`@details`** addition for method selection guidance:
```r
#' @details
#' \strong{Choosing between \code{ieppa} and \code{ieppa_soft}}: Use
#' \code{method="ieppa"} (default for AUTO routing on tight-bounds inputs)
#' when capacity bounds are slack — most cells will not hit \code{max_weight}.
#' Switch to \code{method="ieppa_soft"} when:
#' \itemize{
#'   \item \code{ieppa} returns with \code{max_error > 1e-3} on a feasible problem
#'   \item Many observations cluster near \code{max_weight} (visible in result
#'         diagnostics: \code{n_bounds_clamped > 0.05 * n})
#'   \item You need the calibration to find the best constrained optimum rather
#'         than getting stuck at the boundary
#' }
#' \code{ieppa_soft} is roughly 10-30\% slower than \code{ieppa} due to ALM
#' bookkeeping; final weights respect bounds exactly (hard projection), but
#' \code{sum(weights)} may differ from \code{n} by up to \code{1e-6 * n} on
#' adversarial inputs (reported in \code{alm_sum_drift}).
```

**`@return` itemize update** — add 4 ALM diagnostic fields:
```r
#' \item \code{alm_capacity_mu_final}: final value of the ALM penalty after
#'   adaptive scaling (\code{method="ieppa_soft"} only; \code{0} otherwise).
#' \item \code{alm_n_growth_events}: count of adaptive penalty growth events
#'   during solving. \code{0} = auto value sufficed; high values suggest
#'   manually increasing \code{capacity_penalty}.
#' \item \code{alm_max_dual_norm}: maximum absolute Lagrange dual at exit.
#'   Near the cap (\code{10 * capacity_mu_base * max_weight}) suggests very
#'   binding constraints or infeasibility.
#' \item \code{alm_sum_drift}: \code{|sum(weights) - n|} after final projection.
#'   Bounded by \code{1e-6 * n} on well-conditioned inputs.
```

**Function signature** (insertion point after `max_weight`, before `...`):
```r
harvest <- function(data, target,
                    method = "auto",
                    max_weight = 5,
                    min_weight = 0,
                    capacity_penalty = NULL,    # NEW (rev 3)
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
after `start_weights` (slot 8). All subsequent slots shift by 1.

**Complete slot table after insertion (32 args total)**:

| Slot | Name | Type | Source |
|------|------|------|--------|
| 1 | data | data.frame | df |
| 2 | target | list | targets |
| 3 | min_weight | double | param |
| 4 | max_weight | double | param |
| 5 | method | character | param |
| 6 | verbose | integer | param |
| 7 | max_iterations | integer | param |
| 8 | start_weights | double or NULL | param |
| **9** | **capacity_penalty** | **double or NULL** | **NEW (rev 3)** |
| 10 | tol_abs (legacy) | double | hardcoded |
| 11 | bounds_mode_int | integer | param |
| 12 | homotopy_levels | integer | param |
| 13 | homotopy_start_factor | double | param |
| 14 | homotopy_end_factor | double | param |
| 15 | homotopy_budget_p | double | param |
| 16 | scheduler | character | param |
| 17 | eta_schedule | character | param |
| 18 | eta_start | double | param |
| 19 | eta_end | double | param |
| 20 | eta_schedule_power | double | param |
| 21 | conv_pct_tol | double | parsed conv |
| 22 | conv_absolute_tol | double | parsed conv |
| 23 | conv_metric_int | integer | parsed conv |
| 24 | conv_rule_int | integer | parsed conv |
| 25 | conv_stop_when_int | integer | parsed conv |
| 26 | sor_enabled | integer | param |
| 27 | sor_auto | integer | param |
| 28 | sor_omega_init | double | param |
| 29 | sor_omega_min | double | param |
| 30 | sor_omega_fixed | double | param |
| 31 | sor_burnin | integer | param |
| 32 | accelerate_bool | integer | param |

**Mandatory grep audit** (acceptance criterion 10): `grep -rn 'C_rk_calibrate' R/ tests/ benchmarks/ python/ src/r_bridge.cpp` must show ALL call sites updated. Current sites: `R/harvest.R` (one .Call), `tests/testthat/test-safety.R` (B7 test), `src/r_bridge.cpp` registration table. Any positional caller (none found in benchmarks, but grep is mandatory before merge) must insert NULL or value at slot 9.

**Registration update**: `src/r_bridge.cpp` `R_CallMethodDef` for `C_rk_calibrate` must change arity from 31 to 32. The static_assert tripwire: if R sees a 31-arg call to a 32-arg registration, R errors with "Incorrect number of arguments to .Call". This catches any missed call site at runtime in test execution.

### ALM Diagnostics Exposed

Add to `rk_result_t` (in `src/leafblower.h`):
```c
double alm_capacity_mu_final;   /* final capacity_mu after adaptive scaling; 0 if not ieppa_soft */
int    alm_n_growth_events;     /* count of adaptive growth fires; 0 if not ieppa_soft */
double alm_max_dual_norm;       /* max |lambda_cell[c]| at solver exit; for monitoring */
double alm_sum_drift;           /* |sum(X) - n| after final projection; user-visible bound check */
```

**ABI update**: 4 new fields × {8, 4, 8, 8} bytes = 28B + 4B padding = 32B.
- Old `EXPECTED_RK_RESULT_BYTES = 448`
- New `EXPECTED_RK_RESULT_BYTES = 480` (update in `src/leafblower.h:171` along with comment line 170)

**ABI update for rk_params_t**: 1 new field (8B) + 0 padding = 8B.
- Old `EXPECTED_RK_PARAMS_BYTES = 224`
- New `EXPECTED_RK_PARAMS_BYTES = 232` (update in `src/leafblower.h:189`)

Exposed in `attr(r, "result")` (4 new R-level fields):
- `alm_capacity_mu_final / capacity_penalty_input ≈ 1.0` → adaptive growth not needed
- `alm_capacity_mu_final / capacity_penalty_input ≥ 1000` → hit cap, increase `capacity_penalty` manually by 10×
- `alm_max_dual_norm` near `lambda_cap = 10·capacity_mu_base·max_weight` → constraint very binding, may indicate infeasibility
- `alm_sum_drift` should be < 1e-9·n on well-conditioned inputs; > 1e-6·n indicates infeasible bounds

### Cross-method behavior: capacity_penalty ignored for non-ieppa_soft

If user passes `capacity_penalty` with `method != "ieppa_soft"`, harvest.R emits a warning:
```r
if (!is.null(capacity_penalty) && method != "ieppa_soft") {
  warning("leafblower: capacity_penalty is only used by method='ieppa_soft'; ignored for method='", method, "'",
          call. = FALSE)
}
```

This matches the existing pattern for `accelerate` (harvest.R:206) and gives users explicit feedback rather than silent param drop.

### Files to Modify (with implementation order — Step 0 prerequisite)

| Order | File | Change |
|-------|------|--------|
| **0** | **`tests/testthat/fixtures/ieppa_pre_alm_ref.rds`** | **PREREQUISITE: capture pre-merge ieppa weights via `data-raw/gen_ieppa_pre_alm_ref.R`. Must commit fixture BEFORE step 1.** |
| 1 | `src/types.hpp` | Add `double capacity_mu = 0.0;` to CalibState |
| 2 | `src/cell_table.hpp` | Add `double capacity_mu_auto = 0.0;` to CellTable struct |
| 3 | `src/cell_table.cpp` | `build_cell_table`: compute `capacity_mu_auto = M_cell/n` (with `n>0 && M_cell>0` guard) |
| 4 | `src/leafblower.h` | Add 4 `alm_*` fields to `rk_result_t`; add `capacity_penalty` to `rk_params_t`; bump `EXPECTED_RK_RESULT_BYTES` (448→480) and `EXPECTED_RK_PARAMS_BYTES` (224→232); update tripwire comments |
| 5 | `src/r_bridge.cpp` | Read `capacity_penalty` SEXP with `Rf_isNull`/`LENGTH` guards; populate `p.capacity_penalty`; resolve auto vs manual; update `R_CallMethodDef` arity 31→32; pack 4 ALM diagnostic fields into result |
| 6 | `src/c_api.cpp` | Same routing; add comment block documenting validation contract for direct C callers (R-layer guards do not apply; callers must validate) |
| 7 | `src/ieppa.cpp` | Replace P1.1 capacity block (lines 800-816) with ALM update; declare persistent `capacity_mu_adaptive` / `eta_i_current` / `alm_violation_streak`; add lambda_cell vector (size M_cell); update u[] reset policy (level transitions + linear→log fallback); add adaptive growth; add final clamp+rescale projection; populate 4 result fields |
| 8 | `src/ieppa.hpp` | (no changes — diagnostics carried via `rk_result_t`) |
| 9 | `R/harvest.R` | Add `capacity_penalty = NULL` param at signature position 6; validation block; pass to `.Call` slot 9; warn when used with non-ieppa_soft method; update `@param method` (add ieppa_soft); add `@param capacity_penalty`; add `@details` method-selection guidance; add 4 ALM fields to `@return` itemize |
| 10 | `man/harvest.Rd` | Auto-regenerated by `devtools::document()` |
| 11 | `tests/testthat/test-calibration-solvers.R` | T3 strengthened; T4-T7 new; T9 (backward compat reading Step 0 fixture); T10 (warning) |
| 12 | `tests/testthat/test-safety.R` | Update existing `.Call` test (B7) to pass NULL at slot 9 (capacity_penalty) |
| 13 | `data-raw/gen_ieppa_pre_alm_ref.R` | Script that produced Step 0 fixture (committed alongside fixture for reproducibility) |

**Compile order**: 1→2→3 must compile clean before 4 (cell_table.hpp included by leafblower.h consumers? actually no — types.hpp and cell_table.hpp are independent of leafblower.h). 1, 2, 3, 4 are independent — pick any order, then 5/6 (which depend on 4 + cell_table). Then 7 (depends on 1, 2). Then 9 (depends on R-side ABI being stable in 4). Then 10 auto. Then 11/12 last.

**Step 0 procedure** (must complete before any code change):
```bash
# In a clean checkout of the pre-merge state:
Rscript data-raw/gen_ieppa_pre_alm_ref.R
# Script content:
#   1. Set seed; generate canonical 5000-row 5-category test problem
#   2. Run harvest(df, tgt, method="ieppa", max_weight=1.8, ...)
#   3. saveRDS(list(df=df, tgt=tgt, max_weight=1.8, min_weight=0,
#                    max_iterations=500, convergence=list(...),
#                    weights=as.numeric(r),
#                    result=attr(r,"result")),
#               "tests/testthat/fixtures/ieppa_pre_alm_ref.rds")
git add tests/testthat/fixtures/ieppa_pre_alm_ref.rds data-raw/gen_ieppa_pre_alm_ref.R
git commit -m "test(ieppa): capture pre-ALM reference fixture for backward compat regression"
# Only THEN start steps 1-13.
```

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

### T8 — Wall-time budget (regression guard, interleaved)

```r
test_that("ieppa_soft wall-time within 1.5x of ieppa on stepstone (interleaved)", {
  skip_on_cran()
  skip_if_not_installed("bench")
  if (!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet")) skip("no fixture")
  df <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet"); df$uuid <- NULL
  tgt <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
  tgt <- lapply(tgt, function(t) { t <- unlist(t); t / sum(t) })
  for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])

  # Interleaved comparison via bench::mark — eliminates GC/scheduler variance
  # that would dominate sequential system.time() measurements.
  res <- bench::mark(
    ieppa      = harvest(df, tgt, method="ieppa",      max_iterations=300, attach_weights=FALSE),
    ieppa_soft = harvest(df, tgt, method="ieppa_soft", max_iterations=300, attach_weights=FALSE),
    iterations    = 3,
    check         = FALSE,           # weights differ by design (T3)
    min_iterations = 3,
    filter_gc      = FALSE
  )
  t_hard <- as.numeric(res$median[res$expression == "ieppa"])
  t_soft <- as.numeric(res$median[res$expression == "ieppa_soft"])
  expect_lt(t_soft, 1.5 * t_hard,
    label = sprintf("ieppa_soft median (%.2fs) must be within 1.5x of ieppa median (%.2fs)",
                    t_soft, t_hard))
})
```

### T9 — Backward compatibility for `method="ieppa"` (regression guard)

```r
test_that("method='ieppa' produces bit-identical weights vs pre-ALM reference fixture", {
  # Fixture captured pre-merge as Step 0 (see Implementation Order below).
  # Hash-locked: any change to ieppa hard-clamp behavior fails this test.
  fixture_path <- testthat::test_path("fixtures/ieppa_pre_alm_ref.rds")
  if (!file.exists(fixture_path)) {
    skip("Step 0 fixture not captured — run data-raw/gen_ieppa_pre_alm_ref.R first")
  }
  ref <- readRDS(fixture_path)

  # Reference inputs encoded in the fixture for self-containment:
  df  <- ref$df
  tgt <- ref$tgt

  r <- harvest(df, tgt, method="ieppa",
               max_weight=ref$max_weight, min_weight=ref$min_weight,
               max_iterations=ref$max_iterations,
               convergence=ref$convergence, attach_weights=FALSE)
  w_post <- as.numeric(r)
  expect_equal(w_post, ref$weights, tolerance=1e-12,
               label="method='ieppa' must produce bit-identical weights pre/post ALM merge")
  # Result attribute fields also unchanged (excluding new ALM fields, which
  # should be 0/NA for non-ieppa_soft):
  res_post <- attr(r, "result")
  expect_equal(res_post$status,     ref$result$status)
  expect_equal(res_post$iterations, ref$result$iterations)
  expect_equal(res_post$max_error,  ref$result$max_error, tolerance=1e-12)
})
```

### T10 — `capacity_penalty` warns when used with non-ieppa_soft method

```r
test_that("capacity_penalty emits warning when passed to non-ieppa_soft method", {
  set.seed(10); n <- 100L
  df <- data.frame(v=factor(sample(c("a","b"), n, TRUE)))
  tgt <- list(v=c(a=0.5, b=0.5))
  expect_warning(
    harvest(df, tgt, method="ieppa", capacity_penalty=0.5, attach_weights=FALSE),
    regexp = "capacity_penalty.*ieppa_soft.*ignored"
  )
  expect_warning(
    harvest(df, tgt, method="raking", capacity_penalty=0.5, attach_weights=FALSE),
    regexp = "capacity_penalty.*ieppa_soft.*ignored"
  )
})
```

---

## Acceptance Criteria

1. **T3 GREEN**: `me_soft < me_hard - 1e-6` on tight-bounds problem (strict improvement, not noise).
2. **T4 GREEN**: `capacity_penalty=NULL` routes to auto; validation rejects all invalid inputs (NaN, Inf, 0, negative, vector, string) with descriptive message.
3. **T5 GREEN**: Final weights satisfy bounds exactly via IEEE comparison; sum drift < 1e-6·n on adversarial asymmetric bounds; ieppa_soft weights provably differ from ieppa hard-clamp result on the same input.
4. **T6 GREEN**: On stepstone fulldata fixture, `me_soft ≤ me_hard + 1e-9`. Skip on CRAN. Manual verification required if fixture absent.
5. **T7 GREEN**: Adaptive growth fires on adversarial start (`capacity_penalty=1e-6`, target=99%/1%): `alm_n_growth_events > 0` AND `alm_capacity_mu_final > capacity_penalty_input * 1.5`.
6. **T8 GREEN**: `bench::mark(min_iterations=3)` median wall time of ieppa_soft within 1.5× of ieppa on stepstone.
7. **T9 GREEN (backward compat)**: `method="ieppa"` produces bit-identical weights to pre-ALM reference fixture (Step 0). Tolerance 1e-12.
8. **T10 GREEN (warning)**: `capacity_penalty` passed with `method != "ieppa_soft"` emits warning matching `capacity_penalty.*ieppa_soft.*ignored`.
9. **No regressions**: `devtools::test()` FAIL count exactly equal to pre-merge baseline (currently 3, locked at PR open time, not "≤3").
10. **`.Call` arity audit**: `grep -rn 'C_rk_calibrate' R/ tests/ benchmarks/ python/ src/r_bridge.cpp` shows ALL call sites updated to 32 args; CI run confirms no "incorrect number of arguments" errors.
11. **harvest.Rd content**: `grep -c "capacity_penalty" man/harvest.Rd` returns ≥ 3 (param block + details + return block); `grep -c "ieppa_soft" man/harvest.Rd` returns ≥ 4 (method param + details + each instance in the fix description).
12. **Diagnostics observable**: `attr(r, "result")` exposes `alm_capacity_mu_final`, `alm_n_growth_events`, `alm_max_dual_norm`, `alm_sum_drift` for `method="ieppa_soft"`. Type checked: numeric scalars, finite (or 0/NA for non-ieppa_soft).
13. **ABI tripwires fire correctly**: `EXPECTED_RK_RESULT_BYTES = 480` and `EXPECTED_RK_PARAMS_BYTES = 232` static_asserts compile clean (= correct field count).
14. **Step 0 fixture committed first**: `git log --oneline tests/testthat/fixtures/ieppa_pre_alm_ref.rds` shows commit BEFORE any commit touching `src/ieppa.cpp` ALM code.

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
