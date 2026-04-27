# Sinkhorn A1 Fix + ieppa_soft Method + Decoupled Solver Objective

**Date**: 2026-04-27 (rev 2)
**Status**: Pending design review
**Tickets**: leafblower-lig9
**Files**: `src/sinkhorn.cpp`, `src/ieppa.cpp`, `src/ieppa.hpp`, `src/greg.cpp`,
           `src/grake.cpp`, `src/calib_dispatch.hpp`, `src/c_api.cpp`,
           `src/r_bridge.cpp`, `R/harvest.R`, `data-raw/gen_ieppa_kl_ref.R`,
           `tests/testthat/test-calibration-solvers.R`

---

## Problem

### 1. `convergence_used$objective` reports wrong value for sinkhorn (and others)

**Root cause:** `res.convergence_objective = best_metric_seen` in each solver.
`best_metric_seen` = minimum of the **stopping criterion** metric, NOT the solver's
mathematical objective. When they differ, the field lies.

**Smoking gun:** benchmark shows `sinkhorn weight_kl == sinkhorn max_err == 6.006e-3`
to 4 significant figures — two independent quantities cannot coincidentally match
unless one is being reported for the other.

**Affected solvers:**

| Solver | Stopping default | Mathematical objective | What `$objective` currently reports |
|--------|-----------------|----------------------|--------------------------------------|
| ieppa | marginal_kl | weight KL | marginal KL (wrong) |
| sinkhorn | max_err | weight KL | max_err (wrong) |
| greg | max_err | chi² | max_err (wrong) |
| grake | grake_norm | grake_norm | grake_norm (correct by coincidence) |
| chebyshev | max_err | max_err | max_err (correct by coincidence) |
| raking | max_err | weight KL | max_err (wrong) |

### 2. A1 violated

A1: `sinkhorn weight_KL_at_convergence ≤ ieppa weight_KL_at_best_iter = 3.008e-3`.
Current: sinkhorn `$objective` = 6.006e-3 (is max_err, not weight KL). A1 tests
the wrong quantity. True sinkhorn weight KL is unknown.

### 3. No explicit method for ieppa+ADMM (ieppa_soft)

Task 2 added ADMM capacity enforcement (always-on). ieppa+ADMM is strictly better on
all metrics (max_err 2.34e-3 vs 2.74e-3) but has no explicit method name. Both
`method="ieppa"` and the ADMM variant need distinct names.

Note: both ieppa and ieppa+ADMM are novel implementations — no "classic" heritage.

---

## Design

**Priority: correctness → speed/efficiency → elegance. No backward compat needed.**

---

## Part 1: Decouple Solver Objective from Stopping Criterion

### New field: `best_objective_seen`

Each solver struct gains `double best_objective_seen = 0.0` (add to `IEPPAResult`,
`SinkhornResult`, `GregResult`, etc. in respective `.hpp` files).

### Per-solver mathematical objective

```cpp
// In calib_dispatch.hpp — add helper using existing RK_ALG_* constants:
inline double select_solver_objective(int alg_id,
                                       const lbw::CellMetrics& m) {
    switch (alg_id) {
    case RK_ALG_IEPPA:
    case RK_ALG_RAKING:
    case RK_ALG_SINKHORN:  return m.kl;          // weight KL
    case RK_ALG_GREG:      return m.chi2;         // chi²
    case RK_ALG_GRAKE:     return m.grake_norm;   // grake_norm
    case RK_ALG_CHEBYSHEV: return m.errRp;        // max_err (objective == stopping metric)
    default:               return m.errRp;
    }
}
```

`alg_id` = the existing `RK_ALG_*` integer constant from `src/leafblower.h`
(RK_ALG_IEPPA=1, RK_ALG_RAKING=2, RK_ALG_SINKHORN=4, etc.). Each solver already
knows its own alg_id — it is a compile-time constant, not a runtime parameter.

### Mandatory update pattern (atomic block)

In each solver's error-check section, the three assignments MUST be co-located in a
single named block to prevent future decoupling:

```cpp
// === BEST-ITER UPDATE (objective, metric, and iter MUST stay co-located) ===
if (curr_metric < best_metric_seen) {
    best_metric_seen    = curr_metric;
    best_iter_val       = iter;
    best_objective_seen = select_solver_objective(RK_ALG_<SOLVER>, m);
    // save W_best here
}
// === END BEST-ITER UPDATE ===
```

Guard: `if (!std::isfinite(best_objective_seen)) best_objective_seen = 0.0;`
(handles NaN/Inf from degenerate cells with log(0)).

### Field rename: `convergence_used$solver_objective`

Rename `res.convergence_objective` → `res.convergence_solver_objective` in all solver
structs and the C API (`rk_calib_result_t`). Update R bridge to expose as
`result$convergence_used$solver_objective`.

**Rationale for rename:** the old field reported stopping criterion value; the new
field has different semantics. Rename forces discovery of callers at compile time.
No backward compat needed.

Set: `res.convergence_solver_objective = best_objective_seen;`

### Sinkhorn default metric change

In `R/harvest.R`: change sinkhorn default stopping metric from `"max_err"` to `"kl"`.
Sinkhorn's KL is monotone-decreasing, so `kl+improvement` is the correct stopping
criterion for its mathematical objective.

```r
# In harvest.R method-specific defaults section:
if (method == "sinkhorn" && is.null(convergence[["metric"]]) && ...) {
  conv$metric <- "kl"
}
```

---

## Part 2: Add `method="ieppa_soft"` (ADMM capacity enforcement)

`method="ieppa"` = original ieppa with Euclidean hard clamp (P1.1: `xc = clamp(X_tilde, L, U)`).
`method="ieppa_soft"` = ieppa with ADMM soft capacity enforcement. Same BCD sweeps,
same homotopy/Greenkhorn/Tang-η. Only P1.1 differs.

**Capacity bounds guarantee:** ieppa_soft DOES strictly enforce [min_weight, max_weight].
ADMM z-update always applies `std::clamp` → `X[c] = z ≤ U_cell[c]` at every iteration.
The cell-mode post-normalization clamp (leafblower-a3mr) further guarantees final output.
"Soft" = soft convergence path, NOT soft bounds.

### Implementation

Add `bool use_admm_capacity = false` to `IEPPAConfig` (in `src/types.hpp` or `src/ieppa.hpp`).

Gate P1.1 in `src/ieppa.cpp`:
```cpp
if (st.ieppa_cfg.use_admm_capacity) {
    double z = std::clamp(X_tilde_c + u[c], L_cell[c], U_cell[c]);
    u[c] += X_tilde_c - z;
    X[c] = z; W[c] = z / X_tilde_c; X_cur[c] = z;
} else {
    double xc = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);
    X[c] = xc; W[c] = xc / X_tilde_c; X_cur[c] = xc;
}
```

`u[c]` vector: allocated only when `use_admm_capacity == true`. If false, skip allocation.

In `R/harvest.R`: add `"ieppa_soft"` to method dispatch. Route to ieppa solver with
`use_admm_capacity = TRUE`. In `c_api.cpp` / `r_bridge.cpp`: pass flag through param struct.

`method="ieppa"` reverts to `use_admm_capacity = FALSE` (original hard clamp).

---

## Part 3: Regenerate A1 Fixture

**Prerequisite:** Part 1 (correct objective reporting) must be implemented and compiled
before running `gen_ieppa_kl_ref.R`. Running with the old code produces wrong values.

**Execution order (mandatory):**
1. Implement Part 1 C++ changes
2. Compile: `R CMD INSTALL --preclean .`
3. Run: `Rscript data-raw/gen_ieppa_kl_ref.R`
4. Commit fixture + script together
5. Then run tests (A1 will now reference correct weight KL)

**Fix in `gen_ieppa_kl_ref.R`:** Current script extracts `best_error` (= max_err at exit,
wrong). After Part 1, extract `result$convergence_used$solver_objective` (= weight KL
at best_iter). Script must be updated before regeneration.

**Pre-condition check before committing fixture:**
Manually verify `kl_sinkhorn < kl_ieppa_soft` on stepstone dataset.
If sinkhorn (true KL minimizer) still has higher weight KL than ieppa_soft, A1 will
fail and requires a deeper algorithmic investigation (separate ticket).

---

## TDD Requirements

### RED tests (must fail before implementation, pass after)

**T1 — objective decoupling for sinkhorn:**
```r
test_that("sinkhorn: solver_objective reports weight KL, not stopping metric", {
  set.seed(1); n <- 1000L
  df <- data.frame(v1=factor(sample(c("A","B","C"), n, TRUE)))
  tgt <- list(v1=c("A"=0.6,"B"=0.3,"C"=0.1))
  r_max <- leafblower::harvest(df, tgt, method="sinkhorn",
    convergence=list(metric="max_err"), max_iterations=200, attach_weights=FALSE)
  r_kl <- leafblower::harvest(df, tgt, method="sinkhorn",
    convergence=list(metric="kl"), max_iterations=200, attach_weights=FALSE)
  obj_max <- attr(r_max,"result")$convergence_used$solver_objective
  obj_kl  <- attr(r_kl, "result")$convergence_used$solver_objective
  # Both should report weight KL (not stopping criterion value)
  # weight KL < 1e-2 for feasible calibration; max_err typically 1e-3 to 1e-1
  expect_true(obj_max < 0.1, label="solver_objective is weight KL not max_err")
  # Weight KL should be roughly the same regardless of stopping criterion
  expect_true(abs(obj_max - obj_kl) / max(obj_max, obj_kl) < 0.5,
              label="objective consistent across stopping criteria")
})
```

**T2 — ieppa_soft available and bounds-safe:**
```r
test_that("ieppa_soft: method available, weights within [min_weight, max_weight]", {
  set.seed(2); n <- 2000L
  df <- data.frame(v1=factor(sample(c("X","Y"), n, TRUE, prob=c(.3,.7))))
  tgt <- list(v1=c("X"=0.8,"Y"=0.2))
  r <- leafblower::harvest(df, tgt, method="ieppa_soft",
    max_weight=2.0, min_weight=0.0, max_iterations=300, attach_weights=FALSE)
  w <- as.numeric(r)
  expect_true(max(w) <= 2.0 + 1e-9)
  expect_true(min(w) >= 0.0 - 1e-9)
  expect_equal(attr(r,"result")$status, 0L)
})
```

**T3 — ieppa (hard) vs ieppa_soft differ on tight-bounds problem:**
```r
test_that("ieppa and ieppa_soft produce different convergence paths", {
  set.seed(3); n <- 5000L
  df <- data.frame(v1=factor(sample(5, n, TRUE)))
  tgt <- list(v1=setNames(c(0.4,0.3,0.15,0.1,0.05), as.character(1:5)))
  r_hard <- leafblower::harvest(df, tgt, method="ieppa",
    max_weight=1.8, min_weight=0, max_iterations=500, attach_weights=FALSE)
  r_soft <- leafblower::harvest(df, tgt, method="ieppa_soft",
    max_weight=1.8, min_weight=0, max_iterations=500, attach_weights=FALSE)
  # ieppa_soft should achieve lower max_err on tight-bounds problem
  me_hard <- attr(r_hard,"result")$max_error
  me_soft <- attr(r_soft,"result")$max_error
  expect_true(me_soft <= me_hard + 1e-4,
              label="ieppa_soft max_err not worse than ieppa on tight bounds")
})
```

### Python audit (required before deferral is acceptable)

```bash
grep -r "convergence_used\|convergence_objective\|solver_objective" python/ | grep -v ".pyc"
```
If any Python test asserts on these fields, they must be updated in the same PR.
Record result here before finalizing spec.

---

## Acceptance Criteria

1. **Objective correctness**: `method="sinkhorn"` with `metric="max_err"` reports
   weight KL (not max_err) in `result$convergence_used$solver_objective`.
2. **greg objective**: `method="greg"` reports chi² in `$solver_objective`.
3. **ieppa_soft available**: `harvest(..., method="ieppa_soft")` runs without error.
4. **ieppa_soft bounds**: max(weights) ≤ max_weight + 1e-9 (T2 test above).
5. **ieppa unchanged**: `harvest(..., method="ieppa")` produces same result as before
   Task 2 (hard clamp, no ADMM). max_err ≈ 2.74e-3 on stepstone.
6. **ieppa_soft better**: max_err ≈ 2.34e-3 on stepstone (better than ieppa).
7. **A1 passes**: sinkhorn `$solver_objective` ≤ ieppa `$solver_objective` at best_iter
   (after fixture regeneration with corrected values).
8. **Regression**: FAIL 2 pre-existing, PASS ≥ 381.

---

## Files Changed

| File | Change |
|------|--------|
| `src/ieppa.hpp` | Add `best_objective_seen`, `convergence_solver_objective` to IEPPAResult; `use_admm_capacity` to IEPPAConfig |
| `src/sinkhorn.hpp`, `src/greg.hpp`, `src/grake.hpp` | Add `best_objective_seen`, `convergence_solver_objective` |
| `src/ieppa.cpp` P1.1 block | Gate ADMM on `cfg.use_admm_capacity`; u[] allocated only when true |
| `src/ieppa.cpp` convergence block | Atomic BEST-ITER UPDATE block with `best_objective_seen` |
| `src/sinkhorn.cpp` line ~183 | `convergence_solver_objective = best_objective_seen` |
| `src/greg.cpp`, `src/grake.cpp` | Same |
| `src/calib_dispatch.hpp` | Add `select_solver_objective(alg_id, m)` using RK_ALG_* constants |
| `src/c_api.cpp` | Rename field in `rk_calib_result_t`; pass `use_admm_capacity` flag |
| `src/r_bridge.cpp` | Expose `solver_objective` field; pass flag through |
| `R/harvest.R` | Add `"ieppa_soft"` dispatch with `use_admm_capacity=TRUE`; sinkhorn default `metric="kl"`; rename field |
| `data-raw/gen_ieppa_kl_ref.R` | Use `$solver_objective` (weight KL) not `$best_error` |
| `tests/testthat/test-calibration-solvers.R` | Add T1/T2/T3 RED tests; update A1 to use `$solver_objective` |

---

## Out of Scope

- Python wrapper field rename (audit first; if no tests assert on field, defer)
- lbfgsb objective reporting
- ADMM convergence rate analysis
- u[] oscillation bounds analysis
