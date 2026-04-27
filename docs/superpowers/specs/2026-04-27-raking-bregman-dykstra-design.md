# Raking: Bregman Dykstra + SOR + Greedy A/B Test

**Date**: 2026-04-27 (rev 2 — design-review-gate round-1 fixes)
**Status**: Pending design review
**File**: `src/raking.cpp` only (plus fixture + docstring)

---

## Problem

Raking's bounded calibration uses a **hybrid metric**: IPF marginal updates (multiplicative / KL-Bregman geometry) mixed with Euclidean additive Dykstra box and hyperplane projections. The code documents this limitation at lines 13-16. Consequence: no unified convergence guarantee.

On `stepstone-fulldata` (n=1.58M, K=9): raking max_err=5.44e-3 vs ieppa 2.34e-3 (2.3× gap). `ipfn` (unconstrained, pure multiplicative) achieves 1.59e-3 — bounding alone does not explain the gap.

---

## Fix: Three independent improvements

### 1. Bregman/KL Dykstra (always-on, base fix)

Replaces additive Euclidean Dykstra corrections with multiplicative KL-Dykstra corrections. Bauschke & Lewis (2000) guarantees convergence to the KL minimum if feasible.

**Current data structures** (confirmed from source):
- `p` = `std::vector<double>(ct.M_cell, 0.0)` — per-cell box correction (vector)
- `q_hyp` = `double = 0.0` — scalar hyperplane correction

**Box correction** (`p[c]`, init 0→1, additive→multiplicative ratio):
```cpp
// Before: yc = X[c] + p[c]; Xc = clamp(yc, L, U); p[c] = yc - Xc;
// After:  yc = X[c] * p[c]; Xc = clamp(yc, L, U);
//         if (Xc > 0) p[c] = yc / Xc; else p[c] = 1.0;  // guard: Xc=0 when L=0 and X[c]=0
```

**Guard for Xc=0**: When `min_weight=0` (default) and a cell has zero total weight (structural zero), `yc=X[c]*p[c]=0`, `Xc=clamp(0,0,U)=0`. Without the guard, `p[c]=0/0=NaN` propagates silently. If `Xc ≤ 0`, set `p[c]=1.0` (multiplicative identity = no correction).

**Hyperplane correction** (`q_hyp`, init 1.0, shift→scale):
```cpp
// Before: X[c] += q_hyp; scale to sum=n; q_hyp = -shift;
// After:  X[c] *= q_hyp; scale to sum=n;
//         if (scale > 0) q_hyp = 1.0/scale; else q_hyp = 1.0;  // guard: scale=0 for empty data
```

**Guard for scale=0**: If all cells have zero weight (degenerate input), `scale=0`. Set `q_hyp=1.0` (no correction). Input validation in harvest.R should catch this upstream.

**Post-loop finalizer**: apply same multiplicative pattern (shown in raking.cpp lines 334-345).

**A5 fixture note**: Bregman Dykstra changes the fixed point — cell-table raking no longer matches the obs-level Euclidean reference. `data-raw/gen_raking_obs_ref.R` (confirmed present) must be rerun with the new algorithm after confirming AC4.

### 2. SOR (Successive Over-Relaxation) wiring

Wire existing `CalibState::sor_cfg` (type `CalibSorCfg`, line 95 of types.hpp, accessible from `raking_solve(CalibState& st)`) into raking's IPF marginal step.

Apply under-relaxation to scaling: `X[c] *= pow(T/S, eff_omega)` where `eff_omega ≤ 1`.

**"Bounds active" condition** (triggers SOR): `lo > 0.0 || hi < 1e300`. If bounds are absent (unconstrained IPF), oscillation doesn't occur and `eff_omega=1.0` (effectively disabled).

Update `harvest.R` docstring at line ~48: remove "iEPPA only" from `@param sor` description; add raking.

Same `sor=list(auto=TRUE, omega_min=0.3)` API as ieppa. No new types.

### 3. Greedy margin ordering

Wire existing `CalibState::scheduler` (type `SchedulerConfigLbw`, line 93 of types.hpp) into raking's K-margin sweep at line ~156 (`for (int k=0; k<st.K; k++)`).

**Prerequisite**: Raking currently computes only a scalar max `errRp` (no per-margin breakdown). Greedy requires a per-margin residual vector `errRp_k[K]`. Must add this vector computation before the Greedy sort (analogous to ieppa's `per_margin_err_prev`).

Same `scheduler="greedy"` API as ieppa. No new types.

---

## RED Test (TDD prerequisite)

Before implementing Bregman Dykstra, write a failing test:

```r
test_that("raking-bregman: Euclidean cycling breaks, Bregman converges on tight bounds", {
  # Synthetic: 1 margin, 2 categories, L=U=0.5 (equality bound)
  # Pure Euclidean Dykstra cycles here; Bregman should converge to exact target
  set.seed(1); n <- 100L
  df <- data.frame(a = factor(sample(c("1","2"), n, TRUE)))
  tgt <- list(a = c("1"=0.5, "2"=0.5))
  r <- leafblower::harvest(df, tgt, method="raking",
    min_weight=0.5, max_weight=0.5,  # equality bounds
    max_iterations=200L, attach_weights=FALSE)
  expect_equal(attr(r,"result")$status, 0L)
  expect_lt(attr(r,"result")$max_error, 1e-6)
})
```

With Euclidean Dykstra this may NOCONV (oscillate). With Bregman Dykstra it must converge to `max_err < 1e-6`. This is the RED test.

---

## A/B Test

Two configurations benchmarked on `stepstone-fulldata` (n=1.58M, K=9, max_weight=5):

| Config | Parameters | Label |
|--------|-----------|-------|
| **Approach A** | `method="raking", sor=list(auto=TRUE)` | Bregman + SOR |
| **Approach C** | `method="raking", sor=list(auto=TRUE), scheduler="greedy"` | Bregman + SOR + Greedy |

Metrics: max_err, marg_kl, weight_kl, wall, iters, DEFF, ESS.

Winner = lower max_err AND marg_kl (calibration quality primary). If tied, prefer A (simpler).

---

## Acceptance Criteria

1. **AC1**: `p[c]` init=1.0; box uses ratio `yc/Xc` with `Xc≤0` guard (set `p[c]=1.0`).
2. **AC2**: `q_hyp` init=1.0; hyperplane uses scale with `scale≤0` guard.
3. **AC3**: SOR fires when `sor_cfg.enabled && (lo>0||hi<1e300)` — eff_omega adapts per margin.
4. **AC4 (prerequisite for AC5)**: Winner config achieves max_err < 5.44e-3 on stepstone. **This must be verified BEFORE regenerating the A5 fixture.** If AC4 fails, halt and do not regenerate.
5. **AC5**: After AC4 confirmed, regenerate `raking_obs_reference_stepstone.rds` via `data-raw/gen_raking_obs_ref.R`. Cell-table raking then matches new obs-level reference to 1e-8.
6. **AC6**: Greedy: per-margin `errRp_k[K]` vector computed; margins sorted descending.
7. **AC7**: `harvest.R` `@param sor` docstring updated — no longer says "iEPPA only".
8. **AC8**: RED test passes after implementation; FAIL ≤ 2 (pre-existing, named: test-sor.R/lhs, test-ieppa-nonuniform-d.R).

---

## Files Changed

| File | Change |
|------|--------|
| `src/raking.cpp` | Bregman box+hyperplane+finalizer; `errRp_k[K]` vector; SOR wiring; Greedy wiring |
| `R/harvest.R` | Remove "iEPPA only" from `@param sor` docstring |
| `data-raw/gen_raking_obs_ref.R` | Rerun to regenerate fixture (after AC4 confirmed) |
| `tests/testthat/fixtures/raking_obs_reference_stepstone.rds` | Regenerated |
| `tests/testthat/test-calibration-solvers.R` | RED test (written before implementation) |

---

## Out of Scope

- P-A progressive bound tightening for raking
- Changes to ieppa, sinkhorn, or other solvers
- Python wrapper changes
