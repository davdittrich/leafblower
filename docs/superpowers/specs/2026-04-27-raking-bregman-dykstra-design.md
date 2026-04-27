# Raking: Bregman Dykstra + SOR + Greedy A/B Test

**Date**: 2026-04-27
**Status**: Pending design review
**File**: `src/raking.cpp` only

---

## Problem

Raking's bounded calibration uses a **hybrid metric**: IPF marginal updates (multiplicative / KL-Bregman geometry) mixed with Euclidean additive Dykstra box and hyperplane projections. The code explicitly documents this at lines 13-16:

> "The composition of multiplicative (IPF) and Euclidean (Dykstra) projections is not covered by Boyle-Dykstra's single-metric theorem. Euclidean Dykstra corrections are NOT applied to the IPF step — they diverge against multiplicative updates."

Consequence: no unified convergence guarantee. On `stepstone-fulldata` (n=1.58M, K=9), raking achieves max_err=5.44e-3 vs ieppa's 2.34e-3 — a 2.3× gap. `ipfn` (unconstrained IPF, pure multiplicative) achieves 1.59e-3.

---

## Fix: Three independent improvements

### 1. Bregman/KL Dykstra (always-on, base fix)

Replace additive Euclidean corrections with multiplicative KL-Dykstra corrections. Unifies all projection steps under KL geometry — Bauschke & Lewis (2000) guarantees convergence to the KL minimum if a feasible solution exists.

**Box correction** (`p[c]`, initialized to 0 → 1):
```cpp
// Before (Euclidean):  p[c] = 0.0;  yc = X[c] + p[c]; Xc = clamp(yc, L, U); p[c] = yc - Xc;
// After  (Bregman):    p[c] = 1.0;  yc = X[c] * p[c]; Xc = clamp(yc, L, U); p[c] = yc / Xc;
```

**Hyperplane correction** (`q_hyp`, scalar shift → scale):
```cpp
// Before (Euclidean):  q_hyp = 0.0; X[c] += q_hyp; scale to n; q_hyp = -shift;
// After  (Bregman):    q_hyp = 1.0; X[c] *= q_hyp; scale to n; q_hyp = 1.0 / scale;
```

Post-loop finalizer: same multiplicative pattern.

**Note**: Changes the fixed point — cell-table raking no longer matches the obs-level Euclidean reference exactly. A5 fixture (`raking_obs_reference_stepstone.rds`) must be regenerated with the new Bregman Dykstra algorithm as baseline.

### 2. SOR (Successive Over-Relaxation) wiring

Wire existing `SorConfig` (from `CalibState`, already tested in ieppa) into raking's IPF marginal step. Apply under-relaxation to scaling:
```cpp
// Before: X[c] *= (T / S)
// After:  X[c] *= std::pow(T / S, eff_omega)  where eff_omega <= 1
```

Same `sor=list(auto=TRUE, omega_min=0.3)` API as ieppa. Reuses existing SOR adaptation logic — no new types.

**Condition**: SOR only if bounds are active (tight constraints cause oscillation; unconstrained IPF converges faster without SOR).

### 3. Greedy margin ordering

Wire existing `SchedulerMode::GREEDY` into raking's K-margin sweep. Sorts margins by residual `errRp_k` (descending) each outer iteration. Same `scheduler="greedy"` API as ieppa. O(K log K) overhead per sweep.

---

## A/B Test

Two configurations benchmarked on `stepstone-fulldata` (n=1.58M, K=9, max_weight=5):

| Config | Parameters |
|--------|-----------|
| **Approach A** | `method="raking", sor=list(auto=TRUE)` (Bregman + SOR) |
| **Approach C** | `method="raking", sor=list(auto=TRUE), scheduler="greedy"` (Bregman + SOR + Greedy) |

Metrics: max_err, marg_kl, weight_kl, wall, iters, DEFF, ESS.

Winner = better max_err AND marg_kl (calibration quality primary). If tied, prefer simpler (A).

---

## Acceptance Criteria

1. **AC1**: Bregman Dykstra: `p[c]` init=1.0, box uses ratio `yc/Xc`, hyperplane uses scale not shift.
2. **AC2**: SOR fires when `sor=list(auto=TRUE)` — omega adapts per margin.
3. **AC3**: Greedy fires when `scheduler="greedy"` — margins sorted by residual.
4. **AC4**: Winner config achieves max_err < 5.44e-3 (current raking baseline) on stepstone.
5. **AC5**: A5 fixture regenerated; cell-table raking matches new obs-level reference to 1e-8.
6. **AC6**: `devtools::test()` FAIL ≤ 2 (pre-existing).

---

## Files Changed

| File | Change |
|------|--------|
| `src/raking.cpp` | Bregman Dykstra (box + hyperplane + finalizer); SOR wiring; Greedy wiring |
| `data-raw/gen_raking_obs_ref.R` (or equivalent) | Regenerate A5 obs-level reference with Bregman Dykstra |
| `tests/testthat/fixtures/raking_obs_reference_stepstone.rds` | Regenerated fixture |

No API changes. `sor=` and `scheduler=` already accepted by harvest.R; just not wired to raking.

---

## Out of Scope

- Changes to ieppa, sinkhorn, or other solvers
- P-A progressive bound tightening for raking (separate ticket if needed)
- Python wrapper changes
