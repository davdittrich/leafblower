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

**Box correction** (`p[c]`, init 0→1, additive→multiplicative ratio).
Exact variable names: `L_cell[c]`, `U_cell[c]` (per-cell bounds at lines 99-100):
```cpp
// Before (Euclidean, line 179-181):
//   double yc = X[c] + p[c];
//   double Xc = std::clamp(yc, L_cell[c], U_cell[c]);
//   p[c] = yc - Xc;
// After  (Bregman):
//   double yc = X[c] * p[c];
//   double Xc = std::clamp(yc, L_cell[c], U_cell[c]);
//   p[c] = (Xc > 0.0) ? yc / Xc : 1.0;  // guard: Xc=0 when L_cell=0 and X[c]=0
```

**Guard for Xc=0**: When `min_weight=0` (default) and a cell has zero total weight (structural zero), `yc=X[c]*p[c]=0`, `Xc=clamp(0,0,U)=0`. Without the guard, `p[c]=0/0=NaN` propagates silently. If `Xc ≤ 0`, set `p[c]=1.0` (multiplicative identity = no correction).

**Hyperplane correction** (`q_hyp` scalar at line 105, init 1.0, shift→scale).
Hyperplane lambda `hyperplane_step()` at lines 141-147 must be replaced:
```cpp
// Before (Euclidean, hyperplane_step lines 143-146):
//   for (int c = 0; c < ct.M_cell; c++) { X[c] += q_hyp; s += X[c]; }
//   double shift = (n - s) / M_cell;
//   for (int c = 0; c < ct.M_cell; c++) X[c] += shift;
//   q_hyp = -shift;
// After  (Bregman):
//   for (int c = 0; c < ct.M_cell; c++) { X[c] *= q_hyp; s += X[c]; }
//   double scale = static_cast<double>(st.n) / s;
//   for (int c = 0; c < ct.M_cell; c++) X[c] *= scale;
//   q_hyp = (scale > 0.0) ? 1.0 / scale : 1.0;  // guard: scale=0 for degenerate input
```

**Post-loop finalizer** (lines 335-346) — apply same multiplicative pattern:
```cpp
// Before (Euclidean, lines 338-340):
//   double yc = X[c] + p[c]; Xc = clamp(yc, L_cell, U_cell); p[c] = yc - Xc;
//   + hyperplane_step() (additive)
// After  (Bregman):
//   double yc = X[c] * p[c]; Xc = clamp(yc, L_cell, U_cell);
//   p[c] = (Xc > 0.0) ? yc / Xc : 1.0;
//   + hyperplane_step() (multiplicative, updated above)
```

**Guard for scale=0**: If all cells have zero weight (degenerate input), `scale=0`. Set `q_hyp=1.0` (no correction). Input validation in harvest.R should catch this upstream.

**Post-loop finalizer**: apply same multiplicative pattern (shown in raking.cpp lines 334-345).

**A5 fixture note**: Bregman Dykstra changes the fixed point — cell-table raking no longer matches the obs-level Euclidean reference. `data-raw/gen_raking_obs_ref.R` (confirmed present) must be rerun:
```bash
# Run ONLY after AC4 is verified (max_err < 5.44e-3). If AC4 fails, HALT.
OMP_NUM_THREADS=1 Rscript data-raw/gen_raking_obs_ref.R
```
Verify: the generated fixture's `max_error` field ≤ 5.44e-3 (confirms AC4 holds for obs-level ref too).

### 2. SOR (Successive Over-Relaxation) wiring

Wire existing `CalibState::sor_cfg` (type `CalibSorCfg`, line 95 of types.hpp, accessible from `raking_solve(CalibState& st)`) into raking's IPF marginal step.

Apply under-relaxation to scaling: `X[c] *= pow(T/S, eff_omega)` where `eff_omega ≤ 1`.

**"Bounds active" condition** (triggers SOR): `st.min_weight > 0.0 || hi < 1e300` where local variable `hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300` at raking.cpp line 96. Both are accessible via `CalibState& st`. If bounds are absent (unconstrained IPF), oscillation doesn't occur and `eff_omega=1.0` (effectively disabled).

**S=0 guard in SOR marginal scaling**: In the IPF step, `X[c] *= pow(T/S, eff_omega)` where `S` = sum of weights in margin group. If `S=0` (structurally empty group), skip scaling for that group (`continue`). Add: `if (S <= 0.0) continue;` before the pow() call.

Update `harvest.R` docstring at **line 48** (confirmed exact):
```r
# Before: #' @param sor Named list for SOR adaptive under-relaxation (iEPPA only).
# After:  #' @param sor Named list for SOR adaptive under-relaxation (iEPPA and raking).
```

Same `sor=list(auto=TRUE, omega_min=0.3)` API as ieppa. No new types.

### 3. Greedy margin ordering

Wire existing `CalibState::scheduler` (type `SchedulerConfigLbw`, line 93 of types.hpp) into raking's K-margin sweep at line ~156 (`for (int k=0; k<st.K; k++)`).

**Prerequisite — per-margin residual vector**: Raking's marginal sweep (lines 158-171) already computes per-category `bucket[j]` sums for each margin `k`. Track per-margin `errRp_k[k]` DURING the sweep (not after), where those sums are available:

```cpp
// Declare before the outer iter loop (alongside existing declarations):
std::vector<double> errRp_k(st.K, 1.0 / st.K);  // init uniform; updated each iter
std::vector<int> order(st.K);
std::iota(order.begin(), order.end(), 0);

// Inside the K-margin sweep, replace `for (int k = 0; k < st.K; k++)` with:
// (when scheduler==GREEDY) sort order descending by errRp_k:
if (st.scheduler.mode == SchedulerMode::GREEDY)
    std::sort(order.begin(), order.end(),
              [&](int a, int b){ return errRp_k[a] > errRp_k[b]; });

for (int ki = 0; ki < st.K; ki++) {
    int k = (st.scheduler.mode == SchedulerMode::GREEDY) ? order[ki] : ki;
    // ... existing sweep body for margin k, which fills bucket[j] ...
    // After sweep, update per-margin residual:
    double ek = 0.0, W = 0.0;
    for (int j = 0; j < nj; j++) W += bucket[j];  // total
    for (int j = 0; j < nj; j++) {
        double e = std::fabs(bucket[j]/W - st.targets[k][j]);
        if (e > ek) ek = e;
    }
    errRp_k[k] = ek;
}
```

`bucket[j]` and `nj` are already computed within the sweep body — no extra cell-pass needed.

Same `scheduler="greedy"` API as ieppa. No new types.

---

## RED Test (TDD prerequisite)

The RED test must assert a PROPERTY that changes when Bregman Dykstra is implemented. Single-margin equality bounds are easy for Euclidean Dykstra and may not fail.

**Correct RED test** — on an unconstrained problem (no bounds), pure IPF = Bregman. Before Bregman fix, the Euclidean hyperplane correction (`q_hyp` additive) changes the fixed point vs pure IPF. After fix, raking without bounds must give the SAME weight_kl as ipfn on the same data.

```r
test_that("raking-bregman: unconstrained raking matches ipfn weight_kl (unified KL fixed point)", {
  # Without bounds, raking = pure IPF = same KL fixed point as ipfn.
  # Before Bregman fix: Euclidean hyperplane correction changes fixed point → different weight_kl.
  # After  Bregman fix: multiplicative hyperplane = identity for IPF → same weight_kl as ipfn.
  set.seed(1); n <- 2000L
  df <- data.frame(
    v1 = factor(sample(3L, n, TRUE)),
    v2 = factor(sample(2L, n, TRUE))
  )
  tgt <- list(v1=c("1"=0.5,"2"=0.3,"3"=0.2), v2=c("1"=0.6,"2"=0.4))

  r_raking <- leafblower::harvest(df, tgt, method="raking",
    min_weight=0, max_weight=Inf, max_iterations=500L, attach_weights=FALSE)
  r_ipfn   <- leafblower::harvest(df, tgt, method="ieppa",  # ieppa = pure IPF at cell level
    min_weight=0, max_weight=Inf, max_iterations=500L, attach_weights=FALSE)

  wkl_raking <- attr(r_raking,"result")$convergence_used$solver_objective
  wkl_ipfn   <- attr(r_ipfn, "result")$convergence_used$solver_objective

  # Before fix: wkl_raking ≠ wkl_ipfn (Euclidean hyperplane diverges from KL minimum)
  # After  fix: wkl_raking ≈ wkl_ipfn (same KL fixed point, tolerance 1e-4)
  expect_equal(wkl_raking, wkl_ipfn, tolerance=1e-4,
               label="unconstrained raking must reach same KL minimum as ipfn")
})
```

This is provably RED before implementation: the Euclidean additive hyperplane correction (`q_hyp -= shift`) changes the fixed point — it normalizes by subtraction instead of scaling, which is not the KL projection onto `{sum=n}`. The KL projection is multiplicative (`q_hyp = 1/scale`). Therefore, unconstrained raking with Euclidean Dykstra converges to a different fixed point than ieppa (pure multiplicative IPF), making the comparison fail.

**Field existence confirmed**: `convergence_solver_objective` is populated in `src/raking.cpp` at line 349 and exposed via `R/harvest.R` line 280 as `result$convergence_used$solver_objective`. The test will not ERROR on a NULL field.

**AC4→AC5 halting gate**: The spec mandates that fixture regeneration is a MANUAL step. An implementer must run the A/B benchmark, confirm AC4, then run `gen_raking_obs_ref.R`. This is not automated in CI. Add a comment in the CI test file: `# MANUAL GATE: regenerate fixture only after AC4 benchmark confirms max_err < 5.44e-3`.

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
