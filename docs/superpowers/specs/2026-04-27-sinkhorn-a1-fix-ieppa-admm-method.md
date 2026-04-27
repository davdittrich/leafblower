# Sinkhorn A1 Fix + iEPPA-ADMM Method + Decoupled Solver Objective

**Date**: 2026-04-27
**Status**: Pending design review
**Tickets**: leafblower-lig9
**Files**: `src/sinkhorn.cpp`, `src/ieppa.cpp`, `src/greg.cpp`, `src/grake.cpp`,
           `src/calib_dispatch.hpp`, `R/harvest.R`, `R/harvest_aliases.R` (or inline),
           `data-raw/gen_ieppa_kl_ref.R`

---

## Problem

### 1. A1 violated: sinkhorn weight_kl > ieppa weight_kl

A1 acceptance criterion: `sinkhorn_weight_kl_at_convergence ≤ ieppa_weight_kl_at_best_iter = 3.008e-3`.
Current benchmark: sinkhorn `convergence_used$objective` = 6.006e-3 > ieppa 4.97e-3.

Root cause: `res.convergence_objective = best_metric_seen` (sinkhorn.cpp line 183).
`best_metric_seen` = minimum of the **stopping criterion** (max_err for sinkhorn, since
max_err is the sinkhorn default). The field reports 6.006e-3 = min(max_err), not weight KL.

`weight_kl ≈ max_err = 6.006e-3` (4 sig fig match) is the smoking gun — two independent
quantities cannot coincidentally match unless one is being reported for the other.

This bug exists in all solvers where stopping criterion ≠ mathematical objective:
- ieppa: stopping=marginal_kl, objective=weight KL → reports marginal_kl (wrong)
- sinkhorn: stopping=max_err, objective=weight KL → reports max_err (wrong)
- greg: stopping=max_err(default), objective=chi² → reports max_err (wrong)

### 2. ieppa+ADMM has no explicit method name

ADMM capacity enforcement (Task 2) is always-on in `method="ieppa"`. No way to
explicitly request classic hard-clamp ieppa for comparison or regression testing.

---

## Design

### Principle: Correctness > Elegance > Simplicity. No backward compatibility needed.

---

## Part 1: Decouple Solver Objective from Stopping Criterion

Each solver maintains two independent trackers:

**`best_metric_seen`** (unchanged): minimum of stopping criterion metric over iterations.
Drives best_iter selection and convergence detection. Set by CalibMetric choice.

**`best_objective_seen`** (new): value of the solver's MATHEMATICAL OBJECTIVE at the
iteration selected by `best_metric_seen` (i.e., at `best_iter`). Always computed
regardless of stopping criterion.

```
solver          mathematical objective       field
─────────────── ─────────────────────────── ──────────────────────────────────
ieppa           weight KL                   m.kl (from compute_cell_metrics)
sinkhorn        weight KL                   m.kl
raking          weight KL                   m.kl
greg            chi²                        m.chi2
chebyshev       max_err                     m.errRp
grake           grake_norm                  m.grake_norm
```

**Implementation**: in each solver's error-check block, after computing `m` via
`compute_cell_metrics`, when `best_iter_val` is updated (stopping criterion minimum
achieved), also record:
```cpp
best_objective_seen = select_objective(solver_type, m);
// where select_objective picks m.kl for ieppa/sinkhorn/raking,
// m.chi2 for greg, m.errRp for chebyshev, m.grake_norm for grake
```

Add `best_objective_seen` to each solver result struct. Wire to
`res.convergence_objective = best_objective_seen`.

**sinkhorn.cpp additional change**: change default convergence metric to `"kl"` in
`R/harvest.R`. Sinkhorn's KL is monotone-decreasing, so `kl+improvement` is the
correct stopping criterion for its mathematical objective.

---

## Part 2: Add `method="ieppa_classic"` (pre-ADMM hard clamp)

`method="ieppa"` → ieppa+ADMM (current, best, keeps simple name).
`method="ieppa_classic"` → ieppa with original Euclidean hard clamp (no u[c] dual variable).

**Motivation:**
- `method="ieppa"` is strictly better on all metrics; simple name → best method (correct)
- `method="ieppa_classic"` needed for: regression testing, A1 fixture baseline,
  comparison studies, debugging ADMM convergence
- Consistent with codebase convention: add specific variant names when needed, keep
  generic name for best version

**Implementation**: The `u[c]` ADMM dual variable is always allocated and used in the
current ieppa.cpp. To support `ieppa_classic`, gate the ADMM update on an `IEPPAConfig`
flag:

```cpp
// In IEPPAConfig (types.hpp or ieppa.hpp):
bool use_admm_capacity = true;  // default: ADMM (ieppa); false = classic
```

In P1.1 block:
```cpp
if (cfg.use_admm_capacity) {
    double z = std::clamp(X_tilde_c + u[c], L_cell[c], U_cell[c]);
    u[c] += X_tilde_c - z;
    X[c] = z; W[c] = z / X_tilde_c; X_cur[c] = z;
} else {
    double xc = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);
    X[c] = xc; W[c] = xc / X_tilde_c; X_cur[c] = xc;
}
```

In `R/harvest.R`: add `"ieppa_classic"` to method dispatch, route to ieppa solver with
`use_admm_capacity = FALSE`. In `c_api.cpp` / `r_bridge.cpp`: pass flag through param struct.

`method="ieppa_classic"` uses the same `"marginal_kl"` default stopping criterion as
`method="ieppa"` (same solver, different capacity enforcement).

---

## Part 3: Regenerate A1 Fixture

The A1 fixture `tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds` was generated
pre-ADMM. With ADMM and marginal_kl stopping, ieppa's weight KL at best_iter changes.

`data-raw/gen_ieppa_kl_ref.R` must regenerate using the current `method="ieppa"` (ADMM)
and correctly extract weight KL from `result$convergence_used$objective` (which, after
Part 1 fix, will be weight KL at best_iter).

A1 test becomes:
```r
kl_s <- attr(w_s, "result")$convergence_used$objective  # sinkhorn weight KL at best_iter
kl_i <- ref$kl_at_best_iter                             # ieppa+ADMM weight KL at best_iter
expect_lte(kl_s, kl_i, label="sinkhorn weight KL ≤ ieppa+ADMM weight KL at best_iter")
```

---

## Acceptance Criteria

1. **A1 passes**: `sinkhorn_weight_kl ≤ ieppa_weight_kl_at_best_iter` — both values now
   correctly report weight KL (not stopping criterion values).
2. **Objective correctness**: `method="greg"` reports chi² in `convergence_used$objective`,
   not max_err. Verify for all solvers.
3. **Method dispatch**: `method="ieppa_classic"` runs without ADMM; benchmark shows
   max_err ≈ 2.74e-3 (pre-ADMM baseline). `method="ieppa"` shows max_err ≈ 2.34e-3 (ADMM).
4. **Sinkhorn default**: running `harvest(..., method="sinkhorn")` with no convergence
   arg uses `metric="kl"` internally.
5. **Regression**: All existing tests pass (FAIL 2 pre-existing unchanged).

---

## Files Changed

| File | Change |
|------|--------|
| `src/types.hpp` or `src/ieppa.hpp` | Add `bool use_admm_capacity = true` to IEPPAConfig |
| `src/ieppa.cpp` P1.1 block | Gate ADMM on `cfg.use_admm_capacity` |
| `src/ieppa.cpp` convergence block | Track `best_objective_seen` = weight KL at best_iter |
| `src/sinkhorn.cpp` line 183 | `res.convergence_objective = best_objective_seen` (weight KL) |
| `src/greg.cpp` | `res.convergence_objective = best_objective_seen` (chi²) |
| `src/grake.cpp` | `res.convergence_objective = best_objective_seen` (grake_norm) |
| `src/calib_dispatch.hpp` | Add `select_objective(solver_type, m)` helper |
| `R/harvest.R` | Add `"ieppa_classic"` dispatch; sinkhorn default `metric="kl"` |
| `src/c_api.cpp` / `src/r_bridge.cpp` | Pass `use_admm_capacity` flag |
| `data-raw/gen_ieppa_kl_ref.R` | Regenerate with ieppa+ADMM, extract weight KL |
| `tests/testthat/test-calibration-solvers.R` | Update A1 test to use new objective values |

---

## Out of Scope

- Python wrapper changes (add separately after R is stable)
- lbfgsb objective reporting (minor, separate ticket)
- ADMM convergence rate analysis
