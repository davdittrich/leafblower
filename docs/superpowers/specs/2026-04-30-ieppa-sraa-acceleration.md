# iEPPA + SRAA-m (Anderson Acceleration) — Design Spec (rev 2)

## Overview

Apply Safeguarded Regularized Anderson Acceleration (SRAA-m) to the iEPPA and ieppa_soft solvers using the log-factor (`lf`) vector as the iterate. This extends the existing `accelerate=TRUE` parameter (currently a no-op for ieppa) to ieppa and ieppa_soft, achieving substantial iteration reduction on large-K or tight-bounds problems with minimal wall-time overhead.

---

## Problem Statement

iEPPA's outer loop is a Block Coordinate Descent fixed-point iteration on `lf[k][j]` (log-factors). For 9-margin surveys like Stepstone, it converges in ~50 outer iterations with `accelerate=FALSE`. On tighter problems (K≥6, max_weight<2, or skewed margins) it needs 100–500+ iterations. SRAA-m on raking achieves 2–5× iteration reduction there; ieppa has the same mathematical structure and should benefit similarly.

**Core insight:** iEPPA's true primal variable is `lf[k][j]` — the log-factor for margin k, category j. `lf` has dimension `n_cats_total = Σ_k cat_counts[k]` ≈ 50–500 for typical surveys. `X_cur[c]` is derived: `X_cur[c] = X_init[c] · exp(Σ_k lf[cat_offset[k] + g_per_cell[k][c]])`. Applying SRAA in lf-space means the history matrix is `m × n_cats_total ≈ 2250 doubles` vs `m × M_cell ≈ 5M doubles` in X_cur-space on a 1M-row survey. lf-space also preserves the fixed-point geometry cleanly because lf uniquely parameterizes the active state.

**lf vector layout note:** `lf` is a flat vector of size `n_cats_total_with_na = Σ_k (cat_counts[k]+1)` — each margin has a NA bucket (slot at cat_offset[k]+cat_counts[k]). NA slots are always zero and never updated by the sweep. pack/unpack operate on the full vector including NA slots; the inert zeros participate in the SRAA linear system but have zero contribution (ΔX=0, ΔR=0 for NA slots). This is harmless and avoids index complexity.

---

## Architecture

### Files changed

| File | Change |
|---|---|
| `src/sraa.hpp` | Add `apply_clamp = true` defaulted parameter to `sraa_step`; guard L/U access |
| `src/ieppa.cpp` | Add `pack_lf`, `unpack_lf` helpers; F_eval lambda; SRAA outer loop; `aa_accepted_count` field write |
| `src/ieppa.hpp` | Add `int aa_accepted_count = 0` field to `IEPPAResult` |
| `src/r_bridge.cpp` | Extend `accelerate_bool` to include ieppa/ieppa_soft; expose `aa_accepted_count` |
| `R/harvest.R` | Update `accelerate` docstring: list ieppa/ieppa_soft; note SOR disabled |

### New components (file-local to ieppa.cpp)

**`pack_lf(lf, dst)`** — O(n_cats_total_with_na): `std::copy(lf.begin(), lf.end(), dst.begin())`. Seeds SRAA iterate. `lf` is indexed as `lf[cat_offset[k]+j]` (same layout as f_lin).

**`unpack_lf(src, lf, f_lin, cell_lf, X_cur, ct, X_init, K, cat_offset)`** — O(n_cats_total + K·M_cell):
1. Copy src → lf; derive `f_lin[i] = exp(lf[i])` for i in 0..n_cats_total_with_na-1 (inv_f_old_lin is a scratch buffer re-derived inside the sweep lambdas from f_lin; unpack does NOT set it)
2. Rebuild `cell_lf[c] = Σ_k lf[cat_offset[k] + g_per_cell[k][c]]` for c in 0..M_cell-1 (g<0 skipped)
3. `X_cur[c] = X_init[c] * exp(cell_lf[c])` for c in 0..M_cell-1

Note: `inv_f_old_lin` is an inner-loop scratch buffer computed at the start of each `apply_single_margin_linear` call as `1/f_lin[off+j]`. It does not need to be set by `unpack_lf`.

**`f_eval_lf` lambda** — defined inside homotopy level loop, captures all sweep state by reference. Structure:
```cpp
auto f_eval_lf = [&](std::vector<double>& flat) -> double {
    unpack_lf(flat, lf, f_lin, cell_lf, X_cur, ct, X_init, st.K, cat_offset);
    bool overflow = false;
    for (int k = 0; k < st.K && !overflow; k++)
        overflow = apply_single_margin_linear(k);
    pack_lf(lf, flat);   // always pack — on overflow, lf is mid-sweep (partial)
    if (overflow) return std::numeric_limits<double>::infinity();
    // compute and return errRp from X_cur
    ...
};
```

**Overflow in F_eval:** When `apply_single_margin_linear` returns true (overflow), `lf` has been partially updated (some margins updated, rest not). `pack_lf` writes this partial state into `flat`. The SRAA safeguard compares `err_AA = inf > err_plain` → rejects, reverts `X` (= `lf_flat`) to the pre-step `F_cur` snapshot. This restores the consistent pre-overflow lf state. The next plain step will also hit overflow → log-fallback fires → SRAA history cleared.

**SRAA outer loop** — replaces the `for (iter_in_lvl=1; ...)` loop entirely when `st.accelerate && use_linear`. The for-loop is NOT preserved alongside SRAA; it is an either/or branch:

```cpp
if (st.accelerate && use_linear) {
    // --- SRAA accelerated path ---
    lbw::SRAAState ieppa_sraa;
    ieppa_sraa.init(n_cats_total_with_na, lbw::kSRAAm);
    std::vector<double> lf_flat(n_cats_total_with_na);
    std::vector<double> dummy_L, dummy_U;   // empty; not accessed with apply_clamp=false
    pack_lf(lf, lf_flat);

    int f_evals_used = 0;
    bool converged   = false;
    while (f_evals_used < budget_lvl && !converged) {
        ieppa_sraa.F_cur = lf_flat;   // seed F_cur before each step
        auto r = lbw::sraa_step(f_eval_lf, lf_flat, dummy_L, dummy_U,
                                ieppa_sraa, /*apply_clamp=*/false);
        f_evals_used += r.f_evals;
        res.base.iterations = f_evals_used;   // single write location; no conflict

        // best-iterate tracking (same pattern as raking.cpp:353)
        if (r.err_result < best_metric_seen) {
            best_metric_seen = r.err_result;
            W_best = X_cur;
            best_iter_val    = f_evals_used;
            best_objective_seen = compute_weight_kl(...);
        }
        // outer stall guard (mirrors raking.cpp:358-368)
        if (r.err_result > best_metric_seen * (1.0 + lbw::kSRAAOuterSlack)) {
            if (++sraa_outer_stall_count >= lbw::kSRAAOuterStallWindow) {
                X_cur = W_best;   unpack_lf_state_only(W_best, lf, f_lin, cell_lf);
                ieppa_sraa.clear(); sraa_outer_stall_count = 0;
            }
        } else { sraa_best_errRp = std::min(sraa_best_errRp, r.err_result);
                 sraa_outer_stall_count = 0; }
        // convergence check
        converged = lbw::check_convergence(st.convergence_cfg, r.err_result, ...);
    }
    res.aa_accepted_count = ieppa_sraa.aa_accepted_count;   // cumulative over this level
} else {
    // --- existing non-accelerated for-loop (unchanged) ---
    for (int iter_in_lvl = 1; iter_in_lvl <= budget_lvl; iter_in_lvl++) {
        ...
        res.base.iterations = total_iters + iter_in_lvl;
    }
}
```

`res.base.iterations` is set in exactly one place per path — no collision. When `accelerate=TRUE`, it reports F-evals consumed (AA step = 2, plain step = 1). When `accelerate=FALSE`, it reports outer BCD iterations as before. This matches the raking precedent.

---

## Interface Changes

### `src/sraa.hpp` — `apply_clamp` parameter with guard

```cpp
template<typename FEval>
SRAAStepResult sraa_step(
    FEval& f_eval,
    std::vector<double>& X,
    const std::vector<double>& L_cell,
    const std::vector<double>& U_cell,
    SRAAState& state,
    bool apply_clamp = true)   // NEW; default preserves all existing callers
{
    // ... steps 1-7 unchanged ...

    // Step 8: build extrapolated candidate
    for (int c = 0; c < M; c++) {
        double val = X[c] + Rk_c - corr;
        // CHANGED: guard L/U access behind apply_clamp
        state.scratch[c] = apply_clamp
            ? std::clamp(val, L_cell[c], U_cell[c])
            : val;
    }
    // ... steps 9-10 unchanged ...
}
```

When `apply_clamp=false`, `L_cell` and `U_cell` are never dereferenced — empty vectors are safe to pass.

Backward-compatible: all existing callers (raking.cpp, greenkhorn.cpp) omit the parameter → `apply_clamp=true` → behavior unchanged.

### `src/ieppa.hpp` — new field on IEPPAResult

```cpp
struct IEPPAResult {
    CalibResult base;
    // ... existing ieppa-private fields ...
    int  aa_accepted_count = 0;   // NEW: SRAA-m AA steps accepted this solve (0 if accelerate=FALSE)
    // ...
};
```

### `src/r_bridge.cpp`

```cpp
// Before:
accelerate_bool = isTRUE(accelerate) && method %in% c("raking", "greenkhorn");
// After:
accelerate_bool = isTRUE(accelerate) && method %in% c("raking", "greenkhorn",
                                                       "ieppa", "ieppa_soft");
```

Expose `aa_accepted_count` in the result list returned to R (alongside existing ieppa diagnostics).

### `R/harvest.R` — docstring additions

```r
#' @param accelerate Logical. Enable SRAA-m (Safeguarded Regularized Anderson
#'   Acceleration, window m=5) for \code{method="raking"}, \code{"greenkhorn"},
#'   \code{"ieppa"}, and \code{"ieppa_soft"}. Default \code{FALSE}.
#'
#'   \strong{When accelerate=TRUE for ieppa/ieppa_soft:}
#'   \itemize{
#'     \item SOR adaptive under-relaxation is disabled (omega fixed at omega_init=1.0).
#'       Use accelerate=FALSE if adaptive SOR is required.
#'     \item If scheduler="greedy" is also set, greedy is silently downgraded to
#'       round-robin (logged at verbose>=1). AA requires a deterministic sweep map.
#'     \item res$iterations reports F-evals consumed (AA step = 2, plain step = 1).
#'     \item res$aa_accepted_count reports how many AA steps were accepted.
#'     \item SRAA history is reset at each homotopy level boundary and on
#'       linear→log path fallback.
#'   }
```

---

## Behavioral Constraints

### SOR disabled when accelerate=TRUE (ieppa/ieppa_soft)

```cpp
const bool sor_auto_v = st.sor_cfg.auto_adapt && !st.accelerate;
```

Omega stays at `omega_init` (default 1.0). Documented in R docstring. Non-breaking: `accelerate=FALSE` preserves existing SOR behavior exactly.

### Greedy scheduler downgrades with log message

When `st.accelerate && scheduler == GREEDY`, round-robin is used. A message is emitted at verbose≥1:

```cpp
if (st.accelerate && st.scheduler.mode == SchedulerMode::GREEDY) {
    if (st.verbose >= 1)
        st.log("[ieppa] greedy scheduler disabled under SRAA-m; using round_robin");
    use_greedy_for_sraa = false;
}
```

### Log-path fallback: SRAA exit and result accounting

When linear overflow triggers `use_linear = false`:
1. `ieppa_sraa.clear()` — discard lf-space history
2. `res.aa_accepted_count = ieppa_sraa.aa_accepted_count` — snapshot count before clear (aa_accepted_count is NOT reset by clear(), so this captures the total from the linear phase)
3. Control falls through to the non-accelerated for-loop for the log path
4. `res.base.iterations` on log-path exit reflects the for-loop iter counter (not F-evals)

This is documented in the R docstring as "SRAA history reset on linear→log fallback."

### Alpha damping active inside F_eval

`alpha = compute_alpha()` is part of the sweep's fixed-point map. It remains active. The infeas_streak state is consistent on SRAA reversions because `ieppa_sraa.F_cur` snapshots the lf state before the step, and reverting `X` (= lf_flat) restores the pre-step lf. `infeas_streak` is not reverted — it is monotone and should not be.

### ALM outer loop (ieppa_soft) — combined code block

Both `lambda_cell` reset and SRAA clear must occur together:

```cpp
// When capacity_mu is updated in the ALM violation handler:
if (capacity_mu_updated) {
    // Reset ALM duals first, then reseed SRAA from the updated lf state
    std::fill(lambda_cell.begin(), lambda_cell.end(), 0.0);
    if (st.accelerate) {
        ieppa_sraa.clear();
        pack_lf(lf, lf_flat);   // reseed from current lf (post-dual-reset)
    }
}
```

Order matters: `lambda_cell` reset happens before `pack_lf` so the reseeded lf_flat reflects the state from which the new ALM level will begin.

`aa_accepted_count` for ieppa_soft reflects the cumulative count across all ALM levels (SRAAState.aa_accepted_count is not reset by `clear()`).

### n_cats_total invariant

`n_cats_total_with_na` is immutable for the lifetime of a homotopy level solve. pack_lf, unpack_lf, and SRAAState.init all use the same value. Validated once at SRAA init: `assert(lf.size() == n_cats_total_with_na)`.

---

## SRAA State Sizing

```
SRAAState.init(n_cats_total_with_na, kSRAAm=5)

Stepstone (K=9, avg 17 cats/margin, 1 NA/margin):
  n_cats_total_with_na = 9 × 18 = 162
  History: 2 × 5 × 162 × 8 = 12,960 bytes
  Gram: 5×5 × 8 = 200 bytes

X_cur-space equivalent: 2 × 5 × 800K × 8 ≈ 64 MB
```

---

## Testing

New file: `tests/testthat/test-ieppa-sraa.R`

| Test | Fixture | Assertion |
|---|---|---|
| Convergence parity | stepstone_small, ieppa | `max_error(accel) ≤ max_error(plain) + 1e-3` at same max_iterations |
| Iteration reduction | K=5 tight (max_weight=1.8, 5000 obs) | `res$iterations(accel) < 0.7 * res$iterations(plain)` (≥30% fewer F-evals) |
| ieppa_soft convergence | T5 fixture (5 cats, max_weight=1.8, bounds_mode="unit") | `max_error ≤ 0.01` with accelerate=TRUE |
| Greedy downgrade | any small fixture | no error; `res$status == 0`; verbose=1 message contains "round_robin" |
| Output correlation | well-conditioned K=2 fixture | `cor(w_plain, w_sraa) > 0.9999` |
| aa_accepted_count | tight K=5 fixture, 200 max_iters | `res$aa_accepted_count ≥ 5` |
| Regression: accel not worse at equal budget | tight K=5, max_iterations=50 | `max_error(accel) ≤ 2 × max_error(plain)` at same budget |
| Benchmark assertion | tight K=5 to tol=1e-4 | `f_evals(accel) < 0.7 × iters(plain)` to reach same tol |

The regression test (row 7) uses equal iteration budget; if SRAA diverges and the outer stall guard reverts to W_best, the returned weights must be the best-seen, not a degraded state.

---

## Open Questions / Non-goals

- **Log-path SRAA**: lf-space F_eval for the log path is feasible (same lf vector) but the log sweep has different convergence properties. Deferred to follow-on.
- **Greedy + SRAA**: compatible in theory; deferred after round-robin is validated.
- **SOR + SRAA composition**: possible if omega adaptation is aware of SRAA steps; deferred.
- **Cross-homotopy warm-starting**: SRAA history could be retained across levels with rescaling. Not in scope.
- **Per-ALM-level aa_accepted_count breakdown**: currently cumulative across ALM levels for ieppa_soft. Per-level tracking deferred.
