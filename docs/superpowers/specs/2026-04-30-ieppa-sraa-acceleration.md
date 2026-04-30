# iEPPA + SRAA-m (Anderson Acceleration) — Design Spec

## Overview

Apply Safeguarded Regularized Anderson Acceleration (SRAA-m) to the iEPPA and ieppa_soft solvers using the log-factor (`lf`) vector as the iterate. This extends the existing `accelerate=TRUE` parameter (currently a no-op for ieppa) to ieppa and ieppa_soft, achieving substantial iteration reduction on large-K or tight-bounds problems with minimal wall-time overhead.

---

## Problem Statement

iEPPA's outer loop is a Block Coordinate Descent fixed-point iteration on `lf[k][j]` (log-factors). For 9-margin surveys like Stepstone, it converges in ~50 outer iterations with `accelerate=FALSE`. On tighter problems (K≥6, max_weight<2, or skewed margins) it needs 100–500+ iterations. SRAA-m on the raking solver (already implemented in `sraa.hpp`) achieves 2–5× iteration reduction there. The same mechanism should benefit ieppa.

The key design insight: iEPPA's true primal variable is **`lf[k][j]`** (n_cats_total ≈ 50–500 dimensions for typical surveys), NOT `X_cur` (M_cell ≈ 10K–1M). Applying SRAA in lf-space means the history matrix is `m × n_cats_total ≈ 2250 doubles` vs `m × M_cell ≈ 5M doubles` if applied in X_cur-space. The lf-space approach is also mathematically principled: lf uniquely parameterizes the ieppa state, while X_cur does not.

---

## Architecture

### Files changed

| File | Change |
|---|---|
| `src/sraa.hpp` | Add `apply_clamp = true` defaulted parameter to `sraa_step` |
| `src/ieppa.cpp` | Add `pack_lf`, `unpack_lf` helpers; F_eval lambda; SRAA outer loop |
| `src/r_bridge.cpp` | Extend `accelerate_bool` to include ieppa/ieppa_soft |
| `R/harvest.R` | Update `accelerate` docstring |

### New components

**`pack_lf(lf, dst)`** — O(n_cats_total): copy `lf[cat_offset[k]+j]` → `dst`. Seed SRAA iterate.

**`unpack_lf(src, lf, f_lin, inv_f_old_lin, cell_lf, X_cur, ...)`** — O(n_cats_total + K·M_cell): copy src → lf, derive f_lin/inv_f_old_lin via `exp`, rebuild `cell_lf[c] = Σ_k lf[cat_offset[k]+g_per_cell[k][c]]`, derive `X_cur[c] = X_init[c]·exp(cell_lf[c])`. The K·M_cell pass is same order as one sweep — acceptable.

**`f_eval_lf` lambda** — defined inside homotopy level loop. Calls `unpack_lf`, runs one round-robin linear sweep (`apply_single_margin_linear` for k=0..K-1), calls `pack_lf`, returns errRp from X_cur. On overflow: pack unchanged lf, return `infinity` (SRAA safeguard rejects and reverts to plain).

**SRAA outer loop** — replaces the `for (iter_in_lvl=1; ...)` loop when `st.accelerate && use_linear`:
```
while (f_evals_used < budget_lvl && !converged):
    seed F_cur = lf_flat
    r = sraa_step(f_eval_lf, lf_flat, {}, {}, sraa_state, apply_clamp=false)
    f_evals_used += r.f_evals
    [best-iterate tracking, outer stall guard, convergence check]
```

---

## Interface Changes

### `sraa.hpp`

```cpp
template<typename FEval>
SRAAStepResult sraa_step(
    FEval& f_eval,
    std::vector<double>& X,
    const std::vector<double>& L_cell,
    const std::vector<double>& U_cell,
    SRAAState& state,
    bool apply_clamp = true);   // NEW: false for lf-space (unconstrained iterate)
```

Backward-compatible: all existing callers (raking, greenkhorn) omit the parameter and get `apply_clamp=true`.

When `apply_clamp=false`, L_cell and U_cell may be empty — they are not accessed.

### `r_bridge.cpp`

```cpp
// Before:
accelerate_bool = isTRUE(accelerate) && method %in% c("raking", "greenkhorn");
// After:
accelerate_bool = isTRUE(accelerate) && method %in% c("raking", "greenkhorn",
                                                       "ieppa", "ieppa_soft");
```

---

## Behavioral Constraints

### SOR disabled when accelerate=TRUE

SOR adapts `sor_omega[k]` based on monotone convergence trajectory. SRAA's non-monotone extrapolation confuses the adaptive SOR. When `st.accelerate`:

```cpp
const bool sor_auto_v = st.sor_cfg.auto_adapt && !st.accelerate;
```

Omega stays at `omega_init` (default 1.0). Users requiring SOR should not set `accelerate=TRUE`.

### Greedy scheduler silently downgrades to round-robin

SRAA requires a deterministic fixed-point map F. Greedy margin ordering is state-dependent (it picks argmax errRp_k, which changes with X_cur). When `st.accelerate && scheduler == GREEDY`, the SRAA loop uses round-robin.

No error is emitted; behavior matches the raking precedent at raking.cpp:144-149.

### Alpha damping active inside F_eval

`alpha = compute_alpha()` runs at the top of each sweep, reading `infeas_streak`. It's part of the fixed-point map, not an external accelerator. It remains active. On a reverted (rejected) SRAA step, the post-revert lf state is consistent with the pre-step infeas_streak — alpha is recomputed correctly at the next F_eval.

### Log-path fallback clears SRAA history

When the linear path overflows and `use_linear` is set to false, ieppa switches to the log-space sweep — a different fixed-point map. SRAA history from the linear phase is stale:

```cpp
if (overflow_to_log) {
    ieppa_sraa.clear();
    pack_lf(lf, lf_flat);   // reseed
}
```

The non-accelerated for-loop handles the log path (SRAA only applies to linear path).

### Homotopy level boundary

SRAA history is cleared between homotopy levels (L_cell/U_cell change → different fixed-point). This already happens implicitly since the SRAA state is declared inside the homotopy level loop.

---

## ieppa_soft (ALM) Interaction

ieppa_soft runs inner BCD via the same sweep machinery as ieppa. The ALM outer loop:

1. Runs SRAA-accelerated inner BCD until inner convergence or budget
2. Checks capacity violation
3. If violation > `alm_violation_tol`: updates `capacity_mu` (adaptive) + resets `lambda_cell` duals
4. Go to 1

`capacity_mu` update changes the effective sweep objective → stale SRAA history. Clear on update:

```cpp
if (capacity_mu_updated && st.accelerate) {
    ieppa_sraa.clear();
    pack_lf(lf, lf_flat);
}
```

`f_eval_lf` already incorporates `st.alm_mu` and `lambda_cell` via capture — the ALM penalty is part of F. No change to F_eval.

---

## SRAA State Sizing

```
SRAAState.init(n_cats_total, kSRAAm=5)

n_cats_total = Σ_k cat_counts[k]

Stepstone (K=9, avg 16 cats/margin): n_cats_total ≈ 145
Memory: 2 × 5 × 145 × 8 = 11,600 bytes  (history buffers dX/dR)
Gram:   5×5 × 8 = 200 bytes

Contrast: X_cur space: 2 × 5 × 800K × 8 ≈ 64 MB
```

The lf-space SRAA state fits in a single cache line group vs megabytes for X_cur-space.

---

## Result Struct Updates

`res.iterations` reports F-evals consumed (analogous to raking's `f_evals_used`) when `accelerate=TRUE`. Each AA-accepted step counts 2 F-evals; each plain step counts 1.

Add `res.aa_accepted_count` to `IEPPAResult` (from `ieppa_sraa.aa_accepted_count`) so callers can diagnose SRAA effectiveness.

---

## F_eval Overflow Handling

The linear sweep can overflow (`new_f > kLinearOverflowTrip`). Inside F_eval:

- Pack unchanged lf into flat before returning
- Return `infinity`
- `sraa_step` receives `inf` from the second F_eval (AA attempt) → safeguard rejects → reverts X to pre-step state (which is the pre-overflow lf_flat)
- Next iteration: the plain step (1 F_eval) also overflows → ieppa's existing log-fallback fires → SRAA history cleared → non-accelerated log path takes over

This correctly handles the overflow→fallback transition within the SRAA framework.

---

## Testing

New file: `tests/testthat/test-ieppa-sraa.R`

| Test | What it checks |
|---|---|
| Convergence parity | `accelerate=TRUE` reaches same or better max_error as plain on stepstone_small |
| Iteration reduction | On tight-bounds K=5 problem, F-evals with AA < iters without AA |
| ieppa_soft convergence | T5 fixture (5 cats, max_weight=1.8) converges with ieppa_soft + accelerate |
| Greedy downgrade | scheduler="greedy" + accelerate=TRUE: no error, converges |
| Output correlation | cor(w_plain, w_sraa) > 0.9999 on well-conditioned problem |
| aa_accepted_count | res$aa_accepted_count > 0 when SRAA fires on a non-trivial problem |

---

## Open Questions / Non-goals

- **Log-path SRAA**: lf-space F_eval for the log path is feasible (same lf vector) but the log sweep has different convergence properties. Deferred to follow-on.
- **Greedy + SRAA**: compatible in theory (greedy is still a fixed-point map); deferred to follow-on after round-robin is validated.
- **SOR + SRAA composition**: possible if omega adaptation is aware of SRAA steps; deferred.
- **Cross-homotopy warm-starting**: SRAA history could be retained across levels with rescaling. Not in scope.
