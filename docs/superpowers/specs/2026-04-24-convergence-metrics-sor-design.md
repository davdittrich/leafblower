# Convergence Criterion Reform + SOR Stabilization + Best-Iterate Design

**Date:** 2026-04-24
**Status:** Draft — pending design-review-gate

## Problem

iEPPA solver on stepstone-fulldata (1.58M rows, K=9, max_weight=5) diverges: errRp increases from ~2.22e-3 at iter 50 to ~6.5e-3 at iter 3000. Root causes:

1. **Wrong stopping criterion:** `errRp` (L∞ marginal error) oscillates in the capacity-clamp non-smooth regime. Solver declares NOCONV but reports the oscillated (worse) final iterate rather than the best-seen iterate.
2. **No oscillation damping:** Fixed multiplicative step size ω=1.0 amplifies oscillation across iterations.
3. **Incomplete quality reporting:** Only `max_error` returned; practitioners need `mean_error`, `KL`, `chi²` for cross-solver comparison and survey-statistics compliance.

## Goals

1. **Replace `max_err` as default convergence criterion** with weight-change `pct` — oscillation-resistant, industry-standard (GREG, calibrate, anesrake, autumn).
2. **Allow user-selectable convergence criteria:** `pct`, `max_err`, `mean_err`, `KL`, `chi2`.
3. **Always report all five quality metrics** at exit regardless of stopping criterion.
4. **Auto-enable SOR adaptive under-relaxation** (Thibault et al. 2021, Soma & Uschmajew 2024) to dampen oscillation in the margin update step.
5. **Always track best-iterate** (W at minimum observed errRp) and expose it alongside the final iterate.

## Non-goals

- Do not change the multiplicative IPF / iEPPA core algorithm structure.
- Do not break backward compatibility on well-behaved inputs (ω=1.0 when no oscillation, pct=0.001 declares convergence at same point max_err would on smooth problems).
- Do not implement full Polyak-Ruppert averaging (batched iterate mean) — best-iterate is sufficient and simpler.
- Do not expose SOR ω directly on the C ABI (internal implementation detail).

## Design

### §1 Convergence Criterion System

**New default:** `convergence = list(pct = 0.001)` — stop when `max_i |W_i[t] - W_i[t-1]| / W_i[t-1] < 0.001` (0.1% max proportional weight change).

**Pluggable criterion:**

```r
convergence = list(pct = 0.001)                         # default
convergence = list(absolute = 1e-6)                     # backward compat (max_err)
convergence = list(absolute = 1e-3, criterion = "chi2")
convergence = list(absolute = 1e-3, criterion = "mean_err")
convergence = list(absolute = 1e-3, criterion = "KL")
convergence = list(pct = 0.001, absolute = 1e-6, combine = "or")  # OR logic
```

**`combine`:** relevant only when BOTH `pct` and `absolute` are non-zero. Default `"or"` (stop when EITHER fires). `"and"` available for strict mode (both must be satisfied). When only one threshold is set, `combine` is ignored.

**New `CalibConvergenceCfg` struct fields (added to `CalibState`):**

| field | type | default | meaning |
|---|---|---|---|
| `pct_tol` | double | 0.001 | weight-change pct threshold (0 = off) |
| `absolute_tol` | double | 0.0 | marginal error absolute threshold (0 = off) |
| `criterion` | enum | PCT | PCT / MAX_ERR / MEAN_ERR / KL / CHI2 |
| `combine` | enum | OR | OR / AND |

**pct computation** (inside inner loop, after capacity block, before convergence check):
```cpp
double pct_change = 0.0;
for (int c = 0; c < M_cell; c++) {
    double delta = std::fabs(W[c] - W_prev[c]) / std::max(W_prev[c], 1e-12);
    if (delta > pct_change) pct_change = delta;
}
// Store W_prev = W after check
```

**Alternative criteria:**

- **max_err:** current `errRp` (unchanged formula)
- **mean_err:** `mean_k(max_j |Ŝ_kj/n - T_kj|)` — L1-over-margins, L∞-within-margin
- **KL:** `max_k Σ_j T_kj * log((T_kj + ε) / (Ŝ_kj/n + ε))`, ε = 1e-10
- **chi2:** `Σ_k Σ_j (Ŝ_kj - T_kj*n)² / (T_kj*n + ε)`, ε = 1

All alternative criteria computed at the same `kErrCheckInterval` cadence as `errRp`.

**Backward compatibility:** `convergence = list(absolute = X)` without `criterion` = `{criterion=MAX_ERR, absolute_tol=X, pct_tol=0}` — exact current behavior.

---

### §2 Quality Metric Reporting

Always computed at exit. Added to the `C_rk_calibrate` return list (C ABI + R bridge extension):

| R field | metric | formula |
|---|---|---|
| `max_error` | L∞ marginal | unchanged `errRp` |
| `mean_error` | L1-over-margins | `mean_k(max_j |Ŝ_kj/n - T_kj|)` |
| `kl_divergence` | max KL | `max_k Σ_j T_kj log((T_kj + ε)/(Ŝ_kj/n + ε))` |
| `chi_square` | total χ² | `Σ_k Σ_j (obs-exp)²/(exp + ε)` |
| `pct_change` | weight pct | last iteration's `max_i |ΔW_i/W_i|` |

Computed from final iterate (or best iterate if different — both computed). `DEFF` and `ESS` remain R-side helpers using the returned weights.

---

### §3 SOR Adaptive Under-Relaxation (Rank 1)

**Theory:** Thibault–Chizat–Dossal–Papadakis 2017 (HAL/Algorithms 14(5):143, 2021); Soma–Uschmajew arxiv:2410.14104 (2024). Replace Sinkhorn step `f_k_new[j] = f_k_old[j] * (T_kj / Ŝ_kj)` with:

```cpp
f_k_new[j] = f_k_old[j] * std::pow(T_kj / S_kj, omega[k]);
// omega[k] < 1: under-relaxation (suppresses oscillation)
// omega[k] = 1: standard Sinkhorn (current behavior)
```

**Default: `sor = list(auto = TRUE, omega_min = 0.3)`**

Auto mode mechanics per margin k:
- `omega[k]` initialized to 1.0
- Burn-in: first `sor_burnin = 20` outer iterations, ω unchanged
- After burn-in: track `prev_errRp_k` per margin
  - Sign flip (`errRp_k[t] > errRp_k[t-1]` after prior decrease): `omega[k] = max(omega_min, omega[k] * 0.7)`
  - Monotone decreasing: `omega[k] = min(1.0, omega[k] * 1.05)` (slow recovery)

**Performance:** `std::pow(x, 1.0)` = x (compiler trivially elides). Zero overhead when ω=1.0 (no oscillation). `omega[]` array = K doubles (K=9 for stepstone-fulldata, negligible).

**R API:**

```r
sor = list(auto = TRUE, omega_min = 0.3)   # default
sor = list(omega = 0.5)                    # fixed uniform ω (same value, no per-margin adaptation)
sor = list(enabled = FALSE)                # disable
```

**New `CalibSorCfg` struct:**

| field | type | default |
|---|---|---|
| `auto` | bool | true |
| `omega_init` | double | 1.0 |
| `omega_min` | double | 0.3 |
| `omega_fixed` | double | -1 (auto) |
| `burnin` | int | 20 |

SOR diagnostics reported: `calib_result$sor_min_omega` (lowest ω reached across any margin) and `calib_result$sor_n_damped` (count of margin×iteration pairs where ω < 1).

---

### §4 Best-Iterate Tracking (Rank 3)

Always enabled. Negligible overhead: copies O(M_cell) doubles only when errRp improves (monotone copy-on-improvement).

```cpp
// In inner loop, after errRp computed:
if (errRp < best_errRp) {
    best_errRp = errRp;
    W_best     = W;       // copy M_cell doubles
    best_iter  = iter;
}
// At exit:
// 1. expand W_best to obs-level (same logic as W_final)
// 2. attach as separate output
```

**R API:**

```r
result <- harvest(data, target, ...)
# Final iterate (unchanged default return):
w_final <- if (attach_weights) result$weights else as.numeric(result)
# Best iterate:
w_best  <- attr(result, "best_weights")   # always populated
# Diagnostics:
attr(result, "result")$best_error  # errRp at best iterate
attr(result, "result")$best_iter   # which iteration number
```

For stepstone-fulldata: `best_iter ≈ 50`, `best_error ≈ 2.22e-3` — retrieved without convergence criterion change.

**Note:** When SOR (§3) also active, `best_weights` reflect W at the iteration where SOR + convergence together achieved the minimum errRp.

---

### §5 Files Affected

| File | Change |
|---|---|
| `src/types.hpp` | Add `CalibConvergenceCfg`, `CalibSorCfg`; attach to `CalibState`; add `W_prev` (initialized to start_weights at solver entry), `W_best` (initialized to W at iter=1), `omega[]` (K doubles, stack-allocated in solver, initialized to omega_init) scratch fields |
| `src/leafblower.h` | Extend C ABI: convergence config struct, SOR config struct, new result fields |
| `src/ieppa.cpp` | Add pct computation; add alternative criteria computation; add SOR ω-adaptation; add best-iterate tracking |
| `src/raking.cpp` | Add pct computation + best-iterate tracking (same pattern) |
| `src/c_api.cpp` | Marshal new config structs; return new result fields |
| `src/r_bridge.cpp` | Unpack new args; pack new result fields |
| `R/harvest.R` | Accept `convergence$criterion`, `convergence$combine`, `convergence$pct`, `sor`; attach `best_weights` attribute |
| `python/leafblower/_harvest.py` | Mirror R API changes |
| `python/leafblower/_bindings.cpp` | Mirror C ABI changes |
| `tests/testthat/test-convergence-criteria.R` | New: pct criterion, alternative criteria, pct vs max_err equivalence on smooth inputs |
| `tests/testthat/test-sor.R` | New: SOR auto disabled on smooth inputs (ω stays 1.0); SOR dampens on oscillatory fixture |
| `tests/testthat/test-best-iterate.R` | New: best_weights populated; best_iter correct; best_error ≤ max_error at exit |
| `man/harvest.Rd` | Updated via roxygen |

---

### §6 Acceptance Criteria

**A1:** `convergence = list(absolute = X)` (no criterion field) = exact current behavior. Existing 232 tests green.

**A2:** `pct = 0.001` default: on a smooth synthetic input, solver converges at same (±5%) iteration count as current `absolute = 1e-6` max_err.

**A3:** SOR auto: on the tight-clamp synthetic fixture (max_weight=2, K=5, n=5000), ω drops below 0.9 for at least one margin → confirmed oscillation detected.

**A4:** SOR auto: on a smooth synthetic input (max_weight=10, well-separated targets), `sor_min_omega = 1.0` (no damping triggered). Zero performance regression.

**A5:** `best_weights` populated in all harvest() calls. `best_error ≤ max_error` at exit always.

**A6:** For stepstone-fulldata (from saved benchmark report): `attr(result, "best_weights")` corresponding `best_error ≈ 2.22e-3` (not the diverged 6.5e-3 final iterate).

**A7:** All five quality metrics (`max_error`, `mean_error`, `kl_divergence`, `chi_square`, `pct_change`) present in `attr(result, "result")` for every call.

**A8:** `devtools::test()` FAIL 0. `R CMD check --as-cran` 0 ERROR, 0 WARNING.

---

### §7 Open Questions (resolved by design)

- **pct default value:** 0.001 (0.1% max weight change). Tighter than anesrake (0.01); matches precision-calibration use.
- **chi2 threshold:** User-supplied; no default (must pass `absolute = X` with `criterion = "chi2"`). No implied n-dependence.
- **KL ε:** 1e-10 (numerical floor); prevents log(0) on empty cells; small enough not to distort KL on non-empty cells.
- **Best-iterate storage:** Cell-level (M_cell), not obs-level, during iteration. Expanded at exit. Prevents 1.58M × T memory.
- **SOR + raking:** SOR applied to iEPPA only in first implementation. Raking convergence is already guaranteed (smooth operators); SOR unnecessary there but `pct` criterion + best-iterate tracking apply to raking too.
