# iEPPA Speed, Convergence, and Bounds Hardening Design

**Status:** Draft rev 1

**Author:** Dennis Alexis Valin Dittrich
**Date:** 2026-04-24
**Triggers:**
- Assessment `docs/ieppa_assessmen.md` — identified four deviations from standard raking, two actionable: §1 obs-level bounds, §4 sum-to-N drift under capacity clamping.
- Post-WU-1..4 merge state: kk1204 per-iter ratio 2.17× raking; NOCONV at 500 iter on that regime. Stepstone errRp 2.21e-3 (best of 3 solvers).
- User priority: **C (speed) > B (convergence) > A (bounds)**.

**Baseline commit:** `1b29df1` (WU-4 structural vs transient infeas split).

**Scope.** Three passes in dependency order. Each pass independent-committable; later passes do not presuppose earlier-pass constants but share plumbing:

- **Pass 1 — Speed.** C1 fuse sweep+X_tilde+capacity into one outer-iter pass + C4 bump errRp interval 10→50. Target: kk1204 per-iter ratio ≤ 1.5× (down from 2.17×).
- **Pass 2 — Convergence.** B3 adaptive damping schedule (replace hard alpha=0.5 with `alpha(stress)`) + B1 Anderson(m=5) acceleration on `lf`/`f_lin` updates. Target: kk1204 reaches tol_abs=1e-3 within 500 iter (currently NOCONV).
- **Pass 3 — Bounds.** A1 new `bounds_mode` parameter ∈ {`cell` (default, current), `observation`} with A2 intra-cell water-filling behind `observation` mode. Target: when `bounds_mode="observation"`, every returned weight satisfies `w_i ∈ [min_weight, max_weight]` strictly.

## 1. Problem Statement

**(C) Speed.** Linear-path per-iter work at K=20, M_cell=1M: three O(M_cell) passes per outer iter (sweep sum, X_tilde rebuild, capacity block), each touching `X_cur` or derivatives. Measured 2.17× raking's single sweep. Three passes = 3 memory round-trips per iter; fusing into one pass with local temporaries cuts traffic to 1 round-trip. errRp check runs every 10 iters at O(M_cell), adding ~10% overhead; a 50-iter interval drops this to 2%.

**(B) Convergence.** Classical Sinkhorn is sublinear (rate ≈ iter^−0.5). kk1204 at K=20, dense compression, 500-iter budget reaches errRp ~2e-3 but still NOCONV at tol_abs=1e-3. Published acceleration (Chen-Lin-Jordan 2020; Mishchenko-Defazio 2020; Altschuler-Niles-Weed-Rigollet 2017) via Anderson mixing on the log-factor iterate closes 5–10× of iters-to-convergence on this regime class without changing fixed-point semantics.

**(A) Bounds.** Current code enforces `sum_{i∈c} w_i ∈ [L·|c|, U·|c|]` (cell aggregate). On cells with skewed base weights `d_i`, applying the uniform cell multiplier can produce individual `w_i > max_weight` or `w_i < min_weight`. Documented in spec §8 of prior design; test `test-ieppa-nonuniform-d.R` verifies the cell-aggregate invariant. Users who need strict per-observation bounds (statistical agency convention) have no route today.

## 2. Scope

Three passes, seven work units:

| Pass | WU | Addresses | Files | Priority |
|---|---|---|---|---|
| 1 | P1.1 | Fuse sweep + X_tilde + capacity | `src/ieppa.cpp` | P2 |
| 1 | P1.2 | Bump errRp interval 10 → 50 | `src/ieppa.cpp` | P3 |
| 2 | P2.1 | Adaptive damping schedule (replace hard 0.5) | `src/ieppa.cpp` | P2 |
| 2 | P2.2 | Anderson(m=5) acceleration on log-space path | `src/ieppa.cpp` | P1 |
| 2 | P2.3 | Anderson on linear-space f_lin | `src/ieppa.cpp` | P2 |
| 3 | P3.1 | `bounds_mode` parameter (C API + R/Python wrappers) | `src/ieppa.cpp`, `src/c_api.cpp`, `src/leafblower.h`, `R/harvest.R`, `python/leafblower/_harvest.py` | P2 |
| 3 | P3.2 | Intra-cell water-filling expansion | `src/ieppa.cpp` | P2 |

## 3. Non-scope

- No changes to raking or L-BFGS-B solvers.
- No changes to AUTO routing (`leafblower-3c2`).
- No new benchmark harness.
- Sum-to-N drift (§4 of assessment) is not addressed directly. Post-WU-4 behaviour: wrapper normalization `weights / mean(weights)` guarantees `sum(w) = n` to floating-point. Any drift inside the solver's iterates is absorbed by normalization. Tracking separately if evidence accumulates.

## 4. Pass 1 — Speed

### 4.1 P1.1 Fused outer-iter pass

**Current** (linear path): three passes over `c ∈ [0, M_cell)` per outer iter:
1. Sweep: per bucket (k, j) sum `S_kj = Σ X_cur[c] / f_kj_old`; per bucket update `f_lin[k][j]`; rescale `X_cur` for bucket cells.
2. Post-sweep: compute `X_tilde[c] = X_cur[c] / W[c]` (log path: `X_tilde[c] = exp(log_X_init + Σ lf[m])`).
3. Capacity: `xc = clamp(X_tilde[c], L, U); W[c] = xc/X_tilde[c]; X[c] = xc`; rebuild `X_cur[c] = X_tilde[c] * W[c]_new`.

Step 1 already writes `X_cur[c]` for every touched cell. Steps 2–3 overwrite `X_cur[c]` and `W[c]` again. Memory traffic: 3× `X_cur` reads + 3× writes per iter.

**Fused**: after step 1 (which leaves `X_cur` reflecting the final `f_lin` for this iter and the pre-capacity `W`), a single pass over `c` computes `X_tilde_c = X_cur[c] / W[c]`, clamps to `[L_cell[c], U_cell[c]]`, updates `W[c]` and `X[c]`, then rebuilds `X_cur[c] = X_tilde_c * W[c]` — all with one load of `X_cur[c]` and `W[c]`, and one store of each. Memory traffic: 1× round-trip. No change to math.

Log-path analog: step 2 computes `s = log_X_init[c] + Σ lf[m]` (K-loop, unavoidable); fuse step 3 inline using `X_tilde_c = exp(clip(s, kLogClip))` directly without a separate `X_tilde` store. `X_cur` is not used in log path, so no prefactored update needed. Fusion benefit smaller in log path (K-loop dominates), but non-zero.

Expected impact: 20–40% wall reduction on kk1204.

### 4.2 P1.2 errRp interval

Change `constexpr int kErrCheckInterval = 10` to `50` at default `inner_max_iter=500`. Scale: `max(10, inner_max_iter/10)`. errRp check is O(M_cell + Σ|bucket|) = O(K·M_cell). At 10-iter interval, amortizes to ~10% per-iter overhead on kk1204; at 50-iter interval, 2%.

Guard against over-coarse intervals on small budgets: `const int kErrCheckInterval = std::max(10, st.inner_max_iter / 10);`. Below budget 100, falls back to 10.

## 5. Pass 2 — Convergence

### 5.1 P2.1 Adaptive damping schedule

**Current** (WU-3): when `infeas_streak[idx] ≥ kInfeasPersistence/2 = 2` for any bucket, `alpha` latches at 0.5 for the remainder of the solve. Hard 0.5 halves convergence rate unconditionally after first stress, even if stress subsides.

**New**: `alpha = 1.0 / (1.0 + beta * stress)` where `stress = max_k,j infeas_streak[cat_offset[k]+j]` over all buckets with `target > 0`, `beta = 0.5`. At `stress=0`, `alpha=1.0` (fast path). At `stress=2`, `alpha=0.5`. At `stress=10`, `alpha=0.17`. Monotone, smooth, unlatched (alpha recovers as streaks reset).

Fast-path preserved: `if (stress == 0) { ... naive update ... }` branch avoids per-sweep recomputation of alpha.

WU-3's `force_damping_on` env var remains: when forced, `alpha = 0.5` constant regardless of streak.

Preserves Peyré-Cuturi 2019 Proposition 4.8 convergence guarantee (alpha ∈ (0, 1]). Rate slows by 1/alpha at any iter.

### 5.2 P2.2 Anderson acceleration on log-space `lf`

**Method.** After each full margin sweep (one `k=0..K` cycle), treat `lf` as the iterate in a fixed-point iteration `lf_{t+1} = G(lf_t)`. Anderson(m) mixing: maintain history of last `m+1` iterates `{lf_t, lf_{t-1}, ..., lf_{t-m}}` and residuals `{r_t, r_{t-1}, ..., r_{t-m}}` where `r_t = G(lf_t) - lf_t`. Next iterate: `lf_{t+1}_AA = lf_t + sum_j γ_j (G(lf_{t-j}) - lf_{t-j})` where γ solves a least-squares problem over residual differences.

**Parameter.** `m = 5`. History buffer: 2 × (m+1) × total_cats doubles. At total_cats=2000, ≈ 200KB per solver call — negligible.

**Restart condition.** If the least-squares system is ill-conditioned (condition number > 1e12 via thin QR), skip Anderson for this step and use the plain sweep. Prevents blowups in degenerate regimes.

**Interaction with damping.** Anderson applied AFTER damping: damping modifies the iterate `lf_{t+1}_plain = alpha * G(lf_t) + (1-alpha) * lf_t`; Anderson then mixes `lf_{t+1}_plain` with history. Damped-plain residual `r_t = lf_{t+1}_plain - lf_t = alpha * (G(lf_t) - lf_t)`.

**Engagement gate.** Anderson off for the first 5 iters (classical Sinkhorn preconditioner builds a useful starting point). After iter 5, always on unless condition-number check fails.

**Auto-disable flag.** `LBW_IEPPA_ACCEL ∈ {on, off, unset}` env var. Default `on`. Tests using both paths required.

Expected impact on kk1204: 5× iter reduction. At current per-iter ~40ms (post-P1.1), 500 iter → 100 iter = ~4s for kk1204 (was 20+s NOCONV). Still sublinear but practical.

### 5.3 P2.3 Anderson on linear-space `f_lin`

Parallel to 5.2 but on the linear iterate `f_lin`. Subtle: Anderson in linear space can produce non-positive iterates (coefficients γ can be negative). Guard: clamp `f_lin_AA[k][j] = max(f_lin_AA[k][j], kMinFLin)` where `kMinFLin = 1e-300`; below that threshold, revert to plain update for that cell (no Anderson contribution).

Alternatively: apply Anderson on `log(f_lin)` = `lf` equivalent, then exponentiate. This is numerically identical to P2.2 semantics. Preferred.

Implementation note: if the linear path is active, use `lf = log(f_lin)` as the Anderson state; after Anderson update, compute `f_lin = exp(lf_AA)` for the next sweep. Adds a `log`+`exp` pair per (k,j) per iter when Anderson fires — O(total_cats) = O(K·avg_cats) per iter, negligible vs sweep.

## 6. Pass 3 — Bounds

### 6.1 P3.1 `bounds_mode` parameter

**C API addition.** New field in `rk_params_t`:
```c
typedef enum { RK_BOUNDS_CELL = 0, RK_BOUNDS_OBSERVATION = 1 } rk_bounds_mode_t;
```
Add `rk_bounds_mode_t bounds_mode;` to `rk_params_t`. Default value via existing `rk_params_default()`: `RK_BOUNDS_CELL` (preserves current behaviour exactly).

**R wrapper.** `harvest(..., bounds_mode = c("cell", "observation"))`. Default `"cell"`. Maps to C enum.

**Python wrapper.** `harvest(..., bounds_mode: Literal["cell", "observation"] = "cell")`. Same mapping.

**ABI impact.** `rk_params_t` size grows by `sizeof(int)` = 4 bytes (or 8 after padding). Existing compiled callers that don't set the field get 0 = `RK_BOUNDS_CELL` = current behaviour. Backward compatible since the field is added at the end of the struct and zero-initialization is the current default.

### 6.2 P3.2 Intra-cell water-filling expansion

**Trigger.** Only when `bounds_mode == RK_BOUNDS_OBSERVATION`.

**Current expansion** (`src/ieppa.cpp` end of `ieppa_solve`):
```cpp
for (int c = 0; c < ct.M_cell; c++) mult[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
for (int i = 0; i < st.n; i++) st.weights[i] = st.weights[i] * mult[ct.cell_of[i]];
```

**New expansion** (observation mode). Per cell `c`:
1. Apply uniform multiplier `mult[c] = X[c] / X_init[c]` as before.
2. Clamp each `w_i` in cell `c` to `[min_weight, max_weight]`; accumulate `excess_above`, `excess_below`, total `free_mass` from non-clamped obs.
3. Redistribute `net_excess = excess_above − excess_below` proportionally over non-clamped obs in cell.
4. Repeat clamp+redistribute until no violations (bounded inner iter count = 50; exceeds → emit warning, return weights with remaining out-of-bound values clamped to nearest bound — distortion absorbed into that cell).

Per-cell inner loop converges in O(k) iters where k = number of distinct bound-hit patterns (typically ≤ 3 in practice). Worst case 50 iters * |cell| per cell = O(n * 50) total; n=10^6 → 5×10^7 ops, negligible against solver cost.

**Invariant.** After water-filling: `sum_{i∈c} w_i = X[c]` preserved (redistribution conserves mass), and `w_i ∈ [min_weight, max_weight]` for all i (up to clamp failure in pathological cells, warned).

**Implementation.** Group observations by cell ahead of time (`cells_of_obs[c]` list). Or iterate over `cells_of_obs[c] = [i : ct.cell_of[i] == c]` built in one O(n) pass. This index is built once pre-loop; memory O(n) ints.

**Corner case.** If cell has only one observation and its `w_i` is out of bounds, water-filling cannot redistribute. Fall back: clamp to nearest bound, emit warning, accept the cell-aggregate violation.

### 6.3 Tests

- `test-ieppa-bounds-mode.R`: when `bounds_mode="observation"` on highly skewed `d_i` input, assert `max(w) ≤ max_weight + 1e-9` and `min(w) ≥ min_weight - 1e-9` and `|sum(w) - n| < 1e-6`.
- Existing `test-ieppa-nonuniform-d.R`: unchanged behaviour when `bounds_mode="cell"` (default).

## 7. Testing

### 7.1 Speed regression (P1.1, P1.2)

- kk1204 per-iter ratio ≤ 1.5× (down from 2.17×). Script `/tmp/wu2_kk1204.R` reused.
- Full suite: FAIL=0, PASS ≥ 181 (no regressions).

### 7.2 Convergence (P2.x)

- kk1204 at tol_abs=1e-3, max_iter=500: expect RK_OK. Currently NOCONV.
- Existing dense-equivalence test (test-ieppa-faithful.R:109-136): `|res_lin - res_log| < 1e-8` must still hold. Anderson on both paths yields identical iterates.
- New: Anderson off vs on comparison. With `LBW_IEPPA_ACCEL=off`, iter count matches pre-P2 baseline ± 1. With `on`, at least 2× fewer iters on kk1204.
- Stepstone: errRp ≤ 2.21e-3 preserved (baseline from 1b29df1).

### 7.3 Bounds (P3.x)

- `test-ieppa-bounds-mode.R`: skewed-d input returns strict bounds under `bounds_mode="observation"`.
- Regression: default `bounds_mode="cell"` returns current behaviour (test-ieppa-nonuniform-d.R untouched).
- C API: zero-initialized `rk_params_t` yields `RK_BOUNDS_CELL`.
- Python: `bounds_mode` not passed → default cell behaviour.

## 8. Merge Gate

All of:
- kk1204 per-iter ratio ≤ 1.5× raking (post-P1).
- kk1204 reaches RK_OK at tol_abs=1e-3, max_iter=500 (post-P2).
- Stepstone errRp ≤ 2.21e-3 with no INFEAS (post-P2; must not regress).
- Strict obs-level bounds under `bounds_mode="observation"` on skewed-d input (post-P3).
- Default `bounds_mode="cell"` preserves current weight vectors to 1e-12 tolerance on all existing tests (backward compat).
- Full suite FAIL=0, PASS ≥ 181 at every intermediate commit.

## 9. Deliverables

Seven commits, one per WU (P1.1, P1.2, P2.1, P2.2, P2.3, P3.1, P3.2), in dependency order within each pass. Each commit includes:
- Source edit.
- Test (or test extension).
- All existing tests pass.

## 10. Non-Goals / Out of Scope

- Sum-to-N drift (assessment §4) — handled at wrapper level; not re-opening.
- Configurable `kInfeasPersistence` / exposing streak state to user API.
- Parallel within-margin scatter (B-axis C5) — scope creep; latent ticket only if benchmarks motivate.
- SIMD-specific vectorization (C3) — latent ticket; benchmark-driven only.
- CRAN / PyPI distribution changes.
