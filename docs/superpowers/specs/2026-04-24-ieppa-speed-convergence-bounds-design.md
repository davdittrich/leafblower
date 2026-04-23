# iEPPA Speed, Convergence, and Bounds Hardening Design

**Status:** Draft rev 2 (post iter-1 PM + Architect adversarial review)

**Rev 2 changes:**
- Dropped P1.2 (errRp interval bump) per user directive.
- §4.1: corrected "3 passes → 1" framing to "2 post-sweep passes → 0"; sweep already maintains `X_cur` per WU-2.
- §5.2: Anderson + capacity BCD interaction — added engagement gate `n_cap_active == 0`; mixing lf-residuals across iters with different W breaks the contraction argument otherwise.
- §5.3 P2.3 merged into P2.2 — PM+Architect both flagged redundancy; log-domain Anderson is a dispatch branch inside P2.2, not a separate WU. Total WUs: 4.
- §6.1: added `memset(p, 0, sizeof(*p))` requirement in `rk_params_init()` to prevent uninitialized `bounds_mode` in stack-allocated structs. Added `static_assert` on sizeof(rk_params_t) as documentation anchor for ABI size change.
- §6.1: reframed "statistical agency convention" claim; no concrete citation, softened to "advanced use case."
- §7: per-pass stepstone regression checkpoint added after each commit.

**Author:** Dennis Alexis Valin Dittrich
**Date:** 2026-04-24
**Triggers:**
- Assessment `docs/ieppa_assessmen.md` — identified four deviations from standard raking, two actionable: §1 obs-level bounds, §4 sum-to-N drift under capacity clamping.
- Post-WU-1..4 merge state: kk1204 per-iter ratio 2.17× raking; NOCONV at 500 iter on that regime. Stepstone errRp 2.21e-3 (best of 3 solvers).
- User priority: **C (speed) > B (convergence) > A (bounds)**.

**Baseline commit:** `1b29df1` (WU-4 structural vs transient infeas split).

**Scope.** Three passes in dependency order. Each pass independent-committable; later passes do not presuppose earlier-pass constants but share plumbing:

- **Pass 1 — Speed.** C1 fuse sweep+X_tilde+capacity into one outer-iter pass. Target: kk1204 per-iter ratio ≤ 1.5× (down from 2.17×).
- **Pass 2 — Convergence.** B3 adaptive damping schedule (replace hard alpha=0.5 with `alpha(stress)`) + B1 Anderson(m=5) acceleration on `lf`/`f_lin` updates. Target: kk1204 reaches tol_abs=1e-3 within 500 iter (currently NOCONV).
- **Pass 3 — Bounds.** A1 new `bounds_mode` parameter ∈ {`cell` (default, current), `observation`} with A2 intra-cell water-filling behind `observation` mode. Target: when `bounds_mode="observation"`, every returned weight satisfies `w_i ∈ [min_weight, max_weight]` strictly.

## 1. Problem Statement

**(C) Speed.** Linear-path per-iter work at K=20, M_cell=1M: three O(M_cell) passes per outer iter (sweep sum, X_tilde rebuild, capacity block), each touching `X_cur` or derivatives. Measured 2.17× raking's single sweep. Three passes = 3 memory round-trips per iter; fusing into one pass with local temporaries cuts traffic to 1 round-trip. errRp check runs every 10 iters at O(M_cell), adding ~10% overhead; a 50-iter interval drops this to 2%.

**(B) Convergence.** Classical Sinkhorn is sublinear (rate ≈ iter^−0.5). kk1204 at K=20, dense compression, 500-iter budget reaches errRp ~2e-3 but still NOCONV at tol_abs=1e-3. Published acceleration (Chen-Lin-Jordan 2020; Mishchenko-Defazio 2020; Altschuler-Niles-Weed-Rigollet 2017) via Anderson mixing on the log-factor iterate closes 5–10× of iters-to-convergence on this regime class without changing fixed-point semantics.

**(A) Bounds.** Current code enforces `sum_{i∈c} w_i ∈ [L·|c|, U·|c|]` (cell aggregate). On cells with skewed base weights `d_i`, applying the uniform cell multiplier can produce individual `w_i > max_weight` or `w_i < min_weight`. Documented in spec §8 of prior design; test `test-ieppa-nonuniform-d.R` verifies the cell-aggregate invariant. Users who need strict per-observation bounds (statistical agency convention) have no route today.

## 2. Scope

Three passes, four work units. Within-pass order below is the mandated execution order regardless of priority label.

| Pass | WU | Addresses | Files | Priority |
|---|---|---|---|---|
| 1 | P1.1 | Eliminate post-sweep X_tilde + X_cur rebuild loops; fuse capacity inline | `src/ieppa.cpp` | P2 |
| 2 | P2.1 | Adaptive damping schedule (replace hard 0.5) | `src/ieppa.cpp` | P2 |
| 2 | P2.2 | Anderson(m=5) acceleration on `lf` (both paths via log-domain dispatch) | `src/ieppa.cpp` | P2 |
| 3 | P3.1 | `bounds_mode` parameter + water-filling expansion | `src/ieppa.cpp`, `src/c_api.cpp`, `src/leafblower.h`, `R/harvest.R`, `python/leafblower/_harvest.py` | P2 |

## 3. Non-scope

- No changes to raking or L-BFGS-B solvers.
- No changes to AUTO routing (`leafblower-3c2`).
- No new benchmark harness.
- Sum-to-N drift (§4 of assessment) is not addressed directly. Post-WU-4 behaviour: wrapper normalization `weights / mean(weights)` guarantees `sum(w) = n` to floating-point. Any drift inside the solver's iterates is absorbed by normalization. Tracking separately if evidence accumulates.

## 4. Pass 1 — Speed

### 4.1 P1.1 Eliminate post-sweep X_tilde + X_cur rebuild

**Current linear-path structure** (`src/ieppa.cpp` HEAD 1b29df1). Per outer iter:
1. **Sweep** (two O(M_cell) passes + per-bucket sum/rescale) leaves `X_cur[c] = X_init[c] · W[c] · ∏_m f_lin[m][g_m(c)]` current as of the last f_lin update.
2. **X_tilde compute** (separate O(M_cell) pass): `X_tilde[c] = X_cur[c] / W[c]` (linear) or `exp(log_X_init + Σ lf[m])` (log).
3. **Capacity block** (O(M_cell)): `xc = clamp(X_tilde[c], L, U); W[c] = xc/X_tilde[c]; X[c] = xc`.
4. **X_cur rebuild** (O(M_cell)): `X_cur[c] = X_tilde[c] * W[c]_new` (linear path only).

**Fusion.** Steps 2+3+4 collapse to a single O(M_cell) pass. Invariant: after the sweep, `X_cur[c] / W[c]` equals the correct pre-capacity X_tilde by construction (sweep maintains `X_cur` proportional to `W · ∏f_lin`). So:

```cpp
// Linear path: fused post-sweep block
for (int c = 0; c < ct.M_cell; c++) {
    if (X_init[c] <= 0.0 || W[c] <= 0.0) {
        X[c] = 0.0; X_cur[c] = 0.0; continue;
    }
    double X_tilde_c = X_cur[c] / W[c];          // implicit, no store
    double xc = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);
    X[c] = xc;
    W[c] = xc / X_tilde_c;                       // new W
    X_cur[c] = xc;                               // = X_tilde_c * W_new; one store
    if (xc != X_tilde_c) n_cap++;
}
```

Log-path analog: step 2 computes `s = log_X_init[c] + Σ lf[m]` (K-loop, unchanged); fuse capacity inline using `X_tilde_c = exp(clip(s, kLogClip))` directly without a separate `X_tilde` store. No X_cur in log path.

Memory traffic reduction: linear path 3× (X_cur, W, X_tilde) round-trips → 1×. K-loop in log path unchanged so fusion gain smaller there.

**Framing correction (iter-1 Architect).** Prior rev 1 claimed "3 passes → 1"; the sweep itself is not fused (it remains the dominant pass). What P1.1 actually eliminates is the 2 post-sweep passes (X_tilde + X_cur rebuild) by computing capacity inline using the sweep's final `X_cur`. Framed as "2 post-sweep passes → 0."

Expected impact: 20–40% wall reduction on kk1204 linear path; marginal on stepstone log path (K-loop dominates).

## 5. Pass 2 — Convergence

### 5.1 P2.1 Adaptive damping schedule

**Current** (WU-3): when `infeas_streak[idx] ≥ kInfeasPersistence/2 = 2` for any bucket, `alpha` latches at 0.5 for the remainder of the solve. Hard 0.5 halves convergence rate unconditionally after first stress, even if stress subsides.

**New**: `alpha = 1.0 / (1.0 + beta * stress)` where `stress = max_k,j infeas_streak[cat_offset[k]+j]` over all buckets with `target > 0`, `beta = 0.5`. At `stress=0`, `alpha=1.0` (fast path). At `stress=2`, `alpha=0.5`. At `stress=10`, `alpha=0.17`. Monotone, smooth, unlatched (alpha recovers as streaks reset).

Fast-path preserved: `if (stress == 0) { ... naive update ... }` branch avoids per-sweep recomputation of alpha.

WU-3's `force_damping_on` env var remains: when forced, `alpha = 0.5` constant regardless of streak.

Preserves Peyré-Cuturi 2019 Proposition 4.8 convergence guarantee (alpha ∈ (0, 1]). Rate slows by 1/alpha at any iter.

### 5.2 P2.2 Anderson(m=5) acceleration on `lf` (both paths via log-domain dispatch)

**Method.** After each full margin sweep (one `k=0..K` cycle) followed by the capacity block, treat `lf` as the iterate in a fixed-point iteration `lf_{t+1} = G(lf_t)`. Anderson(m) mixing: maintain history of last `m+1` iterates `{lf_t, lf_{t-1}, ..., lf_{t-m}}` and residuals `{r_t, r_{t-1}, ..., r_{t-m}}` where `r_t = G(lf_t) - lf_t`. Next iterate: `lf_{t+1}_AA = lf_t + sum_j γ_j (G(lf_{t-j}) - lf_{t-j})` where γ solves a least-squares problem over residual differences.

**Both paths via log-domain dispatch.** Anderson always operates on `lf`. In the linear path, maintain a shadow `lf = log(f_lin)` state: after each sweep compute `lf[kj] = log(f_lin[kj])` for updated buckets (O(total_cats)); apply Anderson; then materialize `f_lin[kj] = exp(lf_AA[kj])` (O(total_cats)). The `log+exp` per (k,j) per iter is O(total_cats) = O(K·avg_cats), dominated by the O(K·M_cell) sweep at M_cell ≫ avg_cats. Architect iter-1 flagged this as a possible offset to P1.1 savings — empirical check required (merge-gate budget section §7.1).

**Parameter.** `m = 5`. History buffer: 2 × (m+1) × total_cats doubles. At total_cats=2000, ≈ 200KB per solver call — negligible.

**Capacity-block interaction (iter-1 Architect BLOCKER).** Each outer iter is: sweep → fused capacity block (P1.1). Capacity updates `W[c]`; the next sweep operates with new `W`, so the effective operator `G` is not stationary across iters when `n_cap_active > 0`. Anderson history mixing residuals computed under different `W` breaks the contraction argument and can diverge.

**Engagement gate.** Anderson **only fires on outer iters where `n_cap_active == 0` at the end of the capacity block.** When any cell is capacity-active, skip Anderson for that step: use the plain (damped) iterate, clear Anderson history buffer, and resume history accumulation from the next uncapacitated iter. In the kk1204 regime at `max_weight=3`, capacity tends to be saturated early then release; once all cells are unclamped, Anderson kicks in and carries the iteration. This is a conservative restriction; later ticket can explore joint `(lf, log W)` Anderson if evidence motivates.

**Restart condition.** If the least-squares system is ill-conditioned (condition number > 1e12 via thin QR), skip Anderson for this step, keep history. Prevents blowups in degenerate regimes.

**Interaction with damping.** Anderson applied AFTER damping: damping modifies the iterate `lf_{t+1}_plain = alpha * G(lf_t) + (1-alpha) * lf_t`; Anderson then mixes `lf_{t+1}_plain` with history. Damped-plain residual `r_t = lf_{t+1}_plain - lf_t = alpha * (G(lf_t) - lf_t)`.

**Warmup.** Anderson off for the first 5 iters (Sinkhorn preconditioner builds a useful starting point). After iter 5, engagement subject to the gates above.

**Auto-disable flag.** `LBW_IEPPA_ACCEL ∈ {on, off, unset}` env var. Default `on`. Tests using both paths required.

Expected impact on kk1204: 3–5× iter reduction on uncapacitated iters. With the `n_cap_active == 0` gate, benefit depends on fraction of iters where capacity settles — empirical check required.

## 6. Pass 3 — Bounds

Motivation (iter-1 PM): current cell-aggregate bounding can produce individual `w_i` outside user-stated `[min_weight, max_weight]` when base weights `d_i` are skewed within a cell. This is **documented but surprising** behaviour — an advanced use case for analysts requiring strict per-obs bounds. No authoritative statistical agency citation is claimed; prior rev 1 overstated provenance.

### 6.1 P3.1 `bounds_mode` parameter + water-filling expansion

**C API addition** (`src/leafblower.h`):
```c
typedef enum { RK_BOUNDS_CELL = 0, RK_BOUNDS_OBSERVATION = 1 } rk_bounds_mode_t;
```
Append `rk_bounds_mode_t bounds_mode;` to the end of `rk_params_t`. Default value `RK_BOUNDS_CELL` (preserves current behaviour exactly).

**`rk_params_init()` fix** (`src/c_api.cpp`, iter-1 Architect BLOCKER): the current init function assigns fields individually without zeroing. Stack-allocated callers would get garbage `bounds_mode`. Add `memset(p, 0, sizeof(*p));` as the first line of `rk_params_init()`. Also add `static_assert(sizeof(rk_params_t) == EXPECTED_BYTES, "rk_params_t size changed; check ABI consumers");` as a documentation anchor for ABI size change.

**R wrapper.** `harvest(..., bounds_mode = c("cell", "observation"))`. Default `"cell"`. Maps to C enum.

**Python wrapper.** `harvest(..., bounds_mode: Literal["cell", "observation"] = "cell")`. Same mapping.

**ABI impact.** `rk_params_t` size grows by `sizeof(int)` + padding. Backward compatible for callers using `rk_params_init()` after the memset fix. Raw-struct callers (stack-allocated without `rk_params_init`) relying on prior `sizeof(rk_params_t)` must recompile — the size grows. Documented in spec §11 release notes.

### 6.2 Intra-cell water-filling expansion

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

### 7.1 Speed regression (P1.1)

- kk1204 per-iter ratio ≤ 1.5× (down from 2.17×). Script `/tmp/wu2_kk1204.R` reused.
- Full suite: FAIL=0, PASS ≥ 181 (no regressions).
- Stepstone per-pass regression checkpoint (iter-1 PM): after P1.1 commit, rerun `/tmp/stepstone_2algo.R`; ieppa errRp must stay ≤ 2.21e-3 ± 1e-4 (numerical invariance across fusion).

### 7.2 Convergence (P2.1, P2.2)

- kk1204 at tol_abs=1e-3, max_iter=500: expect RK_OK. Currently NOCONV.
- Existing dense-equivalence test (test-ieppa-faithful.R:109-136): `|res_lin - res_log| < 1e-8` must still hold. Anderson on both paths via shared log-domain dispatch yields identical iterates.
- New: Anderson off vs on comparison. With `LBW_IEPPA_ACCEL=off`, iter count matches pre-P2 baseline ± 1. With `on`, at least 2× fewer iters on kk1204.
- Per-pass stepstone regression (iter-1 PM): after P2.1 and P2.2 commits, rerun stepstone bench; errRp ≤ 2.21e-3 ± 1e-4 preserved.

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

Four commits, one per WU (P1.1, P2.1, P2.2, P3.1), in pass-then-within-pass order. Each commit includes:
- Source edit.
- Test (or test extension).
- All existing tests pass.
- Stepstone per-pass regression check (errRp ≤ 2.21e-3 ± 1e-4).

## 10. Non-Goals / Out of Scope

- Sum-to-N drift (assessment §4) — handled at wrapper level; not re-opening.
- Configurable `kInfeasPersistence` / exposing streak state to user API.
- Parallel within-margin scatter (B-axis C5) — scope creep; latent ticket only if benchmarks motivate.
- SIMD-specific vectorization (C3) — latent ticket; benchmark-driven only.
- CRAN / PyPI distribution changes.
