# iEPPA Speed, Convergence, and Bounds Hardening Design

**Status:** Draft rev 5 (post iter-3 design-review-gate: 5/5 NEEDS_REVISION; gate cap exceeded by user authorization to close residual blockers)

**Rev 5 changes (iter-3 consensus blockers):**

*Field-name copy-paste bug (4-reviewer agreement).* §7.3 cell-mode warning test asserts `result.n_bounds_violated > 0` (was incorrectly `n_bounds_clamped`). Unit-mode test asserts `n_bounds_clamped == 0` on benign, `> 0` on degenerate.

*LAPACK near-singular detection (3-reviewer agreement).* `dgels` alone cannot detect near-singular systems. Added two defensive layers: (a) `m_active = min(m, total_cats - 1)` cap before call (prevents underdetermined branch); (b) post-solve `||γ||∞ > 1e4 → clear history, use plain iterate`. Removed contradictory "condition number > 1e12 via thin QR" sentence — we use dgels, not a separate QR step. See §5.2.

*P1.1 fallback state reset (2-reviewer agreement).* Named the existing WU-2 fallback semantics explicitly in §4.1: on `overflow_detected` mid-loop `break`, the fallback clears `X_cur[]`, `W[]`, `X[]`, `X_tilde[]`, `lf[]`, `f_lin[]`, `infeas_streak[]`, `result.n_xcur_writes_per_iter_linear` before restarting the outer iter loop in log-space.

*P2.1 alpha observability (CTO iter-3).* Added `double min_alpha_seen = 1.0` to `IEPPAResult`. Test asserts `result$min_alpha_seen < 1.0` on stress input + `>= 0.95` after stress subsides — struct-based, no log parsing.

*§8 iter floor (PM iter-3).* `250` iter gate softened to `kk1204 iter-to-RK_OK ≤ 400 with ACCEL_ANDERSON=on AND at least 2× fewer iters than ACCEL_ANDERSON=off baseline`. The "2× fewer" claim is robust to capacity-gate fraction; absolute iter count to 400 is generous (room for capacity to settle before Anderson engages).

*§8 unit-mode bounds gate (PM iter-3).* Added: `n_bounds_clamped == 0` on benign unit-mode input (dense, uniform-d); `n_bounds_clamped < 0.1% · n` on skewed-d stress input (not perfect adherence but bounded).

*R helper pattern (Designer iter-3).* `parse_bounds_mode` returns character, `.Call` site in `R/harvest.R` converts via `as.integer(factor(bounds_mode, levels=c("cell","unit"))) - 1L`. Matches the `as.character(method)` pattern used by `map_method`.

*EXPECTED_BYTES (Security iter-3).* Spec §6.1 now instructs implementer to compute actual `sizeof(rk_params_t)` at implementation time, record the value in a comment adjacent to the static_assert, and commit the value as a source-of-truth tripwire.

*Env var orthogonality claim (Security iter-3).* §7.4 rationale narrowed: `FORCE_PATH × FORCE_DAMPING` is orthogonal-by-construction; `ACCEL_ANDERSON × *` combinations require explicit coverage — T2/T3/T4/T5 provide it.

*Acknowledged and deferred (PM iter-3 B2).* P3.1 ABI growth breaks raw-struct callers. Documented in §10 release notes. No version-bump mechanism added — ABI break is the cost of carrying new API without a major version. If a raw-struct caller exists outside our tree, they recompile.

*Acknowledged and deferred (Architect iter-3 W2 caveat).* Water-fill "typically ≤ 3 iters" claim softened to "empirical observation on test inputs; 50-iter fallback ensures termination in all cases." No proof; acceptable per Security iter-3 confirmation that fallback correctness holds.

**Rev 4 changes (iter-2 consensus blockers):**

*LAPACK route (Architect+Security+CTO convergent blocker):* Replace `dgeqrf + dtrsv` with **`dgels`** (one-call least squares with QR). `dgels` handles rank-deficient cases internally via `dgelsd` fallback pattern, returns condition info via `INFO`. Explicit rank-revealing truncation via pivoted QR (`dgeqp3`) as fallback if `dgels` returns `INFO > 0`. See §5.2.

*`n_bounds_clamped` dual-use (Security iter-2 B6):* Rename cell-mode diagnostic field to `n_bounds_violated` (no action taken). Unit-mode fallback uses `n_bounds_clamped` (clamping action taken). Both fields exposed in `rk_calib_result_t`; wrappers warn on either non-zero. See §6.1, §6.2.

*RED-test mechanics (CTO iter-2 MF1+MF2):* Replace compile-time instrumentation and log-string parsing with result-struct counters. Add to `IEPPAResult`: `int n_xcur_writes_per_iter_linear` (P1.1 structural assertion); `int n_anderson_iters_engaged` (P2.2 engagement gate). Both always computed, always zero when not applicable. Tests assert on these integers directly. See §7.

*Env-var interaction matrix (Security iter-2 B8):* Add §7.4 enumerating 27-combo coverage. Minimum required: `linear × on × on`, `log × on × on`, `linear × on × off`, `log × off × on`, `unset × unset × unset` (auto). Others declared safe-by-construction with rationale.

*P1.1 overflow fallback integration (Architect iter-2 B1):* Fused block detects linear-overflow inline and triggers the existing WU-2 fallback path (reset state, switch to log-space, restart outer loop). Spec now shows this explicitly in §4.1 pseudocode.

*`damped_latched` removal (Architect iter-2 W1):* P2.1 deletes `damped_latched` — new formula is unlatched. `force_damping_on` becomes "force `alpha = 0.5` constant regardless of streak." See §5.1.

*§8 merge gate P2.2 iter floor (PM iter-2 B1):* Added concrete gate: `kk1204 iter-to-RK_OK ≤ 250` with `LBW_IEPPA_ACCEL_ANDERSON=on` (pairs with "at least 2×" from §7.2).

*Scope-table priorities (PM iter-2 B2):* P1.1=P1 (speed is user-top), P2.1/P2.2=P2, P3.1=P3.

*Env var deprecation note (PM iter-2 B3):* `LBW_IEPPA_ACCEL` → `LBW_IEPPA_ACCEL_ANDERSON`. Old name never shipped in a release; no deprecation path needed. Documented in §10.5 note.

*Designer minor revisions:* (1) parse_bounds_mode returns character, converted to int at `.Call` site matching map_method pattern. (2) Roxygen caveat added: "In degenerate single-obs cells, strict bounds cannot always be guaranteed; `result$n_bounds_clamped` reports residual clamp count." (3) §11 → §10.5 numbering reference corrected in preamble (done in rev 3, confirmed).

**Rev 3 changes (kept; superseded items marked):**

**Rev 3 changes (iter-1 design-review-gate):**

*Naming / API ergonomics (Designer):*
- `bounds_mode` values renamed: `"cell"` / `"observation"` → `"cell"` / `"unit"`. C enum: `RK_BOUNDS_CELL` / `RK_BOUNDS_UNIT`. "Unit" matches survey-calibration literature and disambiguates from "observation" (data row vs survey response).
- R wrapper uses `parse_bounds_mode()` helper matching `map_method()` pattern in `R/harvest.R:168`.
- Python wrapper types `bounds_mode: str` (not `Literal`) to match existing `method: str` typing in `_harvest.py:23`; docstring constrains values.
- User-facing copy mandated in §6.1 (roxygen and Python docstring text).
- Cell-mode silent violation: warning emitted when detected (§6.2).
- Env var `LBW_IEPPA_ACCEL` → `LBW_IEPPA_ACCEL_ANDERSON` (clearer). Env var table added in §11.

*Safety / numerics (Security):*
- §4.1 fused block: explicit `W[c] = 1.0` reset when X_tilde == 0 (prevent stale W propagation).
- §5.2 Anderson: `std::isfinite` check on all γ_j after LAPACK solve before writing back to `lf`; clear history on NaN detection.
- §5.2 Anderson: history buffer declared as `std::vector<double>` (RAII-safe on early exit).
- §5.2 Anderson: linear-path shadow-log conversion guarded by `f_lin[kj] > 0`.
- §6.1: `static_assert(RK_ALG_AUTO == 0, ...)` added to prevent memset silently changing default algorithm.
- §6.2: clamp-failure surfaced via new `result.n_bounds_clamped` counter (not silent); R/Python wrappers emit `warning()` / `warnings.warn()` when nonzero.

*Engineering (CTO):*
- §5.2 thin QR implementation specified: LAPACK `dgeqrf` + `dtrsv` via R's `R_ext/Lapack.h` (already a dependency).
- §8 merge gate tolerance 1e-12 → 1e-8 matching coarsest existing test (`test-ieppa-faithful.R:109-136`).
- §7.1 per-WU RED test added for P1.1 (structural assertion: only one post-sweep pass over M_cell in linear path, verifiable via `.wolf` trace or a test-only counter).
- §7.2 per-WU RED test for P2.1 (alpha recovers from 0.5 → 1.0 after stress subsides).
- §7.2 per-WU RED test for P2.2 engagement gate (capacity-toggle scenario).
- §7.3 cross-language ABI consistency test for `bounds_mode` enum mapping.
- §5.2 Anderson 3–5× iter-reduction claim softened to "empirical" with active-iter-fraction diagnostic logged at verbose=2.
- §9: P2.1 commit allowed to leave kk1204 gate RED; interim assertion: `errRp @ 500 iter ≤ 2.5e-3` so P2.1 does not regress beyond current NOCONV plateau.

**Rev 2 changes (kept):**
- Dropped P1.2 (errRp interval bump) per user directive.
- §4.1: corrected "3 passes → 1" framing to "2 post-sweep passes → 0"; sweep already maintains `X_cur` per WU-2.
- §5.2: Anderson + capacity BCD interaction — engagement gate `n_cap_active == 0`.
- §5.3 P2.3 merged into P2.2.
- §6.1: `memset` requirement in `rk_params_init()`.
- §6.1: reframed "statistical agency convention" claim.
- §7: per-pass stepstone regression checkpoint.

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
- **Pass 3 — Bounds.** A1 new `bounds_mode` parameter ∈ {`cell` (default, current), `observation`} with A2 intra-cell water-filling behind `observation` mode. Target: when `bounds_mode="unit"`, every returned weight satisfies `w_i ∈ [min_weight, max_weight]` strictly.

## 1. Problem Statement

**(C) Speed.** Linear-path per-iter work at K=20, M_cell=1M: three O(M_cell) passes per outer iter (sweep sum, X_tilde rebuild, capacity block), each touching `X_cur` or derivatives. Measured 2.17× raking's single sweep. Three passes = 3 memory round-trips per iter; fusing into one pass with local temporaries cuts traffic to 1 round-trip. errRp check runs every 10 iters at O(M_cell), adding ~10% overhead; a 50-iter interval drops this to 2%.

**(B) Convergence.** Classical Sinkhorn is sublinear (rate ≈ iter^−0.5). kk1204 at K=20, dense compression, 500-iter budget reaches errRp ~2e-3 but still NOCONV at tol_abs=1e-3. Published acceleration (Chen-Lin-Jordan 2020; Mishchenko-Defazio 2020; Altschuler-Niles-Weed-Rigollet 2017) via Anderson mixing on the log-factor iterate closes 5–10× of iters-to-convergence on this regime class without changing fixed-point semantics.

**(A) Bounds.** Current code enforces `sum_{i∈c} w_i ∈ [L·|c|, U·|c|]` (cell aggregate). On cells with skewed base weights `d_i`, applying the uniform cell multiplier can produce individual `w_i > max_weight` or `w_i < min_weight`. Documented in spec §8 of prior design; test `test-ieppa-nonuniform-d.R` verifies the cell-aggregate invariant. Users who need strict per-observation bounds (statistical agency convention) have no route today.

## 2. Scope

Three passes, four work units. Priority labels match user-stated C>B>A ordering (speed=top, convergence=middle, bounds=bottom). Within-pass order is the mandated execution order.

| Pass | WU | Addresses | Files | Priority |
|---|---|---|---|---|
| 1 | P1.1 | Eliminate post-sweep X_tilde + X_cur rebuild loops; fuse capacity inline | `src/ieppa.cpp` | **P1** |
| 2 | P2.1 | Adaptive damping schedule (replace hard 0.5 latch with `alpha = 1/(1+β·stress)`) | `src/ieppa.cpp` | **P2** |
| 2 | P2.2 | Anderson(m=5) acceleration on `lf` (both paths via log-domain dispatch) | `src/ieppa.cpp` | **P2** |
| 3 | P3.1 | `bounds_mode` parameter + water-filling expansion | `src/ieppa.cpp`, `src/c_api.cpp`, `src/leafblower.h`, `R/harvest.R`, `python/leafblower/_harvest.py` | **P3** |

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
// Linear path: fused post-sweep block (Security iter-1 B2: stale-W fix;
// Architect iter-2 B1: inline overflow detection + existing WU-2 fallback).
bool overflow_detected = false;
for (int c = 0; c < ct.M_cell; c++) {
    if (X_init[c] <= 0.0 || W[c] <= 0.0) {
        X[c] = 0.0; X_cur[c] = 0.0;
        W[c] = 1.0;
        continue;
    }
    double X_tilde_c = X_cur[c] / W[c];
    if (!std::isfinite(X_tilde_c) || X_tilde_c > kLinearOverflowTrip) {
        overflow_detected = true; break;                      // trigger fallback
    }
    if (X_tilde_c <= 0.0) {                                   // degenerate
        X[c] = 0.0; X_cur[c] = 0.0; W[c] = 1.0; continue;
    }
    double xc = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);
    X[c] = xc;
    W[c] = xc / X_tilde_c;
    X_cur[c] = xc;
    res.n_xcur_writes_per_iter_linear++;                      // CTO iter-2 MF1 counter
    if (xc != X_tilde_c) n_cap++;
}
if (overflow_detected) {
    // Full state reset (Architect+Security iter-3): clear ALL cell-level arrays
    // that may carry partial writes from the broken loop, PLUS the WU-2 fallback
    // invariants. Partial mid-loop writes to W[c], X[c], X_cur[c] are undone.
    std::fill(X_cur.begin(),  X_cur.end(),  0.0);
    std::fill(W.begin(),      W.end(),      1.0);
    std::fill(X.begin(),      X.end(),      0.0);
    std::fill(X_tilde.begin(),X_tilde.end(),0.0);
    std::fill(lf.begin(),     lf.end(),     0.0);
    std::fill(f_lin.begin(),  f_lin.end(),  1.0);
    std::fill(infeas_streak.begin(), infeas_streak.end(), 0);
    res.n_xcur_writes_per_iter_linear = 0;  // counter reset; post-fallback is log-space
    // Flip to log-space and restart outer iter loop.
    use_linear = false;
    linear_fallback_used = true;
    continue;
}
```

Log-path analog: step 2 computes `s = log_X_init[c] + Σ lf[m]` (K-loop, unchanged); fuse capacity inline using `X_tilde_c = exp(clip(s, kLogClip))` directly without a separate `X_tilde` store. No X_cur in log path.

Memory traffic reduction: linear path 3× (X_cur, W, X_tilde) round-trips → 1×. K-loop in log path unchanged so fusion gain smaller there.

**Framing correction (iter-1 Architect).** Prior rev 1 claimed "3 passes → 1"; the sweep itself is not fused (it remains the dominant pass). What P1.1 actually eliminates is the 2 post-sweep passes (X_tilde + X_cur rebuild) by computing capacity inline using the sweep's final `X_cur`. Framed as "2 post-sweep passes → 0."

Expected impact: 20–40% wall reduction on kk1204 linear path; marginal on stepstone log path (K-loop dominates).

## 5. Pass 2 — Convergence

### 5.1 P2.1 Adaptive damping schedule

**Current** (WU-3): when `infeas_streak[idx] ≥ kInfeasPersistence/2 = 2` for any bucket, `alpha` latches at 0.5 for the remainder of the solve. Hard 0.5 halves convergence rate unconditionally after first stress, even if stress subsides.

**New**: `alpha = 1.0 / (1.0 + beta * stress)` where `stress = min(max_k,j infeas_streak[cat_offset[k]+j], kInfeasPersistence)` over buckets with `target > 0`, `beta = 0.5`. At `stress=0`, `alpha=1.0` (fast path). At `stress=2`, `alpha=0.5`. At `stress=5` (cap), `alpha ≈ 0.286`. Monotone, smooth, unlatched (alpha recovers as streaks reset).

**Stress cap rationale (P2.1 landing empirical finding, commit a3a5f77):** without the `min(..., kInfeasPersistence)` cap, stepstone drove stress to 490+ (one slow-settling bucket with 490-iter transient streak). Uncapped formula gave α ≈ 0.004, stalling the solver and regressing errRp from 2.21e-3 → 2.28e-3. With cap, stepstone improves to 2.14e-3. The cap reflects the design invariant that `kInfeasPersistence` is the boundary between "transient recoverable" and "de-facto stuck"; damping stronger than α=0.29 doesn't un-stick a stuck bucket (it only slows the rest of the solver).

An equivalent principled formulation would exclude `infeas_streak[idx] >= kInfeasPersistence` buckets from the stress max (declaring them structural for damping purposes while keeping them outside `structural_infeas_pairs` for the INFEAS status gate). Both formulations produce α ≥ 0.286 in practice; the `min()` cap is simpler.

Fast-path preserved: `if (stress == 0) { ... naive update ... }` branch avoids per-sweep recomputation of alpha.

**`damped_latched` removed** (Architect iter-2 W1). The pre-WU-3 boolean latch is deleted; alpha is recomputed each sweep from current streak state. Matching the new semantics, `LBW_IEPPA_FORCE_DAMPING=on` forces `alpha = 0.5` constant (regardless of streak); `=off` forces `alpha = 1.0` regardless of streak; unset yields the adaptive schedule above.

Preserves Peyré-Cuturi 2019 Proposition 4.8 convergence guarantee (alpha ∈ (0, 1]). Rate slows by 1/alpha at any iter.

### 5.2 P2.2 Anderson(m=5) acceleration on `lf` (both paths via log-domain dispatch)

**Method.** After each full margin sweep (one `k=0..K` cycle) followed by the capacity block, treat `lf` as the iterate in a fixed-point iteration `lf_{t+1} = G(lf_t)`. Anderson(m) mixing: maintain history of last `m+1` iterates `{lf_t, lf_{t-1}, ..., lf_{t-m}}` and residuals `{r_t, r_{t-1}, ..., r_{t-m}}` where `r_t = G(lf_t) - lf_t`. Next iterate: `lf_{t+1}_AA = lf_t + sum_j γ_j (G(lf_{t-j}) - lf_{t-j})` where γ solves a least-squares problem over residual differences.

**Both paths via log-domain dispatch.** Anderson always operates on `lf`. In the linear path, maintain a shadow `lf = log(f_lin)` state: after each sweep compute `lf[kj] = (f_lin[kj] > 0.0) ? log(f_lin[kj]) : -kLogClip` (Security iter-1: `f_lin == 0` guard prevents -inf) for updated buckets (O(total_cats)); apply Anderson; then materialize `f_lin[kj] = exp(lf_AA[kj])` (O(total_cats)). The `log+exp` per (k,j) per iter is O(total_cats) = O(K·avg_cats), dominated by the O(K·M_cell) sweep at M_cell ≫ avg_cats. Empirical check required (merge-gate budget section §7.1).

**Parameter.** `m = 5`. History buffer: 2 × (m+1) × total_cats doubles. At total_cats=2000, ≈ 200KB per solver call — negligible.

**Implementation details (rev-5 fix).** Use LAPACK `dgels` (one-call least-squares solver via QR) from `R_ext/Lapack.h`. Signature:
```c
void F77_CALL(dgels)(const char* trans, int* m, int* n, int* nrhs,
                      double* A, int* lda, double* B, int* ldb,
                      double* work, int* lwork, int* info);
```
Call with `trans="N"`, `A = F` (residual differences, column-major, `total_cats × m_active`), `B = -r_t` (current residual, length `total_cats`). On return: first `m_active` entries of `B` contain γ; `INFO == 0` on success.

**Shape guard (Architect iter-3).** Before the call: `m_active = std::min(m_active, total_cats - 1)`. Prevents `dgels` from entering the underdetermined branch (`m < n`) which silently returns minimum-norm solutions instead of LS. At `total_cats ≤ 5` (e.g., K=1 binary margin), effective history depth is capped. Degenerate `total_cats == 0` or `1` → skip Anderson entirely this step.

**Rank-deficiency / near-singular handling (Architect+Security+CTO iter-3 consensus).** Three defensive layers, applied in order after `dgels` returns:
1. **Exact rank-deficiency:** `INFO > 0` → clear history (reset `m_active = 0`), use plain iterate, increment `result.n_anderson_nan_fallbacks`.
2. **NaN blowup:** any `!std::isfinite(γ_j)` → same action as (1).
3. **Near-singular silent-wrong-answer:** `std::fabs(γ_max) = max_j |γ_j| > kGammaNormMax` (constant: `1e4`) → same action as (1). This closes the near-singular silent-corruption gap without a second LAPACK call (O(m) scan).

`dgels` does not return condition number. The three checks together bound both exact-singular and near-singular failure modes. `dgelsd` (SVD + rank truncation) is NOT used — adds implementation complexity and per-call cost; the γ-norm check is sufficient evidence.

**History matrices.** `F` (residual differences), `X_hist` (iterate differences), both `total_cats × m` column-major, stored as `std::vector<double>` (RAII on all early-exit paths). `m_active = min(m, t - warmup_iters, total_cats - 1)` columns used at iter `t`.

**Removed rev-4 sentence:** the "restart if condition number > 1e12 via thin QR" clause is deleted — contradicted the move to `dgels`.

**Capacity-block interaction (iter-1 Architect BLOCKER; revised post-P2.2 landing for joint iterate).**

Original (single-iterate Anderson, post-P2.2 landing `adba0c7`): Anderson fires only on `n_cap_active == 0` iters. Empirically FAILED on kk1204 K=20 max_weight=3: `n_cap_active` saturates at ~19,485 from iter 3 and never clears. Anderson never engages; kk1204 remains NOCONV.

**Revised (P2.2c APVA — Asymmetric Partial-Variable Anderson):** attempt 1 (P2.2b — joint iterate with W write-back, spec rev post-adba0c7) empirically diverged: Anderson-extrapolated `W_aa` conflicts with the deterministic capacity-block update `W = X/X_tilde`, producing oscillatory divergence. Scale mismatch between `lf` (bounded O(1)) and `log W` (large negative) corrupted the LS γ solution. Spec revised 2026-04-24 after web research (Fu-Liang-Toh 2020 arxiv:1902.04273; Pollock-Rebholz 2021; Scieur-D'Aspremont-Bach 2020 arxiv:2010.03816; Tang et al 2024 arxiv:2401.07127).

**APVA semantics:**

1. **Joint residual** `R_t = [r_lf_t; σ · r_W_t]` of length `total_cats + M_cell`. `r_W_t = log_W_t - log_W_{t-1}` treated as CONDITIONING SIGNAL only. `σ = 1 / max(1.0, std::max(|log_W|))` rescales W-component to match lf-scale. Computed fresh each iter; no hand-tuned constant.
2. **LS solve** via `dgels` on the expanded residual: `F_hist` columns of length `total_cats + M_cell`, solving `min_γ || F_hist · γ - R_t ||` with the 3-layer guard (INFO + isfinite(γ) + |γ|∞ < 1e4).
3. **lf-only update.** `lf_new[kj] = lf[kj] + r_lf_t[kj] - Σ_j γ_j (X_hist_lf[kj, j] + F_hist_lf[kj, j])` where `X_hist_lf` and `F_hist_lf` are the lf-slice (first total_cats rows) of the history buffers.
4. **No W write-back.** W is NOT modified by Anderson. The capacity block continues to set `W = clamp(X_tilde, L, U) / X_tilde` from the Anderson-updated `lf`. This breaks the oscillation loop — W stays consistent with the capacity KKT conditions; Anderson only influences W indirectly via the new lf.

**Engagement gate:** drop `n_cap_active == 0`. Warmup (iter > kAndersonWarmup) + NaN/γ-norm guards remain. Add **safeguarding** (Fu et al 2020): if `||r_lf_aa|| > ||r_lf_plain||` (accelerated step didn't reduce lf residual), reject the Anderson update, clear history, use the plain damped iterate.

**History buffer sizes:**
- `F_hist`, `X_hist`: `(total_cats + M_cell) × m` column-major std::vector<double>. At M_cell=1M, m=5: 2 × 5 × 1M × 8B = 80 MB per solver call.
- `log_W`, `log_W_prev`: length M_cell each. 16 MB at M_cell=1M.
- Allocations as std::vector<double>; RAII on all early-exit paths.

**Non-smoothness mitigation layers** (apply in order):
1. P2.1 `compute_alpha` damping — α ∈ [0.286, 1.0] applied before Anderson mix.
2. APVA — no W write-back avoids clamp-induced discontinuity in the Anderson-modified state.
3. Safeguarding — reject Anderson step that doesn't improve lf residual.
4. NaN / γ-norm guard — hard stop on LS blowup.

**Halpern fallback (follow-up ticket, not in this scope):** if APVA still fails on adversarial inputs, a Halpern-style mixing `lf_{t+1} = lf_0/(t+2) + ((t+1)/(t+2)) · G(lf_t)` provides O(1/k) convergence guarantee immune to oscillation, slower than Anderson's super-linear potential but architecturally stable.

## §5.3 P2.2e — Tang 2024 primal-dual with Newton on constraint duals [SHELVED 2026-04-24]

**Status:** Research note only; NOT an active implementation target. Plan-review-gate iter-1 on this section FAILED unanimously (all 3 reviewers — Feasibility/Completeness/Scope). Root issues:

- **Theorem 1 K-count mismatch.** Paper's K counts INEQUALITY constraints; iEPPA box constraints produce K_ineq = 2·M_cell. At kk1204 this yields `n^(2/(2M_cell+1)) ≈ 1.0`, making the convergence bound trivially large and useless.
- **Δ runtime-unknowable.** Theorem 1 requires `η·Δ/(K+1)` grow large, where Δ is the LP-polytope constant. Paper (p.3): "depending only on the LP" — not computable from inputs at runtime without solving the LP.
- **η annealing required** per paper §3.1 for practical convergence; spec deferred it to follow-up but that schedule is the primary acceleration mechanism.
- **Capacity-block replacement** substantially exceeds the authorized "Anderson on lf" scope.
- **Per-cell 2D Newton decomposition** ≠ Algorithm 1's joint Newton; convergence equivalence unproven.

If Halpern (§5.4, active P2.2d) fails to close the kk1204 gate, this section will be reopened with: runtime Δ estimator, η-annealing schedule, joint-Newton shape analysis, and explicit user re-authorization of the capacity-block replacement.

Historical design content (DO NOT IMPLEMENT as-is):

**Correct reference:** Tang, Rahmanian, Shavlovsky, Thekumparampil, Xiao, Ying (2024), *"A Sinkhorn-type Algorithm for Constrained Optimal Transport"*, arxiv:2403.05054. (Prior spec rev cited wrong arxiv ID; corrected 2026-04-24.)

**Motivation.** APVA (P2.2c, commit 7a837fe) failed kk1204: capacity-clamp non-smoothness rejected 487/500 Anderson steps. Tang 2024 handles box-constraint non-smoothness at the formulation level — lift constraints into dual variables, solve via Newton + Sinkhorn alternation.

**Tang 2024 variational formulation (eq 5–6, algorithm 1).** Constrained entropic OT:
```
min_{P, s_1..s_K}  C·P + (1/η) H(P, s_1..s_K)
s.t.  P·1 = r,  P^T·1 = c,  D_k·P = s_k (k=1..K inequality),  P ≥ 0, s_k ≥ 0
```
H is joint entropy; s_k are slack variables; D_k matrices encode inequality constraints.

Dual problem after eliminating P, s: `max_{x,y,a} f(x, y, a)` where
```
f(x,y,a) = -1/η · Σ_ij exp(η·(-C_ij + Σ_m a_m D_m_ij + x_i + y_j) − 1)
           + Σ_i x_i r_i + Σ_j y_j c_j − 1/η · Σ_k exp(-η·a_k − 1)
```
`x`, `y` = standard row/column Sinkhorn duals (generalize to K margins via existing sweep). `a` = constraint duals, one per inequality constraint.

**Algorithm 1 alternating updates:**
1. **x update** = argmax_x̃ f(x̃, y, a). Closed-form Sinkhorn row scaling.
2. **y update** = argmax_ỹ f(x, ỹ, a). Closed-form column scaling.
3. **(a, t) update** = argmax_{ã, t̃} f(x + t̃·1, y, ã). **Newton's method** on a (with t auxiliary for numerical stability). Gradient `∂_ak f = exp(-η·a_k − 1) − P·D_k` ; Hessian computed analytically. Step direction `(Δa, Δt) = -(∇²f)^-1 ∇f` with backtracking line search.

Intermediate transport plan: `P = exp(η·(-C + Σ_m a_m D_m + x·1^T + 1·y^T) − 1)`.

**Mapping to iEPPA box constraint `L_cell[c] ≤ X[c] ≤ U_cell[c]`:**
- Each cell `c` generates **two inequality constraints**: `-X[c] ≤ -L_cell[c]` and `X[c] ≤ U_cell[c]`. Total `2·M_cell` inequality constraints, hence `2·M_cell` dual variables `a_{c,lo}`, `a_{c,hi}`.
- D matrices are sparse — each constraint touches one cell. Vectorized updates.
- x, y in Tang = lf dual in iEPPA (per-margin-per-category Sinkhorn factor), generalized K times. Current sweep already implements the K-margin equivalent of x, y cycling.

**Algorithm integration with current ieppa_solve:**

```
// Outer iter: (lf-sweep) + (a-Newton) + (errRp check)
For outer_iter = 1..kMaxOuter:
  // Sweep: k=0..K-1 cycles of lf-update (x, y generalized) — UNCHANGED from current
  sweep_lf(lf, a_lo, a_hi, ...);  // but log_S_kj accumulator now includes a-terms

  // Newton on constraint duals (replaces capacity block):
  For c = 0..M_cell-1:
    Compute X[c] from (lf, a_lo, a_hi):
      X[c] = X_init[c] · exp(η·(Σ_m lf[m][g_m(c)] - a_lo[c] + a_hi[c]))
    // Gradient: ∂_a f at c: for a_lo[c]: exp(-η·a_lo[c] - 1) - (-X[c] - (-L_cell[c])) · [indicator];
    //                       for a_hi[c]: exp(-η·a_hi[c] - 1) - (X[c] - U_cell[c]) · [indicator].
    // Hessian scalar ≥ 0; Newton step with backtracking.
    newton_update(a_lo[c], a_hi[c], X[c], L_cell[c], U_cell[c]);

  // errRp check (unchanged)
```

**No capacity block clamp.** The clamp is replaced entirely by the `a` dual variables — KKT complementarity (`a_hi[c] · (X[c] - U_cell[c]) = 0` and `a_lo[c] · (L_cell[c] - X[c]) = 0`) is the optimality condition, not a post-update projection.

**Newton inner iteration:** per cell, scalar 2-D Newton on (a_lo[c], a_hi[c]). Convergence to machine precision in 3-5 Newton iters per outer iter. Backtracking line search bounds step to maintain `a_lo ≥ 0`, `a_hi ≥ 0`. Cost: O(M_cell) per outer iter, negligible vs sweep O(K·M_cell).

**Regularization η.** Tang 2024 proposes η-annealing: start small (strong regularization, smooth problem), increase toward ∞ (weak regularization, closer to original LP solution). For initial implementation: fixed η = 1/tol_abs per user setting (so `exp(-η·Δ)` matches `tol_abs` error bound from Theorem 1). Follow-up ticket: dynamic annealing if fixed-η insufficient.

**Convergence guarantee (Theorem 1 + Theorem 3 of Tang):** sublinear first-order in dual space. Error bound: `||P_η* − P*||_1 ≤ 8n^(2/(K+1))(K+1)·exp(-η·Δ/(K+1))`. At η = 1e3 (tol=1e-3), `K+1=21`: exp bound ~ exp(-47·Δ). For feasible problems Δ > 0 bounded below; error drops as η grows.

**Memory:** `a_lo`, `a_hi` vectors length M_cell. At M_cell=1M: 16 MB. Cheap.

**API:** `LBW_IEPPA_SOLVER ∈ {tang2024, direct, off}`. Default `direct` (current behavior, capacity block) — byte-identical preservation. `tang2024` activates primal-dual Newton path. Targets kk1204-class non-smooth regimes.

**Merge gate:**
- kk1204 SOLVER=tang2024: `status == RK_OK` AND `total outer iter × sweep-work ≤ 400 iter-equivalents`.
- Stepstone SOLVER=tang2024: errRp ≤ 2.15e-3 (same as direct on sparse regimes — Newton on capacity duals is cheap).
- Stepstone SOLVER=direct: byte-identical to current (regression safety).
- Full suite FAIL=0 PASS ≥ 198.

**Deviation from my prior spec:** earlier spec described "augmented Lagrangian with quadratic penalty + ρ annealing". That's a DIFFERENT technique (Birgin-Martínez textbook NLP). Tang 2024 is variational-primal-dual with Newton — simpler inner loop, no ρ annealing, Theorem 1 convergence guarantee for specific bounds.

**Fallback: Halpern mixing (§5.4 follow-up)** if Tang approach fails to close kk1204 gate.

## §5.4 P2.2d — Halpern mixing [ACTIVE IMPLEMENTATION TARGET, post-Tang pivot 2026-04-24]

**Status:** Promoted from fallback to primary P2.2 engineering response after Tang 2024 (§5.3) review-gate failure. Halpern mixing provides O(1/k) convergence guarantee on the existing iEPPA fixed-point operator without algorithmic replacement.

**Empirical status post-P2.2c APVA (commit 7a837fe):** kk1204 ACCEL=on shows engaged=6, nan_fallbacks=487 of 500 iters. Safeguard rejects nearly every Anderson step because capacity-clamp non-smoothness generates extrapolated lf residuals larger than plain damped residuals. Stepstone + 198 tests green — APVA doesn't regress other regimes, just fails on kk1204's structural non-smoothness.

**Halpern mixing (Halpern 1967; rediscovered for fixed-point iteration in Lieder 2021 "Halpern's Iteration in Convex Programming"):**

```
For m = iter - kHalpernAnchor (m = 0 at anchor iter, m = 1 at first post-anchor, ...):
  lf_{k+1}[kj] = (1 / (m + 1)) · lf_anchor[kj] + (m / (m + 1)) · G(lf_k)[kj]
```

- `lf_anchor` = snapshot of `lf` at iter `kHalpernAnchor` (= 5).
- `G(lf_t)` = one full outer iter (sweep + capacity block) applied to `lf_t`.
- At m=0: weight 1.0 to anchor (capture-only; no mix). At m=1: weight 1/2 / 1/2. At m→∞: weight → 0 / 1.0. Lieder 2021 canonical form.

**Convergence guarantee:** for non-expansive G (Sinkhorn + capacity clamp is non-expansive in ℓ∞ norm on the feasible region), Halpern converges with `||r_t|| = O(1/t)`. Not super-linear, but immune to oscillation — no safeguard required.

**Replaces APVA entirely** (P2.2c code removed in P2.2d commit):
- No dgels, no LS solve, no history buffers.
- No `n_anderson_iters_engaged` / `n_anderson_nan_fallbacks` counters (replaced with `n_halpern_iters`).
- Env var renamed: `LBW_IEPPA_ACCEL_ANDERSON` → `LBW_IEPPA_ACCEL` with values `{halpern, off}`. **Default `off`** (conservative opt-in during P2.2d rollout; flip to `halpern` default in a later commit after empirical validation across all regimes).
- Memory: ~80 MB history buffer gone; `lf_anchor` is total_cats doubles (≪ 1 MB at kk1204).

**Compatibility with P2.1 damping:** Halpern applied AFTER damped iterate. Let `lf_damped = α·G_raw(lf_t) + (1-α)·lf_t` (from compute_alpha). Then `lf_halpern = weighted blend of lf_anchor and lf_damped`. Composition preserves contraction since both steps are convex combinations.

**Merge gate:**
- Stepstone errRp ≤ 2.15e-3 preserved (Halpern should not regress sparse regimes)
- kk1204 ACCEL=halpern: `iters ≤ 450` AND `status == RK_OK`. Relaxed from APVA's 400 because Halpern's O(1/t) rate is slower than super-linear Anderson would have been.
- Full suite FAIL=0 PASS ≥ 198

**Fallback to APVA if Halpern fails kk1204:** Halpern's O(1/t) theoretical rate may be too slow to reach tol_abs=1e-3 in 500 iters on n=1M. If so, file ticket for E (Lagrangian relaxation) as follow-up.

**Restart condition.** If the least-squares system is ill-conditioned (condition number > 1e12 via thin QR), skip Anderson for this step, keep history. Prevents blowups in degenerate regimes.

**Interaction with damping.** Anderson applied AFTER damping: damping modifies the iterate `lf_{t+1}_plain = alpha * G(lf_t) + (1-alpha) * lf_t`; Anderson then mixes `lf_{t+1}_plain` with history. Damped-plain residual `r_t = lf_{t+1}_plain - lf_t = alpha * (G(lf_t) - lf_t)`.

**Warmup.** Anderson off for the first 5 iters (Sinkhorn preconditioner builds a useful starting point). After iter 5, engagement subject to the gates above.

**Auto-disable flag.** `LBW_IEPPA_ACCEL_ANDERSON ∈ {on, off, unset}` env var. Default `on`. Tests using both paths required.

Expected impact on kk1204: empirical — claim of 3–5× iter reduction conditional on fraction of iters where capacity is inactive. Diagnostic logged at verbose=2: `"Anderson: <M> of <N> iters engaged"`. Merge gate is kk1204 RK_OK within 500 iter, not a specific speedup multiplier.

## 6. Pass 3 — Bounds

Motivation (iter-1 PM): current cell-aggregate bounding can produce individual `w_i` outside user-stated `[min_weight, max_weight]` when base weights `d_i` are skewed within a cell. This is **documented but surprising** behaviour — an advanced use case for analysts requiring strict per-obs bounds. No authoritative statistical agency citation is claimed; prior rev 1 overstated provenance.

### 6.1 P3.1 `bounds_mode` parameter + water-filling expansion

**C API addition** (`src/leafblower.h`):
```c
typedef enum { RK_BOUNDS_CELL = 0, RK_BOUNDS_UNIT = 1 } rk_bounds_mode_t;
```
Append `rk_bounds_mode_t bounds_mode;` to the end of `rk_params_t`. Default value `RK_BOUNDS_CELL` (preserves current behaviour exactly). "Unit" chosen (Designer iter-1 B1) over "observation" — matches survey-calibration literature; "cell" denotes aggregate-level bounding, "unit" denotes per-respondent bounding.

**`rk_params_init()` fix** (`src/c_api.cpp`, Architect iter-1): current init assigns fields individually without zeroing; stack-allocated callers would get garbage `bounds_mode`. Add `memset(p, 0, sizeof(*p));` as first line. Add two guards in the same header:
```cpp
static_assert(RK_ALG_AUTO == 0, "memset(0) default must equal RK_ALG_AUTO");  // Security iter-1
// EXPECTED_BYTES: compute sizeof(rk_params_t) on the target platform at implementation time
// (likely 56-64 bytes depending on padding), hard-code the observed value, and record it in
// a comment adjacent to this assert as the ABI tripwire. Any future field addition that
// changes the size fails this assert, forcing an intentional re-bless of consumers.
static_assert(sizeof(rk_params_t) == EXPECTED_BYTES, "rk_params_t size changed; check ABI consumers");
```

**R wrapper** (Designer iter-3 fix — character-return matching map_method pattern). Use a helper `parse_bounds_mode(x)` in `R/harvest.R` that returns a validated character string (not an integer):
```r
parse_bounds_mode <- function(x = c("cell", "unit")) {
  match.arg(x)  # returns "cell" or "unit"; errors with a user-friendly message otherwise
}
```
Conversion to the C enum integer happens at the `.Call` bridge site (same pattern as `as.character(method)` used by `map_method`):
```r
bounds_mode_char <- parse_bounds_mode(bounds_mode)
bounds_mode_int  <- match(bounds_mode_char, c("cell", "unit")) - 1L  # 0 or 1
# ... passed through raw() / .Call / as.integer(bounds_mode_int)
```
`harvest(..., bounds_mode = "cell")` is the default.

**Python wrapper** (Designer iter-1 B3). Match existing `method: str` typing in `_harvest.py:23`:
```python
def harvest(..., bounds_mode: str = "cell"):
    if bounds_mode not in ("cell", "unit"):
        raise ValueError(f"bounds_mode must be 'cell' or 'unit', got {bounds_mode!r}")
    params.bounds_mode = {"cell": 0, "unit": 1}[bounds_mode]
```
Docstring constrains values. Future typing overhaul applies `Literal` to both `method` and `bounds_mode` together (out of scope for this spec).

**User-facing copy (Designer iter-1 R3).** R roxygen `@param bounds_mode`:
```
Bounds enforcement mode. "cell" (default): per-cell aggregate bounds —
sum of weights within each cross-classified demographic cell is bounded,
but individual weights may fall outside [min_weight, max_weight] if
initial design weights are skewed within a cell. "unit": per-observation
bounds — every returned weight strictly satisfies [min_weight, max_weight]
via intra-cell water-filling redistribution.
```
Python docstring mirrors same text.

**ABI impact.** `rk_params_t` size grows by `sizeof(int)` + padding. Backward compatible for callers using `rk_params_init()` after the memset fix. Raw-struct callers (stack-allocated without `rk_params_init`) relying on prior `sizeof(rk_params_t)` must recompile — the size grows. Documented in §11 release notes.

### 6.2 Intra-cell water-filling expansion

**Trigger.** Only when `bounds_mode == RK_BOUNDS_UNIT`.

**Cell-mode diagnostic (Designer iter-1 R1; Security iter-2 B6 — split field).** When `bounds_mode == RK_BOUNDS_CELL` (default), scan the final expanded weight vector for bound violations. If any `w_i ∉ [min_weight, max_weight]`, set `result.n_bounds_violated = count_of_violations` (no clamping performed in cell mode — diagnostic only). R/Python wrappers emit a `warning()` / `warnings.warn()` at cell-mode exit when `n_bounds_violated > 0`:
> "cell-mode bounds: N of M weights fell outside [min_weight, max_weight] due to skewed base weights within cells. Consider bounds_mode='unit' for strict per-observation bounds."

**Clamp-failure signal in unit mode (Security iter-1 B6).** Distinct field `int n_bounds_clamped;` (not shared with cell-mode diagnostic). Under `bounds_mode == RK_BOUNDS_UNIT`, if the per-cell water-fill inner loop exhausts its 50-iter budget with residual violations, clamp remaining out-of-bound weights and increment `n_bounds_clamped`. R/Python wrappers emit a distinct warning when nonzero:
> "unit-mode bounds: N weights clamped to nearest bound after water-fill exhausted; degenerate cells could not be fully redistributed."

Both fields added to `IEPPAResult` (`src/ieppa.hpp`) and propagated to `rk_calib_result_t` via `c_api.cpp`. Wrappers fire separate warnings keyed on each field so users can tell diagnosis (cell-mode) from action (unit-mode) apart.

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

Per-cell inner loop typically converges in ≤ 3 iters (few distinct bound-hit patterns). Pathological alternating `max`/`min` patterns exist where the loop does not strictly reduce violators per iter; the 50-iter fallback ensures termination. Early-exit (Security iter-1 #5): `if (!bounds_violated) break` at top of each inner iter. On exhaustion, clamp residuals to nearest bound and increment `result.n_bounds_clamped` (§6.2 clamp-failure signal). Worst case 50 iters × |cell| per cell; n=10^6 → 5×10^7 ops, negligible against solver cost.

**Invariant.** After water-filling: `sum_{i∈c} w_i = X[c]` preserved (redistribution conserves mass), and `w_i ∈ [min_weight, max_weight]` for all i (up to clamp failure in pathological cells, warned).

**Implementation.** Group observations by cell ahead of time (`cells_of_obs[c]` list). Or iterate over `cells_of_obs[c] = [i : ct.cell_of[i] == c]` built in one O(n) pass. This index is built once pre-loop; memory O(n) ints.

**Corner case.** If cell has only one observation and its `w_i` is out of bounds, water-filling cannot redistribute. Fall back: clamp to nearest bound, emit warning, accept the cell-aggregate violation.

### 6.3 Tests

- `test-ieppa-bounds-mode.R`: when `bounds_mode="unit"` on highly skewed `d_i` input, assert `max(w) ≤ max_weight + 1e-9` and `min(w) ≥ min_weight - 1e-9` and `|sum(w) - n| < 1e-6`.
- Existing `test-ieppa-nonuniform-d.R`: unchanged behaviour when `bounds_mode="cell"` (default).

## 7. Testing

### 7.1 Speed regression (P1.1)

- **Per-WU RED test** (CTO iter-2 MF1 fix). Add `int n_xcur_writes_per_iter_linear;` to `IEPPAResult`. The fused block (§4.1) increments it per `X_cur[c]` write. Test asserts `result$n_xcur_writes_per_iter_linear / result$iterations == ct.M_cell` on a linear-path input (not `2 * M_cell` or `3 * M_cell` as pre-P1.1). Field always computed — no compile-time instrumentation. Pre-P1.1 the field exists but is zero (no writes tracked); the test fails with "expected M_cell, got 0" — genuine RED. Post-P1.1 the test passes.
- kk1204 per-iter ratio ≤ 1.5× (down from 2.17×). Script `/tmp/wu2_kk1204.R` reused.
- Full suite: FAIL=0, PASS ≥ 181 (no regressions).
- Stepstone per-pass regression checkpoint (PM iter-1): after P1.1 commit, rerun `/tmp/stepstone_2algo.R`; ieppa errRp must stay ≤ 2.21e-3 ± 1e-4.

### 7.2 Convergence (P2.1, P2.2)

**P2.1 per-WU RED tests** (CTO iter-3 fix — struct-based, no log parsing):
1. Alpha observability. Add `double min_alpha_seen = 1.0;` to `IEPPAResult`. Each sweep: `min_alpha_seen = std::min(min_alpha_seen, alpha)`. Tests assert on this struct field.
   - Stress input (adversarial): assert `result$min_alpha_seen < 1.0` (damping engaged).
   - Benign input (no stress): assert `result$min_alpha_seen == 1.0` (fast path).
   - Stress-then-recovery: run with stress present early; assert `result$min_alpha_seen < 0.9` AND final alpha at exit ≥ 0.95 (recovered). Add `double final_alpha` field for the latter.
2. Interim bound: after P2.1 alone (without P2.2), kk1204 errRp @ 500 iter must be ≤ 2.5e-3 (not worse than current NOCONV plateau at ~2.0e-3; prevents adaptive schedule regressing convergence).

**P2.2 per-WU RED tests** (CTO iter-2 MF2 fix — result-struct counters, not log parsing):
1. Capacity-toggle engagement. Add `int n_anderson_iters_engaged;` to `IEPPAResult` (always computed; zero when Anderson off). Construct an input that alternates `n_cap_active > 0` ↔ `n_cap_active == 0` across iters. Run twice: once with `LBW_IEPPA_ACCEL_ANDERSON=off`, once with `=on`. Assert `result$n_anderson_iters_engaged == 0` under `off`. Under `on`, assert `n_anderson_iters_engaged` equals exactly the count of iters where `n_cap_active == 0` (minus the 5-iter warmup).
2. NaN guard. Add `int n_anderson_nan_fallbacks;` to `IEPPAResult`. Synthetic input that produces rank-deficient residuals; assert `result$n_anderson_nan_fallbacks > 0` AND solver weights are all finite (no NaN leakage).

**Aggregate P2 gates.**
- kk1204 at tol_abs=1e-3, max_iter=500: expect RK_OK after P2.2. Currently NOCONV.
- Existing dense-equivalence test (test-ieppa-faithful.R:109-136): `|res_lin - res_log| < 1e-8` must still hold. Anderson on both paths via shared log-domain dispatch yields identical iterates.
- Anderson off vs on comparison. With `LBW_IEPPA_ACCEL_ANDERSON=off`, iter count matches pre-P2 baseline ± 1. With `on`, at least 2× fewer iters on kk1204.
- Per-pass stepstone regression (PM iter-1): after P2.1 and P2.2 commits, rerun stepstone bench; errRp ≤ 2.21e-3 ± 1e-4 preserved.

### 7.3 Bounds (P3.1)

- `test-ieppa-bounds-mode.R`: skewed-d input returns strict bounds under `bounds_mode="unit"`. Assert `max(w) ≤ max_weight + 1e-9`, `min(w) ≥ min_weight - 1e-9`, `|sum(w) - n| < 1e-6`, `result.n_bounds_clamped == 0` on benign input.
- Cell-mode warning: skewed-d input with cell-mode default assert `result.n_bounds_violated > 0` (NOT `n_bounds_clamped` — cell-mode is diagnostic-only; 4-reviewer iter-3 consensus) and R `expect_warning(harvest(...), "cell-mode bounds")`.
- Clamp-failure path: pathological single-obs cell with out-of-bound weight → `result.n_bounds_clamped > 0` + warning surfaces even in unit mode.
- Regression: default `bounds_mode="cell"` returns current behaviour (test-ieppa-nonuniform-d.R untouched).
- C API: zero-initialized `rk_params_t` yields `RK_BOUNDS_CELL`.
- Python: `bounds_mode` not passed → default cell behaviour.
- Cross-language ABI consistency (CTO iter-1): test file directly invokes `.Call` with `params$bounds_mode = 1L` (raw integer) and asserts unit-mode semantics — validates C enum integer mapping matches the R/Python string→int logic.

### 7.4 Env var interaction matrix (Security iter-2 B8)

Three vars × 3 values each = 27 combinations. Minimum required tests covering all non-trivial interactions:

| Test | `FORCE_PATH` | `FORCE_DAMPING` | `ACCEL_ANDERSON` | Rationale |
|---|---|---|---|---|
| T1 | unset | unset | unset | auto — default production path |
| T2 | linear | on | on | stress all three: damped linear sweep + shadow-lf Anderson |
| T3 | log | on | on | damped log sweep + direct-lf Anderson |
| T4 | linear | off | on | Anderson w/o damping on linear (most common acceleration case) |
| T5 | log | off | on | Anderson w/o damping on log (stepstone-class input) |
| T6 | linear | off | off | bare linear path (WU-2 baseline behaviour) |
| T7 | log | off | off | bare log path (pre-WU-3 baseline behaviour) |

**Orthogonality narrowed** (Designer+Security iter-3). Only `FORCE_PATH × FORCE_DAMPING` is truly orthogonal-by-construction: path-dispatch is pre-loop, damping mutates `alpha` each sweep, no shared state.

`ACCEL_ANDERSON × *` combinations are NOT orthogonal — Anderson reads `alpha` (set by damping) via the damped-residual mixing and maintains shadow `lf` for the linear path. The 7 explicit tests cover the Anderson interactions: T2 (linear+on+on), T3 (log+on+on), T4 (linear+off+on), T5 (log+off+on). Remaining 20 combinations are `ACCEL_ANDERSON=off` or subsets already covered. Rationale per combination documented in the table caption.

## 8. Merge Gate

All of:
- kk1204 per-iter **wall-clock** ratio ≤ 1.5× raking (post-P2.2, cumulative). Measurement: `/tmp/wu2_kk1204.R` reuses existing wall-clock / iter formula. CTO iter-2 SF2: "per-iter" defined as wall-clock seconds ÷ iterations returned by solver. Post-P1.1 alone, per-iter ratio drops ≈5% (2.17× → ~2.09×); the sweep itself is O(K·M_cell) and dominates per-iter cost at kk1204 K=20. The 1.5× target requires iter reduction from P2.1 damping recovery + P2.2 Anderson acceleration, not just fusion. Post-P1.1 alone gate: ratio strictly improved vs 2.17× baseline; no regression.
- kk1204 reaches RK_OK at tol_abs=1e-3, max_iter=500 (post-P2). Iter-count floor softened (PM iter-3 B3): with `LBW_IEPPA_ACCEL_ANDERSON=on` require iter-to-RK_OK ≤ 400 AND at least 2× fewer iters than the `=off` baseline run on same input. The 2× claim is robust to capacity-gate engagement fraction; the absolute 400 cap is generous to allow initial capacity settling before Anderson fires.
- Unit-mode bounds (PM iter-3): `n_bounds_clamped == 0` on benign unit-mode input (uniform-d, dense cells); `n_bounds_clamped < 0.001 · n` on skewed-d stress input.
- Stepstone errRp ≤ 2.21e-3 ± 1e-4 with no INFEAS (post-P2; must not regress).
- Strict obs-level bounds under `bounds_mode="unit"` on skewed-d input (post-P3).
- Default `bounds_mode="cell"` preserves current weight vectors to **1e-8 tolerance** (CTO iter-1: matching coarsest existing test `test-ieppa-faithful.R:109-136`; 1e-12 was un-testable) on all existing tests.
- Full suite FAIL=0, PASS ≥ 181 at every intermediate commit.

## 9. Deliverables

Four commits, one per WU (P1.1, P2.1, P2.2, P3.1), in pass-then-within-pass order. Each commit includes:
- Source edit.
- Test (or test extension).
- All existing tests pass.
- Stepstone per-pass regression check (errRp ≤ 2.21e-3 ± 1e-4).

**P2.1 intermediate state (CTO iter-1):** P2.1 commit is permitted to leave the kk1204 aggregate-convergence gate RED (RK_OK at tol_abs=1e-3 requires P2.2 Anderson). Interim P2.1-alone assertion: `kk1204 errRp @ 500 iter ≤ 2.5e-3` (prevents adaptive damping regressing convergence relative to current NOCONV plateau near 2.0e-3). P2.2 commit must satisfy the full gate.

## 10.5 Env vars (reference)

Consolidated from WU-2, WU-3, and this design:

| Var | Values | Default | Role |
|---|---|---|---|
| `LBW_IEPPA_FORCE_PATH` | `linear`, `log`, unset | unset → auto | Test-only override of linear-vs-log dispatch |
| `LBW_IEPPA_FORCE_DAMPING` | `on`, `off`, unset | unset → auto | Test-only override of adaptive damping |
| `LBW_IEPPA_ACCEL_ANDERSON` | `on`, `off`, unset | unset → on | Anderson acceleration toggle (was `LBW_IEPPA_ACCEL` in rev 2) |

All three read via `std::getenv` once per solver entry. Microsecond cost. Not documented in user-facing man pages; test-oriented.

## 11. Non-Goals / Out of Scope

- Sum-to-N drift (assessment §4) — handled at wrapper level; not re-opening.
- Configurable `kInfeasPersistence` / exposing streak state to user API.
- Parallel within-margin scatter (B-axis C5) — scope creep; latent ticket only if benchmarks motivate.
- SIMD-specific vectorization (C3) — latent ticket; benchmark-driven only.
- CRAN / PyPI distribution changes.
