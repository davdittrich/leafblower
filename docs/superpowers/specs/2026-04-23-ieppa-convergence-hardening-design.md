# iEPPA Convergence Hardening Design

**Status:** Draft rev 3 (post design-review-gate iter 2: 3 APPROVED / 2 NEEDS_REVISION → fixes applied)

**Rev 3 blockers resolved:**
- Security B2: overflow trip accounts for `X_init` base weight — `kLinearOverflowTrip = pow(DBL_MAX / (2.0 * max_X_init), 1.0/K)` where `max_X_init = max_c X_init[c]`. Prevents `factor *= f[m]` overflow when survey cell masses exceed 1e6 at high K.
- Security B3: explicit `if (alpha == 1.0)` fast-path in WU-3 — `std::pow(_, 1.0)` is NOT a no-op, costs transcendental overhead. Branch preserves zero-cost stable mode.
- Security B4: WU-1 uses `== kInfeasPersistence` instead of `>= kInfeasPersistence` for `persistent_infeas_pairs.emplace` — insertion fires exactly once per (k,j); subsequent checks no longer re-enter the O(log N) set path.
- CTO B3: dropped `LBW_TEST_HOOKS` compile-time macro. Env var `LBW_IEPPA_FORCE_PATH` always read (cheap, once per solve); no build-time divergence between R CMD INSTALL and CRAN. NOT_CRAN environment at user site gates which CI runs this test, not which binary is shipped.
- CTO B4: stepstone full data is NOT bundled with R package. `.Rbuildignore` excludes `benchmarks/stepstone_fulldata_*`. Integration test calls `testthat::skip_on_cran()` explicitly; manual-run only.
- CTO B5: WU-3 iteration-count assertion relaxed from "1.5-2× ratio" to strict monotonic `iter_damped > iter_stable` (avoids flakiness across OS/BLAS).

**Rev 2 blockers resolved:**

**Rev 2 blockers resolved:**
- Architect B1: WU-3 damping formula corrected to geometric blend (f_new = f_old^(1-α) · naive^α), matching log-space `lf_new = α·Δ + (1-α)·lf_old` — prior formula dropped f_old and converged to wrong marginal.
- Security B1: kLinearOverflowTrip derived at runtime from `pow(DBL_MAX/2.0, 1.0/K)` — prior 1e100 fixed threshold overflowed at K=64 (1e100^63 ≫ DBL_MAX).
- Designer B1: ieppa_damping is internal-only (not exposed in rk_params_t); WU-3 becomes purely auto-triggered per-call stabilization. No API change.
- Designer B2: R/Python wrappers update error-message strings to distinguish persistent vs structural infeasibility.
- CTO B1: WU-2 adds test-only env var `LBW_IEPPA_FORCE_PATH=linear|log` hook (build-time guarded; ignored in release R CMD INSTALL without explicit define).
- CTO B2: merge-gate success criterion rephrased from "errRp < 1e-3 in <30s" (unsupported; raking also stalls at 8.25e-3 on stepstone) to "iEPPA returns RK_OK with errRp no worse than raking's errRp on the same input" (defensible, falsifiable).

**Author:** Dennis Alexis Valin Dittrich
**Date:** 2026-04-23
**Triggers:**
- `leafblower-bel` (P1) — iEPPA INFEAS on stepstone full data despite feasibility
- `leafblower-g8f` (P2) — iEPPA 22× per-iter overhead at M_cell=n on kk1204 regime
- `leafblower-3c2` (P2) — AUTO routing fallback at M_cell/n > 0.9

**Baseline evidence:**
- Stepstone (n=1.58M, K=9, 836 cats, **M_cell/n=0.0183 = 54× compression**): ieppa returns RK_ERR_INFEAS at max_weight ∈ {5, 10, 50}; raking succeeds in 2.3s; lbfgsb NOCONV in 244s.
- kk1204 (n=1M, K=20, cat=5, M_cell/n=1.00): ieppa 369s NOCONV, raking 16s NOCONV, per-iter ratio 22×.

## 1. Problem Statement

Two orthogonal failure modes:

**(A) False-positive infeasibility detection.** Current code latches `is_infeasible = true` on the first occurrence of an empty (margin, category) bucket during Sinkhorn iteration. Under K≥5 overlapping margins, multiplicative updates transiently push some buckets to near-zero mass before subsequent margins re-populate them. The latch triggers on transient conditions, returning `RK_ERR_INFEAS` on structurally feasible problems.

**(B) Per-iteration overhead at dense compression.** At `M_cell/n ≥ 0.5`, log-space Sinkhorn + log-sum-exp stabilization adds transcendental-function cost per cell per margin per iteration with no corresponding numerical-stability benefit (weights stay moderate under realistic bounds). Raking uses linear-space multiplicative updates and is ~22× faster per iter on this regime.

## 2. Scope

Three work units:

| WU | Addresses | Files | Priority |
|---|---|---|---|
| 1 | (A) persistent-infeas tracker | `src/ieppa.cpp` | P1 |
| 2 | (B) linear-space dispatch | `src/ieppa.cpp` | P2 |
| 3 | Hardening: adaptive damping | `src/ieppa.cpp` | P3 |

All WUs single-file; no API changes; no dispatch changes (user-facing `method="ieppa"` semantics preserved).

## 3. Non-scope

- No changes to raking or L-BFGS-B solvers.
- No changes to AUTO routing (separate ticket `leafblower-3c2`).
- No new benchmark harness (uses existing `benchmarks/ieppa_vs_raking_bench.R` from WU-5 of prior plan).
- No public API additions (internal-only refactor).

## 4. WU-1 Persistent-infeas tracker

**Current behavior** (`src/ieppa.cpp` ~lines 111-118, 153-158):
```cpp
if (cells.empty()) { if (st.targets[k][j] > 0.0) { is_infeasible = true; ... } }
// ... and ...
if (log_S_kj < log_threshold) { if (st.targets[k][j] > 0.0) { is_infeasible = true; ... } }
```
`is_infeasible` is a bool, latched once and never cleared. One transient hit → permanent INFEAS status.

**New behavior:** track per-(k,j) persistence. Latch only when the same bucket has been empty for N consecutive outer iterations.

State additions:
```cpp
constexpr int kInfeasPersistence = 5;
// One counter per (margin, category). Flat layout matching cat_offset.
std::vector<int> infeas_streak(total_cats, 0);
std::set<std::pair<int,int>> persistent_infeas_pairs;  // promoted from transient
```

Update logic (replaces current latch sites):
```cpp
auto record_empty = [&](int k, int j) {
    if (st.targets[k][j] <= 0.0) return;  // zero target → no constraint, ignore
    int idx = cat_offset[k] + j;
    infeas_streak[idx]++;
    // Use == not >= so the set insertion fires exactly once per (k,j).
    // Subsequent sweeps with streak > persistence are no-ops here, avoiding
    // O(log N) tree traversal on the hot path.
    if (infeas_streak[idx] == kInfeasPersistence) {
        persistent_infeas_pairs.emplace(k, j);
    }
};

auto record_nonempty = [&](int k, int j) {
    if (st.targets[k][j] <= 0.0) return;
    int idx = cat_offset[k] + j;
    infeas_streak[idx] = 0;  // reset on any non-empty check
};
```

Margin sweep fires `record_empty` on structural-empty (cells.empty) OR log_S_kj below threshold. Fires `record_nonempty` on successful Sinkhorn update.

On solver exit:
```cpp
bool is_infeasible = !persistent_infeas_pairs.empty();
if (is_infeasible && res.status == RK_ERR_NOCONV)
    res.status = RK_ERR_INFEAS;
```

`persistent_infeas_pairs` replaces current `infeasible_pairs` for verbose enumeration (same structure, stricter population criterion).

**Design choices:**
- `kInfeasPersistence = 5`: at default `inner_max_iter=500` this is 1% of the iteration budget. Under K=9 overlapping margins the worst-case settle time for a transient near-zero bucket is ~2-3 outer iters (measured on stepstone during debugging). 5 is comfortably above settle time, well below budget.
- `targets[k][j] <= 0.0` skip: structurally unconstrained categories never count toward infeasibility.
- Reset on non-empty: a bucket that recovers doesn't carry historical transient counts.

**Wrapper error-message update (per designer review iter 1):** both `R/harvest.R` and `python/leafblower/_harvest.py` hardcode the infeasibility message:
```r
stop("leafblower: infeasible problem — empty cell with positive target.")
```
Rev 2 update (same commit as WU-1 source changes):
```r
stop("leafblower: infeasible problem — persistent empty cell with positive target (detected after N consecutive outer iterations).")
```
And equivalent Python `raise RuntimeError(...)`. Distinguishes persistent-iteration infeasibility from input-validation infeasibility (rk_calibrate pre-check). Both cases use `RK_ERR_INFEAS=2`; the message disambiguates.

**Oscillation edge case (per architect review iter 1):** if a bucket oscillates empty ↔ non-empty such that streak resets before reaching kInfeasPersistence, WU-1 will NOT flag INFEAS. On genuinely infeasible inputs with K ≥ 3 this can happen. Outcome: solver hits `max_iter`, returns `RK_ERR_NOCONV` with high errRp. Acceptable: user sees "did not converge" warning rather than "infeasible." A future WU could promote these via a second-order check (e.g., aggregate `sum(infeas_streak) > K * kInfeasPersistence / 2`), but not in this design — keeps WU-1 focused.

**Why not alternative (a) loosen threshold:** `1e-15 → 1e-12` masks transient AND structural infeasibility. The persistence check differentiates them. Strictly better safety.

## 5. WU-2 Linear-space dispatch at dense compression

**Current behavior:** all iterations run in log-space with log-sum-exp stabilization. Per-iter cost at M_cell=n: 2 × M_cell × K `log/exp` calls plus the LSE max-reduction pass.

**New behavior:** at solver entry, compute compression ratio `r = n / M_cell`. If `r < kLinearSpaceThreshold` (i.e., M_cell/n > 1/kLinearSpaceThreshold), dispatch to a linear-space inner loop. Otherwise keep log-space path unchanged.

**Dispatch threshold:** `kLinearSpaceThreshold = 2.0`, i.e. switch to linear-space when M_cell/n > 0.5 (compression ≤ 2×). Below 2× compression, log-space overhead dominates the compression savings; linear-space wins. Above 2× compression, log-space is fine (per-iter cost dominated by O(K) inner loop, not transcendentals).

**Linear-space inner loop:**
```cpp
// State: X[c], X_init[c], W[c], f[k][j]  (all linear, positive)
for (int k = 0; k < st.K; k++) {
    for (int j = 0; j < st.cat_counts[k]; j++) {
        double S_kj = 0.0;
        for (int c : cells_by_margin_cat[cat_offset[k] + j]) {
            // Current X[c] = X_init[c] * W[c] * prod_{m!=k} f[m][g_m(c)]
            double factor = X_init[c] * W[c];
            for (int m = 0; m < st.K; m++) {
                if (m != k) factor *= f[m][cat_offset[m] + ct.g_per_cell[m][c]];
            }
            S_kj += factor;
        }
        // Empty-bucket check identical to log-space path, in linear terms
        if (S_kj < kEmptyBucketThreshold * ct.W_input) { record_empty(k, j); continue; }
        record_nonempty(k, j);
        f[cat_offset[k] + j] = (st.targets[k][j] * ct.W_input) / S_kj;
    }
}
```

**Overflow handling:** track max `f` per sweep. Trip threshold is derived at runtime, accounting for both K (K-way product accumulation) AND the maximum base weight `X_init[c]` (survey cell masses routinely exceed 1e6):
```cpp
// The hot-loop accumulator is: factor = X_init[c] * W[c] * ∏_m f[m]
// W[c] ≤ 1 always (capacity cap). X_init[c] is bounded by total_weight ≤ n.
// So factor ≤ max_X_init * trip^K. Require factor < DBL_MAX / 2:
const double max_X_init = *std::max_element(X_init.begin(), X_init.end());
const double kLinearOverflowTrip = std::pow(
    std::numeric_limits<double>::max() / (2.0 * std::max(max_X_init, 1.0)),
    1.0 / st.K
);
// At K=20, max_X_init=1e6:  ~1e14
// At K=64, max_X_init=1e6:  ~5.5e4
// Always strictly less than (DBL_MAX / max_X_init)^(1/K); the K-way product
// with X_init prefactor cannot reach DBL_MAX/2 before this trip fires.
```

On trip: abort linear-space path, fully reset solver state (`lf`, `W`, `X_tilde`, `X`, `infeas_streak`, `persistent_infeas_pairs`, `iter` counter), re-initialize from `X_init`, switch to log-space, continue from iter 0. One-shot per solve (a boolean `linear_fallback_used` prevents re-entry). Rare; provides safety without ongoing overhead.

**State-clean fallback checklist** (per architect review iter 1): the reset must clear every piece of solver state set by the linear-space path. Explicit enumeration:
- `lf[]` zeroed (log-space path will reinitialize from `X_init`)
- `W[c]` reset to 1.0
- `X_tilde[c]` / `X[c]` reset to `X_init[c]`
- `infeas_streak[idx]` reset to 0 for all idx
- `persistent_infeas_pairs` cleared
- `iter` counter reset to 0 (outer loop restarts; the `inner_max_iter` budget is re-granted)
- `alpha` reset to 1.0 (stable mode; WU-3 re-triggers if streak re-builds)

**Why not always-linear:** at extreme compression (M_cell ≪ n) with tight bounds, `f[k][j]` can grow unboundedly between sweeps before the capacity block catches up. Log-space handles this gracefully; linear-space overflows. The ratio threshold selects the regime where log-space cost isn't justified.

## 6. WU-3 Adaptive damping (hardening)

**Motivation:** even after WU-1, occasional transient near-zero cells still occur. Damping each Sinkhorn update reduces the amplitude of transients.

**Corrected formula (rev 2, per architect review):** Damping is a geometric blend of the previous factor with the naive full step — NOT `pow(naive, alpha)` alone, which would discard `f_old` and converge to a wrong marginal. Correct linear-space:
```cpp
// Instead of:
// f[k][j] = (target * W) / S_kj
// Use:
double f_old = exp(lf[cat_offset[k] + j]);    // or stored directly in linear path
double naive = (st.targets[k][j] * ct.W_input) / S_kj;
// Geometric blend: at alpha=1 full step (f_new = naive), at alpha=0 no step (f_new = f_old).
double f_new = std::pow(f_old, 1.0 - alpha) * std::pow(naive, alpha);
// Equivalent in log-space (preferred form for the existing log-space path):
// lf_new = (1 - alpha) * lf_old + alpha * (log_target - log_S_kj)
```

Standard Sinkhorn stabilization (Peyré-Cuturi 2019 §4.4 equation 4.52 "softened updates").

**Auto-triggering (no API addition):** ieppa_damping remains a solver-internal constant. Default `alpha = 1.0`. When WU-1's `infeas_streak[k][j]` reaches `kInfeasPersistence / 2` (mid-streak), auto-enable `alpha = 0.5` for subsequent sweeps for THAT solve. Resets between solves. Two states:

- **Stable mode (default):** alpha=1.0, byte-identical to pre-WU-3 behavior. Applied until any bucket hits persistence/2.
- **Damped mode (auto):** alpha=0.5, applies to ALL margin updates (simplest). Triggered when transient stress detected.

Transition latched per-solve; does not revert even if streak resets (prevents oscillation between modes).

**Cost:** two extra `std::pow` calls per (margin, category) when damped. **Stable-mode must explicitly branch on `alpha == 1.0`** — `std::pow(_, 1.0)` does NOT short-circuit and pays full transcendental cost:
```cpp
if (alpha == 1.0) {
    lf[cat_offset[k] + j] = log_target - log_S_kj;  // fast path
} else {
    lf[cat_offset[k] + j] = (1.0 - alpha) * lf[cat_offset[k] + j]
                          + alpha * (log_target - log_S_kj);  // damped
}
```
Convergence slows by ~2× in damped mode; acceptable because damped mode only engages when WU-1 is already headed for NOCONV without it.

**Why separate from WU-1:** independent mechanism. WU-1 fixes the latching logic (what gets reported); WU-3 reduces the frequency of transients in the first place (makes WU-1 less likely to fire). Additive benefits, composed via the `infeas_streak` state WU-1 already maintains.

## 7. Testing

**WU-1 test:** new `tests/testthat/test-ieppa-persistent-infeas.R`:
- Structurally feasible input that triggers current false-positive (stepstone-like: K=5, 2-4 cats each, interaction margin forces transient near-zero)
- Pre-fix: returns RK_ERR_INFEAS
- Post-fix: converges to RK_OK with errRp < tol_abs
- Truly infeasible input (one target cell has no matching observations): returns RK_ERR_INFEAS on both pre- and post-fix (regression guard for genuine detection)

**WU-2 test:** extend `tests/testthat/test-ieppa-faithful.R`:
- Dense regime test (M_cell/n ≈ 1, n=10000): confirm result weights identical (to 1e-8) between log-space and linear-space paths.
  **Force-path mechanism (per CTO review iter 1, refined iter 2):** env var `LBW_IEPPA_FORCE_PATH` read once at solver entry. Values: `linear`, `log`, or unset (auto-dispatch). Setting forces one path regardless of M_cell/n. Always compiled in (no `#ifdef LBW_TEST_HOOKS` — that macro cannot differentiate R CMD INSTALL from CRAN builds cleanly). Cost: one `getenv` call per solve (microseconds; amortized over O(iter·M_cell·K) work). Documented as test-only in the source comment; not mentioned in public docs. CRAN packaging is unaffected — env var remains dormant unless explicitly set.
- Sparse regime test (M_cell/n ≈ 0.01, n=10000): confirm log-space path selected (verify via verbose=1 log line), no regression.
- Overflow synthesis (very tight bounds, high K): confirm linear-space fallback to log-space fires, final result correct. Force `LBW_IEPPA_FORCE_PATH=linear` at start, observe fallback trip + completion.

**WU-3 test:** in `tests/testthat/test-ieppa-faithful.R`:
- Default stable mode (no persistence stress): byte-identical to current behavior (regression check)
- Damped mode (triggered via engineered persistence stress): converges to same solution, strictly slower. **Assertion: `iter_damped > iter_stable`** (monotone, not ratio-bounded — per CTO review iter 2, ratio assertions flake across OS/BLAS).

**Integration test (stepstone full data):** `testthat::skip_on_cran()` required. Data file `benchmarks/stepstone_fulldata_*` excluded from CRAN tarball via `.Rbuildignore`. Test is developer-only; runs locally + in non-CRAN CI.

**Integration test:** rerun stepstone full data (n=1.58M, K=9):
- Pre-fix: iEPPA returns INFEAS
- Post-WU-1: iEPPA converges (target: errRp < 1e-3 at < 30s wall-clock given 54× compression)
- Post-WU-2: either not triggered (sparse regime, log-space stays) or identical result

## 8. Convergence guarantees

- **WU-1**: preserves all existing termination conditions; stricter infeasibility criterion. Never latches on genuinely feasible input. For truly infeasible inputs, latches after exactly `kInfeasPersistence` outer iterations instead of iter 1 — delay is bounded, correctness preserved.
- **WU-2**: linear-space and log-space paths produce mathematically identical iterates modulo floating-point roundoff. Overflow fallback guaranteed to terminate (one-shot). No convergence regression.
- **WU-3**: at `alpha ∈ (0, 1]` Sinkhorn convergence preserved (Peyré-Cuturi 2019 Proposition 4.8); rate slows by 1/alpha factor. Default alpha=1.0 ⇒ byte-identical.

## 9. Error Handling

| Condition | Code | Detection |
|---|---|---|
| Persistent empty bucket (WU-1) | `RK_ERR_INFEAS` | After `kInfeasPersistence` consecutive transient hits |
| Transient empty bucket | non-fatal | Skip update this sweep; reset streak on next non-empty |
| Linear-space overflow (WU-2) | non-fatal | One-shot fallback to log-space; continue |
| Invalid damping (WU-3) | `RK_ERR_BADARG` | `ieppa_damping` outside (0, 1] |

## 10. Non-Goals

- Routing changes (`leafblower-3c2`)
- Full re-benchmark (leafblower-w5d's BLSE already scheduled)
- Raking changes (separate WU list)
- L-BFGS-B changes

## 11. Deliverables

Three commits (one per WU), each with:
- Source edit in `src/ieppa.cpp`
- Unit test in `tests/testthat/`
- All existing tests pass (`FAIL 0`)
- Integration verification on stepstone data (WU-1 only — post-fix must succeed)

Post-merge acceptance (rev 2: rephrased per CTO review iter 1 — prior claim "errRp < 1e-3" unsupported given raking also stalls at 8.25e-3 on stepstone):
- **Stepstone full data**: `method="ieppa"` no longer returns `RK_ERR_INFEAS` (false-positive eliminated). Returned errRp is **no worse than raking's errRp on the same input** (raking: 8.25e-3; iEPPA must be ≤ this, same or stricter). Wall-clock ≤ 2× raking's 2.3s on this data (i.e., ≤5s target).
- **kk1204 regime** (M_cell/n=1.0): `method="ieppa"` per-iter cost within 2× of raking (vs current 22×) when linear-space dispatch fires.
- **All existing test suite**: `[FAIL 0 | PASS ≥ 170]`.
- **No false-positive INFEAS regressions**: on the existing 169-test suite, any test previously returning RK_OK must still return RK_OK (no new RK_ERR_INFEAS from tightened persistence logic).
