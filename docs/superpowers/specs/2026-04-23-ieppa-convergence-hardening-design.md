# iEPPA Convergence Hardening Design

**Status:** Draft (pre design-review-gate)
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
    if (infeas_streak[idx] >= kInfeasPersistence) {
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

**Overflow handling:** track max `f` per sweep. If `max_f > kLinearOverflowTrip = 1e100`, abort linear-space path, re-initialize from `X_init` in log-space, continue in log-space from iter 0. One-shot fallback per solve. Rare; provides safety without ongoing overhead.

**Why not always-linear:** at extreme compression (M_cell ≪ n) with tight bounds, `f[k][j]` can grow unboundedly between sweeps before the capacity block catches up. Log-space handles this gracefully; linear-space overflows. The ratio threshold selects the regime where log-space cost isn't justified.

## 6. WU-3 Adaptive damping (hardening)

**Motivation:** even after WU-1, occasional transient near-zero cells still occur. Damping each Sinkhorn update reduces the amplitude of transients:
```cpp
// Instead of:
// f[k][j] = (target * W) / S_kj
// Use:
const double alpha = 0.8;
double naive = (target * W) / S_kj;
f[k][j] = std::pow(naive, alpha);  // damped step
```

At `alpha = 1.0` equivalent to current; at `alpha = 0.5` a half-step; at `alpha → 0` no update. Standard Sinkhorn stabilization technique (see Peyré-Cuturi 2019 §4.4).

**Cost:** one extra `std::pow` call per (margin, category) per iteration. Cheap. Convergence slows by factor proportional to 1/alpha — at alpha=0.8, ~25% more iterations required to reach same tol.

**Default:** `alpha = 1.0` (current behavior). Enable via `rk_params_t` new field `ieppa_damping` (double, default 1.0). Users who hit edge cases can set 0.5-0.8 for stability.

**Why separate from WU-1:** independent mechanism. WU-1 fixes the latching logic; WU-3 reduces the frequency of transients in the first place. Additive benefits.

## 7. Testing

**WU-1 test:** new `tests/testthat/test-ieppa-persistent-infeas.R`:
- Structurally feasible input that triggers current false-positive (stepstone-like: K=5, 2-4 cats each, interaction margin forces transient near-zero)
- Pre-fix: returns RK_ERR_INFEAS
- Post-fix: converges to RK_OK with errRp < tol_abs
- Truly infeasible input (one target cell has no matching observations): returns RK_ERR_INFEAS on both pre- and post-fix (regression guard for genuine detection)

**WU-2 test:** extend `tests/testthat/test-ieppa-faithful.R`:
- Dense regime test (M_cell/n ≈ 1, n=10000): confirm result weights identical (to 1e-8) between log-space and linear-space paths
- Sparse regime test (M_cell/n ≈ 0.01, n=10000): confirm log-space path selected, no regression
- Overflow synthesis (very tight bounds, high K): confirm linear-space fallback to log-space fires, final result correct

**WU-3 test:** in `tests/testthat/test-ieppa-faithful.R`:
- Default alpha=1.0: byte-identical to current behavior (regression check)
- alpha=0.5: converges to same solution, slower (expect 1.5-2× iteration count)
- alpha=0.0: returns RK_ERR_BADARG (degenerate, no progress possible)

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

Post-merge acceptance:
- Stepstone full data: `method="ieppa"` returns `RK_OK` with `errRp < 1e-3`
- kk1204 regime: `method="ieppa"` per-iter cost within 2× of raking (vs current 22×)
- All existing test suite: `[FAIL 0 | PASS ≥ 170]`
