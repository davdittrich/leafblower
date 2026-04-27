# iEPPA Linear-Space Overflow Fix + Log-Path Acceleration

**Date**: 2026-04-27 (rev 4 — T1.B replaces T1.A; T2.B removed)
**Status**: Pending design review
**Ticket**: leafblower-svbx
**File**: `src/ieppa.cpp` only

---

## Problem

On K=20 / M_cell=1M / max_weight=3 (kk1204 benchmark), iEPPA runs at 0.105s/iter for the first ~45 iterations, then triggers `linear-space overflow trip; fallback to log-space` and slows to ~0.75-1.0s/iter for all remaining iterations.

Root causes:
1. **Linear overflow**: The PRODUCT `Π_k f_lin[k][g_k(c)]` across K=20 margins per cell grows exponentially on oscillating skewed-target problems. Individual f_lin values are small (~5–19 at overflow) but their product reaches `kLinearOverflowTrip ≈ 2.1e15`. One-shot permanent fallback to log-space.
2. **Log-path bottleneck**: X_tilde rebuild requires K=20 simultaneous sequential streams through `g_per_cell[m]` (20 separate 4MB arrays = 80MB total). Hardware prefetcher saturates → DRAM-bound. ~400ms/iter.

Measured: 37s for 80 iters (first 45 fast, remaining 35 slow). P1 gate target: <30s.

### Root cause detail: product overflow, not individual overflow

Empirically confirmed (n=1M, iter ~45): individual f_lin values per margin ≈ 5–19 at overflow — far below `kLinearOverflowTrip ≈ 2.1e15`. The overflow fires because the PRODUCT across K=20 margins per cell reaches the trip. Specifically `X_tilde[c] = X_init[c] × Π_k f_lin[k][g_k(c)] ≈ 2.5e15`.

Prior approach T1.A (per-margin geometric mean trigger) was broken: the average of `log f_lin` per margin stays near 3 while the threshold is 17.7, because one category grows while others shrink. The geometric mean approach is fundamentally incapable of detecting product overflow from per-margin averages.

---

## Fix: Two-Layer Defense

### Layer 1 (primary): T1.B — Linear-path cell_lf with high-water correction

**Goal**: Track `log(Π_k f_lin[k][g_k(c)])` exactly per cell via `cell_lf[c]`, detect product overflow via a high-water mark, and correct before overflow fires.

**Invariant**: `X_cur[c] = X_init[c] × W[c] × Π_k f_lin[k][g_k(c)]`

**Data structures** (declared alongside existing `f_lin` and `X_cur`):
```cpp
std::vector<double> lf;       // shadow: lf[cat_offset[k]+j] = log(f_lin[k][j])
std::vector<double> cell_lf;  // per-cell: cell_lf[c] = Σ_k lf[k][g_k(c)]
double cell_lf_hwm = -std::numeric_limits<double>::infinity();  // high-water mark
```

Allocated at solver start (same site as `f_lin`). `lf` has same size as `f_lin` (total_cats doubles). `cell_lf` has size `ct.M_cell`.

`log_X_init[c]` is already maintained in the existing code (used by the log path). `cell_lf_hwm` tracks `max_c(log_X_init[c] + cell_lf[c])`, i.e., the log of the worst-case `X_tilde` product.

**Integration into `apply_single_margin_linear(k)`**

After updating `f_lin[cat_offset[k] + j]`, capture the log delta and propagate to `cell_lf`:
```cpp
// After: f_lin[cat_offset[k] + j] = new_f_lin  (existing update)
double lf_new = std::log(new_f_lin);
double delta = lf_new - lf[cat_offset[k] + j];
lf[cat_offset[k] + j] = lf_new;
if (std::fabs(delta) > 1e-12) {
    for (int c : cells_by_margin_cat[cat_offset[k] + j]) {
        cell_lf[c] += delta;
        if (delta > 0.0) {
            double val = cell_lf[c] + log_X_init[c];
            if (val > cell_lf_hwm) cell_lf_hwm = val;
        }
    }
}
```

**Post-sweep correction block** (insert AFTER the K-margin BCD loop, BEFORE `if (overflow_trip && !linear_fallback_used)`):
```cpp
// T1.B: correct before X_tilde product overflows.
// cell_lf_hwm = max_c log(X_init[c] * Π_k f_lin[k][g_k(c)]) ≈ log(max X_tilde).
// Fire when hwm reaches log(kLinearOverflowThreshold) = log(sqrt(kLinearOverflowTrip)).
if (!overflow_trip && cell_lf_hwm >= std::log(kLinearOverflowThreshold)) {
    double shift = cell_lf_hwm - std::log(kLinearOverflowThreshold);
    double lf_correction = -shift / static_cast<double>(st.K);
    double x_scale = std::exp(-shift);

    // Distribute shift evenly across all K margins (lf and f_lin)
    for (int k = 0; k < st.K; k++) {
        for (int j = 0; j < st.cat_counts[k]; j++) {
            lf[cat_offset[k] + j] += lf_correction;
            f_lin[cat_offset[k] + j] = std::exp(lf[cat_offset[k] + j]);
        }
    }
    // Update cell_lf and X_cur to maintain invariant X_cur = X_init × W × Π f_lin
    for (int c = 0; c < ct.M_cell; c++) {
        cell_lf[c] -= shift;
        X_cur[c] *= x_scale;
    }
    cell_lf_hwm = std::log(kLinearOverflowThreshold);

    if (st.verbose >= 2) {
        char msg[128];
        std::snprintf(msg, sizeof(msg), "iEPPA T1.B renorm shift=%.2e", shift);
        st.log(msg);
    }
}
```

**Mathematical correctness**:

`cell_lf[c] = Σ_k lf[k][g_k(c)] = Σ_k log(f_lin[k][g_k(c)]) = log(Π_k f_lin[k][g_k(c)])`.

After correction with `shift = cell_lf_hwm - log(threshold)`:
- `lf[k][j] += -shift/K` → `f_lin[k][j] *= exp(-shift/K)` for all k, j
- `Π_k f_lin_new[k][g_k(c)] = Π_k f_lin_old[k][g_k(c)] × exp(-shift)` (K terms, each ×exp(-shift/K))
- `X_cur_new[c] = X_cur_old[c] × exp(-shift)`
- `X_init × W × Π f_lin_new = X_init × W × Π f_lin_old × exp(-shift) = X_cur_old × exp(-shift) = X_cur_new` ✓
- `X_tilde_new = X_cur_new / W = X_tilde_old × exp(-shift)` — decreases to threshold ✓

**kLinearOverflowThreshold**: declare immediately after `kLinearOverflowTrip` (line ~212):
```cpp
const double kLinearOverflowThreshold = std::sqrt(kLinearOverflowTrip);
// Trigger at sqrt(trip): gives ~factor-of-2 headroom in log space for the correction to fire.
```

**Initialization**: At solver start, `lf[k][j] = log(f_lin[k][j])` for all k,j (f_lin initialized to 1 → lf = 0). `cell_lf[c] = 0` for all c (since lf = 0). `cell_lf_hwm = log(max_X_init)` (initial X_tilde before any f_lin growth).

**Overhead**:
- `cell_lf` maintenance: O(M_cell) additions per BCD sweep (bucket-sequential, cache-local). ~1.5ms for n=1M.
- `cell_lf_hwm` tracking: O(1) per bucket update, branch-predicted (delta > 0 fires only when f_lin grows).
- Correction: O(K×J) f_lin/lf updates + O(M_cell) cell_lf+X_cur updates. ~2ms per event. Fires every ~20 iters. Amortized ~0.1ms/iter.
- **Total overhead: ~1.6ms/iter (~1.5% of 105ms linear-path baseline).**

---

### Layer 1b: T4.B — Defer X_tilde allocation to log-path entry

**Goal**: Save 8MB allocation when linear path succeeds.

`std::vector<double> X_tilde(ct.M_cell)` at line ~139 is allocated unconditionally. `apply_single_margin_linear` never reads the `X_tilde` vector (line 677's `X_tilde_c` is a LOCAL scalar). The vector is only used at fallback reset and in the log path.

**Fix**: Change line ~139 from unconditional allocation to deferred:
```cpp
// Replace: std::vector<double> X_tilde(ct.M_cell);
std::vector<double> X_tilde;  // deferred — allocated at first fallback/log-path use
```

Add guard `if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);` at THREE sites:
1. Start of `if (overflow_trip && !linear_fallback_used)` block (~line 594)
2. Start of `if (overflow_detected)` block (~line 711)
3. Before the log-path X_tilde rebuild loop (~line 728)

When T1.B prevents overflow, X_tilde is never allocated. Saves 8MB on every `ieppa_solve` call.

---

### Layer 2 (fallback): T2.A + T3.A — Faster log path

**T2.A is a free side-effect of T1.B.** `cell_lf` is allocated and maintained from iteration 1. When linear→log fallback triggers:

1. Existing fallback code resets `lf[] = 0`.
2. Add: `std::fill(cell_lf.begin(), cell_lf.end(), 0.0); cell_lf_hwm = -inf;` in fallback reset block.
3. `cell_lf[c] = 0` is consistent with `lf = 0` — no row-major initialization needed.
4. `apply_single_margin_log(k)` updates `cell_lf` incrementally (delta capture, same pattern as T1.B).
5. X_tilde rebuild uses single exp pass (eliminates K=20 DRAM streams):
```cpp
// Was: for c: for k: s += lf[g_per_cell[k][c]]  (K=20 simultaneous streams, DRAM-bound)
for (int c = 0; c < ct.M_cell; c++)
    X_tilde[c] = std::exp(log_X_init[c] + cell_lf[c]);  // 1 sequential stream
```

**T2.B is NOT needed.** T2.B was designed to populate `cell_lf` at log-path entry from non-zero `lf` values. At fallback, `lf[]` is always reset to 0, so `cell_lf[]` is trivially reset to 0. T2.B would compute `Σ_k 0 = 0` — a no-op. Removed from scope.

**T3.A: Verify SIMD exp() vectorization**

Check: `objdump -d $(R_HOME)/library/leafblower/libs/leafblower.so | grep -c "_ZGVdN4v_exp"`. If count > 0: already vectorized (GCC -O3 -mavx2 -lmvec auto-vectorizes). If absent: add `#pragma GCC ivdep` before the T2.A exp loop.

---

## Data structures added

- `std::vector<double> lf` (total_cats doubles) — shadow log of f_lin; same lifetime as f_lin.
- `std::vector<double> cell_lf` (M_cell doubles = 8MB) — per-cell log product; allocated at solver start.
- `double cell_lf_hwm` — scalar high-water mark; initialized to `log(max_X_init)`.
- `std::vector<double> X_tilde` — changed to deferred allocation (was unconditional 8MB).
- No new public symbols. No ABI changes.

---

## Performance model

| Scenario | Cost |
|----------|------|
| 80-iter benchmark, T1.B prevents overflow | 80 × 0.107s = **8.6s** ✓ |
| 80-iter benchmark, T1.B fails, T2.A fallback | 45 × 0.107s + 35 × 0.365s = **17.6s** ✓ |
| 80-iter baseline (no fix) | 45 × 0.105s + 35 × 0.75s = **31s** ✗ |
| 460-iter full convergence, T1.B | 460 × 0.107s = **49s** |
| 460-iter full convergence, T2.A only | 45 × 0.105s + 415 × 0.365s = **156s** |

T1.B is 2× faster than T2.A-only for 80 iters, 3× faster for full convergence.

---

## Verification

1. **T1.B fires**: `harvest(K=20, n=1M, skewed, verbose=2)` shows "iEPPA T1.B renorm" messages; NO "linear-space overflow trip".
2. **Benchmark**: wall < 30s on kk1204 (K=20, n=1M, max_weight=3).
3. **Regression test**: K=20, n=100000, max_weight=3, skewed targets, max_iterations=200. Assert: `status==0L`, no overflow message in verbose=1 log. (n=100000 exercises product overflow mechanism without S_lin collapse of n=1000.)
4. **E1/E2 GREEN** (smoke tests for unrelated solvers).

---

## Files changed

- `src/ieppa.cpp` only.
- No ABI changes, no new public symbols.

---

## Out of scope

- T1.A (broken — geometric mean trigger cannot detect product overflow).
- T2.B (zero benefit — lf=0 at log-path entry makes row-major population a no-op).
- L1 (S_kj direct maintenance) — deferred; low ROI for kk1204 (87% cells clamped per iter).
- Changes to raking/sinkhorn/chebyshev.
- Python/R API changes.
