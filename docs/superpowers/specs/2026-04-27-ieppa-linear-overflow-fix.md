# iEPPA Linear-Space Overflow Fix + Log-Path Acceleration

**Date**: 2026-04-27 (rev 5 — design-review-gate round-1 fixes)
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

Empirically confirmed (n=1M, iter ~45): individual f_lin values per margin ≈ 5–19 at overflow — far below `kLinearOverflowTrip ≈ 2.1e15`. Overflow fires because the PRODUCT across K=20 margins per cell reaches the trip. Specifically: `X_tilde[c] = X_init[c] × Π_k f_lin[k][g_k(c)] ≈ 2.5e15`.

Prior approach **T1.A** (per-margin geometric mean trigger) was broken: the average of `log f_lin` per margin stays near 3 while the threshold is 17.7, because one hot category grows while others shrink. **T1.A must be removed before T1.B is inserted** — the two mechanisms are not composable (T1.A applies per-margin X_cur scaling while T1.B applies a global scale; simultaneous execution violates both invariants).

---

## Implementation Step 0 (REQUIRED): Remove existing T1.A block

**Before any other change**, delete the T1.A block in `src/ieppa.cpp`. This block is currently at approximately lines 599–651 (the `if (!overflow_trip)` block containing the two-pass product trigger and per-margin renorm loop). Confirm by searching for:
```cpp
// T1.A: Renormalize f_lin per margin when the worst-case X_tilde product
```
Delete from this comment through the closing `}` of the outer `if (!overflow_trip)` block.

**Compile gate** after removal before proceeding.

---

## Fix: Two-Layer Defense

### Layer 1 (primary): T1.B — Linear-path `lf` + `cell_lf` with high-water correction

**Goal**: Track `log(Π_k f_lin[k][g_k(c)])` exactly per cell via the existing `lf[]` vector and a new `cell_lf[]` vector. Detect product overflow via a high-water mark. Correct before overflow fires.

**Invariant**: `X_cur[c] = X_init[c] × W[c] × Π_k f_lin[k][g_k(c)]`

#### Reuse of existing `lf[]` vector

**`lf[]` already exists** at line ~135 as `std::vector<double> lf(n_cats_total, 0.0)`. In the log path, `lf[k][j]` stores the additive log calibration factor for margin k, category j. In the linear path, `lf[k][j]` is currently unused.

T1.B reuses `lf[]` in the linear path as `log(f_lin[k][j])`. This is semantically correct: both paths store the log of the calibration factor for category j of margin k. No new `lf` vector is declared. **No naming collision** — it is the same variable serving a consistent semantic role.

At solver start: `f_lin[] = 1` → `lf[] = log(1) = 0`. The reset at fallback (line ~658: `std::fill(lf.begin(), lf.end(), 0.0)`) is correct for both paths since the log path starts with zero factors.

**Verification that T2.B removal is safe**: The fallback block at line ~658 includes `std::fill(lf.begin(), lf.end(), 0.0)`. A second fallback path at line ~772 also resets `lf[]`. Every fallback path resets `lf[] = 0`. Therefore `cell_lf[] = 0` after reset (since `cell_lf[c] = Σ_k lf[k][g_k(c)] = 0` when `lf = 0`). T2.B row-major initialization is a no-op and is not needed.

#### New `cell_lf[]` vector

```cpp
// Declare alongside f_lin (same scope, after line ~135):
std::vector<double> cell_lf(n_cells_total, 0.0);  // per-cell Σ_k lf[k][g_k(c)]
double cell_lf_hwm = std::numeric_limits<double>::lowest();  // max_c(log_X_init[c] + cell_lf[c])
```

`cell_lf` has size `ct.M_cell` (= `n_cells_total`). Initialize `cell_lf` to all zeros (consistent with `lf = 0` at start). Initialize `cell_lf_hwm` to `std::numeric_limits<double>::lowest()` (effectively −∞); it will be populated to `max(log_X_init[c])` after the first BCD sweep's delta updates. This is safe: the correction block only fires when `cell_lf_hwm >= log(kLinearOverflowThreshold)`, which cannot happen while hwm is −∞.

**`kLinearOverflowThreshold`** is already declared at line ~218 as `sqrt(kLinearOverflowTrip)`. No new declaration needed.

#### Integration into `apply_single_margin_linear(k)`

`apply_single_margin_linear` is a lambda (defined around line ~385) capturing outer scope by `[&]`, including `lf`, `f_lin`, `cells_by_margin_cat`, `cat_offset`, `log_X_init`. `cell_lf` and `cell_lf_hwm` are also in the outer scope and captured by `[&]`.

Inside the per-j loop, **after** `f_lin[cat_offset[k] + j] = new_f` (currently line ~457) and **before** the rescale step, insert:

```cpp
// T1.B: maintain lf and cell_lf for linear-path overflow detection.
// Guard: new_f must be positive and finite (f_lin[k][j] is always > 0
// by construction — update formula is T_kj * W / S_lin > 0 when S_lin > 0
// and T_kj > 0. Guard below catches degenerate buckets).
if (std::isfinite(new_f) && new_f > 1e-300) {
    double lf_new = std::log(new_f);
    double delta = lf_new - lf[cat_offset[k] + j];
    lf[cat_offset[k] + j] = lf_new;
    if (std::fabs(delta) > 1e-12) {
        for (int c : cells_by_margin_cat[cat_offset[k] + j]) {
            cell_lf[c] += delta;
            // Update high-water mark: only grows on positive delta.
            // HWM is monotone-nondecreasing between corrections.
            // Stale-high between negative deltas is INTENTIONAL:
            // it causes early-but-never-missed correction, trading
            // occasional unnecessary corrections for O(1) maintenance.
            if (delta > 0.0) {
                double val = cell_lf[c] + log_X_init[c];
                if (val > cell_lf_hwm) cell_lf_hwm = val;
            }
        }
    }
} else {
    // Degenerate bucket (zero or non-finite f_lin): hold lf constant.
    // X_cur update for this bucket is handled by rescale_lin below.
}
```

#### Post-sweep correction block

Insert **after** the K-margin BCD loop closes (after existing line ~597 `}`), **before** `if (overflow_trip && !linear_fallback_used)`:

```cpp
// T1.B: correct before X_tilde product overflows.
// cell_lf_hwm ≈ max_c log(X_init[c] * Π_k f_lin[k][g_k(c)]) = max_c log(X_tilde[c]).
// Fire when hwm reaches log(kLinearOverflowThreshold) = log(sqrt(kLinearOverflowTrip)).
// After correction, hwm is reset to log(threshold); it repopulates on next positive delta.
if (!overflow_trip && cell_lf_hwm >= std::log(kLinearOverflowThreshold)) {
    double shift = cell_lf_hwm - std::log(kLinearOverflowThreshold);
    // shift > 0 guaranteed by the >= condition above.
    double lf_correction = -shift / static_cast<double>(st.K);
    double x_scale = std::exp(-shift);  // < 1, bounded: x_scale in (0, 1]

    // Distribute shift evenly across all K margins (lf and f_lin).
    for (int k2 = 0; k2 < st.K; k2++) {
        for (int j = 0; j < st.cat_counts[k2]; j++) {
            lf[cat_offset[k2] + j] += lf_correction;
            f_lin[cat_offset[k2] + j] = std::exp(lf[cat_offset[k2] + j]);
        }
    }
    // Update cell_lf and X_cur to maintain invariant X_cur = X_init × W × Π f_lin.
    // W is NOT updated: the scale change propagates through f_lin and X_cur only.
    // Invariant check: X_cur_new = X_cur_old * exp(-shift)
    //   X_init × W × Π f_lin_new = X_init × W × Π f_lin_old × exp(-shift) = X_cur_old × exp(-shift) ✓
    //   X_tilde_new = X_cur_new / W = X_tilde_old × exp(-shift) → max X_tilde = threshold ✓
    for (int c = 0; c < ct.M_cell; c++) {
        cell_lf[c] -= shift;
        X_cur[c] *= x_scale;
    }
    cell_lf_hwm = std::log(kLinearOverflowThreshold);
    // hwm reset to threshold; will repopulate to true max via next positive delta.

    if (st.verbose >= 2) {
        char msg[128];
        std::snprintf(msg, sizeof(msg), "iEPPA T1.B renorm shift=%.2e", shift);
        st.log(msg);
    }
}
```

**W is intentionally not updated.** The invariant `X_cur = X_init × W × Π f_lin` is maintained by the combined `f_lin *= exp(-shift/K)` (across K margins) and `X_cur *= exp(-shift)`. W is the capacity weight from the PREVIOUS iteration's capacity step; it is unchanged and the invariant holds.

---

### Layer 1b: T4.B — Defer X_tilde allocation to log-path entry

**Goal**: Save 8MB allocation when linear path succeeds (T1.B prevents overflow).

Change line ~139 from unconditional allocation to deferred:
```cpp
// Replace: std::vector<double> X_tilde(ct.M_cell);
std::vector<double> X_tilde;  // deferred — allocated at first fallback/log-path use
```

Add guard `if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);` at THREE sites:
1. Start of `if (overflow_trip && !linear_fallback_used)` block (~line 594)
2. Start of `if (overflow_detected)` block (~line 711)
3. Before the log-path X_tilde rebuild loop (~line 728)

Audit all `X_tilde[...]` read sites to confirm no path accesses the vector before these guards. Grep: `grep -n "X_tilde\[" src/ieppa.cpp` — every occurrence must be downstream of one of the three guard sites or the pre-existing overflow-trip path.

---

### Layer 2 (fallback): T2.A + T3.A — Faster log path

#### T2.A: incremental `cell_lf` in log path (free side-effect of T1.B)

`cell_lf` is allocated at solver start by T1.B. At fallback:
```cpp
// Add to fallback reset block (after lf[] and f_lin[] resets):
std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
cell_lf_hwm = std::numeric_limits<double>::lowest();
```

In `apply_single_margin_log(k)` (around line ~464–520), after the two-branch lf assignment (net_log==1.0 fast path OR blended path), capture delta and update cell_lf. Both branches must be covered:

```cpp
// Replace the two-branch lf assignment:
//   if (net_log == 1.0) { lf[off+j] = log_target - log_S_kj; }
//   else { lf_old = lf[off+j]; lf[off+j] = (1-net)*lf_old + net*(log_target-log_S_kj); }
// With:
{
    double lf_old = lf[cat_offset[k] + j];
    double lf_new = (net_log == 1.0)
        ? (log_target - log_S_kj)
        : ((1.0 - net_log) * lf_old + net_log * (log_target - log_S_kj));
    lf[cat_offset[k] + j] = lf_new;
    double delta = lf_new - lf_old;
    if (std::fabs(delta) > 1e-12) {
        for (int c : cells_by_margin_cat[cat_offset[k] + j])
            cell_lf[c] += delta;
    }
}
```

X_tilde rebuild (log path): replace K-inner loop with single sequential exp pass:
```cpp
// Was: for c: for k: s += lf[g_per_cell[k][c]]  (K=20 simultaneous DRAM streams)
// Now:
for (int c = 0; c < ct.M_cell; c++)
    X_tilde[c] = std::exp(log_X_init[c] + cell_lf[c]);  // 1 sequential stream
```

#### T2.B: NOT NEEDED

At fallback, `lf[] = 0` (reset at lines ~658 and ~772) → `cell_lf[] = 0` (reset in fallback block above). No row-major population required.

#### T3.A: Verify SIMD exp() vectorization

```bash
objdump -d $(R_HOME)/library/leafblower/libs/leafblower.so | grep -c "_ZGVdN4v_exp"
```
If count > 0: already vectorized (GCC -O3 -mavx2 -lmvec). If 0: add `#pragma GCC ivdep` before the T2.A exp loop and rebuild.

---

## Data structures

| Variable | Size | Lifetime | New? |
|----------|------|---------|------|
| `lf[]` | total_cats doubles | Solver lifetime | **Reused** — already at line ~135 |
| `cell_lf[]` | M_cell doubles (8MB) | Solver lifetime | NEW — declared after line ~135 |
| `cell_lf_hwm` | 1 double | Solver lifetime | NEW — declared after `cell_lf` |
| `X_tilde` | M_cell doubles (8MB) | Deferred (log-path only) | **Changed** from unconditional |

No new public symbols. No ABI changes.

---

## Performance model

| Scenario | Cost |
|----------|------|
| 80-iter benchmark, T1.B prevents overflow | 80 × 0.107s = **8.6s** ✓ |
| 80-iter benchmark, fallback, T2.A log path | 45 × 0.107s + 35 × 0.365s = **17.6s** ✓ |
| 80-iter baseline (no fix) | 45 × 0.105s + 35 × 0.75s = **31s** ✗ |
| 460-iter full convergence, T1.B | 460 × 0.107s = **49s** |
| 460-iter full convergence, T2.A only | 45 × 0.105s + 415 × 0.365s = **156s** |

T1.B overhead: ~1.5ms/iter for cell_lf maintenance (O(M_cell) bucket-sequential additions) + ~0.1ms/iter amortized correction. Total ~1.5% of 105ms baseline. Figures reference M_cell (unique cells in cell table), not n (rows); for kk1204 M_cell ≈ n = 1M.

---

## Verification

### Regression tests

**Test 1 (new, replaces T1.A test)**: K=20, n=100000, max_weight=3, min_weight=0.2, skewed targets (0.3/0.175×4), max_iterations=200.

Math: at n=100000 each bucket has ~200000/5=20000 cells — no S_lin collapse. Per-factor growth: each BCD sweep multiplies f_lin for the hot category by ≈0.3/0.175 ÷ current_ratio ≈ 1.5x/sweep when far from target. Product across K=20 margins: (1.5)^(20×N) where N = sweep count. At N=5 sweeps: (1.5)^100 ≈ 4×10^17 >> kLinearOverflowTrip ≈ 2.1e15. Product overflows at ~4–5 sweeps WITHOUT T1.B. T1.B fires at ~3 sweeps (when product reaches sqrt(trip)).

**Before T1.B** (T1.A removed, T1.B not yet added): test FAILS — overflow trip fires, status≠0.
**After T1.B**: test PASSES — renorm messages appear at verbose=2, no overflow trip, status=0.

**Test 2 (existing, keep)**: K=20, n=1000 test. After T1.A removal this test also goes RED (overflow fires via product collapse). After T1.B it goes GREEN. Rename to reference T1.B.

Assert for both: `status == 0L`, no "overflow trip" in verbose=1 log.

### Benchmark
`ieppa n=1M K=20 max_weight=3 skewed verbose=2` → wall < 30s, "T1.B renorm" messages present, no "overflow trip".

### Smoke tests
E1/E2 GREEN (existing solver regression suite).

---

## Files changed

- `src/ieppa.cpp` only.
- `tests/testthat/test-calibration-solvers.R` — update/replace T1.A test.
- No ABI changes, no new public symbols.

---

## Out of scope

- **T1.A** — broken (geometric mean trigger cannot detect product overflow). **Must be removed.**
- **T2.B** — zero benefit (lf=0 at fallback makes row-major population a no-op).
- **L1** (S_kj direct maintenance) — deferred; low ROI at kk1204 (87% cells clamped per iter).
- Changes to raking/sinkhorn/chebyshev.
- Python/R API changes.
