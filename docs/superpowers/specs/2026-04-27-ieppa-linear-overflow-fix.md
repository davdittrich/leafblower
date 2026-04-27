# iEPPA Linear-Space Overflow Fix + Log-Path Acceleration

**Date**: 2026-04-27
**Status**: Approved for implementation
**Ticket**: leafblower-svbx
**File**: `src/ieppa.cpp` only

---

## Problem

On K=20 / M_cell=1M / max_weight=3 (kk1204 benchmark), iEPPA runs at 0.105s/iter for the first ~45 iterations, then triggers `linear-space overflow trip; fallback to log-space` and slows to ~0.75-1.0s/iter for all remaining iterations.

Root causes (two layers):
1. **Linear overflow**: `f_lin[k][j]` grows exponentially on oscillating skewed-target problems. One-shot fallback to log-space is permanent.
2. **Log-path bottleneck**: X_tilde rebuild requires K=20 simultaneous sequential streams through `g_per_cell[m]` (20 separate 4MB arrays = 80MB total). Hardware prefetcher saturates → DRAM-bound. ~400ms/iter just for the rebuild.

Measured: 37s for 80 iters (first 45 fast, remaining 35 slow). P1 gate target: <30s for full convergence.

---

## Fix: Two-Layer Defense

### Layer 1 (primary): T1.A — f_lin renormalization

**Prevent the overflow that triggers the fallback.**

After each full BCD sweep (K margin updates), check for overflow risk in `f_lin`. When detected, renormalize per-margin: divide all `f_lin[k][j]` by the geometric mean of that margin's active factors, and compensate `X_cur` accordingly.

**Trigger condition**: `max(f_lin[k][j]) > kOverflowThreshold` (= 1e7) OR `min(f_lin[k][j]) < 1/kOverflowThreshold`. Check at end of every BCD sweep.

**Renormalization**:
```cpp
// After BCD sweep — check each margin
for (int k = 0; k < st.K; k++) {
    double log_sum = 0.0; int cnt = 0;
    for (int j = 0; j < st.cat_counts[k]; j++) {
        if (bucket_kj > kEmptyBucketThreshold) {
            log_sum += std::log(std::fabs(f_lin[cat_offset[k]+j]));
            cnt++;
        }
    }
    if (cnt == 0) continue;
    double c_k = std::exp(log_sum / cnt);  // geometric mean
    if (c_k < 1.0/kOverflowThreshold || c_k > kOverflowThreshold) {
        for (int j = 0; j < st.cat_counts[k]; j++)
            f_lin[cat_offset[k]+j] /= c_k;
        // Compensate X_cur: multiply all cells by c_k
        for (int c = 0; c < ct.M_cell; c++) X_cur[c] *= c_k;
    }
}
```

**Mathematical correctness**: X_cur = X_init × W × Π_k f_lin[k][g_k(c)]. Dividing one margin's factors by c_k and multiplying X_cur by c_k is identity. Calibration result unchanged.

**Frequency**: Check every BCD sweep (K calls to `apply_single_margin_linear`). Renormalization fires rarely — only when factors approach overflow. When it fires: 2 × O(M_cell) sequential passes = ~5ms, negligible.

**Expected**: linear overflow never occurs → solver stays in linear path → 0.105s/iter permanently.

---

### Layer 1b (secondary): T4.B — skip X_tilde when M_cell=n

**Eliminate redundant 8MB array in linear path.**

When `M_cell == st.n` (every obs is in its own cell, i.e., no compression), `X_tilde[c]` is identical to `X_cur[c]`. The capacity step uses `X[c] = clamp(X_tilde[c] × W[c], L_cell[c], U_cell[c])` — but `X_tilde[c] = X_cur[c]` when `n_per_cell[c] = 1` for all c.

**Change**: In `apply_single_margin_linear(k)`, where `X_tilde` is referenced after the BCD update, replace `X_tilde[c]` with `X_cur[c]` when M_cell == st.n. Eliminate the `std::vector<double> X_tilde(ct.M_cell)` allocation when `M_cell == st.n`. 

Saves 8MB allocation + one sequential O(M_cell) copy per capacity step.

---

### Layer 2 (fallback): T2.A + T2.B + T3.A — faster log path

For edge cases where log fallback still triggers (extreme sparsity, structural infeasibility, very tight bounds):

**T2.A: Incremental `cell_lf` maintenance**

Instead of rebuilding `X_tilde` from K=20 arrays each iteration:
```cpp
// Maintain cell_lf[c] = Σ_k lf[k][g_k(c)] — updated when lf changes
// When lf[k][j] changes by delta: cell_lf[c] += delta for c in bucket(k,j)
// X_tilde rebuild: X_tilde[c] = exp(log_X_init[c] + cell_lf[c])
```

The rebuild is now ONE sequential pass (no K-inner loop), eliminating the K=20-stream DRAM problem. X_tilde rebuild cost: O(M_cell) exp() = ~15ms vs ~400ms current.

Add `std::vector<double> cell_lf(ct.M_cell)` initialized as `Σ_k lf[k][g_k(c)]` at log-path entry. Update incrementally in `apply_single_margin_log(k)`.

**T2.B: Row-major loop for full rebuild**

When full rebuild is needed (initial entry into log path, after homotopy level reset):
```cpp
// Row-major: one sequential pass per margin (not K simultaneous streams)
std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
for (int m = 0; m < st.K; m++)           // sequential per-margin pass
    for (int c = 0; c < ct.M_cell; c++)  // single 4MB array per pass
        cell_lf[c] += lf[cat_offset[m] + ct.g_per_cell[m][c]];
// Then: X_tilde[c] = exp(log_X_init[c] + cell_lf[c])
```

3-5× speedup vs current column-major access.

**T3.A: SIMD exp() verification**

Confirm GCC auto-vectorizes the exp loop with `-O3 -mavx2 -lmvec`. Check via:
```bash
objdump -d $(R_HOME)/library/leafblower/libs/leafblower.so | grep -A5 "call.*exp"
```

If not vectorized: add `#pragma GCC ivdep` before the exp loop over M_cell cells. libmvec (`_ZGVdN4v_exp`) handles 4 doubles simultaneously → 4× on exp overhead.

---

## Data structures added

- `std::vector<double> cell_lf(ct.M_cell)` — only when in log path. Initialized at fallback entry, maintained incrementally.

---

## Verification

1. Verbose log must NOT show "linear-space overflow trip" for kk1204 K=20 after T1.A.
2. Benchmark: `ieppa n=1M K=20 max_weight=3 skewed` → wall < 30s.
3. E1/E2 correctness tests remain GREEN.
4. Per-iter cost stays ~0.105s throughout 500 iters.

---

## Files changed

- `src/ieppa.cpp` only.
- No ABI changes, no new public symbols.

---

## Out of scope

- L1 (S_kj direct maintenance) — deferred due to bounds complexity.
- Changes to raking/sinkhorn/chebyshev.
- Python/R API changes.
