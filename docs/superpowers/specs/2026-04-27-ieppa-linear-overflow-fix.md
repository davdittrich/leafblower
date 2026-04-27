# iEPPA Linear-Space Overflow Fix + Log-Path Acceleration

**Date**: 2026-04-27 (rev 2 — design review fixes)
**Status**: Approved for implementation
**Ticket**: leafblower-svbx
**File**: `src/ieppa.cpp` only

---

## Problem

On K=20 / M_cell=1M / max_weight=3 (kk1204 benchmark), iEPPA runs at 0.105s/iter for the first ~45 iterations, then triggers `linear-space overflow trip; fallback to log-space` and slows to ~0.75-1.0s/iter for all remaining iterations.

Root causes:
1. **Linear overflow**: `f_lin[k][j]` grows exponentially on oscillating skewed-target problems. One-shot permanent fallback to log-space.
2. **Log-path bottleneck**: X_tilde rebuild requires K=20 simultaneous sequential streams through `g_per_cell[m]` (20 separate 4MB arrays = 80MB total). Hardware prefetcher saturates → DRAM-bound. ~400ms/iter.

Measured: 37s for 80 iters (first 45 fast, remaining 35 slow). P1 gate target: <30s.

---

## Fix: Two-Layer Defense

### Layer 1 (primary): T1.A — f_lin renormalization

**Goal**: Prevent the overflow that triggers the permanent fallback.

**Invariant**: `X_cur[c] = X_init[c] × W[c] × Π_k f_lin[k][g_k(c)]`

**Insertion point**: AFTER the linear BCD K-margin loop (line ~593, after `for (int k = 0; k < st.K && !overflow_trip; k++)` closes) and BEFORE the P1.1 fused capacity block (line ~677, `double X_tilde_c = X_cur[c] / W[c]`).

At this point X_cur has absorbed W from the PREVIOUS iteration's capacity step. The renorm identity: dividing `f_lin[k][j]` by `c_k` for ALL j of margin k, and multiplying `X_cur[c]` by `c_k`, is exact per-cell (each cell belongs to exactly one j per margin k, so the factor cancels precisely).

**Trigger**: check after every full BCD sweep (K margin updates).

**Algorithm**:
```cpp
// Insert between line ~593 (end of K-loop) and ~677 (P1.1 capacity block)
// Skip if overflow already tripped (fallback will handle it)
if (!overflow_trip) {
    double log_cumul = 0.0;  // track cumulative renorm magnitude
    for (int k = 0; k < st.K; k++) {
        double log_sum = 0.0; int cnt = 0;
        for (int j = 0; j < st.cat_counts[k]; j++) {
            double fkj = f_lin[cat_offset[k] + j];
            // Guard: skip zero/negative/subnormal f_lin (empty/degenerate buckets)
            if (fkj > 1e-300) {
                log_sum += std::log(fkj);
                cnt++;
            }
        }
        if (cnt == 0) continue;
        double log_c = log_sum / cnt;  // log of geometric mean
        // Only renorm if factors are approaching overflow threshold
        if (std::fabs(log_c) < std::log(kLinearOverflowThreshold)) continue;

        // Cumulative guard: prevent K simultaneous renorms from overflowing X_cur
        log_cumul += log_c;
        if (std::fabs(log_cumul) > 700.0) {
            // Product of all c_k would overflow X_cur — fall back gracefully
            overflow_trip = true; break;
        }

        double c_k = std::exp(log_c);
        // Renorm f_lin[k][j] for all j of this margin
        for (int j = 0; j < st.cat_counts[k]; j++)
            f_lin[cat_offset[k] + j] /= c_k;
        // Compensate X_cur: multiply ALL cells by c_k to restore invariant
        for (int c = 0; c < ct.M_cell; c++) X_cur[c] *= c_k;

        if (st.verbose >= 2) {
            char msg[128];
            std::snprintf(msg, sizeof(msg),
                "iEPPA linear renorm k=%d c_k=%.2e (overflow prevention)", k, c_k);
            st.log(msg);
        }
    }
}
```

**kLinearOverflowThreshold**: declare immediately after `kLinearOverflowTrip` (line ~212):
```cpp
// kLinearOverflowTrip already defined above (line ~212)
const double kLinearOverflowThreshold = std::sqrt(kLinearOverflowTrip);
// √trip ≈ (DBL_MAX/2)^(1/(2K)) — halfway in log space before the overflow trip
```
Guard uses `std::fabs(log_c) >= std::log(kLinearOverflowThreshold)`. Guard for degenerate f_lin: `std::isfinite(fkj) && fkj > 1e-300` (catches NaN upstream).

**Mathematical correctness**: At insertion point, X_cur[c] = X_init[c] × W_prev[c] × Π_k f_lin[k][g_k(c)]. After renorm for margin k: f_lin[k][g_k(c)] → f_lin[k][g_k(c)]/c_k, X_cur[c] → X_cur[c] × c_k. Net: X_cur[c] unchanged. ✓

**Cumulative overflow guard**: If K margins all fire renorm in one sweep with c_k near the threshold, the running log sum `log_cumul` is bounded to 700 (exp(700) ≈ 1e304 < DBL_MAX). If it would exceed 700, fall back to log-space (clean fallback, not the panic fallback).

---

### Layer 1b: T4.B — defer X_tilde allocation to log-path entry

**Goal**: Save 8MB allocation when linear path succeeds.

`std::vector<double> X_tilde(ct.M_cell)` at line 139 is allocated unconditionally. Analysis: `apply_single_margin_linear` never reads the `X_tilde` vector; line 677's `X_tilde_c` is a LOCAL scalar, not the vector. The vector is only used in:
- Line 604: fallback reset (inside `if (overflow_trip && !linear_fallback_used)`)
- Lines 543, 549, 627, 634, 655, 662, 728+: log-path rebuild

**Fix**: Change line 139 from unconditional allocation to optional:
```cpp
// Replace: std::vector<double> X_tilde(ct.M_cell);
// With:
std::vector<double> X_tilde;  // deferred — allocate only at log-path/fallback entry
```

Add the guard `if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);` at THREE sites:
1. **Line ~594** — start of `if (overflow_trip && !linear_fallback_used)` block (before line ~604's write)
2. **Line ~711** — start of `if (overflow_detected)` block (before `std::fill(X_tilde...)`)
3. **Log-path entry** — before lines 728+ (X_tilde rebuild)

Both overflow blocks (BCD-loop trip at line ~594 and P1.1 capacity trip at line ~706) are independent paths that can each be reached with an empty deferred X_tilde.

When T1.A prevents overflow, X_tilde is NEVER allocated. Saves 8MB allocation on every `ieppa_solve` call for large M_cell problems.

---

### Layer 2 (fallback): T2.A + T2.B + T3.A — faster log path

For edge cases where log fallback still triggers despite T1.A (extreme sparsity, very tight bounds, structural infeasibility):

**T2.A: Incremental `cell_lf` maintenance**

Maintain `cell_lf[c] = Σ_k lf[k][g_k(c)]` incrementally. Initialized at log-path entry (after fallback at line ~604):
```cpp
// Initialize cell_lf when entering log path (after fallback reset, lf[] = all zeros):
if (cell_lf.empty()) cell_lf.assign(ct.M_cell, 0.0);
// cell_lf[c] = 0 initially (lf is reset to 0 at fallback line 599)
```

In `apply_single_margin_log(k)`: when `lf[k][j]` changes by `delta`, update:
```cpp
for (int c : cells_by_margin_cat[cat_offset[k] + j]) cell_lf[c] += delta;
```

`cells_by_margin_cat` is exhaustive (covers all M_cell cells across all margins) — verified: built at lines 143-147, every cell added to its category for every margin.

X_tilde rebuild: replace K-inner loop with single sequential exp pass:
```cpp
// Was: for c: for k: s += lf[g_per_cell[k][c]]  (K=20 DRAM streams)
// Now:
for (int c = 0; c < ct.M_cell; c++)
    X_tilde[c] = std::exp(log_X_init[c] + cell_lf[c]);  // 1 stream
```

**T2.B: Row-major initial cell_lf population**

When entering log path and `cell_lf` needs full initialization (not from-zero):
```cpp
std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
for (int m = 0; m < st.K; m++)           // one 4MB array per pass
    for (int c = 0; c < ct.M_cell; c++)  // sequential
        cell_lf[c] += lf[cat_offset[m] + ct.g_per_cell[m][c]];
```
3-5× speedup vs current column-major access. Used only at log-path entry — O(1) amortized.

**T3.A: Verify SIMD exp() vectorization**

Check: `objdump -d $(R_HOME)/library/leafblower/libs/leafblower.so | grep -A3 "_ZGVdN4v_exp"`. If absent: add `#pragma GCC ivdep` before the exp loop. With `-O3 -mavx2 -lmvec`, GCC should auto-vectorize. Expected: 3-4× on exp() overhead.

---

## Data structures added

- `std::vector<double> cell_lf(ct.M_cell)` — allocated at log-path entry only.
- `std::vector<double> X_tilde` — changed to deferred allocation (was unconditional).
- No new public symbols. No ABI changes.

---

## Verification

1. **No overflow message**: `harvest(K=20, skewed, verbose=1)` must NOT show "linear-space overflow trip". T1.A verbose≥2 may show "linear renorm" messages (expected — degradation signal for debugging).
2. **Benchmark**: `ieppa n=1M K=20 max_weight=3 skewed` → wall < 30s.
3. **New regression test**: `method="ieppa"`, K=3, skewed targets (0.6/0.2/0.2), n=1000, max_weight=5, max_iterations=200. Assert: `status==0L`, `max_error < 0.01`, no overflow message in verbose=1 log. This exercises multi-margin f_lin growth that would overflow without T1.A.
4. **E1/E2 GREEN** (smoke test for unrelated solvers; insufficient alone for T1.A regression).

---

## Files changed

- `src/ieppa.cpp` only.
- No ABI changes, no new public symbols.

---

## Out of scope

- L1 (S_kj direct maintenance) — deferred due to bounds complexity.
- Changes to raking/sinkhorn/chebyshev.
- Python/R API changes.
