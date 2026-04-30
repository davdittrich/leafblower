# Code Review: C++ Source Files (cell_table, chebyshev, sinkhorn, grake, greg)

**Date:** 2026-04-28  
**Scope:** src/cell_table.cpp/hpp, chebyshev.cpp (first 200 lines), sinkhorn.cpp, grake.cpp/hpp, greg.cpp  
**Criteria:** Memory allocation, division patterns, dead code, integer overflow, edge cases, status codes

---

## Issues Found (10 total, ranked by severity)

### HIGH SEVERITY

1. **chebyshev.cpp:77 — Division by zero risk (n_d)**
   - **Location:** Line 77: `T_flat[m] = Tgt[m] / n_d`
   - **Problem:** `n_d = static_cast<double>(st.n)` is never checked for zero. If `st.n == 0` (edge case), division is undefined.
   - **Severity:** High
   - **Fix:** Add guard after line 56: `if (st.n <= 0) { res.status = RK_ERR_BADARG; return res; }`
   - **Confidence:** 95 (exact code match)

2. **chebyshev.cpp:175 — Zero-size matrix allocation (nct_red)**
   - **Location:** Line 175: `std::vector<double> N_red((size_t)nct_red * (size_t)nct_red)`
   - **Problem:** If all margins have `cat_counts[k] < 2`, then `nct_red = 0`, allocating a 0×0 matrix. This is valid but wasteful. If code later assumes `nct_red > 0`, it's a logic error.
   - **Severity:** High
   - **Fix:** Add guard after line 85: `if (nct_red == 0) { /* handle single-cat-only case */ }`
   - **Confidence:** 90 (logic inference)

3. **sinkhorn.cpp:204 — Redundant empty check**
   - **Location:** Line 204: `if (std::isfinite(best_metric_seen) && !W_best.empty())`
   - **Problem:** `W_best` allocated unconditionally at line 84 with `ct.M_cell` elements. The `.empty()` check at line 204 is always false — dead condition.
   - **Severity:** High
   - **Fix:** Remove `&& !W_best.empty()` from the condition; rely only on `std::isfinite(best_metric_seen)`
   - **Confidence:** 95 (code trace)

### MEDIUM SEVERITY

4. **greg.cpp:40 — Unbounded integer addition (n_cats_total)**
   - **Location:** Line 40: `for (int k = 0; k < st.K; k++) { ... n_cats_total += st.cat_counts[k]; }`
   - **Problem:** Accumulation `n_cats_total += st.cat_counts[k]` can overflow if sum > INT_MAX (though unlikely with K ≤ 64). The guard at line 42 catches this, but no underflow check.
   - **Severity:** Medium
   - **Fix:** Change to `size_t n_cats_total = 0;` or add explicit overflow check before accumulation.
   - **Confidence:** 70 (defensive coding pattern)

5. **chebyshev.cpp:77 & 71 — Precomputable division**
   - **Location:** Lines 71, 77: `Tgt[...] = st.targets[...] * n_d` and `T_flat[m] = Tgt[m] / n_d`
   - **Problem:** `T_flat[m] = Tgt[m] / n_d` is computed once at init, but could simplify to `T_flat[m] = st.targets[k][j]` directly, avoiding two operations per element.
   - **Severity:** Medium (efficiency, not correctness)
   - **Fix:** Inline: `T_flat[cat_offset[k] + j] = st.targets[k][j]` (no multiplication/division)
   - **Confidence:** 85 (code simplification opportunity)

6. **greg.cpp:57 — Potential N matrix overflow (n_cats_total²)**
   - **Location:** Line 57: `std::vector<double> N(static_cast<size_t>(n_cats_total) * static_cast<size_t>(n_cats_total))`
   - **Problem:** If `n_cats_total = 2048` (limit), `N` size = 4,194,304 doubles = 32MB. Allocation itself is safe due to guard at line 42, but no runtime check on allocation success.
   - **Severity:** Medium
   - **Fix:** Add check: `if (N.capacity() < size_cast) { res.status = RK_ERR_BADARG; return res; }` or trust allocator exceptions.
   - **Confidence:** 60 (defensive pattern, not critical)

7. **cell_table.cpp:51 — Early return suppresses error message**
   - **Location:** Line 51: `if (K <= 0 || n <= 0) return -1;`
   - **Problem:** No error message set. Caller sees `-1` but cannot diagnose K/n issue. Compare to line 90+ which sets `res.message`.
   - **Severity:** Medium (usability)
   - **Fix:** Caller should document that `-1 → bad K or n`. Or set a static error context.
   - **Confidence:** 75 (design pattern issue)

### LOW SEVERITY

8. **grake.cpp:1 — Dead code / incomplete implementation**
   - **Location:** File `grake.cpp` contains only `#include "grake.hpp"` (1 line)
   - **Problem:** No function definitions. If `grake_solve` exists in header, this file is a stub. If no implementation is intended, should be deleted or marked `// TODO`.
   - **Severity:** Low
   - **Fix:** Add implementation or document why file exists as placeholder.
   - **Confidence:** 95 (literal file contents)

9. **sinkhorn.cpp:98 — Unused bucket_tmp allocation**
   - **Location:** Line 98: `std::vector<double> bucket(max_cats);`
   - **Problem:** Pre-allocated at init time. If iterations are low (<= 2), reallocation inside loop would be cheaper. But pre-allocation is safer; this is a minor inefficiency.
   - **Severity:** Low (performance micro-optimization, not a bug)
   - **Fix:** Profile before optimizing; current approach is acceptable.
   - **Confidence:** 50 (depends on workload)

10. **chebyshev.cpp:23 — Default status RK_ERR_NOCONV**
    - **Location:** Line 23: `res.status = RK_ERR_NOCONV;`
    - **Problem:** Default status is "did not converge" rather than "success" or "unknown". If early exit (e.g., build_cell_table fails at line 30), status remains RK_ERR_NOCONV even though it's a fatal error, not a convergence failure.
    - **Severity:** Low (semantic correctness, not functional)
    - **Fix:** Use `RK_ERR_BADARG` for initialization, or set status immediately on entry errors.
    - **Confidence:** 75 (code review pattern)

---

## Summary by Category

| Category | Count | Files |
|----------|-------|-------|
| **Division / Arithmetic** | 2 | chebyshev.cpp, greg.cpp |
| **Memory Allocation** | 4 | chebyshev.cpp, sinkhorn.cpp, greg.cpp |
| **Edge Cases (K=0, n=0)** | 3 | chebyshev.cpp, cell_table.cpp, greg.cpp |
| **Dead Code** | 2 | sinkhorn.cpp, grake.cpp |
| **Status Codes** | 1 | chebyshev.cpp |
| **Integer Overflow** | 1 | greg.cpp |

## Files Status

| File | Issues | Severity |
|------|--------|----------|
| **chebyshev.cpp** | 4 | 1 High, 2 Med, 1 Low |
| **sinkhorn.cpp** | 2 | 1 High, 1 Low |
| **greg.cpp** | 2 | 1 Med, 1 Low |
| **cell_table.cpp** | 1 | 1 Med |
| **grake.cpp** | 1 | 1 Low |

---

## Notes

- **RK_ERR_BUDGET / RK_ERR_STALL**: Not used in these 5 files; only `raking.cpp` uses them.
- **Positive**: Good pre-allocation of work vectors (bucket, bucket_b, etc.) outside main loops.
- **Positive**: Proper guards on division by X_init (lines 215 sinkhorn, 261 sinkhorn) with fallback to 1.0.
- **Positive**: Overflow checks on n_cats_total (greg line 42, chebyshev line 51).
