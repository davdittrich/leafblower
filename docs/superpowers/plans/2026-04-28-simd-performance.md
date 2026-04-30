# SIMD Performance — libmvec Vectorization — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development

**Goal:** Vectorize all hot scalar log/exp loops over M_cell using glibc libmvec (_ZGVdN4v_log, _ZGVdN4v_exp) for 3-4x speedup in math-dominated paths.

**Architecture:** S0 (lbw_math.hpp additions) is prerequisite for S1, S2, S6. S3 adds bulk_exp_clipped to lbw_math.hpp (also a prerequisite). S4 uses existing bulk_scaled_exp. S5 is trivial. All S1-S6 can run after S0.

**Tech Stack:** C++17, AVX2 intrinsics, glibc libmvec _ZGVdN4v_exp/_ZGVdN4v_log, -fopenmp-simd, bench::mark()

---

**Mechanism:** Explicit _ZGVdN4v_log declaration + bulk_log helper; precompute exp_a[] before bisection; vectorize all M_cell-scale exp/log loops
**Forbidden:** Using -ffast-math; vectorizing loops with indirect addressing (gather patterns); removing scalar fallbacks
**Audit:** bench::mark() before/after on stepstone small fixture; devtools::test() must remain identical

---

## Source Ground Truth (read before coding)

Key facts verified from source:

- `src/lbw_math.hpp` already declares `extern "C" __m256d _ZGVdN4v_exp(__m256d)` and implements `bulk_scaled_exp`. No log equivalent exists.
- `src/lbw_config.h`: `LBW_HAS_GLIBC_MVEC 1`, `LBW_HAS_OMP_SIMD 1`, `LBW_HAS_OMP 0`.
- `lbfgsb_solver.cpp` already `#include "lbw_math.hpp"` (line 10). Wolfe logit loops already have `#pragma omp simd` at lines 382, 396, 507, 521. The `std::log` calls inside those loops (lines 402, 527) are the vectorization targets.
- `raking.cpp` `compute_weight_kl` lambda (lines 163-171): `wkl += X[c] * std::log(X[c] / X_init[c]) * inv_n` — no guard beyond `X_init[c] > 0.0 && X[c] > 0.0`.
- `ieppa.cpp` `apply_single_margin_log` (lines 520-560): `std::log(W[c])` called per cell inside a per-category inner loop with indirect `cells[r]` indexing — gather pattern, cannot vectorize.
- `ieppa.cpp` X_tilde rebuild: three loop sites at lines ~731-735, ~754-758, ~833-840. Pattern: `double s = log_X_init[c] + cell_lf[c]; double s_clip = (s > kLogClip) ? kLogClip : s; X_tilde[c] = std::exp(s_clip);`. All are direct (no gather). The ~833 site also tracks `max_log_X_tilde` and `overflow_detected` — these scalar reductions cannot move into bulk_exp_clipped.
- `ieppa.cpp` T1.B f_lin update (line 670): `f_lin[cat_offset[k] + j] = std::exp(lf[cat_offset[k] + j])` inside a double loop over K × cat_counts[k] — total_cats entries, not M_cell. Use `bulk_scaled_exp(1.0, lf.data(), f_lin.data(), total_cats)` only if the index layout is contiguous (it is: both arrays sized total_cats, same offset).
- `sinkhorn.cpp` `bisect_capacity` (lines 14-45): inner eval lambda `f(mu)` computes `X[c] * std::exp(a[c] + mu)` for all c in a tight loop called ~80 times. Pre-computing `exp_a[c] = exp(a[c])` before bisection reduces this to `exp_a[c] * std::exp(mu)` — one scalar exp per bisect call instead of M_cell vector exps.
- `ieppa.cpp` homotopy outer loop begins after `X_cur` init (~line 253). Pre-alloc scratch buffers (`s_buf`, `log_W`, `kl_scratch`, `exp_a`) must go before the homotopy `for` loop, not inside it.

---

## Task 0 (Prerequisite): S0 — Add bulk_log + bulk_scaled_log to lbw_math.hpp

**File:** `src/lbw_math.hpp`

**Ticket:** one ticket — "S0: add bulk_log/bulk_scaled_log to lbw_math.hpp"

**Before** (end of file, after `bulk_scaled_exp` closing brace, before `} // namespace lbw`):
```cpp
} // namespace lbw
```

**After** (note: `extern "C"` is at file scope — OUTSIDE `namespace lbw`, matching the existing `_ZGVdN4v_exp` pattern on line 7 of lbw_math.hpp):
```cpp
// Place this block BEFORE namespace lbw { — at file scope, inside the same
// #if LBW_HAS_GLIBC_MVEC guard that wraps the existing _ZGVdN4v_exp declaration.
// The existing file structure is:
//   #if LBW_HAS_GLIBC_MVEC
//     #include <immintrin.h>
//     extern "C" __m256d _ZGVdN4v_exp(__m256d) ...;   ← line 7, file scope
//   #endif
//   namespace lbw { ... }
//
// Append _ZGVdN4v_log to the SAME #if block, after _ZGVdN4v_exp:
extern "C" __m256d _ZGVdN4v_log(__m256d) __attribute__((target("avx2")));
// Then close #endif, then namespace lbw continues.
```

**Concrete edit to lbw_math.hpp** — replace the existing `#if` block at lines 5-8:

**Before:**
```cpp
#if LBW_HAS_GLIBC_MVEC
#  include <immintrin.h>
extern "C" __m256d _ZGVdN4v_exp(__m256d) __attribute__((target("avx2")));
#endif
```

**After:**
```cpp
#if LBW_HAS_GLIBC_MVEC
#  include <immintrin.h>
extern "C" __m256d _ZGVdN4v_exp(__m256d) __attribute__((target("avx2")));
extern "C" __m256d _ZGVdN4v_log(__m256d) __attribute__((target("avx2")));
#endif
```

Then inside `namespace lbw`, before `} // namespace lbw`, add:

```cpp
// bulk_log: out[i] = log(in[i]) for i in [0, n).
// in[i] must be > 0 (caller's responsibility — mirrors std::log contract).
inline void bulk_log(const double* __restrict__ in,
                     double*       __restrict__ out,
                     int n) {
#if LBW_HAS_GLIBC_MVEC
    int i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256d v = _mm256_loadu_pd(in + i);
        _mm256_storeu_pd(out + i, _ZGVdN4v_log(v));
    }
    for (; i < n; ++i) out[i] = std::log(in[i]);
#else
    for (int i = 0; i < n; ++i) out[i] = std::log(in[i]);
#endif
}

// bulk_scaled_log: out[i] = log(scale * in[i]) for i in [0, n).
// scale > 0, in[i] > 0 required.
inline void bulk_scaled_log(double scale,
                             const double* __restrict__ in,
                             double*       __restrict__ out,
                             int n) {
#if LBW_HAS_GLIBC_MVEC
    const __m256d vs = _mm256_set1_pd(scale);
    int i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256d v = _mm256_mul_pd(_mm256_loadu_pd(in + i), vs);
        _mm256_storeu_pd(out + i, _ZGVdN4v_log(v));
    }
    for (; i < n; ++i) out[i] = std::log(scale * in[i]);
#else
    for (int i = 0; i < n; ++i) out[i] = std::log(scale * in[i]);
#endif
}

} // namespace lbw
```

**Test** (`tests/testthat/test-simd-math.R`, new file):

NOTE: `lbw_has_mvec()` and `.force_scalar = TRUE` do NOT exist in the codebase (confirmed by grep of `src/r_bridge.cpp` and `src/c_api.cpp` — no SIMD capability query is exposed to R). `LBW_HAS_GLIBC_MVEC` is a compile-time constant only. Tests verify correctness by running `harvest()` and comparing against a pre-saved scalar reference RDS generated from a scalar-only build. The SIMD path is always active when present.

```r
test_that("ieppa with bulk_log produces weights matching scalar reference within 1e-12", {
  df <- arrow::read_parquet(
    testthat::test_path("fixtures/stepstone_small.parquet"))
  result <- harvest(df, method = "ieppa")
  ref    <- readRDS(testthat::test_path("task0_ieppa_scalar_ref.rds"))
  expect_equal(result$weights, ref$weights, tolerance = 1e-12)
})
```

The scalar reference RDS (`task0_ieppa_scalar_ref.rds`) must be generated from the unmodified codebase (before this PR) and committed to `tests/testthat/`.

**Benchmark** (R console, wall-clock comparison before/after the code change — no `.force_scalar` parameter):
```r
# Run BEFORE applying bulk_log patch; save timing:
library(bench)
df <- arrow::read_parquet("tests/testthat/fixtures/stepstone_small.parquet")
# Before-patch timing (record result):
t_before <- bench::mark(harvest(df, method = "ieppa"), min_iterations = 10)

# After applying patch and R CMD INSTALL --preclean .:
t_after <- bench::mark(harvest(df, method = "ieppa"), min_iterations = 10)

# Compare:
print(t_before$median)
print(t_after$median)
# Expected: t_after$median < t_before$median (ideally >= 2x on M_cell >= 1000).
```

**Compile gate:** `R CMD INSTALL --preclean . 2>&1 | grep -E "error:|warning:"` — must be clean.

---

## Task 1 (Prerequisite): S3-prep — Add bulk_exp_clipped to lbw_math.hpp

**File:** `src/lbw_math.hpp` — same commit as Task 0.

**Ticket:** one ticket — "S3-prep: add bulk_exp_clipped to lbw_math.hpp"

**Why separate ticket:** Task 3 (S3) depends on it; keep changes auditable.

**Insert** before `} // namespace lbw` (after bulk_scaled_log):

```cpp
// bulk_exp_clipped: out[i] = exp(clamp(in[i], -clip, clip)) for i in [0, n).
// clip must be <= 700.0 (exp(700) < DBL_MAX).
inline void bulk_exp_clipped(const double* __restrict__ in,
                              double*       __restrict__ out,
                              int n, double clip) {
#if LBW_HAS_GLIBC_MVEC
    const __m256d vclip  = _mm256_set1_pd(clip);
    const __m256d vnclip = _mm256_set1_pd(-clip);
    int i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256d v = _mm256_loadu_pd(in + i);
        v = _mm256_min_pd(_mm256_max_pd(v, vnclip), vclip);
        _mm256_storeu_pd(out + i, _ZGVdN4v_exp(v));
    }
    for (; i < n; ++i) out[i] = std::exp(std::clamp(in[i], -clip, clip));
#else
    for (int i = 0; i < n; ++i) out[i] = std::exp(std::clamp(in[i], -clip, clip));
#endif
}
```

**Test** (append to `tests/testthat/test-simd-math.R`):

NOTE: No `lbw_has_mvec()` function exists; `LBW_HAS_GLIBC_MVEC` is build-time only. Test compares against pre-saved reference RDS from the unmodified build.

```r
test_that("ieppa with bulk_exp_clipped X_tilde rebuild matches pre-patch reference within 1e-12", {
  df  <- arrow::read_parquet(
    testthat::test_path("fixtures/stepstone_small.parquet"))
  ref <- readRDS(testthat::test_path("task1_ref.rds"))
  result <- harvest(df, method = "ieppa")
  expect_equal(result$weights, ref$weights, tolerance = 1e-12)
})
```

---

## Task 2: S1 — raking.cpp compute_weight_kl via bulk_scaled_log

**File:** `src/raking.cpp`

**Ticket:** one ticket — "S1: vectorize compute_weight_kl in raking.cpp with bulk_log"

**Prerequisite:** Task 0 merged.

**Target site:** `compute_weight_kl` lambda, lines 163-171. Current pattern:
```cpp
auto compute_weight_kl = [&]() -> double {
    double wkl = 0.0;
    const double inv_n = 1.0 / static_cast<double>(st.n);
    for (int c = 0; c < ct.M_cell; c++) {
        if (X_init[c] > 0.0 && X[c] > 0.0)
            wkl += X[c] * std::log(X[c] / X_init[c]) * inv_n;
    }
    return std::isfinite(wkl) ? wkl : 0.0;
};
```

**Problem:** The guarded loop prevents direct bulk_scaled_log because of the `X_init[c] > 0.0 && X[c] > 0.0` branches. We pre-fill a scratch buffer with valid ratios, then log them in bulk.

**Step 1:** Add `#include "lbw_math.hpp"` at the top of `raking.cpp` (if not already present — grep first).

**Step 2:** Declare scratch before the `while`-loop (after `W_best` declaration, before `compute_weight_kl` lambda):
```cpp
std::vector<double> kl_ratio_scratch(ct.M_cell, 0.0);
std::vector<double> kl_log_scratch(ct.M_cell, 0.0);
```

**Step 3:** Replace compute_weight_kl lambda body:

**Before:**
```cpp
auto compute_weight_kl = [&]() -> double {
    double wkl = 0.0;
    const double inv_n = 1.0 / static_cast<double>(st.n);
    for (int c = 0; c < ct.M_cell; c++) {
        if (X_init[c] > 0.0 && X[c] > 0.0)
            wkl += X[c] * std::log(X[c] / X_init[c]) * inv_n;
    }
    return std::isfinite(wkl) ? wkl : 0.0;
};
```

**After:**
```cpp
auto compute_weight_kl = [&]() -> double {
    const double inv_n = 1.0 / static_cast<double>(st.n);
    // Fill ratio buffer; zero out invalid cells.
    int valid_count = 0;
    for (int c = 0; c < ct.M_cell; c++) {
        if (X_init[c] > 0.0 && X[c] > 0.0) {
            kl_ratio_scratch[valid_count] = X[c] / X_init[c];
            kl_log_scratch[valid_count]   = X[c];  // weight for summation
            valid_count++;
        }
    }
    // Vectorized log over valid ratios only (contiguous, no gather).
    lbw::bulk_log(kl_ratio_scratch.data(), kl_ratio_scratch.data(), valid_count);
    double wkl = 0.0;
    for (int i = 0; i < valid_count; i++)
        wkl += kl_log_scratch[i] * kl_ratio_scratch[i] * inv_n;
    return std::isfinite(wkl) ? wkl : 0.0;
};
```

**Correctness note:** `bulk_log` with in==out (aliased) is safe: the AVX2 loop reads a 4-wide vector, calls `_ZGVdN4v_log`, writes back — no overlap issue because each 4-lane chunk is consumed before the next.

**Verification:** All existing `devtools::test()` raking tests must pass. Explicitly check:
```r
testthat::test_file("tests/testthat/test-raking.R")
```
KL values must match to `tolerance = 1e-12` vs baseline `task1_ref.rds`.

**Benchmark** (wall-clock before/after the code change — no `.force_scalar` parameter exists):
```r
df <- arrow::read_parquet("tests/testthat/fixtures/stepstone_small.parquet")
# Run before applying patch, record; run again after R CMD INSTALL --preclean .:
bench::mark(harvest(df, method = "raking"), min_iterations = 10)
```

---

## Task 3: S2 — ieppa.cpp log(W[c]) precomputation with bulk_log

**File:** `src/ieppa.cpp`

**Ticket:** one ticket — "S2: precompute log_W[] with bulk_log in ieppa apply_single_margin_log"

**Prerequisite:** Task 0 merged.

**Analysis:** `apply_single_margin_log` (lines 520-560) calls `std::log(W[c])` inside a category loop with indirect cell index `cells[r]` — a gather pattern. We cannot vectorize the inner loop. Instead, precompute `log_W[c]` for all c before the homotopy loop and update it after each W[] update in the capacity block.

**Step 1:** Add scratch declaration after `std::vector<double> W(ct.M_cell, 1.0);` (~line 137):
```cpp
std::vector<double> log_W(ct.M_cell, 0.0);  // log_W[c] = log(W[c]); W init=1 → log=0
```

**Step 2:** After the capacity block's W update in the log-path (lines ~861-869, the `W[c] = xc / X_tilde[c]` block):

**Before:**
```cpp
            for (int c = 0; c < ct.M_cell; c++) {
                double xc = std::clamp(X_tilde[c], L_cell[c], U_cell[c]);
                X[c] = xc;
                if (X_tilde[c] > 0.0) {
                    W[c] = xc / X_tilde[c];
                } else {
                    W[c] = 1.0;
                }
                if (W[c] != 1.0) n_cap++;
            }
```

**After:**
```cpp
            for (int c = 0; c < ct.M_cell; c++) {
                double xc = std::clamp(X_tilde[c], L_cell[c], U_cell[c]);
                X[c] = xc;
                if (X_tilde[c] > 0.0) {
                    W[c] = xc / X_tilde[c];
                } else {
                    W[c] = 1.0;
                }
                if (W[c] != 1.0) n_cap++;
            }
            // Refresh log_W after W[] update (log-path only; used by apply_single_margin_log).
            lbw::bulk_log(W.data(), log_W.data(), ct.M_cell);
```

**Step 3:** In `apply_single_margin_log`, replace `std::log(W[c])`:

**Before:**
```cpp
                    double s = log_X_init[c] + std::log(W[c]);
```

**After:**
```cpp
                    double s = log_X_init[c] + log_W[c];
```

**W write-site audit (all sites must be followed by log_W refresh):**

Grep results from `src/ieppa.cpp` for `W[c]`, `W.assign`, `W.resize`, `std::fill.*W`, `W =`:

| Line | Site | Action |
|------|------|--------|
| ~242 | `X_cur[c] = X_init[c] * W[c]` | Read-only (W not written) — no refresh needed |
| ~537 | `std::log(W[c])` in apply_single_margin_log | Read-only — replaced by `log_W[c]` in Step 3 |
| ~700 | `std::fill(W.begin(), W.end(), 1.0)` — linear→log fallback reset | Add `std::fill(log_W.begin(), log_W.end(), 0.0)` immediately after |
| ~781,~792 | `W[c] = 1.0` inside linear-path per-cell guard (X_init<=0 or X_tilde<=0) | Handled by the bulk_log refresh after the entire per-cell loop (Step 2b below) |
| ~800 | `W[c] = z / X_tilde_c` (ADMM branch) | Handled by bulk_log refresh after linear-path loop |
| ~803 | `W[c] = xc / X_tilde_c` (non-ADMM branch) | Handled by bulk_log refresh after linear-path loop |
| ~812 | `std::fill(W.begin(), W.end(), 1.0)` — overflow-detected full reset | Add `std::fill(log_W.begin(), log_W.end(), 0.0)` immediately after |
| ~863,~865 | `W[c] = xc / X_tilde[c]` or `W[c] = 1.0` — log-path capacity block | Add `lbw::bulk_log(W.data(), log_W.data(), ct.M_cell)` after the loop (Step 2, already planned) |

**Step 2a — linear-path W update (new):** After the linear-path per-cell loop (after `res.n_cap_active = n_cap;`, before the `if (overflow_detected)` check), add:
```cpp
        // Refresh log_W after linear-path W[] update (log_W not used in linear path
        // but must be current if fallback to log-space occurs mid-homotopy).
        if (!overflow_detected)
            lbw::bulk_log(W.data(), log_W.data(), ct.M_cell);
```

**Step 2b — overflow-detected full reset (~line 812):** Add immediately after `std::fill(W.begin(), W.end(), 1.0)`:
```cpp
                std::fill(log_W.begin(), log_W.end(), 0.0);
```

**Step 2c — linear→log soft reset (~line 700):** Add immediately after `std::fill(W.begin(), W.end(), 1.0)` in the soft-reset block:
```cpp
            std::fill(log_W.begin(), log_W.end(), 0.0);
```

**Edge case:** At solver entry W[c]=1 so log_W[c]=0 — correct (initialized at Step 1). All three reset paths now zero log_W in sync with W reset.

**Verification:**
```r
testthat::test_file("tests/testthat/test-ieppa.R")
testthat::test_file("tests/testthat/test-ieppa-faithful.R")
```
Weights must match `task2_ieppa_ref.rds` to `tolerance = 1e-12`.

---

## Task 4: S3 — ieppa.cpp X_tilde rebuild via bulk_exp_clipped

**File:** `src/ieppa.cpp`

**Ticket:** one ticket — "S3: vectorize X_tilde rebuild loops with bulk_exp_clipped"

**Prerequisite:** Task 1 (bulk_exp_clipped in lbw_math.hpp) merged.

**Target sites:** Three X_tilde rebuild loops. Two are simple (greedy ~lines 731-735, greedy refresh ~lines 754-758). The third (~lines 833-840) also tracks `max_log_X_tilde` and `overflow_detected` — that scalar logic must remain scalar; only the exp call is vectorized.

**Step 1:** Add `s_buf` scratch before homotopy loop, after the `log_W` declaration added in Task 3:
```cpp
std::vector<double> s_buf(ct.M_cell, 0.0);  // staging for log_X_init[c] + cell_lf[c]
```

**Site A — greedy initial rebuild (~lines 728-735):**

**Before:**
```cpp
                if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);
                for (int c = 0; c < ct.M_cell; c++) {
                    if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
                    // T2.A: single-stream exp via cell_lf (was K=20 DRAM streams)
                    double s = log_X_init[c] + cell_lf[c];
                    double s_clip = (s > kLogClip) ? kLogClip : s;
                    X_tilde[c] = std::exp(s_clip);
                }
```

**After:**
```cpp
                if (X_tilde.empty()) X_tilde.assign(ct.M_cell, 0.0);
                for (int c = 0; c < ct.M_cell; c++) {
                    if (X_init[c] <= 0.0) { s_buf[c] = -kLogClip; continue; }
                    s_buf[c] = log_X_init[c] + cell_lf[c];
                }
                lbw::bulk_exp_clipped(s_buf.data(), X_tilde.data(), ct.M_cell, kLogClip);
                // Zero-out cells with X_init==0 (bulk_exp_clipped produced exp(-kLogClip)≈0 but be exact).
                for (int c = 0; c < ct.M_cell; c++) {
                    if (X_init[c] <= 0.0) X_tilde[c] = 0.0;
                }
```

**Site B — greedy refresh (~lines 752-758):** Apply identical transformation.

**Site C — main log-path rebuild (~lines 833-840):** The overflow tracking loop cannot be fully replaced, but the exp call can be:

**Before:**
```cpp
            for (int c = 0; c < ct.M_cell; c++) {
                if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
                double s = log_X_init[c] + cell_lf[c];
                if (s > max_log_X_tilde) max_log_X_tilde = s;
                double s_clip = (s > kLogClip) ? kLogClip : s;
                if (s > kLogClip && U_cell[c] >= 1e299) {
                    overflow_detected = true;
                }
                X_tilde[c] = std::exp(s_clip);
            }
```

**After:**
```cpp
            // Pass 1: scalar — compute s values, track overflow, fill s_buf.
            for (int c = 0; c < ct.M_cell; c++) {
                if (X_init[c] <= 0.0) { s_buf[c] = -kLogClip; continue; }
                double s = log_X_init[c] + cell_lf[c];
                if (s > max_log_X_tilde) max_log_X_tilde = s;
                if (s > kLogClip && U_cell[c] >= 1e299) overflow_detected = true;
                s_buf[c] = s;
            }
            // Pass 2: vectorized exp with clip.
            lbw::bulk_exp_clipped(s_buf.data(), X_tilde.data(), ct.M_cell, kLogClip);
            for (int c = 0; c < ct.M_cell; c++) {
                if (X_init[c] <= 0.0) X_tilde[c] = 0.0;
            }
```

**Verification:**
```r
testthat::test_file("tests/testthat/test-ieppa.R")
testthat::test_file("tests/testthat/test-homotopy.R")
testthat::test_file("tests/testthat/test-ieppa-bounds-mode.R")
```

**Benchmark** (wall-clock before/after; no `.force_scalar` parameter exists):
```r
df <- arrow::read_parquet("tests/testthat/fixtures/stepstone_small.parquet")
# Run before applying Task 4 patch, record; run again after R CMD INSTALL --preclean .:
bench::mark(harvest(df, method = "ieppa"), min_iterations = 10)
```

---

## Task 5: S4 — sinkhorn.cpp bisect_capacity exp_a precompute

**File:** `src/sinkhorn.cpp`

**Ticket:** one ticket — "S4: precompute exp_a[] before bisect_capacity to eliminate M_cell exp per bisect eval"

**No new helpers needed** — `bulk_scaled_exp(1.0, ...)` already exists.

**Add include** if not present: `#include "lbw_math.hpp"` at top of sinkhorn.cpp.

**Analysis:** `bisect_capacity` is called inside the main sinkhorn sweep loop. Its `f(mu)` lambda computes `X[c] * std::exp(a[c] + mu)` for all M_cell cells, called ~80 times per bisection. Precomputing `exp_a[c] = exp(a[c])` reduces each eval to `exp_a[c] * std::exp(mu)` — one scalar multiply+exp per cell instead of one vector exp.

**Step 1:** Declare `exp_a` scratch before the main sinkhorn sweep loop (after `std::vector<double> a(ct.M_cell, 0.0);`):
```cpp
std::vector<double> exp_a(ct.M_cell, 1.0);  // exp_a[c] = exp(a[c]); a init=0 → exp=1
```

**Step 2:** After each `a[c]` update (after the `X_proj` acceptance block where `a[c] += std::log(X[c]) - std::log(X_proj[c])`):
```cpp
        lbw::bulk_scaled_exp(1.0, a.data(), exp_a.data(), ct.M_cell);
```

**Step 3:** Add a `bisect_capacity_fast` overload (or modify signature) that accepts `exp_a` and replaces `std::exp(a[c] + mu)` with `exp_a[c] * exp_mu`:

The cleanest approach: add an overload `bisect_capacity` taking `const double* exp_a_data`:

```cpp
static bool bisect_capacity(const std::vector<double>& X,
                             const double* exp_a_data,   // exp(a[c]), pre-computed
                             const std::vector<double>& L,
                             const std::vector<double>& U,
                             int M_cell,
                             double target_mass,
                             double& mu_out,
                             std::vector<double>& X_proj)
{
    double sum_L = 0.0, sum_U = 0.0;
    for (int c = 0; c < M_cell; c++) { sum_L += L[c]; sum_U += U[c]; }
    if (sum_L > target_mass + 1e-9 || sum_U < target_mass - 1e-9) return false;

    auto f = [&](double mu) -> double {
        const double exp_mu = std::exp(mu);
        double s = 0.0;
        for (int c = 0; c < M_cell; c++)
            s += std::clamp(X[c] * exp_a_data[c] * exp_mu, L[c], U[c]);
        return s - target_mass;
    };

    double lo = -50.0, hi = 50.0;
    while (lo > -500.0 && f(lo) > 0.0) lo *= 2.0;
    while (hi < 500.0  && f(hi) < 0.0) hi *= 2.0;
    if (f(lo) > 0.0 || f(hi) < 0.0) return false;
    for (int i = 0; i < 80; i++) {
        double mid = 0.5 * (lo + hi);
        if (f(mid) < 0.0) lo = mid; else hi = mid;
        if (hi - lo < 1e-12) break;
    }
    mu_out = 0.5 * (lo + hi);
    const double exp_mu_out = std::exp(mu_out);
    for (int c = 0; c < M_cell; c++)
        X_proj[c] = std::clamp(X[c] * exp_a_data[c] * exp_mu_out, L[c], U[c]);
    return true;
}
```

**Call site count:** `grep -n "bisect_capacity" src/sinkhorn.cpp` yields exactly 2 hits: the function definition at line 14 and one call site at line 143. Only one call site to update.

Update call site (~line 143) to pass `exp_a.data()`.

**Verification:**
```r
testthat::test_file("tests/testthat/test-calibration-solvers.R")
```

---

## Task 6: S5 — ieppa.cpp T1.B f_lin update via bulk_scaled_exp

**File:** `src/ieppa.cpp`

**Ticket:** one ticket — "S5: vectorize T1.B f_lin = exp(lf) with bulk_scaled_exp"

**Same commit as Task 3** (both in ieppa.cpp — batch to minimize recompile cycles).

**Target site:** T1.B correction block (~line 668-671):
```cpp
                for (int k = 0; k < st.K; k++) {
                    for (int j = 0; j <= st.cat_counts[k]; j++) {
                        lf[cat_offset[k] + j] += lf_correction;
                        f_lin[cat_offset[k] + j] = std::exp(lf[cat_offset[k] + j]);
                    }
                }
```

**Step 1:** The double loop over k × (cat_counts[k]+1) touches all `total_cats` entries sequentially (by construction of cat_offset). But it also updates `lf[]` in-place before reading it for exp. Split into two passes:

**After:**
```cpp
                for (int k = 0; k < st.K; k++) {
                    for (int j = 0; j <= st.cat_counts[k]; j++) {
                        lf[cat_offset[k] + j] += lf_correction;
                    }
                }
                lbw::bulk_scaled_exp(1.0, lf.data(), f_lin.data(), total_cats);
```

**Why valid:** `cat_offset` is defined such that `cat_offset[st.K] == total_cats`, and `lf` and `f_lin` are both allocated to `total_cats`. The NA bucket per margin occupies the `j == cat_counts[k]` slot — the existing loop already writes it, and the bulk call covers all total_cats entries including NA buckets, which is identical to the original.

**Verification:**
```r
testthat::test_file("tests/testthat/test-ieppa.R")
testthat::test_file("tests/testthat/test-ieppa-faithful.R")
```

---

## Task 7: S6 — lbfgsb_solver.cpp Wolfe logit log vectorization

**File:** `src/lbfgsb_solver.cpp`

**Ticket:** one ticket — "S6: vectorize Wolfe logit std::log calls in lbfgsb_solver"

**Prerequisite:** Task 0 merged. `lbw_math.hpp` is already included (line 10 confirmed).

**Target sites:** Two structurally identical logit loops at lines ~396-411 and ~521-536. Each computes:
```cpp
double Hi = L * u_work[i] + R * std::log(denom / UmL);
```
inside an `#pragma omp simd` loop. `denom = P + Q * ei` where `ei = e_vec[i]` (pre-computed exp vector). `denom / UmL` is the argument to log.

**Problem:** The loop body mixes exp-derived arithmetic with a `std::log` call. The `#pragma omp simd` is already present but auto-vectorization of `std::log` requires libmvec linkage that the compiler may not inject automatically.

**Step 1:** Compile with `-fopt-info-vec-omp-missed` flag (add to `src/Makevars` temporarily or use `PKG_CXXFLAGS`):
```
PKG_CXXFLAGS="$(PKG_CXXFLAGS) -fopt-info-vec-omp-missed"
```
Run `R CMD INSTALL --preclean .` and grep the output for the logit loop location. If the compiler reports it vectorized with `_ZGVdN4v_log`, no source change needed.

**Step 2 (if not auto-vectorized):** Pre-compute `log_denom_over_UmL[]` into a scratch buffer before the logit loop using `bulk_log`, then substitute:

Declare scratch before the Wolfe line-search function body or as a resized-on-demand member:
```cpp
std::vector<double> log_denom_scratch(st.n);
```

Replace loop body (shown for the first occurrence; apply identically to second):

**Before:**
```cpp
#pragma omp simd reduction(+:phi_acc) reduction(+:slope_acc)
                for (int i = 0; i < st.n; i++) {
                    double ei = e_vec[i];
                    double denom = P + Q * ei;
                    double Fi = (A + B * ei) / denom;
                    double Hi = L * u_work[i] + R * std::log(denom / UmL);
                    phi_acc   += d[i] * Hi;
                    slope_acc += d[i] * Fi * du[i];
                }
```

**After:**
```cpp
                // Pre-compute log(denom/UmL) via bulk_log for vectorization.
                for (int i = 0; i < st.n; i++) {
                    log_denom_scratch[i] = (P + Q * e_vec[i]) / UmL;
                }
                lbw::bulk_log(log_denom_scratch.data(), log_denom_scratch.data(), st.n);
#pragma omp simd reduction(+:phi_acc) reduction(+:slope_acc)
                for (int i = 0; i < st.n; i++) {
                    double ei = e_vec[i];
                    double denom = P + Q * ei;
                    double Fi = (A + B * ei) / denom;
                    double Hi = L * u_work[i] + R * log_denom_scratch[i];
                    phi_acc   += d[i] * Hi;
                    slope_acc += d[i] * Fi * du[i];
                }
```

**Correctness note:** `log_denom_scratch` aliasing is safe — it is fully written before the simd loop reads it.

**Verification:**
```r
testthat::test_file("tests/testthat/test-lbfgsb.R")
testthat::test_file("tests/testthat/test-logit.R")
```

**Benchmark** (wall-clock before/after; no `.force_scalar` parameter exists):
```r
df <- arrow::read_parquet("tests/testthat/fixtures/stepstone_small.parquet")
# Run before applying Task 7 patch, record; run again after R CMD INSTALL --preclean .:
bench::mark(harvest(df, method = "lbfgsb"), min_iterations = 10)
```

---

## Execution Order & Dependencies

```
Task 0 (S0: bulk_log + bulk_scaled_log)  ─┬─→ Task 2 (S1: raking KL)
Task 1 (S3-prep: bulk_exp_clipped)        │ ─→ Task 3 (S2: ieppa log_W)
  [same commit as Task 0]                 │ ─→ Task 4 (S3: X_tilde) [needs Task 1]
                                           │ ─→ Task 6 (S5: T1.B f_lin) [same commit as Task 4]
                                           │ ─→ Task 7 (S6: lbfgsb logit)
                                           └─→ Task 5 (S4: sinkhorn exp_a) [no new helper needed]
```

Tasks 2, 3, 5, 7 are independent after Task 0. Tasks 4 and 6 are independent after Task 1. Tasks 4 and 6 should land in the same commit (both ieppa.cpp).

---

## Global Verification Protocol

After all tasks are committed:

```r
# 1. Full test suite
devtools::test()

# 2. Comprehensive benchmark (stepstone_small fixture)
df <- arrow::read_parquet("tests/testthat/fixtures/stepstone_small.parquet")
bench::mark(
  raking  = harvest(df, method = "raking"),
  ieppa   = harvest(df, method = "ieppa"),
  sinkhorn = harvest(df, method = "sinkhorn"),
  lbfgsb  = harvest(df, method = "lbfgsb"),
  min_iterations = 5
)

# 3. Baseline reference comparison
ref_ieppa <- readRDS("tests/testthat/task2_ieppa_ref.rds")
new_ieppa <- harvest(df, method = "ieppa")
stopifnot(max(abs(new_ieppa$weights - ref_ieppa$weights)) < 1e-12)
```

**Pass criteria:** All tests green; no weight regression > 1e-12; ieppa log-path wall time reduced ≥ 2x on M_cell ≥ 1000.

---

## Self-Review Checklist

- [x] No `-ffast-math` introduced — all changes use explicit libmvec calls, not compiler auto-flags.
- [x] No gather-pattern loops vectorized — `apply_single_margin_log` inner loop (indirect `cells[r]`) is not touched; only scalar `std::log(W[c])` is replaced by a precomputed lookup.
- [x] Scalar fallbacks present in every helper — `#else` branch in all `#if LBW_HAS_GLIBC_MVEC` blocks.
- [x] All scratch buffers pre-allocated before homotopy loop — no per-iteration `std::vector` construction.
- [x] T1.B f_lin patch (Task 6) correctly splits add-then-exp into two passes — avoids reading `lf` before its update within the same loop.
- [x] Site C X_tilde (Task 4) retains scalar `max_log_X_tilde` and `overflow_detected` tracking — correctness-critical.
- [x] `log_W` reset to 0 on linear→log fallback — matches `W[]` reset to 1.0.
- [x] `bulk_log` with aliased in==out is safe for AVX2 non-temporal 32-byte stores — verified: `_mm256_storeu_pd` writes 32 bytes at `out+i` which are already consumed in the same iteration.
- [x] `sinkhorn.cpp` overload does not break the `a`-vector `std::clamp` guard (lines ~51, `kAmax=30.0`) — `exp_a` is updated after each `a[]` update, which happens after clamping.
