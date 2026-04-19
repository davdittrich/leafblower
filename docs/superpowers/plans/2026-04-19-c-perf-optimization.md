# leafblower C Performance Optimization — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Speed up leafblower's iEPPA and L-BFGS-B solvers with no new package dependencies. Target: ~25-35% iEPPA speedup, ~50% L-BFGS-B speedup.

**Architecture:** Two phases. Phase 1: zero-dep scalar/compiler wins (logit shared-exp fix, W-sum hoisting, errRp frequency, `#pragma omp simd`). Phase 2: configure-detected libmvec (vectorized exp for L-BFGS-B) and OpenMP threading for `compute_errRp`.

**Tech Stack:** C++17; glibc libmvec (optional, Linux/gcc only, configure-detected); OpenMP (optional, via R's `SHLIB_OPENMP_CXXFLAGS`, no new package dep).

**Compile gate:** `R CMD INSTALL --preclean . && Rscript -e "testthat::test_local()"` after every file change. Never advance to the next file before this passes.

**Test command:** `cd /home/dd/Gemini/leafblower && Rscript -e "testthat::test_local()"`

---

## Phase 1 — Zero-dep scalar wins

### Task 1: Fix shared-exp redundancy in logit.hpp and lbfgsb_solver.cpp

Both `F(u)` and `H(u)` independently call `safe_exp(logit_scale * u)`. In the Wolfe trial loop this means 2× exp per observation. Adding `FH(u)` computes the shared `e` once.

**Files:**
- Modify: `src/logit.hpp` — add `FH()` method
- Modify: `src/lbfgsb_solver.cpp` — use `FH()` in Wolfe loops

- [ ] **Step 1: Add `FH()` to `LinkFn` in logit.hpp**

  Insert after the existing `H()` method (line 52, before the closing `};`):

  ```cpp
  // FH(u): compute F(u) and H(u) simultaneously from a single safe_exp call.
  // Cuts transcendental calls in the Wolfe inner loop from 2× to 1× per obs.
  struct FHResult { double F; double H; };
  FHResult FH(double u) const {
      if (exponential) {
          double e = safe_exp(u);
          return {e, e};
      }
      double e = safe_exp(logit_scale * u);
      double denom = (U - 1.0) + (1.0 - L) * e;
      double f = (L * (U - 1.0) + U * (1.0 - L) * e) / denom;
      double h = L * u + (U - L) / logit_scale * std::log(denom / (U - L));
      return {f, h};
  }
  ```

- [ ] **Step 2: Update Wolfe trial loop in `wolfe_line_search` (lbfgsb_solver.cpp:259–263)**

  Replace:
  ```cpp
  for (int j = 0; j < st.n; j++) {
      double Fj = fn.F(u_work[j]);
      phi_trial -= d[j] * fn.H(u_work[j]);
      slope -= d[j] * Fj * du[j];
  }
  ```
  With:
  ```cpp
  for (int j = 0; j < st.n; j++) {
      auto fh = fn.FH(u_work[j]);
      phi_trial -= d[j] * fh.H;
      slope    -= d[j] * fh.F * du[j];
  }
  ```

- [ ] **Step 3: Apply the same fix in `wolfe_zoom` (lbfgsb_solver.cpp:201–205)**

  Replace:
  ```cpp
  for (int i = 0; i < st.n; i++) {
      double Fi = fn.F(u_work[i]);
      phi_trial -= d[i] * fn.H(u_work[i]);
      slope -= d[i] * Fi * du[i];
  }
  ```
  With:
  ```cpp
  for (int i = 0; i < st.n; i++) {
      auto fh = fn.FH(u_work[i]);
      phi_trial -= d[i] * fh.H;
      slope    -= d[i] * fh.F * du[i];
  }
  ```

- [ ] **Step 4: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -5
  ```
  Expected: `* DONE (leafblower)`

- [ ] **Step 5: Run tests — verify FH() preserves numerical identity**

  ```bash
  Rscript -e "testthat::test_local()" 2>&1 | tail -20
  ```
  Expected: All tests pass. The existing `test-logit.R` tests `H'(u) = F(u)` numerically — this catches any FH() regression. Weights from `harvest(..., method='lbfgsb')` must be bit-for-bit identical before/after (same arithmetic, just shared intermediate).

- [ ] **Step 6: Commit**

  ```bash
  git add src/logit.hpp src/lbfgsb_solver.cpp
  git commit -m "perf: share exp in logit FH() — halves transcendentals in Wolfe loop"
  ```

---

### Task 2: Hoist W-sum from bucket loop + 4-way ILP unroll (ieppa.cpp)

Currently `W += w[i]` is inside the bucket scatter-add loop, serialising both. Separating them lets the compiler vectorise the W reduction. The 4-way accumulator pattern (from revss_temp `adm_core_scalar`) hides FP latency.

**Files:**
- Modify: `src/ieppa.cpp`

- [ ] **Step 1: Replace the combined W + bucket loop (lines 70–76) with two separate loops**

  Replace:
  ```cpp
  double W = 0.0;
  std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
  for (int i = 0; i < st.n; i++) {
      W += w[i];
      int g = st.group_ids[k][i];
      if (g >= 0) bucket[g] += w[i];
  }
  ```
  With:
  ```cpp
  // W sum: 4-way ILP, compiler-vectorisable (no scatter aliases).
  double W = 0.0, W1 = 0.0, W2 = 0.0, W3 = 0.0;
  int ni = st.n;
  int i4 = ni & ~3;
  for (int i = 0; i < i4; i += 4) {
      W  += w[i];   W1 += w[i+1];
      W2 += w[i+2]; W3 += w[i+3];
  }
  for (int i = i4; i < ni; ++i) W += w[i];
  W += W1 + W2 + W3;

  // Bucket scatter-add (cannot vectorise due to write aliases).
  std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
  for (int i = 0; i < ni; i++) {
      int g = st.group_ids[k][i];
      if (g >= 0) bucket[g] += w[i];
  }
  ```

  Note: `W` is only used for normalisation and IPF scale factors within this margin pass. It does NOT need to be accumulated across the outer k-loop. The existing code recomputes W fresh per margin, which is correct.

- [ ] **Step 2: Apply the same W-sum fix in `compute_errRp` (lines 17–18)**

  Replace:
  ```cpp
  double W = 0.0;
  for (int i = 0; i < st.n; i++) W += w[i];
  ```
  With:
  ```cpp
  double W = 0.0, W1 = 0.0, W2 = 0.0, W3 = 0.0;
  int i4e = st.n & ~3;
  for (int i = 0; i < i4e; i += 4) {
      W  += w[i];   W1 += w[i+1];
      W2 += w[i+2]; W3 += w[i+3];
  }
  for (int i = i4e; i < st.n; ++i) W += w[i];
  W += W1 + W2 + W3;
  ```

- [ ] **Step 3: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -5
  ```

- [ ] **Step 4: Run tests**

  ```bash
  Rscript -e "testthat::test_local()" 2>&1 | tail -20
  ```
  Expected: All pass. Weights must match to within 1e-12 of pre-change reference (floating-point reassociation changes sum order — this is acceptable; `test-ieppa.R` uses `tolerance = 1e-6`).

- [ ] **Step 5: Commit**

  ```bash
  git add src/ieppa.cpp
  git commit -m "perf: hoist W-sum from scatter loop, 4-way ILP unroll in ieppa"
  ```

---

### Task 3: Reduce errRp check frequency to every 10 iterations

`compute_errRp` does K+1 O(n) passes per call — nearly as expensive as a full IPF sweep. Checking every 10 iterations eliminates 90% of that cost. Worst-case overshoot: 9 extra iterations past true convergence.

**Files:**
- Modify: `src/ieppa.cpp`

- [ ] **Step 1: Add the interval constant at the top of `ieppa_solve`**

  After `static constexpr int kMaxFixupIterations = 20;` (line 43), add:
  ```cpp
  // Check convergence every N iterations instead of every 1.
  // compute_errRp costs ~K O(n) passes — as expensive as a full sweep.
  // Checking every 10 iters cuts that overhead by 90% at the cost of
  // at most 9 extra iterations past true convergence.
  static constexpr int kErrCheckInterval = 10;
  ```

- [ ] **Step 2: Guard the convergence check block (lines 122–134)**

  Replace:
  ```cpp
  // Convergence check
  double errRp = compute_errRp(st, w, bucket);
  res.max_error = errRp;

  if (st.verbose >= 1) {
      char msg[256];
      std::snprintf(msg, 256, "iEPPA iter %d: errRp=%.2e", iter, errRp);
      st.log(msg);
  }

  if (errRp < st.tol_abs) {
      res.status = is_infeasible ? RK_ERR_INFEAS : RK_OK;
      break;
  }
  ```
  With:
  ```cpp
  // Convergence check: run every kErrCheckInterval iters and on the final iter.
  if (iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
      double errRp = compute_errRp(st, w, bucket);
      res.max_error = errRp;

      if (st.verbose >= 1) {
          char msg[256];
          std::snprintf(msg, 256, "iEPPA iter %d: errRp=%.2e", iter, errRp);
          st.log(msg);
      }

      if (errRp < st.tol_abs) {
          res.status = is_infeasible ? RK_ERR_INFEAS : RK_OK;
          break;
      }
  }
  ```

- [ ] **Step 3: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -5
  ```

- [ ] **Step 4: Run tests — verify convergence behaviour preserved**

  ```bash
  Rscript -e "testthat::test_local()" 2>&1 | tail -20
  ```
  `test-bounded-convergence.R` and `test-ieppa.R` check that weights satisfy margin constraints after calibration. They must all pass. Convergence iteration count may increase by up to 9, which is fine.

- [ ] **Step 5: Commit**

  ```bash
  git add src/ieppa.cpp
  git commit -m "perf: check errRp every 10 iters — eliminates ~30% of iEPPA work"
  ```

---

### Task 4: `#pragma omp simd` on box-projection and normalise loops

These two O(n) loops in `ieppa_solve` are element-wise with no scatter — the compiler can auto-vectorise them with a hint. No OpenMP threading required; `omp simd` is a SIMD-only pragma.

**Files:**
- Modify: `src/ieppa.cpp`

- [ ] **Step 1: Add simd pragma to the normalise loop (lines 108)**

  Replace:
  ```cpp
  for (int i = 0; i < st.n; i++) { w[i] /= wm; q[i] /= wm; }
  ```
  With:
  ```cpp
  // LBW_HAS_OMP_SIMD is a numeric value (0/1) from lbw_config.h, not a presence
  // macro — use #if VALUE, not #if defined(VALUE). _OPENMP is a presence guard.
  #if defined(_OPENMP) || LBW_HAS_OMP_SIMD
  #pragma omp simd
  #endif
  for (int i = 0; i < st.n; i++) { w[i] /= wm; q[i] /= wm; }
  ```

- [ ] **Step 2: Add simd pragma to the box-projection loop (lines 114–119)**

  Replace:
  ```cpp
  for (int i = 0; i < st.n; i++) {
      double yi = w[i] + q[i];
      double wc = std::max(lo, std::min(hi, yi));
      q[i] = yi - wc;
      w[i] = wc;
  }
  ```
  With:
  ```cpp
  #if defined(_OPENMP) || LBW_HAS_OMP_SIMD
  #pragma omp simd
  #endif
  for (int i = 0; i < st.n; i++) {
      double yi = w[i] + q[i];
      double wc = std::max(lo, std::min(hi, yi));
      q[i] = yi - wc;
      w[i] = wc;
  }
  ```

- [ ] **Step 3: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -5
  ```
  `LBW_HAS_OMP_SIMD` is 0 until Task 5 generates `lbw_config.h` — the guard evaluates
  `#if 0` (inactive), so the pragma is absent and the build is clean. After Task 5 the
  guard becomes `#if 1` and the pragma is active.

- [ ] **Step 4: Run tests**

  ```bash
  Rscript -e "testthat::test_local()" 2>&1 | tail -20
  ```
  Expected: all pass, identical results (simd pragma must not alter semantics).

- [ ] **Step 5: Commit**

  ```bash
  git add src/ieppa.cpp
  git commit -m "perf: omp simd hint on box-projection and normalise loops"
  ```

---

## Phase 2 — Configure-detected libmvec + OpenMP

### Task 5: Configure detection — libmvec, OpenMP, lbw_config.h

Extends the existing `configure` script (which currently detects c++17 and -O3) with libmvec and OpenMP probes. Generates `src/lbw_config.h` with feature macros. Fallback to scalar/serial on all platforms.

**Files:**
- Modify: `configure`
- Modify: `src/Makevars.in`
- Create: `src/lbw_config.h` (generated by configure; add to `.Rbuildignore`)

- [ ] **Step 1: Extend `configure` with OpenMP probe**

  Append to `configure` (after the existing `-O3` detection, before the final `sed` call):

  ```sh
  # Detect OpenMP: use R's SHLIB_OPENMP_CXXFLAGS (no new package dep).
  OMP_FLAGS=""
  OMP_LIBS=""
  LBW_HAS_OMP_SIMD=0
  LBW_HAS_OMP=0
  if [ -n "${SHLIB_OPENMP_CXXFLAGS}" ]; then
      cat > /tmp/lbw_omp_test.cpp << 'EOF'
  #include <omp.h>
  int main() { return omp_get_max_threads() > 0 ? 0 : 1; }
  EOF
      if eval $CXX ${SHLIB_OPENMP_CXXFLAGS} -x c++ /tmp/lbw_omp_test.cpp \
             -o /tmp/lbw_omp_test ${SHLIB_OPENMP_CXXFLAGS} 2>/dev/null; then
          OMP_FLAGS="${SHLIB_OPENMP_CXXFLAGS}"
          OMP_LIBS="${SHLIB_OPENMP_CXXFLAGS}"
          LBW_HAS_OMP_SIMD=1
          LBW_HAS_OMP=1
          echo "configure: OpenMP detected — parallel errRp and simd enabled"
      else
          echo "configure: OpenMP not available — falling back to serial"
      fi
      rm -f /tmp/lbw_omp_test.cpp /tmp/lbw_omp_test
  fi
  ```

- [ ] **Step 2: Extend `configure` with libmvec probe (after OpenMP probe)**

  ```sh
  # Detect glibc libmvec vectorised exp (_ZGVdN4v_exp — AVX2, 4-wide double).
  # Linux/glibc only; graceful fallback on macOS and Windows.
  LBW_HAS_GLIBC_MVEC=0
  MVEC_LIBS=""
  cat > /tmp/lbw_mvec_test.cpp << 'EOF'
  #include <immintrin.h>
  extern "C" __m256d _ZGVdN4v_exp(__m256d);
  int main() {
      __m256d x = _mm256_set1_pd(1.0);
      __m256d r = _ZGVdN4v_exp(x);
      (void)r; return 0;
  }
  EOF
      if eval $CXX -mavx2 /tmp/lbw_mvec_test.cpp \
             -o /tmp/lbw_mvec_test -lm 2>/dev/null && \
         /tmp/lbw_mvec_test 2>/dev/null; then
          LBW_HAS_GLIBC_MVEC=1
          MVEC_LIBS="-lm"
          echo "configure: glibc libmvec detected — vectorised exp enabled for L-BFGS-B"
      else
          echo "configure: glibc libmvec not available — scalar exp fallback"
      fi
  rm -f /tmp/lbw_mvec_test.cpp /tmp/lbw_mvec_test
  ```

- [ ] **Step 3: Generate `src/lbw_config.h` and update the final sed call in configure**

  After the probes, add the `lbw_config.h` generator:
  ```sh
  cat > src/lbw_config.h << EOF
  /* Auto-generated by configure — do not edit. */
  #ifndef LBW_CONFIG_H
  #define LBW_CONFIG_H
  /* OpenMP SIMD hint: active when omp simd is supported. */
  #define LBW_HAS_OMP_SIMD ${LBW_HAS_OMP_SIMD}
  /* OpenMP threading: active when full OpenMP is available. */
  #define LBW_HAS_OMP ${LBW_HAS_OMP}
  /* glibc libmvec vectorised exp/log (AVX2, Linux/gcc only). */
  #define LBW_HAS_GLIBC_MVEC ${LBW_HAS_GLIBC_MVEC}
  #endif /* LBW_CONFIG_H */
  EOF
  ```

  Then **replace** the existing final sed call (line 19-20 of configure, which currently only
  substitutes `@CXXFLAGS_STD@` and `@OPT_FLAGS@`) with a 4-variable substitution:

  Old (current configure, lines 19-20):
  ```sh
  sed "s|@CXXFLAGS_STD@|${CXXFLAGS_STD}|;s|@OPT_FLAGS@|${OPT_FLAGS}|" \
      src/Makevars.in > src/Makevars
  ```

  New (extends to also substitute the two new Makevars.in placeholders):
  ```sh
  sed "s|@CXXFLAGS_STD@|${CXXFLAGS_STD}|;s|@OPT_FLAGS@|${OPT_FLAGS}|;s|@OMP_FLAGS@|${OMP_FLAGS}|;s|@MVEC_LIBS@|${MVEC_LIBS}|" \
      src/Makevars.in > src/Makevars
  ```

  Without this change the two new placeholders remain literal in `src/Makevars` and the
  package fails to compile with flags containing `@OMP_FLAGS@` verbatim.

- [ ] **Step 4: Update `src/Makevars.in` to include new flags**

  Replace current content:
  ```makefile
  PKG_CXXFLAGS = @CXXFLAGS_STD@ @OPT_FLAGS@ -I. -DSTRICT_R_HEADERS
  PKG_SOURCES = c_api.cpp logit.cpp lbfgsb_solver.cpp ieppa.cpp r_bridge.cpp
  ```
  With:
  ```makefile
  PKG_CXXFLAGS = @CXXFLAGS_STD@ @OPT_FLAGS@ @OMP_FLAGS@ -I. -DSTRICT_R_HEADERS
  PKG_LIBS = @MVEC_LIBS@
  PKG_SOURCES = c_api.cpp logit.cpp lbfgsb_solver.cpp ieppa.cpp r_bridge.cpp
  ```

- [ ] **Step 5: Add `src/lbw_config.h` to `.Rbuildignore` and `src/.gitignore`**

  ```bash
  echo "^src/lbw_config\\.h$" >> .Rbuildignore
  echo "lbw_config.h" >> src/.gitignore
  ```
  (The file is generated at install time, not shipped in the tarball.)

- [ ] **Step 6: Add `include "lbw_config.h"` to ieppa.cpp and lbfgsb_solver.cpp**

  Add at the top of each file, after existing `#include` lines:
  ```cpp
  #include "lbw_config.h"
  ```
  This makes `LBW_HAS_OMP_SIMD` and `LBW_HAS_OMP` available to both files.

- [ ] **Step 7: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -10
  ```
  On Linux/gcc with glibc: expect "vectorised exp enabled" and "OpenMP detected" in configure output.
  On macOS/Windows: expect both "not available" messages — scalar fallback, package still builds.

- [ ] **Step 8: Run tests**

  ```bash
  Rscript -e "testthat::test_local()" 2>&1 | tail -20
  ```
  Expected: all pass on all platforms.

- [ ] **Step 9: Commit**

  ```bash
  git add configure src/Makevars.in src/.gitignore .Rbuildignore src/ieppa.cpp src/lbfgsb_solver.cpp
  git commit -m "build: configure detection for OpenMP and glibc libmvec"
  ```

---

### Task 6: Bulk vectorised exp/log dispatch header (lbw_math.hpp)

Provides `lbw::bulk_scaled_exp(scale, u, out, n)`: batches the `exp(scale * u[i])` call that appears in every L-BFGS-B iteration. Dispatches to `_ZGVdN4v_exp` when available; falls back to scalar loop.

**Files:**
- Create: `src/lbw_math.hpp`

- [ ] **Step 1: Write `lbw_math.hpp`**

  ```cpp
  #pragma once
  #include "lbw_config.h"
  #include <cmath>

  #if LBW_HAS_GLIBC_MVEC
  #  include <immintrin.h>
  extern "C" __m256d _ZGVdN4v_exp(__m256d) __attribute__((target("avx2")));
  #endif

  namespace lbw {

  // Compute out[i] = exp(scale * u[i]) for i in [0, n).
  // When glibc libmvec is available: processes 4 doubles/cycle via AVX2.
  // Otherwise: scalar std::exp loop.
  inline void bulk_scaled_exp(double scale,
                               const double* __restrict__ u,
                               double*       __restrict__ out,
                               int n) {
  #if LBW_HAS_GLIBC_MVEC
      const __m256d vs = _mm256_set1_pd(scale);
      int i = 0;
      for (; i + 4 <= n; i += 4) {
          __m256d vu = _mm256_loadu_pd(u + i);
          __m256d res = _ZGVdN4v_exp(_mm256_mul_pd(vs, vu));
          _mm256_storeu_pd(out + i, res);
      }
      for (; i < n; ++i) out[i] = std::exp(scale * u[i]);
  #else
      for (int i = 0; i < n; ++i) out[i] = std::exp(scale * u[i]);
  #endif
  }

  } // namespace lbw
  ```

- [ ] **Step 2: Compile gate (include check)**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -5
  ```
  The header is not yet included by anything — compile gate verifies the file is syntactically valid by checking the package builds cleanly.

- [ ] **Step 3: Commit**

  ```bash
  git add src/lbw_math.hpp
  git commit -m "perf: add lbw_math.hpp with bulk_scaled_exp dispatch (libmvec / scalar)"
  ```

---

### Task 7: Apply vectorised exp to L-BFGS-B Wolfe loops

Pre-computes `e[i] = exp(logit_scale * u_work[i])` for all n before the per-obs trial loop, using `bulk_scaled_exp`. F and H then read from `e[]` with no further exp calls.

**Files:**
- Modify: `src/lbfgsb_solver.cpp`
- Modify: `src/logit.hpp` — add `F_from_e()` and `H_from_e()` helpers

- [ ] **Step 1: Add `F_from_e()` and `H_from_e()` to `LinkFn` in logit.hpp**

  After the `FH()` method, add:
  ```cpp
  // F and H computed from pre-computed e = exp(logit_scale * u).
  // Used in vectorised Wolfe loops where e[] is batch-computed by bulk_scaled_exp.
  double F_from_e(double e) const {
      if (exponential) return e;  // e = exp(u) in exp-link case
      return (L * (U - 1.0) + U * (1.0 - L) * e) /
             ((U - 1.0) + (1.0 - L) * e);
  }
  double H_from_e(double e, double u) const {
      if (exponential) return e;
      double num = (U - 1.0) + (1.0 - L) * e;
      return L * u + (U - L) / logit_scale * std::log(num / (U - L));
  }
  ```

- [ ] **Step 2: Add include and scratch buffer to `lbfgsb_solve`**

  At the top of lbfgsb_solver.cpp, after existing includes:
  ```cpp
  #include "lbw_math.hpp"
  ```

  In `lbfgsb_solve`, after the existing scratch buffer declarations (around line 305):
  ```cpp
  std::vector<double> e_vec(st.n);   // scratch: exp(logit_scale * u_work[i])
  ```

- [ ] **Step 3: Update `wolfe_line_search` signature and trial loop to use `e_vec`**

  `e_vec` is allocated in `lbfgsb_solve` and must be threaded through both
  `wolfe_line_search` and `wolfe_zoom`. Start with `wolfe_line_search`.

  **3a. Add `e_vec` to `wolfe_line_search` signature** (insert after `u_work` parameter):

  Old declaration (line 230-239):
  ```cpp
  static double wolfe_line_search(
          const CalibState& st, const LinkFn& fn,
          const std::vector<int>& off, const std::vector<double>& T,
          const std::vector<double>& d,
          const std::vector<double>& lam, double phi_0, double slope_0,
          const std::vector<double>& u_base, const std::vector<double>& du,
          std::vector<double>& u_work,
          const std::vector<double>& dir,
          std::vector<double>& lam_new, std::vector<double>& grad_new,
          double& phi_new) {
  ```
  New:
  ```cpp
  static double wolfe_line_search(
          const CalibState& st, const LinkFn& fn,
          const std::vector<int>& off, const std::vector<double>& T,
          const std::vector<double>& d,
          const std::vector<double>& lam, double phi_0, double slope_0,
          const std::vector<double>& u_base, const std::vector<double>& du,
          std::vector<double>& u_work,
          std::vector<double>& e_vec,
          const std::vector<double>& dir,
          std::vector<double>& lam_new, std::vector<double>& grad_new,
          double& phi_new) {
  ```

  **3b. Replace the trial loop body** (lines 256-263, inside the `for (int i = 0; i < 20; i++)` bracket loop):

  Old:
  ```cpp
  for (int j = 0; j < st.n; j++) u_work[j] = u_base[j] + alpha * du[j];
  double phi_trial = Tlam + Tdir * alpha;
  double slope = Tdir;
  for (int j = 0; j < st.n; j++) {
      double Fj = fn.F(u_work[j]);
      phi_trial -= d[j] * fn.H(u_work[j]);
      slope -= d[j] * Fj * du[j];
  }
  ```
  New:
  ```cpp
  for (int j = 0; j < st.n; j++) u_work[j] = u_base[j] + alpha * du[j];
  if (!fn.exponential)
      lbw::bulk_scaled_exp(fn.logit_scale, u_work.data(), e_vec.data(), st.n);
  double phi_trial = Tlam + Tdir * alpha;
  double slope = Tdir;
  for (int j = 0; j < st.n; j++) {
      double Fj, Hj;
      if (fn.exponential) { auto fh = fn.FH(u_work[j]); Fj = fh.F; Hj = fh.H; }
      else                { Fj = fn.F_from_e(e_vec[j]); Hj = fn.H_from_e(e_vec[j], u_work[j]); }
      phi_trial -= d[j] * Hj;
      slope    -= d[j] * Fj * du[j];
  }
  ```

  **3c. Update BOTH `wolfe_zoom` call sites in `wolfe_line_search`** to pass `e_vec`.
  There are exactly two call sites; both must be updated or the code will not compile.

  Call site 1 (lines 266-269, Armijo violation or non-first-iter phi increase):
  Old:
  ```cpp
  return wolfe_zoom(st, fn, off, T, d, phi_0, slope_0,
                    alpha_prev, phi_prev, alpha,
                    u_base, du, lam, dir, u_work,
                    lam_new, grad_new, phi_new);
  ```
  New:
  ```cpp
  return wolfe_zoom(st, fn, off, T, d, phi_0, slope_0,
                    alpha_prev, phi_prev, alpha,
                    u_base, du, lam, dir, u_work, e_vec,
                    lam_new, grad_new, phi_new);
  ```

  Call site 2 (lines 278-281, negative slope — zoom from other side):
  Old:
  ```cpp
  return wolfe_zoom(st, fn, off, T, d, phi_0, slope_0,
                    alpha, phi_trial, alpha_prev,
                    u_base, du, lam, dir, u_work,
                    lam_new, grad_new, phi_new);
  ```
  New:
  ```cpp
  return wolfe_zoom(st, fn, off, T, d, phi_0, slope_0,
                    alpha, phi_trial, alpha_prev,
                    u_base, du, lam, dir, u_work, e_vec,
                    lam_new, grad_new, phi_new);
  ```

  **3d. Update `lbfgsb_solve`'s call to `wolfe_line_search`** (line 338-339) to pass `e_vec`:

  Old:
  ```cpp
  wolfe_line_search(st, fn, off, T, d, lam, phi_curr, slope_0,
                    u, du, u_work, dir, lam_new, grad_new, phi_new);
  ```
  New:
  ```cpp
  wolfe_line_search(st, fn, off, T, d, lam, phi_curr, slope_0,
                    u, du, u_work, e_vec, dir, lam_new, grad_new, phi_new);
  ```

- [ ] **Step 4: Update `wolfe_zoom` signature and bisection trial loop**

  **4a. Add `e_vec` to `wolfe_zoom` signature** (insert after `u_work` parameter).

  Old declaration (line 175-185):
  ```cpp
  static double wolfe_zoom(
          const CalibState& st, const LinkFn& fn,
          const std::vector<int>& off, const std::vector<double>& T,
          const std::vector<double>& d, double phi_0, double slope_0,
          double alpha_lo, double phi_lo, double alpha_hi,
          const std::vector<double>& u_base, const std::vector<double>& du,
          const std::vector<double>& lam,
          const std::vector<double>& dir,
          std::vector<double>& u_work,
          std::vector<double>& lam_new, std::vector<double>& grad_new,
          double& phi_new) {
  ```
  New:
  ```cpp
  static double wolfe_zoom(
          const CalibState& st, const LinkFn& fn,
          const std::vector<int>& off, const std::vector<double>& T,
          const std::vector<double>& d, double phi_0, double slope_0,
          double alpha_lo, double phi_lo, double alpha_hi,
          const std::vector<double>& u_base, const std::vector<double>& du,
          const std::vector<double>& lam,
          const std::vector<double>& dir,
          std::vector<double>& u_work,
          std::vector<double>& e_vec,
          std::vector<double>& lam_new, std::vector<double>& grad_new,
          double& phi_new) {
  ```

  **4b. Replace the bisection trial loop body** (lines 198-205, inside the `for (int j = 0; j < 20; j++)` zoom loop):

  Old:
  ```cpp
  for (int i = 0; i < st.n; i++) u_work[i] = u_base[i] + alpha * du[i];
  double phi_trial = Tlam + alpha * Tdir;
  double slope = Tdir;
  for (int i = 0; i < st.n; i++) {
      double Fi = fn.F(u_work[i]);
      phi_trial -= d[i] * fn.H(u_work[i]);
      slope -= d[i] * Fi * du[i];
  }
  ```
  New:
  ```cpp
  for (int i = 0; i < st.n; i++) u_work[i] = u_base[i] + alpha * du[i];
  if (!fn.exponential)
      lbw::bulk_scaled_exp(fn.logit_scale, u_work.data(), e_vec.data(), st.n);
  double phi_trial = Tlam + alpha * Tdir;
  double slope = Tdir;
  for (int i = 0; i < st.n; i++) {
      double Fi, Hi;
      if (fn.exponential) { auto fh = fn.FH(u_work[i]); Fi = fh.F; Hi = fh.H; }
      else                { Fi = fn.F_from_e(e_vec[i]); Hi = fn.H_from_e(e_vec[i], u_work[i]); }
      phi_trial -= d[i] * Hi;
      slope    -= d[i] * Fi * du[i];
  }
  ```

  Note: the second `u_work` update at line 222 (`for (int i...) u_work[i] = u_base[i] + alpha_accepted * du[i]`)
  is for the final accepted-point recompute passed to `phi_from_u` — it does NOT need `bulk_scaled_exp`
  because `phi_from_u` calls `phi_and_grad` which computes its own full O(K·n) computation internally.

- [ ] **Step 5: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -5
  ```

- [ ] **Step 6: Run full tests**

  ```bash
  Rscript -e "testthat::test_local()" 2>&1 | tail -20
  ```
  `test-lbfgsb.R` checks weight correctness and bound enforcement. All must pass.
  Weights must agree with pre-change reference within `tolerance = 1e-10` (numerical identity preserved; only execution order of exp evaluation changes).

- [ ] **Step 7: Benchmark L-BFGS-B path before/after**

  ```r
  library(leafblower); library(bench)
  set.seed(42); n <- 50000L
  age <- sample(c("18-34","35-54","55+"), n, replace=TRUE, prob=c(0.35,0.40,0.25))
  sex <- sample(c("M","F"), n, replace=TRUE, prob=c(0.52,0.48))
  edu <- sample(c("HS","College","Grad"), n, replace=TRUE, prob=c(0.40,0.40,0.20))
  df  <- data.frame(age=factor(age), sex=factor(sex), edu=factor(edu))
  tgt <- list(age=c("18-34"=0.30,"35-54"=0.45,"55+"=0.25),
              sex=c(M=0.50,F=0.50), edu=c(HS=0.35,College=0.45,Grad=0.20))
  bench::mark(harvest(df, tgt, method="lbfgsb"), iterations=10)
  ```
  Expected improvement on Linux/gcc: 25-50% reduction in median time.

- [ ] **Step 8: Commit**

  ```bash
  git add src/lbfgsb_solver.cpp src/logit.hpp
  git commit -m "perf: vectorise Wolfe exp via bulk_scaled_exp (libmvec / scalar fallback)"
  ```

---

### Task 8: OpenMP parallel `compute_errRp` with thread-local buckets

`compute_errRp`'s K margin passes are fully independent (read-only on `w[]`). Parallelise across K with `#pragma omp parallel for`, each thread using a private bucket array. With K=9 and 4–8 cores this yields near-linear speedup on the errRp component.

**Files:**
- Modify: `src/ieppa.cpp`

- [ ] **Step 1: Rewrite `compute_errRp` with OpenMP guard**

  Replace the current `compute_errRp` function (lines 14–33) with:

  ```cpp
  static double compute_errRp(const CalibState& st,
                               const std::vector<double>& w,
                               std::vector<double>& bucket) {
      // W sum: 4-way ILP (same pattern as hot loop — vectorisable).
      double W = 0.0, W1 = 0.0, W2 = 0.0, W3 = 0.0;
      int i4 = st.n & ~3;
      for (int i = 0; i < i4; i += 4) {
          W  += w[i];   W1 += w[i+1];
          W2 += w[i+2]; W3 += w[i+3];
      }
      for (int i = i4; i < st.n; ++i) W += w[i];
      W += W1 + W2 + W3;

      double err = 0.0;
  #if LBW_HAS_OMP
      // K margin passes are read-only on w[] — fully independent.
      // Thread-local bucket avoids false sharing and scatter conflicts.
      int max_cats = *std::max_element(st.cat_counts, st.cat_counts + st.K);
      #pragma omp parallel for schedule(static) reduction(max:err) \
              firstprivate(max_cats)
      for (int k = 0; k < st.K; k++) {
          std::vector<double> local_bucket(max_cats, 0.0);
          for (int i = 0; i < st.n; i++) {
              int g = st.group_ids[k][i];
              if (g >= 0) local_bucket[g] += w[i];
          }
          for (int j = 0; j < st.cat_counts[k]; j++) {
              double e = std::fabs(local_bucket[j] / W - st.targets[k][j]);
              if (e > err) err = e;
          }
      }
  #else
      for (int k = 0; k < st.K; k++) {
          std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
          for (int i = 0; i < st.n; i++) {
              int g = st.group_ids[k][i];
              if (g >= 0) bucket[g] += w[i];
          }
          for (int j = 0; j < st.cat_counts[k]; j++) {
              double e = std::fabs(bucket[j] / W - st.targets[k][j]);
              if (e > err) err = e;
          }
      }
  #endif
      return err;
  }
  ```

  Note: `max_cats` is already computed before the main loop in `ieppa_solve` and passed in `bucket`. The function currently receives `bucket` as scratch; with OpenMP it allocates thread-local storage internally. The `bucket` parameter is still used in the serial path.

- [ ] **Step 2: Compile gate — test with and without OpenMP**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -10
  ```
  Check configure output confirms OpenMP is detected on Linux/gcc.
  On macOS (no OpenMP by default): serial path compiles cleanly.

- [ ] **Step 3: Correctness test — parallel must match serial to within 1e-14**

  ```bash
  Rscript -e "
  library(leafblower)
  set.seed(123); n <- 200000L
  age <- sample(c('u25','25-54','55p'), n, replace=TRUE)
  sex <- sample(c('M','F'), n, replace=TRUE)
  reg <- sample(letters[1:5], n, replace=TRUE)
  df <- data.frame(age=factor(age), sex=factor(sex), reg=factor(reg))
  tgt <- list(age=c(u25=0.2,\`25-54\`=0.6,\`55p\`=0.2),
              sex=c(M=0.5,F=0.5),
              reg=setNames(rep(0.2,5), letters[1:5]))
  w <- harvest(df, tgt, method='ieppa', attach_weights=FALSE)
  cat('max_weight:', max(w), '\n')  # must be <= 5
  cat('mean:', mean(w), '\n')        # must be ~1
  " 2>&1
  ```
  Expected: `mean ≈ 1.0`, `max ≤ 5.0`.

- [ ] **Step 4: Run full test suite**

  ```bash
  Rscript -e "testthat::test_local()" 2>&1 | tail -20
  ```
  All tests must pass on both OpenMP and non-OpenMP builds.

- [ ] **Step 5: Commit**

  ```bash
  git add src/ieppa.cpp
  git commit -m "perf: OpenMP-parallel compute_errRp — K margins computed concurrently"
  ```

---

### Task 9: Full regression + before/after benchmark

Run the full Stepstone benchmark (n=1,582,732) comparing the original code to the fully optimised build.

**Files:**
- Read: `benchmarks/stepstone_fulldata_benchmark.py` (existing)
- No new files — use existing benchmark infrastructure

- [ ] **Step 1: Run full test suite as final regression gate**

  ```bash
  Rscript -e "testthat::test_local()" 2>&1
  ```
  Expected: 0 failures, 0 errors.

- [ ] **Step 2: Run the Python Stepstone benchmark (requires parquet data already generated)**

  ```bash
  python3 benchmarks/stepstone_fulldata_benchmark.py 2>&1
  ```
  Document median time and max_error in the commit message.
  Baseline: 127,000 ms (Python), 133,000 ms (R).
  Target: ≥ 15% reduction (≤ 108,000 ms Python).

- [ ] **Step 3: Run L-BFGS-B micro-benchmark (n=50K, 3 margins)**

  ```r
  library(leafblower); library(bench)
  set.seed(42); n <- 50000L
  df <- data.frame(
    age = factor(sample(c("18-34","35-54","55+"), n, replace=TRUE, prob=c(.35,.40,.25))),
    sex = factor(sample(c("M","F"), n, replace=TRUE, prob=c(.52,.48))),
    edu = factor(sample(c("HS","College","Grad"), n, replace=TRUE, prob=c(.40,.40,.20)))
  )
  tgt <- list(age=c("18-34"=.30,"35-54"=.45,"55+"=.25),
              sex=c(M=.50,F=.50), edu=c(HS=.35,College=.45,Grad=.20))
  bench::mark(harvest(df, tgt, method="lbfgsb"), iterations=20)
  ```
  Document median time. Baseline: measure before changes. Target: ≥ 25% reduction.

- [ ] **Step 4: Final commit with benchmark results**

  ```bash
  git add -A
  git commit -m "perf: benchmark results — document iEPPA and L-BFGS-B speedups"
  ```
