# Code Review Fixes — ieppa.cpp + configure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix all blocking and required items from the post-merge code review of the C++ performance optimisation branch.

**Architecture:** Six self-contained tasks, each touching one file. All changes compile-gated after every file. No new public APIs introduced.

**Tech Stack:** C++17, POSIX sh (`configure`), R package (`R CMD INSTALL`), testthat

---

## File Map

| File | Tasks |
|------|-------|
| `src/ieppa.cpp` | Task 1 (include order), Task 2 (convergence gate), Task 3 (W-sum helper) |
| `configure` | Task 4 (OMP simd probe, run-gate, flag split) |
| `src/Makevars.in` | Task 4 (add `@SIMD_FLAGS@` placeholder) |
| `src/lbw_config.h.in` | Task 4 (no change — generated file) |
| `src/ieppa.cpp` (box loop) | Task 5 (std::clamp) |

---

## Task 1: Fix `lbw_config.h` include order in `ieppa.cpp`

**Files:** `src/ieppa.cpp` only — move one line.

`lbw_config.h` is a configure-generated feature-detection header. It must appear before any other project header so that `LBW_HAS_*` macros are defined before anything that might transitively use them. `lbfgsb_solver.cpp` already has it first. `ieppa.cpp` does not. This is a latent bug that becomes active the moment any dependency header checks `LBW_HAS_*`.

**DoD:**
- [ ] `lbw_config.h` is the first `#include` in `ieppa.cpp`
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)`
- [ ] `Rscript -e "testthat::test_local()" 2>&1 | tail -3` = `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 49 ]`
- [ ] Committed: `fix: lbw_config.h must be first include in ieppa.cpp`

- [ ] **Step 1: Read and confirm current include order**

  ```bash
  head -8 src/ieppa.cpp
  ```

  Expected (current wrong state):
  ```cpp
  #include "ieppa.hpp"
  #include "leafblower.h"
  #include "lbw_config.h"   // ← line 3, wrong
  ```

- [ ] **Step 2: Move `lbw_config.h` to line 1**

  New top of file:
  ```cpp
  #include "lbw_config.h"
  #include "ieppa.hpp"
  #include "leafblower.h"
  #include <cmath>
  #include <cstdio>
  #include <algorithm>
  #include <vector>
  ```

- [ ] **Step 3: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 4: Test suite**

  ```bash
  Rscript -e "testthat::test_local()" 2>&1 | tail -3
  ```

- [ ] **Step 5: Commit**

  ```bash
  git add src/ieppa.cpp
  git commit -m "fix: lbw_config.h must be first include in ieppa.cpp"
  ```

---

## Task 2: Fix convergence gate — check on iteration 1

**Files:** `src/ieppa.cpp` only — one condition change.

The loop runs `iter = 1..inner_max_iter`. The check `iter % kErrCheckInterval == 0` fires at 10, 20, 30, … — never at 1–9. A problem that converges in fewer than 10 iterations returns `RK_ERR_NOCONV` with `res.max_error` reflecting the initial point (1.0), not the converged iterate. Adding `iter == 1` ensures the first check occurs at iter 1, catching fast-converging problems correctly and keeping subsequent checks at the 10-iteration cadence.

The comment on `kErrCheckInterval` must be updated to document the iter-1 exception so future readers don't remove it.

**DoD:**
- [ ] `iter == 1` is part of the convergence-check condition
- [ ] `kErrCheckInterval` comment documents the iter-1 exception
- [ ] Compile gate passes
- [ ] Test suite: 49 PASS, 0 FAIL
- [ ] Committed: `fix: check convergence on iter 1 to catch fast-converging problems`

- [ ] **Step 1: Read the current convergence check block**

  Find lines around:
  ```cpp
  if (iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
  ```

- [ ] **Step 2: Update the condition and constant comment**

  Old constant comment:
  ```cpp
  static constexpr int kErrCheckInterval = 10;  // Check convergence every N iterations instead of every 1.
                                                 // compute_errRp costs K O(n) passes — nearly as expensive as a full sweep.
                                                 // Every-10 reduces that overhead by 90% at the cost of ≤9 extra IPF iters.
  ```

  New:
  ```cpp
  static constexpr int kErrCheckInterval = 10;  // Check convergence every N inner iterations.
                                                 // compute_errRp costs K O(n) passes — nearly as expensive as a full sweep.
                                                 // Every-10 reduces that overhead by 90% at the cost of ≤9 extra IPF iters.
                                                 // Exception: always check on iter 1 to catch problems that converge immediately.
  ```

  Old condition:
  ```cpp
  if (iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
  ```

  New:
  ```cpp
  if (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
  ```

- [ ] **Step 3: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 4: Test suite**

  ```bash
  Rscript -e "testthat::test_local()" 2>&1 | tail -3
  ```

- [ ] **Step 5: Commit**

  ```bash
  git add src/ieppa.cpp
  git commit -m "fix: check convergence on iter 1 to catch fast-converging problems"
  ```

---

## Task 3: Extract `sum_weights_ilp` and unify the two ILP W-sum sites

**Files:** `src/ieppa.cpp` only — new static helper, two call sites.

The 4-way ILP W-sum unroll appears at exactly two sites: `compute_errRp` (lines 18–25, variables `W/W1/W2/W3/i4e`) and the margin loop (lines 81–88, variables `W/W1/W2/W3/i4/ni`). These two sites are structurally identical but use different variable names, which is a maintenance hazard — a future change to one will silently miss the other. The normalize block (line 119) and post-loop fixup (line 180) use plain scalar `Wsum += w[i]` loops — a different and correct pattern; those are left unchanged.

Extract a single `sum_weights_ilp` static helper that owns the ILP pattern. Both ILP sites call the helper. This eliminates the naming divergence (`i4e` vs `i4`, `st.n` vs `ni`).

**Trace:** For w = [1.0, 2.0, 3.0, 4.0, 5.0] (n=5):
- `n4 = 5 & ~3 = 4`
- Unrolled: W=1, W1=2, W2=3, W3=4 → W += 5 (tail) → result = 1+5+2+3+4 = 15 ✓
- Scalar: 1+2+3+4+5 = 15 ✓

**DoD:**
- [ ] `sum_weights_ilp` is a file-scope `static` helper declared before `compute_errRp`
- [ ] Both ILP W-sum sites (`compute_errRp` and margin loop) call `sum_weights_ilp(w, st.n)`
- [ ] No `W1`/`W2`/`W3`/`i4e`/`i4`/`ni` variables exist outside the helper
- [ ] Normalize block (`Wsum` scalar loop) and post-loop fixup (`Wsum` scalar loop) are **unchanged**
- [ ] Compile gate passes
- [ ] Weight regression: `max(abs(w_after - w_ref)) < 1e-12` (iEPPA path)
- [ ] Test suite: 49 PASS, 0 FAIL
- [ ] Committed: `refactor: extract sum_weights_ilp helper for the two ILP W-sum sites in ieppa.cpp`

- [ ] **Step 1: Read ieppa.cpp** — confirm exact text of the two ILP unroll blocks (lines 18–25 and 81–88)

- [ ] **Step 2: Add `sum_weights_ilp` before `compute_errRp`**

  Insert immediately before `static double compute_errRp(...)`:

  ```cpp
  // Sum weights with 4-way ILP unroll.
  // Separate accumulators break the loop-carried dependency chain, letting
  // the compiler pipeline four additions in parallel. The tail loop handles n % 4.
  static double sum_weights_ilp(const std::vector<double>& w, int n) {
      double W = 0.0, W1 = 0.0, W2 = 0.0, W3 = 0.0;
      const int n4 = n & ~3;
      for (int i = 0; i < n4; i += 4) {
          W  += w[i];   W1 += w[i+1];
          W2 += w[i+2]; W3 += w[i+3];
      }
      for (int i = n4; i < n; ++i) W += w[i];
      return W + W1 + W2 + W3;
  }
  ```

- [ ] **Step 3: Replace ILP site 1 — `compute_errRp`**

  Remove the inline unroll block (lines 18–25), replace with:
  ```cpp
  double W = sum_weights_ilp(w, st.n);
  ```

- [ ] **Step 4: Replace ILP site 2 — margin loop**

  Remove the inline unroll block (lines 81–88, including `int ni = st.n, i4 = ...`), replace with:
  ```cpp
  double W = sum_weights_ilp(w, st.n);
  ```
  Replace all remaining `ni` references in that scope with `st.n`.

- [ ] **Step 5: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 6: Weight regression (iEPPA)**

  `tests/testthat/task2_ieppa_ref.rds` is a frozen reference computed from the iEPPA solver output during the prior performance optimisation branch. None of Tasks 1–5 change the numerical weight computation (include reordering, convergence-check timing for problems that converge in <10 iters, helper extraction, configure/Makevars edits, and std::clamp are all algebraically neutral for n=10000 with 2 margins). The regression asserts this is still true post-refactor.

  ```bash
  Rscript -e "
  library(leafblower)
  set.seed(42); n <- 10000L
  df  <- data.frame(
    age = factor(sample(c('18-34','35-54','55+'), n, replace=TRUE)),
    sex = factor(sample(c('M','F'), n, replace=TRUE))
  )
  tgt <- list(age=c('18-34'=0.30,'35-54'=0.45,'55+'=0.25), sex=c(M=0.50,F=0.50))
  w_after <- harvest(df, tgt, method='ieppa', attach_weights=FALSE)
  w_ref   <- readRDS('tests/testthat/task2_ieppa_ref.rds')
  delta   <- max(abs(w_after - w_ref))
  cat('delta:', delta, '\n')
  stopifnot(delta < 1e-12)
  cat('PASS\n')
  "
  ```

- [ ] **Step 7: Test suite**

  ```bash
  Rscript -e "testthat::test_local()" 2>&1 | tail -3
  ```

- [ ] **Step 8: Commit**

  ```bash
  git add src/ieppa.cpp
  git commit -m "refactor: extract sum_weights_ilp helper for the two ILP W-sum sites in ieppa.cpp"
  ```

---

## Task 4: Fix configure — separate OMP-simd probe, run binary, split flags

**Files:** `configure` + `src/Makevars.in`

Three issues in configure, fixed together as one atomic configure rewrite:

1. **`LBW_HAS_OMP_SIMD` and `LBW_HAS_OMP` always identical** — the simd-only path (`-fopenmp-simd`) is never probed. A user with `-fopenmp-simd` but no full OpenMP runtime gets no SIMD hints even though their compiler supports them.
2. **OMP probe never runs the binary** — unlike the libmvec probe, compile+link success is the only gate. A build where OpenMP compiles but the runtime is broken would silently set `LBW_HAS_OMP=1`.
3. **`SHLIB_OPENMP_CXXFLAGS` appears twice** in the compile+link command — cargo-cult duplication; any compiler-only flag in that variable would cause a linker error.

Fix strategy:
- Probe `LBW_HAS_OMP_SIMD` first via `-fopenmp-simd` (compile+link+run a `#pragma omp simd` loop).
- Probe `LBW_HAS_OMP` second via `SHLIB_OPENMP_CXXFLAGS` with separate compile and link steps, then run.
- `SIMD_FLAGS` = `-fopenmp-simd` only when simd works but full OMP does not (avoids redundant flag when `-fopenmp` already covers simd).
- Add `@SIMD_FLAGS@` to `Makevars.in` so the flag reaches the compiler when simd-only mode is active.

**Concrete probe logic:**

```sh
# Step A: probe simd-only (-fopenmp-simd, no threading runtime required)
LBW_HAS_OMP_SIMD=0
SIMD_FLAGS=""
cat > ${_lbw_tmp}/simd_test.cpp << 'EOF'
int main() {
    double a[8] = {};
    #pragma omp simd
    for (int i = 0; i < 8; i++) a[i] = i * 2.0;
    return a[0] != 0.0;   /* prevents dead-code elim */
}
EOF
if eval $CXX -fopenmp-simd -x c++ ${_lbw_tmp}/simd_test.cpp \
       -o ${_lbw_tmp}/simd_test 2>/dev/null && \
   ${_lbw_tmp}/simd_test 2>/dev/null; then
    LBW_HAS_OMP_SIMD=1
    SIMD_FLAGS="-fopenmp-simd"   # provisional; cleared if full OMP takes over below
    echo "configure: -fopenmp-simd detected — omp simd hints enabled"
fi

# Step B: probe full OpenMP (threading runtime required for LBW_HAS_OMP)
LBW_HAS_OMP=0
OMP_FLAGS=""
if [ -n "${SHLIB_OPENMP_CXXFLAGS}" ]; then
    cat > ${_lbw_tmp}/omp_test.cpp << 'EOF'
#include <omp.h>
int main() { return omp_get_max_threads() > 0 ? 0 : 1; }
EOF
    # Compile and link in separate steps to avoid flag duplication
    if eval $CXX ${SHLIB_OPENMP_CXXFLAGS} -x c++ -c ${_lbw_tmp}/omp_test.cpp \
           -o ${_lbw_tmp}/omp_test.o 2>/dev/null && \
       eval $CXX ${SHLIB_OPENMP_CXXFLAGS} ${_lbw_tmp}/omp_test.o \
           -o ${_lbw_tmp}/omp_test 2>/dev/null && \
       ${_lbw_tmp}/omp_test 2>/dev/null; then
        OMP_FLAGS="${SHLIB_OPENMP_CXXFLAGS}"
        LBW_HAS_OMP=1
        LBW_HAS_OMP_SIMD=1   # full OMP implies simd
        SIMD_FLAGS=""          # -fopenmp already covers simd; no duplicate flag needed
        echo "configure: OpenMP detected — threading and simd enabled"
    else
        echo "configure: OpenMP not available — serial / simd-only fallback"
    fi
fi
```

**DoD:**
- [ ] `LBW_HAS_OMP_SIMD` is probed independently via `-fopenmp-simd` before the full OMP probe
- [ ] Full OMP probe: compile step and link step are separate `eval $CXX` calls; binary is run
- [ ] `SHLIB_OPENMP_CXXFLAGS` appears at most once per `eval $CXX` invocation
- [ ] `SIMD_FLAGS` is cleared to `""` when full OMP is detected (no redundant `-fopenmp-simd`)
- [ ] `src/Makevars.in` has `@SIMD_FLAGS@` between `@OMP_FLAGS@` and `@MAVX2_FLAG@`
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)`
- [ ] `Rscript -e "testthat::test_local()" 2>&1 | tail -3` = 49 PASS
- [ ] On this machine (no full OMP, has libmvec): `grep LBW_HAS src/lbw_config.h` shows `LBW_HAS_OMP_SIMD 0`, `LBW_HAS_OMP 0`, `LBW_HAS_GLIBC_MVEC 1`
  (OpenMP and fopenmp-simd are not available on this machine per earlier configure output)
- [ ] Committed: `fix: configure — separate simd probe, run OMP binary, split compile/link`

- [ ] **Step 1: Read current configure** — identify exact block to replace (the OpenMP probe section)

- [ ] **Step 2: Replace the OpenMP probe section**

  Replace everything from `# Detect OpenMP:` through `fi` with the two-step probe shown above.

- [ ] **Step 3: Update `src/Makevars.in`**

  Old:
  ```makefile
  PKG_CXXFLAGS = @CXXFLAGS_STD@ @OPT_FLAGS@ @OMP_FLAGS@ @MAVX2_FLAG@ -I. -DSTRICT_R_HEADERS
  ```
  New:
  ```makefile
  PKG_CXXFLAGS = @CXXFLAGS_STD@ @OPT_FLAGS@ @OMP_FLAGS@ @SIMD_FLAGS@ @MAVX2_FLAG@ -I. -DSTRICT_R_HEADERS
  ```

- [ ] **Step 4: Update the sed command at bottom of configure**

  Old sed has 5 substitutions. Add `s|@SIMD_FLAGS@|${SIMD_FLAGS}|`:
  ```sh
  sed "s|@CXXFLAGS_STD@|${CXXFLAGS_STD}|;\
       s|@OPT_FLAGS@|${OPT_FLAGS}|;\
       s|@OMP_FLAGS@|${OMP_FLAGS}|;\
       s|@SIMD_FLAGS@|${SIMD_FLAGS}|;\
       s|@MAVX2_FLAG@|${MAVX2_FLAG}|;\
       s|@MVEC_LIBS@|${MVEC_LIBS}|" \
      src/Makevars.in > src/Makevars
  ```

- [ ] **Step 5: Compile gate and substitution verification**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```
  Verify `src/lbw_config.h` values are correct for this machine:
  ```bash
  grep "LBW_HAS" src/lbw_config.h
  ```
  Verify `@SIMD_FLAGS@` was substituted into the generated `src/Makevars`:
  ```bash
  grep "^PKG_CXXFLAGS" src/Makevars
  ```
  The `@SIMD_FLAGS@` placeholder must not appear literally — it must be gone (replaced by either `-fopenmp-simd` or the empty string). On this machine (no -fopenmp-simd support per prior configure output) it will be an empty substitution, so the line will show `@OMP_FLAGS@` and `@MAVX2_FLAG@` regions also collapsed to their values.

- [ ] **Step 6: Test suite**

  ```bash
  Rscript -e "testthat::test_local()" 2>&1 | tail -3
  ```

- [ ] **Step 7: Commit**

  ```bash
  git add configure src/Makevars.in
  git commit -m "fix: configure — separate simd probe, run OMP binary, split compile/link"
  ```

---

## Task 5: Replace `std::max`/`std::min` with `std::clamp` in SIMD loops

**Files:** `src/ieppa.cpp` only — two loop bodies.

`#pragma omp simd` instructs the compiler to vectorise the loop. Some compilers (MSVC, older GCC) do not vectorise through `std::max(std::min(...))` call chains. `std::clamp` (C++17, already required by this package) produces the same result from a single function that compilers pattern-match to conditional moves (`CMOV` / `vmaxpd`/`vminpd`), giving the vectoriser a cleaner signal.

**Trace for yi = 0.3, lo = 0.1, hi = 5.0:**
- `std::clamp(0.3, 0.1, 5.0)` = `0.3` ✓
- `std::max(0.1, std::min(5.0, 0.3))` = `std::max(0.1, 0.3)` = `0.3` ✓ (identical)

**Trace for yi = -0.5:**
- `std::clamp(-0.5, 0.1, 5.0)` = `0.1` ✓

**Trace for yi = 7.0:**
- `std::clamp(7.0, 0.1, 5.0)` = `5.0` ✓

**DoD:**
- [ ] Both `std::max(lo, std::min(hi, ...))` expressions in `ieppa.cpp` replaced with `std::clamp`
- [ ] `#include <algorithm>` present (already is — `std::clamp` is in `<algorithm>`)
- [ ] Compile gate passes
- [ ] Weight regression `< 1e-12`
- [ ] Test suite: 49 PASS
- [ ] Committed: `refactor: std::clamp in box-projection loop for cleaner SIMD hint`

- [ ] **Step 1: Find both clamp sites**

  ```bash
  grep -n "std::max.*std::min\|std::min.*std::max" src/ieppa.cpp
  ```
  Expect two hits: one in the box-projection loop (main iteration) and one in the post-loop fixup.

- [ ] **Step 2: Replace box-projection loop body**

  Old:
  ```cpp
  double wc = std::max(lo, std::min(hi, yi));
  ```
  New:
  ```cpp
  double wc = std::clamp(yi, lo, hi);
  ```

- [ ] **Step 3: Replace post-loop fixup body**

  Old:
  ```cpp
  double wc = std::max(lo, std::min(hi, w[i]));
  ```
  New:
  ```cpp
  double wc = std::clamp(w[i], lo, hi);
  ```

- [ ] **Step 4: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 5: Weight regression**

  ```bash
  Rscript -e "
  library(leafblower)
  set.seed(42); n <- 10000L
  df  <- data.frame(
    age = factor(sample(c('18-34','35-54','55+'), n, replace=TRUE)),
    sex = factor(sample(c('M','F'), n, replace=TRUE))
  )
  tgt <- list(age=c('18-34'=0.30,'35-54'=0.45,'55+'=0.25), sex=c(M=0.50,F=0.50))
  w_after <- harvest(df, tgt, method='ieppa', attach_weights=FALSE)
  w_ref   <- readRDS('tests/testthat/task2_ieppa_ref.rds')
  stopifnot(max(abs(w_after - w_ref)) < 1e-12)
  cat('PASS\n')
  "
  ```

- [ ] **Step 6: Test suite**

  ```bash
  Rscript -e "testthat::test_local()" 2>&1 | tail -3
  ```

- [ ] **Step 7: Commit**

  ```bash
  git add src/ieppa.cpp
  git commit -m "refactor: std::clamp in box-projection loop for cleaner SIMD hint"
  ```
