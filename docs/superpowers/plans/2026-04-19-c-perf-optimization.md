# leafblower C Performance Optimization — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Speed up leafblower's iEPPA and L-BFGS-B solvers with no new package dependencies. Target: ~25-35% iEPPA speedup, ~50% L-BFGS-B speedup.

**Architecture:** Two phases. Phase 1: zero-dep scalar/compiler wins (logit shared-exp fix, W-sum hoisting, errRp frequency, `#pragma omp simd`). Phase 2: configure-detected libmvec (vectorized exp for L-BFGS-B) and OpenMP threading for `compute_errRp`.

**Tech Stack:** C++17; glibc libmvec (optional, Linux/gcc only, configure-detected); OpenMP (optional, via R's `SHLIB_OPENMP_CXXFLAGS`, no new package dep).

## TDD Convention

Every task that introduces a **new function** has a RED step before code exists:
- Write the test → confirm it fails (function not found or assertion fails)
- Implement the function
- Run the test → GREEN

Every task that **modifies existing behavior** has a reference snapshot before the edit:
- Save reference weights → implement → assert `max(abs(w_new - w_ref)) < tol`

**Compile gate (mandatory after every single file change):**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -1
```
Expected: `* DONE (leafblower)`. Never advance to the next file before this passes.

**Test gate (mandatory before every commit):**
```bash
Rscript -e "testthat::test_local()" 2>&1 | grep -E "^(FAIL|ERROR)" | wc -l
```
Expected: `0`. Any non-zero count = STOP, diagnose, fix.

---

## Phase 1 — Zero-dep scalar wins

---

### Task 1a: Add FH() to LinkFn in logit.hpp

Both `F(u)` and `H(u)` call `safe_exp(logit_scale * u)` independently. `FH(u)` computes one shared `e` and returns both, halving transcendental calls in Wolfe loops.

**Files:** `src/logit.hpp` only — no other file touched.

**DoD:**
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)`
- [ ] `Rscript -e "testthat::test_local()" 2>&1 | grep -c 'FAIL\|ERROR'` = `0`
- [ ] FH formula verified: `FH(1.0).F == F(1.0)` and `FH(1.0).H == H(1.0)` (traced in Step 1)
- [ ] Committed with message `perf: add FH() to LinkFn — single exp for Wolfe loop`

- [ ] **Step 1 (RED): Capture timing baseline and save reference weights — BEFORE any code change**

  This baseline is referenced in Task 9 Step 4 DoD ("≥25% reduction"). It must be
  measured BEFORE Task 1a modifies any source file. Save it now.

  ```bash
  Rscript -e "
  library(leafblower); library(bench)
  # ── lbfgsb timing baseline ──────────────────────────────────────────────────
  set.seed(42); n <- 50000L
  df_bench <- data.frame(
    age = factor(sample(c('18-34','35-54','55+'), n, replace=TRUE, prob=c(.35,.40,.25))),
    sex = factor(sample(c('M','F'), n, replace=TRUE, prob=c(.52,.48))),
    edu = factor(sample(c('HS','College','Grad'), n, replace=TRUE, prob=c(.40,.40,.20)))
  )
  tgt_bench <- list(age=c('18-34'=.30,'35-54'=.45,'55+'=.25),
                    sex=c(M=.50,F=.50), edu=c(HS=.35,College=.45,Grad=.20))
  bm <- bench::mark(harvest(df_bench, tgt_bench, method='lbfgsb'), iterations=10)
  cat('lbfgsb baseline median ms:', as.numeric(bm\$median) * 1000, '\n')
  saveRDS(bm\$median, 'tests/testthat/lbfgsb_baseline_time.rds')

  # ── regression weight reference ─────────────────────────────────────────────
  set.seed(7); n <- 5000L
  df  <- data.frame(x = factor(sample(c('A','B','C'), n, replace=TRUE)))
  tgt <- list(x = c(A=0.4, B=0.35, C=0.25))
  w_before <- harvest(df, tgt, method='lbfgsb', attach_weights=FALSE)
  saveRDS(w_before, 'tests/testthat/task1_ref.rds')
  cat('ref saved. mean:', mean(w_before), 'max:', max(w_before), '\n')
  "
  ```

- [ ] **Step 2 (GREEN): Add FH() to logit.hpp after the H() method (line 52, before `};`)**

  ```cpp
  // FH(u): F and H from a single safe_exp call.
  // Halves transcendental evaluations in the Wolfe inner loop.
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

  Trace for u=1.0, L=0.2, U=5.0, logit_scale=1.5:
  - `e = exp(1.5 * 1.0) = exp(1.5) ≈ 4.4817`
  - `denom = (U-1) + (1-L)*e = 4 + 0.8*4.4817 ≈ 7.5854`
  - `f = (L*(U-1) + U*(1-L)*e) / denom = (0.8 + 5*0.8*4.4817) / 7.5854 ≈ 18.726/7.5854 ≈ 2.469`
  - `h = L*u + (U-L)/logit_scale * log(denom/(U-L)) = 0.2 + (4.8/1.5)*log(7.5854/4.8) ≈ 0.2 + 3.2*0.456 ≈ 1.659`
  - These must equal `F(1.0) ≈ 2.469` and `H(1.0) ≈ 1.659` respectively (verify against the actual logit.hpp F() and H() formulas before implementing FH()).

- [ ] **Step 3: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 4 (GREEN check): Assert FH() matches F()+H() and no weight regression**

  ```bash
  Rscript -e "
  library(leafblower)
  set.seed(7); n <- 5000L
  df  <- data.frame(x = factor(sample(c('A','B','C'), n, replace=TRUE)))
  tgt <- list(x = c(A=0.4, B=0.35, C=0.25))
  w_after <- harvest(df, tgt, method='lbfgsb', attach_weights=FALSE)
  w_ref   <- readRDS('tests/testthat/task1_ref.rds')
  delta   <- max(abs(w_after - w_ref))
  cat('delta:', delta, '\n')
  stopifnot(delta < 1e-12)
  cat('PASS\n')
  "
  ```
  (FH() is not yet called anywhere — this confirms logit.hpp change doesn't break the build.)

- [ ] **Step 5: Commit**

  ```bash
  git add src/logit.hpp tests/testthat/task1_ref.rds tests/testthat/lbfgsb_baseline_time.rds
  git commit -m "perf: add FH() to LinkFn — single exp for Wolfe loop"
  ```

---

### Task 1b: Use FH() in wolfe_line_search trial loop

Replace the two separate `fn.F()` / `fn.H()` calls in the Wolfe bracketing trial loop with one `fn.FH()` call.

**Files:** `src/lbfgsb_solver.cpp` only — `wolfe_line_search` function, lines 259–263.

**DoD:**
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)`
- [ ] `max(abs(w_after - w_ref)) < 1e-12` (same `task1_ref.rds` from Task 1a)
- [ ] Committed with message `perf: use FH() in wolfe_line_search trial loop`

- [ ] **Step 1: Replace the trial loop body in wolfe_line_search (lines 259–263)**

  Old (lines 259–263 inside `for (int i = 0; i < 20; i++)` loop):
  ```cpp
  for (int j = 0; j < st.n; j++) {
      double Fj = fn.F(u_work[j]);
      phi_trial -= d[j] * fn.H(u_work[j]);
      slope -= d[j] * Fj * du[j];
  }
  ```
  New:
  ```cpp
  for (int j = 0; j < st.n; j++) {
      auto fh = fn.FH(u_work[j]);
      phi_trial -= d[j] * fh.H;
      slope    -= d[j] * fh.F * du[j];
  }
  ```

- [ ] **Step 2: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 3 (GREEN): Assert weights unchanged**

  ```bash
  Rscript -e "
  library(leafblower)
  set.seed(7); n <- 5000L
  df  <- data.frame(x = factor(sample(c('A','B','C'), n, replace=TRUE)))
  tgt <- list(x = c(A=0.4, B=0.35, C=0.25))
  w_after <- harvest(df, tgt, method='lbfgsb', attach_weights=FALSE)
  w_ref   <- readRDS('tests/testthat/task1_ref.rds')
  stopifnot(max(abs(w_after - w_ref)) < 1e-12)
  cat('PASS\n')
  "
  ```

- [ ] **Step 4: Commit**

  ```bash
  git add src/lbfgsb_solver.cpp
  git commit -m "perf: use FH() in wolfe_line_search trial loop"
  ```

---

### Task 1c: Use FH() in wolfe_zoom bisection trial loop

Same FH() substitution in `wolfe_zoom`, which is called when Wolfe bracketing finds an interval to zoom into.

**Files:** `src/lbfgsb_solver.cpp` only — `wolfe_zoom` function, lines 201–205.

**DoD:**
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)`
- [ ] `max(abs(w_after - w_ref)) < 1e-12` (same `task1_ref.rds`)
- [ ] `Rscript -e "testthat::test_local()" 2>&1 | grep -c 'FAIL\|ERROR'` = `0`
- [ ] Committed with message `perf: use FH() in wolfe_zoom bisection loop`

- [ ] **Step 1: Replace the trial loop body in wolfe_zoom (lines 201–205)**

  Old:
  ```cpp
  for (int i = 0; i < st.n; i++) {
      double Fi = fn.F(u_work[i]);
      phi_trial -= d[i] * fn.H(u_work[i]);
      slope -= d[i] * Fi * du[i];
  }
  ```
  New:
  ```cpp
  for (int i = 0; i < st.n; i++) {
      auto fh = fn.FH(u_work[i]);
      phi_trial -= d[i] * fh.H;
      slope    -= d[i] * fh.F * du[i];
  }
  ```

- [ ] **Step 2: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 3 (GREEN): Assert weights unchanged + full test suite**

  ```bash
  Rscript -e "
  library(leafblower)
  set.seed(7); n <- 5000L
  df  <- data.frame(x = factor(sample(c('A','B','C'), n, replace=TRUE)))
  tgt <- list(x = c(A=0.4, B=0.35, C=0.25))
  w_after <- harvest(df, tgt, method='lbfgsb', attach_weights=FALSE)
  w_ref   <- readRDS('tests/testthat/task1_ref.rds')
  stopifnot(max(abs(w_after - w_ref)) < 1e-12)
  cat('PASS\n')
  " && Rscript -e "testthat::test_local()" 2>&1 | tail -5
  ```

- [ ] **Step 4: Commit**

  ```bash
  git add src/lbfgsb_solver.cpp
  git commit -m "perf: use FH() in wolfe_zoom bisection loop"
  ```

---

### Task 2: Hoist W-sum from bucket loop + 4-way ILP unroll (ieppa.cpp)

`W += w[i]` is serialised inside the scatter-add loop. Separating into a 4-way ILP reduction lets the compiler vectorise the sum. Two sites: the IPF margin loop (lines 70–76) and `compute_errRp` (lines 17–18).

**Files:** `src/ieppa.cpp` only.

**DoD:**
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)`
- [ ] `max(abs(w_after - w_ref)) < 1e-6` (acceptable: FP reassociation changes sum order)
- [ ] `Rscript -e "testthat::test_local()" 2>&1 | grep -c 'FAIL\|ERROR'` = `0`
- [ ] Committed with message `perf: hoist W-sum from scatter loop, 4-way ILP unroll`

- [ ] **Step 1 (RED): Save iEPPA reference weights before any change**

  ```bash
  Rscript -e "
  library(leafblower)
  set.seed(1); n <- 10000L
  df  <- data.frame(
    age = factor(sample(c('A','B','C'), n, replace=TRUE)),
    sex = factor(sample(c('M','F'), n, replace=TRUE))
  )
  tgt <- list(age=c(A=0.4,B=0.35,C=0.25), sex=c(M=0.5,F=0.5))
  w_ref <- harvest(df, tgt, method='ieppa', attach_weights=FALSE)
  saveRDS(w_ref, 'tests/testthat/task2_ieppa_ref.rds')
  cat('ref saved. mean:', mean(w_ref), '\n')
  "
  ```

- [ ] **Step 2: Replace the combined W + bucket loop (lines 70–76) with two separate loops**

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
  // W sum separated from scatter-add so the compiler can vectorise it.
  double W = 0.0, W1 = 0.0, W2 = 0.0, W3 = 0.0;
  int ni = st.n, i4 = ni & ~3;
  for (int i = 0; i < i4; i += 4) {
      W  += w[i];   W1 += w[i+1];
      W2 += w[i+2]; W3 += w[i+3];
  }
  for (int i = i4; i < ni; ++i) W += w[i];
  W += W1 + W2 + W3;

  // Bucket scatter-add: write aliases prevent vectorisation.
  std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
  for (int i = 0; i < ni; i++) {
      int g = st.group_ids[k][i];
      if (g >= 0) bucket[g] += w[i];
  }
  ```

- [ ] **Step 3: Apply the same ILP pattern to compute_errRp W-sum (lines 17–18)**

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

- [ ] **Step 4: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 5 (GREEN): Assert weights within tolerance + test suite**

  ```bash
  Rscript -e "
  library(leafblower)
  set.seed(1); n <- 10000L
  df  <- data.frame(
    age = factor(sample(c('A','B','C'), n, replace=TRUE)),
    sex = factor(sample(c('M','F'), n, replace=TRUE))
  )
  tgt <- list(age=c(A=0.4,B=0.35,C=0.25), sex=c(M=0.5,F=0.5))
  w_after <- harvest(df, tgt, method='ieppa', attach_weights=FALSE)
  w_ref   <- readRDS('tests/testthat/task2_ieppa_ref.rds')
  delta   <- max(abs(w_after - w_ref))
  cat('delta:', delta, '\n')   # FP reorder ok; tolerance is 1e-6
  stopifnot(delta < 1e-6)
  cat('PASS\n')
  " && Rscript -e "testthat::test_local()" 2>&1 | tail -5
  ```

- [ ] **Step 6: Commit**

  ```bash
  git add src/ieppa.cpp tests/testthat/task2_ieppa_ref.rds
  git commit -m "perf: hoist W-sum from scatter loop, 4-way ILP unroll"
  ```

---

### Task 3: Reduce errRp check frequency to every 10 iterations (ieppa.cpp)

`compute_errRp` costs ~K O(n) passes per call — nearly as expensive as a full IPF sweep. Checking every 10 iterations cuts that overhead by 90%. Worst-case overshoot: 9 extra iterations past true convergence.

**Files:** `src/ieppa.cpp` only.

**DoD:**
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)`
- [ ] `Rscript -e "testthat::test_local()" 2>&1 | grep -c 'FAIL\|ERROR'` = `0`
- [ ] Margin error ≤ tol_abs after calibration (confirmed by existing `test-bounded-convergence.R`)
- [ ] Committed with message `perf: check errRp every 10 iters — ~30% iEPPA speedup`

- [ ] **Step 1: Add the interval constant after `static constexpr int kMaxFixupIterations = 20;` (line 43)**

  ```cpp
  // Check convergence every N iterations instead of every 1.
  // compute_errRp costs K O(n) passes — nearly as expensive as a full sweep.
  // Every-10 reduces that overhead by 90% at the cost of ≤9 extra IPF iters.
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
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 4 (GREEN): Verify convergence preserved**

  ```bash
  Rscript -e "
  library(leafblower)
  set.seed(1); n <- 10000L
  df  <- data.frame(
    age = factor(sample(c('A','B','C'), n, replace=TRUE)),
    sex = factor(sample(c('M','F'), n, replace=TRUE))
  )
  tgt <- list(age=c(A=0.4,B=0.35,C=0.25), sex=c(M=0.5,F=0.5))
  r   <- harvest(df, tgt, method='ieppa')
  cat('max_error:', r\$max_error, '\n')  # must be < 1e-3 (default tol)
  stopifnot(r\$max_error < 1e-3)
  cat('PASS\n')
  " && Rscript -e "testthat::test_local()" 2>&1 | tail -5
  ```

- [ ] **Step 5: Commit**

  ```bash
  git add src/ieppa.cpp
  git commit -m "perf: check errRp every 10 iters — ~30% iEPPA speedup"
  ```

---

### Task 4a: Create stub lbw_config.h and wire include into ieppa.cpp

Required before Task 4b/4c. The `#pragma omp simd` guards reference `LBW_HAS_OMP_SIMD`, which lives in `lbw_config.h`. Without a pre-existing header, the `#if LBW_HAS_OMP_SIMD` evaluates an undefined identifier — a `-Wundef` error under CRAN checks. The real values are filled in by Task 5.

**Files:** Create `src/lbw_config.h`; modify `src/ieppa.cpp` (add include only).

**DoD:**
- [ ] `src/lbw_config.h` exists with all features set to 0 (generated locally; NOT committed)
- [ ] `src/ieppa.cpp` has `#include "lbw_config.h"` after existing includes
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)`
- [ ] `lbw_config.h` in `src/.gitignore` (so it is never accidentally staged); `^src/lbw_config\\.h$` in `.Rbuildignore`
- [ ] Committed with message `build: add stub lbw_config.h with all features disabled` — but WITHOUT `src/lbw_config.h` itself (it is gitignored and configure-generated)

- [ ] **Step 1: Create stub src/lbw_config.h**

  ```bash
  cat > src/lbw_config.h << 'EOF'
  /* Stub — overwritten by configure (Task 5). All features disabled. */
  #ifndef LBW_CONFIG_H
  #define LBW_CONFIG_H
  #define LBW_HAS_OMP_SIMD 0
  #define LBW_HAS_OMP 0
  #define LBW_HAS_GLIBC_MVEC 0
  #endif
  EOF
  ```

- [ ] **Step 2: Add include to ieppa.cpp (after existing #include lines)**

  Add:
  ```cpp
  #include "lbw_config.h"
  ```

- [ ] **Step 3: Add to .Rbuildignore and src/.gitignore BEFORE creating the stub**

  Write the gitignore entry FIRST so the stub file is never accidentally staged:
  ```bash
  echo "^src/lbw_config\\.h$" >> .Rbuildignore
  echo "lbw_config.h" >> src/.gitignore
  ```

  Then create the stub (Step 1 above should run after Step 3's gitignore entries are written):
  ```bash
  cat > src/lbw_config.h << 'EOF'
  /* Stub — overwritten by configure (Task 5). All features disabled. */
  #ifndef LBW_CONFIG_H
  #define LBW_CONFIG_H
  #define LBW_HAS_OMP_SIMD 0
  #define LBW_HAS_OMP 0
  #define LBW_HAS_GLIBC_MVEC 0
  #endif
  EOF
  ```

  Verify the stub is correctly gitignored (git status must NOT show it):
  ```bash
  git status src/lbw_config.h  # expected: "Ignored: src/lbw_config.h" or absent from output
  ```

- [ ] **Step 4: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 5: Commit (WITHOUT src/lbw_config.h — it is gitignored and configure-generated)**

  ```bash
  # NOTE: src/lbw_config.h is intentionally NOT staged here.
  # configure regenerates it at install time from the real probe values (Task 5).
  git add src/ieppa.cpp .Rbuildignore src/.gitignore
  git commit -m "build: add stub lbw_config.h with all features disabled"
  ```

---

### Task 4b: Add #pragma omp simd to normalise loop (ieppa.cpp)

The normalise loop `w[i] /= wm; q[i] /= wm` is element-wise with no scatter — the compiler can auto-vectorise with a hint. The guard falls back to no-op when OpenMP is absent.

**Files:** `src/ieppa.cpp` only — normalise loop at line 108.

**DoD:**
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)`
- [ ] `Rscript -e "testthat::test_local()" 2>&1 | grep -c 'FAIL\|ERROR'` = `0`
- [ ] Committed with message `perf: omp simd hint on normalise loop`

- [ ] **Step 1: Wrap the normalise loop with the simd pragma**

  Replace:
  ```cpp
  for (int i = 0; i < st.n; i++) { w[i] /= wm; q[i] /= wm; }
  ```
  With:
  ```cpp
  // LBW_HAS_OMP_SIMD is numeric (0/1) from lbw_config.h — use #if VALUE not #if defined(VALUE).
  // _OPENMP is a presence macro (undef or defined) — use #if defined(_OPENMP).
  #if defined(_OPENMP) || LBW_HAS_OMP_SIMD
  #pragma omp simd
  #endif
  for (int i = 0; i < st.n; i++) { w[i] /= wm; q[i] /= wm; }
  ```

- [ ] **Step 2: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 3 (GREEN): Test suite**

  ```bash
  Rscript -e "testthat::test_local()" 2>&1 | tail -5
  ```

- [ ] **Step 4: Commit**

  ```bash
  git add src/ieppa.cpp
  git commit -m "perf: omp simd hint on normalise loop"
  ```

---

### Task 4c: Add #pragma omp simd to box-projection loop (ieppa.cpp)

Same SIMD hint for the Dykstra box-projection loop — element-wise clamp with no dependencies.

**Files:** `src/ieppa.cpp` only — box-projection loop at lines 114–119.

**DoD:**
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)`
- [ ] `max(abs(w_after - w_ref)) < 1e-12` (box semantics must be bit-identical)
- [ ] `Rscript -e "testthat::test_local()" 2>&1 | grep -c 'FAIL\|ERROR'` = `0`
- [ ] Committed with message `perf: omp simd hint on box-projection loop`

- [ ] **Step 1: Wrap the box-projection loop**

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

- [ ] **Step 2: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 3 (GREEN): Assert box semantics preserved + test suite**

  ```bash
  Rscript -e "
  library(leafblower)
  set.seed(1); n <- 10000L
  df  <- data.frame(
    age = factor(sample(c('A','B','C'), n, replace=TRUE)),
    sex = factor(sample(c('M','F'), n, replace=TRUE))
  )
  tgt <- list(age=c(A=0.4,B=0.35,C=0.25), sex=c(M=0.5,F=0.5))
  w_after <- harvest(df, tgt, method='ieppa', max_weight=3, attach_weights=FALSE)
  w_ref   <- readRDS('tests/testthat/task2_ieppa_ref.rds')
  cat('max_weight:', max(w_after), '\n')  # must be <= 3
  stopifnot(max(w_after) <= 3.0 + 1e-10)
  cat('PASS\n')
  " && Rscript -e "testthat::test_local()" 2>&1 | tail -5
  ```

- [ ] **Step 4: Commit**

  ```bash
  git add src/ieppa.cpp
  git commit -m "perf: omp simd hint on box-projection loop"
  ```

---

## Phase 2 — Configure-detected libmvec + OpenMP

---

### Task 5: Configure detection — OpenMP, libmvec, lbw_config.h

Extends the existing `configure` script (currently detects C++17 and -O3) with libmvec and OpenMP probes. Overwrites the stub `src/lbw_config.h` from Task 4a with real values. Adds `@OMP_FLAGS@`, `@MAVX2_FLAG@`, and `@MVEC_LIBS@` to `src/Makevars.in`.

**Files:** `configure`, `src/Makevars.in`. (`src/lbw_config.h` is overwritten by configure at install time — not committed.)

**DoD:**
- [ ] `R CMD INSTALL --preclean . 2>&1 | grep -E "(OpenMP|libmvec|DONE)"` shows detection result and `DONE`
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)` on Linux AND macOS paths
- [ ] `Rscript -e "testthat::test_local()" 2>&1 | grep -c 'FAIL\|ERROR'` = `0`
- [ ] Committed with message `build: configure detection for OpenMP, -mavx2, and glibc libmvec`

- [ ] **Step 1: Append OpenMP probe to configure (after the -O3 detection, before the final sed call)**

  ```sh
  # Detect OpenMP: use R's SHLIB_OPENMP_CXXFLAGS (no new package dep).
  OMP_FLAGS=""
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
          LBW_HAS_OMP_SIMD=1
          LBW_HAS_OMP=1
          echo "configure: OpenMP detected — parallel errRp and simd enabled"
      else
          echo "configure: OpenMP not available — falling back to serial"
      fi
      rm -f /tmp/lbw_omp_test.cpp /tmp/lbw_omp_test
  fi
  ```

- [ ] **Step 2: Append libmvec probe (after OpenMP probe, before sed call)**

  ```sh
  # Detect glibc libmvec (_ZGVdN4v_exp — AVX2, 4-wide double, glibc >= 2.22).
  # _ZGVdN4v_exp is in libmvec.so, NOT libm.so. Linker flag must be -lmvec.
  # -mavx2 is required in PKG_CXXFLAGS so __m256d intrinsics compile.
  LBW_HAS_GLIBC_MVEC=0
  MVEC_LIBS=""
  MAVX2_FLAG=""
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
         -o /tmp/lbw_mvec_test -lm -lmvec 2>/dev/null && \
     /tmp/lbw_mvec_test 2>/dev/null; then
      LBW_HAS_GLIBC_MVEC=1
      MVEC_LIBS="-lm -lmvec"
      MAVX2_FLAG="-mavx2"
      echo "configure: glibc libmvec detected — vectorised exp enabled for L-BFGS-B"
  else
      echo "configure: glibc libmvec not available — scalar exp fallback"
  fi
  rm -f /tmp/lbw_mvec_test.cpp /tmp/lbw_mvec_test
  ```

- [ ] **Step 3: Overwrite lbw_config.h with real probe values and update the sed call**

  After the probes, add:
  ```sh
  cat > src/lbw_config.h << EOF
  /* Auto-generated by configure — do not edit. */
  #ifndef LBW_CONFIG_H
  #define LBW_CONFIG_H
  #define LBW_HAS_OMP_SIMD ${LBW_HAS_OMP_SIMD}
  #define LBW_HAS_OMP ${LBW_HAS_OMP}
  #define LBW_HAS_GLIBC_MVEC ${LBW_HAS_GLIBC_MVEC}
  #endif /* LBW_CONFIG_H */
  EOF
  ```

  Replace the existing final sed call (current line 19–20 of configure, 2 substitutions) with 5:
  ```sh
  sed "s|@CXXFLAGS_STD@|${CXXFLAGS_STD}|;s|@OPT_FLAGS@|${OPT_FLAGS}|;s|@OMP_FLAGS@|${OMP_FLAGS}|;s|@MAVX2_FLAG@|${MAVX2_FLAG}|;s|@MVEC_LIBS@|${MVEC_LIBS}|" \
      src/Makevars.in > src/Makevars
  ```

- [ ] **Step 4: Update src/Makevars.in**

  Replace current content:
  ```makefile
  PKG_CXXFLAGS = @CXXFLAGS_STD@ @OPT_FLAGS@ -I. -DSTRICT_R_HEADERS
  PKG_SOURCES = c_api.cpp logit.cpp lbfgsb_solver.cpp ieppa.cpp r_bridge.cpp
  ```
  With:
  ```makefile
  PKG_CXXFLAGS = @CXXFLAGS_STD@ @OPT_FLAGS@ @OMP_FLAGS@ @MAVX2_FLAG@ -I. -DSTRICT_R_HEADERS
  PKG_LIBS = @MVEC_LIBS@
  PKG_SOURCES = c_api.cpp logit.cpp lbfgsb_solver.cpp ieppa.cpp r_bridge.cpp
  ```

- [ ] **Step 5: Add lbw_config.h include to lbfgsb_solver.cpp (after existing includes)**

  ```cpp
  #include "lbw_config.h"
  ```
  (`ieppa.cpp` already has this from Task 4a.)

- [ ] **Step 6: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | grep -E "(OpenMP|libmvec|DONE)"
  ```
  On Linux/gcc: expect "OpenMP detected" and "libmvec detected" and "DONE".
  On macOS: expect both "not available" messages and "DONE".

- [ ] **Step 7: Test suite**

  ```bash
  Rscript -e "testthat::test_local()" 2>&1 | tail -5
  ```

- [ ] **Step 8: Commit**

  ```bash
  git add configure src/Makevars.in src/lbfgsb_solver.cpp
  git commit -m "build: configure detection for OpenMP, -mavx2, and glibc libmvec"
  ```

---

### Task 6: Create lbw_math.hpp + save lbfgsb reference weights

Adds `lbw::bulk_scaled_exp()` — dispatches to `_ZGVdN4v_exp` (glibc libmvec AVX2) when available, scalar fallback otherwise. Also saves reference lbfgsb weights NOW (before the AVX2 path is wired in Task 7) to enable numerical regression testing.

**Files:** Create `src/lbw_math.hpp` only.

**DoD:**
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)` (header not yet included)
- [ ] `tests/testthat/lbfgsb_ref_weights.rds` exists and committed
- [ ] `max(w_ref) <= max_weight` and `mean(w_ref) ≈ 1.0` (confirmed in Step 3)
- [ ] Committed with message `perf: add lbw_math.hpp + save lbfgsb reference weights`

- [ ] **Step 1 (RED): Save lbfgsb reference weights before vectorised path is wired in**

  ```bash
  Rscript -e "
  library(leafblower)
  set.seed(42); n <- 20000L
  df  <- data.frame(
    age = factor(sample(c('18-34','35-54','55+'), n, replace=TRUE)),
    sex = factor(sample(c('M','F'), n, replace=TRUE))
  )
  tgt <- list(age=c('18-34'=0.30,'35-54'=0.45,'55+'=0.25), sex=c(M=0.50,F=0.50))
  w_ref <- harvest(df, tgt, method='lbfgsb', attach_weights=FALSE)
  saveRDS(w_ref, 'tests/testthat/lbfgsb_ref_weights.rds')
  cat('ref saved. max:', max(w_ref), 'mean:', mean(w_ref), '\n')
  stopifnot(abs(mean(w_ref) - 1.0) < 0.01)
  cat('PASS\n')
  "
  ```

- [ ] **Step 2: Write src/lbw_math.hpp**

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
  // With glibc libmvec: 4 doubles/cycle via AVX2 _ZGVdN4v_exp.
  // Otherwise: scalar std::exp.
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

- [ ] **Step 3: Compile gate (header not yet included — tests package builds)**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 4: Commit**

  ```bash
  git add src/lbw_math.hpp tests/testthat/lbfgsb_ref_weights.rds
  git commit -m "perf: add lbw_math.hpp + save lbfgsb reference weights"
  ```

---

### Task 7a: Add F_from_e() and H_from_e() to logit.hpp

These read from a pre-computed `e = exp(logit_scale * u)` array. Used in Tasks 7c/7d where `bulk_scaled_exp` pre-fills the array before the per-obs loop.

**Files:** `src/logit.hpp` only.

**DoD:**
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)`
- [ ] `F_from_e(e) == F(u)` for 3 test values (traced in Step 1)
- [ ] `H_from_e(e, u) == H(u)` for 3 test values (traced in Step 1)
- [ ] `max(abs(w_after - w_ref)) < 1e-12` (functions not yet called — no change expected)
- [ ] Committed with message `perf: add F_from_e/H_from_e to LinkFn for pre-computed exp`

- [ ] **Step 1 (RED): Trace F_from_e and H_from_e on concrete values before implementing**

  For u=0.5, L=0.2, U=5.0, logit_scale=1.5:
  - `e = exp(1.5 * 0.5) = exp(0.75) ≈ 2.117`
  - `F_from_e(2.117)`: `num = 0.2*4 + 5*0.8*2.117 = 9.268`, `denom = 4 + 0.8*2.117 = 5.694` → `f = 9.268/5.694 ≈ 1.628`
  - `F(0.5)` must equal `1.628` — verify against existing F() before implementing F_from_e
  - `H_from_e(2.117, 0.5)`: `num = 4 + 0.8*2.117 = 5.694` → `h = 0.2*0.5 + (4.8/1.5)*ln(5.694/4.8) ≈ 0.395`
  - `H(0.5)` must equal `0.395`

- [ ] **Step 2 (GREEN): Add F_from_e() and H_from_e() after FH() in logit.hpp**

  ```cpp
  // F and H from pre-computed e = exp(logit_scale * u).
  // Precondition: e was computed as exp(logit_scale * u) — not for exponential link.
  // Callers must branch on fn.exponential and fall back to FH() when true.
  double F_from_e(double e) const {
      if (exponential) return e;
      return (L * (U - 1.0) + U * (1.0 - L) * e) /
             ((U - 1.0) + (1.0 - L) * e);
  }
  double H_from_e(double e, double u) const {
      if (exponential) return e;
      double num = (U - 1.0) + (1.0 - L) * e;
      return L * u + (U - L) / logit_scale * std::log(num / (U - L));
  }
  ```

- [ ] **Step 3: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 4 (GREEN): Assert no weight regression (functions unused yet)**

  ```bash
  Rscript -e "
  library(leafblower)
  set.seed(42); n <- 20000L
  df  <- data.frame(
    age = factor(sample(c('18-34','35-54','55+'), n, replace=TRUE)),
    sex = factor(sample(c('M','F'), n, replace=TRUE))
  )
  tgt <- list(age=c('18-34'=0.30,'35-54'=0.45,'55+'=0.25), sex=c(M=0.50,F=0.50))
  w_after <- harvest(df, tgt, method='lbfgsb', attach_weights=FALSE)
  w_ref   <- readRDS('tests/testthat/lbfgsb_ref_weights.rds')
  stopifnot(max(abs(w_after - w_ref)) < 1e-12)
  cat('PASS\n')
  "
  ```

- [ ] **Step 5: Commit**

  ```bash
  git add src/logit.hpp
  git commit -m "perf: add F_from_e/H_from_e to LinkFn for pre-computed exp"
  ```

---

### Task 7b: Add lbw_math.hpp include and e_vec scratch buffer to lbfgsb_solve

Wire `lbw_math.hpp` into `lbfgsb_solver.cpp` and allocate the `e_vec` scratch buffer in `lbfgsb_solve`. No loop bodies changed yet — this is a pure plumbing step.

**Files:** `src/lbfgsb_solver.cpp` only — two additions: one include, one vector declaration.

**DoD:**
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)`
- [ ] `max(abs(w_after - w_ref)) < 1e-12` (e_vec unused — no change expected)
- [ ] Committed with message `perf: wire lbw_math.hpp and e_vec scratch into lbfgsb_solve`

- [ ] **Step 1: Add #include "lbw_math.hpp" to lbfgsb_solver.cpp (after existing includes)**

  ```cpp
  #include "lbw_math.hpp"
  ```

- [ ] **Step 2: Add e_vec declaration in lbfgsb_solve after the existing scratch buffers (~line 306)**

  After `std::vector<double> u_work(st.n);`:
  ```cpp
  std::vector<double> e_vec(st.n);   // scratch: exp(logit_scale * u_work[i]) per trial
  ```

- [ ] **Step 3: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 4 (GREEN): Assert no weight regression**

  ```bash
  Rscript -e "
  library(leafblower)
  set.seed(42); n <- 20000L
  df  <- data.frame(
    age = factor(sample(c('18-34','35-54','55+'), n, replace=TRUE)),
    sex = factor(sample(c('M','F'), n, replace=TRUE))
  )
  tgt <- list(age=c('18-34'=0.30,'35-54'=0.45,'55+'=0.25), sex=c(M=0.50,F=0.50))
  w_after <- harvest(df, tgt, method='lbfgsb', attach_weights=FALSE)
  w_ref   <- readRDS('tests/testthat/lbfgsb_ref_weights.rds')
  stopifnot(max(abs(w_after - w_ref)) < 1e-12)
  cat('PASS\n')
  "
  ```

- [ ] **Step 5: Commit**

  ```bash
  git add src/lbfgsb_solver.cpp
  git commit -m "perf: wire lbw_math.hpp and e_vec scratch into lbfgsb_solve"
  ```

---

### Task 7c: Update wolfe_line_search to use bulk_scaled_exp + e_vec

Adds `e_vec` as a parameter to `wolfe_line_search`, calls `bulk_scaled_exp` before the trial loop, and replaces FH() with F_from_e/H_from_e per observation. Also updates both `wolfe_zoom` call sites in this function to pass `e_vec`.

**Files:** `src/lbfgsb_solver.cpp` only — `wolfe_line_search` function (lines 230–291) and its two `wolfe_zoom` call sites (lines 266, 278). Also update `lbfgsb_solve`'s call to `wolfe_line_search` (line 338).

**DoD:**
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)`
- [ ] `max(abs(w_after - w_ref)) < 1e-10` (AVX2 vs scalar exp may differ by ≤ ULP)
- [ ] `Rscript -e "testthat::test_local()" 2>&1 | grep -c 'FAIL\|ERROR'` = `0`
- [ ] Committed with message `perf: bulk_scaled_exp in wolfe_line_search trial loop`

- [ ] **Step 1: Add e_vec to wolfe_line_search signature (after u_work)**

  Old signature (lines 230–239):
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
  New (insert `std::vector<double>& e_vec,` after `u_work`):
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

- [ ] **Step 2: Replace trial loop body in wolfe_line_search (lines 256–263)**

  Old (lines 256–263 inside the `for (int i = 0; i < 20; i++)` bracket loop):
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

- [ ] **Step 3: Update BOTH wolfe_zoom call sites to pass e_vec**

  Call site 1 (lines 266–269, Armijo-fail or phi-increase):
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

  Call site 2 (lines 278–281, negative slope):
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

- [ ] **Step 4: Update lbfgsb_solve's call to wolfe_line_search (line 338)**

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

- [ ] **Step 5: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 6 (GREEN): Assert weights within AVX2/scalar ULP tolerance**

  ```bash
  Rscript -e "
  library(leafblower)
  set.seed(42); n <- 20000L
  df  <- data.frame(
    age = factor(sample(c('18-34','35-54','55+'), n, replace=TRUE)),
    sex = factor(sample(c('M','F'), n, replace=TRUE))
  )
  tgt <- list(age=c('18-34'=0.30,'35-54'=0.45,'55+'=0.25), sex=c(M=0.50,F=0.50))
  w_after <- harvest(df, tgt, method='lbfgsb', attach_weights=FALSE)
  w_ref   <- readRDS('tests/testthat/lbfgsb_ref_weights.rds')
  delta   <- max(abs(w_after - w_ref))
  cat('delta:', delta, '\n')
  stopifnot(delta < 1e-10)
  cat('PASS\n')
  " && Rscript -e "testthat::test_local()" 2>&1 | tail -5
  ```

- [ ] **Step 7: Commit**

  ```bash
  git add src/lbfgsb_solver.cpp
  git commit -m "perf: bulk_scaled_exp in wolfe_line_search trial loop"
  ```

---

### Task 7d: Update wolfe_zoom to use bulk_scaled_exp + e_vec

Same vectorised exp pattern for the zoom bisection loop. `wolfe_zoom` is called with `u_work` already passed by reference; `e_vec` needs to be added to its signature.

**Files:** `src/lbfgsb_solver.cpp` only — `wolfe_zoom` function (lines 175–225).

**DoD:**
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -1` = `* DONE (leafblower)`
- [ ] `max(abs(w_after - w_ref)) < 1e-10` (same lbfgsb_ref_weights.rds)
- [ ] `Rscript -e "testthat::test_local()" 2>&1 | grep -c 'FAIL\|ERROR'` = `0`
- [ ] Committed with message `perf: bulk_scaled_exp in wolfe_zoom bisection loop`

- [ ] **Step 1: Add e_vec to wolfe_zoom signature (after u_work)**

  Old signature (lines 175–185):
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
  New (insert `std::vector<double>& e_vec,` after `u_work`):
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

- [ ] **Step 2: Replace bisection trial loop body in wolfe_zoom (lines 198–205)**

  Old (inside `for (int j = 0; j < 20; j++)` zoom loop):
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

  Note: the second `u_work` update at line 222 (final accepted-point recompute before
  `phi_from_u`) does NOT need `bulk_scaled_exp` — `phi_from_u` runs its own full O(K·n)
  computation internally.

- [ ] **Step 3: Compile gate**

  ```bash
  R CMD INSTALL --preclean . 2>&1 | tail -1
  ```

- [ ] **Step 4 (GREEN): Assert weights + full test suite**

  ```bash
  Rscript -e "
  library(leafblower)
  set.seed(42); n <- 20000L
  df  <- data.frame(
    age = factor(sample(c('18-34','35-54','55+'), n, replace=TRUE)),
    sex = factor(sample(c('M','F'), n, replace=TRUE))
  )
  tgt <- list(age=c('18-34'=0.30,'35-54'=0.45,'55+'=0.25), sex=c(M=0.50,F=0.50))
  w_after <- harvest(df, tgt, method='lbfgsb', attach_weights=FALSE)
  w_ref   <- readRDS('tests/testthat/lbfgsb_ref_weights.rds')
  delta   <- max(abs(w_after - w_ref))
  cat('delta:', delta, '\n')
  stopifnot(delta < 1e-10)
  cat('PASS\n')
  " && Rscript -e "testthat::test_local()" 2>&1 | tail -5
  ```

- [ ] **Step 5: Commit**

  ```bash
  git add src/lbfgsb_solver.cpp
  git commit -m "perf: bulk_scaled_exp in wolfe_zoom bisection loop"
  ```

---

### Task 8 (removed): OpenMP-parallel compute_errRp

> **Deferred.** After Task 3 reduces errRp frequency to every-10 iterations,
> `compute_errRp` is ~0.7% of total runtime. K=9 OpenMP adds Windows/Rtools
> risk for negligible gain. Tasks 1–7 cover the remaining ~97% of runtime.

---

### Task 8: Final regression, CRAN check, and benchmark

**Files:** `DESCRIPTION` (update SystemRequirements). No source changes.

**DoD:**
- [ ] `R CMD check --as-cran leafblower_*.tar.gz` produces 0 ERRORs, 0 WARNINGs
- [ ] `Rscript -e "testthat::test_local()"` produces 0 failures
- [ ] Benchmark median time documented in commit message
- [ ] Committed with message `perf: benchmark results + CRAN-ready`

- [ ] **Step 1: Update DESCRIPTION SystemRequirements**

  In `DESCRIPTION`, find:
  ```
  SystemRequirements: C++17 compiler
  ```
  Replace with:
  ```
  SystemRequirements: C++17 compiler, OpenMP (optional, for parallel errRp)
  ```

  Commit:
  ```bash
  git add DESCRIPTION
  git commit -m "chore: add OpenMP to SystemRequirements per CRAN policy"
  ```

- [ ] **Step 2: Full test suite (final regression gate)**

  ```bash
  Rscript -e "testthat::test_local()" 2>&1
  ```
  Expected: 0 failures, 0 errors.

- [ ] **Step 3: R CMD check --as-cran**

  ```bash
  R CMD build . && R CMD check --as-cran leafblower_*.tar.gz 2>&1 \
    | grep -E "^(ERROR|WARNING|NOTE|checking)"
  ```
  Expected: 0 ERRORs, 0 WARNINGs. Address any NOTEs before proceeding.

- [ ] **Step 4: L-BFGS-B micro-benchmark (n=50K, 3 margins) — compare to baseline**

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
  bm_after   <- bench::mark(harvest(df, tgt, method="lbfgsb"), iterations=20)
  bm_before  <- readRDS("tests/testthat/lbfgsb_baseline_time.rds")
  pct_faster <- (1 - as.numeric(bm_after$median) / as.numeric(bm_before)) * 100
  cat(sprintf("before: %.0f ms  after: %.0f ms  speedup: %.1f%%\n",
              as.numeric(bm_before)*1000,
              as.numeric(bm_after$median)*1000,
              pct_faster))
  stopifnot(pct_faster >= 20)  # 20% floor; target is ≥25%
  ```
  Target: ≥ 25% reduction on Linux/gcc (libmvec active).

- [ ] **Step 5: Stepstone full-data benchmark (n=1,582,732)**

  ```bash
  python3 benchmarks/stepstone_fulldata_benchmark.py 2>&1
  ```
  Baseline: 127,000 ms (Python), 133,000 ms (R).
  Target: ≥ 15% reduction (≤ 108,000 ms Python).

- [ ] **Step 6: Final commit with benchmark results**

  ```bash
  git add DESCRIPTION
  git commit -m "perf: benchmark results — <X>ms lbfgsb (<Y>% speedup), <Z>ms stepstone"
  ```
