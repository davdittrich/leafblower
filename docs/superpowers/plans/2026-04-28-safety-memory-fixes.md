# Critical Safety & Memory Safety — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate all memory-safety bugs and infrastructure hazards in src/r_bridge.cpp and src/Makevars.

**Architecture:** Four independent fixes: try/catch around solver dispatch; start_weights OOB guard; Makevars gitignored; n*sizeof size_t casts. Each fix is atomic — one commit per task.

**Tech Stack:** C++17, R API (PROTECT/UNPROTECT, Rf_error, LENGTH), git, R CMD INSTALL --preclean .

---

**Mechanism:** R API boundary safety (try/catch, PROTECT counting, LENGTH validation)
**Forbidden:** bare catch(...), Rf_error inside try block, amending prior commits
**Audit:** devtools::test() before and after each task; R CMD INSTALL --preclean . compile gate

---

## Baseline

Before starting any task, record the baseline test count:

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

All tasks must leave this count unchanged or reduced.

---

## Task 1: B4 — Wrap solver dispatch in try/catch

**Files:**
- Modify: `src/r_bridge.cpp:366-520`
- Test: `tests/testthat/test-safety.R` (new file, shared across tasks)

**Context:** `C_rk_calibrate` makes zero PROTECT calls before the solver dispatch block at line 366. The first PROTECT in the function is `wts` at line 533, which is after the dispatch. Therefore N=0 — no UNPROTECT is needed in the catch block for PROTECTs made before the try. The try block must not contain any `Rf_error` calls (they longjmp); call `Rf_error` only after the try/catch in the catch handler.

- [ ] **Step 1: Write the static-analysis check + regression test**

**Note on RED/GREEN:** B4 is defensive hardening with no observable behavioral change for
correct inputs — triggering `std::bad_alloc` in a test is not feasible. The strategy is:
(a) a grep-based static assertion that `catch` wraps the dispatch block, and
(b) a regression guard test that confirms no crash on valid inputs. There is no RED phase.

Create `tests/testthat/test-safety.R`:

```r
# tests/testthat/test-safety.R
# Safety & memory-safety regression tests (B4, B7, R6)

library(leafblower)

# Helper: minimal valid harvest() call
.mini_harvest <- function(n = 20, method = "raking") {
  set.seed(1)
  df <- data.frame(
    sex  = sample(c("M", "F"), n, replace = TRUE),
    wt   = rep(1.0, n)
  )
  targets <- list(sex = c(M = 0.5, F = 0.5))
  harvest(df, targets = targets, weight_col = "wt", method = method)
}

# B4: solver dispatch must not crash R on a valid call (regression guard)
test_that("B4: solver dispatch returns without crashing for each method", {
  for (m in c("raking", "lbfgsb", "ieppa", "auto")) {
    expect_no_error(.mini_harvest(method = m),
                    label = paste("method =", m))
  }
})
```

- [ ] **Step 1b: Verify try/catch presence via grep (static assertion)**

```bash
grep -c "catch (const std::exception" src/r_bridge.cpp
```

Expected: `1` (exactly one catch block wrapping the dispatch). If the output is `0`, the
fix in Step 3 was not applied — halt and re-check.

- [ ] **Step 2: Run test to verify it passes (baseline)**

```bash
Rscript -e "devtools::test(filter='safety')"
```

Expected: PASS (baseline — confirms the test harness works before the C++ change)

- [ ] **Step 3: Implement the fix**

In `src/r_bridge.cpp`, wrap lines 366–520 (the entire solver dispatch block, from the first `if (strcmp(method_str, "lbfgsb")` through the closing `}` at line 520) in a try/catch. Replace:

```cpp
    if (strcmp(method_str, "lbfgsb") == 0) {
        // ... all dispatch branches ...
    }
```

with:

```cpp
    try {
        if (strcmp(method_str, "lbfgsb") == 0) {
            auto res = lbw::lbfgsb_solve(st);
            res_status     = res.status;
            res_iterations = res.iterations;
            res_max_error  = res.max_error;
            res_alg_used   = (int)RK_ALG_LBFGSB;
            pack_solver_result(res);
            res_best_weights = std::move(res.best_weights);
        } else if (strcmp(method_str, "raking") == 0) {
            auto res = lbw::raking_solve(st);
            res_status     = res.status;
            res_iterations = res.iterations;
            res_max_error  = res.max_error;
            res_alg_used   = (int)RK_ALG_RAKING;
            pack_solver_result(res);
            res_best_weights = std::move(res.best_weights);
        } else if (strcmp(method_str, "auto") == 0) {
            // ... keep existing auto block verbatim ...
        } else if (strcmp(method_str, "sinkhorn") == 0) {
            // ... keep existing sinkhorn block verbatim ...
        } else if (strcmp(method_str, "greg") == 0) {
            // ... keep existing greg block verbatim ...
        } else {
            // ... keep existing else block verbatim ...
        }
    } catch (const std::exception& e) {
        // N=0: no PROTECT calls were made before this try block.
        Rf_error("leafblower: internal solver error — %s", e.what());
    }
```

**Critical constraint:** `Rf_error` is called in the catch handler, outside the try block. This is correct — `Rf_error` longjmps and must not be inside a try block whose stack frames hold C++ objects with destructors.

- [ ] **Step 4: Compile**

```bash
R CMD INSTALL --preclean .
```

Expected: no errors, no warnings about exception specifications.

- [ ] **Step 5: Run tests**

```bash
Rscript -e "devtools::test(filter='safety')"
```

Expected: PASS

```bash
Rscript -e "devtools::test()"
```

Expected: same pass count as baseline.

- [ ] **Step 6: Commit**

```bash
git add src/r_bridge.cpp tests/testthat/test-safety.R
git commit -m "$(cat <<'EOF'
fix(r_bridge): wrap solver dispatch in try/catch to prevent C++ exception escape

C++ exceptions thrown inside lbw::*_solve() would bypass R's longjmp-based
error handling and terminate the R session. The catch handler calls Rf_error
outside the try block (required: Rf_error longjmps, cannot be in a try that
owns C++ RAII objects). N=0 PROTECTs before the dispatch, so no UNPROTECT
needed in the catch path. Fixes B4.
EOF
)"
```

---

## Task 2: B7 — Validate start_weights length

**Files:**
- Modify: `src/r_bridge.cpp:177-182`
- Test: `tests/testthat/test-safety.R` (append)

**Context:** Lines 177–182 unconditionally call `REAL(start_weights_sexp)` and iterate `n` times. If the caller passes a vector shorter than `n`, this reads past the end of an R-allocated buffer — undefined behavior. The guard must be inserted after the `Rf_isNull` check and before `REAL(...)` is called.

**R-layer note:** `R/harvest.R:normalize_start_weights` (line 536–548) already validates
`length(start_weights) != n` and calls `stop("start_weights length must equal nrow(data)")`.
This means `harvest()` intercepts the bad input BEFORE reaching `.Call`. The C-layer guard is
still necessary (defense-in-depth for any caller that bypasses R), but the test MUST use
`.Call("C_rk_calibrate", ...)` directly to bypass the R wrapper and reach the C guard. The
test verifies the C-layer check, not the R-layer check.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-safety.R`:

```r
# B7: start_weights of wrong length must error at the C layer
# (uses .Call directly to bypass R/harvest.R's normalize_start_weights check)
test_that("B7: C_rk_calibrate start_weights length mismatch produces error", {
  set.seed(2)
  n <- 30L
  df <- data.frame(
    sex = sample(c("M", "F"), n, replace = TRUE),
    wt  = rep(1.0, n)
  )
  targets <- list(sex = c(M = 0.5, F = 0.5))
  bad_sw  <- rep(1.0, n - 5L)   # wrong length: 25 instead of 30

  # All 31 args to C_rk_calibrate, matching the .Call in R/harvest.R lines 241-276.
  # Arg 8 (sw_vec/start_weights) is bad_sw — wrong length — to trigger the C-layer guard.
  # Non-critical args use safe scalar defaults; the error must fire before the solver runs.
  expect_error(
    .Call("C_rk_calibrate",
          df,                          # 1:  data
          targets,                     # 2:  target
          as.double(0),                # 3:  min_weight
          as.double(1e6),              # 4:  max_weight
          as.character("raking"),      # 5:  method
          as.integer(0L),              # 6:  verbose
          as.integer(50L),             # 7:  max_iterations
          bad_sw,                      # 8:  start_weights (wrong length — the probe)
          as.double(1e-6),             # 9:  tol_abs (legacy slot)
          as.integer(0L),              # 10: bounds_mode_int (0 = "cell")
          as.integer(1L),              # 11: homotopy_levels
          as.double(1.0),              # 12: homotopy_start_factor
          as.double(1.0),              # 13: homotopy_end_factor
          as.double(1.0),              # 14: homotopy_budget_p
          as.character("round_robin"), # 15: scheduler
          as.character("none"),        # 16: eta_schedule
          as.double(0.5),              # 17: eta_start
          as.double(0.1),              # 18: eta_end
          as.double(1.0),              # 19: eta_schedule_power
          as.double(1e-4),             # 20: conv pct_tol
          as.double(0.0),              # 21: conv absolute_tol
          as.integer(4L),              # 22: conv metric_int (4 = grake_norm)
          as.integer(0L),              # 23: conv rule_int (0 = threshold)
          as.integer(0L),              # 24: conv stop_when_int (0 = any)
          as.integer(0L),              # 25: sor enabled
          as.integer(0L),              # 26: sor auto
          as.double(1.0),              # 27: sor omega_init
          as.double(0.5),              # 28: sor omega_min
          as.double(1.0),              # 29: sor omega_fixed
          as.integer(0L),              # 30: sor burnin
          as.integer(0L)               # 31: accelerate_bool
    ),
    regexp = "start_weights length"
  )
})

- [ ] **Step 2: Run test to verify it fails**

```bash
Rscript -e "devtools::test(filter='safety')"
```

Expected: the B7 test FAILS (no LENGTH check in C layer yet — `.Call` either silently reads
garbage or crashes before producing the expected error message).

- [ ] **Step 3: Implement the fix**

In `src/r_bridge.cpp`, replace lines 177–182:

```cpp
    if (Rf_isNull(start_weights_sexp)) {
        for (int i = 0; i < n; i++) weights[i] = 1.0;
    } else {
        const double* sw = REAL(start_weights_sexp);
        for (int i = 0; i < n; i++) weights[i] = sw[i];
    }
```

with:

```cpp
    if (Rf_isNull(start_weights_sexp)) {
        for (int i = 0; i < n; i++) weights[i] = 1.0;
    } else {
        if (LENGTH(start_weights_sexp) != n)
            Rf_error("leafblower: start_weights length %d != n=%d",
                     (int)LENGTH(start_weights_sexp), n);
        const double* sw = REAL(start_weights_sexp);
        for (int i = 0; i < n; i++) weights[i] = sw[i];
    }
```

- [ ] **Step 4: Compile**

```bash
R CMD INSTALL --preclean .
```

Expected: no errors.

- [ ] **Step 5: Run tests**

```bash
Rscript -e "devtools::test(filter='safety')"
```

Expected: B7 test now PASSES.

```bash
Rscript -e "devtools::test()"
```

Expected: same pass count as baseline.

- [ ] **Step 6: Commit**

```bash
git add src/r_bridge.cpp tests/testthat/test-safety.R
git commit -m "$(cat <<'EOF'
fix(r_bridge): guard start_weights length before REAL() dereference

Reading REAL(start_weights_sexp)[i] for i >= LENGTH(...) is UB. Added
LENGTH check before the copy loop and Rf_error on mismatch. Fixes B7.
EOF
)"
```

---

## Task 3: B8 — Remove src/Makevars from VCS

**Files:**
- Remove from VCS: `src/Makevars`
- Modify: `src/.gitignore`
- Verify: `src/Makevars.in` already has `@MVEC_LIBS@` placeholder (confirmed: line 2)
- Verify: `configure` line 103–104 runs `sed ... src/Makevars.in > src/Makevars`

**Context:** `src/Makevars` is a build artifact generated by `configure` from `src/Makevars.in`. Tracking it in git causes spurious diffs and means developers on systems without libmvec get the wrong flags committed. `configure` runs automatically during `R CMD INSTALL`. `src/.gitignore` currently only ignores `lbw_config.h`; `Makevars` must be added.

No test file changes — this task is infrastructure only.

- [ ] **Step 1: Verify the current tracked state**

```bash
git ls-files src/Makevars src/Makevars.bak
```

Expected: both `src/Makevars` and `src/Makevars.bak` are printed (both tracked).

- [ ] **Step 2: Remove both from VCS (keep files on disk)**

```bash
git rm --cached src/Makevars src/Makevars.bak
```

Expected:
```
rm 'src/Makevars'
rm 'src/Makevars.bak'
```

- [ ] **Step 3: Add both to .gitignore**

Current `src/.gitignore` content: `lbw_config.h`

Edit `src/.gitignore` — append two new lines:

```
Makevars
Makevars.bak
```

Final `src/.gitignore`:
```
lbw_config.h
Makevars
Makevars.bak
```

- [ ] **Step 4: Verify configure regenerates Makevars**

```bash
R CMD INSTALL --preclean .
```

Expected: configure runs, produces `src/Makevars`, compilation succeeds. The file must exist after install:

```bash
ls -la src/Makevars
```

Expected: file present with today's timestamp.

- [ ] **Step 5: Verify git no longer tracks it**

```bash
git status src/Makevars
```

Expected: either not listed, or listed under "Untracked files" — NOT under "Changes to be committed" or "Changes not staged".

- [ ] **Step 6: Run tests**

```bash
Rscript -e "devtools::test()"
```

Expected: same pass count as baseline.

- [ ] **Step 7: Commit**

```bash
git add src/.gitignore
git commit -m "$(cat <<'EOF'
fix(build): remove src/Makevars from VCS; it is a configure-generated artifact

src/Makevars is produced by configure via sed substitution from Makevars.in.
Committing it causes spurious diffs when building on systems with or without
libmvec. Added to src/.gitignore. R CMD INSTALL --preclean . regenerates it
automatically. Fixes B8.
EOF
)"
```

---

## Task 4: R6 — Cast n to size_t in memcpy calls

**Files:**
- Modify: `src/r_bridge.cpp:534, 653, 656`
- Test: `tests/testthat/test-safety.R` (no new tests needed — compile gate + existing tests suffice)

**Context:** Three `memcpy`/`std::memcpy` calls multiply `n` (type `int`) by `sizeof(T)`. On platforms where `int` is 32-bit and `size_t` is 64-bit, `n * sizeof(T)` promotes to `ptrdiff_t` (signed 64-bit) rather than `size_t`, which is harmless at small n but is technically implementation-defined for the multiplication. Casting `n` to `size_t` before the multiply makes the arithmetic well-defined and suppresses `-Wsign-conversion` warnings. Confirmed locations:

- Line 534: `memcpy(REAL(wts), weights.data(), (size_t)n * sizeof(double));` — already correct (cast present)
- Line 653: `std::memcpy(INTEGER(cell_of_sexp), ct.cell_of.data(), n * sizeof(int));` — missing cast
- Line 656: `std::memcpy(INTEGER(npc), ct.n_per_cell.data(), ct.M_cell * sizeof(int));` — `ct.M_cell` is int, missing cast

Grep all memcpy calls to confirm no others exist:

```bash
grep -n "memcpy" src/r_bridge.cpp
```

- [ ] **Step 1: Grep to find all memcpy sites**

```bash
grep -n "memcpy" /home/dd/Gemini/leafblower/src/r_bridge.cpp
```

Expected output (3 lines):
```
534:    memcpy(REAL(wts), weights.data(), (size_t)n * sizeof(double));
653:    std::memcpy(INTEGER(cell_of_sexp), ct.cell_of.data(), n * sizeof(int));
656:    std::memcpy(INTEGER(npc), ct.n_per_cell.data(), ct.M_cell * sizeof(int));
```

Line 534 already has `(size_t)n`. Fix only lines 653 and 656.

- [ ] **Step 2: Implement the fix**

In `src/r_bridge.cpp`, replace line 653:

```cpp
    std::memcpy(INTEGER(cell_of_sexp), ct.cell_of.data(), n * sizeof(int));
```

with:

```cpp
    std::memcpy(INTEGER(cell_of_sexp), ct.cell_of.data(), (size_t)n * sizeof(int));
```

Replace line 656:

```cpp
    std::memcpy(INTEGER(npc), ct.n_per_cell.data(), ct.M_cell * sizeof(int));
```

with:

```cpp
    std::memcpy(INTEGER(npc), ct.n_per_cell.data(), (size_t)ct.M_cell * sizeof(int));
```

- [ ] **Step 3: Compile**

```bash
R CMD INSTALL --preclean .
```

Expected: no errors, no `-Wsign-conversion` warnings on the patched lines.

- [ ] **Step 4: Run tests**

```bash
Rscript -e "devtools::test(filter='cell-table')"
```

Expected: PASS (cell table probe exercises the patched memcpy sites).

```bash
Rscript -e "devtools::test()"
```

Expected: same pass count as baseline.

- [ ] **Step 5: Commit**

```bash
git add src/r_bridge.cpp
git commit -m "$(cat <<'EOF'
fix(r_bridge): cast n/M_cell to size_t in memcpy size arguments

int * sizeof(T) is technically implementation-defined for the size_t
promotion. Cast to (size_t) before multiply on lines 653 and 656 to
match the existing pattern at line 534. Fixes R6.
EOF
)"
```

---

---

## Task 5: R6-py — Cast n to size_t in Python bindings memcpy

**Files:**
- Modify: `python/leafblower/_bindings.cpp:131`

**Context:** `python/leafblower/_bindings.cpp` line 131 contains:
```cpp
std::memcpy(weights_out.mutable_data(), weights_copy.data(), n * sizeof(double));
```
`n` is `int`; multiplying by `sizeof(double)` (type `size_t`) technically promotes through
`ptrdiff_t` (signed) rather than `size_t` on 64-bit platforms — implementation-defined for the
promotion. Cast `n` to `size_t` for well-defined arithmetic and to suppress `-Wsign-conversion`.

Grep all memcpy calls in the Python bindings to confirm no other sites:

```bash
grep -n "memcpy" python/leafblower/_bindings.cpp
```

- [ ] **Step 1: Grep to find all memcpy sites**

```bash
grep -n "memcpy" python/leafblower/_bindings.cpp
```

Expected output (1 line):
```
131:            std::memcpy(weights_out.mutable_data(), weights_copy.data(), n * sizeof(double));
```

- [ ] **Step 2: Implement the fix**

In `python/leafblower/_bindings.cpp`, replace line 131:

```cpp
            std::memcpy(weights_out.mutable_data(), weights_copy.data(), n * sizeof(double));
```

with:

```cpp
            std::memcpy(weights_out.mutable_data(), weights_copy.data(), (size_t)n * sizeof(double));
```

- [ ] **Step 3: Compile**

```bash
cd python && pip install -e . --no-build-isolation 2>&1 | tail -10
```

Expected: no errors, no `-Wsign-conversion` warning on line 131.

- [ ] **Step 4: Run Python tests**

```bash
cd python && python -m pytest python/leafblower/test_python.py -v 2>&1 | tail -20
```

Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add python/leafblower/_bindings.cpp
git commit -m "$(cat <<'EOF'
fix(python): cast n to size_t in _bindings.cpp memcpy size argument

int * sizeof(double) is implementation-defined for the size_t promotion.
Cast to (size_t) before multiply at line 131 to match the R-bridge pattern.
Fixes R6-py.
EOF
)"
```

---

## Completion Checklist

- [ ] All 5 tasks committed (5 separate commits, no bundling)
- [ ] `git log --oneline -5` shows 5 fix commits
- [ ] `Rscript -e "devtools::test()"` pass count equals or exceeds baseline
- [ ] `git ls-files src/Makevars` returns empty (not tracked)
- [ ] `git ls-files src/Makevars.bak` returns empty (not tracked)
- [ ] `src/.gitignore` contains `Makevars` and `Makevars.bak`
- [ ] `tests/testthat/test-safety.R` exists with B4 and B7 tests
- [ ] `grep -c "catch (const std::exception" src/r_bridge.cpp` returns `1`
- [ ] `grep -n "memcpy" python/leafblower/_bindings.cpp` shows `(size_t)n` cast
