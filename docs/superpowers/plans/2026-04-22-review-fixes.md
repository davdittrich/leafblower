# Plan: Post-P3 Code Review Fixes

**Date:** 2026-04-22  
**Scope:** 12 items from critical-code-reviewer of P1–P3 calibration refactor  
**Files touched:** `src/lbfgsb_solver.cpp`, `src/logit.hpp`, `src/ieppa.cpp`,
`src/types.hpp`, `R/harvest.R`, `tests/testthat/test-lbfgsb.R`,
`python/leafblower/_harvest.py`

---

## Work Units

### WU-1 — Fix misleading "Do NOT normalize" comment (Review item 1, Blocking)
**File:** `src/lbfgsb_solver.cpp`  
**Location:** Inside `compute_final_weights_and_error`, comment at lines ~169–171:
```cpp
// Do NOT normalize: bridge normalizes start_weights to mean=1 before
// rk_calibrate(); re-normalizing after clamping invalidates calibration constraints.
```
**Problem:** `harvest.R:99` executes `weights <- weights / mean(weights)` unconditionally after `rk_calibrate()`. The C++ comment claims this "invalidates calibration constraints" — false. Constraints are proportional so scaling preserves them. Contradicts the R code; will mislead future maintainers.  
**Fix:** Replace with:
```cpp
// C returns raw weights. Caller (harvest.R) normalises to mean=1 post-call;
// safe because calibration constraints are proportional (scaling w preserves margins).
```

---

### WU-2 — Fix overconfident "guarantees sum(w)=n" comment (Review item 7, Required)
**File:** `src/lbfgsb_solver.cpp`  
**Location:** Inside `lbfgsb_solve`:
```cpp
// ALM fields inactive: dual calibration guarantees sum(w)=n at convergence.
```
**Problem:** The guarantee has unchecked preconditions: targets must sum to exactly 1 per margin, no observation may have NA across all calibration variables. The actual safety net is `harvest.R`'s normalization.  
**Fix:**
```cpp
// ALM inactive. sum(w)≈n holds at convergence when targets sum to 1 per margin;
// harvest.R normalises post-call as a safety net for non-converged iterates.
```

---

### WU-3 — Fix `dF` recomputing `logit_scale` inline (Review item 4, Required)
**File:** `src/logit.hpp`  
**Location:** `LinkFn::dF` body, the local `ls` variable:
```cpp
double ls = (U - L) / ((U - 1.0) * (1.0 - L));
return ls * (f - L) * (U - f) / (U - L);
```
**Problem:** Struct member `logit_scale` stores this identical formula. Two independent copies will diverge silently if the formula is updated.  
**Fix:** Delete the local `ls` declaration; use the field:
```cpp
return logit_scale * (f - L) * (U - f) / (U - L);
```

---

### WU-4 — Fix `box_ok` wrong comparator in iEPPA post-loop finalizer (Review item 5, Required)
**File:** `src/ieppa.cpp`  
**Location:** Post-loop Dykstra finalizer; the `box_ok` check inside `for (int fixup ...)`:
```cpp
if (wc != w[i]) box_ok = false;
```
**Problem:** `wc = clamp(yi, lo, hi)` where `yi = w[i] + q[i]`. When no clamping occurs, `wc = yi != w[i]` whenever `q[i] != 0`. This fires `box_ok = false` during residual-correction drain — extra unnecessary finalizer iterations. Comment "terminates in 1 iteration at true convergence" is false under this condition.  
**Fix:** Detect actual clamp events:
```cpp
if (yi != wc) box_ok = false;
```
`yi` is already declared in the same loop body, two lines above this check.  
**Verification:** The `box_ok` fix is performance-only — correctness (box constraints satisfied on output) is already verified by `tests/testthat/test-ieppa.R:55–56`:
```r
expect_true(max(res) <= 2.0 + 1e-10)
expect_true(min(res) >= 0.2 - 1e-10)
```
These bounds assertions fail if the finalizer terminates with out-of-range weights. No new test required.

---

### WU-5 — Rename `cres` to `calib_result` in `harvest.R` (Review item 6, Required)
**File:** `R/harvest.R`  
**Location:** Line 87 and 4 downstream references (`cres$status`, `cres$max_error`, `cres$algorithm_used`, `cres$message`):
```r
cres <- raw$result
```
**Fix:**
```r
calib_result <- raw$result
```
Update all 4 downstream references. No downstream files import this variable.

---

### WU-6 — Fix `alg_names` dead "auto" slot (Review item 10, Suggestion)
**File:** `R/harvest.R`  
**Location:** Line ~106:
```r
alg_names <- c("auto", "ieppa", "lbfgsb")  # index 0 (auto) unreachable after match.arg; kept for 1-indexed enum alignment
```
**Problem:** If C returns `algorithm_used=0` via a non-R path, `alg_names[0+1L]` silently mislabels the algorithm as "auto".  
**Fix:**
```r
alg_names <- c("", "ieppa", "lbfgsb")  # index 0 (auto) removed from user API
```
`""` is immediately visible in output; `"auto"` would falsely imply success.

---

### WU-7 — Rename 3 misleading ALM test names (Review item 2, Blocking)
**File:** `tests/testthat/test-lbfgsb.R`  
**Location:** `test_that(...)` descriptions at lines 52, 64, 77  
**Problem:** ALM outer loop was removed. These tests verify L-BFGS-B bound/mean properties with no ALM involvement.  
**Fix:**
- Line 52: `"L-BFGS-B converges with tight bounds (max=1.5, min=0.2)"`
- Line 64: `"L-BFGS-B stable near infeasibility boundary (90/10 split, tight bounds)"`
- Line 77: `"L-BFGS-B mean=1 after normalization, loose bounds"`

---

### WU-8 — Replace vacuous mean=1 test with calibration quality check (Review item 3, Blocking)
**File:** `tests/testthat/test-lbfgsb.R`  
**Location:** Line 16 inside the "L-BFGS-B converges on 3-margin no-bounds case" test  
**Current:**
```r
expect_lt(abs(mean(result$weights) - 1.0), 1e-8)
```
**Problem:** `harvest.R` always normalises to mean=1. This assertion is true by construction; tests the normalization arithmetic, not L-BFGS-B convergence quality.  
**Fix:** Replace with calibration error check. The test at line 1–17 has `df` (data frame with `age`, `sex`, `edu` columns) and `tgt` (calibration targets) already in scope — verified at lines 7–12:
```r
diag <- diagnose_weights(df, tgt, result$weights)
expect_true(all(abs(diag$error_weighted) < 1e-4))
```
`diagnose_weights` is exported from leafblower (`NAMESPACE` export confirmed). `df` is the original data frame (not `result`). Tolerance 1e-4 is appropriate for n=50k with 3 margins.

---

### WU-9 — Add normalization parity to Python `_harvest.py` (Review item 8, Required)
**File:** `python/leafblower/_harvest.py`  
**Location:** After the status-check block (lines ~157–165), before `if not attach_weights`  
**Problem:** R always does `weights / mean(weights)`. Python returns raw C output. Non-converged iterates have `mean != 1` — behavioral divergence.  
**Fix:**
```python
# Match R behaviour: normalise to mean=1 (preserves proportional constraints).
w_mean = weights_out.mean()
if w_mean > 0:
    weights_out = weights_out / w_mean
```

---

### WU-10 — Fix Python silent failure on bad input type (Review item 9, Required)
**File:** `python/leafblower/_harvest.py`  
**Location:** Lines ~73–74 (current guard), replacing with the following structure after the dict-to-DataFrame conversion block  
**Current:**
```python
if _PANDAS_AVAILABLE and not isinstance(data, pd.DataFrame):
    raise TypeError("data must be a pd.DataFrame or dict")
```
**Problem:** When `_PANDAS_AVAILABLE = False`, the condition short-circuits to `False` for all inputs. A numpy array or list passes silently and hits `data[varname]` with an opaque `AttributeError`.  
**Fix:** Replace the entire validation block (the single `if` guard at line ~73) with:
```python
if isinstance(data, dict):
    if not _PANDAS_AVAILABLE:
        raise ImportError("pandas required to use dict input; install with pip install pandas")
    data = pd.DataFrame(data)
elif _PANDAS_AVAILABLE:
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pd.DataFrame or dict")
else:
    raise TypeError("data must be a dict when pandas is not installed")
```
This replaces both the existing dict-to-DataFrame block (lines ~68–71) AND the existing type guard (line ~73) as a single consolidated structure. Remove the original two separate blocks.

---

### WU-11 — Simplify ALM guard to single `alm_mu > 0.0` check (Review item 12, Suggestion)
**Files:** `src/lbfgsb_solver.cpp` (4 guard sites), `src/types.hpp` (comment on `alm_lambda`)  
**Current at each of 4 sites:**
```cpp
if (st.alm_mu > 0.0 || st.alm_lambda != 0.0)
```
**Problem:** `alm_lambda != 0.0` is a float equality comparison. `alm_lambda` is only meaningful when `alm_mu > 0` — allowing lambda active with mu=0 is an undocumented half-state.  
**Fix in `lbfgsb_solver.cpp`:** Change all 4 occurrences to `if (st.alm_mu > 0.0)`.  
**Fix in `src/types.hpp`:** Add to the `alm_lambda` comment:
```cpp
double alm_lambda = 0.0;  // dual variable for sum(w)=n; only read when alm_mu > 0
double alm_mu     = 0.0;  // penalty coefficient; 0.0 = ALM inactive
```

---

### WU-12 — Document `lbfgsb_solve_inner` extraction purpose (Review item 11, Suggestion)
**File:** `src/lbfgsb_solver.cpp`  
**Location:** Above the `static LBFGSResult lbfgsb_solve_inner(...)` declaration  
**Review finding:** The extraction has no current justification — function has one caller and was scaffolded for a removed ALM outer loop. The review explicitly offered two resolutions: document OR inline.  
**Chosen resolution:** Document (not inline). Inlining 80 lines into `lbfgsb_solve` produces an unreadable monolith with no separation between setup (offsets, targets, design weights) and optimization (L-BFGS-B iteration). The extraction is clean and keeps future ALM work feasible.  
**Fix:** Add above the function:
```cpp
// Extracted from lbfgsb_solve to isolate the L-BFGS-B iteration kernel.
// If an ALM outer loop is added, lbfgsb_solve calls this per outer iteration.
```

---

## Order of Execution

1. **WU-1, WU-2, WU-12** simultaneously — comment-only C++ changes, no recompile needed
2. **WU-3** — `logit.hpp`, recompile gate after
3. **WU-4** — `ieppa.cpp`, recompile gate after
4. **WU-11** — `lbfgsb_solver.cpp` (4 sites) + `types.hpp`, recompile gate after
5. **WU-5, WU-6** simultaneously — R naming, no compile needed
6. **WU-7, WU-8** simultaneously — test updates, no compile needed
7. **WU-9, WU-10** simultaneously — Python, no compile needed
8. `R CMD INSTALL --preclean .` — full compile gate
9. `devtools::test()` — must show `FAIL 0`
10. `python -c "from leafblower import harvest; print('ok')"` — Python smoke test
11. Single commit

---

## Verification Checklist
- [ ] `R CMD INSTALL --preclean .` succeeds, zero new compiler warnings
- [ ] `devtools::test()` FAIL 0
- [ ] Python smoke test passes
- [ ] Review item cross-reference (WU → review item mapping):
  - Blocking: items 1→WU-1, 2→WU-7, 3→WU-8
  - Required: items 4→WU-3, 5→WU-4, 6→WU-5, 7→WU-2, 8→WU-9, 9→WU-10
  - Suggestions: items 10→WU-6, 11→WU-12, 12→WU-11
- [ ] No functional change to calibration algorithms (comment/test/name changes only in C++/R algo files)
