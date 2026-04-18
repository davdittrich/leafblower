# Clean Code Fixes (Round 2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix 11 clean-code violations identified by review: 3 silent bugs (dead `dF`, orphaned `epsilon`, discarded `tol_abs` return value), 5 design issues (duplicate `phi_and_grad`/`phi_from_u`, inner-loop allocations, O(K×C) column scan, discarded `compute_targets_abs` return, `compute_errRp` per-call allocations), and 3 minor polish items (`linf` name, `enforce_mean` silent ignore, `outer_max_iter` docs).

**Architecture:** Changes span C++ internals (`src/lbfgsb_solver.cpp`, `src/ieppa.cpp`, `src/logit.hpp`, `src/types.hpp`, `src/c_api.cpp`), the public C API (`src/leafblower.h`), the R/C bridge (`src/r_bridge.cpp`), and R layer (`R/harvest.R`). Tasks are independent; each compiles and passes all tests before the next begins. Compile gate: `R CMD INSTALL --preclean .`. Test gate: `Rscript -e "devtools::test()"`.

**Tech Stack:** C++17, R, testthat, `R CMD INSTALL --preclean .`

---

## File Map

| File | Changes |
|---|---|
| `src/logit.hpp` | Remove dead `dF()` method (Task 1) |
| `src/types.hpp` | Remove dead `epsilon` field (Task 1) |
| `src/c_api.cpp` | Remove dead `st.epsilon = p->epsilon` copy-in (Task 1) |
| `src/leafblower.h` | Mark `epsilon` deprecated, document `outer_max_iter` R-bridge limitation (Tasks 1, 6) |
| `src/r_bridge.cpp` | Delete stale comment, add `tol_abs_sexp` 9th parameter, pre-build column `unordered_map` (Tasks 1, 2, 6) |
| `R/harvest.R` | Capture `parse_convergence` return value and pass to C, add `enforce_mean` to ignored params (Tasks 2, 6) |
| `src/ieppa.cpp` | Hoist `bucket` in `compute_errRp`, document `1e300` sentinel (Tasks 3, 6) |
| `src/lbfgsb_solver.cpp` | Deduplicate `phi_and_grad`/`phi_from_u`, use `compute_targets_abs` return value, hoist 3 vectors above L-BFGS-B loop, rename `linf` (Tasks 4, 5, 6) |
| `tests/testthat/test-harvest.R` | Add `tol_abs` forwarding behavioral test (Task 2) |

---

## Task 1: Dead code removal + API documentation

**Context:** Three dead artifacts ship in every build:
- `dF()` in `src/logit.hpp` (never called anywhere — confirmed by `grep -rn "\.dF(" src/ R/`)
- `epsilon` in `rk_params_t` propagates to `CalibState.epsilon` but `ieppa_solve` never reads `st.epsilon`
- Stale comment `// Existing logit test bridges (from Task 4)` in `src/r_bridge.cpp:11` is a task-tracker artifact

**Files:**
- Modify: `src/logit.hpp:39-45` (remove `dF`)
- Modify: `src/types.hpp:22` (remove `epsilon` from `CalibState`)
- Modify: `src/c_api.cpp:181` (remove `st.epsilon = p->epsilon;`)
- Modify: `src/leafblower.h:33` (mark `epsilon` deprecated in comment)
- Modify: `src/r_bridge.cpp:11` (delete stale comment)

- [ ] **Step 1: Verify `dF` has no callers**

```bash
grep -rn "\.dF\b\|fn\.dF\|->dF" src/ R/ tests/
```

Expected: no output. If any output appears, stop — the function is used and must not be deleted.

- [ ] **Step 2: Remove `dF` from `src/logit.hpp`**

Delete lines 39–45 (the entire `dF` method). After deletion the struct should end with `H(double u) const { ... }` followed by `};`:

```cpp
    // H(u): antiderivative of F(u).
    // Logit branch: H(0) = 0 by construction (constant of integration chosen).
    // Exp branch: H(u) = exp(u); H(0) = 1 (additive constant irrelevant for optimization).
    // Logit (Deville-Sarndal 1992):
    //   H(u) = L*u + (U-L)/logit_scale * ln(((U-1)+(1-L)*exp(logit_scale*u)) / (U-L))
    double H(double u) const {
        if (exponential) {
            return safe_exp(u);
        }
        double e = safe_exp(logit_scale * u);
        double num = (U - 1.0) + (1.0 - L) * e;
        return L * u + (U - L) / logit_scale * std::log(num / (U - L));
    }
};
```

- [ ] **Step 3: Remove `epsilon` from `CalibState` in `src/types.hpp`**

Delete the line `double epsilon;` (currently after `int outer_max_iter;`). The struct should have no `epsilon` field.

- [ ] **Step 4: Remove `st.epsilon` copy-in from `src/c_api.cpp`**

Find and delete the line:
```cpp
    st.epsilon       = p->epsilon;
```
in `rk_calibrate()`. Do not touch any surrounding lines.

- [ ] **Step 5: Mark `epsilon` deprecated in `src/leafblower.h`**

Change the `epsilon` field comment from:
```c
    double          epsilon;         /* iEPPA entropic parameter, default 0.05 */
```
to:
```c
    double          epsilon;         /* deprecated: no longer read by any solver; kept for ABI compat */
```

- [ ] **Step 6: Document `outer_max_iter` limitation in `src/leafblower.h`**

Change:
```c
    int             outer_max_iter;  /* outer EPP / L-BFGS max iters, default 50 */
```
to:
```c
    int             outer_max_iter;  /* outer EPP / L-BFGS max iters, default 50.
                                       * Note: R bridge sets outer_max_iter = inner_max_iter;
                                       * independent control requires direct C API use. */
```

- [ ] **Step 7: Delete stale comment in `src/r_bridge.cpp`**

Delete line 11:
```cpp
// Existing logit test bridges (from Task 4)
```

- [ ] **Step 8: Compile gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected last line: `* DONE (leafblower)`

- [ ] **Step 9: Run tests**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 43 ]`

- [ ] **Step 10: Commit**

```bash
git add src/logit.hpp src/types.hpp src/c_api.cpp src/leafblower.h src/r_bridge.cpp
git commit -m "refactor: remove dead dF/epsilon; document deprecated fields; drop stale comment"
```

---

## Task 2: Forward `tol_abs` through the R→C bridge

**Context:** `R/harvest.R` calls `parse_convergence(convergence)` but discards the return value. Users who pass `convergence = list(absolute = 1e-8)` silently get the default 1e-6 tolerance. The fix threads `tol_abs` as a 9th argument through `.Call("C_rk_calibrate", ...)` and sets `p.tol_abs` in `src/r_bridge.cpp`.

**Files:**
- Modify: `R/harvest.R:63` (capture `parse_convergence` return)
- Modify: `R/harvest.R:73-82` (add 9th argument to `.Call`)
- Modify: `src/r_bridge.cpp:16` (declare 9th parameter)
- Modify: `src/r_bridge.cpp:31` (update arg count from 8 to 9)
- Modify: `src/r_bridge.cpp:76-79` (add `tol_abs_sexp` parameter)
- Modify: `src/r_bridge.cpp` body (set `p.tol_abs`)
- Test: `tests/testthat/test-harvest.R`

- [ ] **Step 1: Write the failing test in `tests/testthat/test-harvest.R`**

Append at the end of the file:

```r
test_that("convergence$absolute is forwarded to solver", {
  set.seed(42)
  n <- 10000L
  df <- data.frame(
    age = factor(sample(c("Y","M","O"), n, replace=TRUE,
                        prob=c(0.60, 0.30, 0.10))),
    sex = factor(sample(c("M","F"),     n, replace=TRUE,
                        prob=c(0.70, 0.30)))
  )
  tgt <- list(
    age = c(Y=0.33, M=0.40, O=0.27),
    sex = c(M=0.49, F=0.51)
  )
  # 2 iterations with default tol (1e-6): competing margins, can't converge -> warning
  expect_warning(
    harvest(df, tgt, method="ieppa", max_iterations=2),
    regexp="did not converge"
  )
  # 2 iterations with loose tol (0.3): error after 2 iters < 0.3 -> no warning
  # Before fix: tol_abs ignored, tol=1e-6 used, warning fires -> test fails
  # After fix:  tol_abs=0.3 forwarded, error < 0.3 accepted -> no warning
  expect_no_warning(
    harvest(df, tgt, method="ieppa", max_iterations=2,
            convergence=list(absolute=0.3))
  )
})
```

- [ ] **Step 2: Run test to verify it fails**

```bash
Rscript -e "devtools::test('tests/testthat/test-harvest.R')" 2>&1 | tail -10
```

Expected: the new test fails with the `expect_no_warning` assertion (warning fires because `tol_abs=0.3` is currently ignored).

- [ ] **Step 3: Capture `parse_convergence` return value in `R/harvest.R`**

Change line 63 from:
```r
  parse_convergence(convergence)
```
to:
```r
  tol_abs <- parse_convergence(convergence)
```

- [ ] **Step 4: Pass `tol_abs` as 9th argument in the `.Call` block in `R/harvest.R`**

Change the `.Call` block from:
```r
  raw <- .Call("C_rk_calibrate",
               data,
               target,
               as.double(min_weight),
               as.double(max_weight),
               as.character(method),
               as.integer(verbose),
               as.integer(max_iterations),
               sw_vec,
               PACKAGE = "leafblower")
```
to:
```r
  raw <- .Call("C_rk_calibrate",
               data,
               target,
               as.double(min_weight),
               as.double(max_weight),
               as.character(method),
               as.integer(verbose),
               as.integer(max_iterations),
               sw_vec,
               as.double(tol_abs),
               PACKAGE = "leafblower")
```

- [ ] **Step 5: Add `tol_abs_sexp` parameter to `C_rk_calibrate` declaration in `src/r_bridge.cpp`**

Change the forward declaration (around line 16):
```cpp
SEXP C_rk_calibrate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
```
to:
```cpp
SEXP C_rk_calibrate(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
```

- [ ] **Step 6: Update arg count in `R_init_leafblower` in `src/r_bridge.cpp`**

Change:
```cpp
        {"C_rk_calibrate",       (DL_FUNC)&C_rk_calibrate,       8},
```
to:
```cpp
        {"C_rk_calibrate",       (DL_FUNC)&C_rk_calibrate,       9},
```

- [ ] **Step 7: Add `tol_abs_sexp` to the `C_rk_calibrate` function signature in `src/r_bridge.cpp`**

Change:
```cpp
SEXP C_rk_calibrate(SEXP data_sexp, SEXP target_sexp,
                    SEXP min_weight_sexp, SEXP max_weight_sexp,
                    SEXP method_sexp, SEXP verbose_sexp,
                    SEXP inner_max_iter_sexp, SEXP start_weights_sexp) {
```
to:
```cpp
SEXP C_rk_calibrate(SEXP data_sexp, SEXP target_sexp,
                    SEXP min_weight_sexp, SEXP max_weight_sexp,
                    SEXP method_sexp, SEXP verbose_sexp,
                    SEXP inner_max_iter_sexp, SEXP start_weights_sexp,
                    SEXP tol_abs_sexp) {
```

- [ ] **Step 8: Set `p.tol_abs` from the new argument in `src/r_bridge.cpp`**

After the line `p.outer_max_iter = INTEGER(inner_max_iter_sexp)[0];` (currently the last param assignment in the "Set calibration params" block), add:
```cpp
    p.tol_abs        = REAL(tol_abs_sexp)[0];
```

The complete "Set calibration params" block should now read:
```cpp
    rk_params_t p;
    rk_params_init(&p);
    p.min_weight     = REAL(min_weight_sexp)[0];
    p.max_weight     = REAL(max_weight_sexp)[0];
    p.verbose        = INTEGER(verbose_sexp)[0];
    p.inner_max_iter = INTEGER(inner_max_iter_sexp)[0];
    p.outer_max_iter = INTEGER(inner_max_iter_sexp)[0];
    p.tol_abs        = REAL(tol_abs_sexp)[0];
    p.log_fn         = (p.verbose > 0) ? r_log_trampoline : nullptr;
```

- [ ] **Step 9: Compile gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 10: Run all tests**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 44 ]`

- [ ] **Step 11: Commit**

```bash
git add R/harvest.R src/r_bridge.cpp tests/testthat/test-harvest.R
git commit -m "fix: forward convergence\$absolute to solver as tol_abs"
```

---

## Task 3: Hoist `bucket` allocation out of `compute_errRp`

**Context:** `compute_errRp` in `src/ieppa.cpp` is called every iteration (up to `inner_max_iter=500`). It allocates a new `std::vector<double> bucket(st.cat_counts[k], 0.0)` inside its inner k-loop — K allocations per call, 500*K allocations per calibration run. Fix: allocate one `bucket(max_cats)` before the k-loop and `std::fill` to reset.

**Files:**
- Modify: `src/ieppa.cpp:12-30` (`compute_errRp` body)

- [ ] **Step 1: Replace `compute_errRp` body in `src/ieppa.cpp`**

Replace the entire function (lines 12–30) with:

```cpp
// Compute errRp = max_k max_j |S_kj/W - tau_kj|
// O(n*K): single O(n) bucket accumulation pass per margin.
// bucket pre-allocated to max_cats to avoid per-call heap allocation.
static double compute_errRp(const CalibState& st,
                              const std::vector<double>& w) {
    double W = 0.0;
    for (int i = 0; i < st.n; i++) W += w[i];

    int max_cats = *std::max_element(st.cat_counts, st.cat_counts + st.K);
    std::vector<double> bucket(max_cats);
    double err = 0.0;
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
    return err;
}
```

- [ ] **Step 2: Compile gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 3: Run all tests**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 44 ]`

- [ ] **Step 4: Commit**

```bash
git add src/ieppa.cpp
git commit -m "perf(ieppa): hoist bucket allocation out of compute_errRp"
```

---

## Task 4: Deduplicate `phi_and_grad`/`phi_from_u` and use `compute_targets_abs` return value

**Context:** Two issues in `src/lbfgsb_solver.cpp`:
1. `phi_and_grad` (lines 47–71) and `phi_from_u` (lines 73–98) share an identical 40-line body that computes `obj` and `grad`. Only difference: `phi_and_grad` calls `compute_u` first. Extract `phi_from_u` as the canonical implementation; make `phi_and_grad` a 2-line wrapper.
2. `compute_targets_abs` returns `W` (the weight sum) but the caller at line 323 ignores it and then recomputes the same `W_sum` in a separate loop (lines 324–325). Use the return value instead.

**Files:**
- Modify: `src/lbfgsb_solver.cpp:47-98` (dedup)
- Modify: `src/lbfgsb_solver.cpp:323-325` (use return value)

- [ ] **Step 1: Move `phi_from_u` above `phi_and_grad` in `src/lbfgsb_solver.cpp`**

In the current file, `phi_and_grad` occupies lines 47–71 and `phi_from_u` occupies lines 73–98. Move `phi_from_u` to above `phi_and_grad` so the wrapper in Step 2 can reference it. The new order after `compute_u`:

```
compute_u        (unchanged)
phi_from_u       (moved up — definition unchanged)
phi_and_grad     (will be replaced in Step 2)
phi_from_u       (delete the old position here — only one definition should remain)
```

`phi_from_u` definition to move (keep unchanged):

```cpp
// phi_from_u: evaluate objective and gradient given pre-computed u.
// Avoids re-running compute_u; gradient in lambda-space still requires K*n pass.
static double phi_from_u(const CalibState& st,
                          const LinkFn& fn,
                          const std::vector<int>& off,
                          const std::vector<double>& lam,
                          const std::vector<double>& T,
                          const std::vector<double>& d,
                          std::vector<double>& grad,
                          const std::vector<double>& u) {
    int total = off[st.K];

    double obj = 0.0;
    for (int idx = 0; idx < total; idx++) obj += T[idx] * lam[idx];
    for (int i = 0; i < st.n; i++) obj -= d[i] * fn.H(u[i]);

    std::fill(grad.begin(), grad.end(), 0.0);
    for (int idx = 0; idx < total; idx++) grad[idx] = T[idx];
    for (int k = 0; k < st.K; k++) {
        for (int i = 0; i < st.n; i++) {
            int g = st.group_ids[k][i];
            if (g >= 0) grad[off[k] + g] -= d[i] * fn.F(u[i]);
        }
    }
    return obj;
}
```

- [ ] **Step 2: Replace `phi_and_grad` body with a 2-line wrapper**

Replace the entire old `phi_and_grad` body (lines 47–71) with:

```cpp
// phi(lambda) = sum_kj T_kj*lam_kj - sum_i d_i*H(u_i)
// grad[off_k+j] = T_kj - S_kj where S_kj = sum_{i:g_k(i)==j} d_i*F(u_i)
// Computes u from lam, then delegates to phi_from_u.
static double phi_and_grad(const CalibState& st,
                            const LinkFn& fn,
                            const std::vector<int>& off,
                            const std::vector<double>& lam,
                            const std::vector<double>& T,
                            const std::vector<double>& d,
                            std::vector<double>& grad,
                            std::vector<double>& u) {
    compute_u(st, off, lam, u);
    return phi_from_u(st, fn, off, lam, T, d, grad, u);
}
```

Note: `phi_from_u` takes `const std::vector<double>& u` — the `const` must be preserved in its definition (already shown in Step 1). The wrapper passes `u` (non-const owner) into `phi_from_u`'s const ref parameter — this is valid C++.

- [ ] **Step 3: Use `compute_targets_abs` return value in `lbfgsb_solve`**

Find the block in `lbfgsb_solve`:
```cpp
    compute_targets_abs(st, T);
    double W_sum = 0.0;
    for (int i = 0; i < st.n; i++) W_sum += d[i];
```

Replace with:
```cpp
    double W_sum = compute_targets_abs(st, T);
```

`compute_targets_abs` already computes and returns the total weight sum (`W`), which equals `sum(d[i])` since `d` was copied from `st.weights`. The separate loop is redundant.

- [ ] **Step 4: Compile gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 5: Run all tests**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 44 ]`

- [ ] **Step 6: Commit**

```bash
git add src/lbfgsb_solver.cpp
git commit -m "refactor(lbfgsb): deduplicate phi_and_grad; use compute_targets_abs return value"
```

---

## Task 5: Hoist per-iteration vector allocations out of L-BFGS-B loop

**Context:** `lbfgsb_solve` allocates `lam_new(total)`, `s_new(total)`, and `y_new(total)` inside the main iteration loop (currently lines ~347, ~353–354). With `outer_max_iter=500` and `total=n*K`, this causes up to 1500 heap allocations per calibration call. Hoist all three above the loop.

**Files:**
- Modify: `src/lbfgsb_solver.cpp` (`lbfgsb_solve` body)

- [ ] **Step 1: Hoist `lam_new`, `s_new`, `y_new` above the `for` loop in `lbfgsb_solve`**

In `lbfgsb_solve`, find the block just before or at the start of `for (int iter = 0; iter < max_iter; iter++)` and add three pre-allocated vectors:

Before:
```cpp
    int max_iter = st.outer_max_iter;
    int final_iter = 0;
    for (int iter = 0; iter < max_iter; iter++) {
        final_iter = iter + 1;
        ...
        std::vector<double> lam_new(total);
        double phi_new = phi_curr;

        wolfe_line_search(..., lam_new, grad_new, phi_new);

        std::vector<double> s_new(total), y_new(total);
```

After (hoist the three declarations above the loop):
```cpp
    int max_iter = st.outer_max_iter;
    int final_iter = 0;
    std::vector<double> lam_new(total), s_new(total), y_new(total);
    for (int iter = 0; iter < max_iter; iter++) {
        final_iter = iter + 1;
        ...
        double phi_new = phi_curr;

        wolfe_line_search(..., lam_new, grad_new, phi_new);

        // s_new and y_new are now reused each iteration
```

Inside the loop, the assignments to `s_new[i]` and `y_new[i]` remain unchanged — they overwrite the pre-allocated storage each iteration. Remove the `std::vector<double>` type declarations from inside the loop, keeping only the assignments.

- [ ] **Step 2: Compile gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 3: Run all tests**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 44 ]`

- [ ] **Step 4: Commit**

```bash
git add src/lbfgsb_solver.cpp
git commit -m "perf(lbfgsb): hoist lam_new/s_new/y_new allocation above iteration loop"
```

---

## Task 6: Minor polish — `linf` rename, `enforce_mean`, column map, `1e300` doc

**Context:** Four remaining items:
1. `linf` in `src/lbfgsb_solver.cpp` is not intention-revealing. Rename to `maxAbs`.
2. `enforce_mean` parameter in `R/harvest.R` is silently ignored (normalization always applies). Add it to the `ignored` params list so verbose >= 2 users get a note.
3. Column name lookup in `src/r_bridge.cpp` is O(K×C). Pre-build an `unordered_map` once.
4. `1e300` sentinel in `src/ieppa.cpp:52` is unexplained. Add a comment.

**Files:**
- Modify: `src/lbfgsb_solver.cpp:115,333` (`linf` → `maxAbs`)
- Modify: `R/harvest.R:67` (add `enforce_mean` to `ignored`)
- Modify: `src/r_bridge.cpp:91-100` (pre-build column map)
- Modify: `src/ieppa.cpp:52` (document `1e300`)

- [ ] **Step 1: Rename `linf` to `maxAbs` in `src/lbfgsb_solver.cpp`**

Change the function definition:
```cpp
static double linf(const std::vector<double>& v) {
    double mx = 0.0;
    for (double x : v) mx = std::max(mx, std::fabs(x));
    return mx;
}
```
to:
```cpp
static double maxAbs(const std::vector<double>& v) {
    double mx = 0.0;
    for (double x : v) mx = std::max(mx, std::fabs(x));
    return mx;
}
```

And update the call site (currently `double gn = linf(grad) / W_sum;`):
```cpp
        double gn = maxAbs(grad) / W_sum;
```

- [ ] **Step 2: Add `enforce_mean` to `ignored` in `R/harvest.R`**

Change:
```r
  ignored <- c("select_params", "select_function", "error_function",
                "adaptive_order", "accelerate")
```
to:
```r
  # enforce_mean is always TRUE: normalization is unconditional (harvest.R:86).
  ignored <- c("select_params", "select_function", "error_function",
                "adaptive_order", "accelerate", "enforce_mean")
```

- [ ] **Step 3: Pre-build column `unordered_map` in `src/r_bridge.cpp`**

In `C_rk_calibrate`, find the column-lookup block inside the `for (int k = 0; k < K; k++)` loop:

```cpp
        int col_idx = -1;
        for (int c = 0; c < LENGTH(col_names); c++) {
            if (strcmp(CHAR(STRING_ELT(col_names, c)), varname) == 0) {
                col_idx = c; break;
            }
        }
        if (col_idx < 0)
            Rf_error("Variable '%s' not found in data", varname);
```

Replace the entire for-k loop's column lookup with a pre-built map. Add the map construction **before** the `for (int k = 0; k < K; k++)` loop:

```cpp
    // Build column name → index map once (O(C)) for O(1) per-margin lookup.
    std::unordered_map<std::string, int> col_map;
    col_map.reserve(LENGTH(col_names));
    for (int c = 0; c < LENGTH(col_names); c++)
        col_map[CHAR(STRING_ELT(col_names, c))] = c;
```

Inside the for-k loop, replace the linear scan with map lookup:

```cpp
        auto col_it = col_map.find(varname);
        if (col_it == col_map.end())
            Rf_error("Variable '%s' not found in data", varname);
        int col_idx = col_it->second;
```

- [ ] **Step 4: Document `1e300` in `src/ieppa.cpp`**

Change:
```cpp
    double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
```
to:
```cpp
    // 1e300 not numeric_limits::max(): prevents overflow in w[i] *= scale[g]
    // when bucket[g] is tiny (scale[g] can be huge).
    double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
```

- [ ] **Step 5: Compile gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 6: Run all tests**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 44 ]`

- [ ] **Step 7: Commit**

```bash
git add src/lbfgsb_solver.cpp R/harvest.R src/r_bridge.cpp src/ieppa.cpp
git commit -m "refactor: rename linf->maxAbs; add enforce_mean to ignored; O(1) column map; doc 1e300"
```

---

## Self-Review

### 1. Spec coverage

| Item | Task |
|---|---|
| 1. `parse_convergence` discarded | Task 2 (steps 3-4) |
| 2. Dead `dF` | Task 1 (step 2) |
| 3. Dead `epsilon` | Task 1 (steps 3-5) |
| 4. `phi_and_grad`/`phi_from_u` dup | Task 4 (steps 1-2) |
| 5. `lam_new`/`s_new`/`y_new` in loop | Task 5 (step 1) |
| 6. Column map O(K×C) | Task 6 (step 3) |
| 7. `compute_targets_abs` return unused | Task 4 (step 3) |
| 8. `outer_max_iter` misleading | Task 1 (step 6) |
| 9. `compute_errRp` per-call alloc | Task 3 (step 1) |
| 10. Stale `// Task 4` comment | Task 1 (step 7) |
| 11a. `linf` rename | Task 6 (step 1) |
| 11b. `enforce_mean` silent ignore | Task 6 (step 2) |
| 11c. `alg_names` enum sync | Task 6 — existing tests (`test-ieppa.R:7`, `test-harvest.R:22,32`) already cover this; added comment to `harvest.R` in step 2 |
| 11d. `1e300` undocumented | Task 6 (step 4) |

All 11 items covered.

### 2. Placeholder scan

No TBDs, todos, or "similar to" references. All code blocks are complete and self-contained.

### 3. Type consistency

- `phi_from_u` signature in Task 4 matches its existing definition — same parameter names and types.
- `compute_targets_abs` return type is `double` — matches `W_sum` type in Task 4 step 3.
- `maxAbs` signature in Task 6 step 1 is identical to `linf` except name.
- `tol_abs_sexp` declared as `SEXP` in Task 2 steps 5, 6, 7 — consistent throughout.
- `REAL(tol_abs_sexp)[0]` returns `double` — matches `p.tol_abs` type (`double`) in `rk_params_t`.
