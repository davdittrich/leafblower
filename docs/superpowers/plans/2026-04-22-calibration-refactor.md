# Calibration Refactor Implementation Plan

> **For agentic workers:** Use `superpowers:subagent-driven-development` or `superpowers:executing-plans` to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Remove `method="auto"` from the R/Python API, make iEPPA the unconditional default, make L-BFGS-B a first-class explicit option, and make both solvers enforce `sum(w) = n` (equivalently `mean(w) = 1`) as a hard constraint so bounds are honored without post-hoc fixup loops.

**Architecture:**
- **Phase 1 (P1)**: API surface only. Remove `method="auto"` from `R/harvest.R`, `python/leafblower/_harvest.py`, and `src/r_bridge.cpp`. Remove `select_algorithm()` and `kComplexityThreshold` from `src/c_api.cpp`. Delete obsolete tests. Ship user-visible API change.
- **Phase 2 (P2)**: Remove in-loop normalization (`ieppa.cpp:116-133`) from iEPPA; replace with Dykstra hyperplane projection onto `{w : sum(w) = n}` with correction vector `q_hyp[]`. Remove post-solve fixup loop (`ieppa.cpp:173-190`).
- **Phase 3 (P3)**: Refactor `lbfgsb_solve` to extract `lbfgsb_solve_inner`; wrap in ALM outer loop enforcing `sum(w) = n` with correct gradient in λ_kj dual space. Remove the unconditional post-normalization in `R/harvest.R:97` once both solvers guarantee mean=1.

**Tech stack:** C++17, Rcpp (r_bridge), R 4.x, testthat 3.x, Python 3.x, `R CMD INSTALL --preclean .` as the per-file compile gate.

---

## Phase 1 — Remove AUTO routing (user-visible API change)

### Task 1: Write failing test — default method routes to iEPPA

**Files:**
- Modify: `tests/testthat/test-algo-selection.R` (append)

- [ ] **Step 1: Add test**

```r
test_that("default method (no method arg) routes to ieppa", {
  set.seed(7L)
  n  <- 30000L
  df <- data.frame(
    m1 = factor(sample(paste0("c", 1:4), n, replace = TRUE)),
    m2 = factor(sample(paste0("c", 1:4), n, replace = TRUE))
  )
  tgt <- list(
    m1 = c(c1 = 0.25, c2 = 0.25, c3 = 0.25, c4 = 0.25),
    m2 = c(c1 = 0.25, c2 = 0.25, c3 = 0.25, c4 = 0.25)
  )
  res <- leafblower::harvest(df, tgt,
                              max_weight = Inf, min_weight = 0,
                              convergence = list(absolute = 1e-3))
  expect_equal(attr(res, "algorithm"), "ieppa")
})
```

- [ ] **Step 2: Run to verify FAIL (RED)**

`Rscript -e "devtools::test(filter='algo-selection')"`

Expected: FAIL with `"lbfgsb" != "ieppa"`. Current `harvest.R:32` default is `method="auto"`; `c_api.cpp:141-144` routes unconstrained complexity=240K < 500K threshold to LBFGSB.

---

### Task 2: Write + run tests — explicit lbfgsb routing and caps

**Files:**
- Modify: `tests/testthat/test-algo-selection.R` (append)

- [ ] **Step 1: Add routing test**

```r
test_that("method='lbfgsb' with max_weight routes to lbfgsb", {
  set.seed(1L)
  n  <- 500L
  df <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  res <- leafblower::harvest(df, tgt, method = "lbfgsb", max_weight = 5)
  expect_equal(attr(res, "algorithm"), "lbfgsb")
})
```

- [ ] **Step 2: Add caps-honored test (will fail until P3; document skip)**

```r
test_that("method='lbfgsb' output weights satisfy max_weight/min_weight", {
  skip("Enforced only after Phase 3 ALM lands (beads: P3)")
  set.seed(2L)
  n   <- 500L
  df  <- data.frame(x = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.3, c = 0.2))
  res <- leafblower::harvest(df, tgt, method = "lbfgsb",
                              max_weight = 1.5, min_weight = 0.2)
  expect_true(max(res) <= 1.5 + 1e-6)
  expect_true(min(res) >= 0.2 - 1e-6)
})
```

- [ ] **Step 3: Run tests — routing must PASS, caps test SKIPPED**

`Rscript -e "devtools::test(filter='algo-selection')"`

---

### Task 3: Remove select_algorithm() from c_api.cpp

**Files:**
- Modify: `src/c_api.cpp:128-145` (delete block), `src/c_api.cpp:162` (fix init), `src/c_api.cpp:169-172` (replace routing), `src/c_api.cpp:193-203` (delete verbose block)

- [ ] **Step 0: Verify RK_ALG_AUTO enum ABI**

Grep for `RK_ALG_AUTO` in all C headers to confirm the enum value stays defined:

```bash
grep -rn 'RK_ALG_AUTO' src/ inst/
```

`RK_ALG_AUTO=0` must remain in the enum definition (for ABI stability — direct C callers may pass it). The Task 3 Step 2 ternary `(p->algorithm == RK_ALG_LBFGSB) ? RK_ALG_LBFGSB : RK_ALG_IEPPA` routes `RK_ALG_AUTO=0` silently to iEPPA — this is the intended behavior. Confirm the enum definition line is NOT deleted by this task (only the `select_algorithm()` function and `kComplexityThreshold` constant are deleted).

- [ ] **Step 1: Delete kComplexityThreshold (line 133) and select_algorithm() (lines 128-145)**

Remove the entire block:
```cpp
// Routing note: std::isfinite(max_weight) fires for any finite upper bound ...
// kComplexityThreshold only differentiates ...
static constexpr int64_t kComplexityThreshold = 500000L;

static rk_algorithm_t select_algorithm(int n, int K,
                                        const int* cat_counts,
                                        const rk_params_t* p,
                                        int64_t& complexity_out) {
    complexity_out = INT64_C(0);
    for (int k = 0; k < K; k++) complexity_out += (int64_t)n * cat_counts[k];
    if (p->algorithm != RK_ALG_AUTO) return p->algorithm;
    if (complexity_out > kComplexityThreshold || std::isfinite(p->max_weight) || p->min_weight > 0.0)
        return RK_ALG_IEPPA;
    return RK_ALG_LBFGSB;
}
```

- [ ] **Step 2: Replace routing block at lines 169-172**

The null-guard `(cat_counts && K > 0 && n > 0)` must be preserved — it protects against degenerate inputs before accessing cat_counts[]. Replace:

Before:
```cpp
int64_t complexity = INT64_C(0);
rk_algorithm_t alg = (cat_counts && K > 0 && n > 0)
    ? select_algorithm(n, K, cat_counts, p, complexity)
    : p->algorithm;
```

After:
```cpp
rk_algorithm_t alg = (cat_counts && K > 0 && n > 0)
    ? ((p->algorithm == RK_ALG_LBFGSB) ? RK_ALG_LBFGSB : RK_ALG_IEPPA)
    : p->algorithm;
```

- [ ] **Step 3: Fix result init at line 162**

Before:
```cpp
result->algorithm_used = RK_ALG_AUTO;
```

After:
```cpp
result->algorithm_used = RK_ALG_IEPPA;
```

- [ ] **Step 4: Delete verbose AUTO routing block (lines 193-203)**

Remove the `if (p->verbose >= 1 && p->algorithm == RK_ALG_AUTO) { ... }` block entirely.

- [ ] **Step 5: Compile gate**

`R CMD INSTALL --preclean . 2>&1 | tail -10`

Expected: `* DONE (leafblower)` and 0 errors.

---

### Task 4: Remove "auto" from R API (harvest.R and Python bindings)

**Files:**
- Modify: `R/harvest.R:7` (roxygen), `R/harvest.R:32` (default), `R/harvest.R:104` (comment), `R/harvest.R:150` (match.arg)
- Modify: `python/leafblower/_harvest.py:23` (default), `python/leafblower/_harvest.py:41` (docstring), `python/leafblower/_harvest.py:87` (alg_map)
- Regenerate: `man/harvest.Rd`

#### harvest.R

- [ ] **Step 1: Change default at line 32**

Before:
```r
method = "auto",
```

After:
```r
method = "ieppa",
```

- [ ] **Step 2: Change match.arg at line 150**

Before:
```r
match.arg(method, c("auto", "ieppa", "lbfgsb"))
```

After:
```r
match.arg(method, c("ieppa", "lbfgsb"))
```

- [ ] **Step 3: Update roxygen at line 7**

Before:
```r
#' @param method One of "auto", "ieppa", "lbfgsb", "rake", "nr". Default "auto".
```

After (verify `rake`/`nr` are still supported by checking match.arg before editing):
```r
#' @param method One of "ieppa", "lbfgsb". Default "ieppa".
```

- [ ] **Step 4: Add defensive comment to alg_names at line 104**

`alg_names` stays as-is — the enum index 0 (auto) is dead after the match.arg change above, but the vector is kept for 1-indexed enum alignment. Add a one-line comment:

Before:
```r
alg_names <- c("auto", "ieppa", "lbfgsb")
```

After:
```r
alg_names <- c("auto", "ieppa", "lbfgsb")  # index 0 (auto) unreachable after match.arg; kept for 1-indexed enum alignment
```

#### Python bindings

- [ ] **Step 5: Update python/leafblower/_harvest.py**

Read the file first to get exact line content, then:

At line 23, change default:
```python
# Before:
method: str = "auto"
# After:
method: str = "ieppa"
```

At line 41 (docstring), update the method description to remove "auto" from the list of valid choices and update the default description.

At line 87, update alg_map:
```python
# Before:
alg_map = {"auto": 0, "ieppa": 1, "lbfgsb": 2}
# After:
alg_map = {"ieppa": 1, "lbfgsb": 2}  # "auto" (0) removed from user API
```

#### Finalize

- [ ] **Step 6: Verify Python binding end-to-end**

```bash
python -c "
import pandas as pd, numpy as np
try:
    import leafblower
    df = pd.DataFrame({'x': np.random.choice(['a','b'], 100)})
    res = leafblower.harvest(df, {'x': {'a': 0.5, 'b': 0.5}})
    print('default method: PASS, algo =', res.attrs.get('algorithm'))
    res2 = leafblower.harvest(df, {'x': {'a': 0.5, 'b': 0.5}}, method='lbfgsb')
    print('lbfgsb method: PASS, algo =', res2.attrs.get('algorithm'))
    try:
        leafblower.harvest(df, {'x': {'a': 0.5, 'b': 0.5}}, method='auto')
        print('auto not rejected: FAIL')
    except Exception as e:
        print('auto rejected: PASS —', str(e)[:60])
except Exception as e:
    print('FAIL:', e)
"
```

Expected: default routes to ieppa, lbfgsb routes to lbfgsb, `method='auto'` raises an error.

- [ ] **Step 7: Regenerate docs**

`Rscript -e "devtools::document()"`

- [ ] **Step 9: Verify match.arg rejects "auto"**

```r
Rscript -e '
library(leafblower)
tryCatch(
  harvest(data.frame(x=factor(c("a","b"))), list(x=c(a=0.5,b=0.5)), method="auto"),
  error = function(e) cat("Error correctly raised:", conditionMessage(e), "\n")
)
'
```

Expected: `Error correctly raised: 'arg' should be one of "ieppa", "lbfgsb"`

---

### Task 5: Update r_bridge.cpp fallback

**Files:**
- Modify: `src/r_bridge.cpp:175`

- [ ] **Step 1: Change fallback**

Before:
```cpp
else                                          p.algorithm = RK_ALG_AUTO;
```

After:
```cpp
else                                          p.algorithm = RK_ALG_IEPPA;
```

- [ ] **Step 2: Compile gate**

`R CMD INSTALL --preclean . 2>&1 | tail -10`

Expected: 0 errors. After Task 4 match.arg filters "auto" out at R level, this branch is unreachable from R but kept as a safe default for direct C callers.

---

### Task 6: Fix existing tests using method="auto"

**Files:**
- Modify: `tests/testthat/test-harvest.R` (DELETE lines 25-33; fix any remaining method="auto")
- Modify: `tests/testthat/test-lbfgsb.R` (DELETE sub-block at lines 47-50 only)
- Modify: `tests/testthat/test-algo-selection.R` (remove method="auto"; delete obsolete tests)

- [ ] **Step 1: Audit call sites**

```bash
grep -rn 'method[[:space:]]*=[[:space:]]*"auto"' tests/
```

Record all line numbers. Each must be addressed.

- [ ] **Step 2: Handle method="auto" blocks in test-harvest.R**

There are two `method="auto"` call sites in `test-harvest.R`:

**Lines 15-23** (`test_that` asserting `algorithm == "ieppa"` for iEPPA routing): Remove the `method = "auto"` argument from the `harvest()` call only. The assertion `algorithm == "ieppa"` will still pass after the default changes to `"ieppa"`. Do NOT delete this test.

**Lines 25-34** (`test_that("auto-routing selects lbfgsb for small unconstrained problems", ...)` asserting `algorithm == "lbfgsb"`): This block must be **deleted entirely**. Removing the `method` arg would make the assertion wrong (`"ieppa"` not `"lbfgsb"`), so patching is not an option. Delete from the opening `test_that(` to the closing `})` (inclusive of line 34).

- [ ] **Step 3: DELETE sub-block at test-lbfgsb.R lines 47-50 ONLY**

The outer `test_that` at lines 29-56 ("near-one max_weight") must **survive intact**. Only lines 47-50 (the sub-block with `# auto:` comment that tests `method="auto"` fallback to iEPPA at the logit singularity boundary) must be deleted.

Read lines 29-56 before editing to confirm the exact sub-block boundaries. Delete only the `# auto:` comment line and the associated `expect_no_error(...)` call. Compile + run the enclosing test to confirm the outer test still passes.

- [ ] **Step 4: Fix test-algo-selection.R**

- Line 178 ("constrained max_weight=5 → ieppa"): change `method = "auto"` → remove the `method` arg.
- Delete lines 182-200 ("unconstrained large complexity routes to L-BFGS-B"): this test was the regression guard for the old AUTO→LBFGSB path which no longer exists.
- Delete lines 202-214 (Case B skip block): obsolete.

- [ ] **Step 5: Fix any remaining method="auto" call sites found in Step 1**

Remove the `method` argument from each remaining call (uses new default `"ieppa"`). Keep all other arguments intact.

- [ ] **Step 6: Run full test suite**

`Rscript -e "devtools::test()"`

Expected: 0 FAIL, all Task 1/2 tests pass, Task 2 caps-test SKIPPED with phase-3 reason.

---

### Task 7a: Commit Phase 1 code changes

**Files:** all source + test files modified above.

- [ ] **Step 1: Stage code files only**

```bash
git add src/c_api.cpp src/r_bridge.cpp R/harvest.R man/harvest.Rd \
        python/leafblower/_harvest.py \
        tests/testthat/test-algo-selection.R \
        tests/testthat/test-harvest.R \
        tests/testthat/test-lbfgsb.R
```

- [ ] **Step 2: Commit**

```bash
git commit -m "$(cat <<'EOF'
refactor(api): remove method="auto", make iEPPA the default

Benchmark (2026-04-22) showed iEPPA dominates L-BFGS-B across the full
tested (log_complexity, log_tol) space with no 1.2x contour found. The
AUTO routing logic and kComplexityThreshold constant are removed:

- R: method="auto" rejected at match.arg; default is "ieppa".
- Python: method="auto" default changed to "ieppa"; "auto" removed from alg_map.
- C: select_algorithm() and kComplexityThreshold removed; routing is
  method=="lbfgsb" ? LBFGSB : IEPPA with the null-guard preserved.
- L-BFGS-B caps honoring tracked as Phase 3 (ALM outer loop).
EOF
)"
```

---

### Task 7b: Commit benchmark artifacts separately

**Files:** benchmark outputs from the 2026-04-22 LSE benchmark run.

- [ ] **Step 1: Stage benchmark artifacts**

```bash
git add benchmarks/algo_selection_results.rds \
        benchmarks/algo_selection_contour.pdf \
        benchmarks/algo_selection_uncertainty.pdf \
        benchmarks/algo_selection_k_stability.pdf
```

- [ ] **Step 2: Commit**

```bash
git commit -m "$(cat <<'EOF'
chore(bench): add 2026-04-22 Bayesian LSE benchmark outputs

8 LHC evals, classified=100%, 0 adaptive evals needed. All y << log(1.2).
No 1.2x contour found: iEPPA dominates across the full tested space.
These artifacts justify the AUTO routing removal in the preceding commit.
EOF
)"
```

---

## Phase 2 — iEPPA sum(w)=n Dykstra projection

### Task 8: Write regression guard — iEPPA output weights have mean=1

**Important note on RED/GREEN state:** The existing post-solve fixup loop at `ieppa.cpp:173-190` already produces `mean(w) = 1` to machine epsilon at termination. When the fixup loop exits with `!changed`, the final normalize step sets `w[i] /= wm` where `wm = Wsum/n`, so `Σ(w[i]/wm) = n` exactly (by definition of wm). A test on `mean(res) == 1.0` with `tolerance = 1e-10` will therefore PASS on current code and cannot serve as a RED test.

This test is therefore written as a **regression guard** (GREEN from the start). Its role is to catch if the P2 refactor accidentally breaks the mean=1 guarantee.

**Files:**
- Modify: `tests/testthat/test-ieppa.R` (or create if absent)

- [ ] **Step 1: Add regression guard**

```r
test_that("iEPPA output weights have mean=1 and respect bounds", {
  set.seed(5L)
  n   <- 1000L
  df  <- data.frame(x = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.3, c = 0.2))
  res <- leafblower::harvest(df, tgt, method = "ieppa",
                              max_weight = 2.0, min_weight = 0.2)
  # mean=1 is guaranteed by both the old fixup loop and the new Dykstra projection;
  # this test guards against regressions in the P2 refactor.
  expect_equal(mean(res), 1.0, tolerance = 1e-10)
  expect_true(max(res) <= 2.0 + 1e-10)
  expect_true(min(res) >= 0.2 - 1e-10)
})
```

- [ ] **Step 2: Run — must PASS (regression guard, not RED test)**

`Rscript -e "devtools::test(filter='ieppa')"`

---

### Task 9: Replace in-loop normalization with Dykstra hyperplane projection

**Files:**
- Modify: `src/ieppa.cpp`

**What to remove — in-loop normalization block (lines 116-133):**

The current code normalizes `w` and rescales `q[]` inside the main loop:
```cpp
// Normalize to mean=1 so box bounds match the scale R returns.
{
    double Wsum = 0.0;
    for (int i = 0; i < st.n; i++) Wsum += w[i];
    double wm = Wsum / st.n;
    if (wm > kWeightCollapseThreshold) {
#if defined(_OPENMP) || LBW_HAS_OMP_SIMD
#pragma omp simd
#endif
        for (int i = 0; i < st.n; i++) { w[i] /= wm; q[i] /= wm; }
    }
}
```

This block must be **deleted entirely**. The `q[i] /= wm` rescaling violates the Dykstra fixed-point invariant: the correction `q[i]` accumulates "overshoot" from the box projection and must not be rescaled — rescaling shifts the invariant point of the box projection. Removing this block is the core of P2.

**What to add — Dykstra hyperplane projection:**

After the box projection (currently lines 137-145), add a Dykstra projection onto the hyperplane `{w : sum(w) = n}`:

```cpp
// Dykstra hyperplane projection: {w : sum(w) = n}
// q_hyp[i] accumulates overshoot from previous hyperplane projections.
// Dykstra invariant: y_i = w[i] + q_hyp[i] is the pre-projection iterate.
{
    for (int i = 0; i < st.n; i++) w[i] += q_hyp[i];
    double s = 0.0;
    for (int i = 0; i < st.n; i++) s += w[i];
    double shift = (static_cast<double>(st.n) - s) / static_cast<double>(st.n);
    for (int i = 0; i < st.n; i++) {
        double w_proj = w[i] + shift;
        q_hyp[i] = w[i] - w_proj;   // Dykstra correction = pre - post projection
        w[i] = w_proj;
    }
}
```

`q_hyp` must be declared alongside `q`:
```cpp
std::vector<double> q(st.n, 0.0);
std::vector<double> q_hyp(st.n, 0.0);  // ADD: Dykstra correction for hyperplane
```

**Verify Dykstra correction sign consistency** by comparing with the existing box correction in lines 140-145:
```cpp
// Box step (existing):
double yi = w[i] + q[i];
double wc = std::clamp(yi, lo, hi);
q[i] = yi - wc;    // yi = pre-projection; wc = post-projection; q = pre - post
w[i] = wc;
```
The hyperplane correction above uses the same sign convention: `q_hyp[i] = w[i] - w_proj` where `w[i]` is pre-projection and `w_proj` is post-projection.

**What to remove — post-solve fixup loop (lines 173-190):**

```cpp
// Final normalization-and-clamp fixup: ...
bool fixup_converged = false;
for (int fixup = 0; fixup < kMaxFixupIterations; fixup++) {
    ...
}
if (!fixup_converged)
    st.log("iEPPA: fixup loop did not reach fixed point ...");
```

Delete this entire block. The Dykstra hyperplane projection inside the main loop now guarantees `sum(w) = n` at every convergence check; the post-hoc fixup is no longer needed or correct.

Also delete the now-unused constants:
```cpp
static constexpr int    kMaxFixupIterations      = 20;
static constexpr double kWeightCollapseThreshold = 1e-300;
```

**Revised iteration order (after P2):**

```
for each iteration:
  1. IPF margins (multiplicative, K steps)      [unchanged]
  2. Box projection [lo,hi]^n + Dykstra q[]     [unchanged]
  3. Hyperplane projection sum(w)=n + Dykstra q_hyp[]  [NEW, replaces renorm]
  4. Convergence check every kErrCheckInterval   [unchanged]
```

- [ ] **Step 1: Read src/ieppa.cpp in full before editing**

Verify exact line numbers of lines 116-133 (in-loop normalization) and 173-190 (fixup loop) against current file.

- [ ] **Step 2: Delete in-loop normalization block (lines 116-133)**

- [ ] **Step 3: Add q_hyp[] declaration alongside q[] (at approx line 70)**

- [ ] **Step 4: Add Dykstra hyperplane projection block after box projection**

- [ ] **Step 5: Delete post-solve fixup loop (lines 173-190) and unused constants**

- [ ] **Step 6: Compile gate**

`R CMD INSTALL --preclean . 2>&1 | tail -10`

Expected: 0 errors.

- [ ] **Step 7: Run Task 8 regression guard — must still PASS**

`Rscript -e "devtools::test(filter='ieppa')"`

Expected: Task 8 PASS. Full suite 0 FAIL.

- [ ] **Step 8: Run existing calibration benchmark to confirm no regression**

```bash
Rscript benchmarks/stepstone_benchmark.R 2>&1 | tail -20
```

Expected: median harvest time within 10% of pre-change baseline.

---

### Task 10: Commit Phase 2

- [ ] **Step 1: Stage and commit**

```bash
git add src/ieppa.cpp tests/testthat/test-ieppa.R
git commit -m "$(cat <<'EOF'
refactor(ieppa): replace in-loop normalization with Dykstra hyperplane projection

Removes the in-loop w[i]/=wm; q[i]/=wm normalization (ieppa.cpp:116-133)
which violated the Dykstra fixed-point invariant: rescaling q[] shifts the
invariant point of the box projection.

Replaces it with a first-class Dykstra projection onto {w : sum(w) = n}
with correction vector q_hyp[]. This fits cleanly into the alternating
Dykstra cycle (box projection + hyperplane projection) with correct
correction signs.

Also removes the post-solve renormalize-reclamp fixup loop (ieppa.cpp:173-190)
which was a non-standard heuristic compensating for the above bug.
Mean(w) = 1 is now guaranteed by the Dykstra cycle, not post-hoc patching.
EOF
)"
```

---

## Phase 3 — L-BFGS-B sum(w)=n via augmented Lagrangian (ALM)

### Task 11: Un-skip Task 2's caps test; verify it fails (RED for P3)

**Files:**
- Modify: `tests/testthat/test-algo-selection.R` (remove the skip from the caps test)

- [ ] **Step 1: Remove `skip("Enforced only after Phase 3 ALM lands")`**

- [ ] **Step 2: Run — must FAIL**

`Rscript -e "devtools::test(filter='algo-selection')"`

Expected: FAIL on `max(res) <= 1.5 + 1e-6`.

---

### Task 12: Refactor lbfgsb_solve — extract lbfgsb_solve_inner

**Files:**
- Modify: `src/lbfgsb_solver.cpp`

This task is a prerequisite for Task 13. The current `lbfgsb_solve` is monolithic. The ALM outer loop must call an inner solver repeatedly with different `lambda_eq` and `mu` values; this requires the inner L-BFGS-B call to be extractable.

- [ ] **Step 0: Read src/lbfgsb_solver.cpp in full**

Identify the body of `lbfgsb_solve`: the L-BFGS-B call, objective/gradient computation function (e.g. `phi_from_u` or equivalent), and the local variable names for design weights (`d[]`), log-dual iterate (`u[]`), gradient array (`grad_lam[]`), category offsets (`off[]`), and the link function instance (`fn`). These names are needed by Task 13 to write correct ALM gradient code.

- [ ] **Step 1: Grep callers of lbfgsb_solve to verify signature**

```bash
grep -rn 'lbfgsb_solve' src/
```

Confirm the current signature is `LBFGSResult lbfgsb_solve(CalibState& st)` and all call sites pass a single `CalibState&`. The refactor must preserve this exact signature for `lbfgsb_solve` (the ALM wrapper). `lbfgsb_solve_inner` is a new internal function, not exposed.

- [ ] **Step 2: Add ALM fields to CalibState**

`CalibState` is defined in `src/types.hpp:11`. Add:

```cpp
// ALM equality constraint state (used by lbfgsb_solver)
double alm_lambda = 0.0;  // dual variable for sum(w)=n
double alm_mu     = 0.0;  // penalty coefficient; 0.0 = no ALM term
```

When `alm_mu = 0.0` and `alm_lambda = 0.0` (the default), the ALM gradient term `(alm_lambda + alm_mu * residual) * dphi_kj = 0` — no division by zero, augmented objective term is zero. This is correct.

- [ ] **Step 3: Extract lbfgsb_solve_inner**

Create `LBFGSResult lbfgsb_solve_inner(CalibState& st)` containing the existing L-BFGS-B solve body (objective + gradient + L-BFGS-B call). The function reads `st.alm_lambda` and `st.alm_mu` to include the augmented Lagrangian term in the gradient.

- [ ] **Step 4: Compile gate after extraction**

`R CMD INSTALL --preclean . 2>&1 | tail -10`

Expected: 0 errors. Run full test suite: 0 FAIL (pure refactor, no behavior change — alm_lambda=0, alm_mu=0 means augmented term = 0, identical output to pre-refactor).

---

### Task 13: Add F' derivative to LinkFn; add ALM gradient to lbfgsb_solve_inner

**Files:**
- Modify: `src/logit.hpp` (add dF method)
- Modify: `src/lbfgsb_solver.cpp` (add ALM gradient)

**Gradient derivation (dual space):**

The L-BFGS-B optimizer works in dual `λ_kj` space (dimension `Σ cat_counts[k]`). The per-observation weight is `w_i = d_i · F(u_i)` where `u_i = Σ_{k,j: g_k(i)=j} λ_kj` and F is the logit link function. The equality constraint is `h(λ) = Σ_i w_i - n = 0`.

Gradient of ALM penalty term `λ_eq·h + (μ/2)·h²` w.r.t. `λ_kj` (using chain rule through `u_i`):
```
∂/∂λ_kj [λ_eq·Σw + (μ/2)(Σw-n)²] = (λ_eq + μ·(Σw-n)) · Σ_{i: g_k(i)=j} d_i·F'(u_i)
```

This is a per-(k,j)-pair sum, computed via a single O(K·n) scatter-add pass — not a nested O(n·Σcat_counts) loop.

**F'(u) for the logit link:**

`F(u) = L + (U-L)/(1+exp(-u))` → `F'(u) = (F(u)-L)·(U-F(u))/(U-L)`

where L=`st.min_weight`, U=`st.max_weight`. For the exp link (U=Inf, L=0): `F'(u) = F(u)`.

- [ ] **Step 0: Add `dF(double u) const` to LinkFn in `src/logit.hpp`**

Read `src/logit.hpp` in full. The Deville-Sarndal logit parameterization at lines 33-39 defines `F(u)` using `logit_scale = (U-L)/((U-1)*(1-L))` and `e = exp(logit_scale * u)`. Fields are `L` and `U` (not `lo_`/`hi_`).

The derivative `F'(u) = dF/du` via chain rule:
- Logit branch (both U and L finite, U≠1, L≠1): `F'(u) = logit_scale · (F(u)-L) · (U-F(u)) / (U-L)`
- Exp branch (U = ∞): `F'(u) = F(u)` (since exp link F'(u) = exp(u))

Read the actual branch condition in `F(u)` and mirror it exactly in `dF`. The formula must use `L` and `U` (actual field names):

```cpp
double dF(double u) const {
    double f = F(u);
    if (/* COPY the same if-condition from F(u) verbatim — e.g. if (exponential) */) {
        // Exp branch: F'(u) = F(u)
        return f;
    } else {
        // Logit branch: F'(u) = logit_scale * (f-L) * (U-f) / (U-L)
        double ls = (U - L) / ((U - 1.0) * (1.0 - L));  // logit_scale
        return ls * (f - L) * (U - f) / (U - L);
    }
}
```

Branch order matches `F(u)`'s ordering (exp in `if`, logit in `else`). Confirm by reading `F(u)` before writing `dF`.

Compile gate after adding: `R CMD INSTALL --preclean . 2>&1 | tail -10`

- [ ] **Step 1: Identify the gradient/objective computation function in lbfgsb_solver.cpp**

From Task 12 Step 0, identify the function that computes BOTH the objective value AND the gradient (e.g., `phi_from_u`, `phi_and_grad`, or similar). Read that function's signature and local variable names:
- What parameter/variable holds the design weights (the reviewer identified `d[]` at lbfgsb_solver.cpp:310)?
- What variable holds the log-dual iterate (`u[]` at line 316)?
- What variable holds the gradient array (actual name is `grad`, not `grad_lam`)?
- What variable holds the objective accumulator (actual name is `obj`, returned by value)?
- What variable holds category offsets (`off` at line 307)?

The ALM block must be inserted INSIDE this gradient/objective function, after the existing contribution to the gradient is computed.

- [ ] **Step 2: Add ALM gradient + objective inside the gradient computation function**

Using the actual variable names from Step 1, insert this block (adapt names to match the codebase):

```cpp
// ALM: equality constraint sum(w) = n in dual space.
if (st.alm_mu > 0.0 || st.alm_lambda != 0.0) {
    // sum_w = Σ d[i]*F(u[i]) — iterate weights; d[] and u[] are already in scope here
    double sum_w = 0.0;
    for (int i = 0; i < st.n; i++) sum_w += d[i] * fn.F(u[i]);
    double residual = sum_w - static_cast<double>(st.n);
    double alm_scale = st.alm_lambda + st.alm_mu * residual;

    // Scatter-add ALM gradient contribution: O(K*n) single pass
    // Use actual gradient array name (grad, not grad_lam)
    for (int k = 0; k < st.K; k++) {
        for (int i = 0; i < st.n; i++) {
            int g = st.group_ids[k][i];
            if (g >= 0)
                grad[off[k] + g] += alm_scale * d[i] * fn.dF(u[i]);
        }
    }

    // ALM objective term — add to the existing objective accumulator (obj, not total_f)
    obj += st.alm_lambda * residual + (st.alm_mu / 2.0) * residual * residual;
}
```

Adapt: `grad` → actual gradient array name; `obj` → actual objective accumulator name.

- [ ] **Step 3: Compile gate**

`R CMD INSTALL --preclean . 2>&1 | tail -10`

---

### Task 14: Add ALM outer loop to lbfgsb_solve

**Files:**
- Modify: `src/lbfgsb_solver.cpp` only (no new header — ALMParams is a local struct)

- [ ] **Step 1: Replace lbfgsb_solve body with ALM outer loop**

Replace the existing `lbfgsb_solve(CalibState& st)` body with the ALM wrapper. ALMParams is a local struct (no `src/alm.h`):

```cpp
LBFGSResult lbfgsb_solve(CalibState& st) {
    struct ALMParams {
        double mu0       = 10.0;
        double rho       = 10.0;   // penalty growth factor; empirically validated in Task 15
        int    max_outer = 10;
        double tol_eq    = 1e-8;
    } alm;

    double lambda_eq = 0.0;
    double mu = alm.mu0;
    LBFGSResult res;

    for (int outer = 0; outer < alm.max_outer; outer++) {
        st.alm_lambda = lambda_eq;
        st.alm_mu     = mu;
        res = lbfgsb_solve_inner(st);

        // sum_w from output weights written by lbfgsb_solve_inner
        double sum_w = 0.0;
        for (int i = 0; i < st.n; i++) sum_w += st.weights[i];
        double residual = sum_w - static_cast<double>(st.n);

        if (std::fabs(residual) < alm.tol_eq) break;
        lambda_eq += mu * residual;
        mu *= alm.rho;
    }
    return res;
}
```

Note: `st.weights[i]` here is safe because `lbfgsb_solve_inner` calls `compute_final_weights_and_error` before returning, which writes the output weights to `st.weights`. This is distinct from the inside-gradient use in Task 13 where `st.weights` still holds design weights.

- [ ] **Step 2: Compile gate**

`R CMD INSTALL --preclean . 2>&1 | tail -10`

- [ ] **Step 3: Remove single-shot clamp at lbfgsb_solver.cpp:150**

Before:
```cpp
wi = std::max(st.min_weight, std::min(st.max_weight, wi));
```

The ALM enforces `sum(w) = n`; the logit link (`src/lbfgsb_solver.cpp:306`) with `LinkFn(st.min_weight, st.max_weight)` already enforces the box constraints `[lo, hi]` by construction. The single-shot post-hoc clamp is no longer needed. Delete this line.

If the comment at lines 153-154 ("Do NOT normalize: bridge normalizes start_weights to mean=1 before rk_calibrate()") needs updating after P3 changes, update it.

- [ ] **Step 4: Compile gate**

`R CMD INSTALL --preclean . 2>&1 | tail -10`

---

### Task 15: Verify ALM convergence on known cases

**Files:**
- Modify: `tests/testthat/test-lbfgsb.R` (add ALM convergence tests)

- [ ] **Step 1: Add tests for ALM convergence**

```r
test_that("L-BFGS-B ALM converges to sum(w)=n with tight bounds", {
  set.seed(11L)
  n   <- 500L
  df  <- data.frame(x = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.3, c = 0.2))
  res <- leafblower::harvest(df, tgt, method = "lbfgsb",
                              max_weight = 1.5, min_weight = 0.2)
  expect_equal(mean(res), 1.0, tolerance = 1e-6)
  expect_true(max(res) <= 1.5 + 1e-6)
  expect_true(min(res) >= 0.2 - 1e-6)
})

test_that("L-BFGS-B ALM mean=1 without tight bounds", {
  set.seed(12L)
  n   <- 1000L
  df  <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  res <- leafblower::harvest(df, tgt, method = "lbfgsb",
                              max_weight = 100, min_weight = 0)
  expect_equal(mean(res), 1.0, tolerance = 1e-6)
})
```

- [ ] **Step 2: Add ALM stress test near infeasibility boundary**

This validates that `rho=10.0` geometric growth does not destabilize L-BFGS-B when bounds are tight. Choose a case where nearly all observations must hit `max_weight` (highly skewed targets):

```r
test_that("L-BFGS-B ALM stable near infeasibility boundary (tight upper bound)", {
  set.seed(99L)
  n   <- 300L
  df  <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE,
                                prob = c(0.9, 0.1))))
  # Target forces rare category to nearly 50% share; very few observations must
  # carry large weight → many weights at max_weight=3.0
  tgt <- list(x = c(a = 0.5, b = 0.5))
  res <- leafblower::harvest(df, tgt, method = "lbfgsb",
                              max_weight = 3.0, min_weight = 0.1)
  expect_equal(mean(res), 1.0, tolerance = 1e-5)
  expect_true(max(res) <= 3.0 + 1e-5)
  expect_true(min(res) >= 0.1 - 1e-5)
})
```

If this test FAILS or produces convergence warnings, the ALM parameters (`rho`, `max_outer`) need empirical tuning before the Phase 3 commit.

- [ ] **Step 3: Run Task 2 caps test — must PASS now**

`Rscript -e "devtools::test(filter='algo-selection')"`

---

### Task 16: Remove unconditional normalization from harvest.R

**Files:**
- Modify: `R/harvest.R:97`

- [ ] **Step 1: Replace unconditional normalize with assertion**

Before:
```r
# Normalize to mean=1 (preserves calibration constraints which are proportional)
weights <- weights / mean(weights)
```

After:
```r
# Both solvers guarantee mean(w) == 1 (iEPPA: Dykstra hyperplane; L-BFGS-B: ALM).
stopifnot(abs(mean(weights) - 1) < 1e-6)
```

- [ ] **Step 2: Run full test suite**

`Rscript -e "devtools::test()"`

Expected: 0 FAIL.

---

### Task 17: Commit Phase 3

- [ ] **Step 1: Stage and commit**

```bash
git add src/lbfgsb_solver.cpp src/logit.hpp src/types.hpp \
        R/harvest.R \
        tests/testthat/test-lbfgsb.R \
        tests/testthat/test-algo-selection.R
git commit -m "$(cat <<'EOF'
refactor(lbfgsb): ALM for sum(w)=n; remove harvest.R post-normalization

Wraps L-BFGS-B in an augmented Lagrangian outer loop enforcing
sum(w) = n as a hard equality constraint. ALM gradient is derived in
dual lambda_kj space: (lambda_eq + mu*(Σw-n)) * Σ_{i: g_k(i)=j} d_i*F'(u_i)
for each (k,j) pair — NOT a per-observation u_i gradient.

logit link (LinkFn) already enforces box bounds [lo,hi] by construction;
single-shot post-hoc clamp at lbfgsb_solver.cpp:150 removed.

harvest.R:97 unconditional normalization replaced by a sanity-check
assertion. Both solvers now guarantee mean(w) = 1 internally.
EOF
)"
```

---

## Self-Review Checklist

- [ ] All file paths verified to exist
- [ ] Each task has explicit code/command (no "TBD", no "implement similar to")
- [ ] ALMParams is a LOCAL struct inside lbfgsb_solve — no src/alm.h header file
- [ ] CalibState fields (alm_lambda, alm_mu) added to src/types.hpp (not src/leafblower.h)
- [ ] `dF(double u)` added to LinkFn in src/logit.hpp BEFORE Task 13 gradient code
- [ ] dF formula uses field names `L`/`U` (not `lo_`/`hi_`) and includes `logit_scale` factor
- [ ] dF logit branch: `logit_scale * (f-L) * (U-f) / (U-L)` where `ls = (U-L)/((U-1)*(1-L))`
- [ ] ALM gradient uses local `d[]` vector NOT `st.design_weights[]` (which doesn't exist)
- [ ] ALM gradient sum_w computed from `d[i]*fn.F(u[i])` NOT `st.weights[i]` inside gradient
- [ ] ALM gradient uses actual gradient array name (`grad`, not `grad_lam`)
- [ ] ALM objective uses actual accumulator name (`obj`, not `total_f`)
- [ ] ALM block inserted INSIDE gradient/objective function (e.g. phi_from_u), not in lbfgsb_solve
- [ ] ALM outer loop in Task 14 uses `st.weights[i]` (safe: post-compute_final_weights_and_error)
- [ ] ALM gradient O(K·n) scatter-add pass, not O(n·Σcat_counts) nested loop
- [ ] ALM gradient formula verified to be in λ_kj dual space (Task 13 Step 2)
- [ ] `lbfgsb_solve` caller grep done in Task 12 Step 1 BEFORE refactor
- [ ] `lbfgsb_solve_inner` created in Task 12 BEFORE it is called in Task 14
- [ ] `lbfgsb_solve` external signature `LBFGSResult lbfgsb_solve(CalibState&)` preserved
- [ ] RK_ALG_AUTO=0 enum definition verified to STAY in C header (Task 3 Step 0)
- [ ] Test-harvest.R:25-33 DELETE (not patch) — old assertion expects "lbfgsb", wrong after routing change
- [ ] Test-lbfgsb.R:47-50 sub-block DELETE — outer test_that (lines 29-56) must survive
- [ ] Python bindings updated in Task 4
- [ ] Benchmark artifacts in SEPARATE commit (Task 7b)
- [ ] alg_names[0] comment added (Task 4 Step 4)
- [ ] Task 8 regression guard passes on current code (not a RED test)
- [ ] Compile gate after every C++ file change
- [ ] Phase 1 ships independently (user-visible API change)
- [ ] Phase 2 independent of Phase 3 (ordering: P1 → P2 → P3)
- [ ] Phase 3 depends on Phase 2 being merged (both solvers must guarantee mean=1 before harvest.R:97 normalization is removed)
- [ ] Benchmark regression check after iEPPA change (Task 9 Step 8)

---

## Execution Notes

- **Budget**: P1 ≈ 3-4 hours, P2 ≈ 4-6 hours (Dykstra math verification), P3 ≈ 1-2 days (ALM parameter tuning + gradient verification)
- **Subagent dispatch**: Delegate to Gemini per user preference. One implementer per task; spec + code reviewer after each major task (per Task 3, 9, 12-14, 16).
- **Risk profile**: P1 low risk (mechanical), P2 moderate (math + algorithmic), P3 high (numerical tuning + gradient derivation verification).
- **ALM parameter sensitivity**: `mu0=10.0`, `rho=10.0`, `tol_eq=1e-8` are initial values; empirical tuning may be needed. Benchmark against iEPPA on the stepstone dataset to verify no accuracy regression.
