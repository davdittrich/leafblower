# Clean Code Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix one correctness bug and seven code smells identified in the clean-code review of the bounded-convergence fix.

**Architecture:** All changes are internal to three files. No API changes, no new files. Task 1 restructures execution order in `rk_calibrate` (algorithm selection before validation) to give the guard algorithm-level precision. Tasks 2–4 are pure refactors with no behaviour change — tests must stay green.

**Tech Stack:** C++17, R 4.x, testthat 3.x. Compile gate: `R CMD INSTALL --preclean .`. Test gate: `Rscript -e "devtools::test()"`.

---

## File Map

| File | What changes |
|------|-------------|
| `src/c_api.cpp` | Task 1: reorder select_algorithm/validate_inputs; update validate_inputs signature; fix guard condition. Task 4: kComplexityThreshold constant; routing verbose message fix |
| `src/ieppa.cpp` | Task 2: combine normalization loops; rename Wf→Wsum; rename infeas_flag→is_infeasible. Task 3: hoist bucket/scale out of inner loop. Task 4: named constants |
| `tests/testthat/test-lbfgsb.R` | Task 1: update singularity guard test to cover ieppa acceptance. Task 4: remove stale comment |

---

## Task 1: Fix Singularity Guard Scope (Bug 1)

**The bug:** `validate_inputs` fires the logit singularity guard for ALL calls with finite `max_weight`, including `method="ieppa"`. iEPPA uses no logit link — it uses direct multiplicative scaling. `harvest(df, tgt, method="ieppa", max_weight=1.0+5e-7)` currently returns `RK_ERR_BADARG` with "logit link undefined". It should succeed.

**Root cause:** `validate_inputs` runs before `select_algorithm`. The guard predicate `std::isfinite(p->max_weight)` conflates "max_weight is finite" with "L-BFGS-B will be used". These are different after the routing fix.

**Fix:** Call `select_algorithm` before `validate_inputs`. Pass the resolved algorithm to `validate_inputs`. Guard fires only when `alg == RK_ALG_LBFGSB`.

**Files:**
- Modify: `src/c_api.cpp:35-123` (`validate_inputs`), `src/c_api.cpp:137-227` (`rk_calibrate`)
- Modify: `tests/testthat/test-lbfgsb.R:30-45`

- [ ] **Step 1: Write the failing test**

Add the following test to `tests/testthat/test-lbfgsb.R` (replace the existing `"near-one max_weight rejected with informative error"` block at line 30):

```r
test_that("near-one max_weight rejected for lbfgsb, accepted for ieppa", {
  set.seed(1)
  df <- data.frame(x = factor(sample(c("a","b"), 200, replace=TRUE)))
  tgt <- list(x = c(a=0.5, b=0.5))

  # lbfgsb: exact and near-1 max_weight → logit singularity error
  expect_error(harvest(df, tgt, method="lbfgsb", max_weight=1.0),
               regexp="logit link undefined")
  expect_error(harvest(df, tgt, method="lbfgsb", max_weight=1.0 + 5e-7),
               regexp="logit link undefined")
  expect_error(harvest(df, tgt, method="lbfgsb", max_weight=1.0 - 5e-7),
               regexp="logit link undefined")

  # ieppa: near-1 max_weight → valid (no logit link); may or may not converge
  expect_no_error(suppressWarnings(
    harvest(df, tgt, method="ieppa", max_weight=1.0 + 5e-7)
  ))

  # auto: near-1 max_weight → valid (routes to ieppa)
  expect_no_error(suppressWarnings(
    harvest(df, tgt, method="auto", max_weight=1.0 + 5e-7)
  ))

  # lbfgsb: max_weight=2.0 (well outside eps) → valid
  expect_no_error(suppressWarnings(harvest(df, tgt, method="lbfgsb", max_weight=2.0)))
})
```

- [ ] **Step 2: Run test to verify the new ieppa/auto assertions currently fail**

```bash
Rscript -e "devtools::test(filter='lbfgsb')" 2>&1 | tail -10
```

Expected: `FAIL 1` (the two new `expect_no_error` assertions for ieppa/auto fail with "logit link undefined").

- [ ] **Step 3: Update `validate_inputs` signature in `src/c_api.cpp`**

Add `rk_algorithm_t alg` as the last parameter to `validate_inputs` (line 35). Replace the existing logit singularity block (lines 60–71):

```cpp
static int validate_inputs(int n, int K,
                            const double* weights,
                            const int32_t** group_ids,
                            const int* cat_counts,
                            const double** targets,
                            const rk_params_t* p,
                            rk_result_t* result,
                            rk_algorithm_t alg) {
```

Replace the singularity guard block (currently lines 60–71) with:

```cpp
    // Logit singularity guard — only applies to L-BFGS-B.
    // iEPPA uses multiplicative IPF scaling; no logit link, no singularity.
    // logit_scale = (U-L)/((U-1)*(1-L)); denominator → 0 when L or U near 1,
    // producing logit_scale ~1/eps → overflow/cancellation in L-BFGS-B.
    if (alg == RK_ALG_LBFGSB) {
        const double kSingularityEps = 1e-6;
        if (std::fabs(p->min_weight - 1.0) < kSingularityEps)
            return err("logit link undefined: min_weight near 1 makes denominator (1-L)~0");
        if (std::fabs(p->max_weight - 1.0) < kSingularityEps)
            return err("logit link undefined: max_weight near 1 makes denominator (U-1)~0");
    }
```

- [ ] **Step 4: Reorder select_algorithm / validate_inputs in `rk_calibrate`**

In `rk_calibrate` (currently line 156 calls `validate_inputs`, line 160 calls `select_algorithm`), restructure so algorithm is resolved first:

Replace lines 156–160:
```cpp
    int rc = validate_inputs(n, K, weights, group_ids, cat_counts, targets, p, result);
    if (rc != RK_OK) return rc;

    int64_t complexity = INT64_C(0);
    rk_algorithm_t alg = select_algorithm(n, K, cat_counts, p, complexity);
```

With:
```cpp
    // Resolve algorithm before validation so the singularity guard knows which
    // link function will be used. Guard against null cat_counts or invalid K/n
    // — validate_inputs will reject those cases with a proper error message.
    int64_t complexity = INT64_C(0);
    rk_algorithm_t alg = (cat_counts && K > 0 && n > 0)
        ? select_algorithm(n, K, cat_counts, p, complexity)
        : p->algorithm;

    int rc = validate_inputs(n, K, weights, group_ids, cat_counts, targets, p, result, alg);
    if (rc != RK_OK) return rc;
```

- [ ] **Step 5: Compile**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected last line: `* DONE (leafblower)`

- [ ] **Step 6: Run tests to verify new assertions pass and nothing regressed**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 42 ]` (41 previous + 2 new ieppa/auto assertions)

- [ ] **Step 7: Commit**

```bash
git add src/c_api.cpp tests/testthat/test-lbfgsb.R
git commit -m "fix(c_api): scope logit singularity guard to L-BFGS-B only

validate_inputs fired the guard for all finite max_weight, rejecting
valid iEPPA calls with a spurious 'logit link undefined' error. iEPPA
uses multiplicative IPF scaling, not a logit link.

Fix: resolve algorithm via select_algorithm before calling
validate_inputs; pass alg parameter; guard fires only for
RK_ALG_LBFGSB."
```

---

## Task 2: ieppa.cpp — Naming and Logic Cleanup

**Items covered:** Smell 2 (combine normalization loops), Smell 3 (Wf→Wsum), Smell 8 (infeas_flag→is_infeasible). Pure rename/refactor — no behaviour change.

**Files:**
- Modify: `src/ieppa.cpp`

- [ ] **Step 1: Rename `infeas_flag` → `is_infeasible` throughout `ieppa_solve`**

In `ieppa_solve`, replace all occurrences of `infeas_flag`:

```cpp
// Line ~49: declaration
bool is_infeasible = false;

// Line ~71: set flag
if (Tkj > 0.0) is_infeasible = true;

// Line ~115: status assignment
res.status = is_infeasible ? RK_ERR_INFEAS : RK_OK;
```

- [ ] **Step 2: Rename `Wf` → `Wsum` in the fixup loop**

In the post-convergence fixup loop (around line 127):

```cpp
// Before:
double Wf = 0.0;
for (int i = 0; i < st.n; i++) Wf += w[i];
double wm = (Wf > 1e-300) ? Wf / st.n : 1.0;

// After:
double Wsum = 0.0;
for (int i = 0; i < st.n; i++) Wsum += w[i];
double wm = (Wsum > 1e-300) ? Wsum / st.n : 1.0;
```

- [ ] **Step 3: Combine the two normalization loops into one**

In the main iteration loop, the normalization block (around lines 88–93):

```cpp
// Before:
if (wm > 1e-300) for (int i = 0; i < st.n; i++) w[i] /= wm;
// Scale the Dykstra box correction to stay consistent after renormalization.
if (wm > 1e-300) for (int i = 0; i < st.n; i++) q[i] /= wm;

// After:
if (wm > 1e-300) {
    for (int i = 0; i < st.n; i++) { w[i] /= wm; q[i] /= wm; }
}
```

- [ ] **Step 4: Compile**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 5: Run tests to verify no regression**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 42 ]`

- [ ] **Step 6: Commit**

```bash
git add src/ieppa.cpp
git commit -m "refactor(ieppa): rename infeas_flag, Wf; combine normalization loops

- is_infeasible replaces infeas_flag (_flag suffix is noise)
- Wsum replaces Wf (Wf is meaningless; matches main-loop naming)
- Two sequential w[]/q[] normalization loops merged into one pass"
```

---

## Task 3: ieppa.cpp — Hoist Inner-Loop Allocations

**Item covered:** Smell 4. `bucket` and `scale` are `std::vector` allocated on the heap inside the IPF loop — O(K × max_iter) heap allocations. Hoist to outer scope; use `std::fill` to reset per-k. No behaviour change.

**Files:**
- Modify: `src/ieppa.cpp`

- [ ] **Step 1: Hoist `bucket` and `scale` above the iteration loop**

`max_cats` is the maximum of `st.cat_counts[0..K-1]`. Allocate once with that capacity; reset with `std::fill` inside the k loop.

In `ieppa_solve`, before the `for (int iter = 1; ...)` loop, add:

```cpp
    int max_cats = *std::max_element(st.cat_counts, st.cat_counts + st.K);
    std::vector<double> bucket(max_cats), scale(max_cats);
```

Inside the k loop, replace:

```cpp
    // Before (allocates on every k iteration):
    std::vector<double> bucket(st.cat_counts[k], 0.0);
    ...
    std::vector<double> scale(st.cat_counts[k], 1.0);
```

With:

```cpp
    // After (reset pre-allocated buffer):
    std::fill(bucket.begin(), bucket.begin() + st.cat_counts[k], 0.0);
    ...
    std::fill(scale.begin(), scale.begin() + st.cat_counts[k], 1.0);
```

The rest of the bucket accumulation and scale application code is unchanged — indexing `bucket[g]` and `scale[j]` still use `g < cat_counts[k]` and `j < cat_counts[k]`, so the oversized buffer causes no issue.

- [ ] **Step 2: Compile**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 3: Run tests to verify no regression**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 42 ]`

- [ ] **Step 4: Commit**

```bash
git add src/ieppa.cpp
git commit -m "perf(ieppa): hoist bucket/scale allocation outside IPF loop

Eliminated O(K * max_iter) heap allocations by pre-allocating bucket
and scale with max_cats capacity and resetting with std::fill per
margin. No behaviour change."
```

---

## Task 4: Named Constants + Minor Cleanup

**Items covered:** Smell 5 (magic numbers in ieppa.cpp and c_api.cpp), Smell 6 (verbose routing message fires for explicit selection), Smell 7 (stale TDD comment in test).

**Files:**
- Modify: `src/ieppa.cpp`
- Modify: `src/c_api.cpp`
- Modify: `tests/testthat/test-lbfgsb.R`

- [ ] **Step 1: Name magic constants in `ieppa_solve`**

At the top of `ieppa_solve` (before any local variables), add:

```cpp
    static constexpr double kEmptyBucketThreshold  = 1e-15;  // bucket near-zero: skip scale
    static constexpr double kWeightCollapseThreshold = 1e-300; // weights collapsed: skip norm
    static constexpr int    kMaxFixupIterations     = 20;    // post-convergence fixup cap
```

Replace usages:
- `1e-15` (line ~70) → `kEmptyBucketThreshold`
- `1e-300` (lines ~91, ~129) → `kWeightCollapseThreshold`
- `20` (line ~126) → `kMaxFixupIterations`

- [ ] **Step 2: Name complexity threshold in `c_api.cpp`**

Before `select_algorithm`, add a file-scope constant:

```cpp
static constexpr int64_t kComplexityThreshold = 500000L;
```

Replace `500000L` in `select_algorithm` (line 132):

```cpp
    if (complexity_out > kComplexityThreshold || std::isfinite(p->max_weight) || p->min_weight > 0.0)
```

- [ ] **Step 3: Fix verbose routing message to only fire for auto-routing**

In `rk_calibrate`, the verbose block (around lines 183–192):

```cpp
// Before:
if (p->verbose >= 1) {
    char msg[256];
    if (alg == RK_ALG_IEPPA)
        snprintf(msg, 256, "Auto-selected iEPPA: complexity=%lld, max_weight=%.2f, min_weight=%.2f",
                 (long long)complexity, p->max_weight, p->min_weight);
    else
        snprintf(msg, 256, "Auto-selected L-BFGS-B: complexity=%lld <= 500000, max_weight=%.2f, min_weight=%.2f",
                 (long long)complexity, p->max_weight, p->min_weight);
    st.log(msg);
}

// After:
if (p->verbose >= 1 && p->algorithm == RK_ALG_AUTO) {
    char msg[256];
    if (alg == RK_ALG_IEPPA)
        snprintf(msg, 256, "Auto-selected iEPPA: complexity=%lld, max_weight=%.2f, min_weight=%.2f",
                 (long long)complexity, p->max_weight, p->min_weight);
    else
        snprintf(msg, 256, "Auto-selected L-BFGS-B: complexity=%lld <= 500000, max_weight=%.2f, min_weight=%.2f",
                 (long long)complexity, p->max_weight, p->min_weight);
    st.log(msg);
}
```

- [ ] **Step 4: Remove stale TDD comment from `tests/testthat/test-lbfgsb.R`**

Delete line 13:
```r
  # RED: harvest() not yet implemented
```

The test block `"L-BFGS-B converges on 3-margin no-bounds case"` (lines 1–18) after removal:

```r
test_that("L-BFGS-B converges on 3-margin no-bounds case", {
  set.seed(42)
  n <- 50000L
  age <- sample(c("18-34","35-54","55+"), n, replace=TRUE, prob=c(0.35,0.40,0.25))
  sex <- sample(c("M","F"), n, replace=TRUE, prob=c(0.52,0.48))
  edu <- sample(c("HS","College","Grad"), n, replace=TRUE, prob=c(0.40,0.40,0.20))
  df  <- data.frame(age=factor(age), sex=factor(sex), edu=factor(edu))
  tgt <- list(
    age = c("18-34"=0.30, "35-54"=0.45, "55+"=0.25),
    sex = c(M=0.50, F=0.50),
    edu = c(HS=0.35, College=0.45, Grad=0.20)
  )
  result <- harvest(df, tgt, method="lbfgsb")
  expect_s3_class(result, "data.frame")
  expect_true("weights" %in% names(result))
  expect_lt(abs(mean(result$weights) - 1.0), 1e-8)
})
```

- [ ] **Step 5: Compile**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Expected: `* DONE (leafblower)`

- [ ] **Step 6: Run tests to verify no regression**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 42 ]`

- [ ] **Step 7: Commit**

```bash
git add src/ieppa.cpp src/c_api.cpp tests/testthat/test-lbfgsb.R
git commit -m "refactor: name magic constants; fix routing verbose; drop stale comment

- ieppa.cpp: kEmptyBucketThreshold, kWeightCollapseThreshold,
  kMaxFixupIterations replace bare literals
- c_api.cpp: kComplexityThreshold replaces 500000L; routing verbose
  message now guarded by p->algorithm==RK_ALG_AUTO so it doesn't
  fire when the user explicitly selects an algorithm
- test-lbfgsb.R: remove stale '# RED: harvest() not yet implemented'"
```

---

## Self-Review

**Spec coverage:**
- Bug 1 (guard too broad) → Task 1 ✓
- Smell 2 (double loop) → Task 2 ✓
- Smell 3 (Wf) → Task 2 ✓
- Smell 4 (inner alloc) → Task 3 ✓
- Smell 5 (magic numbers, both files) → Task 4 ✓
- Smell 6 (routing message) → Task 4 ✓
- Smell 7 (stale comment) → Task 4 ✓
- Smell 8 (infeas_flag) → Task 2 ✓

**Placeholder scan:** No TBDs, no vague steps. Every code block is complete.

**Type consistency:**
- `is_infeasible` used consistently in Task 2 (declaration, set, read).
- `Wsum` used consistently in Task 2 (declaration, accumulation, condition).
- `kEmptyBucketThreshold`, `kWeightCollapseThreshold`, `kMaxFixupIterations`, `kComplexityThreshold` each defined once and referenced in exactly one location.
- `bucket` and `scale` hoisted in Task 3; subsequent tasks reference them as already-declared.
- Task 1 `alg` parameter: added to signature and call site consistently.
