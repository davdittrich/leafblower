# Clean Code Fixes (Round 4) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix one confirmed correctness bug (infeasibility status never returned by iEPPA solver) and clean up three minor documentation gaps identified by round-4 code review.

**Architecture:** Changes span one C++ solver (`src/ieppa.cpp`), one C API init function (`src/c_api.cpp`), and one R helper (`R/design_effect.R`). Tasks are independent; each compiles and passes all tests before the next begins.

**Tech Stack:** C++17, R, testthat, `R CMD INSTALL --preclean .`.

**Beads issues:** leafblower-9o9 (T1), leafblower-dm4 (T2).

---

## File Map

| File | Changes |
|---|---|
| `tests/testthat/test-harvest.R` | Add infeasibility test (Task 1) |
| `src/ieppa.cpp` | Fix infeasibility status after main loop (Task 1); improve kEmptyBucketThreshold comment (Task 2) |
| `src/c_api.cpp` | Fix epsilon initialization 0.05 → 0.0 (Task 2) |
| `R/design_effect.R` | Add denominator comment (Task 2) |

---

## Task 1: Fix iEPPA infeasibility status

**Beads:** leafblower-9o9

**Context:** `ieppa_solve` initializes `res.status = RK_ERR_NOCONV` (line 46). It sets `is_infeasible = true` when `bucket[j] < kEmptyBucketThreshold * W` and `target[j] > 0` (line 83). `RK_ERR_INFEAS` is returned ONLY from line 132 inside `if (errRp < st.tol_abs)`. But infeasible problems can never satisfy `errRp < tol_abs` — the empty bucket always contributes `|0/W − τ_kj| = τ_kj > 0` to errRp. Therefore the convergence break is unreachable for truly infeasible inputs, and `is_infeasible` can never cause `RK_ERR_INFEAS` to be returned. All infeasible problems silently report `RK_ERR_NOCONV`, and the R layer emits a "did not converge" warning instead of the correct "infeasible problem" stop.

**Fix:** After the main loop and before the fixup loop, check `is_infeasible` and override the status.

**Files:**
- Modify: `tests/testthat/test-harvest.R` (new failing test)
- Modify: `src/ieppa.cpp:135` (add status override after main loop)

- [ ] **Step 1: Write the failing test in `tests/testthat/test-harvest.R`**

Append to `tests/testthat/test-harvest.R`:

```r
test_that("infeasible problem (empty cell with positive target) stops with infeasible error", {
  # All 3 observations are in group "a"; no observation is in group "b".
  # Target for "b" = 0.4 > 0 → empty cell with positive target → infeasible.
  df  <- data.frame(x = factor(c("a", "a", "a"), levels = c("a", "b")))
  tgt <- list(x = c(a = 0.6, b = 0.4))
  expect_error(
    harvest(df, tgt, method = "ieppa"),
    regexp = "infeasible problem"
  )
})
```

Note: `factor(..., levels = c("a", "b"))` ensures both levels are present in the factor so R encodes the column correctly; the data has no "b" observations, creating the empty bucket.

- [ ] **Step 2: Run test to verify it fails**

```bash
Rscript -e "devtools::test()" 2>&1 | grep -E "FAIL|infeasible|PASS" | tail -10
```

Expected: The new test FAILS — currently returns a "did not converge" warning (status 1) instead of an "infeasible problem" stop (status 2).

- [ ] **Step 3: Fix `src/ieppa.cpp`**

Read `src/ieppa.cpp`. Find the closing brace of the main loop at line 135 (the `}` that ends `for (int iter = 1; iter <= st.inner_max_iter; iter++)`). Insert this block immediately after:

```cpp
    // Infeasibility detected during iteration: override NOCONV with INFEAS.
    // Truly infeasible problems (empty bucket with positive target) can never
    // converge — the empty bucket always contributes τ_kj > 0 to errRp,
    // so errRp never drops below tol_abs and the convergence break never fires.
    // Check the flag here and return the correct status code.
    if (is_infeasible && res.status == RK_ERR_NOCONV)
        res.status = RK_ERR_INFEAS;
```

The insertion point is between the closing `}` of the main loop and the comment `// Final normalization-and-clamp fixup:`. The fixup loop runs regardless of infeasibility (it normalizes weights that were produced before infeasibility was detected).

After the insertion, the block order is:
1. Main iteration loop → `break` sets `RK_OK` or `RK_ERR_INFEAS`; loop exhaustion leaves `RK_ERR_NOCONV`
2. **NEW:** Infeasibility override: if `is_infeasible && NOCONV` → set `RK_ERR_INFEAS`
3. Fixup loop
4. Copy weights to `st.weights`
5. Return `res`

- [ ] **Step 4: Compile gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

Expected: `* DONE (leafblower)`

- [ ] **Step 5: Run all tests**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 49 ]` (48 existing + 1 new)

- [ ] **Step 6: Commit**

```bash
git add src/ieppa.cpp tests/testthat/test-harvest.R
git commit -m "fix(ieppa): return RK_ERR_INFEAS for infeasible problems (empty cell + positive target)"
```

---

## Task 2: Documentation and naming cleanup

**Beads:** leafblower-dm4

**Context:** Three minor documentation gaps:

1. `kEmptyBucketThreshold = 1e-15` in `src/ieppa.cpp:41`. Comment says "bucket near-zero: skip scale" but the threshold is relative to W (`bucket[j] < 1e-15 * W`), not absolute. Maintainers reading the constant may expect it to be an absolute threshold.

2. `p->epsilon = 0.05` in `src/c_api.cpp:29`. Field is deprecated ("no longer read by any solver; kept for ABI compat" per `leafblower.h:37`) but is initialized to 0.05 rather than 0.0. A deprecated no-op field should have a neutral default.

3. `var(outcome)` in `R/design_effect.R:25`. Uses R's `var()` which divides by n−1 (sample variance). The Henry-Valliant calibration design effect formula compares weighted variance to unweighted variance; using n−1 in both is the correct survey-methodology convention (both are unbiased estimators of the same population variance), but no comment documents this choice.

**Files:**
- Modify: `src/ieppa.cpp:41` (kEmptyBucketThreshold comment)
- Modify: `src/c_api.cpp:29` (epsilon init)
- Modify: `R/design_effect.R:25` (add denominator comment)

- [ ] **Step 1: Improve `kEmptyBucketThreshold` comment in `src/ieppa.cpp`**

Read `src/ieppa.cpp`. Change line 41 from:
```cpp
    static constexpr double kEmptyBucketThreshold   = 1e-15;   // bucket near-zero: skip scale
```
to:
```cpp
    static constexpr double kEmptyBucketThreshold   = 1e-15;   // relative threshold: bucket[j] < 1e-15*W → treat as empty, skip IPF scale
```

- [ ] **Step 2: Fix `epsilon` initialization in `src/c_api.cpp`**

Read `src/c_api.cpp`. Change line 29 from:
```cpp
    p->epsilon       = 0.05;
```
to:
```cpp
    p->epsilon       = 0.0;   /* deprecated: not read by any solver; see leafblower.h */
```

- [ ] **Step 3: Add denominator comment to `R/design_effect.R`**

Read `R/design_effect.R`. Change line 25 from:
```r
  var_u   <- var(outcome)
```
to:
```r
  var_u   <- var(outcome)  # n-1 denominator (sample variance): correct per H&V (both num/denom use n-1 → ratio is unbiased)
```

- [ ] **Step 4: Compile gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

Expected: `* DONE (leafblower)`

- [ ] **Step 5: Run all tests**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 49 ]`

- [ ] **Step 6: Commit**

```bash
git add src/ieppa.cpp src/c_api.cpp R/design_effect.R
git commit -m "docs: clarify kEmptyBucketThreshold is relative; fix epsilon default to 0.0; document design_effect n-1 denominator"
```

---

## Self-Review

### 1. Spec coverage

| Review finding | Task |
|---|---|
| Infeasibility status bug (is_infeasible never triggers RK_ERR_INFEAS) | Task 1 |
| No infeasibility test | Task 1 |
| kEmptyBucketThreshold comment misleading (relative vs absolute) | Task 2 |
| epsilon initialized to 0.05 (deprecated, not read) | Task 2 |
| design_effect n-1 denominator undocumented | Task 2 |

All 5 findings addressed.

### 2. Placeholder scan

No TBDs, vague steps, or "similar to" references. All code blocks are exact.

### 3. Type consistency

- Task 1: `res.status` is `int` (from IEPPAResult); `RK_ERR_INFEAS` is `int` constant. Assignment type-safe.
- Task 2: Comment-only and scalar default change — no type concerns.

### 4. Correctness of the infeasibility fix

The fix adds a post-loop check that only fires when both conditions hold:
- `is_infeasible`: empty cell with positive target detected during iteration
- `res.status == RK_ERR_NOCONV`: loop exited without converging (as expected for truly infeasible problems)

If a problem happens to converge before the flag is set (impossible for truly infeasible problems, but defensive), the convergence break at line 132 would already have set the status to `RK_OK` or `RK_ERR_INFEAS` (impossible to reach `RK_OK` if `is_infeasible` is true) — the post-loop check would not fire since `res.status != RK_ERR_NOCONV`. The fix is safe and conservative.
