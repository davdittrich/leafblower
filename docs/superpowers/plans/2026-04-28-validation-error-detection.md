# Validation & Error Detection Gaps — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Close all pre-flight validation holes so invalid inputs produce RK_ERR_BADARG/RK_ERR_INFEAS immediately.

**Architecture:** Four fixes across validation.hpp, calib_validate.cpp, and cell_table.cpp. B5 (iEPPA/raking validation) deferred to existing ticket leafblower-0991. B6 closed as phantom — already handled (see note below).

**Tech Stack:** C++17, R, devtools::test(), R CMD INSTALL --preclean .

---

**Mechanism:** Unified tolerance + pre-flight checks + R-layer input validation
**Forbidden:** Changing calib_validate.cpp tolerance (it's the reference), adding checks inside hot solver loops
**Audit:** Each bug: RED test demonstrating wrong behavior, then GREEN after fix

---

## Task 1 — B3: Unify target-sum tolerance to 1e-6

**File:** `src/validation.hpp` line 84

**Bug mechanism:** `validate_calibrate_inputs` rejects targets summing to 1.0 ± (1e-9 … 1e-8) with `RK_ERR_BADARG`, but `calib_validate_preentry` in `calib_validate.cpp` line 81 uses `1e-6` as the authoritative tolerance. Any target vector that passes `calib_validate_preentry` but was first validated through `validate_calibrate_inputs` at the bridge would be rejected at the tighter 1e-8 gate even though the solver-level check would accept it.

**Exact current code (validation.hpp line 84):**
```cpp
if (std::fabs(sum - 1.0) > 1e-8)
    return err("targets[k] does not sum to 1 (within 1e-8)");
```

**Exact fix (validation.hpp line 84):**
```cpp
if (std::fabs(sum - 1.0) > 1e-6)
    return err("targets[k] does not sum to 1 (within 1e-6)");
```

**Unchanged:** Everything else in `validate_calibrate_inputs`. The `calib_validate.cpp` tolerance stays at `1e-6` (reference).

**Regressions prevented:** No solver sees a target vector with sum error > 1e-6 — the `calib_validate_preentry` check remains as the inner guard.

---

### RED test (add to `tests/testthat/test-harvest.R`)

```r
test_that("B3: target sum 1.0 + 5e-7 is accepted after tolerance unification", {
  # 5e-7 is between old gate (1e-8) and new gate (1e-6).
  # Before fix: harvest() errors with "does not sum to 1 (within 1e-8)".
  # After fix:  harvest() accepts and converges normally.
  set.seed(1)
  n  <- 500L
  df <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE)))
  # Deliberately add 5e-7 off-sum to hit the gap zone.
  tgt <- list(x = c(a = 0.5 + 5e-7, b = 0.5))  # sum = 1.0000005
  # Before fix this errors; after fix it should succeed (converge or BUDGET/STALL).
  expect_no_error(
    harvest(df, tgt, method = "raking", convergence = list(absolute = 1e-4))
  )
})

test_that("B3: target sum 1.0 + 2e-6 is still rejected (outside 1e-6 gate)", {
  # 2e-6 > 1e-6 — must be rejected both before and after the fix.
  set.seed(1)
  n  <- 500L
  df <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5 + 2e-6, b = 0.5))  # sum = 1.000002
  expect_error(
    harvest(df, tgt, method = "raking", convergence = list(absolute = 1e-4)),
    regexp = "does not sum to 1"
  )
})
```

**Build & verify:**
```bash
R CMD INSTALL --preclean /home/dd/Gemini/leafblower && \
  Rscript -e "devtools::test('/home/dd/Gemini/leafblower', filter='harvest')"
```

---

## B6 — CLOSED: already handled by `parse_convergence()`

**Finding (feasibility review):** `parse_convergence()` in `R/harvest.R` lines 483–509 calls `match.arg()` for all three fields before returning `conv`:

- `metric` — line 483: `match.arg(metric_raw, c("max_err", "mean_err", "kl", "chi2", "grake_norm", "l1_weight", "marginal_kl", "pct"))`
- `rule` — line 494: `match.arg(rule_raw, c("threshold", "improvement", "plateau"))`
- `stop_when` — line 509: `match.arg(convergence[["stop_when"]] %||% "any", c("any", "all"))`

`match.arg()` raises a clear R error on any unknown name before the `[[]]` lookup at lines 264–266 is ever reached. The NULL→NA→`NA_integer_`→C UB path does not exist. Calling `harvest(df, tgt, convergence=list(metric="bad_metric"))` already produces `"'arg' should be one of ..."` — no fix required.

**No task. No test. No code change.**

---

## Task 2 — B13: NA-only category INFEAS detection

**File:** `src/validation.hpp` — append new check inside `validate_calibrate_inputs`, after the existing group_ids range loop (after line 98 of the current file).

**Bug mechanism:** If all observations for category `j` of margin `k` have `group_ids[k][i] == -1` (NA), the calibration target `targets[k][j] > 0` is structurally infeasible — no observation can be assigned to that category. The current validation only checks `group_ids[k][i] >= cat_counts[k]` (out-of-range), not zero-population cells with positive targets. The solver silently diverges or hits STALL/INFEAS later with no diagnostic.

**Exact new code** — insert immediately after the closing `}` of the group_ids range loop (after `if (g >= cat_counts[k]) return err(...)`):

```cpp
    // B13: NA-only category with positive target → structural INFEAS.
    // For each (k, j): if targets[k][j] > 0 but count of obs with group_ids[k][i]==j
    // is zero, no solution exists. O(n*K) — acceptable at pre-flight time only.
    for (int k = 0; k < K; k++) {
        for (int j = 0; j < cat_counts[k]; j++) {
            if (targets[k][j] <= 0.0) continue;  // zero target: ignore
            int count = 0;
            for (int i = 0; i < n; i++)
                if (group_ids[k][i] == j) count++;
            if (count == 0) {
                char msg[256];
                std::snprintf(msg, sizeof(msg),
                    "margin %d, category %d: target=%.6g but 0 observations assigned "
                    "(all NA or missing); problem is structurally infeasible",
                    k, j, targets[k][j]);
                if (result) {
                    result->status = RK_ERR_INFEAS;
                    std::snprintf(result->message, 256, "%s", msg);
                }
                return RK_ERR_INFEAS;
            }
        }
    }
```

**Unchanged:** All existing checks in `validate_calibrate_inputs`. `calib_validate.cpp` untouched.

**Regressions prevented:** The check only fires when `targets[k][j] > 0.0` and count is exactly zero; valid sparse data with some NA observations but at least one non-NA per targeted category is unaffected.

---

### RED test (add to `tests/testthat/test-harvest.R`)

```r
test_that("B13: NA-only category with positive target returns INFEAS error", {
  # 2-margin problem. All obs for category 'a' of margin x have x=NA (-1).
  # targets[x]['a'] = 0.3 > 0 but count('a') == 0 → structurally infeasible.
  # Before fix: solver runs, may STALL or give wrong weights, no clear error.
  # After fix:  harvest() stops immediately with an error mentioning "infeasible".
  n  <- 100L
  df <- data.frame(
    x = factor(c(rep(NA,  50L), rep("b", 50L)), levels = c("a", "b")),
    y = factor(sample(c("p", "q"), n, replace = TRUE))
  )
  tgt <- list(
    x = c(a = 0.3, b = 0.7),   # 'a' has 0 obs with x='a'
    y = c(p = 0.5, q = 0.5)
  )
  expect_error(
    harvest(df, tgt, method = "raking", convergence = list(absolute = 1e-4)),
    regexp = "infeasible|INFEAS",
    ignore.case = TRUE
  )
})

test_that("B13: partial NA (some obs assigned) does not trigger INFEAS", {
  # Only 10 obs assigned to 'a', but > 0, so no INFEAS.
  n  <- 100L
  df <- data.frame(
    x = factor(c(rep("a", 10L), rep(NA, 40L), rep("b", 50L)), levels = c("a", "b"))
  )
  tgt <- list(x = c(a = 0.3, b = 0.7))
  expect_no_error(
    harvest(df, tgt, method = "raking", convergence = list(absolute = 1e-3))
  )
})
```

**Build & verify:**
```bash
R CMD INSTALL --preclean /home/dd/Gemini/leafblower && \
  Rscript -e "devtools::test('/home/dd/Gemini/leafblower', filter='harvest')"
```

---

## Task 3 — R4: estimate_M_cell K > 8 fix

**File:** `src/cell_table.cpp` lines 142–147

**Bug mechanism:** When `K > 8`, `estimate_M_cell` returns `n` unconditionally (line 146: `if (K > 8) return n;`). For `K = 9` margins with 2 categories each, the true upper bound on distinct cells is `2^9 = 512`. Returning `n = 10000` makes `M_cell_est * 10 = 100000 > n * 9 = 90000` → `use_raking = true`, routing AUTO to raking. The correct product estimate is `512`, so `512 * 10 = 5120 < 90000` → `use_raking = false` → iEPPA. The current code always picks raking for K > 8.

**Exact current code (cell_table.cpp lines 142–148):**
```cpp
int estimate_M_cell(int n, int K,
                    const int32_t* const* group_ids,
                    const int* cat_counts) {
    if (K <= 0 || n <= 0) return 0;
    if (K > 8) return n;  // can't pack; assume incompressible
    if (!pack_key_fits(K, cat_counts)) return n;
```

**Exact fix:** Replace lines 146–147 with an int64_t product-of-cat-counts estimate capped at n:

```cpp
int estimate_M_cell(int n, int K,
                    const int32_t* const* group_ids,
                    const int* cat_counts) {
    if (K <= 0 || n <= 0) return 0;
    if (K > 8) {
        // Bit-packing unavailable, but the cross-product upper bound is still
        // a better estimate than n. Use int64_t to avoid 32-bit overflow.
        int64_t prod = 1;
        for (int k = 0; k < K; k++) {
            prod *= static_cast<int64_t>(cat_counts[k]);
            if (prod >= static_cast<int64_t>(n)) return n;  // saturate early
        }
        return static_cast<int>(prod);
    }
    if (!pack_key_fits(K, cat_counts)) return n;
```

**Unchanged:** The `pack_key_fits` path for K ≤ 8, the hash-set enumeration, and all callers.

**Regressions prevented:** The early `if (prod >= n) return n` saturation preserves the existing behavior for any K > 8 problem where the true cross-product exceeds n (e.g., K=9 with 4 cats each: 4^9 = 262144 > any reasonable n).

---

### RED test (add to `tests/testthat/test-harvest.R`)

```r
test_that("R4: K=9 with 2 cats each routes to iEPPA via AUTO, not raking", {
  # K=9, cat_counts all 2 → cross-product = 512.
  # n = 10000, threshold: 512/10000 = 0.0512 < 0.9 → iEPPA.
  # Before fix: estimate returns n=10000, ratio=1.0 > 0.9 → raking.
  # After fix:  estimate returns 512, ratio=0.0512 < 0.9 → iEPPA.
  set.seed(42)
  n <- 10000L
  cats <- paste0("m", 1:9)
  df_list <- lapply(cats, function(nm)
    factor(sample(c("a", "b"), n, replace = TRUE)))
  names(df_list) <- cats
  df  <- as.data.frame(df_list)
  tgt <- setNames(
    lapply(cats, function(nm) c(a = 0.5, b = 0.5)),
    cats
  )
  result <- harvest(df, tgt, method = "auto",
                    convergence = list(absolute = 1e-3))
  alg <- attr(result, "algorithm")
  expect_identical(alg, "ieppa",
    info = paste("Expected ieppa but got", alg,
                 "— estimate_M_cell K>8 fix may be missing"))
})

test_that("R4: K=9 with 4 cats each saturates to n (product > n)", {
  # 4^9 = 262144 > n=10000; estimate should saturate to n=10000 → raking.
  set.seed(42)
  n <- 10000L
  cats <- paste0("m", 1:9)
  df_list <- lapply(cats, function(nm)
    factor(sample(letters[1:4], n, replace = TRUE)))
  names(df_list) <- cats
  df  <- as.data.frame(df_list)
  tgt <- setNames(
    lapply(cats, function(nm) setNames(rep(0.25, 4), letters[1:4])),
    cats
  )
  result <- harvest(df, tgt, method = "auto",
                    convergence = list(absolute = 1e-3))
  alg <- attr(result, "algorithm")
  expect_identical(alg, "raking",
    info = paste("Expected raking (product saturates) but got", alg))
})
```

**Build & verify:**
```bash
R CMD INSTALL --preclean /home/dd/Gemini/leafblower && \
  Rscript -e "devtools::test('/home/dd/Gemini/leafblower', filter='harvest')"
```

---

## Task 4 — R9: harvest.R @return documentation update

**File:** `R/harvest.R` — Roxygen `@return` block

**Current (line 94):**
```r
#'         \item \code{status}: 0=converged, 1=max_iter hit, 2=infeasible, 3=bad args.
```

**Bug:** Codes 4 (BUDGET) and 5 (STALL) are real return values from the C solver (defined in `leafblower.h` lines 36–37) that are surfaced to R callers via the result list. Code 1 description says "max_iter hit" but `RK_ERR_NOCONV=1` means "failed to converge" (not necessarily max_iter — the solver may exit early). BUDGET (4) is already exercised in existing tests (`test-harvest.R`: "budget exhausted" warning). STALL (5) is a valid non-error terminal state.

**Exact fix — replace the @return status line:**
```r
#'         \item \code{status}: integer status code. 0=converged (RK_OK);
#'           1=did not converge (RK_ERR_NOCONV); 2=infeasible (RK_ERR_INFEAS);
#'           3=bad argument (RK_ERR_BADARG); 4=budget exhausted — loss still
#'           decreasing, increase \code{max_iterations} (RK_ERR_BUDGET);
#'           5=loss plateau at constrained optimum, weights are valid (RK_ERR_STALL).
```

**No TDD for this task** — documentation only.

**Verify:**
```bash
Rscript -e "devtools::document('/home/dd/Gemini/leafblower')"
grep -A 10 'status.*integer' /home/dd/Gemini/leafblower/man/harvest.Rd
```

Expected: `man/harvest.Rd` contains "budget exhausted" and "RK_ERR_STALL" in the `\item{status}` entry.

---

## Execution Order

Tasks 1, 2, 3 are independent and can be implemented in parallel (they touch different files). Task 4 is documentation-only and can run last. B6 is closed — no work required.

Recommended sequence for a single-agent pass:

1. Write RED tests for T1, T2, T3 → confirm they all fail → `R CMD INSTALL --preclean .` → `devtools::test()`.
2. Apply fix for T1 (validation.hpp line 84). Rebuild. T1 RED→GREEN; others still RED.
3. Apply fix for T2 (validation.hpp new loop). Rebuild. T2 RED→GREEN.
4. Apply fix for T3 (cell_table.cpp K>8 block). Rebuild. T3 RED→GREEN.
5. Apply T4 doc fix. `devtools::document()`. Check harvest.Rd.
6. Run full test suite: `devtools::test()`. Confirm no regressions.

---

## Self-Review

- [x] Task 1: exact file/line identified (`validation.hpp` line 84), exact old/new code given, test exercises the 5e-7 gap zone and the still-rejected 2e-6 zone.
- [x] B6 CLOSED: `parse_convergence()` calls `match.arg()` for metric (line 483), rule (line 494), and stop_when (line 509) — the NA→UB path cannot be reached. No code change, no test. Confirmed from source.
- [x] Task 2: new code inserted after existing group_ids range loop, uses `RK_ERR_INFEAS` not `RK_ERR_BADARG`, O(n*K) only at validation time, partial-NA case tested.
- [x] Task 3: `int64_t` used throughout, early-saturation guard prevents overflow, two test scenarios (product < n and product ≥ n).
- [x] Task 4: no placeholder text, exact diff shown, verification command given.
- [x] No placeholder text anywhere.
- [x] All file paths are absolute-compatible (relative to package root as required by R conventions).
- [x] B5 (iEPPA/raking validation) explicitly excluded per spec.
- [x] calib_validate.cpp tolerance (1e-6) not changed — only validation.hpp loosened to match it.
