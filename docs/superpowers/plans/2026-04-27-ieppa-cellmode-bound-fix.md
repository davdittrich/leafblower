# iEPPA Cell-Mode Post-Normalization Bound Fix

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix iEPPA cell-mode weight output so max(weights) ≤ max_weight always (currently violated by up to 0.3% after normalization).

**Architecture:** Single 3-line change in `src/ieppa.cpp` line ~1193. The post-normalization cell-mode diagnostic block counts violations but never clamps. Add a clamp before counting.

**Tech Stack:** C++17, `src/ieppa.cpp` only. Build: `R CMD INSTALL --preclean . 2>&1 | tail -3`.

---

## Root Cause

```
P1.1 clamp:   X[c] = clamp(X_tilde_c, L_cell[c], U_cell[c])  ← within bounds ✓
Expansion:    w[i] = initial_w[i] * X[c] / X_init[c]          ← within bounds ✓
Normalization: w[i] *= n / total_w  (norm > 1 when total_w < n) ← can exceed max ✗
Cell mode:    counts violations, NO clamp                        ← violation persists ✗
```

Unit mode runs water-fill after normalization which enforces bounds. Cell mode does not.

---

## File Map

| File | Role |
|------|------|
| `src/ieppa.cpp:1193-1201` | Cell-mode diagnostic block — add clamp |
| `tests/testthat/test-calibration-solvers.R` | New regression test |

---

## Task 1: TDD RED — write failing bound test (leafblower-a3mr)

**Files:** Modify `tests/testthat/test-calibration-solvers.R` (append)

- [ ] **Step 1.1: Append test**

```r
test_that("ieppa: cell-mode weights respect max_weight hard cap", {
  # Reliable trigger: K=1, max_weight=1.5, targets=(0.9,0.1), uniform data.
  # Ideal weight for cat "1" ≈ 1.8 > cap 1.5 → cells clamped → W_total < n
  # → norm = n/W_total > 1 → post-norm wmax = 1.5 * norm > 1.5.  Bug fires.
  set.seed(1); n <- 1000L
  df <- data.frame(v1 = factor(sample(c("A","B"), n, replace=TRUE, prob=c(.5,.5))))
  tgt <- list(v1 = c("A"=0.9, "B"=0.1))

  r <- leafblower::harvest(df, tgt, method = "ieppa",
    max_weight = 1.5, min_weight = 0.0,
    bounds_mode = "cell", max_iterations = 500L,
    attach_weights = FALSE, verbose = 0)
  w <- as.numeric(r)

  # Before fix: post-norm wmax ≈ 1.5 * (1000/850) ≈ 1.76 > 1.5 → FAIL
  # After fix:  cell-mode clamp applied → wmax ≤ 1.5 → PASS
  expect_true(max(w) <= 1.5 + 1e-9,
    label = sprintf("max weight %.6f exceeds cap 1.5", max(w)))
  expect_true(min(w) >= 0.0 - 1e-9,
    label = sprintf("min weight %.6f below floor 0.0", min(w)))
})
```

- [ ] **Step 1.2: Verify RED**

```bash
cd /home/dd/Gemini/leafblower
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "FAIL|PASS|bound|cap" | head -5
```

Expected: FAIL (max weight exceeds 3.0).

If PASS (bounds already respected): this benchmark uses a different data distribution than stepstone. Try max_weight=2.5, n=20000 with highly skewed targets to force normalization to exceed bounds. Report NEEDS_CONTEXT.

- [ ] **Step 1.3: Commit**

```bash
git add tests/testthat/test-calibration-solvers.R
git commit -m "test(ieppa): RED — cell-mode weights must respect max_weight cap"
```

---

## Task 2: Fix — clamp in cell-mode diagnostic block (leafblower-a3mr)

**Files:** Modify `src/ieppa.cpp:1193-1201`

- [ ] **Step 2.1: Read lines 1193-1205**

Confirm current cell-mode block:
```cpp
    if (st.bounds_mode == RK_BOUNDS_CELL) {
        // Diagnostic scan: count violations without action.
        int violations = 0;
        for (int i = 0; i < st.n; i++) {
            if (st.weights[i] > st.max_weight || st.weights[i] < st.min_weight) {
                violations++;
            }
        }
        res.n_bounds_violated = violations;
    } else {
```

- [ ] **Step 2.2: Replace with clamp + count**

Replace the cell-mode block with:
```cpp
    if (st.bounds_mode == RK_BOUNDS_CELL) {
        // Clamp post-normalization violations: normalization (w *= n/total_w)
        // can push weights above max_weight when cells were clamped (W_total<n).
        // Cell mode previously only counted — now clamps too, mirroring unit mode.
        // n_bounds_clamped = n_bounds_violated: every violation is clamped.
        int violations = 0;
        for (int i = 0; i < st.n; i++) {
            if (st.weights[i] > st.max_weight) {
                st.weights[i] = st.max_weight;
                violations++;
            } else if (st.weights[i] < st.min_weight) {
                st.weights[i] = st.min_weight;
                violations++;
            }
        }
        res.n_bounds_violated = violations;
        res.n_bounds_clamped  = violations;   // every violation was clamped
    } else {
```

- [ ] **Step 2.3: Compile**

```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
```

- [ ] **Step 2.4: Verify GREEN**

```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "FAIL|PASS|bound|cap" | head -5
```

Expected: PASS.

- [ ] **Step 2.5: Full test suite**

```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```

Expected: FAIL 0 (or same pre-existing 1), PASS ≥ 332.

- [ ] **Step 2.6: Spot-check benchmark wmax**

```bash
OMP_NUM_THREADS=1 Rscript -e '
suppressPackageStartupMessages({library(arrow); library(leafblower)})
df  <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
df$uuid <- NULL
library(jsonlite)
tgt <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
tgt <- lapply(tgt, function(t) { t <- unlist(t); t/sum(t) })
for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
r <- suppressWarnings(leafblower::harvest(df, tgt, method="ieppa",
  max_weight=5, min_weight=0, max_iterations=500L, attach_weights=FALSE))
w <- as.numeric(r)
cat(sprintf("wmax=%.6f wmin=%.6f (cap: 5.0 / 0.0)\n", max(w), min(w)))
' 2>&1 | grep wmax
```

Expected: `wmax <= 5.000000`.

- [ ] **Step 2.7: Commit**

```bash
git add src/ieppa.cpp
git commit -m "fix(ieppa): clamp cell-mode weights after normalization

Post-expansion normalization (w *= n/total_w) can push individual weights
above max_weight when total calibrated weight < n. Cell mode previously
only counted violations without clamping. Add clamp matching unit-mode
semantics. n_bounds_violated now counts post-clamp residuals (should be 0).

Closes: leafblower-a3mr"
```

---

## Self-Review

**Spec coverage:**

| Item | Task |
|------|------|
| ROOT: normalization inflates cell-mode weights | documented in root cause |
| FIX: clamp in cell-mode block (ieppa.cpp:1193) | Task 2 Step 2.2 |
| RED test before fix | Task 1 |
| GREEN after fix | Task 2 Step 2.4 |
| Benchmark wmax ≤ 5.0 confirmed | Task 2 Step 2.6 |

**Placeholder scan:** None.

**Type consistency:** `st.max_weight`, `st.min_weight` (double), `st.weights[i]` (double). Consistent with surrounding code.
