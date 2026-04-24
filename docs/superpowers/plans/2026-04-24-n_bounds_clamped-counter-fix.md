# Fix ieppa n_bounds_clamped Counter Undercount (rev 2)

**Bug:** `leafblower-kssd`. `n_bounds_clamped=0` reported for stepstone unit-mode despite ~36,745 obs pinned to `max_weight=5.0` (pre-water-fill `w_max=5.0079`, post `5.0000`).

**Root cause:** `src/ieppa.cpp:577-619` water-fill counter (`total_clamped`) only increments on two pathological fallback paths (lines 600/601, 615/616). Normal convergent clamps at lines 585 and 589 never bump the counter.

**Impact:** Diagnostic only. Correctness unaffected. Strict bounds still enforced.

## Rev history
- **rev 1 → rev 2** (plan-review-gate iter 1):
  - Scope: switched from post-loop scan to running-counter (2-line increment at 585/589). Semantically exact — post-scan would OVERcount FP-drift-near-bound obs that were never explicitly pinned.
  - Feasibility: verified running-counter has no double-count risk (once pinned, obs equals bound exactly; strict-inequality at 583 excludes; strict-inequality at 608 excludes from redistribute).
  - Completeness GAP-1: `test-ieppa-bounds-mode.R:60` asserts `expect_lt(n_bounds_clamped, 0.001*n)` — skewed-d fixture may now report real clamps under accurate counter. Relax to bounds-holding + nonzero-allowed.
  - Completeness GAP-2: R/harvest.R:141 and python/_harvest.py:181 warning text says "after water-fill exhausted" — update to match new semantic.
  - Completeness GAP-3: test-ieppa-bounds-mode.R:83 benign uniform-d test stays (no clamps fire on feasible benign data).
  - Completeness GAP-3b: new unit test synthetic data must exercise normal-path, not pathological `n_free==0`. Redesigned data.
  - Completeness GAP-4: added `bd close leafblower-kssd` to commit task.
  - Completeness GAP-5: stepstone claim relaxed to "> 0" not exact 36,745.

## Fix design

Running counter — increment `total_clamped` at each clamp site in the normal-path scan. Once pinned to `max_weight`/`min_weight` exactly, obs stays there (line 583 strict `>`; line 608 strict `<`). No double-count risk.

Pathological-path increments at 600-601, 615-616 are REDUNDANT (obs already counted at 585/589 in the same iter that detected the violation). Remove them.

## Scope

One file: `src/ieppa.cpp`. ~4-line change (2 adds, 4 removes).
Two docstring updates: `R/harvest.R:141`, `python/leafblower/_harvest.py:181`.
One test relax: `tests/testthat/test-ieppa-bounds-mode.R:60`.
One test add: unit test covering normal-path clamp counting.

## Pre-flight

- [ ] **Step P.1:** `git status --short` → baseline clean (plus pre-existing `.wolf/anatomy.md`).
- [ ] **Step P.2:** Baseline test count: 205 pass, 0 fail (post-rev-3 normalization commit `bac2877`).

## Task 1: Counter logic in src/ieppa.cpp

- [ ] **Step 1.1:** At line 585 add increment inside the clamp block:
```cpp
if (st.weights[i] > st.max_weight) {
    excess += st.weights[i] - st.max_weight;
    st.weights[i] = st.max_weight;
    any_violation = true;
    total_clamped++;   // NEW: count this clamp event
}
```
- [ ] **Step 1.2:** Same at line 589 for min bound:
```cpp
else if (st.weights[i] < st.min_weight) {
    excess -= st.min_weight - st.weights[i];
    st.weights[i] = st.min_weight;
    any_violation = true;
    total_clamped++;   // NEW: count this clamp event
}
```
- [ ] **Step 1.3:** Remove the pathological-path `total_clamped++` at lines 600, 601, 615, 616. The clamp actions (`st.weights[i] = st.max_weight` etc.) STAY — they enforce bounds — only the counter increments go. Rationale: those paths re-clamp obs that were ALREADY clamped (and counted) in the same iter's scan at lines 585/589.

- [ ] **Step 1.4:** Add a one-line comment above line 585 clamp block:
```cpp
// Running counter: increment at each normal-path clamp (strict-inequality
// above/below bound). Pinned weights stay at bound exactly, so next-iter
// scan (line 583 uses strict >, line 608 strict <) excludes them — no
// double count. Pathological re-clamps (n_free==0, budget exhausted) are
// redundant with these increments and do not re-count.
```

- [ ] **Step 1.5:** Build gate: `R CMD INSTALL --preclean .` clean.

## Task 2: Warning text updates

- [ ] **Step 2.1:** `R/harvest.R:140-142`:
```r
if (!is.null(calib_result$n_bounds_clamped) && calib_result$n_bounds_clamped > 0) {
  warning(sprintf(
    "unit-mode bounds: %d weights clamped to [%.3f, %.3f] during per-cell water-filling.",
    calib_result$n_bounds_clamped, min_weight, max_weight))
}
```

- [ ] **Step 2.2:** `python/leafblower/_harvest.py:179-181`:
```python
if result_dict.get("n_bounds_clamped", 0) > 0:
    warnings.warn(
        f"unit-mode bounds: {result_dict['n_bounds_clamped']} weights clamped "
        f"to [{params['min_weight']}, {params['max_weight']}] during per-cell water-filling.",
        ...)
```
(Preserve existing warning category and stacklevel args.)

## Task 3: Test fixture relax + new unit test

- [ ] **Step 3.1:** `tests/testthat/test-ieppa-bounds-mode.R:60` — existing assertion was written against broken counter (expected near-zero). Under accurate counter, skewed-d fixture may legitimately clamp. The real invariant is lines 57-58 (`max ≤ 3 + 1e-9`, `min ≥ 0.3 - 1e-9`). Relax line 60 to:
```r
# Accurate counter reports real clamps; strict bounds (lines 57-58) are the
# correctness check. Counter non-negative by construction.
expect_gte(info$n_bounds_clamped, 0L)
```

- [ ] **Step 3.2:** `tests/testthat/test-ieppa-bounds-mode.R:83` (benign uniform-d, spec §8) — KEEP unchanged. Benign data doesn't require clamping; running counter stays 0. Verified: uniform d + feasible targets → expansion weights naturally inside `[min, max]` → no clamp event at line 585/589.

- [ ] **Step 3.3:** Append new unit test that EXERCISES normal redistribution path (not pathological). Data design: 3 cells, minority cell has some obs that exceed cap but MAJORITY obs stay within bounds — water-fill redistributes cell's excess to the cell's free obs (normal convergent path):
```r
test_that("n_bounds_clamped counts normal-path water-fill clamps accurately", {
  # 3-cell synthetic: cell "a" has 20 obs at d=5 (will exceed cap), 80 obs
  # at d=0.8 (will stay free). Cells "b","c" match targets cleanly.
  # Within cell "a": 20 violators + 80 free → normal redistribute, not
  # pathological n_free==0.
  set.seed(17L)
  n <- 300L
  cat_a <- c(rep("a", 100), rep("b", 100), rep("c", 100))
  design <- c(rep(5.0, 20), rep(0.8, 80), rep(1.0, 200))
  df <- data.frame(a = factor(cat_a))
  tgt <- list(a = c(a = 1/3, b = 1/3, c = 1/3))
  res <- harvest(df, tgt, method = "ieppa",
                 max_weight = 1.5, min_weight = 0.3,
                 design_weights = design,
                 bounds_mode = "unit",
                 max_iterations = 500L,
                 convergence = list(absolute = 1e-6),
                 attach_weights = FALSE)
  info <- attr(res, "result")
  w <- as.numeric(res)
  expect_lte(max(w), 1.5 + 1e-9)
  expect_gte(min(w), 0.3 - 1e-9)
  # Under the old broken counter this returned 0 on normal-path clamps;
  # under the fix it must be > 0 because cell "a" has 20 cap violators.
  expect_gt(info$n_bounds_clamped, 0L)
})
```

- [ ] **Step 3.4:** Run full suite: `[FAIL 0 | PASS ≥ 206]` (205 baseline + 1 new; line 60 relax doesn't change count).

## Task 4: Stepstone diagnostic check (quick)

- [ ] **Step 4.1:** Re-run `benchmarks/stepstone_compare_current.R`. Assertions:
  - unit-mode `n_bounds_clamped > 0` (was 0 pre-fix; exact count not pinned — expected ≈36,745 but count depends on water-fill iteration semantics and may differ slightly)
  - unit-mode `max(w) <= 5.0 + 1e-9` still holds (correctness unchanged)
  - cell-mode values identical to pre-fix (diagnostic-only change doesn't touch cell-mode)

## Task 5: Commit

- [ ] **Step 5.1:**
```bash
git add src/ieppa.cpp R/harvest.R python/leafblower/_harvest.py \
        tests/testthat/test-ieppa-bounds-mode.R \
        docs/superpowers/plans/2026-04-24-n_bounds_clamped-counter-fix.md
git commit -m "$(cat <<'EOF'
fix(ieppa): n_bounds_clamped counter reflects every clamp, not only pathological paths

Prior: counter only incremented in two fallback paths (n_free==0
no-redistribute; 50-iter budget exhausted). Normal convergent water-fill
(the common case) clamped weights at st.max_weight / st.min_weight at
lines 585/589 without bumping the counter. On stepstone this produced
n_bounds_clamped=0 despite thousands of obs actually pinned to
max_weight=5.0.

Fix: inline increment at the two normal-path clamp sites. No double-count
risk: once pinned, obs equals the bound exactly; the per-iter scan uses
strict inequality (> max / < min) so pinned obs never re-triggers the
violation branch, and the redistribute check (lines 608) uses strict
inequality so pinned obs is never multiplied by factor. The pathological-
path increments are redundant (they re-clamp obs already counted in the
same iter's scan) and are removed.

Updated warning text in R and Python wrappers to match the broadened
semantic ("clamped to [lo, hi] during per-cell water-filling" vs the
prior "after water-fill exhausted" which was only accurate for the
pathological case).

Relaxed tests/testthat/test-ieppa-bounds-mode.R:60 which asserted near-
zero clamps under the broken counter; the real strict-bounds check at
lines 57-58 is the correctness invariant. Added a new test covering the
normal-redistribution path explicitly.

Stepstone unit-mode now reports non-zero n_bounds_clamped; weight values
unchanged (diagnostic-only fix).

Closes leafblower-kssd.
EOF
)"
```

- [ ] **Step 5.2:** `bd close leafblower-kssd --reason="Counter now reflects normal-path clamps; see commit"`

## Self-review

1. Running-counter is 2 increments + 4 removals + 2 docstring edits + 1 test relax + 1 test add. Net small.
2. No behavior change in weights. Diagnostic-only.
3. Double-count invariant proven by strict-inequality chain at 583, 591, 608.
4. Existing tests: line 83 (benign) unchanged; line 60 (skewed-d) relaxed because it asserted broken behavior.
5. New test exercises normal path via carefully constructed `design_weights` (20 violators + 80 free within cell "a").
