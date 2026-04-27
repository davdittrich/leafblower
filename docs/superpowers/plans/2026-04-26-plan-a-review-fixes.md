# Plan A Review Fixes

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix six issues from the /critical-code-reviewer pass on Plan A: two c_api.cpp field aliasing bugs, one Inf init bug, two test quality gaps, and one r_bridge documentation fix.

**Architecture:** All fixes are local and independent. No new test files. c_api.cpp changes are verified by inspection and build gate (no R-level test can reach the C ABI init path — r_bridge.cpp bypasses rk_result_init defaults entirely). Two commits: one for C++ (c_api.cpp), one for tests+docs (r_bridge.cpp, test-calibration-solvers.R).

**Honesty note on testability:** `r_bridge.cpp` calls the solver directly and reads `res.convergence_objective` via `pack_solver_result` — it never touches `rk_result_t` fields set by `rk_result_init`. The c_api.cpp fixes affect the C ABI consumer path only. R-level tests cannot distinguish the two paths. The fixes are verified by inspection (grep confirms old patterns gone) and by the existing test suite remaining green.

**Tech Stack:** C++17 (src/c_api.cpp, src/r_bridge.cpp), R (tests/testthat/test-calibration-solvers.R), testthat 3.

**Baseline:** FAIL 0 | PASS 346 | SKIP 4

---

## File Structure

| File | Action | Tickets |
|---|---|---|
| `src/c_api.cpp` | Modify | leafblower-y8g7, leafblower-9uhn |
| `src/r_bridge.cpp` | Modify | leafblower-1xa6, leafblower-aiil |
| `tests/testthat/test-calibration-solvers.R` | Modify | leafblower-56d4, leafblower-4dtl |

---

## Task 1 — c_api.cpp: fix field aliasing in all three solver branches (leafblower-y8g7)

**File:** `src/c_api.cpp`

No new test — the fix is verified by inspection and the existing test suite. The r_bridge path already correctly uses `res.convergence_objective` (line 345 of r_bridge.cpp); this fix makes the C ABI path consistent.

- [ ] **Step 1.1: Verify current code pattern**

```bash
grep -n "convergence_objective\s*=\s*res\.best_error\|convergence_minimized_metric\s*=\s*res\.convergence_metric" src/c_api.cpp
```
Expected: 6 lines — 3 for each pattern (lbfgsb, raking, ieppa branches).

- [ ] **Step 1.2: Fix field aliasing — replace all 6 occurrences**

In each of the three solver branches (~lines 208-209, 231-232, 263-264), replace:
```cpp
// OLD:
result->convergence_objective           = res.best_error;
result->convergence_minimized_metric    = res.convergence_metric;

// NEW:
result->convergence_objective           = res.convergence_objective;
result->convergence_minimized_metric    = res.convergence_minimized_metric;
```

Apply to ALL three branches. Do not change any other lines.

- [ ] **Step 1.3: Verify fix applied and r_bridge.cpp path is undisturbed**

```bash
# Must produce no output (old pattern gone):
grep -n "convergence_objective\s*=\s*res\.best_error" src/c_api.cpp
grep -n "convergence_minimized_metric\s*=\s*res\.convergence_metric" src/c_api.cpp

# r_bridge.cpp must still use res.convergence_objective in pack_solver_result (unchanged):
grep -n "res\.convergence_objective\|res\.convergence_minimized_metric" src/r_bridge.cpp
```
Expected: r_bridge.cpp lines 345-346 unchanged.

- [ ] **Step 1.4: Build gate**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```
Expected: `* DONE (leafblower)`

---

## Task 2 — c_api.cpp: initialize convergence_objective to Inf (leafblower-9uhn)

**File:** `src/c_api.cpp`

No new test — this fixes early-exit C ABI behavior not reachable from R tests.

- [ ] **Step 2.1: Locate initialization**

```bash
grep -n "convergence_objective\|best_error.*infinity\|numeric_limits" src/c_api.cpp | head -10
```
Expected: line ~64 sets `r->best_error = std::numeric_limits<double>::infinity()`. Line ~69 sets `r->convergence_objective = 0.0;` (from memset default).

- [ ] **Step 2.2: Fix initialization**

In `rk_result_init`, immediately after the `r->best_error = std::numeric_limits<double>::infinity();` line, add:
```cpp
r->convergence_objective = std::numeric_limits<double>::infinity();
```

`<limits>` is already included in c_api.cpp — no new include needed.

- [ ] **Step 2.3: Build gate + regression**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: `* DONE (leafblower)`, FAIL 0, PASS ≥ 346.

- [ ] **Step 2.4: Commit Tasks 1+2**
```bash
git add src/c_api.cpp
git commit -m "$(cat <<'EOF'
fix(c_api): use dedicated convergence fields + Inf init

result->convergence_objective now reads res.convergence_objective (not
res.best_error); result->convergence_minimized_metric reads
res.convergence_minimized_metric (not res.convergence_metric). Semantically
equivalent today but semantically correct for future divergence.
rk_result_init: convergence_objective initialized to Inf (consistent with
best_error) so C ABI early-exit callers get Inf sentinel, not 0.0.
EOF
)"
```

---

## Task 3 — T1a: add regexp to expect_error (leafblower-56d4)

**File:** `tests/testthat/test-calibration-solvers.R`

- [ ] **Step 3.1: Find the T1a test**

```bash
grep -n "T1a\|expect_error" tests/testthat/test-calibration-solvers.R | head -10
```

- [ ] **Step 3.2: Add regexp**

Inside the `for (m in c("sinkhorn", ...))` loop, change:
```r
# OLD:
expect_error(
  leafblower::harvest(data, target, max_weight=3, method=m, attach_weights=FALSE),
  info = paste("method", m, "should error not crash")
)

# NEW:
expect_error(
  leafblower::harvest(data, target, max_weight=3, method=m, attach_weights=FALSE),
  regexp = "not yet implemented",
  info = paste("method", m, "should error with 'not yet implemented'")
)
```

- [ ] **Step 3.3: Confirm T1a still passes**

```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "T1a|FAIL|PASS" | head -5
```
Expected: T1a passes (r_bridge.cpp's stub message contains "not yet implemented").

---

## Task 4 — A1 test: remove library() calls (leafblower-4dtl)

**File:** `tests/testthat/test-calibration-solvers.R`

- [ ] **Step 4.1: Find library() calls in A1 test**

```bash
grep -n "library(arrow)\|library(jsonlite)" tests/testthat/test-calibration-solvers.R
```

- [ ] **Step 4.2: Replace with skip guards and remove redundant calls**

The A1 test already uses `arrow::read_parquet()` and `jsonlite::fromJSON()` with package namespace. The bare `library()` calls are redundant AND unsafe. Replace:

```r
# OLD:
library(arrow); library(jsonlite)
data <- arrow::read_parquet(...)
tgt_raw <- jsonlite::fromJSON(...)

# NEW (remove library() lines entirely — namespace calls suffice):
skip_if_not_installed("arrow")
skip_if_not_installed("jsonlite")
data <- arrow::read_parquet(...)
tgt_raw <- jsonlite::fromJSON(...)
```

---

## Task 5 — r_bridge.cpp: annotate stub check + fix stale comment (leafblower-1xa6, leafblower-aiil)

**File:** `src/r_bridge.cpp`

- [ ] **Step 5.1: Annotate stub check for Plan C**

```bash
grep -n "Stub methods\|sinkhorn\|chebyshev\|grake\|greg" src/r_bridge.cpp | head -10
```

Find the early-exit block for stub methods. Add comment immediately before the `if (strcmp(method_str, "sinkhorn") == 0 ...)` block:
```cpp
/* Stub methods: error immediately rather than falling through to "unknown method".
 * Plan C: remove this block and wire p.algorithm = RK_ALG_SINKHORN/CHEBYSHEV/GREG/GRAKE
 * in the if/else chain below. */
```

- [ ] **Step 5.2: Fix stale count comment**

Find the `Rf_allocVector(VECSXP, 30)` line. Update comment:
```cpp
// OLD: 14 prior + 8 scalars + best_weights + 5 convergence fields + 2 new
// NEW: 14 prior + 8 scalars + best_weights + 7 convergence fields
```

- [ ] **Step 5.3: Build gate**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```
Expected: `* DONE (leafblower)`

- [ ] **Step 5.4: Commit Tasks 3+4+5**
```bash
git add tests/testthat/test-calibration-solvers.R src/r_bridge.cpp
git commit -m "$(cat <<'EOF'
fix(tests+docs): T1a regexp, A1 skip guards, r_bridge annotations

T1a: expect_error gains regexp='not yet implemented' — passes on any error
was too weak. A1: library() removed (redundant with :: calls); replaced with
skip_if_not_installed() guards. r_bridge: stub check annotated for Plan C
removal; stale comment updated from '5+2' to '7 convergence fields'.
EOF
)"
```

---

## Final Verification

- [ ] `grep "convergence_objective\s*=\s*res\.best_error" src/c_api.cpp` → no output
- [ ] `grep "convergence_minimized_metric\s*=\s*res\.convergence_metric" src/c_api.cpp` → no output
- [ ] `grep "res\.convergence_objective" src/r_bridge.cpp` → line 345 present (r_bridge path undisturbed)
- [ ] `grep "convergence_objective.*infinity\|infinity.*convergence_objective" src/c_api.cpp` → init line present
- [ ] `Rscript -e 'devtools::test()' 2>&1 | tail -3` → FAIL 0, PASS ≥ 346
