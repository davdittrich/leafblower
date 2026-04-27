# Convergence Reform Post-Review Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix all issues found in the post-merge critical code review of the a2p2 convergence reform (blocking bugs, required changes, suggestions) plus earlier follow-up tickets.

**Architecture:** Twelve atomic tasks ordered by severity. Tasks 1–3 fix blocking bugs. Task 4 completes the P1 validation consolidation. Tasks 5–8 fix required changes. Tasks 9–12 fix suggestions and docs. Each task is one beads ticket, one atomic commit. TDD throughout — failing test first for every behavioral change.

**Tech Stack:** C++17 (ieppa.cpp, lbfgsb_solver.cpp, r_bridge.cpp, c_api.cpp), R (harvest.R), Python (_harvest.py), testthat 3, pytest.

**Baseline:** `devtools::test()` → FAIL 0 | WARN 33 | SKIP 3 | PASS 285

**Status constants:** `RK_OK=0`, `RK_ERR_NOCONV=1`, `RK_ERR_INFEAS=2`, `RK_ERR_BADARG=3` (from `src/leafblower.h:32-34`)

---

## File Structure

| File | Tasks | Change |
|---|---|---|
| `src/ieppa.cpp` | 1, 5, 7 | W_best fix; stall warning; metrics gate |
| `src/c_api.cpp` + new `src/validation.hpp` | 4 | Extract validate_inputs to shared header |
| `src/r_bridge.cpp` | 2, 4 | Immediate enum guards; then replace with shared validator |
| `src/lbfgsb_solver.cpp` | 6 | Remove dead provisional capture |
| `R/harvest.R` | 3, 9, 10, 11 | Unknown-key validation; docs |
| `python/leafblower/_harvest.py` | 12 | attach_weights=False diagnostics |
| `tests/testthat/test-best-iterate.R` | 1 | W_best homotopy regression test |
| `tests/testthat/test-convergence-criteria.R` | 3, 5 | New tests |

---

## Task 1 — Fix W_best snapshot across homotopy levels (leafblower-z8wx)

**Files:** `src/ieppa.cpp:801`, `tests/testthat/test-best-iterate.R`

**Root cause:** Line 801 captures `W_best = W` where `W[c]` is the per-level capacity multiplier. Across homotopy levels `W[c]` resets at each level boundary; the level-0 contribution is folded into `X[c]`, not reflected in `W[c]`. The correct cumulative multiplier is `X[c] / X_init[c]` — the same quantity line 1036 uses for final weight expansion.

**Note on X_init[c] == 0:** Cells with `X_init[c] == 0` have zero observations (`n_per_cell[c] == 0`). No observation `i` maps to such a cell via `cell_of[i]`, so `W_best[c] = 0` for those cells never contributes to the expansion at line 1020. The `sum(best_weights) == n` guarantee is safe.

- [ ] **Step 1.1: Write failing regression test**

Append to `tests/testthat/test-best-iterate.R`:

```r
test_that("z8wx: best_weights sum=n and best_error<=max_error with homotopy_levels=3", {
  set.seed(77)
  n <- 800
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2","3","4"), n, replace = TRUE))
  )
  target <- list(
    a = c("1"=0.4, "2"=0.4, "3"=0.2),
    b = c("1"=0.25, "2"=0.25, "3"=0.25, "4"=0.25)
  )
  w <- leafblower::harvest(data, target, max_weight = 2, method = "ieppa",
                           max_iterations = 300,
                           convergence = list(absolute = 1e-6),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_equal(sum(result$best_weights), n, tolerance = 1e-6)
  expect_lte(result$best_error, result$max_error)
  expect_true(all(result$best_weights >= 0))
})
```

Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-best-iterate.R")' 2>&1 | tail -5`

With single-level default (no explicit homotopy setting), this test may currently pass. To force multi-level: check what `homotopy_levels` defaults to. If the default is 1, modify the call to add explicit homotopy settings when the API supports it, OR accept that this test will pass on the unfixed code and serve as a non-regression guard going forward (the bug manifests only when `homotopy_levels > 1` AND best iterate is on a non-final level).

- [ ] **Step 1.2: Fix W_best capture in `src/ieppa.cpp`**

Change lines 797–802 from:

```cpp
// WU-E: update best-iterate snapshot (tracks errRp regardless of active criterion).
if (errRp < best_errRp_seen) {
    best_errRp_seen = errRp;
    best_iter_val   = iter;
    W_best          = W;  // cell-level capacity multiplier snapshot
}
```

To:

```cpp
// WU-E: update best-iterate snapshot. Use X[c]/X_init[c] (cumulative
// multiplier from initial state) not W[c] (per-level factor only).
// W[c] resets at each homotopy level; X[c]/X_init[c] does not.
if (errRp < best_errRp_seen) {
    best_errRp_seen = errRp;
    best_iter_val   = iter;
    for (int c = 0; c < ct.M_cell; c++) {
        W_best[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
    }
}
```

- [ ] **Step 1.3: Build gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```
Expected: `* DONE (leafblower)`

- [ ] **Step 1.4: Run test**

```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-best-iterate.R")' 2>&1 | tail -5
```
Expected: all tests PASS.

- [ ] **Step 1.5: Full regression**

```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 286.

- [ ] **Step 1.6: Commit**

```bash
git add src/ieppa.cpp tests/testthat/test-best-iterate.R
git commit -m "$(cat <<'EOF'
fix(ieppa): snapshot X[c]/X_init[c] not W[c] for best_weights

W[c] is the per-level capacity multiplier; it resets at each homotopy
level boundary. X[c]/X_init[c] is the cumulative multiplier from the
initial state — the same quantity used for final weight expansion at
line 1036. Bug was dormant on single-level runs (W == X/X_init when
only one level). Regression test covers multi-level case.
EOF
)"
```

- [ ] **Step 1.7: Close ticket**

```bash
bd close leafblower-z8wx
```

---

## Task 2 — Immediate enum range guards in r_bridge.cpp (new ticket)

**Files:** `src/r_bridge.cpp:247-248`

Add the two missing guards now. The full validation consolidation (sharing validate_inputs across c_api and r_bridge) is Task 4. This task gets its own ticket since it is a distinct atomic change.

- [ ] **Step 2.0: Create ticket**

```bash
bd create --title "fix(r_bridge): immediate criterion/stop_when enum range guards" --type bug --priority 1 --description "r_bridge bypasses rk_calibrate() and must independently guard criterion in [0,4] and stop_when in [0,1] to prevent UB in solver dispatch switch. Partial fix for p7ry before full validation consolidation." 2>&1 | tail -2
# Note the created ticket ID — use it in Step 2.3
```

- [ ] **Step 2.1: Insert guards in `src/r_bridge.cpp`**

After line 246 (closing `}` of logit singularity block) and before line 248 (`// WU-E: call C++ solvers...`), insert:

```cpp
/* Criterion and stop_when must be in valid enum range.
   R callers are guarded by match.arg; direct C/Python callers are not. */
if (p.criterion < 0 || p.criterion > 4)
    Rf_error("leafblower: invalid arguments — criterion out of range [0,4]"
             " (0=PCT 1=MAX_ERR 2=MEAN_ERR 3=KL 4=CHI2)");
if (p.stop_when < 0 || p.stop_when > 1)
    Rf_error("leafblower: invalid arguments — stop_when out of range [0,1]"
             " (0=ANY 1=ALL)");
```

- [ ] **Step 2.2: Build gate + full regression**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: build clean, FAIL 0, PASS unchanged.

- [ ] **Step 2.3: Commit + close ticket**

```bash
git add src/r_bridge.cpp
git commit -m "$(cat <<'EOF'
fix(r_bridge): add criterion/stop_when enum range guards

r_bridge bypasses rk_calibrate() and must guard against out-of-range
enum casts that reach static_cast<lbw::CalibCriterion>(n). R callers
are protected by match.arg; direct C/Python callers are not.
Full validation consolidation follows in Task 4 (p7ry).
EOF
)"
# Close the ticket created in Step 2.0:
bd close <ticket-id-from-step-2.0>
```

---

## Task 3 — Unknown-key validation in parse_convergence and parse_sor (leafblower-azra)

**Files:** `R/harvest.R:286-318`, `tests/testthat/test-convergence-criteria.R`

**Note on NULL/non-list inputs:** Add a top-level type guard to `parse_convergence`. `names(NULL) = NULL`, `setdiff(NULL, valid) = character(0)` — safe. But `convergence = 0.001` passes the setdiff check and crashes later on `convergence[["pct"]]`. Guard explicitly.

- [ ] **Step 3.1: Write failing tests**

Append to `tests/testthat/test-convergence-criteria.R`:

```r
test_that("parse_convergence rejects unknown keys", {
  expect_error(
    leafblower::harvest(
      data.frame(a = factor(c("1","2"))),
      list(a = c("1"=0.5, "2"=0.5)),
      max_weight = 3, method = "ieppa",
      convergence = list(pct_tol = 0.001),
      attach_weights = FALSE
    ),
    regexp = "Unknown convergence key"
  )
})

test_that("parse_convergence rejects non-list convergence", {
  expect_error(
    leafblower::harvest(
      data.frame(a = factor(c("1","2"))),
      list(a = c("1"=0.5, "2"=0.5)),
      max_weight = 3, method = "ieppa",
      convergence = 1e-6,
      attach_weights = FALSE
    ),
    regexp = "convergence must be a named list"
  )
})

test_that("parse_sor rejects unknown keys", {
  expect_error(
    leafblower::harvest(
      data.frame(a = factor(c("1","2"))),
      list(a = c("1"=0.5, "2"=0.5)),
      max_weight = 3, method = "ieppa",
      sor = list(omega_minimum = 0.3),
      attach_weights = FALSE
    ),
    regexp = "Unknown sor key"
  )
})
```

Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")' 2>&1 | tail -5`
Expected: FAIL on all 3 new tests.

- [ ] **Step 3.2: Add validation to `R/harvest.R`**

Replace the opening of `parse_convergence` (line 286):

```r
parse_convergence <- function(convergence) {
  if (!is.null(convergence) && !is.list(convergence))
    stop("convergence must be a named list or NULL (e.g. list(pct = 0.001))")
  valid_keys <- c("pct", "absolute", "criterion", "stop_when")
  bad <- setdiff(names(convergence), valid_keys)
  if (length(bad))
    stop(sprintf("Unknown convergence key(s): %s. Valid keys: %s",
                 paste(bad, collapse = ", "),
                 paste(valid_keys, collapse = ", ")))
  `%||%` <- function(a, b) if (is.null(a)) b else a
  # ... rest of function unchanged ...
```

Add to `parse_sor` after the NULL early-return (before the list() return):

```r
  valid_keys <- c("auto", "omega_min", "omega", "omega_init", "burnin")
  bad <- setdiff(names(sor), valid_keys)
  if (length(bad))
    stop(sprintf("Unknown sor key(s): %s. Valid keys: %s",
                 paste(bad, collapse = ", "),
                 paste(valid_keys, collapse = ", ")))
```

- [ ] **Step 3.3: Build + run tests**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")' 2>&1 | tail -5
```
Expected: all 3 new tests PASS.

- [ ] **Step 3.4: Full regression**

```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 288.

- [ ] **Step 3.5: Commit**

```bash
git add R/harvest.R tests/testthat/test-convergence-criteria.R
git commit -m "$(cat <<'EOF'
fix(harvest.R): reject unknown + non-list keys in parse_convergence/parse_sor

list(pct_tol=0.001) silently used default; now raises error with
valid-key list. Non-list convergence arg raises immediately. Python
already raised ValueError — R now matches that behavior.
EOF
)"
```

- [ ] **Step 3.6: Close ticket**

```bash
bd close leafblower-azra
```

---

## Task 4 — Full validation consolidation: shared validate_inputs (leafblower-p7ry)

**Files:** new `src/validation.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`

**Goal:** Extract `validate_inputs()` from `c_api.cpp` into a shared inline header `src/validation.hpp`. Both `c_api.cpp` and `r_bridge.cpp` call the same function. Eliminates divergence risk.

**Confirmed: c_api.cpp already contains the enum guards** (lines 157-161):
```cpp
/* Enum range guards: prevent silent UB from out-of-range static_cast */
if (p->criterion < 0 || p->criterion > 4)
    return err("criterion out of range [0,4]: ...");
if (p->stop_when < 0 || p->stop_when > 1)
    return err("stop_when out of range [0,1]: ...");
```
The validation.hpp verbatim copy will include them. No new guards are added — this is extraction only.

The current `validate_inputs` signature (c_api.cpp:66-73):
```cpp
static int validate_inputs(int n, int K,
    const double* weights, const int32_t** group_ids,
    const int* cat_counts, const double** targets,
    const rk_params_t* p, rk_result_t* result, rk_algorithm_t alg);
```

- [ ] **Step 4.1: Create `src/validation.hpp`**

```cpp
#pragma once
#include "leafblower.h"
#include <cmath>
#include <cstddef>
#include <cstring>

namespace lbw {

// Validate inputs for rk_calibrate. Returns RK_OK or RK_ERR_BADARG.
// On RK_ERR_BADARG, sets result->status and result->message if result != nullptr.
inline int validate_calibrate_inputs(int n, int K,
    const double* weights, const int32_t** group_ids,
    const int* cat_counts, const double** targets,
    const rk_params_t* p, rk_result_t* result, rk_algorithm_t alg)
{
    auto err = [&](const char* msg) -> int {
        if (result) {
            result->status = RK_ERR_BADARG;
            snprintf(result->message, 256, "%s", msg);
        }
        return RK_ERR_BADARG;
    };
    // ── copy the full body of validate_inputs from c_api.cpp verbatim ──
    // (the enum guards are already at the end: criterion∈[0,4], stop_when∈[0,1])
    if (!weights)    return err("weights pointer is NULL");
    if (!group_ids)  return err("group_ids pointer is NULL");
    if (!cat_counts) return err("cat_counts pointer is NULL");
    if (!targets)    return err("targets pointer is NULL");
    if (n <= 0)      return err("n must be > 0");
    if (K <= 0)      return err("K must be > 0");
    if (K > 64)      return err("K exceeds maximum (64); too many margin columns");
    if (p->min_weight >= p->max_weight)
        return err("min_weight must be strictly less than max_weight");
    if (alg == RK_ALG_LBFGSB) {
        const double kSingularityEps = 1e-6;
        if (std::fabs(p->min_weight - 1.0) < kSingularityEps)
            return err("logit link undefined: min_weight near 1 makes denominator (1-L)~0");
        if (std::fabs(p->max_weight - 1.0) < kSingularityEps)
            return err("logit link undefined: max_weight near 1 makes denominator (U-1)~0");
    }
    size_t total_cats = 0;
    for (int k = 0; k < K; k++) {
        if (cat_counts[k] <= 0)  return err("cat_counts[k] must be > 0 for all k");
        if (cat_counts[k] > n)   return err("cat_counts[k] > n: more categories than observations");
        total_cats += (size_t)cat_counts[k];
    }
    if ((size_t)n * total_cats > SIZE_MAX / 2)
        return err("problem too large for platform size_t");
    double total_w = 0.0;
    for (int i = 0; i < n; i++) {
        if (!std::isfinite(weights[i])) return err("NaN or Inf in initial weights[]");
        total_w += weights[i];
    }
    if (total_w < 1e-15) return err("total weight is zero or negative");
    for (int k = 0; k < K; k++) {
        if (!targets[k]) return err("targets[k] is NULL");
        double sum = 0.0;
        for (int j = 0; j < cat_counts[k]; j++) {
            if (!std::isfinite(targets[k][j])) return err("NaN or Inf in targets[]");
            if (targets[k][j] < 0.0) return err("targets[k][j] < 0");
            sum += targets[k][j];
        }
        if (std::fabs(sum - 1.0) > 1e-8) return err("targets[k] does not sum to 1 (within 1e-8)");
    }
    for (int k = 0; k < K; k++) {
        if (!group_ids[k]) return err("group_ids[k] is NULL");
        for (int i = 0; i < n; i++) {
            int g = group_ids[k][i];
            if (g < -1)            return err("group_ids[k][i] < -1: only -1 (NA) is valid");
            if (g >= cat_counts[k]) return err("group_ids[k][i] >= cat_counts[k]");
        }
    }
    if (!std::isfinite(p->tol_abs) || p->tol_abs <= 0.0)
        return err("tol_abs must be finite and positive");
    if (p->criterion < 0 || p->criterion > 4)
        return err("criterion out of range [0,4]: 0=PCT 1=MAX_ERR 2=MEAN_ERR 3=KL 4=CHI2");
    if (p->stop_when < 0 || p->stop_when > 1)
        return err("stop_when out of range [0,1]: 0=ANY 1=ALL");
    return RK_OK;
}

} // namespace lbw
```

- [ ] **Step 4.2: Replace `validate_inputs` in `src/c_api.cpp`**

**Prerequisite:** Confirm `src/validation.hpp` from Step 4.1 contains all 13+ checks before replacing. Run: `grep -c "return err" /home/dd/Gemini/leafblower/src/validation.hpp` — expected: ≥ 13.

Replace the entire body of `static int validate_inputs(...)` (lines 66–164) with a one-liner that delegates to the shared header:

```cpp
#include "validation.hpp"

static int validate_inputs(int n, int K,
    const double* weights, const int32_t** group_ids,
    const int* cat_counts, const double** targets,
    const rk_params_t* p, rk_result_t* result, rk_algorithm_t alg) {
    return lbw::validate_calibrate_inputs(
        n, K, weights, group_ids, cat_counts, targets, p, result, alg);
}
```

- [ ] **Step 4.3: Replace partial validation in `src/r_bridge.cpp` with shared validator**

Find the validation block at lines ~229–246 (K>64, cat_counts, min/max, logit, criterion, stop_when guards). Replace with a `rk_result_t` call to the shared validator:

```cpp
// Full validation via shared validator (same as c_api.cpp path).
{
    rk_result_t validation_result;
    rk_result_init(&validation_result);
    rk_algorithm_t alg_for_validation =
        (strcmp(method_str, "lbfgsb") == 0) ? RK_ALG_LBFGSB :
        (strcmp(method_str, "raking") == 0) ? RK_ALG_RAKING : RK_ALG_IEPPA;
    int vrc = lbw::validate_calibrate_inputs(
        n, K, weights.data(),
        const_cast<const int32_t**>(group_ids.data()),
        cat_counts.data(),
        const_cast<const double**>(targets.data()),
        &p, &validation_result, alg_for_validation);
    if (vrc != RK_OK)
        Rf_error("leafblower: invalid arguments — %s", validation_result.message);
}
```

Remove the old per-check `Rf_error` calls that the shared validator now replaces.

- [ ] **Step 4.4: Build gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```
Expected: `* DONE (leafblower)` — no compile errors.

- [ ] **Step 4.5: Full regression**

```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 288.

- [ ] **Step 4.6: Commit**

```bash
git add src/validation.hpp src/c_api.cpp src/r_bridge.cpp
git commit -m "$(cat <<'EOF'
refactor(validation): extract shared validate_calibrate_inputs to validation.hpp

c_api.cpp and r_bridge.cpp now both call lbw::validate_calibrate_inputs,
eliminating the validation divergence risk. Future validation changes
apply to both call paths automatically.
EOF
)"
```

- [ ] **Step 4.7: Close ticket**

```bash
bd close leafblower-p7ry
```

---

## Task 5 — iEPPA stalled-infeasible PCT warning (leafblower-hawe)

**Files:** `src/ieppa.cpp` (exit block), `tests/testthat/test-convergence-criteria.R`

**Design decision:** Emit a log-channel warning only — do NOT change status to `RK_ERR_NOCONV`. Rationale: changing status would break backward-compatible callers who test `status == 0`. The warning tells users to check `max_error`. This is sufficient to address the ticket without introducing a backward-compat break. The status change option is left open for a future major version.

**Threshold derivation:** In a well-posed calibration problem, `errRp` and `pct_change` are coupled: as weights converge toward targets, both decrease together. The ratio `errRp / pct_change` stays in the range 1–5× at convergence on well-posed benchmarks (stepstone small, kk1204). When the problem is infeasible, weights stall (pct_change → 0) while errRp stays large — decoupling the ratio. A ratio of 10× (one order of magnitude) is the minimum gap that clearly separates stall-on-infeasible from normal convergence: `errRp > 10 * pct_tol`. This avoids false positives on problems where pct convergence fires slightly ahead of marginal-error convergence (ratio 1–5×) while reliably catching infeasible stalls (ratio 100×+). Exact threshold: `errRp > 10.0 * cfg.pct_tol`.

- [ ] **Step 5.1: Read ieppa.cpp exit block to find status assignment**

```bash
grep -n "RK_ERR_NOCONV\|RK_OK\|res\.status\|status = RK" /home/dd/Gemini/leafblower/src/ieppa.cpp | tail -20
```

Find the line where `res.status` is set to `RK_OK` after convergence (not the NOCONV/INFEAS assignments). Identify the line number.

- [ ] **Step 5.2: Write failing test**

Append to `tests/testthat/test-convergence-criteria.R`:

```r
test_that("hawe: iEPPA warns on PCT stall with large max_error", {
  # Contradictory: all of var1 is "A" maps to "2" in var2, but targets
  # demand opposite distributions — infeasible at max_weight=1.5
  n <- 400
  var1 <- factor(rep(c("A","B"), each = n/2))
  var2 <- factor(rep(c("2","1"), each = n/2))  # A↔2, B↔1
  data <- data.frame(var1 = var1, var2 = var2)
  target <- list(
    var1 = c(A = 0.95, B = 0.05),
    var2 = c("1" = 0.95, "2" = 0.05)
  )
  # With max_weight=1.5, solver cannot satisfy both margins — stall expected
  expect_warning(
    leafblower::harvest(data, target, max_weight = 1.5, method = "ieppa",
                        max_iterations = 300,
                        convergence = list(pct = 0.001),
                        attach_weights = FALSE),
    regexp = "PCT convergence stall"
  )
})
```

Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")' 2>&1 | tail -5`
Expected: FAIL (no stall warning emitted yet).

- [ ] **Step 5.3: Add stall detection to `src/ieppa.cpp`**

Find the block near the end of `ieppa_solve` where `res.status` gets its final value (after the main loop, before return). Insert after the status is confirmed:

```cpp
// PCT stall detection: if PCT convergence fired but max_error >> pct_tol,
// the problem is likely infeasible (weights stalled far from targets).
// Threshold: errRp > 10 * pct_tol (one order of magnitude gap between
// weight-change step size and marginal calibration error = stall signal).
// Emits warning only — status unchanged for backward compatibility.
{
    const auto& cfg = st.convergence_cfg;
    if (cfg.pct_tol > 0.0 &&
        res.max_error > 10.0 * cfg.pct_tol &&
        st.log_fn != nullptr) {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
            "PCT convergence stall: pct_change < %.3g but max_error=%.3g "
            "(%.0fx pct_tol). Possible contradictory or infeasible targets.",
            cfg.pct_tol, res.max_error,
            res.max_error / cfg.pct_tol);
        st.log_fn(msg, st.log_ctx);
    }
}
```

This warning fires only when `verbose > 0` because `st.log_fn` is only set when `p.verbose > 0` (r_bridge.cpp:187).

**harvest.R already calls warning() on non-convergence via the `msg` field in the result. Check if stall needs an explicit `warning()` in R-land:** Read `R/harvest.R` around line 200 (result handling). If the R bridge only warns on `status != RK_OK`, the stall won't surface unless we also emit an R-level warning. Add to `R/harvest.R` after unpacking the result:

Insert after line 213 in `R/harvest.R` (immediately after the existing `if (calib_result$status == 1L) warning(...)` block):

```r
# Stall detection: PCT converged (status=0) but max_error >> pct_tol
# signals infeasible problem. Threshold: 10x pct_tol (one order of magnitude gap).
if (calib_result$status == 0L &&
    !is.null(conv$pct_tol) && conv$pct_tol > 0 &&
    !is.null(calib_result$max_error) &&
    calib_result$max_error > 10 * conv$pct_tol) {
  warning(sprintf(
    "leafblower: PCT convergence stall: max_error=%.3g >> pct_tol=%.3g (%.0fx). ",
    calib_result$max_error, conv$pct_tol,
    calib_result$max_error / conv$pct_tol),
    "Possible contradictory or infeasible targets.",
    call. = FALSE)
}
```

- [ ] **Step 5.4: Build + run test**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")' 2>&1 | tail -5
```
Expected: stall test PASS. If `expect_warning(regexp="PCT convergence stall")` doesn't fire (infeasible fixture not stalling in 300 iters), increase `max_iterations` to 1000 or loosen the targets further.

- [ ] **Step 5.5: Full regression**

```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0. If a test expects NO warnings and the stall fires: that test's input is actually infeasible — pin it with `convergence = list(absolute = 1e-6)` or adjust `max_weight` to make it feasible.

- [ ] **Step 5.6: Commit**

```bash
git add src/ieppa.cpp R/harvest.R tests/testthat/test-convergence-criteria.R
git commit -m "$(cat <<'EOF'
fix(ieppa): warn on PCT convergence stall with large max_error

When pct_change < pct_tol but max_error > 10*pct_tol (one order of
magnitude gap: weight steps << calibration error), emit a warning.
Threshold 10x: well-posed problems have errRp/pct_change ratio 1-5x;
infeasible stalls show 100x+; 10x cleanly separates the two regimes.
Status unchanged for backward compat. R-level warning via harvest.R.
EOF
)"
```

- [ ] **Step 5.7: Close ticket**

```bash
bd close leafblower-hawe
```

---

## Task 6 — Verify best_iter is cumulative (leafblower-qbsf)

**Files:** `src/ieppa.cpp`

The ticket claims best_iter is per-level. Verify by reading the code and adding a test.

- [ ] **Step 6.1: Verify cumulative semantics**

```bash
grep -n "total_iters\|best_iter_val\|const int iter =" /home/dd/Gemini/leafblower/src/ieppa.cpp
```

Expected output includes:
- `int total_iters = 0;` (declared before the homotopy loop)
- `const int iter = total_iters + iter_in_lvl;` (cumulative: prior levels' iter count + current level's position)
- `total_iters = res.iterations;` (updated at end of each level)

This confirms `iter` IS cumulative across homotopy levels.

- [ ] **Step 6.2: Add confirming regression test**

Append to `tests/testthat/test-best-iterate.R`:

```r
test_that("qbsf: best_iter is positive (cumulative counter, not per-level reset)", {
  set.seed(88)
  n <- 500
  data <- data.frame(a = factor(sample(c("1","2","3"), n, replace = TRUE)))
  target <- list(a = c("1"=0.4, "2"=0.4, "3"=0.2))
  w <- leafblower::harvest(data, target, max_weight = 2, method = "ieppa",
                           max_iterations = 100,
                           convergence = list(absolute = 1e-10),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_gte(result$best_iter, 1L)
  expect_true(is.integer(result$best_iter) || is.numeric(result$best_iter))
})
```

Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-best-iterate.R")' 2>&1 | tail -5`
Expected: PASS (confirms existing implementation is correct).

- [ ] **Step 6.3: Build + full regression**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 1 higher.

- [ ] **Step 6.4: Commit + close ticket**

```bash
git add tests/testthat/test-best-iterate.R
git commit -m "test(ieppa): confirm best_iter is cumulative across homotopy levels"
bd close leafblower-qbsf --reason "iter = total_iters + iter_in_lvl where total_iters accumulates across levels. IS cumulative. Confirmed by test."
```

---

## Task 7 — Remove dead provisional best_weights capture in lbfgsb (leafblower-wmj2)

**Files:** `src/lbfgsb_solver.cpp:264-268`

- [ ] **Step 7.1: Remove dead code**

Lines 264–268 in `src/lbfgsb_solver.cpp` contain a provisional capture comment + 3 assignment lines that are immediately overwritten at line 683 post-normalization. Remove them:

```cpp
// DELETE these 5 lines (264-268):
// WU-E: L-BFGS-B is a batch solver with a single errRp evaluation at exit.
// best_error == max_error.  best_weights is a provisional capture here;
// lbfgsb_solve() re-assigns it after the Σ=n normalization to satisfy spec §4.
res.best_error   = max_err;
res.best_iter    = iterations;
res.best_weights = std::vector<double>(st.weights, st.weights + st.n);
```

Replace with (keep only the two scalar assignments; best_weights is set later):

```cpp
res.best_error = max_err;
res.best_iter  = iterations;
// best_weights assigned post-normalization in lbfgsb_solve() below.
```

- [ ] **Step 7.2: Build + full regression**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS unchanged.

- [ ] **Step 7.3: Commit + close**

```bash
git add src/lbfgsb_solver.cpp
git commit -m "cleanup(lbfgsb): remove dead provisional best_weights capture"
bd close leafblower-wmj2
```

---

## Task 8 — Gate extra metrics on active criterion (leafblower-a0gk)

**Files:** `src/ieppa.cpp`, `src/raking.cpp`

**Context:** mean_err, kl, chi2 computed at every kErrCheckInterval regardless of active criterion. For K=20, M_cell=1M this is O(20M) extra scatter-adds per check. The A7 test in `tests/testthat/test-quality-metrics.R` verifies metrics are populated at exit — this test is the regression guard for this task.

- [ ] **Step 8.1: Confirm A7 non-regression test behavior**

**Note:** Task 8 is a pure performance optimization with no user-visible behavior change. A test that fails BEFORE gating and passes AFTER cannot be written at the functional level — only benchmarks measure the perf gain. The correct TDD approach for pure performance changes is a non-regression guard: verify the behavior we're preserving (metrics present at exit) before AND after the change.

Verify `tests/testthat/test-quality-metrics.R` contains the A7 tests:

```bash
grep -n "A7" /home/dd/Gemini/leafblower/tests/testthat/test-quality-metrics.R
```

Expected: at least 3 test blocks present (iEPPA, raking, lbfgsb). These are the non-regression guards for Task 8.

Now add one more non-regression guard for the MAX_ERR criterion path (the gated path after implementation):

Append to `tests/testthat/test-quality-metrics.R`:

```r
test_that("a0gk: metrics non-zero at exit with MAX_ERR criterion (gated path)", {
  # After gating, mean_err/kl/chi2 are only computed on the final iter of each
  # level when criterion=MAX_ERR. This verifies exit-path metrics are populated.
  set.seed(42)
  n <- 500
  data <- data.frame(a = factor(sample(c("1","2"), n, replace = TRUE)))
  target <- list(a = c("1"=0.5, "2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                           max_iterations = 50,
                           convergence = list(absolute = 1e-6),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_true(is.finite(result$mean_error) && result$mean_error > 0)
  expect_true(is.finite(result$kl) && result$kl >= 0)
  expect_true(is.finite(result$chi2) && result$chi2 >= 0)
})
```

Run BEFORE implementation — expected PASS (metrics currently always computed). Run AGAIN after — still expected PASS (non-regression confirmed).

- [ ] **Step 8.2: Read ieppa.cpp metrics block**

```bash
grep -n "mean_err\|kl_max\|chi2_total\|need_extra" /home/dd/Gemini/leafblower/src/ieppa.cpp | head -20
```

Find the exact line range of the mean_err/kl/chi2 computation loop (part of the kErrCheckInterval block).

- [ ] **Step 8.3: Add `need_extra_metrics` gate in `src/ieppa.cpp`**

The kErrCheckInterval block in ieppa.cpp is triggered at line 763:
`if (iter == 1 || iter % kErrCheckInterval == 0 || iter_in_lvl == budget_lvl)`
where `budget_lvl` is the per-level iteration budget (line 327) and `iter_in_lvl` is the within-level counter (line 367). The `iter_in_lvl == budget_lvl` condition is the "final iteration of this level" gate.

Before the metrics accumulation loop (after `W_total` and `pct_change` computation), insert:

```cpp
// Gate: skip mean_err/kl/chi2 when they're not the active stopping criterion.
// Always compute on the final iteration of each level for exit-path reporting.
const bool need_extra_metrics =
    (cfg.criterion == lbw::CalibCriterion::MEAN_ERR ||
     cfg.criterion == lbw::CalibCriterion::KL       ||
     cfg.criterion == lbw::CalibCriterion::CHI2     ||
     iter_in_lvl == budget_lvl);  // final iter: always populate for calib_result
```

Wrap the mean_err/kl/chi2 loop with `if (need_extra_metrics) { ... }`.

Keep `res.mean_error = mean_err;` etc. unconditional after the gated block (variables init to 0.0 when not computed — that's correct for intermediate checks; the final-iter gate ensures non-zero at exit).

- [ ] **Step 8.4: Mirror in `src/raking.cpp`**

Same pattern. Find the metrics loop and add `need_extra_metrics` gate including `iter == st.inner_max_iter` (or whatever the final-iter sentinel is in raking.cpp) as an always-compute condition.

- [ ] **Step 8.5: Build + run A7 regression test**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
Rscript -e 'devtools::test_active_file("tests/testthat/test-quality-metrics.R")' 2>&1 | tail -5
```
Expected: all A7 tests PASS including new a0gk test.

If any A7 test fails with `mean_error == 0` at exit: the `iter_in_lvl == budget_lvl` final-iter gate is not firing. Verify `iter_in_lvl` is in scope at the insertion point (it is — declared at line 367, inside the same for-loop body). Also verify the gate was placed INSIDE the `if (iter == 1 || iter % kErrCheckInterval == 0 || iter_in_lvl == budget_lvl)` block, not outside it. Do NOT proceed to full regression until A7 is green.

- [ ] **Step 8.6: Full regression**

```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0.

- [ ] **Step 8.7: Commit + close**

```bash
git add src/ieppa.cpp src/raking.cpp tests/testthat/test-quality-metrics.R
git commit -m "$(cat <<'EOF'
perf(solvers): skip mean_err/kl/chi2 when not active stopping criterion

O(K*M_cell) scatter-adds skipped for unused metrics at intermediate
kErrCheckInterval. Gate: need_extra_metrics = (criterion is MEAN_ERR|KL|CHI2)
OR (final iteration of level). Exit-path reporting preserved: A7 tests pass.
EOF
)"
bd close leafblower-a0gk
```

---

## Task 9 — Docs: lbfgsb pct_change semantics (leafblower-q8pu)

**Files:** `R/harvest.R` (roxygen), then `Rscript -e 'devtools::document()'`

- [ ] **Step 9.1: Add lbfgsb note to `@param convergence` in `R/harvest.R`**

In the `@param convergence` block, after the `criterion` bullet, add:

```r
#'   \strong{Note for \code{method = "lbfgsb"}:} \code{pct_change} in the
#'   result measures the start-to-final weight shift (batch solver, single
#'   pass), not iteration-to-iteration shift as in iEPPA and raking. A
#'   \code{pct} threshold tuned for iEPPA will behave differently with lbfgsb.
```

- [ ] **Step 9.2: Regenerate Rd + regression**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS unchanged.

- [ ] **Step 9.3: Commit + close**

```bash
git add R/harvest.R man/harvest.Rd
git commit -m "docs(harvest.R): lbfgsb pct_change measures start->final not iter->iter"
bd close leafblower-q8pu
```

---

## Task 10 — Docs: chi2 W_total not comparable across solvers (leafblower-6pqx)

**Files:** `R/harvest.R` (roxygen)

Option (b) (normalize W_total=n in all solvers) is deferred — it changes metric values and requires a separate spec/plan. Option (a) (document) is the correct scope here.

- [ ] **Step 10.1: Add chi2 note to `@param convergence` criterion bullet**

In the `criterion` bullet documentation, add after the chi2 note about n-scaling:

```r
#'   \strong{chi2 cross-solver note:} chi2 is not directly comparable across
#'   methods. iEPPA uses unnormalized cell mass as \code{W_total}; raking and
#'   lbfgsb use \code{n}. Use chi2 as a convergence criterion within a single
#'   method, not to compare across methods.
```

- [ ] **Step 10.2: Regenerate Rd + regression**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'devtools::test()' 2>&1 | tail -3
```

- [ ] **Step 10.3: Commit + close**

```bash
git add R/harvest.R man/harvest.Rd
git commit -m "docs(harvest.R): chi2 W_total scale differs across solvers — document"
bd close leafblower-6pqx
```

---

## Task 11 — Docs: best_weights=zeros guard (leafblower-dguw)

**Files:** `R/harvest.R` (roxygen @return section)

- [ ] **Step 11.1: Add best_weights guard note to `@return`**

In the `@return` section of `harvest()`, add to the result list description:

```r
#'   \item{\code{best_weights}}{Numeric vector of length \code{n}, sum
#'     normalized to \code{n}. Weights at the iterate with minimum observed
#'     marginal error. When \code{best_error} is \code{Inf} (solver exited
#'     before first convergence check), \code{best_weights} is all-zero.
#'     Guard: \code{if (is.finite(attr(r, "result")$best_error)) ...} before use.}
```

- [ ] **Step 11.2: Regenerate Rd + regression**

```bash
Rscript -e 'devtools::document()'
Rscript -e 'devtools::test()' 2>&1 | tail -3
```

- [ ] **Step 11.3: Commit + close**

```bash
git add R/harvest.R man/harvest.Rd
git commit -m "docs(harvest.R): document best_weights=zeros when best_error=Inf"
bd close leafblower-dguw
```

---

## Task 12 — Python attach_weights=False drops result diagnostics (leafblower-j3v7)

**Files:** `python/leafblower/_harvest.py`, `python/leafblower/test_python.py`

- [ ] **Step 12.1: Read `python/leafblower/_harvest.py`**

Find the `attach_weights=False` return path. Understand the return type and where `result_dict` is currently dropped.

- [ ] **Step 12.2: Write failing test**

Append to `python/leafblower/test_python.py`:

```python
def test_result_accessible_when_attach_weights_false():
    """attach_weights=False must still expose result diagnostics."""
    data, target = _make_fixture(n=500)
    res = leafblower.harvest(data, target, max_weight=5, method="ieppa",
                             attach_weights=False)
    # Result diagnostics must survive when weights are not attached to df
    assert hasattr(res, "attrs") or hasattr(res, "result"), \
        "result diagnostics lost with attach_weights=False"
    r = res.attrs.get("result") if hasattr(res, "attrs") else getattr(res, "result", None)
    assert r is not None, "result dict is None with attach_weights=False"
    assert "max_error" in r
    assert "pct_change" in r
```

Run: `pytest python/leafblower/test_python.py::test_result_accessible_when_attach_weights_false -v`
Expected: FAIL.

- [ ] **Step 12.3: Implement WeightsResult wrapper in `python/leafblower/_harvest.py`**

After reading the file, locate the `attach_weights=False` path. Add a `WeightsResult` class or similar thin wrapper that carries `.attrs["result"]`:

```python
class WeightsResult:
    """Weights array that preserves result diagnostics (attach_weights=False)."""
    def __init__(self, weights_array, attrs):
        self._arr = weights_array
        self.attrs = attrs

    def __array__(self, dtype=None):
        return self._arr if dtype is None else self._arr.astype(dtype)

    def __len__(self):
        return len(self._arr)

    def __getitem__(self, key):
        return self._arr[key]

    def __repr__(self):
        return f"WeightsResult(n={len(self._arr)}, best_error={self.attrs.get('result', {}).get('best_error')})"
```

Return `WeightsResult(weights_array, {"result": result_dict})` when `attach_weights=False`.

- [ ] **Step 12.4: Run tests**

```bash
pip install -e . && pytest python/leafblower/test_python.py -v 2>&1 | tail -15
```
Expected: all tests PASS including new test.

- [ ] **Step 12.5: Commit + close**

```bash
git add python/leafblower/_harvest.py python/leafblower/test_python.py
git commit -m "$(cat <<'EOF'
fix(python): preserve result diagnostics when attach_weights=False

WeightsResult wrapper carries .attrs["result"] so callers retain
pct_change, best_error, sor, etc. regardless of attach_weights setting.
EOF
)"
bd close leafblower-j3v7
```

---

## Final Verification

- [ ] `devtools::test()` → FAIL 0, PASS ≥ 292
- [ ] `pytest python/leafblower/` → all green
- [ ] `R CMD check --as-cran` → 0 ERROR, 0 WARNING
- [ ] All 11 beads tickets closed (z8wx, azra, p7ry, hawe, qbsf, wmj2, a0gk, q8pu, 6pqx, dguw, j3v7)
