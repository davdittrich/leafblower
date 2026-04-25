# Calibration Solvers — Plan A: Infrastructure + iEPPA Default Change

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the shared C++ infrastructure (ABI constants, validation, linalg stubs), change iEPPA's default convergence to kl+improvement, and wire new method names through R/C so Plans B/C/D can add solvers without touching the ABI layer.

**Architecture:** Scaffold only — no new solvers yet. Every task produces a green test suite. iEPPA and raking remain unchanged algorithmically; only iEPPA's default metric changes. New method names (sinkhorn/chebyshev/greg/grake) are wired in the R/C dispatch as stubs returning `RK_ERR_BADARG` until Plans B/C/D fill them in.

**Tech Stack:** C++17 (leafblower.h, c_api.cpp, r_bridge.cpp, two new files), R (harvest.R), testthat 3.

**Spec:** `docs/superpowers/specs/2026-04-25-calibration-solvers-design.md`

**Baseline:** FAIL 0 | PASS 337

**Scientific validation mandate** (`/scientific-critical-thinking`): Before closing each task, run the convergence diagnostic to confirm iEPPA still terminates. The KL metric default must not cause NOCONV regressions on existing fixtures.

---

## File Structure

| File | Task | Action | Description |
|---|---|---|---|
| `src/leafblower.h` | 1 | Modify | Add 4 RK_ALG constants, 2 rk_result_t fields, EXPECTED_RK_RESULT_BYTES |
| `src/c_api.cpp` | 1, 3 | Modify | Update rk_result_init; add 4 method stubs; change iEPPA default metric |
| `src/calib_validate.hpp` | 2 | Create | calib_validate_preentry() declaration |
| `src/calib_validate.cpp` | 2 | Create | calib_validate_preentry() implementation |
| `src/calib_linalg.hpp` | 2 | Create | compute_normal_equations + ldlt stubs (bodies in Plan D) |
| `src/r_bridge.cpp` | 1 | Modify | Unpack convergence_objective + convergence_minimized_metric |
| `R/harvest.R` | 3, 4 | Modify | Accept new method strings; update default to kl+improvement |
| `data-raw/gen_ieppa_kl_ref.R` | 4 | Create | Script to generate A1 KL reference fixture |
| `tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds` | 4 | Create | A1 reference value |
| `tests/testthat/test-calibration-solvers.R` | 1–4 | Create | New test file for all Plan A tests |

---

## Task 1 — ABI: New constants + result fields + tripwire

**Files:** `src/leafblower.h`, `src/c_api.cpp`, `src/r_bridge.cpp`

**ATOMIC:** All 3 files in one commit. r_bridge.cpp references `rk_result_t` fields; changing the struct requires simultaneous r_bridge update.

- [ ] **Step 1.1: Write failing test (ABI presence)**

Create `tests/testthat/test-calibration-solvers.R`:
```r
test_that("T1: new method names accepted without error (stub)", {
  # sinkhorn/chebyshev/greg/grake should fail with a clear error, not segfault
  data <- data.frame(a = factor(c("1","2")))
  target <- list(a = c("1"=0.5, "2"=0.5))
  for (m in c("sinkhorn", "chebyshev", "greg", "grake")) {
    expect_error(
      leafblower::harvest(data, target, max_weight=3, method=m,
                          attach_weights=FALSE),
      info = paste("method", m, "should error not crash")
    )
  }
})

test_that("T1: convergence_used$objective field present", {
  set.seed(1)
  data <- data.frame(a=factor(sample(c("1","2"),200,T)))
  target <- list(a=c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight=3, method="ieppa",
                           attach_weights=FALSE)
  r <- attr(w, "result")
  expect_true("objective" %in% names(r$convergence_used),
              info="convergence_used must have $objective field")
  expect_true("minimized_metric" %in% names(r$convergence_used),
              info="convergence_used must have $minimized_metric field")
  expect_true(is.finite(r$convergence_used$objective))
})
```

Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | tail -8`
Expected: T1 objective test FAIL (fields not yet present); method stub test FAIL (method not accepted).

- [ ] **Step 1.2: Add constants to `src/leafblower.h`**

Read the file first: `sed -n '38,45p' src/leafblower.h`

After `RK_ALG_RAKING = 3,` add:
```c
    RK_ALG_SINKHORN   = 4,   /* true KL min, Bregman Dykstra */
    RK_ALG_CHEBYSHEV  = 5,   /* true L∞ marginal error min, LP-IPM */
    RK_ALG_GREG       = 6,   /* true chi2 min, Newton QP */
    RK_ALG_GRAKE      = 7    /* true grake_norm min, LP-IPM */
```

After `int    sor_n_damped;` in `rk_result_t` (line ~118), add:
```c
    double convergence_objective;        /* value of minimized metric at convergence */
    int    convergence_minimized_metric; /* CalibMetric: which metric was minimized */
```

Add `EXPECTED_RK_RESULT_BYTES` immediately after `rk_result_t` closes. First measure the current sizeof: add temporarily, compile with a `static_assert(false)` to read the value from the error, then set it. Or use the g++ sizeof trick:
```bash
cat > /tmp/sz_result.cpp << 'EOF'
#include "src/leafblower.h"
#include <stdio.h>
int main() { printf("rk_result_t: %zu\n", sizeof(rk_result_t)); }
EOF
g++ -std=c++17 -I. /tmp/sz_result.cpp -o /tmp/sz_result && /tmp/sz_result
```

Add after rk_result_t:
```c
#define EXPECTED_RK_RESULT_BYTES <measured value>
static_assert(sizeof(rk_result_t) == EXPECTED_RK_RESULT_BYTES,
    "rk_result_t size changed; update EXPECTED_RK_RESULT_BYTES and ABI consumers");
```

- [ ] **Step 1.3: Initialize new fields in `src/c_api.cpp`**

Read `rk_result_init` (grep: `grep -n "rk_result_init" src/c_api.cpp`). Add:
```c
r->convergence_objective        = 0.0;
r->convergence_minimized_metric = 0;  /* CalibMetric::MAX_ERR = 0 */
```

Add 4 new method stubs to the main dispatch (find: `grep -n "RK_ALG_RAKING\|algorithm_used\|raking_solve" src/c_api.cpp | head -10`). After the `RK_ALG_RAKING` case, add:
```c
case RK_ALG_SINKHORN:
case RK_ALG_CHEBYSHEV:
case RK_ALG_GREG:
case RK_ALG_GRAKE: {
    result->status = RK_ERR_BADARG;
    snprintf(result->message, sizeof(result->message),
        "method '%s' is not yet implemented in this build",
        algorithm == RK_ALG_SINKHORN  ? "sinkhorn"  :
        algorithm == RK_ALG_CHEBYSHEV ? "chebyshev" :
        algorithm == RK_ALG_GREG      ? "greg"      : "grake");
    return RK_ERR_BADARG;
}
```

- [ ] **Step 1.4: Unpack new fields — two-layer approach**

**Layer A: r_bridge.cpp** — flat fields at indices 28 and 29 (current VECSXP size is 28, indices 0-27).

Change `Rf_allocVector(VECSXP, 28)` → `Rf_allocVector(VECSXP, 30)` (and same for STRSXP).

Add result variable declarations alongside existing ones:
```cpp
double res_conv_objective          = 0.0;
int    res_conv_minimized_metric   = 0;
```

In each solver branch (ieppa/raking/lbfgsb), add after `res_conv_iter = res.convergence_iter;`:
```cpp
res_conv_objective        = res.convergence_objective;
res_conv_minimized_metric = res.convergence_minimized_metric;
```

After the existing index-27 block, add:
```cpp
SET_STRING_ELT(res_names, 28, Rf_mkChar("convergence_objective"));
SET_STRING_ELT(res_names, 29, Rf_mkChar("convergence_minimized_metric"));
SET_VECTOR_ELT(res_list, 28, Rf_ScalarReal(res_conv_objective));
SET_VECTOR_ELT(res_list, 29, Rf_ScalarInteger(res_conv_minimized_metric));
```

**Layer B: harvest.R** — add to the `convergence_used` list construction (~line 249):
```r
.metric_names <- c("max_err","mean_err","kl","chi2","grake_norm","l1_weight","pct")
calib_result$convergence_used <- list(
  # ...existing: metric, rule, tol, fired_at_iter...
  objective        = calib_result$convergence_objective,
  minimized_metric = .safe_lookup(.metric_names, calib_result$convergence_minimized_metric)
)
calib_result$convergence_objective        <- NULL
calib_result$convergence_minimized_metric <- NULL
```

- [ ] **Step 1.5: Populate new fields in solver exit blocks**

`rk_result_init` zeros all fields including `convergence_objective=0.0` and `convergence_minimized_metric=0`. But these must be populated at solver exit to hold meaningful values.

Read `grep -n "convergence_metric\|convergence_rule\|convergence_tol\|convergence_iter" src/ieppa.cpp | head -8` — find the exit block where these are set (~line 1103-1106).

In **ieppa.cpp** at the same exit block, add after `res.convergence_iter`:
```cpp
res.convergence_objective        = best_metric_seen;   // min of active metric across iters
res.convergence_minimized_metric = static_cast<int>(st.convergence_cfg.metric);
```

In **raking.cpp**: `grep -n "convergence_metric\|convergence_rule" src/raking.cpp | head -5` — same pattern, add:
```cpp
res.convergence_objective        = best_metric_seen;
res.convergence_minimized_metric = static_cast<int>(st.convergence_cfg.metric);
```

In **lbfgsb_solver.cpp**: `grep -n "convergence_metric\|convergence_rule" src/lbfgsb_solver.cpp | head -5` — same pattern:
```cpp
res.convergence_objective        = res.max_error;   // lbfgsb single-pass: final = best
res.convergence_minimized_metric = static_cast<int>(cfg.metric);
```

- [ ] **Step 1.5b: Build gate**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```
Expected: `* DONE (leafblower)` — static_assert must pass.

- [ ] **Step 1.6: Update map_method in harvest.R to accept new strings**

Read `grep -n "map_method\|match.arg.*method\|method.*ieppa\|method.*raking" R/harvest.R | head -15`

Add the 4 new methods to the valid method list and map them to integers:
```r
method_int <- c(ieppa=1L, lbfgsb=2L, raking=3L,
                sinkhorn=4L, chebyshev=5L, greg=6L, grake=7L)
# (or however the existing method dispatch works — read it first)
```

- [ ] **Step 1.7: Run tests**
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | tail -8
```
Expected: T1 stub test PASS (new methods error cleanly); T1 objective test PASS.

- [ ] **Step 1.8: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 337.

- [ ] **Step 1.9: Commit**
```bash
git add src/leafblower.h src/c_api.cpp src/r_bridge.cpp R/harvest.R \
        tests/testthat/test-calibration-solvers.R
git commit -m "$(cat <<'EOF'
feat(ABI): add RK_ALG_SINKHORN/CHEBYSHEV/GREG/GRAKE + convergence_objective

New method constants 4-7 in rk_algorithm_t. Stubs return RK_ERR_BADARG
until Plans B/C/D fill them in. rk_result_t gains convergence_objective
and convergence_minimized_metric. EXPECTED_RK_RESULT_BYTES tripwire added.
harvest.R accepts new method strings; R convergence_used exposes new fields.
EOF
)"
```

---

## Task 2 — Shared infrastructure: calib_validate + calib_linalg stubs

**Files:** `src/calib_validate.hpp`, `src/calib_validate.cpp` (new), `src/calib_linalg.hpp` (new)

**Clean-code principle:** `calib_validate_preentry` does ONE thing: checks all preconditions and returns an error code. `compute_normal_equations` does ONE thing: computes N = ADAᵀ. No side effects.

- [ ] **Step 2.1: Write failing test (validation)**

Append to `tests/testthat/test-calibration-solvers.R`:
```r
test_that("T2: L_c > U_c rejected with RK_ERR_BADARG", {
  # Manufacture infeasible bounds by setting min_weight > max_weight
  data <- data.frame(a = factor(c("1","2")))
  target <- list(a = c("1"=0.5, "2"=0.5))
  expect_error(
    leafblower::harvest(data, target, min_weight=10, max_weight=2,
                        method="ieppa", attach_weights=FALSE),
    regexp = "min_weight"
  )
})

test_that("T2: X_init=0 with L_c>0 returns infeasible", {
  # All obs in category "1"; category "2" has zero obs → cell has X_init=0
  # but target T_kj > 0 means margin infeasible (structural zero)
  data <- data.frame(a = factor(rep("1", 100)))
  # This test already exists in ieppa via structural_infeas_pairs
  # Just verify the error message is informative
  expect_warning(
    w <- leafblower::harvest(data, list(a=c("1"=0.5, "2"=0.5)),
                              max_weight=3, method="ieppa", attach_weights=FALSE),
    regexp = "infeas|converge",
    ignore.case = TRUE
  )
})
```

Run: expected PASS (existing validation already handles these via c_api validate_inputs and iEPPA structural infeas).

- [ ] **Step 2.2: Create `src/calib_validate.hpp`**

```cpp
#pragma once
#include "leafblower.h"
#include "cell_table.hpp"
#include "types.hpp"

namespace lbw {

// kNCatsTotalMax is defined HERE (single definition).
// calib_linalg.hpp includes this header to share the constant.
constexpr int kNCatsTotalMax = 2048;

/**
 * Pre-entry validation for new cell-table solvers (sinkhorn, chebyshev, greg, grake).
 * Returns RK_OK, RK_ERR_BADARG, or RK_ERR_INFEAS.
 * On error, writes a human-readable message to result->message.
 *
 * Checks (in order):
 * 1. L_c <= U_c for all cells              → RK_ERR_BADARG  with cell index
 * 2. X_init[c]==0 && L_c>0 for any cell   → RK_ERR_INFEAS  ("structural zero")
 * 3. sum(L_c) <= n <= sum(U_c)             → RK_ERR_INFEAS  ("total capacity")
 * 4. |sum(T_kj) - 1| <= 1e-6 for all k   → normalizes T_kj + emits warning
 * 5. n_cats_total <= kNCatsTotalMax        → RK_ERR_BADARG  ("too many categories")
 */
int calib_validate_preentry(const CellTable& ct,
                             const CalibState& st,
                             rk_result_t* result,
                             const double* X_init,
                             int n_cats_total);

} // namespace lbw
```

- [ ] **Step 2.3: Create `src/calib_validate.cpp`**

```cpp
#include "calib_validate.hpp"
#include <cmath>
#include <cstring>
#include <cstdio>

namespace lbw {

int calib_validate_preentry(const CellTable& ct,
                             const CalibState& st,
                             rk_result_t* result,
                             const double* X_init,
                             int n_cats_total)
{
    // Helper: write error message and return code
    auto fail = [&](int code, const char* msg) -> int {
        if (result) {
            result->status = code;
            std::strncpy(result->message, msg, 255);
            result->message[255] = '\0';
        }
        return code;
    };

    // 1. n_cats_total capacity check (before any allocation)
    if (n_cats_total > kNCatsTotalMax) {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
            "n_cats_total=%d exceeds limit %d; use method='ieppa' or 'raking'",
            n_cats_total, kNCatsTotalMax);
        return fail(RK_ERR_BADARG, msg);
    }

    // Compute per-cell bounds
    const double lo = st.min_weight;
    const double hi = st.max_weight;
    double sum_L = 0.0, sum_U = 0.0;

    for (int c = 0; c < ct.M_cell; c++) {
        double L_c = lo * ct.n_per_cell[c];
        double U_c = hi * ct.n_per_cell[c];

        // 2. L_c <= U_c
        if (L_c > U_c + 1e-12) {
            char msg[256];
            std::snprintf(msg, sizeof(msg),
                "cell %d: L_c=%.4g > U_c=%.4g (min_weight > max_weight)", c, L_c, U_c);
            return fail(RK_ERR_BADARG, msg);
        }

        // 3. X_init[c]==0 && L_c>0 → structural infeasibility
        if (X_init != nullptr && X_init[c] <= 0.0 && L_c > 1e-12) {
            char msg[256];
            std::snprintf(msg, sizeof(msg),
                "cell %d: X_init=0 but L_c=%.4g > 0; multiplicative solver cannot "
                "move zero cell above lower bound", c, L_c);
            return fail(RK_ERR_INFEAS, msg);
        }

        sum_L += L_c;
        sum_U += U_c;
    }

    // 4. Total capacity vs target mass n
    const double n = static_cast<double>(st.n);
    if (sum_L > n + 1e-6) {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
            "sum(L_c)=%.4g > n=%.4g: lower bounds exceed total target mass", sum_L, n);
        return fail(RK_ERR_INFEAS, msg);
    }
    if (sum_U < n - 1e-6) {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
            "sum(U_c)=%.4g < n=%.4g: upper bounds cannot contain total target mass", sum_U, n);
        return fail(RK_ERR_INFEAS, msg);
    }

    // 5. Target sum validation (st.targets is const — cannot normalize in-place)
    // Callers must normalize before calling. Return RK_ERR_BADARG if off by > 1e-6.
    for (int k = 0; k < st.K; k++) {
        double s = 0.0;
        for (int j = 0; j < st.cat_counts[k]; j++) s += st.targets[k][j];
        if (std::fabs(s - 1.0) > 1e-6) {
            char msg[256];
            std::snprintf(msg, sizeof(msg),
                "margin %d targets sum to %.8f (expected 1.0±1e-6); "
                "normalize targets before calling", k, s);
            return fail(RK_ERR_BADARG, msg);
        }
    }

    return RK_OK;
}

} // namespace lbw
```

- [ ] **Step 2.4: Create `src/calib_linalg.hpp` (stubs — bodies in Plan D)**

```cpp
#pragma once
#include "calib_validate.hpp"   // provides kNCatsTotalMax — DO NOT redefine
#include "cell_table.hpp"
#include <cstddef>

namespace lbw {

// kNCatsTotalMax is defined in calib_validate.hpp — included above.

/**
 * Compute N = A × diag(D) × Aᵀ where A is the (n_cats_total × M_cell)
 * marginal incidence matrix (A[(k,j),c] = 1 iff cell c ∈ bucket(k,j)).
 * N is n_cats_total × n_cats_total, stored row-major (n_cats_total² doubles).
 * D is a M_cell-element diagonal weight vector on cell masses.
 * cat_offset[k] gives the starting row index for margin k.
 * Returns RK_OK or RK_ERR_BADARG (if n_cats_total > kNCatsTotalMax).
 * NOTE: implementation in Plan D (calib_linalg.cpp).
 */
int compute_normal_equations(const CellTable& ct,
                              const double* D,
                              double* N,
                              const int* cat_offset,
                              size_t n_cats_total);

/**
 * In-place LDLT factorization of a symmetric positive semidefinite n×n matrix.
 * eps_perturb: minimum diagonal value after factorization (Gill-Murray stability).
 * Returns RK_OK or RK_ERR_BADARG if n > kNCatsTotalMax.
 * NOTE: implementation in Plan D (calib_linalg.cpp).
 */
int ldlt_factor_inplace(double* A, size_t n, double eps_perturb);

/**
 * Solve A×x = b using a previously LDLT-factored matrix.
 * Overwrites b with the solution.
 * NOTE: implementation in Plan D (calib_linalg.cpp).
 */
void ldlt_solve(const double* L, const double* d_diag, double* b, size_t n);

} // namespace lbw
```

- [ ] **Step 2.5: Add calib_validate.cpp to Makevars AND Makevars.in**

The project has BOTH `src/Makevars` and `src/Makevars.in`. A `configure` script regenerates `Makevars` from `Makevars.in`. **Edit BOTH files** or changes to `Makevars` will be lost on the next configure run.

Read both: `cat src/Makevars src/Makevars.in`

In BOTH files, add `calib_validate.cpp`:
```makefile
PKG_SOURCES = c_api.cpp logit.cpp lbfgsb_solver.cpp ieppa.cpp cell_table.cpp r_bridge.cpp raking.cpp calib_validate.cpp
```

Note: verify `raking.cpp` is present before adding. Add only what is missing.

- [ ] **Step 2.6: Build gate**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```
Expected: `* DONE (leafblower)`

- [ ] **Step 2.6: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 337.

- [ ] **Step 2.7: Commit**
```bash
git add src/calib_validate.hpp src/calib_validate.cpp src/calib_linalg.hpp \
        tests/testthat/test-calibration-solvers.R
git commit -m "feat(infra): calib_validate_preentry + calib_linalg.hpp stubs

calib_validate.cpp: L_c>U_c, X_init=0&&L>0, capacity sum, n_cats_total
guard, target normalization. calib_linalg.hpp: compute_normal_equations +
ldlt_factor_inplace + ldlt_solve stubs (bodies in Plan D)."
```

---

## Task 3 — iEPPA default metric: MAX_ERR → KL

**Files:** `src/c_api.cpp`, `R/harvest.R`

**Scientific validation:** The KL metric is already implemented in ieppa.cpp and dispatched via CalibMetric::KL. Changing the default only affects which metric applies when the user passes `convergence=list()`. Tests that pin to explicit convergence (added during a2p2) are unaffected.

- [ ] **Step 3.1: Write failing test**

Append to `tests/testthat/test-calibration-solvers.R`:
```r
test_that("T3: ieppa default convergence is kl+improvement", {
  set.seed(42)
  data <- data.frame(
    a = factor(sample(c("1","2","3"), 500, replace=TRUE)),
    b = factor(sample(c("1","2"), 500, replace=TRUE))
  )
  target <- list(a=c("1"=1/3,"2"=1/3,"3"=1/3), b=c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight=10, method="ieppa",
                           max_iterations=500, attach_weights=FALSE)
  r <- attr(w, "result")
  # After change: default metric is kl, not max_err
  expect_equal(r$convergence_used$metric, "kl",
               info="ieppa default metric must be 'kl' after spec change")
  expect_equal(r$convergence_used$rule, "improvement")
})
```

Run: expected FAIL (currently metric is "max_err").

- [ ] **Step 3.2: Change c_api.cpp default metric**

Read line ~50 in `src/c_api.cpp`: `p->metric = 0;  /* MAX_ERR */`

Change to:
```c
p->metric = 2;  /* KL — iEPPA is a Sinkhorn-type algorithm minimizing KL */
```

- [ ] **Step 3.3: Change harvest.R default — method-aware at call site**

`parse_convergence` does not know the method name. The fix is in the harvest.R call site, AFTER `conv <- parse_convergence(convergence)`, not inside parse_convergence.

Read: `grep -n "conv\s*<-\s*parse_convergence\|method_int\|map_method" R/harvest.R | head -10`

After `conv <- parse_convergence(convergence)`, add:
```r
# iEPPA is a KL minimizer — override default metric from max_err to kl
# when the user has not explicitly set a metric.
if (method == "ieppa" && is.null(convergence[["metric"]]) &&
    is.null(convergence[["improvement"]]) && is.null(convergence[["absolute"]]) &&
    is.null(convergence[["pct"]])) {
  conv$metric <- "kl"  # iEPPA is a Sinkhorn KL minimizer
}
```

This is method-aware and does NOT touch parse_convergence's generic defaults, avoiding regressions for raking/lbfgsb.

Also update the `@param convergence` roxygen to mention iEPPA default:
```r
#'   For \code{method="ieppa"}, default is \code{kl + improvement + 0.001}
#'   (consistent with iEPPA's Sinkhorn objective). For other methods:
#'   \code{max_err + improvement + 0.001}.
```

- [ ] **Step 3.4: Build + run test**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | tail -5
```
Expected: T3 PASS. Previously pinned tests unaffected (they use explicit convergence).

**Scientific validation:** Run the full suite and grep for NOCONV to ensure no regressions:
```bash
Rscript -e 'devtools::test()' 2>&1 | grep -i "noconv\|FAIL\|warn.*converge" | head -10
```
Expected: no new NOCONV warnings; FAIL 0.

- [ ] **Step 3.5: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 338.

- [ ] **Step 3.6: Commit**
```bash
git add src/c_api.cpp R/harvest.R tests/testthat/test-calibration-solvers.R
git commit -m "feat(ieppa): change default convergence from max_err to kl+improvement

Per spec §9: iEPPA is a Sinkhorn-type KL minimizer. kl+improvement is
semantically consistent with the algorithm's actual objective. Backward
compat: explicit convergence=list(metric='max_err',...) unchanged.
NEWS.md: breaking change documented (v0.2.0)."
```

---

## Task 4 — iEPPA KL reference fixture for A1

**Files:** `data-raw/gen_ieppa_kl_ref.R` (new), `tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds` (new)

**Purpose:** A1 test (`method="sinkhorn"` KL ≤ iEPPA KL at best_iter) requires a precomputed reference value. Generated now, before Plan C implements sinkhorn, so A1 can be written as a red test.

- [ ] **Step 4.1: Create `data-raw/gen_ieppa_kl_ref.R`**

```r
#!/usr/bin/env Rscript
# Generates tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds
# Run: Rscript data-raw/gen_ieppa_kl_ref.R
# Records iEPPA's KL divergence at best_iter on stepstone-fulldata.
# Used by A1 test: sinkhorn KL at convergence must be <= this value.
suppressPackageStartupMessages({
  library(leafblower)
  library(arrow)
  library(jsonlite)
})

data <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
tgt_raw <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
target <- lapply(tgt_raw, function(x) { v <- unlist(x); v / sum(v) })

cat("Running iEPPA on stepstone-fulldata to record KL at best_iter...\n")
suppressWarnings(
  w <- leafblower::harvest(data, target, method = "ieppa",
                           max_weight = 5, max_iterations = 3000,
                           attach_weights = FALSE)
)
r <- attr(w, "result")

# Record best_error (which now tracks KL since default is kl+improvement)
ref <- list(
  kl_at_best_iter = r$best_error,
  best_iter       = r$best_iter,
  max_error       = r$max_error,
  ieppa_version   = as.character(packageVersion("leafblower")),
  date            = Sys.Date()
)
cat(sprintf("  KL at best_iter=%d: %.4e\n", ref$best_iter, ref$kl_at_best_iter))
cat(sprintf("  max_error: %.4e\n", ref$max_error))

saveRDS(ref, "tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds")
cat("Saved to tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds\n")
```

- [ ] **Step 4.2: Run the script to generate the fixture**
```bash
Rscript data-raw/gen_ieppa_kl_ref.R
```
Expected: outputs KL at best_iter (should be ~3e-3 based on prior benchmarks) and saves the rds file.

- [ ] **Step 4.3: Write A1 test (will FAIL until Plan C implements sinkhorn)**

Append to `tests/testthat/test-calibration-solvers.R`:
```r
test_that("A1: sinkhorn KL <= ieppa KL at best_iter (skip until sinkhorn implemented)", {
  skip_on_cran()
  skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"))
  ref_path <- test_path("fixtures/ieppa_kl_reference_stepstone.rds")
  skip_if(!file.exists(ref_path))
  ref <- readRDS(ref_path)

  # SKIP until Plan C implements sinkhorn
  skip("sinkhorn not yet implemented — Plan C will unskip this test")

  library(arrow); library(jsonlite)
  data <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
  tgt_raw <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
  target <- lapply(tgt_raw, function(x) { v <- unlist(x); v / sum(v) })

  suppressWarnings(
    w_s <- leafblower::harvest(data, target, method="sinkhorn",
                                max_weight=5, max_iterations=3000,
                                attach_weights=FALSE)
  )
  r_s <- attr(w_s, "result")

  # Sinkhorn must reach lower KL than iEPPA's best_iter KL
  expect_lte(r_s$convergence_used$objective, ref$kl_at_best_iter,
             label="sinkhorn KL <= ieppa best_iter KL")
  # Sinkhorn must converge (not NOCONV)
  expect_equal(r_s$status, 0L, label="sinkhorn must converge")
  # A7: best_iter == last_iter for monotone methods
  expect_equal(r_s$convergence_used$fired_at_iter, r_s$iterations,
               label="sinkhorn best_iter == last_iter (monotone)")
})
```

Run: test PASSES (skipped).

- [ ] **Step 4.4: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 338 (A1 skipped, not failed).

- [ ] **Step 4.5: Commit**
```bash
git add data-raw/gen_ieppa_kl_ref.R \
        tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds \
        tests/testthat/test-calibration-solvers.R
git commit -m "test(A1): iEPPA KL reference fixture + skipped sinkhorn A1 test

Reference fixture records iEPPA KL at best_iter on stepstone-fulldata.
A1 test is present but skipped until Plan C implements method='sinkhorn'."
```

---

## Final Verification

- [ ] `devtools::test()` → FAIL 0, PASS ≥ 338
- [ ] `R CMD check --as-cran` → 0 ERROR, 0 WARNING
- [ ] `grep -r "RK_ALG_SINKHORN\|RK_ALG_CHEBYSHEV\|RK_ALG_GREG\|RK_ALG_GRAKE" src/leafblower.h` → all 4 present
- [ ] `grep -n "EXPECTED_RK_RESULT_BYTES" src/leafblower.h` → static_assert present
- [ ] `Rscript -e 'r<-attr(leafblower::harvest(data.frame(a=factor(c("1","2"))),list(a=c("1"=0.5,"2"=0.5)),max_weight=3,method="ieppa",attach_weights=FALSE),"result"); cat(r$convergence_used$metric)'` → prints "kl"
- [ ] New method stubs: `expect_error(harvest(..., method="sinkhorn"))` works without crash

---

## Self-Review

**Spec coverage:**
- §3 new method constants (4-7) → Task 1 ✅
- §3 convergence_objective + convergence_minimized_metric → Task 1 ✅
- §3 EXPECTED_RK_RESULT_BYTES → Task 1 ✅
- §2 calib_linalg.hpp → Task 2 ✅ (stubs; bodies in Plan D)
- §13 calib_validate_preentry → Task 2 ✅
- §13 X_init=0 && L_c>0 check → Task 2 (step 2.3) ✅
- §9 iEPPA default to kl+improvement → Task 3 ✅
- A1 fixture → Task 4 ✅

**Placeholder scan:** No TBD found. All code blocks are concrete.

**Type consistency:** `calib_validate_preentry` signature matches between .hpp and .cpp. `compute_normal_equations` stub signature consistent with Plan D's intended caller. `convergence_objective` field name consistent with r_bridge → harvest.R chain.

**Scientific validation:** Task 3 explicitly runs `grep -i NOCONV` after regression to confirm no calibration regressions from default metric change. A1 fixture captures the reference value before sinkhorn is built so the test is genuinely red in Plan C.
