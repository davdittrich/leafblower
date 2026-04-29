# Chebyshev IPM Rewrite + GREG Quality Warning

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix chebyshev/grake NaN convergence failures and greg 11% error on K=9 overlapping margins. Chebyshev gets ieppa warm-start + Mehrotra predictor-corrector + Jacobi preconditioning. Greg gets quality warning.

**Architecture:** Three-layer chebyshev fix: (1) obs-level ieppa warm-start passed into chebyshev_ipm via w_warm_obs parameter, aggregated to cell masses internally after CellTable build (avoids double build); (2) Mehrotra predictor-corrector replacing fixed-barrier-reduction with near-quadratic convergence; (3) Jacobi diagonal preconditioning on normal equations for K=9 ill-conditioning. Greg fix: single harvest.R warning when max_err > 0.05.

**Tech Stack:** C++17 (chebyshev.cpp/hpp, r_bridge.cpp), R (harvest.R), testthat3.

**Spec:** `docs/superpowers/specs/2026-04-29-chebyshev-greg-fix.md`

**Failing baselines (must verify in Task 1):**
- `greg` on stepstone K=9 → max_err = 1.10e-01 (wrong solution, status=0)
- `chebyshev` on stepstone K=9 → NaN, status=RK_ERR_NOCONV (500-iter cap hit)
- `grake` on stepstone K=9 → NaN, status=RK_ERR_NOCONV (chebyshev wrapper)

**Forbidden alternatives:**
- Do NOT replace LDLT with iterative CG (Mehrotra needs exact factor reuse).
- Do NOT add a homotopy-on-mu schedule on top of Mehrotra (μ becomes adaptive — schedule destroys it).
- Do NOT silently fall back from chebyshev to raking when warm-start fails (return finite cold-start result).
- Do NOT remove the existing schur_nu diagnostic (verbose>=2 logging stays).
- Do NOT mutate `st.weights` in r_bridge.cpp pre-solve dispatch (use weights_copy).

**Audit (spy/mock strategy):**
- Verify ieppa pre-solve runs by snapshot-checking `verbose>=2` log line "chebyshev: warm-start from ieppa max_err=...".
- Verify Mehrotra two-phase by counting LDLT factorizations: exactly 1 per outer iter (not 2). Add a `n_factorizations` counter to ChebyshevResult.diagnostics.
- Verify Jacobi by checking `cond(N_scaled) < cond(N_raw)` at iter 0 (verbose>=3 dumps both).

---

## Task 1: RED tests — failing-first baseline

**Mechanism:** testthat3 expectations on greg quality warning + chebyshev finite convergence.
**Forbidden:** Skipping the RED confirmation step. Writing GREEN tests before implementation.
**Audit:** `devtools::test(filter="chebyshev")` output shows 3 FAIL/ERROR before Task 2.

**Files:**
- NEW: `tests/testthat/test-chebyshev.R`

**Steps:**

- [ ] 1.1 Read `tests/testthat/test-harvest.R` (for fixture conventions) and `tests/testthat/test-calibration-solvers.R` (for solver-call patterns); list the existing `expect_*` idioms used for `harvest()` calls.

- [ ] 1.2 Write `tests/testthat/test-chebyshev.R` with three tests:

```r
# tests/testthat/test-chebyshev.R
# RED tests for chebyshev IPM rewrite + greg quality warning.
# T_greg_warn  : greg on tight bounds emits "greg.*unreliable" warning
# T_cheby_warm : chebyshev K=3 returns finite weights with status=0 (warm-start)
# T_cheby_warm_fallback : chebyshev still returns finite result if warm-start degenerates

library(testthat)

# --- Tiny synthetic fixture: K=3 overlapping margins, n=200, tight bounds ---
.make_tight_fixture <- function(n = 200L, seed = 42L) {
  set.seed(seed)
  cat1 <- sample.int(4L, n, replace = TRUE)
  cat2 <- sample.int(3L, n, replace = TRUE)
  cat3 <- sample.int(5L, n, replace = TRUE)
  design <- data.frame(cat1 = factor(cat1), cat2 = factor(cat2), cat3 = factor(cat3))
  d <- rep(1.0, n)
  # Targets perturb uniform proportions enough to force tight-bounds activation
  t1 <- as.numeric(table(cat1)) * 1.20
  t2 <- as.numeric(table(cat2)) * 0.85
  t3 <- as.numeric(table(cat3)) * 1.10
  target <- list(cat1 = t1, cat2 = t2, cat3 = t3)
  list(design = design, d = d, target = target)
}

test_that("T_greg_warn: greg emits unreliability warning on tight bounds", {
  skip_if_not_installed("leafblower")
  fx <- .make_tight_fixture()
  expect_warning(
    leafblower::harvest(
      design = fx$design, d = fx$d, target = fx$target,
      method = "greg",
      max_weight = 1.5,   # tight bounds: forces greg into >5% regime
      max_iterations = 100L, verbose = 0L
    ),
    regexp = "greg.*unreliable",
    fixed = FALSE
  )
})

test_that("T_cheby_warm: chebyshev returns finite weights with status=0", {
  skip_if_not_installed("leafblower")
  fx <- .make_tight_fixture()
  res <- leafblower::harvest(
    design = fx$design, d = fx$d, target = fx$target,
    method = "chebyshev",
    max_weight = 5.0,
    max_iterations = 200L, verbose = 0L
  )
  expect_true(all(is.finite(res$weights)),
              info = "chebyshev returned non-finite weights (NaN/Inf)")
  expect_equal(res$status, 0L,
               info = sprintf("chebyshev status=%d (expect 0); message=%s",
                              res$status %||% -1L, res$message %||% "<NA>"))
  # Quality bound: chebyshev should beat or tie raking on max_err
  ref <- leafblower::harvest(
    design = fx$design, d = fx$d, target = fx$target,
    method = "raking", max_weight = 5.0,
    max_iterations = 500L, verbose = 0L
  )
  expect_lte(res$max_error, ref$max_error * 1.05,
             label = "chebyshev max_err should be <= raking max_err * 1.05")
})

test_that("T_cheby_warm_fallback: chebyshev finite even with weak warm-start", {
  skip_if_not_installed("leafblower")
  fx <- .make_tight_fixture(n = 80L, seed = 7L)  # smaller n: ieppa warm-start degenerates faster
  res <- leafblower::harvest(
    design = fx$design, d = fx$d, target = fx$target,
    method = "chebyshev",
    max_weight = 5.0,
    max_iterations = 50L,    # short budget: ieppa pre-solve will produce poor warm-start
    verbose = 0L
  )
  expect_true(all(is.finite(res$weights)),
              info = "chebyshev returned non-finite even with degraded warm-start")
})

# Optional K=9 stepstone test — skip if fixture absent
test_that("T_cheby_K9: chebyshev finite on stepstone K=9 (optional)", {
  skip_if_not_installed("leafblower")
  fx_path <- system.file("fixtures", "stepstone_k9.rds", package = "leafblower")
  skip_if(!nzchar(fx_path) || !file.exists(fx_path), "stepstone K=9 fixture absent")
  fx <- readRDS(fx_path)
  res <- leafblower::harvest(
    design = fx$design, d = fx$d, target = fx$target,
    method = "chebyshev", max_weight = 5.0,
    max_iterations = 200L, verbose = 0L
  )
  expect_true(all(is.finite(res$weights)))
  expect_equal(res$status, 0L)
})
```

- [ ] 1.3 Build current code (no implementation changes yet):
```bash
R CMD INSTALL --preclean .
```
Expected: clean build, no warnings, exit 0.

- [ ] 1.4 Run RED suite:
```bash
Rscript -e 'devtools::test(filter = "chebyshev")'
```
Expected outputs:
- `T_greg_warn` → FAIL (no warning emitted, current greg returns silently).
- `T_cheby_warm` → FAIL (`status != 0` because chebyshev hits NaN / RK_ERR_NOCONV).
- `T_cheby_warm_fallback` → FAIL (likely non-finite from cold-start divergence).
- `T_cheby_K9` → SKIPPED (fixture absent in package install).

- [ ] 1.5 Verify RED state. If any test passes accidentally, halt and diagnose — the test does not pin failing behavior. Do not proceed until 3 of 3 (T_greg_warn, T_cheby_warm, T_cheby_warm_fallback) FAIL.

- [ ] 1.6 Commit:
```bash
git add tests/testthat/test-chebyshev.R
git commit -m "test(chebyshev): add RED suite for warm-start + greg warning"
```

**Done when:** 3 tests FAIL with informative messages; 1 SKIPPED; build clean.

---

## Task 2: GREG quality warning in `R/harvest.R`

**Mechanism:** Single conditional `warning()` after greg dispatch when `max_error > 0.05`.
**Forbidden:** Erroring out, silently rerunning with another method, modifying greg C++ code.
**Audit:** Spy: assert exact regex `greg.*unreliable` matches; assert `result$weights` still returned (warn-and-continue, not abort).

**Files:**
- `R/harvest.R`

**Steps:**

- [ ] 2.1 Read `R/harvest.R` lines 350-410 (where `calib_result` is assembled post-C++ dispatch) to find the exact return-site conventions and the variable name holding the chosen method (line 530 has `match.arg`).

- [ ] 2.2 Locate the post-dispatch block where `calib_result$max_error` is finalized (after the `convergence_used` list is attached but before `return(calib_result)`).

- [ ] 2.3 Add the warning block immediately before `return(calib_result)`:

```r
# greg quality check: warn (do not abort) when max_err > 5%.
# greg's regression model is unreliable for K>=5 overlapping margins or tight bounds
# because the design Gram matrix becomes near-singular. We keep the result and warn.
if (identical(method, "greg") &&
    !is.null(calib_result$max_error) &&
    is.finite(calib_result$max_error) &&
    calib_result$max_error > 0.05) {
  warning(
    sprintf(
      paste0("greg converged but max_err=%.4g (>5%%). ",
             "greg may be unreliable for K=%d margins or tight bounds ",
             "(max_weight=%.4g). Consider method='raking' or 'ieppa'."),
      calib_result$max_error,
      length(target),
      max_weight
    ),
    call. = FALSE
  )
}
```

- [ ] 2.4 Compile + reinstall:
```bash
R CMD INSTALL --preclean .
```
Expected: clean build.

- [ ] 2.5 Run targeted test:
```bash
Rscript -e 'devtools::test(filter = "chebyshev")'
```
Expected:
- `T_greg_warn` → PASS.
- `T_cheby_warm`, `T_cheby_warm_fallback` → still FAIL (chebyshev untouched).

- [ ] 2.6 Run full harvest suite to catch regressions:
```bash
Rscript -e 'devtools::test(filter = "harvest")'
```
Expected: no new failures.

- [ ] 2.7 Commit:
```bash
git add R/harvest.R
git commit -m "feat(greg): warn when max_err > 5% (likely unreliable solution)"
```

**Done when:** `T_greg_warn` GREEN; chebyshev tests still RED; harvest suite no regressions.

---

## Task 3: Warm-start infrastructure (signature + dispatch plumbing)

**Mechanism:** `r_bridge.cpp` runs ieppa pre-solve into a deep-copied weights buffer; passes obs-level result + delta_warm to chebyshev_ipm via new optional parameters.
**Forbidden:** Mutating `st.weights` directly. Calling `ieppa_solve(st)` without a copy. Building CellTable twice (aggregation must happen inside chebyshev after its own build).
**Audit:** Compile gate only — initialization logic in Task 4. Verify by inspection: `weights_copy` lifetime ≥ ieppa_solve scope; `st.weights` pointer unchanged after pre-solve block.

**Files:**
- `src/chebyshev.hpp`
- `src/chebyshev.cpp` (signature update only, no logic yet)
- `src/r_bridge.cpp`

**Steps:**

- [ ] 3.1 Read `src/chebyshev.hpp` (43 lines, full file) and `src/r_bridge.cpp` lines 480-510 (chebyshev/grake dispatch site near line 488 where `greg_solve` is called — chebyshev dispatch is in the same conditional ladder).

- [ ] 3.2 Read `src/ieppa.hpp` lines 1-60 to confirm `IEPPAResult { double max_error; std::vector<double> best_weights; ... }` and `IEPPAResult ieppa_solve(CalibState&)`.

- [ ] 3.3 Update `src/chebyshev.hpp`:

```cpp
// chebyshev.hpp — UPDATED signature
#pragma once
#include "calib_state.hpp"
#include <vector>

namespace lbw {

enum class LpVariant { CHEBYSHEV, GRAKE };

struct ChebyshevResult {
    int    status;
    int    iterations;
    double max_error;
    int    n_factorizations;   // Mehrotra audit counter (one LDLT per outer iter)
    char   message[256];
    std::vector<double> best_weights;   // obs-level (matches existing ChebyshevResult field name)
};

ChebyshevResult chebyshev_ipm(
    CalibState& st,
    LpVariant   variant,
    const std::vector<double>& w_warm_obs = {},   // obs-level warm weights (size n);
                                                  //   empty -> cold start (uniform).
                                                  //   chebyshev aggregates obs->cell internally
                                                  //   AFTER its own build_cell_table.
    double      delta_warm = -1.0                 // initial delta_0; -1 -> default 1.0
);

inline ChebyshevResult chebyshev_solve(CalibState& st) {
    return chebyshev_ipm(st, LpVariant::CHEBYSHEV);
}
inline ChebyshevResult chebyshev_solve_warm(
    CalibState& st,
    const std::vector<double>& w_warm_obs,
    double delta_warm)
{
    return chebyshev_ipm(st, LpVariant::CHEBYSHEV, w_warm_obs, delta_warm);
}
inline ChebyshevResult grake_solve(CalibState& st) {
    return chebyshev_ipm(st, LpVariant::GRAKE);
}
inline ChebyshevResult grake_solve_warm(
    CalibState& st,
    const std::vector<double>& w_warm_obs,
    double delta_warm)
{
    return chebyshev_ipm(st, LpVariant::GRAKE, w_warm_obs, delta_warm);
}

} // namespace lbw
```

- [ ] 3.4 Update `src/chebyshev.cpp` function definition signature only (no init logic — added in Task 4):

```cpp
ChebyshevResult chebyshev_ipm(
    CalibState& st,
    LpVariant variant,
    const std::vector<double>& w_warm_obs,
    double delta_warm)
{
    // ... existing body unchanged for now ...
    // Initialize Mehrotra audit counter:
    ChebyshevResult res;
    res.n_factorizations = 0;
    // ... rest of existing body ...
}
```

Add `#include "ieppa.hpp"` only if not already present (it is NOT needed in chebyshev.cpp itself — pre-solve happens in r_bridge.cpp).

- [ ] 3.5 Update `src/r_bridge.cpp` chebyshev/grake dispatch. **CRITICAL**: The actual code does NOT use bare `else if (strcmp...)` branches. It uses a `dispatch_cheb` lambda at ~line 526:
  ```cpp
  auto dispatch_cheb = [&](lbw::LpVariant variant, int alg_code) {
      auto res = lbw::chebyshev_ipm(st, variant);   // line 527 — this is what to replace
      ...
  };
  if (strcmp(method_str, "chebyshev") == 0)      dispatch_cheb(LpVariant::CHEBYSHEV, ...);
  else if (strcmp(method_str, "grake") == 0)     dispatch_cheb(LpVariant::GRAKE, ...);
  ```
  Verify with: `grep -n "dispatch_cheb\|LpVariant::CHEBYSHEV\|LpVariant::GRAKE" src/r_bridge.cpp | head -10`
  
  Modify `dispatch_cheb` to accept `w_warm_obs` and `delta_warm`, or replace it with an inline block. The warm-start pre-solve runs BEFORE both chebyshev and grake dispatch (both use the same ieppa pre-solve). Restructure as:

```cpp
} else if (strcmp(method_str, "chebyshev") == 0 ||
           strcmp(method_str, "grake")     == 0) {
    const lbw::LpVariant variant =
        (strcmp(method_str, "chebyshev") == 0) ? lbw::LpVariant::CHEBYSHEV
                                               : lbw::LpVariant::GRAKE;

    // Warm-start from ieppa pre-solve.
    // CRITICAL: ieppa_solve(st) MUTATES st.weights in-place.
    // We MUST hand it a copy so the caller's CalibState is untouched.
    std::vector<double> w_warm_obs;     // empty => cold start
    double delta_warm = -1.0;           // -1 => use default delta_0=1.0
    {
        std::vector<double> weights_copy(st.weights, st.weights + st.n);
        CalibState st_warm = st;                              // shallow copy of struct
        st_warm.weights = weights_copy.data();                // redirect pointer to copy
        // Cap pre-solve budget: 10% of total, clamped to [5, 100]
        st_warm.inner_max_iter = std::max(5,
                                  std::min(100, st.inner_max_iter / 10));
        auto ieppa_res = lbw::ieppa_solve(st_warm);
        // SAFETY: st_warm holds dangling pointer once weights_copy goes out of scope.
        //         st_warm MUST NOT escape this block. weights_copy still alive
        //         only until block close. We move best_weights out before that.
        if (!ieppa_res.best_weights.empty() &&
            static_cast<int>(ieppa_res.best_weights.size()) == st.n &&
            std::isfinite(ieppa_res.max_error)) {
            w_warm_obs = std::move(ieppa_res.best_weights);   // size n
            delta_warm = ieppa_res.max_error * 1.5;           // start delta slightly above ieppa quality
        }
        // weights_copy and st_warm destroyed here.
        // st.weights pointer + st.* unchanged.
    }

    auto res = lbw::chebyshev_ipm(st, variant, w_warm_obs, delta_warm);
    // ... existing result-handling code (status copy, weights copy, message, etc.) ...
}
```

NOTE: preserve all surrounding diagnostics/return-payload code that was already in the chebyshev/grake branches. The only change is wrapping with the warm-start block and switching from `chebyshev_solve(st)` / `grake_solve(st)` to the unified `chebyshev_ipm(st, variant, w_warm_obs, delta_warm)` call.

- [ ] 3.6 Compile gate:
```bash
R CMD INSTALL --preclean .
```
Expected: clean build, no `-Wmaybe-uninitialized` or pointer-lifetime warnings on the weights_copy block. If clang/gcc emits a dangling-pointer warning, halt and reread step 3.5 — `st_warm` must NOT be passed downstream.

- [ ] 3.7 Run tests (signature wired but warm-start init still missing):
```bash
Rscript -e 'devtools::test(filter = "chebyshev")'
```
Expected:
- `T_greg_warn` → PASS (Task 2 stays GREEN).
- `T_cheby_warm` → still FAIL (chebyshev_ipm receives w_warm_obs but doesn't use it yet).
- `T_cheby_warm_fallback` → still FAIL.
- Full ieppa suite must NOT regress:
```bash
Rscript -e 'devtools::test(filter = "ieppa")'
```
Expected: no new failures (st.weights untouched proves the copy works).

- [ ] 3.8 Commit:
```bash
git add src/chebyshev.hpp src/chebyshev.cpp src/r_bridge.cpp
git commit -m "feat(chebyshev): plumb warm-start params (no-op until Task 4)

- Add w_warm_obs and delta_warm optional params to chebyshev_ipm.
- r_bridge.cpp pre-solves ieppa into a weights_copy (st.weights untouched).
- Wire grake to same path via LpVariant.
- Signature change only; init logic in Task 4."
```

**Done when:** Build clean; ieppa tests still GREEN; chebyshev tests still RED.

---

## Task 4: Warm-start initialization inside `chebyshev_ipm`

**Mechanism:** Aggregate obs-level `w_warm_obs` to cell masses using the existing `ct.cell_of[]` (built by `build_cell_table` already inside chebyshev_ipm), then mass-preserving clamp to `[L_cell, U_cell]`. Use `delta_warm` for `delta_0`.
**Forbidden:** Building a second CellTable. Skipping the mass-preservation rescale (would violate Σ X = N constraint). Allowing X[c] outside [L_cell, U_cell] after init.
**Audit:** Verify post-init: (a) `Σ X_primal == Σ st.d` ± 1e-9·n; (b) `L_cell[c] ≤ X[c] ≤ U_cell[c]` for all c; (c) at verbose>=2, log "chebyshev: warm-start total_pre=%g total_post=%g".

**Files:**
- `src/chebyshev.cpp`

**Steps:**

- [ ] 4.1 Read `src/chebyshev.cpp` lines 30-135 (initialization region: `build_cell_table` at line 32, `L_cell/U_cell` at line 41, existing `X_init` cold-start at line 109, `X` clamped init at line 125).

- [ ] 4.2 Locate the cold-start block (lines ~108-132). The current code does:
```cpp
std::vector<double> X_init(ct.M_cell, 0.0);
for (int i = 0; i < st.n; i++) X_init[ct.cell_of[i]] += st.weights[i];
// ... then clamps X[c] = clamp(X_init[c], L_cell+eps, U_cell-eps) ...
```

- [ ] 4.3 Replace the X_init build with warm-start-aware version:

```cpp
// --- Warm-start aggregation: obs-level w_warm_obs -> cell masses ---
// If w_warm_obs is empty or wrong size, fall back to st.weights cold start.
std::vector<double> X_init(ct.M_cell, 0.0);
const bool have_warm = (!w_warm_obs.empty() &&
                        static_cast<int>(w_warm_obs.size()) == st.n);
const double* w_src = have_warm ? w_warm_obs.data() : st.weights;
for (int i = 0; i < st.n; i++) X_init[ct.cell_of[i]] += w_src[i];

if (have_warm) {
    // Mass-preserving clamp:
    //   1) record total_pre = Σ X_init (== Σ w_warm_obs == Σ st.d, by ieppa invariant)
    //   2) clamp each cell to [L_cell, U_cell]
    //   3) rescale so Σ X == total_pre (preserve total mass)
    //   4) clamp again (rescale can push out of bounds)
    double total_pre = 0.0;
    for (int c = 0; c < ct.M_cell; c++) total_pre += X_init[c];
    for (int c = 0; c < ct.M_cell; c++)
        X_init[c] = std::clamp(X_init[c], L_cell[c], U_cell[c]);
    double total_post = 0.0;
    for (int c = 0; c < ct.M_cell; c++) total_post += X_init[c];
    if (total_post > 1e-12 && total_pre > 1e-12) {
        const double scale = total_pre / total_post;
        for (int c = 0; c < ct.M_cell; c++)
            X_init[c] = std::clamp(X_init[c] * scale, L_cell[c], U_cell[c]);
    }
    if (st.verbose >= 2) {
        char msg[160];
        double total_final = 0.0;
        for (int c = 0; c < ct.M_cell; c++) total_final += X_init[c];
        std::snprintf(msg, sizeof(msg),
            "chebyshev: warm-start total_pre=%.6e total_post=%.6e total_final=%.6e",
            total_pre, total_post, total_final);
        st.log_fn(msg);
    }
}
// (Cold-start path falls through unchanged: X_init holds Σ st.weights per cell.)

// Existing eps-shifted clamp into X (unchanged):
//   double eps_shift = ...;
//   for (int c = 0; c < ct.M_cell; c++) {
//       double gap = U_cell[c] - L_cell[c];
//       if (gap < 2.0*eps_shift) X[c] = 0.5*(L_cell[c] + U_cell[c]);
//       else X[c] = std::clamp(X_init[c], L_cell[c]+eps_shift, U_cell[c]-eps_shift);
//   }
```

- [ ] 4.4 Replace `delta_0 = 1.0` initialization with the warm-aware version. Locate the line that sets the initial delta (look for `mu = 1.0` near line 159 — delta init is in the same neighborhood; if delta is stored in a scalar `delta` variable, find it via `grep -n "delta_0\|delta =" src/chebyshev.cpp`):

```cpp
// delta initialization: warm if positive, else cold (1.0)
const double delta_0 = (delta_warm > 0.0 && std::isfinite(delta_warm))
                       ? std::clamp(delta_warm, 1e-6, 1.0)
                       : 1.0;
double delta = delta_0;
```

- [ ] 4.5 Compile gate:
```bash
R CMD INSTALL --preclean .
```
Expected: clean build.

- [ ] 4.6 Run RED tests:
```bash
Rscript -e 'devtools::test(filter = "chebyshev")'
```
Expected:
- `T_greg_warn` → PASS.
- `T_cheby_warm_fallback` → PASS (warm-start with weak ieppa now gives a finite starting X; single-step IPM may still hit budget but does not produce NaN — finiteness gate satisfied).
- `T_cheby_warm` → may still FAIL on `status==0` and `max_err <= raking*1.05` (single-step barrier still too slow on K=3 tight) — Task 5 fixes this. Confirm the failure is now "status=RK_ERR_BUDGET" or similar non-NaN failure, not a NaN.

- [ ] 4.7 Run regression suite:
```bash
Rscript -e 'devtools::test()'
```
Expected: no new failures outside chebyshev tests.

- [ ] 4.8 Commit:
```bash
git add src/chebyshev.cpp
git commit -m "feat(chebyshev): aggregate ieppa warm-start into cell masses

- obs-level w_warm_obs -> X_init via existing ct.cell_of[]
- mass-preserving clamp to [L_cell, U_cell] (clamp -> rescale -> reclamp)
- delta_0 = ieppa_max_err * 1.5 (clamped to [1e-6, 1.0]) when warm
- verbose>=2 logs total_pre/post/final for audit"
```

**Done when:** `T_cheby_warm_fallback` GREEN; `T_cheby_warm` failure mode shifts from NaN to non-convergence within budget; build clean.

---

## Task 5: Mehrotra predictor-corrector + Jacobi preconditioning

**Mechanism:** Replace single-step barrier reduction (μ *= 0.1) with two-phase Mehrotra: (Phase A) affine predictor with Jacobi-scaled LDLT factor; (Phase B) corrector reusing same LDLT with second-order RHS. Adaptive μ via complementarity arithmetic — no schedule.
**Forbidden:** Refactoring N twice per iteration. Applying Jacobi to Phase A then unscaling and re-Jacobi for Phase B (must scale once, reuse factor). Using fixed σ. Allowing σ ≤ 0 or σ > 1. Skipping the m==0 guard (degenerate LP must fall through to plain barrier μ *= 0.1).
**Audit:** ChebyshevResult.n_factorizations counter must equal `iterations` (one factor per outer iter, two solves). Verbose>=3 dumps `cond(N_raw)` and `cond(N_scaled)` at iter 0 — Jacobi must reduce condition number.

**Files:**
- `src/chebyshev.cpp`

**Steps:**

- [ ] 5.1 Read `src/chebyshev.cpp` lines 150-410 (the existing IPM loop body: `mu = 1.0` init line 159, sigma_mu line 246, rhs assembly lines 269-289, Sherman-Morrison line 289, schur_nu diagnostic lines 322-356, residual computation lines 372-401).

- [ ] 5.2 Identify the LDLT factor + solve site. Currently the code factorizes a `nct_red × nct_red` matrix (call it `N_red`) once per outer iter then solves once. We will keep the factor-once pattern but:
  (a) fold Jacobi scaling into the matrix BEFORE factorization,
  (b) build a separate Phase-A RHS,
  (c) solve A → unscale → compute α_aff, μ_aff, σ,
  (d) build Phase-B RHS (with second-order term),
  (e) solve B reusing the same factor → unscale → take corrected step.

- [ ] 5.3 Add the Mehrotra + Jacobi block. Replace the single-step IPM iteration body. The skeleton (preserve existing variable names where they exist; `nct_red`, `red_to_full`, `is_ref`, `s_up`, `s_dn`, `s_lo`, `s_hi`, `y_up`, `y_dn`, `y_lo`, `y_hi`, `D_marg`, `w_kj`, etc. all stay):

```cpp
// === Mehrotra predictor-corrector with Jacobi preconditioning =================
// Variables in scope from earlier in chebyshev_ipm:
//   nct_red, red_to_full, is_ref, full_to_red
//   X (size M_cell), delta (scalar)
//   s_up, s_dn, s_lo, s_hi, y_up, y_dn, y_lo, y_hi, s_delta, y_delta
//   w_kj, T_flat, S, D_marg, etc.
// Existing helper builders (build_N_red, build_rhs, ldlt_factor, ldlt_solve, ...)
// must already exist or be lifted from current code; no new linear-algebra deps.

// --- m: number of complementarity pairs ---
//   2*nct (s_up,s_dn) + 2*M_cell (s_lo,s_hi) + 1 (s_delta)
const int m_pairs = 2*nct + 2*ct.M_cell + 1;

// --- compute current barrier parameter mu (adaptive, replaces fixed schedule) ---
auto compute_mu = [&]() {
    double xs = 0.0;
    for (int j = 0; j < nct; j++) {
        xs += s_up[j]*y_up[j] + s_dn[j]*y_dn[j];
    }
    for (int c = 0; c < ct.M_cell; c++) {
        xs += s_lo[c]*y_lo[c] + s_hi[c]*y_hi[c];
    }
    xs += s_delta * y_delta;
    return xs / static_cast<double>(std::max(m_pairs, 1));
};

double mu = compute_mu();

// --- Jacobi diagonal scaling buffers, sized once ---
std::vector<double> D_jac(nct_red);

for (int iter = 0; iter < max_ipm; iter++) {

    // Convergence check (existing residual-based criterion, unchanged) ---
    // if (residuals satisfy tol) break;

    // === DEGENERATE GUARD: m_pairs == 0 fallback ===
    if (m_pairs == 0) {
        // No complementarity pairs => degenerate LP. Fall back to plain barrier
        // step with mu *= 0.1 (existing single-step path). This branch is
        // unreachable for non-empty K (kept for safety).
        mu *= 0.1;
        continue;
    }

    // === Phase A — Affine predictor (sigma=0 implicit; barrier centering = 0) ===

    // 1. Build N_red = A_red * D * A_red^T using current s_*, y_*
    //    (existing build_N_red helper; produces nct_red x nct_red dense matrix
    //    in row-major std::vector<double> N(nct_red*nct_red).)
    std::vector<double> N(nct_red * nct_red, 0.0);
    build_N_red(/* in: s_up, s_dn, y_up, y_dn, s_lo, s_hi, y_lo, y_hi, ...,
                   out: */ N, nct_red, /* ... */);

    // 2. Compute Jacobi diagonal: D_jac[j] = 1 / sqrt(max(N[j,j], 1e-12))
    for (int j = 0; j < nct_red; j++) {
        const double djj = N[j*nct_red + j];
        D_jac[j] = 1.0 / std::sqrt(std::max(djj, 1e-12));
    }

    // 3. Scale N in-place: N_scaled[i,j] = D_jac[i] * N[i,j] * D_jac[j]
    for (int i = 0; i < nct_red; i++) {
        const double di = D_jac[i];
        for (int j = 0; j < nct_red; j++) {
            N[i*nct_red + j] *= di * D_jac[j];
        }
    }

    // 4. (verbose>=3 audit) dump cond estimate of raw vs scaled N
    if (st.verbose >= 3 && iter == 0) {
        // diag-only condition estimator: (max diag) / (min diag)
        double dmax = 0.0, dmin = 1e300;
        // We already destroyed the raw diagonal — recover via 1/D_jac^2.
        for (int j = 0; j < nct_red; j++) {
            const double raw_djj = 1.0 / (D_jac[j]*D_jac[j]);
            dmax = std::max(dmax, raw_djj);
            dmin = std::min(dmin, raw_djj);
        }
        char msg[160];
        std::snprintf(msg, sizeof(msg),
            "chebyshev[mehrotra]: raw_diag_ratio=%.3e (Jacobi scales to ~1)",
            dmax/std::max(dmin, 1e-300));
        st.log_fn(msg);
    }

    // 5. Build Phase-A RHS (affine, centering term = 0)
    //    rhs_A[m] = -(S[m] - T_flat[m]*W) + D_marg[m] * (-y_up[m] + y_dn[m])
    //                                      [centering rmu_*/s_* terms with sigma=0]
    std::vector<double> rhs_A_full(nct, 0.0);
    for (int m = 0; m < nct; m++) {
        // sigma = 0 in Phase A => rmu_up = -s_up*y_up; rmu_dn = -s_dn*y_dn
        // term -y_up + y_dn comes from rmu_up/s_up + rmu_dn/s_dn telescoping
        rhs_A_full[m] = -(S[m] - T_flat[m]*W)
                      + D_marg[m] * (-y_up[m] + y_dn[m]);
    }
    // Reduce to nct_red via existing reduction routine (skip is_ref entries)
    std::vector<double> rhs_A_red(nct_red);
    for (int r = 0; r < nct_red; r++) rhs_A_red[r] = rhs_A_full[red_to_full[r]];

    // 6. Apply Jacobi to RHS: rhs_A_scaled[j] = D_jac[j] * rhs_A_red[j]
    for (int j = 0; j < nct_red; j++) rhs_A_red[j] *= D_jac[j];

    // 7. Factorize scaled N once (one LDLT per outer iteration)
    LdltFactor F;
    if (ldlt_factor(N, nct_red, F) != 0) {
        // Fall back: retry with stronger Tikhonov regularization on N
        for (int j = 0; j < nct_red; j++) N[j*nct_red + j] += 1e-8;
        if (ldlt_factor(N, nct_red, F) != 0) {
            res.status  = RK_ERR_NOCONV;
            std::snprintf(res.message, sizeof(res.message),
                          "chebyshev: LDLT factorization failed (Jacobi+Tikhonov)");
            return res;
        }
    }
    res.n_factorizations++;

    // 8. Solve Phase A: F · Δλ_A_scaled = rhs_A_scaled
    std::vector<double> dlam_A_red(nct_red);
    ldlt_solve(F, rhs_A_red.data(), dlam_A_red.data(), nct_red);

    // 9. Unscale: Δλ_A_red[j] = D_jac[j] * Δλ_A_red_scaled[j]
    for (int j = 0; j < nct_red; j++) dlam_A_red[j] *= D_jac[j];

    // 10. Expand to full Δλ_A_full (zero on is_ref entries)
    std::vector<double> dlam_A_full(nct, 0.0);
    for (int r = 0; r < nct_red; r++) dlam_A_full[red_to_full[r]] = dlam_A_red[r];

    // 11. Reconstruct Δs_aff, Δy_aff, ΔX_aff, Δδ_aff via existing recovery code
    //     (Sherman-Morrison for δ-row, etc.). Use sigma=0 in centering terms.
    std::vector<double> dS_up_aff(nct), dS_dn_aff(nct);
    std::vector<double> dY_up_aff(nct), dY_dn_aff(nct);
    std::vector<double> dX_aff(ct.M_cell);
    std::vector<double> dY_lo_aff(ct.M_cell), dY_hi_aff(ct.M_cell);
    double dDelta_aff = 0.0, dY_delta_aff = 0.0;
    recover_directions(/* sigma=*/ 0.0, dlam_A_full,
                       dS_up_aff, dS_dn_aff, dY_up_aff, dY_dn_aff,
                       dX_aff, dY_lo_aff, dY_hi_aff,
                       dDelta_aff, dY_delta_aff);

    // 12. Affine step length α_aff (max α s.t. all primal/dual slacks ≥ 0)
    auto step_max = [&](double s, double ds) -> double {
        return (ds < 0.0) ? std::min(1.0, -s / ds) : 1.0;
    };
    double alpha_p_aff = 1.0, alpha_d_aff = 1.0;
    for (int j = 0; j < nct; j++) {
        alpha_p_aff = std::min(alpha_p_aff, step_max(s_up[j], dS_up_aff[j]));
        alpha_p_aff = std::min(alpha_p_aff, step_max(s_dn[j], dS_dn_aff[j]));
        alpha_d_aff = std::min(alpha_d_aff, step_max(y_up[j], dY_up_aff[j]));
        alpha_d_aff = std::min(alpha_d_aff, step_max(y_dn[j], dY_dn_aff[j]));
    }
    for (int c = 0; c < ct.M_cell; c++) {
        alpha_p_aff = std::min(alpha_p_aff, step_max(s_lo[c], dX_aff[c]));
        alpha_p_aff = std::min(alpha_p_aff, step_max(s_hi[c], -dX_aff[c]));
        alpha_d_aff = std::min(alpha_d_aff, step_max(y_lo[c], dY_lo_aff[c]));
        alpha_d_aff = std::min(alpha_d_aff, step_max(y_hi[c], dY_hi_aff[c]));
    }
    alpha_p_aff = std::min(alpha_p_aff, step_max(s_delta, dDelta_aff));
    alpha_d_aff = std::min(alpha_d_aff, step_max(y_delta, dY_delta_aff));
    const double alpha_aff = std::min(alpha_p_aff, alpha_d_aff);

    // 13. Affine mu: mu_aff = (x+α*Δx)·(s+α*Δs) / (2*m_pairs), clamped
    auto mu_aff_term = [&](double s, double ds, double y, double dy, double a) {
        return (s + a*ds) * (y + a*dy);
    };
    double mu_aff_num = 0.0;
    for (int j = 0; j < nct; j++) {
        mu_aff_num += mu_aff_term(s_up[j], dS_up_aff[j], y_up[j], dY_up_aff[j], alpha_aff);
        mu_aff_num += mu_aff_term(s_dn[j], dS_dn_aff[j], y_dn[j], dY_dn_aff[j], alpha_aff);
    }
    for (int c = 0; c < ct.M_cell; c++) {
        mu_aff_num += mu_aff_term(s_lo[c], dX_aff[c], y_lo[c], dY_lo_aff[c], alpha_aff);
        mu_aff_num += mu_aff_term(s_hi[c], -dX_aff[c], y_hi[c], dY_hi_aff[c], alpha_aff);
    }
    mu_aff_num += mu_aff_term(s_delta, dDelta_aff, y_delta, dY_delta_aff, alpha_aff);
    double mu_aff = mu_aff_num / static_cast<double>(2 * m_pairs);
    mu_aff = std::clamp(mu_aff, 0.0, mu * 100.0);   // guard

    // 14. Centering parameter sigma = clamp((mu_aff/mu)^3, 1e-8, 1.0)
    double sigma = (mu > 1e-300)
                   ? std::pow(mu_aff / mu, 3.0)
                   : 1.0;
    sigma = std::clamp(sigma, 1e-8, 1.0);

    // === Phase B — Corrector (centering + second-order) ===

    // 15. Build Phase-B RHS:
    //     rhs_B[m] = rhs_A[m] + (centering = sigma*mu / barrier-structured terms)
    //                         + second-order = -(Δs_aff * Δy_aff)/s contribution
    //     Detailed assembly mirrors Phase A but with rmu_* = sigma*mu - s*y
    //     AND subtracts Δs_aff*Δy_aff from the per-pair rmu term.
    std::vector<double> rhs_B_full(nct, 0.0);
    const double sigma_mu = sigma * mu;
    for (int m = 0; m < nct; m++) {
        const double rmu_up = sigma_mu - s_up[m]*y_up[m] - dS_up_aff[m]*dY_up_aff[m];
        const double rmu_dn = sigma_mu - s_dn[m]*y_dn[m] - dS_dn_aff[m]*dY_dn_aff[m];
        rhs_B_full[m] = -(S[m] - T_flat[m]*W)
                      + D_marg[m] * (rmu_up/s_up[m] - rmu_dn/s_dn[m]);
    }
    std::vector<double> rhs_B_red(nct_red);
    for (int r = 0; r < nct_red; r++) rhs_B_red[r] = rhs_B_full[red_to_full[r]];

    // 16. Scale RHS_B with same Jacobi diagonals (factor reused, NOT re-scaled)
    for (int j = 0; j < nct_red; j++) rhs_B_red[j] *= D_jac[j];

    // 17. Solve REUSING factor F (no refactor). Counter NOT incremented.
    std::vector<double> dlam_B_red(nct_red);
    ldlt_solve(F, rhs_B_red.data(), dlam_B_red.data(), nct_red);

    // 18. Unscale: Δλ_B_red[j] = D_jac[j] * Δλ_B_red_scaled[j]
    for (int j = 0; j < nct_red; j++) dlam_B_red[j] *= D_jac[j];

    // 19. Expand to full Δλ and recover corrected directions
    std::vector<double> dlam_B_full(nct, 0.0);
    for (int r = 0; r < nct_red; r++) dlam_B_full[red_to_full[r]] = dlam_B_red[r];
    std::vector<double> dS_up(nct), dS_dn(nct), dY_up(nct), dY_dn(nct);
    std::vector<double> dX(ct.M_cell), dY_lo(ct.M_cell), dY_hi(ct.M_cell);
    double dDelta = 0.0, dY_delta = 0.0;
    recover_directions(sigma, dlam_B_full,
                       dS_up, dS_dn, dY_up, dY_dn,
                       dX, dY_lo, dY_hi,
                       dDelta, dY_delta);
    // Add Mehrotra second-order correction on slacks (subtract dS_aff*dY_aff/s)
    for (int j = 0; j < nct; j++) {
        dY_up[j] -= dS_up_aff[j]*dY_up_aff[j] / s_up[j];
        dY_dn[j] -= dS_dn_aff[j]*dY_dn_aff[j] / s_dn[j];
    }

    // 20. Step lengths with 0.99 damping (not 1.0 — keep strict interior)
    double alpha_p = 1.0, alpha_d = 1.0;
    for (int j = 0; j < nct; j++) {
        alpha_p = std::min(alpha_p, step_max(s_up[j], dS_up[j]));
        alpha_p = std::min(alpha_p, step_max(s_dn[j], dS_dn[j]));
        alpha_d = std::min(alpha_d, step_max(y_up[j], dY_up[j]));
        alpha_d = std::min(alpha_d, step_max(y_dn[j], dY_dn[j]));
    }
    for (int c = 0; c < ct.M_cell; c++) {
        alpha_p = std::min(alpha_p, step_max(s_lo[c], dX[c]));
        alpha_p = std::min(alpha_p, step_max(s_hi[c], -dX[c]));
        alpha_d = std::min(alpha_d, step_max(y_lo[c], dY_lo[c]));
        alpha_d = std::min(alpha_d, step_max(y_hi[c], dY_hi[c]));
    }
    alpha_p = std::min(alpha_p, step_max(s_delta, dDelta));
    alpha_d = std::min(alpha_d, step_max(y_delta, dY_delta));
    constexpr double kDamp = 0.99;
    alpha_p *= kDamp;
    alpha_d *= kDamp;

    // 21. Update primal+dual+slack
    for (int c = 0; c < ct.M_cell; c++) X[c]    += alpha_p * dX[c];
    delta    += alpha_p * dDelta;
    for (int j = 0; j < nct; j++) {
        s_up[j] += alpha_p * dS_up[j]; s_dn[j] += alpha_p * dS_dn[j];
        y_up[j] += alpha_d * dY_up[j]; y_dn[j] += alpha_d * dY_dn[j];
    }
    for (int c = 0; c < ct.M_cell; c++) {
        s_lo[c] = std::max(X[c]-L_cell[c], kEps);
        s_hi[c] = std::max(U_cell[c]-X[c], kEps);
        y_lo[c] += alpha_d * dY_lo[c]; y_hi[c] += alpha_d * dY_hi[c];
    }
    s_delta = std::max(s_delta + alpha_p*dDelta, kEps);
    y_delta += alpha_d * dY_delta;

    // 22. Recompute mu adaptively (no schedule — Mehrotra-driven)
    mu = compute_mu();
}
// === end Mehrotra loop ========================================================
```

Notes:
- `recover_directions(sigma, dlam_full, ...)`: **does NOT exist as a named function** — it must be WRITTEN AS A NEW LAMBDA inside `chebyshev_ipm`, extracting the per-pair recovery logic from lines ~390-401 of chebyshev.cpp. Read those lines first, then write the lambda so Phase A (sigma=0) and Phase B (sigma>0) share the formula.
- `build_N_red`: **also does NOT exist as a named helper** — `N_red` is built inline at lines ~175-178. Extract this into a lambda if the Mehrotra loop needs to call it twice; otherwise, build it once and reuse the factored form across Phase A and Phase B. The key: LDLT-factor N once per outer iteration, solve twice (Phase A RHS, Phase B RHS). Do NOT refactor N between phases.
- `ldlt_factor`/`ldlt_solve` USE the existing `lbw::ldlt_factor_inplace` and `lbw::ldlt_solve` from `calib_linalg.hpp` — these DO exist. No new helpers needed for these.
- Before writing ANY Phase A/B code, read `src/chebyshev.cpp` lines 83-420 fully to understand the existing IPM structure. The reviewer confirmed `N_red` exists (lines ~175-178) but `build_N_red` as a callable does not.
- The schur_nu diagnostic at lines 322-356 still runs once at iter 0 against the RAW (pre-Jacobi) N — preserve it before step 3 of the new loop.

- [ ] 5.4 Set `res.n_factorizations = iter_count` at return so the audit counter is exposed in the C++ result. Pass through to R via `r_bridge.cpp` if a return-list field already aggregates diagnostics; otherwise add `n_factorizations` to the SEXP list returned (named "n_factorizations").

- [ ] 5.5 Compile gate:
```bash
R CMD INSTALL --preclean .
```
Expected: clean build. If LDLT helpers were inline-only, the lambda-lift will surface — fix and rebuild.

- [ ] 5.6 Run target tests:
```bash
Rscript -e 'devtools::test(filter = "chebyshev")'
```
Expected:
- `T_greg_warn` → PASS (untouched).
- `T_cheby_warm_fallback` → PASS (warm-start finite + Mehrotra converges).
- `T_cheby_warm` → PASS (status=0, max_err ≤ raking*1.05).
- `T_cheby_K9` → PASS if fixture present, else SKIPPED.

- [ ] 5.7 Run full test suite:
```bash
Rscript -e 'devtools::test()'
```
Expected: no regressions across ieppa, raking, sinkhorn, lbfgsb, harvest, calibration-solvers.

- [ ] 5.8 Verification before completion. Confirm via spy logs (verbose=3) on the K=3 fixture:
```bash
Rscript -e '
fx <- {
  set.seed(42L); n <- 200L
  c1 <- factor(sample.int(4L,n,replace=TRUE))
  c2 <- factor(sample.int(3L,n,replace=TRUE))
  c3 <- factor(sample.int(5L,n,replace=TRUE))
  list(design=data.frame(cat1=c1,cat2=c2,cat3=c3), d=rep(1.0,n),
       target=list(cat1=as.numeric(table(c1))*1.20,
                   cat2=as.numeric(table(c2))*0.85,
                   cat3=as.numeric(table(c3))*1.10))
}
r <- leafblower::harvest(fx$design, fx$d, fx$target,
                         method="chebyshev", max_weight=5,
                         max_iterations=200L, verbose=3L)
cat("status=",r$status," iter=",r$iterations,
    " max_err=",r$max_error," n_factorizations=",r$n_factorizations,"\n")
stopifnot(r$status==0L, r$n_factorizations==r$iterations)
'
```
Expected:
- `status=0`.
- `iter` (Mehrotra) significantly fewer than the old 500-cap (target: ≤ 50 for K=3 tight).
- `n_factorizations == iter` (proves one factor per outer iter, two solves).
- Verbose log contains "chebyshev[mehrotra]: raw_diag_ratio=…".

- [ ] 5.9 Commit:
```bash
git add src/chebyshev.cpp src/chebyshev.hpp src/r_bridge.cpp
git commit -m "feat(chebyshev): Mehrotra predictor-corrector + Jacobi preconditioning

- Replace fixed mu *= 0.1 schedule with adaptive Mehrotra two-phase IPM.
- Phase A (affine, sigma=0) builds N=A·D·A^T, applies Jacobi diagonal scaling,
  factors LDLT once, solves, unscales -> Δλ_A.
- Compute alpha_aff, mu_aff = clamp(.,0,mu*100), sigma = clamp((mu_aff/mu)^3, 1e-8, 1).
- Phase B (corrector) reuses the same LDLT factor; RHS adds centering and
  second-order Mehrotra term -Δs_aff*Δy_aff. One factor, two solves per iter.
- m_pairs==0 guard falls back to plain barrier (degenerate LP safety).
- 0.99 step damping keeps strict interior.
- Expose n_factorizations diagnostic counter for audit (==iterations invariant).

Fixes K=9 stepstone NaN convergence. Used together with ieppa warm-start
(Tasks 3-4), chebyshev now converges in O(20-50) iter on K=9 vs. the prior
500-iter hard cap that produced NaN."
```

**Done when:** All chebyshev tests GREEN; full suite GREEN; `n_factorizations == iterations` audit invariant holds; Mehrotra iter count << 500.

---

## Final Verification Checklist (before declaring plan COMPLETE)

- [ ] All five Tasks committed individually (5 commits, no squashing).
- [ ] `R CMD INSTALL --preclean .` clean from a fresh checkout of HEAD.
- [ ] `Rscript -e 'devtools::test()'` returns 0 with no new failures vs. the baseline before Task 1.
- [ ] `tests/testthat/test-chebyshev.R` shows 3 PASS + 1 SKIP (or 4 PASS if stepstone fixture present).
- [ ] `n_factorizations == iterations` invariant holds in the verbose=3 spy log.
- [ ] verbose>=2 log on chebyshev path emits the warm-start total_pre/post line.
- [ ] greg on tight-bounds emits "greg.*unreliable" warning AND still returns weights (warn-and-continue).
- [ ] `git log --oneline` shows 5 conventional-commit-style commits, no AI attribution / emojis / Co-Authored-By trailers.
- [ ] Anatomy / memory updated: append entries to `.wolf/memory.md` and `.wolf/anatomy.md` for `tests/testthat/test-chebyshev.R` (new file) + signature change in `src/chebyshev.hpp`.

PLAN_COMPLETE
