# Sinkhorn A1 Fix + ieppa_soft Method Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix `convergence_used$solver_objective` to report each solver's mathematical objective (not stopping criterion value). Add `method="ieppa_soft"` (ADMM capacity). Revert `method="ieppa"` to hard clamp.

**Architecture:**
- Task 1 (leafblower-4ijo): RED tests — must ERROR before any code change.
- Task 2 (leafblower-3mq8): Decouple solver_objective — `calib_dispatch.hpp` + all structs + convergence blocks + field rename + sinkhorn default.
- Task 3 (leafblower-i3tj): `method="ieppa_soft"` — new RK_ALG constant + `use_admm_capacity` flag + P1.1 gate + ieppa hard-clamp revert.
- Task 4 (leafblower-zee9): Fixture regeneration + A1 test update.
- Task 5 (leafblower-lig9): Final benchmark + verification.

**Tech Stack:** C++17, R. Build: `R CMD INSTALL --preclean . 2>&1 | tail -3`.

**Key invariants:**
- `rk_calib_result_t` field rename: `convergence_objective` → `convergence_solver_objective` (same position/type, no ABI layout change).
- `RK_ALG_IEPPA_SOFT = 8` — new constant added to `leafblower.h` enum (no struct size change).
- ieppa weight KL = `Σ_c X[c] * log(X[c] / X_init[c]) / n` (same as `m.kl` from compute_cell_metrics but computed inline).

---

## File Map

| File | Tasks | Change |
|------|-------|--------|
| `tests/testthat/test-calibration-solvers.R` | 1, 4 | RED tests T1/T2/T3; update A1 |
| `src/calib_dispatch.hpp` | 2 | Add `select_solver_objective(alg_id, m)` |
| `src/ieppa.hpp` | 2 | Add `best_objective_seen`; rename `convergence_objective` → `convergence_solver_objective` |
| `src/sinkhorn.hpp` (or `.cpp`) | 2 | Same |
| `src/greg.hpp` / `src/grake.hpp` | 2 | Same |
| `src/leafblower.h` | 2, 3 | Rename C API field; add `RK_ALG_IEPPA_SOFT = 8` |
| `src/ieppa.cpp` | 2, 3 | Atomic BEST-ITER UPDATE blocks; gate ADMM; revert default |
| `src/sinkhorn.cpp` | 2 | Atomic BEST-ITER UPDATE; `convergence_solver_objective = best_objective_seen` |
| `src/raking.cpp` | 2 | Same (raking minimizes weight KL — missing from original plan, added) |
| `src/greg.cpp` | 2 | Same |
| `src/grake.cpp` | 2 | Same |
| `src/c_api.cpp` | 2, 3 | Field rename; `ieppa_soft` dispatch + `use_admm_capacity` |
| `src/r_bridge.cpp` | 2, 3 | Field rename exposure; pass flag |
| `src/types.hpp` | 3 | Add `bool use_admm_capacity = false` to CalibState |
| `R/harvest.R` | 2, 3 | Field rename; sinkhorn default `metric="kl"`; add `"ieppa_soft"` dispatch |
| `data-raw/gen_ieppa_kl_ref.R` | 4 | Use `$solver_objective` not `$best_error` |

---

## Task 1: RED tests T1/T2/T3 (leafblower-4ijo)

**Files:** `tests/testthat/test-calibration-solvers.R` (append)

Expected RED mode: T1 = ERROR (field `$solver_objective` not found). T2/T3 = ERROR (`method="ieppa_soft"` unknown method). ERROR = valid RED.

- [ ] **Step 1.1: Append three tests**

```r
# ── T1: solver_objective field existence (RED: field not found before Task 2) ──
test_that("T1: sinkhorn convergence_used$solver_objective exists and is weight KL", {
  set.seed(1L); n <- 800L
  df  <- data.frame(v1 = factor(sample(c("A","B","C"), n, TRUE)))
  tgt <- list(v1 = c("A"=0.6, "B"=0.3, "C"=0.1))
  r_mx <- leafblower::harvest(df, tgt, method="sinkhorn",
    convergence=list(metric="max_err"), max_iterations=200L, attach_weights=FALSE)
  r_kl <- leafblower::harvest(df, tgt, method="sinkhorn",
    convergence=list(metric="kl"),     max_iterations=200L, attach_weights=FALSE)
  obj_mx <- attr(r_mx, "result")$convergence_used$solver_objective
  obj_kl <- attr(r_kl, "result")$convergence_used$solver_objective
  # Both stopping criteria → same mathematical objective (weight KL, not stopping value)
  expect_false(is.null(obj_mx),    label="solver_objective field exists")
  expect_true(is.finite(obj_mx),  label="solver_objective is finite")
  expect_true(obj_mx < 0.5,       label="solver_objective is weight KL (not max_err ~0.05)")
  expect_true(abs(obj_mx - obj_kl) / max(obj_mx, obj_kl) < 0.5,
              label="objective consistent across stopping criteria")
})

# ── T2: ieppa_soft available (RED: unknown method before Task 3) ──
test_that("T2: ieppa_soft method exists and respects max_weight", {
  set.seed(2L); n <- 2000L
  df  <- data.frame(v1 = factor(sample(c("X","Y"), n, TRUE, prob=c(.3,.7))))
  tgt <- list(v1 = c("X"=0.8, "Y"=0.2))
  r <- leafblower::harvest(df, tgt, method="ieppa_soft",
    max_weight=2.0, min_weight=0.0, max_iterations=300L, attach_weights=FALSE)
  w <- as.numeric(r)
  expect_true(max(w) <= 2.0 + 1e-9, label="ieppa_soft wmax ≤ max_weight")
  expect_true(min(w) >= 0.0 - 1e-9, label="ieppa_soft wmin ≥ min_weight")
  expect_equal(attr(r, "result")$status, 0L)
})

# ── T3: ieppa_soft strictly better than ieppa on tight-bounds problem ──
test_that("T3: ieppa_soft max_err strictly < ieppa max_err on tight bounds", {
  set.seed(3L); n <- 5000L
  df  <- data.frame(v1 = factor(sample(5L, n, TRUE)))
  tgt <- list(v1 = setNames(c(0.4, 0.3, 0.15, 0.1, 0.05), as.character(1:5)))
  r_hard <- leafblower::harvest(df, tgt, method="ieppa",
    max_weight=1.8, min_weight=0, max_iterations=500L, attach_weights=FALSE)
  r_soft <- leafblower::harvest(df, tgt, method="ieppa_soft",
    max_weight=1.8, min_weight=0, max_iterations=500L, attach_weights=FALSE)
  me_hard <- attr(r_hard, "result")$max_error
  me_soft <- attr(r_soft, "result")$max_error
  expect_true(me_soft < me_hard,
              label="ieppa_soft max_err strictly better on tight bounds")
})
```

- [ ] **Step 1.2: Verify all three ERROR (not PASS)**
```bash
cd /home/dd/Gemini/leafblower
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "T1|T2|T3|Error|error|FAIL|PASS" | head -10
```
Expected: T1/T2/T3 show ERROR (not PASS). Any PASS is a false positive.

- [ ] **Step 1.3: Commit RED tests**
```bash
git add tests/testthat/test-calibration-solvers.R
git commit -m "test(sinkhorn): RED — T1 solver_objective field, T2/T3 ieppa_soft method

T1: ERROR before Task 2 (field not found).
T2/T3: ERROR before Task 3 (unknown method).

Closes: leafblower-4ijo"
```

---

## Task 2: Decouple solver_objective + rename + sinkhorn default (leafblower-3mq8)

**Files:** `src/calib_dispatch.hpp`, `src/ieppa.hpp`, `src/leafblower.h`, `src/ieppa.cpp`, `src/sinkhorn.cpp`, `src/raking.cpp`, `src/greg.cpp`, `src/grake.cpp`, `src/c_api.cpp`, `src/r_bridge.cpp`, `R/harvest.R`

### Step 2.1: Add `select_solver_objective` to calib_dispatch.hpp

**Critical distinction:** `m.kl` from `compute_cell_metrics` is **marginal KL** (max over margins
of Σ_j T_kj log(T_kj/achieved_kj)), NOT weight-space KL. Sinkhorn and ieppa minimize
weight-space KL = Σ_i w_i log(w_i/d_i). These are different quantities — using m.kl for
solver_objective of KL-minimizing solvers would make A1 compare apples to oranges.

For KL-minimizing solvers (ieppa, sinkhorn, raking), weight-space KL is computed inline
via a `compute_weight_kl` lambda (see Step 2.4). `select_solver_objective` handles only
the non-KL solvers where m.chi2/m.grake_norm/m.errRp are the correct objectives:

Read `src/calib_dispatch.hpp`. After the existing `select_metric` function, add:
```cpp
// Solver mathematical objective for NON-KL solvers only.
// KL-minimizing solvers (ieppa, sinkhorn, raking) compute weight-space KL inline
// via compute_weight_kl lambda — DO NOT use this for those solvers.
// chebyshev: objective == stopping metric by design (documented).
inline double select_solver_objective(int alg_id, const lbw::CellMetrics& m) {
    switch (alg_id) {
    case RK_ALG_GREG:      return m.chi2;
    case RK_ALG_GRAKE:     return m.grake_norm;
    case RK_ALG_CHEBYSHEV: return m.errRp;   // objective == stopping metric
    default:               return m.errRp;
    }
}
```
(RK_ALG_IEPPA_SOFT is added in Task 3 — it is a KL solver, so not in this switch.)

### Step 2.2: Rename field in `src/leafblower.h`

Find `double convergence_objective;` (line ~123). Replace with:
```cpp
double convergence_solver_objective;  /* solver's mathematical objective at best_iter */
```
Also find `EXPECTED_RK_RESULT_BYTES`. Verify size unchanged (same type, same position — no change needed). If the assert fires, investigate before proceeding.

### Step 2.3: Rename field in solver result structs

In `src/ieppa.hpp`: find `double convergence_objective = 0.0;` (line ~40). Replace with:
```cpp
double best_objective_seen         = 0.0;   // internal: objective at best_iter (ieppa weight KL)
double convergence_solver_objective = 0.0;  // exposed: same as best_objective_seen
```

In `src/sinkhorn.hpp` (or at the top of `sinkhorn.cpp`): find the SinkhornResult struct. Add same two fields.

For greg/grake: read `src/greg.hpp` and `src/grake.hpp`, add same two fields to each result struct.

### Step 2.4: Update ieppa.cpp convergence blocks — atomic BEST-ITER UPDATE

Read lines 899-918 of `src/ieppa.cpp` (BLOCK 1 and BLOCK 1b).

For each best-iter update block, add `best_objective_seen` computation. ieppa's mathematical objective = weight KL = `Σ_c X[c] * log(X[c] / X_init[c]) / n`.

The three update sites (BLOCK 1 at ~902, BLOCK 1b at ~912, need_extra_metrics at ~1035) must each update `best_objective_seen`. Use this helper lambda (add once before the convergence block):

```cpp
// Weight KL objective for ieppa — computed independently of stopping criterion.
// Must remain co-located with best_iter update (atomic block invariant).
auto compute_weight_kl = [&]() -> double {
    double wkl = 0.0;
    const double inv_n = 1.0 / static_cast<double>(st.n);
    for (int c = 0; c < ct.M_cell; c++) {
        if (X_init[c] > 0.0 && X[c] > 0.0)
            wkl += X[c] * std::log(X[c] / X_init[c]) * inv_n;
    }
    return std::isfinite(wkl) ? wkl : 0.0;
};
```

Update BLOCK 1 (MAX_ERR, ~line 901):
```cpp
            if (st.convergence_cfg.metric == lbw::CalibMetric::MAX_ERR) {
                if (errRp < best_metric_seen) {
                    // === BEST-ITER UPDATE (metric, iter, objective MUST stay co-located) ===
                    best_metric_seen    = errRp;
                    best_iter_val       = iter;
                    best_objective_seen = compute_weight_kl();
                    for (int c = 0; c < ct.M_cell; c++)
                        W_best[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
                    // === END BEST-ITER UPDATE ===
                }
            }
```

Update BLOCK 1b (MARGINAL_KL, ~line 911): same pattern.

Update the `need_extra_metrics` block (around line 1035): replace
```cpp
if (std::isfinite(curr_best) && curr_best < best_metric_seen) {
    best_metric_seen = curr_best;
    best_iter_val    = iter;
    ...
}
```
with the same atomic block pattern, adding `best_objective_seen = compute_weight_kl();` inside.

At line ~1186 (`res.convergence_objective = best_metric_seen;`), change to:
```cpp
res.convergence_solver_objective = best_objective_seen;
```

Also at line ~1193 (`res.best_error = best_metric_seen;`) — leave unchanged (best_error tracks the stopping metric, correct).

### Step 2.5: Update sinkhorn.cpp (src/sinkhorn.hpp)

Read lines 163-185 of `src/sinkhorn.hpp`. Add `compute_weight_kl` lambda before the
error-check loop (same lambda as in ieppa.cpp, Step 2.4). Declare
`double best_objective_seen = 0.0;` near `best_metric_seen` (~line 81).

The best-iter update pattern:
```cpp
            if (std::isfinite(curr_best) && curr_best < best_metric_seen) {
                // === BEST-ITER UPDATE ===
                best_metric_seen    = curr_best;
                best_iter_val       = iter;
                best_objective_seen = compute_weight_kl();  // weight-space KL, not m.kl
                W_best              = X;
                // === END BEST-ITER UPDATE ===
            }
```

Also declare `double best_objective_seen = 0.0;` near `best_metric_seen` declaration (line ~81).

At line 183 (`res.convergence_objective = best_metric_seen;`), change to:
```cpp
res.convergence_solver_objective = best_objective_seen;
```

### Step 2.6: Update raking.cpp, greg.cpp, and grake.cpp

**raking.cpp** also has `res.convergence_objective = best_metric_seen` (line 330) and
best_iter tracking. Raking minimizes weight KL — add `compute_weight_kl` lambda and
atomic BEST-ITER UPDATE blocks, same pattern as sinkhorn (Steps 2.4-2.5).
Change `res.convergence_objective = best_metric_seen;` → `res.convergence_solver_objective = best_objective_seen;`.

For **greg and grake**, the mathematical objective IS in `m` (chi² and grake_norm respectively),
so `select_solver_objective(alg_id, m)` can be used:
```cpp
best_objective_seen = select_solver_objective(RK_ALG_GREG, m);   // m.chi2
// or
best_objective_seen = select_solver_objective(RK_ALG_GRAKE, m);  // m.grake_norm
```

Declare `double best_objective_seen = 0.0;` and add the atomic BEST-ITER UPDATE block
in each solver's best_iter update location. Change
`res.convergence_objective = best_metric_seen;` to
`res.convergence_solver_objective = best_objective_seen;` in each.

### Step 2.7: Update c_api.cpp

Find all occurrences of `result->convergence_objective` (or `res.convergence_objective`). Replace with `result->convergence_solver_objective` (or `res.convergence_solver_objective`).

### Step 2.8: Update r_bridge.cpp

Find the line that exposes `convergence_objective` (search: `"convergence_objective"` — this
is the literal string in the R element name, at line ~575):
```cpp
SET_STRING_ELT(res_names, 28, Rf_mkChar("convergence_objective"));
```
Replace with:
```cpp
SET_STRING_ELT(res_names, 28, Rf_mkChar("solver_objective"));
```
Update the corresponding `SET_VECTOR_ELT` line to use `res.convergence_solver_objective`.

### Step 2.9: Update R/harvest.R — sinkhorn default + field reference

Add sinkhorn default metric:
```r
  if (method == "sinkhorn" &&
      is.null(convergence[["metric"]]) &&
      is.null(convergence[["criterion"]]) &&
      is.null(convergence[["improvement"]]) &&
      is.null(convergence[["pct"]]) &&
      is.null(convergence[["absolute"]])) {
    conv$metric <- "kl"
  }
```

Search for any R code accessing `result$convergence_used$objective` and update to `$solver_objective`.

### Step 2.10: Compile gate
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
```

### Step 2.11: Verify T1 goes GREEN, T2/T3 still ERROR
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "T1|T2|T3|FAIL|PASS|Error" | head -10
```
Expected: T1 PASS, T2/T3 ERROR.

### Step 2.12: Full test suite
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 2 (pre-existing), PASS ≥ 381.

### Step 2.13: Commit Task 2
```bash
git add src/calib_dispatch.hpp src/ieppa.hpp src/sinkhorn.hpp src/greg.hpp src/grake.hpp \
        src/leafblower.h src/ieppa.cpp src/sinkhorn.cpp src/raking.cpp src/greg.cpp src/grake.cpp \
        src/c_api.cpp src/r_bridge.cpp R/harvest.R
git commit -m "feat: decouple solver_objective from stopping criterion

Add select_solver_objective(alg_id, m) to calib_dispatch.hpp.
Track best_objective_seen in all solver convergence blocks (atomic).
Rename convergence_objective -> convergence_solver_objective everywhere.
sinkhorn default metric: max_err -> kl (its mathematical objective).
ieppa weight KL computed inline via compute_weight_kl() lambda.

Closes: leafblower-3mq8"
```

---

## Task 3: ieppa_soft method + use_admm_capacity gate (leafblower-i3tj)

**Files:** `src/leafblower.h`, `src/types.hpp`, `src/ieppa.cpp`, `src/c_api.cpp`, `src/r_bridge.cpp`, `R/harvest.R`

### Step 3.1: Add RK_ALG_IEPPA_SOFT to leafblower.h

Find the `rk_algorithm_t` enum (lines ~39-47). Add after `RK_ALG_GRAKE = 7`:
```c
RK_ALG_IEPPA_SOFT = 8    /* ieppa + ADMM soft capacity enforcement */
```

Do NOT add RK_ALG_IEPPA_SOFT to `select_solver_objective` — ieppa_soft is a KL-minimizing
solver and uses `compute_weight_kl` inline (same as ieppa), not m.kl which is marginal KL.
The only change to calib_dispatch.hpp in Task 3 is the RK_ALG_IEPPA_SOFT constant itself.

### Step 3.2: Add `use_admm_capacity` to CalibState

Read `src/types.hpp` struct CalibState. Add:
```cpp
bool use_admm_capacity = false;   // ieppa_soft: ADMM P1.1; default false = hard clamp
```

### Step 3.3: Gate ADMM P1.1 in ieppa.cpp + revert default + conditional u[] allocation

**Change u[] declaration** (line ~151). Current:
```cpp
    std::vector<double> u(ct.M_cell, 0.0);
```
Change to:
```cpp
    std::vector<double> u;  // allocated only for ieppa_soft (ADMM)
    if (st.use_admm_capacity) u.assign(ct.M_cell, 0.0);
```

**Gate ADMM in P1.1 block** (lines ~777-784). Current:
```cpp
                double z = std::clamp(X_tilde_c + u[c], L_cell[c], U_cell[c]);
                u[c] += X_tilde_c - z;
                X[c] = z; W[c] = z / X_tilde_c; X_cur[c] = z;
```
Change to:
```cpp
                if (st.use_admm_capacity) {
                    double z = std::clamp(X_tilde_c + u[c], L_cell[c], U_cell[c]);
                    u[c] += X_tilde_c - z;
                    X[c] = z; W[c] = z / X_tilde_c; X_cur[c] = z;
                } else {
                    double xc = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);
                    X[c] = xc; W[c] = xc / X_tilde_c; X_cur[c] = xc;
                }
```

**Gate u resets in fallback blocks** (lines ~666 and ~682). Current:
```cpp
                std::fill(u.begin(), u.end(), 0.0);
```
Change each to:
```cpp
                if (st.use_admm_capacity) std::fill(u.begin(), u.end(), 0.0);
```

### Step 3.4: Add ieppa_soft dispatch in c_api.cpp

Find the algorithm dispatch switch (search for `RK_ALG_IEPPA:` case). Add:
```cpp
case RK_ALG_IEPPA_SOFT:
    st.use_admm_capacity = true;
    // fall through to ieppa_solve
    // (RK_ALG_IEPPA_SOFT uses the same solver with ADMM capacity enabled)
    ieppa_solve(st, res);
    result->algorithm_used = RK_ALG_IEPPA_SOFT;
    break;
```
(Adapt to the exact dispatch pattern in c_api.cpp — do NOT fall through silently; add explicit call.)

Also add `"ieppa_soft"` to the algorithm string mapping where `"ieppa"` is mapped to `RK_ALG_IEPPA`.

### Step 3.5: Update r_bridge.cpp

Add `"ieppa_soft"` to the method string→int mapping alongside `"ieppa"`.

### Step 3.6: Update R/harvest.R

Add `"ieppa_soft"` to the `match.arg(method, c("auto", "ieppa", ...))` vector and to the internal method dispatch. No other changes — ieppa_soft uses the same harvest.R path as ieppa, routed to the new C algorithm constant.

### Step 3.7: Compile gate
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
```

### Step 3.8: Verify T2 and T3 go GREEN
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "T1|T2|T3|FAIL|PASS" | head -5
```
Expected: T1, T2, T3 all PASS.

### Step 3.9: Full test suite
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```

### Step 3.10: Quick benchmark sanity check
```bash
OMP_NUM_THREADS=1 Rscript -e '
suppressPackageStartupMessages({library(arrow); library(leafblower)})
df <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
df$uuid <- NULL
tgt <- lapply(jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json"),
  function(t) { t <- unlist(t); t/sum(t) })
for(nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
for (m in c("ieppa","ieppa_soft")) {
  t0 <- proc.time()["elapsed"]
  r <- suppressWarnings(leafblower::harvest(df, tgt, method=m,
    max_weight=5, min_weight=0, max_iterations=500L, attach_weights=FALSE))
  wall <- proc.time()["elapsed"] - t0
  res <- attr(r,"result")
  cat(sprintf("%-12s wall=%.1fs iters=%d max_err=%.4e solver_obj=%.4e\n",
    m, wall, res$iterations, res$max_error,
    res$convergence_used$solver_objective))
}
' 2>&1 | grep -E "ieppa"
```
Expected:
- `ieppa`: max_err ≈ 2.74e-3 (hard clamp, pre-ADMM behavior restored)
- `ieppa_soft`: max_err ≈ 2.34e-3 (ADMM, better)

### Step 3.11: Commit Task 3
```bash
git add src/leafblower.h src/types.hpp src/ieppa.cpp src/calib_dispatch.hpp \
        src/c_api.cpp src/r_bridge.cpp R/harvest.R
git commit -m "feat(ieppa_soft): ADMM method + revert ieppa to hard clamp

Add RK_ALG_IEPPA_SOFT=8. Add CalibState::use_admm_capacity (default false).
Gate P1.1 ADMM on flag; allocate u[] only for ieppa_soft.
method='ieppa' reverts to hard clamp (max_err ~2.74e-3 on stepstone).
method='ieppa_soft' uses ADMM (max_err ~2.34e-3, +14% improvement).

Closes: leafblower-i3tj"
```

---

## Task 4: Regenerate A1 fixture + update A1 test (leafblower-zee9)

**PREREQUISITE:** Tasks 2 and 3 must be compiled and tested before this step.

### Step 4.1: Audit gen_ieppa_kl_ref.R

Read `data-raw/gen_ieppa_kl_ref.R`. Find where it extracts the ieppa KL value. Current (wrong):
```r
kl_at_best_iter = res$best_error  # wrong: best_error = max_err, not weight KL
```
Change to:
```r
kl_at_best_iter = res$convergence_used$solver_objective  # correct: weight KL at best_iter
```
Also update `method` to `"ieppa"` (hard clamp, not ADMM — the A1 baseline should be the hard-clamp ieppa since that is the standard version).

### Step 4.2: Run regeneration script

**MANDATORY:** Must compile (Task 3 complete) BEFORE running this script.
```bash
cd /home/dd/Gemini/leafblower && OMP_NUM_THREADS=1 Rscript data-raw/gen_ieppa_kl_ref.R 2>&1 | tail -5
```

### Step 4.3: Verify pre-condition: sinkhorn weight KL < new ieppa weight KL

```bash
OMP_NUM_THREADS=1 Rscript -e '
suppressPackageStartupMessages({library(arrow); library(jsonlite); library(leafblower)})
ref <- readRDS("tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds")
cat("New ieppa weight KL at best_iter:", ref$kl_at_best_iter, "\n")
df <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
df$uuid <- NULL
tgt <- lapply(fromJSON("benchmarks/stepstone_fulldata_bench_targets.json"),
  function(t) { t <- unlist(t); t/sum(t) })
for(nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
r_s <- suppressWarnings(leafblower::harvest(df, tgt, method="sinkhorn",
  max_weight=5, min_weight=0, max_iterations=5000L, attach_weights=FALSE))
kl_s <- attr(r_s,"result")$convergence_used$solver_objective
cat("Sinkhorn weight KL at convergence:", kl_s, "\n")
cat("A1 satisfied (kl_s < kl_i)?", kl_s < ref$kl_at_best_iter, "\n")
' 2>&1 | grep -E "ieppa|Sinkhorn|A1"
```

If A1 is NOT satisfied (sinkhorn weight KL > ieppa weight KL), do NOT commit — file a follow-up ticket for deeper sinkhorn algorithmic investigation.

### Step 4.4: Update A1 test to use `$solver_objective`

Read `tests/testthat/test-calibration-solvers.R` lines 39-66 (the A1 test). Find:
```r
  expect_lte(r_s$convergence_used$objective, ref$kl_at_best_iter, ...)
```
Change to:
```r
  expect_lte(attr(w_s,"result")$convergence_used$solver_objective, ref$kl_at_best_iter,
             label="sinkhorn weight KL ≤ ieppa weight KL at best_iter")
```

### Step 4.5: Run A1 test
```bash
Rscript -e '
  devtools::test_active_file("tests/testthat/test-calibration-solvers.R")
' 2>&1 | grep -E "A1|T1|T2|T3|PASS|FAIL" | head -5
```
Expected: A1 PASS, T1 PASS, T2 PASS, T3 PASS.

### Step 4.6: Full suite
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```

### Step 4.7: Commit fixture + script together
```bash
git add data-raw/gen_ieppa_kl_ref.R tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds \
        tests/testthat/test-calibration-solvers.R
git commit -m "fix: regenerate A1 fixture with correct weight KL + update A1 test

gen_ieppa_kl_ref.R now extracts \$solver_objective (weight KL at best_iter)
not \$best_error (max_err). A1 test updated to use \$solver_objective.
Fixture regenerated with method='ieppa' (hard clamp, standard baseline).
Pre-condition verified: kl_sinkhorn < kl_ieppa.

Closes: leafblower-zee9 leafblower-lig9"
```

---

## Task 5: Final benchmark + verification (leafblower-lig9)

- [ ] **Step 5.1: Run full method comparison**
```bash
OMP_NUM_THREADS=1 Rscript benchmarks/stepstone_all_methods.R 2>&1 | grep -E "^(===|ieppa|sinkhorn|autumn)"
```

- [ ] **Step 5.2: Verify acceptance criteria**

| AC | Expected | Actual |
|----|----------|--------|
| AC1: sinkhorn solver_objective ≠ max_err | obj ≠ max_err | |
| AC2: greg solver_objective = chi² | chi² value | |
| AC3: ieppa_soft available | no error | |
| AC4: ieppa_soft wmax ≤ max_weight | ≤ 5.000 | |
| AC5: ieppa max_err ≈ 2.74e-3 | ~2.74e-3 | |
| AC6: ieppa_soft max_err ≈ 2.34e-3 | ~2.34e-3 | |
| AC7: A1 PASS | PASS | |
| AC8: FAIL 2 pre-existing | FAIL 2 | |

- [ ] **Step 5.3: Python audit**
```bash
grep -r "convergence_used\|convergence_objective\|solver_objective" python/ 2>/dev/null | grep -v ".pyc" | head -10
```
If any assertions found, create separate ticket for Python update.

- [ ] **Step 5.4: Close tickets**
```bash
bd close leafblower-lig9 --reason="A1 fixed, ieppa_soft added, solver_objective decoupled. All ACs verified."
```

---

## Self-Review

**Spec coverage:**

| Spec item | Task |
|-----------|------|
| `select_solver_objective(alg_id, m)` with RK_ALG_* | 2.1 |
| `best_objective_seen` atomic update in all solvers | 2.4-2.6 |
| Field rename `convergence_solver_objective` everywhere | 2.2/2.3/2.7/2.8/2.9 |
| Sinkhorn default `metric="kl"` | 2.9 |
| `RK_ALG_IEPPA_SOFT = 8` | 3.1 |
| `CalibState::use_admm_capacity = false` | 3.2 |
| ADMM P1.1 gated on flag | 3.3 |
| `u[]` conditional allocation | 3.3 |
| `method="ieppa_soft"` dispatch | 3.4-3.6 |
| `method="ieppa"` reverts to hard clamp | 3.3 |
| gen_ieppa_kl_ref.R uses `$solver_objective` | 4.1 |
| Fixture regenerated AFTER Part 1+2 compiled | 4.2 (order explicit) |
| A1 pre-condition check before commit | 4.3 |
| T1 GREEN after Task 2, T2/T3 GREEN after Task 3 | 2.11, 3.8 |
| Python audit | 5.3 |

**Placeholder scan:** None. All code blocks complete.

**Type consistency:** `best_objective_seen` as `double`. `use_admm_capacity` as `bool`. `RK_ALG_IEPPA_SOFT` as int constant. All consistent throughout.
