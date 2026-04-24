# Convergence Reform + SOR + Best-Iterate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the design approved in `docs/superpowers/specs/2026-04-24-convergence-metrics-sor-design.md` (rev3, commit `e69fec3`): replace the default convergence criterion with `pct = 0.001`, add pluggable alternative criteria (max_err / mean_err / kl / chi2), report all 5 quality metrics at exit, add SOR adaptive under-relaxation (iEPPA only), and add best-iterate tracking across all three solvers.

**Architecture:** Six atomic work units (WU-A through WU-F) matching the spec §9 decomposition. Each WU is one beads ticket, one atomic commit. Scaffolding first (no behavioral change), then convergence criteria, then quality metrics, then SOR, then best-iterate, then Python parity + docs. Build gates (`R CMD INSTALL --preclean .`) after every C++ edit. TDD: failing test first for every WU.

**Tech Stack:** C++17 solvers (iEPPA/raking/lbfgsb), R wrapper (`R/harvest.R`), C ABI (`src/leafblower.h`), Python bindings, testthat 3, pytest.

**Baseline:** commit `e69fec3` on branch `feat/ieppa-homotopy-greenkhorn`. `devtools::test()` reports FAIL 0 | PASS 232 | SKIP 2.

**R API ground truth** (from `R/harvest.R:37–62`):
- `target` (singular), `method`, `max_iterations`, `convergence = list(...)`, `bounds_mode`.
- Overlay args from prior plan already wired: `homotopy_levels`, `homotopy_start_factor`, `homotopy_end_factor`, `homotopy_budget_p`, `scheduler`, `eta_schedule`, `eta_start`, `eta_end`, `eta_schedule_power`.

**Memory rules (non-negotiable):**
- Named target vectors only (bridge crashes on unnamed — named like `c("1"=0.1,"2"=0.9)`).
- `harvest()` defaults to `attach_weights=TRUE` (data.frame return). Tests use `attach_weights=FALSE` + `as.numeric()`.
- `capture.output(type='output')` for `Rprintf` logs.
- Atomic per-WU commits. Never amend.
- `R CMD INSTALL --preclean .` gates every C++ edit.
- One beads ticket per WU.

---

## File Structure

| File | WUs touching | Responsibility |
|---|---|---|
| `src/types.hpp` | A | Add `CalibConvergenceCfg`, `CalibSorCfg` to `CalibState` |
| `src/leafblower.h` | A | Extend `rk_params_t`, `rk_result_t`; update `EXPECTED_RK_PARAMS_BYTES` |
| `src/ieppa.hpp` | A, E | `IEPPAResult` extended (best_error, best_iter, sor_min_omega, sor_n_damped) |
| `src/raking.hpp` | A | `RakingResult` extended (best_error, best_iter, 4 new metrics) |
| `src/lbfgsb_solver.hpp` | A | `LBFGSBResult` extended (best_error, best_iter, 4 new metrics) |
| `src/c_api.cpp` | A | Marshal new config fields; zero-init `rk_result_t`; pack new result fields |
| `src/r_bridge.cpp` | A, E | Extend `res_list` from 14 → 15 elements; add `best_weights` REALSXP |
| `src/ieppa.cpp` | B, C, D, E | pct + criteria + SOR + best-iterate |
| `src/raking.cpp` | B, C, E | pct + criteria + best-iterate (no SOR) |
| `src/lbfgsb_solver.cpp` | B, C, E | pct + criteria + best-iterate (no SOR) |
| `R/harvest.R` | A, D, E | Rewrite `parse_convergence()`; add `sor` arg; expose `$sor` diagnostics |
| `python/leafblower/_bindings.cpp` | F | Mirror C ABI changes |
| `python/leafblower/_harvest.py` | F | Mirror R API changes |
| `python/leafblower/test_python.py` | F | Python parity tests |
| `tests/testthat/test-convergence-criteria.R` | B | New — A1, A2, A8 |
| `tests/testthat/test-quality-metrics.R` | C | New — A7 |
| `tests/testthat/test-sor.R` | D | New — A3, A4 |
| `tests/testthat/test-best-iterate.R` | E | New — A5, A6 |
| `tests/testthat/fixtures/stepstone_best_error_ref.rds` | E | New — A6 reference |
| `NEWS.md` | F | Default criterion change + new field list |
| `man/harvest.Rd` | F | Regenerated via `devtools::document()` |

---

## Work Unit A — Scaffold (no behavioral change)

**Beads:** `bd create --title "WU-A: convergence/SOR config scaffolding" --type task --priority 2`

**Files:** `src/types.hpp`, `src/leafblower.h`, `src/ieppa.hpp`, `src/raking.hpp`, `src/lbfgsb_solver.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`, `R/harvest.R`, `tests/testthat/test-convergence-criteria.R` (initial skeleton).

- [ ] **Step A1: Write failing default-off test.**

Create `tests/testthat/test-convergence-criteria.R`:

```r
test_that("convergence=list(absolute=1e-6) is backward compat (max_err criterion)", {
  set.seed(101)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(
    a = c("1"=0.3, "2"=0.5, "3"=0.2),
    b = c("1"=0.6, "2"=0.4)
  )
  # Must accept explicit absolute without any new fields — no error.
  w <- leafblower::harvest(
    data, target, max_weight = 3, method = "ieppa",
    convergence = list(absolute = 1e-6),
    attach_weights = FALSE
  )
  expect_true(is.numeric(as.numeric(w)))
  expect_length(as.numeric(w), n)
})

test_that("harvest accepts sor argument without error", {
  set.seed(101)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.5, "2"=0.5))
  # sor = NULL is disable; sor = list(auto=TRUE) enables auto mode.
  w1 <- leafblower::harvest(data, target, max_weight = 3, method = "ieppa",
                            sor = NULL, attach_weights = FALSE)
  w2 <- leafblower::harvest(data, target, max_weight = 3, method = "ieppa",
                            sor = list(auto = TRUE, omega_min = 0.3),
                            attach_weights = FALSE)
  expect_length(as.numeric(w1), n)
  expect_length(as.numeric(w2), n)
})
```

Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")'`
Expected: FAIL (new `sor` arg not accepted).

- [ ] **Step A2: Extend `src/types.hpp`.**

Append inside `namespace lbw { }` after existing WU-1..WU-5 overlay structs:

```cpp
enum class CalibCriterion : int {
    PCT      = 0,
    MAX_ERR  = 1,
    MEAN_ERR = 2,
    KL       = 3,
    CHI2     = 4
};
enum class CalibStopWhen : int { ANY = 0, ALL = 1 };

struct CalibConvergenceCfg {
    double        pct_tol      = 0.001;
    double        absolute_tol = 0.0;
    CalibCriterion criterion    = CalibCriterion::PCT;
    CalibStopWhen  stop_when    = CalibStopWhen::ANY;
};

struct CalibSorCfg {
    bool   enabled       = true;
    bool   auto_adapt    = true;
    double omega_init    = 1.0;
    double omega_min     = 0.3;
    double omega_fixed   = -1.0;  // sentinel: use auto
    int    burnin        = 20;
};
```

Attach to `CalibState`:

```cpp
CalibConvergenceCfg convergence_cfg;
CalibSorCfg         sor_cfg;
```

- [ ] **Step A3: Extend `src/leafblower.h`.**

Add enums and struct fields. Append after existing overlay fields in `rk_params_t`:

```c
/* Convergence config */
double pct_tol;          /* default 0.001 */
double absolute_tol;     /* default 0.0 */
int    criterion;        /* 0=PCT 1=MAX_ERR 2=MEAN_ERR 3=KL 4=CHI2 */
int    stop_when;        /* 0=ANY 1=ALL */

/* SOR config (iEPPA only; ignored by raking/lbfgsb) */
int    sor_enabled;
int    sor_auto;
double sor_omega_init;
double sor_omega_min;
double sor_omega_fixed;  /* -1.0 = use auto */
int    sor_burnin;
```

Append new fields to `rk_result_t`:

```c
double mean_error;       /* L1-over-margins */
double kl;               /* max KL across margins */
double chi2;             /* total chi-square */
double pct_change;       /* final iter max proportional weight change */
double best_error;       /* errRp at best iterate */
int    best_iter;
double sor_min_omega;    /* iEPPA only; non-iEPPA = 0.0 (zero-init by c_api.cpp) */
int    sor_n_damped;     /* iEPPA only; non-iEPPA = 0 (zero-init by c_api.cpp) */
```

Update `EXPECTED_RK_PARAMS_BYTES`: measure `sizeof(rk_params_t)` after adding 10 new fields (~80B depending on alignment), set the constant to the measured value. Update the layout comment block. The existing `static_assert` fails on first compile if forgotten — this is the safety net.

- [ ] **Step A4: Initialize new fields in `rk_params_init` and `rk_result_init`.**

In `src/c_api.cpp`, inside `rk_params_init`:

```c
p->pct_tol          = 0.001;
p->absolute_tol     = 0.0;
p->criterion        = 0;  // PCT
p->stop_when        = 0;  // ANY
p->sor_enabled      = 1;
p->sor_auto         = 1;
p->sor_omega_init   = 1.0;
p->sor_omega_min    = 0.3;
p->sor_omega_fixed  = -1.0;
p->sor_burnin       = 20;
```

In `rk_result_init`:

```c
r->mean_error     = 0.0;
r->kl             = 0.0;
r->chi2           = 0.0;
r->pct_change     = 0.0;
r->best_error     = std::numeric_limits<double>::infinity();
r->best_iter      = 0;
r->sor_min_omega  = 1.0;
r->sor_n_damped   = 0;
```

(Use explicit zero-init via memset at struct allocation if needed — confirm current style in existing `rk_result_init`.)

- [ ] **Step A5: Thread config fields through `src/c_api.cpp` into `CalibState`.**

In `rk_calibrate` (or equivalent entry point), copy:

```c
state.convergence_cfg.pct_tol      = params->pct_tol;
state.convergence_cfg.absolute_tol = params->absolute_tol;
state.convergence_cfg.criterion    = static_cast<CalibCriterion>(params->criterion);
state.convergence_cfg.stop_when    = static_cast<CalibStopWhen>(params->stop_when);
state.sor_cfg.enabled              = (params->sor_enabled != 0);
state.sor_cfg.auto_adapt           = (params->sor_auto != 0);
state.sor_cfg.omega_init           = params->sor_omega_init;
state.sor_cfg.omega_min             = params->sor_omega_min;
state.sor_cfg.omega_fixed           = params->sor_omega_fixed;
state.sor_cfg.burnin                = params->sor_burnin;
```

- [ ] **Step A6: Extend solver result structs.**

`src/ieppa.hpp` — extend `IEPPAResult`:

```cpp
struct IEPPAResult {
    // ... existing fields ...
    double mean_error     = 0.0;
    double kl             = 0.0;
    double chi2           = 0.0;
    double pct_change     = 0.0;
    double best_error     = std::numeric_limits<double>::infinity();
    int    best_iter      = 0;
    double sor_min_omega  = 1.0;
    int    sor_n_damped   = 0;
};
```

Do the same for `RakingResult` in `src/raking.hpp` and `LBFGSBResult` in `src/lbfgsb_solver.hpp`, minus the two `sor_*` fields (raking/lbfgsb never use SOR).

- [ ] **Step A7: Pack new result fields in `src/c_api.cpp`.**

After the solver returns, copy from the internal solver Result into `rk_result_t`:

```c
result->mean_error    = res.mean_error;
result->kl            = res.kl;
result->chi2          = res.chi2;
result->pct_change    = res.pct_change;
result->best_error    = res.best_error;
result->best_iter     = res.best_iter;
// For iEPPA only:
result->sor_min_omega = res.sor_min_omega;
result->sor_n_damped  = res.sor_n_damped;
// For raking/lbfgsb: skip sor_* (already zero-init'd)
```

- [ ] **Step A8: Extend `src/r_bridge.cpp` to unpack new config + expose new result fields.**

Currently `res_list` is a 14-element `VECSXP`. Extend to add:

```c
// New scalar result fields — append to the existing result list:
SET_VECTOR_ELT(result_list, 14, Rf_ScalarReal(r.mean_error));     // mean_error
SET_VECTOR_ELT(result_list, 15, Rf_ScalarReal(r.kl));             // kl
SET_VECTOR_ELT(result_list, 16, Rf_ScalarReal(r.chi2));           // chi2
SET_VECTOR_ELT(result_list, 17, Rf_ScalarReal(r.pct_change));     // pct_change
SET_VECTOR_ELT(result_list, 18, Rf_ScalarReal(r.best_error));     // best_error
SET_VECTOR_ELT(result_list, 19, Rf_ScalarInteger(r.best_iter));   // best_iter
SET_VECTOR_ELT(result_list, 20, Rf_ScalarReal(r.sor_min_omega));  // sor_min_omega
SET_VECTOR_ELT(result_list, 21, Rf_ScalarInteger(r.sor_n_damped));// sor_n_damped
// best_weights REALSXP is added in WU-E Step E6; placeholder now = NA vector length 0.
```

Adjust `Rf_allocVector(VECSXP, 22)` count. Update names vector. Also unpack new input args from the R-side argument list (the `convergence` list and `sor` list parsed by R are flattened into `rk_params_t` fields).

- [ ] **Step A9: Rewrite `parse_convergence()` in `R/harvest.R`.**

Replace lines 243–247:

```r
# Rewritten: pct is now a primary key. Deprecation warning removed.
# Returns list: pct_tol, absolute_tol, criterion, stop_when.
parse_convergence <- function(convergence) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  explicit_pct <- !is.null(convergence[["pct"]])
  explicit_abs <- !is.null(convergence[["absolute"]])
  pct_tol <- if (explicit_pct) convergence[["pct"]]
             else if (!explicit_abs) 0.001
             else 0.0
  absolute_tol <- convergence[["absolute"]] %||% 0.0
  criterion <- match.arg(
    convergence[["criterion"]] %||%
      (if (explicit_pct || !explicit_abs) "pct" else "max_err"),
    c("pct", "max_err", "mean_err", "kl", "chi2")
  )
  stop_when <- match.arg(convergence[["stop_when"]] %||% "any", c("any", "all"))
  list(pct_tol = pct_tol, absolute_tol = absolute_tol,
       criterion = criterion, stop_when = stop_when)
}
```

Update the call site at line 100 to consume the list:

```r
conv <- parse_convergence(convergence)
# ... forward through .Call as 4 separate args, not a single tol_abs ...
```

- [ ] **Step A10: Add `sor` argument to `harvest()` signature.**

Insert `sor = list(auto = TRUE, omega_min = 0.3)` as a new named argument in the signature (alongside the existing overlay args from prior plan). Add:

```r
parse_sor <- function(sor) {
  if (is.null(sor)) {
    return(list(enabled = 0L, auto = 0L, omega_init = 1.0,
                omega_min = 0.3, omega_fixed = -1.0, burnin = 20L))
  }
  list(
    enabled     = 1L,
    auto        = if (isTRUE(sor[["auto"]])) 1L else 0L,
    omega_init  = as.double(sor[["omega_init"]] %||% 1.0),
    omega_min   = as.double(sor[["omega_min"]] %||% 0.3),
    omega_fixed = as.double(sor[["omega"]] %||% -1.0),
    burnin      = as.integer(sor[["burnin"]] %||% 20L)
  )
}
```

Thread the parsed config through the `.Call` invocation.

- [ ] **Step A11: Build + run the failing tests from Step A1.**

```bash
R CMD INSTALL --preclean .
Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")'
```
Expected: both tests PASS.

- [ ] **Step A12: Full regression.**

`Rscript -e 'devtools::test()'`
Expected: `FAIL 0 | PASS ≥ 234` (baseline 232 + 2 new).

If any existing test fails due to the default criterion change (pct instead of max_err), pin that test's call by adding `convergence = list(absolute = 1e-6)`. The spec §Migration lists the 6 candidate files.

- [ ] **Step A13: Atomic commit.**

```bash
git add src/types.hpp src/leafblower.h src/ieppa.hpp src/raking.hpp \
        src/lbfgsb_solver.hpp src/c_api.cpp src/r_bridge.cpp R/harvest.R \
        tests/testthat/test-convergence-criteria.R
git commit -m "$(cat <<'EOF'
feat(WU-A): scaffold convergence/SOR config + extended result fields

Adds CalibConvergenceCfg, CalibSorCfg, CalibCriterion enum, CalibStopWhen
enum to types.hpp. Extends rk_params_t with pct_tol/absolute_tol/criterion/
stop_when and 6 sor_* fields. Extends rk_result_t with mean_error/kl/chi2/
pct_change/best_error/best_iter/sor_min_omega/sor_n_damped. Updates
EXPECTED_RK_PARAMS_BYTES. Rewrites parse_convergence() to remove deprecated
pct warning and use pct as primary key. Adds sor= argument to harvest().
No behavioral change — fields are plumbed but not yet consumed by solvers.
EOF
)"
```

- [ ] **Step A14: Close WU-A ticket.** `bd close <WU-A ticket id>`.

---

## Work Unit B — Convergence Criteria

**Beads:** `bd create --title "WU-B: pct + alternative criteria in all solvers" --type task --priority 2`

**Files:** `src/ieppa.cpp`, `src/raking.cpp`, `src/lbfgsb_solver.cpp`, `tests/testthat/test-convergence-criteria.R`.

- [ ] **Step B1: Write failing A2 test (pct criterion fires on smooth input).**

Append to `tests/testthat/test-convergence-criteria.R`:

```r
test_that("A2: pct=0.001 default converges on smooth synthetic", {
  set.seed(42)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE)),
    c = factor(sample(c("1","2","3","4"), n, replace = TRUE))
  )
  target <- list(
    a = c("1"=1/3, "2"=1/3, "3"=1/3),
    b = c("1"=0.5, "2"=0.5),
    c = c("1"=0.25, "2"=0.25, "3"=0.25, "4"=0.25)
  )
  w <- leafblower::harvest(
    data, target, max_weight = 10, method = "ieppa",
    max_iterations = 500,
    attach_weights = FALSE
  )
  result <- attr(w, "result")
  # pct criterion must have fired (pct_change below threshold at exit)
  expect_lt(result$pct_change, 0.001 * 1.5)  # allow 1.5x slack for inter-iter granularity
  expect_lt(result$max_error, 1e-3)
  # Status 0 = OK (not NOCONV=1)
  expect_equal(result$status, 0L)
})
```

- [ ] **Step B2: Write failing A8 tests (alternative criteria).**

Append to the same file:

```r
test_that("A8a: criterion='mean_err' actively stops solver", {
  set.seed(43)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
    convergence = list(absolute = 1e-4, criterion = "mean_err"),
    max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  expect_lt(result$mean_error, 1e-4)
  expect_equal(result$status, 0L)
})

test_that("A8b: criterion='kl' actively stops solver", {
  set.seed(44)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
    convergence = list(absolute = 1e-6, criterion = "kl"),
    max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  expect_lt(result$kl, 1e-6)
  expect_equal(result$status, 0L)
})

test_that("A8c: criterion='chi2' actively stops solver", {
  set.seed(45)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  # chi2 scales with n; threshold must be n-scaled (here ~1e-3 * 2000 = 2)
  w <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
    convergence = list(absolute = 2.0, criterion = "chi2"),
    max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  expect_lt(result$chi2, 2.0)
  expect_equal(result$status, 0L)
})
```

Run: FAIL expected (criteria not yet implemented).

- [ ] **Step B3: Implement pct + alternative criteria in `src/ieppa.cpp`.**

Inside `inner_solve_one_level` (or `ieppa_solve` depending on WU-3 structure):

Declare scratch vectors before the main loop:
```cpp
std::vector<double> W_prev(ct.M_cell);
// Initialize W_prev = W at entry (after initial W is set up from start_weights expansion):
W_prev = W;
```

Inside the `kErrCheckInterval` convergence check block (after `errRp` computed):

```cpp
// Compute pct_change
double pct_change = 0.0;
for (int c = 0; c < ct.M_cell; c++) {
    double rel = std::fabs(W[c] - W_prev[c]) / std::max(W_prev[c], 1e-12);
    if (rel > pct_change) pct_change = rel;
}

// Compute alternative metrics (all computed regardless of active criterion):
double W_total = 0.0;
for (int c = 0; c < ct.M_cell; c++) W_total += X[c];

double mean_err_sum = 0.0;
double kl_max       = 0.0;
double chi2_total   = 0.0;
constexpr double kMetricEps = 1e-10;
constexpr double kChi2Floor = 1.0;

for (int k = 0; k < st.K; k++) {
    const int nj  = st.cat_counts[k];
    const int off = cat_offset[k];
    std::fill(S_lin.begin(), S_lin.begin() + nj, 0.0);
    const int* gk = ct.g_per_cell[k].data();
    for (int c = 0; c < ct.M_cell; c++) {
        int j = gk[c];
        if (j >= 0 && j < nj) S_lin[j] += X[c];
    }
    double max_k = 0.0;
    double kl_k  = 0.0;
    for (int j = 0; j < nj; j++) {
        double S_p = S_lin[j] / W_total;
        double T   = st.targets[k][j];
        double err = std::fabs(S_p - T);
        if (err > max_k) max_k = err;
        if (T > 0.0) {
            kl_k += T * std::log((T + kMetricEps) / (S_p + kMetricEps));
        }
        double obs = S_lin[j];
        double exp = T * W_total;
        chi2_total += (obs - exp) * (obs - exp) / (exp + kChi2Floor);
    }
    mean_err_sum += max_k;
    if (kl_k > kl_max) kl_max = kl_k;
}
double mean_err = mean_err_sum / st.K;

// Active criterion dispatch:
bool converged_by_pct = (st.convergence_cfg.pct_tol > 0.0) &&
                        (pct_change < st.convergence_cfg.pct_tol);
double active_err = 0.0;
bool converged_by_absolute = false;
if (st.convergence_cfg.absolute_tol > 0.0) {
    switch (st.convergence_cfg.criterion) {
        case CalibCriterion::MAX_ERR:  active_err = errRp;    break;
        case CalibCriterion::MEAN_ERR: active_err = mean_err; break;
        case CalibCriterion::KL:       active_err = kl_max;   break;
        case CalibCriterion::CHI2:     active_err = chi2_total; break;
        case CalibCriterion::PCT:      active_err = pct_change; break;
    }
    converged_by_absolute = (active_err < st.convergence_cfg.absolute_tol);
}

bool converged = false;
bool have_pct  = (st.convergence_cfg.pct_tol > 0.0);
bool have_abs  = (st.convergence_cfg.absolute_tol > 0.0);
if (have_pct && have_abs) {
    converged = (st.convergence_cfg.stop_when == CalibStopWhen::ALL)
                ? (converged_by_pct && converged_by_absolute)
                : (converged_by_pct || converged_by_absolute);
} else if (have_pct) {
    converged = converged_by_pct;
} else if (have_abs) {
    converged = converged_by_absolute;
}

// Store last-seen metrics for §2 reporting:
res.mean_error = mean_err;
res.kl         = kl_max;
res.chi2       = chi2_total;
res.pct_change = pct_change;

W_prev = W;   // Update AFTER computing pct_change

if (converged) {
    res.status = RK_OK;
    break;
}
```

- [ ] **Step B4: Mirror pct + criteria in `src/raking.cpp`.**

Raking operates on obs-level weights. Add:

```cpp
std::vector<double> W_prev(st.n);
// Initialize W_prev = st.weights (post sw_vec normalization):
for (int i = 0; i < st.n; i++) W_prev[i] = st.weights[i];

// Inside each iteration's convergence check:
double pct_change = 0.0;
for (int i = 0; i < st.n; i++) {
    double rel = std::fabs(st.weights[i] - W_prev[i]) / std::max(W_prev[i], 1e-12);
    if (rel > pct_change) pct_change = rel;
}
// Same alternative-criteria computation as iEPPA (accumulate S_kj over obs-level).
// Same dispatch logic.
W_prev.assign(st.weights, st.weights + st.n);
```

Accumulating margins for raking requires iterating observations into per-(margin,category) sums (raking already does this for its errRp computation — extract as a helper or reuse).

- [ ] **Step B5: Mirror pct + criteria in `src/lbfgsb_solver.cpp`.**

L-BFGS-B computes weights via `fn.F(u[i])` from the dual. Compute pct across obs-level weights between successive outer iterations. Alternative criteria via the same margin-accumulation helper.

- [ ] **Step B6: Build + run WU-B tests.**

```bash
R CMD INSTALL --preclean .
Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")'
```
Expected: all 5 tests PASS (2 from WU-A + 3 new from A8 + 1 from A2).

- [ ] **Step B7: Audit 6 prior test files (A1).**

For each of `test-ieppa.R`, `test-compare.R`, `test-config-defaults.R`, `test-ieppa-nonuniform-d.R`, `test-compat.R`, `test-bounded-convergence.R`:
- grep for `harvest(` calls that do NOT pass `convergence = list(...)`.
- If the test asserts on `calib_result$max_error` or iteration counts, add `convergence = list(absolute = 1e-6)` to that call.
- If the test only asserts calibration quality (e.g., `diag$error_weighted < 1e-6`), no change needed.

Run `Rscript -e 'devtools::test()'` after each file's audit. Expected: no net FAIL regression.

- [ ] **Step B8: Full regression.** FAIL 0, PASS ≥ 237 (232 + 2 WU-A + 1 A2 + 3 A8 − 1 if any test obsolete).

- [ ] **Step B9: Atomic commit.**

```bash
git add src/ieppa.cpp src/raking.cpp src/lbfgsb_solver.cpp \
        tests/testthat/test-convergence-criteria.R \
        tests/testthat/test-ieppa.R tests/testthat/test-compare.R \
        tests/testthat/test-config-defaults.R \
        tests/testthat/test-ieppa-nonuniform-d.R tests/testthat/test-compat.R \
        tests/testthat/test-bounded-convergence.R
git commit -m "feat(WU-B): pct + alternative convergence criteria in all solvers

Implements pct_change computation + mean_err/kl/chi2 metrics + active
criterion dispatch + stop_when ANY/ALL logic in ieppa.cpp, raking.cpp,
lbfgsb_solver.cpp. Pinned 6 prior test files to convergence=list(absolute=1e-6)
where assertions depend on max_err criterion. A1 (backward compat) + A2
(pct quality on smooth) + A8a/b/c (mean_err/kl/chi2 active stop) tests pass."
```

- [ ] **Step B10: Close WU-B ticket.**

---

## Work Unit C — Quality Metrics Reporting

**Beads:** `bd create --title "WU-C: 5-metric reporting at exit in all solvers" --type task --priority 2`

**Files:** `src/ieppa.cpp`, `src/raking.cpp`, `src/lbfgsb_solver.cpp`, `tests/testthat/test-quality-metrics.R`.

- [ ] **Step C1: Write failing A7 test.**

Create `tests/testthat/test-quality-metrics.R`:

```r
test_that("A7: all 5 quality metrics present in calib_result for iEPPA", {
  set.seed(42)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2"), n, replace = TRUE)),
    b = factor(sample(c("1","2","3"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.5,"2"=0.5), b = c("1"=1/3,"2"=1/3,"3"=1/3))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                           max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  expected_names <- c("max_error", "mean_error", "kl", "chi2", "pct_change")
  for (nm in expected_names)
    expect_true(nm %in% names(result),
                info = sprintf("metric '%s' missing from calib_result", nm))
  expect_true(is.finite(result$max_error))
  expect_true(is.finite(result$mean_error))
  expect_true(is.finite(result$kl))
  expect_true(is.finite(result$chi2) || is.infinite(result$chi2))
  expect_true(is.finite(result$pct_change))
})

test_that("A7: all 5 quality metrics present in calib_result for raking", {
  set.seed(43)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2"), n, replace = TRUE)),
    b = factor(sample(c("1","2","3"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.5,"2"=0.5), b = c("1"=1/3,"2"=1/3,"3"=1/3))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "raking",
                           max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  for (nm in c("max_error", "mean_error", "kl", "chi2", "pct_change"))
    expect_true(nm %in% names(result))
})

test_that("A7: all 5 quality metrics present in calib_result for lbfgsb", {
  set.seed(44)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2"), n, replace = TRUE)),
    b = factor(sample(c("1","2","3"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.5,"2"=0.5), b = c("1"=1/3,"2"=1/3,"3"=1/3))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "lbfgsb",
                           max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  for (nm in c("max_error", "mean_error", "kl", "chi2", "pct_change"))
    expect_true(nm %in% names(result))
})
```

- [ ] **Step C2: Ensure metrics computed at exit in all three solvers.**

WU-B already computes mean_err, kl, chi2, pct_change at every `kErrCheckInterval`. WU-C guarantees these are computed at the FINAL iteration too (exit path) and packed into the result struct.

In `ieppa.cpp`, after the main loop breaks (whether via convergence or max_iter), ensure the metrics reflect the final state. If the loop exits before hitting a `kErrCheckInterval`, compute the metrics inline before the final return.

Same pattern in `raking.cpp` and `lbfgsb_solver.cpp`.

- [ ] **Step C3: Build + run A7 tests.**

```bash
R CMD INSTALL --preclean .
Rscript -e 'devtools::test_active_file("tests/testthat/test-quality-metrics.R")'
```
Expected: 3 tests PASS.

- [ ] **Step C4: Full regression.** FAIL 0, PASS ≥ 240.

- [ ] **Step C5: Atomic commit.**

```bash
git add src/ieppa.cpp src/raking.cpp src/lbfgsb_solver.cpp \
        tests/testthat/test-quality-metrics.R
git commit -m "feat(WU-C): compute all 5 quality metrics at exit across all solvers

Guarantees max_error, mean_error, kl, chi2, pct_change are populated in the
result struct for iEPPA, raking, and lbfgsb regardless of whether exit is via
convergence or max_iterations. A7 test verifies presence for all three methods."
```

- [ ] **Step C6: Close WU-C ticket.**

---

## Work Unit D — SOR Adaptive Under-Relaxation (iEPPA only)

**Beads:** `bd create --title "WU-D: SOR adaptive under-relaxation in iEPPA" --type task --priority 2`

**Files:** `src/ieppa.cpp`, `R/harvest.R` (expose `$sor` nested list), `tests/testthat/test-sor.R`.

- [ ] **Step D1: Write failing A3 + A4 tests.**

Create `tests/testthat/test-sor.R`:

```r
test_that("A3: SOR auto triggers on oscillatory tight-clamp synthetic", {
  set.seed(31415)
  n <- 5000
  data <- data.frame(
    v1 = factor(sample(c("A","B","C","D"), n, replace = TRUE)),
    v2 = factor(sample(c("X","Y","Z"), n, replace = TRUE)),
    v3 = factor(sample(c("1","2","3","4","5"), n, replace = TRUE)),
    v4 = factor(sample(c("p","q"), n, replace = TRUE)),
    v5 = factor(sample(c("a","b","c","d","e","f"), n, replace = TRUE))
  )
  target <- list(
    v1 = c(A=0.1, B=0.4, C=0.4, D=0.1),
    v2 = c(X=0.5, Y=0.3, Z=0.2),
    v3 = c("1"=0.1,"2"=0.1,"3"=0.4,"4"=0.3,"5"=0.1),
    v4 = c(p=0.7, q=0.3),
    v5 = c(a=0.05,b=0.05,c=0.5,d=0.2,e=0.15,f=0.05)
  )
  w <- leafblower::harvest(data, target, max_weight = 2, method = "ieppa",
                           max_iterations = 500,
                           sor = list(auto = TRUE, omega_min = 0.3),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  # SOR must have damped at least one margin below 0.9:
  expect_lt(result$sor_min_omega, 0.9)
})

test_that("A4: SOR silent on smooth input — no damping, no measurable overhead", {
  set.seed(202)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
                           max_iterations = 500,
                           sor = list(auto = TRUE, omega_min = 0.3),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_equal(result$sor_min_omega, 1.0)
  expect_equal(result$sor_n_damped, 0L)
})
```

- [ ] **Step D2: Implement SOR in `ieppa.cpp`.**

Inside `inner_solve_one_level` (or `ieppa_solve`) — add SOR state locals:

```cpp
std::vector<double> omega(st.K, st.sor_cfg.omega_init);
std::vector<double> prev_errRp_k(st.K, std::numeric_limits<double>::infinity());
std::vector<bool>   prev_decreasing(st.K, false);
double sor_min_omega = 1.0;
int    sor_n_damped  = 0;
const bool sor_active = st.sor_cfg.enabled;
const bool sor_auto   = st.sor_cfg.auto_adapt;
const double omega_fixed_val = st.sor_cfg.omega_fixed;
```

In the margin-sweep inner body where `f_new[j] = f_old[j] * (T_kj / S_kj)`:

```cpp
if (S_lin[j] < kEmptyBucketThreshold * ct.W_input) {
    // skip (existing empty-bucket path)
} else {
    const double ratio = (st.targets[k][j] * ct.W_input) / S_lin[j];
    double effective_omega;
    if (!sor_active) {
        effective_omega = 1.0;
    } else if (!sor_auto && omega_fixed_val > 0.0) {
        effective_omega = omega_fixed_val;
    } else {
        effective_omega = omega[k];
    }
    const double r = (effective_omega == 1.0) ? ratio
                                              : std::pow(ratio, effective_omega);
    f_lin[off + j] *= r;
}
```

After the per-margin sweep completes, if `sor_active && sor_auto && iter >= st.sor_cfg.burnin`, update omega[k] based on `errRp_k` sign-flip detection. Compute `errRp_k` per margin first.

**Scope note:** `compute_margin_errRp_linear` is defined as a lambda at `ieppa.cpp:435` inside the `ieppa_solve` function scope (from prior WU-4). Insert the SOR update block INSIDE the same scope (anywhere after the lambda is defined, inside `ieppa_solve` or `inner_solve_one_level` whichever contains the iteration loop) so the lambda capture-by-reference is valid. If the SOR block must live outside that scope for any reason, promote `compute_margin_errRp_linear` to a named static helper with explicit parameters; but the default path keeps it as an in-scope lambda call.

```cpp
if (sor_active && sor_auto && iter >= st.sor_cfg.burnin) {
    for (int k = 0; k < st.K; k++) {
        double errRp_k = compute_margin_errRp_linear(k);  // in-scope lambda from WU-4
        bool decreasing = (errRp_k < prev_errRp_k[k]);
        bool sign_flip  = !decreasing && prev_decreasing[k];
        if (sign_flip) {
            omega[k] = std::max(st.sor_cfg.omega_min, omega[k] * 0.7);
        } else if (decreasing) {
            omega[k] = std::min(1.0, omega[k] * 1.05);
        }
        prev_decreasing[k] = decreasing;
        prev_errRp_k[k]    = errRp_k;
        if (omega[k] < sor_min_omega) sor_min_omega = omega[k];
        if (omega[k] < 1.0) sor_n_damped++;
    }
}
```

At solver exit:

```cpp
res.sor_min_omega = sor_min_omega;
res.sor_n_damped  = sor_n_damped;
```

- [ ] **Step D3: Expose `calib_result$sor` nested list in R.**

In `R/harvest.R`, after receiving the raw result from `.Call`, reshape:

```r
calib_result <- raw$result
# Nest SOR diagnostics under $sor for namespace hygiene:
calib_result$sor <- list(
  min_omega = calib_result$sor_min_omega,
  n_damped  = calib_result$sor_n_damped
)
calib_result$sor_min_omega <- NULL
calib_result$sor_n_damped  <- NULL
```

Update tests to use `result$sor$min_omega` if preferred — but the test above accesses `result$sor_min_omega` directly. Choose one and keep consistent. Spec §3 says nested `$sor` — update test to match:

```r
expect_lt(result$sor$min_omega, 0.9)
expect_equal(result$sor$min_omega, 1.0)
expect_equal(result$sor$n_damped, 0L)
```

- [ ] **Step D4: Build + run WU-D tests.**

```bash
R CMD INSTALL --preclean .
Rscript -e 'devtools::test_active_file("tests/testthat/test-sor.R")'
```
Expected: PASS.

- [ ] **Step D5: Full regression.** FAIL 0, PASS ≥ 242.

- [ ] **Step D6: Atomic commit.**

```bash
git add src/ieppa.cpp R/harvest.R tests/testthat/test-sor.R
git commit -m "feat(WU-D): SOR adaptive under-relaxation (iEPPA only)

Per-margin omega adaptation with sign-flip oscillation detection. Burn-in 20
iters, floor omega_min=0.3 empirical, recovery factor 1.05 per monotone
decrease, damping factor 0.7 on sign flip. std::pow guard via
kEmptyBucketThreshold. Diagnostics exposed as calib_result\$sor\$min_omega and
calib_result\$sor\$n_damped. A3 (dampens on tight-clamp synthetic) + A4 (silent
on smooth) tests pass."
```

- [ ] **Step D7: Close WU-D ticket.**

---

## Work Unit E — Best-Iterate Tracking

**Beads:** `bd create --title "WU-E: best-iterate tracking in all solvers" --type task --priority 2`

**Files:** `src/ieppa.cpp`, `src/raking.cpp`, `src/lbfgsb_solver.cpp`, `src/r_bridge.cpp`, `R/harvest.R`, `tests/testthat/test-best-iterate.R`, `tests/testthat/fixtures/stepstone_best_error_ref.rds`.

- [ ] **Step E1: Write failing A5 + A6 tests.**

Create `tests/testthat/test-best-iterate.R`:

```r
test_that("A5: best_error <= max_error for iEPPA", {
  set.seed(42)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.5,"2"=0.5), b = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                           max_iterations = 200, attach_weights = FALSE)
  result <- attr(w, "result")
  expect_lte(result$best_error, result$max_error)
  expect_true(is.finite(result$best_error))
  expect_true(result$best_iter > 0)
})

test_that("A5: best_error <= max_error for raking", {
  set.seed(43)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.5,"2"=0.5), b = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "raking",
                           max_iterations = 200, attach_weights = FALSE)
  result <- attr(w, "result")
  expect_lte(result$best_error, result$max_error)
})

test_that("A5: best_weights is obs-level, length n, post-normalized", {
  set.seed(44)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_length(result$best_weights, n)
  # sum(best_weights) normalized to n:
  expect_equal(sum(result$best_weights), n, tolerance = 1e-6)
})

test_that("A6: stepstone best-iterate within 5% of reference", {
  skip_on_cran()
  ref_path <- test_path("fixtures/stepstone_best_error_ref.rds")
  skip_if(!file.exists(ref_path))
  ref <- readRDS(ref_path)
  # Regenerated by WU-E Step E7: ref is a numeric scalar (best_error from the
  # approved baseline run on stepstone-fulldata).
  fx <- "benchmarks/stepstone_fulldata_bench_data.parquet"
  tg <- "benchmarks/stepstone_fulldata_bench_targets.json"
  skip_if(!file.exists(fx) || !file.exists(tg))
  data <- arrow::read_parquet(fx)
  tgt_raw <- jsonlite::fromJSON(tg)
  target <- lapply(tgt_raw, function(x) { v <- unlist(x); v / sum(v) })
  w <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                           max_iterations = 3000,
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_lte(result$best_error, ref * 1.05)
})
```

- [ ] **Step E2: Implement best-iterate in `src/ieppa.cpp`.**

Inside `ieppa_solve` (not the inner loop — accumulator survives across homotopy levels):

```cpp
double best_errRp = std::numeric_limits<double>::infinity();
int    best_iter  = 0;
std::vector<double> W_best(ct.M_cell, 0.0);
```

After `errRp` computed at each `kErrCheckInterval`:

```cpp
if (errRp < best_errRp) {
    best_errRp = errRp;
    best_iter  = iter_total;  // across homotopy levels
    W_best     = W;
}
```

At exit, expand W_best to obs-level WITHOUT water-fill/bounds-clamping:

```cpp
std::vector<double> best_weights(st.n);
for (int i = 0; i < st.n; i++) {
    best_weights[i] = st.weights_initial[i] * W_best[ct.cell_of[i]];
}
double s = 0.0;
for (double w : best_weights) s += w;
if (s > 0.0) {
    const double scale = static_cast<double>(st.n) / s;
    for (double &w : best_weights) w *= scale;
}
// Store into res:
res.best_error = best_errRp;
res.best_iter  = best_iter;
res.best_weights = std::move(best_weights);
```

(Add `std::vector<double> best_weights` field to `IEPPAResult` in ieppa.hpp.)

- [ ] **Step E3: Mirror best-iterate in `raking.cpp` and `lbfgsb_solver.cpp`.**

For raking: W_best is n-element (obs-level) since raking operates directly on obs weights. No cell expansion; only sum-normalization at exit.

```cpp
std::vector<double> W_best(st.n, 0.0);
double best_errRp = std::numeric_limits<double>::infinity();
int    best_iter  = 0;

// At each check:
if (errRp < best_errRp) {
    best_errRp = errRp;
    best_iter  = iter;
    for (int i = 0; i < st.n; i++) W_best[i] = st.weights[i];
}

// At exit: sum-normalize
double s = 0.0;
for (double w : W_best) s += w;
if (s > 0.0) {
    const double scale = static_cast<double>(st.n) / s;
    for (double &w : W_best) w *= scale;
}
res.best_error   = best_errRp;
res.best_iter    = best_iter;
res.best_weights = std::move(W_best);
```

Same pattern in `lbfgsb_solver.cpp`.

- [ ] **Step E4: Extend `src/r_bridge.cpp` to return `best_weights` as REALSXP in result list.**

Update `res_list` element 22 to be `best_weights` REALSXP (length n, allocated + PROTECT'd):

```c
SEXP best_weights_sxp = PROTECT(Rf_allocVector(REALSXP, n));
double* bw = REAL(best_weights_sxp);
for (int i = 0; i < n; i++) bw[i] = r.best_weights[i];
SET_VECTOR_ELT(result_list, 22, best_weights_sxp);
```

Adjust `VECSXP` size to 23. Update names vector. Balance PROTECT/UNPROTECT.

- [ ] **Step E5: Ensure `attr(result, "result")$best_weights` is accessible.**

No R-side change needed — the result list already propagates via `attr(., "result")`.

- [ ] **Step E6: Generate stepstone_best_error_ref.rds.**

Run the benchmark once with WU-E code, capture `best_error`, save:

```r
data <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
tgt_raw <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
target <- lapply(tgt_raw, function(x) { v <- unlist(x); v / sum(v) })
w <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                         max_iterations = 3000, attach_weights = FALSE)
ref <- attr(w, "result")$best_error
saveRDS(ref, "tests/testthat/fixtures/stepstone_best_error_ref.rds")
cat(sprintf("stepstone best_error reference: %.3e\n", ref))
```

Expected: ref ≈ 2.22e-3.

- [ ] **Step E7: Build + run WU-E tests.**

```bash
R CMD INSTALL --preclean .
Rscript -e 'devtools::test_active_file("tests/testthat/test-best-iterate.R")'
```
Expected: PASS.

- [ ] **Step E8: Full regression.** FAIL 0, PASS ≥ 246.

- [ ] **Step E9: Atomic commit.**

```bash
git add src/ieppa.cpp src/ieppa.hpp src/raking.cpp src/raking.hpp \
        src/lbfgsb_solver.cpp src/lbfgsb_solver.hpp \
        src/r_bridge.cpp tests/testthat/test-best-iterate.R \
        tests/testthat/fixtures/stepstone_best_error_ref.rds
git commit -m "feat(WU-E): best-iterate tracking across all three solvers

Tracks W at minimum observed errRp; expanded + normalized to obs-level
best_weights at exit. No water-fill or bounds-clamping applied to W_best
(mid-loop snapshot). calib_result\$best_weights is element 22 of res_list.
calib_result\$best_error <= max_error always. A5/A6 tests pass."
```

- [ ] **Step E10: Close WU-E ticket.**

---

## Work Unit F — Python Parity + Documentation

**Beads:** `bd create --title "WU-F: Python parity + NEWS.md + roxygen" --type task --priority 2`

**Files:** `python/leafblower/_bindings.cpp`, `python/leafblower/_harvest.py`, `python/leafblower/test_python.py`, `NEWS.md`, `R/harvest.R` (roxygen), `man/harvest.Rd`.

- [ ] **Step F1: Extend `python/leafblower/_bindings.cpp`.**

Mirror the same new config fields (pct_tol, absolute_tol, criterion, stop_when, 6 sor_*) + result fields (mean_error, kl, chi2, pct_change, best_error, best_iter, sor_min_omega, sor_n_damped, best_weights) through pybind11.

- [ ] **Step F2: Extend `python/leafblower/_harvest.py`.**

Add Python-side parameters matching R defaults:
```python
def harvest(data, target, *,
    max_weight=5,
    method="ieppa",
    max_iterations=500,
    convergence=None,   # dict-like with keys pct, absolute, criterion, stop_when
    sor=None,           # dict with auto, omega_min, omega_fixed, burnin
    attach_weights=True,
    ...
):
    # Parse convergence mirror of parse_convergence() in R
    # Parse sor mirror of parse_sor() in R
    # Forward to pybind11 extension
```

- [ ] **Step F3: Write Python parity test.**

Extend `python/leafblower/test_python.py`:

```python
def test_convergence_criterion_default_is_pct():
    data = _make_fixture(n=1000, K=2)
    target = _target_fixture()
    res = leafblower.harvest(data, target, max_weight=5, method="ieppa",
                             attach_weights=False)
    r = res.attrs["result"]
    assert r["pct_change"] <= 0.001 * 1.5
    assert r["max_error"] < 1e-3

def test_all_5_metrics_present_in_python():
    data = _make_fixture(n=1000, K=2)
    target = _target_fixture()
    res = leafblower.harvest(data, target, max_weight=5, method="ieppa",
                             attach_weights=False)
    r = res.attrs["result"]
    for key in ("max_error", "mean_error", "kl", "chi2", "pct_change"):
        assert key in r

def test_best_weights_accessible_in_python():
    data = _make_fixture(n=1000, K=2)
    target = _target_fixture()
    res = leafblower.harvest(data, target, max_weight=5, method="ieppa",
                             attach_weights=False)
    r = res.attrs["result"]
    assert "best_weights" in r
    assert len(r["best_weights"]) == len(data)
```

Run: `pip install -e . && pytest python/leafblower/test_python.py -v`
Expected: PASS.

- [ ] **Step F4: Update R roxygen for `harvest()`.**

In `R/harvest.R`, update `@param convergence`:

```r
#' @param convergence Named list controlling the stopping criterion.
#'   Accepted keys:
#'   \itemize{
#'     \item \code{pct}: proportional weight-change threshold.
#'       Default \code{0.001} (0.1 percent) when neither \code{pct} nor
#'       \code{absolute} is provided. Set to 0 to disable.
#'     \item \code{absolute}: absolute threshold for the active criterion.
#'       Default \code{0} (disabled). When set without \code{pct},
#'       \code{pct} is disabled and \code{absolute} is applied to the
#'       \code{criterion} metric.
#'     \item \code{criterion}: one of \code{"pct"} (default),
#'       \code{"max_err"}, \code{"mean_err"}, \code{"kl"}, \code{"chi2"}.
#'       Note: \code{chi2} scales with sample size; supply an n-scaled
#'       \code{absolute} threshold.
#'     \item \code{stop_when}: \code{"any"} (default) or \code{"all"} —
#'       applied when both \code{pct} and \code{absolute} are non-zero.
#'   }
#'   Backward compat: \code{list(absolute = 1e-6)} preserves the
#'   pre-v0.2 behavior (max_error criterion).
```

Add `@param sor`:

```r
#' @param sor Named list controlling SOR adaptive under-relaxation
#'   (iEPPA only). \code{NULL} disables SOR. Accepted keys:
#'   \itemize{
#'     \item \code{auto}: logical, default TRUE. Enables per-margin
#'       adaptive omega.
#'     \item \code{omega_min}: lower bound on omega. Default 0.3.
#'     \item \code{omega}: fixed omega (disables auto-adapt).
#'     \item \code{burnin}: iterations before adaptation begins. Default 20.
#'   }
```

Run `Rscript -e 'devtools::document()'` to regenerate `man/harvest.Rd`.

- [ ] **Step F5: Create `NEWS.md`.**

```markdown
# leafblower NEWS

## leafblower 0.2.0 (unreleased)

### Breaking changes

- **Default convergence criterion changed** from `absolute = 1e-6`
  (max marginal error) to `pct = 0.001` (max proportional weight change).
  To preserve prior behavior: `convergence = list(absolute = 1e-6)`.
- `convergence[["pct"]]` was previously a deprecated key that triggered a
  warning and was discarded. It is now the primary convergence threshold.
  Users who suppressed the deprecation warning and set `pct = X`:
  your value will now be honored.

### New features

- Pluggable convergence criteria: `convergence = list(criterion = "kl")`
  accepts `"pct"` (default), `"max_err"`, `"mean_err"`, `"kl"`, `"chi2"`.
- Five quality metrics always reported: `max_error`, `mean_error`, `kl`,
  `chi2`, `pct_change` in `attr(result, "result")`.
- SOR adaptive under-relaxation for iEPPA via new `sor` argument
  (default: auto-enabled).
- Best-iterate tracking: `attr(result, "result")$best_weights` is always
  populated with weights at the iteration with minimum observed errRp.

### Python

Python `harvest()` API is extended with matching keyword arguments.
See `python/leafblower/_harvest.py` docstring.
```

- [ ] **Step F6: Final CRAN gate + pytest.**

```bash
R CMD build .
R CMD check --as-cran leafblower_*.tar.gz
pip install -e .
pytest python/leafblower/ -v
```
Expected: 0 ERROR, 0 WARNING, NOTEs ≤ baseline. pytest PASS.

- [ ] **Step F7: Full regression.** FAIL 0, PASS ≥ 246.

- [ ] **Step F8: Atomic commit.**

```bash
git add python/leafblower/_bindings.cpp python/leafblower/_harvest.py \
        python/leafblower/test_python.py NEWS.md R/harvest.R man/
git commit -m "docs+bindings(WU-F): Python parity + NEWS.md + roxygen

Python harvest() accepts matching convergence and sor arguments; all 5
quality metrics and best_weights exposed through pybind11. Parity tests
mirror test-convergence-criteria.R / test-quality-metrics.R /
test-best-iterate.R. NEWS.md documents the pct=0.001 default change,
the pct deprecation-reversal migration note, and the new result fields.
Roxygen @param convergence enumerates all keys with defaults; @param sor
added. R CMD check --as-cran clean."
```

- [ ] **Step F9: Close WU-F ticket.**

---

## Final Verification

- [ ] `devtools::test()` → FAIL 0, PASS ≥ 246
- [ ] `pytest python/leafblower/` → green
- [ ] `R CMD check --as-cran` → 0 ERROR, 0 WARNING, NOTEs ≤ baseline
- [ ] 6 beads tickets (WU-A..WU-F) all closed
- [ ] `benchmarks/stepstone_fulldata_homotopy.R` re-run shows `attr(result,"result")$best_error ≈ 2.22e-3` (matches stepstone_best_error_ref.rds)
- [ ] NEWS.md updated

## Self-review

- **Spec coverage:** every spec section §1–§10 mapped to at least one WU task. §1 criteria → WU-B; §2 metrics → WU-C; §3 SOR → WU-D; §4 best-iterate → WU-E; §5 C ABI → WU-A; §6 files → distributed across WUs; §7 AC → A1–A9 mapped to WU-B/WU-C/WU-D/WU-E/WU-F tests.
- **Placeholders:** none. Every step has concrete code/commands. No "TBD", "similar to Task N", "add appropriate".
- **Type consistency:** `CalibConvergenceCfg`, `CalibSorCfg`, `CalibCriterion`, `CalibStopWhen` used consistently. `pct_tol`, `absolute_tol`, `criterion`, `stop_when` field names match spec §5 C ABI table. R-side arg names (`convergence$pct`, `convergence$absolute`, `convergence$criterion`, `convergence$stop_when`, `sor$auto`, `sor$omega_min`, `sor$omega`) match the ground-truth `parse_convergence()` and `parse_sor()` contracts.
- **R API names ground-truthed:** `target`, `method`, `max_iterations`, `convergence = list(...)` — verified against spec §1 edge case table. No `targets`/`algorithm`/`max_iter`/bare `absolute=` anywhere in test code.
- **Result field names:** `max_error`, `mean_error`, `kl`, `chi2`, `pct_change`, `best_error`, `best_iter`, `best_weights`, `sor$min_omega`, `sor$n_damped` — consistent with spec §2 + §3.
- **Backward compat coverage:** 6 prior test files audited in WU-B Step B7 per spec Migration section.
