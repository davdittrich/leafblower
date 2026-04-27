# Improvement-Based Convergence Criterion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `criterion = "improvement"` — stop when the loss function (errRp) stops improving between checks, i.e., `(errRp[t-1] - errRp[t]) / errRp[t-1] < pct_tol`. Make this the new default, replacing `criterion = "pct"` (weight-change stopping). Natural plateau detection fires on both slow convergence AND oscillation; combined with best-iterate tracking, this returns the best-seen weights when the solver stalls.

**Mathematical basis:** Convergence in the optimization sense = the objective function no longer decreases meaningfully. `criterion="pct"` (weight change) is a proxy that can misfire: weights can stop moving while calibration error remains large (infeasible stall). `criterion="improvement"` measures what we actually care about: is the calibration error still getting better?

**Behavior on oscillating problems (stepstone):** When errRp diverges (goes up), relative improvement is negative → fires immediately → best_iterate returns the best-seen weights. This is the correct behavior: stop wasting iterations once the solver is diverging.

**Architecture:** Five atomic work units. TDD throughout. Backward compat preserved: `convergence = list(absolute = 1e-6)` still works; `convergence = list(pct = 0.001)` still means weight-change stopping. `convergence = list(improvement = 0.001)` is the new default.

**Tech Stack:** C++17 (types.hpp, leafblower.h, ieppa.cpp, raking.cpp, lbfgsb_solver.cpp, c_api.cpp), R (harvest.R), Python (_harvest.py), testthat 3, pytest.

**Baseline:** `devtools::test()` → FAIL 0 | PASS 298

**Constants:**
- `RK_OK=0`, `RK_ERR_NOCONV=1` (leafblower.h:32-33)
- Current `CalibCriterion`: PCT=0, MAX_ERR=1, MEAN_ERR=2, KL=3, CHI2=4 (types.hpp:37-41)
- New: `IMPROVEMENT=5`

---

## File Structure

| File | WU | Change |
|---|---|---|
| `src/types.hpp` | A | Add `CalibCriterion::IMPROVEMENT = 5` |
| `src/leafblower.h` | A | Add `criterion=5` comment; update rk_params_init default |
| `src/c_api.cpp` | A | Change `rk_params_init` default criterion from 0→5 |
| `src/validation.hpp` | A | Update valid range comment (0-5 not 0-4); update error message |
| `src/r_bridge.cpp` | A | Update criterion range guard (0-4 → 0-5) |
| `src/ieppa.cpp` | B | Add `prev_errRp_improvement` tracking + IMPROVEMENT dispatch |
| `src/raking.cpp` | C | Mirror ieppa IMPROVEMENT tracking |
| `src/lbfgsb_solver.cpp` | D | Add IMPROVEMENT (start→final improvement, batch semantics) |
| `R/harvest.R` | E | Add "improvement" to match.arg; change default; update docs |
| `python/leafblower/_harvest.py` | E | Mirror criterion map |
| `NEWS.md` | E | Document new default |

---

## Implementation Notes

### IMPROVEMENT criterion semantics

```
relative_improvement = (errRp_prev - errRp_curr) / errRp_prev
```

- `errRp_prev` initialized to `+∞`
- First kErrCheckInterval: `improvement = (∞ - errRp) / ∞ = 1.0` → NOT triggered (1.0 >> 0.001)
- Subsequent checks:
  - Decreasing errRp: `improvement ∈ (0, 1)` → fire when `< pct_tol`
  - Plateau: `improvement ≈ 0` < pct_tol → CONVERGE
  - Diverging: `improvement < 0` < pct_tol → CONVERGE (then best_iterate returns minimum)
- Zero errRp: guard with `if (errRp_prev == 0) improvement = 1.0` (already converged last check)

### Dispatch integration

IMPROVEMENT uses `pct_tol` (same field as PCT). In the convergence dispatch:

```cpp
bool converged_pct = (cfg.pct_tol > 0.0) && (
    (cfg.criterion == CalibCriterion::PCT         && pct_change < cfg.pct_tol) ||
    (cfg.criterion == CalibCriterion::IMPROVEMENT && relative_improvement < cfg.pct_tol)
);
```

The `absolute_tol` / `active_val` switch: add `case CalibCriterion::IMPROVEMENT: active_val = relative_improvement; break;` for the `stop_when=ALL` use case.

### Variable: `prev_errRp_improvement`

Solver-local variable, separate from:
- `X_prev` / `W_prev` (used for pct_change = weight tracking)
- `prev_errRp_k` (per-margin, used by SOR)

Initialized to `std::numeric_limits<double>::infinity()`. Updated AFTER computing `relative_improvement` (same pattern as X_prev update order for pct_change).

---

## Work Unit A — Scaffold: Add IMPROVEMENT to C ABI + types

**Beads ticket:** `bd create --title "feat: add CalibCriterion::IMPROVEMENT=5 to C ABI and types" --type feature --priority 1`

**Files:** `src/types.hpp`, `src/leafblower.h`, `src/c_api.cpp`, `src/validation.hpp`, `src/r_bridge.cpp`

- [ ] **Step A0: Create + claim ticket**
```bash
bd create --title "feat: add CalibCriterion::IMPROVEMENT=5 to C ABI and types" \
  --type feature --priority 1 \
  --description "Add IMPROVEMENT=5 to CalibCriterion enum. Update C ABI range guard from [0,4] to [0,5]. Change rk_params_init default criterion from 0 (PCT) to 5 (IMPROVEMENT). No solver behavior change yet — solvers fall through to existing PCT logic until WU-B/C/D." 2>&1 | tail -2
# Claim the returned ticket ID
```

- [ ] **Step A1: Extend `src/types.hpp`**

Add `IMPROVEMENT = 5` to `CalibCriterion`:
```cpp
enum class CalibCriterion : int {
    PCT         = 0,
    MAX_ERR     = 1,
    MEAN_ERR    = 2,
    KL          = 3,
    CHI2        = 4,
    IMPROVEMENT = 5   // relative improvement in errRp between checks
};
```

Change default in `CalibConvergenceCfg`:
```cpp
CalibCriterion criterion = CalibCriterion::IMPROVEMENT;  // was PCT
```

- [ ] **Step A2: Update `src/leafblower.h`**

In `rk_params_t` comment block, update criterion documentation:
```c
int criterion;  /* 0=PCT 1=MAX_ERR 2=MEAN_ERR 3=KL 4=CHI2 5=IMPROVEMENT */
```

- [ ] **Step A3: Update `src/c_api.cpp` rk_params_init**

Change:
```c
p->criterion = 0;  /* PCT */
```
To:
```c
p->criterion = 5;  /* IMPROVEMENT — mathematical convergence: loss stops improving */
```

- [ ] **Step A4: Update `src/validation.hpp` range guard**

Change:
```cpp
if (p->criterion < 0 || p->criterion > 4)
    return err("criterion out of range [0,4]: 0=PCT 1=MAX_ERR 2=MEAN_ERR 3=KL 4=CHI2");
```
To:
```cpp
if (p->criterion < 0 || p->criterion > 5)
    return err("criterion out of range [0,5]: 0=PCT 1=MAX_ERR 2=MEAN_ERR 3=KL 4=CHI2 5=IMPROVEMENT");
```

- [ ] **Step A5: Update `src/r_bridge.cpp` range guard**

Change:
```cpp
if (p.criterion < 0 || p.criterion > 4)
    Rf_error("leafblower: invalid arguments — criterion out of range [0,4]...");
```
To:
```cpp
if (p.criterion < 0 || p.criterion > 5)
    Rf_error("leafblower: invalid arguments — criterion out of range [0,5]"
             " (0=PCT 1=MAX_ERR 2=MEAN_ERR 3=KL 4=CHI2 5=IMPROVEMENT)");
```

- [ ] **Step A6: Build gate**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```
Expected: `* DONE (leafblower)`. The new default (IMPROVEMENT=5) falls through to PCT behavior in solvers because IMPROVEMENT has no dispatch branch yet — all solvers will behave as if neither pct_tol nor absolute_tol is satisfied and exit at max_iter. Regression tests may change! See Step A7.

- [ ] **Step A7: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
If tests fail because the default criterion changed to IMPROVEMENT and it's unimplemented (falls through to never-converge): add `convergence = list(absolute = 1e-6)` to affected tests OR run after WU-B implements IMPROVEMENT (recommended — defer full regression to after WU-B).

- [ ] **Step A8: Commit**
```bash
git add src/types.hpp src/leafblower.h src/c_api.cpp src/validation.hpp src/r_bridge.cpp
git commit -m "$(cat <<'EOF'
feat(types): add CalibCriterion::IMPROVEMENT=5 to C ABI

Mathematical convergence criterion: stop when errRp stops improving
(relative_improvement < pct_tol). Replaces PCT (weight-change) as
default in rk_params_init and CalibConvergenceCfg. Solvers dispatch
added in subsequent commits. C ABI range guard updated [0,4]→[0,5].
EOF
)"
```

- [ ] **Step A9: Close ticket** `bd close <ticket-id>`

---

## Work Unit B — ieppa.cpp IMPROVEMENT dispatch

**Beads ticket:** `bd create --title "feat(ieppa): IMPROVEMENT criterion — relative errRp improvement tracking" --type feature --priority 1`

**Files:** `src/ieppa.cpp`

- [ ] **Step B0: Create + claim ticket**

- [ ] **Step B1: Write failing test (TDD)**

Append to `tests/testthat/test-convergence-criteria.R`:

```r
test_that("B_improvement: iEPPA converges on smooth input with improvement criterion", {
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
    convergence = list(improvement = 0.001),
    attach_weights = FALSE
  )
  result <- attr(w, "result")
  # Must converge (not NOCONV) and max_error must be reasonable
  expect_equal(result$status, 0L)
  expect_lt(result$max_error, 1e-3)
})

test_that("B_improvement: iEPPA stops early on oscillating input (no waste)", {
  # Oscillating problem: improvement criterion fires when errRp stops improving
  # Returns best_error < max_error (best iterate used)
  set.seed(31415)
  n <- 2000
  data <- data.frame(
    v1 = factor(sample(c("A","B","C","D"), n, replace = TRUE)),
    v2 = factor(sample(c("X","Y","Z"), n, replace = TRUE))
  )
  target <- list(
    v1 = c(A=0.1, B=0.4, C=0.4, D=0.1),
    v2 = c(X=0.5, Y=0.3, Z=0.2)
  )
  w_imp <- leafblower::harvest(
    data, target, max_weight = 2, method = "ieppa",
    max_iterations = 500,
    convergence = list(improvement = 0.001),
    attach_weights = FALSE
  )
  res_imp <- attr(w_imp, "result")
  # Improvement criterion should stop before 500 iters once errRp plateaus
  # best_error <= max_error by definition
  expect_lte(res_imp$best_error, res_imp$max_error)
  # Should not exhaust all 500 iterations (plateaus earlier)
  expect_lt(res_imp$iterations, 500L)
})
```

Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")' 2>&1 | tail -8`
Expected: FAIL (IMPROVEMENT not yet dispatched in ieppa).

- [ ] **Step B2: Add `prev_errRp_improvement` tracking in `src/ieppa.cpp`**

**Read ieppa.cpp first.** Find where solver-local variables are declared before the main loop (around lines 270-295, after `X_prev` and `W_best` declarations).

Add after existing declarations:
```cpp
// IMPROVEMENT criterion: track previous errRp to measure relative improvement.
// Initialized to ∞ so first check (improvement=1.0) never triggers convergence.
double prev_errRp_improvement = std::numeric_limits<double>::infinity();
```

- [ ] **Step B3: Add IMPROVEMENT dispatch in `src/ieppa.cpp`**

Find the convergence dispatch block (around lines 940-980, the `converged_pct` computation). Currently:
```cpp
bool converged_pct = (cfg.pct_tol > 0.0) && (pct_change < cfg.pct_tol);
```

Replace with:
```cpp
// Compute relative improvement in errRp for IMPROVEMENT criterion.
double relative_improvement = 1.0;  // default: always "improving" until second check
if (prev_errRp_improvement < std::numeric_limits<double>::infinity()) {
    if (prev_errRp_improvement <= 0.0) {
        relative_improvement = 1.0;  // already at zero — converged, handled by errRp check
    } else {
        relative_improvement = (prev_errRp_improvement - errRp) / prev_errRp_improvement;
    }
}
prev_errRp_improvement = errRp;  // update AFTER computing (same pattern as X_prev)

bool converged_pct = (cfg.pct_tol > 0.0) && (
    (cfg.criterion == lbw::CalibCriterion::PCT         && pct_change < cfg.pct_tol) ||
    (cfg.criterion == lbw::CalibCriterion::IMPROVEMENT && relative_improvement < cfg.pct_tol)
);
```

Also add IMPROVEMENT to the `active_val` switch (for `stop_when=ALL` support):
```cpp
case lbw::CalibCriterion::IMPROVEMENT: active_val = relative_improvement; break;
```

**IMPORTANT:** `prev_errRp_improvement` must be RESET at the start of each homotopy level (same pattern as `X_prev` reset at line ~361). Find the homotopy level reset block and add:
```cpp
prev_errRp_improvement = std::numeric_limits<double>::infinity();
```

- [ ] **Step B4: Build gate**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

- [ ] **Step B5: Run test**
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")' 2>&1 | tail -8
```
Expected: B_improvement tests PASS. Second test (oscillating, stops early) — verify `result$iterations < 500`.

- [ ] **Step B6: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 300. If tests fail due to the new default (criterion=IMPROVEMENT), they had harvest() calls without pinned convergence that now behave differently. For tests that assert `status == 0` on inputs that converge smoothly: add `convergence = list(improvement = 0.001)` (new default) — they should pass. For tests that pinned `convergence = list(absolute = 1e-6)`: unchanged.

- [ ] **Step B7: Commit + close**
```bash
git add src/ieppa.cpp tests/testthat/test-convergence-criteria.R
git commit -m "$(cat <<'EOF'
feat(ieppa): IMPROVEMENT criterion — relative errRp improvement tracking

Stop when (errRp_prev - errRp_curr) / errRp_prev < pct_tol.
Fires on plateau (improvement near 0) AND on divergence (negative
improvement = errRp going up). Combined with best-iterate, returns
best-seen weights when solver oscillates. prev_errRp_improvement
reset at homotopy level boundaries.
EOF
)"
bd close <ticket-id>
```

---

## Work Unit C — raking.cpp IMPROVEMENT dispatch

**Beads ticket:** `bd create --title "feat(raking): IMPROVEMENT criterion dispatch" --type feature --priority 1`

**Files:** `src/raking.cpp`

- [ ] **Step C0: Create + claim ticket**

- [ ] **Step C1: Write failing test**

Append to `tests/testthat/test-convergence-criteria.R`:

```r
test_that("C_improvement: raking converges with improvement criterion", {
  set.seed(43)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight=10, method="raking",
                           max_iterations=500,
                           convergence = list(improvement = 0.001),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_equal(result$status, 0L)
  expect_lt(result$max_error, 1e-3)
})
```

Run: FAIL expected.

- [ ] **Step C2: Mirror in `src/raking.cpp`**

Read raking.cpp. Find the convergence dispatch block. Add `prev_errRp_improvement` (same as ieppa: initialized to ∞, reset at... raking has no homotopy levels, so just initialize before the main loop).

Same dispatch logic:
```cpp
double prev_errRp_improvement = std::numeric_limits<double>::infinity();
// ... inside kErrCheckInterval block:
double relative_improvement = 1.0;
if (prev_errRp_improvement < std::numeric_limits<double>::infinity()) {
    relative_improvement = (prev_errRp_improvement > 0.0)
        ? (prev_errRp_improvement - errRp) / prev_errRp_improvement
        : 1.0;
}
prev_errRp_improvement = errRp;

bool converged_pct = (cfg.pct_tol > 0.0) && (
    (cfg.criterion == lbw::CalibCriterion::PCT         && pct_change < cfg.pct_tol) ||
    (cfg.criterion == lbw::CalibCriterion::IMPROVEMENT && relative_improvement < cfg.pct_tol)
);
```

- [ ] **Step C3: Build + test + regression**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")' 2>&1 | tail -5
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: C_improvement test PASS; FAIL 0.

- [ ] **Step C4: Commit + close**
```bash
git add src/raking.cpp tests/testthat/test-convergence-criteria.R
git commit -m "feat(raking): IMPROVEMENT criterion dispatch"
bd close <ticket-id>
```

---

## Work Unit D — lbfgsb_solver.cpp IMPROVEMENT dispatch

**Beads ticket:** `bd create --title "feat(lbfgsb): IMPROVEMENT criterion (start→final improvement)" --type feature --priority 1`

**Files:** `src/lbfgsb_solver.cpp`

**Note:** lbfgsb is a batch solver — single pass, no inner iteration loop. IMPROVEMENT measures start→final errRp improvement: `(errRp_start - errRp_final) / errRp_start`. If this is < pct_tol, the solver made negligible progress — return as-is (best_weights = final weights, which is the only iterate).

- [ ] **Step D0: Create + claim ticket**

- [ ] **Step D1: Implement in `src/lbfgsb_solver.cpp`**

Read the file to find `compute_final_weights_and_error` and `lbfgsb_solve`.

In `lbfgsb_solve`, before the solve call, record:
```cpp
const double errRp_start = res.max_error;  // initial errRp before optimization
```

After `compute_final_weights_and_error` returns (which sets `res.max_error`):
```cpp
// IMPROVEMENT criterion for lbfgsb: start→final errRp improvement (batch semantics).
// Unlike iEPPA/raking, lbfgsb has no inner iteration loop — measure total improvement.
double relative_improvement_lbfgsb = 1.0;
if (errRp_start > 0.0 && std::isfinite(errRp_start)) {
    relative_improvement_lbfgsb =
        (errRp_start - res.max_error) / errRp_start;
}
```

Add to convergence dispatch (lbfgsb has its own dispatch logic — read and extend):
```cpp
bool converged_pct = (cfg.pct_tol > 0.0) && (
    (cfg.criterion == lbw::CalibCriterion::PCT         && pct_change < cfg.pct_tol) ||
    (cfg.criterion == lbw::CalibCriterion::IMPROVEMENT && relative_improvement_lbfgsb < cfg.pct_tol)
);
```

- [ ] **Step D2: Build + test + regression**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
Rscript -e 'devtools::test()' 2>&1 | tail -3
```

- [ ] **Step D3: Commit + close**
```bash
git add src/lbfgsb_solver.cpp
git commit -m "feat(lbfgsb): IMPROVEMENT criterion (start->final errRp improvement)"
bd close <ticket-id>
```

---

## Work Unit E — R/harvest.R + Python + NEWS.md

**Beads ticket:** `bd create --title "feat(harvest.R): improvement criterion as default + docs + Python mirror" --type feature --priority 1`

**Files:** `R/harvest.R`, `python/leafblower/_harvest.py`, `NEWS.md`

- [ ] **Step E0: Create + claim ticket**

- [ ] **Step E1: Update `R/harvest.R` parse_convergence**

**Read R/harvest.R first.** Find `parse_convergence` (line ~286).

1. Add `"improvement"` to `valid_keys`:
```r
valid_keys <- c("pct", "absolute", "criterion", "stop_when", "improvement")
```

Wait — `improvement` is not a separate key in the convergence list. Looking at the API design:

The user passes: `convergence = list(improvement = 0.001)`. This means:
- `convergence[["improvement"]]` sets the threshold
- `criterion` is implicitly "improvement"

BUT the current design uses `pct = X` with `criterion = "pct"`. To match this pattern:
- `convergence = list(improvement = 0.001)` → `pct_tol = 0.001`, `criterion = "improvement"`
- `convergence = list(pct = 0.001)` → `pct_tol = 0.001`, `criterion = "pct"` (unchanged)

This means `improvement` is a NEW key in the convergence list (like `pct` and `absolute`).

Update `parse_convergence`:
```r
parse_convergence <- function(convergence) {
  if (!is.null(convergence) && !is.list(convergence))
    stop("convergence must be a named list or NULL")
  valid_keys <- c("pct", "absolute", "criterion", "stop_when", "improvement")
  bad <- setdiff(names(convergence), valid_keys)
  if (length(bad))
    stop(sprintf("Unknown convergence key(s): %s. Valid keys: %s",
                 paste(bad, collapse=", "), paste(valid_keys, collapse=", ")))
  `%||%` <- function(a, b) if (is.null(a)) b else a
  explicit_pct  <- !is.null(convergence[["pct"]])
  explicit_abs  <- !is.null(convergence[["absolute"]])
  explicit_impr <- !is.null(convergence[["improvement"]])

  # pct_tol: threshold for either pct (weight change) or improvement criterion
  pct_tol <- if (explicit_impr) convergence[["improvement"]]
             else if (explicit_pct) convergence[["pct"]]
             else if (!explicit_abs) 0.001   # default improvement threshold
             else 0.0
  absolute_tol <- convergence[["absolute"]] %||% 0.0

  # default criterion:
  #   "improvement" when neither pct nor abs is explicit (new default)
  #   "pct"         when user explicitly sets pct=
  #   "max_err"     when user sets absolute= without pct/improvement
  criterion <- match.arg(
    convergence[["criterion"]] %||%
      (if (explicit_impr) "improvement"
       else if (explicit_pct) "pct"
       else if (!explicit_abs) "improvement"
       else "max_err"),
    c("pct", "max_err", "mean_err", "kl", "chi2", "improvement")
  )
  stop_when <- match.arg(convergence[["stop_when"]] %||% "any", c("any", "all"))
  list(pct_tol = pct_tol, absolute_tol = absolute_tol,
       criterion = criterion, stop_when = stop_when)
}
```

2. Update criterion integer map to include `improvement = 5L`:
```r
criterion_int <- c(pct = 0L, max_err = 1L, mean_err = 2L, kl = 3L, chi2 = 4L,
                   improvement = 5L)
```

3. Update `@param convergence` roxygen:
```r
#'   \item{\code{improvement}}: relative improvement threshold in marginal
#'     calibration error (errRp) between checks (default \code{0.001}).
#'     Stop when \code{(errRp_prev - errRp_curr) / errRp_prev < 0.001}.
#'     Fires on plateau (errRp barely moves) AND oscillation (errRp increases).
#'     Combined with \code{best_weights}, returns minimum-error iterate.
#'     \strong{This is the default criterion.}
```

- [ ] **Step E2: Verify backward compatibility**

Edge cases in `parse_convergence`:
- `list()` → pct_tol=0.001, criterion="improvement" ✅ (new default)
- `list(absolute=1e-6)` → pct_tol=0.0, criterion="max_err" ✅ (unchanged)
- `list(pct=0.001)` → pct_tol=0.001, criterion="pct" ✅ (unchanged)
- `list(improvement=0.001)` → pct_tol=0.001, criterion="improvement" ✅

- [ ] **Step E3: Update Python `_harvest.py`**

Add `improvement` to criterion map:
```python
criterion_map = {"pct": 0, "max_err": 1, "mean_err": 2, "kl": 3, "chi2": 4, "improvement": 5}
```

Update `_parse_convergence` default:
```python
def _parse_convergence(conv):
    if conv is None:
        conv = {}
    explicit_impr = "improvement" in conv
    explicit_pct  = "pct" in conv
    explicit_abs  = "absolute" in conv
    pct_tol = conv.get("improvement", conv.get("pct",
              0.001 if not explicit_abs else 0.0))
    absolute_tol = conv.get("absolute", 0.0)
    if explicit_impr:
        criterion_str = "improvement"
    elif explicit_pct:
        criterion_str = "pct"
    elif not explicit_abs:
        criterion_str = "improvement"  # new default
    else:
        criterion_str = "max_err"
    criterion_str = conv.get("criterion", criterion_str)
    criterion = criterion_map.get(criterion_str, 5)  # default 5=IMPROVEMENT
    stop_when = {"any": 0, "all": 1}.get(conv.get("stop_when", "any"), 0)
    return pct_tol, absolute_tol, criterion, stop_when
```

- [ ] **Step E4: Update NEWS.md**

Add to the NEWS.md top:
```markdown
## leafblower 0.2.1 (unreleased)

### Breaking changes (convergence default)

- **Default convergence criterion changed** from `pct = 0.001` (max proportional
  weight change) to `improvement = 0.001` (relative improvement in marginal
  calibration error between checks). Mathematical convergence: stop when the
  loss function stops improving, not when weights stop moving.
  Old behavior: `convergence = list(pct = 0.001)`.
  New behavior: `convergence = list()` or `convergence = list(improvement = 0.001)`.

### New features

- New `criterion = "improvement"`: stop when `(errRp[t-1] - errRp[t]) / errRp[t-1] < tol`.
  Fires on both plateau (convergence) and oscillation (divergence). Combined with
  `best_weights`, returns minimum-error iterate when solver cannot converge.
```

- [ ] **Step E5: Regenerate Rd**
```bash
Rscript -e 'devtools::document()'
```

- [ ] **Step E6: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
pip install -e . && pytest python/leafblower/ -v 2>&1 | tail -10
```
Expected: FAIL 0. Fix any test that fails because the default criterion changed: if a test uses `harvest()` without pinned convergence and asserts on convergence-sensitive fields, add `convergence = list(improvement = 0.001)` (or keep existing pin).

- [ ] **Step E7: R CMD check**
```bash
R CMD build . && R CMD check --as-cran leafblower_*.tar.gz 2>&1 | grep -E "^(ERROR|WARNING)"
```
Expected: 0 ERROR, 0 WARNING.

- [ ] **Step E8: Commit + close**
```bash
git add R/harvest.R man/harvest.Rd python/leafblower/_harvest.py NEWS.md
git commit -m "$(cat <<'EOF'
feat(harvest.R): improvement criterion as default convergence

convergence=list() now uses criterion='improvement' (relative errRp
improvement < 0.001) instead of criterion='pct' (weight change < 0.001).
Mathematical convergence: stop when loss stops improving, not when
weights stop moving. Backward compat: list(pct=X) and list(absolute=X)
unchanged. Python criterion_map updated; NEWS.md documents the change.
EOF
)"
bd close <ticket-id>
```

---

## Final Verification

- [ ] `devtools::test()` → FAIL 0, PASS ≥ 302
- [ ] `pytest python/leafblower/` → all green
- [ ] `R CMD check --as-cran` → 0 ERROR, 0 WARNING
- [ ] `list()` default → criterion="improvement", pct_tol=0.001
- [ ] `list(absolute=1e-6)` → criterion="max_err", backward compat preserved
- [ ] `list(pct=0.001)` → criterion="pct" (weight-change), preserved
- [ ] `list(improvement=0.001)` → criterion="improvement", works explicitly
- [ ] iEPPA oscillating test: `iterations < 500` (stops early)
- [ ] All 5 beads tickets closed
- [ ] NEWS.md updated
