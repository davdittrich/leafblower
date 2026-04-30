# Status Code Propagation Fixes — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development

**Goal:** Ensure every solver emits correct RK_ERR_BUDGET(4)/RK_ERR_STALL(5); eliminate dead is_infeasible in raking; fix lbfgsb convergence_rule.

**Architecture:** Three independent fixes in lbfgsb_solver.cpp and raking.cpp. No inter-dependencies except harvest.R doc update comes last.

**Tech Stack:** C++17, R, status codes RK_OK=0 NOCONV=1 INFEAS=2 BADARG=3 BUDGET=4 STALL=5

---

**Mechanism:** Exit-path classification (iter==max_iter→BUDGET, gradient-plateau→STALL, water_fill infeasible→INFEAS)
**Forbidden:** Emitting RK_ERR_NOCONV(1) from any solver — it is a legacy alias; harvest.R already warns "legacy" when status==1L
**Audit:** RED tests asserting specific status codes before fix; GREEN after

---

## Source Ground Truth (read before editing)

### Status codes — `src/leafblower.h` lines 32–37
```c
#define RK_OK           0  /* Converged */
#define RK_ERR_NOCONV   1  /* Legacy alias — no longer emitted by new solvers */
#define RK_ERR_INFEAS   2  /* Infeasible: empty cell with positive target */
#define RK_ERR_BADARG   3  /* Invalid argument */
#define RK_ERR_BUDGET   4  /* Budget exhausted while loss still decreasing */
#define RK_ERR_STALL    5  /* Loss function plateau */
```

### CalibRule enum — `src/types.hpp` lines 47–50
```cpp
enum class CalibRule : int {
    THRESHOLD   = 0,
    IMPROVEMENT = 1,
    PLATEAU     = 2
};
```

### Current lbfgsb defect — `src/lbfgsb_solver.cpp` lines 295–312
`compute_final_weights_and_error` sets `res.status = converged ? RK_OK : RK_ERR_NOCONV` (line 305) with no BUDGET/STALL distinction. It also hardcodes `res.convergence_rule = static_cast<int>(lbw::CalibRule::THRESHOLD)` (line 309) regardless of `cfg.rule`.

The inner loop `lbfgsb_solve_inner` (lines 602–664):
- breaks early when `gn < st.tol_abs` (line 606) — gradient-norm convergence
- otherwise runs to `iter == max_iter - 1` and falls through to `compute_final_weights_and_error`
- `final_iter` is set to `iter + 1` each iteration; after loop exhaustion `final_iter == max_iter`

### Current raking defect — `src/raking.cpp`
- `is_infeasible` initialized `false` at line 117 (correct)
- `is_infeasible` is set `true` inside `water_fill_cat` lambda (lines 183, 197, 287) on bound-infeasibility events
- After loop: code comments say "INFEAS only overrides on stall (post-loop check below)" — but the post-loop override **does not exist** (lines 629–665 have no `is_infeasible` check)
- If loop stalls: `res.status = RK_ERR_STALL` is set (line 620), which is correct but INFEAS is silently lost

---

## Task 1 — B2: lbfgsb BUDGET/STALL status codes

**File:** `src/lbfgsb_solver.cpp`
**Scope:** `compute_final_weights_and_error` (called from `lbfgsb_solve_inner`) and the main loop in `lbfgsb_solve_inner`.

### Mechanism

`compute_final_weights_and_error` receives `iterations` (= `final_iter` from the loop). The loop runs `for (int iter = 0; iter < max_iter; iter++)` and breaks when `gn < st.tol_abs`. Two distinct exits:

1. **Gradient-norm early break** (`gn < st.tol_abs`): `final_iter < max_iter`. This is a genuine gradient plateau — emit STALL only if not threshold-converged.
2. **Loop exhausted** (`final_iter == max_iter` and `gn >= st.tol_abs` on last iteration): emit BUDGET.
3. **Threshold convergence** (`converged == true`): emit RK_OK regardless.

STALL detection for lbfgsb uses the gradient-norm exit. A separate per-iteration gradient plateau tracker (tracking `prev_gn` over 5 iters) is the more principled approach requested in the plan header. However, the simpler and provably correct approach is: if `final_iter < max_iter` and not `converged` → STALL (gradient hit tol_abs); if `final_iter == max_iter` and not `converged` → BUDGET.

`compute_final_weights_and_error` must receive `max_iter` to distinguish these cases. Current signature: `compute_final_weights_and_error(st, fn, d, off, u, final_iter)`. Add `max_iter` as a parameter.

### Step 1.1 — Write RED test (file: `tests/testthat/test-lbfgsb.R`, append)

```r
test_that("lbfgsb emits BUDGET(4) when max_iterations exhausted before convergence", {
  # 2-category problem with target far from sample proportions.
  # max_iterations=3 is far below what convergence requires (~50+).
  # Expected: status == 4L (BUDGET). Currently returns status == 1L (NOCONV).
  set.seed(42)
  n  <- 5000L
  df <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE,
                                     prob = c(0.9, 0.1))))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  res_raw <- suppressWarnings(
    harvest(df, tgt, method = "lbfgsb",
            max_iterations = 3L,
            max_weight = 10,
            convergence = list(pct = 1e-4),
            attach_weights = FALSE)
  )
  result_attr <- attr(
    harvest(df, tgt, method = "lbfgsb",
            max_iterations = 3L,
            max_weight = 10,
            convergence = list(pct = 1e-4)),
    "result"
  )
  # Before fix: result_attr$status == 1L (NOCONV)
  # After fix:  result_attr$status == 4L (BUDGET)
  expect_equal(result_attr$status, 4L,
               info = paste("got status:", result_attr$status,
                            "expected BUDGET=4"))
})

test_that("lbfgsb emits RK_OK(0) or STALL(5) on already-converged input, never NOCONV", {
  # Well-conditioned problem run to convergence (max_iterations=1000).
  # Should return 0 (converged) or 5 (stall at optimum), never 1 (NOCONV).
  set.seed(7)
  n   <- 10000L
  df  <- data.frame(
    age = factor(sample(c("Y","M","O"), n, replace=TRUE, prob=c(0.33,0.34,0.33))),
    sex = factor(sample(c("M","F"),     n, replace=TRUE, prob=c(0.50,0.50)))
  )
  tgt <- list(age = c(Y=0.33, M=0.34, O=0.33),
              sex = c(M=0.50, F=0.50))
  result <- suppressWarnings(
    harvest(df, tgt, method = "lbfgsb",
            max_iterations = 1000L,
            convergence = list(pct = 1e-4))
  )
  s <- attr(result, "result")$status
  expect_true(s %in% c(0L, 5L),
              info = paste("expected status 0 or 5, got:", s))
})
```

### Step 1.2 — Run RED, confirm failure

```bash
cd /home/dd/Gemini/leafblower
Rscript -e "devtools::test(filter='lbfgsb')" 2>&1 | grep -E 'FAIL|PASS|Error|status'
```

Expected output includes: `Failure ... expected status: 4, got: 1`

### Step 1.3 — Implement fix in `src/lbfgsb_solver.cpp`

**Change 1:** Update `compute_final_weights_and_error` signature to accept `max_iter`.

Find the function definition (around line 200; exact line confirmed by reading):
```cpp
static LBFGSResult compute_final_weights_and_error(CalibState& st,
                                                    const LinkFn& fn,
                                                    const std::vector<double>& d,
                                                    const std::vector<int>& off,
                                                    const std::vector<double>& u,
                                                    int iterations) {
```
Change to:
```cpp
static LBFGSResult compute_final_weights_and_error(CalibState& st,
                                                    const LinkFn& fn,
                                                    const std::vector<double>& d,
                                                    const std::vector<int>& off,
                                                    const std::vector<double>& u,
                                                    int iterations,
                                                    int max_iter) {
```

**Change 2:** Replace the status-setting block (current lines 305–309) inside `compute_final_weights_and_error`:

Current:
```cpp
        res.status             = converged ? RK_OK : RK_ERR_NOCONV;
        res.convergence_metric = static_cast<int>(cfg.metric);
        // lbfgsb is a batch solver: single optimization pass regardless of rule requested.
        // Report THRESHOLD as the applied rule for accurate diagnostic output.
        res.convergence_rule   = static_cast<int>(lbw::CalibRule::THRESHOLD);
```

Replacement:
```cpp
        if (converged) {
            res.status = RK_OK;
        } else if (iterations >= max_iter) {
            // Loop exhausted with gradient still above tol_abs — budget spent.
            res.status = RK_ERR_BUDGET;
        } else {
            // Early break via gn < st.tol_abs — gradient plateau at constrained optimum.
            res.status = RK_ERR_STALL;
        }
        res.convergence_metric = static_cast<int>(cfg.metric);
        // Report THRESHOLD as the applied rule: lbfgsb is a batch solver;
        // IMPROVEMENT/PLATEAU rules have no per-iteration baseline here.
        res.convergence_rule   = static_cast<int>(lbw::CalibRule::THRESHOLD);
```

**Change 3:** Update the call site in `lbfgsb_solve_inner` (around line 664):

Current:
```cpp
    return compute_final_weights_and_error(st, fn, d, off, u, final_iter);
```
Replacement:
```cpp
    return compute_final_weights_and_error(st, fn, d, off, u, final_iter, max_iter);
```

### Step 1.4 — Compile

```bash
cd /home/dd/Gemini/leafblower
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Must show: `* DONE (leafblower)`. Fix any compile errors before proceeding.

### Step 1.5 — Run GREEN test

```bash
Rscript -e "devtools::test(filter='lbfgsb')" 2>&1 | grep -E 'FAIL|PASS|Error|status'
```

Expected: all lbfgsb tests PASS including the two new ones.

### Step 1.6 — Full test suite regression check

```bash
Rscript -e "devtools::test()" 2>&1 | tail -20
```

No new failures. If any existing test checks `status == 1L` for lbfgsb, that test is itself wrong (NOCONV is forbidden); update it to the correct code.

---

## Task 2 — B9: raking is_infeasible → RK_ERR_INFEAS

**File:** `src/raking.cpp`
**Scope:** Post-loop block after `// end else flat loop` comment (around line 627).

### Mechanism

`is_infeasible` is correctly tracked during `water_fill_cat` calls. The post-loop block normalizes weights but never consults `is_infeasible`. The fix: if the loop exited with STALL or BUDGET **and** `is_infeasible` is true, override status to INFEAS. Do not override RK_OK — the convergence comment at line 593 explains why: transient infeasibility during convergence must not corrupt a successful result.

### Step 2.1 — Write RED test (file: `tests/testthat/test-raking.R`, append)

```r
test_that("raking returns INFEAS(2) when bounds make targets structurally unreachable", {
  # Construct a provably infeasible problem:
  # n=10, single margin with 2 categories, min_weight=0.8.
  # Category "a" has 2 obs → cell mass L_cell = 0.8*2 = 1.6 (cell aggregate lower bound).
  # Category "b" has 8 obs → cell mass L_cell = 0.8*8 = 6.4.
  # Target for "a" = 0.8 * n = 8.0. But max achievable in "a" = U_cell = max_weight*2.
  # Set max_weight=1.5: U_cell_a = 1.5*2 = 3.0 < 8.0 → structurally infeasible.
  # Without fix: returns STALL(5). With fix: returns INFEAS(2).
  n  <- 20L
  df <- data.frame(x = factor(c(rep("a", 2L), rep("b", 18L))))
  # Target: 90% in "a" — but only 2/20 obs are "a", max_weight=1.5 caps cell at 3.
  # 3 < 0.9 * 20 = 18 → cannot reach target mass for "a".
  tgt <- list(x = c(a = 0.9, b = 0.1))
  # harvest() stops with an error on status==2L (infeasible hard stop in harvest.R).
  # We catch that error and verify the status propagated correctly.
  expect_error(
    suppressWarnings(
      harvest(df, tgt, method = "raking",
              min_weight = 0.5, max_weight = 1.5,
              max_iterations = 200L,
              convergence = list(pct = 1e-4))
    ),
    regexp = "infeasible",
    info   = "expected infeasible hard-stop from harvest.R; got something else"
  )
})

test_that("raking returns INFEAS(2) — status attribute before harvest.R error intercept", {
  # Access raw C result before harvest.R raises the error.
  # Use tryCatch to capture the condition and inspect its structure.
  # harvest.R raises stop("leafblower: infeasible...") after setting convergence_reason="infeasible".
  n  <- 20L
  df <- data.frame(x = factor(c(rep("a", 2L), rep("b", 18L))))
  tgt <- list(x = c(a = 0.9, b = 0.1))
  err <- tryCatch(
    harvest(df, tgt, method = "raking",
            min_weight = 0.5, max_weight = 1.5,
            max_iterations = 200L,
            convergence = list(pct = 1e-4)),
    error = function(e) e
  )
  expect_true(inherits(err, "error"))
  expect_match(conditionMessage(err), "infeasible",
               info = paste("message was:", conditionMessage(err)))
})
```

### Step 2.2 — Run RED, confirm failure

```bash
Rscript -e "devtools::test(filter='raking')" 2>&1 | grep -E 'FAIL|PASS|Error|infeasible'
```

Expected: new infeasible test fails (no error thrown, or wrong status).

### Step 2.3 — Implement fix in `src/raking.cpp`

Locate the post-loop normalization block (lines ~629–637). Insert after `}  // end else flat loop` and before the post-loop normalization:

Current (around line 629):
```cpp
    // For stalled iterations: return STALL (status=5) + best weights rather than
    // hard-erroring — caller can use the best achievable calibration.

    // Post-loop: normalize sum to n (water-filling already enforces bounds)
    {
        double s_post = 0.0;
```

Insert after `// hard-erroring — caller can use the best achievable calibration.`:
```cpp
    // Post-loop infeasibility promotion: if loop exited without converging AND
    // water_fill_cat detected a structural bound conflict (is_infeasible==true),
    // escalate STALL/BUDGET to INFEAS. Never override RK_OK — transient
    // infeasibility during a successful convergence run must not corrupt status.
    if (is_infeasible && res.status != RK_OK && res.status != RK_ERR_BADARG) {
        res.status = RK_ERR_INFEAS;
    }
```

### Step 2.4 — Compile

```bash
cd /home/dd/Gemini/leafblower
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Must show: `* DONE (leafblower)`.

### Step 2.5 — Run GREEN test

```bash
Rscript -e "devtools::test(filter='raking')" 2>&1 | grep -E 'FAIL|PASS|Error|infeasible'
```

Expected: all raking tests PASS including the two new ones.

### Step 2.6 — Full test suite regression check

```bash
Rscript -e "devtools::test()" 2>&1 | tail -20
```

Verify no regressions. The existing stall test (`descent monitor aborts early`) must still pass — that test uses an imbalanced sample that stalls but is NOT structurally infeasible, so `is_infeasible` stays false and STALL is preserved.

---

## Task 3 — R5: lbfgsb convergence_rule reports as-requested

**File:** `src/lbfgsb_solver.cpp`
**Scope:** `compute_final_weights_and_error`, the `convergence_rule` assignment (line 309, already touched in Task 1).

### Mechanism

Currently line 309 hardcodes `lbw::CalibRule::THRESHOLD`. The fix: emit `cfg.rule` (the rule the caller requested) so that `harvest.R`'s `convergence_used$rule` faithfully reflects the caller's intent. Add a log message noting that lbfgsb reduced all rules to a threshold check internally. This is a diagnostic transparency fix only — no algorithm change.

### Step 3.1 — Write RED test (file: `tests/testthat/test-lbfgsb.R`, append)

```r
test_that("lbfgsb convergence_used$rule reflects requested rule", {
  # Request rule="improvement" (CalibRule::IMPROVEMENT = 1 → "improvement" in .rule_names).
  # Before fix: convergence_used$rule == "threshold" (hardcoded).
  # After fix:  convergence_used$rule == "improvement".
  set.seed(42)
  n   <- 5000L
  df  <- data.frame(x = factor(sample(c("a","b","c"), n, replace=TRUE)))
  tgt <- list(x = c(a=1/3, b=1/3, c=1/3))
  result <- suppressWarnings(
    harvest(df, tgt, method = "lbfgsb",
            max_weight = 5,
            convergence = list(rule = "improvement", pct = 1e-4))
  )
  cu <- attr(result, "result")$convergence_used
  expect_equal(cu$rule, "improvement",
               info = paste("expected 'improvement', got:", cu$rule))
})
```

### Step 3.2 — Run RED, confirm failure

```bash
Rscript -e "devtools::test(filter='lbfgsb')" 2>&1 | grep -E 'FAIL|PASS|convergence_used|rule'
```

Expected: `Failure ... expected 'improvement', got: 'threshold'`

### Step 3.3 — Implement fix in `src/lbfgsb_solver.cpp`

Inside `compute_final_weights_and_error`, **after** the status block from Task 1, update the `convergence_rule` line:

Current (now inside the if/else block from Task 1):
```cpp
        // Report THRESHOLD as the applied rule: lbfgsb is a batch solver;
        // IMPROVEMENT/PLATEAU rules have no per-iteration baseline here.
        res.convergence_rule   = static_cast<int>(lbw::CalibRule::THRESHOLD);
```

Replacement:
```cpp
        // Report the rule as requested by the caller. Internally lbfgsb reduces
        // IMPROVEMENT/PLATEAU to a threshold check on the batch metric (no
        // per-iteration baseline exists). The as-requested rule is more useful
        // for downstream diagnostics than a hardcoded THRESHOLD.
        if (cfg.rule != lbw::CalibRule::THRESHOLD) {
            char msg[128];
            std::snprintf(msg, sizeof(msg),
                "[lbfgsb] rule=%d applied as threshold (batch solver); "
                "reporting rule as-requested for diagnostics",
                static_cast<int>(cfg.rule));
            st.log(msg);
        }
        res.convergence_rule = static_cast<int>(cfg.rule);
```

### Step 3.4 — Compile

```bash
cd /home/dd/Gemini/leafblower
R CMD INSTALL --preclean . 2>&1 | tail -5
```

Must show: `* DONE (leafblower)`.

### Step 3.5 — Run GREEN test

```bash
Rscript -e "devtools::test(filter='lbfgsb')" 2>&1 | grep -E 'FAIL|PASS|rule'
```

Expected: all lbfgsb tests PASS.

### Step 3.6 — Full test suite regression check

```bash
Rscript -e "devtools::test()" 2>&1 | tail -20
```

No regressions. Note: `harvest.R`'s `.rule_names` vector is `c("threshold", "improvement", "plateau")` using 0-based index lookup via `.safe_lookup`; `CalibRule::IMPROVEMENT=1` → index 1 → `"improvement"`. No R changes needed.

---

## Execution Order and Dependency Notes

Tasks 1, 2, and 3 are independent — no shared state. However:
- Tasks 1 and 3 both touch `compute_final_weights_and_error` in `lbfgsb_solver.cpp`. Execute them in sequence (Task 1 first, Task 3 second) on the same file state to avoid conflicts.
- Task 2 touches only `raking.cpp`. Can be done in any order relative to Tasks 1/3.
- Recommended sequence: Task 1 → Task 3 (same file, same compile) → Task 2 → full suite.

## Compile-once shortcut for Tasks 1+3

After applying both `lbfgsb_solver.cpp` edits, compile once:
```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
Rscript -e "devtools::test(filter='lbfgsb')" 2>&1 | grep -E 'FAIL|PASS'
```

## Post-completion verification

```bash
# All tests must pass:
Rscript -e "devtools::test()" 2>&1 | grep -E 'FAIL|ERROR|OK'

# Status codes confirmed absent from lbfgsb path:
grep -n 'RK_ERR_NOCONV' /home/dd/Gemini/leafblower/src/lbfgsb_solver.cpp
# Expected: 0 matches (or only in comments/old-code remarks)

# INFEAS guard confirmed present in raking:
grep -n 'is_infeasible.*RK_ERR_INFEAS\|RK_ERR_INFEAS.*is_infeasible' \
  /home/dd/Gemini/leafblower/src/raking.cpp
# Expected: 1 match (the new post-loop guard)

# convergence_rule no longer hardcoded:
grep -n 'CalibRule::THRESHOLD' /home/dd/Gemini/leafblower/src/lbfgsb_solver.cpp
# Expected: 0 matches (condition check only, not assignment)
```

## Files Modified

| File | Task | Change |
|------|------|--------|
| `src/lbfgsb_solver.cpp` | T1, T3 | `compute_final_weights_and_error`: add `max_iter` param, replace NOCONV with BUDGET/STALL dispatch, emit as-requested `convergence_rule` |
| `src/raking.cpp` | T2 | Post-loop: promote STALL/BUDGET to INFEAS when `is_infeasible == true` |
| `tests/testthat/test-lbfgsb.R` | T1, T3 | Append 3 new tests |
| `tests/testthat/test-raking.R` | T2 | Append 2 new tests |

`harvest.R` documentation update (line 94 comment `1=max_iter hit`) should reflect `4=budget/5=stall` — update in the same commit as the C++ changes.
