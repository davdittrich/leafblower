# iEPPA & Raking Algorithm Correctness — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development

**Goal:** Fix seven correctness bugs in ieppa.cpp, raking.cpp, and calib_dispatch.hpp that silently corrupt convergence metrics, stall detection, and log-path computations.

**Architecture:** Seven atomic fixes. B11/B12 in ieppa.cpp; B16/R2/R8 in raking.cpp; B17 in calib_dispatch.hpp; R3 in ieppa.cpp (perf). Each fix independent — one commit per task.

**Tech Stack:** C++17, R CMD INSTALL --preclean ., devtools::test()

---

**Mechanism:** Initialization ordering (B11), sentinel values (B12), snapshot timing (B16), convergence guard (B17), DRY (R2), hoisting (R3), mutual exclusion (R8)
**Forbidden:** Changing algorithm logic beyond the specific fix; bundling multiple bugs in one commit
**Audit:** For each bug: capture concrete evidence of wrong behavior before fix; verify correct behavior after

---

## Pre-Flight: Evidence Map

All line numbers verified against current source. Confidence scores per CLAUDE.md.

| Bug | File | Exact Line(s) | Confidence |
|-----|------|---------------|------------|
| B11 | src/ieppa.cpp | 401 (loop in lvl=0 block) | 97 |
| B12 | src/ieppa.cpp | 589, 606 | 97 |
| B16 | src/raking.cpp | 495 | 97 |
| B17 | src/calib_dispatch.hpp | 93 | 97 |
| R2  | src/raking.cpp | 532-572 | 95 |
| R3  | src/ieppa.cpp | 553 (in-loop) | 97 |
| R8  | src/raking.cpp | 159 (use_greedy decl) | 95 |

---

## Task 1 — B11: ieppa X_prev zero-init clobber

**Ticket:** `Task [B11: guard X_prev assignment at lvl=0/iter=0] ! [restructuring outer loop, touching homotopy budget logic]`

### Bug Mechanism (confidence: 97)

At `src/ieppa.cpp:290-291`, X_prev is correctly initialised from X_init:
```cpp
std::vector<double> X_prev(ct.M_cell);
for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X_init[c];
```

Then at line 401 — inside the `for (int lvl = 0; lvl < N_levels ...)` loop, unconditionally before `iter_in_lvl` starts — the comment says "reset X_prev at the start of each homotopy level":
```cpp
for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];
```

At lvl=0, iter_in_lvl=0, `X` is still all-zeros (it is populated by the first sweep). This wipes the X_init baseline. The l1_weight computed at the first convergence check (iter=1) compares X[c] (after one sweep) against X_prev[c]=0, yielding l1_weight ≈ W_total/n regardless of how close X is to X_init. This inflates early convergence checks and corrupts improvement/plateau baselines for the first level.

### Fix

In `src/ieppa.cpp` at line 401, guard the reset so it is skipped when entering lvl=0 for the first time:

**Before:**
```cpp
for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];
```

**After:**
```cpp
if (lvl > 0) {
    // Reset X_prev at the start of each subsequent homotopy level so that
    // pct_change measures iteration-to-iteration shift within a level.
    // At lvl=0, X_prev retains the X_init baseline set at declaration.
    for (int c = 0; c < ct.M_cell; c++) X_prev[c] = X[c];
}
```

**Unchanged components:** X_prev initialisation at lines 290-291 stays. WU-C reset of `prev_metric_for_rule` at line 403 is untouched. Inner iteration logic is untouched.

**Regressions prevented:** Multi-level homotopy (N_levels>1) still resets X_prev correctly at lvl=1,2,…. Single-level (N_levels=1) no longer corrupts the first-iteration l1_weight reading.

### TDD

**RED test** — write before touching source:
```r
# tests/testthat/test-b11-xprev-init.R
test_that("B11: l1_weight_change at iter=1 reflects X_init distance, not zeros", {
  # Construct a simple 1-margin, 2-category case with uniform design.
  # With X_init uniform and targets [0.4, 0.6], the first sweep moves X
  # by a predictable amount. Before the fix, l1_weight = W_total/n (comparing
  # post-sweep X to zeros). After the fix, l1_weight = real shift from X_init.
  result <- harvest(
    data         = data.frame(x = c("A","A","B","B","B"), w = 1),
    targets      = list(x = c(A=0.4, B=0.6)),
    algorithm    = "ieppa",
    max_iter     = 1,
    rule         = "improvement",
    verbose      = 0
  )
  # After fix: l1_weight_change measures X[iter=1] vs X_init, not vs 0.
  # l1_weight_change must be << 1.0 for a near-feasible problem.
  # Before fix: l1_weight_change ≈ 1.0 (comparing post-sweep to all-zeros).
  expect_lt(result$diagnostics$l1_weight_change[1], 0.5)
})
```

**Compile gate:** `R CMD INSTALL --preclean . 2>&1 | tail -5`

**GREEN verify:** Run test; confirm l1_weight_change < 0.5 and not ≈ 1.0.

**Commit message:**
```
fix(ieppa): guard X_prev reset at lvl=0 to preserve X_init baseline

B11: unconditional `X_prev = X` at homotopy level start clobbered the
X_init initialization when lvl=0 and X was still all-zeros. Adds `if
(lvl > 0)` guard so the lvl=0 entry retains the correct X_init baseline
for l1_weight_change and improvement/plateau convergence checks.
```

---

## Task 2 — B12: compute_margin_errRp returns 0.0 on empty weight sum

**Ticket:** `Task [B12: return infinity on W_total<=0 in errRp functions] ! [changing convergence tolerance values, touching non-empty paths]`

### Bug Mechanism (confidence: 97)

Two lambdas in `src/ieppa.cpp` share the same flaw:

`compute_margin_errRp_linear` at line 589:
```cpp
if (W_total <= 0.0) return 0.0;
```

`compute_margin_errRp_log` at line 606:
```cpp
if (W_total <= 0.0) return 0.0;
```

The greedy scheduler (when active) uses these per-margin residuals to prioritise margins. Returning `0.0` signals perfect convergence on that margin. A degenerate initialisation (all weights zero, e.g., infeasibility early in homotopy level 0) causes the greedy scheduler to de-prioritise every margin and makes convergence checking falsely report zero error — masking a completely un-calibrated state.

The correct sentinel is `+∞`: the margin is maximally far from target when we cannot evaluate it, and the scheduler should treat it with highest priority.

### Fix

**src/ieppa.cpp line 589** (`compute_margin_errRp_linear`):
```cpp
// Before:
if (W_total <= 0.0) return 0.0;
// After:
if (W_total <= 0.0) return std::numeric_limits<double>::infinity();
```

**src/ieppa.cpp line 606** (`compute_margin_errRp_log`):
```cpp
// Before:
if (W_total <= 0.0) return 0.0;
// After:
if (W_total <= 0.0) return std::numeric_limits<double>::infinity();
```

**Unchanged components:** All non-zero-weight code paths are identical. The `X_tilde.empty()` guard at line 603 already returns `infinity` for the log path when X_tilde is empty — this fix makes the W_total guard consistent with that existing sentinel.

**Regressions prevented:** Normal operation (W_total > 0) is unaffected. Greedy scheduler now correctly prioritises all margins during initialisation phases where weights are zero.

### TDD

**RED test** — write before touching source. Two complementary tests; at least one must be RED before implementation:

**Test A — indirect scheduler test (preferred):**
```r
# tests/testthat/test-b12-empty-weight-sentinel.R
test_that("B12: greedy scheduler selects non-zero margin when a structural zero exists", {
  # Set up a 2-margin problem where margin 1 has zero weight for category A
  # (structural zero, e.g. no observations in that cell) and margin 2 is normal.
  # With the bug: compute_margin_errRp_linear returns 0.0 for the zero-weight
  # margin → greedy treats it as perfectly converged and never prioritises it.
  # With the fix: returns Inf → greedy correctly prioritises the structural-zero
  # margin (maximum residual) above the normal margin.
  #
  # Verifiable by: the algorithm must iterate at least 2 full sweeps
  # (it cannot exit at iter=1 claiming convergence when a margin has Inf errRp).
  result <- harvest(
    data    = data.frame(
      x = c("A", "B", "B", "B", "B"),   # category A has only 1 observation
      y = c("P", "P", "Q", "Q", "Q"),
      w = c(1,   1,   1,   1,   1)
    ),
    targets   = list(
      x = c(A=0.4, B=0.6),
      y = c(P=0.4, Q=0.6)
    ),
    algorithm = "ieppa",
    scheduler = "greedy",
    max_iter  = 10
  )
  # If spurious convergence fires at iter=1, iterations will be 1.
  # After fix, at least 2 iterations must occur on this non-trivial problem.
  expect_gt(result$iterations, 1L)
  # Additionally, any reported error must be finite (not NaN from 0.0 sentinel confusion)
  expect_true(is.finite(result$max_error))
})
```

**Test B — direct unit test for the sentinel value:**
```r
test_that("B12: compute_margin_errRp_linear returns Inf when W_total=0", {
  # This test directly exercises the guard. Because the lambda is internal,
  # we drive it through a harvest() call designed to hit W_total=0 on first
  # entry: a 1-margin, 1-category problem with all-zero initial weights.
  # After fix, the function must return Inf (not 0.0).
  #
  # If compute_margin_errRp_linear is accessible via Rcpp test harness:
  #   expect_equal(compute_margin_errRp_linear_test(W_total=0), Inf)
  # Otherwise, guard via observable: algorithm must NOT claim errRp=0 before sweep.
  result <- tryCatch(
    harvest(
      data    = data.frame(x = c("A","B"), w = c(0, 0)),  # zero weights
      targets = list(x = c(A=0.5, B=0.5)),
      algorithm = "ieppa",
      scheduler = "greedy",
      max_iter  = 5
    ),
    error = function(e) NULL
  )
  # If the solver handles zero-weight input at all, it must not report
  # max_error=0.0 (which would indicate false perfect convergence from the 0.0 sentinel).
  if (!is.null(result)) {
    expect_false(isTRUE(all.equal(result$max_error, 0.0, tolerance=1e-15)))
  }
})
```

**Compile gate:** `R CMD INSTALL --preclean . 2>&1 | tail -5`

**Commit message:**
```
fix(ieppa): return infinity from errRp lambdas when W_total <= 0

B12: both compute_margin_errRp_linear and compute_margin_errRp_log
returned 0.0 when W_total<=0, signalling false convergence to the
greedy scheduler on empty-weight initialisation. Replace with
std::numeric_limits<double>::infinity() to match the existing
X_tilde.empty() sentinel and drive correct greedy prioritisation.
```

---

## Task 3 — B16: SQUAREM stale X_prev_sq after fallback

**Ticket:** `Task [B16: always update X_prev_sq after SQUAREM step] ! [changing SQUAREM step geometry, touching alpha/beta computation]`

### Bug Mechanism (confidence: 97)

In `src/raking.cpp` at line 495, inside the SQUAREM accelerated loop:
```cpp
if (!fell_back) X_prev_sq = X;
```

The comment at line 454 says: "Skip snapshot update on fell_back (prevents wchange=0 spurious stall)."

However this reasoning is backward. After a fallback, `X` holds the plain gradient step (not the super-step), so it *is* the correct new iterate. By skipping `X_prev_sq = X` on fallback, `X_prev_sq` remains two iterations stale. On the *next* iteration, the weight-change stall detector computes:
```cpp
wchange += std::fabs(X[c] - X_prev_sq[c]) * inv_n_per_cell[c];
```
...against a snapshot that is two steps old, producing an artificially large `wchange` that suppresses stall detection. Conversely, on the iteration *after* a fallback that happened to move X very little, `wchange` can be near-zero against the two-steps-old snapshot, triggering a spurious stall. The correct behaviour: always snapshot X after any accepted step.

### Fix

**src/raking.cpp line 495**:

**Before:**
```cpp
if (!fell_back) X_prev_sq = X;
```

**After:**
```cpp
// Always update X_prev_sq to the current accepted iterate (plain step on
// fallback, super-step otherwise). The two-steps-stale snapshot on fallback
// corrupts wchange and triggers both false stalls and false non-stalls.
X_prev_sq = X;
```

Remove or update the now-incorrect comment at line 454 ("Skip snapshot update on fell_back").

**Unchanged components:** `fell_back` flag logic, super-step computation, all `errRp` tracking, the stall counter `n_no_improve` accumulation logic.

**Regressions prevented:** wchange now always reflects one-step movement, making the stall window `kMaxNoImprove` function as designed.

### TDD

**RED test:**
```r
# tests/testthat/test-b16-squarem-xprev-staleness.R
test_that("B16: SQUAREM stall detection does not fire spuriously after fallback", {
  # Construct a well-conditioned calibration problem where SQUAREM should
  # converge smoothly without stalling. With the bug, a fallback step
  # can cause the stall detector to misfire.
  result <- harvest(
    data      = data.frame(x = rep(c("A","B"), 50), w = 1),
    targets   = list(x = c(A=0.5, B=0.5)),
    algorithm = "raking",
    accelerate = TRUE,
    max_iter  = 200
  )
  # A uniform problem with target=empirical proportion converges in <5 iters.
  # A spurious stall (B16) would cause status=RK_ERR_STALL at iter<10.
  expect_false(grepl("stall", result$status, ignore.case=TRUE))
  expect_lt(result$iterations, 20)
})
```

**Compile gate:** `R CMD INSTALL --preclean . 2>&1 | tail -5`

**Commit message:**
```
fix(raking): always update X_prev_sq after SQUAREM step including fallback

B16: skipping X_prev_sq update on fallback left the weight-change stall
snapshot two iterations stale, corrupting wchange computation. The plain
gradient step on fallback IS the accepted iterate and must be snapshotted.
Removes erroneous "prevents wchange=0 spurious stall" comment.
```

---

## Task 4 — B17: PLATEAU false convergence when prev=0

**Ticket:** `Task [B17: guard PLATEAU rule against prev=0 in calib_dispatch.hpp] ! [changing IMPROVEMENT rule logic, touching THRESHOLD branch]`

### Bug Mechanism (confidence: 97)

In `src/calib_dispatch.hpp` at line 93:
```cpp
case CalibRule::PLATEAU:
    // Converge when curr did NOT drop by at least tol fraction vs prev.
    // Skip on first check (prev==inf); once prev is finite always evaluate.
    if (std::isfinite(prev)) {
        converged = !(curr < prev * (1.0 - tol));
    }
    break;
```

`prev` is initialised to `std::numeric_limits<double>::infinity()` (confirmed at ieppa.cpp line 295 and raking.cpp equivalent). The comment says "skip on first check (prev==inf)". However, `prev` is unconditionally written to `curr` at line 99 after every call:
```cpp
prev = curr;
```

If on any call `curr == 0.0` (e.g., metric reaches machine zero), then `prev` is set to `0.0`. On the *next* call with any finite `curr > 0`, the guard `std::isfinite(prev)` fires (0.0 is finite), and the test `!(curr < 0.0 * (1 - tol))` = `!(curr < 0.0)` = `!false` = `true`, triggering false convergence regardless of the actual metric value.

### Fix

**src/calib_dispatch.hpp line 93**:

**Before:**
```cpp
if (std::isfinite(prev)) {
    converged = !(curr < prev * (1.0 - tol));
}
```

**After:**
```cpp
if (std::isfinite(prev) && prev > 0.0) {
    converged = !(curr < prev * (1.0 - tol));
}
```

**Unchanged components:** THRESHOLD and IMPROVEMENT branches. The `prev = curr` assignment at line 99. All callers (`check_convergence`).

**Regressions prevented:** PLATEAU now correctly skips evaluation when prev=0 (metric was previously at machine-zero), requiring a genuine non-zero baseline before firing.

### TDD

**RED test** — write before touching source:
```r
# tests/testthat/test-b17-plateau-prev-zero.R
test_that("B17: PLATEAU rule does not converge when prev=0", {
  # Indirect test: run ieppa on an already-feasible problem (targets match
  # empirical proportions exactly). At iter=1, errRp=0.0, so prev becomes 0.
  # At iter=2, if errRp is any positive value (numerical noise), PLATEAU
  # must NOT fire. With the bug, !(any_positive < 0*(1-tol)) = TRUE → converges.
  result <- harvest(
    data      = data.frame(x = c("A","A","B","B","B"), w = c(2,2,3,3,3)),
    targets   = list(x = c(A=0.4, B=0.6)),  # exact match to empirical
    algorithm = "ieppa",
    rule      = "plateau",
    pct_tol   = 0.01,
    max_iter  = 50
  )
  # The problem converges at iteration 1 via THRESHOLD (errRp=0),
  # not at iteration 2 via a false PLATEAU. The distinction:
  # - false PLATEAU: converges at iter=2 with errRp>0 reported
  # - correct: converges at iter=1 with errRp=0
  expect_equal(result$max_error, 0.0, tolerance=1e-10)
})

# Direct unit test for apply_rule (if exposed via Rcpp test harness):
# apply_rule(CalibRule::PLATEAU, curr=0.5, prev=0.0, tol=0.01) must be FALSE.
```

**Compile gate:** `R CMD INSTALL --preclean . 2>&1 | tail -5`

**Commit message:**
```
fix(calib_dispatch): PLATEAU rule must not fire when prev=0

B17: std::isfinite(0.0) is true, so when the metric previously reached
zero, any subsequent positive metric value triggered false PLATEAU
convergence via !(curr < 0*(1-tol)). Add prev>0.0 guard to require a
positive baseline before evaluating fractional improvement.
```

---

## Task 5 — R2: Remove duplicate metric aggregation in raking flat loop

**Ticket:** `Task [R2: replace inline metric aggregation with compute_cell_metrics call] ! [changing compute_cell_metrics signature, altering metric semantics]`

### Bug Mechanism (confidence: 95)

`src/raking.cpp` lines 532-572 re-implement bucket aggregation, kl_max, chi2_total, grake_norm, and mean_err inside the flat loop's `need_extra` block. This is a ~40-line verbatim duplicate of the logic in `lbw::compute_cell_metrics`. Any divergence between the two implementations (tolerances, edge cases) silently produces inconsistent metric values depending on which solver path is active.

**Exact signature of `lbw::compute_cell_metrics` (verified from `src/calib_dispatch.hpp` lines 146-149):**
```cpp
inline CellMetrics compute_cell_metrics(
    const CalibState& st, const CellTable& ct,
    const std::vector<double>& X, double W,
    std::vector<double>& bucket) noexcept
```
- `W` is the pre-computed weight sum (pass `W_total` from the flat loop)
- `bucket` is a pre-allocated scratch vector of size `>= max(st.cat_counts[k])` (already available in the flat loop scope)
- Returns `lbw::CellMetrics` with fields: `errRp`, `mean_err`, `kl`, `chi2`, `grake_norm`, `l1`

### Fix

Replace lines 532-572 in `src/raking.cpp` with a single call to `lbw::compute_cell_metrics`:

```cpp
// Before (lines 532-572): ~40 lines of inline bucket aggregation
if (need_extra) {
    double W_tot2 = 0.0;
    for (int c = 0; c < ct.M_cell; c++) W_tot2 += X[c];
    // ... 35 more lines ...
}

// After — exact replacement call using verified signature:
if (need_extra) {
    const lbw::CellMetrics m = lbw::compute_cell_metrics(st, ct, X, W_total, bucket);
    mean_err   = m.mean_err;
    kl_max     = m.kl;
    chi2_total = m.chi2;
    grake_norm = m.grake_norm;
    // m.errRp is the per-iteration max_err; update if the flat loop also tracks it:
    // max_err = m.errRp;  — verify whether max_err is a local in this scope
}
```

**Pre-edit verification:** Read lines 532-572 of `src/raking.cpp` in full before editing to confirm the local variable names (`W_total`, `kl_max`, `chi2_total`, `grake_norm`, `mean_err`) match the assignments above. If names differ, adjust accordingly.

**Verification requirement:** Run `devtools::test()` before and after. All metric values in test output must be bit-identical (or within floating-point associativity tolerance if accumulation order differs).

**Compile gate:** `R CMD INSTALL --preclean . 2>&1 | tail -5`

**Commit message:**
```
refactor(raking): replace duplicate metric aggregation with compute_cell_metrics

R2: flat loop re-implemented bucket aggregation, kl, chi2, grake_norm in
~40 lines — exact duplicate of lbw::compute_cell_metrics. Replace with a
single function call to eliminate the divergence risk. No behavioral change.
```

---

## Task 6 — R3: Hoist log_threshold out of the hot cell loop

**Ticket:** `Task [R3: hoist log_threshold constant above homotopy loop] ! [changing threshold value, touching non-log paths]`

### Bug Mechanism (confidence: 97)

In `src/ieppa.cpp` at line 553, inside the per-(k,j) inner loop of the log-path sweep:
```cpp
double log_threshold = std::log(kEmptyBucketThreshold * ct.W_input);
```

`kEmptyBucketThreshold` is a compile-time constant (line 66: `constexpr double kEmptyBucketThreshold = 1e-15;`) and `ct.W_input` is read-only after construction. This recomputes an identical `std::log(...)` call on every (k,j) pair in every iteration of every homotopy level. Pure waste — no behavioral change from hoisting.

### Fix

**Step 1:** Locate the homotopy loop start (approximately line 357) and the log-path section entry. The hoist point is just before the `for (int k = 0; k < st.K; k++)` log-path loop, after `ct.W_input` is known to be stable.

**Step 2:** Add before the log-path loop:
```cpp
// Hoisted: kEmptyBucketThreshold and ct.W_input are loop-invariant.
const double log_empty_threshold = std::log(kEmptyBucketThreshold * ct.W_input);
```

**Step 3:** Replace line 553:
```cpp
// Before:
double log_threshold = std::log(kEmptyBucketThreshold * ct.W_input);
if (!std::isfinite(log_S_kj) || log_S_kj < log_threshold) {
// After:
if (!std::isfinite(log_S_kj) || log_S_kj < log_empty_threshold) {
```

**Step 4:** If the log-path lambda captures `log_threshold` by reference, update the capture list to capture `log_empty_threshold` instead.

**Unchanged components:** `kEmptyBucketThreshold` constant, `ct.W_input`, all threshold semantics. This is pure constant-folding — identical results.

**Compile gate:** `R CMD INSTALL --preclean . 2>&1 | tail -5`
**Test:** `devtools::test()` — no behavioral change means all existing tests must pass unchanged.

**Commit message:**
```
perf(ieppa): hoist log_empty_threshold out of log-path cell loop

R3: std::log(kEmptyBucketThreshold * ct.W_input) was recomputed on every
(k,j) pair in every iteration despite being loop-invariant. Hoist to a
const before the margin loop. Pure performance — no behavioral change.
```

---

## Task 7 — R8: Enforce greedy+SQUAREM mutual exclusion

**Ticket:** `Task [R8: disable greedy scheduler when SQUAREM accelerate=TRUE] ! [changing SQUAREM step logic, touching non-greedy paths]`

### Bug Mechanism (confidence: 95)

At `src/raking.cpp` line 159, `use_greedy` is declared as `const bool`:
```cpp
const bool use_greedy = (st.scheduler.mode == SchedulerMode::GREEDY);
```

`F_eval` is a `[&]`-capturing lambda defined at line 251. It captures `use_greedy` by reference from the enclosing scope. Inside F_eval, `use_greedy` appears at **two sites**:
- Line 255: `if (use_greedy)` — controls the greedy sort
- Line 260: `const int k = use_greedy ? margin_order[ki] : ki;` — controls margin selection

At line 354, inside the SQUAREM branch:
```cpp
(void)errRp_w1;  // errRp_k updated inside F_eval but not consumed (use_greedy=false when accelerate=true)
```

This comment documents that the two modes are mutually exclusive — but the exclusion is enforced only by documentation, not by code. If a caller sets both `accelerate=TRUE` and `scheduler="greedy"`, `use_greedy=true` is captured by `F_eval` and used inside SQUAREM sub-steps (w1, w2), producing undefined interaction: the greedy priority list was designed for sequential raking, not for the SQUAREM fixed-point iteration which requires a consistent sweep order.

### Fix

**Step 1:** After line 159, introduce a non-const shadow variable:

```cpp
const bool use_greedy = (st.scheduler.mode == SchedulerMode::GREEDY);

// R8: greedy scheduler and SQUAREM acceleration are mutually exclusive.
// Greedy priority reordering is designed for sequential raking; using it
// inside SQUAREM sub-steps (w1, w2) breaks the fixed-point geometry.
// Silently demote to round_robin and log the override.
bool use_greedy_effective = use_greedy;
if (st.accelerate && use_greedy_effective) {
    use_greedy_effective = false;
    st.log("[raking] greedy scheduler disabled under SQUAREM acceleration; using round_robin");
}
```

**Step 2:** Update the `F_eval` lambda's two `use_greedy` use sites to `use_greedy_effective`. Because `F_eval` uses `[&]` capture, the reference will automatically bind to `use_greedy_effective` once you rename the variable at those two sites:

```cpp
// Before (line 255):
if (use_greedy)
// After:
if (use_greedy_effective)

// Before (line 260):
const int k = use_greedy ? margin_order[ki] : ki;
// After:
const int k = use_greedy_effective ? margin_order[ki] : ki;
```

The `[&]` capture list does **not** need to be modified — `F_eval` already captures all enclosing variables by reference; the rename at the use sites is sufficient.

**Step 3:** In the flat (non-SQUAREM) loop branch, replace its own `use_greedy` references with `use_greedy_effective` for consistency (the value is identical on the non-SQUAREM path, so no behavioral change there).

**All use_greedy occurrences in raking.cpp (verified by grep):**
- Line 159: declaration (`const bool use_greedy = ...`) — keep, add shadow after
- Line 255: inside F_eval — rename to `use_greedy_effective`
- Line 260: inside F_eval — rename to `use_greedy_effective`
- Line 354: comment only — update text to reflect the guard

**Unchanged components:** Non-SQUAREM path with greedy scheduler. SQUAREM sub-step geometry. `F_eval` `[&]` capture list (no change needed).

**Regressions prevented:** `accelerate=TRUE, scheduler="round_robin"` path is unchanged. `accelerate=FALSE, scheduler="greedy"` path is unchanged.

### TDD

**RED test:**
```r
# tests/testthat/test-r8-greedy-squarem-exclusion.R
test_that("R8: accelerate=TRUE with greedy scheduler runs without error", {
  # Before fix: undefined behavior — greedy reordering inside SQUAREM sub-steps.
  # After fix: silently uses round_robin, completes normally.
  expect_no_error({
    result <- harvest(
      data      = data.frame(x = rep(c("A","B"), 25), w = 1),
      targets   = list(x = c(A=0.5, B=0.5)),
      algorithm = "raking",
      accelerate = TRUE,
      scheduler  = "greedy",
      max_iter  = 50
    )
  })
  # Confirm it actually ran (not zero iterations)
  expect_gt(result$iterations, 0)
})
```

**Compile gate:** `R CMD INSTALL --preclean . 2>&1 | tail -5`

**Commit message:**
```
fix(raking): enforce mutual exclusion between greedy scheduler and SQUAREM

R8: greedy margin reordering is incompatible with SQUAREM sub-step geometry.
Prior code relied on undocumented caller convention; add explicit runtime
guard that demotes to round_robin when accelerate=TRUE and logs the override.
```

---

## Execution Checklist

Each task is fully independent. Execute sequentially or in parallel per superpowers:subagent-driven-development.

**Per-task gate:** After every fix:
1. `R CMD INSTALL --preclean . 2>&1 | tail -5` — must show `* DONE (leafblower)`
2. `devtools::test()` — must show 0 failures
3. Commit with the message specified in the task

**Global order recommendation:** B11 → B17 → B12 → B16 → R8 → R2 → R3
(correctness bugs before DRY/perf; convergence machinery before scheduler fixes)

**No bundling:** Each bug gets exactly one commit. Do not combine even if two edits touch the same file in the same session.
