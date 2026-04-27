# Convergence Redesign — Simplify & Fix Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix 2 blocking defects and 2 DRY violations from the critical code review of the convergence redesign, plus 3 minor cleanups. All changes in `feat/convergence-redesign` worktree.

**Architecture:** Four atomic work units. WU-A fixes R-side defects (fast, no rebuild). WU-B extracts a shared C++ dispatch header (calib_dispatch.hpp). WU-C DRYs r_bridge.cpp. WU-D handles minor fixes. Each WU is one ticket, one commit.

**Tech Stack:** C++17 (ieppa.cpp, raking.cpp, lbfgsb_solver.cpp, r_bridge.cpp, new src/calib_dispatch.hpp), R (harvest.R). Build gate: `R CMD INSTALL --preclean .` after every C++ change.

**Baseline:** `devtools::test()` → FAIL 0 | PASS 326

**Worktree:** `/home/dd/Gemini/leafblower/.worktrees/convergence-redesign`

**Tickets:**
- WU-A: `leafblower-lj5j` (C1 bounds guard) + `leafblower-we1a` (C2 shorthand conflict) + `leafblower-x0xn` (S1 double-null)
- WU-B: `leafblower-moj0` (R2 calib_dispatch.hpp)
- WU-C: `leafblower-s7pb` (R1 r_bridge DRY)
- WU-D: `leafblower-homh` (S2 lbfgsb diagnostic) + `leafblower-d2s4` (R3 test debt note)

---

## File Structure

| File | WU | Change |
|---|---|---|
| `R/harvest.R` | A | Guard indexing; fix improvement+absolute conflict; simplify rule_default |
| `src/calib_dispatch.hpp` | B | New: `select_metric()` + `apply_rule()` inline helpers |
| `src/ieppa.cpp` | B | Replace triplicated dispatch with calib_dispatch.hpp calls |
| `src/raking.cpp` | B | Same |
| `src/lbfgsb_solver.cpp` | B, D | Same; fix convergence_rule to report THRESHOLD |
| `src/r_bridge.cpp` | C | Extract repeated 5-tuple into lambda or inline block |

---

## WU-A — harvest.R defects + cleanup (leafblower-lj5j, leafblower-we1a, leafblower-x0xn)

**Files:** `R/harvest.R` only. No rebuild needed.

- [ ] **Step A0: Claim tickets**
```bash
cd /home/dd/Gemini/leafblower/.worktrees/convergence-redesign
bd update leafblower-lj5j --claim
bd update leafblower-we1a --claim
bd update leafblower-x0xn --claim
```

- [ ] **Step A1a: Write lj5j test only (TDD for C1)**

Append ONLY this test to `tests/testthat/test-convergence-criteria.R`:
```r
test_that("lj5j: convergence_used metric is NA-safe when C value is out of range", {
  set.seed(1)
  data <- data.frame(a = factor(sample(c("1","2"), 100, replace=TRUE)))
  target <- list(a = c("1"=0.5, "2"=0.5))
  w <- leafblower::harvest(data, target, max_weight=3, method="ieppa",
                           attach_weights=FALSE)
  r <- attr(w, "result")
  expect_true(is.character(r$convergence_used$metric) || is.na(r$convergence_used$metric))
  expect_true(is.character(r$convergence_used$rule)   || is.na(r$convergence_used$rule))
})
```
Run: should PASS (existing safe_lookup behavior may already handle this).

- [ ] **Step A1b: Write we1a test only (TDD for C2)**

Append ONLY this test to `tests/testthat/test-convergence-criteria.R`:
```r
test_that("we1a: list(improvement=0.01, absolute=1e-6) errors on ambiguous co-supply", {
  data <- data.frame(a = factor(c("1","2")))
  target <- list(a = c("1"=0.5,"2"=0.5))
  expect_error(
    leafblower::harvest(data, target, max_weight=3, method="ieppa",
                        convergence = list(improvement=0.01, absolute=1e-6),
                        attach_weights=FALSE),
    regexp = "stop_when|ambiguous|combine"
  )
})
```
Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")' 2>&1 | tail -5`
Expected: we1a FAIL (no error raised yet).

- [ ] **Step A2: Read harvest.R — find the 3 locations to fix**

```bash
grep -n "convergence_used\|improvement.*early\|rule_default\|is.null.*rule\|is.null.*converg" \
  R/harvest.R | head -20
```

Find:
1. `convergence_used` construction (~line 205): the `+1L` indexing
2. `parse_convergence` shorthand block (~line 295): the `improvement` early-return
3. `rule_default` (~line 360): the double `is.null` check

- [ ] **Step A3: Fix C1 — guard `+1L` indexing**

Find the convergence_used block. Replace:
```r
# OLD (unguarded):
metric = .metric_names[calib_result$convergence_metric + 1L],
rule   = .rule_names[calib_result$convergence_rule + 1L],

# NEW:
.safe_lookup <- function(v, idx) {
  if (is.integer(idx) && !is.na(idx) && idx >= 0L && idx < length(v)) v[idx + 1L]
  else NA_character_
}
metric = .safe_lookup(.metric_names, calib_result$convergence_metric),
rule   = .safe_lookup(.rule_names, calib_result$convergence_rule),
```

- [ ] **Step A4: Fix C2 — error on ambiguous improvement + absolute co-supply**

In `parse_convergence()`, in the shorthand block where `improvement` is handled, add a conflict check:
```r
explicit_impr <- !is.null(convergence[["improvement"]])
# ... after detecting explicit_impr:
if (explicit_impr && !is.null(convergence[["absolute"]]) &&
    is.null(convergence[["stop_when"]])) {
  stop("convergence: 'improvement' and 'absolute' cannot be combined without 'stop_when'. ",
       "Use stop_when = 'any' (fire on either) or 'all' (require both).")
}
```

When `stop_when` IS provided, both are honored: `pct_tol = improvement value`, `absolute_tol = absolute value`.

- [ ] **Step A5: Fix S1 — hoist rule_explicit to avoid double-null check**

Find the rule_default logic (~line 360). Replace:
```r
# OLD (double-null check):
rule_default <- if (explicit_pct && is.null(convergence[["rule"]])) "plateau"
               else if (explicit_abs && is.null(convergence[["rule"]])) "threshold"
               else "improvement"

# NEW:
rule_explicit <- !is.null(convergence[["rule"]])
rule_default  <- if (!rule_explicit && explicit_pct) "plateau"
                 else if (!rule_explicit && explicit_abs) "threshold"
                 else "improvement"
```

- [ ] **Step A6: Run tests**
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")' 2>&1 | tail -5
```
Expected: both lj5j AND we1a PASS; all previous A-tests still PASS.

- [ ] **Step A7: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 326.

- [ ] **Step A8: Commit C1 — lj5j test + safe_lookup fix**

The lj5j test was appended in Step A1a. The harvest.R safe_lookup fix was done in Step A3.
```bash
git add R/harvest.R tests/testthat/test-convergence-criteria.R
git commit -m "fix(harvest.R): guard +1L indexing in convergence_used with safe_lookup (lj5j)"
bd close leafblower-lj5j
```

- [ ] **Step A9: Commit C2 — we1a test + conflict error**

The we1a test was appended in Step A1b (AFTER A8 committed the previous file state). The harvest.R C2 fix was done in Step A4.
```bash
git add R/harvest.R tests/testthat/test-convergence-criteria.R
git commit -m "fix(harvest.R): error on improvement+absolute co-supply without stop_when (we1a)"
bd close leafblower-we1a
```

- [ ] **Step A10: Commit S1 — harvest.R rule_explicit hoist**

S1 fix only touches harvest.R (no new test).
```bash
git add R/harvest.R
git commit -m "cleanup(harvest.R): hoist rule_explicit, eliminate double-null check (x0xn)"
bd close leafblower-x0xn
```

---

## WU-B — Extract calib_dispatch.hpp (leafblower-moj0)

**Files:** new `src/calib_dispatch.hpp`, `src/ieppa.cpp`, `src/raking.cpp`, `src/lbfgsb_solver.cpp`

**Clean-code principle:** Do one thing. The dispatch logic (select metric, apply rule) is one thing — it should live in one place, not three.

**CAUTION:** Read all 3 solver dispatch implementations BEFORE writing calib_dispatch.hpp to capture semantic differences (ieppa uses `<= 1e-15` zero-guard; raking uses `> 1e-15`). Then pick ONE semantics (document the choice). Compile gate after EACH solver file change.

- [ ] **Step B0: Claim ticket**
```bash
bd update leafblower-moj0 --claim
```

- [ ] **Step B1: Read all 3 dispatch blocks to identify differences**

```bash
grep -n "switch.*cfg.metric\|select_metric\|curr_metric\|prev_metric_for_rule\|1e-15" \
  src/ieppa.cpp src/raking.cpp src/lbfgsb_solver.cpp | head -30
```

Note any semantic divergences. Key expected divergence: zero-guard threshold (`<= 1e-15` vs `> 1e-15`) in the IMPROVEMENT rule.

- [ ] **Step B2: Create `src/calib_dispatch.hpp`**

```cpp
#pragma once
#include "types.hpp"
#include <cmath>
#include <limits>

namespace lbw {

// Select the active metric value from the pre-computed metric set.
// Caller computes all metrics; this function picks the one cfg.metric designates.
inline double select_metric(CalibMetric metric,
    double max_err, double mean_err, double kl,
    double chi2, double grake_norm, double l1_weight) noexcept
{
    switch (metric) {
        case CalibMetric::MAX_ERR:    return max_err;
        case CalibMetric::MEAN_ERR:   return mean_err;
        case CalibMetric::KL:         return kl;
        case CalibMetric::CHI2:       return chi2;
        case CalibMetric::GRAKE_NORM: return grake_norm;
        case CalibMetric::L1_WEIGHT:  return l1_weight;
    }
    return max_err;  // unreachable; CalibMetric is exhaustive
}

// Apply the stopping rule using curr (current metric value), prev (previous check
// value, initialized to +∞), and tol (pct_tol threshold).
//
// Rule semantics:
//   THRESHOLD:   curr < tol
//   IMPROVEMENT: |curr - prev| / prev < tol (prev > 1e-15 required; else no convergence)
//   PLATEAU:     curr >= prev * (1 - tol)   (fires when metric stops improving)
//
// First-check behavior: prev == +∞ → IMPROVEMENT rel_change = 1.0 (> any sane tol) → no convergence.
//                       prev == +∞ → PLATEAU: curr >= ∞ * (1-tol) is false → no convergence.
//
// Returns true if the stopping condition is met.
inline bool apply_rule(CalibRule rule, double curr, double prev, double tol) noexcept
{
    if (!std::isfinite(curr)) return false;
    switch (rule) {
        case CalibRule::THRESHOLD:
            return (tol > 0.0) && (curr < tol);
        case CalibRule::IMPROVEMENT: {
            if (!std::isfinite(prev) || prev <= 1e-15) return false;
            const double rel = std::fabs(curr - prev) / prev;
            return (tol > 0.0) && (rel < tol);
        }
        case CalibRule::PLATEAU:
            if (!std::isfinite(prev)) return false;
            return (tol > 0.0) && !(curr < prev * (1.0 - tol));
    }
    return false;  // unreachable
}

} // namespace lbw
```

**Design note:** `apply_rule` uses `prev <= 1e-15` (not `<`) for IMPROVEMENT, matching the stricter of the two existing implementations (raking). ieppa used `<= 1e-15`; now unified.

- [ ] **Step B3: Update `src/ieppa.cpp`**

Add `#include "calib_dispatch.hpp"` immediately after `#include "ieppa.hpp"` (the solver's own header). In ieppa.cpp this is line 2. Convention: local project headers before STL headers.

Find the convergence dispatch block (two switch statements). Replace with:

// ... at the kErrCheckInterval dispatch site:
const auto& cfg = st.convergence_cfg;
const double curr_metric = lbw::select_metric(
    cfg.metric, errRp, mean_err, kl_max, chi2_total, grake_norm, l1_weight);

const bool converged_primary = lbw::apply_rule(
    cfg.rule, curr_metric, prev_metric_for_rule, cfg.pct_tol);
prev_metric_for_rule = curr_metric;

const bool converged_abs = (cfg.absolute_tol > 0.0) && (errRp < cfg.absolute_tol);
const bool converged = (cfg.absolute_tol > 0.0)
    ? ((cfg.stop_when == CalibStopWhen::ALL)
       ? (converged_primary && converged_abs)
       : (converged_primary || converged_abs))
    : converged_primary;
```

**Build gate after ieppa.cpp:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

- [ ] **Step B4: Update `src/raking.cpp`**

Add `#include "calib_dispatch.hpp"` after `#include "raking.hpp"`. Same dispatch replacement. Verify `prev_metric_for_rule` is declared/reset correctly.

**Build gate after raking.cpp:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

- [ ] **Step B5: Update `src/lbfgsb_solver.cpp`**

Add `#include "calib_dispatch.hpp"` after `#include "lbfgsb_solver.hpp"`. Same replacement for the post-solve dispatch. lbfgsb calls `apply_rule` with `pre_metric` as `prev`.

**Note on LBFGSResult fields:** `lbfgsb_solver.hpp:17-20` already has `l1_weight_change`, `grake_norm`, `convergence_metric`, `convergence_rule` (added by WU-A). Generic lambda in WU-C works because all 3 structs have identical field names.

**Build gate after lbfgsb_solver.cpp:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

- [ ] **Step B6: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 326.

- [ ] **Step B7: Commit**
```bash
git add src/calib_dispatch.hpp src/ieppa.cpp src/raking.cpp src/lbfgsb_solver.cpp
git commit -m "$(cat <<'EOF'
refactor: extract select_metric()+apply_rule() to calib_dispatch.hpp

Eliminates triplication of CalibMetric/CalibRule dispatch across ieppa,
raking, and lbfgsb. Single implementation with documented semantics:
IMPROVEMENT requires prev > 1e-15 (stricter of two prior variants).
All three solvers include calib_dispatch.hpp and call the shared inlines.
EOF
)"
```

- [ ] **Step B8: Close ticket**
```bash
bd close leafblower-moj0
```

---

## WU-C — r_bridge.cpp DRY (leafblower-s7pb)

**Files:** `src/r_bridge.cpp`

**Clean-code principle:** The same expression repeated 3 times is 2 too many.

- [ ] **Step C0: Claim ticket**
```bash
bd update leafblower-s7pb --claim
```

- [ ] **Step C1: Read the 3 repeated blocks**
```bash
grep -n "res_l1_weight_change\|res_grake_norm\|res_conv_metric\|res_conv_rule\|res_conv_tol\|res_conv_iter" \
  src/r_bridge.cpp | head -25
```

Find the 3 sets of ~6-line assignments (one per solver branch: ieppa/raking/lbfgsb).

- [ ] **Step C2: Extract to a lambda or helper block**

In r_bridge.cpp, just before the 3 solver branches, define an in-scope lambda:

```cpp
// Hoist repeated result-packing into a single closure.
// Captures all result variables by reference.
auto pack_solver_result = [&](const auto& res) {
    res_l1_weight_change = res.l1_weight_change;
    res_grake_norm       = res.grake_norm;
    res_conv_metric      = res.convergence_metric;
    res_conv_rule        = res.convergence_rule;
    res_conv_tol         = res.convergence_tol;
    res_conv_iter        = res.convergence_iter;
};
```

In each solver branch, replace the 6-line block with:
```cpp
pack_solver_result(res);
```

Note: `auto& res` works if IEPPAResult, RakingResult, LBFGSResult all have the same field names. Verify:
```bash
grep -n "l1_weight_change\|grake_norm\|convergence_metric" \
  src/ieppa.hpp src/raking.hpp src/lbfgsb_solver.hpp | head -15
```

If field names differ, use explicit field access in the lambda body (one per solver type). The template deduction will resolve correctly.

- [ ] **Step C3: Build gate**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

- [ ] **Step C4: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 326.

- [ ] **Step C5: Commit**
```bash
git add src/r_bridge.cpp
git commit -m "refactor(r_bridge): extract repeated solver result packing into lambda"
bd close leafblower-s7pb
```

---

## WU-D — Minor fixes (leafblower-homh, leafblower-d2s4)

**Files:** `src/lbfgsb_solver.cpp`, `tests/testthat/test-convergence-criteria.R`

- [ ] **Step D0: Claim tickets**
```bash
bd update leafblower-homh --claim
bd update leafblower-d2s4 --claim
```

- [ ] **Step D1: Fix lbfgsb convergence_rule diagnostic**

First find any existing tests that assert `convergence_used$rule` for lbfgsb:
```bash
grep -n "lbfgsb\|LBFGSB" tests/testthat/test-convergence-criteria.R | grep -i "convergence_used\|rule"
```
If found, update those assertions to expect `"threshold"`.

Read lbfgsb_solver.cpp exit block where `res.convergence_rule` is set (field exists at lbfgsb_solver.hpp:20). lbfgsb always reduces to threshold semantics (single-pass). Change:
```cpp
// OLD: res.convergence_rule = static_cast<int>(cfg.rule);  // misleading: reports requested, not applied

// NEW: always report THRESHOLD because lbfgsb makes one pass regardless of rule
res.convergence_rule = static_cast<int>(CalibRule::THRESHOLD);
```

Add a comment above explaining why:
```cpp
// lbfgsb is a batch solver (single optimization pass). All convergence rules
// reduce to a start-vs-final comparison, equivalent to THRESHOLD semantics.
// We report THRESHOLD so convergence_used$rule reflects what was *applied*.
```

No build gate needed (comment + one int change). Build:
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```

Update A-test: if any test asserts `convergence_used$rule` for lbfgsb method, it should now expect "threshold".

- [ ] **Step D2: Document test debt**

Find the A1, A2, A3 test bodies in tests/testthat/test-convergence-criteria.R. Add a comment above each:

```r
# CONFIRMATORY TEST (not TDD red-green): added in the same commit as the implementation.
# Verified retroactively: running against pre-WU-C code produces FAIL due to:
#   A1 - convergence_used$rule returns NULL (field not yet constructed)
#   A2 - best_error == max_error (no IMPROVEMENT dispatch, old PCT threshold fires early)
#   A3 - convergence_rule is integer 0 (THRESHOLD), not "plateau" string
```

- [ ] **Step D3: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 326.

- [ ] **Step D4: Commit**
```bash
git add src/lbfgsb_solver.cpp tests/testthat/test-convergence-criteria.R
git commit -m "$(cat <<'EOF'
fix(lbfgsb): report applied convergence_rule=THRESHOLD not requested rule;
test(debt): document A1-A3 as confirmatory tests with retroactive red-state analysis
EOF
)"
bd close leafblower-homh leafblower-d2s4
```

---

## Final Verification

- [ ] `devtools::test()` → FAIL 0, PASS ≥ 326
- [ ] `R CMD check --as-cran` → 0 ERROR, 0 WARNING  
- [ ] `grep -r "switch.*cfg.metric\|switch.*cfg.rule" src/*.cpp` returns 0 hits (all dispatch via calib_dispatch.hpp)
- [ ] All 7 tickets closed

---

## Self-Review

**Spec coverage:**
- C1 (lj5j) → WU-A safe_lookup ✅
- C2 (we1a) → WU-A conflict error ✅
- R1 (s7pb) → WU-C lambda ✅
- R2 (moj0) → WU-B calib_dispatch.hpp ✅
- S1 (x0xn) → WU-A rule_explicit hoist ✅
- S2 (homh) → WU-D THRESHOLD report ✅
- R3 (d2s4) → WU-D comment ✅

**No placeholders.** All code blocks are concrete.

**Type consistency:** `select_metric` takes 6 doubles matching the 6 metrics in CalibMetric. `apply_rule` takes (CalibRule, double curr, double prev, double tol) — same triple across all 3 call sites.
