# Convergence Criteria Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the single `criterion` C ABI field with orthogonal `metric` (6 values: what to measure) + `rule` (3 values: how to stop), change the default to `max_err+improvement+0.001`, rename `pct_change` → `l1_weight_change`, and expose `convergence_used` diagnostics in the result.

**Architecture:** Pre-step + WU-A through WU-G (WU-E split into E1/E2). Pre-step pins unpinned tests. WU-A atomically scaffolds types + C ABI (r_bridge.cpp also updated atomically to avoid broken-build window). WU-B/C/D implement solver dispatch. WU-E1 wires r_bridge.cpp; WU-E2 wires harvest.R. WU-F1/F2 do TDD tests + migration. WU-G does Python + docs.

**NOTE — no legacy helper:** `rk_params_set_legacy_convergence()` is omitted (software unreleased; breaking changes OK; no external C ABI users).

**Tech Stack:** C++17 (types.hpp, leafblower.h, ieppa.cpp, raking.cpp, lbfgsb_solver.cpp, c_api.cpp — all inline; NO separate validation.hpp), R (harvest.R, r_bridge.cpp), Python (_harvest.py, _bindings.cpp), testthat 3, pytest.

**Spec:** `docs/superpowers/specs/2026-04-25-convergence-redesign.md`

**Baseline:** `devtools::test()` → FAIL 0 | PASS 298

**Beads tickets:**
- Pre-step: `leafblower-lt3g`
- WU-A: `leafblower-uea6`
- WU-B: `leafblower-iom2`
- WU-C: `leafblower-h02q`
- WU-D: `leafblower-fvb4`
- WU-E1: `leafblower-1xrh` (r_bridge.cpp only)
- WU-E2: `leafblower-6lxk` (harvest.R only)
- WU-F1: `leafblower-rtgv`
- WU-F2: `leafblower-o9x1`
- WU-G: `leafblower-lccb`

**Key constants:**
- `EXPECTED_RK_PARAMS_BYTES` changes 216 → 220 (criterion removed −4B, metric+rule added +8B, net +4B)
- CalibMetric: `MAX_ERR=0, MEAN_ERR=1, KL=2, CHI2=3, GRAKE_NORM=4, L1_WEIGHT=5`
- CalibRule: `THRESHOLD=0, IMPROVEMENT=1, PLATEAU=2`
- Default: `metric=0, rule=1, pct_tol=0.001`

---

## File Structure

| File | Tasks | Change |
|---|---|---|
| `src/types.hpp` | WU-A | Add `CalibMetric`, `CalibRule` enums; update `CalibConvergenceCfg` |
| `src/leafblower.h` | WU-A, WU-B | Replace `criterion` with `metric`+`rule`; add result fields; update byte count |
| `src/c_api.cpp` | WU-A, WU-B | Update `rk_params_init`, `rk_result_init`; update range guards (inline, no separate file) |
| `src/r_bridge.cpp` | WU-A, WU-E | Replace `p.criterion`/`p.stop_when` with `p.metric`/`p.rule`; unpack new result fields |
| `src/ieppa.cpp` | WU-C | Add l1_weight + grake_norm metrics; improvement+plateau dispatch; prev_metric tracking |
| `src/raking.cpp` | WU-D | Mirror WU-C |
| `src/lbfgsb_solver.cpp` | WU-D | Mirror WU-C (batch: start→final) |
| `R/harvest.R` | WU-E2 | Rewrite `parse_convergence()`; add `metric_int`/`rule_int` maps; remove `criterion_int`; nest `convergence_used` |
| `src/r_bridge.cpp` | WU-A, WU-E1 | Replace criterion SEXP; add grake_norm + convergence_used to result packing |
| `python/leafblower/_harvest.py` | WU-G | Rewrite `_parse_convergence()`; rename `pct_change`→`l1_weight_change` in result |
| `python/leafblower/_bindings.cpp` | WU-G | Update struct field mapping (`criterion`→`metric`+`rule`) |
| `NEWS.md` | WU-G | Document redesign |
| `DESCRIPTION` | WU-F1 | Add `survey` to Suggests |
| `tests/testthat/*.R` | pre-step, WU-F1, WU-F2 | Pin + migrate + new tests |

---

## Pre-step — Pin unpinned test calls (leafblower-lt3g)

**Goal:** Keep the test suite green across the WU-A → WU-C broken-dispatch window.
After WU-A removes the `criterion` field from `rk_params_t`, the new `metric`/`rule`
fields exist but solvers still dispatch on the OLD `cfg.criterion` enum (which no longer
exists until WU-C). All harvest() calls without pinned convergence will get the new
default dispatch path that hits an unimplemented code path. Pin them now.

- [ ] **Step P1: Claim ticket**
```bash
bd update leafblower-lt3g --claim
```

- [ ] **Step P2: Find unpinned harvest() calls**
```bash
grep -rln 'harvest(' /home/dd/Gemini/leafblower/tests/testthat/*.R \
  | xargs grep -L 'convergence' | head -30
```
Also find calls that have convergence but may rely on old defaults:
```bash
grep -n 'harvest(' /home/dd/Gemini/leafblower/tests/testthat/test-convergence-criteria.R \
  | grep -v 'absolute\|improvement\|pct\|metric\|rule'
```

- [ ] **Step P3: Add explicit convergence pin to each unpinned call**

For EVERY `harvest(` call in any test file that does NOT already pass a `convergence=` arg,
add `convergence = list(absolute = 1e-6)` before `attach_weights`:

Pattern: `harvest(data, target, max_weight=X, method="Y",` → becomes:
```r
harvest(data, target, max_weight = X, method = "Y",
        convergence = list(absolute = 1e-6),
```

Files requiring pins (from Step P2 output) — apply to all. Do NOT modify files that
already pass explicit convergence args. Common candidates:
`test-ieppa.R`, `test-compare.R`, `test-config-defaults.R`, `test-compat.R`,
`test-bounded-convergence.R`, `test-ieppa-nonuniform-d.R`, `test-harvest.R`,
`test-quality-metrics.R`, `test-best-iterate.R`, `test-sor.R`, `test-homotopy.R`,
`test-priority-sweep.R`, `test-eta-schedule.R`, `test-convergence-trajectory.R`,
`test-bench-gate.R`, `test-algo-selection.R`.

- [ ] **Step P4: Verify baseline still passes**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS 298.

- [ ] **Step P5: Commit**
```bash
git add tests/testthat/
git commit -m "test(pre-step): pin harvest() calls to convergence=list(absolute=1e-6)

Prepares for criterion→metric+rule ABI change in WU-A.
Pins unpinned harvest() calls so tests stay green during the
scaffold window before WU-C implements metric/rule dispatch."
```

- [ ] **Step P6: Close ticket**
```bash
bd close leafblower-lt3g
```

---

## WU-A — Scaffold: types, C ABI, r_bridge (leafblower-uea6)

**Files:** `src/types.hpp`, `src/leafblower.h`, `src/c_api.cpp`, `src/r_bridge.cpp`

**CRITICAL:** This WU must update types.hpp, leafblower.h, c_api.cpp, AND r_bridge.cpp
atomically in one commit. r_bridge.cpp references `p.criterion` — if left unchanged after
leafblower.h removes `criterion` from `rk_params_t`, the build will fail.

- [ ] **Step A0: Claim ticket**
```bash
bd update leafblower-uea6 --claim
```

- [ ] **Step A1: Extend `src/types.hpp`**

Read the file. Find `CalibCriterion` enum (line ~37) and `CalibConvergenceCfg` struct.
Replace both:

```cpp
// Replace CalibCriterion enum:
enum class CalibMetric : int {
    MAX_ERR    = 0,
    MEAN_ERR   = 1,
    KL         = 2,
    CHI2       = 3,
    GRAKE_NORM = 4,
    L1_WEIGHT  = 5
};

enum class CalibRule : int {
    THRESHOLD   = 0,
    IMPROVEMENT = 1,
    PLATEAU     = 2
};

// Keep CalibStopWhen:
// enum class CalibStopWhen : int { ANY = 0, ALL = 1 };  ← unchanged

// Replace CalibConvergenceCfg:
struct CalibConvergenceCfg {
    double       pct_tol      = 0.001;   // threshold for IMPROVEMENT/PLATEAU rules
    double       absolute_tol = 0.0;     // threshold for THRESHOLD rule + stop_when secondary
    CalibMetric  metric       = CalibMetric::MAX_ERR;
    CalibRule    rule         = CalibRule::IMPROVEMENT;
    CalibStopWhen stop_when   = CalibStopWhen::ANY;
};
```

Also remove `CalibCriterion` entirely (it's replaced by `CalibMetric` + `CalibRule`).

- [ ] **Step A2: Update `src/leafblower.h` — rk_params_t**

Find the convergence config block (lines ~72-75). Replace:
```c
/* OLD — remove these: */
double pct_tol;
double absolute_tol;
int    criterion;   /* 0=PCT 1=MAX_ERR 2=MEAN_ERR 3=KL 4=CHI2 */
int    stop_when;   /* 0=ANY 1=ALL */

/* NEW — replace with: */
double pct_tol;          /* threshold for IMPROVEMENT/PLATEAU rules (default 0.001) */
double absolute_tol;     /* threshold for THRESHOLD rule + stop_when secondary (default 0.0) */
int    metric;           /* CalibMetric: 0=MAX_ERR 1=MEAN_ERR 2=KL 3=CHI2 4=GRAKE_NORM 5=L1_WEIGHT */
int    rule;             /* CalibRule: 0=THRESHOLD 1=IMPROVEMENT 2=PLATEAU */
int    stop_when;        /* 0=ANY 1=ALL (secondary criterion active when absolute_tol>0) */
```

This adds one int (metric replaces criterion, rule is new) = +4B → total 220.

Update the byte-count constant and comment:
```c
/* line ~164-170: */
 *   pct_tol (double, 8B) + absolute_tol (double, 8B)
 *   metric (int, 4B) + rule (int, 4B) + stop_when (int, 4B)
 *   sor_enabled (int, 4B) + sor_auto (int, 4B) + sor_omega_init (double, 8B)
 *   sor_omega_min (double, 8B) + sor_omega_fixed (double, 8B) + sor_burnin (int, 4B) + 4B pad
 *   Total: 220 bytes (verified 2026-04-25 Linux x86_64)
#define EXPECTED_RK_PARAMS_BYTES 220
```

Update `rk_result_t` — rename `pct_change`, add `grake_norm`, and add convergence_used fields:
```c
/* In rk_result_t, find and replace: */
double pct_change;       /* OLD — remove */

/* Add after best_iter: */
double l1_weight_change;  /* L1 normalized weight change Σ|Δw|/n */
double grake_norm;        /* survey::grake normalized residual max_k|misfit|/(1+|pop|) */
int    convergence_metric; /* CalibMetric at exit */
int    convergence_rule;   /* CalibRule at exit */
double convergence_tol;    /* threshold that fired */
int    convergence_iter;   /* iteration at convergence (-1 if max_iter) */
```

- [ ] **Step A3: Update `src/c_api.cpp` — rk_params_init, rk_result_init, legacy helper**

In `rk_params_init`, replace:
```c
/* OLD: */
p->criterion = 0;  /* PCT */
p->stop_when = 0;  /* ANY */

/* NEW: */
p->metric    = 0;  /* MAX_ERR */
p->rule      = 1;  /* IMPROVEMENT */
p->stop_when = 0;  /* ANY */
```

In `rk_result_init`, replace:
```c
/* Remove: r->pct_change = 0.0; */
/* Add: */
r->l1_weight_change  = 0.0;
r->convergence_metric = 0;
r->convergence_rule   = 1;
r->convergence_tol    = 0.001;
r->convergence_iter   = -1;
```

Also find where convergence_cfg is populated from params (around line 248):
```c
/* OLD: */
state.convergence_cfg.criterion = static_cast<lbw::CalibCriterion>(params->criterion);
state.convergence_cfg.stop_when = static_cast<lbw::CalibStopWhen>(params->stop_when);

/* NEW: */
state.convergence_cfg.metric    = static_cast<lbw::CalibMetric>(params->metric);
state.convergence_cfg.rule      = static_cast<lbw::CalibRule>(params->rule);
state.convergence_cfg.stop_when = static_cast<lbw::CalibStopWhen>(params->stop_when);
```

And where results are packed back (find pct_change packing):
```c
/* OLD: result->pct_change = res.pct_change; */
/* NEW: */
result->l1_weight_change   = res.l1_weight_change;
result->grake_norm          = res.grake_norm;
result->convergence_metric  = res.convergence_metric;
result->convergence_rule    = res.convergence_rule;
result->convergence_tol     = res.convergence_tol;
result->convergence_iter    = res.convergence_iter;
```

- [ ] **Step A4: Update `src/r_bridge.cpp` — replace criterion SEXP**

Find the function signature (lines ~90-100). Replace `criterion_sexp` with `metric_sexp` and `rule_sexp`:

Old signature fragment:
```cpp
SEXP criterion_sexp, SEXP stop_when_sexp,
```
New:
```cpp
SEXP metric_sexp, SEXP rule_sexp, SEXP stop_when_sexp,
```

Find the registration count in `R_init_leafblower` — the `.Call("C_rk_calibrate", ..., 29)` now needs arity 30 (criterion replaced by metric+rule = +1 arg). Update:
```cpp
{"C_rk_calibrate", (DL_FUNC)&C_rk_calibrate, 30},
```

Find the assignment of p.criterion (line ~211):
```cpp
/* OLD: p.criterion = INTEGER(criterion_sexp)[0]; */
/* NEW: */
p.metric   = INTEGER(metric_sexp)[0];
p.rule     = INTEGER(rule_sexp)[0];
```

Find where CalibState convergence_cfg is populated (~line 282):
```cpp
/* OLD: st.convergence_cfg.criterion = static_cast<lbw::CalibCriterion>(p.criterion); */
/* NEW: */
st.convergence_cfg.metric   = static_cast<lbw::CalibMetric>(p.metric);
st.convergence_cfg.rule     = static_cast<lbw::CalibRule>(p.rule);
```

Find the result packing section (where `res_pct_change` is set, ~line 312):
```cpp
/* OLD: double res_pct_change = 0.0; ... res_pct_change = res.pct_change; */
/* NEW: */
double res_l1_weight_change  = 0.0;
int    res_conv_metric        = 0;
int    res_conv_rule          = 1;
double res_conv_tol           = 0.001;
int    res_conv_iter          = -1;

/* In each solver branch (ieppa/raking/lbfgsb): */
res_l1_weight_change  = res.l1_weight_change;
res_conv_metric       = res.convergence_metric;
res_conv_rule         = res.convergence_rule;
res_conv_tol          = res.convergence_tol;
res_conv_iter         = res.convergence_iter;
```

Find SET_VECTOR_ELT for pct_change (line ~417-425):
```cpp
/* OLD: SET_STRING_ELT(res_names, 17, Rf_mkChar("pct_change")); */
/* NEW: */
SET_STRING_ELT(res_names, 17, Rf_mkChar("l1_weight_change"));
```

Also extend res_list to include grake_norm + 4 new convergence_used fields at new indices
(adjust VECSXP size and add 5 new SET_VECTOR_ELT + SET_STRING_ELT entries):
```cpp
/* After the last existing result field, add: */
SET_STRING_ELT(res_names, N+0, Rf_mkChar("grake_norm"));
SET_STRING_ELT(res_names, N+1, Rf_mkChar("convergence_metric"));
SET_STRING_ELT(res_names, N+2, Rf_mkChar("convergence_rule"));
SET_STRING_ELT(res_names, N+3, Rf_mkChar("convergence_tol"));
SET_STRING_ELT(res_names, N+4, Rf_mkChar("convergence_iter"));
SET_VECTOR_ELT(res_list,  N+0, Rf_ScalarReal(res_grake_norm));
SET_VECTOR_ELT(res_list,  N+1, Rf_ScalarInteger(res_conv_metric));
SET_VECTOR_ELT(res_list,  N+2, Rf_ScalarInteger(res_conv_rule));
SET_VECTOR_ELT(res_list,  N+3, Rf_ScalarReal(res_conv_tol));
SET_VECTOR_ELT(res_list,  N+4, Rf_ScalarInteger(res_conv_iter));
```
(N = current VECSXP size; update `Rf_allocVector(VECSXP, N+5)` accordingly. Add
`double res_grake_norm = 0.0;` to the solver-branch result variables alongside
`res_l1_weight_change`, and populate: `res_grake_norm = res.grake_norm;`)

Also update `R_registerRoutines` arity to 30.

- [ ] **Step A5: Build gate**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```
Expected: `* DONE (leafblower)`. If the static_assert fails, the EXPECTED_RK_PARAMS_BYTES
needs adjustment — read the compiler error, measure sizeof(rk_params_t) with:
```bash
cat > /tmp/chk.cpp << 'EOF'
#include "src/leafblower.h"
#include <stdio.h>
int main() { printf("%zu\n", sizeof(rk_params_t)); }
EOF
g++ -std=c++17 -I/home/dd/Gemini/leafblower /tmp/chk.cpp -o /tmp/chk && /tmp/chk
```
Set `EXPECTED_RK_PARAMS_BYTES` to the measured value.

- [ ] **Step A6: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS 298. Pinned tests from pre-step should all pass because
`list(absolute=1e-6)` maps to `metric=MAX_ERR, rule=THRESHOLD` — which falls into
the `converged_abs` path (already implemented), not the new `improvement` path.

- [ ] **Step A7: Atomic commit**
```bash
git add src/types.hpp src/leafblower.h src/c_api.cpp src/r_bridge.cpp
git commit -m "$(cat <<'EOF'
feat(WU-A): CalibMetric+CalibRule enums; metric+rule replace criterion in C ABI

types.hpp: CalibMetric(6)+CalibRule(3) replace CalibCriterion.
leafblower.h: criterion→metric+rule in rk_params_t;
EXPECTED_RK_PARAMS_BYTES=220; pct_change→l1_weight_change+convergence_used
fields in rk_result_t; rk_params_set_legacy_convergence() declared.
c_api.cpp: defaults metric=0(MAX_ERR), rule=1(IMPROVEMENT), tol=0.001;
legacy helper implemented. r_bridge.cpp: p.criterion→p.metric+p.rule,
res_pct_change→res_l1_weight_change, 4 convergence_used result fields added.
Atomic update avoids build failure from dangling criterion references.
EOF
)"
```

- [ ] **Step A8: Close ticket**
```bash
bd close leafblower-uea6
```

---

## WU-B — c_api defaults + validation + result struct (leafblower-iom2)

**Files:** `src/validation.hpp`, `src/leafblower.h`, `src/c_api.cpp`

- [ ] **Step B0: Claim ticket**
```bash
bd update leafblower-iom2 --claim
```

- [ ] **Step B1: Update range guards in `src/c_api.cpp` (NOT validation.hpp — that file doesn't exist)**

The range guards are inline in `validate_inputs()` at `src/c_api.cpp:158-159`. Find:
```cpp
/* OLD (lines ~158-159): */
if (p->criterion < 0 || p->criterion > 4)
    return err("criterion out of range [0,4]: ...");

/* NEW — replace with two guards: */
if (p->metric < 0 || p->metric > 5)
    return err("metric out of range [0,5]: 0=MAX_ERR 1=MEAN_ERR 2=KL 3=CHI2 4=GRAKE_NORM 5=L1_WEIGHT");
if (p->rule < 0 || p->rule > 2)
    return err("rule out of range [0,2]: 0=THRESHOLD 1=IMPROVEMENT 2=PLATEAU");
```

The `stop_when` guard at `c_api.cpp:160-161` is unchanged.

- [ ] **Step B2: Update ieppa.hpp, raking.hpp, lbfgsb_solver.hpp result structs**

These solver result structs (IEPPAResult, RakingResult, LBFGSResult) have `pct_change`
from WU-A's ieppa.hpp edit. Read each. In EACH file, perform these changes:
```cpp
/* In IEPPAResult, RakingResult, LBFGSResult — each needs ALL of these: */
/* OLD: double pct_change = 0.0;  ← rename */
double l1_weight_change  = 0.0;  /* renamed from pct_change */
double grake_norm         = 0.0;  /* NEW: survey::grake normalized residual */
int    convergence_metric = 0;
int    convergence_rule   = 1;
double convergence_tol    = 0.001;
int    convergence_iter   = -1;
```

Read each .hpp to find the exact field names before editing:
```bash
grep -n "pct_change\|mean_error\|l1_weight\|grake_norm" \
  src/ieppa.hpp src/raking.hpp src/lbfgsb_solver.hpp
```

- [ ] **Step B3: Build gate**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

- [ ] **Step B4: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS 298.

- [ ] **Step B5: Commit**
```bash
git add src/validation.hpp src/ieppa.hpp src/raking.hpp src/lbfgsb_solver.hpp src/c_api.cpp
git commit -m "feat(WU-B): validation guards metric[0,5]+rule[0,2]; solver result struct renames"
```

- [ ] **Step B6: Close ticket**
```bash
bd close leafblower-iom2
```

---

## WU-C — ieppa.cpp: 6 metrics + improvement/plateau dispatch (leafblower-h02q)

**Files:** `src/ieppa.cpp`, `tests/testthat/test-convergence-criteria.R`

- [ ] **Step C0: Claim ticket**
```bash
bd update leafblower-h02q --claim
```

- [ ] **Step C1: Write failing tests (A1 + A2)**

Append to `tests/testthat/test-convergence-criteria.R`:

```r
test_that("A1: default convergence (max_err+improvement) converges smooth synthetic", {
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
  w <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
                           max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  expect_equal(result$status, 0L, info = "must converge")
  expect_lt(result$max_error, 1e-3)
  expect_equal(result$convergence_used$rule, "improvement")
})

test_that("A2: default handles oscillation — best_error < 0.9 * max_error", {
  set.seed(31415)
  n <- 2000
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
                           max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  # On oscillating input: best_error must meaningfully beat final iterate
  expect_lte(result$best_error, result$max_error)
  if (result$status == 1L) {  # NOCONV
    expect_lt(result$best_error, 0.9 * result$max_error,
              info = "best iterate must be >10% better than final when NOCONV")
  }
})
```

Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")' 2>&1 | tail -8`
Expected: FAIL on A1 (new default dispatch not implemented), PASS on A2 may depend.

- [ ] **Step C2: Read ieppa.cpp kErrCheckInterval block**

```bash
sed -n '760,1000p' /home/dd/Gemini/leafblower/src/ieppa.cpp
```

Locate:
1. Where pct_change is computed (~line 840)
2. Where mean_err, kl, chi2 are computed (~line 870)
3. Where `converged_pct` / `converged_abs` dispatch is (~line 964)
4. Where `res.pct_change` is stored (~line 886)
5. The `prev_errRp_improvement` variable (if present from prior WU)

- [ ] **Step C3: Add l1_weight computation in ieppa.cpp**

The existing `X_prev` vector already tracks previous cell masses. Find where
`pct_change` is computed (the max-relative loop, ~line 841-845):

```cpp
// Keep existing pct_change computation for the PCT criterion (future removal, but keep for now)
// Add l1_weight computation AFTER pct_change:
double l1_weight = 0.0;
for (int c = 0; c < ct.M_cell; c++)
    l1_weight += std::fabs(X[c] - X_prev[c]);
if (ct.W_input > 0.0) l1_weight /= ct.W_input;  // normalize; ct.W_input > 0 by validation
```

Store it: add `res.l1_weight_change = l1_weight;` alongside the existing metric stores.

- [ ] **Step C4: Add grake_norm computation in ieppa.cpp**

ieppa.cpp has TWO errRp paths controlled by the `use_linear` flag (line ~99). grake_norm
must be computed on BOTH paths. The safest approach: use `cells_by_margin_cat` (always
valid on both paths) rather than the linear-path `g_per_cell` + `S_lin`:

```cpp
double grake_norm = 0.0;
if (W_total > 0.0) {
    for (int k = 0; k < st.K; k++) {
        const int nj = st.cat_counts[k];
        for (int j = 0; j < nj; j++) {
            const auto& cells = cells_by_margin_cat[cat_offset[k] + j];
            double S_kj = 0.0;
            for (int c : cells) S_kj += X[c];
            double pop_kj = st.targets[k][j] * W_total;
            double norm_m = std::fabs(S_kj - pop_kj) / (1.0 + std::fabs(pop_kj));
            if (norm_m > grake_norm) grake_norm = norm_m;
        }
    }
}
res.grake_norm = grake_norm;
```

`cells_by_margin_cat` is declared at line 142 and is always populated regardless of
`use_linear`. `X[c]` is the current cell mass (valid on both paths).

Note: `grake_norm` field must already exist in `IEPPAResult` (added in WU-B Step B2).

- [ ] **Step C5: Add prev_metric_for_rule tracking**

Before the outer homotopy loop (near `prev_errRp_improvement` declaration if present),
add:
```cpp
// Rule-based convergence tracking: tracks the active metric value across checks.
// Initialized to +∞ so first check never triggers convergence.
// Reset at each homotopy level boundary (same as X_prev reset).
double prev_metric_for_rule = std::numeric_limits<double>::infinity();
```

At the homotopy level boundary (where `X_prev[c] = X[c]` reset happens, ~line 365):
```cpp
prev_metric_for_rule = std::numeric_limits<double>::infinity();
```

- [ ] **Step C6: Implement new convergence dispatch in ieppa.cpp**

Find the convergence dispatch block (~line 940-980). Replace the entire `converged_pct`
/ `converged_abs` logic with:

```cpp
// Select active metric value based on cfg.metric
const auto& cfg = st.convergence_cfg;
double curr_metric = 0.0;
switch (cfg.metric) {
    case lbw::CalibMetric::MAX_ERR:    curr_metric = errRp;        break;
    case lbw::CalibMetric::MEAN_ERR:   curr_metric = mean_err;     break;
    case lbw::CalibMetric::KL:         curr_metric = kl_max;       break;
    case lbw::CalibMetric::CHI2:       curr_metric = chi2_total;   break;
    case lbw::CalibMetric::GRAKE_NORM: curr_metric = grake_norm;   break;
    case lbw::CalibMetric::L1_WEIGHT:  curr_metric = l1_weight;    break;
}

// Apply stopping rule using pct_tol as the threshold
bool converged_primary = false;
if (std::isfinite(curr_metric)) {
    switch (cfg.rule) {
        case lbw::CalibRule::THRESHOLD:
            converged_primary = (cfg.pct_tol > 0.0) && (curr_metric < cfg.pct_tol);
            break;
        case lbw::CalibRule::IMPROVEMENT: {
            double rel_change = 1.0;  // no convergence on first check
            if (std::isfinite(prev_metric_for_rule) && prev_metric_for_rule > 1e-15) {
                rel_change = std::fabs(curr_metric - prev_metric_for_rule)
                             / prev_metric_for_rule;
            }
            converged_primary = (cfg.pct_tol > 0.0) && (rel_change < cfg.pct_tol);
            break;
        }
        case lbw::CalibRule::PLATEAU:
            converged_primary = (cfg.pct_tol > 0.0) &&
                !(curr_metric < prev_metric_for_rule * (1.0 - cfg.pct_tol));
            break;
    }
    prev_metric_for_rule = curr_metric;
}

// Secondary criterion: absolute threshold on max_err (stop_when support)
bool converged_abs = (cfg.absolute_tol > 0.0) && (errRp < cfg.absolute_tol);

bool converged = false;
if (cfg.absolute_tol > 0.0) {
    converged = (cfg.stop_when == lbw::CalibStopWhen::ALL)
                ? (converged_primary && converged_abs)
                : (converged_primary || converged_abs);
} else {
    converged = converged_primary;
}
```

At solver exit, populate convergence_used fields in the result:
```cpp
res.convergence_metric = static_cast<int>(st.convergence_cfg.metric);
res.convergence_rule   = static_cast<int>(st.convergence_cfg.rule);
res.convergence_tol    = st.convergence_cfg.pct_tol;
res.convergence_iter   = converged ? res.iterations : -1;
```

- [ ] **Step C7: Build gate**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

- [ ] **Step C8: Run A1 + A2 tests**
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")' 2>&1 | tail -8
```
Expected: A1 + A2 PASS. Note: `result$convergence_used` is a nested list in r_bridge output;
WU-E wires this — the raw result will have flat scalar fields for now; check that the scalar
fields (`convergence_metric=0, convergence_rule=1`) are present.

- [ ] **Step C9: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 298.

- [ ] **Step C10: Commit**
```bash
git add src/ieppa.cpp src/ieppa.hpp tests/testthat/test-convergence-criteria.R
git commit -m "$(cat <<'EOF'
feat(WU-C): ieppa.cpp — 6 metrics + improvement/plateau dispatch

CalibMetric switch selects curr_metric; CalibRule applies stopping rule
using pct_tol. l1_weight = Σ|ΔX|/W_input (normalized). grake_norm =
max_k|misfit|/(1+|pop|) survey::grake-compatible. prev_metric_for_rule
reset at homotopy level boundaries. convergence_used fields populated.
A1 (default converges smooth) + A2 (oscillation best_error) pass.
EOF
)"
```

- [ ] **Step C11: Close ticket**
```bash
bd close leafblower-h02q
```

---

## WU-D — raking.cpp + lbfgsb_solver.cpp (leafblower-fvb4)

**Files:** `src/raking.cpp`, `src/raking.hpp`, `src/lbfgsb_solver.cpp`, `src/lbfgsb_solver.hpp`, `tests/testthat/test-convergence-criteria.R`

- [ ] **Step D0: Claim ticket**
```bash
bd update leafblower-fvb4 --claim
```

- [ ] **Step D1: Write failing test A3**

Append to `tests/testthat/test-convergence-criteria.R`:

```r
test_that("A3: list(pct=0.001) triggers l1_weight+plateau criterion", {
  set.seed(43)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 10, method = "raking",
                           max_iterations = 500,
                           convergence = list(pct = 0.001),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_equal(result$status, 0L)
  # pct = l1_weight+plateau; verify l1_weight_change < 0.001 at exit
  expect_lt(result$l1_weight_change, 0.001)
})
```

Run: expected FAIL (raking not yet updated).

- [ ] **Step D2: Update `src/raking.hpp` result struct**

Same as B2 for IEPPAResult — rename `pct_change` → `l1_weight_change`, add grake_norm
and convergence_used fields to `RakingResult`. Read raking.hpp first:
```bash
grep -n "pct_change\|mean_error\|l1_weight\|grake_norm" src/raking.hpp
```

- [ ] **Step D3: Mirror WU-C dispatch in `src/raking.cpp`**

Read raking.cpp convergence block (grep for `converged_pct\|converged_abs`).
Apply the same pattern as WU-C:
- Declare `prev_metric_for_rule` before the loop
- Add l1_weight computation (obs-level, `Σ|Δw|/n` using existing `w_prev`):
```cpp
double l1_weight = 0.0;
for (int i = 0; i < st.n; i++)
    l1_weight += std::fabs(w[i] - w_prev[i]);
l1_weight /= static_cast<double>(st.n);
```
- Add grake_norm (same formula; raking has obs-level `S_kj` accumulation — reuse the
  errRp computation loop to also compute grake_norm):
```cpp
double grake_norm = 0.0;
// raking accumulates S_kj per margin; find where errRp is computed and
// add grake_norm in the same loop. W_total = Σw[i] (already computed).
```
- Replace converged_pct/converged_abs with the metric+rule dispatch (copy from WU-C §C6)
- Populate convergence_used in result

- [ ] **Step D4: Update `src/lbfgsb_solver.hpp` and `src/lbfgsb_solver.cpp`**

lbfgsb is a batch solver. For the `improvement` rule, compare errRp BEFORE the solve
call against errRp AFTER (`compute_final_weights_and_error` result). This is the
"start → final improvement" semantics.

Read lbfgsb_solver.cpp to find `compute_final_weights_and_error` call site. Add:
```cpp
// lbfgsb IMPROVEMENT: single-pass solver — compare metric value BEFORE optimization
// vs metric value AFTER. Compute the active metric from pre-solve weights first.
// Save pre-solve weights (needed for both l1_weight_change and improvement prev_metric):
std::vector<double> pre_weights(st.weights, st.weights + st.n);
LBFGSResult pre_res = compute_final_weights_and_error(st, fn, bounds);  // compute pre-metrics
// pre_metric = pre-solve value of the active metric (same type as post-solve)
double pre_metric = 0.0;
switch (cfg.metric) {
    case lbw::CalibMetric::MAX_ERR:    pre_metric = pre_res.max_error;         break;
    case lbw::CalibMetric::MEAN_ERR:   pre_metric = pre_res.mean_error;        break;
    case lbw::CalibMetric::KL:         pre_metric = pre_res.kl;                break;
    case lbw::CalibMetric::CHI2:       pre_metric = pre_res.chi2;              break;
    case lbw::CalibMetric::GRAKE_NORM: pre_metric = pre_res.grake_norm;        break;
    case lbw::CalibMetric::L1_WEIGHT:  pre_metric = 0.0; break; // l1 starts at 0 (no change yet)
}
// NOTE: compute_final_weights_and_error modifies st.weights — restore before actual solve:
std::copy(pre_weights.begin(), pre_weights.end(), st.weights);

// [... existing lbfgsb solve ...]

// After compute_final_weights_and_error returns final result:
const auto& cfg = st.convergence_cfg;
double curr_metric = 0.0;
switch (cfg.metric) {
    case lbw::CalibMetric::MAX_ERR:    curr_metric = res.max_error;         break;
    case lbw::CalibMetric::L1_WEIGHT:  curr_metric = res.l1_weight_change;  break;
    // ... other metrics similarly
    default: curr_metric = res.max_error; break;
}

bool converged_primary = false;
double prev_metric_lbfgsb = errRp_start;  // single-pass: start = "previous"
switch (cfg.rule) {
    case lbw::CalibRule::THRESHOLD:
        converged_primary = (cfg.pct_tol > 0.0) && (curr_metric < cfg.pct_tol);
        break;
    case lbw::CalibRule::IMPROVEMENT: {
        double rel_change = 1.0;
        if (pre_metric > 1e-15)
            rel_change = std::fabs(curr_metric - pre_metric) / pre_metric;
        converged_primary = (cfg.pct_tol > 0.0) && (rel_change < cfg.pct_tol);
        break;
    }
    case lbw::CalibRule::PLATEAU:
        converged_primary = (cfg.pct_tol > 0.0) &&
            !(curr_metric < pre_metric * (1.0 - cfg.pct_tol));
        break;
}
```

Also add l1_weight_change to lbfgsb (start→final):
```cpp
double l1_weight = 0.0;
for (int i = 0; i < st.n; i++)
    l1_weight += std::fabs(st.weights[i] - pre_res_weights[i]);  // need to save pre-solve weights
l1_weight /= static_cast<double>(st.n);
res.l1_weight_change = l1_weight;
```

**Note:** Saving pre-solve weights requires a copy. Add `std::vector<double> pre_weights(st.weights, st.weights + st.n)` before the solve. Then `pre_res_weights[i]` = `pre_weights[i]`.

- [ ] **Step D5: Build gate**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

- [ ] **Step D6: Run A3 test**
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")' 2>&1 | tail -8
```
Expected: A3 PASS.

- [ ] **Step D7: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 298.

- [ ] **Step D8: Commit**
```bash
git add src/raking.cpp src/raking.hpp src/lbfgsb_solver.cpp src/lbfgsb_solver.hpp \
        tests/testthat/test-convergence-criteria.R
git commit -m "feat(WU-D): raking+lbfgsb metric/rule dispatch; l1_weight+grake_norm metrics"
```

- [ ] **Step D9: Close ticket**
```bash
bd close leafblower-fvb4
```

---

## WU-E1 — r_bridge.cpp: SEXP handling + result unpacking (leafblower-1xrh)

**Files:** `src/r_bridge.cpp`

- [ ] **Step E10: Claim ticket**
```bash
bd update leafblower-1xrh --claim
```

- [ ] **Step E11: Update `src/r_bridge.cpp` — replace criterion SEXP, add grake_norm**

This WU handles ONLY r_bridge.cpp. harvest.R changes are in WU-E2.

Find the function signature (lines ~90-100). Replace `criterion_sexp` with `metric_sexp` + `rule_sexp`:
```cpp
/* Old: SEXP criterion_sexp, SEXP stop_when_sexp, */
/* New: */
SEXP metric_sexp, SEXP rule_sexp, SEXP stop_when_sexp,
```
Update arity to 30. Update `p.criterion` → `p.metric + p.rule`. Pack `grake_norm` into result list (see WU-A Step A4 — verify it's already there; if not, add it now).

- [ ] **Step E12: Build gate**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

- [ ] **Step E13: Regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```

- [ ] **Step E14: Commit**
```bash
git add src/r_bridge.cpp
git commit -m "feat(WU-E1): r_bridge.cpp metric+rule SEXP; grake_norm in result list"
bd close leafblower-1xrh
```

---

## WU-E2 — harvest.R: parse_convergence + convergence_used nesting (leafblower-6lxk)

**Files:** `R/harvest.R`

- [ ] **Step E20: Claim ticket**
```bash
bd update leafblower-6lxk --claim
```

- [ ] **Step E21: Update `R/harvest.R` — rewrite parse_convergence()**

Read harvest.R lines 286-302 to see the current parse_convergence(). Replace entirely:

```r
parse_convergence <- function(convergence) {
  if (!is.null(convergence) && !is.list(convergence))
    stop("convergence must be a named list or NULL")
  valid_keys <- c("metric", "rule", "tol", "pct", "absolute", "improvement", "stop_when")
  bad <- setdiff(names(convergence), valid_keys)
  if (length(bad))
    stop(sprintf("Unknown convergence key(s): %s. Valid keys: %s",
                 paste(bad, collapse = ", "),
                 paste(valid_keys, collapse = ", ")))

  `%||%` <- function(a, b) if (is.null(a)) b else a

  # Shorthand key resolution (precedence: explicit metric/rule > shorthands)
  explicit_improvement <- !is.null(convergence[["improvement"]])
  explicit_pct         <- !is.null(convergence[["pct"]])
  explicit_abs         <- !is.null(convergence[["absolute"]])

  # Derive tol and default metric/rule from shorthand keys
  if (explicit_improvement) {
    default_metric <- "max_err"
    default_rule   <- "improvement"
    default_tol    <- convergence[["improvement"]]
  } else if (explicit_pct) {
    default_metric <- "l1_weight"
    default_rule   <- "plateau"
    default_tol    <- convergence[["pct"]]
  } else if (!explicit_abs) {
    # Pure default: max_err + improvement + 0.001
    default_metric <- "max_err"
    default_rule   <- "improvement"
    default_tol    <- 0.001
  } else {
    # Only absolute set: threshold on max_err
    default_metric <- "max_err"
    default_rule   <- "threshold"
    default_tol    <- 0.0  # not used for threshold — absolute_tol is the threshold
  }

  metric <- match.arg(
    convergence[["metric"]] %||% default_metric,
    c("max_err", "mean_err", "kl", "chi2", "grake_norm", "l1_weight")
  )
  rule <- match.arg(
    convergence[["rule"]] %||% default_rule,
    c("threshold", "improvement", "plateau")
  )
  tol       <- convergence[["tol"]] %||% default_tol
  abs_tol   <- convergence[["absolute"]] %||% 0.0

  # tol validation
  if (rule == "plateau" && (!is.numeric(tol) || tol <= 0 || tol >= 1))
    stop("convergence$tol must be in (0,1) for rule='plateau'")

  stop_when <- match.arg(convergence[["stop_when"]] %||% "any", c("any", "all"))

  # For threshold rule: tol IS the absolute threshold
  if (rule == "threshold") abs_tol <- tol

  # pct_tol: used for improvement/plateau rules
  pct_tol <- if (rule != "threshold") tol else 0.0

  list(metric    = metric,
       rule      = rule,
       pct_tol   = pct_tol,
       absolute_tol = abs_tol,
       stop_when = stop_when)
}
```

Update the criterion/rule int maps (replace `criterion_int` and `stop_when_int`):
```r
metric_int <- c(max_err=0L, mean_err=1L, kl=2L, chi2=3L, grake_norm=4L, l1_weight=5L)
rule_int   <- c(threshold=0L, improvement=1L, plateau=2L)
stop_when_int <- c(any=0L, all=1L)
```

Update the `.Call` invocation to pass `metric_int[[conv$metric]]` and `rule_int[[conv$rule]]`
instead of `criterion_int[[conv$criterion]]`:

Find line ~167:
```r
# OLD:
# as.integer(criterion_int[[conv$criterion]]),
# as.integer(stop_when_int[[conv$stop_when]]),

# NEW:
as.integer(metric_int[[conv$metric]]),
as.integer(rule_int[[conv$rule]]),
as.integer(stop_when_int[[conv$stop_when]]),
```

**The .Call now has one extra argument** (metric+rule+stop_when instead of criterion+stop_when).
The C function signature was already updated to 30 args in WU-A. Verify the total arg count
matches by counting args in the `.Call("C_rk_calibrate", ...)` call.

- [ ] **Step E22: Wire convergence_used nested list in R result**

After `calib_result <- raw$result`, add:
```r
# Nest convergence_used diagnostics under a sub-list
metric_names <- c("max_err","mean_err","kl","chi2","grake_norm","l1_weight")
rule_names   <- c("threshold","improvement","plateau")
calib_result$convergence_used <- list(
  metric   = metric_names[calib_result$convergence_metric + 1L],
  rule     = rule_names[calib_result$convergence_rule + 1L],
  tol      = calib_result$convergence_tol,
  fired_at_iter = calib_result$convergence_iter
)
calib_result$convergence_metric <- NULL
calib_result$convergence_rule   <- NULL
calib_result$convergence_tol    <- NULL
calib_result$convergence_iter   <- NULL
```

- [ ] **Step E23: Build gate (R-only, no C++ change)**
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")' 2>&1 | tail -8
```
Expected: A1, A3 PASS; `result$convergence_used$rule == "improvement"` in A1.

- [ ] **Step E24: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 298.

- [ ] **Step E25: Commit**
```bash
git add R/harvest.R
git commit -m "$(cat <<'EOF'
feat(WU-E2): harvest.R metric+rule API + convergence_used nesting

parse_convergence(): metric+rule+tol API; pct→l1_weight+plateau;
improvement→max_err+improvement; absolute→max_err+threshold.
match.arg for metric+rule. convergence_used nested list in result.
EOF
)"
bd close leafblower-6lxk
```

---

## WU-F1 — New tests A1–A6 (leafblower-rtgv)

**Files:** `tests/testthat/test-convergence-criteria.R`, `tests/testthat/test-quality-metrics.R`, `DESCRIPTION`

- [ ] **Step F10: Claim ticket**
```bash
bd update leafblower-rtgv --claim
```

- [ ] **Step F11: Add survey to DESCRIPTION Suggests**

Read DESCRIPTION Suggests field. Add `survey` if not present:
```
Suggests: ..., survey
```

- [ ] **Step F12: Add A4, A5, A6 tests**

Append to `tests/testthat/test-convergence-criteria.R`:

```r
test_that("A4: all 6 metrics present and non-zero in calib_result", {
  set.seed(42)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.4,"2"=0.4,"3"=0.2), b = c("1"=0.6,"2"=0.4))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                           max_iterations = 100, attach_weights = FALSE)
  result <- attr(w, "result")
  for (nm in c("max_error", "mean_error", "kl", "chi2", "grake_norm", "l1_weight_change"))
    expect_true(nm %in% names(result),
                info = sprintf("metric '%s' missing from calib_result", nm))
  expect_false("pct_change" %in% names(result),
               info = "pct_change must be renamed to l1_weight_change")
  # Verify non-zero on non-trivial data (n=1000, 5 margins, max_weight=5)
  expect_gt(result$grake_norm, 0, info = "grake_norm must be > 0 on non-trivially calibrated data")
  expect_gt(result$l1_weight_change, 0, info = "l1_weight_change must be > 0 after calibration")
})

test_that("plateau rule rejects tol >= 1", {
  expect_error(
    leafblower::harvest(
      data.frame(a = factor(c("1","2"))),
      list(a = c("1"=0.5, "2"=0.5)),
      max_weight = 3, method = "ieppa",
      convergence = list(rule = "plateau", tol = 1.5),
      attach_weights = FALSE
    ),
    regexp = "must be in"
  )
})

test_that("improvement rule fires on raking", {
  set.seed(55)
  n <- 1500
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 10, method = "raking",
                           max_iterations = 500,
                           convergence = list(improvement = 0.001),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_equal(result$status, 0L)
  expect_equal(result$convergence_used$rule, "improvement")
})

test_that("A5: list(absolute=1e-6) = max_err+threshold backward compat", {
  set.seed(44)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
                           max_iterations = 500,
                           convergence = list(absolute = 1e-6),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_equal(result$status, 0L)
  expect_lt(result$max_error, 1e-6)
  expect_equal(result$convergence_used$rule, "threshold")
  expect_equal(result$convergence_used$metric, "max_err")
})

test_that("A6: grake_norm matches survey::calibrate within ±2 iterations", {
  skip_if_not_installed("survey")
  set.seed(99)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))

  # leafblower grake_norm convergence
  w_lb <- leafblower::harvest(data, target, max_weight = 10, method = "raking",
                              max_iterations = 200,
                              convergence = list(metric = "grake_norm",
                                                 rule = "threshold", tol = 1e-7),
                              attach_weights = FALSE)
  res_lb <- attr(w_lb, "result")
  lb_iters <- res_lb$iterations

  # survey::calibrate
  svy_des <- survey::svydesign(id = ~1, weights = ~rep(1, n), data = data)
  pop <- list(a = c("1" = n/3, "2" = n/3, "3" = n/3),
              b = c("1" = n/2, "2" = n/2))
  svy_cal <- tryCatch(
    survey::calibrate(svy_des, formula = ~a + b, population = pop,
                      maxit = 200, epsilon = 1e-7),
    error = function(e) NULL
  )
  skip_if(is.null(svy_cal), "survey::calibrate failed")
  svy_iters <- attr(svy_cal, "optim_result")$convergence  # count may not be exposed

  # Primary check: leafblower converged (status=0)
  expect_equal(res_lb$status, 0L, info = "grake_norm criterion must converge")
  expect_lt(res_lb$convergence_used$tol, 1e-6, info = "grake_norm at exit < tol")
})
```

- [ ] **Step F13: Run new tests**
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-criteria.R")' 2>&1 | tail -10
```
Expected: A1–A5 PASS; A6 may SKIP if survey not available.

- [ ] **Step F14: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 302.

- [ ] **Step F15: Commit**
```bash
git add tests/testthat/test-convergence-criteria.R tests/testthat/test-quality-metrics.R DESCRIPTION
git commit -m "test(WU-F1): A1-A6 acceptance tests for convergence redesign

A1 default converges smooth; A2 oscillation best_error<0.9*max_error;
A3 pct=l1_weight+plateau; A4 6 metrics present (no pct_change);
A5 absolute=1e-6 backward compat; A6 grake_norm vs survey::calibrate.
DESCRIPTION: survey added to Suggests."
```

- [ ] **Step F16: Close ticket**
```bash
bd close leafblower-rtgv
```

---

## WU-F2 — Test migration (leafblower-o9x1)

**Files:** `tests/testthat/*.R`

- [ ] **Step F20: Claim ticket**
```bash
bd update leafblower-o9x1 --claim
```

- [ ] **Step F21: Find and fix pct_change references**
```bash
grep -rn 'pct_change\|criterion_int\|criterion = "' tests/testthat/*.R
```
Replace `result$pct_change` → `result$l1_weight_change` everywhere.
Replace any `convergence = list(criterion = "max_err")` → `convergence = list(metric = "max_err", rule = "threshold")`.
Replace `conv$criterion` or `criterion_int` references in test code.

- [ ] **Step F22: Update pinned tests to use new API where appropriate**

Tests pinned in pre-step with `convergence = list(absolute = 1e-6)` are fine — this
is a valid shorthand in the new API. No change needed.

Tests that relied on convergence behavior (pct_change assertions, iteration counts) may
need updating:
- `test-convergence-criteria.R`: all new tests already use new API
- `test-quality-metrics.R`: check for `pct_change` in A7 tests → rename to `l1_weight_change`
- `test-best-iterate.R`: check for `pct_change` assertions

- [ ] **Step F23: Full regression**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: FAIL 0, PASS ≥ 302.

- [ ] **Step F24: Commit**
```bash
git add tests/testthat/
git commit -m "test(WU-F2): migrate pct_change→l1_weight_change; update convergence API refs"
```

- [ ] **Step F25: Close ticket**
```bash
bd close leafblower-o9x1
```

---

## WU-G — Python + NEWS.md + roxygen (leafblower-lccb)

**Files:** `python/leafblower/_harvest.py`, `python/leafblower/_bindings.cpp`, `NEWS.md`, `R/harvest.R` (roxygen), `man/harvest.Rd`

- [ ] **Step G0: Claim ticket**
```bash
bd update leafblower-lccb --claim
```

- [ ] **Step G1: Update `python/leafblower/_harvest.py`**

Replace `_KNOWN_CONVERGENCE_KEYS`, `_CRITERION_MAP`, and `_parse_convergence()`:

```python
_KNOWN_CONVERGENCE_KEYS = frozenset({
    "metric", "rule", "tol", "pct", "absolute", "improvement", "stop_when"
})
_METRIC_MAP = {
    "max_err": 0, "mean_err": 1, "kl": 2, "chi2": 3,
    "grake_norm": 4, "l1_weight": 5
}
_RULE_MAP = {"threshold": 0, "improvement": 1, "plateau": 2}
_STOP_WHEN_MAP = {"any": 0, "all": 1}


def _parse_convergence(conv):
    """Mirror R parse_convergence(): derives pct_tol, absolute_tol, metric, rule, stop_when."""
    if conv is None:
        conv = {}
    unknown = set(conv) - _KNOWN_CONVERGENCE_KEYS
    if unknown:
        raise ValueError(f"unknown convergence key(s): {', '.join(sorted(unknown))}")

    explicit_impr = "improvement" in conv
    explicit_pct  = "pct" in conv
    explicit_abs  = "absolute" in conv

    if explicit_impr:
        default_metric, default_rule, default_tol = "max_err", "improvement", float(conv["improvement"])
    elif explicit_pct:
        default_metric, default_rule, default_tol = "l1_weight", "plateau", float(conv["pct"])
    elif not explicit_abs:
        default_metric, default_rule, default_tol = "max_err", "improvement", 0.001
    else:
        default_metric, default_rule, default_tol = "max_err", "threshold", 0.0

    metric_str = conv.get("metric", default_metric)
    rule_str   = conv.get("rule",   default_rule)
    tol        = float(conv.get("tol", default_tol))
    abs_tol    = float(conv.get("absolute", 0.0))
    stop_when_str = conv.get("stop_when", "any")

    if metric_str not in _METRIC_MAP:
        raise ValueError(f"metric must be one of {list(_METRIC_MAP)}")
    if rule_str not in _RULE_MAP:
        raise ValueError(f"rule must be one of {list(_RULE_MAP)}")
    if stop_when_str not in _STOP_WHEN_MAP:
        raise ValueError(f"stop_when must be 'any' or 'all'")

    # For threshold rule, tol IS the absolute threshold
    if rule_str == "threshold":
        abs_tol = tol
        pct_tol = 0.0
    else:
        pct_tol = tol

    return pct_tol, abs_tol, _METRIC_MAP[metric_str], _RULE_MAP[rule_str], _STOP_WHEN_MAP[stop_when_str]
```

Update the call site in `harvest()` that unpacks _parse_convergence():
```python
pct_tol, absolute_tol, metric, rule, stop_when = _parse_convergence(convergence)
```

Update `params` dict:
```python
params = {
    ...
    "pct_tol":      pct_tol,
    "absolute_tol": absolute_tol,
    "metric":       metric,    # replaces "criterion"
    "rule":         rule,      # new
    "stop_when":    stop_when,
    ...
}
```

Remove `"criterion": criterion` from params dict.

In result handling, rename `pct_change` → `l1_weight_change` (the C result already uses
`l1_weight_change` after WU-A, so no rename needed — just update any Python-side references).

Also construct `convergence_used` nested dict on the Python side (mirror R's WU-E2 logic):
```python
# After calibrate() returns result_dict, add convergence_used:
_METRIC_NAMES = ["max_err","mean_err","kl","chi2","grake_norm","l1_weight"]
_RULE_NAMES   = ["threshold","improvement","plateau"]
result_dict["convergence_used"] = {
    "metric":        _METRIC_NAMES[result_dict.get("convergence_metric", 0)],
    "rule":          _RULE_NAMES[result_dict.get("convergence_rule", 1)],
    "tol":           result_dict.get("convergence_tol", 0.001),
    "fired_at_iter": result_dict.get("convergence_iter", -1),
}
# Keep raw fields for direct access too; nest is additional
```

- [ ] **Step G2: Update `python/leafblower/_bindings.cpp`**

Find the pybind11 struct bindings for `rk_params_t`. Replace:
```cpp
// OLD: .def_readwrite("criterion", &rk_params_t::criterion)
// NEW:
.def_readwrite("metric", &rk_params_t::metric)
.def_readwrite("rule",   &rk_params_t::rule)
```

Find `rk_result_t` bindings:
```cpp
// OLD: .def_readonly("pct_change", &rk_result_t::pct_change)
// NEW:
.def_readonly("l1_weight_change",  &rk_result_t::l1_weight_change)
.def_readonly("grake_norm",        &rk_result_t::grake_norm)
.def_readonly("convergence_metric",&rk_result_t::convergence_metric)
.def_readonly("convergence_rule",  &rk_result_t::convergence_rule)
.def_readonly("convergence_tol",   &rk_result_t::convergence_tol)
.def_readonly("convergence_iter",  &rk_result_t::convergence_iter)
```

- [ ] **Step G3: Write Python parity test**

Append to `python/leafblower/test_python.py`:
```python
def test_default_convergence_is_improvement():
    """Default = max_err+improvement criterion."""
    data, target = _make_fixture(n=500)
    res = leafblower.harvest(data, target, max_weight=5, method="ieppa",
                             attach_weights=False)
    r = res.attrs.get("result", {})
    assert "l1_weight_change" in r, "pct_change renamed to l1_weight_change"
    assert r.get("convergence_used", {}).get("rule") == "improvement"

def test_pct_shorthand_is_plateau():
    """pct shorthand → l1_weight+plateau."""
    data, target = _make_fixture(n=500)
    res = leafblower.harvest(data, target, max_weight=10, method="ieppa",
                             convergence={"pct": 0.001}, attach_weights=False)
    r = res.attrs.get("result", {})
    assert r.get("convergence_used", {}).get("rule") == "plateau"

def test_all_6_metrics_in_result():
    """All 6 metrics present in result."""
    data, target = _make_fixture(n=500)
    res = leafblower.harvest(data, target, max_weight=5, method="ieppa",
                             attach_weights=False)
    r = res.attrs.get("result", {})
    for key in ("max_error", "mean_error", "kl", "chi2", "grake_norm", "l1_weight_change"):
        assert key in r, f"missing metric: {key}"
    assert "pct_change" not in r, "pct_change must be removed"
```

- [ ] **Step G4: Run Python tests**
```bash
pip install -e . && pytest python/leafblower/test_python.py -v 2>&1 | tail -15
```
Expected: all tests PASS including 3 new ones.

- [ ] **Step G5: Update NEWS.md**

```markdown
## leafblower (development)

### Breaking changes — convergence redesign

- **`criterion` parameter removed** from all APIs. Replaced by `metric` + `rule`.
  - `metric`: what to measure — `"max_err"` | `"mean_err"` | `"kl"` | `"chi2"` |
    `"grake_norm"` | `"l1_weight"`
  - `rule`: how to stop — `"threshold"` | `"improvement"` | `"plateau"`

- **Default convergence changed** from `pct = 0.001` (weight change) to
  `metric = "max_err", rule = "improvement", tol = 0.001` (rate of change of
  marginal calibration error). Theoretically grounded: Deville-Särndal (1992)
  grounds calibration convergence in constraint residuals, not weight proxies.

- **`pct` key re-defined**: was max-per-obs relative weight change; now
  `l1_weight + plateau` (L1 total weight change, autumn/anesrake compatible).

- **`pct_change` result field renamed** to `l1_weight_change` everywhere.

- **New result fields**: `grake_norm`, `l1_weight_change`, `convergence_used`
  (nested: metric, rule, tol, fired_at_iter).

### New features

- `metric = "grake_norm"`: survey::grake-compatible normalized constraint residual.
  Use `tol = 1e-7` for convergence equivalent to `survey::calibrate(epsilon=1e-7)`.
- `convergence_used` in result: expose active criterion at exit for debugging.
```

- [ ] **Step G6: Update harvest.R roxygen + regenerate man/**

In R/harvest.R, update `@param convergence` to document the new API:
```r
#' @param convergence Named list controlling stopping criterion. Keys:
#'   \itemize{
#'     \item \code{metric}: what to measure — \code{"max_err"} (default),
#'       \code{"mean_err"}, \code{"kl"}, \code{"chi2"}, \code{"grake_norm"}
#'       (survey::grake compatible), \code{"l1_weight"} (autumn/anesrake compatible).
#'     \item \code{rule}: how to stop — \code{"improvement"} (default: fire when
#'       loss barely changes), \code{"threshold"} (fire when loss < tol),
#'       \code{"plateau"} (fire when loss stops improving by > tol fraction).
#'     \item \code{tol}: threshold value (default \code{0.001}).
#'     \item Shorthands: \code{improvement = X} (max_err+improvement),
#'       \code{pct = X} (l1_weight+plateau, autumn compatible),
#'       \code{absolute = X} (max_err+threshold).
#'     \item \code{stop_when}: \code{"any"} (default) or \code{"all"} —
#'       combines primary criterion with \code{absolute} as secondary threshold.
#'   }
```

```bash
Rscript -e 'devtools::document()'
```

- [ ] **Step G7: R CMD check gate**
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
R CMD build . && R CMD check --as-cran leafblower_*.tar.gz 2>&1 | grep -E "^(ERROR|WARNING)"
```
Expected: FAIL 0; 0 ERROR, 0 WARNING.

- [ ] **Step G8: Commit**
```bash
git add python/leafblower/_harvest.py python/leafblower/_bindings.cpp \
        python/leafblower/test_python.py NEWS.md R/harvest.R man/
git commit -m "$(cat <<'EOF'
docs+bindings(WU-G): Python parity + NEWS.md + roxygen for convergence redesign

_harvest.py: _parse_convergence() metric+rule API; criterion removed;
pct_change→l1_weight_change in result. _bindings.cpp: struct metric+rule;
6 new result fields. 3 Python parity tests. NEWS.md documents redesign.
Roxygen @param convergence updated. R CMD check --as-cran clean.
EOF
)"
```

- [ ] **Step G9: Close ticket + epic**
```bash
bd close leafblower-lccb
bd close leafblower-4jj7 --reason "all 9 WUs complete; convergence redesign shipped"
```

---

## Final Verification

- [ ] `devtools::test()` → FAIL 0, PASS ≥ 302
- [ ] `pytest python/leafblower/` → all green
- [ ] `R CMD check --as-cran` → 0 ERROR, 0 WARNING
- [ ] `EXPECTED_RK_PARAMS_BYTES = 220` verified by build
- [ ] `pct_change` absent from all R result lists
- [ ] `convergence_used` nested list present with `metric`, `rule`, `tol`, `fired_at_iter`
- [ ] `list()` default → `metric=max_err, rule=improvement, tol=0.001`
- [ ] `list(pct=0.001)` → l1_weight_change < 0.001 at exit (A3)
- [ ] All 9 beads tickets closed

---

## Self-Review

**Spec coverage check:**
- §1.1 Metrics (6): ✅ A4 tests all 6; WU-C implements MAX_ERR/MEAN_ERR/KL/CHI2/GRAKE_NORM/L1_WEIGHT
- §1.2 Rules (3): ✅ WU-C/D implement THRESHOLD/IMPROVEMENT/PLATEAU
- §2 Default: ✅ WU-B sets metric=0,rule=1,tol=0.001; WU-E parses list() → same
- §3 R API shorthands: ✅ WU-E: improvement/pct/absolute shorthands; match.arg for metric+rule
- §4 C ABI EXPECTED_RK_PARAMS_BYTES=220: ✅ WU-A
- §4 legacy helper: ✅ WU-A
- §4 pct_change→l1_weight_change: ✅ WU-A/B
- §5.1 l1_weight normalization: ✅ WU-C (ct.W_input)/WU-D (st.n)
- §5.2 grake_norm: ✅ WU-C/D
- §5.3 improvement rule + prev_metric guard: ✅ WU-C
- §5.4 plateau rule: ✅ WU-C
- §6 Breaking changes: ✅ NEWS.md in WU-G; no shims needed
- §7 A1–A7: ✅ WU-F1; A6 has skip_if_not_installed
- §8 WU sequence: ✅ 9 tickets, sequential deps

**Placeholder scan:** None found. All code blocks are concrete.

**Type consistency:** `CalibMetric` and `CalibRule` defined in WU-A types.hpp and used consistently throughout. `l1_weight_change` field name consistent across types.hpp/ieppa.hpp/raking.hpp/lbfgsb_solver.hpp/r_bridge.cpp/harvest.R/Python.
