# Bounded Convergence Fix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix two convergence bugs that prevent leafblower from reaching 1e-6 tolerance on problems where `max_weight` binds: iEPPA's BCD-with-clamp replaced with Dykstra's alternating projections, and L-BFGS-B's misrouted link function corrected to logit for finite `max_weight`.

**Architecture:** iEPPA gets a full solver-body replacement — Dykstra's alternating projections with K+1 correction vectors replaces `bcd_sweep` + EPPA outer loop + post-convergence projection. L-BFGS-B gets a one-line fix in the link-function selector and a matching guard update in the C API validator. `compute_errRp` is retained unchanged.

**Tech Stack:** C++17, R testthat v3, devtools, `R CMD INSTALL --preclean .` as the compile gate.

---

## File Map

| File | Action | What changes |
|------|--------|-------------|
| `tests/testthat/test-bounded-convergence.R` | **Create** | Two RED tests written first (TDD) |
| `src/logit.hpp` | **Modify** | Line 19: `exponential = !std::isfinite(U)` |
| `src/c_api.cpp` | **Modify** | Line 61: `bool use_logit = std::isfinite(p->max_weight);` |
| `src/ieppa.cpp` | **Modify** | Delete `bregman_dist` + `bcd_sweep`; replace `ieppa_solve` body with Dykstra's |

Existing tests: 6 test files across `test-compat.R`, `test-design.R`, `test-harvest.R`, `test-ieppa.R`, `test-lbfgsb.R`, `test-logit.R`. All must remain GREEN after every task. Do not rely on any hardcoded test count — use `devtools::test()` output as the authority.

---

### Task 1: Write RED Tests (TDD)

**Files:**
- Create: `tests/testthat/test-bounded-convergence.R`

- [ ] **Step 1: Create beads tracking issue**

```bash
bd create \
  --title="bounded-convergence: write RED tests" \
  --description="TDD: write failing tests before fixing bugs. Both tests must FAIL before implementation starts." \
  --type=task --priority=1
# Note the issue ID printed; claim it:
bd update <id> --claim
```

- [ ] **Step 2: Create the test file**

Create `tests/testthat/test-bounded-convergence.R` with this exact content:

```r
test_that("iEPPA converges on bound-hitting problem (max_weight=5, skewed sample)", {
  set.seed(42); n <- 10000L
  df <- data.frame(
    age = factor(sample(c("18-34","35-54","55+"), n, replace=TRUE,
                        prob=c(0.60,0.30,0.10))),
    sex = factor(sample(c("M","F"), n, replace=TRUE, prob=c(0.70,0.30))),
    edu = factor(sample(c("HS","Some","BA","Grad"), n, replace=TRUE,
                        prob=c(0.50,0.25,0.15,0.10)))
  )
  tgt <- list(age=c("18-34"=0.33,"35-54"=0.40,"55+"=0.27),
              sex=c(M=0.49,F=0.51),
              edu=c(HS=0.28,Some=0.30,BA=0.27,Grad=0.15))
  result <- harvest(df, tgt, method="ieppa", max_weight=5)
  expect_true(max(result$weights) <= 5.0 + 1e-10)
  diag <- diagnose_weights(result, tgt, result$weights)
  expect_true(all(abs(diag$error_weighted) < 1e-6))
})

test_that("L-BFGS-B converges on bound-hitting problem (max_weight=5, skewed sample)", {
  set.seed(42); n <- 2000L
  df <- data.frame(
    age = factor(sample(c("18-34","35-54","55+"), n, replace=TRUE,
                        prob=c(0.60,0.30,0.10))),
    sex = factor(sample(c("M","F"), n, replace=TRUE, prob=c(0.70,0.30)))
  )
  tgt <- list(age=c("18-34"=0.33,"35-54"=0.40,"55+"=0.27),
              sex=c(M=0.49,F=0.51))
  expect_no_warning(
    result <- harvest(df, tgt, method="lbfgsb", max_weight=5)
  )
  diag <- diagnose_weights(result, tgt, result$weights)
  expect_true(all(abs(diag$error_weighted) < 1e-6))
})
```

- [ ] **Step 3: Run tests — verify they FAIL**

```bash
cd /home/dd/Gemini/leafblower
Rscript -e "devtools::test(filter='bounded-convergence')"
```

Expected: **Both tests FAIL.** The iEPPA test fails on `all(abs(diag$error_weighted) < 1e-6)` (actual ~2.3e-3). The L-BFGS-B test either warns or fails on the convergence check (actual ~1.7e-2).

**If either test PASSES at this stage, do not proceed.** The bug is not reproducible with the exact seed/n — investigate before continuing.

- [ ] **Step 4: Commit RED tests**

```bash
git add tests/testthat/test-bounded-convergence.R
git commit -m "test(ieppa,lbfgsb): add RED bounded-convergence tests

Two failing tests for bound-hitting convergence at max_weight=5 with skewed
sample (60/30/10 age, 70/30 sex). Both fail until Dykstra's iEPPA fix and
L-BFGS-B logit-link fix are applied."
bd close <id>
```

---

### Task 2: Fix L-BFGS-B — Logit Link for Finite max_weight

**Files:**
- Modify: `src/logit.hpp:19`
- Modify: `src/c_api.cpp:61`

**Why this fixes the bug:** `LinkFn` selects the exp link whenever `L == 0.0` (the default `min_weight`), even when `max_weight` is finite. Exp link maps the dual variable to `(0, ∞)`, so `max_weight` is applied as a post-hoc clamp that breaks dual-primal correspondence — the gradient at the clamped point is not zero, so the solver never converges. With `L=0, U=5`, the logit link gives `logit_scale = (5-0)/((5-1)*(1-0)) = 1.25` and `F(u) = 5·exp(1.25u)/(4+exp(1.25u))`, which maps the dual exactly to `(0, 5)` — making the clamp a no-op and restoring gradient correctness.

- [ ] **Step 1: Create beads tracking issue**

```bash
bd create \
  --title="lbfgsb: fix logit link selection for finite max_weight" \
  --description="Fix LinkFn::exponential condition: L==0.0 should not force exp link when U is finite." \
  --type=bug --priority=1
bd update <id> --claim
```

- [ ] **Step 2: Fix src/logit.hpp line 19**

Current line 19 (inside `explicit LinkFn(double min_weight, double max_weight)` constructor):
```cpp
        exponential = (L == 0.0 || !std::isfinite(U));
```

Replace with:
```cpp
        exponential = !std::isfinite(U);
```

The full constructor after the fix (lines 15–21 of logit.hpp):
```cpp
    explicit LinkFn(double min_weight, double max_weight)
        : L(min_weight), U(max_weight)
    {
        // Use exponential link only when max_weight is infinite
        exponential = !std::isfinite(U);
        logit_scale = exponential ? 0.0 : (U - L) / ((U - 1.0) * (1.0 - L));
    }
```

- [ ] **Step 3: Fix src/c_api.cpp line 61**

Current line 61 (inside `validate_inputs`, under `// Logit singularity guard`):
```cpp
    bool use_logit = (p->min_weight > 0.0) && std::isfinite(p->max_weight);
```

Replace with:
```cpp
    bool use_logit = std::isfinite(p->max_weight);
```

The `L=0, U=1.0` edge case (logit_scale denominator `(U-1)=0`) is already caught by the existing guard on lines 65–66:
```cpp
        if (p->max_weight == 1.0)
            return err("logit link undefined: max_weight=1 makes denominator (U-1)=0");
```
No other change needed in that block.

- [ ] **Step 4: Compile gate**

```bash
cd /home/dd/Gemini/leafblower
R CMD INSTALL --preclean .
```

Expected: install completes with no errors. **If compilation fails, do not proceed.** Read the exact error message and fix the root cause.

- [ ] **Step 5: Run bounded-convergence tests**

```bash
Rscript -e "devtools::test(filter='bounded-convergence')"
```

Expected: **L-BFGS-B test GREEN; iEPPA test still RED** (iEPPA not yet fixed). If L-BFGS-B test still fails: verify the edit was applied (`grep -n 'exponential' src/logit.hpp` should show `exponential = !std::isfinite(U);`).

- [ ] **Step 6: Run full test suite — verify no regressions**

```bash
Rscript -e "devtools::test()"
```

Expected: All existing test blocks GREEN (1 new iEPPA block still RED is acceptable). Any newly failing existing test is a regression — investigate before continuing.

- [ ] **Step 7: Commit**

```bash
git add src/logit.hpp src/c_api.cpp
git commit -m "fix(lbfgsb): use logit link for finite max_weight when min_weight=0

LinkFn selected exp link whenever L==0.0, even with finite U. Exp link maps
the dual to (0,inf), so max_weight was applied as a post-hoc clamp that broke
dual-primal correspondence and prevented convergence.

Fix: exponential = !std::isfinite(U). With L=0, U=5: logit_scale=1.25 and
F(u)=5*exp(1.25u)/(4+exp(1.25u)) maps exactly to (0,5); clamp becomes no-op.

Matching use_logit guard in c_api.cpp updated to std::isfinite(p->max_weight);
existing singularity checks for max_weight==1 remain correct."
bd close <id>
```

---

### Task 3: Fix iEPPA — Dykstra's Alternating Projections

**Files:**
- Modify: `src/ieppa.cpp` (full rewrite of solver body; `compute_errRp` kept verbatim)

**Why Dykstra's, not BCD continuation:** BCD-with-clamp cycles near the constraint boundary — the clamp after each IPF step shifts the iterate off the marginal constraint set, so the next margin starts from a biased point. Dykstra's alternating projections maintains K+1 correction vectors (`p[k]` per margin, `q` for the box) that accumulate the "overshoot" from each projection. This provably forces convergence to the exact intersection of all constraint sets, with monotone error decrease from iteration 1. No outer EPPA loop needed.

**What is deleted:**
- `bregman_dist` function (entire)
- `bcd_sweep` function (entire)
- The outer EPPA loop (`for (int outer = 1; ...)`)
- The post-convergence projection loop (`for (int proj = 0; proj < 200; ...)`)

`inner_max_iter` (default 500 in `rk_params_init`) becomes the single iteration budget. `outer_max_iter` is preserved in `CalibState` for API compatibility but is not read inside `ieppa_solve`.

**Memory:** K+1 extra vectors of size n. At n=10000, K=3: ~320 KB. Acceptable.

- [ ] **Step 1: Create beads tracking issue**

```bash
bd create \
  --title="ieppa: replace BCD-with-clamp with Dykstra's alternating projections" \
  --description="Delete bregman_dist, bcd_sweep, EPPA outer loop, post-projection loop. Replace ieppa_solve body with Dykstra's. Keep compute_errRp verbatim." \
  --type=bug --priority=1
bd update <id> --claim
```

- [ ] **Step 2: Replace src/ieppa.cpp**

Write the following as the complete contents of `src/ieppa.cpp`. This replaces everything except the header includes and the `compute_errRp` function, which is kept verbatim.

```cpp
#include "ieppa.hpp"
#include "leafblower.h"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <vector>

namespace lbw {

// Compute errRp = max_k max_j |S_kj/W - tau_kj|
// O(n*K): single O(n) bucket accumulation pass per margin.
static double compute_errRp(const CalibState& st,
                              const std::vector<double>& w) {
    double W = 0.0;
    for (int i = 0; i < st.n; i++) W += w[i];

    double err = 0.0;
    for (int k = 0; k < st.K; k++) {
        std::vector<double> bucket(st.cat_counts[k], 0.0);
        for (int i = 0; i < st.n; i++) {
            int g = st.group_ids[k][i];
            if (g >= 0) bucket[g] += w[i];
        }
        for (int j = 0; j < st.cat_counts[k]; j++) {
            double e = std::fabs(bucket[j] / W - st.targets[k][j]);
            if (e > err) err = e;
        }
    }
    return err;
}

// Dykstra's alternating projections for constrained raking.
// Finds w in the intersection of K marginal constraint sets {C_k} and box [lo,hi]^n.
// Correction vectors p[k] (marginal) and q (box) prevent cycling at constraint boundary.
// inner_max_iter is the single iteration budget; outer_max_iter is unused.
IEPPAResult ieppa_solve(CalibState& st) {
    IEPPAResult res;
    res.status = RK_ERR_NOCONV;
    res.iterations = 0;
    res.max_error = 1.0;

    std::vector<double> w(st.weights, st.weights + st.n);

    // K marginal correction vectors + 1 box correction vector (all zero-initialized)
    std::vector<std::vector<double>> p(st.K, std::vector<double>(st.n, 0.0));
    std::vector<double> q(st.n, 0.0);

    double lo = st.min_weight;
    double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    bool infeas_flag = false;

    for (int iter = 1; iter <= st.inner_max_iter; iter++) {
        res.iterations = iter;

        // Marginal projections: IPF step with no clamp, Dykstra correction applied
        for (int k = 0; k < st.K; k++) {
            // Apply correction to get corrected input y
            std::vector<double> y(st.n);
            for (int i = 0; i < st.n; i++) y[i] = w[i] + p[k][i];

            // Bucket accumulation for IPF scale computation
            double W = 0.0;
            std::vector<double> bucket(st.cat_counts[k], 0.0);
            for (int i = 0; i < st.n; i++) {
                W += y[i];
                int g = st.group_ids[k][i];
                if (g >= 0) bucket[g] += y[i];
            }

            // IPF scale factors
            std::vector<double> scale(st.cat_counts[k], 1.0);
            for (int j = 0; j < st.cat_counts[k]; j++) {
                double Tkj = st.targets[k][j] * W;
                if (bucket[j] < 1e-15 * W) {
                    if (Tkj > 0.0) infeas_flag = true;
                } else {
                    scale[j] = Tkj / bucket[j];
                }
            }

            // IPF step — NO CLAMP. g==-1 (NA) entries pass through unchanged.
            // Update Dykstra correction: p[k][i] = y[i] - projection(y[i])
            for (int i = 0; i < st.n; i++) {
                int g = st.group_ids[k][i];
                double w_new = (g >= 0) ? y[i] * scale[g] : y[i];
                p[k][i] = y[i] - w_new;
                w[i] = w_new;
            }
        }

        // Box projection [lo, hi]^n with Dykstra correction
        for (int i = 0; i < st.n; i++) {
            double yi = w[i] + q[i];
            double wc = std::max(lo, std::min(hi, yi));
            q[i] = yi - wc;
            w[i] = wc;
        }

        // Convergence check
        double errRp = compute_errRp(st, w);
        res.max_error = errRp;

        if (st.verbose >= 1) {
            char msg[256];
            std::snprintf(msg, 256, "Dykstra iter %d: errRp=%.2e", iter, errRp);
            st.log(msg);
        }

        if (errRp < st.tol_abs) {
            res.status = infeas_flag ? RK_ERR_INFEAS : RK_OK;
            break;
        }
    }

    if (infeas_flag && res.status == RK_OK) res.status = RK_ERR_INFEAS;

    for (int i = 0; i < st.n; i++) st.weights[i] = w[i];

    return res;
}

} // namespace lbw
```

- [ ] **Step 3: Compile gate**

```bash
cd /home/dd/Gemini/leafblower
R CMD INSTALL --preclean .
```

Expected: successful install. If it fails with any undefined reference or type error, read the exact compiler output. Do not proceed with a broken build.

- [ ] **Step 4: Run bounded-convergence tests**

```bash
Rscript -e "devtools::test(filter='bounded-convergence')"
```

Expected: **Both tests GREEN.** Both `max(result$weights) <= 5.0 + 1e-10` and `all(abs(diag$error_weighted) < 1e-6)` pass for both iEPPA and L-BFGS-B.

If iEPPA test still fails: run with `verbose=1` to trace the errRp trajectory:
```bash
Rscript -e "
library(leafblower)
set.seed(42); n <- 10000L
df <- data.frame(
  age=factor(sample(c('18-34','35-54','55+'),n,replace=TRUE,prob=c(.60,.30,.10))),
  sex=factor(sample(c('M','F'),n,replace=TRUE,prob=c(.70,.30))),
  edu=factor(sample(c('HS','Some','BA','Grad'),n,replace=TRUE,prob=c(.50,.25,.15,.10)))
)
tgt <- list(age=c('18-34'=.33,'35-54'=.40,'55+'=.27),
            sex=c(M=.49,F=.51),
            edu=c(HS=.28,Some=.30,BA=.27,Grad=.15))
r <- harvest(df, tgt, method='ieppa', max_weight=5, verbose=1)
cat('max_error:', attr(r, 'max_error'), '\n')
"
```
If errRp decreases monotonically but stops above 1e-6, increase `max_iterations` to diagnose whether it's a budget issue.

- [ ] **Step 5: Run full test suite**

```bash
Rscript -e "devtools::test()"
```

Expected: All test blocks GREEN (existing + 2 new). Pay particular attention to:
- `test-ieppa.R` "iEPPA respects max_weight=2 on tight bounds" — exercises box constraint on different data
- `test-ieppa.R` "iEPPA respects min_weight=0.5" — exercises lo > 0 in box projection
- `test-logit.R` — unaffected (logit math unchanged); must remain GREEN

- [ ] **Step 6: Commit**

```bash
git add src/ieppa.cpp
git commit -m "fix(ieppa): replace BCD-with-clamp with Dykstra's alternating projections

bcd_sweep clamped weights during IPF scaling, violating the Sinkhorn invariant
and causing cycling near constraint boundaries (~2.3e-3 residual after 2000 iters).
The post-convergence projection loop (formerly lines 132-153) degraded errRp
after the convergence criterion fired.

Dykstra's alternating projections uses K+1 correction vectors (p[k] per margin,
q for box) to provably converge to the intersection of all constraint sets.
IPF steps are clamp-free; box constraint is a separate Dykstra projection step.

Deleted: bcd_sweep, bregman_dist, EPPA outer loop, post-projection loop.
inner_max_iter (default 500) is the single iteration budget.
outer_max_iter field preserved in CalibState for API compatibility; unused."
bd close <id>
```

---

### Task 4: Final Verification Gate

**Files:** None (read-only verification)

- [ ] **Step 1: Create beads tracking issue**

```bash
bd create \
  --title="bounded-convergence: final verification gate" \
  --description="Run full test suite and benchmark to confirm both bugs fixed, no regressions, success criteria met." \
  --type=task --priority=1
bd update <id> --claim
```

- [ ] **Step 2: Run full test suite**

```bash
cd /home/dd/Gemini/leafblower
Rscript -e "devtools::test()"
```

Expected: All test blocks GREEN (existing + 2 new). Record the output. Any failure is a blocker.

- [ ] **Step 3: Verify convergence on the benchmark that previously failed (warnings asserted)**

```bash
Rscript -e "
library(leafblower)
# options(warn=2) converts any warning to an error — failure here means a warning was emitted
options(warn=2)
set.seed(42); n <- 10000L
df <- data.frame(
  age = factor(sample(c('18-34','35-54','55+'), n, replace=TRUE, prob=c(0.60,0.30,0.10))),
  sex = factor(sample(c('M','F'), n, replace=TRUE, prob=c(0.70,0.30))),
  edu = factor(sample(c('HS','Some','BA','Grad'), n, replace=TRUE, prob=c(0.50,0.25,0.15,0.10)))
)
tgt <- list(age=c('18-34'=0.33,'35-54'=0.40,'55+'=0.27),
            sex=c(M=0.49,F=0.51),
            edu=c(HS=0.28,Some=0.30,BA=0.27,Grad=0.15))
r_ep <- harvest(df, tgt, method='ieppa', max_weight=5)
r_lb <- harvest(df, tgt, method='lbfgsb', max_weight=5)
diag_ep <- diagnose_weights(r_ep, tgt, r_ep\$weights)
diag_lb <- diagnose_weights(r_lb, tgt, r_lb\$weights)
ep_err <- max(abs(diag_ep\$error_weighted))
lb_err <- max(abs(diag_lb\$error_weighted))
cat('iEPPA    max|error_weighted|:', ep_err, '\n')
cat('L-BFGS-B max|error_weighted|:', lb_err, '\n')
stopifnot(ep_err < 1e-6)
stopifnot(lb_err < 1e-6)
cat('All success criteria met.\n')
"
```

Expected output ends with `All success criteria met.` Any warning emitted by `harvest()` will raise an error (via `options(warn=2)`) before `stopifnot` is reached. If the script exits non-zero, check the error message.

Success criteria (from spec):
- Script exits 0 with `All success criteria met.` printed
- `iEPPA max|error_weighted|` < 1e-6 (enforced by `stopifnot`)
- `L-BFGS-B max|error_weighted|` < 1e-6 (enforced by `stopifnot`)
- No convergence warnings emitted (enforced by `options(warn=2)`)

- [ ] **Step 4: Close tracking**

```bash
bd close <id>
```
