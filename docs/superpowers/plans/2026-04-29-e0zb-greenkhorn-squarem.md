# Fix leafblower-e0zb: greenkhorn accelerate=TRUE silently ignores SQUAREM

> **For agentic workers:** Use superpowers:subagent-driven-development to execute.

**Goal:** Allow `harvest(method="greenkhorn", accelerate=TRUE)` to pass SQUAREM through to greenkhorn_solve, and scope harvest.R's R8 raking/SQUAREM mutual-exclusion guard so it does not apply to greenkhorn.

**Mechanism:** Extend `accelerate_bool` from `method == "raking"` to `method %in% c("raking","greenkhorn")`; scope R8 greedy guard to `method == "raking"` only; rewrite T_acc to verify SQUAREM fires using K*3 divisibility.

**Forbidden:** Touching greenkhorn.cpp; touching logit; adding new harvest() params; amending prior commits.

**Audit:** `grep accelerate_bool R/harvest.R` → must show `%in% c("raking","greenkhorn")`; T_acc checks `iters_acc %% 6 == 0` (K=2×3=6 per super-step proves SQUAREM fired).

---

## Pre-flight

Read `R/harvest.R` lines 275-300 to confirm exact current state:

- Line 277-280: accelerate_bool gated to `method == "raking"` → CONFIRMED
- Line 297-298: R8 guard `if (accelerate_bool && scheduler == "greedy")` → CONFIRMED

Baseline: `Rscript -e "devtools::test()" 2>&1 | tail -3` → FAIL 3 | PASS 479

---

## Task 1: R/harvest.R — extend accelerate-allowed list + scope R8 guard

**Step 1: Replace lines 277-280**

Before (current):
```r
if (isTRUE(accelerate) && method != "raking")
  warning("accelerate=TRUE is only supported for method='raking'; ignoring for method='",
          method, "'")
accelerate_bool <- isTRUE(accelerate) && method == "raking"
```

After:
```r
if (isTRUE(accelerate) && !method %in% c("raking", "greenkhorn"))
  warning("accelerate=TRUE is only supported for method='raking' or 'greenkhorn'; ignoring for method='",
          method, "'", call. = FALSE)
accelerate_bool <- isTRUE(accelerate) && method %in% c("raking", "greenkhorn")
```

**Step 2: Replace lines 297-298**

Before (current):
```r
# Greedy inside SQUAREM changes the fixed point each F-call, degrading CBB accuracy.
if (accelerate_bool && scheduler == "greedy") scheduler <- "round_robin"
```

After:
```r
# R8: raking-only guard. Greedy re-sorts within each F_eval for raking (non-stationary).
# greenkhorn sorts ONCE at F_eval entry (stationary); scheduler param irrelevant to it.
if (accelerate_bool && method == "raking" && scheduler == "greedy")
  scheduler <- "round_robin"
```

**Step 3: Compile gate**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
# Must: * DONE (leafblower)
```

**Step 4: Smoke test — SQUAREM fires**

```bash
Rscript -e "
library(leafblower)
set.seed(99); n <- 2000L
df  <- data.frame(x=factor(sample(letters[1:3],n,TRUE)),
                  y=factor(sample(c('M','F'),n,TRUE)))
tgt <- list(x=c(a=0.3,b=0.4,c=0.3), y=c(M=0.5,F=0.5))
r_acc   <- suppressWarnings(harvest(df,tgt,method='greenkhorn',accelerate=TRUE,  max_iterations=500L,attach_weights=FALSE))
r_plain <- suppressWarnings(harvest(df,tgt,method='greenkhorn',accelerate=FALSE, max_iterations=500L,attach_weights=FALSE))
i_acc   <- attr(r_acc,  'result')\$iterations
i_plain <- attr(r_plain,'result')\$iterations
cat('plain iters:', i_plain, '  squarem iters:', i_acc, '  squarem%%6=', i_acc%%6L, '\n')
" 2>&1 | grep "plain iters"
# Expected: squarem iters divisible by 6 (K=2 × 3 F_eval = 6 per super-step)
# Expected: plain iters NOT divisible by 6
```

**Step 5: Full test suite**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
# Must: FAIL 3 (unchanged)
```

**Step 6: Commit**

```bash
git add R/harvest.R
git commit -m "fix(harvest): allow accelerate=TRUE for greenkhorn; scope R8 guard to raking

harvest.R unconditionally set accelerate_bool=FALSE for any method other
than 'raking', silently ignoring SQUAREM for greenkhorn despite its solver
having a complete SQUAREM implementation (F_eval = K greedy steps sorted
once at entry — stationary operator, valid for CBB extrapolation).

Two targeted changes:
1. Extend accelerate-allowed list: method %in% c('raking','greenkhorn')
2. Scope R8 greedy/SQUAREM guard to method=='raking' only:
   greenkhorn has its own internal argmax selection; harvest.R scheduler
   param is irrelevant to it; R8 guard must not apply.

Fixes: leafblower-e0zb"
```

---

## Task 2: tests/testthat/test-calibration-solvers.R — rewrite T_acc

**Current T_acc problem:** `expect_lt(iters_acc, iters_plain * 2L)` is semantically broken after the fix. SQUAREM counts K*3 per super-step (K=2 → 6 per iter); plain counts 1. Cross-mode count comparison is meaningless and will produce misleading failures or spurious passes.

**Step 1: Find the Tacc test** (search for `test_that("Tacc:`)

**Step 2: Replace it with:**

```r
test_that("Tacc: greenkhorn accelerate=TRUE fires SQUAREM and converges", {
  set.seed(99); n <- 2000L
  K_exp <- 2L  # K=2 margins in this problem
  df  <- data.frame(x=factor(sample(letters[1:3],n,TRUE)),
                    y=factor(sample(c("M","F"),n,TRUE)))
  tgt <- list(x=c(a=0.3,b=0.4,c=0.3), y=c(M=0.5,F=0.5))

  r_acc   <- suppressWarnings(
    harvest(df, tgt, method="greenkhorn", accelerate=TRUE,  max_iterations=500L,
            attach_weights=FALSE))
  r_plain <- suppressWarnings(
    harvest(df, tgt, method="greenkhorn", accelerate=FALSE, max_iterations=500L,
            attach_weights=FALSE))

  me_acc   <- attr(r_acc,   "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  iters_acc   <- attr(r_acc,   "result")$iterations
  iters_plain <- attr(r_plain, "result")$iterations

  # (1) Both must converge to good quality
  expect_lt(me_acc, 1e-3,
    label=sprintf("greenkhorn+SQUAREM must converge: got max_err=%.2e", me_acc))
  expect_lt(me_plain, 1e-3,
    label=sprintf("greenkhorn plain must converge: got max_err=%.2e", me_plain))

  # (2) SQUAREM fired: iters divisible by K*3 = 6 (K=2 margins x 3 F_eval/super-step)
  # Plain path counts 1 per main loop iter; SQUAREM counts K*3. Divisibility proves
  # the SQUAREM branch executed, not the plain single-step branch.
  squarem_stride <- K_exp * 3L
  expect_equal(iters_acc %% squarem_stride, 0L,
    label=sprintf(
      "greenkhorn+SQUAREM iters (%d) must be divisible by K*3=%d (confirms SQUAREM fired)",
      iters_acc, squarem_stride))

  # (3) Accelerated and plain produce qualitatively similar results
  expect_lt(abs(me_acc - me_plain), 5e-4,
    label=sprintf("SQUAREM and plain must converge similarly (acc=%.2e, plain=%.2e)",
                  me_acc, me_plain))
})
```

**Step 3: Verify**

```bash
Rscript -e "devtools::test(filter='calibration-solvers')" 2>&1 | grep -E "Tacc|PASS|FAIL" | head -5
Rscript -e "devtools::test()" 2>&1 | tail -3
# Must: FAIL 3
```

**Step 4: Commit**

```bash
git add tests/testthat/test-calibration-solvers.R
git commit -m "test(calibration-solvers): rewrite Tacc to verify SQUAREM fires via K*3 divisibility

Prior Tacc compared iters_acc < iters_plain * 2 — semantically broken after
harvest.R fix because SQUAREM path counts K*3 per super-step while plain
counts 1. Cross-mode count comparison is meaningless.

New assertions:
1. Both accelerated and plain converge: max_error < 1e-3
2. SQUAREM fired: iters_acc %% K*3 == 0 (K=2 -> stride=6)
3. Quality within 5e-4 of each other"
```

---

## Acceptance Criteria

| # | Criterion | Verify |
|---|-----------|--------|
| AC1 | `accelerate_bool` allows greenkhorn | `grep "accelerate_bool" R/harvest.R` shows `%in% c("raking","greenkhorn")` |
| AC2 | R8 guard scoped to raking | `grep "scheduler.*greedy" R/harvest.R` shows `method == "raking"` in condition |
| AC3 | T_acc PASS: `iters_acc %% 6 == 0` | `devtools::test(filter='calibration-solvers')` |
| AC4 | No regression: FAIL==3 | `devtools::test()` |
| AC5 | Stepstone: squarem iters %% 27 == 0 | Manual benchmark run (K=9 × 3 = 27) |
