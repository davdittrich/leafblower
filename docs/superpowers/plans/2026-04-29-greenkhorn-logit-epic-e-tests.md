# Epic E: R API + Tests — Greenkhorn & Logit

**Mechanism:** Roxygen2 `@param`/`@return` edits + verbatim TDD tests appended to existing test file  
**Forbidden:** Modifying test code from spec; introducing new R functions; touching C++ files  
**Audit:** `grep -c` assertions on generated Rd; `devtools::test()` FAIL count invariant (==3)

**Branch:** `fix/correctness-performance-2026-04-28`  
**Date:** 2026-04-29

---

## Task E1: harvest.R Roxygen updates

**Ticket:** One bead per task — do not bundle with E2.

### Steps

**1. Read R/harvest.R @param method block.**

Target lines: the `@param method` block starting at harvest.R line ~17.
Current text ends at `\code{"grake"} (normalized Chebyshev via IPM).`.

**2. Append to @param method — add greenkhorn and logit entries.**

Insert immediately after `\code{"grake"} (normalized Chebyshev via IPM).`:

```r
#'   \code{"greenkhorn"} (Greedy coordinate-descent IPF — Altschuler-Weed-Rigollet 2017;
#'   picks the single hardest margin per step. Use when K>=2 margins have very unequal
#'   difficulty. Distinct from \code{scheduler="greedy"} on raking, which sorts all K
#'   margins but still sweeps all per round),
#'   \code{"logit"} (Deville-Sarndal 1992 logit-distance Newton calibration; bounds
#'   enforced analytically — no clamping. Use when raking stalls on tight
#'   \code{max_weight} bounds; typically converges in 10-20 Newton steps),
```

**3. Update @param scheduler to clarify Greenkhorn distinction.**

The existing `@param scheduler` reads:
> `(Greenkhorn priority).`

Append a clarifying note:
> `Note: \code{scheduler="greedy"} is NOT the same as \code{method="greenkhorn"} — greedy-scheduler raking still sweeps all K margins per round in residual-sorted order; \code{method="greenkhorn"} performs pure single-margin coordinate descent (one margin per step).`

**4. Update @return — add alm_* field notes for non-ieppa_soft methods.**

The four `alm_*` fields (`alm_capacity_mu_final`, `alm_n_growth_events`,
`alm_max_dual_norm`, `alm_sum_drift`) already exist in the `@return` block.
Append to each: `(\code{0} for \code{method} other than \code{"ieppa_soft"}).`

**5. Regenerate man/harvest.Rd.**

```bash
Rscript -e "devtools::document('/home/dd/Gemini/leafblower')"
```

**6. Verify Rd content.**

```bash
grep -c "greenkhorn" man/harvest.Rd   # must be >= 2
grep -c "logit"      man/harvest.Rd   # must be >= 3
```

**7. Verify test suite FAIL count unchanged.**

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
# Expected: FAIL == 3 (pre-existing failures, unchanged)
```

**8. Commit.**

```
docs(harvest): document greenkhorn and logit methods; API Notes on greedy distinction
```

---

## Task E2: Tests T1-T8

**Ticket:** Separate from E1.

### Prerequisite

E1 committed (harvest.Rd regenerated). Tests are RED before implementation (Epic A–D).

### Step 1: Read end of test file

Read `/home/dd/Gemini/leafblower/tests/testthat/test-calibration-solvers.R` from line 950+
to confirm the exact append point (after T10 block, which ends at the last `})` on line 967).

### Step 2: Append T1–T8 verbatim

Append the following block exactly as written to `tests/testthat/test-calibration-solvers.R`:

```r

# ══════════════════════════════════════════════════════════════════════════════
# Greenkhorn tests (T1–T4)
# ══════════════════════════════════════════════════════════════════════════════

# ── T1: Greenkhorn available, calibrates, returns correct algorithm_used ──────

test_that("T1: greenkhorn available and calibrates", {
  set.seed(1); n <- 1000L
  df  <- data.frame(sex=factor(sample(c("M","F"),n,TRUE)),
                    age=factor(sample(c("Y","O"),n,TRUE)))
  tgt <- list(sex=c(M=0.5,F=0.5), age=c(Y=0.6,O=0.4))
  r   <- harvest(df, tgt, method="greenkhorn", max_iterations=500L)
  expect_lt(attr(r,"result")$max_error, 1e-3)
  expect_equal(attr(r,"result")$algorithm_used, "greenkhorn",
    label="algorithm_used must be 'greenkhorn' (requires alg_names extended to index 9)")
})

# ── T2: Machine precision on trivial 1-margin problem ────────────────────────

test_that("T2: greenkhorn reaches near-machine precision on 1-margin 2-cat", {
  # Greenkhorn on a single 2-cat margin with no bounds converges in 1 step
  # (one exact scale factor). Use absolute convergence to not stop early.
  set.seed(2); n <- 5000L
  df  <- data.frame(g=factor(sample(c("A","B"), n, TRUE, prob=c(0.3,0.7))))
  tgt <- list(g=c(A=0.5, B=0.5))
  r   <- harvest(df, tgt, method="greenkhorn",
                 max_iterations=100L,
                 convergence=list(absolute=1e-12))  # REQUIRED: default pct=1e-4 stops at 1e-4
  me  <- attr(r,"result")$max_error
  # 1 step → exact scale → machine precision
  expect_lt(me, 1e-10,
    label=sprintf("greenkhorn should reach ~machine precision on 1-margin (got %.2e)", me))
})

# ── T3: Quality within 2× of raking on 3-margin problem ─────────────────────

test_that("T3: greenkhorn max_err within 2x of raking", {
  set.seed(42); n <- 10000L
  df <- data.frame(
    a=factor(sample(letters[1:3],n,TRUE)),
    b=factor(sample(LETTERS[1:4],n,TRUE)),
    c=factor(sample(c("x","y"),n,TRUE))
  )
  tgt <- list(a=c(a=0.3,b=0.4,c=0.3),
              b=c(A=0.25,B=0.25,C=0.25,D=0.25),
              c=c(x=0.6,y=0.4))
  r_rk  <- harvest(df, tgt, method="raking",
                   convergence=list(absolute=1e-6))
  r_grk <- harvest(df, tgt, method="greenkhorn",
                   convergence=list(absolute=1e-6))
  me_rk  <- attr(r_rk,  "result")$max_error
  me_grk <- attr(r_grk, "result")$max_error
  expect_lt(me_grk, 2.0 * me_rk + 1e-6)
})

# ── T4: Bounds respected ──────────────────────────────────────────────────────

test_that("T4: greenkhorn respects max_weight and min_weight exactly", {
  set.seed(5); n <- 2000L
  df  <- data.frame(v=factor(sample(5, n, TRUE)))
  tgt <- list(v=setNames(c(0.4,0.3,0.15,0.1,0.05), as.character(1:5)))
  r   <- harvest(df, tgt, method="greenkhorn", max_weight=2.0, min_weight=0.1)
  w   <- r$weights
  expect_true(max(w) <= 2.0 + 1e-9)
  expect_true(min(w) >= 0.1 - 1e-9)
})

# ══════════════════════════════════════════════════════════════════════════════
# Logit calibration tests (T5–T8)
# ══════════════════════════════════════════════════════════════════════════════

# ── T5: Logit available and calibrates ───────────────────────────────────────

test_that("T5: logit available and calibrates", {
  set.seed(5); n <- 1000L
  df  <- data.frame(sex=factor(sample(c("M","F"),n,TRUE)),
                    age=factor(sample(c("Y","O"),n,TRUE)))
  tgt <- list(sex=c(M=0.5,F=0.5), age=c(Y=0.6,O=0.4))
  r   <- harvest(df, tgt, method="logit", max_iterations=50L)
  expect_lt(attr(r,"result")$max_error, 1e-3)
  expect_equal(attr(r,"result")$algorithm_used, "logit")
})

# ── T6: Logit bounds by construction (Newton steps < raking rounds on K=2) ───

test_that("T6: logit bounds by construction (Newton steps < raking rounds on K=2 tight problem)", {
  set.seed(6); n <- 5000L
  df  <- data.frame(
    v=factor(sample(5, n, TRUE)),
    g=factor(sample(c("M","F"), n, TRUE))  # second margin
  )
  tgt <- list(v=setNames(c(0.4,0.3,0.15,0.1,0.05),as.character(1:5)),
              g=c(M=0.55, F=0.45))
  r   <- harvest(df, tgt, method="logit", max_weight=1.5, min_weight=0.1)
  w   <- r$weights
  expect_true(max(w) <= 1.5 + 1e-9)
  expect_true(min(w) >= 0.1 - 1e-9)
  n_iters <- attr(r,"result")$iterations
  expect_lt(n_iters, 50L)
  # K=2 means raking needs >=2 rounds/sweeps; logit Newton should converge faster:
  r_rk <- harvest(df, tgt, method="raking", max_weight=1.5, min_weight=0.1,
                  convergence=list(absolute=1e-6))
  n_rk <- attr(r_rk,"result")$iterations
  expect_lt(n_iters, n_rk,
    label=sprintf("logit Newton (%d steps) < raking (%d rounds) on K=2 tight problem", n_iters, n_rk))
})

# ── T7: Logit matches raking on unconstrained problem (max_err < 1e-4) ───────

test_that("T7: logit and raking reach same calibration target (max_err < 1e-4)", {
  # Without tight bounds, logit and raking should both converge to same weights
  set.seed(7); n <- 5000L
  df  <- data.frame(
    a=factor(sample(letters[1:3],n,TRUE)),
    b=factor(sample(LETTERS[1:4],n,TRUE))
  )
  tgt <- list(a=c(a=0.3,b=0.4,c=0.3), b=c(A=0.25,B=0.25,C=0.25,D=0.25))
  r_logit <- harvest(df, tgt, method="logit",
                     convergence=list(absolute=1e-6))
  me_logit <- attr(r_logit,"result")$max_error
  expect_lt(me_logit, 1e-4)
})

# ── T8: Logit vs raking on tight-bounds problem (logit must be competitive) ───

test_that("T8: logit max_err within 2x of raking on tight-bounds problem", {
  set.seed(8); n <- 5000L
  df  <- data.frame(v=factor(sample(5, n, TRUE)))
  tgt <- list(v=setNames(c(0.4,0.3,0.15,0.1,0.05), as.character(1:5)))
  r_rk    <- harvest(df, tgt, method="raking", max_weight=1.8, min_weight=0,
                     convergence=list(absolute=1e-6))
  r_logit <- harvest(df, tgt, method="logit", max_weight=1.8, min_weight=0,
                     convergence=list(absolute=1e-6))
  me_rk    <- attr(r_rk,    "result")$max_error
  me_logit <- attr(r_logit, "result")$max_error
  expect_lt(me_logit, 2.0 * me_rk + 1e-6)
})

# ── T_acc: Greenkhorn with accelerate=TRUE (SQUAREM round-level acceleration) ─

test_that("Tacc: greenkhorn with accelerate=TRUE runs without error and converges", {
  set.seed(99); n <- 2000L
  df  <- data.frame(x=factor(sample(letters[1:3],n,TRUE)),
                    y=factor(sample(c("M","F"),n,TRUE)))
  tgt <- list(x=c(a=0.3,b=0.4,c=0.3), y=c(M=0.5,F=0.5))
  r_acc   <- harvest(df, tgt, method="greenkhorn", accelerate=TRUE,  max_iterations=500L)
  r_plain <- harvest(df, tgt, method="greenkhorn", accelerate=FALSE, max_iterations=500L)
  me_acc   <- attr(r_acc,   "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  # Both must converge to good quality
  expect_lt(me_acc, 1e-3)
  # Accelerated must not massively regress in iterations vs plain
  iters_acc   <- attr(r_acc,   "result")$iterations
  iters_plain <- attr(r_plain, "result")$iterations
  expect_lt(iters_acc, iters_plain * 2L,
    label=sprintf("accelerate=TRUE used %d rounds vs plain %d; must not be >2x worse",
                  iters_acc, iters_plain))
})
```

### Step 3: Verify tests RED (before implementation)

```bash
Rscript -e "devtools::test(filter='calibration-solvers')" 2>&1 | tail -20
# T1-T8 must ERROR (method not found) or FAIL — not PASS
```

### Step 4: Post-implementation — verify all 8 PASS

After Epics A–D are implemented and built:

```bash
R CMD INSTALL --preclean /home/dd/Gemini/leafblower
Rscript -e "devtools::test(filter='calibration-solvers')" 2>&1 | tail -20
# T1 PASS: max_error < 1e-3; algorithm_used == "greenkhorn"
# T2 PASS: max_error < 1e-10 (convergence=list(absolute=1e-12))
# T3 PASS: greenkhorn max_err < 2 * raking max_err + 1e-6
# T4 PASS: max(w) <= 2.0 + 1e-9; min(w) >= 0.1 - 1e-9
# T5 PASS: max_error < 1e-3; algorithm_used == "logit"
# T6 PASS: bounds respected; n_iters < 50 AND < n_rk
# T7 PASS: logit max_err < 1e-4
# T8 PASS: logit max_err < 2 * raking max_err + 1e-6
# Tacc PASS: greenkhorn accelerate=TRUE max_err < 1e-3; iters < 2 * plain iters
```

### Step 5: Verify FAIL count invariant

```bash
Rscript -e "devtools::test()" 2>&1 | tail -3
# FAIL == 3 (pre-existing failures unchanged)
```

### Step 6: Commit

```
test: add T1-T8 Greenkhorn and Logit calibration tests from spec
```

---

## Critical invariants

| Invariant | Why |
|-----------|-----|
| T2 MUST use `convergence=list(absolute=1e-12)` | Default `pct=1e-4` plateau-stops before reaching 1e-10 |
| T6 MUST be K=2 (v + g margin) | Ensures `n_rk > n_iters_logit` is reliable; K=1 raking also converges in 1 step |
| `alg_names` in harvest.R MUST extend to 11 elements | Missing index 9/10 → `algorithm_used` returns `""` → T1/T5 fail |
| `harvest.R stop()` for `status==2` MUST use `res$message` | Hardcoded override discards logit's "singular normal equations" message |
| `map_method()` MUST include `"greenkhorn"` and `"logit"` | Missing → harvest() rejects method before C++ dispatch → all T1/T5 fail |

## Acceptance criteria

| # | Criterion | Gate |
|---|-----------|------|
| AC-E1a | `grep -c "greenkhorn" man/harvest.Rd` >= 2 | E1 commit |
| AC-E1b | `grep -c "logit" man/harvest.Rd` >= 3 | E1 commit |
| AC-E1c | `devtools::test()` FAIL == 3 | E1 commit |
| AC-E2a | T1–T4 PASS (greenkhorn) | Post A–D build |
| AC-E2b | T5–T8 PASS (logit) | Post A–D build |
| AC-E2c | T_acc PASS (greenkhorn accelerate=TRUE) | Post A–D build |
| AC-E2d | `devtools::test()` FAIL == 3 | Post A–D build |
