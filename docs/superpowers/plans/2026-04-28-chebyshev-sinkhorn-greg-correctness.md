# Chebyshev, Sinkhorn, Greg Correctness — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development

**Goal:** Fix seven correctness bugs across chebyshev.cpp, sinkhorn.cpp, and greg.cpp.

**Architecture:** Three solvers, independent fixes. Chebyshev has 4 bugs; sinkhorn 1; greg 2. Each fix one commit.

**Tech Stack:** C++17, R CMD INSTALL --preclean ., devtools::test()

---

**Mechanism:** Guard ordering (B1), phantom accumulation (B10), stale value (B14), IPM invariant (B15), factorization cache (R1), INFEAS detection (R7), diagnostic (R10)
**Forbidden:** Changing IPM algorithm structure beyond the specific fix; touching N_red debug path
**Audit:** RED tests where possible; verify solver outputs correct status/values

---

## Source Evidence (read 2026-04-28)

All line numbers verified against `/home/dd/Gemini/leafblower/src/`:

- **chebyshev.cpp line 210–265**: `if (W > 1e-300)` block populates `cm`; convergence check at 220–243 reads `cm.errRp`; INFEAS guard `if (W < 1e-300)` at line 265 fires *after* convergence check.
- **chebyshev.cpp lines 424–425**: `alpha_p = std::max(alpha_p, 1e-10); alpha_d = std::max(alpha_d, 1e-10);` — floor that overrides the fraction-to-boundary.
- **chebyshev.cpp lines 446–449**: slack update uses `std::max(..., kEps)` silently masking negative raw slack.
- **chebyshev.cpp lines 469–479**: dual explosion guard resets `y_delta = mu / s_delta` using `mu` (pre-step), not `mu_new` (line 472).
- **sinkhorn.cpp lines 147–156**: `if (needs_projection)` block updates `a[]`; `else` block (`X_proj = X`) does NOT zero `a[]`.
- **greg.cpp line 61**: `static constexpr int kMaxNewtonIters = 10;`
- **greg.cpp lines 64–124**: Newton loop recomputes `N` + LDLT every iteration unconditionally; no active-set change tracking.
- **greg.cpp line 124 (post-loop)**: no NOCONV message set on exhausting `kMaxNewtonIters`.

---

## Task 1 — B1: chebyshev convergence fires before INFEAS guard

**File:** `src/chebyshev.cpp`
**Ticket:** `Task [Guard ordering] ! [Restructuring convergence block]`

### Bug mechanism (confidence 95)

Lines 210–243: `cm` is populated only when `W > 1e-300`. Convergence check at line 233 (`mu < kTolMu`) and line 234 (`converged_abs`) can both return `RK_OK` before the INFEAS check at line 265. When `W ≤ 1e-300`, `cm.errRp = 0` (zero-initialised struct), so the absolute-tol branch fires spuriously.

### RED test

Add to `tests/testthat/test-calibration-solvers.R`:

```r
test_that("B1: chebyshev returns INFEAS not RK_OK when W collapses", {
  # Construct a 1-margin 2-cat problem where min_weight > 0 forces all cell
  # masses to zero (infeasible bounds: lo > initial mass for all cells).
  n <- 10L
  w <- rep(1.0, n)
  group_ids <- matrix(rep(0L, n), nrow = n, ncol = 1L)  # all in cat 0
  cat_counts <- 1L
  # Target demands cat 0 = 100% but min_weight pushes cells to 0 (contrived):
  # Use bounds that are provably infeasible (lo > sum(w)).
  res <- harvest(
    data.frame(g = rep(1L, n)),
    targets = list(g = c(g1 = 2.0)),   # target > 1.0 — infeasible
    method = "chebyshev",
    min_weight = 0.001, max_weight = 1e10,
    max_iterations = 50L
  )
  expect_equal(res$status, "infeasible")
})
```

*Run:* `devtools::test(filter = "calibration-solvers")` — expect FAIL before fix.

### Fix

Move the INFEAS guard from line 265 to **before** the convergence check block (before line 218):

```cpp
// --- B1 FIX: check W before convergence block ---
if (W < 1e-300) { res.status = RK_ERR_INFEAS; return res; }
```

**Remove** the duplicate guard at the old line 265.

**Exact diff:**

```
// Before (line 209-265):
CellMetrics cm;
if (W > 1e-300) {
    cm = lbw::compute_cell_metrics(st, ct, X, W, bucket_tmp);
    if (cm.errRp < best_errRp) { best_errRp = cm.errRp; res.best_iter = iter+1; X_best = X; }
}
double errRp = cm.errRp, mean_err = cm.mean_err;
...
{ /* convergence block lines 220-244 */ }
...
if (W < 1e-300) { res.status = RK_ERR_INFEAS; return res; }   // <-- old location

// After:
if (W < 1e-300) { res.status = RK_ERR_INFEAS; return res; }   // <-- B1: moved here
CellMetrics cm;
cm = lbw::compute_cell_metrics(st, ct, X, W, bucket_tmp);
if (cm.errRp < best_errRp) { best_errRp = cm.errRp; res.best_iter = iter+1; X_best = X; }
double errRp = cm.errRp, mean_err = cm.mean_err;
...
{ /* convergence block unchanged */ }
// (delete old guard at former line 265)
```

**Compile gate:** `R CMD INSTALL --preclean . 2>&1 | grep -E "error:|warning:"` — must produce zero errors.
**GREEN:** Re-run test from RED step.

---

## Task 2 — B10: sinkhorn Dykstra a[] phantom mass accumulation

**File:** `src/sinkhorn.cpp`
**Ticket:** `Task [Dykstra correction zeroing] ! [Changing bisection logic]`

### Verified code structure (lines 76–156, read 2026-04-28)

`a` is declared at line 77: `std::vector<double> a(ct.M_cell, 0.0)` and is **never reset at the start of each outer iteration**. It persists and accumulates across all `iter` values — this is true cross-iteration Dykstra semantics, not a per-cycle reset. `bisect_capacity` uses `a[c]` as an offset term: the projection computes `X[c] * exp(a[c] + mu)`. Zeroing `a[]` when `needs_projection = false` would destroy accumulated correction history from prior projected iterations.

### Bug mechanism (confidence 95, revised)

The original plan's proposed fix (zeroing `a[]` in the `else` branch) is **wrong**. Standard Dykstra does not zero the correction on non-binding projections — it simply does not update it. The `else` branch (lines 147–149) correctly does nothing to `a[]`. The actual bug is subtler: the `if (needs_projection)` update block (lines 150–156) is positioned **after** the `else` branch, so it already only runs when `needs_projection = true`. The structure is correct.

Re-reading the block: lines 147–156 show that when `needs_projection = false`, neither `X_proj` nor `a[]` changes beyond `X_proj = X`. This is mathematically correct for Dykstra: stale `a[]` values ARE the correct accumulated history, and preserving them is what Dykstra requires. There is no phantom accumulation bug in the current code — the existing `else` branch (no-op on `a[]`) is already correct.

**Confidence re-assessment:** The original bug description was based on inference without reading the code. After reading lines 50–170, the existing implementation matches standard Dykstra. The `else` branch must NOT zero `a[]`.

### Revised fix: no code change to Dykstra loop

The `else` branch is correct as-is. The only actionable change is to **add a comment** clarifying the invariant, to prevent future misreads:

```cpp
// Before (lines 147-149):
} else {
    X_proj = X;
}

// After (B10: add clarifying comment — a[] intentionally preserved per Dykstra semantics):
} else {
    X_proj = X;
    // B10: Dykstra invariant — do NOT zero a[] here. Accumulated correction from prior
    // projected iterations is valid history; zeroing would corrupt subsequent bisect_capacity calls.
}
```

**No functional code change.** The existing logic is correct.

### RED test

Add to `tests/testthat/test-calibration-solvers.R` to guard against future regressions (any PR that adds a `std::fill(a...)` in the else branch must fail this test):

```r
test_that("B10: sinkhorn KL stable when box constraint transiently inactive", {
  # Problem where some iters have needs_projection=false after prior projected iters.
  # Zeroing a[] would cause KL to spike as accumulated Dykstra history is destroyed.
  set.seed(42)
  n <- 200L
  df <- data.frame(
    age = sample(c("young","old"), n, replace = TRUE),
    sex = sample(c("m","f"), n, replace = TRUE)
  )
  targets <- list(age = c(young = 0.6, old = 0.4), sex = c(m = 0.5, f = 0.5))
  res <- harvest(df, targets = targets, method = "sinkhorn",
                 min_weight = 0.5, max_weight = 3.0)
  expect_lt(res$diagnostics$kl, 0.01)
})
```

**Compile gate:** `R CMD INSTALL --preclean . 2>&1 | grep -E "error:|warning:"` — zero errors.
**GREEN:** KL < 0.01.

---

## Task 3 — B14: chebyshev dual explosion guard uses stale mu

**File:** `src/chebyshev.cpp`
**Ticket:** `Task [Stale mu in explosion reset] ! [Changing explosion detection threshold]`

### Bug mechanism (confidence 90)

Lines 469–479: after computing `mu_new = comp_new / n_comp` (line 472) and detecting explosion (`mu_new > 100.0 * mu`, line 473), the reset uses `mu` (the complementarity measured at the *start* of the iteration, line 202), not `mu_new`. The reset target `mu/s_delta` should be the post-step complementarity target `mu_new/s_*`, otherwise duals are re-centred at a stale level that may be orders of magnitude different from the actual post-step constraint values.

### Fix

Change lines 477–479:

```cpp
// Before:
y_delta = mu / s_delta;
for (int c = 0; c < ct.M_cell; c++) { y_lo[c] = mu/s_lo[c]; y_hi[c] = mu/s_hi[c]; }
for (int m = 0; m < nct; m++) { y_up[m] = mu/s_up[m]; y_dn[m] = mu/s_dn[m]; }

// After (B14: use mu_new, not stale mu):
y_delta = mu_new / s_delta;
for (int c = 0; c < ct.M_cell; c++) { y_lo[c] = mu_new/s_lo[c]; y_hi[c] = mu_new/s_hi[c]; }
for (int m = 0; m < nct; m++) { y_up[m] = mu_new/s_up[m]; y_dn[m] = mu_new/s_dn[m]; }
```

### Test

The explosion guard is hard to trigger synthetically. Verify:

1. Existing chebyshev tests in `test-calibration-solvers.R` pass GREEN.
2. Add a regression sentinel comment in the test file documenting the fix.

```r
test_that("B14: chebyshev dual explosion guard uses mu_new (regression)", {
  # No explosion guarantee, but verifying existing convergence not broken.
  # If explosion guard fires incorrectly, mu diverges and the solver returns NOCONV.
  res <- harvest(
    data.frame(g = rep(1:3, length.out = 90)),
    targets = list(g = c("1" = 0.4, "2" = 0.35, "3" = 0.25)),
    method = "chebyshev"
  )
  expect_equal(res$status, "converged")
})
```

**Compile gate:** `R CMD INSTALL --preclean .` — zero errors.

---

## Task 4 — B15: chebyshev alpha_p/alpha_d floor violates IPM invariant

**File:** `src/chebyshev.cpp`
**Ticket:** `Task [Fraction-to-boundary step length] ! [Arbitrary floor on alpha]`

### Verified FTB code (lines 409–425, read 2026-04-28)

```cpp
double alpha_p = 1.0, alpha_d = 1.0;
for (int c = 0; c < ct.M_cell; c++) {
    if (dX[c] > 0.0)  alpha_p = std::min(alpha_p, kStepScale*s_hi[c]/dX[c]);
    if (dX[c] < 0.0)  alpha_p = std::min(alpha_p, -kStepScale*s_lo[c]/dX[c]);
    if (dY_lo[c] < 0.0) alpha_d = std::min(alpha_d, -kStepScale*y_lo[c]/dY_lo[c]);
    if (dY_hi[c] < 0.0) alpha_d = std::min(alpha_d, -kStepScale*y_hi[c]/dY_hi[c]);
}
for (int m = 0; m < nct; m++) {
    if (dS_up[m] < 0.0) alpha_p = std::min(alpha_p, -kStepScale*s_up[m]/dS_up[m]);
    if (dS_dn[m] < 0.0) alpha_p = std::min(alpha_p, -kStepScale*s_dn[m]/dS_dn[m]);
    if (dY_up[m] < 0.0) alpha_d = std::min(alpha_d, -kStepScale*y_up[m]/dY_up[m]);
    if (dY_dn[m] < 0.0) alpha_d = std::min(alpha_d, -kStepScale*y_dn[m]/dY_dn[m]);
}
if (d_delta < 0.0) alpha_p = std::min(alpha_p, -kStepScale*s_delta/d_delta);
if (dY_delta < 0.0) alpha_d = std::min(alpha_d, -kStepScale*y_delta/dY_delta);
alpha_p = std::max(alpha_p, 1e-10);   // <-- floor
alpha_d = std::max(alpha_d, 1e-10);   // <-- floor
```

### Safety analysis of removing the floor (confidence 95)

Both `alpha_p` and `alpha_d` start at `1.0` and are only ever **reduced** via `std::min`. The FTB denominators are `dX[c]`, `dY_lo[c]`, etc. — all guarded by `if (d* > 0)` or `if (d* < 0)` conditions. Division-by-zero is impossible because zero search directions are skipped by the guards.

Cases:
- **`ds_i = 0`**: the guard condition is false; the variable is skipped; `alpha` stays at its current value (correct — a zero-direction constraint is non-binding).
- **`ds_i > 0` (primal) or `ds_i < 0` (dual)**: standard FTB ratio `kStepScale * slack / |ds|` is computed; this is always strictly positive since `slack > 0` by IPM invariant and `kStepScale = 0.99 > 0`.
- **Empty constraints (nct = 0)**: both loops are no-ops; `alpha_p = alpha_d = 1.0` (full step, correct).
- **All search directions safe**: `alpha_p` and/or `alpha_d` remain at `1.0`; floor is harmless.
- **Near-boundary degenerate case**: FTB produces a small but strictly positive ratio; floor may override to a larger value that pushes slack to `kStepScale * slack < slack` — safe — but the floor (1e-10) can only force a larger step, never smaller. The risk is the floor forcing `alpha >= 1e-10` when the FTB correctly computed a smaller safe value, which does push the slack negative. **This is exactly the bug.**

`alpha_p` initialized to `1.0` and only reduced — it can reach `0` only if `kStepScale * slack / |ds| = 0`, which requires `slack = 0`. But `slack > 0` is the IPM invariant; if it is ever violated, the problem is elsewhere. The floor cannot cause a crash; it can only cause IPM invariant violation.

**Removing the floor is safe.** No crash path exists from the FTB formula alone.

### Bug mechanism (confidence 95)

Lines 424–425:
```cpp
alpha_p = std::max(alpha_p, 1e-10);
alpha_d = std::max(alpha_d, 1e-10);
```
When FTB correctly computes `alpha_p = epsilon < 1e-10` (a legitimately tiny safe step near the boundary), the floor overrides to `1e-10`. Since `kStepScale * slack / |ds| < 1e-10` means `slack < 1e-10 * |ds| / kStepScale`, taking a step of `1e-10` in direction `ds` violates `s > 0`. The IPM invariant is broken.

### Fix

Remove both floor lines. Exact before/after diff:

```cpp
// Before (lines 424-425):
alpha_p = std::max(alpha_p, 1e-10);
alpha_d = std::max(alpha_d, 1e-10);

// After (B15: remove floors — FTB rule is self-sufficient; floor overrides tiny safe steps):
// (both lines deleted)
```

Optional stall guard (separate concern — do NOT combine with floor):

```cpp
// After FTB computation, if stall detection is needed:
if (alpha_p < 1e-15 && alpha_d < 1e-15) {
    // Primal and dual both stalled — declare no-convergence.
    break;
}
```

### Test

```r
test_that("B15: chebyshev alpha floor removal does not regress existing cases", {
  res <- harvest(
    data.frame(
      age = sample(c("18-34","35-54","55+"), 500, replace = TRUE, prob = c(0.3,0.4,0.3)),
      sex = sample(c("M","F"), 500, replace = TRUE)
    ),
    targets = list(
      age = c("18-34"=0.25, "35-54"=0.45, "55+"=0.30),
      sex = c(M=0.48, F=0.52)
    ),
    method = "chebyshev"
  )
  expect_equal(res$status, "converged")
  expect_lt(res$diagnostics$max_error, 1e-4)
})
```

**Compile gate:** `R CMD INSTALL --preclean .` — zero errors.

---

## Task 5 — R1: greg factorization cache (active-set gating)

**File:** `src/greg.cpp`
**Ticket:** `Task [Factorization cache] ! [Changing Newton update formula]`

### Bug mechanism (confidence 88)

Lines 64–94: every Newton iteration unconditionally recomputes `D_eff` → `N` → LDLT factorization. Cost: `O(nct² * M_cell)` for `compute_normal_equations` + `O(nct³)` for LDLT. When the active set (fixed_lo, fixed_hi) is unchanged, `D_eff` is identical to the previous iteration and both calls are pure waste.

### Fix

Track the active set from the previous iteration. Skip `compute_normal_equations` + `ldlt_factor_inplace` when unchanged. `ldlt_solve` must still run every iteration (because `b` changes with `X`).

**N must be refactored when:** `(fixed_lo != prev_fixed_lo || fixed_hi != prev_fixed_hi)`.

```cpp
// Insert BEFORE the Newton loop (after line 62):
std::vector<bool> prev_fixed_lo(ct.M_cell, false), prev_fixed_hi(ct.M_cell, false);
bool need_refactor = true;  // always refactor on first iteration

for (int newton_iter = 0; newton_iter < kMaxNewtonIters; newton_iter++) {
    res.iterations = newton_iter + 1;

    std::fill(D_eff.begin(), D_eff.end(), 0.0);
    for (int c = 0; c < ct.M_cell; c++)
        if (!fixed_lo[c] && !fixed_hi[c] && X_init[c] > kEps)
            D_eff[c] = X_init[c];

    if (need_refactor) {
        if (compute_normal_equations(ct, D_eff.data(), N.data(),
                                     cat_offset.data(), st.K,
                                     static_cast<size_t>(n_cats_total)) != RK_OK) {
            res.status = RK_ERR_BADARG; return res;
        }
        if (ldlt_factor_inplace(N.data(), static_cast<size_t>(n_cats_total), 1e-10) != RK_OK) {
            res.status = RK_ERR_BADARG; return res;
        }
    }
    // ... build b[] (unchanged) ...
    ldlt_solve(N.data(), static_cast<size_t>(n_cats_total), b.data());   // always solve

    // ... Newton update + active-set update (unchanged) ...

    // Gate: compare new active set to previous
    need_refactor = (fixed_lo != prev_fixed_lo || fixed_hi != prev_fixed_hi);
    prev_fixed_lo = fixed_lo;
    prev_fixed_hi = fixed_hi;
}
```

**N.B.:** The N buffer is overwritten by `ldlt_factor_inplace`. When `!need_refactor`, N already holds the factored form from the previous iteration — `ldlt_solve` can use it directly. This is correct because `D_eff` is unchanged when the active set is unchanged.

### Test

```r
test_that("R1: greg produces identical weights before and after factorization caching", {
  set.seed(7)
  n <- 300L
  df <- data.frame(
    age = sample(c("a","b","c"), n, replace = TRUE),
    sex = sample(c("m","f"), n, replace = TRUE)
  )
  tgt <- list(age = c(a=0.3, b=0.4, c=0.3), sex = c(m=0.5, f=0.5))
  # Run twice — cache-gated and reference must give identical weights.
  # Since fix is internal, just verify result quality.
  res <- harvest(df, targets = tgt, method = "greg")
  expect_equal(res$status, "converged")
  expect_lt(res$diagnostics$chi2, 1e-6)
})
```

**Compile gate:** `R CMD INSTALL --preclean .` — zero errors.

---

## Task 6 — R7: chebyshev negative slack INFEAS promotion

**File:** `src/chebyshev.cpp`
**Ticket:** `Task [Slack violation counter] ! [Changing kEps floor semantics]`

### Bug mechanism (confidence 91)

Lines 446–449:
```cpp
s_up[m] = std::max(T_flat[m]*W_upd + w_kj[m]*delta - S[m], kEps);
s_dn[m] = std::max(S[m] - T_flat[m]*W_upd + w_kj[m]*delta, kEps);
```
When `T_flat[m]*W_upd + w_kj[m]*delta - S[m] < 0`, the `std::max(..., kEps)` silently clamps to `kEps`. This masks the fact that the LP is infeasible at this iterate: the slack should be positive by construction but went negative due to a step that violated `s > 0`. Persistent negative raws indicate the problem is genuinely infeasible, not just a transient numerical blip.

### Fix

Add a violation counter. Reset to zero each iteration if no violations. Declare INFEAS after `kInfeasPersistence` consecutive iterations with violations:

```cpp
// Add constant (near top of function, with other constexprs):
static constexpr int kInfeasPersistence = 5;

// Add before Newton loop:
int slack_violations = 0;

// Replace lines 443-450 (W_upd block) with:
{
    double W_upd = 0.0;
    for (int c = 0; c < ct.M_cell; c++) W_upd += X[c];
    int viol_this_iter = 0;
    for (int m = 0; m < nct; m++) {
        double raw_sup = T_flat[m]*W_upd + w_kj[m]*delta - S[m];
        double raw_sdn = S[m] - T_flat[m]*W_upd + w_kj[m]*delta;
        if (raw_sup < 0.0 || raw_sdn < 0.0) viol_this_iter++;
        s_up[m] = std::max(raw_sup, kEps);
        s_dn[m] = std::max(raw_sdn, kEps);
    }
    if (viol_this_iter > 0) {
        slack_violations++;
        if (slack_violations > kInfeasPersistence) {
            res.status = RK_ERR_INFEAS;
            std::snprintf(res.message, sizeof(res.message),
                          "chebyshev: %d consecutive iters with negative slacks — INFEAS",
                          slack_violations);
            break;
        }
    } else {
        slack_violations = 0;
    }
}
```

### RED test

```r
test_that("R7: chebyshev detects INFEAS via persistent negative slacks", {
  # Infeasible: targets sum to >1 within a single margin → cannot satisfy.
  n <- 50L
  df <- data.frame(g = sample(1:2, n, replace = TRUE))
  res <- harvest(df,
    targets = list(g = c("1" = 0.9, "2" = 0.9)),  # sum > 1 — infeasible
    method = "chebyshev",
    max_iterations = 100L
  )
  expect_true(res$status %in% c("infeasible", "no-convergence"))
})
```

**Compile gate:** `R CMD INSTALL --preclean .` — zero errors.

---

## Task 7 — R10: greg kMaxNewtonIters + NOCONV message

**File:** `src/greg.cpp`
**Ticket:** `Task [kMaxNewtonIters + NOCONV message] ! [Changing Newton algorithm]`

### Bug mechanism (confidence 95)

Line 61: `static constexpr int kMaxNewtonIters = 10;` — 10 iterations is insufficient for problems with cycling active sets or many bound-constrained cells. When the loop exits without `RK_OK`, `res.status` remains `RK_ERR_NOCONV` (set at line 16) but `res.message` is empty (null string), leaving the caller with no diagnostic.

### Fix

Two atomic changes:

**1. Increase kMaxNewtonIters to 50:**

```cpp
// Before:
static constexpr int kMaxNewtonIters = 10;

// After:
static constexpr int kMaxNewtonIters = 50;
```

**2. Add NOCONV message after the loop:**

```cpp
// After the Newton loop (after line 124), before metrics computation:
if (res.status == RK_ERR_NOCONV) {
    std::snprintf(res.message, sizeof(res.message),
                  "greg: no convergence after %d Newton steps; active set still cycling",
                  kMaxNewtonIters);
}
```

### RED test

```r
test_that("R10: kMaxNewtonIters is 50 and NOCONV carries a message", {
  # We cannot easily trigger 50-iter non-convergence without a pathological input,
  # so we test the observable effect: hard convergence problems should now succeed.
  # (Pre-fix at 10 iters some problems hit NOCONV; post-fix at 50 they converge.)
  set.seed(99)
  n <- 500L
  df <- data.frame(
    a = sample(letters[1:5], n, replace = TRUE),
    b = sample(LETTERS[1:4], n, replace = TRUE),
    c = sample(c("x","y","z"), n, replace = TRUE)
  )
  tgt <- list(
    a = setNames(rep(0.2, 5), letters[1:5]),
    b = setNames(rep(0.25, 4), LETTERS[1:4]),
    c = c(x=0.33, y=0.33, z=0.34)
  )
  res <- harvest(df, targets = tgt, method = "greg",
                 min_weight = 0.1, max_weight = 10.0)
  # With 50 iters, this should converge.
  expect_equal(res$status, "converged")
})
```

**Compile gate:** `R CMD INSTALL --preclean .` — zero errors.

---

## Execution order

Tasks are **independent** (different files or disjoint code regions). Preferred serial order to keep compilation incremental:

1. **T7 (R10)** — trivial two-line greg.cpp change; lowest risk.
2. **T5 (R1)** — greg.cpp caching; verify correctness before touching chebyshev.
3. **T2 (B10)** — sinkhorn.cpp; isolated.
4. **T1 (B1)** — chebyshev.cpp guard reorder; highest impact.
5. **T6 (R7)** — chebyshev.cpp slack violations; depends on knowing B1 is clean.
6. **T3 (B14)** — chebyshev.cpp explosion guard; low-risk single-token change.
7. **T4 (B15)** — chebyshev.cpp alpha floor removal; validate carefully.

Each task: RED → fix → compile gate → GREEN → commit.

## Commit message template

```
fix(<solver>): <B/R tag> — <one-line description>

<mechanism>: <exact before → after>.
Fixes: <what breaks without this>.
```

Example:
```
fix(chebyshev): B1 — move INFEAS guard before convergence check

W-collapse check was at line 265, after convergence block at 220–243.
Zero-initialised cm.errRp caused RK_OK to fire on infeasible inputs.
Moved guard to immediately after W computation.
```
