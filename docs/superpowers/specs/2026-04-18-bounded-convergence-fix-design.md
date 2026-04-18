# Bounded Convergence Fix — Design Spec

## Problem

Two bugs prevent leafblower from converging on problems where `max_weight` binds:

**Bug 1 — iEPPA:** `bcd_sweep` clamps weights during each IPF scaling step, violating
the Sinkhorn invariant. The next margin's bucket accumulation starts from a biased point,
causing cycling near the constraint boundary. The outer EPPA loop exhausts all iterations
without converging. The post-convergence projection loop (lines 138–153 of `ieppa.cpp`)
compounds the damage: it modifies weights after the convergence criterion fires, so
`res.max_error` reflects the post-projection degradation, not the solver's actual output.
Observed: `max_error = 2.3e-3` after 2000 outer iterations on a 3-margin, n=10000 problem
with `max_weight=5`.

**Bug 2 — L-BFGS-B:** `LinkFn` uses the exp link whenever `L == 0.0` (the default),
even when `max_weight` is finite. Exp link maps the dual to `(0, ∞)`, so `max_weight`
is applied as a post-hoc clamp that breaks dual-primal correspondence. The gradient at
the clamped point is not zero; the solver can never declare convergence.
Observed: `max_error = 1.7e-2` after 2000 iterations on the same problem.

Benchmark comparison with `autumn 0.2.0` (unconstrained): autumn converges to `1.9e-12`;
leafblower fails to converge on the constrained problem.

---

## Architecture

Three files change. No new files. No API changes.

| File | Change |
|------|--------|
| `src/ieppa.cpp` | Replace `bcd_sweep` + EPPA outer loop + post-projection with `dykstra_solve` |
| `src/logit.hpp` | Change link selection: `exponential = !std::isfinite(U)` |
| `src/c_api.cpp` | Update `use_logit` guard to match new `LinkFn` condition |
| `tests/testthat/test-bounded-convergence.R` | New file: RED tests written first (TDD) |

---

## Fix 1: iEPPA — Dykstra's Alternating Projections

### Why Dykstra's, not "BCD continuation"

BCD-with-clamp (current approach) does not converge to the true constrained optimum —
it cycles near the constraint boundary. Dykstra's alternating projections maintains
per-constraint correction vectors that provably force convergence to the intersection
of all constraint sets. It is also faster: no wasted EPPA outer iterations, no cleanup
phase, monotone convergence from iteration 1.

### Algorithm

The constrained raking problem is: find `w` in the intersection of K marginal constraint
sets `{C_k}` and one box constraint `Box = [L, U]^n`.

```
// Correction vectors (Dykstra)
p[k][i] = 0.0   for k=0..K-1, i=0..n-1    // marginal corrections
q[i]    = 0.0   for i=0..n-1               // box correction

for iter = 1..max_iter:
  // Marginal projections
  for k = 0..K-1:
    y[i]    = w[i] + p[k][i]               // apply correction
    W       = sum(y)
    bucket[j] = sum_{g_k(i)==j} y[i]
    if bucket[j] < 1e-15*W and tau[k][j]>0: infeas_flag = true
    scale[j]  = tau[k][j] * W / max(bucket[j], 1e-15*W)
    w_new[i]  = y[i] * scale[g_k(i)]       // IPF step — NO CLAMP
    p[k][i]   = y[i] - w_new[i]            // update correction
    w = w_new

  // Box projection
  y[i]        = w[i] + q[i]                // apply correction
  w_box[i]    = clamp(y[i], L, U)
  q[i]        = y[i] - w_box[i]            // update correction
  w = w_box

  errRp = compute_errRp(w)
  if errRp < tol_abs: status = RK_OK; break

if infeas_flag and status == RK_OK: status = RK_ERR_INFEAS
```

### What is deleted

- `bcd_sweep` function (entire)
- `bregman_dist` function (entire)
- The outer EPPA loop (`for outer = 1..outer_max_iter`)
- The post-convergence projection loop (lines 138–153)

`compute_errRp` is retained unchanged.

### Parameter mapping

`inner_max_iter` becomes the single iteration budget (renamed conceptually to
`max_iter` internally). `outer_max_iter` is no longer used by iEPPA internally;
the struct field is preserved for API compatibility and ignored.

For backward compatibility, the effective budget is `inner_max_iter` iterations
(default 500). Users who previously relied on `outer_max_iter * inner_max_iter`
total work should increase `max_iterations` in the R call.

### Memory

K+1 extra vectors of size n. At n=10000, K=3: 320 KB. At n=1M, K=20: 168 MB.
Acceptable — comparable to the L-BFGS-B s/y deque.

---

## Fix 2: L-BFGS-B — Logit Link for Finite max_weight

### `src/logit.hpp` line 19

```cpp
// Before (broken):
exponential = (L == 0.0 || !std::isfinite(U));

// After (fixed):
exponential = !std::isfinite(U);
```

With `L=0, U=5`:
- `logit_scale = (5-0)/((5-1)*(1-0)) = 1.25`
- `F(u) = 5·exp(1.25u) / (4 + exp(1.25u))` maps exactly to `(0, 5)`
- The post-hoc clamp in `compute_final_weights_and_error` becomes a no-op

### `src/c_api.cpp` — singularity guard

```cpp
// Before:
bool use_logit = (p->min_weight > 0.0) && std::isfinite(p->max_weight);

// After:
bool use_logit = std::isfinite(p->max_weight);
```

Existing singularity checks (`min_weight == 1.0`, `max_weight == 1.0`) remain
correct. The `L=0, U=1.0` edge case is already caught by the `max_weight == 1.0`
guard (produces `(U-1)=0` in `logit_scale`).

### No routing changes

Auto-routing thresholds (`complexity > 500000`, `max_weight < 3.0`, `min_weight > 0`)
are unchanged. Both fixes are algorithm-internal.

---

## Tests (TDD)

New file: `tests/testthat/test-bounded-convergence.R`

Written RED (before C++ changes). Two tests:

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

**Regression guard:** All existing tests in `test-ieppa.R`, `test-harvest.R`,
`test-lbfgsb.R` (if present) must remain GREEN after both fixes.

---

## Success Criteria

1. Both new tests GREEN at `1e-6` tolerance
2. All existing 32 tests remain GREEN
3. `R CMD INSTALL --preclean .` succeeds (compilation gate)
4. No convergence warning emitted on the bound-hitting benchmark that previously failed
