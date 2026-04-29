# i0am: Safeguarded Regularized Anderson Acceleration (SRAA-m)

**Date**: 2026-04-29
**Status**: Approved for implementation
**Ticket**: leafblower-i0am
**Files**: `src/sraa.hpp` (new), `src/greenkhorn.cpp`, `src/raking.cpp`,
           `tests/testthat/test-calibration-solvers.R`

---

## Problem

`greenkhorn+squarem` converges to a **worse** fixed point than plain greenkhorn:

| Method | max_err | iters | vs plain |
|---|---|---|---|
| greenkhorn plain | 1.57e-3 | 1030 | baseline |
| greenkhorn+squarem | 2.12e-3 | 540 | **+35% worse** |
| raking plain | 1.60e-3 | 50 | baseline |
| raking+squarem | 1.75e-3 | 125 | +9% worse |

**Root cause**: The greenkhorn SQUAREM has no backtracking. Its acceptance criterion
`err_star ≤ err_w2 * 1.01` compares against the second F_eval result (not the starting
point) with 1% slack — consistently accepting steps that are slightly worse, accumulating
across 540 super-steps. Raking's step-halving backtracking partially mitigates this but
still overshoots.

The CBB step size α is calibrated for smooth operators. At bounds, the fixed-point
operator has a kink (non-smooth). The scalar α overestimates the useful step length,
extrapolating past the bounded optimum.

**Anderson Acceleration (AA-m)** replaces the scalar CBB extrapolation with an m-point
least-squares combination of past iterates. It naturally handles non-smooth bounded
operators and uses 2 F_evals per accepted step vs SQUAREM's 3.

---

## Design

### Architecture

**New file: `src/sraa.hpp`** — header-only template (~200 lines)

```
src/
  sraa.hpp          ← new: SRAAState + sraa_step<FEval>() template
  greenkhorn.cpp    ← remove CBB-SQUAREM block (lines ~159-210), add sraa_step call
  raking.cpp        ← remove step-halving SQUAREM block (lines ~360-515), add sraa_step call
tests/testthat/test-calibration-solvers.R  ← update T_acc
```

**`SRAAState`** — pre-allocated flat circular buffer, zero heap allocation per step:

```cpp
struct SRAAState {
    int M;            // M_cell
    int m;            // window size (default 5)
    int head;         // circular buffer write index (0..m-1)
    int count;        // valid entries (0..m)
    bool has_prev;

    // Column-major layout: buf[i * M + c] = slot i, cell c
    // Column-major → dot products read two contiguous M-length arrays → SIMD-vectorizable
    std::vector<double> dX_buf;     // m × M: ΔX_i = X_i - X_{i-1}
    std::vector<double> dR_buf;     // m × M: ΔR_i = R_i - R_{i-1}
    std::vector<double> R_prev;     // M: R_{k-1}
    std::vector<double> X_prev;     // M: X_{k-1}
    std::vector<double> F_cur;      // M: F(X_k) result
    std::vector<double> scratch;    // M: AA extrapolation staging

    // Stack-allocated m×m LS (lives in L1 cache, m ≤ 10)
    double gram[10 * 10];           // ΔR^T ΔR + δI
    double rhs[10];                 // ΔR^T R_k
    double gamma_[10];              // LS solution γ

    double prev_resid_norm;

    void init(int M_cell, int window);
    void clear() { head = 0; count = 0; has_prev = false; }
};
```

Memory at stepstone (m=5, M=1.58M):
- dX_buf + dR_buf: 2 × 5 × 12.6 MB = 126 MB
- R_prev + X_prev + F_cur + scratch: 4 × 12.6 MB = 50 MB
- gram/rhs/gamma (stack): 880 bytes
- **Total: ~176 MB** (vs SQUAREM's 9 × 12.6 MB = 113 MB scratch → net +63 MB)

**`sraa_step` signature:**

```cpp
struct SRAAStepResult { bool aa_accepted; int f_evals; double err_result; };

template<typename FEval>
SRAAStepResult sraa_step(
    FEval& f_eval,                        // (std::vector<double>&) -> double errRp
    std::vector<double>& X,               // current iterate, modified in-place
    const std::vector<double>& L_cell,
    const std::vector<double>& U_cell,
    SRAAState& state
);
```

The `f_eval` callable is self-contained: recomputes S and W from X before running K
steps. Required because SRAA calls F at arbitrary extrapolated X_AA values.

---

### Algorithm: SRAA-m Step (Type II Anderson Acceleration)

```
Inputs:  X_k, L_cell, U_cell, f_eval, SRAAState
Outputs: X updated in-place, SRAAStepResult

1. state.F_cur = f_eval(X_k)              [1 F_eval, costs K iters]
   err_plain   = return value of f_eval
   R_k[c]      = F_cur[c] - X_k[c]  ∀c

2. Append to circular buffer (if has_prev):
     slot = state.head
     dX_buf[slot * M + c] = X_k[c] - X_prev[c]  ∀c
     dR_buf[slot * M + c] = R_k[c]  - R_prev[c]  ∀c
     state.head  = (state.head + 1) % m
     state.count = min(state.count + 1, m)

3. Restart check:
     norm_k = ||R_k||²
     if norm_k > kSRAARestartGamma² * state.prev_resid_norm:
       state.clear()
       goto ACCEPT_PLAIN        [diverging — restart from fresh history]

4. state.X_prev = X_k
   state.R_prev = R_k
   state.prev_resid_norm = norm_k
   state.has_prev = true

5. if state.count < kSRAAMinCount: goto ACCEPT_PLAIN   [not enough history yet]

6. Build m×m Gram matrix (symmetric, sequential dot products over M):
     G[i][j] = dR[i] · dR[j] = Σ_c dR_buf[i*M+c] * dR_buf[j*M+c]
   Regularize: G[i][i] += kSRAAdeltaReg * max_j(G[j][j])
   RHS: rhs[i] = dR[i] · R_k

7. Solve (G + δI)γ = rhs  via ldlt_factor_inplace + ldlt_solve  [reuse calib_linalg.hpp]
   If ldlt_factor_inplace returns !RK_OK: state.clear(); goto ACCEPT_PLAIN

8. Extrapolate + clamp (vectorizable; m=5 unrolls over inner loop):
     scratch[c] = clamp(X_k[c] + R_k[c] - Σ_i γ[i]*(dX[i*M+c] + dR[i*M+c]),
                        L_cell[c], U_cell[c])  ∀c

9. err_AA = f_eval(state.scratch)           [1 F_eval, costs K iters]

10. Safeguard:
    if err_AA ≤ err_plain:
      X = state.scratch                     [accept AA — better than plain step]
      return {aa_accepted=true, f_evals=2, err_result=err_AA}
    else:
      state.clear()                         [history caused bad step; reset]
      goto ACCEPT_PLAIN

ACCEPT_PLAIN:
    X = state.F_cur                         [reuse already-computed F(X_k)]
    return {aa_accepted=false, f_evals=1, err_result=err_plain}
```

**Key invariants:**
- ACCEPT_PLAIN costs only 1 F_eval (F_cur already computed in step 1)
- AA path costs exactly 2 F_evals regardless of acceptance
- History cleared on: restart trigger, safeguard rejection, LDLT failure
- Safeguard compares `err_AA ≤ err_plain` — strictly better than plain F_eval (no slack)

---

### Parameters

| Constant | Value | Derivation |
|---|---|---|
| `kSRAAm` | 5 | 5 history pairs; 176 MB at stepstone; empirically optimal for survey calibration |
| `kSRAAMinCount` | 2 | Need ≥2 ΔX/ΔR pairs for LS to be non-trivial |
| `kSRAAdeltaReg` | 1e-10 | Relative Tikhonov; activates only when G is near-singular late in convergence |
| `kSRAARestartGamma` | 2.0 | Restart when `‖R_k‖² > 4 × prev_norm` (residual doubled = diverging) |

All constexpr, defined at top of `sraa.hpp`.

---

### Integration: greenkhorn.cpp

**Remove** (~9 scratch vectors + CBB block, lines ~139–210):
```cpp
// DELETE: sq_X_snap, sq_w1, sq_w2, sq_S1, sq_S2, sq_X_star, sq_Sstar, sq_W1, sq_W2, sq_Wstar
// DELETE: CBB alpha computation, extrapolation, accept/reject block
```

**Add** at solver entry (before iteration loop):
```cpp
SRAAState grk_sraa;
if (accelerate) grk_sraa.init(M, kSRAAm);
```

**Self-contained f_eval_sraa** (recomputes S/W from Xv):
```cpp
std::vector<double> Sv_sraa(K * S_stride, 0.0);  // scratch for SRAA F_eval
auto f_eval_sraa = [&](std::vector<double>& Xv) -> double {
    double Wv = 0.0;
    std::fill(Sv_sraa.begin(), Sv_sraa.end(), 0.0);
    for (int c = 0; c < M; c++) {
        Wv += Xv[c];
        for (int k = 0; k < K; k++) {
            int g = ct.g_per_cell[k][c];
            if (g >= 0 && g < st.cat_counts[k])
                Sv_sraa[k * S_stride + g] += Xv[c];
        }
    }
    std::swap(X, Xv); std::swap(S_flat, Sv_sraa); std::swap(W, Wv);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
        [&](int a, int b){ return errRp[a] > errRp[b]; });
    for (int ki = 0; ki < K; ki++) greenkhorn_step(order[ki]);
    std::swap(X, Xv); std::swap(S_flat, Sv_sraa); std::swap(W, Wv);
    return *std::max_element(errRp.begin(), errRp.end());
};
```

**Replace** the `if (accelerate)` block:
```cpp
if (accelerate) {
    auto r = sraa_step(f_eval_sraa, X, L_cell, U_cell, grk_sraa);
    res.iterations += K * (r.aa_accepted ? 2 : 1);
} else {
    int k_star = (int)(std::max_element(errRp.begin(), errRp.end()) - errRp.begin());
    greenkhorn_step(k_star);
    res.iterations += 1;
}
```

---

### Integration: raking.cpp

Same pattern: remove step-halving SQUAREM block (~lines 360–515), add `SRAAState rk_sraa`,
wrap raking's existing F_eval as `f_eval_sraa`.

**Preserve unchanged** (in outer loop, not in SQUAREM block):
- Infeasibility tracking (`is_infeasible`)
- Weight-change stall detection (`X_prev_sq`, wchange convergence check)
- `errRp_new` update: set from `sraa_step` return value `r.err_result`

---

### TDD Requirements

#### T_sraa_grk — AA fires and improves greenkhorn quality

```r
test_that("T_sraa_grk: greenkhorn+AA max_err <= plain and converges faster", {
  set.seed(99); n <- 2000L
  df  <- data.frame(x=factor(sample(letters[1:3],n,TRUE)),
                    y=factor(sample(c("M","F"),n,TRUE)))
  tgt <- list(x=c(a=0.3,b=0.4,c=0.3), y=c(M=0.5,F=0.5))
  r_aa    <- suppressWarnings(
    harvest(df, tgt, method="greenkhorn", accelerate=TRUE,  max_iterations=500L,
            attach_weights=FALSE))
  r_plain <- suppressWarnings(
    harvest(df, tgt, method="greenkhorn", accelerate=FALSE, max_iterations=500L,
            attach_weights=FALSE))
  me_aa    <- attr(r_aa,    "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  iters_aa    <- attr(r_aa,    "result")$iterations
  iters_plain <- attr(r_plain, "result")$iterations
  K_exp <- 2L

  # (1) AA quality must match or beat plain
  expect_lte(me_aa, me_plain * 1.001,
    label=sprintf("AA (%.2e) must not be worse than plain (%.2e)", me_aa, me_plain))

  # (2) AA fired: iters divisible by K*2 (all steps accepted on this well-conditioned problem)
  expect_equal(iters_aa %% (K_exp * 2L), 0L,
    label=sprintf("AA iters (%d) divisible by K*2=%d proves AA fired",
                  iters_aa, K_exp * 2L))

  # (3) AA faster than plain
  expect_lt(iters_aa, iters_plain,
    label=sprintf("AA (%d iters) must be faster than plain (%d iters)", iters_aa, iters_plain))
})
```

#### T_sraa_rk — raking+AA quality matches or beats raking+squarem

```r
test_that("T_sraa_rk: raking+AA max_err <= raking plain", {
  set.seed(99); n <- 2000L
  df  <- data.frame(x=factor(sample(letters[1:3],n,TRUE)),
                    y=factor(sample(c("M","F"),n,TRUE)))
  tgt <- list(x=c(a=0.3,b=0.4,c=0.3), y=c(M=0.5,F=0.5))
  r_aa    <- suppressWarnings(
    harvest(df, tgt, method="raking", accelerate=TRUE,  max_iterations=500L,
            attach_weights=FALSE))
  r_plain <- suppressWarnings(
    harvest(df, tgt, method="raking", accelerate=FALSE, max_iterations=500L,
            attach_weights=FALSE))
  me_aa    <- attr(r_aa,    "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  expect_lte(me_aa, me_plain * 1.001,
    label=sprintf("raking+AA (%.2e) must not be worse than plain (%.2e)", me_aa, me_plain))
})
```

#### T_acc update (replaces K*3 divisibility check)

Replace existing T_acc body with `T_sraa_grk` assertions above. The K*3 divisibility
check is SQUAREM-specific and no longer valid after AA replaces CBB.

#### Regression: T5–T8, T_logit_armijo, T_logit_init must remain GREEN

---

## Acceptance Criteria

| # | Criterion | Verify |
|---|-----------|--------|
| AC1 | greenkhorn+AA max_err ≤ greenkhorn plain | `devtools::test(filter="calibration-solvers")` |
| AC2 | raking+AA max_err ≤ raking plain | Same |
| AC3 | T_sraa_grk GREEN (AA fires, faster) | Same |
| AC4 | T_sraa_rk GREEN | Same |
| AC5 | FAIL count unchanged (3) | `devtools::test()` |
| AC6 (benchmark) | stepstone greenkhorn+AA max_err ≤ 1.57e-3 | `Rscript benchmarks/stepstone_all_methods.R` |
| AC7 (benchmark) | stepstone raking+AA max_err ≤ 1.60e-3 | Same |

---

## Out of Scope

- Exposing `m`, `delta_reg`, `restart_gamma` to R user (internal; hardcoded is simpler)
- AA for `method="sinkhorn"` (different algorithm structure)
- Adaptive m (start small, grow) — follow-on if m=5 quality warrants
- Per-obs weight-level AA (cell-aggregate is sufficient and cheaper)
- AA for logit (uses Newton, not fixed-point iteration)
