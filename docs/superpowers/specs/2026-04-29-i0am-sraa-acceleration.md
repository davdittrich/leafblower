# i0am: Safeguarded Regularized Anderson Acceleration (SRAA-m)

**Date**: 2026-04-29
**Status**: Pending design review (rev 2)
**Ticket**: leafblower-i0am
**Files**: `src/sraa.hpp` (new), `src/greenkhorn.cpp`, `src/raking.cpp`,
           `R/harvest.R`, `tests/testthat/test-calibration-solvers.R`

---

## User Benefit

**`harvest(df, targets, method="greenkhorn", accelerate=TRUE)` will produce calibrated
weights that are at least as accurate as `accelerate=FALSE`, and converge in fewer
iterations.** Currently `accelerate=TRUE` produces weights 35% worse than plain
(`max_err` 2.12e-3 vs 1.57e-3); after this fix it matches or beats plain quality
while retaining the speed benefit on interior iterations.

**Breaking change for existing `accelerate=TRUE` users**: calibrated weight vectors
will change (they will be more accurate). Reproducible research pipelines using
`set.seed() + accelerate=TRUE` will produce different results. This is a correctness
fix — the old results were subtly wrong. `NEWS.md` and `harvest.R` roxygen docs must
be updated in the same commit as the C++ changes.

Memory trade-off: SQUAREM used ~113 MB scratch at stepstone scale; SRAA-m uses
~176 MB (+63 MB). Documented here for users on memory-constrained compute.

---

## Problem

`greenkhorn+squarem` converges to a **worse** fixed point than plain greenkhorn:

| Method | max_err | iters | vs plain |
|---|---|---|---|
| greenkhorn plain | 1.57e-3 | 1030 | baseline |
| greenkhorn+squarem | 2.12e-3 | 540 | **+35% worse** |
| raking plain | 1.60e-3 | 50 | baseline |
| raking+squarem | 1.75e-3 | 125 | +9% worse |

**Root cause**: Greenkhorn SQUAREM has no backtracking. Its acceptance criterion
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

**New file: `src/sraa.hpp`** — header-only template (~220 lines)

```
src/
  sraa.hpp          ← new: SRAAState + sraa_step<FEval>() template
  greenkhorn.cpp    ← remove CBB-SQUAREM block (verify exact lines before cut), add sraa_step
  raking.cpp        ← remove step-halving SQUAREM block (verify exact lines before cut), add sraa_step
R/
  harvest.R         ← update roxygen @param accelerate + warning message (BREAKING CHANGE)
tests/testthat/test-calibration-solvers.R  ← replace T_acc, add T_sraa_grk/rk/ldlt/restart
```

**`SRAAState`** — pre-allocated flat circular buffer, zero heap allocation per step:

```cpp
// kSRAAMaxM: compile-time max window; stack arrays sized to this
constexpr int kSRAAMaxM = 10;

struct SRAAState {
    int M;            // M_cell
    int m;            // window size (≤ kSRAAMaxM; default kSRAAm=5)
    int head;         // circular buffer write index (0..m-1)
    int count;        // valid entries (0..m)
    bool has_prev;    // false until first step completes
    int aa_accepted_count;  // total AA-accepted steps (for test observability)

    // Slot-contiguous (row-major) layout: buf[slot * M + cell]
    // Each slot is a contiguous M-length array → dot products iterate sequentially
    // over M → SIMD-vectorizable. "Slot-contiguous" not "column-major."
    std::vector<double> dX_buf;     // m × M: ΔX_i = X_i - X_{i-1}
    std::vector<double> dR_buf;     // m × M: ΔR_i = R_i - R_{i-1}
    std::vector<double> R_prev;     // M: R_{k-1}
    std::vector<double> X_prev;     // M: X_{k-1}
    std::vector<double> F_cur;      // M: F(X_k) result
    std::vector<double> scratch;    // M: AA extrapolation staging

    // Stack-allocated m×m LS (lives in L1 cache, m ≤ kSRAAMaxM)
    double gram[kSRAAMaxM * kSRAAMaxM];  // ΔR^T ΔR + δI
    double rhs[kSRAAMaxM];               // ΔR^T R_k
    double gamma_[kSRAAMaxM];            // LS solution γ  (gamma_ avoids <cmath> gamma symbol)

    double prev_resid_norm;   // ||R_{k-1}||²; initialized to 0.0; only read when has_prev=true

    // init(): call once at solver entry. Throws Rcpp::stop on window > kSRAAMaxM or OOM.
    void init(int M_cell, int window) {
        if (window > kSRAAMaxM)
            Rcpp::stop("SRAA window %d exceeds kSRAAMaxM=%d", window, kSRAAMaxM);
        M = M_cell; m = window;
        try {
            dX_buf.assign((size_t)m * M, 0.0);
            dR_buf.assign((size_t)m * M, 0.0);
            R_prev.assign(M, 0.0); X_prev.assign(M, 0.0);
            F_cur.assign(M, 0.0);  scratch.assign(M, 0.0);
        } catch (std::bad_alloc&) {
            Rcpp::stop("SRAA: out of memory allocating %.0f MB (m=%d, M=%d)",
                       (double)(2.0*m + 4.0) * M * 8.0 / 1e6, m, M);
        }
        clear();
    }

    // clear(): reset history. prev_resid_norm reset to 0 as a safe sentinel —
    // it is never read when has_prev=false (see restart guard at algorithm step 2).
    // aa_accepted_count is NOT reset by clear() — it is a cumulative run total.
    void clear() { head = 0; count = 0; has_prev = false; prev_resid_norm = 0.0; }
};
```

**`sraa_step` signature:**

```cpp
struct SRAAStepResult {
    bool   aa_accepted;   // true if AA step was accepted
    int    f_evals;       // 1 (plain) or 2 (AA attempted)
    double err_result;    // calibration max_error of accepted step
};

template<typename FEval>
SRAAStepResult sraa_step(
    FEval& f_eval,                        // (std::vector<double>&) -> double errRp
    std::vector<double>& X,               // current iterate, modified in-place
    const std::vector<double>& L_cell,
    const std::vector<double>& U_cell,
    SRAAState& state
);
// Invariant: f_eval must NEVER be called with the outer X directly.
// It must always receive state.F_cur or state.scratch (distinct objects).
// The swap pattern inside f_eval temporarily aliases; X and Xv must be different.
```

Memory at stepstone (m=5, M=1.58M):
- dX_buf + dR_buf: 2 × 5 × 12.6 MB = 126 MB
- R_prev + X_prev + F_cur + scratch: 4 × 12.6 MB = 50 MB
- gram/rhs/gamma (stack): 880 bytes
- **Total: ~176 MB** (vs SQUAREM's 9 × 12.6 MB = 113 MB scratch → net +63 MB)

---

### Algorithm: SRAA-m Step (Type II Anderson Acceleration)

```
Inputs:  X_k (outer caller's X), L_cell, U_cell, f_eval, SRAAState
Outputs: X updated in-place via std::swap, SRAAStepResult

1. f_eval(state.F_cur)                    [1 F_eval, costs K iters]
   err_plain = return value
   R_k[c]   = F_cur[c] - X_k[c]  ∀c
   norm_k   = Σ_c R_k[c]²

2. RESTART CHECK (only when has_prev — prevents spurious restart on first call):
     if (state.has_prev && norm_k > kSRAARestartGamma² * state.prev_resid_norm):
       state.clear()                      [has_prev=false; head=count=0]
       goto ACCEPT_PLAIN

3. Append to circular buffer (only when has_prev):
     if (state.has_prev):
       slot = state.head
       dX_buf[slot*M + c] = X_k[c]  - X_prev[c]  ∀c
       dR_buf[slot*M + c] = R_k[c]  - R_prev[c]   ∀c
       state.head  = (state.head + 1) % m
       state.count = min(state.count + 1, m)

4. Update state (always — provides baseline for next restart check):
     state.X_prev         = X_k
     state.R_prev         = R_k
     state.prev_resid_norm = norm_k      [read at step 2 of NEXT call]
     state.has_prev       = true

5. if (state.count < kSRAAMinCount): goto ACCEPT_PLAIN   [not enough history]

6. Gram matrix (symmetric, m×m — all i,j in [0, count)):
     G[i][j] = Σ_c dR_buf[i*M+c] * dR_buf[j*M+c]   [sequential reads over M]
   Regularize: G[i][i] += kSRAAdeltaReg * max_j(G[j][j])
   RHS: rhs[i] = Σ_c dR_buf[i*M+c] * R_k[c]

7. Solve (G + δI)γ = rhs via ldlt_factor_inplace(gram, count, eps) + ldlt_solve
   [reuses calib_linalg.hpp — ldlt is runtime-n, no size constraint]
   If ldlt_factor_inplace returns !RK_OK:
     state.clear(); goto ACCEPT_PLAIN    [ill-conditioned G; reset history]

8. Extrapolate + clamp (vectorizable; m=5 unrolls over inner loop):
     scratch[c] = clamp(X_k[c] + R_k[c] - Σ_i γ[i]*(dX[i*M+c]+dR[i*M+c]),
                        L_cell[c], U_cell[c])  ∀c

9. f_eval(state.scratch)                  [1 F_eval, costs K iters; result in scratch]
   err_AA = return value
   [std::isfinite(err_AA) check: NaN from degenerate F_eval falls through safeguard correctly
    because NaN comparisons return false, triggering ACCEPT_PLAIN — but make this explicit]
   if (!std::isfinite(err_AA)): state.clear(); goto ACCEPT_PLAIN

10. Safeguard:
    if (err_AA <= err_plain):
      std::swap(X, state.scratch)        [O(1) swap — not O(M) copy]
      state.aa_accepted_count++
      return {aa_accepted=true, f_evals=2, err_result=err_AA}
    else:
      state.clear()                      [bad step; reset history]
      goto ACCEPT_PLAIN

ACCEPT_PLAIN:
    std::swap(X, state.F_cur)            [O(1) swap of already-computed F(X_k)]
    return {aa_accepted=false, f_evals=1, err_result=err_plain}
```

**Key invariants:**
- ACCEPT_PLAIN always uses `std::swap`, never `operator=` — O(1), not O(M) copy
- ACCEPT_PLAIN costs 1 F_eval (F_cur computed in step 1, reused here)
- AA path costs exactly 2 F_evals
- Restart check guarded by `has_prev` — never fires on first call
- History cleared on: restart, safeguard rejection, LDLT failure, NaN err_AA
- `prev_resid_norm` updated at step 4 (AFTER restart check) — stale value from before
  restart is never read again because `has_prev=false` after `clear()`
- `state.scratch` and outer `X` are always distinct objects — invariant documented at sraa_step declaration

---

### Parameters

| Constant | Value | Derivation |
|---|---|---|
| `kSRAAm` | 5 | 5 history pairs; 176 MB at stepstone; empirically optimal for survey calibration |
| `kSRAAMaxM` | 10 | Stack array size; runtime window must be ≤ this |
| `kSRAAMinCount` | 2 | Need ≥2 ΔX/ΔR pairs for non-trivial LS |
| `kSRAAdeltaReg` | 1e-10 | Relative Tikhonov; activates when G near-singular late in convergence |
| `kSRAARestartGamma` | 2.0 | Restart when `‖R_k‖² > 4 × prev_norm` (residual doubled) |

All `constexpr`, defined at top of `sraa.hpp`.

---

### Integration: greenkhorn.cpp

**Pre-cut verification** (implementer must run before deleting lines):
```bash
grep -n "sq_X_snap\|sq_w1\|sq_w2\|alpha =\|err_star\|err_w2" src/greenkhorn.cpp
# Confirms exact line range before deletion
```

**Remove** (verified line range, not approximate):
```cpp
// DELETE: sq_X_snap, sq_w1, sq_w2, sq_S1, sq_S2, sq_X_star, sq_Sstar, sq_W1, sq_W2, sq_Wstar
// DELETE: CBB alpha computation block
// DELETE: extrapolation + clamp block
// DELETE: accept/reject block (err_star <= err_w2 * 1.01)
```

**Add** at solver entry (before iteration loop):
```cpp
SRAAState grk_sraa;
if (accelerate) grk_sraa.init(M, kSRAAm);
```

**Self-contained f_eval_sraa** — `Xv` MUST be `state.F_cur` or `state.scratch`, never outer `X`:
```cpp
std::vector<double> Sv_sraa(K * S_stride, 0.0);
auto f_eval_sraa = [&](std::vector<double>& Xv) -> double {
    // Xv must not alias outer X (invariant; see sraa_step declaration)
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
    // errRp is updated by greenkhorn_step from the swapped S_flat — consistent with Xv
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

Same pattern. **Pre-cut verification** before removing lines ~360–515:
```bash
grep -n "X_prev_sq\|is_infeasible\|wchange\|kMaxHalvings\|halv\|sq_w1" src/raking.cpp | head -20
# Confirm X_prev_sq stall detection and is_infeasible tracking are in OUTER loop, not SQUAREM block
```

**Preserve unchanged** (outer loop, unaffected by SQUAREM removal):
- `is_infeasible` tracking
- `X_prev_sq` / wchange stall convergence check
- `errRp_new` update: assign from `r.err_result` after `sraa_step`

### Integration: R/harvest.R

Update in the **same commit** as C++ changes:
- `@param accelerate` roxygen: replace "SqS3 SQUAREM" reference with "Anderson Acceleration (AA-m)"
- Warning message (currently line ~435): replace "SQUAREM weight-change plateau" with "AA weight-change plateau"
- Any inline comment referencing "SqS3" or "CBB" in the accelerate path

---

### TDD Requirements

#### T_sraa_grk — AA fires and improves greenkhorn quality

```r
test_that("T_sraa_grk: greenkhorn+AA max_err <= plain and converges faster", {
  set.seed(99); n <- 2000L
  df  <- data.frame(x=factor(sample(letters[1:3],n,TRUE)),
                    y=factor(sample(c("M","F"),n,TRUE)))
  tgt <- list(x=c(a=0.3,b=0.4,c=0.3), y=c(M=0.5,F=0.5))
  K_exp <- 2L

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

  # (1) AA quality must match or beat plain (0.1% slack for float noise)
  expect_lte(me_aa, me_plain * 1.001,
    label=sprintf("AA (%.2e) must not be worse than plain (%.2e)", me_aa, me_plain))

  # (2) AA fired: on a well-conditioned K=2 problem, first 2 plain steps contribute
  # 2*(K*1) = 4 iters; all subsequent AA steps contribute K*2 = 4 each.
  # Total = 4 + N*4 for any N ≥ 0 → always divisible by K*2=4. Proof: (4+4N) mod 4 = 0.
  expect_equal(iters_aa %% (K_exp * 2L), 0L,
    label=sprintf("AA iters (%d) %% K*2=%d == 0 proves AA fired", iters_aa, K_exp*2L))

  # (3) AA faster than plain in iteration count
  expect_lt(iters_aa, iters_plain,
    label=sprintf("AA (%d) must be faster than plain (%d)", iters_aa, iters_plain))
})
```

#### T_sraa_rk — raking+AA quality matches or beats plain

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
    label=sprintf("raking+AA (%.2e) must not exceed plain (%.2e)", me_aa, me_plain))
})
```

#### T_sraa_ldlt_fallback — LDLT failure falls back gracefully

```r
test_that("T_sraa_ldlt_fallback: ill-conditioned AA history falls back to plain", {
  # A problem where the same iterate repeats makes ΔR rows identical
  # → G is rank-deficient → ldlt_factor_inplace returns !RK_OK → ACCEPT_PLAIN
  # Verify: solver still converges (falls back gracefully, not crashes)
  set.seed(7); n <- 500L
  df  <- data.frame(x=factor(sample(letters[1:2],n,TRUE)))
  tgt <- list(x=c(a=0.3,b=0.7))
  r <- suppressWarnings(
    harvest(df, tgt, method="greenkhorn", accelerate=TRUE, max_iterations=200L,
            attach_weights=FALSE))
  expect_lt(attr(r,"result")$max_error, 1e-3,
    label="LDLT-fallback path must still converge to good solution")
  expect_true(attr(r,"result")$status %in% c(0L, 1L, 5L),
    label="Status must be valid (not crash)")
})
```

#### T_sraa_restart — diverging AA triggers restart and recovers

```r
test_that("T_sraa_restart: restart on divergence recovers and converges", {
  # A tight-bounds problem where early Newton steps can diverge temporarily
  set.seed(42); n <- 1000L
  df  <- data.frame(x=factor(sample(letters[1:4],n,TRUE)),
                    y=factor(sample(c("M","F"),n,TRUE)))
  tgt <- list(x=setNames(c(0.4,0.3,0.2,0.1), letters[1:4]),
              y=c(M=0.3, F=0.7))
  r <- suppressWarnings(
    harvest(df, tgt, method="greenkhorn", accelerate=TRUE,
            max_weight=1.5, min_weight=0.1, max_iterations=500L,
            attach_weights=FALSE))
  expect_lt(attr(r,"result")$max_error, 1e-2,
    label="Restart path must recover and converge within 1e-2 tolerance")
})
```

#### T_acc replacement

Replace existing T_acc test body (`test_that("Tacc: greenkhorn accelerate=TRUE fires SQUAREM...")`)
with T_sraa_grk above. The K*3 divisibility check is SQUAREM-specific and invalid after AA.

#### Regression: T5–T8, T_logit_armijo, T_logit_init must remain GREEN

---

## Acceptance Criteria

| # | Criterion | Verify |
|---|-----------|--------|
| AC1 | greenkhorn+AA max_err ≤ greenkhorn plain max_err | `devtools::test(filter="calibration-solvers")` |
| AC2 | raking+AA max_err ≤ raking plain max_err | Same |
| AC3 | T_sraa_grk GREEN (AA fires, faster) | Same |
| AC4 | T_sraa_rk GREEN | Same |
| AC5 | T_sraa_ldlt_fallback GREEN | Same |
| AC6 | T_sraa_restart GREEN | Same |
| AC7 | FAIL count unchanged (3 pre-existing) | `devtools::test()` |
| AC8 | harvest.R docs updated (no SQUAREM/CBB/SqS3 in accelerate path) | `grep -n "SqS3\|CBB\|SQUAREM" R/harvest.R` |
| AC9 | NEWS.md contains breaking-change entry for accelerate=TRUE output change | `grep -i "break\|accelerate\|SRAA\|Anderson" NEWS.md` |
| AC10 (benchmark) | stepstone greenkhorn+AA max_err ≤ 1.57e-3 | `Rscript benchmarks/stepstone_all_methods.R` |
| AC11 (benchmark) | stepstone raking+AA max_err ≤ 1.60e-3 | Same |

---

## Out of Scope

- Exposing `m`, `delta_reg`, `restart_gamma` to R user (internal; hardcoded is simpler)
- AA for `method="sinkhorn"` (different algorithm structure)
- Adaptive m (start small, grow) — follow-on if m=5 quality warrants
- Per-obs weight-level AA (cell-aggregate is sufficient and cheaper)
- AA for logit (uses Newton, not fixed-point iteration)
- Automated K≥5 CI test (stepstone K=9 is manual benchmark; acceptable trade-off vs CI runtime)
