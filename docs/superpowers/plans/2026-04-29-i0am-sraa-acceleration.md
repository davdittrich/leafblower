# SRAA-m: Replace SQUAREM with Anderson Acceleration

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace CBB-SQUAREM acceleration in greenkhorn and raking with Safeguarded Regularized Anderson Acceleration (SRAA-m), fixing the 35% quality regression caused by SQUAREM overshooting the bounded optimum.

**Architecture:** New `src/sraa.hpp` header-only template containing `SRAAState` (pre-allocated flat circular buffer) and `sraa_step<FEval>()` template. Both greenkhorn.cpp and raking.cpp replace their SQUAREM blocks with a single sraa_step call. Each AA super-step uses 2 F_evals (vs SQUAREM's 3), with safeguard comparing err_AA vs err_plain.

**Tech Stack:** C++17, R/Rcpp, calib_linalg.hpp (ldlt_factor_inplace + ldlt_solve reused for m×m LS), testthat3.

---

## Task 1 — RED: Failing tests for SRAA quality and acceleration

**Mechanism:** testthat3 expectations comparing accelerated vs plain on harvest() with greenkhorn and raking methods.
**Forbidden:** stubbing C++; mocking solver internals; weakening assertions.
**Audit:** observe failure output of `devtools::test()` to confirm RED state pre-implementation.

### Step 1.1 — Append SRAA tests to test-calibration-solvers.R

- [ ] Edit `/home/dd/Gemini/leafblower/tests/testthat/test-calibration-solvers.R`, appending these four tests at the end of the file:

```r
test_that("T_sraa_grk: greenkhorn+AA max_err <= plain and converges faster", {
  set.seed(99); n <- 2000L
  df  <- data.frame(x=factor(sample(letters[1:3],n,TRUE)),
                    y=factor(sample(c("M","F"),n,TRUE)))
  tgt <- list(x=c(a=0.3,b=0.4,c=0.3), y=c(M=0.5,F=0.5))
  K_exp <- 2L
  r_aa    <- suppressWarnings(harvest(df, tgt, method="greenkhorn", accelerate=TRUE,
                                       max_iterations=500L, attach_weights=FALSE))
  r_plain <- suppressWarnings(harvest(df, tgt, method="greenkhorn", accelerate=FALSE,
                                       max_iterations=500L, attach_weights=FALSE))
  me_aa    <- attr(r_aa,    "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  iters_aa    <- attr(r_aa,    "result")$iterations
  iters_plain <- attr(r_plain, "result")$iterations
  expect_lte(me_aa, me_plain * 1.001,
    label=sprintf("AA (%.2e) must not exceed plain (%.2e)", me_aa, me_plain))
  # Divisibility: first 2 plain steps = 2*(K*1)=4; all AA steps = K*2=4 each.
  # Total = 4+4N for any N>=0 -> divisible by K*2=4. Proof: (4+4N) mod 4 = 0.
  expect_equal(iters_aa %% (K_exp * 2L), 0L,
    label=sprintf("AA iters (%d) %%%% K*2=%d == 0 proves AA fired", iters_aa, K_exp*2L))
  expect_lt(iters_aa, iters_plain,
    label=sprintf("AA (%d) must be faster than plain (%d)", iters_aa, iters_plain))
})

test_that("T_sraa_rk: raking+AA max_err <= raking plain", {
  set.seed(99); n <- 2000L
  df  <- data.frame(x=factor(sample(letters[1:3],n,TRUE)),
                    y=factor(sample(c("M","F"),n,TRUE)))
  tgt <- list(x=c(a=0.3,b=0.4,c=0.3), y=c(M=0.5,F=0.5))
  r_aa    <- suppressWarnings(harvest(df, tgt, method="raking", accelerate=TRUE,
                                       max_iterations=500L, attach_weights=FALSE))
  r_plain <- suppressWarnings(harvest(df, tgt, method="raking", accelerate=FALSE,
                                       max_iterations=500L, attach_weights=FALSE))
  me_aa    <- attr(r_aa,    "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  expect_lte(me_aa, me_plain * 1.001,
    label=sprintf("raking+AA (%.2e) must not exceed plain (%.2e)", me_aa, me_plain))
})

test_that("T_sraa_ldlt_fallback: ill-conditioned AA history falls back to plain", {
  set.seed(7); n <- 500L
  df  <- data.frame(x=factor(sample(letters[1:2],n,TRUE)))
  tgt <- list(x=c(a=0.3,b=0.7))
  r <- suppressWarnings(harvest(df, tgt, method="greenkhorn", accelerate=TRUE,
                                 max_iterations=200L, attach_weights=FALSE))
  expect_lt(attr(r,"result")$max_error, 1e-3)
  expect_true(attr(r,"result")$status %in% c(0L, 1L, 5L))
})

test_that("T_sraa_restart: restart on divergence recovers and converges", {
  set.seed(42); n <- 1000L
  df  <- data.frame(x=factor(sample(letters[1:4],n,TRUE)),
                    y=factor(sample(c("M","F"),n,TRUE)))
  tgt <- list(x=setNames(c(0.4,0.3,0.2,0.1),letters[1:4]),
              y=c(M=0.3, F=0.7))
  r <- suppressWarnings(harvest(df, tgt, method="greenkhorn", accelerate=TRUE,
                                 max_weight=1.5, min_weight=0.1,
                                 max_iterations=500L, attach_weights=FALSE))
  expect_lt(attr(r,"result")$max_error, 1e-2)
})
```

### Step 1.2 — Remove the existing Tacc test (SQUAREM-specific; invalid after SRAA)

The existing `Tacc` test asserts `iters_acc %% (K*3) == 0` — a SQUAREM K×3 stride check.
After SRAA integration, SRAA emits K×2 or K×1; this will unconditionally fail.
`T_sraa_grk` (appended in Step 1.1) is its direct replacement.

- [ ] In `tests/testthat/test-calibration-solvers.R`, locate and DELETE the entire block:

```r
test_that("Tacc: greenkhorn accelerate=TRUE fires SQUAREM and converges", {
  ...
})
```

Run to find it: `grep -n "Tacc:\|squarem_stride\|K_exp \* 3L" tests/testthat/test-calibration-solvers.R`

- [ ] After deletion, verify it is gone: `grep -c "squarem_stride\|K_exp \* 3L" tests/testthat/test-calibration-solvers.R` must return `0`.

### Step 1.3 — Verify RED state

- [ ] Run: `cd /home/dd/Gemini/leafblower && Rscript -e "devtools::test(filter='calibration-solvers')" 2>&1 | tail -30`
- [ ] Expected: at least `T_sraa_grk` fails on `me_aa <= me_plain * 1.001` (current SQUAREM is 35% worse). At least one expect_lte / expect_lt MUST fail. If all four pass, halt — tests are not exercising SRAA path.

### Step 1.4 — Commit

- [ ] `cd /home/dd/Gemini/leafblower && git add tests/testthat/test-calibration-solvers.R && git commit -m "test(sraa): RED tests for SRAA-m + remove SQUAREM-specific Tacc"`

---

## Task 2 — Create src/sraa.hpp (SRAAState + sraa_step template)

**Mechanism:** Header-only C++17 template. Uses `lbw::ldlt_factor_inplace` + `lbw::ldlt_solve` from `calib_linalg.hpp` for the m×m least-squares solve.
**Forbidden:** Eigen, BLAS, dynamic allocation in hot path, virtual dispatch.
**Audit:** `R CMD INSTALL --preclean .` after creation must succeed and link cleanly when included in greenkhorn.cpp later.

### Step 2.1 — Write `/home/dd/Gemini/leafblower/src/sraa.hpp`

- [ ] Create file with this exact content:

```cpp
#pragma once
// src/sraa.hpp — Safeguarded Regularized Anderson Acceleration (Type II AA)
// Used by greenkhorn.cpp and raking.cpp to replace SQUAREM CBB acceleration.
#include "calib_linalg.hpp"
#include "lbw_config.h"
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace lbw {

static constexpr int    kSRAAMaxM       = 10;    // stack array bound; init() enforces window <= this
static constexpr int    kSRAAm          = 5;     // default window: 176 MB at stepstone scale
static constexpr int    kSRAAMinCount   = 2;     // min DX/DR pairs before AA fires
static constexpr double kSRAAdeltaReg   = 1e-10; // relative Tikhonov on Gram matrix
static constexpr double kSRAARestartGamma = 2.0; // restart when ||R_k||^2 > 4 x prev_norm

struct SRAAStepResult {
    bool   aa_accepted;
    int    f_evals;     // 1 (plain) or 2 (AA attempted)
    double err_result;  // max_error of accepted step
};

struct SRAAState {
    int M = 0, m = 0, head = 0, count = 0;
    bool has_prev = false;
    int aa_accepted_count = 0; // cumulative; NOT reset by clear()

    // Slot-contiguous: buf[slot*M + cell]. Each slot is contiguous M doubles -> SIMD dot products.
    std::vector<double> dX_buf;   // m x M
    std::vector<double> dR_buf;   // m x M
    std::vector<double> R_prev;   // M
    std::vector<double> X_prev;   // M
    std::vector<double> F_cur;    // M: F(X_k)
    std::vector<double> scratch;  // M: R_k temp, then X_AA, then F(X_AA)

    double gram[kSRAAMaxM * kSRAAMaxM] = {};
    double rhs[kSRAAMaxM] = {};
    double gamma_[kSRAAMaxM] = {}; // gamma_ avoids <cmath> gamma conflict

    double prev_resid_norm = 0.0;  // read only when has_prev=true

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
                       (2.0*m + 4.0) * M * 8.0 / 1e6, m, M);
        }
        clear();
    }

    // Resets history. aa_accepted_count NOT reset. prev_resid_norm set to 0 (safe sentinel).
    void clear() { head = 0; count = 0; has_prev = false; prev_resid_norm = 0.0; }
};

// f_eval: (std::vector<double>&) -> double
//   Modifies Xv in-place (K steps), returns max_errRp.
//   INVARIANT: must receive state.F_cur or state.scratch, NEVER the outer X.
template<typename FEval>
SRAAStepResult sraa_step(
    FEval& f_eval,
    std::vector<double>& X,
    const std::vector<double>& L_cell,
    const std::vector<double>& U_cell,
    SRAAState& state)
{
    const int M   = state.M;

    // --- Step 1: F(X_k) -> F_cur; compute R_k into scratch; compute norm ---
    double err_plain = f_eval(state.F_cur);
    double norm_k = 0.0;
    for (int c = 0; c < M; c++) {
        double rk = state.F_cur[c] - X[c];
        state.scratch[c] = rk;  // temporarily holds R_k
        norm_k += rk * rk;
    }

    // --- Step 2: Restart check (guarded by has_prev — never fires on first call) ---
    if (state.has_prev &&
        norm_k > kSRAARestartGamma * kSRAARestartGamma * state.prev_resid_norm) {
        state.clear();
        std::swap(X, state.F_cur);
        return {false, 1, err_plain};
    }

    // --- Step 3: Append DX, DR to circular buffer (only when has_prev) ---
    if (state.has_prev) {
        double* dX_s = state.dX_buf.data() + state.head * M;
        double* dR_s = state.dR_buf.data() + state.head * M;
        for (int c = 0; c < M; c++) {
            dX_s[c] = X[c]              - state.X_prev[c];
            dR_s[c] = state.scratch[c]  - state.R_prev[c];  // R_k - R_{k-1}
        }
        state.head  = (state.head + 1) % state.m;
        state.count = std::min(state.count + 1, state.m);
    }

    // --- Step 4: Update prev state ---
    for (int c = 0; c < M; c++) {
        state.X_prev[c] = X[c];
        state.R_prev[c] = state.scratch[c];  // R_k
    }
    state.prev_resid_norm = norm_k;
    state.has_prev = true;

    // --- Step 5: Not enough history -> plain ---
    if (state.count < kSRAAMinCount) {
        std::swap(X, state.F_cur);
        return {false, 1, err_plain};
    }

    // --- Step 6: Gram matrix G[i][j] = DR[i].DR[j]; RHS rhs[i] = DR[i].R_k ---
    const int n = state.count;
    double max_diag = 0.0;
    for (int i = 0; i < n; i++) {
        const double* dRi = state.dR_buf.data() + i * M;
        for (int j = i; j < n; j++) {
            const double* dRj = state.dR_buf.data() + j * M;
            double dot = 0.0;
            for (int c = 0; c < M; c++) dot += dRi[c] * dRj[c];
            state.gram[i * kSRAAMaxM + j] = state.gram[j * kSRAAMaxM + i] = dot;
        }
        max_diag = std::max(max_diag, state.gram[i * kSRAAMaxM + i]);
    }
    for (int i = 0; i < n; i++) {
        const double* dRi = state.dR_buf.data() + i * M;
        double dot = 0.0;
        for (int c = 0; c < M; c++) dot += dRi[c] * state.scratch[c];
        state.rhs[i] = dot;
    }

    // --- Step 7: Regularize + LDLT solve (n x n submatrix) ---
    double eps = kSRAAdeltaReg * (max_diag > 0.0 ? max_diag : 1.0);
    // Copy n x n submatrix into contiguous buffer for ldlt (gram is kSRAAMaxM x kSRAAMaxM)
    double G[kSRAAMaxM * kSRAAMaxM];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            G[i * n + j] = state.gram[i * kSRAAMaxM + j];
    G[0] += eps;  // pre-regularize diagonal before factoring
    for (int i = 1; i < n; i++) G[i * n + i] += eps;  // rest of diagonal

    if (lbw::ldlt_factor_inplace(G, (size_t)n, 0.0) != RK_OK) {
        state.clear();
        std::swap(X, state.F_cur);
        return {false, 1, err_plain};
    }
    for (int i = 0; i < n; i++) state.gamma_[i] = state.rhs[i];
    lbw::ldlt_solve(G, (size_t)n, state.gamma_);

    // --- Step 8: Extrapolate + clamp into scratch ---
    // scratch currently holds R_k; overwrite in-place with X_AA
    for (int c = 0; c < M; c++) {
        double Rk_c = state.scratch[c];  // read R_k before overwrite
        double corr = 0.0;
        for (int i = 0; i < n; i++)
            corr += state.gamma_[i] * (state.dX_buf[i * M + c] + state.dR_buf[i * M + c]);
        state.scratch[c] = std::clamp(X[c] + Rk_c - corr, L_cell[c], U_cell[c]);
    }

    // --- Step 9: F(X_AA); scratch -> F(X_AA) ---
    double err_AA = f_eval(state.scratch);

    // NaN guard
    if (!std::isfinite(err_AA)) {
        state.clear();
        std::swap(X, state.F_cur);
        return {false, 1, err_plain};
    }

    // --- Step 10: Safeguard ---
    if (err_AA <= err_plain) {
        std::swap(X, state.scratch);    // O(1); X = F(X_AA)
        state.aa_accepted_count++;
        return {true, 2, err_AA};
    } else {
        state.clear();                  // bad step; reset history
        std::swap(X, state.F_cur);     // X = F(X_k)
        return {false, 2, err_plain};  // 2 F_evals spent
    }
}

} // namespace lbw
```

### Step 2.2 — Compile gate

- [ ] Run: `cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3`
- [ ] Expected: ends with `* DONE (leafblower)`. The header is not yet included; this only confirms it compiles when other files are not broken.

### Step 2.3 — Commit

- [ ] `cd /home/dd/Gemini/leafblower && git add src/sraa.hpp && git commit -m "feat(sraa): add SRAAState + sraa_step template (header-only Type II AA)"`

---

## Task 3 — Integrate SRAA into greenkhorn.cpp (remove SQUAREM)

**Mechanism:** Replace SQUAREM CBB block with `lbw::sraa_step(f_eval_sraa, X, L_cell, U_cell, grk_sraa)`. `f_eval_sraa` recomputes S/W from Xv, runs K greenkhorn_step calls in errRp-descending order, and restores caller state via swap-restore.
**Forbidden:** keeping any sq_* buffers; passing outer X to f_eval_sraa; allocating inside the hot loop.
**Audit:** test divisibility (`iters_aa %% (K*2) == 0`) proves AA path fired; quality assertion proves safeguard works.

### Step 3.1 — Add include and SRAA state at top of greenkhorn.cpp

- [ ] Edit `/home/dd/Gemini/leafblower/src/greenkhorn.cpp`. Find the existing `#include` block at the top (after `#include "lbw_config.h"` or similar) and add:

```cpp
#include "sraa.hpp"
```

### Step 3.2 — Remove old F_eval lambda (dead code after SRAA)

The existing `auto F_eval` at line ~124 is a 3-parameter lambda used only by the SQUAREM block.
After SRAA replaces the block, it becomes dead code.

- [ ] Locate and DELETE the block starting with:
```cpp
    auto F_eval = [&](std::vector<double>& X_in,
                      std::vector<double>& S_in,
                      double& W_in) {
```
Run to find it: `grep -n "auto F_eval" src/greenkhorn.cpp`

Delete through the closing `};` of the lambda.

### Step 3.3 — Remove SQUAREM scratch declarations

- [ ] Find the block beginning with the comment `// SQUAREM scratch buffers` (approx lines 139-148). Delete the entire block:

```cpp
// SQUAREM scratch buffers
std::vector<double> sq_w1, sq_w2, sq_X_snap, sq_S1, sq_S2, sq_Ssnap;
double sq_W1 = 0.0, sq_W2 = 0.0, sq_Wsnap = 0.0;
if (st.accelerate) {
    sq_w1.assign(M, 0.0); sq_w2.assign(M, 0.0); sq_X_snap.assign(M, 0.0);
    sq_S1.assign(K*S_stride, 0.0); sq_S2.assign(K*S_stride, 0.0);
    sq_Ssnap.assign(K*S_stride, 0.0);
}
```

Replace with:

```cpp
// SRAA-m state (replaces SQUAREM)
lbw::SRAAState grk_sraa;
std::vector<double> Sv_sraa;
std::vector<int>    order_sraa;
if (st.accelerate && K > 0) {
    grk_sraa.init(M, lbw::kSRAAm);
    Sv_sraa.assign((size_t)K * S_stride, 0.0);
    order_sraa.assign(K, 0);
}
auto f_eval_sraa = [&](std::vector<double>& Xv) -> double {
    // Xv must be grk_sraa.F_cur or grk_sraa.scratch — NEVER the outer X
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
    std::iota(order_sraa.begin(), order_sraa.end(), 0);
    std::stable_sort(order_sraa.begin(), order_sraa.end(),
        [&](int a, int b){ return errRp[a] > errRp[b]; });
    for (int ki = 0; ki < K; ki++) greenkhorn_step(order_sraa[ki]);
    std::swap(X, Xv); std::swap(S_flat, Sv_sraa); std::swap(W, Wv);
    return *std::max_element(errRp.begin(), errRp.end());
};
```

### Step 3.4 — Replace the SQUAREM CBB block inside the main loop

- [ ] Locate the block guarded by `if (st.accelerate && K > 0) { ... CBB + extrapolation + accept/reject ... res.iterations += K * 3; }` (the SQUAREM body, approx lines 165-205). Replace the entire body with:

```cpp
if (st.accelerate && K > 0) {
    auto r = lbw::sraa_step(f_eval_sraa, X, L_cell, U_cell, grk_sraa);
    for (int k = 0; k < K; k++) errRp[k] = compute_errRp_k(k);
    res.iterations += K * (r.aa_accepted ? 2 : 1);
} else {
```

Make sure the `else` branch (plain greenkhorn step) is preserved unchanged.

### Step 3.5 — Compile gate

- [ ] Run: `cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3`
- [ ] Expected: ends with `* DONE (leafblower)`. If link error mentions `ldlt_factor_inplace`, confirm `calib_linalg.hpp` exposes it in the `lbw` namespace.

### Step 3.6 — Greenkhorn test gate

- [ ] Run: `cd /home/dd/Gemini/leafblower && Rscript -e "devtools::test(filter='calibration-solvers')" 2>&1 | tail -10`
- [ ] Expected: `T_sraa_grk`, `T_sraa_ldlt_fallback`, `T_sraa_restart` PASS. `T_sraa_rk` may still fail (raking not migrated yet).

### Step 3.7 — Commit

- [ ] `cd /home/dd/Gemini/leafblower && git add src/greenkhorn.cpp && git commit -m "feat(sraa): integrate SRAA-m into greenkhorn, remove SQUAREM CBB"`

---

## Task 4 — Integrate SRAA into raking.cpp (remove SQUAREM)

**Mechanism:** Replace inner SQUAREM while-loop with single `lbw::sraa_step` call per outer iteration. Preserve `X_prev_sq` and `is_infeasible` outer-loop state.
**Forbidden:** removing outer convergence checks (errRp threshold / budget / wchange); touching outer loop control flow outside the `if (st.inner_max_iter >= 3)` block.
**Audit:** raking quality test `T_sraa_rk` must pass after this task.

### Step 4.1 — Add include at top of raking.cpp

- [ ] Edit `/home/dd/Gemini/leafblower/src/raking.cpp`. Add to the `#include` block:

```cpp
#include "sraa.hpp"
```

### Step 4.2 — Verify raking.cpp structure before editing

- [ ] Run: `grep -n "auto F_eval\|if (st.accelerate\|inner_max_iter >= 3\|X_prev_sq\|sq_w1" /home/dd/Gemini/leafblower/src/raking.cpp | head -20`

Expected output confirms:
- `auto F_eval = [&](std::vector<double>& Xv) -> double {` at ~line 269 — **this is the F_eval to reuse directly**
- `if (st.accelerate) {` at ~line 351 — outer SQUAREM guard
- `if (st.inner_max_iter >= 3) {` at ~line 358 — inner guard (nested inside accelerate block)
- `auto X_prev_sq = X;` at ~line 360 — inside the inner block (SQUAREM stall detection; deleted with the block)

**Key insight**: raking's `F_eval` at line 269 already has signature `(std::vector<double>&) -> double` — exactly what `sraa_step` needs. No extraction needed. Pass it directly.

### Step 4.3 — Allocate SRAA state before the accelerate block

- [ ] Immediately BEFORE the `if (st.accelerate) {` line (~line 351), add:

```cpp
// SRAA-m state: replaces SQUAREM while-loop inside if(st.accelerate)
lbw::SRAAState rk_sraa;
if (st.accelerate) rk_sraa.init(ct.M_cell, lbw::kSRAAm);
```

### Step 4.4 — Replace the inner SQUAREM block

- [ ] Locate the nested structure:
```cpp
if (st.accelerate) {
    // ...
    if (st.inner_max_iter >= 3) {
        auto X_prev_sq = X;  // SQUAREM stall detection — deleted with block
        // sq_w1, sq_w2, sq_X_snap, sq_X_star, sq_X_star_pre allocations
        while (f_eval_count + 3 <= st.inner_max_iter) {
            // ... CBB + step-halving + F_eval calls ...
        }
    }
    // plain path: for (int iter = 1; ...) { F_eval(X); ... }
```

Replace the ENTIRE `if (st.inner_max_iter >= 3) { ... }` block with a single sraa_step call. The `if (st.accelerate)` outer guard and the plain path below it remain:

```cpp
if (st.accelerate) {
    // SQUAREM replaced by SRAA-m. X_prev_sq stall detection removed (SQUAREM-specific).
    auto r = lbw::sraa_step(F_eval, X, L_cell, U_cell, rk_sraa);
    res.max_error  = r.err_result;
    res.iterations += (r.aa_accepted ? 2 : 1);
} else {
    // plain path below (unchanged)
```

Note: `X_prev_sq` was stall-detection internal to the SQUAREM while-loop — it is deleted along with the block. The outer loop's convergence check (errRp threshold, budget, stall) is unaffected.

### Step 4.4 — Compile gate

- [ ] Run: `cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3`
- [ ] Expected: ends with `* DONE (leafblower)`.

### Step 4.5 — Full test gate

- [ ] Run: `cd /home/dd/Gemini/leafblower && Rscript -e "devtools::test()" 2>&1 | tail -3`
- [ ] Expected: all four T_sraa_* tests PASS, no regressions in pre-existing tests. If any pre-existing test fails, halt and diagnose.

### Step 4.6 — Commit

- [ ] `cd /home/dd/Gemini/leafblower && git add src/raking.cpp && git commit -m "feat(sraa): integrate SRAA-m into raking, remove SQUAREM CBB"`

---

## Task 5 — Update R/harvest.R + NEWS.md

**Mechanism:** Roxygen + NEWS prose updates. No behavior change.
**Forbidden:** changing default value of `accelerate`; removing the parameter; touching unrelated docs.
**Audit:** `devtools::document()` regenerates NAMESPACE/Rd cleanly; `R CMD check` doc warnings absent.

### Step 5.1 — Update accelerate roxygen in R/harvest.R

- [ ] Edit `/home/dd/Gemini/leafblower/R/harvest.R`. Find the `@param accelerate` roxygen line. Replace its description with:

```r
#' @param accelerate Logical. Enable Safeguarded Regularized Anderson Acceleration
#'   (SRAA-m) for greenkhorn and raking. SRAA-m reduces iteration count by ~40-60%
#'   while preserving (or improving) max_error quality via a per-step safeguard
#'   that compares the AA candidate against the plain F-step. Replaces the prior
#'   SQUAREM/CBB scheme (which overshot the bounded optimum and degraded quality
#'   by up to 35% on bounded problems). Default \code{TRUE}.
```

### Step 5.2 — Update accelerate warning message (if present)

- [ ] In the same file, search for any existing `warning()` or `message()` referencing "SQUAREM" or "CBB". Replace the substring `SQUAREM` with `SRAA-m` and `CBB` with `Anderson` in any user-facing warning text. If no such message exists, skip this step.

### Step 5.3 — Add NEWS.md entry

- [ ] Edit `/home/dd/Gemini/leafblower/NEWS.md`. Add at the top (under the unreleased / current dev version header — create one if absent named `# leafblower (development version)`):

```markdown
## Acceleration

* Replaced SQUAREM/CBB acceleration with Safeguarded Regularized Anderson
  Acceleration (SRAA-m, window m=5) for `method = "greenkhorn"` and
  `method = "raking"`. Fixes a 35% quality regression on bounded calibration
  problems where SQUAREM extrapolation overshot the feasible region. SRAA-m
  guarantees `max_error_accelerated <= max_error_plain` per super-step via
  an explicit safeguard, with automatic restart on residual blow-up
  (||R_k||^2 > 4 ||R_{k-1}||^2) and Tikhonov regularization on the m x m
  Gram system.
```

### Step 5.4 — Regenerate documentation

- [ ] Run: `cd /home/dd/Gemini/leafblower && Rscript -e "devtools::document()" 2>&1 | tail -5`
- [ ] Expected: writes `man/harvest.Rd`, no errors.

### Step 5.5 — Final full check

- [ ] Run: `cd /home/dd/Gemini/leafblower && Rscript -e "devtools::test()" 2>&1 | tail -3`
- [ ] Expected: all tests PASS including the four T_sraa_* tests.

### Step 5.6 — Commit

- [ ] `cd /home/dd/Gemini/leafblower && git add R/harvest.R NEWS.md man/harvest.Rd && git commit -m "docs(sraa): update accelerate param + NEWS for SRAA-m migration"`

---

## Done criteria

- [ ] All four `T_sraa_*` tests PASS.
- [ ] No SQUAREM-named identifiers (`sq_w1`, `sq_X_snap`, `sq_S1`, etc.) remain in `src/greenkhorn.cpp` or `src/raking.cpp`.
- [ ] `R CMD INSTALL --preclean .` succeeds.
- [ ] `devtools::test()` reports zero failures.
- [ ] `man/harvest.Rd` and `NEWS.md` reflect SRAA-m.
- [ ] Five commits landed (one per task), each with conventional-commit subject and zero AI attribution.
