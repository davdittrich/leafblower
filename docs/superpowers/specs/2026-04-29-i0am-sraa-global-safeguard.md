# i0am-B: SRAA-m Global Safeguard + Revert-to-Best

**Date**: 2026-04-29
**Status**: Approved for implementation
**Ticket**: leafblower-i0am (continued)
**Files**: `src/sraa.hpp`, `tests/testthat/test-calibration-solvers.R`

---

## Problem

Current SRAA-m uses a **per-step local safeguard**: accept AA step if `err_AA ≤ err_plain`.
This prevents single-step degradation but not trajectory drift to a suboptimal fixed point.

On K=9 stepstone: SRAA finds a KL-optimal but max_err-suboptimal basin (max_err 2.12e-3 vs
plain 1.57e-3, marg_kl 5.93e-4 vs 8.87e-3). Each AA step is locally acceptable, but the
cumulative trajectory escapes the max_err-optimal basin that plain greenkhorn's greedy-sort
strategy finds.

**Root cause**: `err_plain` (one plain step ahead) is a loose baseline — it doesn't reflect the
quality of the globally optimal trajectory.

**Transfer from revss_temp**: That project's Aitken guard was over-rejecting valid steps
(criterion too tight in one direction). Our safeguard is under-rejecting invalid steps
(criterion too loose in the other direction). Both stem from using a wrong LOCAL criterion.
Fix: track the best quality seen and gate acceptance against it — analogous to accepting
Aitken candidates unconditionally if they improve on the best known value.

---

## Design

### Architecture

Two additions to `src/sraa.hpp` only. No R changes, no C++ solver changes.

```
src/sraa.hpp    ← add best_err_seen, X_best, stall_count to SRAAState;
                  update init(), clear();
                  change safeguard in sraa_step; add best-tracking + revert logic
tests/testthat/test-calibration-solvers.R  ← add T_sraa_global
```

---

### Changes to SRAAState

**New fields:**
```cpp
// Global quality floor — NOT reset by clear() — only improves across entire solver run
double best_err_seen = std::numeric_limits<double>::infinity();
std::vector<double> X_best;  // iterate achieving best_err_seen; pre-allocated in init()
int stall_count = 0;          // steps since last improvement in best_err_seen
```

**New constants:**
```cpp
static constexpr double kSRAAGlobalEps   = 1e-3;   // 0.1% slack on global safeguard
static constexpr int    kSRAAStallWindow = 15;     // stall steps before revert+restart
```

**`init()` additions:**
```cpp
X_best.assign(M, 0.0);
best_err_seen = std::numeric_limits<double>::infinity();
stall_count   = 0;
```

**`clear()` updated** — resets history + stall counter but preserves quality floor:
```cpp
void clear() {
    head = 0; count = 0; has_prev = false; prev_resid_norm = 0.0;
    stall_count = 0;
    // best_err_seen and X_best NOT reset — quality floor persists across all restarts.
    // After clear(), the first plain steps must beat best_err_seen*(1+kSRAAGlobalEps)
    // before AA can fire again — prevents regression to a worse basin.
}
```

`aa_accepted_count` also NOT reset (unchanged from current).

Memory delta: +1 × M doubles for `X_best` = **+12.6 MB at stepstone** (total ~189 MB).

---

### Changes to `sraa_step`

**Change 1 — Replace safeguard comparison (Step 10):**

```cpp
// OLD:
if (err_AA <= err_plain)

// NEW:
if (err_AA <= state.best_err_seen * (1.0 + kSRAAGlobalEps))
```

`err_plain` is no longer the acceptance criterion. AA is accepted only when it can match
the best quality seen across the entire solver run (with 0.1% slack). Plain steps always
execute unconditionally and update `best_err_seen`.

Note: `err_plain` is still computed (needed for NaN guard and ACCEPT_PLAIN path) but
removed from the safeguard comparison.

**Change 2 — Best-tracking + revert-to-best (after every accepted step):**

After the `std::swap` that sets X to the accepted iterate — applies to ALL exit paths
(ACCEPT_PLAIN, restart, LDLT failure, AA accept, AA reject):

```cpp
// Track best quality. Applies after X is set to the accepted iterate.
double accepted_err = r.err_result;  // the quality of X after acceptance
if (accepted_err < state.best_err_seen) {
    state.best_err_seen = accepted_err;
    state.X_best        = X;    // O(M) copy — new best position
    state.stall_count   = 0;
} else {
    state.stall_count++;
}

// Revert-to-best when stalled: no improvement for kSRAAStallWindow steps.
// X_best is the best iterate seen across the entire run.
if (state.stall_count >= kSRAAStallWindow &&
    state.best_err_seen < std::numeric_limits<double>::infinity()) {
    std::swap(X, state.X_best);  // O(1): X = best position, X_best = stale old X
    state.clear();               // reset AA history + stall_count; preserves best_err_seen
}
```

Revert fires at most once per `kSRAAStallWindow` steps (clear() resets stall_count to 0).
After revert, next iteration starts from X_best position with fresh AA history.

**Callers (greenkhorn.cpp, raking.cpp) unchanged.** The `grk_sraa.F_cur = X` seeding
before each sraa_step call is still required and sufficient.

---

### Semantics of best_err_seen

`best_err_seen` is **monotonically non-increasing** across the entire solver run:
- Initialized to `+∞` by `init()` (cold start)
- Never reset by `clear()` (quality floor persists across AA restarts)
- Updated on every accepted step (both plain and AA)

Effect: `best_err_seen` after the first few plain steps reflects the quality of the
max_err-optimal trajectory. AA can only fire if it stays in or returns to that trajectory's
quality basin. When AA extrapolates to a KL-better but max_err-worse basin, `err_AA`
exceeds `best_err_seen * 1.001` → rejected → plain steps continue → basin recovered.

---

## TDD Requirements

### Step 1 — RED gate: verify T_sraa_global is RED with current sraa.hpp

Run BEFORE implementing: `Rscript -e "devtools::test(filter='calibration-solvers')" 2>&1 | grep T_sraa_global`

Expected: T_sraa_global FAILS (old safeguard overshoots on K=4).

### T_sraa_global — New test

```r
test_that("T_sraa_global: greenkhorn+SRAA max_err <= plain on K=4 overlapping-margin problem", {
  # K=4 designed to reproduce the basin-escape failure seen at K=9 (stepstone).
  # Old local safeguard: AA overshoots max_err-optimal basin on multi-margin problems.
  # New global safeguard: AA stays in or returns to max_err-optimal basin.
  set.seed(5); n <- 3000L
  df <- data.frame(
    a = factor(sample(letters[1:4], n, TRUE)),
    b = factor(sample(LETTERS[1:3], n, TRUE)),
    c = factor(sample(c("x","y"),   n, TRUE)),
    d = factor(sample(c("M","F"),   n, TRUE))
  )
  tgt <- list(
    a = setNames(c(0.4, 0.3, 0.2, 0.1), letters[1:4]),
    b = setNames(c(0.5, 0.3, 0.2),      LETTERS[1:3]),
    c = c(x=0.6, y=0.4),
    d = c(M=0.45, F=0.55)
  )
  r_aa    <- suppressWarnings(harvest(df, tgt, method="greenkhorn", accelerate=TRUE,
                                       max_iterations=500L, attach_weights=FALSE))
  r_plain <- suppressWarnings(harvest(df, tgt, method="greenkhorn", accelerate=FALSE,
                                       max_iterations=500L, attach_weights=FALSE))
  me_aa    <- attr(r_aa,    "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  expect_lte(me_aa, me_plain * 1.001 + 1e-10,
    label=sprintf("SRAA K=4 (%.2e) must not exceed plain (%.2e)", me_aa, me_plain))
})
```

### Existing tests must remain GREEN

T_sraa_grk, T_sraa_rk, T_sraa_ldlt_fallback, T_sraa_restart, T5–T8, T_logit_*.

---

## Acceptance Criteria

| # | Criterion | Verify |
|---|-----------|--------|
| AC1 | T_sraa_grk GREEN | `devtools::test(filter="calibration-solvers")` |
| AC2 | T_sraa_rk GREEN | Same |
| AC3 | T_sraa_global GREEN (new K=4) | Same |
| AC4 | FAIL count = 3 (unchanged) | `devtools::test()` |
| AC5 (benchmark) | stepstone greenkhorn+SRAA max_err ≤ 1.57e-3 | `Rscript benchmarks/stepstone_all_methods.R` |
| AC6 (benchmark) | stepstone raking+SRAA max_err ≤ 1.60e-3 | Same |

---

## Calibration of constants

| Constant | Value | Derivation |
|---|---|---|
| `kSRAAGlobalEps` | 1e-3 | 0.1% slack; prevents over-rejection of equivalent-quality steps (revss lesson) |
| `kSRAAStallWindow` | 15 | ~5× kSRAAMinCount=2 warmup; enough plain steps for convergence signal before forcing revert |

---

## Out of Scope

- Exposing `best_err_seen`, `stall_count`, `X_best` to R user
- Adaptive `kSRAAStallWindow`
- Sharing `best_err_seen` across greenkhorn+raking multi-method runs
- Applying global safeguard to raking separately (same `sraa.hpp` change covers both)
