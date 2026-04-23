# Plan: L-BFGS-B Audit Fixes (2026-04-23) — rev 3

## Context

Opus audit of `src/lbfgsb_solver.cpp` found 6 actionable items plus 3
supplementary per-iteration copies (lines 381, 413, 414) and an alpha
allocation. The code is mathematically correct and conforming to
Deville-Sarndal (1992) logit calibration. The Armijo sign claim was
re-verified as **not a bug** under the maximization convention — action
there is a clarifying comment only.

Rev 3 addresses iter-2 plan-review-gate blockers: (a) WU-5 stall guard
dropped (scope creep), (b) WU-6 uses `st.log` not REprintf (preserves
verbose-gate contract, no include-chain risk), (c) WU-2 comment reworded
(Armijo failure → bracket narrow), (d) WU-3c aliasing proof added.

## Scope

Primary file: `src/lbfgsb_solver.cpp` (441 lines). Read-only: `src/logit.hpp`,
`src/lbw_math.hpp`, `src/lbw_config.h`, `src/types.hpp`, `src/leafblower.h`.
No API/ABI changes, no header changes, no new dependencies, no new tests.

## Non-scope

- No change to iEPPA, R/Python bindings, link-function math, ALM activation.
- No new tests.
- No signature changes on static helpers.
- No new features (no stall guard, no REprintf bypass).

## Verification gate (all WUs)

- `R CMD INSTALL --preclean .` succeeds after EACH WU.
- `devtools::test()` ends with `[ FAIL 0 | WARN 2 | SKIP 0 | PASS 100 ]`.
- Before final commit: `Rscript benchmarks/autumn_nr_benchmark.R`; wall-clock
  delta vs `master@dcd804e` recorded; no WU accepted if slowdown > 2 %.

---

## WU-1: Precise docstring for sum(w) at convergence

**Target.** `src/lbfgsb_solver.cpp`, lines 435-438 (comment + two assignments).
Replacement preserves assignments verbatim; only comment text changes.

**Current (435-438):**
```cpp
    // ALM inactive. sum(w)≈n holds at convergence when targets sum to 1 per margin;
    // harvest.R normalises post-call as a safety net for non-converged iterates.
    st.alm_lambda = 0.0;
    st.alm_mu     = 0.0;
```

**New (435-438):**
```cpp
    // ALM inactive. Solver-level invariant: at dual convergence ∇phi=0 ⟹
    // for each margin k, Σ_j S_kj = Σ_j T_kj = W (given targets sum to 1),
    // and summing observations once gives sum(w) = W = Σ d_i. Caller contract
    // (harvest.R, _harvest.py): input weights are normalised so Σ d_i = n,
    // therefore sum(w) = n downstream. The /mean(weights) step in the caller
    // is a safety net for non-converged iterates; a no-op at true convergence.
    st.alm_lambda = 0.0;
    st.alm_mu     = 0.0;
```

**Acceptance.** `grep -c "Σ d_i" src/lbfgsb_solver.cpp` returns ≥ 1. Tests pass.

---

## WU-2: Clarifying comment on Armijo maximization convention

**Target.** `src/lbfgsb_solver.cpp`, above lines 243 and 319. Each guard tests
the Armijo condition; on failure (step insufficient) the control flow
narrows the Wolfe bracket (line 244: `alpha_hi = alpha`) or enters the zoom
phase (line 320: `return wolfe_zoom(...)`). The comment describes the test,
not the branch.

**Insertion above line 243 (one line, no logic change):**
```cpp
        // Armijo test (maximization: slope_0 > 0). Fail ⟹ phi_trial did NOT
        // exceed phi_0 by at least c1*α*slope_0 ⟹ step too long ⟹ narrow bracket.
```

**Insertion above line 319 (one line, no logic change):**
```cpp
        // Armijo test (maximization: slope_0 > 0). Fail ⟹ phi_trial did NOT
        // exceed phi_0 by at least c1*α*slope_0 ⟹ step too long ⟹ enter zoom.
```

**Acceptance.** `grep -c "Armijo test" src/lbfgsb_solver.cpp` returns 2.
Tests pass.

---

## WU-3: Eliminate per-iteration vector copies

Three candidate sites in `lbfgsb_solve_inner`. Action differs per site.

### 3a — History push at 407-408: **skipped** (with correct justification)

Sites 407-408 deep-copy `s_new` and `y_new` (`total`-sized) into the deque.
On typical workloads (K = 1-10 margins, avg cat_count 5-20), `total ≈ 5-200`.
Cost per iteration: 2 × `total` doubles = 80-3200 bytes copied. This is
small vs. the O(K·n) gradient pass that dominates each outer step.
A `std::move` + `assign(total, 0.0)` pattern still requires an allocation
on the next iteration (moved-from capacity is unspecified; `assign` to a
moved-from vector generally reallocates). Net: the simplest correct
alternative does not save an allocation. **Skip.**

### 3b — State rotation at 413-414: **`std::swap`** (O(1))

**Current:**
```cpp
        lam = lam_new;
        grad = grad_new;
```

**New:**
```cpp
        // O(1) pointer swap. lam/grad take the new values; lam_new/grad_new
        // retain the old buffers and are fully overwritten on the next
        // Wolfe call (lam_new at 257, 327, 341; grad_new via phi_from_u at
        // 259, 328, 343) before any further read.
        std::swap(lam, lam_new);
        std::swap(grad, grad_new);
```

**Safety trace.** `lam_new` and `grad_new` are read at lines 396-400 (to
compute `s_new`, `y_new`, and `sy`) which execute BEFORE the swap at 413.
After swap, the previous-iterate data sits in `lam_new`/`grad_new` buffers;
they are overwritten on the next iteration's Wolfe call / `phi_from_u`
before any read. Verified across all three Wolfe paths (accepted,
zoom, fallback).

### 3c — Initial direction at 381: **skipped** (with aliasing proof)

**Context.** Line 381 `dir = grad;` executes only inside
`if (svec.empty())` — iteration 0, or whenever the curvature guard has
rejected every prior update (cold start). Not a hot path.

**Aliasing proof (added per reviewer concern):** `dir` is declared
`std::vector<double> dir(total);` at line 363 as a value (not reference or
pointer). `dir = grad;` at line 381 performs a value-copy; after the
assignment `dir` owns an independent buffer with the same contents as
`grad`. WU-3b's `std::swap(grad, grad_new)` at line 413 swaps the internal
pointers of `grad` and `grad_new`; it does not affect `dir`'s buffer in
any way. No aliasing hazard. Leaving `dir = grad;` as-is is safe.

**Acceptance.**
- `grep -n "std::swap(lam" src/lbfgsb_solver.cpp` returns one hit near line 413.
- `grep -n "^        lam = lam_new;" src/lbfgsb_solver.cpp` returns 0.
- Tests pass. Benchmark no worse than master.

---

## WU-4: Hoist `alpha` scratch via `thread_local`

**Target.** `src/lbfgsb_solver.cpp` line 145 inside `lbfgs_direction`.

**Current:**
```cpp
    std::vector<double> alpha(m);
```

**New:**
```cpp
    // Per-thread scratch, grown once per thread. m ≤ st.lbfgs_m (default 10).
    // thread_local std::vector preferred over std::array<double, N>: lbfgs_m
    // is runtime-configurable (leafblower.h:38), so a hard-coded array size
    // would either truncate or require an assert. thread_local also avoids
    // stack growth on recursion (hypothetical future ALM outer loop).
    thread_local std::vector<double> alpha;
    if ((int)alpha.size() < m) alpha.resize(m);
```

**Acceptance.**
- `grep -c "thread_local std::vector<double> alpha" src/lbfgsb_solver.cpp` returns 1.
- `grep -c "std::vector<double> alpha(m)" src/lbfgsb_solver.cpp` returns 0.
- No call-site changes. Tests pass.

---

## WU-5: Relative curvature gate (threshold only, no stall guard)

**Target.** `src/lbfgsb_solver.cpp`, lines 402-403.

**Justification.** Relative curvature is the canonical L-BFGS safeguard:
- Liu & Nocedal (1989), *Math. Prog.* 45:503-528, §3:
  `s · y > ε · ‖s‖ · ‖y‖` with ε ∈ [1e-8, 1e-6].
- Nocedal & Wright (2006), *Numerical Optimization* 2e, §7.2, same.
- Current `sy > 1e-20` is absolute, non-standard, and effectively disabled
  (denormal-adjacent).

**Current:**
```cpp
        constexpr double kCurvMin = 1e-20;
        if (sy > kCurvMin && yy > kCurvMin) {
```

**New:**
```cpp
        // Relative curvature gate (Liu & Nocedal 1989 §3; Nocedal-Wright 2e §7.2).
        // Squared form avoids sqrt: (sy)² > ε²·‖s‖²·‖y‖². On rejection, fall
        // through to steepest-ascent on this iteration (history unchanged).
        constexpr double kCurvRel = 1e-8;
        double s_norm2 = dot(s_new, s_new);
        bool curv_ok = (sy > 0.0) && (sy * sy > kCurvRel * kCurvRel * s_norm2 * yy);
        if (curv_ok) {
```

No stall guard. If the relative gate repeatedly rejects, the solver
degrades to pure steepest ascent on a concave objective — convergence is
slower but still guaranteed (by concavity + Wolfe). If empirical behaviour
shows problematic stalls, file a follow-up issue with evidence.

**Acceptance.**
- `grep -c "kCurvRel" src/lbfgsb_solver.cpp` returns ≥ 1.
- `grep -c "kCurvMin" src/lbfgsb_solver.cpp` returns 0.
- Tests pass. Iteration counts on `test-lbfgsb.R` do not regress by > 20 %
  (record before/after in commit body). If they do, revert and reopen.
- **Mandatory before commit:** file a bd issue titled
  "L-BFGS-B stall risk on ill-conditioned / near-saturated inputs" tracking
  the scenario where `kCurvRel=1e-8` rejects every update and solver
  degrades to steepest-ascent. Commit body must reference that issue ID.

---

## WU-6: Saturation diagnostic via `st.log` (verbose-gated)

**Target.** `src/lbfgsb_solver.cpp`, inside `compute_final_weights_and_error`,
lines 158-188. Emits via `st.log` — same channel as iEPPA's per-iter progress
log. At `verbose=0` the message is suppressed (standard contract); at
`verbose≥1` it augments the existing NOCONV signal with actionable
saturation counts.

**Insertion** after the `max_err` computation loop (around line 184),
before `res.max_error = max_err;`:

```cpp
    if (max_err >= st.tol_abs && !fn.exponential) {
        int n_sat_low = 0, n_sat_high = 0;
        const double kSatTol = 1e-9;
        const double bw = fn.U - fn.L;
        for (int i = 0; i < st.n; i++) {
            if      (st.weights[i] <= fn.L + kSatTol * bw) n_sat_low++;
            else if (st.weights[i] >= fn.U - kSatTol * bw) n_sat_high++;
        }
        double sat_frac = static_cast<double>(n_sat_low + n_sat_high) / st.n;
        if (sat_frac > 0.5) {
            char msg[256];
            std::snprintf(msg, 256,
                "L-BFGS-B: %.1f%% of weights pinned at [L,U] bounds "
                "(low=%d high=%d); targets may be infeasible given "
                "min/max_weight.",
                100.0 * sat_frac, n_sat_low, n_sat_high);
            st.log(msg);
        }
    }
```

**Heuristic limitation (acknowledged in commit body).** The proxy
`st.weights[i] ≤ fn.L + tol` tests `d[i]·F(u[i]) ≤ L + tol` rather than
`F(u[i]) ≤ L + tol`. When input design weights `d[i]` deviate substantially
from 1, saturated observations with `d[i] ≫ 1` are false-negatives
(weight = `d[i]·L` > `L`). Acceptable for a diagnostic; caller contract
(harvest.R, _harvest.py) normalises `Σd = n`, so most `d[i]` cluster near 1.

**Include addition.** Add `#include <cstdio>` at the top of
`src/lbfgsb_solver.cpp` for `std::snprintf` (do not rely on transitive
inclusion via `<R_ext/Print.h>`). `st.log` is an inline member of
`CalibState` in `types.hpp`; available wherever `types.hpp` is transitively
visible. `lbfgsb_solver.cpp` line 2 includes `lbfgsb_solver.hpp`, which
line 2 includes `types.hpp`. ✓

**No stderr contamination.** Tests run with default `verbose=0`; `st.log`
at line 32 of `types.hpp` returns early if `verbose <= 0`. Zero impact on
the 100-test regression suite. `log_fn` path similarly gated.

**Acceptance.**
- `grep -c "pinned at \[L,U\]" src/lbfgsb_solver.cpp` returns 1.
- All 100 tests pass.
- Manual smoke: R session, infeasible problem (target cat A = 0.9, sample
  5 % A, `max_weight=1.2`, `verbose=1`) — warning appears on stderr.
  Recorded in commit body.

---

## WU-7: SIMD hints on dot/axpy (optional, benchmark-gated)

**Target.** `src/lbfgsb_solver.cpp` lines 127-132 (`dot`) and 148, 150, 153
(inner loops of `lbfgs_direction`).

**Fix.**
```cpp
static double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    int n = (int)a.size();
#if LBW_HAS_OMP_SIMD
#pragma omp simd reduction(+:s)
#endif
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}
```
And `#pragma omp simd` (no reduction) above the three loops at 148, 150, 153.

**Gate.** Run `Rscript benchmarks/autumn_nr_benchmark.R` before and after.
Keep WU-7 only if L-BFGS-B path shows ≥ 3 % speedup on n ≥ 100k. Revert
otherwise. Advisory, not required.

**Acceptance.**
- Compiles clean (no new warnings).
- Benchmark delta recorded. Merge only if ≥ 3 %.

---

## Ordering

1. **Commit 1** (docs, parallel-safe): WU-1 + WU-2.
2. **Commit 2** (mechanical): WU-3b (swap at 413-414) + WU-4 (thread_local).
3. **Commit 3** (robustness): WU-5 (curvature gate only, no stall guard).
4. **Commit 4** (diagnostic): WU-6 (st.log verbose-gated).
5. **Commit 5** (optional, benchmark-gated): WU-7. Revert if < 3 %.

Each commit independently revertable; each passes the full 100-test regression.

## Rollback

Per-WU `git revert <sha>`. No schema, ABI, or API changes.
