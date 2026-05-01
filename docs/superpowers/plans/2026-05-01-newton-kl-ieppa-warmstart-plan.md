# Newton-KL IEPPA Warm-Start — Plan (rev 3)

**Epic:** `leafblower-usg8` (Epic-C; follow-up to Epic-B BLOCKED)
**Spec reference:** `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` (a new "IEPPA warm-start" subsection lands in WI-0)
**Predecessors:** epic `leafblower-5k08` (LM rev 2; landed); `leafblower-91u7` (Epic-B target homotopy; BLOCKED — see `docs/investigations/2026-05-01-newton-kl-homotopy-result.md`)
**Tasks (sequential):** WI-0 (spec) → WI-0b (basin kill-switch) → WI-1 (ieppa-inner extraction) → WI-2 (warm-start wiring) → WI-1c (r_bridge SEXP surfacing) → WI-3 (T6 + T7 + T8 + T9 + T9b + T10 tests) → WI-4 (verify) → WI-5a (bench) → WI-5b (verdict + investigation doc)
**Date:** 2026-05-01

---

## Rev 3 Changelog (plan-review-gate iteration 1 — REVISION NEEDED: Feasibility FAIL ×3, Scope FAIL ×2)

Applied after plan-review-gate rev 2 returned REVISION NEEDED. Fixes are lettered P–U.

| Fix | Section(s) affected | Summary |
|-----|---------------------|---------|
| P | §Forbidden | Remove `src/r_bridge.cpp` from the untouched-files list (narrowly: SEXP-pack of two new scalar fields only — `n_warmstart_iters_used` and `warmstart_max_err_at_handoff`). Does NOT surface `lam`. Resolves Feasibility FAIL: WI-3 tests read `res$n_warmstart_iters_used` and `res$warmstart_max_err_at_handoff` from R, which requires r_bridge.cpp to be touched. |
| Q | §Plan Steps | Insert new ticket WI-1c (r_bridge SEXP surfacing) between WI-2 and WI-3. Touches `src/r_bridge.cpp` only. Adds two scalar entries to `pack_solver_result` on the newton_kl path: `n_warmstart_iters_used` (integer) and `warmstart_max_err_at_handoff` (double). DOES NOT surface `lam`. Atomic ≤30-line change. WI-3 now depends on WI-1c, not WI-2. |
| R | §Plan Steps | Insert new ticket WI-0b (basin-overlap kill-switch) between WI-0 and WI-1. Pure analysis, no new code. One-shot R script `benchmarks/basin_overlap_kill_switch.R` comparing converged IEPPA weights vs Newton weights on stepstone via `max(|log(w_ieppa / w_newton)|)`. PROCEED iff `max(|log_ratio|) >= 1e-3`. Supersedes Fix M (WI-4 basin-overlap experiment), which is removed from WI-4. |
| S | §Cost (new §Epic-C Value Justification) | Add explicit value statement: IEPPA+SRAA standalone achieves 1.13e-4 on stepstone — MISSES the 1e-4 gate by 13%. Newton polish closes that gap. On kk1204 severe-skew, master Newton-KL DIVERGES; IEPPA-warm-started Newton-KL is the only path that converges to <1e-4 within 5s wall. PARTIAL (wall ≥3s) is strict improvement over master's diverging behavior. |
| T | §WI-2 deliverables, §WI-1 | Add `#include "ieppa.hpp"` at top of `newton_calib.cpp` to WI-2 required touches (and note it must exist in WI-1 if the shim is defined in newton_calib.cpp). The static helper `newton_kl_ieppa_warmstart` calls `ieppa_solve` from this include. |
| U | §Plan Steps sequence | New 9-ticket sequence: WI-0 → WI-0b → WI-1 → WI-2 → WI-1c → WI-3 → WI-4 → WI-5a → WI-5b. Header line updated. WI-1 now depends on WI-0b; WI-3 now depends on WI-1c. |

---

## Rev 2 Changelog (plan-review-gate NEEDS_REVISION — 3 reviewers)

Applied after plan-review-gate rev 1 returned NEEDS_REVISION. Fixes are lettered A–O.

| Fix | Section(s) affected | Summary |
|-----|---------------------|---------|
| A | §Mechanism conversion, WI-1 ticket | Explicitly distinguish `cat_offset` (IEPPA, size sum(cat_counts)+K, NA included) vs `lam_off` (Newton, size sum(cat_counts-1), no NA). Conversion reads `lf[cat_offset_ieppa[k] + j]` for j ∈ [0, cat_counts[k]-1]; NEVER accesses NA slot at `cat_counts[k]`. |
| B | §WI-1 / Q1 | Capture `lf_best` not raw `lf`: mirror `ieppa.cpp:386` `W_best`. lf_capture writes ONCE at iter where `best_iter_val` is recorded, not at every exit path with current lf. |
| C | §WI-1 / Q1 | Exit-path enumeration vs RAII: plan mandates RAII guard struct (simpler, less error-prone). Documents the choice. |
| D | §WI-3 T7 | Monotonicity assertion: `expect_lte(res$max_error, res$warmstart_max_err_at_handoff*1.01 + 1e-10)`. Requires `warmstart_max_err_at_handoff` field added to WI-2 deliverables (already present; assertion is new). |
| E | §WI-5a | Bench script must include `harvest(method="ieppa", accelerate=TRUE)` row on same kk1204 fixture for direct comparison baseline. |
| F | §WI-3 | New test T9: synthetic fixture with `max_weight=1` + skewed targets triggers IEPPA infeasibility. Assert `status==0 AND n_warmstart_iters_used==0` (cold-start fallback fired). |
| G | §WI-1 | Empirical K_warm=8 validation: one-shot R script after shim wired. Run IEPPA on stepstone for K_warm ∈ {1,2,4,8,16,32}, log max_err per K_warm. Capture in WI-1 commit body. ABORT if K_warm=8 doesn't reach max_err ≤ 1e-2. |
| H | §WI-3 | New test T9b: K=20 moderate-skew (targets [0.4,0.25,0.15,0.12,0.08]) — assert `system.time(harvest(...)) < 1.5 × baseline_pre_warmstart_seconds`. |
| I | §WI-3 | New test T10: stepstone K=9 with K_warm=0 (via env-var override, Fix L) — assert `max_err ≤ 2.8e-4` (TH-1a baseline preserved). |
| J | §WI-3 T6, §WI-2 | T6: expose test-only C entry point `_newton_kl_warmstart_diagnostic` returning `(lf, λ, max(|u_newton - u_ieppa - C|))`. Asserts conversion math directly. Entry point added in WI-2. |
| K | §WI-2 ticket | Clarify iter accumulation: `res.base.iterations` = Newton-only count; warm-start iters live only in `n_warmstart_iters_used`. |
| L | §WI-2, Q3 | Env var `LBW_NEWTON_KL_WARMSTART_ITERS` (integer, optional). Overrides hardcoded K_warm=8. Used by T6/T10. |
| M | §WI-4 | Basin-overlap experiment in WI-4: one-shot R script logs `max(|λ_warm_final - λ_cold_final|)` on stepstone. If < 1e-3, warm-start is functionally a no-op → escalate WI-5b BLOCKED. |
| N | §Forbidden | Add rationale: "AUTO routing unchanged — warm-start is internal to newton_kl; AUTO selection independent of inner-solver mechanics." |
| O | §Q2 / §Q6 | Clarify `accelerate=FALSE` handling: warm-start always uses SRAA regardless of caller's `accelerate`; user `accelerate` controls IEPPA main loop only when `method='ieppa'`. |

---

## Mechanism Header

- **Mechanism:** Run IEPPA's coordinate-descent inner sweep for a small fixed budget (default `K_warm = 8` outer sweeps) on the *original* targets `T`, convert IEPPA's log-factor vector `lf` to Newton's reference-eliminated `λ` parameterisation, then enter the existing LM-Newton inner with that warm `λ`. Newton polishes from inside the convergence basin.
- **Forbidden:** Target homotopy of any kind (Epic-B BLOCKED). Touching the IEPPA coordinate-descent loop body. Loosening any test gate. Disabling SRAA on the IEPPA warm-start (default ON; rationale in §Mechanism). Touching `harvest()` public signature. Per-method internal homotopy-bounds reuse.
- **Audit:** Two new `NewtonCalibResult` fields — `n_warmstart_iters_used` (count of IEPPA outer sweeps actually run, including any early-exit-on-converged) and `warmstart_max_err_at_handoff` (the IEPPA-side `max_err` measured at handoff). T2 and T8 read both via testthat snapshots; T7 asserts `n_warmstart_iters_used >= 1 && status == 0` (loose, not brittle).

---

## Context

### What landed before this epic

- **Epic 5k08 (LM rev 2; `5355f69`):** Newton-KL solver with LM damping (gain-ratio adaptive `lm_mu`, scale-invariant Hessian damping with additive floor, Armijo line search, best-iterate restoration). Converges on stepstone K=9 to a *natural floor* `gap ≈ 1.4e-4` → `max_err ≈ 2.8e-4` (vs `<1e-4` gate). The floor is the rank-deficient-Hessian drift wall: Newton steps past optimum into unbounded-dual regions where best-iterate rolls back.
- **Epic-B (`91u7`; `5355f69` TH-1a refactor only):** target homotopy regressed T2 from 2.8e-4 to 2.29e-3 (8× worse). Mechanism failure: per-margin-uniform `T_eps` for the stepstone overlap structure (margins `rk_age10_gender`, `rk_gender_time`, `rk_i_loc_time_gender`) lands `λ*(T_eps→0+)` in a *structurally distinct basin* from `λ*(T_0)`. Warm-starting Newton from there is worse than starting from `λ=0`. Implementation discarded; master at TH-1a state.
- **TH-1a refactor (kept):** `src/newton_calib.cpp:139` defines `auto run_newton_inner = [&](const std::vector<double>& T_eps, int max_iter_inner) -> bool` capturing all outer-scope state by reference (`lam`, `lm_mu`, `Z_curr`, `u_max_curr`, `res`, `H`, `G`, `delta`, etc.). Single top-level call at line 371: `run_newton_inner(T, max_iter);`. `NewtonCalibResult` carries `n_homotopy_levels_used = 0` (stale) plus `lm_mu_final`, `dual_gap`, `step_norm`, `line_alpha`. This refactor is the affordance Epic-C consumes — no further changes needed in the Newton inner.

### Why IEPPA warm-start is the right next step

IEPPA (faithful coordinate descent in log-factor space; `src/ieppa.cpp:127` `ieppa_solve`) is *empirically robust on overlapping margins*. On the same stepstone fixture where target homotopy regressed by 8×, ieppa+sraa already converges to `max_err = 1.13e-4` (better than the LM baseline 2.8e-4 and approaching the `<1e-4` gate by itself). On kk1204 severe-skew (n=1M, K=20, nj=5, max_weight=3, OMP=1), ieppa+sraa is the published top performer at 3.7s / max_err 2.4e-14 (spec table). IEPPA's per-margin-per-category coordinate updates cannot overshoot across margin boundaries — a sweep is a sequence of one-dimensional moves, each on the linearly-separable problem `α* = log(T_kj · Z_total / S_kj)`. The dual basin is reached by construction, not by Hessian luck.

The *combined* picture: IEPPA brings us inside Newton's quadratic basin in O(K) outer sweeps (linear convergence rate); LM-Newton then polishes from inside the basin at the rate that originally motivated the Newton-KL solver. This is the textbook coordinate-descent → second-order combination for overlapping-margin KL minimisation.

### Why this is mathematically a no-op handoff (correctness argument)

- IEPPA carries `lf[cat_offset[k] + j]` for all `j = 0..cat_counts[k]-1` (no reference elimination). Each obs has `u_i_ieppa = Σ_k lf[cat_offset[k] + j_k(i)]` (positive `j` only; `j<0` is the NA bucket and is skipped — same convention as the Newton inner at `src/newton_calib.cpp:92-93`).
- Newton carries `λ[lam_off[k] + j - 1]` for `j = 1..cat_counts[k]-1` (reference category 0 fixed at `λ_{k,0} = 0`). Each obs has `u_i_newton = Σ_{k: j_k(i)>0} λ[lam_off[k] + j_k(i) - 1]`.
- **Conversion:** for each margin `k` and each `j ≥ 1`, set
  ```
  λ[lam_off[k] + j - 1] = lf[cat_offset[k] + j] - lf[cat_offset[k] + 0]
  ```
- Then `u_i_newton(λ) = u_i_ieppa(lf) - Σ_k lf[cat_offset[k] + 0] = u_i_ieppa(lf) - C` where `C` is a `λ`-independent constant.
- Newton's stable LSE form (`src/newton_calib.cpp:106-118`) computes `Z_stable = Σ d_i exp(u_i - u_max)`, `g(λ) = u_max + log(Z_stable) - T·λ`. Both `Z_stable` and `u_max` shift by the same constant under `u → u - C`; their cancellation in `log(Z_stable) - max` is exact. The recovered weights `w_i ∝ d_i · exp(u_i)` are scale-equivalent under any constant shift of `u`. Hence the converted `λ` produces *bit-identical primal weights* to IEPPA's `lf`.
- **Important caveat — NA buckets (Fix A).** IEPPA stores `lf` for `j = cat_counts[k]` (the NA bucket at `cat_offset_ieppa[k] + cat_counts[k]`). Note: `cat_offset_ieppa[k+1] = cat_offset_ieppa[k] + cat_counts[k] + 1` (the `+1` is this NA bucket). Newton has no NA-bucket dual; `lam_off[k+1] = lam_off[k] + cat_counts[k] - 1` (reference-eliminated, no NA slot). Conversion loop reads `lf[cat_offset_ieppa[k] + j]` for `j ∈ [0, cat_counts[k]-1]` strictly — the loop bound `j < cat_counts[k]` ensures the NA slot at `j = cat_counts[k]` is NEVER accessed. The NA `lf` slot is silently dropped; this is correct because no obs with `j_k(i) = NA` enters Newton's dual sums. Confirmed against `src/ieppa.cpp:649` (`s += lf[cat_offset[m] + gm]` only fires when `gm >= 0`) and `src/newton_calib.cpp:91-93`.

### What the IEPPA solver gives us, structurally

- `ieppa_solve(CalibState&)` (`src/ieppa.cpp:127`, 1924 lines) is a single monolithic function. Its inner data — `lf` (`src/ieppa.cpp:202`, size `cat_offset[st.K]`), `cell_lf` (incremental), `X_cur`, SRAA-m state, SOR-omega state, bounds machinery — lives on its stack frame. The faithful inner loop body is the per-`(k, j)` sweep starting around `src/ieppa.cpp:637` (the `for (k) { for (j) { compute α* and update lf } }` block).
- IEPPA already exposes a *sweep budget* via `st.outer_max_iter`. We will **not** carve a new inner-loop helper out of `ieppa.cpp`. Instead, WI-1 introduces a thin orchestration shim that calls `ieppa_solve` against a duplicated `CalibState` configured with the warm-start budget, then extracts `lf` via a single new accessor. Rationale: extracting the per-sweep loop body cleanly would require touching ~50 captures and risk behaviour drift on every existing IEPPA caller. The shim is an additive change.

---

## Mechanism

### Conversion lf → λ (exact)

**Index space distinction (Fix A):**
- **IEPPA `cat_offset`**: total size `sum(cat_counts) + K` (one NA bucket per margin at index `cat_counts[k]`). `cat_offset[k+1] = cat_offset[k] + cat_counts[k] + 1`. Valid category indices for margin k: `j ∈ [0, cat_counts[k]-1]`. Index `j == cat_counts[k]` is the NA bucket — NEVER read in conversion.
- **Newton `lam_off`**: total size `sum(cat_counts - 1)` (reference category j=0 eliminated per margin). `lam_off[k+1] = lam_off[k] + (cat_counts[k]-1)`. Free dual indices: `j ∈ [1, cat_counts[k]-1]`.

```cpp
// Inputs:
//   lf            : IEPPA lf_best (captured at best_iter_val, not at exit) [Fix B]
//                   size cat_offset_ieppa[K], stride per margin = cat_counts[k]+1 (NA bucket)
//                   lf[cat_offset_ieppa[k] + j] for j=0..cat_counts[k]-1 ONLY
//                   j == cat_counts[k] is the NA bucket — NEVER access [Fix A]
//   cat_offset_ieppa : ieppa-side cat_offset (size K+1;
//                      cat_offset_ieppa[k+1] = cat_offset_ieppa[k] + cat_counts[k] + 1)
//   lam_off       : newton-side lam_off (size K+1;
//                      lam_off[k+1] = lam_off[k] + cat_counts[k]-1)
// Output:
//   lam           : size lam_off[K], lam[lam_off[k] + j - 1] for j=1..cat_counts[k]-1
//
// Math: λ[k,j-1] = lf[k,j] - lf[k,0]. Constant column shift cancels in LSE.
for (int k = 0; k < K; ++k) {
    const double lf0 = lf[cat_offset_ieppa[k] + 0];         // reference category
    for (int j = 1; j < cat_counts[k]; ++j)                 // STOP before cat_counts[k] (NA)
        lam[lam_off[k] + j - 1] = lf[cat_offset_ieppa[k] + j] - lf0;
}
// INVARIANT: no access to lf[cat_offset_ieppa[k] + cat_counts[k]] (NA bucket)
```

Verification before commit: a unit test asserts that for a tiny fixture (n=8, K=2, nj=3), running 1 IEPPA sweep + conversion produces `λ` whose `compute_u(λ, i)` differs from IEPPA's `cell_lf[c_of(i)]` by a single per-iter constant — i.e., `var(u_newton - u_ieppa) < 1e-12`. This is the structural correctness gate (T6, see §Plan Steps).

### IEPPA-iter-budget reasoning

IEPPA on stepstone-class problems converges at a *linear rate* with contraction factor empirically `≈ 0.5–0.7` per outer sweep (see `docs/investigations/2026-04-23-kk1204-convergence.md` figures). Starting from `lf = 0` (uniform), after `K_warm` sweeps the residual is `≈ ρ^K_warm × max_err_initial`. Stepstone initial `max_err` at `lf=0` is `O(1)`; for `ρ = 0.6` and target handoff residual `1e-2` (well inside Newton's basin given the LM solver's basin radius observed empirically as `gap < 1e-1`), `K_warm = log(1e-2) / log(0.6) ≈ 9`. Round to **8** as default; it provides a 5-10× margin within Newton's basin while costing one sub-second IEPPA pass on K=9 / one ~3s pass on K=20.

**Decision:** `K_warm = 8` is hard-coded at top of `newton_calibrate`, with the comment citing the contraction-rate calculation. Not exposed via `CalibState` in this epic — keep the API surface flat. The constant is `constexpr int kNewtonKLWarmstartIters = 8;`. Tunability is a follow-up if benchmarks show wall-time savings from `K_warm = 4` (early handoff) or quality gains from `K_warm = 16` (deeper warm-up).

**Adaptive early-exit (kept simple):** if IEPPA's internal convergence (`st.tol_abs`) trips inside the budget, IEPPA returns early. Newton sees `n_warmstart_iters_used < kNewtonKLWarmstartIters` and proceeds. No special-cased "skip Newton on early IEPPA convergence" path — see Q4 below.

### Q1: IEPPA refactor — minimal extraction shim

Instead of refactoring `ieppa_solve` into `ieppa_inner(CalibState&, int max_iters) -> lf`, introduce a thin **orchestration shim** in a new translation unit `src/newton_calib_warmstart.cpp` (or inline at top of `newton_calib.cpp` — implementer chooses):

```cpp
// Run a warm-start IEPPA pass and return the converted λ.
// Configures a budget-limited IEPPA call against a copied CalibState
// (so we don't mutate caller's max_iter / accelerate / convergence_cfg).
// Returns:
//   lam_out  : size lam_off[K], converted from IEPPA's final lf
//   n_iters  : number of IEPPA outer sweeps actually run
//   max_err  : IEPPA-side max_err at handoff (for diagnostics)
struct WarmstartResult {
    std::vector<double> lam;
    int n_iters;
    double max_err;
    int status;  // RK_OK or first IEPPA error code
};
WarmstartResult newton_kl_ieppa_warmstart(const CalibState& st, int K_warm);
```

The shim:
1. Copies `st` into a local `CalibState st_warm` (struct-copy is cheap; pointer fields are caller-owned and shared safely).
2. Sets `st_warm.outer_max_iter = K_warm` and `st_warm.accelerate = true` (SRAA on; see Q2).
3. Calls `ieppa_solve(st_warm, &lf_out)` and captures the returned `IEPPAResult.base.max_error`, `iterations`, and `status`.
4. **lf-capture implementation — RAII guard, not exit-path enumeration (Fix C):** Introduce a RAII guard struct inside `ieppa_solve` that writes `*lf_capture = lf_best` on destruction. This fires at ALL exits (normal return, exception, early-exit) without enumerating every `return` statement. The guard captures `lf_capture` and `lf_best` by pointer/reference. **`lf_best`, not `lf` (Fix B):** mirror `ieppa.cpp:386` `W_best` — IEPPA already maintains a best-iterate `lf_best` (updated at each iter where the objective improves). The capture writes `lf_best` once, at the iter where `best_iter_val` is recorded, NOT the raw `lf` at function exit (which may be past-best due to SRAA overstep). This ensures the converted `λ` reflects IEPPA's quality-checked best point.

   ```cpp
   // RAII lf-capture guard (Fix B + Fix C)
   struct LfCaptureGuard {
       std::vector<double>* out;
       const std::vector<double>& lf_best;
       ~LfCaptureGuard() { if (out) *out = lf_best; }
   } lf_guard{lf_capture, lf_best};
   // No exit-path enumeration needed — destructor fires at ALL exits.
   ```

5. Performs the conversion above into `lam_out` using `cat_offset_ieppa` and stopping at `cat_counts[k]` (Fix A).

**K_warm empirical validation (Fix G):** After WI-1 shim is wired, before closing the ticket, the implementer MUST run a one-shot K_warm sweep on stepstone:

```r
# Run after shim is wired (can be a scratch script, not committed)
for (K_warm in c(1, 2, 4, 8, 16, 32)) {
  Sys.setenv(LBW_NEWTON_KL_WARMSTART_ITERS = K_warm)
  res <- harvest(method = "newton_kl", ...)   # stepstone fixture
  cat(sprintf("K_warm=%d  max_err=%.2e  n_iters_used=%d\n",
              K_warm, res$warmstart_max_err_at_handoff, res$n_warmstart_iters_used))
}
```

**ABORT criterion (Fix G):** If K_warm=8 does NOT produce `warmstart_max_err_at_handoff ≤ 1e-2` on stepstone, halt and escalate before proceeding to WI-2. Capture the per-K_warm table in the WI-1 commit body.

**This is the minimal refactor.** It touches one new function declaration in `src/ieppa.hpp`, one RAII guard + `lf_best` capture in `src/ieppa.cpp`, and adds the shim. Existing callers (`harvest(method="ieppa")`, AUTO routing, `r_bridge.cpp`) are bit-identical. WI-1 ticket implements exactly this and only this.

### Q2: SRAA on or off for warm-start? — **ON**

**`accelerate=FALSE` handling (Fix O):** The warm-start shim always sets `st_warm.accelerate = true` (SRAA on), regardless of the *caller's* `accelerate` parameter. The user's `accelerate` flag controls the IEPPA main loop only when `method='ieppa'`. When `method='newton_kl'`, SRAA inside the warm-start shim is an internal optimisation detail invisible to the public API.

`accelerate = true` (SRAA-m) on the IEPPA warm-start. Justification:
- SRAA-m converges 3-5× faster on stepstone (per `docs/superpowers/specs/2026-04-30-ieppa-sraa-acceleration.md`). The whole point of an 8-iter budget is to land deep in Newton's basin cheaply. Plain IEPPA at 8 iters gives roughly a 10× residual contraction; SRAA-m gives roughly 30-100×. Newton sees a tighter handoff for the same wall cost.
- SRAA-m's failure mode is super-step rejection (falls back to plain IEPPA on that sweep), not divergence. There is no scenario where SRAA-on hurts.
- The IEPPA+SRAA path is the published top performer on kk1204 (3.7s / max_err 2.4e-14). Reusing the same code path is the lowest-risk choice.
- SRAA off (plain IEPPA) is the *fallback* in WI-2's status-1 recovery path: if SRAA produces NaN or status≠0 inside the warm-start, the shim re-runs once with `accelerate = false` before bailing to "skip warm-start, run Newton from λ=0".

### Q3: Warm-start budget — **`K_warm = 8`, fixed**

See the contraction-rate derivation above. Not adaptive on gap because (a) the gap-to-Newton-basin-radius mapping isn't a clean threshold (basin depends on local `H` curvature, which Newton itself measures), (b) every dynamic-budget scheme adds branching that needs its own tests, and (c) 8 iters of IEPPA+SRAA on K=9 is `<0.5s` (per spec table extrapolation; ieppa+sraa needs 10 iters for full convergence at 3.7s on K=20 with n=1M; per-iter cost there is ~0.4s; on K=9 n=200K it is ~0.05s). The budget is dwarfed by Newton's own per-iter cost on these sizes.

**Env-var override for testing (Fix L):** The env var `LBW_NEWTON_KL_WARMSTART_ITERS` (integer, optional) overrides `kNewtonKLWarmstartIters` at runtime. Read via `std::getenv` at the top of `newton_calibrate`. If set and parseable as a non-negative integer, use that value instead of 8. If `0`, warm-start is skipped entirely (K_warm=0 path). Used by T6, T10, and the empirical K_warm sweep (Fix G). Not exposed in user docs — test/debug knob only.

### Q4: What if IEPPA already converges (max_err < 1e-4) inside K_warm?

**Decision: always run Newton afterwards (option (a)).** Justification:
1. Newton's polishing cost from a converged warm start is ~1-3 iters (basin entry → quadratic floor in 1-2 steps). Per-iter cost on stepstone is ~0.05s; total Newton tail is ~0.1s. Negligible.
2. Newton from a converged IEPPA `λ` is a *correctness check by construction*: if Newton's first-iter gradient is `<tol_abs`, it exits at iter 0 with `status = RK_OK` and the same primal weights IEPPA produced. The user gets a quality-checked Newton-KL result with no measurable cost. If somehow Newton's gradient is *not* tiny (e.g., target rounding mismatch between ieppa's per-cat residual and newton's per-(k,j-1) gradient), Newton polishes; the user still wins.
3. Option (b) "adaptive handoff: skip Newton if IEPPA converged" branches the API surface for negligible savings and creates two distinct user-facing convergence paths to test. Avoid.

### Q5: Failure modes

| Failure | Detection | Recovery |
|---|---|---|
| IEPPA SRAA produces NaN / `status = RK_ERR_*` inside warm-start | `IEPPAResult.base.status != RK_OK` from shim | Retry once with `accelerate = false` (plain IEPPA, same K_warm). If that also fails, **skip warm-start entirely**: log a one-line warning, set `n_warmstart_iters_used = 0`, run Newton from `λ = 0`. Newton's existing best-iterate-fallback handles the rest. |
| IEPPA bound infeasibility (max_weight=3 + skew, no feasible primal exists) | `IEPPAResult.base.status` returns the bounds error code | Same as above: skip warm-start, run Newton from `λ=0`. Newton itself doesn't enforce bounds (it's the dual optimiser only); the bounds-infeasibility surface is the IEPPA cell-mode invariant, not Newton's concern. The user sees a Newton result with the unmodified `max_err` characterisation. |
| Conversion produces non-finite `λ` (lf overflow during warm-start before SRAA quenches) | `std::isfinite(lam[a])` check on conversion output | Skip warm-start, log warning, Newton from `λ=0`. |
| Newton from warm-started `λ` has `status != RK_OK` after `max_iter` | Caller-side: existing post-Newton path handles this exactly as today. | No change — Newton's own NOCONV/best-iterate logic handles it. |

In all four failure modes, the behaviour reduces to "Newton from `λ = 0`", which is exactly today's master behaviour. The warm-start is *strictly additive*: it can only improve quality, never regress it (modulo the warm-start's own wall cost).

### Q6: Public API impact — none

`harvest(method="newton_kl", ...)` signature is unchanged. AUTO routing in `src/c_api.cpp:182` is unchanged. The only user-visible changes are:
1. Two new fields on `NewtonCalibResult` (`n_warmstart_iters_used`, `warmstart_max_err_at_handoff`) — additive.
2. A NEWS.md bullet (mandatory in same commit as WI-2 wiring per discipline §4).

### Q7: Field naming — rename `n_homotopy_levels_used` → `n_warmstart_iters_used`

**Decision: rename in place (option (a)).** Justification:
- The field was added by Epic-B's TH-1a (`5355f69`); Epic-B is BLOCKED and its semantics did not ship to users. The R-side surfacing (LM-1c) was *deliberately deferred* per the homotopy investigation's "Outstanding scope" section: *"`lm_mu_final` and `n_homotopy_levels_used` R-side surfacing deferred until a method actually needs it"*. So no R-side consumer exists; no `r_bridge` SEXP-pack code currently reads the field.
- Grep for `n_homotopy_levels_used` in the worktree (verified during planning) finds it only in `src/newton_calib.cpp` (zero writes — the homotopy code that wrote it was discarded), `src/newton_calib.hpp` (the field declaration), the spec doc, the homotopy plan, and the homotopy investigation. No live reader exists outside the source files this epic touches anyway.
- Option (b) "deprecate-in-place + add new field" leaves a permanently-zero stale field on a publicly-allocated struct. Worse hygiene; same risk profile.

The rename is a single-commit edit on `src/newton_calib.hpp` (rename declaration), `src/newton_calib.cpp` (rename in `NewtonCalibResult res;` initialisation paths). Done in WI-2 alongside the warm-start wiring (algorithmic-with-code colocation per discipline §4).

Add `double warmstart_max_err_at_handoff = 0.0;` as a second new field. Rationale: it gives reviewers and users a one-glance signal of "did the warm-start actually pre-converge or not?" — useful for debugging without re-running with `verbose`.

---

## Forbidden

Carried forward from rev 2 homotopy plan where applicable, plus IEPPA-specific entries. Italicised items are *not* applicable to this epic and are listed only to pre-empt review re-raising.

- **No target homotopy of any form.** Epic-B BLOCKED on this; do not reintroduce.
- **No touching the IEPPA coordinate-descent loop body** (`src/ieppa.cpp:637-690` and surrounding sweep machinery). The shim only changes IEPPA's *external* affordances (one new lf-capture parameter on `ieppa_solve`, defaulted to `nullptr` so existing callers are bit-identical).
- **No SRAA disabled by default** on the warm-start. SRAA off is the recovery path only.
- **No loosening of any test gate.** T2 stays at `<1e-4`. T8 stays at `<1e-3`. T7 (new) lands at `<1e-4`. NO relaxation to mask convergence shortfall.
- **No exposing `K_warm` via `CalibState` or `harvest()`** in this epic. Hard-coded constant. Tunability is a follow-up.
- **No "skip Newton if IEPPA converged" branch.** Always run Newton after warm-start.
- **No homotopy-bounds reuse.** `ieppa.homotopy_levels` is a different concept (bounds homotopy on `max_weight`); leave it alone. The warm-start runs at the user's `max_weight` directly.
- **No public API changes.** `harvest()` signature, AUTO routing in `src/c_api.cpp:182` and `src/r_bridge.cpp:425` — all unchanged. **Rationale (Fix N):** AUTO routing is unchanged because warm-start is internal to `newton_kl`; AUTO's method-selection logic is independent of the inner-solver mechanics. AUTO dispatches to `newton_kl` or `ieppa` based on fixture characteristics (K, skew, cell density); the warm-start does not change which method is selected.
- **No NEWS.md update deferred to a later ticket.** Algorithmic change colocates with code in WI-2 per discipline §4.
- **No skipping pre-commit hooks.**
- **No T2 amendment.** WI-3 explicitly forbids touching T2; only ADDs T7 + T8 (T8 is the K=20 severe-skew unit test mirrored from Epic-B's design).
- **No bundling tickets.** WI-0, WI-1, WI-2, WI-3, WI-4, WI-5a, WI-5b are each one bead. WI-2 is "warm-start wiring + NEWS.md + field rename" — atomic only because the rename is meaningless without the wiring.
- **No touching files outside:** `src/newton_calib.cpp`, `src/newton_calib.hpp`, `src/ieppa.cpp` (one new lf-capture parameter on `ieppa_solve`), `src/ieppa.hpp` (one new function signature), `src/r_bridge.cpp` **(Fix P — narrowly in scope for WI-1c only: SEXP-pack of `n_warmstart_iters_used` (integer) and `warmstart_max_err_at_handoff` (double) on the newton_kl result path. Does NOT surface `lam`. ≤30 lines. No other changes to r_bridge.cpp.)**, `tests/testthat/test-newton-kl.R`, `benchmarks/newton_kl_bench.R`, `benchmarks/results/`, `NEWS.md`, `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` (WI-0 only), `benchmarks/basin_overlap_kill_switch.R` (WI-0b only). NO touches to `src/c_api.cpp`, `R/*.R`.
- *Not applicable: `T_eps` schedules, `lm_mu` per-ε reset semantics, `T_eps` feasibility guard, ε-schedule tail-jump risk, `n_homotopy_levels_used == N_EPS` brittle assertion.* All Epic-B concerns; this epic does no homotopy.
- *Not applicable: "rank deficiency in `H` is fixed by warm start".* It isn't. Warm start lands `λ` inside the basin; rank deficiency can still strike if the LM solver oversteps from the warm-start. The LM-Newton's own best-iterate-restoration handles that, exactly as today.

---

## Plan Steps

Each ticket is independently revertible. The sequence is strict (9 tickets): WI-0 → WI-0b → WI-1 → WI-2 → WI-1c → WI-3 → WI-4 → WI-5a → WI-5b. No parallel paths in this epic.

---

### WI-0 — Spec amendment

- **File:** `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md`
- **Objective:** Add an "IEPPA warm-start" subsection after the "Target homotopy" subsection (which Epic-B's TH-0 already landed). Document: rationale (rank deficiency robust to coordinate descent, basin entry by construction), `K_warm = 8` derivation, lf → λ conversion math, SRAA default, failure-fallback behaviour, the two new `NewtonCalibResult` fields, AUTO routing impact (none).
- **Constraints:** Spec amendment only. No code touched. No NEWS.md. Conventional commit `docs(spec): Newton-KL — add IEPPA warm-start section`.
- **Body sketch (skeletal — controller expands at ticket creation):**
  - Subsection header.
  - Math block: lf → λ conversion (copied verbatim from §Mechanism).
  - Justification block: linear-rate-contraction argument for K_warm = 8.
  - Field documentation: `n_warmstart_iters_used`, `warmstart_max_err_at_handoff`, deprecation of `n_homotopy_levels_used`.
  - Cross-reference to the homotopy investigation as the "why not target homotopy" explanation.

---

### WI-0b — Basin-overlap kill-switch (Fix R)

- **Files:** `benchmarks/basin_overlap_kill_switch.R` (new, committed). `benchmarks/results/basin_overlap_killswitch.csv` (output).
- **Blocked by:** WI-0 (spec must exist before referencing it in the script header).
- **Objective:** Pure analysis ticket. NO new source code. Empirically determine whether IEPPA and Newton-KL converge to the same primal basin on stepstone, BEFORE any warm-start code is written. This is the kill-switch: if they already converge to the same basin, warm-start is a functional no-op and the epic should be aborted.
- **Method:** Uses CURRENT MASTER (TH-1a state, no warm-start code). One-shot R script:
  1. `harvest(method="ieppa", accelerate=TRUE)` on stepstone → `converged_ieppa_weights`.
  2. `harvest(method="newton_kl")` on stepstone (cold-start, current master path, drift to gap=1.4e-4 then best-iterate restoration) → `newton_basin_weights`.
  3. Compute per-obs log-weight-ratio: `log_ratio[i] = log(w_ieppa[i] / w_newton[i])` for all observations.
  4. Compute `max_log_ratio = max(abs(log_ratio))`.
  5. Save full per-obs CSV + summary row to `benchmarks/results/basin_overlap_killswitch.csv`.
- **Decision criterion:**
  - **PROCEED** iff `max_log_ratio >= 1e-3` (warm-start would land Newton in a basin meaningfully different from cold start → the warm-start can change behavior).
  - **ABORT** iff `max_log_ratio < 1e-3` (both methods already converge to the same primal basin → warm-start is a no-op → halt epic, file BLOCKED on leafblower-usg8 with evidence).
- **Supersedes:** Fix M (basin-overlap experiment in WI-4). Fix M is removed from WI-4 (see below). Running this earlier means we abort before writing any code if the fundamental hypothesis is wrong.
- **Script structure:**
  ```r
  # benchmarks/basin_overlap_kill_switch.R
  # Basin-overlap kill-switch for Epic-C (Newton-KL IEPPA warm-start)
  # Uses CURRENT MASTER — no warm-start code must be present.
  # Spec: docs/superpowers/plans/2026-05-01-newton-kl-ieppa-warmstart-plan.md §WI-0b

  library(leafblower)
  library(arrow)

  # Load stepstone fixture
  stepstone <- read_parquet("tests/testthat/fixtures/stepstone_small.parquet")
  targets   <- readRDS("tests/testthat/fixtures/stepstone_small_targets.rds")

  # Run IEPPA
  res_ieppa <- harvest(stepstone, targets, method = "ieppa", accelerate = TRUE)
  w_ieppa   <- res_ieppa$weights

  # Run Newton-KL (cold start, current master)
  res_newton <- harvest(stepstone, targets, method = "newton_kl")
  w_newton   <- res_newton$weights

  # Compute log-ratio
  log_ratio <- log(w_ieppa / w_newton)
  max_lr    <- max(abs(log_ratio))

  cat(sprintf("max|log(w_ieppa / w_newton)| = %.4e\n", max_lr))
  cat(sprintf("Decision: %s\n", if (max_lr >= 1e-3) "PROCEED" else "ABORT EPIC"))

  # Save results
  dir.create("benchmarks/results", showWarnings = FALSE)
  write.csv(
    data.frame(obs = seq_along(log_ratio), log_ratio = log_ratio),
    "benchmarks/results/basin_overlap_killswitch.csv", row.names = FALSE
  )
  ```
- **Constraints:** Commit the script (not scratch). `R CMD INSTALL --preclean .` before running. The script itself does NOT modify any source. Commit: `bench(epic-c): basin-overlap kill-switch experiment`.
- **Halt criterion:** If `max_log_ratio < 1e-3`, output ABORT EPIC, mark `leafblower-usg8` BLOCKED, stop. Do NOT proceed to WI-1.

---

### WI-1 — IEPPA lf-capture orchestration shim (no Newton wiring yet)

- **Files:** `src/ieppa.hpp`, `src/ieppa.cpp`, `src/newton_calib.cpp` (new helper function — does NOT call into Newton inner), `src/newton_calib.hpp` (forward declare warmstart helper if exposed at TU level; otherwise file-static).
- **Blocked by:** WI-0b (basin kill-switch must PROCEED before any warm-start code is written).
- **Fix T — `#include "ieppa.hpp"` in newton_calib.cpp:** If the `newton_kl_ieppa_warmstart` shim is defined in `newton_calib.cpp`, it calls `ieppa_solve`. That requires `#include "ieppa.hpp"` at the top of `newton_calib.cpp`. Add this include in this ticket.
- **Objective:** Land the orchestration shim (`newton_kl_ieppa_warmstart`) and the IEPPA lf-capture affordance. **No behavioural change to the Newton solver in this ticket** — `run_newton_inner(T, max_iter)` is still called with `lam` initialised to zero. The shim is dead code at end of WI-1 (called by no one). Build clean. All existing tests pass bit-identically.
- **Constraints:** `ieppa_solve` signature change must default-nullptr the new lf-capture parameter so all existing callers compile and behave bit-identically. Confirm this with grep before commit: *all* call sites of `ieppa_solve` (the `harvest` path through `r_bridge.cpp`, AUTO routing through `c_api.cpp`, any internal callers) take the no-arg default. Build clean (`R CMD INSTALL --preclean .`); zero new warnings.
- **Body sketch:**
  1. Edit `src/ieppa.hpp`: change `IEPPAResult ieppa_solve(CalibState& state);` to `IEPPAResult ieppa_solve(CalibState& state, std::vector<double>* lf_capture = nullptr);`.
  2. Edit `src/ieppa.cpp`: at the function's exit paths (convergence return + recovery returns + early-exit returns), if `lf_capture != nullptr`, write `*lf_capture = lf;` (copy current lf vector).
  3. Edit `src/newton_calib.cpp`: add `static WarmstartResult newton_kl_ieppa_warmstart(const CalibState& st, int K_warm);` definition with body per §Mechanism Q1. Include the lf → λ conversion. Include the SRAA-fallback recovery logic per Q5.
  4. Add a TU-local unit-test affordance: a small inline assertion (gated on `#ifndef NDEBUG`) that verifies `var(u_newton(lam) - u_ieppa(lf)) < 1e-12` on the shim's output, run on a tiny synthetic fixture (n=8, K=2, nj=3) at TU init time. Or — cleaner — drop this in favour of a proper testthat unit test as part of WI-3 (T6). Implementer's choice; the testthat path is recommended.
  5. Conventional commit: `feat(newton_kl,ieppa): add IEPPA lf-capture shim for Newton-KL warm-start (no wiring yet)`.

---

### WI-2 — Warm-start wiring + field rename + NEWS.md + diagnostic entry point

- **Files:** `src/newton_calib.cpp`, `src/newton_calib.hpp`, `NEWS.md`. Single atomic commit.
- **Fix T — `#include "ieppa.hpp"` required:** `src/newton_calib.cpp` must include `#include "ieppa.hpp"` at the top (or wherever system includes are placed) so that `newton_kl_ieppa_warmstart` can call `ieppa_solve`. This include is a required touch in this ticket. List it explicitly in §IV step 2 of the WI-2 ticket.
- **Objective:** Rename `n_homotopy_levels_used` → `n_warmstart_iters_used`; add `warmstart_max_err_at_handoff`; wire the warm-start shim before the Newton inner call at `src/newton_calib.cpp:371`; add env-var override `LBW_NEWTON_KL_WARMSTART_ITERS` (Fix L); add test-only diagnostic entry point (Fix J); clarify iter accumulation (Fix K). Update NEWS.md.
- **Constraints:** Bit-identical Newton inner behaviour for `K_warm = 0`. Build clean.
- **Fix J — Diagnostic entry point:** Add a C-linkage function (gated on `#ifdef LEAFBLOWER_TESTING`) callable from testthat via `.Call`:
  ```cpp
  // Test-only: returns list(lf=lf_best, lam=converted_lam, max_conv_residual=double)
  // where max_conv_residual = max|u_newton(lam) - u_ieppa(lf_best) - C|
  SEXP _newton_kl_warmstart_diagnostic(SEXP r_calib_state);
  ```
  This entry point runs `newton_kl_ieppa_warmstart` on the provided CalibState, computes both `u_newton(λ)` and `u_ieppa(lf_best)` for each obs, computes C = mean(u_ieppa - u_newton), and returns `max(|u_newton + C - u_ieppa|)`. Used by T6 to assert conversion residual < 1e-12.
- **Fix K — Iter accumulation:** `res.base.iterations` counts Newton-only iterations (unchanged from existing semantics). `res.n_warmstart_iters_used` counts IEPPA outer sweeps in the warm-start shim. These are never summed. Do NOT add them together in any field or log message.
- **Body sketch:**
  1. `src/newton_calib.hpp`: rename `int n_homotopy_levels_used = 0;` to `int n_warmstart_iters_used = 0;`. Add `double warmstart_max_err_at_handoff = 0.0;`.
  2. `src/newton_calib.cpp`: above the lambda definition, add `constexpr int kNewtonKLWarmstartIters = 8;`. Comment cites contraction-rate derivation and pins `8` as default.
  3. Replace the single `run_newton_inner(T, max_iter);` at line 371 with:
     ```cpp
     // Warm-start via IEPPA coordinate descent (8 sweeps, SRAA on).
     // Lands λ inside Newton's quadratic basin on overlapping-margin fixtures
     // where bare LM-Newton drifts (e.g., stepstone K=9). Strictly additive:
     // failure modes fall back to λ=0 with no behaviour change. See spec.
     const auto warm = newton_kl_ieppa_warmstart(st, kNewtonKLWarmstartIters);
     if (warm.status == RK_OK && warm.lam.size() == static_cast<size_t>(n_lam)) {
         std::copy(warm.lam.begin(), warm.lam.end(), lam.begin());
         res.n_warmstart_iters_used = warm.n_iters;
         res.warmstart_max_err_at_handoff = warm.max_err;
         // Refresh the Newton inner's local g_curr / Z_curr / u_max_curr
         // so the first iter's convergence check sees the warm λ.
         g_curr = eval_dual(lam, Z_curr, u_max_curr);
     }
     // else: lam is still zero, n_warmstart_iters_used = 0; Newton runs cold.
     run_newton_inner(T, max_iter);
     ```
  4. NEWS.md (same commit): add bullet under `## Newton-KL calibration`:
     > * `method="newton_kl"` now warm-starts via 8 iterations of IEPPA coordinate descent (with SRAA-m acceleration) on the original targets before entering the LM-damped Newton inner. Lands `λ` inside Newton's quadratic basin on overlapping-margin fixtures (e.g., stepstone K=9) where bare LM-Newton drifts past the optimum. Strictly additive: warm-start failure falls back to cold start with no behaviour change. New result fields: `n_warmstart_iters_used`, `warmstart_max_err_at_handoff`. The Epic-B `n_homotopy_levels_used` field is renamed to `n_warmstart_iters_used` (no R-side consumer existed).
  5. Build clean.
  6. Conventional commit: `feat(newton_kl): IEPPA warm-start before LM-Newton inner; field rename`.

---

### WI-1c — r_bridge SEXP surfacing of warm-start fields (Fix Q)

- **Files:** `src/r_bridge.cpp` only. Single atomic commit. ≤30 lines changed.
- **Blocked by:** WI-2 (`NewtonCalibResult` must carry `n_warmstart_iters_used` and `warmstart_max_err_at_handoff` before r_bridge.cpp can pack them).
- **Objective:** Surface the two new scalar fields from `NewtonCalibResult` into the R list returned by `harvest(method="newton_kl")`. This is the prerequisite for WI-3 tests that read `res$n_warmstart_iters_used` and `res$warmstart_max_err_at_handoff`.
- **Scope (strictly bounded):**
  - ADD `n_warmstart_iters_used` as an integer scalar to the SEXP list in `pack_solver_result` (or equivalent newton_kl result-packing block in r_bridge.cpp).
  - ADD `warmstart_max_err_at_handoff` as a double scalar.
  - DO NOT surface `lam` — kept internal to C++, not exposed to R.
  - DO NOT change any other field in the result list.
  - DO NOT touch any code path outside the newton_kl result-packing block.
- **Constraints:** Build clean (`R CMD INSTALL --preclean .`). All existing tests pass unchanged. `testthat` suite unchanged (no new test here — WI-3 adds the tests). Atomic.
- **Body sketch:**
  1. Read `src/r_bridge.cpp` — locate the newton_kl result packing block (search for `NewtonCalibResult` or `newton_calibrate`).
  2. In the SEXP list construction (e.g., near `pack_solver_result` or equivalent), add:
     ```cpp
     SET_VECTOR_ELT(result_list, idx++, Rf_ScalarInteger(newton_res.n_warmstart_iters_used));
     SET_VECTOR_ELT(result_list, idx++, Rf_ScalarReal(newton_res.warmstart_max_err_at_handoff));
     ```
     with corresponding `SET_STRING_ELT` name assignments.
  3. Build clean.
  4. Conventional commit: `feat(r_bridge): surface n_warmstart_iters_used + warmstart_max_err_at_handoff for newton_kl`.

---

### WI-3 — T6 + T7 + T8 + T9 + T9b + T10 regression tests

- **Files:** `tests/testthat/test-newton-kl.R`.
- **Blocked by:** WI-1c. R-side field names `res$n_warmstart_iters_used` and `res$warmstart_max_err_at_handoff` must exist in the harvest() result before these tests can be written.
- **Objective:** Six new tests. T6 = lf → λ conversion correctness (structural gate via diagnostic entry point). T7 = stepstone K=9 + warm-start audit + monotonicity. T8 = K=20 severe-skew unit test. T9 = IEPPA-infeasibility cold-start fallback. T9b = K=20 moderate-skew wall regression. T10 = cold-start non-regression.
- **Constraints:** **DO NOT modify T2.** T2's `<1e-4` gate stays. Only ADD T6–T10.
- **WI-2 prerequisite (Fix J, Fix K):** WI-2 must expose test-only C entry point `_newton_kl_warmstart_diagnostic(...)` returning `(lf, λ, max(|u_newton - u_ieppa - C|))` on a tiny fixture. Also clarify that `res.base.iterations` = Newton-only count; warm-start iters live only in `n_warmstart_iters_used` (Fix K).
- **Body sketch:**
  - **T6 (lf → λ conversion — Fix J).** Tiny synthetic n=8, K=2, nj=3 fixture. Call `_newton_kl_warmstart_diagnostic` (test-only C entry point added in WI-2). Assert `max(|u_newton - u_ieppa - C|) < 1e-12` where C is the constant column shift. This directly verifies the conversion math, not just the downstream weights. Label: `test_that("T6: lf->lambda conversion residual < 1e-12 on tiny fixture", {...})`.
  - **T7 (stepstone K=9 + audit + monotonicity — Fix D).** Existing stepstone fixture from T2. Run `harvest(method="newton_kl")`. Assert:
    - `res$max_error < 1e-4`
    - `res$status == 0`
    - `res$n_warmstart_iters_used >= 1`
    - `res$warmstart_max_err_at_handoff < 1e-2`
    - **Monotonicity (Fix D):** `expect_lte(res$max_error, res$warmstart_max_err_at_handoff * 1.01 + 1e-10)` — Newton must not increase max_error beyond the handoff value (1% slack for floating-point rounding; +1e-10 floor for near-zero handoff values).
    Label: `test_that("T7: stepstone K=9 warm-start audit + Newton monotonicity", {...})`.
  - **T8 (K=20 severe-skew unit).** n=10000, K=20, nj=5, target {0.6, 0.2, 0.1, 0.07, 0.03}, `max_weight = 3`. Assert `status %in% c(0L, 1L)` (empirically determined) AND `max_err < 1e-3`. Catches regressions where warm-start fails to break the K=20 drift. Label: `test_that("T8: K=20 severe-skew newton_kl max_err < 1e-3", {...})`.
  - **T9 (IEPPA infeasibility fallback — Fix F).** Synthetic fixture with `max_weight=1` + highly skewed targets (e.g., one category at 0.95, rest near-zero) that triggers IEPPA infeasibility (status ≠ RK_OK from warm-start). Assert `result$status == 0 AND result$n_warmstart_iters_used == 0` — cold-start fallback fired and Newton succeeded from λ=0. Label: `test_that("T9: IEPPA infeasibility triggers cold-start fallback", {...})`.
  - **T9b (K=20 moderate-skew wall regression — Fix H).** n=10000, K=20, nj=5, targets {0.4, 0.25, 0.15, 0.12, 0.08}, `max_weight=3`. Measure `system.time(harvest(method="newton_kl", ...))["elapsed"]`. Assert `< 1.5 × baseline_pre_warmstart_seconds` where baseline is established by running the same fixture with `LBW_NEWTON_KL_WARMSTART_ITERS=0` in the same test. Label: `test_that("T9b: K=20 moderate-skew wall < 1.5x cold baseline", {...})`.
  - **T10 (cold-start non-regression — Fix I).** Stepstone K=9 with `LBW_NEWTON_KL_WARMSTART_ITERS=0` (env-var override per Fix L). Assert `max_err ≤ 2.8e-4` (TH-1a baseline preserved — warm-start off must not regress Newton's cold-start quality). Label: `test_that("T10: cold-start (K_warm=0) non-regression <= 2.8e-4 on stepstone", {...})`.
  - Conventional commit: `test(newton_kl): T6–T10 warm-start regression tests`.

---

### WI-4 — Verify test suite

- **Files:** none (verification only). No commit.
- **Objective:** Run the full testthat suite and confirm zero failures. Specific gate: T2 (stepstone) passes at `<1e-4`. The basin-overlap kill-switch (Fix R) was already run in WI-0b BEFORE any code was written; no re-run needed here.
- **Body:**
  - `Rscript -e "testthat::test_file('tests/testthat/test-newton-kl.R')"` — T1, T2, T3, T4, T5, T6, T7, T8, T9, T9b, T10 all PASS.
  - `Rscript -e "testthat::test_local('.')"` — FAIL=0.
  - **Halt criterion (SPEC_FAILURE):** if T2 still fails (warm-start did not break the basin floor), output `SPEC_FAILURE` and halt. Do NOT proceed to WI-5a.
  - No commit on this ticket — verification only.

---

### WI-5a — kk1204 benchmark (newton_kl + ieppa+sraa baseline)

- **Files:** `benchmarks/results/newton_kl_kk1204.csv`, `benchmarks/newton_kl_bench.R` (if edit needed).
- **Objective:** Run kk1204 severe-skew benchmark with `OMP_NUM_THREADS=1`. Capture wall, iters (Newton + warm-start separately), `max_err`, `n_warmstart_iters_used`, `warmstart_max_err_at_handoff`, `lm_mu_final`, status. **Also capture ieppa+sraa baseline on same fixture (Fix E)** for direct comparison.
- **Body:**
  - Edit `benchmarks/newton_kl_bench.R` to include a `harvest(method="ieppa", accelerate=TRUE)` row on the same kk1204 fixture (Fix E). The CSV must contain one row per method: `newton_kl` (warm-started) and `ieppa_sraa` (baseline). This is the only way to determine whether `newton_kl` + warm-start beats, matches, or loses to plain `ieppa+sraa` on wall time and quality.
  - `OMP_NUM_THREADS=1 Rscript benchmarks/newton_kl_bench.R`
  - Save CSV to `benchmarks/results/newton_kl_kk1204.csv`.
  - Conventional commit: `bench(newton_kl): kk1204 newton_kl warm-start vs ieppa+sraa baseline`.

---

### WI-5b — Verdict + investigation doc + close epic

- **Files:** `docs/investigations/2026-05-01-newton-kl-ieppa-warmstart-result.md` (new); possibly `NEWS.md` (verdict bullet only).
- **Objective:** Write the result document. Decide GATE_MET / PARTIAL / BLOCKED. Close epic appropriately.
- **Decision rules:**
  - **GATE_MET** (kk1204 severe-skew wall `<3s` ∧ `max_err <1e-4` ∧ stepstone `<1e-4`): commit investigation doc + verdict NEWS.md bullet, close Epic-C as resolved. *Honest assessment: GATE_MET is plausible on quality (warm-start solves the structural problem) but unlikely on wall (see §Cost). Expected outcome is PARTIAL.*
  - **PARTIAL** (quality OK on T2/T7/T8 but wall ≥3s on kk1204, OR severe-skew `max_err ∈ [1e-4, 1e-3]`): commit `fix(newton_kl): IEPPA warm-start partial — wall budget unmet` analogue. **Close Epic-C** with the warm-start as the *correct* algorithm; file Epic-D ("Newton-KL wall-time tuning" or "AUTO routing on quality verdict") as a follow-up. PARTIAL is the planned successful closure mode.
  - **BLOCKED** (T2 stays `>1e-4` even with warm-start): the basin-warm-start hypothesis itself is wrong. The investigation doc explains why; epic closes BLOCKED. The follow-up at this point is Epic-E: ESCAPE-only AUTO routing (route severe-skew K≥5 to `ieppa+sraa` directly, leaving Newton-KL for moderate-skew K≥5). This is the rev 2 homotopy plan's contingency, copy-pasted forward.

---

## Epic-C Value Justification (Fix S)

**Why not just use IEPPA+SRAA alone?**

IEPPA+SRAA standalone achieves `1.13e-4` on stepstone — MISSES the `<1e-4` gate by 13%. The Epic-C value-add is the Newton polish that closes that 13% gap. IEPPA alone cannot reach the gate by itself on this fixture.

On kk1204 severe-skew (n=1M, K=20, nj=5, max_weight=3, OMP=1): master Newton-KL DIVERGES (status=1, drift past optimum). IEPPA-warm-started Newton-KL is the only available path that converges to `<1e-4` within a 5s wall budget. The PARTIAL verdict (wall ≥3s) is still a strict improvement over master's diverging behavior (status=1, `max_err >> 1e-4`).

**Summary:** IEPPA alone misses quality gate by 13%. Newton alone diverges on severe-skew. IEPPA-warm-started Newton is the combination that achieves both convergence and quality.

---

## Cost Estimate

Honest, source-derived. Numbers cited from the spec table at `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md:9-21` and the LM rev 2 measurements.

| Component | kk1204 severe-skew (K=20, n=1M) | stepstone (K=9, n=200K) |
|---|---|---|
| IEPPA+SRAA warm-start (8 sweeps) | ≈ 8/10 × 3.7s = **3.0s** | ≈ 8/10 × ~0.3s = **0.24s** |
| LM-Newton tail (1-3 iters from inside basin) | ≈ 1-3 × 0.6s = **0.6 — 1.8s** | ≈ 1-3 × 0.05s = **0.05 — 0.15s** |
| **Total** | **3.6 — 4.8s** | **≈ 0.3 — 0.4s** |

**Honest gate:** kk1204 < 3s. **Realistic outcome:** PARTIAL on wall (3.6-4.8s vs 3s gate). The warm-start is paying for IEPPA's already-3.7s cost and adding Newton on top.

The original Newton-KL spec gate (<2s on kk1204) is dead under any approach that uses IEPPA as a warm-start. The path to <2s would require either (a) IEPPA with fewer sweeps that still land Newton in the basin (need to test K_warm = 4 in a follow-up), or (b) a different first-order pre-conditioner (out of scope; the only one we have that works on these fixtures is IEPPA).

**Quality gate:** stepstone `<1e-4`. **Realistic outcome:** GATE_MET. IEPPA+SRAA already hits 1.13e-4 on stepstone alone; warm-start to Newton from there is mathematically a polish step. Expected stepstone result: `max_err < 1e-7` (Newton converges to its tol_abs from inside the basin).

**Quality gate (T8):** K=20 severe-skew unit, `<1e-3`. **Realistic outcome:** GATE_MET. The K=20 severe-skew failure under bare Newton is a basin-escape problem; warm-start prevents it.

---

## Risks & Open Questions

- **Risk: warm-start lands `λ` inside the basin but Newton's first step still drifts.** Mitigation: existing best-iterate-restoration in `run_newton_inner` (pre-existing TH-1a code at `src/newton_calib.cpp:357-365`) catches this. Worst case: Newton returns its warm-started `λ` unchanged, with quality identical to IEPPA-alone (1.13e-4 on stepstone) — still better than master's 2.8e-4.
- **Risk: IEPPA+SRAA fails on a problem where Newton-KL would have worked.** Mitigation: failure-fallback recovery (Q5) reverts to `λ = 0` and runs Newton as today. Strictly additive.
- **Risk: lf-capture parameter on `ieppa_solve` lands a subtle copy-cost regression.** Mitigation: `lf_capture` defaults to nullptr; the cost is one if-check at function exit. Confirm via WI-1's existing-tests-pass-bit-identically gate — if any IEPPA bench regresses by even 1%, dig into it.
- **Risk: `kNewtonKLWarmstartIters = 8` is wrong.** Mitigation: the constant is derived from a contraction-rate calculation, not pulled from a hat. WI-5a's benchmark CSV captures `n_warmstart_iters_used`; if early-exit fires at 3-4 sweeps consistently, future Epic-D can lower the budget. If post-warm-start Newton consistently runs >5 iters (bad warm-up), Epic-D raises it. The follow-up tuning is small-scope.
- **Open: should T6 use the public R API or a TU-local C++ test harness?** Recommended path: testthat-only via the public API (run IEPPA via `harvest(method="ieppa")` + Newton via `harvest(method="newton_kl")` and compare recovered weights). This is the cleanest scope but doesn't directly probe the lf → λ conversion. Implementer chooses; if testthat-only is too coarse, a tiny standalone C++ unit-test file `tests/testcpp/test_newton_kl_warmstart.cpp` would be the alternative.
- **Open: AUTO routing impact.** This epic does NOT change AUTO routing. If WI-5b lands GATE_MET, AUTO's `K ≥ 5 ∧ M_cell/n ≥ 0.9 → newton_kl` rule becomes more aggressive in practice (newton_kl is now reliable); no code change needed. If WI-5b lands PARTIAL on wall, AUTO routing logic stays as-is (newton_kl is still correct, just not faster than ieppa+sraa on K=20). If BLOCKED, the rev 2 homotopy plan's ESCAPE patch (severe-skew gate routing K=20 to ieppa+sraa) becomes the next epic.

---

## Reviewer Concerns from Epic-B That DO NOT Apply

Pre-empting the plan-review-gate (rev 2 homotopy reviewers caught these; they should not re-raise on Epic-C):

| Concern | Why not applicable here |
|---|---|
| C1 (mechanism handwaving — claims without proof) | §Mechanism includes the exact lf → λ conversion, the algebraic argument that recovered weights are bit-identical (constant column shift cancels in LSE), and the contraction-rate justification for `K_warm = 8`. No "this should work because basin" hand-waving — IEPPA is *empirically known* to converge on overlapping margins (1.13e-4 on stepstone today). |
| C3 (T2 gate loosening) | §Forbidden explicitly forbids it; WI-3 explicitly forbids touching T2. |
| C4 (iter-accumulation contract) | Not relevant — we don't accumulate iter counts across an outer loop. The Newton inner runs once with its full budget; the warm-start uses IEPPA's own counter (one number). |
| C7 (NEWS.md update timing) | WI-2 colocates NEWS.md update with the wiring commit per discipline §4. |
| C10 (per-method internal homotopy) | We don't use per-method homotopy at all. SRAA is acceleration (different concept). `ieppa.homotopy_levels` is bounds-homotopy (different concept again). Both are inert from this epic's perspective. |
| C12 (stagnation criterion plan/ticket disagreement) | No stagnation criterion exists in this plan. Newton uses its existing `tol_abs + max_iter` + best-iterate-restoration. |
| F2 (ε schedule tail jump) | No ε schedule. |
| F3 (`lm_mu` per-ε reset) | No per-ε anything. `lm_mu` is initialised to `1e-3` exactly as today (no change). |
| F4 (T_eps feasibility guard) | No T_eps. |
| F7 (`n_homotopy_levels_used == N` brittle assertion) | T7 uses `>= 1`, not `== K_warm`, exactly because of this concern. |
| F8 (target-skew gate in c_api.cpp) | This epic doesn't touch `c_api.cpp`. AUTO-routing gate is deferred to Epic-D/E if WI-5b lands BLOCKED. |
| S6 (`c_api.cpp:182` ternary, not fall-through) | Out of scope; `c_api.cpp` is in §Forbidden's untouched-files list. |

---

## Out of Scope

- AUTO routing changes (deferred to Epic-D/E based on WI-5b verdict).
- R-side surfacing of `n_warmstart_iters_used`, `warmstart_max_err_at_handoff`, `lm_mu_final` via `r_bridge.cpp` (the prior LM-1c parallel ticket). Filed as Epic-D scope only if WI-5b GATE_MET / PARTIAL and a user-facing tuning loop ever needs it.
- `K_warm` parameter exposed via `CalibState` or `harvest()`. Tunability is a follow-up.
- SIMD / multithreading optimisation of either IEPPA's sweep or Newton's accumulation.
- Different first-order pre-conditioners (raking-warmstart, sinkhorn-warmstart, etc.). IEPPA is chosen because (a) it works empirically on overlapping margins, (b) it has SRAA acceleration already, (c) it has the same `CalibState` ABI as Newton.
- Different K_warm values per fixture / dynamic adjustment based on K. Hard-coded constant; revisit only on bench-driven evidence in a follow-up.
- Touching the IEPPA coordinate-descent loop body for any reason (cosmetic, performance, or otherwise). Off-limits in this epic per §Forbidden.
