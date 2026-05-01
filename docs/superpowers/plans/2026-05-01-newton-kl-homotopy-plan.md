# Newton-KL Target Homotopy — Plan (rev 2)

## Rev 2 changelog (post plan-review-gate)

3 reviewers NEEDS_REVISION on rev 1. Critical fixes:
- **§Context mechanism reframed**: warm-start basin trick, not rank fix. T_eps blends targets, NOT sample structure; rank deficiency persists, but Newton starting near `λ*(T_0)` from `λ*(T_eps→0)` stays inside the convergence basin instead of drifting (F1).
- **6-level ε schedule** `{0.5, 0.1, 0.03, 0.01, 0.003, 0}` (was 4-level; closes 100× jump 0.01→0) (F2).
- **`lm_mu = max(prev/3, 1e-6)` per ε** (was hard reset to 1e-3) — keeps recovered damping while avoiding saturation lock-in (F3).
- **Stagnation criterion removed** (plan/ticket disagreement); rely on `tol_abs` + `max_iter_inner=20` (C12).
- **TH-1 split** → TH-1a (refactor inner) + TH-1b-impl (homotopy wrap). Each atomic, independent revertible (S4).
- **NEWS.md update lives in TH-1b-impl commit** (algorithmic change colocates with code, per discipline §4) (C7/S3).
- **T2 gate loosening explicitly forbidden** in TH-1b-impl Forbidden block (C3/S7).
- **K=20 severe-skew T8 added** to TH-1c (formerly TH-1b) regression tests (C1/F5).
- **Iter-accumulation contract** explicit in TH-1b-impl (inner writes local count to `res.base.iterations`; outer reads-then-accumulates `total_iters` *before* next inner call overwrites) (C4).
- **`T_eps` feasibility check**: at start of each ε, abort with `RK_ERR_BADARG` if no obs has positive sample-frequency for any cat with `T_eps[k][j] > 0` (F4).
- **n_homotopy_levels_used semantics**: `==` count of ε levels with ≥1 inner iter; T8 asserts `>=1 && status==0` not `==4` (F7).
- **ESCAPE patch**: `src/c_api.cpp:182` is a ternary, NOT a fall-through chain. Rewrite to add `target_skew = max_T/min_T` gate that selects `RK_ALG_IEPPA` when severe (S6). Homotopy code STAYS in for moderate-skew K=20 use (C10).
- **PARTIAL closes Epic-B + files Epic-C** (IEPPA warm-start) explicit (S8/C11).
- **Honest gate documented**: <3s aspirational; <10s realistic; PARTIAL is the expected closure mode (F7 elaboration).

**Epic:** `leafblower-91u7` (Epic-B follow-up, post-NEEDS_HOMOTOPY verdict on epic `leafblower-5k08`)
**Spec reference:** `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` (LM section, plus appendix to add in TH-0)
**Tasks (sequential):** TH-0 (spec) → TH-1 (impl) → TH-1b (μ-reset test) → TH-2 (verify) → TH-3a (bench) → TH-3b (verdict+commit) + parallel LM-1c (r_bridge surfacing of `lm_mu_final`)
**Date:** 2026-05-01

## Context (rev 2)

LM-2 verdict on epic 5k08 was NEEDS_HOMOTOPY. Newton-KL on stepstone K=9 converges to a natural basin floor at `gap ≈ 1.4e-4` (from there λ drifts past optimum into unbounded-dual regions; best-iterate fallback rolls back to the floor; ref-cat amplification × ≈ 2 → `max_err ≈ 2.8e-4` vs `<1e-4` gate).

**Mechanism (corrected per F1):** the rank deficiency in `H_pre` comes from the *sample's overlapping-margin support structure* — empty `(j_{k1}, j_{k2})` cells across margin pairs. **`T_eps` does not fix this** (it doesn't change which obs share which (j_k) patterns). What `T_eps` *does* is shift `λ*(T_eps) → λ*(T_0)` continuously: each ε's optimum is a warm start for the next ε's Newton. The rank deficiency is still present at every ε, but Newton starting near a warm start with smaller `||λ_init − λ*||` stays inside the basin where its quadratic model is locally accurate. Without homotopy, Newton at `λ=0` is far enough from `λ*` that the model misleads it past optimum into rank-collapsed drift territory.

This is the standard interior-point / continuation-method argument: each problem on the homotopy path is *easier than the next* in basin-of-attraction terms, not in rank-of-Hessian terms. Final ε=0 step inherits a `λ_init` close enough to `λ*(T_0)` that even the rank-deficient Newton stays inside its quadratic basin and converges to higher precision than from `λ=0`.

## Mechanism (rev 2)

**Outer ε schedule** (6 levels):
```
ε ∈ {0.5, 0.1, 0.03, 0.01, 0.003, 0.0}
```
Geometric ratio ~3× (was 5× / 10× / 100× — the 100× tail jump 0.01→0 was risky). Each ε's `λ*` is closer to the next ε's `λ_init` than under the 4-level schedule.

**Inner solve per ε:**
- Build `T_eps[k][j] = (1-ε)·T[k][j] + ε·(1/cat_counts[k])` (per-margin uniform; NOT joint uniform).
- **Feasibility guard** (F4): for each margin k and category j with `T_eps[k][j] > 0`, verify there exists at least one obs with `cat_k(i) == j AND d_i > 0`. If any (k,j) has zero sample support but positive target, the dual is unbounded for that direction; abort with `RK_ERR_BADARG` and `lm_mu_final` left at the entry value.
- **`lm_mu` reset semantics** (F3): `lm_mu_init_eps = max(lm_mu_final_prev / 3.0, 1e-6)` (i.e., un-clamp from saturation but keep recovered damping). At ε=ε[0] (first), use `lm_mu = 1e-3`.
- Existing LM-damped Newton iter loop with max 20 iters per ε.
- Per-ε convergence criterion: `||∇g_eps||_∞ < tol_abs` (same as outer). NO stagnation criterion (revert to plain max_iter + tol).
- λ carries across ε boundaries (warm start).

**Final ε=0:** runs same inner. From `λ*(0.003) ≈ λ*(0)` warm start, expected ~3-5 Newton iters. Best-iterate fallback still active (rollback within inner if ε=0 still drifts past optimum).

**Recovery:** weights from final ε=0's `λ`. Unbiased.

**Audit:**
- `n_homotopy_levels_used` field on `NewtonCalibResult`: count of ε levels with ≥1 inner Newton iter (NOT just "attempted"). T8 asserts `>= 1 && status == 0`.
- `dual_gap`, `lm_mu_final`, `step_norm` reflect final ε=0 state.
- Per-ε gap trajectory in debug-gated print (`NEWTON_KL_DEBUG=1`).
- Total iters: `res.base.iterations = sum_over_ε(inner_iters_at_that_ε)`. Inner writes local iter count to `res.base.iterations` per call; outer reads it BEFORE the next inner call overwrites and accumulates into `total_iters` local. After loop: `res.base.iterations = total_iters`.

## Forbidden (rev 2)

- **`exp(u/τ)` recovery-only smoothing** (the user's brainstorm-#2 *literal* reading) — equivalent to a coordinate change λ' = λ/τ, no actual smoothing.
- **Joint-uniform smoothing** (`T_eps[k]` uniform over joint cells) — explosion of cell count.
- **Loosening any test gate.** T2 `<1e-4` stays. T8 (new) gate `<1e-4`. NO relaxation to mask convergence shortfall.
- **Hard `lm_mu` reset to `1e-3`** per ε (rev 1 had this; reviewer F3 caught: leaves recovered damping behind). Use `max(prev_final / 3, 1e-6)` instead.
- **4-level ε schedule with 100× tail jump** (rev 1 used `{0.5, 0.1, 0.01, 0}`). Use 6-level `{0.5, 0.1, 0.03, 0.01, 0.003, 0}`.
- **Stagnation criterion** as additional inner stop. Plan rev 1 had this informally; rev 2 removes (kept simple: `tol_abs + max_iter_inner=20`).
- **Skipping the final ε=0 step.** Must always hit ε=0.
- **Skipping `T_eps` feasibility guard.**
- **Per-method internal homotopy reuse confusion.** `ieppa.homotopy_levels` is bounds-homotopy (different concept). DO NOT couple to it. `harvest()` API unchanged.
- **`n_homotopy_levels_used == 4` brittle assertion** (rev 1 had this). Use `>= 1 && status == 0` semantics.
- **Test for K=20 severe-skew omitted.** T8 mandatory in TH-1c.
- **NEWS.md update deferred to TH-3b.** Algorithmic change colocates with code in TH-1b-impl per discipline §4.
- **Touching files outside:** `src/newton_calib.cpp`, `src/newton_calib.hpp`, `src/r_bridge.cpp` (LM-1c only), `tests/testthat/test-newton-kl.R`, `benchmarks/newton_kl_bench.R`, `benchmarks/results/`, `NEWS.md`, the spec file (TH-0 only). On ESCAPE only, also `src/c_api.cpp`.
- **Skipping pre-commit hooks.**

## Plan Steps (rev 2)

Sequenced: TH-0 → TH-1a → TH-1b-impl → TH-1c (was TH-1b) → TH-2 → TH-3a → TH-3b. LM-1c parallel after TH-1b-impl.

### TH-0 — Spec amendment
- Append "Target homotopy" subsection to `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` after the LM section: ε schedule, inner solve, μ-reset semantics, n_homotopy_levels_used field, AUTO-routing impact (none; homotopy is internal).
- Conventional commit `docs(spec): Newton-KL — add target homotopy section`.

### TH-1a — Refactor inner loop (no behavioral change)
- Edit `src/newton_calib.cpp`: extract iter loop body into a callable lambda `run_newton_inner(const std::vector<double>& T_eps, int max_iter_inner) -> bool`.
  - Inner contract: takes target vector + iter cap; mutates outer-scope `lam`, `lm_mu`, `Z_curr`, `u_max_curr`, `res` (status, dual_gap, step_norm, line_alpha, iterations, lm_mu_final). Returns `true` if inner exited via converged path (status=0 or step-stall); `false` if BADARG/NOCONV (caller halts outer loop).
  - Inner exit paths to handle: convergence-break (gap<tol), step_norm<1e-14 (mark_converged), LDLT 3× fail (RK_ERR_BADARG, return false), 3 consecutive line-search failures (RK_ERR_NOCONV, return false). Best-iterate restoration runs INSIDE inner before return.
  - Inner writes `res.base.iterations` to its LOCAL iter count (NOT cumulative).
- Add `int n_homotopy_levels_used = 0;` to `NewtonCalibResult` in `src/newton_calib.hpp` (default 0; TH-1b-impl will set it).
- Top-level call site: invoke once with original `T` and `max_iter` from CalibState. Behavior must be bit-identical to current master.
- Build clean: `R CMD INSTALL --preclean .` zero new warnings.
- All existing tests (T1-T5) still pass.
- Conventional commit `refactor(newton_kl): extract iter loop to run_newton_inner; add n_homotopy_levels_used field`.

### TH-1b-impl — Add target-homotopy outer loop (was TH-1)
- Edit `src/newton_calib.cpp`:
  - Above `run_newton_inner` definition: add `const double EPS_SCHEDULE[] = {0.5, 0.1, 0.03, 0.01, 0.003, 0.0};` and `constexpr int N_EPS = 6`.
  - Replace top-level single inner-call with outer ε loop:
    ```cpp
    int total_iters = 0;
    res.n_homotopy_levels_used = 0;
    std::vector<double> T_eps(n_lam, 0.0);
    for (int e = 0; e < N_EPS; ++e) {
        const double eps = EPS_SCHEDULE[e];
        // Build T_eps: per-margin uniform.
        for (int k = 0; k < K; ++k)
            for (int j = 1; j < st.cat_counts[k]; ++j)
                T_eps[lam_off[k] + j - 1] = (1.0 - eps) * st.targets[k][j]
                                             + eps / st.cat_counts[k];
        // Feasibility guard: every (k,j) with T_eps>0 must have ≥1 obs.
        // (Cheap O(n·K) check on first ε; subsequent ε share support.)
        if (e == 0 && /* feasibility violated */) {
            res.base.status = RK_ERR_BADARG;
            std::snprintf(res.message, sizeof(res.message),
                "newton_kl: T_eps infeasible — empty (margin,cat) with positive smoothed target");
            res.lm_mu_final = lm_mu;
            return res;
        }
        // lm_mu init for this ε.
        lm_mu = (e == 0) ? 1e-3 : std::max(res.lm_mu_final / 3.0, 1e-6);
        // Run inner; iter count captured PER call.
        const bool ok = run_newton_inner(T_eps, /*max_iter_inner=*/20);
        const int iters_this_eps = res.base.iterations;
        total_iters += iters_this_eps;
        if (iters_this_eps > 0) res.n_homotopy_levels_used = e + 1;
        if (!ok) break;  // BADARG / NOCONV from inner — halt outer.
    }
    res.base.iterations = total_iters;
    ```
- NEWS.md (in same commit, per discipline §4 / reviewer C7): add bullet under `## Newton-KL calibration`:
  > * `method="newton_kl"` now wraps the LM-damped Newton solver in a 6-level target-homotopy outer loop (ε ∈ {0.5, 0.1, 0.03, 0.01, 0.003, 0}). Smooths targets toward per-margin uniform and anneals to original; warm-starts Newton across ε boundaries to escape rank-deficient drift on overlapping-margin fixtures (e.g., stepstone K=9). New result field: `n_homotopy_levels_used`.
- Build clean.
- Conventional commit `feat(newton_kl): target-homotopy outer loop with 6 ε levels`.

### TH-1c — Add T7 + T8 regression tests (was TH-1b)
- Append to `tests/testthat/test-newton-kl.R`:
  - **T7: stepstone K=9 + n_homotopy_levels_used.** As in rev 1's T7 but assertion changed to `n_homotopy_levels_used >= 1 && status == 0` (NOT `==4`).
  - **T8: kk1204-style K=20 severe-skew unit test.** n=10000 (small enough for testthat), K=20, nj=5, target {0.6,0.2,0.1,0.07,0.03}, max_weight=3. Assert `status == 0` (or `1` if bounds violation expected; TBD by implementer based on actual bounds-violation rate at this n) AND `max_err < 1e-3` (looser than stepstone's 1e-4 because severe skew at K=20 is a harder regime). T8 catches regressions where homotopy doesn't break the K=20 drift.
  - **NO amendment to T2.** T2's `<1e-4` gate stays. Forbidden block in this ticket explicitly: "DO NOT modify T2; only ADD T7 + T8."

### LM-1c (parallel after TH-1b-impl) — Surface lm_mu_final + n_homotopy_levels_used in r_bridge
- Same as rev 1.

### TH-2 — Verify test suite
- `Rscript -e "testthat::test_file('tests/testthat/test-newton-kl.R')"` ⇒ T1-T7 all PASS.
- `Rscript -e "testthat::test_local('.')"` ⇒ FAIL=0.
- If T2 still fails (stepstone homotopy did not break the basin floor), STOP — escalate to BLOCKED.

### TH-3a — Run kk1204 benchmark
- Same as LM-3a: `OMP_NUM_THREADS=1 Rscript benchmarks/newton_kl_bench.R`.
- CSV at `benchmarks/results/newton_kl_kk1204.csv`.

### TH-3b — Verdict + investigation + AUTO routing (rev 2)

Decision rules:
- **GATE_MET** (severe-skew wall<3s ∧ max_err<1e-4): commit investigation doc + verdict bullet in NEWS.md, close epic. *Realistically unreachable on wall — see §Cost.*
- **PARTIAL** (severe converges with `max_err < 1e-4` but wall ≥3s, OR `max_err ∈ [1e-4, 1e-3]`): commit `fix(newton_kl): target homotopy partial — wall budget unmet` (or analogous). **Close Epic-B AND file Epic-C** (IEPPA warm-start) explicit follow-up. PARTIAL is the planned successful closure mode.
- **ESCAPE** (severe diverges or `lm_mu_final` saturated): commit `fix(newton_kl,auto): homotopy ESCAPE — AUTO routing target-skew gate`. Patch:
  - **`src/c_api.cpp:182`** is currently `alg = (K >= 5) ? RK_ALG_NEWTON_KL : RK_ALG_RAKING;` (a ternary, NOT a fall-through chain — reviewer S6). **Rewrite** to:
    ```cpp
    // Severe target skew: newton_kl drifts (post-homotopy verdict ESCAPE).
    // Route to ieppa+sraa for K≥5 severe-skew problems.
    double max_T = 0.0, min_T = 1.0;
    for (int k = 0; k < K; ++k)
        for (int j = 0; j < cat_counts[k]; ++j) {
            const double t = targets[k][j];
            if (t > max_T) max_T = t;
            if (t > 1e-12 && t < min_T) min_T = t;
        }
    const bool severe_skew = (max_T / std::max(min_T, 1e-12) > 5.0);
    if (K >= 5 && !severe_skew) alg = RK_ALG_NEWTON_KL;
    else if (K >= 5 && severe_skew) alg = RK_ALG_IEPPA;  // skew → ieppa+sraa
    else alg = RK_ALG_RAKING;
    ```
  - Mirror equivalent change at `src/r_bridge.cpp:425`.
  - **Homotopy code STAYS** in `newton_calib.cpp` (still useful for moderate-skew K=20 calls when AUTO selects newton_kl, plus explicit `method="newton_kl"` calls from users).

### LM-1c (parallel) — Surface `lm_mu_final` + `n_homotopy_levels_used` in r_bridge
- Edit `src/r_bridge.cpp` SEXP-packing for `NewtonCalibResult`: add `lm_mu_final` and `n_homotopy_levels_used` fields to the result list.
- Required by TH-1b's μ-reset test which asserts `n_homotopy_levels_used == 4`.
- Atomic single-file change.

## Cost Estimate

- 4 ε-levels × ~5 inner Newton iters each = 20 inner iters.
- Per-iter cost on K=20 n=1M ≈ 0.6-1.0s (per LM-2 measurements).
- Total wall: ~12-20s on kk1204 severe-skew. **Above the <3s honest gate.** GATE_MET unlikely on wall basis even if max_err passes.
- Stepstone K=9 (n=200K): per-iter ~0.05s × 20 = 1s. Within budget.

If wall budget is the binding constraint, TH-3b verdict will be PARTIAL (correct algorithm, too slow).

## Risks & Open Questions

- **Risk: homotopy doesn't push past stepstone basin floor.** Mitigation: TH-2 STOP-on-failure path; if so, the rank-deficiency is structural (not landscape-curvature) and homotopy can't fix it. Escalate to Epic-C (IEPPA warm-start).
- **Risk: per-ε `lm_mu=1e-3` reset is wrong direction.** Some homotopy methods carry forward warm μ. Standard practice for interior-point methods is reset; for trust-region-Newton-LM, mixed. We document the choice and tag it as tunable.
- **Risk: ε schedule {0.5, 0.1, 0.01, 0} too coarse.** Geometric (5-level): {0.5, 0.158, 0.05, 0.0158, 0} smoother. Defer tuning to TH-2 results.
- **Open: should `tol_abs` tighten across ε steps?** Standard practice: looser early, tighter for final. We use the same tol throughout to keep it simple; revisit if convergence stalls at intermediate ε.

## Out of Scope

- IEPPA warm-start (= Epic-C, alternative path).
- AUTO routing changes beyond ESCAPE-only safety guard (= follow-up).
- SIMD / multithreading optimization of inner loop.
- Different smoothing kernels (e.g., Dirichlet prior with non-uniform `α_kj`).
