# Newton-KL Target Homotopy — Plan

**Epic:** `leafblower-91u7` (Epic-B follow-up, post-NEEDS_HOMOTOPY verdict on epic `leafblower-5k08`)
**Spec reference:** `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` (LM section, plus appendix to add in TH-0)
**Tasks (sequential):** TH-0 (spec) → TH-1 (impl) → TH-1b (μ-reset test) → TH-2 (verify) → TH-3a (bench) → TH-3b (verdict+commit) + parallel LM-1c (r_bridge surfacing of `lm_mu_final`)
**Date:** 2026-05-01

## Context

LM-2 verdict on epic 5k08 was NEEDS_HOMOTOPY. Newton-KL converges to a natural basin floor of `gap ≈ 1.4e-4` on stepstone K=9 (rank-deficient `H_pre` from overlapping margins; ref-cat amplification × ≈ 2 → `max_err ≈ 2.8e-4` vs `<1e-4` gate). The basin floor is intrinsic to the dual-landscape combination of sample × target; LM, LSE, and best-iterate fallback together stabilize convergence to that floor but cannot push past it. To reach `<1e-4`, deform the landscape: target homotopy `T_τ = (1-ε)·T + ε·uniform` with ε annealed `0.5 → 0.1 → 0.01 → 0`. Each intermediate problem has a strictly-feasible interior optimum where Newton converges fully; final ε=0 step starts from a near-optimal warm start, side-stepping the rank-deficient region.

## Mechanism

**Outer ε schedule** (4 levels):
```
ε ∈ {0.5, 0.1, 0.01, 0.0}
```
Geometric-ish; final ε=0 gives the original problem. Smaller ε leaves more "smoothing"; ε=0.5 is the most-smoothed (mostly uniform target); ε=0 is exact.

**Inner solve per ε:**
- Build `T_eps[k] = (1-ε)·T[k] + ε·(1/cat_counts[k])` (uniform per margin, NOT joint uniform).
- Reset `lm_mu = 1e-3` at start of each ε (don't carry warm-start damping that may be stale).
- Use existing LM-damped Newton iter loop with max 20 iters per ε (4 ε × 20 = 80 iter budget total; current `max_iterations` (50) is already in the right ballpark for total budget).
- Per-ε convergence criterion: same `tol_abs` as outer (default 1e-6); but allow ε-stages to stop early if `||∇g_eps||_∞` stops improving for 3 iters in a row (saves budget for ε=0 where precision matters).
- Carry λ across ε boundaries (warm start).

**Final ε=0:** runs the same iter budget; expected to converge in ~3-5 iters from a near-optimal warm start. Best-iterate fallback still active inside the inner loop.

**Recovery:** weights computed from final ε=0's converged λ (no homotopy bias).

**Audit:**
- New diagnostic field `n_homotopy_levels_used` on `NewtonCalibResult` (records actual ε steps that fired before exit).
- Existing `dual_gap` and `lm_mu_final` reflect final state after ε=0.
- Inner gap trajectory captured per ε stage in debug-gated print (`NEWTON_KL_DEBUG=1`).

## Forbidden

- **`exp(u/τ)` recovery-only smoothing** (the user's brainstorm-#2 *literal* reading) — equivalent to a coordinate change, no actual smoothing.
- **Joint-uniform smoothing** (`T_eps[k] uniform over joint cells`) — explosion of cell count, not what's needed.
- **Loosening test gates** to mask incomplete homotopy convergence.
- **Per-ε `lm_mu` carry-over** — fresh damping per ε is the standard interior-point homotopy pattern; carry-over risks ε=0 starting at lm_mu=1e-12 (essentially undamped) on a possibly ill-conditioned problem.
- **Skipping the final ε=0 step** (some homotopy implementations stop at small ε; we MUST hit ε=0 for unbiased weights).
- **Touching files outside this list:** `src/newton_calib.cpp`, `src/newton_calib.hpp`, `src/r_bridge.cpp` (LM-1c only — surface `lm_mu_final` and `n_homotopy_levels_used`), `tests/testthat/test-newton-kl.R`, `benchmarks/newton_kl_bench.R`, `benchmarks/results/`, `NEWS.md`, the spec file (TH-0 only).
- **Skipping pre-commit hooks.**

## Plan Steps

### TH-0 — Spec amendment
- Append "Target homotopy" subsection to `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` after the LM section: ε schedule, inner solve, μ-reset semantics, n_homotopy_levels_used field, AUTO-routing impact (none; homotopy is internal).
- Conventional commit `docs(spec): Newton-KL — add target homotopy section`.

### TH-1 — Implement target-homotopy outer loop
- Edit `src/newton_calib.cpp`:
  - Refactor existing iter loop into a callable inner: `run_newton_inner(T_eps, max_iter_inner, lam, &lam_best, &best_gap, ...)`.
  - Add outer `for (ε in {0.5, 0.1, 0.01, 0.0})` loop that builds `T_eps`, resets `lm_mu=1e-3`, calls inner.
  - Track total iters across ε stages → set `res.base.iterations = sum`.
  - Track final ε's `dual_gap`, `lm_mu`, `step_norm` in result struct.
  - Add `int n_homotopy_levels_used = 0;` to `NewtonCalibResult` (touch `newton_calib.hpp`).
- DO NOT change recovery (step 6/7); it consumes final λ and Z_curr from ε=0.
- Build clean: `R CMD INSTALL --preclean .` zero new warnings.

### TH-1b — μ-reset + ε-step regression test
- Add T7 to `tests/testthat/test-newton-kl.R`:
  - Synthetic K=3 well-conditioned: assert `n_homotopy_levels_used == 4` (all ε steps fire).
  - Stepstone K=9 stress: assert `max_err < 1e-4` after homotopy (the gate the original epic missed).
- Optionally amend T2 expectations to reflect the now-passing fixture.

### TH-2 — Verify test suite
- `Rscript -e "testthat::test_file('tests/testthat/test-newton-kl.R')"` ⇒ T1-T7 all PASS.
- `Rscript -e "testthat::test_local('.')"` ⇒ FAIL=0.
- If T2 still fails (stepstone homotopy did not break the basin floor), STOP — escalate to BLOCKED.

### TH-3a — Run kk1204 benchmark
- Same as LM-3a: `OMP_NUM_THREADS=1 Rscript benchmarks/newton_kl_bench.R`.
- CSV at `benchmarks/results/newton_kl_kk1204.csv`.

### TH-3b — Verdict + commit + AUTO routing
- **GATE_MET** (severe-skew wall<3s ∧ max_err<1e-4): commit final implementation + investigation doc, close epic.
- **PARTIAL** (severe converges, but max_err in [1e-4, 1e-3]): document, ship as best-effort, draft Epic-C (IEPPA warm-start) follow-up.
- **ESCAPE** (severe still diverges with homotopy): ship AUTO-routing safety guard at `src/c_api.cpp` line ~171 — guard `K≥5 && M_cell/n≥0.9 → newton_kl` rule with target-skew check (`max_T/min_T > 5 ⇒ fall through to ieppa+sraa`).

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
