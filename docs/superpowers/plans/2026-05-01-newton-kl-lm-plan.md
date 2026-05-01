# Newton-KL LM-Damping Convergence Fix — Plan

**Epic:** `leafblower-5k08`
**Tasks:** `leafblower-5k08.1` (LM-1) → `leafblower-5k08.2` (LM-2) → `leafblower-5k08.3` (LM-3)
**Spec reference:** `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md`
**Date:** 2026-05-01

---

## Mechanism

**Target pattern:** Levenberg-Marquardt damped Newton with scale-invariant damping
$$(H + \mu \cdot \text{diag}(H)) \cdot \delta = \nabla g, \quad \lambda \mathrel{-}= \alpha \cdot \delta$$

**Adaptive μ:** init `μ = 1e-3`; on Armijo full-step accept (α=1) `μ ← max(μ/3, 1e-12)`; on Armijo failure (line search exhausted) `μ ← min(μ·10, 1e12)`.

**Audit:** spy on `dual_gap` trajectory across iters; in-place inspection of `μ` evolution via debug-gated print (env var `NEWTON_KL_DEBUG`); existing testthat T1-T5 cover correctness; new bench script verifies kk1204 wall+max_err.

## Forbidden

- **µ·I** (spherical) damping — must be `µ·diag(H)`.
- **Trust-region clip on δ** — LM replaces it; both together over-regularize.
- **Larger LDLT pivot floor (>1e-12)** — LM makes the matrix well-conditioned by construction; pivot bump is symptom-suppression.
- **Coordinate change tricks** (logit reparameterization, IPF/Sinkhorn substitution) — keep the spec algorithm structure.
- **Loosening test gates to mask convergence failure** — if T2 still fails, escalate to Epic-B (target homotopy), do not loosen.
- **Fixture changes to bench** — kk1204 severe-skew params (0.6/0.2/0.1/0.07/0.03 + max_weight=3 + n=1M + K=20) are spec-fixed.
- **Skipping pre-commit hooks**.
- **Touching files outside `src/newton_calib.cpp`, `benchmarks/newton_kl_bench.R`, `benchmarks/results/`, and one new investigation doc.**

## Diagnosis (Verified 2026-05-01)

1. **Original Newton-KL** produced NaN weights on K=20 severe-skew via exp(u_i) overflow → `Z=0` → false convergence (gap=0 from NaN comparisons).
2. **A1 — LSE stabilization** (already shipped in current `newton_calib.cpp`): pre-compute `u_max(λ)`, accumulate `f_i = d_i · exp(u_i − u_max)`. Eliminates overflow. **Keep all of it.**
3. **Bandaids** (currently in code, to be removed): trust-region clip κ=1, pivot floor 1e-8/1e-4. Together produced regression on stepstone K=9 (T2 max_err = 3.85e-4 vs gate 1e-4) by over-regularizing well-conditioned problems.
4. **Root cause identified:** sample Hessian rank-collapses dynamically as λ drifts into low-density exp(u) regions. The few obs with large `f_i` dominate; H eigenvalues span many orders of magnitude. Trust-region clipping cripples direction in collapsed modes; pivot floor doesn't repair direction.

## Why LM Scale-Invariant

`diag(H)_a = p_a(1−p_a)` reflects sampled cell mass for category `a`, which varies by orders of magnitude under severe skew. Adding `μ·diag(H)` (rather than `μ·I`) interpolates Newton ↔ scaled gradient descent **per-coordinate**, without imposing a uniform identity floor that mismatches H's scale. Standard textbook fix for the exact diagnosed disease.

Adaptive μ via Armijo outcome:
- Full-step accept ⇒ μ ↘ /3 (recover Newton near optimum).
- Line search exhausted ⇒ μ ↗ ×10 (back off toward gradient descent in collapsed-rank directions).

## Plan Steps

1. **LM-1 (`leafblower-5k08.1`)** — Edit `src/newton_calib.cpp`:
   - Confirm `double mu = 1e-3; mu_min=1e-12; mu_max=1e12;` declared above iter loop (already added in prior turn).
   - Replace trust-region block (lines ~184-208) with damp-then-LDLT: `for a: H[a*n+a] *= (1+μ)`, then `ldlt_factor_inplace(H, n_lam, 1e-12)`.
   - Add μ-update after Armijo: `if (accepted && α≥0.999) μ /= 3 else if (!accepted) μ *= 10`.
   - Remove `trust_radius`, `delta_max`, `clip_scale` references entirely.
   - Verify build: `R CMD INSTALL --preclean . | tail -3` ⇒ `* DONE`.
2. **LM-2 (`leafblower-5k08.2`)** — Run testthat:
   - `Rscript -e "testthat::test_file('tests/testthat/test-newton-kl.R')"` ⇒ T1-T5 all PASS.
   - `Rscript -e "testthat::test_local('.')"` ⇒ FAIL=0.
   - If T2 still FAIL with same gradient-stuck symptom ⇒ verdict=NEEDS_HOMOTOPY, halt epic, draft Epic-B ticket.
3. **LM-3 (`leafblower-5k08.3`)** — Run `benchmarks/newton_kl_bench.R`:
   - Capture severe + moderate skew rows.
   - Decide verdict per rules:
     * **GATE_MET**: severe wall<2s ∧ max_err<1e-4 ⇒ commit `feat(newton_kl): LM solves kk1204 severe-skew in <2s`, close epic.
     * **NEEDS_HOMOTOPY**: severe converges but misses gate ⇒ commit partial result, draft Epic-B (target homotopy `T_τ = (1−ε)T + ε·uniform`, ε: 0.5→0.01→0).
     * **ESCAPE**: severe diverges or max_err stuck even with LM ⇒ draft Epic-D (drop newton_kl from AUTO K=20 path; route to ieppa+sraa).
   - Write `docs/investigations/2026-05-01-newton-kl-lm-result.md` with numbers + verdict + reasoning.

## Test & Validation Strategy

- **Unit:** `tests/testthat/test-newton-kl.R` T1-T5 (already authored; T2 currently failing). Must pass without bound loosening.
- **Integration:** `testthat::test_local('.')` full suite — FAIL=0 baseline.
- **Benchmark:** `benchmarks/newton_kl_bench.R` — kk1204 severe + moderate.
- **Comparison:** newton_kl vs ieppa+sraa on identical fixtures, single-thread (`OMP_NUM_THREADS=1`), 2 bench iterations.
- **No mocks:** problem is numerical correctness; mocks would mask failure modes.

## Risks & Open Questions

- **Risk: LM alone insufficient for severe-skew K=20.** Mitigation: explicit verdict path NEEDS_HOMOTOPY → Epic-B. Honest failure reporting required.
- **Risk: μ schedule (÷3 / ×10) not optimal.** Standard LM literature uses ÷10 / ×10 or ÷2 / ×2. Choice of /3 / ×10 is hand-tuned for our problem. If LM-2 reveals slow convergence, consider tuning in follow-up.
- **Open: cost overhead.** LM adds one diagonal scan per iter (~80 ops, negligible). LDLT cost unchanged. Per-iter wall ≈ existing.
- **Open: should μ be returned as part of `NewtonCalibResult`?** Defer — keep struct unchanged this iteration; add `lm_mu_final` field in a follow-up if useful.

## Out of Scope

- Target homotopy / temperature smoothing (= Epic-B, conditional follow-up).
- IEPPA warm-start (= Epic-C, alternative path).
- AUTO routing changes (= Epic-D11, escape path).
- Public API changes (`R/harvest.R`, `c_api.cpp`, `r_bridge.cpp`, `newton_calib.hpp`) — none.
- Multithreading / SIMD (separate optimization concern).
