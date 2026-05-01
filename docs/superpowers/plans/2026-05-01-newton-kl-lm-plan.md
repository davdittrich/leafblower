# Newton-KL LM-Damping Convergence Fix — Plan (rev 2)

**Epic:** `leafblower-5k08`
**Tasks:** LM-0 (spec amend) → LM-1 (impl) → LM-2 (test) → LM-3a (bench) → LM-3b (verdict+commit)
**Spec reference:** `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md`
**Date:** 2026-05-01

## Rev 2 changelog (post plan-review-gate)

3 reviewers flagged NEEDS_REVISION. Critical fixes incorporated:
- **F1 gain-ratio μ-update** replaces naive `α=1` test (Marquardt criterion ρ>0.75).
- **F2 Armijo 3-way bookkeeping** (FULL_ACCEPT / BACKTRACKED / FAILED); on FAILED retry with new μ, do not take tiny step.
- **F5 additive damping floor** prevents `μ·0=0` where rank collapses: `H[a,a] = max(H[a,a]·(1+μ), μ·d_floor)`, `d_floor = mean(diag(H))`.
- **C7 μ-trajectory observable**: `lm_mu_final` field added to `NewtonCalibResult`.
- **S1 spec amendment** as LM-0 (atomic, before code).
- **S2** rename `mu` → `lm_mu` (disambiguates `st.alm.mu` in lbfgsb_solver).
- **S3 NEWS.md** entry in same commit.
- **S4 ESCAPE safety patch** moves into LM-3b (not deferred).
- **F7 honest gate**: <3s not <2s (per-iter ~1s × ≥3 iters at K=20 n=1M is physically lower bound).
- **F4 μ_init=1.0** (more conservative; iter-0 starts near gradient descent without trust region).

---

## Mechanism

**Target pattern:** Levenberg-Marquardt damped Newton with scale-invariant damping + additive floor:
$$H_{\text{damped}}[a,a] = \max\!\big(H[a,a]\cdot(1+\mu),\; \mu \cdot d_{\text{floor}}\big), \quad d_{\text{floor}} = \overline{\mathrm{diag}(H)}$$
$$H_{\text{damped}} \cdot \delta = \nabla g, \quad \lambda \mathrel{-}= \alpha \cdot \delta$$

The additive floor ensures damping is non-zero even when `H[a,a]→0` (rank-collapsed direction — exactly the failure mode). Multiplicative-only `(1+μ)` was insufficient because it preserves zero diagonals.

**Adaptive μ via Marquardt gain ratio:**
$$\rho = \frac{g_{\text{curr}} - g_{\text{trial}}}{\alpha (G\cdot\delta) - \tfrac12 \alpha^2 \delta^\top H \delta}$$

Schedule:
- `ρ > 0.75` AND Armijo accepted ⇒ `lm_mu ← max(lm_mu/3, 1e-12)` (good prediction → trust Newton).
- `ρ < 0.25` OR Armijo exhausted ⇒ `lm_mu ← min(lm_mu·10, 1e12)` (poor prediction → back off toward grad descent).
- `0.25 ≤ ρ ≤ 0.75` ⇒ keep `lm_mu` unchanged.

**Armijo three-way outcome:**
- `FULL_ACCEPT` (α=1, accepted): proceed; ρ-test decides μ.
- `BACKTRACKED` (α<1, accepted): proceed; μ unchanged.
- `FAILED` (line search exhausted): **do NOT take a tiny step**. Increase `lm_mu ×= 10`, **redo** the LDLT solve in the same iter using new μ. Cap retries at 3 to prevent infinite loop. If still failing after 3 retries, return `RK_ERR_NOCONV`.

**Init:** `lm_mu = 1.0` (more conservative than 1.0e-3 — iter 0 starts near gradient descent, the /3 schedule walks down to Newton as iterations succeed). Removes need for trust region.

**Audit:**
- `lm_mu_final` recorded in `NewtonCalibResult` (new field) for LM-3 verdict diagnosis (saturation at 1e12 → use Epic-B; bounded → other pathology).
- `dual_gap` trajectory across iters (existing).
- T6 new unit test: synthetic 2-iter problem asserts μ-schedule correctness.
- T1 amended: assert exit_reason ≠ "step_norm" (must exit via gradient criterion, not false step-norm convergence).
- testthat T1-T5 carry correctness; bench script verifies kk1204 wall+max_err.

## Forbidden

- **µ·I** (spherical) damping — must be `µ·diag(H)` with additive floor.
- **Multiplicative-only damping** (`H[a,a] *= (1+μ)` without floor) — leaves zero diagonals undamped exactly where rank collapses.
- **Trust-region clip on δ** — LM replaces it; both together over-regularize.
- **Naive "α=1 ⇒ μ/=3" schedule** — must use Marquardt gain ratio ρ; α=1 alone proves nothing.
- **Tiny-step Armijo fallback** (`α=1/2^30`) — silently advances λ with bogus step. Replace with FAILED → retry with larger μ.
- **Coordinate change tricks** (logit reparameterization, IPF/Sinkhorn substitution) — keep the spec algorithm structure.
- **NaN-tolerant comparisons** that allow `Z=0` to pass convergence (the original bug).
- **`exp(u)` recovery without `u_max` stabilization** (the LSE invariant; already enforced).
- **λ-norm clipping** as substitute for damping (was the bandaid).
- **Renaming `mu` to anything other than `lm_mu`** (must disambiguate from `st.alm.mu` in lbfgsb_solver.cpp).
- **Loosening test gates to mask convergence failure** — if T2 still fails, escalate to Epic-B (target homotopy), do not loosen.
- **Fixture changes to bench** — kk1204 severe-skew params (0.6/0.2/0.1/0.07/0.03 + max_weight=3 + n=1M + K=20) are spec-fixed.
- **Skipping pre-commit hooks**.
- **Touching files outside this list:** `src/newton_calib.cpp`, `src/newton_calib.hpp`, `tests/testthat/test-newton-kl.R`, `benchmarks/newton_kl_bench.R`, `benchmarks/results/`, `NEWS.md`, `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` (LM-0 only), one new investigation doc, and (only on ESCAPE verdict) `src/c_api.cpp` + `src/r_bridge.cpp` for AUTO-routing safety patch.

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

### LM-0 — Spec amendment (new ticket: `leafblower-5k08.0`)
Add LM-damping subsection to `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md`:
- Append "Levenberg-Marquardt damping" subsection under "Single-Pass Per-Step Algorithm".
- Update line 76 ("Newton step: δ = H⁻¹ ∇g") to reference the damped form.
- Document gain-ratio μ schedule.
- Commit: `docs(spec): Newton-KL — add LM scale-invariant damping section`.
**Sequenced first** (atomic, code-free).

### LM-1 (`leafblower-5k08.1`) — Implement LM solver
- Edit `src/newton_calib.cpp`:
  - Rename `mu` → `lm_mu` (and `mu_min/mu_max` → `lm_mu_min/lm_mu_max`).
  - Init `lm_mu = 1.0` (not 1e-3).
  - Compute `d_floor = mean(diag(H))` after Schur complement (line ~163).
  - Damp: `H[a*n+a] = max(H[a,a]*(1+lm_mu), lm_mu*d_floor)` (additive floor).
  - LDLT pivot stays at 1e-12 (damping by construction sufficient).
  - Replace trust-region clip block entirely.
  - Compute `predicted = α·(G·δ) - 0.5·α²·(δᵀHδ)` and `ρ = (g_curr - g_trial) / predicted`.
  - Three-way Armijo outcome:
    * `FULL_ACCEPT && ρ>0.75` ⇒ `lm_mu /= 3`.
    * `accepted` (any α<1, or α=1 with 0.25<ρ<0.75) ⇒ `lm_mu` unchanged.
    * `ρ<0.25 || !accepted` ⇒ `lm_mu *= 10`, redo LDLT solve in same iter (max 3 retries; if still failing, return `RK_ERR_NOCONV`).
- Edit `src/newton_calib.hpp`: add `double lm_mu_final = 0.0;` field to `NewtonCalibResult`.
- Set `res.lm_mu_final = lm_mu;` before each return path in `newton_calib.cpp`.
- Verify build: `R CMD INSTALL --preclean . 2>&1 | grep -iE "warning|error"` ⇒ no new warnings.

### LM-1b (`leafblower-5k08.1b`) — μ-schedule unit test (new ticket)
- Add T6 to `tests/testthat/test-newton-kl.R`: synthetic 2-iter problem, expose `lm_mu_final`, assert ρ-driven decrease on accept-with-good-prediction; assert increase on FAILED line search.
- Amend T1 to assert exit reason is gradient-based (not `step_norm < 1e-14`).

### LM-2 (`leafblower-5k08.2`) — Test suite
- `Rscript -e "testthat::test_file('tests/testthat/test-newton-kl.R')"` ⇒ T1-T6 all PASS.
- `Rscript -e "testthat::test_local('.')"` ⇒ FAIL=0.
- Verdict: PASS → LM-3a; FAIL on T2 → NEEDS_HOMOTOPY (halt, draft Epic-B); other failure → STOP, escalate.

### LM-3a (`leafblower-5k08.3a`) — Run benchmark
- Pre-step: re-baseline ieppa+sraa on kk1204 severe at HEAD (verify spec-table 3.7s/2.4e-14 still reproduces). If ieppa baseline shifted, capture and compare against new baseline.
- `OMP_NUM_THREADS=1 Rscript benchmarks/newton_kl_bench.R 2>&1 | tail -20`.
- Capture severe + moderate skew rows in `benchmarks/results/newton_kl_kk1204.csv`.
- bench::mark with `iterations=2`, report median.

### LM-3b (`leafblower-5k08.3b`) — Verdict + commit + ESCAPE safety
Decision rules (using **<3s honest gate**, not <2s):
- **GATE_MET**: severe wall<3s ∧ max_err<1e-4 ⇒ write `docs/investigations/2026-05-01-newton-kl-lm-result.md` + `NEWS.md` entry + commit `feat(newton_kl): LM-damped Newton — kk1204 severe-skew GATE_MET`. Close epic.
- **NEEDS_HOMOTOPY**: severe converges (`lm_mu_final` bounded, `dual_gap` decreased monotonically) but misses gate ⇒ write investigation doc + `NEWS.md` partial-result entry + commit `fix(newton_kl): LM partial — NEEDS_HOMOTOPY`. Draft Epic-B stub ticket (target homotopy `T_τ = (1-ε)T + ε·uniform`, ε: 0.5→0.1→0.01→0; reset `lm_mu` at each ε step).
- **ESCAPE**: severe diverges or `lm_mu_final` saturated at 1e12 ⇒ **immediate AUTO-routing safety patch** in `src/c_api.cpp` (line ~171) and `src/r_bridge.cpp` (line ~425): the `K≥5 && M_cell/n≥0.9 → newton_kl` rule must be guarded with an explicit detection (e.g., target-skew metric `max_target/min_target > 5`) or fallback, ensuring `harvest()` calls do not silently route to a diverging solver. Commit `fix(newton_kl,auto): LM ESCAPE — AUTO routing safety guard`. Draft Epic-D stub ticket.

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
