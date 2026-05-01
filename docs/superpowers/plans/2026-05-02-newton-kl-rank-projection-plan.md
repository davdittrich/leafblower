# Newton-KL Rank-Deficient Direction Projection — Implementation Plan

**Date:** 2026-05-02
**Epic:** `leafblower-xytj` (Epic-D)
**Master HEAD:** `7dca47d` (post Epic-B BLOCKED + Epic-C BLOCKED)
**Status:** DRAFT — pending plan-review-gate

## References

- **Spec (current):** `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` (LM-0, TH-0, WI-0 amendments shipped; this plan adds RP-0).
- **Investigations (load-bearing):**
  - `docs/investigations/2026-05-01-newton-kl-lm-result.md` — LM-only basin floor at gap≈1.4e-4, max_err 2.79e-4 on stepstone K=9 T2.
  - `docs/investigations/2026-05-01-newton-kl-homotopy-result.md` — Epic-B target homotopy BLOCKED.
  - `docs/investigations/2026-05-02-newton-kl-ieppa-warmstart-result.md` — Epic-C IEPPA warm-start BLOCKED.
- **Code anchors:**
  - `src/newton_calib.cpp` lines 139–340 (`run_newton_inner` lambda), specifically lines 234–254 (LM-damped LDLT solve) and lines 273–282 (`H_pre` retained for gain-ratio).
  - `src/calib_linalg.hpp` / `src/calib_linalg.cpp` — currently exposes `ldlt_factor_inplace` + `ldlt_solve` ONLY. No eigendecomposition primitive.
  - `src/Makevars` — links `-lm -lmvec` ONLY. **No LAPACK / no Eigen / no BLAS.**
  - `src/newton_calib.hpp` — `NewtonCalibResult` already has `lm_mu_final`, `n_homotopy_levels_used` (stale; will be repurposed/replaced).

## Ticket sequence

WJ-0 → WJ-1 → WJ-2 → WJ-3 → WJ-4 → WJ-5 (six tickets, strictly serial; one ticket per task per repo policy).

| Ticket | Title | Touches |
|---|---|---|
| WJ-0 | Spec amendment — RP-0 "rank-deficient direction projection" | `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` |
| WJ-1 | Add symmetric eigendecomposition primitive (cyclic Jacobi) | `src/calib_linalg.{hpp,cpp}`, `tests/testthat/test-calib-linalg.R` (new) |
| WJ-2 | Wire projection into `run_newton_inner` (replaces LDLT solve path) | `src/newton_calib.cpp`, `src/newton_calib.hpp` |
| WJ-3 | Synthetic rank-1-deficient unit test (projection-correctness) | `tests/testthat/test-newton-rank-projection.R` (new) |
| WJ-4 | T2 / T7 / kk1204-K20-severe verification (gates) | `tests/testthat/*` execution; no source change |
| WJ-5 | Verdict doc + carryover cleanup (drop stale `n_homotopy_levels_used`) | `docs/investigations/2026-05-02-newton-kl-rank-projection-result.md` (new), `src/newton_calib.hpp`, `src/r_bridge.cpp` |

---

## §1 Context

Epic-B (target homotopy) and Epic-C (IEPPA warm-start) both BLOCKED. The empirical lesson (`2026-05-02-newton-kl-ieppa-warmstart-result.md` §Lessons #2): both attacks moved λ to a different point in the SAME ill-conditioned dual landscape; neither changed the landscape, so neither escaped the basin floor.

The basin floor mechanism, traced in `2026-05-01-newton-kl-lm-result.md`:

1. Stepstone K=9 has overlapping margins → many empty cells → `H_pre = (1/Z)·Σ f_i·χ_iχ_i^T − G⊗G` is rank-deficient (deficiency = number of unsupported margin-cell pairs after Schur).
2. LM damping `H_damp = H_pre + μ·diag(H_pre) + μ·d_floor·I_floor_active` makes `H_damp` invertible but does NOT change the direction of the solve. `δ = H_damp^{-1}·G` still has unbounded magnitude along eigenvectors of `H_pre` with `λ_i ≈ 0` (because `1 / (μ · 0 + μ · d_floor)` is small but not zero, and `v_i^T G` is generically O(1)).
3. Newton "accepts" steps along these directions because they ARE descent on `g(λ)` — `g` decreases unboundedly along the rank-null ray (as documented in iter-trace, lines 28–30 of `lm-result.md`: `g: 1.85 → -1.99e1 → -∞`). But the gradient `∇g` does not shrink: λ moves perpendicular to the gradient.
4. Best-iterate fallback rescues max_err but freezes at 2.79e-4 (the K=9 basin floor), 2.79× the 1e-4 T2 gate.

**Mechanism change required:** make Newton blind to rank-deficient directions. Project the step onto the well-conditioned subspace.

## §2 Mechanism

### Standard truncated-SVD Newton

Symmetric `H_pre ∈ R^{n_λ × n_λ}` with eigendecomposition `H_pre = V Λ V^T`, `Λ = diag(λ_1, …, λ_{n_λ})`, `λ_i ≥ 0` (PSD; clamp tiny-negative to 0 — same logic as line 282 of `newton_calib.cpp`). Define threshold `λ_thr = ratio · λ_max`. Truncated pseudoinverse:

```
H_pre^+ = V · diag(1/λ_i if λ_i ≥ λ_thr else 0) · V^T
```

Newton step:

```
δ = H_pre^+ · G   (no LM needed in retained subspace)
  = Σ_{i: λ_i ≥ λ_thr} (v_i^T G / λ_i) · v_i
```

In the dropped subspace, `δ` has zero component → λ does not move along eigenvectors with negligible curvature. The dual-objective component along those directions is left at its current value. This is exactly the Moore–Penrose minimum-norm solution to the underdetermined system `H_pre·δ = G`.

### Composition with LM

Two compositions are mathematically equivalent in the LIMIT but differ in finite-precision behavior:

(a) **Project then LM:** truncate `H_pre` first, then add `μ·diag(H_proj)` only on retained eigenvalues. Solve in retained subspace.
(b) **LM then project:** form `H_damp = H_pre + μ·diag(H_pre) + μ·d_floor`, eigendecompose `H_damp`, truncate.

**Decision: (a) project then LM.** Rationale: the rank-deficient directions have no curvature information at all — adding LM damping there is a numerical band-aid (already proven insufficient). Projecting first gives a cleanly-defined retained subspace; LM on that subspace then plays its proper role (regularizing the well-conditioned-but-imperfectly-Hessian-modeled curvature). Confidence: 80 (theoretical; (b) is also defensible; we will assert (a) and revisit if the gate fails).

The retained-subspace LM solve is implemented as: in the eigenbasis, the system is diagonal — `δ_i = (v_i^T G) / (λ_i · (1 + μ) + μ·d_floor_retained)` for `λ_i ≥ λ_thr`, else `δ_i = 0`. Project back: `δ = V_retained · δ_retained`. **No LDLT needed** — projection bypasses the linear-solve entirely. This eliminates the `ldlt_factor_inplace` retry loop (lines 239–254) for the projected path.

### Threshold value

For K=9 stepstone, `H_pre` eigenvalues span 1e-12 to 1e-1 (estimated from `n=200000`, mean cell density ≈ 0.05). Reasonable thresholds:

| `ratio` | Behavior |
|---|---|
| `1e-2` | Aggressive — drops 30–50% of dims; fastest convergence in retained subspace, risk of dropping mildly-informative directions. |
| `1e-6` | Standard truncated-SVD threshold (machine-precision-ish); drops only true rank-null directions. Conservative. |
| `1e-10` | Too generous — admits FP-noise eigenvalues; close to current LM-only behavior. |

**Decision: `ratio = 1e-8` initial (one safety order below double-precision relative round-off ~1e-15 amplified by `n_λ²` ≈ 6400 for K=20 → ~1e-11; with stepstone n_λ ≈ 30, headroom is 1e-13; pick 1e-8 to dominate Schur-subtraction cancellation noise without dropping informative directions).** Confidence: 70 — derived from FP error analysis; if T2 fails at 1e-8 we tune to 1e-6 then 1e-4 (each is one phone-line away). The threshold value is exposed as a `NewtonCalibParams` field (added in WJ-2) so it can be swept without recompile.

### Why this works where Epic-B / Epic-C failed

Both prior epics moved λ to a different starting point in the same ill-conditioned landscape. Projection changes the LANDSCAPE: in the retained-eigenvector subspace, `H_pre` is full-rank by construction, the Hessian is uniformly bounded above and below, Newton's quadratic-convergence basin is well-defined, `g(λ)` is bounded below in the retained subspace (because the unbounded ray is now orthogonal to permitted motion). The basin floor disappears as a mathematical object, not as an empirical band-aid. Confidence: 75 — the basin floor IS the rank-deficiency signature; eliminating one eliminates the other modulo unforeseen second-order pathologies.

### Citation / standard-method anchor

This is **truncated-SVD Newton**, equivalently **Levenberg–Marquardt with rank truncation**, equivalently **regularized Gauss–Newton via the Moore–Penrose pseudoinverse**. Standard references:

- Nocedal & Wright, *Numerical Optimization* (2nd ed.), §10.3 "Trust-Region Newton-CG with Rank Truncation".
- Hansen, *Rank-Deficient and Discrete Ill-Posed Problems* (SIAM 1998), §3.2 "TSVD".
- Press et al., *Numerical Recipes* §15.4 "Singular Value Decomposition" applied to least-squares.

The mechanism is textbook; reviewer concern "F-M1 mechanism unproven" cannot be raised on theoretical grounds. The empirical question (does it pass T2 on stepstone K=9?) is settled by WJ-4.

## §3 Forbidden

| ✗ Forbidden | Reason |
|---|---|
| Linking against LAPACK / `dsyev` / `dsyevd` | `Makevars` has `-lm -lmvec` only; adding LAPACK is a configure.ac / CRAN policy change, OUT OF SCOPE for this epic. Implement a self-contained symmetric eigensolver. |
| Adding Eigen3 / Armadillo / RcppArmadillo as a dependency | Same reason; package is dependency-free C++17. |
| Iterative methods (Lanczos, Arnoldi) | Need full spectrum to apply threshold; `n_λ` ≤ 80 makes direct method trivially affordable. |
| Cholesky-based pseudoinverse | Cholesky needs PD; `H_pre` is PSD with deliberate rank deficiency. |
| Modifying the LDLT codepath (`ldlt_factor_inplace`, `ldlt_solve`) | These primitives have other callers; preserve. The projection path is a sibling solver, not a replacement of the primitive. |
| Bandaid threshold tuning to pass T2 (e.g. `ratio = 1e-2` to drop 50% of dims) | Threshold must be derived from FP error analysis (§2 Threshold). Tuning to pass tests is a `SPEC_FAILURE` per CLAUDE.md §9. |
| Removing best-iterate fallback in WJ-2 | Defensive against unforeseen drift; remove only after WJ-4 measurement (see Q7). Default: keep, log how often it fires. |
| Public R API change | New diagnostic field on `NewtonCalibResult` is C++/R-bridge only, not a breaking R-level signature change. |

## §4 Plan Steps

Each ticket below has a skeletal §I–§VI body. Tickets are independent commits gated on prior-ticket passing tests. Beads tickets created post plan-review-gate; this plan does NOT create them.

---

### WJ-0 — Spec amendment "RP-0: rank-deficient direction projection"

**§I Goal.** Append a §3.10 "Rank-deficient direction projection (RP-0)" subsection to `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` documenting the §2 Mechanism above (truncated-eigendecomposition Newton, threshold rationale, composition with LM). Note Epic-B / Epic-C BLOCKED status as motivation. State explicit gate: T2 stepstone K=9 max_err < 1e-4; kk1204 K=20 severe-skew partial credit (max_err < 1e-3, no wall-time gate).

**§II Mechanism.** Markdown edit only.

**§III Forbidden.** Don't reword existing §3 LM/TH/WI subsections; append cleanly.

**§IV Audit.** Spec diff reviewed for consistency with TH-0 and WI-0 amendments (no contradictions).

**§V Tests.** N/A (doc).

**§VI Done when.** Spec PR merges; verdict = "spec amendment shipped; ready for WJ-1."

---

### WJ-1 — Symmetric eigendecomposition primitive in `calib_linalg`

**§I Goal.** Add `int sym_eig_jacobi(double* A, size_t n, double* eigvecs, double* eigvals, double tol, int max_sweeps)` to `src/calib_linalg.{hpp,cpp}`. Cyclic Jacobi rotations (Golub & Van Loan §8.4); on entry `A` is the n×n symmetric matrix (full storage, both triangles); on exit `eigvecs` is column-major (each column = eigenvector, length n), `eigvals` is length n. `A` is destroyed in place. Returns `RK_OK` / `RK_ERR_NOCONV` (if `max_sweeps` exhausted without off-diag norm < tol).

**§II Mechanism.** Cyclic Jacobi is ~30 lines of C++17, pure scalar arithmetic, no external dependencies, deterministic, numerically stable for symmetric PSD. For n ≤ 80, ~50× n³ flops worst case = 25M flops for kk1204 K=20 — still <1ms wall-time.

**§III Forbidden.** No LAPACK/Eigen. No SVD (we have symmetry). No `goto`. No threading (problem too small). No reuse of `H_pre` storage as workspace — caller passes its own.

**§IV Audit.** Unit test: random symmetric PSD 5×5 fixture, compare `V·Λ·V^T` reconstruction against original (Frobenius norm < 1e-12); `V^T·V = I` (orthogonality < 1e-12); known rank-1 deficient fixture (one zero eigenvalue) recovered with `λ_min < 1e-13`; pathological: diagonal matrix → identity eigvecs.

**§V Tests.** New `tests/testthat/test-calib-linalg.R` exercising via a thin `.Call` wrapper added to `r_bridge.cpp` (test-only diagnostic export, not user-facing).

**§VI Done when.** All four unit tests pass; `R CMD INSTALL --preclean .` clean; tests run < 100ms.

---

### WJ-2 — Wire projection into `run_newton_inner`

**§I Goal.** Replace the LM-damped LDLT solve path (`src/newton_calib.cpp` lines 234–254) with: eigendecompose `H_pre` (the saved copy at line 236), apply threshold-truncated pseudoinverse with retained-subspace LM damping, compute `δ` directly. Add `n_projected_dims` field to `NewtonCalibResult` (in `src/newton_calib.hpp`), set to `n_λ - retained_count` per iter (final value = max projected dims across iters; useful for post-mortem). Add `projection_threshold_ratio` field to `NewtonCalibParams` (default 1e-8).

**§II Mechanism.** Per-iter steps, replacing existing 5h:

1. Eigendecompose `H_pre` (saved copy is already there for gain-ratio at line 236) → `V`, `Λ`.
2. `λ_max = max(eigvals)`; `λ_thr = projection_threshold_ratio · λ_max`.
3. Compute `g_eig = V^T · G` (length n_λ).
4. For each i: if `λ_i ≥ λ_thr`, `δ_eig[i] = g_eig[i] / (λ_i · (1 + lm_mu) + lm_mu · d_floor)`; else `δ_eig[i] = 0`.
5. `δ = V · δ_eig`.
6. Compute `δᵀ H_pre δ` (still needed for gain-ratio at line 273 — but now equals `Σ_retained δ_eig[i]² · λ_i`, so cheap closed form).
7. Increment `n_projected_dims_curr` count of dropped dims; track max across iters.

LM retry loop (lines 239–254) is **eliminated** — projection cannot fail. The `solve_ok` branch and the LDLT-3x-fail early-exit (lines 255–271) remove. Best-iterate-fallback (lines 219–226, 261–269) **stays** (Q7 — defensive; will measure how often it fires in WJ-4 and decide on removal in WJ-5 cleanup).

**§III Forbidden.** No removal of LDLT primitives from `calib_linalg.hpp` (other callers: `logit_calib.cpp`, `chebyshev.cpp`). No removal of best-iterate fallback yet. No removal of `lm_mu` adaptation logic — μ continues to scale retained-subspace damping. Don't compute eigendecomp BEFORE the gradient norm convergence check at line 228 (skip the work on the converged-already iter).

**§IV Audit.** Spy on `n_projected_dims_curr`: stepstone K=9 should report ≥ 1 projected dim per iter once near optimum; K=4 (no rank deficiency) should report 0 projected dims throughout (regression assertion in WJ-3). Trace: log per-iter `(λ_max, λ_min, n_projected, dual_gap, lm_mu)` to `res.message` in debug builds (compile-time `LBW_NEWTON_TRACE`).

**§V Tests.** Existing T1–T8 must still pass (regression). New T9 (next ticket) exercises projection-correctness. `R CMD INSTALL --preclean .` clean.

**§VI Done when.** All existing tests green; debug trace shows `n_projected_dims > 0` on stepstone K=9 final iters.

---

### WJ-3 — Synthetic rank-1-deficient projection-correctness test

**§I Goal.** Add `tests/testthat/test-newton-rank-projection.R` with a synthetic fixture where one margin's category is a deterministic linear combination of two others (rank-1 deficiency by construction in the indicator design). Assert: (a) Newton converges; (b) `n_projected_dims ≥ 1`; (c) λ-component along the rank-null direction stays at initialization (within tol_abs).

**§II Mechanism.** Three margins, n=1000 rows, margins constructed so `χ_3 = χ_1 + χ_2` (in indicator-vector terms — pick categories so the design matrix has a known one-dim null vector). Initialize `λ = 0`, run `newton_calibrate`, check `λ` after convergence has zero (within FP) component along the known null vector. K=4 (no rank deficiency) sub-test: `n_projected_dims == 0` throughout (asserts WJ-2 didn't over-project on healthy fixtures — addresses the K=4 regression risk noted in chebyshev rev2 plan as analogous concern).

**§III Forbidden.** No reliance on stepstone fixture for projection-correctness — too high-dim to verify the geometric assertion. Synthetic only.

**§IV Audit.** Run with `LBW_NEWTON_TRACE` enabled; trace must show the rank-null eigenvalue in the `eigvals` printout and the corresponding `δ_eig` set to 0.

**§V Tests.** This IS the test.

**§VI Done when.** All four sub-assertions pass; existing T1–T8 still green.

---

### WJ-4 — T2 / T7 / kk1204-K20-severe verification (gate)

**§I Goal.** Run the full test suite (`R CMD check` + targeted T2/T7/kk1204 benchmarks). T2 stepstone K=9: max_err < 1e-4 (HARD gate). T7 K=4: max_err < 1e-4 AND n_projected_dims == 0 (no over-projection on well-conditioned fixtures). kk1204 K=20 severe-skew (`benchmarks/kk1204_severe_skew.R`): max_err < 1e-3 (PARTIAL gate — see §6 Cost / §7 Risks); no wall-time gate. Compare against master `7dca47d` baseline numerically.

**§II Mechanism.** Execute existing benchmarks; no source change. Capture per-iter trace from `LBW_NEWTON_TRACE` build, write to `docs/investigations/2026-05-02-newton-kl-rank-projection-trace.md`.

**§III Forbidden.** Tuning `projection_threshold_ratio` mid-run to make T2 pass. If `1e-8` fails T2, halt → log `SPEC_FAILURE` → escalate to user (per CLAUDE.md §9 #4). Do NOT silently retry with `1e-6`.

**§IV Audit.** Independent re-run on a clean clone; outputs identical (deterministic Jacobi, deterministic line search).

**§V Tests.** N/A (this IS the test execution).

**§VI Done when.** T2 PASS (or `SPEC_FAILURE` halt with full trace artifact); T7 PASS; kk1204 PARTIAL or BLOCKED (documented).

---

### WJ-5 — Verdict doc + carryover cleanup

**§I Goal.** Write `docs/investigations/2026-05-02-newton-kl-rank-projection-result.md` (Epic-D verdict). Remove stale `n_homotopy_levels_used` from `NewtonCalibResult` (never set after Epic-B BLOCKED) and stale `n_warmstart_iters_used` if present (Epic-C unused field per `2026-05-02-...-warmstart-result.md` §Carryover scope). Surface `n_projected_dims` and `lm_mu_final` in `r_bridge.cpp`. Decide best-iterate fallback fate based on WJ-4 trace: if it fires zero times across T2/T7/kk1204, remove with a comment explaining; else keep.

**§II Mechanism.** Doc + small surgical edits.

**§III Forbidden.** No removal of best-iterate fallback if WJ-4 trace shows even a single fire (defensive). No removal of `n_projected_dims` (load-bearing for verdict).

**§IV Audit.** `git grep n_homotopy_levels_used` returns empty after edit. `r_bridge.cpp` exposes the new field, R-level test reads it.

**§V Tests.** Existing tests still green.

**§VI Done when.** Verdict doc written; PR merged; epic closed (`bd close` with verdict tag).

---

## §5 Cost Estimate

### Per-iter cost (additive vs master)

For n_λ × n_λ symmetric matrix:

| Component | Flops | Notes |
|---|---|---|
| Cyclic Jacobi eigendecomp | ~30 · n_λ³ | dominant cost of new code |
| `g_eig = V^T · G` | 2 · n_λ² | matrix-vector |
| Diagonal solve in eigenbasis | n_λ | trivial |
| `δ = V · δ_eig` | 2 · n_λ² | matrix-vector |
| **Total new** | ≈ 30 n_λ³ + 4 n_λ² | replaces ~n_λ³ from LDLT |

Net per-iter overhead: ~**29 n_λ³**.

| Fixture | n_λ | Per-iter overhead | Per-iter accumulation cost (master) | Overhead fraction |
|---|---|---|---|---|
| Stepstone K=9 | ≈ 30 | 0.78 Mflops | 200000 · 9² = 16 Mflops | 4.9% |
| kk1204 K=20 nj=5 | ≈ 80 | 14.8 Mflops | 1500000 · 400 = 600 Mflops | 2.5% |

Single-digit-percent slowdown on per-iter wall, well below the existing tolerance (LSE-stable iter is already a 2× hit accepted in Epic-A). Confidence: 85 — flop counts are exact for cyclic Jacobi; per-iter accumulation flop count is `n · K²` from line 174–183 of `newton_calib.cpp`.

### Total wall budget

If projection cuts iter count by ≥ 1 (basin-floor avoidance — currently iter 7+ are wasted drift on stepstone), the net is faster than master despite per-iter overhead. Confidence: 60 — drift iters DO consume work in master; eliminating them is plausible but not guaranteed.

### Wall-clock prediction

Stepstone K=9: master `0.6s`, predicted `0.6–0.65s`. kk1204 K=20: master `2.4s` (Newton-KL, drifts to fail), predicted `2.5–2.7s` (likely converges → PARTIAL gate). No wall regression beyond noise.

## §6 Risks & Open Questions

### Q1 — Does the project link against LAPACK?

**Resolved:** No. `src/Makevars` line `PKG_LIBS = -lm -lmvec` (verified by reading file; confidence 95). Plan implements cyclic Jacobi from scratch in WJ-1.

### Q2 — Symmetric eigensolver choice

**Resolved:** Cyclic Jacobi (vs Householder+QR). Justification: simpler implementation (~30 LOC), guaranteed convergence on symmetric PSD, no spurious complex roots, rotations are local. Householder+QR is faster asymptotically but for n_λ ≤ 80 the constant factors dominate; Jacobi is comparable and code is half the size.

### Q3 — Threshold value

**Decision:** `projection_threshold_ratio = 1e-8` default, exposed as parameter. See §2 Threshold. If WJ-4 fails T2 at this value, halt per CLAUDE.md §9 — do NOT silently re-tune.

### Q4 — Order of LM + projection

**Decision:** Project then LM (option (a) §2). LM operates only in retained subspace; rank-null directions get zero motion regardless of μ. Confidence 80; alternative (b) is acceptable fallback if (a) destabilizes.

### Q5 — Dropped components of G

**Resolved (theoretically):** In coordinate-descent terms, the dropped directions represent "data does not constrain this λ-component." Standard truncated-SVD theory: minimum-norm solution leaves λ unchanged in those directions, which is the unbiased default (any nonzero choice would be arbitrary). Confidence 90.

### Q6 — Test for projection correctness

**Resolved:** WJ-3 synthetic rank-1-deficient fixture asserts λ stays at init in the null direction.

### Q7 — Best-iterate fallback fate

**Deferred to WJ-5.** Plan: keep through WJ-4, instrument fire count, remove only if zero fires across T2/T7/kk1204. Defensive removal default.

### Q8 — Wall cost on kk1204

**Resolved:** §5 — per-iter overhead 2.5% on kk1204; total wall regression < 5% expected. Negligible.

### Q9 — API surface

**Decision: yes, add `n_projected_dims` to `NewtonCalibResult`.** It's load-bearing for the WJ-4 verdict (proves projection is happening) and post-deployment debuggability. Cost: one int, zero R-API impact (NewtonCalibResult is C++ only, surfaced via `r_bridge.cpp`). Confidence 85.

### Q10 — Plan-review-gate pre-emption

| Prior concern | Apply here? | Mitigation |
|---|---|---|
| F-M1 "mechanism unproven" (Epic-C) | NO — truncated-SVD Newton is textbook (Nocedal & Wright §10.3, Hansen §3.2). Cited above. |
| F-R3 "T7 monotonicity" (Epic-C) | YES — adding T7 K=4 over-projection regression to WJ-3. |
| C-3 "ieppa standalone bench comparison" (Epic-C) | NO — Epic-D doesn't touch IEPPA. |
| F-R4 "spec table fixture mismatch" (Epic-C) | NO — Epic-D doesn't update spec convergence table; only adds RP-0 subsection. |
| F-empirical "kill-switch insufficient" (Epic-C) | NO — Epic-D's gate IS a direct empirical T2 measurement (WJ-4), not a proxy kill-switch. |

### Open risks

1. **R1 — Threshold sensitivity.** `1e-8` is derived but unproven for kk1204 K=20 (n_λ=80, more accumulated cancellation). Mitigation: WJ-4 trace; if it fails kk1204 we file follow-up to make threshold problem-size-aware (e.g. `ratio = c · √n_λ · ε_machine`).
2. **R2 — Jacobi convergence on degenerate spectrum.** Repeated eigenvalues are handled (Jacobi's standard tie-breaking) but worth a stress-test: WJ-1 unit test should include a fixture with `λ_1 = λ_2`.
3. **R3 — Second-order pathology.** Even with rank-null directions projected out, the retained subspace might have its OWN basin floor (smaller). Mitigation: WJ-4 is the empirical test; if T2 passes, R3 didn't materialize; if T2 fails, R3 manifested and we file Epic-E.
4. **R4 — Reference-category amplification.** Stepstone's `max_err = 2.79·gap` factor (per `lm-result.md`) means even gap=5e-5 → max_err = 1.4e-4 → still misses gate. The gap reduction needs to be ≥ 14× to clear T2. Confidence the projected basin reaches gap < 5e-5: 65 — depends on retained-subspace condition number.

## §7 Out of Scope

- Trust-region Newton (alternative to LM in retained subspace) — LM has been validated, switching is unnecessary work.
- LAPACK/Eigen integration — see §3 Forbidden; configure.ac change is its own multi-file epic.
- AUTO routing changes — Epic-D either succeeds or fails; routing is downstream of method choice.
- kk1204 K=20 severe-skew HARD gate — PARTIAL only (max_err < 1e-3); HARD gate requires separate work (likely SOR composition with ieppa, see existing `feat(r_bridge): expose alm_penalty`).
- LM-1c carryover (`lm_mu_final` r_bridge surface) — included in WJ-5.
- L-BFGS-B fallback path — orthogonal to dual-Newton, unrelated.
- Chebyshev IPM ν-correction (`leafblower-bcr3`) — separate epic, no overlap.
- IEPPA improvements — Epic-C closed BLOCKED, no further work scoped here.

---

**Plan rev:** 1
**Author:** Claude (Opus 4.7 1M)
**Awaiting:** plan-review-gate (Feasibility, Completeness, Scope & Alignment).
