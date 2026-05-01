# Epic-Dβ: Truncated-SVD + Steihaug-CG Trust-Region Newton (LAPACK-backed)

## 1. Header

- **Epic ID:** Epic-Dβ
- **Plan ID:** `leafblower-wkmq` / Epic-Dβ rev1
- **Date:** 2026-05-02
- **Author:** controller (Opus, planning subagent)
- **Refs:**
  - Epic-B BLOCKED (KL-Newton μ-floor sweep, plateau 2.79e-04)
  - Epic-C BLOCKED (homotopy continuation, plateau 2.79e-04)
  - Epic-Dα ESCALATE (LM pivot sweep, plateau 2.79e-04)
  - Master HEAD: `benchmarks/results/lm_pivot_sweep.csv` (Epic-Dα CSV committed)
  - Source: `src/newton_calib.cpp` (`run_newton_inner` lambda, lines ~190-225 LM damping; LDLT line ~248)
  - Build: `src/Makevars` currently `-lm -lmvec` only — **LAPACK NOT linked**
- **Ticket sequence:** WL-0 → WL-1 → WL-2 → WL-3 → WL-4 → WL-5 → WL-6 (7 atomic tickets, no bundling per project policy)

## 2. Context

Three consecutive empirical investigations have failed to move the KL-Newton residual gap below **2.79e-04** on the K=4 / K=20 ill-conditioned stepstone:

| Epic | Hypothesis | Sweep | Result | CSV |
|------|------------|-------|--------|-----|
| B | μ_min too low → spurious damping starvation | `lm_mu_min ∈ {1e-12, 1e-10, 1e-8}` | Plateau 2.79e-04 across all | (Epic-B benchmark, master) |
| C | Continuation needed — homotopy on conditioning | K=4 stepstone with continuation steps | Plateau 2.79e-04 | (Epic-C benchmark, master) |
| Dα | LDLT pivot threshold too aggressive | `pivot ∈ {1e-12, 1e-10, 1e-8}` | Plateau 2.79e-04 across all | `benchmarks/results/lm_pivot_sweep.csv` |

**Conclusion:** the bottleneck is **structural Newton-step direction**, not damping magnitude or factorization conditioning. The LDLT step is solving the wrong problem: in unidentified directions of `H_pre`, gradient components are amplified inversely to vanishing eigenvalues, producing a step that is locally valid (LM line search accepts it) but globally orthogonal to the true descent direction in the identified subspace.

User has explicitly authorized: **composed pipeline truncated-SVD + Steihaug-CG trust-region Newton, with LAPACK linkage explicitly OK.**

## 3. Mechanism

### Target pattern
Composed three-stage pipeline replacing the LDLT solve in `run_newton_inner`:

1. **Eigendecompose** `H_pre = V Λ V^T` via LAPACK `dsyevd_` (divide-and-conquer, all eigenvalues).
2. **Truncated-SVD projection.** Define retained set `R = {i : λ_i ≥ ratio · λ_max}` with `ratio = 1e-8` (initial, configurable).
3. **Steihaug-CG trust-region step in retained subspace** (with adaptive Δ via Marquardt-style ρ ratio).
4. **Compose with LM damping in the eigenbasis** (choice (b) per question 6).

### Algorithm (per Newton inner iteration)

```
# Inputs: H_pre (n_λ × n_λ symmetric PSD), G (n_λ gradient), μ (LM damping), Δ (trust radius)
# Constants: ratio_tsvd = 1e-8 (default), tol_cg = 1e-10, max_cg_iter = 2 * |R|

1. H_copy ← H_pre                                    # copy, dsyevd overwrites
2. dsyevd_(H_copy, Λ[1..n_λ], V ← H_copy)            # V is column-eigvecs, Λ ascending
3. λ_max ← Λ[n_λ];  thresh ← ratio_tsvd · λ_max
4. R ← {i : Λ[i] ≥ thresh};  n_keep ← |R|
   IF n_keep == 0: FALLBACK to LDLT (numerical safety)
5. # Project gradient into retained eigenbasis
   FOR i IN R:  g_keep[i] ← Σ_j V[j,i] · G[j]
6. # Compose LM damping in eigenbasis (choice b)
   d_floor ← mean(Λ[R])
   FOR i IN R:  Λ_damped[i] ← Λ[i] · (1 + μ) + μ · d_floor
7. # Decide: direct pseudoinverse vs Steihaug-CG
   FOR i IN R:  δ_keep_pinv[i] ← g_keep[i] / Λ_damped[i]
   IF ‖δ_keep_pinv‖₂ ≤ Δ:
       δ_keep ← δ_keep_pinv                          # trust radius not binding
   ELSE:
       δ_keep ← STEIHAUG_CG(Λ_damped, g_keep, Δ)     # see below
8. # Project step back to original basis
   FOR j IN 1..n_λ:  δ[j] ← Σ_{i ∈ R} V[j,i] · δ_keep[i]
9. # Diagnostic
   result.n_projected_dims ← n_λ - n_keep
```

### Steihaug-CG subroutine (operating on diagonal `Λ_damped` in eigenbasis)

Because the projected Hessian is diagonal `diag(Λ_damped)` with all entries strictly positive, there is **no negative curvature**. Steihaug-CG reduces to standard CG with trust-radius truncation. Body ~30-40 LOC. Exit conditions: (a) ‖p‖₂ exceeds Δ → solve quadratic for τ along last direction, return boundary point; (b) ‖r‖₂ < tol_cg; (c) iter ≥ max_cg_iter. Termination at boundary uses standard formula τ such that ‖p_{k} + τ d_k‖₂ = Δ.

### LM-composition choice (b) — justification (question 6)

| Option | Mechanism | Drawback |
|--------|-----------|----------|
| (a) Damp before eigendecomp | `H_damped = H_pre + μ·diag(H_pre)`; eigendecomp damped matrix | Eigenvalues of damped matrix ≠ eigenvalues of `H_pre`; truncation threshold semantics drifts with μ |
| (b) **Damp in eigenbasis** | Eigendecomp `H_pre`; project; `Λ_damped[i] = Λ[i]·(1+μ) + μ·d_floor` | None — preserves "LM scales eigenvalues" semantic, `d_floor` natural as mean of retained spectrum |
| (c) Damp full diagonal first, then eigendecomp+truncate | Same eigvec structure as (a) | Same issue as (a) |

**Decision: (b).** Confidence 90 — derives from standard truncated-SVD-pseudoinverse formula `(VΛV^T + μI)^+ = V (Λ + μI)^+ V^T` with the per-eigenvalue floor `d_floor` substituted for `I` to preserve scale-invariance. Floor restricted to retained directions only (`d_floor_retained = mean(Λ_keep)`).

### Adaptive trust radius (Marquardt-style ρ ratio, question 5)

```
ρ ← (f(λ) - f(λ + δ)) / m_predicted_decrease
IF ρ > 0.75:    Δ ← min(2 · Δ, 100)
IF ρ < 0.25:    Δ ← Δ / 4
# else unchanged
Δ_init ← 1.0
```

This piggybacks on the existing LM ρ-ratio acceptance test in `run_newton_inner` (no new code path).

## 4. Forbidden

- **NO** AUTO-routing changes. Truncated-SVD-CG is invoked from `run_newton_inner` only when `algo == "newton"` (or AUTO selects newton). Document explicitly: AUTO routing impact = none.
- **NO** Eigen3 dense eigendecomposition (introduces a heavy dep, 2-5x slower than LAPACK `dsyevd` per published benchmarks). Use LAPACK directly.
- **NO** μ_min / pivot retuning workarounds — refuted by `benchmarks/results/lm_pivot_sweep.csv`. Plateau is structural.
- **NO** removal of LM damping (it composes with truncated-SVD, does not replace it).
- **NO** removal of best-iterate fallback in this epic (defensive backstop; revisit in WL-6 if it never fires).
- **NO** algorithmic substitution without interleaved before/after benchmark within a single `bench::mark()` call (per `~/.agents/skills/CLAUDE.md` §2).
- **NO** SPEC pivot. If truncated-SVD-CG fails T2 hard gate, output `SPEC_FAILURE` and halt; do not silently bump ratio or fall back to LM-only.
- **NO** Eigen3 `selfadjointEigenSolver` even if present transitively via RcppEigen — call LAPACK directly to control link surface.

## 5. Audit Strategy

- **Spies:** Add a per-iteration trace vector `proj_trace` in `NewtonCalibResult` (debug build only, gated by `LEAFBLOWER_NEWTON_TRACE` macro). Each entry: `{iter, n_keep, lambda_min_kept, lambda_max, mu, delta, cg_iters_used, rho_ratio, exit_reason}`.
- **Mock strategy:** Unit test injects synthetic `H_pre` with known rank deficiency (e.g. block diag with one zero eigenvalue). Verify `n_projected_dims == 1` and δ has zero component in null direction (within 1e-12).
- **Test how, not if:** Tests assert on (a) `n_projected_dims` per benchmark scenario, (b) `cg_iters_used` distribution (should be ≪ n_keep when well-conditioned in retained subspace), (c) trust radius adaptation (Δ should grow on K=4 well-conditioned, contract on kk1204 K=20).

## 6. Plan Steps (atomic tickets)

Each ticket gets its own bead per project policy (no bundling). Ticket bodies follow hermetic §I-VI skeleton.

---

### **WL-0: LAPACK linkage + dsyevd_ wrapper header**

**§I Goal.** Link LAPACK into the package, add an `extern "C"` declaration for `dsyevd_` accessible to `newton_calib.cpp`. Confidence 95.

**§II Mechanism.**
- Edit `src/Makevars`: change `PKG_LIBS = -lm -lmvec` → `PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) -lm -lmvec`.
- Edit `src/Makevars.in` identically (template kept in sync per R-pkg idiom).
- Create `src/lapack_wrappers.hpp` with:
  ```cpp
  #pragma once
  extern "C" {
    void dsyevd_(const char* JOBZ, const char* UPLO, const int* N,
                 double* A, const int* LDA, double* W,
                 double* WORK, const int* LWORK,
                 int* IWORK, const int* LIWORK, int* INFO);
  }
  ```
- **Audit:** dummy compile-only smoke test in `src/newton_calib.cpp` (call `dsyevd_` on a 2×2 identity, assert eigenvalues = {1,1}, gated `#ifdef LEAFBLOWER_LAPACK_SMOKE`).

**§III Acceptance.** `R CMD INSTALL --preclean .` succeeds (mandatory R compile gate per `~/.agents/skills/CLAUDE.md` §3). `nm src/leafblower.so | grep dsyevd_` returns a defined symbol (or runtime resolved via R's libRblas).

**§IV Risks.** R's LAPACK macro names. **Mitigation:** WebSearch "R-exts LAPACK_LIBS Makevars" before edit to confirm the exact variable names and ordering (`$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)`). Per "Writing R Extensions" §1.2.1 these are correct.

**§V Forbidden.** Hardcoding `-llapack -lblas`. Linking against system LAPACK directly (must use R's resolution to match R's BLAS, including OpenBLAS / MKL / Accelerate variants). Adding RcppEigen.

**§VI Files.** `src/Makevars`, `src/Makevars.in`, `src/lapack_wrappers.hpp` (new).

---

### **WL-1: Eigendecomposition + truncated-SVD projection (no CG, direct pseudoinverse)**

**§I Goal.** Replace the LDLT solve in `run_newton_inner` with eigendecomp + truncated-SVD-pseudoinverse + LM in eigenbasis. No trust region yet (assume Δ = ∞). Confidence 85.

**§II Mechanism.** Implement steps 1-6 and 8 of the algorithm above; use `δ_keep_pinv` directly for `δ_keep`. Add `n_projected_dims` field to `NewtonCalibResult` struct. Default `ratio_tsvd = 1e-8` exposed as `newton_tsvd_ratio` parameter (R-side default).

**§III Acceptance.**
- T0 unit test: synthetic rank-deficient `H_pre` → `n_projected_dims > 0`, δ has zero component in null direction (1e-12).
- T1 baseline: K=4 well-conditioned bench unchanged within 1e-6 (no truncation should fire when spectrum is uniform).
- Compile + tests pass.

**§IV Risks.** Eigvec sign convention from `dsyevd`. **Mitigation:** signs cancel in `V·δ_keep` since they appear in both projection and back-projection.

**§V Forbidden.** Steihaug-CG (deferred to WL-2). Ratio retuning. Modifying LM ρ-ratio acceptance.

**§VI Files.** `src/newton_calib.cpp`, `src/newton_calib.hpp` (struct field), `R/newton.R` (parameter plumbing), `tests/testthat/test-newton-tsvd-projection.R` (new).

---

### **WL-2: T2 hard gate — stepstone <1e-4**

**§I Goal.** Run K=4 stepstone benchmark with WL-1 implementation. **HARD gate:** residual gap < **1e-4** (down from 2.79e-04 plateau). Confidence 70.

**§II Mechanism.** Reuse `benchmarks/lm_pivot_sweep.R` harness, swap algorithm flag. Single-row CSV `benchmarks/results/tsvd_pinv_T2.csv` with `gap`, `n_projected_dims`, `n_iter`, `wall_time_ms`.

**§III Acceptance.** **HARD:** `gap < 1e-4`. **HALT criterion:** if gap ≥ 1e-4, output `SPEC_FAILURE` and stop. Do NOT pivot to ratio retuning. Re-evaluate plan.

**§IV Risks.** Truncation alone may not be sufficient — CG might be needed even without trust radius binding (e.g. if d_floor injection is too aggressive or too loose). **Mitigation:** if ratio at default 1e-8 fails, allow ONE escalation to 1e-6 or 1e-10 with explicit user approval, then halt.

**§V Forbidden.** Tuning ratio_tsvd silently to pass gate. Disabling LM. Increasing `max_iter` beyond master baseline.

**§VI Files.** `benchmarks/tsvd_pinv_T2.R` (new), `benchmarks/results/tsvd_pinv_T2.csv` (new).

---

### **WL-3: Steihaug-CG trust-region in retained subspace**

**§I Goal.** Add Steihaug-CG path for the case `‖δ_keep_pinv‖₂ > Δ`. Adaptive Δ via existing LM ρ-ratio. Confidence 80.

**§II Mechanism.** Implement step 7 alternative branch + Marquardt-style Δ adaptation. CG body ~30-40 LOC, operating on diagonal `Λ_damped`. Initial `Δ = 1.0`. Trust radius adapted post-step: ρ>0.75 → Δ←min(2Δ,100); ρ<0.25 → Δ←Δ/4.

**§III Acceptance.**
- Unit test: synthetic `H` with `‖δ_pinv‖₂ > Δ` → CG returns boundary step with `‖δ‖₂ = Δ` (1e-10).
- T1 stepstone result no worse than WL-2 (gap ≤ WL-2 gap + 1e-6).
- Δ trace shows growth on well-conditioned, contraction on ill-conditioned.

**§IV Risks.** CG numerical drift on ill-scaled diagonal. **Mitigation:** since Λ_damped is diagonal and PSD, CG converges in ≤ n_keep iters analytically; tolerance test redundant but kept defensively.

**§V Forbidden.** Removing direct-pseudoinverse fast path (must keep for the common case `‖δ_pinv‖₂ ≤ Δ`). Curvature checks (provably unnecessary for diagonal PSD).

**§VI Files.** `src/newton_calib.cpp` (CG body), `tests/testthat/test-newton-steihaug-cg.R` (new).

---

### **WL-4: T7 hard gate — K=4 over-projection regression**

**§I Goal.** Verify on a well-conditioned K=4 problem (uniform spectrum, no rank deficiency) that **no eigenvalues are dropped** and the gap matches WL-2 baseline within 1e-4. Confidence 90.

**§II Mechanism.** Construct a K=4 well-conditioned scenario where `λ_max / λ_min < 1e6`. Run WL-3 implementation. Assert `n_projected_dims == 0` and `gap < 1e-4`.

**§III Acceptance.** **HARD:** `n_projected_dims == 0` AND `gap < 1e-4`. Failure → over-projection bug; halt and diagnose.

**§IV Risks.** Default ratio 1e-8 is too aggressive for borderline-conditioned problems. **Mitigation:** scenario is engineered well-conditioned; if it triggers, ratio default needs revision (escalate, do not silently retune).

**§V Forbidden.** Changing ratio_tsvd to make the gate pass. Loosening assertion to `n_projected_dims ≤ 1`.

**§VI Files.** `benchmarks/tsvd_T7_well_conditioned.R` (new), `benchmarks/results/tsvd_T7_K4.csv` (new).

---

### **WL-5: kk1204 K=20 severe-skew partial gate <1e-3**

**§I Goal.** Stretch goal — kk1204 K=20 ill-conditioned problem to gap < **1e-3** (PARTIAL gate, not HARD). Confidence 50.

**§II Mechanism.** Run WL-3 implementation on kk1204 K=20 stepstone. Single-row CSV result.

**§III Acceptance.** **PARTIAL:** `gap < 1e-3` desirable. Not a halt criterion. Documents reach of method.

**§IV Risks.** kk1204 may genuinely require continuation/multi-start; truncated-SVD-CG addresses direction, not initialization. **Mitigation:** failure documented as Out-of-Scope for Epic-Dβ; future Epic-E candidate.

**§V Forbidden.** Engineering a special initial point to pass this gate. Continuation steps (Epic-C territory).

**§VI Files.** `benchmarks/tsvd_kk1204_K20.R` (new), `benchmarks/results/tsvd_kk1204_K20.csv` (new).

---

### **WL-6: Best-iterate fallback dead-code safety check + cleanup ticket filing**

**§I Goal.** Verify whether best-iterate fallback ever fires under WL-3 implementation. If never fires across {T2, T7, kk1204, full master test suite}, file an Epic-Dβ-cleanup bead ticket (do NOT delete in this epic per `~/.agents/skills/CLAUDE.md` §2 — clean up only your own newly created orphans). Confidence 60.

**§II Mechanism.** Add a debug counter `best_iterate_fired_count` (gated by trace macro). Run the four scenarios above. Tabulate.

**§III Acceptance.** Decision made and recorded:
- If counter > 0 anywhere → **KEEP** best-iterate, document scenarios.
- If counter == 0 across all scenarios → **FILE** `leafblower-wkmq-cleanup-best-iterate` ticket for separate epic; do NOT remove in this PR.

**§IV Risks.** Counter macro discipline. **Mitigation:** trace gated under `LEAFBLOWER_NEWTON_TRACE`; release builds unaffected.

**§V Forbidden.** Deleting best-iterate code in this epic. Conflating WL-6 with refactoring of unrelated newton-calib code.

**§VI Files.** `src/newton_calib.cpp` (counter only), bead ticket file (created via `bd create`, not in this PR's diff).

---

## 7. Cost Estimate

| Ticket | LOC ±  | New files | Tests | Bench | Effort |
|--------|--------|-----------|-------|-------|--------|
| WL-0   | +20 / -2 | 1 hpp | 0 | 0 | S (Makevars + smoke) |
| WL-1   | +120 / -25 | 0 | 1 | 0 | M (eigendecomp + projection) |
| WL-2   | +30 / 0 | 1 R, 1 csv | 0 | 1 | S (run + record) |
| WL-3   | +60 / 0 | 0 | 1 | 0 | M (CG + Δ adaptation) |
| WL-4   | +25 / 0 | 1 R, 1 csv | 0 | 1 | S (regression bench) |
| WL-5   | +20 / 0 | 1 R, 1 csv | 0 | 1 | S (stretch bench) |
| WL-6   | +10 / 0 | 0 | 0 | 0 | S (counter audit) |
| **Total** | **~+285 / -27** | **5** | **2** | **3** | **~1.5 dev-day** |

## 8. Risks & Open Questions

1. **R LAPACK macro semantics.** Mitigated WL-0 — verify with WebSearch before edit per `~/.agents/skills/CLAUDE.md` §2 (mandatory for declarative configs).
2. **dsyevd LWORK query.** Standard idiom: call once with `LWORK=-1` for workspace size, allocate, call again. Documented in WL-1.
3. **`d_floor_retained = mean(Λ_keep)` semantic.** Confidence 75. Alternative: `min(Λ_keep)` (more aggressive) or `geometric_mean` (Tikhonov-flavored). Pick mean for symmetry with classical Levenberg-Marquardt and revisit only if T2 fails.
4. **Trust radius initial value.** `Δ=1.0` is heuristic. If WL-3 unit tests show pathological adaptation (oscillation), reconsider — but standard Nocedal-Wright recommendation supports unit init for normalized problems.
5. **Null-space gradient interpretation (question 7).** Convergence check `‖G‖_∞ < tol_abs` will include components in dropped directions where Newton no-ops. **Acceptable** because: (a) those components don't move under any algorithm starting from the same point, (b) the test of "have we found the calibration" is on residual gap (primary metric), not gradient norm. **Add NEWS entry** in WL-1 documenting this.
6. **dsyevr partial spectrum.** Could compute only top-k eigenvalues for n_λ > 100. Question 2 disposition: **stick with dsyevd** for n_λ ≤ 80 (typical leafblower scale); revisit only if profiling shows eigendecomp dominates wall time.
7. **Best-iterate may mask bugs.** WL-6 counter test detects this.

## 9. Reviewer Concerns from Prior Cycles (explicit address)

Three concerns raised during Epic-D plan-review-gate iter 1:

**(1) Scope reviewer: "truncated-SVD ≠ trust-region — pick one."**
Addressed: this plan implements **both** as a composed pipeline (eigendecomp+truncate is direction selection; Steihaug-CG is step length under trust radius). They are orthogonal mechanisms; truncated-SVD addresses unidentified directions, trust-region addresses step magnitude when the projected pseudoinverse step is too large. Confidence 90.

**(2) Feasibility reviewer: "cheap alternatives (μ_min, pivot retune) not yet refuted."**
Addressed: explicitly refuted by `benchmarks/results/lm_pivot_sweep.csv` (Epic-Dα): plateau at 2.79e-04 across `pivot ∈ {1e-12, 1e-10, 1e-8}`. Combined with Epic-B's `lm_mu_min ∈ {1e-12, 1e-10, 1e-8}` plateau, the cheap-first hypothesis is **empirically falsified**. Confidence 95.

**(3) Completeness reviewer: "what if the new method also stagnates?"**
Addressed: best-iterate fallback **kept** (WL-6 verifies it doesn't silently mask a bug). If T2 hard gate fails, plan halts with `SPEC_FAILURE` per Compliance Ultimatum §9.4 — no anti-pivot to "just bump ratio." Re-planning required. Confidence 90.

## 10. Out of Scope

- Continuation / homotopy methods (Epic-C territory; if needed, future Epic-E).
- Multi-start for kk1204 K=20 (Epic-F candidate).
- AUTO routing changes (none — explicit by §4 Forbidden).
- Removal of best-iterate code (deferred per WL-6 disposition; cleanup ticket only).
- Eigen3 / RcppEigen migration (not warranted at current n_λ scale).
- dsyevr partial-spectrum optimization (only if profiling later motivates).
- Removal / refactoring of the existing LM ρ-ratio acceptance test (composed with, not replaced by, this plan).
- Stale `n_homotopy_levels_used` field cleanup (deferred to separate Epic-C cleanup ticket).
