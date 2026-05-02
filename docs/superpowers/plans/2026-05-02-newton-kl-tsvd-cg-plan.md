# Epic-Dβ: Truncated-SVD + Steihaug-CG Trust-Region Newton (LAPACK-backed)

## Rev 2 Changelog (post plan-review-gate iter 1)

Plan-review-gate iter 1 verdict: REVISION NEEDED — Feasibility PASS w/ must-fix, Completeness FAIL ×7, Scope PASS w/ minor. Rev 2 fixes:

- **A — WL-2 gate semantics.** WL-2 changed from epic-halt SPEC_FAILURE on pinv-only failure to a *measurement* gate. CG (WL-3) exists precisely to handle pinv-insufficient cases. Halt only on T1/T4/T5 regression or build failure. (§6 WL-2)
- **B — r_bridge.cpp SEXP-pack update.** WL-1 now expressly updates `pack_solver_result` (or sibling block) in `src/r_bridge.cpp` to surface the new `n_projected_dims` integer scalar. `n_projected_dims` lives on `NewtonCalibResult` (not on `.base`), so it is plumbed via the newton-specific scalar block alongside `lm_mu_final`. `r_bridge.cpp` added to WL-1 §VI Files. (§6 WL-1)
- **C — T7 K=4 fixture concrete DGP.** Plan §6 WL-4 now specifies exact DGP: `n=2000`, `K=4`, `nj_per_margin = (3,4,3,2)`, uniform sample, mild target skew (0.4/0.3/0.2/0.1 for nj≥3, 0.5/0.5 for nj=2). Asserts `λ_max(H_pre)/λ_min(H_pre) < 1e6` AND `n_projected_dims == 0` AND `max_err < 1e-4`. (§6 WL-4)
- **D — NEWS.md timing.** NEWS.md edits colocated with code per discipline §4: WL-1 (truncated-SVD pseudoinverse), WL-3 (Steihaug-CG + Marquardt Δ), WL-6 (verdict bullet GATE_MET / PARTIAL / BLOCKED). (§6 WL-1, WL-3, WL-6)
- **E — Edge cases enumerated.** Plan §3 algorithm explicitly handles `λ_max == 0`, all-zero `H_pre`, `n_keep == 0` (LDLT fallback path = master baseline), `n_keep == n_λ` (no projection = master baseline). `mean(Λ_keep)` NaN-guard via the `n_keep==0` fallback. (§3)
- **F — WL-6 best-iterate counter site specified.** Increment site = the existing best-iterate-restoration block in `run_newton_inner` (the `best_iter_id = iter` block at `src/newton_calib.cpp:224`). Aggregate counter across T1/T2/T4/T5/T7/T8 invocations; print final per-test count. Decision: zero across all → file cleanup ticket. (§6 WL-6)
- **G — `n_homotopy_levels_used` deferral noted.** Plan §10 now explicitly cites discipline §2 (clean up only own newly-created orphans) — Epic-Dβ does NOT delete this Epic-C-deferred field; cleanup ticket post-WL-6.
- **H — Trust radius cap removed.** Per discipline §3 (derive numeric parameters from first principles), the unjustified `min(2·Δ, 100)` upper cap is removed. ρ>0.75 → Δ ← 2·Δ. Cited Nocedal-Wright Algorithm 4.1. (§3, §6 WL-3)
- **I — ρ predicted-decrease formula made explicit.** `m_predicted_decrease = -G^T·δ - 0.5·δ^T·H·δ` where H is the LM-damped Hessian; in retained eigenbasis: `-g_keep^T·δ_keep - 0.5·δ_keep^T·diag(Λ_damped)·δ_keep`. (§3)
- **J — `R_ext/Lapack.h` over bespoke extern "C".** WL-0 uses `#include <R_ext/Lapack.h>` and `F77_CALL(dsyevd)` (declared in `/usr/include/R/R_ext/Lapack.h`) rather than custom `dsyevd_` extern. Avoids signature drift. `src/lapack_wrappers.hpp` becomes a thin convenience header (or removed entirely; WL-0 picks one at edit time). (§1, §6 WL-0)
- **K — Edit `src/Makevars.in` (template), not `src/Makevars` (configure-generated).** Configure substitutes `@MVEC_LIBS@`. `Makevars` regenerates from `Makevars.in` on `./configure`. WL-0 edits the template. (§1, §6 WL-0)
- **L — `PKG_LIBS` ordering.** Standard R-pkg convention: package libs first, LAPACK suffix: `PKG_LIBS = @MVEC_LIBS@ $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` (i.e. retain `@MVEC_LIBS@` substitution and append LAPACK trio). (§6 WL-0)

## 1. Header

- **Epic ID:** Epic-Dβ
- **Plan ID:** `leafblower-wkmq` / Epic-Dβ rev2
- **Date:** 2026-05-02
- **Author:** controller (Opus, planning subagent)
- **Refs:**
  - Epic-B BLOCKED (KL-Newton μ-floor sweep, plateau 2.79e-04)
  - Epic-C BLOCKED (homotopy continuation, plateau 2.79e-04)
  - Epic-Dα ESCALATE (LM pivot sweep, plateau 2.79e-04)
  - Master HEAD: `benchmarks/results/lm_pivot_sweep.csv` (Epic-Dα CSV committed)
  - Source: `src/newton_calib.cpp` (`run_newton_inner` lambda, lines ~190-225 LM damping; LDLT line ~248); `src/r_bridge.cpp` `pack_solver_result` lambda at L394 (newton-specific scalar pack block follows)
  - Build: `src/Makevars` currently `-lm -lmvec` only — **LAPACK NOT linked**. Note: `src/Makevars` is configure-generated from `src/Makevars.in` (substitutes `@MVEC_LIBS@`); edit `Makevars.in` (template).
  - LAPACK header: `<R_ext/Lapack.h>` (system path `/usr/include/R/R_ext/Lapack.h`) declares `F77_NAME(dsyevd)` — use this rather than bespoke `extern "C"` per fix J.
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
2. F77_CALL(dsyevd)(H_copy, Λ[1..n_λ], V ← H_copy)   # V is column-eigvecs, Λ ascending
3. λ_max ← Λ[n_λ]
   # ── EDGE CASE A: degenerate input ──
   IF λ_max == 0 (or all-zero H_pre): δ ← 0; status ← (converged-or-NOCONV depending on ‖G‖); RETURN
4. thresh ← ratio_tsvd · λ_max
   R ← {i : Λ[i] ≥ thresh};  n_keep ← |R|
   # ── EDGE CASE B: n_keep == 0 ──
   IF n_keep == 0: SKIP projection, FALLBACK to existing LM-LDLT path (master baseline). Implicitly NaN-guards mean(Λ_keep).
   # ── EDGE CASE C: n_keep == n_λ ──
   # No projection; algorithm reduces to LM-LDLT in eigenbasis (semantically equivalent to master baseline LDLT solve).
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

### ρ-ratio predicted decrease (fix I)

```
m_predicted_decrease = -G^T·δ - 0.5·δ^T·H·δ
                     # H is the LM-DAMPED Hessian in the original basis.
                     # In retained eigenbasis (since H_proj diagonal):
                     # m_predicted_decrease = -g_keep^T·δ_keep - 0.5·δ_keep^T·diag(Λ_damped)·δ_keep
ρ = (f(λ) - f(λ + δ)) / m_predicted_decrease
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
ρ ← (f(λ) - f(λ + δ)) / m_predicted_decrease   # see ρ-ratio predicted decrease above
IF ρ > 0.75:    Δ ← 2 · Δ                       # NO upper cap (fix H, Nocedal-Wright Alg 4.1)
IF ρ < 0.25:    Δ ← Δ / 4
# else unchanged
Δ_init ← 1.0
```

This piggybacks on the existing LM ρ-ratio acceptance test in `run_newton_inner` (no new code path).

**Rationale (fix H):** the prior `min(2·Δ, 100)` upper cap was unjustified per discipline §3 (numeric parameters from first principles). Removed in rev2. Δ grows under successful steps; pathological growth would itself be diagnostic of a different issue (handled by ρ collapse on the next bad step). Cited Nocedal-Wright "Numerical Optimization" (2nd ed.) Algorithm 4.1, which has no upper cap.

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

### **WL-0: LAPACK linkage + R_ext/Lapack.h adoption**

**§I Goal.** Link LAPACK into the package and make `F77_CALL(dsyevd)` callable from `newton_calib.cpp` via the canonical `<R_ext/Lapack.h>` header (fix J). Confidence 95.

**§II Mechanism.**
- Edit `src/Makevars.in` (template; configure-generated `src/Makevars` regenerates on `./configure` — fix K). Append LAPACK trio after the existing `@MVEC_LIBS@` substitution:
  - Before: `PKG_LIBS = @MVEC_LIBS@`
  - After:  `PKG_LIBS = @MVEC_LIBS@ $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)`
  - Ordering rationale (fix L): package libs first, LAPACK suffix; standard R-exts §1.2.1 idiom. `$(FLIBS)` last so Fortran runtime symbols resolve after BLAS/LAPACK references them.
- Use `<R_ext/Lapack.h>` directly in `src/newton_calib.cpp`:
  ```cpp
  #include <R_ext/Lapack.h>
  // call: F77_CALL(dsyevd)(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, IWORK, &LIWORK, &INFO);
  ```
  R's `R_ext/Lapack.h` (system path `/usr/include/R/R_ext/Lapack.h`) declares `F77_NAME(dsyevd)` with the canonical signature. This avoids ABI/signature drift vs hand-rolled `extern "C"`.
- `src/lapack_wrappers.hpp` is **NOT** created in rev2 (fix J). If a thin convenience header is later needed (e.g. RAII workspace wrapper), file it as a separate WU; not required for WL-1.
- **Audit:** dummy compile-only smoke test in `src/newton_calib.cpp` (call `F77_CALL(dsyevd)` on a 2×2 identity, assert eigenvalues = {1,1}, gated `#ifdef LEAFBLOWER_LAPACK_SMOKE`).

**§III Acceptance.** `R CMD INSTALL --preclean .` succeeds (mandatory R compile gate per `~/.agents/skills/CLAUDE.md` §3). After install, `nm src/leafblower.so | grep -E 'dsyevd|F77'` shows symbol resolved (defined in libRlapack via dynamic link, or runtime-resolved). Re-running `./configure && make -C src` (where applicable) regenerates `Makevars` consistently with the edited `Makevars.in`.

**§IV Risks.** R's LAPACK macro names. **Mitigation:** verified against R-exts §1.2.1 + the live `/usr/include/R/R_ext/Lapack.h` declaring `F77_NAME(dsyevd)`. WebSearch redundant for rev2 since header location is confirmed.

**§V Forbidden.** Hardcoding `-llapack -lblas`. Linking against system LAPACK directly (must use R's resolution to match R's BLAS — OpenBLAS / MKL / Accelerate variants). Adding RcppEigen. Editing `src/Makevars` directly (configure-generated; edit `Makevars.in` only — fix K). Bespoke `extern "C"` re-declaration of `dsyevd_` (fix J).

**§VI Files.** `src/Makevars.in` (edit). `src/newton_calib.cpp` (add `#include <R_ext/Lapack.h>`; smoke test gated by macro). NOT touched: `src/Makevars` (regenerates from template). NOT created: `src/lapack_wrappers.hpp` (fix J).

---

### **WL-1: Eigendecomposition + truncated-SVD projection (no CG, direct pseudoinverse)**

**§I Goal.** Replace the LDLT solve in `run_newton_inner` with eigendecomp + truncated-SVD-pseudoinverse + LM in eigenbasis. No trust region yet (assume Δ = ∞). Confidence 85.

**§II Mechanism.** Implement steps 1-6 and 8 of the algorithm above; use `δ_keep_pinv` directly for `δ_keep`. Add `n_projected_dims` field to `NewtonCalibResult` struct. Default `ratio_tsvd = 1e-8` exposed as `newton_tsvd_ratio` parameter (R-side default). **Edge cases (fix E):** algorithm explicitly handles `λ_max == 0` (δ ← 0), all-zero `H_pre` (δ ← 0), `n_keep == 0` (LDLT fallback path = master baseline; NaN-guards `mean(Λ_keep)`), `n_keep == n_λ` (no projection; equivalent to existing LM-LDLT solve, semantic master baseline). **r_bridge.cpp SEXP-pack update (fix B):** `n_projected_dims` lives on `NewtonCalibResult` directly (NOT on `.base`), so it is plumbed via the newton-specific scalar block in `src/r_bridge.cpp` (alongside how `lm_mu_final` is currently surfaced — see grep evidence in plan §1 Refs). Sub-step: locate the newton-dispatch site, pack `n_projected_dims` as a new integer scalar field in the result list, extend the `res_names` / `res_list` SET_VECTOR_ELT block (mirror the `best_iter` slot pattern at `src/r_bridge.cpp:699`/`:707`). **NEWS.md (fix D):** add bullet "Newton-KL now applies truncated-SVD pseudoinverse to handle rank-deficient Hessian directions on overlapping-margin fixtures." Colocate edit with WL-1 commit.

**§III Acceptance.**
- T0 unit test: synthetic rank-deficient `H_pre` → `n_projected_dims > 0`, δ has zero component in null direction (1e-12).
- T1 baseline: K=4 well-conditioned bench unchanged within 1e-6 (no truncation should fire when spectrum is uniform).
- R-side: `n_projected_dims` round-trips from C++ to R via the result list (sanity check: read back from a small calibration call).
- NEWS.md updated in the same commit.
- Compile + tests pass.

**§IV Risks.** Eigvec sign convention from `dsyevd`. **Mitigation:** signs cancel in `V·δ_keep` since they appear in both projection and back-projection. **r_bridge.cpp drift risk (fix B):** failure to plumb `n_projected_dims` would silently lose the diagnostic. Mitigated by explicit sub-step.

**§V Forbidden.** Steihaug-CG (deferred to WL-3). Ratio retuning. Modifying LM ρ-ratio acceptance. Adding `n_projected_dims` to `.base` (it is a newton-specific diagnostic; goes on `NewtonCalibResult` directly).

**§VI Files.** `src/newton_calib.cpp`, `src/newton_calib.hpp` (struct field), `src/r_bridge.cpp` (newton scalar SEXP-pack — fix B), `R/newton.R` (parameter plumbing), `tests/testthat/test-newton-tsvd-projection.R` (new), `NEWS.md` (fix D).

---

### **WL-2: T2 stepstone — pinv-only measurement gate (fix A)**

**§I Goal.** Run K=4 stepstone benchmark with WL-1 implementation (pinv-only, no CG yet). **Measurement gate**, not epic-halt: report whether truncated-SVD pseudoinverse alone breaks the 2.79e-04 plateau. Confidence 70.

**§II Mechanism.** Reuse `benchmarks/lm_pivot_sweep.R` harness, swap algorithm flag. Single-row CSV `benchmarks/results/tsvd_pinv_T2.csv` with `gap`, `n_projected_dims`, `n_iter`, `wall_time_ms`.

**§III Acceptance (rev2 — fix A: NOT an epic-halt SPEC_FAILURE).**
- **PASS path:** if `gap < 1e-4` → great. Note in WL-2 result: "pinv alone clears T2; CG (WL-3) may be unnecessary — defer KEEP/REMOVE decision to WL-6."
- **MEASUREMENT path:** if `gap ≥ 1e-4` → record outcome as "pinv insufficient → CG required by WL-3"; PROCEED to WL-3. WL-3 is precisely the mechanism designed to handle pinv-insufficient cases (trust-region bounds step magnitude when projected pseudoinverse step is too large or numerically dominated). This is NOT a SPEC_FAILURE.
- **HALT criterion (rev2):** halt only on (a) T1 regression (well-conditioned baseline broken — caught here if applicable), (b) build failure, (c) downstream T4/T5 regression. WL-2 itself does NOT halt the epic on the pinv-only number.

**§IV Risks.** Truncation alone may be insufficient — CG was always planned to handle this (WL-3 exists for this reason). **Mitigation:** see Acceptance — proceed to WL-3 as planned.

**§V Forbidden.** Tuning ratio_tsvd silently to pass the (now soft) gate. Disabling LM. Increasing `max_iter` beyond master baseline. Skipping WL-3 even if WL-2 passes (CG mechanism remains in plan; WL-6 decides on retention).

**§VI Files.** `benchmarks/tsvd_pinv_T2.R` (new), `benchmarks/results/tsvd_pinv_T2.csv` (new).

---

### **WL-3: Steihaug-CG trust-region in retained subspace**

**§I Goal.** Add Steihaug-CG path for the case `‖δ_keep_pinv‖₂ > Δ`. Adaptive Δ via existing LM ρ-ratio. Confidence 80.

**§II Mechanism.** Implement step 7 alternative branch + Marquardt-style Δ adaptation. CG body ~30-40 LOC, operating on diagonal `Λ_damped`. Initial `Δ = 1.0`. Trust radius adapted post-step: ρ>0.75 → Δ ← 2·Δ (NO upper cap, fix H); ρ<0.25 → Δ ← Δ/4. Predicted decrease formula explicit (fix I): `m_predicted_decrease = -g_keep^T·δ_keep - 0.5·δ_keep^T·diag(Λ_damped)·δ_keep` (diagonal in retained eigenbasis). **NEWS.md (fix D):** add bullet "Newton-KL adds Steihaug-CG trust-region step in retained subspace; trust radius adapted via Marquardt gain ratio." Colocate with WL-3 commit.

**§III Acceptance.**
- Unit test: synthetic `H` with `‖δ_pinv‖₂ > Δ` → CG returns boundary step with `‖δ‖₂ = Δ` (1e-10).
- T1 stepstone result no worse than WL-2 (gap ≤ WL-2 gap + 1e-6).
- Δ trace shows growth on well-conditioned, contraction on ill-conditioned.

**§IV Risks.** CG numerical drift on ill-scaled diagonal. **Mitigation:** since Λ_damped is diagonal and PSD, CG converges in ≤ n_keep iters analytically; tolerance test redundant but kept defensively.

**§V Forbidden.** Removing direct-pseudoinverse fast path (must keep for the common case `‖δ_pinv‖₂ ≤ Δ`). Curvature checks (provably unnecessary for diagonal PSD).

**§VI Files.** `src/newton_calib.cpp` (CG body), `tests/testthat/test-newton-steihaug-cg.R` (new), `NEWS.md` (fix D bullet).

---

### **WL-4: T7 hard gate — K=4 over-projection regression**

**§I Goal.** Verify on a well-conditioned K=4 problem (uniform spectrum, no rank deficiency) that **no eigenvalues are dropped** and the gap matches WL-2 baseline within 1e-4. Confidence 90.

**§II Mechanism (concrete DGP — fix C).** Construct K=4 well-conditioned fixture per Epic-A's pattern:
- `n = 2000` observations.
- `K = 4` margins.
- `nj_per_margin = (3, 4, 3, 2)` (i.e. margin 1 has 3 categories, margin 2 has 4, margin 3 has 3, margin 4 has 2).
- Sample observations uniformly across the joint cell space.
- Mild target skew: for each margin with `nj ≥ 3`, target proportions `(0.4, 0.3, 0.2, 0.1)` truncated/normalized to fit `nj`; for `nj == 2`, target `(0.5, 0.5)`.
- Concretely: margin 1 (nj=3) targets `(0.4, 0.3, 0.3)` (renormalized); margin 2 (nj=4) targets `(0.4, 0.3, 0.2, 0.1)`; margin 3 (nj=3) targets `(0.4, 0.3, 0.3)`; margin 4 (nj=2) targets `(0.5, 0.5)`.

Run WL-3 implementation. Compute and assert:
1. `λ_max(H_pre) / λ_min(H_pre) < 1e6` (well-conditioned guarantee — derived from `H_pre` at iter 0).
2. `n_projected_dims == 0` (no truncation should fire at default ratio_tsvd = 1e-8 with a 6-orders-of-magnitude condition number).
3. `max_err < 1e-4`.

**§III Acceptance.** **HARD (all three):** `λ_max/λ_min < 1e6` AND `n_projected_dims == 0` AND `max_err < 1e-4`. Failure → over-projection bug or DGP miscalibration; halt and diagnose. Output single-row CSV with the three measured quantities.

**§IV Risks.** Default ratio 1e-8 is too aggressive for borderline-conditioned problems. **Mitigation:** the engineered DGP is structurally well-conditioned (uniform sample, mild skew); if it still triggers truncation, ratio default needs revision (escalate, do not silently retune).

**§V Forbidden.** Changing ratio_tsvd to make the gate pass. Loosening assertion to `n_projected_dims ≤ 1`. Modifying the DGP to game the condition number.

**§VI Files.** `benchmarks/tsvd_T7_well_conditioned.R` (new — implements DGP above), `benchmarks/results/tsvd_T7_K4.csv` (new — columns: `lambda_ratio`, `n_projected_dims`, `max_err`, `n_iter`, `wall_time_ms`).

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

**§I Goal.** Verify whether best-iterate fallback ever fires under WL-3 implementation. If never fires across {T1, T2, T4, T5, T7, T8 / full master test suite}, file an Epic-Dβ-cleanup bead ticket (do NOT delete in this epic per `~/.agents/skills/CLAUDE.md` §2 — clean up only your own newly created orphans — fix G). Confidence 60.

**§II Mechanism (concrete site — fix F).**
- **Increment site:** in `run_newton_inner` at the existing best-iterate-restoration block. Grep evidence: `src/newton_calib.cpp:224` contains `best_iter_id = iter;` inside a `current_gap < best_gap` comparison block (per Epic-A §5l pattern). Add:
  ```cpp
  static int best_iterate_fired_count = 0;  // gated by LEAFBLOWER_NEWTON_TRACE
  // ... inside the existing `if (current_gap < best_gap)` block:
  ++best_iterate_fired_count;
  ```
  Gate the counter declaration and increment under `#ifdef LEAFBLOWER_NEWTON_TRACE` so release builds are unaffected.
- **Aggregation:** count across T1, T2, T4, T5, T7, T8 (each test invocation re-enters `run_newton_inner` with the static counter accumulating across test runs in the same process). After each test, log the running count via the trace macro and reset (or read & reset) at the test boundary so per-test counts are recoverable. Print final per-test count and grand total.
- **Decision (fix F):** if grand total fires == 0 across all six tests → file cleanup ticket `leafblower-wkmq-cleanup-best-iterate` to remove best-iterate code in a follow-up epic. If any test shows fires > 0 → KEEP and document scenarios in NEWS.

**§III Acceptance.** Decision made and recorded:
- If counter > 0 anywhere → **KEEP** best-iterate, document scenarios in NEWS.
- If counter == 0 across all scenarios → **FILE** `leafblower-wkmq-cleanup-best-iterate` ticket for separate epic; do NOT remove in this PR.
- **NEWS.md (fix D):** add verdict-narrative bullet — one of `GATE_MET` (T2 + T7 hard gates passed), `PARTIAL` (T7 passed, T2 measurement non-trivial improvement, T5 partial), or `BLOCKED` (T1/T4/T5 regression or build failure). Colocate with WL-6 commit.

**§IV Risks.** Counter macro discipline. **Mitigation:** trace gated under `LEAFBLOWER_NEWTON_TRACE`; release builds unaffected. Static-counter aggregation across tests in the same process: acceptable since trace builds are diagnostic; if a fresh process boundary per test is required, fall back to per-test reset.

**§V Forbidden.** Deleting best-iterate code in this epic. Conflating WL-6 with refactoring of unrelated newton-calib code. Deleting the stale `n_homotopy_levels_used` field (Epic-C deferred — fix G; cleanup ticket only).

**§VI Files.** `src/newton_calib.cpp` (counter only, gated), `NEWS.md` (verdict bullet — fix D), bead ticket file (created via `bd create`, not in this PR's diff).

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
- Stale `n_homotopy_levels_used` field cleanup (Epic-C deferred). **Fix G justification:** per `~/.agents/skills/CLAUDE.md` §2 "clean up only your own newly created orphans; mention and create ticket, but do not delete pre-existing dead code unless explicitly asked." This field predates Epic-Dβ; deletion is filed as a separate cleanup ticket post-WL-6, NOT removed in this epic.
