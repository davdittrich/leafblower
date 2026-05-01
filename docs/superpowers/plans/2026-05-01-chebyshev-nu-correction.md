# Chebyshev IPM ν Dual Correction — Implementation Plan

**Date**: 2026-05-01
**Ticket**: leafblower-gvza — research(chebyshev): ν dual correction — enabling dnu without causing INFEAS
**File scope**: `src/chebyshev.cpp` only (uses existing `compute_normal_equations_reduced` from `src/calib_linalg.{hpp,cpp}`)
**Status**: Draft — for review before implementation

---

## Mechanism

Add the normalization-dual ν via reference-category elimination: drop the last category per multi-category margin from the Newton system to break the algebraic degeneracy `schur_nu = D_ν − e^T N_eff^{-1} e ≈ 0`. Once `schur_nu` is non-degenerate, perform a third LDLT back-solve to compute `Δν` via Schur-complement block elimination, applied each Mehrotra phase (predictor + corrector). Confidence: 90 (matches the prior 2026-04-26 spec, infrastructure already in place at `src/chebyshev.cpp:84–97,322–351` and `src/calib_linalg.hpp:28–36`).

## Forbidden

- No renorm workaround / post-hoc rescaling of `X` after each Newton step (superseded approach in the now-defunct `cheb-renorm` branch).
- No use of full-nct `e_vec` SM correction without reference elimination (yields `schur_nu ≈ 0` for sum-to-1 targets → divide-by-near-zero → INFEAS, the failure mode named in the ticket title).
- No new public API in `chebyshev.hpp`. No change to `LpVariant` enum. No change to `ChebyshevResult` struct layout.
- No silent guard `if (schur_nu < eps) dnu = 0` applied to the **full-nct** system (that disables the fix on the very systems we need it for). The `kSchurNuMin = 1e-8` guard in T2/T3 is applied to the **reduced** system after reference elimination — this is intentional and safe because reference elimination makes `schur_nu_red` generically non-zero for well-conditioned problems (the problematic near-zero only arises in the full-nct system due to sum-to-1 degeneracy, which the elimination resolves). If `schur_nu_red < kSchurNuMin` fires in the reduced system, it signals genuine ill-conditioning, not the sum-to-1 degeneracy; `dnu=0` is the correct fallback in that case.
- No bypass of `R CMD INSTALL --preclean .` between tasks.

## Audit

- **Spy on Newton step sum-preservation**: with `verbose >= 2` log `r_nu = W − n_d` and `Δν` per iteration; assert `|W_after − n_d| < 1e-6 · n_d` after the first 5 iterations on K=4 synthetic.
- **Spy on `schur_nu`**: existing diagnostic at `src/chebyshev.cpp:324–351` already logs `schur_nu` at iter 0 with `verbose >= 2`; extend to log `schur_nu` and `dnu` per iteration during T1 diagnostic.
- **Wiring spy**: add `expect_gt(schur_nu, 1e-6)` in a new `tests/testthat/test-calib-linalg.R` fixture that captures the verbose log on a K=2 case.
- **Mehrotra symmetry**: confirm Phase A and Phase B both apply the ν correction with the same `e_red`, `w_e_red`, `D_nu`; otherwise the predictor and corrector live on different affine subspaces and `mu_aff` is meaningless.

---

## Background — what is currently in `src/chebyshev.cpp`

Read confirms (file: `src/chebyshev.cpp`):

- **Calibration LP**: lines 14–18 (entry), 62–78 (`Tgt`, `T_flat`, `w_kj = n_d`).
- **Reference-elimination scaffolding already present** (built but only used for the `schur_nu` diagnostic):
  - `nct_red_count`, `nct_red`, `full_to_red[]` constructed at lines 83–97.
  - Diagnostic-only block at lines 322–351 calls `compute_normal_equations_reduced` and `ldlt_factor_inplace` on `N_red`, computes `schur_nu = D_nu − e_red^T · w_e_red` once at iter 0 with `verbose >= 2`, then discards everything.
- **Normal equations (full nct)**: built and Jacobi-scaled at lines 301–313 into `N0` (size `nct × nct`), factored in-place at line 316.
- **Mehrotra Phase A (predictor, σ=0)**: lines 492–567. RHS at 496–498. δ-direction SM correction at 514–524. `dlambda_A` is full-size `nct`; no ν step.
- **Mehrotra Phase B (corrector, second-order)**: lines 589–700. RHS at 612–619. δ-direction SM correction at 640–644. `dlambda_B` is full-size `nct`; no ν step.
- **Sherman-Morrison for δ (existing)**: `u_vec[m] = w_kj[m]` at line 298; `Theta`, `alpha_sm` at 294–297; first SM in both phases at 516–524 (Phase A) and 640–644 (Phase B). The plan's ν correction is a **second** Sherman-Morrison, layered on top of the δ one.
- **`dX` reconstruction from `dlambda`**: Phase A at 527–534, Phase B at 647–654. Both apply `dX[c] = D_eff[c] · Σ_{k,g} dlambda[cat_offset[k]+g]`. The ν fix appends `dX[c] += D_eff[c] · Δν`.
- **Slack/dual updates** (unchanged): Phase A at 405–419 (degenerate branch; never executed since `n_comp > 0`), Phase B at 683–700, line search at 703–717.

Read confirms (file: `src/calib_linalg.hpp:28–36`): `compute_normal_equations_reduced` already exists with the exact signature this plan needs. **No changes to `calib_linalg.{hpp,cpp}` are required.**

Read confirms (file: `tests/testthat/test-chebyshev.R`):
- `T_cheby_warm` (K=3, n=5000) currently asserts `status==0`, `max_err <= raking * 1.001 + 1e-10`, `max_err < 1e-3`.
- `T_cheby_K9` (Stepstone K=9) skips when parquet not present; asserts `chebyshev max_err <= greenkhorn`.

Confidence on the read: 95 (file content quoted directly).

---

## Math — derivation of the ν correction (T2)

The Chebyshev calibration LP is:

```
min_{X, δ}  δ
s.t.   −w_kj[m]·δ ≤ Σ_c A[m,c]·X[c] − T_flat[m]·W ≤ w_kj[m]·δ   ∀ m            (margins, 2·nct)
       L_cell[c] ≤ X[c] ≤ U_cell[c]                              ∀ c            (box, 2·M_cell)
       Σ_c X[c] = n_d                                                            (normalization)
       δ ≥ 0
```

with `W := Σ_c X[c]`. The unique source of the ν fix is the equality `Σ_c X[c] = n_d`. Its Lagrange multiplier ν enters the KKT stationarity for `X[c]`:

```
0 = Σ_m A[m,c]·λ[m]   −   y_lo[c] + y_hi[c]   −   ν                            (∀ c)
```

(λ collects the two-sided margin duals via the standard `λ = y_dn − y_up` reduction.) Combining with the reduced Newton system on the margin block gives the **augmented normal equations**:

```
[ N_eff      e ] [ Δλ ]   [ rhs_margin ]
[  e^T     D_ν ] [ Δν ] = [   −r_nu    ]
```

where:
- `N_eff = N_0 + α_δ · u·u^T`, the δ-Sherman-Morrison-corrected normal matrix (already built).
- `e[m] = Σ_{c ∈ margin m} D_eff[c]` — column sums of `D_eff` aggregated into margin space.
- `D_ν = Σ_c D_eff[c]` — total effective barrier weight.
- `r_nu = W − n_d` — primal normalization residual (drift).

**Block (Schur) elimination** of Δν gives:

```
Δν · ( D_ν − e^T · N_eff^{-1} · e )   =   − r_nu − e^T · (N_eff^{-1} · rhs_margin)
                                        =   − r_nu − e^T · Δλ_uncorrected
```

Define `schur_nu := D_ν − e^T · N_eff^{-1} · e` (Schur complement). Then:

```
Δν       =  ( − r_nu − e^T · Δλ_uncorrected ) / schur_nu
Δλ_final =  Δλ_uncorrected − Δν · w_e        where  w_e := N_eff^{-1} · e
ΔX[c]   +=  D_ν · Δν / D_ν · D_eff[c]  →  ΔX[c] += D_eff[c] · Δν      (∀ c)
```

The third back-solve (`w_e = N_eff^{-1} · e`) costs one extra `ldlt_solve` per Mehrotra phase, plus a δ-SM correction on `w_e`. Total: from 2 back-solves/phase → 3 back-solves/phase.

### The degeneracy (root cause of past INFEAS)

For any margin k with sum-to-1 targets (the standard case), `Σ_{j=0}^{J_k−1} A[cat_offset[k]+j, c] = 1` for all c, so the rows of A for margin k sum to the all-ones row. This puts `e` in the row span of `N_0`, giving:

```
e^T · N_0^{-1} · e  =  D_ν   (algebraically, modulo numerical noise)
schur_nu            ≈  0      (degenerate denominator → divide-by-zero → INFEAS)
```

The full-nct ν correction in the abandoned `2026-04-26-chebyshev-nu-fix.md` plan suffered exactly this: dividing by `schur_nu ≈ 0` injects garbage Δν, blows up the duals, and triggers `kInfeasPersistence` after ~5 iterations. **This is the failure mode the ticket title names: "enabling dnu without causing INFEAS".**

### The fix — reference-category elimination

For each margin k with `cat_counts[k] >= 2`, drop the last category's row from A. The reduced incidence matrix `A_red` (size `nct_red × M_cell`) no longer has the Σ-rows-equals-ones property because each cell now contributes to only `K' − 1` rows for any margin k' that had multiple categories. Therefore `e_red` is generically NOT in the row span of `N_red = A_red · D_eff · A_red^T`, and `schur_nu = D_ν − e_red^T · N_red^{-1} · e_red > 0` strictly (modulo numerical eps).

Reference-category slack/dual variables `s_up[m_ref], s_dn[m_ref], y_up[m_ref], y_dn[m_ref]` remain in the LP at full size; only the Newton system shrinks. The reference-margin calibration error is enforced **by difference**: once `Σ X[c] = n_d` holds, the reference category's marginal is determined by subtraction from `W = n_d` and the other categories' marginals.

Rule (already encoded at lines 84–86): `cat_counts[k] >= 2` is the gate. Single-category margins contribute no degenerate row and need no reference elimination.

Confidence on the math: 90 (matches Mehrotra-Wright Ch. 11 augmented-system presentation; the reference-elimination identification is the standard ANOVA dummy-variable trick).

---

## T1 — Diagnostic: locate the K threshold

**Ticket**: file as `chebyshev: T1 diagnostic — find K-threshold for current INFEAS/NOCONV` (one ticket per task per global rule).

### Goal

Empirically determine the smallest K at which the **current** (no-ν) Chebyshev IPM stops converging on synthetic problems with sum-to-1 targets. This bounds the regression risk window for T3 and gives a reproducible test corpus for T4.

### Method

Run a sweep over `K ∈ {1,2,3,4,5,6,9}` with fixed `n=2000`, modest `J_k ∈ {2,3,4}`, sum-to-1 targets. Capture: `status`, `iterations`, `max_error`, and (with `verbose=2`) the printed `schur_nu` from the iter-0 diagnostic at `src/chebyshev.cpp:324–351`.

Driver (place in scratch, **do not commit**):

```r
# scratch/diag-chebyshev-K-sweep.R
suppressPackageStartupMessages({library(leafblower)})
set.seed(2026); n <- 2000L
sweep_K <- function(K) {
  Js <- pmin(2L + (seq_len(K) %% 3L), 4L)
  df <- as.data.frame(lapply(Js, function(J) factor(sample(seq_len(J), n, TRUE))))
  names(df) <- paste0("v", seq_len(K))
  tgt <- setNames(lapply(Js, function(J) {
    p <- runif(J); p <- p/sum(p); setNames(p, as.character(seq_len(J)))
  }), names(df))
  t0 <- proc.time()["elapsed"]
  r <- suppressWarnings(harvest(df, tgt, method="chebyshev",
                                 max_iterations=300L, attach_weights=FALSE,
                                 verbose=2L))
  res <- attr(r, "result")
  data.frame(K=K, status=res$status, iters=res$iterations,
             max_err=res$max_error, time=proc.time()["elapsed"]-t0)
}
do.call(rbind, lapply(c(1,2,3,4,5,6,9), sweep_K))
```

### Compile gate
```bash
R CMD INSTALL --preclean .
```
Expected: `* DONE (leafblower)`.

### Test gate (existing chebyshev tests, MUST stay GREEN — this task changes no source)
```bash
Rscript -e 'testthat::test_dir("tests/testthat", filter="chebyshev")'
```
Expected: 4 PASS (T_greg_warn, T_cheby_warm, T_cheby_warm_fallback; T_cheby_K9 skipped without parquet — use `cd benchmarks` parquet if local).

### Diagnostic gate
```bash
Rscript scratch/diag-chebyshev-K-sweep.R 2>&1 | tee /tmp/cheb-K-sweep.txt
grep "schur_nu" /tmp/cheb-K-sweep.txt
```

### Acceptance

- Record the lowest K with `status != 0`. Document the recorded `schur_nu` magnitude at iter 0 for each K.
- Expected (from prior K=9 Stepstone failure): `K=1,2,3` converge cleanly; `schur_nu` measured at iter 0 is `O(1e-12 · D_ν)` or smaller across all K, confirming the degeneracy is K-independent in nature but compounds in step magnitude.
- Output goes into the T1 ticket as a comment; **no source edits in T1**. Confidence: 80 (depends on stochastic seed; the existing `2026-04-28-chebyshev-sinkhorn-greg-correctness.md` notes K=9 NOCONV with max_err≈0.139 as the canonical failure).

---

## T2 — Math derivation in code-ready form

**Ticket**: file as `chebyshev: T2 math — formalize ν Schur-complement block elim`.

### Goal

Document the exact code-level formulas to drop into Phase A and Phase B of the Mehrotra loop. **No code changes.** Output is the math block recorded as a ticket comment that T3 references verbatim.

### Formulas (verbatim, used by T3)

Indices: `m ∈ [0, nct)` (full margin space); `nr ∈ [0, nct_red)` (reduced space, reference categories dropped); `c ∈ [0, M_cell)`.

**Once per IPM iteration (after `D_eff[]` is computed at `src/chebyshev.cpp:271–281`, before any back-solve):**

```cpp
// Reduced-space e_red[nr] = Σ_{c ∈ margin red_to_full[nr]} D_eff[c]
// D_nu = Σ_c D_eff[c]   (independent of reference elim)
std::fill(e_red.begin(), e_red.end(), 0.0);
double D_nu = 0.0;
for (int c = 0; c < ct.M_cell; c++) {
    D_nu += D_eff[c];
    for (int k = 0; k < st.K; k++) {
        int g = ct.g_per_cell[k][c];
        if (g < 0 || g >= st.cat_counts[k]) continue;
        int nr = full_to_red[cat_offset[k] + g];
        if (nr >= 0) e_red[nr] += D_eff[c];   // skip references (nr == -1)
    }
}
const double r_nu = W - n_d;   // W already computed above at chebyshev.cpp:222–223
```

**Per phase** (Phase A: `σ=0` predictor; Phase B: corrector with second-order term), after the existing δ-SM correction has been applied to `dlambda_red`:

```cpp
// Third back-solve: w_e_red = N_red^{-1} · e_red  (jacobi-scaled in/out)
for (int nr = 0; nr < nct_red; nr++) tmp_red[nr] = D_jac_red[nr] * e_red[nr];
ldlt_solve(N_red.data(), static_cast<size_t>(nct_red), tmp_red.data());
for (int nr = 0; nr < nct_red; nr++) w_e_red[nr] = D_jac_red[nr] * tmp_red[nr];

// δ-SM correction on w_e_red  (same alpha_sm, sm_denom, v_red as in δ correction)
double ute = 0.0;
for (int nr = 0; nr < nct_red; nr++) ute += u_red[nr] * w_e_red[nr];
double sm_coeff_e = (std::fabs(sm_denom) > 1e-300) ? (alpha_sm * ute / sm_denom) : 0.0;
for (int nr = 0; nr < nct_red; nr++) w_e_red[nr] -= sm_coeff_e * v_red[nr];

// Schur denominator and Δν
double eTw_e = 0.0, eTdlambda = 0.0;
for (int nr = 0; nr < nct_red; nr++) {
    eTw_e     += e_red[nr] * w_e_red[nr];
    eTdlambda += e_red[nr] * dlambda_red[nr];
}
const double schur_nu = D_nu - eTw_e;
const double dnu = (schur_nu > kSchurNuMin) ? (-r_nu - eTdlambda) / schur_nu : 0.0;

// Apply correction to dlambda_red and to dX[c]
for (int nr = 0; nr < nct_red; nr++) dlambda_red[nr] -= dnu * w_e_red[nr];
for (int c = 0;  c < ct.M_cell; c++) dX[c] += D_eff[c] * dnu;

// d_delta uses w_dot_dlambda computed AFTER the ν correction
double w_dot_dlambda = 0.0;
for (int nr = 0; nr < nct_red; nr++)
    w_dot_dlambda += w_kj[red_to_full[nr]] * dlambda_red[nr];
```

Constant: `kSchurNuMin = 1e-8` (with `D_nu` of order `n_d · M_cell`, `1e-8` is ~14 orders below typical magnitudes; firmer than the spec's `1e-14`).

Confidence on the formulas: 90 (block elimination is textbook; reference elimination is the only non-obvious piece and matches the prior approved spec).

### T2 baseline gate (no code changes — verify existing tests still pass before T3 begins)
```bash
R CMD INSTALL --preclean .
Rscript -e 'testthat::test_dir("tests/testthat", filter="chebyshev")'
```
Expected: `* DONE (leafblower)` and 4/4 chebyshev tests GREEN. If any test is red before T3 starts, **halt** — do not proceed.

---

## T3 — Implementation in `src/chebyshev.cpp`

**Ticket**: file as `chebyshev: T3 impl — ν correction in both Mehrotra phases`.

### Files touched

`src/chebyshev.cpp` only. **No header change.** No `calib_linalg.{hpp,cpp}` change (`compute_normal_equations_reduced` already exists).

### Step 3.1 — Hoist new work vectors

**Before** (lines 180–197 of `src/chebyshev.cpp`, current):
```cpp
const int max_cats = lbw::max_cats_count(st.K, st.cat_counts);
std::vector<double> D_eff(ct.M_cell), D_marg(nct);
// Full nct system (for δ Sherman-Morrison — must use full space to preserve E1/E2)
std::vector<double> N0((size_t)nct * (size_t)nct);
std::vector<double> u_vec(nct), v_vec(nct);
std::vector<double> dlambda(nct);
std::vector<double> dX(ct.M_cell);
std::vector<double> dS_up(nct), dS_dn(nct);
std::vector<double> dY_lo(ct.M_cell), dY_hi(ct.M_cell), dY_up(nct), dY_dn(nct);
std::vector<double> delta_S(nct);
std::vector<double> bucket_tmp(max_cats);
// Mehrotra predictor-corrector workspace
std::vector<double> D_jac(nct);
std::vector<double> rhs_A(nct), rhs_B(nct);
std::vector<double> dlambda_A(nct), dlambda_B(nct);
std::vector<double> dX_A(ct.M_cell), dX_B(ct.M_cell);
std::vector<double> dS_up_A(nct), dS_dn_A(nct);
double best_delta = delta;
```

**After** — replace the full-nct Newton vectors with reduced-space versions; keep duals at full-nct:
```cpp
const int max_cats = lbw::max_cats_count(st.K, st.cat_counts);
std::vector<double> D_eff(ct.M_cell), D_marg(nct);

// Reduced normal equations + reduced Newton vectors (reference-cat elim → schur_nu > 0)
std::vector<double> N_red((size_t)nct_red * (size_t)nct_red);
std::vector<double> u_red(nct_red), v_red(nct_red);
std::vector<double> e_red(nct_red), w_e_red(nct_red), tmp_red(nct_red);
std::vector<double> rhs_A_red(nct_red), rhs_B_red(nct_red);
std::vector<double> dlambda_A_red(nct_red), dlambda_B_red(nct_red);
std::vector<double> dlambda_red(nct_red);            // SM-corrected, post-ν
std::vector<double> D_jac_red(nct_red);

// red_to_full inverse map (built once; mirrors full_to_red built at L88–97)
std::vector<int> red_to_full(nct_red);
{
    int nr = 0;
    for (int m = 0; m < nct; m++) if (full_to_red[m] >= 0) red_to_full[nr++] = m;
}

// Pre-fill u_red[nr] = w_kj[red_to_full[nr]] (constant across iterations: w_kj[m] = n_d)
for (int nr = 0; nr < nct_red; nr++) u_red[nr] = w_kj[red_to_full[nr]];

// Full-nct outputs (unchanged: duals, slacks, primal step on cells/δ)
std::vector<double> dX(ct.M_cell);
std::vector<double> dS_up(nct), dS_dn(nct);
std::vector<double> dY_lo(ct.M_cell), dY_hi(ct.M_cell), dY_up(nct), dY_dn(nct);
std::vector<double> delta_S(nct);
std::vector<double> bucket_tmp(max_cats);
std::vector<double> dX_A(ct.M_cell), dX_B(ct.M_cell);
std::vector<double> dS_up_A(nct), dS_dn_A(nct);
double best_delta = delta;

static constexpr double kSchurNuMin = 1e-8;          // tight; D_nu = O(n_d·M_cell)
```

**Drop**: `N0`, `u_vec`, `v_vec`, `dlambda`, `dlambda_A`, `dlambda_B`, `D_jac`, `rhs_A`, `rhs_B`. They are replaced by their `_red` counterparts.

### Step 3.2 — Replace `compute_normal_equations` with `_reduced` variant

**Before** (lines 300–313):
```cpp
if (lbw::compute_normal_equations(ct, D_eff.data(), N0.data(),
                                  cat_offset.data(), st.K,
                                  static_cast<size_t>(nct)) != RK_OK) {
    res.base.status = RK_ERR_BADARG; return res;
}
for (int j = 0; j < nct; j++)
    D_jac[j] = 1.0 / std::sqrt(std::max(N0[(size_t)j*nct+j], 1e-12));
for (int i = 0; i < nct; i++)
    for (int j = 0; j < nct; j++)
        N0[(size_t)i*nct+j] *= D_jac[i] * D_jac[j];
if (ldlt_factor_inplace(N0.data(), static_cast<size_t>(nct), kEpsLdlt) != RK_OK) {
    res.base.status = RK_ERR_BADARG; return res;
}
res.n_factorizations++;
```

**After** — build `N_red` directly, factor in reduced space:
```cpp
if (lbw::compute_normal_equations_reduced(ct, D_eff.data(), N_red.data(),
                                          cat_offset.data(), st.K,
                                          static_cast<size_t>(nct_red),
                                          full_to_red.data()) != RK_OK) {
    res.base.status = RK_ERR_BADARG; return res;
}
for (int j = 0; j < nct_red; j++)
    D_jac_red[j] = 1.0 / std::sqrt(std::max(N_red[(size_t)j*nct_red+j], 1e-12));
for (int i = 0; i < nct_red; i++)
    for (int j = 0; j < nct_red; j++)
        N_red[(size_t)i*nct_red+j] *= D_jac_red[i] * D_jac_red[j];
if (ldlt_factor_inplace(N_red.data(), static_cast<size_t>(nct_red), kEpsLdlt) != RK_OK) {
    res.base.status = RK_ERR_BADARG; return res;
}
res.n_factorizations++;
```

### Step 3.3 — Compute `e_red`, `D_nu`, `r_nu` once per iteration

**Insert** immediately after the LDLT factor of `N_red` (i.e., just before the diagnostic block that currently spans lines 322–351; the diagnostic block can be **removed** — its entire purpose was to log `schur_nu` for inspection, and T3 now wires `schur_nu` into the actual solve):

```cpp
// Once-per-iteration ν workspace
std::fill(e_red.begin(), e_red.end(), 0.0);
double D_nu = 0.0;
for (int c = 0; c < ct.M_cell; c++) {
    D_nu += D_eff[c];
    for (int k = 0; k < st.K; k++) {
        int g = ct.g_per_cell[k][c];
        if (g < 0 || g >= st.cat_counts[k]) continue;
        int nr = full_to_red[cat_offset[k] + g];
        if (nr >= 0) e_red[nr] += D_eff[c];
    }
}
const double r_nu = W - n_d;   // W from L222–223
```

**Optional preserved diagnostic**: keep the `verbose >= 2 && iter == 0` log line that prints `schur_nu`, but compute it from the values already in scope rather than re-factoring a second matrix:
```cpp
if (st.verbose >= 2 && iter == 0) {
    // schur_nu computed from already-factored N_red below in Phase A
    // (deferred — see Step 3.4 logging line)
}
```

### Step 3.4 — Phase A (predictor) reduced solve + ν correction

**Before** (lines 492–567 — Phase A predictor):
```cpp
for (int m = 0; m < nct; m++)
    rhs_A[m] = -(S[m] - T_flat[m]*W) + D_marg[m] * (-y_up[m] + y_dn[m]);
double margin_delta_center_A = 0.0;
for (int m = 0; m < nct; m++)
    margin_delta_center_A += w_kj[m]*(-y_up[m] - y_dn[m]);
double rmu_delta_A = -s_delta*y_delta;

for (int m = 0; m < nct; m++) rhs_A[m] *= D_jac[m];
ldlt_solve(N0.data(), static_cast<size_t>(nct), rhs_A.data());
for (int m = 0; m < nct; m++) dlambda_A[m] = D_jac[m] * rhs_A[m];

// δ-SM correction (Phase A)
for (int m = 0; m < nct; m++) v_vec[m] = D_jac[m] * u_vec[m];
ldlt_solve(N0.data(), static_cast<size_t>(nct), v_vec.data());
for (int m = 0; m < nct; m++) v_vec[m] *= D_jac[m];
double utv = 0.0, utw_A = 0.0;
for (int m = 0; m < nct; m++) { utv += u_vec[m]*v_vec[m]; utw_A += u_vec[m]*dlambda_A[m]; }
double sm_denom = 1.0 + alpha_sm*utv;
double sm_coeff_A = (std::fabs(sm_denom) > 1e-300) ? (alpha_sm*utw_A/sm_denom) : 0.0;
for (int m = 0; m < nct; m++) dlambda_A[m] -= sm_coeff_A * v_vec[m];

// dx_A from dlambda_A
std::fill(dX_A.begin(), dX_A.end(), 0.0);
for (int k = 0; k < st.K; k++)
    for (int c = 0; c < ct.M_cell; c++) {
        int g = ct.g_per_cell[k][c];
        if (g >= 0 && g < st.cat_counts[k])
            dX_A[c] += dlambda_A[cat_offset[k]+g];
    }
for (int c = 0; c < ct.M_cell; c++) dX_A[c] = D_eff[c] * dX_A[c];

double w_dot_dlambda_A = 0.0;
for (int m = 0; m < nct; m++) w_dot_dlambda_A += w_kj[m]*dlambda_A[m];
double d_delta_A = alpha_sm * (
    rmu_delta_A/s_delta + margin_delta_center_A
    - r_delta_stat - (y_delta/s_delta)*w_dot_dlambda_A );
```

**After** — same structure, reduced-space, with ν correction inserted between SM and `dX_A` reconstruction:

```cpp
// ── Phase A: affine predictor (σ = 0, no centering) ──────────────────
// RHS in reduced space
for (int nr = 0; nr < nct_red; nr++) {
    int m = red_to_full[nr];
    rhs_A_red[nr] = -(S[m] - T_flat[m]*W) + D_marg[m] * (-y_up[m] + y_dn[m]);
}
// margin_delta_center_A still sums over full nct (δ stationarity uses all margins)
double margin_delta_center_A = 0.0;
for (int m = 0; m < nct; m++)
    margin_delta_center_A += w_kj[m]*(-y_up[m] - y_dn[m]);
double rmu_delta_A = -s_delta*y_delta;

// Solve N_red · dlambda_A_red = rhs_A_red (jacobi-scaled)
for (int nr = 0; nr < nct_red; nr++) rhs_A_red[nr] *= D_jac_red[nr];
ldlt_solve(N_red.data(), static_cast<size_t>(nct_red), rhs_A_red.data());
for (int nr = 0; nr < nct_red; nr++) dlambda_A_red[nr] = D_jac_red[nr] * rhs_A_red[nr];

// δ-SM correction on dlambda_A_red (first SM)
for (int nr = 0; nr < nct_red; nr++) tmp_red[nr] = D_jac_red[nr] * u_red[nr];
ldlt_solve(N_red.data(), static_cast<size_t>(nct_red), tmp_red.data());
for (int nr = 0; nr < nct_red; nr++) v_red[nr] = D_jac_red[nr] * tmp_red[nr];
double utv = 0.0, utw_A = 0.0;
for (int nr = 0; nr < nct_red; nr++) {
    utv   += u_red[nr]*v_red[nr];
    utw_A += u_red[nr]*dlambda_A_red[nr];
}
double sm_denom = 1.0 + alpha_sm*utv;
double sm_coeff_A = (std::fabs(sm_denom) > 1e-300) ? (alpha_sm*utw_A/sm_denom) : 0.0;
for (int nr = 0; nr < nct_red; nr++) dlambda_A_red[nr] -= sm_coeff_A * v_red[nr];

// ── ν correction (Phase A): third back-solve + Schur Δν ─────────────
for (int nr = 0; nr < nct_red; nr++) tmp_red[nr] = D_jac_red[nr] * e_red[nr];
ldlt_solve(N_red.data(), static_cast<size_t>(nct_red), tmp_red.data());
for (int nr = 0; nr < nct_red; nr++) w_e_red[nr] = D_jac_red[nr] * tmp_red[nr];
double ute_A = 0.0;
for (int nr = 0; nr < nct_red; nr++) ute_A += u_red[nr] * w_e_red[nr];
double sm_coeff_e_A = (std::fabs(sm_denom) > 1e-300) ? (alpha_sm * ute_A / sm_denom) : 0.0;
for (int nr = 0; nr < nct_red; nr++) w_e_red[nr] -= sm_coeff_e_A * v_red[nr];

double eTw_e_A = 0.0, eTdl_A = 0.0;
for (int nr = 0; nr < nct_red; nr++) {
    eTw_e_A += e_red[nr] * w_e_red[nr];
    eTdl_A  += e_red[nr] * dlambda_A_red[nr];
}
const double schur_nu_A = D_nu - eTw_e_A;
const double dnu_A = (schur_nu_A > kSchurNuMin) ? (-r_nu - eTdl_A) / schur_nu_A : 0.0;
for (int nr = 0; nr < nct_red; nr++) dlambda_A_red[nr] -= dnu_A * w_e_red[nr];

// Verbose iter-0 diagnostic (replaces removed L324–351 block)
if (st.verbose >= 2 && iter == 0) {
    char msg[160];
    std::snprintf(msg, sizeof(msg),
        "chebyshev: schur_nu=%.4e dnu=%.4e r_nu=%.4e (iter 0, Phase A)",
        schur_nu_A, dnu_A, r_nu);
    st.log(msg);
}

// dX_A from dlambda_A_red + Δν · D_eff
std::fill(dX_A.begin(), dX_A.end(), 0.0);
for (int c = 0; c < ct.M_cell; c++) {
    double sum_dlam = 0.0;
    for (int k = 0; k < st.K; k++) {
        int g = ct.g_per_cell[k][c];
        if (g < 0 || g >= st.cat_counts[k]) continue;
        int nr = full_to_red[cat_offset[k] + g];
        if (nr >= 0) sum_dlam += dlambda_A_red[nr];
    }
    dX_A[c] = D_eff[c] * (sum_dlam + dnu_A);
}

// w_dot_dlambda_A AFTER ν correction (Phase A)
double w_dot_dlambda_A = 0.0;
for (int nr = 0; nr < nct_red; nr++)
    w_dot_dlambda_A += w_kj[red_to_full[nr]] * dlambda_A_red[nr];
double d_delta_A = alpha_sm * (
    rmu_delta_A/s_delta + margin_delta_center_A
    - r_delta_stat - (y_delta/s_delta)*w_dot_dlambda_A );
```

The `delta_S`, `dS_up_A`, `dS_dn_A`, `alpha_aff`, `mu_aff`, σ computation downstream of Phase A (lines 545–587 of the current file) are **unchanged** — they consume `dX_A` and `d_delta_A`, not `dlambda_A`.

### Step 3.5 — Phase B (corrector) reduced solve + ν correction

**Before** (lines 612–700 of current file, Phase B):
```cpp
for (int m = 0; m < nct; m++) {
    double corr_up = y_up[m]*dS_up_A[m] + y_up[m]*dS_up_A[m]*dS_up_A[m]/s_up[m];
    double corr_dn = y_dn[m]*dS_dn_A[m] + y_dn[m]*dS_dn_A[m]*dS_dn_A[m]/s_dn[m];
    double rmu_up_B = sigma_mu - s_up[m]*y_up[m] + corr_up;
    double rmu_dn_B = sigma_mu - s_dn[m]*y_dn[m] + corr_dn;
    rhs_B[m] = -(S[m] - T_flat[m]*W) + D_marg[m] * (rmu_up_B/s_up[m] - rmu_dn_B/s_dn[m]);
}
double corr_delta = y_delta*d_delta_A + y_delta*d_delta_A*d_delta_A/s_delta;
double rmu_delta_B = sigma_mu - s_delta*y_delta + corr_delta;
double margin_delta_center_B = 0.0;
for (int m = 0; m < nct; m++) {
    /* … same corr_up/corr_dn/rmu_*_B re-derivation … */
    margin_delta_center_B += w_kj[m]*(rmu_up_B/s_up[m] + rmu_dn_B/s_dn[m]);
}
for (int m = 0; m < nct; m++) rhs_B[m] *= D_jac[m];
ldlt_solve(N0.data(), static_cast<size_t>(nct), rhs_B.data());
for (int m = 0; m < nct; m++) dlambda_B[m] = D_jac[m] * rhs_B[m];
double utw_B = 0.0;
for (int m = 0; m < nct; m++) utw_B += u_vec[m]*dlambda_B[m];
double sm_coeff_B = (std::fabs(sm_denom) > 1e-300) ? (alpha_sm*utw_B/sm_denom) : 0.0;
for (int m = 0; m < nct; m++) dlambda_B[m] -= sm_coeff_B * v_vec[m];

std::fill(dX_B.begin(), dX_B.end(), 0.0);
for (int k = 0; k < st.K; k++)
    for (int c = 0; c < ct.M_cell; c++) {
        int g = ct.g_per_cell[k][c];
        if (g >= 0 && g < st.cat_counts[k])
            dX_B[c] += dlambda_B[cat_offset[k]+g];
    }
for (int c = 0; c < ct.M_cell; c++) dX_B[c] = D_eff[c] * dX_B[c];

double w_dot_dlambda_B = 0.0;
for (int m = 0; m < nct; m++) w_dot_dlambda_B += w_kj[m]*dlambda_B[m];
double d_delta_B = alpha_sm * (
    rmu_delta_B/s_delta + margin_delta_center_B
    - r_delta_stat - (y_delta/s_delta)*w_dot_dlambda_B );

for (int m = 0; m < nct; m++) dlambda[m] = dlambda_B[m];
for (int c = 0; c < ct.M_cell; c++) dX[c] = dX_B[c];
```

**After** — same structure, reduced-space, with ν correction:
```cpp
// RHS_B (margin) — only the reduced rows enter the Newton system
for (int nr = 0; nr < nct_red; nr++) {
    int m = red_to_full[nr];
    double corr_up = y_up[m]*dS_up_A[m] + y_up[m]*dS_up_A[m]*dS_up_A[m]/s_up[m];
    double corr_dn = y_dn[m]*dS_dn_A[m] + y_dn[m]*dS_dn_A[m]*dS_dn_A[m]/s_dn[m];
    double rmu_up_B = sigma_mu - s_up[m]*y_up[m] + corr_up;
    double rmu_dn_B = sigma_mu - s_dn[m]*y_dn[m] + corr_dn;
    rhs_B_red[nr] = -(S[m] - T_flat[m]*W) + D_marg[m] * (rmu_up_B/s_up[m] - rmu_dn_B/s_dn[m]);
}
// δ centering (full nct; δ stationarity row spans ALL margins, including reference)
double corr_delta = y_delta*d_delta_A + y_delta*d_delta_A*d_delta_A/s_delta;
double rmu_delta_B = sigma_mu - s_delta*y_delta + corr_delta;
double margin_delta_center_B = 0.0;
for (int m = 0; m < nct; m++) {
    double corr_up = y_up[m]*dS_up_A[m] + y_up[m]*dS_up_A[m]*dS_up_A[m]/s_up[m];
    double corr_dn = y_dn[m]*dS_dn_A[m] + y_dn[m]*dS_dn_A[m]*dS_dn_A[m]/s_dn[m];
    double rmu_up_B = sigma_mu - s_up[m]*y_up[m] + corr_up;
    double rmu_dn_B = sigma_mu - s_dn[m]*y_dn[m] + corr_dn;
    margin_delta_center_B += w_kj[m]*(rmu_up_B/s_up[m] + rmu_dn_B/s_dn[m]);
}

// Solve in reduced space (REUSE factored N_red — no refactor)
for (int nr = 0; nr < nct_red; nr++) rhs_B_red[nr] *= D_jac_red[nr];
ldlt_solve(N_red.data(), static_cast<size_t>(nct_red), rhs_B_red.data());
for (int nr = 0; nr < nct_red; nr++) dlambda_B_red[nr] = D_jac_red[nr] * rhs_B_red[nr];
double utw_B = 0.0;
for (int nr = 0; nr < nct_red; nr++) utw_B += u_red[nr]*dlambda_B_red[nr];
double sm_coeff_B = (std::fabs(sm_denom) > 1e-300) ? (alpha_sm*utw_B/sm_denom) : 0.0;
for (int nr = 0; nr < nct_red; nr++) dlambda_B_red[nr] -= sm_coeff_B * v_red[nr];

// ── ν correction (Phase B): reuse w_e_red, sm_denom from Phase A; refresh ute ─
double ute_B = 0.0;
for (int nr = 0; nr < nct_red; nr++) ute_B += u_red[nr] * w_e_red[nr];
// sm_coeff_e is already applied to w_e_red in Phase A; skip if reusing;
// simpler: re-solve w_e_red here for symmetry with Phase A and to stay numerically clean.
//
// PRACTICAL: e_red is fixed within an iteration → w_e_red identical between phases.
// Reuse the Phase-A w_e_red (already SM-corrected). schur_nu identical too.
const double schur_nu_B = schur_nu_A;
double eTdl_B = 0.0;
for (int nr = 0; nr < nct_red; nr++) eTdl_B += e_red[nr] * dlambda_B_red[nr];
const double dnu_B = (schur_nu_B > kSchurNuMin) ? (-r_nu - eTdl_B) / schur_nu_B : 0.0;
for (int nr = 0; nr < nct_red; nr++) dlambda_B_red[nr] -= dnu_B * w_e_red[nr];

// dX_B from dlambda_B_red + Δν · D_eff
std::fill(dX_B.begin(), dX_B.end(), 0.0);
for (int c = 0; c < ct.M_cell; c++) {
    double sum_dlam = 0.0;
    for (int k = 0; k < st.K; k++) {
        int g = ct.g_per_cell[k][c];
        if (g < 0 || g >= st.cat_counts[k]) continue;
        int nr = full_to_red[cat_offset[k] + g];
        if (nr >= 0) sum_dlam += dlambda_B_red[nr];
    }
    dX_B[c] = D_eff[c] * (sum_dlam + dnu_B);
}

double w_dot_dlambda_B = 0.0;
for (int nr = 0; nr < nct_red; nr++)
    w_dot_dlambda_B += w_kj[red_to_full[nr]] * dlambda_B_red[nr];
double d_delta_B = alpha_sm * (
    rmu_delta_B/s_delta + margin_delta_center_B
    - r_delta_stat - (y_delta/s_delta)*w_dot_dlambda_B );

// Alias for downstream slack/dual update code (line search, dY_*, X update)
for (int c = 0; c < ct.M_cell; c++) dX[c] = dX_B[c];
```

The `delta_S`, `dS_up`, `dS_dn`, `dY_up`, `dY_dn`, `dY_lo`, `dY_hi`, `dY_delta`, line search, primal/dual update blocks at current lines 670–778 are **unchanged** — they all read from `dX[]`, `delta_S[]`, slacks, and duals (full-nct). The reduced Newton system is invisible to them.

### Step 3.6 — Remove the dead diagnostic block

Delete lines 322–351 of the current file (the `verbose >= 2 && iter == 0 && nct_red > 0` block that builds a separate `N_red` just for the diagnostic). The new `schur_nu`/`dnu` log line in Step 3.4 supersedes it and prints richer information from values already in scope.

### Step 3.7 — Remove the `n_comp == 0` degenerate branch

**Safety verification (mandatory before deletion):**
```bash
grep -n "n_comp" src/chebyshev.cpp
```
Confirm: (a) `n_comp` is assigned exactly once as `n_comp = 2*ct.M_cell + 2*nct + 1`, and (b) `ct.M_cell >= 1` and `nct >= 1` always (they fail build_cell_table earlier if 0), so `n_comp >= 5`. The `if (n_comp == 0)` branch is unreachable. The existing comment at lines 354–355 confirms this. If the grep shows any other write to `n_comp` — **halt and report BLOCKED.**

Lines 356–490 of the current file are dead code (`n_comp = 2*M_cell + 2*nct + 1 >= 1` always). The current comment at 354–355 acknowledges this. Drop the `if (n_comp == 0) { … } else { …Phase A/B… }` wrapper; keep the Phase A/B body unindented. Cuts ~135 lines of unreachable code, simplifies T3 review.

### Compile gate
```bash
R CMD INSTALL --preclean .
```
Expected: `* DONE (leafblower)`. If compile fails, **halt and read errors**; do not retry blindly (per global rule).

### Test gate (focused)
```bash
Rscript -e 'testthat::test_dir("tests/testthat", filter="chebyshev")'
```
Expected: 4 PASS. `T_cheby_warm` must remain GREEN (K=3 max_err <= raking + margin).

### Acceptance for T3

- Source compiles clean.
- `grep -c "dnu_A\|dnu_B\|schur_nu_A\|w_e_red\|e_red\[" src/chebyshev.cpp` ≥ 12 (correction code present in both phases).
- `grep -c "N0\|dlambda_A\[" src/chebyshev.cpp` == 0 (full-nct Newton vectors removed).
- All existing tests in `test-chebyshev.R` GREEN.
- Confidence: 80 (the SM stacking is correct on paper; numerical hazard is the ν step combining with σ — line search may shrink hard early. Mitigation: T4 monitors α_p, α_d magnitudes.)

---

## T4 — Verification: Stepstone K=9 convergence

**Ticket**: file as `chebyshev: T4 verify — Stepstone K=9 converges with ν fix`.

### Goal

Demonstrate the canonical Stepstone K=9 failure now converges (`status==0`) with `max_err` competitive against raking and greenkhorn.

### Step 4.0 — Committed K=4 testthat test (CI-safe, no parquet)

Add to `tests/testthat/test-chebyshev.R` (committed before running other steps):
```r
test_that("chebyshev converges on K=4 overlapping-margin problem (ν fix)", {
  set.seed(7L); n <- 2000L
  df <- data.frame(
    a = factor(sample(3L, n, replace = TRUE)),
    b = factor(sample(4L, n, replace = TRUE)),
    c = factor(sample(3L, n, replace = TRUE)),
    d = factor(sample(2L, n, replace = TRUE))
  )
  tgt <- list(
    a = setNames(rep(1/3, 3L), as.character(1:3)),
    b = setNames(rep(1/4, 4L), as.character(1:4)),
    c = setNames(rep(1/3, 3L), as.character(1:3)),
    d = setNames(rep(1/2, 2L), as.character(1:2))
  )
  r <- suppressWarnings(
    harvest(df, tgt, method = "chebyshev", max_iterations = 300L,
            attach_weights = FALSE, verbose = 0L)
  )
  res <- attr(r, "result")
  expect_equal(res$status, 0L,
    label = sprintf("K=4 chebyshev: status=%d max_err=%.2e", res$status, res$max_error))
  expect_lt(res$max_error, 1e-3,
    label = sprintf("K=4 chebyshev: max_err=%.2e", res$max_error))
})
```
This test fails with the current code (NOCONV on K=4) and passes after T3. Commit it as part of T4 Step 4.0 so CI tracks the fix.

### Step 4.1 — K=2 sanity (must not regress)

```bash
Rscript -e '
library(leafblower)
set.seed(42); n <- 1000L
df  <- data.frame(v1=factor(sample(5,n,TRUE)), v2=factor(sample(4,n,TRUE)))
tgt <- list(v1=setNames(c(0.4,0.3,0.15,0.1,0.05),as.character(1:5)),
            v2=setNames(c(0.4,0.3,0.2,0.1),as.character(1:4)))
r <- harvest(df, tgt, method="chebyshev", max_iterations=200L, attach_weights=FALSE)
res <- attr(r,"result")
cat("K=2 n=1000: max_err=",res$max_error,"iters=",res$iterations,"status=",res$status,"\n")
stopifnot(res$status == 0L, res$max_error < 1e-6)
'
```
Expected: `status= 0`, `max_err < 1e-6`, iters ≤ 50.

### Step 4.2 — Stepstone K=9

```bash
OMP_NUM_THREADS=1 Rscript -e '
suppressPackageStartupMessages({library(arrow);library(jsonlite);library(leafblower)})
df  <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet"); df$uuid <- NULL
tgt <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
tgt <- lapply(tgt, function(t){ t <- unlist(t); t/sum(t) })
for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
t0 <- proc.time()["elapsed"]
r  <- suppressWarnings(leafblower::harvest(df, tgt, method="chebyshev",
                                             max_iterations=500L, attach_weights=FALSE, verbose=0))
res <- attr(r,"result")
cat("chebyshev:", proc.time()["elapsed"]-t0, "s  iters=", res$iterations,
    "  max_err=", res$max_error, "  status=", res$status, "\n")
stopifnot(res$status == 0L)
stopifnot(res$max_error < 1e-3)
' 2>&1 | grep "chebyshev:"
```
Expected: `status= 0`, `max_err < 1e-3` (vs current 0.139 NOCONV).

### Step 4.3 — Test gate (full chebyshev suite)

```bash
Rscript -e 'testthat::test_dir("tests/testthat", filter="chebyshev")'
```
Expected: 4/4 PASS, including `T_cheby_K9` if Stepstone parquet present (asserts `chebyshev max_err <= greenkhorn`).

### Step 4.4 — Full regression

```bash
Rscript -e 'devtools::test()' 2>&1 | tail -3
```
Expected: `FAIL 0 | PASS >= 374`.

### Step 4.5 — Compile gate (final)

```bash
R CMD INSTALL --preclean .
```
Expected: `* DONE (leafblower)`.

### Step 4.6 — Sum-preservation spy (manual)

With `verbose=2`, run K=4 synthetic and confirm the per-iteration log shows `r_nu → 0` within ~5 iterations. If `r_nu` does **not** decay, the ν step is wired wrong; halt per `SPEC_FAILURE` rule, do not pivot.

### Acceptance for T4

- Stepstone K=9 `status==0`, `max_err < 1e-3`.
- K=2 small case `status==0`, `max_err < 1e-6`.
- All 4 chebyshev tests GREEN.
- Full test suite: 0 FAIL.
- `r_nu` log shows decay to `< 1e-6 · n_d` by iter 5 on K=4 synthetic.
- Confidence: 70 (Stepstone is the hardest case in the suite; if iter cap insufficient, increase to 800 — but if convergence rate is slow, that flags a deeper issue, not a cap issue. Halt and re-derive rather than bumping `kMaxIpm`.)

---

## Out of scope

- No change to `harvest()` R-side dispatch.
- No change to `LpVariant` enum, `ChebyshevResult` struct, or any header.
- No change to `compute_normal_equations` (full-nct variant) — kept for any other caller; the chebyshev caller switches to `_reduced`.
- No removal of the `cheb-renorm` branch (out-of-scope cleanup; file ticket if branch still exists).
- Algorithm dispatch / GP routing — unaffected.

## Risks and reservations

1. **Numerical conditioning at small `schur_nu`**: even with reference elimination, when most cells lie on box bounds (late iterations), `D_eff[c]` shrinks for active cells and `schur_nu` can collapse. Mitigation: `kSchurNuMin = 1e-8` gate falls back to `dnu = 0` (i.e., disable ν correction for that step), which is the current behavior anyway; no INFEAS spike. Confidence: 75 (untested in late-iteration boundary regime).
2. **Single-margin (K=1) edge case**: `nct_red = nct − 1`. `e_red` rows now sum to a non-degenerate vector. Should work but lacks coverage in the existing test corpus. Add a K=1 test to `test-chebyshev.R` if T1 reveals issues.
3. **Tight max-weight constraints (E1/E2 corpus)**: prior plan's E1/E2 assertions (`chebyshev max_err <= raking`) currently hold in `T_cheby_warm`. The fix tightens convergence; expect E1/E2 to still hold or improve. If E1/E2 *regress*, the ν correction is over-shooting on min-weight=0.2 / max-weight=5 boxes — investigate `dnu` magnitudes vs `α_p`.
4. **Kept-out: `kMaxIpm = 500`**: the comment at line 20 acknowledges K=9 NOCONV is iter-cap-independent; T4 verifies that the **fixed** code reaches `status==0` well below this cap. If T4 needs >500 iters, the fix is incomplete — do not raise the cap.

## Files changed

- `src/chebyshev.cpp` — main and only source change. Net delta likely −80 to −120 lines (Phase A/B reduced + dead `n_comp==0` branch + dead diagnostic removed) plus +60 lines of ν correction. Estimated final file: ~750 lines (vs current 824).

## Tickets (one per task — no bundling)

- T1 → ticket `chebyshev: T1 diagnostic — find K-threshold for current INFEAS/NOCONV`
- T2 → ticket `chebyshev: T2 math — formalize ν Schur-complement block elim`
- T3 → ticket `chebyshev: T3 impl — ν correction in both Mehrotra phases`
- T4 → ticket `chebyshev: T4 verify — Stepstone K=9 converges with ν fix`

Each ticket carries the bead-lock pattern from §9 of global rules: `Task [ν Schur block elim, reduced-space Newton] ! [renorm hack, schur_nu<eps silent guard, full-nct e_vec SM]`.
