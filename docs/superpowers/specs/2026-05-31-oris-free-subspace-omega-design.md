# ORIS Free-Subspace θ₂ for Box-Constrained Over-Relaxation — Design

**Date:** 2026-05-31
**Status:** Design (awaiting implementation plan)
**Predecessor:** leafblower-mj1p.2 (spectral optimal-ω) closed NO-GO — Lehmann's global θ₂ estimator fails on bounded ORIS.
**Research notebook:** NotebookLM `1e3036a1-fcbb-4d05-bc8e-820854f59d8e` (do not delete).

---

## Mechanism
Free-subspace spectral estimate of the over-relaxation factor ω: estimate the local
convergence rate from error reduction measured over **free (unclamped) cells only**,
feed it into Lehmann's `ω = 2/(1+√(1−θ₂))`.

## Forbidden
- Global/all-cell residual ratio driving θ₂ (the mj1p.2 bug).
- Touching the fixed point: ω changes only the iteration path, never the converged weights/objective.
- Generic LCP solver substitution (Siconos/PATH) — wrong layer, replaces ORIS.
- SRAA adapt path in this cycle (deferred; confounds with Anderson — see Out of Scope).
- ω ≥ 2 (Thibault strict-`<2` boundary; clamp at 1.99).

## Audit
- SageMath artifact computes the true `ρ(M_II)` on a synthetic and is the GO gate.
- Independent re-run of all benchmark numbers by the orchestrator (subagent numbers not trusted).

---

## 1. Problem & root cause

The Lehmann–von Renesse–Sambale–Uschmajew (2022, [2012.12562]) optimal over-relaxation
factor is Young's SOR result:

```
ω_opt = 2 / (1 + √(1 − ρ(M)²))
```

where `ρ(M)` is the second-largest eigenvalue of the linearized Sinkhorn/IPF Jacobi
iteration matrix `M` at the fixed point (the largest is identically 1 from scaling
indeterminacy).

**Box constraints break the global estimate.** With per-cell capacity clamping
`L_c ≤ X[c] ≤ U_c`, the literature consensus (NotebookLM synthesis *Symbolic Derivation
and Numerical Verification of Optimal Over-Relaxation in Box-Constrained Matrix Scaling*;
Lehmann 2022) is:

- **No global closed-form ω_opt exists** for the box-constrained case.
- Error propagates **only through the inactive (free) coordinate block** `I`. The active
  set (cells pinned at `L_c`/`U_c`) is a hard wall.
- The correct, active-set-dependent optimum is:

  ```
  ω_opt(I) = 2 / (1 + √(1 − ρ(M_II)²))
  ```

  where `M_II` is the principal submatrix of `M` on the free coordinates.

ORIS (mj1p.2) fed the **global** residual ratio into the formula. Pinned cells contribute
zero error reduction → the global ratio → 1 → θ₂ → 1 → ω → 1.9 → over-relaxes the free
coordinates → oscillation damp fires → slower than fixed ω=1.5. Confirmed empirically:
stepstone mw=5 gave spectral 40 iters vs fixed 30.

## 2. ORIS mapping (where the active set lives)

ORIS operates on **cells**, not per-obs weights (`src/oris.cpp`):

- Per-cell capacity multiplier `W[c]`, bounded `[L_c/X̃, U_c/X̃]` (`oris.cpp:27`).
- In-loop mass-preserving water-fill clamp, `oris.cpp:997–1038`:
  `X[c] = clamp(X̃[c]·W[c], L_cell[c], U_cell[c])`, one capacity BCD block per outer iter.
- `L_cell[c] = min_weight · n_per_cell[c]`; `U_cell[c]` recomputed per homotopy level from
  `current_max_weight` (`oris.cpp:158–161, 516–519`).

**Active set definition:** cell `c` is *pinned* when `X[c]` sits within `kPinTol` of
`L_cell[c]` or `U_cell[c]` after the clamp block. All other cells are **free**.
`kPinTol` derived from the clamp's own bisection tolerance (relative, ~1e-9), not a magic
number.

**Current residual that drives θ₂ (flat path):** `errRp_k = max_j |S_lin[j]/W_total −
target[k][j]|` (`oris.cpp:1069`), per-margin, computed over all cells.

## 3. The change (flat adapt path only)

Three edits, all in the flat (non-SRAA) adapt block (`oris.cpp:~1556–1569`) plus the
free-cell residual computation:

1. **Free-cell residual.** Alongside the existing per-margin `errRp_k`, compute a
   free-subspace error measure `errF_k` = the same residual restricted to free cells
   (cells contributing to margin `k`'s category `j` that are not pinned). Pinned cells
   excluded from the max/aggregation.

2. **Lag-2 EMA ratio.** Estimate `β̂² = √(errF_p / errF_{p−2})` (lag-2, the operator-
   Sinkhorn standard per Soma–Uschmajew, cancels the alternating-sign transient that makes
   lag-1 noisy), then EMA-smooth (reuse the existing `sor_theta2_ema` machinery with the
   lag-2 numerator). Requires storing residual history at lag 2 (`errF_prev`,
   `errF_prev2`).

3. **ω from free θ₂.** `ω = omega_from_theta2(θ₂_ema, kSorSpectralCeiling=1.99)`. The
   sign-flip oscillation damp (`kSorOscillationDamp`) still fires on divergence and pulls ω
   down. If a margin has *no* free cells (fully pinned), θ₂ is uninformative → fall back to
   ω=1 (plain IPF) for that margin, not the 1.99 ceiling.

`omega_mode_id=2` (spectral, already wired in mj1p.2) gains the free-subspace behavior.
No new ABI field — `sor_theta2_ema` and the lag-2 history are solve-local vectors.

## 4. SageMath verification artifact (the GO gate)

`docs/superpowers/derivations/2026-05-31-free-subspace-omega.sage` (or `.ipynb`):

1. Build a small synthetic box-constrained IPF (3×3 seed, prescribed margins, one cell
   forced to clamp at `U`).
2. Form the unconstrained Jacobi `M` symbolically; extract the free submatrix `M_II`;
   compute `ρ(M)` and `ρ(M_II)` exactly.
3. Iterate the constrained scaling. Confirm:
   - free-coordinate error-reduction ratio → `ρ(M_II)²`,
   - **global** ratio → 1 (reproduces the mj1p.2 bug),
   - `ω_opt(I)=2/(1+√(1−ρ(M_II)²))` beats fixed ω on iteration count.
4. Emit a short markdown result table committed alongside the script.

**GO** = SageMath confirms free-subspace ratio tracks `ρ(M_II)²` AND ω_opt(I) beats fixed
on the synthetic. **NO-GO** = it does not → stop, document, do not touch ORIS.

## 5. Validation (after implementation, only if SageMath GO)

| Fixture | Purpose | Pass criterion |
|---------|---------|----------------|
| SageMath synthetic | spectral match | free ratio → ρ(M_II)² |
| Slow **unconstrained** problem (small spectral gap; NOT stepstone mw≥50 which ties at 20 iters) | spectral beats fixed where ω has room | spectral iters < fixed iters |
| stepstone mw=5 (the NO-GO case) | bounded regression | spectral iters ≤ fixed (no longer 40 vs 30) |
| Fixed-point check, all fixtures | correctness | weights 1e-8, marginal_kl 1e-9 vs baseline |
| R `devtools::test()` + Python parity | no regression | 0 failures; rtol=1e-6 |

A construction step is needed: the **slow unconstrained fixture** does not exist yet
(stepstone is too easy unbounded). The plan must build one — a high-K problem with
near-conflicting margins and loose `max_weight` so unconstrained IPF takes ≫50 iters.

## 6. Software

| Tool | Role | Why |
|------|------|-----|
| **SageMath** | symbolic `M_II` spectral derivation + numerical verify | one tool for both; free; SymPy/NumPy-backed (synthesis report top pick) |
| NumPy/SciPy | iteration simulation inside the artifact | already a Python dep |
| Julia PATHSolver / Complementarity.jl | *optional* ground-truth active set | only if the synthetic's active set is ambiguous |
| Lean/Coq/Isabelle | — | rejected: convergence-rate proof is overkill here |

## 7. Out of scope (this cycle)

- **SRAA adapt path** (`oris.cpp:~857–898`): uses global `err_rp` (index 0) and interacts
  with Anderson acceleration — confounds the free-subspace test. Follow-up ticket; relates
  to leafblower-e65t.1 (SRAA×over-relaxation interplay).
- **Adaptive PSOR (Wolfe-condition step adaptation)** — viable fallback if free-subspace
  spectral NO-GOs again, but needs per-step gradient evaluations ORIS doesn't compute.
  Documented, deferred.
- `omega_max` / fixed-mode removal — only after a successful spectral gate (user directive).

## 8. Risks

- **Free-cell residual cost.** Must be O(M_cell) amortized into the existing residual pass,
  not a new sweep. If it adds a pass, the iteration-count win can be eaten by per-iter cost.
- **Active-set churn.** If the active set changes every iteration, θ₂(I) is non-stationary
  and the EMA lags. Mitigation: the EMA + lag-2 already smooth; if churn dominates, that is
  itself a NO-GO signal to document.
- **Fully-pinned margins.** Handled by the ω=1 fallback (§3.3).
