# Sinkhorn Adaptive Optimal Omega — Theory Derivation

**Ticket:** leafblower-e65t.4.1 (theory/derivation gate)
**Date:** 2026-06-03
**Status:** GO (see Section 6)

---

## 1. Power-Law ↔ Log-SOR Equivalence

### 1.1 Standard Sinkhorn update

The alternating scaling step on row scaling vector `u` is:

```
u^(l+1) = a / (K v^(l))          # full-step Sinkhorn
```

In log-domain with dual potential `α = log u`:

```
α^(l+1) = log a − log(K v^(l))   # full-step in log domain
```

### 1.2 Power-law relaxation

Applying the over-relaxed scaling:

```
u^(l+1) = u^(l) * (a / (K v^(l)) / u^(l))^ω
         = u^(l)^(1−ω) * (a / (K v^(l)))^ω
```

Taking logs:

```
α^(l+1) = (1−ω) α^(l) + ω · [log a − log(K v^(l))]
         = (1−ω) α^(l) + ω · ᾱ^(l)
```

where `ᾱ^(l) = log a − log(K v^(l))` is the full-step update.

**Conclusion:** `scale = ratio^ω` in the primal domain is **exactly** log-domain SOR on dual
potentials `α`. The exponent law in the probability/mass domain maps to a convex combination
in the log domain. The equivalence is exact, not approximate.

Formalized by Thibault et al. (2021) as over-relaxed Bregman projection:

```
log P^ω_{C_k}(γ) = (1−ω) log γ + ω log P_{C_k}(γ)
```

which they prove equals a standard SOR update on the dual objective:

```
α^(l+1) = (1−ω) α^(l) + ω · argmax_α E(α, β^(l))
```

---

## 2. Optimal ω Formula

### 2.1 Formula

```
ω_opt = 2 / (1 + √(1 − θ₂))
```

where `θ₂ = ρ_fullstep` is the **second-largest eigenvalue of the full-step Sinkhorn
iteration matrix** M (i.e., the asymptotic convergence rate of one alternating scaling
pass).

### 2.2 Definition of M

Let `M = diag(1/a) P* diag(1/b) P*ᵀ` (Lehmann 2022, §3). The spectral radius of M minus
the leading eigenvalue 1 gives `θ₂ = λ₂(M)`. This is the Gauss-Seidel full-step rate:

```
θ₂ = ρ_GS = ρ_J²
```

where `ρ_J` is the per-sweep Jacobi spectral radius and `ρ_GS` is the full alternating
Gauss-Seidel sweep rate. **θ₂ is NOT ρ_J itself** — it is the squared Jacobi rate.

### 2.3 Derivation (Gauss-Seidel vs. Jacobi)

Classical Young/SOR theory uses the Jacobi rate `ρ_J`:

```
ω_opt = 2 / (1 + √(1 − ρ_J²))
```

Sinkhorn alternating scaling is a **Gauss-Seidel** method, so its full-step rate satisfies
`ρ_GS = ρ_J²`. Substituting `ρ_J² = ρ_GS = θ₂`:

```
ω_opt = 2 / (1 + √(1 − θ₂))     [CORRECT]
```

The formula `2/(1+√(1−θ₂²))` in Lehmann 2022 Eq. 5 squares the full-step rate a second
time, yielding a sub-optimal ω. Both Thibault et al. (Eq. with `η = 1−ρ`) and
Soma-Uschmajew (Eq. with `β² = ρ`) confirm the un-squared formula.

**References:**
- Thibault, Chizat, Dossal, Papadakis (2021). *Overrelaxed Sinkhorn-Knopp Algorithm for
  Regularized Optimal Transport.* Algorithms 14(5):143. Eq. θ* = 2/(1+√η), η = 1−ρ.
- Lehmann, von Renesse, Sambale, Uschmajew (2022). *A Note on Overrelaxation in
  the Sinkhorn Algorithm.* Optimization Letters. arXiv:2012.12562. Eq. 5.
- Soma, Uschmajew (2024). ω_opt = 2/(1+√(1−β²)), ρ = β². arXiv:2410.14104.

---

## 3. θ₂ Estimator: Iterate-Change Ratio with Block Root

### 3.1 Scalar estimator

Every `I = kErrCheckInterval` sweeps, in mode 2, after burnin:

```
S_dX(m)  = Σ_{c: !is_pinned[c]} (X[c] − X_snapshot[c])²
ratio    = S_dX(m) / S_dX(m−1)          → ρ(M_II)^(2·I)
theta2   = ratio^(1/I)                   → ρ(M_II)²
omega    = 2 / (1 + sqrt(max(0, 1 − clamp(theta2, 0, 1−1e-9))))
X_snapshot ← X
```

**State required:** `X_snapshot[M_cell]` (same size as X), `S_dX_prev` (scalar, init +∞).

### 3.2 Justification

Near a fixed point, the error satisfies the linear recurrence:

```
e^(l+1)_I ≈ (1−ω) e^(l)_I + ω M_II e^(l)_I
```

For the unrelaxed case (ω=1), `e^(l+1)_I = M_II e^(l)_I`, so consecutive squared-norm
ratios converge to `ρ(M_II)²` per sweep. After `I` sweeps: ratio → `ρ(M_II)^(2·I)`.
Taking the `I`-th root recovers `ρ(M_II)²`.

### 3.3 ORIS reference implementation

This estimator is already implemented and SHIP-gated in ORIS (`oris.cpp`, omega_mode_id=2,
commit e18t.9 / `8a30e3f`). The Sinkhorn adaptation transfers the same S_dX mechanism to
the Sinkhorn sweep loop.

**Validation results (from e18t.7 derivation / ORIS):**
- Feasible 3×3: |iterate-change ratio − ρ²| < 0.05 ✓
- Infeasible 4×4: marginal plateaus at 1.25 (ratio → 1 → ω ≈ 2.0); iterate-change ratio → ρ² ✓
- Cadence I=10: block-root recovers ρ² (±0.05) ✓

---

## 4. Bounded Case: Free-Subspace Argument

### 4.1 Variable partition

For box-constrained Sinkhorn (capacity bounds, lower and upper), at a constrained fixed point the
cells partition into:

```
Active   A : cells clamped to lower or upper bound (is_pinned[c] = true)
Inactive I : free cells strictly in the interior    (is_pinned[c] = false)
```

### 4.2 Error propagation

Under local perturbations, the projection acts as a hard clamp on A, eliminating error
there instantaneously: `e^(l+1)_A = 0`. Error propagation is restricted to the free
subspace, governed by the principal submatrix `M_II`:

```
e^(l+1)_I ≈ (1−ω) e^(l)_I + ω M_II e^(l)_I
```

Because the free-cell error is decoupled from clamped cells once the active set stabilizes,
`ρ(M_II) < 1` governs convergence independently of the values at clamped boundaries.

### 4.3 Why marginal-residual fails on infeasible problems

When column targets are unreachable (column mass cap < column target), the marginal
residual of the constrained problem converges to a nonzero **infeasibility floor**:

```
R²_free(m) → C_infeasible > 0   for all m → ∞
```

Consequently:

```
ratio = R²_free(m) / R²_free(m−1) → 1   (plateau)
theta2 → 1
omega  → 2   (SOR divergence / oscillation)
```

This is the **e18t.3 / mj1p.2 failure mode**: stepstone mw=5, 500 iters vs 300 baseline
(stall). The iterate-change estimator is immune because `S_dX → 0` regardless of
feasibility once the free cells have converged.

### 4.4 Transfer to Sinkhorn

The same argument applies to the Sinkhorn iterate `X[c]` (joint distribution cell values).
Free cells (not at capacity bounds) converge at `ρ(M_II)`, and `S_dX` over free cells
tracks this rate faithfully even when global margin feasibility is not achieved.

---

## 5. NotebookLM Corroboration

**Notebook:** 1e3036a1 ("ORIS box-constrained optimal-omega (Lehmann SOR)")
**Conversation:** 81c6b45d-7d33-45f5-a324-1d853b60184e
**Date:** 2026-06-03

### Query 1: Power-law ≡ log-SOR + ω formula

> "Confirm: does applying scale = ratio^omega in Sinkhorn marginal scaling equal log-domain
> SOR on dual potentials? What is the exact optimal omega formula — is it
> 2/(1+sqrt(1-theta2)) or 2/(1+sqrt(1-theta2^2))? What is theta2 in terms of the full-step
> Sinkhorn iteration rate rho?"

**Key quotes from response:**

> "Yes, applying the power-law scale `(ratio)^omega` in the primal domain is exactly
> equivalent to log-domain Successive Over-Relaxation (SOR) on the dual potentials."

> "Thibault et al. formally define this over-relaxed Bregman projection in the log-domain
> as `log P^ω_{C_k}(γ) = (1−ω)log γ + ω log P_{C_k}(γ)`. They explicitly show that this
> corresponds exactly to a standard linear SOR algorithm on the dual objective."

> "[θ₂] is identically the full-step Sinkhorn rate ρ (θ₂ = ρ). ... Since θ₂ is the
> full-step Gauss-Seidel rate (θ₂ = ρ = ρ_J²), substituting it into the classical theory
> yields ω_opt = 2/(1+√(1−θ₂))."

> "By writing ω_opt = 2/(1+√(1−θ₂²)) in their Equation 5, Lehmann et al. erroneously
> square the full-step rate. Both Thibault et al. and Soma-Uschmajew confirm the correct
> un-squared relationship when using the full-step rate."

**Verdict on Query 1:** Power-law ≡ log-SOR CONFIRMED. Correct formula: `2/(1+√(1−θ₂))`
(not squared). Lehmann Eq. 5 has a known sign/squaring discrepancy; Thibault + Soma-Uschmajew agree.

### Query 2: Bounded free-subspace feasibility-agnostic

> "For box-constrained Sinkhorn (cells have upper bounds), does restricting the theta2
> estimation to free cells (not at lower or upper bounds) make the iterate-change ratio
> estimator feasibility-agnostic? Does the error on free cells still converge at the same
> rate rho(M_II) regardless of whether clamped cells are at infeasible margins?"

**Key quotes from response:**

> "Yes, restricting the estimation to the free cells makes the estimator feasibility-agnostic,
> because the error on the free cells will indeed continue to converge at the rate
> ρ(M_II) regardless of whether the clamped cells force an infeasible margin."

> "Under local perturbations, the projection operator acts as a hard clamp on the active set,
> which ensures that the local error is eliminated (e^(t+1)_i = 0) for all clamped variables.
> Consequently, the local propagation of error is restricted entirely to the free variables."

> "the free cells continue to converge at the exact rate ρ(M_{II}) independently of the
> actual values at the clamped boundaries, confirming that tracking the iterate-change ratio
> strictly over the free coordinates correctly isolates the convergence rate even under
> global constraint infeasibility."

**Verdict on Query 2:** Free-subspace θ₂ estimator feasibility-agnostic CONFIRMED.
Error on free cells governed exclusively by `ρ(M_II)` once active set stabilizes.

---

## 6. GO/NO-GO Verdict

| Criterion | Evidence | Result |
|---|---|---|
| Power-law scale = ratio^ω equals log-domain SOR | Thibault 2021 + NotebookLM Q1 | **CONFIRMED** |
| Correct formula: `2/(1+√(1−θ₂))` (not squared) | Thibault + Soma-Uschmajew + NotebookLM Q1 | **CONFIRMED** |
| θ₂ = ρ_fullstep (second eigenvalue of full-step M) | Lehmann 2022 §3 + NotebookLM Q1 | **CONFIRMED** |
| Free-subspace iterate-change feasibility-agnostic | NotebookLM Q2 | **CONFIRMED** |
| Marginal-residual failure mode documented | ORIS e18t.3/mj1p.2 empirical + NotebookLM Q2 | **CONFIRMED** |
| Published counterexample to iterate-change estimator | None found | — |
| Published counterexample to power-law ≡ log-SOR | None found | — |

**DECISION: GO**

All transfer theorems valid. No blocking counterexample found. The Sinkhorn adaptive omega
implementation (e65t.4.2 through e65t.4.6) may proceed.

**Note on Lehmann Eq. 5:** The squaring discrepancy (θ₂ vs θ₂²) noted by NotebookLM means
the ORIS implementation (which uses `2/(1+√(1−theta2))` with `theta2 = ρ_fullstep`) is
correct and consistent with Thibault 2021. The Sinkhorn implementation must use the same
un-squared formula.

---

## 6. POST-HOC CORRECTION (e65t.4.5 empirical re-gate, 2026-06-03)

The GO verdict above was **empirically refuted** for Sinkhorn. The estimator is sound
(it ships in ORIS, omega_mode_id=2, e18t.9) but does NOT transfer to leafblower's
Sinkhorn solver. Verified NO-SHIP on ORIS's own shipping fixtures.

### 6.1 What was verified (self-run, benchmarks/e65t4_stepstone_results.txt)

| fixture (ORIS shipping conditions) | baseline | fixed_1.4 | adaptive |
|---|---|---|---|
| stepstone mw5 (the fixture ORIS shipped on) | 60 | 60 | **220 (+267%)** |
| adversarial unbounded (n=5000,K=8 inverted) | 100 | 50 | 220 (+120%), 1 NOCONV flip |

ORIS adaptive on stepstone mw5: 140→**50** (win). Sinkhorn adaptive: 60→**220** (loss).
Same fixture, same estimator, opposite outcome.

### 6.2 Why §1 (power-law ≡ log-SOR) does not hold for THIS Sinkhorn

§1 proves the equivalence for the **unprojected alternating-scaling skeleton**.
leafblower's `sinkhorn_solve` is scaling **+ capacity bisection + Dykstra correction
vectors `a[c]`** carried across iterations. Under active projection/correction (the
bounded AND adversarial regimes), the equivalence breaks.

### 6.3 Mechanism (LBW_SK_OMEGA_TRACE, stepstone, n_pinned=3/5980 — UNBOUNDED)

theta2_ema never settles: reads 0.96–0.99 → prescribes ω≈1.88 (clamped to 1.8 ceiling)
→ power-law over-relaxation overshoots → S_dX grows → gate-6/gate-9 reset to ω=1.0 →
brief convergence → theta2 high again → overshoot. **Limit cycle ω∈{1.0,1.8}, never
settles.** ORIS settles at a stable omega_mean=1.45.

**Root cause:** the iterate-change observable S_dX conflates (1) genuine SOR slowness
(ω accelerates) with (2) Dykstra/bisection correction grind (ω destabilizes). ORIS
water-fill is stateless → only source (1). Sinkhorn has both → estimator reads combined
slowness as θ₂≈0.99, prescribes ceiling ω, overshoots the part ω cannot accelerate.
n_pinned=3/5980 proves this is the projection/correction dynamics, NOT the bound clamp.

### 6.4 Also corrected

§3.3's "SHIP-gated in ORIS (e18t.9)" is **accurate** (verified harvest.R:970 — e18t.5
NO-SHIP was rescinded, v2 iterate-change ships as omega_mode_id=2 default). The transfer
failure is Sinkhorn-specific, not an ORIS-estimator defect.

**REVISED DECISION: NO-SHIP.** Sinkhorn is the wrong host operator for the iterate-change
estimator. Closed-form θ₂ remains unobtainable (circular unconstrained / ill-posed
piecewise-linear bounded); online estimation works only for stateless-projection solvers
(ORIS), not for Dykstra-corrected Sinkhorn. Fixed ω likewise does not generalize past
near-identity targets (e65t.3 GO is a toy-regime artifact — see e65t.3 flag ticket).
