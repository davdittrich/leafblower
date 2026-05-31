# ORIS Free-Subspace θ₂ for Box-Constrained Over-Relaxation — Design

**Date:** 2026-05-31
**Status:** Design (rev 2 — addresses design-review-gate iteration 1)
**Predecessor:** leafblower-mj1p.2 (spectral optimal-ω) closed NO-GO — Lehmann's *global* θ₂ estimator was slower (40 iters) than fixed ω=1.5 (30 iters) on bounded stepstone.
**Research notebook:** NotebookLM `1e3036a1-fcbb-4d05-bc8e-820854f59d8e` (do not delete; 191 sources).

---

## 0. Deliverable framing (PRIMARY vs CONDITIONAL)

The user's request is *"derive the Lehmann formula for box-constraint cases."* The derivation is
therefore the **primary deliverable**, and SOR is off-by-default in this library
(`sor=NULL`), so a shipped solver feature is **not** assumed. Two phases:

- **Phase 1 (PRIMARY, always done): SageMath derivation + documented finding.** A standalone
  symbolic+numerical artifact proving (or refuting) that the box-constrained optimal ω is
  governed by the free-subspace spectral radius `ρ(M_II)`. Output: committed `.sage`/`.ipynb`
  + a result note. This is the user's actual ask and is expected-value-positive regardless of
  whether ORIS ever changes.
- **Phase 2 (CONDITIONAL): ORIS flat-path implementation.** Only attempted if Phase 1 GOs AND
  a buildable slow-unconstrained fixture (§5) shows the estimator beats fixed ω. If either
  gate fails, the deliverable is the documented derivation + a recorded negative; ORIS is left
  unchanged and `omega_mode_id` default reverts to fixed (mode 1).

This reframing answers PM-BLOCKING-2 / CTO-note-4: value is the derivation; the code change
must earn its way in.

## Mechanism
Free-subspace spectral estimate of ω: estimate the local convergence rate from error
reduction over **free (unclamped) cells only**, feed `ω = 2/(1+√(1−θ₂))`.

## Forbidden
- Global/all-cell residual ratio driving θ₂ (the mj1p.2 bug).
- Changing the fixed point: ω changes only the iteration path.
- Generic LCP solver substitution (Siconos/PATH) — wrong layer.
- SRAA adapt path (deferred; confounds with Anderson — §7).
- ω ≥ 2, or `√` of a negative (Thibault strict `<2`; clamp before the formula — §3.4).

## Audit
- SageMath artifact computes the true `ρ(M_II)` and is the Phase-1 GO gate.
- Orchestrator independently re-runs every benchmark number (subagent numbers not trusted).

---

## 1. Problem & root cause (literature-settled)

Lehmann–von Renesse–Sambale–Uschmajew (2022, arXiv:2012.12562) give Young's optimal SOR factor
`ω_opt = 2/(1+√(1−ρ(M)²))`, where `ρ(M)` is the second-largest eigenvalue of the linearized
Sinkhorn/IPF Jacobi iteration matrix `M` at the fixed point (largest is 1 from scaling
indeterminacy).

**Settled literature finding** (NotebookLM synthesis *Symbolic Derivation and Numerical
Verification of Optimal Over-Relaxation in Box-Constrained Matrix Scaling*; direct query of the
191-source notebook): **no closed-form box-constrained ω_opt exists.** With per-cell clamping
`L_c ≤ X[c] ≤ U_c`, error propagates only through the **inactive (free) coordinate block** `I`;
the active set is a hard wall. The active-set-dependent optimum is:

```
ω_opt(I) = 2 / (1 + √(1 − ρ(M_II)²))
```

`M_II` = principal submatrix of `M` on the free coordinates. Soma–Uschmajew 2024
(arXiv:2410.14104) and the 2024/26 "Numerically stable variants of overrelaxation for operator
Sinkhorn" paper were checked: **neither addresses box/inequality constraints** — so this
derivation is genuinely new work, not a re-import.

mj1p.2 fed the *global* residual ratio into the formula; pinned cells give zero error reduction
→ global ratio → 1 → θ₂ → 1 → ω → 1.9 → over-relaxes the free coords → oscillation → slower
than fixed. (This paragraph makes the finding self-contained per PM-note-3.)

## 2. ORIS code reality (corrected per Architect iteration 1)

**Verified against `src/oris.cpp` (the rev-1 citations were wrong; these are corrected):**

- The two adapt paths are **separate branches**, not one body:
  - **SRAA branch** opens at `oris.cpp:832` (`if (sraa_active_lvl)`). The mass-preserving
    water-fill clamp lives **inside it** at `oris.cpp:997–1043`. This is the rev-1 mis-citation:
    that clamp is SRAA-only and does **not** run on the flat path.
  - **Flat (non-SRAA) branch** is the `else`, the flat adapt block at `oris.cpp:1599–1612`,
    using an **element-wise** clamp, not the SRAA water-fill.
- Per-margin residual `errRp_k = max_j |S_lin[j]/W_total − target[k][j]|` at `oris.cpp:1075`.
- Existing solve-local SOR state (`sor_prev_errRp`, `sor_prev_decreasing`, `sor_theta2_ema`)
  declared at `oris.cpp:~402–403`.
- Per-cell bounds `L_cell[c] = min_weight·n_per_cell[c]` (`:158–161`); `U_cell[c]` recomputed
  per homotopy level (`:516–519`).

**Consequence:** the active-set identification must be **re-derived for the flat path**; the
SRAA water-fill mechanism is NOT reusable. (Architect-BLOCKING-1/2.)

## 3. Phase-2 change (flat adapt path only) — fully specified

### 3.1 Active set is NEW computation (not a read of existing state)
There is no persistent `is_pinned[]` vector today. Add an O(M_cell) determination of
`is_pinned[c]`, folded into the existing hot cell loop (no separate sweep):
`is_pinned[c] = (X[c] ≥ U_cell[c]·(1−kPinTol)) || (X[c] ≤ L_cell[c]·(1+kPinTol))`.
(Architect-BLOCKING-3.)

### 3.2 `kPinTol` is an explicit relative tolerance
The flat path is element-wise clamp (no bisection), so "derive from bisection tol" does not
apply. Define `kPinTol = 1e-9` (relative), applied as above against `U_cell`/`L_cell`.
(Architect-BLOCKING-5.)

### 3.3 Free-cell residual via parallel accumulator (no second pass)
In the same hot loop that builds `S_lin[j]`, maintain a parallel `S_lin_free[j]` accumulator
gated on `!is_pinned[c]`, plus a per-margin free-cell count `n_free_k`. The free residual is
`errF_k = max_j |S_lin_free[j]/W_total_free − target[k][j]|` over categories with free cells.
A **second pass is forbidden**; the cost stays O(M_cell). (Architect-BLOCKING-4, risk §8.)

### 3.4 Estimator + all guards (Security iteration 1)
Per margin `k`, with lag-2 history `errF_prev[k]`, `errF_prev2[k]` (solve-local per-margin
`std::vector<double>(st.K)` alongside `sor_prev_errRp` at `oris.cpp:~402`):

Order of operations (each guard is a precondition, not a post-fix):
1. **Free-set gate:** if `n_free_k == 0` → margin fully pinned → `ω_k = 1`, skip the rest.
   (Security-BLOCKING-3.)
2. **Warm-up gate:** if `iter < sor_burnin` (default 20) or lag-2 history not yet filled →
   `ω_k = 1`. (Security-note-5.)
3. **Converged-denominator gate:** if `errF_prev2[k] < kResidFloor` (= 1e-12) → margin
   effectively converged → `ω_k = 1`, do not divide. (Security-BLOCKING-1.)
4. **Lag-2 ratio:** `β² = √(errF_p / errF_prev2[k])` (operator-Sinkhorn standard; lag-2 cancels
   the alternating-sign transient). EMA-smooth via existing `sor_theta2_ema[k]`.
5. **Clamp before sqrt:** `θ₂ = clamp(sor_theta2_ema[k], 0, 1−1e-9)` — prevents `√(neg)`/NaN
   when the residual increased (`β²≥1`). (Security-BLOCKING-2.)
6. **Formula:** `ω_k = 2/(1+√(1−θ₂))`, then `ω_k = min(ω_k, kSorSpectralCeiling=1.99)`.
7. **Oscillation damp** (`kSorOscillationDamp`) still fires on sign-flip and must read the
   **free-subspace residual** `errF`, not the global `errRp`, so damp and estimator agree on
   "diverging". (Security-note-5.)

### 3.5 Degenerate-state → ω table (complete API contract; Designer-note)
| State | ω used |
|-------|--------|
| Fully-pinned margin (`n_free_k=0`) | 1 |
| Cold start / before `sor_burnin` / lag-2 unfilled | 1 |
| `errF_prev2 < kResidFloor` (converged) | 1 |
| `β² ≥ 1` (residual grew) | clamped → ω near 1 (via θ₂→0 after clamp) |
| Normal | `2/(1+√(1−θ₂))`, capped 1.99 |

### 3.6 Convergence safety under active-set-dependent ω (Security-BLOCKING-4)
ω now depends on `I`, which depends on the weights — a feedback loop a static ω does not have.
The oscillation damp is the safety net, but the spec does not *assume* it dominates. Phase-2
validation (§5) MUST include a regression test: every fixture that CONVERGED under fixed ω must
still converge (not NOCONV/oscillate) under free-subspace ω. If any converts to NOCONV → Phase-2
NO-GO.

### 3.7 No new ABI; binding parity asserted (Designer-BLOCKING-3)
`is_pinned`, `S_lin_free`, `n_free_k`, lag-2 buffers are all solve-local; `CalibSorCfg` gains no
field. The free-subspace logic lives entirely in `oris.cpp` core, reached **identically** by
`r_bridge.cpp` and `c_api.cpp` (no per-binding wiring — the implementer MUST verify this, since
mj1p.1 had a c_api.cpp wiring miss). `omega_mode_id=2` is reused.

### 3.8 Observability diagnostic (Designer-BLOCKING-2)
Expose in the result struct (visible via `attr(w,"result")` in R / result object in Python):
`sor_omega_mean` (mean realized ω over adapted steps) and `sor_n_omega1_fallback` (count of
margin-steps that hit an ω=1 fallback gate). Without this, a silent all-fallback "third NO-GO"
would be invisible.

### 3.9 mode-2 default + docs/migration (Designer-BLOCKING-1, CTO-BLOCKING-3)
mj1p.2 shipped `omega_mode_id=2` (free-subspace's predecessor, "spectral global") as the
DEFAULT. Decision: **the redefined mode 2 ships as default ONLY if regression-neutral** on the
full suite (686 R + 8 Python). If it shifts any existing fixture, the default reverts to mode 1
(fixed) and free-subspace stays opt-in (`omega_mode_id=2L`). Either way, update in the SAME
commit: roxygen `@param` in `R/harvest.R`, the Python binding docstring, `docs/methods/oris.md`
(the guarantees-table NO-GO row becomes "superseded by free-subspace, §…"), and a NEWS entry
noting mode 2's behavioral change.

## 4. Phase-1 SageMath artifact (the PRIMARY deliverable + GO gate)

`docs/superpowers/derivations/2026-05-31-free-subspace-omega.sage` (+ committed result table):
1. Build a small synthetic box-constrained IPF (3×3 seed, prescribed margins, one cell forced
   to clamp at `U`).
2. Form unconstrained Jacobi `M` symbolically; extract free submatrix `M_II`; compute `ρ(M)`,
   `ρ(M_II)` exactly.
3. Iterate the constrained scaling; confirm: free-coordinate error ratio → `ρ(M_II)²`; global
   ratio → 1 (reproduces the bug); `ω_opt(I)` beats fixed ω on iteration count.

**Phase-1 GO** = all three confirmed. **NO-GO** = stop; commit the artifact + a documented
negative; do not start Phase 2.

## 5. Phase-2 fixtures + buildable generator (CTO-BLOCKING-1/2)

The load-bearing fixture is a **slow unconstrained** problem (stepstone ties at 20 iters when
loose — useless). Deterministic generator (so the test is falsifiable):

```
set.seed(20260531)
n = 5000; K = 8 margins, each binary (cardinality 2)
crossed factor structure with deliberately near-conflicting margins:
  - margins 1-4: probs c(0.85, 0.15)
  - margins 5-8: the SAME underlying latent flipped, probs c(0.15, 0.85)
  targets set to the opposite skew of the sample (forces many IPF sweeps)
  max_weight = 1000 (effectively unbounded — 0 cells pinned), min_weight = 0
EXPECTED: unconstrained ORIS (no SOR) takes ≫50 iters to the marginal_kl criterion.
```
The plan MUST first VERIFY this generator actually yields ≫50 unconstrained iters; if it does
not, tune skew/K until it does, and record the achieved baseline iteration count as the test's
asserted reference.

| Fixture | Purpose | Pass criterion |
|---------|---------|----------------|
| SageMath synthetic (Phase 1) | theory gate | free ratio → ρ(M_II)² |
| Slow unconstrained (generator above) | **SHIP gate** — ω has room | spectral iters < fixed iters |
| stepstone mw=5 (the NO-GO case) | bounded **no-regression** only | spectral iters ≤ fixed (NOT a ship gate; demoted per PM-BLOCKING-1) |
| All converged fixtures (§3.6) | convergence safety | none flip to NOCONV |
| Fixed-point check, all fixtures | correctness | weights 1e-8, marginal_kl 1e-9 vs baseline |
| R `devtools::test()` + Python parity | no regression | 0 failures; rtol=1e-6 |

**Ship / no-ship (distinct from the Phase-1 theory gate; CTO-BLOCKING-2):**
- **SHIP** = Phase-1 GO AND spectral beats fixed on the slow-unconstrained fixture AND mw=5 no
  regression AND no converged fixture flips to NOCONV.
- **NO-SHIP** = any of those fail → revert ORIS, set `omega_mode_id` default to mode 1, commit
  the SageMath derivation + a documented second negative. This is an acceptable, expected
  outcome, not a failure of the effort.

## 6. Software (NotebookLM-selected)

| Tool | Role | Why |
|------|------|-----|
| **SageMath** | symbolic `M_II` spectral derivation + numerical verify | one free tool for both; SymPy/NumPy-backed (synthesis report top pick) |
| NumPy/SciPy | iteration simulation in the artifact | already a Python dep |
| Julia PATHSolver / Complementarity.jl | *optional* ground-truth active set | only if the synthetic's active set is ambiguous |
| Lean/Coq/Isabelle | — | rejected: convergence-rate proof is overkill |

## 7. Out of scope (this cycle)
- **SRAA adapt path** (`oris.cpp:864–905`): global `err_rp` + Anderson interaction confounds the
  test. Follow-up; relates to leafblower-e65t.1 (SRAA×over-relaxation).
- **Adaptive PSOR (Wolfe-condition step adaptation,** nonneg-QP literature, c₁=0.89/c₂=0.95/
  λ₁=1.15/λ₂=1.4/ρ=0.85): viable fallback if free-subspace NO-GOs, but needs per-step gradient
  evals ORIS lacks. Documented, deferred.
- `omega_max`/fixed-mode removal — only after a successful spectral ship (user directive).

## 8. Risks
- **Free-cell residual cost:** must be the §3.3 parallel accumulator (O(M_cell)); a second pass
  would let per-iter cost eat the iteration win.
- **Active-set churn:** if `I` changes every iteration, θ₂(I) is non-stationary and the EMA lags;
  §3.6 test catches the dangerous case (converged→NOCONV); persistent churn that merely fails to
  beat fixed is a clean NO-SHIP, not a hazard.
- **Third negative:** this may, like mj1p.2, not beat fixed in any shippable regime. §0/§5 make
  that a documented, expected-value-positive outcome (the derivation stands).

## 9. Commit coherence (CTO-note-5)
- Phase-1: SageMath artifact + result table commit together as gate evidence.
- Phase-2 (if reached): ORIS change + new fixture + observability fields + roxygen/Python
  docstring/oris.md/NEWS updates commit together (CLAUDE.md: docs/tests with code).
