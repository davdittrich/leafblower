# Weighting Methods — Comparative Overview

Eight distinct calibration ("weighting") solvers live in the C++ core (`src/`). Names in the codebase are *mostly* faithful, but the mathematical identity below is read from the **actual update rules**, not the labels. `RK_ALG_AUTO = 0` is a **dispatcher**, not a method; enum slot `2` is removed (was LBFGSB).

All eight solve the same shape of problem — **adjust design weights so every marginal total hits its target** — and differ in (a) the *distance/objective* they keep weights close under, (b) the *algorithm class*, and (c) *how hard per-cell bounds are enforced*.

## What they all share

- **Same input**: design masses `X_init`, margin incidence `A`, targets `T_kj`, optional per-cell box `[L_c, U_c]`.
- **Same exit contract**: `finalize_weights` (`calib_dispatch.hpp`) enforces **Σw = n first, then `bounds_mode`** dispatch (`cell` = count violations only; `unit` = per-cell water-fill). Order is normalize→bounds; renormalising *after* water-fill is forbidden.
- **Same fixed point** *(KL/IPF family)*: IEPPA, Raking, Sinkhorn, Greenkhorn, and Newton-KL all converge to the maximum-entropy (KL) calibration solution when feasible — they are different *algorithms for the same optimum*. GREG (χ²), Logit (logit distance), and Chebyshev (minimax) target **different** optima.
- **Shared infrastructure**: greedy/round-robin scheduler, optional SRAA-m Anderson acceleration (stateless solvers only), optional homotopy multi-level cascade over `max_weight`, `CellMetrics` (errRp/KL/χ²), `BestIterTracker`.

## Where they differ

| Method | Enum | Algorithm class | Distance / objective | Order | Bounds handling | Per-iter cost | Stateless? (SRAA-able) | Convergence test | Guarantee grade |
|--------|------|-----------------|----------------------|-------|-----------------|---------------|------------------------|------------------|-----------------|
| **IEPPA** | 1 | IPF / Sinkhorn + SOR | KL (I-divergence) | 1st | finalize (SOFT: ADMM penalty) | cheap | ✅ | metric improvement / plateau | Moderate |
| **Raking** | 3 | Cyclic IPF + water-fill | bounded KL | 1st | inline water-fill (hard, each step) | cheap-medium | ✅ | metric improvement / plateau | Strong (conv.) |
| **Sinkhorn** | 4 | Sinkhorn + Dykstra | KL on margins+box | 1st | Dykstra log-correction + bisection | medium | ❌ (stateful) | fixed point (no improvement rule) | Strong |
| **Greenkhorn** | 9 | Greedy coordinate IPF | max marginal residual (errRp) | 1st | inline clamp each step | cheap/step | ✅ | metric / max-residual plateau | Moderate→Strong |
| **Chebyshev** | 5 | Mehrotra interior-point LP | minimax (L∞) margin error | 2nd | LP box inequalities | heavy | ❌ | complementarity μ < 1e-6 | Strong theory / Moderate practice |
| **GREG** | 6 | Active-set Newton (QP) | χ² (Deville–Särndal linear) | 2nd | active set (clamp + KKT release) | medium-heavy | ❌ | no active-set change (KKT) | Strong |
| **Newton-KL** | 11 | TSVD trust-region Newton (dual) | KL via smooth dual `log Z − T·λ` | 2nd | finalize | heavy (eig/iter) | ❌ | dual gap ‖∇g‖∞ < tol | Strong (in regime) |
| **Logit** | 10 | Newton + Armijo (GLM) | logit distance (Deville–Särndal) | 2nd | **by construction** (logistic link) | medium-heavy | ❌ | residual plateau | Strong est. / Moderate conv. |

### Reading the table

- **Order**: 1st = first-order (linear convergence, cheap steps); 2nd = second-order (super-linear/quadratic near optimum, expensive steps). The 1st-order block (IEPPA/Raking/Sinkhorn/Greenkhorn) trades iteration count for per-step cheapness; the 2nd-order block (Chebyshev/GREG/Newton-KL/Logit) trades per-step linear algebra for few iterations.
- **Objective families**:
  - **KL / maximum entropy** → IEPPA, Raking, Sinkhorn, Greenkhorn, Newton-KL (positive weights "for free").
  - **χ² (quadratic)** → GREG (interpretable linear weights, but can go negative without the box).
  - **logit distance** → Logit (bounds for free via the link).
  - **minimax (L∞)** → Chebyshev (minimises the *worst* margin error, ignores closeness-to-design).
- **Bounds enforcement** is the sharpest design axis: *finalize-only* (IEPPA, Newton-KL) → *inline projection* (Raking water-fill, Sinkhorn Dykstra, Greenkhorn clamp) → *active set* (GREG) → *parametrised away* (Logit logistic link) → *LP constraint* (Chebyshev).
- **Stateless ⇒ SRAA-able**: only the solvers whose sweep is a pure fixed-point map (no carried correction vectors) can be wrapped by SRAA-m Anderson acceleration: IEPPA, Raking, Greenkhorn. Sinkhorn (Dykstra state), the Newton/IPM methods (iterate-dependent factorisations) cannot.

## Selection guidance (practical)

| If you need… | Prefer | Why |
|--------------|--------|-----|
| Robust default, severe skew | **IEPPA** | Cheap, log-stable, converges where Newton/IPM stall. |
| Hard per-cell bounds, KL-optimal | **Raking** or **Sinkhorn** | Inline bounded-KL projection; Raking is SRAA-able, Sinkhorn is Dykstra-exact. |
| Minimise the *worst* margin error | **Chebyshev** | Only true minimax solver. |
| Classical interpretable linear weights | **GREG** | Deville–Särndal χ²; easy variance estimation. |
| Strict bounds with smooth weights | **Logit** | Bounds by construction; no clamping kinks. |
| Fast second-order on well-conditioned KL | **Newton-KL** | Quadratic convergence in the non-saturating regime. |
| Best worst-case iteration complexity (OT scaling) | **Greenkhorn** | Greedy Õ(n/ε²) advantage. |

## Critical-appraisal summary (GRADE-style)

- **Strongest guarantees**: Sinkhorn, Raking, GREG — rest on textbook convex/projection/QP theorems faithfully implemented and (for Sinkhorn's bracketing) re-proved in source.
- **Strong-but-scoped**: Newton-KL (proven for the non-saturating regime), Logit (estimator + bound-feasibility proven; convergence conditional on feasibility), Chebyshev (LP/IPM proven, but source documents non-convergence on dense K≥9 overlapping margins).
- **Moderate**: IEPPA, Greenkhorn — base algorithm proven, but production embellishments (SOR auto-adaptation, α-damping, greedy tie-handling, Anderson) are empirically guarded rather than formally proven.

> **Honest caveat (no method has a universal convergence proof here).** Every solver carries engineering guards (stall reverts, dual-explosion guards, saturation early-exits, iteration caps) precisely because the clean theorems assume feasibility and conditioning the production inputs do not always meet. The guarantees above are for the *idealised* problem; the guards are what make the *real* inputs terminate.

## Per-method documents

| Method | Doc |
|--------|-----|
| IEPPA (+ SOFT) | [ieppa.md](ieppa.md) |
| Raking | [raking.md](raking.md) |
| Sinkhorn | [sinkhorn.md](sinkhorn.md) |
| Greenkhorn | [greenkhorn.md](greenkhorn.md) |
| Chebyshev | [chebyshev.md](chebyshev.md) |
| GREG | [greg.md](greg.md) |
| Newton-KL | [newton_kl.md](newton_kl.md) |
| Logit | [logit.md](logit.md) |

## Bibliography

All citations in the per-method docs resolve against [`references.bib`](references.bib) (BibTeX, ~45 entries). Each method doc carries its own `## References` section listing the keys it cites; shared seminal sources (Csiszár 1975, Sinkhorn–Knopp 1967, Deville–Särndal 1992, …) are single entries reused across methods.

## Enum reference (`src/leafblower.h`)

```c
RK_ALG_AUTO       = 0   // dispatcher (resolves to a concrete method)
RK_ALG_IEPPA      = 1   // IPF/Sinkhorn + SOR  (KL)
// 2 = removed (was LBFGSB)
RK_ALG_RAKING     = 3   // cyclic IPF + water-fill (bounded KL)
RK_ALG_SINKHORN   = 4   // Sinkhorn + Dykstra (KL on margins+box)
RK_ALG_CHEBYSHEV  = 5   // Mehrotra interior-point (minimax LP)
RK_ALG_GREG       = 6   // active-set Newton (χ², Deville–Särndal)
RK_ALG_IEPPA_SOFT = 8   // IEPPA + ADMM soft capacity
RK_ALG_GREENKHORN = 9   // greedy coordinate-descent IPF (max residual)
RK_ALG_LOGIT      = 10  // logit-link Newton (Deville–Särndal 1992)
RK_ALG_NEWTON_KL  = 11  // TSVD trust-region Newton on dual KL
```
