# SINKHORN — Log-Domain Sinkhorn with Dykstra Capacity Projection

> Enum: `RK_ALG_SINKHORN = 4`
> Source: `src/sinkhorn.cpp`, `src/sinkhorn.hpp`

## Overview

Entropy-regularised optimal-transport-style scaling: alternating **Sinkhorn** marginal rescales combined with **Dykstra's algorithm** for the per-cell capacity box. The distinguishing feature versus RAKING is *how* the box is enforced: Sinkhorn keeps a **log-domain correction vector** `a[c]` that accumulates the residual of each capacity projection (Dykstra correction), rather than projecting category-by-category inside the sweep.

The capacity step is solved by **bisection on a single log-scale multiplier** `μ`.

> **vs. IEPPA**: both reach the same KL/IPF fixed point on the margins; the difference is *how the box is enforced* (Sinkhorn = Dykstra correction every iteration, stateful, not SRAA-able; IEPPA = bounds at finalize, SOR + infeasibility damping, SRAA-able) — see the comparison table in [ieppa.md](ieppa.md#ieppa-vs-sinkhorn-dykstra). The source paper behind IEPPA (arXiv:2011.14312) names *DyKL — Dykstra-with-KL — as the baseline it improves on*, and DyKL ≡ the dual BCD that paper uses internally; this Sinkhorn+Dykstra solver is structurally that baseline.

## Mathematical formulation

### Objective

Two-constraint KL/entropy problem — match all margins *and* keep total mass / cell masses inside `[L_c, U_c]`:

```
min_X   KL(X ‖ X_init)
s.t.    Σ_{c∈(k,j)} X[c] = T_kj   ∀ k,j      (margin constraints)
        L_c ≤ X[c] ≤ U_c          ∀ c        (capacity box, via Dykstra)
        Σ_c X[c] = n                          (preserved by construction)
```

### Update rules (verified)

**Marginal step** (per margin k, category j):

```
scale[j] = T_kj · W / S_kj ,      X[c] *= scale[ cat(c) ]
```

**Capacity projection** — `bisect_capacity`: find `μ` such that

```
Σ_c clamp( X[c]·exp(a[c] + μ), L_c, U_c ) = target_mass
```

`f(μ) = Σ clamp(X·exp(a+μ),L,U) − target` is **strictly increasing**, so bisection is valid; the bracket `[lo,hi]` is derived from per-cell clamp thresholds so that `f(lo) ≤ 0 ≤ f(hi)` is guaranteed. A fast variant `bisect_capacity_fast` precomputes `exp(a[c])` to use one scalar `exp` per bisection step.

**Dykstra correction**: after projecting, the log-residual is accumulated into `a[c]` so the next pass corrects the previous projection's distortion — this is what makes the combined iteration converge to the *joint* KL projection rather than oscillating.

`Σ X[c] = n` is preserved by the Sinkhorn + bisection structure, so **no final renormalisation is needed** (source comment, `sinkhorn.cpp`).

## Architecture

```mermaid
flowchart TD
    A[c_api] --> B[sinkhorn solve]
    B --> C[Sinkhorn marginal sweep]
    C --> D[bisect_capacity: solve mu]
    D --> E[Dykstra: a[c] += log residual]
    E --> F{converged?}
    F -- no --> C
    F -- yes --> G[finalize_weights]
```

- **Marginal sweep + capacity bisection** alternate each iteration.
- **Dykstra state** `a[c]` is carried across iterations (this solver is *not* stateless → no SRAA-m wrapper, unlike RAKING/IEPPA/Greenkhorn).
- `convergence_rule = 0` — sinkhorn uses no improvement criterion; it iterates to the fixed point.
- **Finalize**: Σw=n already holds; `bounds_mode` dispatch still runs for per-obs unit mode.

## Advantages

- **Joint** margin+capacity KL projection via Dykstra — theoretically the cleanest bounded-KL treatment of the box.
- Mass `Σ=n` preserved exactly throughout — no normalize/bounds ordering hazard.
- Bisection capacity step has a **guaranteed valid bracket** and monotone objective → no failure-to-bracket.
- Log domain → numerically stable under extreme scaling.

## Drawbacks

- **Stateful** (Dykstra vectors) → cannot use the Anderson/SRAA acceleration the other IPF-family solvers enjoy.
- Bisection inner loop adds cost per iteration (mitigated by the `_fast` variant).
- First-order convergence; Dykstra can be slow when many cells are simultaneously bound.
- More moving parts than plain raking for the same KL objective.

## Mathematical guarantees and proofs

| Claim | Status | Basis |
|-------|--------|-------|
| Sinkhorn marginal scaling converges (positive kernel) | **Strong** | Sinkhorn–Knopp (1967); Franklin–Lorenz. |
| Dykstra converges to the projection onto the intersection of convex sets | **Strong** | Dykstra (1983); Bauschke–Lewis for Bregman/KL projections. The capacity box and marginal affine sets are convex. |
| Bisection always brackets a root (`f` monotone, valid `[lo,hi]`) | **Strong** | Proven in source by the clamp-threshold derivation (`sinkhorn.cpp` comments); `f` strictly increasing in μ. |
| `Σ X = n` invariant | **Strong** | Asserted and structurally preserved (source comment, verified by fixed-point trace, B11). |
| Convergence rate | **Weak** | No rate bound in-repo; Dykstra rate is problem-dependent. |

**Appraisal: Strong.** Both pillars (Sinkhorn scaling, Dykstra projection) are textbook-proven convex results, and the one numerical sub-step that *could* fail (bracketing) is proven correct in the source. Only the rate is unquantified.

## Code references

| Component | File | Key symbols |
|-----------|------|-------------|
| Capacity projection | `src/sinkhorn.cpp` | `bisect_capacity`, `bisect_capacity_fast`, `mu_out`, `X_proj` |
| Dykstra correction | `src/sinkhorn.cpp` | `a[c]` accumulation |
| Marginal sweep | `src/sinkhorn.cpp` | `scale[j]`, `S_kj` |
| Invariant note | `src/sinkhorn.cpp` | "sum(X)=n preserved … no normalization needed" |
| Finalize | `src/calib_dispatch.hpp` | `finalize_weights` |

---

## History & seminal sources

The iterative scaling heuristic predates its formal proof by decades. Deming & Stephan (1940) applied it as **Iterative Proportional Fitting** (IPFP) in contingency-table adjustment; economists knew it as the **RAS method** [Bacharach 1965].

**Sinkhorn (1964)** [sinkhorn1964relationship] supplied the first formal convergence proof: every strictly positive square matrix can be written as $D_1 A D_2$ where $D_1, D_2$ are diagonal and $D_1 A D_2$ is doubly stochastic. Alternating row/column normalisation converges to this form.

**Sinkhorn & Knopp (1967)** [sinkhorn1967diagonal] extended the result to rectangular matrices and gave tighter conditions.

**Dykstra (1983)** [dykstra1983] proposed his alternating-projection algorithm for finding the nearest point in the intersection of multiple closed convex sets. The original setting was Euclidean (least-squares). The correction vectors (dual variables) that distinguish Dykstra from plain alternating projections ensure convergence to the true intersection even for non-affine constraint sets.

**Cuturi (2013)** [cuturi2013sinkhorn] recast discrete optimal transport as an entropically regularised linear programme and showed that the dual iterates are exactly Sinkhorn row/column scalings on the Gibbs kernel $K_{ij} = \exp(-C_{ij}/\varepsilon)$. This made Wasserstein distances tractable for machine-learning loss functions (handwritten-digit recognition benchmark improved by orders of magnitude). The paper that established OT as a practical ML tool.

**Benamou, Carlier, Cuturi, Nenna & Peyré (2015)** [benamou2015iterative] generalised the Cuturi framework to a broad class of constrained OT problems—barycenters, multi-marginal transport, capacity-constrained transport, partial transport—by interpreting entropic regularisation as KL-projection onto convex constraint intersections. For affine constraints plain iterative Bregman projections suffice; for non-affine constraints (inequalities, capacity boxes) they required the **Bregman–Dykstra** extension of [dykstra1983] to KL divergence. Published in *SIAM Journal on Scientific Computing* 37(2):A1111–A1138, DOI 10.1137/141000439.

**Peyré & Cuturi (2019)** [peyrecuturi2019computational] is the standard textbook reference for computational OT, covering Sinkhorn complexity ($O(n^2/\varepsilon^2)$–$O(n^2/\varepsilon^3)$), log-domain stabilisation, multiscale schedules, and Wasserstein barycenters. *Foundations and Trends in Machine Learning* 11(5-6):355–607.

---

## Practitioner implementations & use cases

| Library | Language | Backend | Primary use |
|---------|----------|---------|------------|
| **POT** [flamary2021pot] | Python | NumPy / CuPy | General-purpose: standard, partial, unbalanced, semi-relaxed OT; Frank-Wolfe conditional-gradient for semi-relaxed problems |
| **OTT-JAX** | Python/JAX | XLA GPU | Differentiable OT in JAX pipelines; Sinkhorn divergences with JIT + vmap |
| **GeomLoss** (KeOps) | Python/PyTorch | KeOps CUDA | OT as differentiable loss in neural networks; memory-efficient GPU Sinkhorn divergences via online log-sum-exp |

**Dominant use cases (post-2013):** domain adaptation, Wasserstein barycenters, generative-model training (Wasserstein GANs), point-cloud registration, colour transfer, statistical inference with distributional robustness, and assignment / matching problems in combinatorial optimisation.

---

## Known caveats & concerns (literature)

### ε → 0 numerical instability

The Gibbs kernel $K_{ij} = \exp(-C_{ij}/\varepsilon)$ spans extreme ranges as $\varepsilon \to 0$. With costs bounded by 10 and $\varepsilon = 0.01$, entries reach $e^{-1000}$, far below IEEE 754 double precision floor. Catastrophic underflow zeros the kernel rows, causing division-by-zero in the scaling update [peyrecuturi2019computational].

### Entropic blur & regularisation bias

Entropy forces the transport plan to have diffuse support; the regularised plan resembles a mixture of near-isotropic Gaussians. Suboptimality vs. the unregularised plan is bounded by $\varepsilon \log(nm)$ [cuturi2013sinkhorn, peyrecuturi2019computational]. Reducing bias requires smaller $\varepsilon$, which worsens instability.

### Log-domain stabilisation cost

Working with dual potentials $f, g$ and the log-sum-exp trick eliminates underflow. Cost: every update requires per-entry transcendental (`exp`, `log`) evaluations instead of batched matrix-vector products—significantly slower on CPU; GPU implementations compensate with warp-level shuffle reductions [arXiv:2605.00837].

### Convergence rate and ε-dependence

Sinkhorn enjoys geometric (linear) convergence with Birkhoff contraction ratio $\lambda \to 1$ as $\varepsilon \to 0$. The rate degrades exponentially: $1 - \lambda \sim e^{-1/\varepsilon}$ in hard cases [peyrecuturi2019computational]. Overall complexity to approximate unregularised OT within accuracy $\delta$: $\tilde{O}(n^2/\delta^3)$ (Altschuler et al. 2017), improved to $\tilde{O}(n^2/\delta^2)$ for Greenkhorn.

### Dykstra accumulation for non-affine constraints

Standard alternating KL projections converge to the intersection only for affine (equality) constraints. Capacity *boxes* (inequality constraints) require Dykstra correction vectors; omitting them yields incorrect projections [benamou2015iterative].

---

## How leafblower deviates

Leafblower's Sinkhorn is **not** a general entropic OT solver and does not operate on a cost matrix or Gibbs kernel. It is a **survey calibration** engine:

- **No cost matrix / Gibbs kernel.** The starting point is cell-level design weights $X[c]$ (not a joint coupling matrix). There is no $C_{ij}$ or regularisation parameter $\varepsilon$; the KL geometry is implicit in the Sinkhorn marginal sweeps.
- **Capacity projection replaces the column marginal step.** Instead of normalising columns to a target marginal, leafblower projects $X$ onto a capacity box $[L[c], U[c]]$ per cell using `bisect_capacity_fast` — a $\mu$-bisection on the log-scale multiplier. This is the Dykstra-projected step from [benamou2015iterative] applied per-cell rather than per row/column of a coupling matrix.
- **Dykstra accumulator `a[c]`.** After each capacity projection, the log-residual $\log X_\text{proj}[c] - \log X[c]$ is added to `a[c]`. On the next bisection call, $a[c]$ shifts the effective scale so the projection remains consistent with all prior constraint resolutions — the KL analogue of Dykstra's correction vectors [dykstra1983, benamou2015iterative].
- **Marginal sweep = Sinkhorn row scaling.** The inner `scale[j]` loop normalises the marginal distribution of $X$ across cells within each category, exactly as Sinkhorn row/column scaling normalises a matrix to its marginals.
- **`at_lower` freeze.** Cells permanently clamped to their lower bound are frozen (`at_lower` flag) and excluded from future Dykstra corrections to prevent unbounded `a[c]` accumulation — a numerical guard not present in the standard textbook algorithm.
- **No entropic bias.** Because there is no $\varepsilon$ parameter, leafblower does not suffer entropic blur. The KL penalty is inherited from the Sinkhorn fixed-point structure, not injected as a regulariser. The target is exact marginal balance, not an approximate transport plan.
- **Output:** calibrated survey weights summing to $n$ (via `finalize_weights`), not a probability coupling matrix.

---

## References

[sinkhorn1964relationship] R. Sinkhorn. "A relationship between arbitrary positive matrices and doubly stochastic matrices." *Annals of Mathematical Statistics*, 35(2):876–879, 1964. DOI: 10.1214/aoms/1177703591.

[sinkhorn1967diagonal] R. Sinkhorn and P. Knopp. "Concerning nonnegative matrices and doubly stochastic matrices." *Pacific Journal of Mathematics*, 21(2):343–348, 1967. DOI: 10.2140/pjm.1967.21.343.

[dykstra1983] R. L. Dykstra. "An algorithm for restricted least squares regression." *The American Statistician*, 37(4):384–386, 1983. DOI: 10.1080/00031305.1983.10483138.

[cuturi2013sinkhorn] M. Cuturi. "Sinkhorn distances: Lightspeed computation of optimal transport." *Advances in Neural Information Processing Systems* 26 (NeurIPS 2013). URL: http://papers.neurips.cc/paper/4927-sinkhorn-distances-lightspeed-computation-of-optimal-transport.

[benamou2015iterative] J.-D. Benamou, G. Carlier, M. Cuturi, L. Nenna, and G. Peyré. "Iterative Bregman projections for regularized transportation problems." *SIAM Journal on Scientific Computing*, 37(2):A1111–A1138, 2015. DOI: 10.1137/141000439.

[peyrecuturi2019computational] G. Peyré and M. Cuturi. "Computational optimal transport: With applications to data sciences." *Foundations and Trends in Machine Learning*, 11(5–6):355–607, 2019. DOI: 10.1561/2200000073. arXiv: 1803.00567.

[flamary2021pot] R. Flamary, N. Courty, A. Gramfort, et al. "POT: Python Optimal Transport." *Journal of Machine Learning Research*, 22(78):1–8, 2021. URL: https://jmlr.org/papers/v22/20-1286.html.
