# NEWTON-KL — Truncated-SVD Trust-Region Newton on the Dual KL Objective

> Enum: `RK_ALG_NEWTON_KL = 11`
> Source: `src/newton_calib.cpp`, `src/newton_calib.hpp`

## Overview

A **second-order (Newton) solver for the same maximum-entropy / KL calibration problem that IEPPA solves by IPF** — but attacked through its **smooth dual**. It minimises the convex dual gap `g(λ) = log Z(λ) − T·λ` directly, using an eigendecomposition of the Hessian with **truncated-SVD (TSVD) regularisation**, **Levenberg–Marquardt damping**, and a **Steihaug–CG trust region**. Designed for the "zero-compression regime" where weights are non-saturating and the dual is well-behaved.

## Mathematical formulation

### Dual objective (verified, `newton_calib.cpp`)

The KL calibration dual (log-linear / exponential family):

```
g(λ) = log Z(λ) − T·λ ,    Z(λ) = Σ_i d_i · exp(u_i),   u_i = (Aᵀλ)_i
weights:  w_i = d_i·exp(u_i) / Z · n
```

with **log-sum-exp stabilisation** `log Z = u_max + log Σ d_i exp(u_i − u_max)` to avoid overflow / false convergence. Gradient `∇g = (margin sums) − T`; Hessian `H = ` weighted covariance of the auxiliaries (PSD).

### Newton step (verified)

```
H = V Λ Vᵀ                                   (eigendecomposition, dsyevd)
retain R = { i : Λ_i ≥ ratio·λ_max }          (TSVD truncation, ratio default 1e-8)
Λ_damp[i] = Λ_i·(1+μ_LM) + μ_LM·floor          (Levenberg–Marquardt damping)
g_keep = Vᵀ∇g  (eigenbasis)
if ‖g_keep/Λ_damp‖ ≤ Δ:  δ_keep = g_keep/Λ_damp        (pseudoinverse fast path)
else:                    δ_keep = Steihaug-CG(g_keep, Λ_damp, Δ)   (trust boundary)
λ ← λ − α·V δ_keep                             (Armijo line search)
```

- **Adaptive μ_LM** via Marquardt gain ratio ρ: ρ>0.75 → μ/=3, ρ<0.25 → μ×10.
- **Adaptive trust radius** Δ (no upper cap, plan rev-2 fix H).
- Fallback to full-Hessian LDLᵀ if EVD fails or no eigenvalues retained.

## Architecture

```mermaid
flowchart TD
    A[c_api] --> B[newton_kl solve]
    B --> C[eval_dual: g, Z, u_max  (log-sum-exp)]
    C --> D[eigendecompose H = V Λ Vᵀ]
    D --> E[TSVD retain Λ_i ≥ ratio·λ_max]
    E --> F[LM damping Λ_damp]
    F --> G{‖step‖ ≤ Δ ?}
    G -- yes --> H[pseudoinverse step]
    G -- no --> I[Steihaug-CG to trust boundary]
    H & I --> J[Armijo line search on λ]
    J --> K{dual gap < tol?}
    K -- no --> C
    K -- yes --> L[weights = exp(u)/Z·n -> finalize]
```

- **Best-iterate** tracked in dual space (`lam_best`, `best_gap`) separately from the weight-space `BestIterTracker`.
- **TSVD ratio** configurable (`st.newton_tsvd_ratio`, default 1e-8).
- Convergence: `‖∇g‖_∞ < tol` (dual-gap optimality).

## Advantages

- **Quadratic local convergence** — far fewer iterations than IPF near the optimum when the dual is well-conditioned.
- **TSVD + LM + trust region** make it robust to rank-deficient / near-singular Hessians (collinear margins) where naïve Newton blows up.
- Log-sum-exp stabilisation prevents the classic overflow false-convergence.
- Operates on the *smooth convex dual* — globally convergent with the trust region.

## Drawbacks

- **Eigendecomposition every iteration** (`dsyevd`) — O(K³) in the number of free duals; expensive for many categories.
- Confined to the **"zero-compression regime"**; saturating/bounded weights are not its target (those go to logit/sinkhorn/raking).
- Many hyperparameters (TSVD ratio, LM init, trust radius policy) — more tuning surface.
- Hard bounds not part of the dual objective here; box enforcement is left to finalize.

## Mathematical guarantees and proofs

| Claim | Status | Basis |
|-------|--------|-------|
| Dual `g(λ)` is convex; its minimiser gives the primal KL solution | **Strong** | Convex duality of max-entropy / exponential-family calibration (standard). |
| Trust-region Newton (Steihaug-CG) is globally convergent to a stationary point | **Strong** | Steihaug (1983); Conn–Gould–Toint, *Trust-Region Methods*. |
| Local quadratic convergence near the optimum | **Strong** | Standard Newton theory under a positive-definite Hessian at the solution. |
| TSVD truncation preserves descent on rank-deficient H | **Moderate** | Pseudoinverse step is a descent direction on the retained subspace; the dropped directions are (near-)null — standard but rate impact not quantified in-repo. |
| Robustness across all inputs | **Moderate** | LDLᵀ fallback + LM damping are engineering guards; convergence proven for the well-conditioned regime it targets. |

**Appraisal: Strong (within its regime).** Convex duality + trust-region Newton are rigorously proven; the TSVD/LM machinery is standard regularisation. The honest caveat is scope — it is designed and proven for the non-saturating regime, not as a universal solver.

## Code references

| Component | File | Key symbols |
|-----------|------|-------------|
| Dual eval | `src/newton_calib.cpp` | `eval_dual`, `g(λ)=log Z − T·λ`, log-sum-exp |
| Eigen + TSVD | `src/newton_calib.cpp` | `dsyevd`, retain `Λ_i ≥ ratio·λ_max` |
| Trust region | `src/newton_calib.cpp` | Steihaug-CG, `Δ`, LM `μ` gain-ratio |
| Dual best-iterate | `src/newton_calib.cpp` | `lam_best`, `best_gap` |
| Config | `src/newton_calib.cpp` | `st.newton_tsvd_ratio` |

---

## History & seminal sources

The problem leafblower's Newton-KL solver addresses has three intellectual lineages:

**Maximum-entropy / KL calibration distance.**  
Ireland and Kullback [irelandkullback1968] provided the information-theoretic foundation: minimising the directed Kullback-Leibler (backward KL) divergence from uniform weights subject to fixed marginal constraints is equivalent to the Iterative Proportional Fitting (raking) procedure introduced by Deming and Stephan (1940). They proved the resulting estimator is consistent, asymptotically normal, and efficient.

**Survey calibration as unified convex optimisation.**  
Deville and Särndal [devillesarndal1992] cast calibration as a single constrained convex optimisation over weight-distance functions, proved that the exponential (raking-ratio) method corresponds exactly to KL divergence, and showed that optimal weights are obtained by solving a *dual* system of nonlinear equations for Lagrange multipliers — the form leafblower attacks directly with Newton steps on `g(λ)`. They also proved all such calibration estimators are asymptotically equivalent to the GREG estimator.  
Deville, Särndal and Sautory [devillesarndalsautory1993] extended the framework to generalised raking on marginal totals and variance estimation; they proved that exponential calibration on two categorical margins recovers classical raking exactly.

**Newton-type dual solvers for calibration.**  
Algorithm 1 of Deville and Särndal [devillesarndal1992] already proposed Newton-Raphson on the dual, with each step being the update `λ ← λ + (XᵀTX)⁻¹(t − t̂)`. This is the direct ancestor of plain `survey::calibrate` (Newton-Raphson path) and all subsequent practitioner implementations. The dual is strictly concave when the Hessian `H = XᵀWX` is positive definite, which makes Newton the natural solver and motivates the regularity conditions explored in [devaudtille2019].

**Trust-region machinery.**  
The theoretical guarantee that the Steihaug truncated-CG subproblem solver finds a good approximate trust-region step — stopping at the boundary when negative curvature is encountered or the CG iterate crosses the trust radius — is due to Steihaug [steihaug1983] (SIAM J Numer Anal, 20:626–637, 1983). Toint (1981) independently proposed a similar approach, leading to the combined name Steihaug-Toint. The definitive algorithmic treatment, including gain-ratio update rules, convergence theory, and Levenberg-Marquardt interpretation, is in Conn, Gould and Toint [conngouldtoint2000].

---

## Practitioner implementations & use cases

All major survey calibration packages implement Newton-Raphson on the dual. **None implements a full trust-region Newton on the dual.** The dominant reason: when the calibration distance function is strictly convex the dual is strictly concave and has a positive-definite Hessian under linear independence of the auxiliary variables; trust-region boundary projection and negative-curvature handling are structurally unnecessary. Plain Newton-Raphson with step-halving (backtracking) when domain constraints are hit is sufficient and universal in practice.

| Implementation | Notes |
|---|---|
| R `survey::calibrate` [lumley2004survey] | Newton-Raphson (all bounded/nonlinear methods); Deville–Särndal (1993) algorithmic basis; drops collinear columns automatically. Reference implementation for the field. |
| R `ReGenesees::e.calibrate` (Zardetto, ~2015) | Newton-Raphson with `maxit=50` default; `partition` argument factorises the problem by sub-population domain to manage memory at Istat scale. |
| R `GECal::GEcalib` (Kwon, Kim & Qiu, 2024) | Generalised entropy calibration; delegates to `nleqslv` (Newton or Broyden); handles ML-derived auxiliary predictors and debiasing constraints. |
| SAS `%SurveyCalibrate` (An, 2020) | Newton-Raphson in SAS/IML; exponential/logit/truncated-linear methods; auto-generates replicate weights. |
| CALMAR/CALMAR2 (Sautory, INSEE) | Original production macro; iterative Newton on margins; widely used in French official statistics and copied internationally. |
| Python `GECal` / `pycalibrate` | Wrappers calling `scipy.optimize` Newton/quasi-Newton; no trust-region subproblem. |

**Use cases:** post-stratification weight adjustment, nonresponse correction, frame undercoverage adjustment, domain estimation, calibration to census totals, entropy balancing for causal inference (Hainmueller 2012).

---

## Known caveats & concerns (literature)

The following failure modes are documented in the survey calibration literature [devillesarndal1992; devaudtille2019; lumley2004survey]:

1. **Rank deficiency / collinear margins.** The Newton Hessian `H = XᵀWX` is singular when the sample design matrix is rank-deficient (empty cells, perfectly collinear auxiliary variables, over-determined margin system). Standard implementations either fail or drop columns; `survey::calibrate` explicitly warns "Calibration with bounds, or on highly collinear data, may fail." [lumley2004survey] Increasing the number of calibration constraints strictly reduces the set of feasible solutions; parsimonious variable selection is therefore required.

2. **Infeasibility under tight bounds.** When box constraints (lower/upper weight bounds) are added, the feasible region shrinks. If bounds are too tight relative to the sample and control totals, no solution exists and Newton diverges [devillesarndal1992; devaudtille2019].

3. **Domain restriction violation ("improper" calibration functions).** Raking/exponential calibration is globally defined (exponential never goes negative), but empirical-likelihood and hyperbolic distance functions have domain restrictions. A full Newton step can cross a singularity; step-halving / backtracking is mandatory [devillesarndalsautory1993].

4. **TSVD truncation bias.** Hard spectral thresholding (TSVD) discards entire singular components, introducing bias in directions associated with small-but-nonzero singular values. In ill-conditioned regimes this trades ill-conditioning for systematic underfitting: the retained dual subspace does not capture the full margin discrepancy. Tikhonov regularisation is a smooth alternative but attenuates *all* directions; both induce regularisation bias [conngouldtoint2000].

5. **Non-convergence in near-saturating / high-dimensional regimes.** When calibration constraints saturate the degrees of freedom (many fine-grained interactions, near-empty cells), weights become extreme, the Hessian eigenspectrum becomes ill-conditioned, and convergence stalls. Newton can diverge or cycle in this regime; IPF / raking or iEPPA are better suited [devillesarndal1992].

6. **Quadratic convergence only near the solution.** Newton guarantees quadratic local convergence under a positive-definite Hessian at the optimum; globally the Hessian can be near-singular during early iterations, making the pure Newton step unreliable without damping.

---

## How leafblower deviates

Standard practitioner implementations (e.g., `survey::calibrate`) apply plain Newton-Raphson:

```
λ ← λ + (XᵀWX)⁻¹ (t − t̂)
```

with a simple direct solve and optional step-halving. Leafblower's Newton-KL solver deviates in four compounding ways:

| Deviation | Leafblower | Practitioner standard | Better when | Worse when |
|---|---|---|---|---|
| **Hessian eigensolver (TSVD)** | Full `dsyevd` eigen-decomposition; retain modes `Λᵢ ≥ ratio·λ_max` | Cholesky / LU direct solve | Hessian is near-singular; avoids numerical blow-up | Well-conditioned H: O(K³) eigen cost unnecessary vs O(K³/3) Cholesky |
| **Trust-region step control** | Steihaug-CG subproblem with adaptive radius Δ; gain-ratio update [steihaug1983; conngouldtoint2000] | None (full Newton step or step-halving) | Early iterations far from optimum; prevents over-shooting | Near the solution: trust region is inactive, adds overhead with no benefit |
| **Levenberg-Marquardt damping** | Adaptive μ on H diagonal; μ relaxed when gain ρ > 0.75, tightened when ρ < 0.25 | No LM; step-halving only | Hessian ill-conditioned; μ stabilises the inversion | Adds a hyper-parameter (μ₀, μ_min, μ_max) that can slow convergence if μ fails to relax |
| **Dual best-iterate tracking** | `lam_best` / `best_gap` retained separately from current iterate | Final iterate only | Oscillatory convergence or early stopping | No overhead concern; purely beneficial |

**Net assessment:** The TSVD + LM + trust-region combination makes leafblower more robust than plain Newton-Raphson when the Hessian is ill-conditioned (rank-deficient margin design, collinear controls, early iterations from a poor warm-start). The cost is O(K³) eigendecomposition every iteration and three extra hyper-parameters. For well-conditioned problems where plain Newton converges in 5–10 iterations, `survey::calibrate` is cheaper. Leafblower's approach is justified for the regime it targets: moderate K, potentially ill-conditioned margins (e.g., stepstone with K=9 interacting categorical margins), and the requirement of a principled, regularisation-aware step rather than a heuristic step-halve.

The approach has no published precedent in the survey calibration literature (notebooklm confirmed no prior trust-region Newton dual implementation exists in the field). It applies a standard optimisation engineering technique (trust-region with Steihaug-CG) to a domain where the well-conditioned case does not need it but the ill-conditioned case can benefit significantly. This is an engineering deviation, not a methodological one.

---

## References

[devillesarndal1992] Deville, J.-C. and Särndal, C.-E. (1992). Calibration estimators in survey sampling. *Journal of the American Statistical Association*, 87(418), 376–382. doi: 10.1080/01621459.1992.10475217

[devillesarndalsautory1993] Deville, J.-C., Särndal, C.-E. and Sautory, O. (1993). Generalized raking procedures in survey sampling. *Journal of the American Statistical Association*, 88(423), 1013–1020.

[steihaug1983] Steihaug, T. (1983). The conjugate gradient method and trust regions in large scale optimization. *SIAM Journal on Numerical Analysis*, 20(3), 626–637. doi: 10.1137/0720042

[conngouldtoint2000] Conn, A. R., Gould, N. I. M. and Toint, P. L. (2000). *Trust-Region Methods*. MPS-SIAM Series on Optimization. Society for Industrial and Applied Mathematics, Philadelphia. doi: 10.1137/1.9780898719857

[lumley2004survey] Lumley, T. (2010). *Complex Surveys: A Guide to Analysis Using R*. John Wiley & Sons, Hoboken, NJ. (R package: `survey`, CRAN, doi: 10.32614/CRAN.package.survey)

[irelandkullback1968] Ireland, C. T. and Kullback, S. (1968). Contingency tables with given marginals. *Biometrika*, 55(1), 179–188.

[devaudtille2019] Devaud, D. and Tillé, Y. (2019). Deville and Särndal's calibration: revisiting a 25-years-old successful optimization problem. *TEST*, 28(4), 1033–1065. doi: 10.1007/s11749-019-00681-3
