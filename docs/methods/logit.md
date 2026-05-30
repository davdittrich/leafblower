# LOGIT — Bounded Logit-Link Newton Calibration (Deville–Särndal 1992)

> Enum: `RK_ALG_LOGIT = 10`
> Source: `src/logit_calib.cpp`, `src/logit_calib.hpp`

## Overview

The **logit distance function** from Deville & Särndal (1992) — the standard *bounded* calibration estimator. Each cell mass is forced into `[L_c, U_c]` by a **logistic link**, so bounds are satisfied *by construction* (no clamping, no active set). The dual multipliers are found by **Newton on the normal equations with Armijo line search**. This is the "autumn::calibrate" style. Compared with GREG (χ², bounds via active set) and raking (KL, bounds via water-fill), LOGIT bakes the box into the parametrisation itself.

## Mathematical formulation

### Link and objective

```
z[c] = Σ_k λ[k, g_k(c)]                        (linear predictor from duals)
σ[c] = 1 / (1 + exp(−z[c]))                     (logistic)
w[c] = L_c + (U_c − L_c)·σ[c]                   (bounded weight — always in [L,U])
```

This is the calibration weight implied by the Deville–Särndal logit distance; minimising that distance s.t. the margins gives the Newton system below. Bounds `L_c ≤ w[c] ≤ U_c` hold identically because `σ ∈ (0,1)`.

### Newton step (verified, `logit_calib.cpp`)

```
D_eff[c] = max( kDeffFloor·range_c , range_c·σ[c]·(1−σ[c]) )      range_c = U_c − L_c
N = Aᵀ D_eff A
b[m] = T[m]·n − Σ_c w[c]                          (margin defect, scaled by n)
solve  LDLᵀ(N) δλ = b
λ ← λ + α·δλ                                       (Armijo, sufficient-decrease c=0.01)
```

- **D_eff floor** `kDeffFloor = 1e-6·range` prevents the Newton step `b/D_eff` blowing up when σ saturates (σ→0 or 1). The inverse-scaling alternative was tested and **rejected** (source comment).
- **Armijo norm guard**: step capped so no `z` coordinate shifts > ~2 (prevents divergence into saturation).
- **Warm start (Layer 2)**: solve `(AAᵀ)λ₀ = A·z_target` with `z_target = logit(σ_target)`, `σ_target = clip((X_init−L)/range, ε, 1−ε)` — places λ₀ in the convergence basin. Rejected if any component > 10 (would saturate).
- **Saturation guard**: early-exit if > 50% of cells saturate (`|z| > 650`).

## Architecture

```mermaid
flowchart TD
    A[c_api] --> B[logit solve]
    B --> C[Layer 2 warm start: λ₀ from logit-inverse]
    C --> D[compute w = L+(U-L)σ(z), D_eff]
    D --> E[normal equations N, defect b]
    E --> F[LDLT solve δλ]
    F --> G[Armijo line search + norm guard]
    G --> H{residual plateau?}
    H -- no, still improving --> D
    H -- budget/stall --> I[finalize_weights]
```

- **Convergence**: residual-shrinkage plateau detection (`best_resid < init·0.999` → budget vs stall).
- **D_eff cache**: `w` and `D_eff` recomputed post-step so the convergence check sees post-step weights.
- **stall_kind**: logit-specific classification on `RK_ERR_STALL`.

## Advantages

- **Bounds satisfied exactly by construction** — no clamping, no active set, no water-fill; the iterate is always strictly inside `(L,U)`.
- **Second-order** (Newton) with Armijo globalisation → fast and safe when feasible.
- Smooth everywhere (logistic) → well-defined Hessian, unlike the kinked active-set/water-fill methods.
- Canonical survey-calibration estimator with established inferential properties.

## Drawbacks

- **Saturation pathology**: as cells approach a bound, `σ(1−σ)→0`, `D_eff→0`, and Newton steps explode — needs the floor + norm guard + 50%-saturation early-exit (substantial guard machinery).
- On infeasible / near-infeasible tight-bound configs it can hit the iteration ceiling regardless of tuning (source profiling note).
- Requires forming/factorising normal equations each iteration — O(K²·M)-ish.
- Convergence basin sensitivity → needs the Layer-2 warm start to behave.

## Mathematical guarantees and proofs

| Claim | Status | Basis |
|-------|--------|-------|
| Logit distance yields bounded calibration weights with the stated form | **Strong** | Deville & Särndal (1992), *JASA*; the logit `G`-function is one of their canonical distances. |
| Weights stay in `[L,U]` for all iterates | **Strong** | Algebraic: `σ∈(0,1) ⇒ w∈(L,U)` identically (verified in update rule). |
| Newton + Armijo converges (sufficient decrease) | **Moderate→Strong** | Armijo line-search Newton converges to a stationary point of a smooth convex calibration objective; standard, but the objective convexity in the bounded regime is assumed not re-proved in-repo. |
| Convergence on infeasible tight-bound problems | **Weak** | Explicitly not guaranteed — hits iteration ceiling (source comment); feasibility is a precondition. |
| Step-explosion control at saturation | **Moderate** | Floor + norm guard are tuned engineering safeguards (profiled across seeds), not a proof. |

**Appraisal: Strong for the estimator and bound-feasibility; Moderate for convergence.** The Deville–Särndal logit method is a rigorously-established calibration estimator and the box-feasibility is algebraically guaranteed. The practical risk is the saturation regime, mitigated by guards rather than theorems.

## Code references

| Component | File | Key symbols |
|-----------|------|-------------|
| Link + Newton weight | `src/logit_calib.cpp` | `σ`, `w = L+(U−L)σ`, `D_eff` |
| D_eff floor | `src/logit_calib.cpp` | `kDeffFloor`, rejected inverse-scaling note |
| Warm start | `src/logit_calib.cpp` | Layer-2 `(AAᵀ)λ₀ = A z_target` |
| Armijo guard | `src/logit_calib.cpp` | `kArmijoC`, `kMaxHalvings`, norm guard |
| Saturation guard | `src/logit_calib.cpp` | `|z| > 650`, 50% early-exit |

---

## History & seminal sources

Deville and Särndal (1992) [devillesarndal1992] introduced a unified framework of calibration estimators parametrised by a distance function G(w, d). The logit distance was presented as "case 6" in their taxonomy: it implies g-weights g_k = w_k/d_k constrained to a user-specified interval [L, U] via an inverse-logit (logistic) link on the dual parameter. The explicit purpose was to eliminate two pathologies of earlier methods: the unbounded linear method (GREG) can produce negative weights or catastrophically large positive ones, and the exponential/raking method prevents negatives but has no upper barrier, leaving estimates exposed to outlier-driven inflation.

Deville, Särndal, and Sautory (1993) [devillesarndalsautory1993] followed immediately with a companion paper focusing on marginal-total calibration (raking on cross-classified tables) and variance estimation, and described the CALMAR SAS macro [sautory1993calmar]. The two papers together — universally cited as "Deville and Särndal 1992" and "Deville et al. 1993" — established the theoretical canon for modern calibration estimation. The 1992 paper has accumulated ~2,000 citations.

Deville was head of the Unit of Statistical Methodology at INSEE (Institut National de la Statistique et des Études Économiques), where the French school of official statistics had been developing reweighting ideas since at least the 1970s (Thionet, Deville 1988, Deville and Särndal 1990). The 1992 paper consolidated and formalised this tradition; Sautory's CALMAR macro was the direct operational implementation.

A 25-year retrospective by Devaud and Tillé (2019) [devaudtille2019] provides a rigorous formalisation of the calibration problem as a convex optimisation problem, proving existence conditions for bounded calibration solutions and comparing the distance functions analytically and in simulation.

---

## Practitioner implementations & use cases

**CALMAR / CALMAR 2 (SAS macro, INSEE).**  
The canonical reference implementation. Developed by Olivier Sautory at INSEE in 1993 [sautory1993calmar], CALMAR (CALage sur MARges) implements linear, raking, and logit distance functions via Newton–Raphson on the dual system. CALMAR 2 extended it with multi-step and domain calibration and improved handling of rank-deficient design matrices via Moore–Penrose pseudoinverse. INSEE used it operationally for virtually all French official household surveys.

**`survey::calibrate(calfun = "logit", bounds = c(L, U))` — R, Lumley [lumley2004survey].**  
The most widely used open-source implementation. `calfun = "logit"` is one of three built-in distance functions (`"linear"`, `"raking"`, `"logit"`); the `bounds` argument is mandatory for logit. The Newton–Raphson solver is exposed as the internal `grake` function; the documentation attributes the algorithm to Deville et al. (1993). Convergence caveats — "calibration with bounds, or on highly collinear data, may fail" — are documented explicitly; `force = TRUE` returns the partially-calibrated design for diagnosis.

**`sampling::calib(method = "logit", bounds = c(low, upp))` — R, Tillé and Matei [tille2016sampling].**  
The `sampling` package implements the same four distance functions as `calib()` (linear, raking, truncated, logit). Bounds default to `c(low = 0, upp = 10)` for the logit and truncated methods. References cite Deville and Särndal (1992) and Deville et al. (1993) directly.

**ReGenesees (`e.calibrate(calfun = "logit")`) — R, Istat / Zardetto.**  
R Evolved Generalized Software for Sampling Estimates and Errors in Surveys, developed by Diego Zardetto at Istat (Italy). Supports linear, raking, and logit distance functions with bounded and partitioned (domain-level) calibration. The `partition` argument decomposes a global calibration problem into independent sub-problems, making logit feasible for very large surveys. Used operationally for the Italian Labour Force Survey (*Indagine Forze di Lavoro*).

**`%SurveyCalibrate` (SAS/STAT macro, SAS Institute).**  
Commercial implementation. Notable for automated bounds search: if L/U are unspecified, the macro minimises U and maximises L iteratively until feasibility is achieved, avoiding the need for manual tuning.

**Bascula (Statistics Netherlands), CLAN (Statistics Sweden), gCALIB (Statistics Belgium).**  
NSI-internal tools all implementing logit as one available calibration function. Adoption is particularly strong in Europe for population surveys requiring non-negative, range-restricted weights.

**Operational use cases.**  
Logit calibration (or bounded variants) is used wherever a practitioner needs a single unified weight vector that is simultaneously (a) always positive, (b) capped to prevent extreme leverage, and (c) exactly matches demographic or administrative margins. Canonical applications:

- **Household and social surveys with nonresponse correction** — SHARE, ESS, GGP and similar panels apply logit calibration to combat unit nonresponse and panel attrition while keeping weights in practitioner-acceptable ranges.
- **Labour force surveys** — Statistics Canada LFS, Istat IFL; bounded weights prevent single respondents from representing implausibly large population fractions.
- **Establishment and business surveys** — BLS and Statistics Canada applications; extreme design weights are endemic in business frames (heavy-tailed size distributions), making upper bounds essential.
- **Population synthesis / synthetic population generation** — spatial microsimulation uses logit-type calibration to fit synthetic populations to census constraints; range restriction prevents degenerate household replication counts.

---

## Known caveats & concerns (literature)

**1. Convergence failure near bounds (saturation).**  
The Newton step scales with the inverse of the logistic derivative σ(1−σ). As a cell weight approaches L or U, σ → 0 or 1 and σ(1−σ) → 0, making the effective Hessian singular and the step size diverge. Deville and Särndal (1992) acknowledge this as a structural limitation: the logit barrier function is asymptotically repulsive but gradient-based methods approaching it lose second-order information. The literature (Devaud & Tillé 2019 [devaudtille2019]; Espuny-Pujol et al. 2018 [espunypujol2018]) treats this as a known numerical pathology requiring safeguards.

**2. Infeasibility with tight bounds.**  
If the interval [L, U] is too narrow relative to the spread required by the margin constraints, no λ satisfying all constraints exists. Devaud and Tillé (2019) provide a rigorous discussion of existence conditions. The practical symptom is algorithm failure before convergence; the `survey::calibrate` documentation states this explicitly. The logit distance acts as a barrier only within the feasible region — if the feasible region is empty, the barrier has no effect and the algorithm diverges.

**3. Choice of L and U.**  
No principled closed-form criterion exists; practitioners use heuristics. The most common approach (documented in Sautory's CALMAR manual and confirmed by the statistical computing literature): (a) first run an unbounded linear calibration to observe the natural g-weight spread; (b) use that spread as a lower bound on what [L, U] must accommodate; (c) then tighten iteratively from a wide starting interval. Automated bound-search (SAS %SurveyCalibrate) minimises U and maximises L algorithmically. The choice is ultimately a policy decision balancing representativeness against weight extreme-value tolerance; no universally optimal selection rule is known.

**4. Comparison with alternatives.**  
When logit is infeasible:
- *Truncated linear*: uses a modified chi-square distance with infinite penalty outside [L, U]. Shares the infeasibility pathology with identical algebraic structure near bounds.
- *Weight trimming (clipping)*: sets a looser feasible interval, solves exactly, then clips extremes and redistributes residual mass. Guaranteed to terminate but breaks exact margin constraints — the calibrated estimates no longer reproduce population totals.
- *Ridge/relaxed calibration*: penalises deviation from margins rather than enforcing hard constraints (Rao et al. 2002; Lesage et al. 2024 [relaxedcalib2024]). Guarantees a bounded, feasible solution at the cost of approximate margin matching; the bias introduced is controlled by the ridge parameter.
- *Global optimisation / always-convergent approaches*: Espuny-Pujol et al. (2018) [espunypujol2018] propose an ℓ₁-minimisation feasibility check followed by an alternate optimisation that degrades gracefully when [L, U] cannot be fully satisfied, avoiding hard failure.

---

## How leafblower deviates

Standard CALMAR-style logit calibration (and its `survey`/`sampling` implementations) operates on observation-level g-weights g_k = w_k/d_k with globally shared [L, U]. Leafblower's implementation diverges in several ways:

**Cell-level bounds.**  
Leafblower applies per-cell bounds `[L_c, U_c]` derived from the `CalibState` (set by the caller per cell). This is strictly more general than a global [L, U]: it subsumes the standard case when all cells share the same bounds but enables heterogeneous per-cell box constraints that CALMAR and `survey::calibrate` do not natively expose (those tools apply bounds multiplicatively to the individual design weight d_k, which is equivalent only if all d_k are equal).

**Logistic link on the cell weight directly, not the g-weight.**  
Standard CALMAR applies the logit link to the *ratio* g_k = w_k/d_k (so the distance function is G(g_k, 1) in Deville-Särndal notation, applied to the ratio relative to the design weight). Leafblower applies the logistic link to `w_c = L_c + (U_c − L_c)·σ(z_c)` — the absolute weight, not the ratio. The two parametrisations are equivalent when d_k = 1 or when the design weights are absorbed into L_c/U_c by the caller, but differ numerically when per-cell bounds encode heterogeneous design-weight information. This is a deliberate choice: it makes the bound constraint exact in the weight scale, not in the g-weight scale.

**D_eff floor and saturation guards (better near-bound stability).**  
Standard implementations (survey's `grake.R`, CALMAR Newton loop) do not document explicit handling of the σ → 0/1 saturation pathology beyond iteration budgets. Leafblower adds:
- `kDeffFloor = 1e-6 · range_c` prevents the Hessian from going singular by flooring the logistic derivative;
- `kMaxDeltaZ` caps the dual step to prevent runaway z outside the linear regime;
- 50%-saturation early-exit aborts and reports `RK_ERR_NOCONV` rather than silently diverging.

These guards make leafblower more robust on near-infeasible tight-bound configurations at the cost of potentially missing solutions that a bare Newton method (with sufficient iterations and step damping) could find.

**Armijo line search (vs. Newton with fixed step).**  
Standard CALMAR uses a Newton method with unit step. `survey::calibrate` via `grake.R` uses a similar update. Leafblower uses Armijo sufficient-decrease backtracking: the step α is halved until the residual norm decreases by at least `kArmijoC = 0.01` of the expected decrease. This globalises convergence (guarantees monotone residual decrease for any starting point with a sufficient-decrease step) at the cost of more residual evaluations per iteration, and is more conservative near bounds.

**When these choices help vs. hurt:**
- *Help*: near-infeasible configurations, tight per-cell bounds, heterogeneous design weights. The saturation guards and Armijo line search prevent silent divergence.
- *Hurt*: well-conditioned problems far from bounds where standard Newton converges in 5–10 steps — leafblower's guards add overhead (floor, step clipping, backtracking evaluations). In that regime, the standard `grake` approach is faster and numerically equivalent.

---

## References

[devillesarndal1992] Deville, J.-C. and Särndal, C.-E. (1992). Calibration estimators in survey sampling. *Journal of the American Statistical Association*, **87**(418):376–382. DOI: 10.1080/01621459.1992.10475217

[devillesarndalsautory1993] Deville, J.-C., Särndal, C.-E., and Sautory, O. (1993). Generalized raking procedures in survey sampling. *Journal of the American Statistical Association*, **88**(423):1013–1020. DOI: 10.1080/01621459.1993.10476369

[sautory1993calmar] Sautory, O. (1993). La macro CALMAR: redressement d'un échantillon par calage sur marges. Document de travail, Institut National de la Statistique et des Études Économiques (INSEE), Paris. (Internal technical report; described in Deville et al. 1993.)

[lumley2004survey] Lumley, T. (2004). Analysis of complex survey samples. *Journal of Statistical Software*, **9**(8):1–19. DOI: 10.18637/jss.v009.i08

[tille2016sampling] Tillé, Y. and Matei, A. (2021). *sampling: Survey Sampling*. R package version 2.10. URL: https://cran.r-project.org/package=sampling. Methods based on Tillé, Y. (2020). *Sampling and Estimation from Finite Populations*. Wiley.

[devaudtille2019] Devaud, D. and Tillé, Y. (2019). Deville and Särndal's calibration: revisiting a 25-years-old successful optimisation problem. *TEST*, **28**(4):1033–1065. DOI: 10.1007/s11749-019-00681-3

[espunypujol2018] Espuny-Pujol, F., Morrissey, K., and Williamson, P. (2018). A global optimisation approach to range-restricted survey calibration. *Statistics and Computing*, **28**(2):427–439. DOI: 10.1007/s11222-017-9739-5

[relaxedcalib2024] Lesage, É. et al. (2024). Relaxed calibration of survey weights. *Survey Methodology*, **50**(2). URL: https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X202400200012
