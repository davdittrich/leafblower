# GREG — Generalized Regression Estimator (Active-Set Newton, χ² / linear method)

> Enum: `RK_ALG_GREG = 6`
> Source: `src/greg.cpp`, `src/greg.hpp`

## Overview

The **linear (chi-square distance) calibration estimator of Deville & Särndal (1992)** — the classical GREG / generalized regression estimator — solved by **active-set Newton** with hard box constraints. The "distance" minimised is quadratic (χ²) rather than KL, which gives the well-known linear calibration formula `w = d·(1 + xᵀλ)`. Bounds are handled by an active set: cells that would leave `[L_c, U_c]` are fixed at the bound and the normal equations are refactorised only on active-set changes.

## Mathematical formulation

### Objective

```
min_X   Σ_c ( X[c] − X_init[c] )² / X_init[c]            (chi-square / GREG distance)
s.t.    Σ_{c∈(k,j)} X[c] = T_kj    ∀ k,j
        L_c ≤ X[c] ≤ U_c           ∀ c
```

### Newton / KKT solve (verified, `greg.cpp`)

For free cells `D_eff[c] = X_init[c]`. Build and factor the normal equations, solve for the per-category duals `λ`, then the linear calibration update:

```
N = Aᵀ D_eff A           (normal matrix; A = margin incidence)
b = Aᵀ D_eff (T − S)     (RHS)
solve  LDLᵀ(N) λ = b
X_new[c] = X_init[c] · ( 1 + Σ_k λ[k, g_k(c)] )           # Deville–Särndal linear form
```

Active-set handling:

```
if X_new[c] < L_c:  X[c]=L_c, fixed_lo[c]=true, mark_refactor
if X_new[c] > U_c:  X[c]=U_c, fixed_hi[c]=true, mark_refactor
KKT release: if a fixed cell's multiplier has the wrong sign → unfix, refactor
```

Optional **Tikhonov ridge** `st.ridge_lambda` adds `τI` to `N` before Cholesky for ill-conditioned designs.

## Architecture

```mermaid
flowchart TD
    A[c_api] --> B[greg solve]
    B --> C[D_eff = X_init for free cells]
    C --> D[compute_normal_equations N,b]
    D --> E[optional Tikhonov ridge]
    E --> F[LDLT factor + solve lambda]
    F --> G[X = X_init*(1+sum lambda)]
    G --> H[clamp to box -> set fixed_lo/hi]
    H --> I[KKT release pass]
    I --> J{active set changed?}
    J -- yes --> D
    J -- no --> K[finalize_weights]
```

- **Refactor cache (R1)**: `compute_normal_equations` written directly into the factored buffer; refactor only when the active set changes.
- **Objective**: `select_solver_objective → m.chi2`; `best_error = chi2`.
- **Convergence**: no active-set change in an iteration ⇒ KKT optimum.

## Advantages

- **Closed-form-ish linear update** (`w = d(1+xᵀλ)`) — the canonical, interpretable calibration estimator with a long survey-statistics pedigree.
- **Second-order**: Newton on a quadratic ⇒ exact in one step *absent bounds*; active-set adds only a few refactors.
- **Refactor caching** makes the constrained solve cheap once the active set settles.
- χ² distance yields weights that are an explicit linear function of the auxiliaries — easy variance estimation.

## Drawbacks

- χ² weights **can go negative / explode** without bounds — the box + active set are essential, and tight boxes cause many refactors.
- Quadratic distance is less robust to extreme adjustments than KL (KL keeps weights positive automatically; χ² does not).
- Forms/factorises the normal equations — O(K²·M) per refactor; expensive for many categories.
- Ridge `τ` (if used) biases the estimator — a bias/conditioning trade-off.

## Mathematical guarantees and proofs

| Claim | Status | Basis |
|-------|--------|-------|
| Unconstrained χ² calibration has the closed-form GREG solution | **Strong** | Deville & Särndal (1992), *JASA* — the foundational calibration paper. |
| Active-set Newton converges to the constrained KKT optimum | **Strong** | The problem is a convex QP; active-set methods converge finitely for strictly convex QPs (Nocedal–Wright). |
| Termination (finite active-set changes) | **Moderate** | True for non-degenerate QPs; cycling possible under degeneracy (not specifically anti-cycling-guarded in source). |
| Positivity of weights | **Weak** | Not guaranteed by χ² itself — only the explicit `[L,U]` box enforces it. |

**Appraisal: Strong.** GREG is the textbook calibration estimator with a rigorous closed-form and convex-QP convergence theory. The only caveats are degenerate-cycling (mild) and the intrinsic non-positivity of χ² distance (handled by the box).

## Code references

| Component | File | Key symbols |
|-----------|------|-------------|
| Newton/active-set loop | `src/greg.cpp` | `D_eff`, `fixed_lo`, `fixed_hi`, `need_refactor` |
| Normal equations | `src/greg.cpp` / `src/calib_linalg.cpp` | `compute_normal_equations`, `N_factored` |
| Linear update | `src/greg.cpp` | `X = X_init*(1 + Σ lambda)` |
| Ridge | `src/greg.cpp` | `st.ridge_lambda` |
| Objective | `src/calib_dispatch.hpp` | `select_solver_objective` → `m.chi2` |

---

## History & seminal sources

The calibration estimator framework was established by Deville & Särndal (1992)
[devillesarndal1992], who showed that a large class of survey estimators can be
derived by minimising a distance function between design weights and calibrated
weights subject to linear calibration constraints. The chi-square (linear) distance
`G(r) = (r−1)²/2` yields the closed-form GREG update `w_k = d_k(1 + xₖᵀλ)`,
where `λ` solves the normal equations. This is the same linear weighting formula
that had appeared as the GREG estimator in Cassel, Särndal & Wretman (1976) and
was systematised in the textbook [sarndalswenssonwretman1992], but Deville &
Särndal's 1992 paper embedded it inside a unified distance-minimisation
framework that admits raking, logit, and other calibration functions as special
cases.

The theoretical context is Deville & Särndal's *Newton-Raphson* approach: for a
general distance function the calibration equations are non-linear in `λ` and
require iterative Newton-Raphson; for the linear (chi-square) distance the system
collapses to a single linear solve, making the unconstrained problem exact in one
step. The survey-statistics community built on this foundation through the 1990s
and 2000s; Särndal (2007) [sarndal2007calibration] provides the definitive
retrospective review, discussing advances in variance estimation, nonresponse
adjustment, multi-phase designs, and methods for controlling extreme weights.

**Prior art.** The calibration idea extends Deming-Stephan (1940) iterative
proportional fitting and the post-stratification estimator. The generalized
regression estimator (GREG) name predates 1992 (Fuller & Isaki 1981; Bethlehem
& Keller 1987); Deville & Särndal's contribution was the unifying distance
framework and the Newton-Raphson proof of equivalence.

**First software implementation.** INSEE's CALMAR SAS macro (Sautory, 1993) was
the first dedicated implementation. Statistics Canada's GES, Statistics Sweden's
CLAN, CBS Netherlands' BASCULA, and ISTAT's GENESEES (predecessor to ReGenesees)
followed in the mid-1990s as national statistical institutes (NSIs) operationalised
calibration for large-scale household and business surveys.

---

## Practitioner implementations & use cases

The following table summarises current open-source implementations
[lumley2004survey, tille2016sampling, zardetto2015regenesees]:

| Package | Maintainer | Algorithm (linear/chi² case) | Key limitation |
|---------|-----------|-------------------------------|---------------|
| `survey::calibrate` | T. Lumley | Unbounded linear: direct linear solve; bounded: Newton-Raphson | Convergence failures with tight bounds or collinear aux; trimming via `trimWeights` violates calibration equations |
| `sampling::calib` | Y. Tillé & A. Matei | Newton-Raphson with Moore-Penrose generalised inverse inside solver | Can fail with very many constraints or overly strict box constraints |
| `ReGenesees::e.calibrate` | D. Zardetto (ISTAT) | Same as `survey`; native partitioned calibration for sub-populations | Taylor-linearised variance only; no custom distance functions |

**Primary use cases in official statistics:**

1. **Census-consistency:** forcing survey estimates to match administrative
   register totals (age × sex × region cells) — the original motivation in
   [devillesarndal1992].
2. **Nonresponse & undercoverage adjustment:** calibrating respondent weights
   to external demographic benchmarks to reduce nonresponse bias
   [sarndal2007calibration].
3. **Business surveys:** calibrating establishment weights to frame totals
   (employment size, industry) — common at Eurostat, ONS, INE.
4. **Multi-purpose weighting:** a single weight vector satisfying several margin
   constraints simultaneously, reusable across many survey variables
   [sarndalswenssonwretman1992].

The chi-square (linear) distance is the default in CALMAR, GES, and
`sampling::calib`. The logit distance is preferred when strictly positive weights
are required without an explicit box constraint.

---

## Known caveats & concerns (literature)

The following problems are documented in the survey-statistics literature and
apply to the linear/chi-square calibration method specifically.

### 1. Negative and extreme weights

The chi-square distance function is defined over all reals; without bounds,
calibrated weights can be negative or very large [devillesarndal1992,
sarndal2007calibration]. Särndal (2007) notes explicitly: *"A few of the weights
computed according to (3.2) can turn out to be quite large or negative."*
Causes include:

- Many simultaneous benchmark constraints (over-determined margins).
- High collinearity among auxiliary variables.
- Leverage outliers in the auxiliary data.

### 2. Variance inflation

Negative weights must be offset by large positive weights to preserve calibration
totals, inflating standard errors — especially for small domains and survey
variables weakly correlated with auxiliaries. Extreme g-weight variation directly
inflates linearised variance estimates. Treating calibrated weights as ordinary
sampling weights in variance calculations (without accounting for the calibration
step) generally overestimates standard errors [sarndal2007calibration].

### 3. Feasibility under box constraints

Imposing a box `[L, U]` on chi-square calibration is not always feasible: no
solution may exist that simultaneously satisfies the calibration equations and the
bounds. Särndal (2007) warns that *"imposing strict intervals on the weights means
a mathematical solution is no longer guaranteed."* This is the fundamental tension
between range restriction and exact calibration.

### 4. Variance estimation with bounded calibration

If calibration fails in a bootstrap or jackknife replicate (a calibration equation
cannot be solved), retaining these replicates biases variance estimates upward.
`survey::calibrate` warns but returns the object; `ReGenesees` is stricter.

### 5. Ridge calibration as remedy

Adding a Tikhonov penalty `τ‖λ‖²` to the normal equations (ridge calibration,
also called penalised calibration) relaxes the exact-match constraint, trades
some calibration accuracy for bounded weights, and regularises ill-conditioned
designs. The ridge GREG weight is a closed-form variant of the standard GREG
formula. Tuning `τ` is a bias-variance trade-off with no universally optimal
choice.

---

## How leafblower deviates

Standard implementations (`survey`, `sampling`, `ReGenesees`) solve the
unconstrained or bounded linear calibration problem with one of two approaches:

- **Unbounded linear:** direct single-step linear solve (one matrix inversion).
- **Bounded (any distance):** Newton-Raphson iterating on the dual `λ` until
  the calibration residual is below tolerance.

**leafblower uses active-set Newton instead of Newton-Raphson on the dual.**
Concretely:

| Property | Standard Newton-Raphson (survey/sampling) | leafblower active-set |
|---------|------------------------------------------|----------------------|
| Solve target | Dual variables `λ` via iterative NR on `∇ℓ(λ) = 0` | Primal cell weights `X[c]` directly |
| Bound handling | Projection or clipping after each NR step; may oscillate | Hard active set: clamped cells removed from factorisation |
| Factorisation reuse | Full refactor every NR iteration | Refactor only on active-set change (R1 cache) |
| Convergence criterion | Calibration-residual norm below tol | No active-set change in one full sweep = KKT optimum |
| Degeneracy / cycling | No anti-cycling guard in NR; can oscillate near bounds | No anti-cycling guard; cycling theoretically possible for degenerate QPs |
| Ridge | Not standard; implemented ad hoc in some packages | `st.ridge_lambda` adds `τI` to `N` before Cholesky |

**Where active-set is better:** the refactor cache makes repeated solves cheap
once the active set stabilises — typical in practice when few cells are at bounds.
The KKT termination criterion gives an exact convex-QP certificate rather than a
residual-norm approximation.

**Where active-set is worse:** active-set Newton works at the cell (aggregate)
level, not the observation level. The number of cells `M_cell` is the problem
dimension, not `n`. For designs with very many fine cells and many cells at
bounds simultaneously, the repeated LDLᵀ factorisations of the normal matrix
(O(K²·M) per refactor) are expensive. Standard NR on the dual is O(K³) per
iteration independent of `M`, which is cheaper when `M ≫ K`.

The ridge term (`st.ridge_lambda`) introduces estimator bias; the caller is
responsible for choosing `τ` appropriately. This parallels the known bias/variance
trade-off in ridge calibration discussed in the literature.

---

## References

[devillesarndal1992] Deville, J.-C. and Särndal, C.-E. (1992). Calibration
estimators in survey sampling. *Journal of the American Statistical Association*,
**87**(418), 376–382.
doi:[10.1080/01621459.1992.10475217](https://doi.org/10.1080/01621459.1992.10475217)

[sarndalswenssonwretman1992] Särndal, C.-E., Swensson, B. and Wretman, J. (1992).
*Model Assisted Survey Sampling*. Springer, New York.
isbn:9780387406206
doi:[10.1007/978-1-4612-4378-6](https://doi.org/10.1007/978-1-4612-4378-6)

[sarndal2007calibration] Särndal, C.-E. (2007). The calibration approach in
survey theory and practice. *Survey Methodology*, **33**(2), 99–119.
Statistics Canada, Catalogue no. 12-001-X.
url:[https://www150.statcan.gc.ca/n1/pub/12-001-x/2007002/article/10488-eng.pdf](https://www150.statcan.gc.ca/n1/pub/12-001-x/2007002/article/10488-eng.pdf)

[lumley2004survey] Lumley, T. (2004). Analysis of complex survey samples.
*Journal of Statistical Software*, **9**(8), 1–19.
doi:[10.18637/jss.v009.i08](https://doi.org/10.18637/jss.v009.i08)

[tille2016sampling] Tillé, Y. and Matei, A. (2024). *sampling: Survey Sampling*.
R package. CRAN.
doi:[10.32614/CRAN.package.sampling](https://doi.org/10.32614/CRAN.package.sampling)
url:[https://CRAN.R-project.org/package=sampling](https://CRAN.R-project.org/package=sampling)

[zardetto2015regenesees] Zardetto, D. (2015). ReGenesees: an advanced R system
for calibration, estimation and sampling error assessment in complex sample
surveys. *Journal of Official Statistics*, **31**(2), 177–203.
doi:[10.1515/jos-2015-0013](https://doi.org/10.1515/jos-2015-0013)
