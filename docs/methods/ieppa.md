# IEPPA — Iterative Entropy-Penalized Proportional Adjustment

> Enum: `RK_ALG_IEPPA = 1` (+ variant `RK_ALG_IEPPA_SOFT = 8`)
> Source: `src/ieppa.cpp`, `src/ieppa_internal.hpp`, `src/ieppa_finalize.cpp`, `src/ieppa_trajectory.cpp`

## Overview

IEPPA is **iterative proportional fitting (IPF / RAS / Sinkhorn–Knopp) on log-Sinkhorn factors, with successive over-relaxation (SOR) stepping and an infeasibility damping factor**. It is the project's default solver. It adjusts cell masses `X[c]` so every marginal sum matches its target while staying as close as possible (in Kullback–Leibler divergence) to the design weights `X_init`.

The solver keeps one multiplicative factor per `(margin k, category j)` — the Sinkhorn dual scaling — and sweeps the margins repeatedly. It runs in **linear space** when the compression ratio is small and switches to **log space** to avoid overflow when factors grow extreme (cutoff `kLinearSpaceThreshold = 2.0`).

## Mathematical formulation

### Objective

Minimise the I-divergence (generalized KL) of the calibrated masses from the design masses subject to all marginal constraints:

```
min_X   Σ_c X[c] · log( X[c] / X_init[c] ) − X[c] + X_init[c]
s.t.    Σ_{c ∈ (k,j)} X[c] = T_kj    for every margin k, category j
```

This is the standard raking / maximum-entropy calibration problem. Its dual has one Lagrange multiplier per category; the multiplicative factor is `f_kj = exp(λ_kj)`.

### Update rule (verified in `ieppa_solve`)

Per category, the *naïve* full Sinkhorn ratio is the target over the current marginal sum:

```
naive = T_kj · W_input / S_kj          where S_kj = Σ_{c∈(k,j)} X_cur[c] / f_old
```

The applied step is a power step in linear space, equivalently a fractional step in log space:

```
linear:  f_new = f_old^(1−α·ω) · naive^(α·ω)
log:     lf_new = (1−α·ω)·lf_old + α·ω·(log T_kj·W − log S_kj)
```

- `ω` (`eff_omega`) is the SOR relaxation parameter, auto-adapted per margin after a burn-in; `ω = 1` recovers plain IPF (fast path, no `pow()`).
- `α` is the **infeasibility-streak damping factor** `α = 1/(1 + β·stress)` (`compute_alpha`), where `stress` = longest consecutive infeasible-bucket streak and `β = kAlphaBeta = 0.5`; no stress ⇒ `α = 1` (fast path). It shrinks the step when a margin cannot be satisfied (keeps `α ∈ (0,1]`, preserving Peyré–Cuturi §4.4 convergence).
- `β` is **only** the constant inside that damping map (the η-schedule scales it per homotopy level, `β = 0.5·η`). It is **not** a proximal/entropic term — the core update carries no proximal term.
- The net exponent on the Sinkhorn ratio is `α·ω`; margins are swept **Gauss–Seidel (BCD-style)**.
- `ω = 1, α = 1` ⇒ exact Sinkhorn–Knopp.

### IEPPA_SOFT variant (enum 8)

Adds an **augmented-Lagrangian / ADMM** soft-capacity term: per-cell capacity bounds are not hard-clamped each sweep but enforced through a penalty `μ` driven up across outer iterations, with the KL Newton step `X̃(1−λ+μz)/(1+ρ)` for the un-normalized-KL generator. (Do **not** "correct" this formula — it is right for this generator; see `CLAUDE.md`.)

## Architecture

```mermaid
flowchart TD
    A[harvest.R / c_api] --> B[ieppa_solve]
    B --> C{compression ratio < 2.0?}
    C -- yes --> D[linear-space sweep]
    C -- no --> E[log-space sweep]
    D & E --> F[SOR omega auto-adapt]
    F --> G{SRAA accel on?}
    G -- yes --> H[Anderson acceleration sraa.hpp]
    G -- no --> I[flat fixed-point loop]
    H & I --> J[ieppa_finalize: obs expand + bounds]
    J --> K[finalize_weights: Σw=n then bounds_mode]
```

- **Core loop**: `ieppa_solve` in `src/ieppa.cpp` (kept in one TU — it is the hot path; cold code lives in `ieppa_finalize.cpp` / `ieppa_trajectory.cpp`, see no-LTO note in `CLAUDE.md`).
- **Scheduler**: round-robin or greedy (largest per-margin residual first).
- **Acceleration**: optional SRAA-m (Anderson) via `sraa.hpp`, opt-in through `st.accelerate`.
- **Homotopy**: optional multi-level cascade over `max_weight` multipliers (`rk_homotopy_cfg_t`); the outer driver iterates these levels.
- **Finalize**: `finalize_weights` (`calib_dispatch.hpp`) enforces Σw = n, then `bounds_mode` dispatch (cell = count only; unit = per-cell water-fill).

## IEPPA vs. Sinkhorn (+Dykstra)

Both solve the **same KL/IPF fixed point on the margins** — at `α = ω = 1` the IEPPA marginal step *is* a Sinkhorn step. They diverge on:

| Axis | IEPPA (this method) | Sinkhorn (`sinkhorn.md`) |
|------|---------------------|--------------------------|
| **Box constraint `[L_c,U_c]`** | Deferred to `finalize_weights` (or, in IEPPA_SOFT, an **ADMM/ALM** penalty `μ` driven up across levels). The *core* sweep is box-free. | Enforced **jointly each iteration** by **Dykstra** correction vectors `a[c]` + a `μ`-bisection that projects onto the capacity box in KL geometry. |
| **State carried between iterations** | Only the Sinkhorn log-factors `lf` (a pure fixed-point map) ⇒ **SRAA-m Anderson acceleration applies**. | The Dykstra correction `a[c]` is iterate-history ⇒ **stateful, not SRAA-able**. |
| **Step control** | SOR over-relaxation `ω` + infeasibility damping `α = 1/(1+β·stress)` (net exponent `α·ω`); homotopy on bounds. **No proximal term.** | Plain Sinkhorn step + Dykstra; no over-relaxation. |
| **Mass `Σw`** | Re-scaled to `n` at finalize (normalize→bounds order). | `Σ = n` preserved *by construction* every iteration (bisection targets `n`). |
| **Stopping** | metric-improvement / plateau. | fixed point, no improvement rule (`convergence_rule = 0`). |

**One-line summary**: IEPPA = *over-relaxed Sinkhorn with bounds pushed to finalize and Anderson acceleration available*; Sinkhorn = *plain Sinkhorn with the bounds folded into every iteration via Dykstra*. IEPPA trades exact per-iterate box-feasibility for cheaper, accelerable sweeps; Sinkhorn trades acceleration for a joint, always-feasible KL projection.

## Relationship to the source paper (arXiv:2011.14312)

> Chu, Liang, Toh & Yang, *"An efficient implementable inexact entropic proximal point algorithm for a class of linear programming problems."* — `docs/iEPPA/arxiv.2011.14312/`. Design contract: `docs/superpowers/specs/2026-04-23-ieppa-faithful-design.md`.

The codebase implements the paper's **inner algBCD specialized to `C = 0`**. The paper's headline *outer* inexact-proximal-point loop is **not** part of this solver — at `C = 0` it is mathematically inert — so the operative guarantee is Csiszár I-projection, not the paper's outer-PPA theorem.

### What the paper's iEPPA is

The paper solves a **linear program** — capacity-constrained multi-marginal optimal transport (CMOT):

```
min ⟨C, X⟩   s.t.   A⁽ⁱ⁾(X) = b⁽ⁱ⁾  (i=1..N marginals),   0 ≤ X ≤ U
```

with an **outer inexact Entropic Proximal Point Algorithm**: at each outer step `k` it *approximately* solves a subproblem with a **re-centered** Bregman (Boltzmann–Shannon) proximal term

```
X^{k+1} ≈ argmin_{X∈Ω°}  ⟨C, X⟩ + ε_k · D_φ(X, X^k)        (re-anchored at X^k)
```

Three load-bearing pieces:
1. **Moving anchor**: `D_φ(X, X^k)` is re-centered on `X^k` each outer iteration, so `ε_k` stays **fixed and moderate** (experiments use `ε = 0.05`) while iterates converge to the **exact** LP optimum — bypassing the `ε → 0` numerical instability that plagues Sinkhorn / DyKL (Dykstra-with-KL).
2. **Inner dual Block Coordinate Descent (BCD)**: multiplicative update `u⁽ⁱ⁾_j ← b⁽ⁱ⁾_j / (A⁽ⁱ⁾(X̄))_j · u⁽ⁱ⁾_j` then capacity clip `X = min(X̄, U)` — a Sinkhorn/IPF sweep with a `min`-clip. (Paper notes DyKL ≡ this dual BCD, per Tibshirani 2017.)
3. A **checkable inexact stopping condition** with a **global-convergence theorem** (iterates reach an exact optimum under summable `ν_k, η_k, μ_k`).

### What the codebase implements

`C = 0` is the correct specialization: survey calibration has **no transport cost**, only closeness-to-prior, so the paper's `min⟨C,X⟩ + ε·D_φ` becomes (spec §2.1)

```
min_X  KL(X ‖ X_init)   s.t.  A(X) = T,  0 ≤ X ≤ U
```

- **No outer proximal-point loop.** At `C = 0` the `ε_k·KL(X‖X^k)` term scales the objective but does **not move the argmin**, so it is mathematically inert (spec §2.1, §9, §11). The calibration solution is the single I-projection of `X_init` onto {margins ∩ box}.
- **Inner dual algBCD**, over the K margin blocks + a capacity block; each block is a closed-form multiplicative Sinkhorn step, the capacity block being the paper's clip `Γ = min{U/M, 1}` ⇒ in code `X[c] = clamp(X̃[c]·W[c], L_c, U_c)`. The core update is a pure Sinkhorn ratio `f[k][j] = τ_kj·W_total / S_kj` — no `ε`/`β` proximal term.
- **Engineering modifiers** on top of the core: SOR `ω`, infeasibility damping `α = 1/(1+β·stress)`, homotopy over `max_weight`, SRAA/Anderson acceleration, and ALM capacity (SOFT). All are step-size or scheduling devices; none introduces a proximal bias.

Mapping paper → code:

| Paper component | In `ieppa_solve` | Evidence |
|-----------------|------------------|----------|
| `C = 0` reduction to KL-to-prior calibration | ✅ | spec §2.1 |
| Inner dual BCD (margin blocks + capacity block, Sinkhorn-like) | ✅ | `ieppa.cpp` algBCD; capacity clip `Γ = min{U/M,1}` |
| Entropic proximal term | ❌ none — `f = τ·W/S` is a pure Sinkhorn ratio | spec §2.3 |
| SOR `ω` + damping `α = 1/(1+β·stress)` | ✅ step-size relaxations (`β` is the damping constant, not a proximal term) | `ieppa.cpp` `compute_alpha`, `sor_omega`, `kAlphaBeta = 0.5` |
| Outer moving-anchor PPA loop `KL(X‖X^k)` | ❌ inert at `C = 0`; the outer driver is homotopy over `max_weight` | spec §2.1/§9/§11; `ieppa.cpp` homotopy driver |
| Capacity via ALM/ADMM | ✅ IEPPA_SOFT only | `use_admm_capacity`, `ALMConfig.capacity_mu` |
| Convergence guarantee | Csiszár (1975) / Csiszár–Tusnády (1984) cyclic I-projection, linear under Slater — not the paper's outer-PPA theorem | spec §9 |

### Appraisal

- **Operative guarantee is Csiszár / Csiszár–Tusnády cyclic I-projection** onto affine sets ∩ log-convex box (linear rate under Slater, spec §9) — the paper's outer-PPA theorem does not apply because that loop is inert at `C = 0`.
- **The core has no entropic-proximal bias.** The update `new_f = f_old^(1−α·ω)·naive^(α·ω)` uses two step-size relaxations (`ω`, `α`); neither moves the fixed point — at a fixed point `naive = 1` regardless of `α`, `ω`, so margins are satisfied exactly. With `ω = 1, α = 1` the solver is plain Sinkhorn/IPF (Sinkhorn–Knopp convergence applies).

*(Confidence: 90 — facts stated in the design spec §2.1/§2.3/§9/§11, source headers, and verified update lines in `ieppa.cpp`.)*

## Advantages

- **Simple, cheap iterations**: each sweep is O(non-zero cell incidences); no linear solve.
- **Robust default**: converges on severely skewed dual landscapes where interior-point/Newton methods stall (see comment in `c_api.cpp`).
- **Log-space fallback** prevents overflow/underflow on extreme compression.
- **SOR + Anderson** give super-linear-ish empirical speed-up while keeping the cheap-iteration structure.

## Drawbacks

- **Linear (first-order) convergence** in the unaccelerated regime — slow when margins are highly correlated (ill-conditioned).
- SOR relaxation `ω` is heuristic and auto-tuned; a bad `ω` can oscillate (mitigated by burn-in + adaptation).
- Hard bounds are handled at finalize, not inside the core loop (the SOFT variant addresses this but adds penalty-tuning complexity).
- The infeasibility damping `α` is an engineering safeguard without a clean optimality interpretation.

## Mathematical guarantees and proofs

| Claim | Status | Basis |
|-------|--------|-------|
| Converges to the unique KL-projection when a feasible interior point exists and ω=1 | **Strong** | Classical Sinkhorn–Knopp / IPF convergence (Csiszár 1975; Sinkhorn–Knopp 1967); the plain step is exactly IPF |
| Bounded-KL (margins ∩ box) convergence | **Strong** | Csiszár (1975) / Csiszár–Tusnády (1984) cyclic I-projection, linear under Slater (spec §9) |
| Geometric (linear) convergence rate | **Moderate** | Holds for IPF under positivity; rate depends on the contraction modulus of the margin coupling. Not re-derived in-repo |
| SOR / damping preserve the fixed point | **Strong** | `α·ω`-step has the same fixed point as IPF (`f_new = f_old` ⇔ `naive = 1` ⇔ marginal satisfied), independent of `α`, `ω` |
| Convergence *with* adaptive SOR | **Weak** | No formal proof; SOR can diverge for ω outside (0,2) in general; mitigated by burn-in/adaptation |
| SRAA/Anderson acceleration convergence | **Weak/Moderate** | Anderson acceleration has local convergence theory under contraction; here wrapped with a stall-revert guard rather than proven |

**Appraisal (GRADE-style): Strong for the base method, Moderate overall.** The core rests on well-proven Sinkhorn/IPF and cyclic I-projection theory. The acceleration/damping layers are step-size relaxations that preserve the fixed point but are empirically guarded, not formally proven, for rate.

## Code references

| Component | File | Key symbols |
|-----------|------|-------------|
| Core solve | `src/ieppa.cpp` | `ieppa_solve`, `eff_omega`, `compute_alpha`, `naive`, `f_lin`, `S_lin` |
| Log/linear switch | `src/ieppa.cpp` | `kLinearSpaceThreshold` |
| Finalize / bounds / obs-expand | `src/ieppa_finalize.cpp` | unit/cell branch |
| Trajectory diagnostics | `src/ieppa_trajectory.cpp` | — |
| Acceleration | `src/sraa.hpp` | SRAA-m Anderson |
| Shared finalize | `src/calib_dispatch.hpp` | `finalize_weights`, `finalize_weights_buf` |
| Source paper + design contract | `docs/iEPPA/arxiv.2011.14312/`, `docs/superpowers/specs/2026-04-23-ieppa-faithful-design.md` | Chu–Liang–Toh–Yang (arXiv:2011.14312); reference MATLAB in `docs/iEPPA/code/` |

## History & seminal sources

The algorithm at the core of this solver — plain Sinkhorn–Knopp / IPF — has been independently discovered across at least three disciplines.

**Earliest recorded use.** J. Kruithof (1937) used what he called the "double factor method" to predict telephone traffic between switching stations, scaling an origin–destination matrix row-by-row then column-by-column [kruithof1937]. G.V. Sheleikhovskii applied a functionally identical procedure to traffic-network analysis around the same period.

**Statistical canonisation.** W.E. Deming and F.F. Stephan (1940) formalized the algorithm to adjust sampled cross-tabulations from the US Census of Population so that internal cells matched updated marginal totals from the full census while preserving the sample correlation structure [demingstephan1940]. Ironically their derivation sought to minimize a Pearson chi-square; the algorithm they constructed actually converges to a maximum-entropy (minimum KL) solution, not chi-square. This paper is conventionally treated as the procedure's birth in the statistical literature.

**Matrix-theoretic analysis.** Richard Sinkhorn (1964) proved that every positive square matrix is diagonally equivalent to a doubly stochastic matrix; Sinkhorn and Knopp (1967) extended this to non-negative matrices and doubly stochastic targets [sinkhorn1967diagonal]. In applied mathematics the procedure is therefore called the *Sinkhorn–Knopp algorithm* or *diagonal scaling*.

**I-divergence geometry.** L.M. Bregman (1967) gave the first proof of necessary and sufficient conditions for convergence of IPF on matrices with zeros using an L1 approach. I. Csiszár (1975) independently established the same conditions via information geometry, showing that IPF computes successive I-projections (KL-projections) onto affine marginal constraint sets [csiszar1975idivergence]. Csiszár and Tusnády (1984) then proved convergence of alternating minimization (the general class containing IPF) from an information-geometry standpoint and derived the linear convergence rate [csiszartusnady1984].

**Optimal-transport context.** Cuturi (2013) recast Sinkhorn iteration as the inner loop for entropy-regularized optimal transport, igniting the modern machine-learning interest in the algorithm and enabling GPU-parallel computation [cuturi2013sinkhorn]. The Chu–Liang–Toh–Yang iEPPA paper (2022) [chu2022ieppa] sits in this lineage, wrapping Sinkhorn/algBCD in an outer proximal-point driver to solve the full capacity-constrained OT LP without the instability of taking the regularisation parameter to zero.

## Practitioner implementations & use cases

IPF/raking is the dominant post-survey adjustment method in official statistics and academic survey research. Implementations across ecosystems:

| Package | Language | Function | Algorithm | Bounds support | Citation |
|---------|----------|----------|-----------|----------------|---------|
| `survey` | R | `rake()` / `calibrate()` | Raking (multiplicative IPF) + GREG + logit | Yes — `bounds=` arg; logit mandatory | [lumley2010survey] |
| `ipfp` | R | `ipfp()` | Classical IPF (C implementation) | No | [blocker2022ipfp] |
| `mipfp` | R | `Ipfp()` | IPFP + ML + chi-square + LSQ on N-dimensional arrays | No explicit bounds | [barthelemy2018mipfp] |
| `icarus` | R | `calibration()` | GREG + raking + logit (Calmar-inspired) | Yes — simplex / bisection | [rebecq2017icarus] |
| `ReGenesees` | R | `e.calibrate()` | GREG + raking + logit via Newton–Raphson | Yes — mandatory for logit | [zardetto2015regenesees] |
| `ipfn` | Python | `ipfn()` | Classical multi-dimensional IPFP | No | — |
| `ipfraking` | Stata | `ipfraking` | Raking (outer-inner cycle) | Yes — hard trimming | [kolenikov2014ipfraking] |
| CALMAR/CALMAR2 | SAS macro | `%CALMAR` | GREG + raking + logit (Deville–Särndal) | Yes — `LO=` / `UP=` args | [sautory1993calmar] |

**Use cases.** Three domains account for the vast majority of deployed IPF:

1. **Survey calibration / raking** — adjusting panel or telephone-poll weights to match known demographic margins (age, sex, education, region). Deville and Särndal (1992) [devillesarndal1992] unified the family of calibration estimators, showing that raking corresponds to the multiplicative distance function minimised subject to marginal constraints. This is the exact problem our solver targets.

2. **Population synthesis for agent-based models** — constructing a synthetic individual-level population from aggregate census cross-tabulations so that simulated agents reproduce observed marginal distributions. IPF is the workhorse here; Barthélemy and Toint (2013) provide a practical review [barthelemy2013popsynth].

3. **Input–output table updating (RAS)** — national accounts statisticians use RAS (a synonym for two-dimensional IPF) to update supply/use tables when only new row and column totals are available [holy2019ras].

## Known caveats & concerns (literature)

The following failures are documented across the sources above and in the practitioner literature.

**Zero-cell / structural-zero failures.** IPF preserves the exact zero pattern of the seed matrix by induction: entry $(i,j)$ is zero after any number of iterations iff it was zero initially [csiszar1975idivergence]. If a target margin is strictly positive but the corresponding sample cell is empty (a *sampling zero*), the multiplicative update divides by zero and the algorithm fails. The necessary and sufficient condition for convergence (Bregman 1967; Csiszár 1975) is the existence of a joint distribution that satisfies all margins *and* has the same zero pattern as the seed. If no such distribution exists, the IPF sequence does not converge to a single matrix.

**Oscillation between two accumulation points.** When margins are structurally incompatible (the intersection of the two constraint sets is empty), IPF does not diverge to infinity — it oscillates. Pukelsheim (2012–2014) [pukelsheim2014biproportional] proved that the even-step and odd-step subsequences each converge to their own distinct accumulation points, giving the algorithm a predictable two-cycle behavior. This is the worst case; if only sampling zeros cause incompatibility, convergence is still lost but the oscillation is less structured.

**Bounds and convergence interference.** Hard weight trimming (clipping individual weights to $[L, U]$) restricts the feasible region. If the trimmed feasible set does not intersect the marginal-constraint set, the bounded algorithm cannot converge to a calibrated solution — it will produce a solution that satisfies bounds but misses margins, or iterate without converging [battaglia2009practical]. This is the central trade-off practitioners face: tighter bounds reduce extreme-weight variance inflation but increase the risk of divergence or residual margin error.

**Margin inconsistency across sources.** If control totals are drawn from different surveys, their sums may not agree (e.g., gender margins sum to N₁ while age margins sum to N₂). IPF has no mechanism to detect or resolve such inconsistency; it will cycle forever or fail silently. Battaglia et al. (2009) [battaglia2009practical] flag this as a frequent real-world failure mode.

**Category sparseness.** Categories comprising fewer than ~1–2% of the sample create near-zero denominators in the Sinkhorn ratio, slowing convergence and amplifying rounding errors. Combining sparse categories before fitting is the standard recommendation.

**Variance estimation after raking.** Raking converts fixed design weights into random variables, invalidating naive standard-error formulas. Taylor linearization requires access to the unadjusted design weights and second-order selection probabilities, which are rarely available in public-use files. Replication methods (jackknife, bootstrap) are the practical fallback, but if calibration fails to converge within a replicate subsample the replicate weight is undefined, biasing the variance estimator upward [battaglia2009practical].

**Margin ordering effects.** Because IPF sweeps margins cyclically, the order in which margins are updated affects transient behavior. For decomposable graphical models, Darroch, Lauritzen and Speed showed that an optimal ordering exists under which IPF converges in one or two full sweeps; in general there is no universal ordering rule. Empirically, margins with many fine-grained categories placed last tend to give worse model fit when they have low correlation with other controls.

## How leafblower deviates

The table below maps standard IPF idioms to leafblower's design choices. All claims are grounded in `src/ieppa.cpp` and the design spec.

| Axis | Standard IPF / common packages | leafblower (this solver) | Better for / Worse for |
|------|-------------------------------|--------------------------|------------------------|
| **Iteration unit** | Observation-level: scale each row/column of the raw weight matrix | **Cell-level**: compress all observations sharing the same covariate pattern into a single mass `X[c]`; only `M_cell ≪ n` cells iterated | Better: O(M_cell·K) vs O(n·K) per sweep — 10–100× faster at high compression. Worse: all obs in a cell receive identical relative adjustment (design weights preserved within cell, not freed) |
| **Step size** | Plain multiplicative ratio (`naive = T / S`) — SOR not used in `survey::rake`, `ipfp`, `mipfp` | SOR over-relaxation `ω` auto-adapted per margin + infeasibility damping `α = 1/(1+β·stress)` | Better: faster convergence in practice (super-linear empirically on well-conditioned problems). Worse: `ω > 1` can oscillate if burn-in estimate is wrong; harder to reason about formally |
| **Bounds handling** | Either no bounds (`ipfp`, `mipfp`), or bounds folded into every iteration via bisection or logit transform (`survey`, `icarus`, `ReGenesees`) | Bounds deferred to `finalize_weights` after convergence (IEPPA_SOFT uses ALM/ADMM to enforce bounds inside the loop) | Better (base): cheaper iterations, SRAA-accelerable; bound-violation is only transient. Worse: core loop is not box-feasible — final weights may differ from the KL minimizer subject to bounds if the box-free projection lands outside `[L,U]` |
| **Acceleration** | None in most packages; SQUAREM in some R wrappers | SRAA-m Anderson acceleration (optional, `st.accelerate`) | Better: empirically reduces iteration counts by 30–80% on well-conditioned cases. Worse: no formal convergence proof for the accelerated sequence beyond the local contraction regime; stall-revert guard required |
| **Zero-cell semantics** | Division-by-zero halts or NaN-propagates | `X_init[c] = 0` cells are skipped; if the corresponding target is positive, `RK_ERR_INFEAS` is latched after the sweep | Better: explicit infeasibility reporting with (k,j) enumeration instead of silent NaN. Same fundamental limitation: a zero seed cell cannot be filled |
| **Convergence criterion** | Max absolute margin error < tol | Metric-improvement + plateau detection (`kErrCheckInterval`, best-iterate tracking via SRAA) | Better: detects stalls that look like slow convergence; NOCONV signalled promptly. Worse: slightly more complex stopping logic; can terminate early on a metric plateau that would eventually resolve |
| **Variance estimation** | Replication or Taylor linearization (design-weight aware) | Not provided — leafblower is a weight-production library, not a variance estimator | Out of scope by design; downstream tools must account for calibration |

**Summary.** leafblower's cell-compression + SOR + deferred-bounds design trades exact per-iterate box-feasibility and formal convergence proof for the accelerated path for significantly cheaper iteration cost and empirical speed. The trade is favourable when `M_cell ≪ n` (high compression) and when `max_weight` bounds are loose enough that the box-free fixed point is already inside `[L,U]` — the common case in production survey calibration. With tight bounds (`max_weight` ≈ 1.1–1.3) IEPPA_SOFT's ALM layer is the recommended variant.

## References

- [demingstephan1940] W.E. Deming and F.F. Stephan, "On a Least Squares Adjustment of a Sampled Frequency Table When the Expected Marginal Totals are Known," *Annals of Mathematical Statistics*, 11(4):427–444, 1940. doi:10.1214/aoms/1177731829
- [sinkhorn1967diagonal] R. Sinkhorn and P. Knopp, "Concerning nonnegative matrices and doubly stochastic matrices," *Pacific Journal of Mathematics*, 21(2):343–348, 1967.
- [csiszar1975idivergence] I. Csiszár, "I-Divergence Geometry of Probability Distributions and Minimization Problems," *Annals of Probability*, 3(1):146–158, 1975. doi:10.1214/aop/1176996454
- [csiszartusnady1984] I. Csiszár and G. Tusnády, "Information Geometry and Alternating Minimization Procedures," *Statistics & Decisions*, Supplement 1:205–237, 1984.
- [chu2022ieppa] H.T.M. Chu, L. Liang, K.-C. Toh, and L. Yang, "An efficient implementable inexact entropic proximal point algorithm for a class of linear programming problems," arXiv:2011.14312v3, 2022. doi:10.48550/arXiv.2011.14312
- [cuturi2013sinkhorn] M. Cuturi, "Sinkhorn Distances: Lightspeed Computation of Optimal Transport Distances," *Advances in Neural Information Processing Systems* 26, 2013.
- [devillesarndal1992] J.-C. Deville and C.-E. Särndal, "Calibration Estimators in Survey Sampling," *Journal of the American Statistical Association*, 87(418):376–382, 1992. doi:10.1080/01621459.1992.10475217
- [lumley2010survey] T. Lumley, *Complex Surveys: A Guide to Analysis Using R*, Wiley, 2010.
- [blocker2022ipfp] A.W. Blocker, *ipfp: Fast Implementation of the Iterative Proportional Fitting Procedure in C*, CRAN, 2022. doi:10.32614/CRAN.package.ipfp
- [barthelemy2018mipfp] J. Barthélemy and T. Suesse, "mipfp: An R Package for Multidimensional Array Fitting and Simulating Multivariate Bernoulli Distributions," *Journal of Statistical Software, Code Snippets*, 86(2):1–20, 2018. doi:10.18637/jss.v086.c02
- [rebecq2017icarus] A. Rebecq, *icarus: Calibrates and Reweights Units in Samples*, CRAN, 2017.
- [zardetto2015regenesees] D. Zardetto, "ReGenesees: an Advanced R System for Calibration, Estimation and Sampling Error Assessment in Complex Sample Surveys," *Journal of Official Statistics*, 31(2):177–203, 2015.
- [kolenikov2014ipfraking] S. Kolenikov, "Calibrating survey data using iterative proportional fitting (raking)," *Stata Journal*, 14(1):22–59, 2014.
- [sautory1993calmar] O. Sautory, "La macro CALMAR: redressement d'un échantillon par calage sur marges," *Documents de Travail*, INSEE, F9310, 1993.
- [barthelemy2013popsynth] J. Barthélemy and P.L. Toint, "Synthetic Population Generation Without a Sample," *Transportation Science*, 47(2):266–279, 2013.
- [holy2019ras] V. Holý and K. Šafr, "Disaggregating Input-Output Tables by the Multidimensional RAS Method," arXiv:1704.07814, 2019.
- [battaglia2009practical] M.P. Battaglia, D.C. Hoaglin, and M.R. Frankel, "Practical Considerations in Raking Survey Data," *Survey Practice*, 2(5), 2009. doi:10.29115/SP-2009-0019
- [pukelsheim2014biproportional] F. Pukelsheim, "Biproportional Scaling of Matrices and the Iterative Proportional Fitting Procedure," *Annals of Operations Research*, 215(1):269–283, 2014.
- [kruithof1937] J. Kruithof, "Telefoonverkeersrekening," *De Ingenieur*, 52:E15–E25, 1937.
