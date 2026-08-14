# Context (DOC-derived)

Running notes from the 11 `DOC`-classified documents. Topic-keyed, attributed, and kept
close to the source wording. These are explanatory/method documents — they describe what
the shipped solvers ARE, and are the best statement of current behaviour in the corpus.
Where a DOC contradicts a SPEC, the conflict is recorded in `../INGEST-CONFLICTS.md`
rather than silently reconciled here.

---

## Solver inventory and enum map
- source: docs/methods/00-overview.md
- Eight solvers ship. `RK_ALG_AUTO = 0` is a dispatcher, not a method; slot 2 is removed
  (was LBFGSB). Live enum: ORIS 1, RAKING 3, SINKHORN 4, CHEBYSHEV 5, GREG 6,
  ORIS_SOFT 8, GREENKHORN 9, LOGIT 10, NEWTON_KL 11. Slot 7 is absent from the listing.
- Slot 7's absence is explained by `docs/superpowers/specs/2026-08-14-removed-solver-slots-supersession.md`
  (SPEC, higher precedence): `RK_ALG_GRAKE = 7` was removed pre-release in commit
  `9a67891` and the slot is permanently reserved, as is slot 2. The live `grake_norm`
  convergence metric is a different thing and is unaffected.
- The names in the codebase are described as "mostly faithful", with the mathematical
  identity read from the ACTUAL update rules rather than the labels.

## What all solvers share
- source: docs/methods/00-overview.md
- Same input (design masses `X_init`, margin incidence `A`, targets `T_kj`, optional
  per-cell box), same exit contract (`finalize_weights` enforces Σw = n FIRST then
  `bounds_mode` dispatch; renormalising after water-fill is forbidden), and shared
  infrastructure (greedy/round-robin scheduler, optional SRAA-m, optional homotopy
  cascade over `max_weight`, `CellMetrics`, `BestIterTracker`).
- ORIS, Raking, Sinkhorn, Greenkhorn and Newton-KL all converge to the SAME maximum-entropy
  (KL) solution when feasible — they are different algorithms for one optimum. GREG (χ²),
  Logit (logit distance) and Chebyshev (minimax) target DIFFERENT optima.

## Bounds-enforcement taxonomy
- source: docs/methods/00-overview.md
- The sharpest design axis, ordered: finalize-only (ORIS, Newton-KL) → inline projection
  (Raking water-fill, Sinkhorn Dykstra, Greenkhorn clamp) → active set (GREG) →
  parametrised away (Logit logistic link) → LP constraint (Chebyshev).

## Selection guidance published to users
- source: docs/methods/00-overview.md
- Robust default / severe skew → ORIS. Hard per-cell bounds with KL optimality → Raking
  (SRAA-able) or Sinkhorn (Dykstra-exact). Minimise the WORST margin error → Chebyshev.
  Classical interpretable linear weights → GREG. Strict bounds with smooth weights →
  Logit. Fast second-order on well-conditioned KL → Newton-KL. Best worst-case iteration
  complexity → Greenkhorn.

## Honest caveat on guarantees
- source: docs/methods/00-overview.md
- Verbatim position: "no method has a universal convergence proof here". Every solver
  carries engineering guards (stall reverts, dual-explosion guards, saturation early-exits,
  iteration caps) precisely because the clean theorems assume feasibility and conditioning
  that production inputs do not always meet. GRADE-style summary: strongest guarantees for
  Sinkhorn, Raking, GREG; strong-but-scoped for Newton-KL, Logit, Chebyshev; moderate for
  ORIS and Greenkhorn, whose production embellishments are empirically guarded rather than
  formally proven.

## ORIS identity and lineage
- source: docs/methods/oris.md
- ORIS is IPF/RAS/Sinkhorn–Knopp on log-Sinkhorn factors with SOR stepping and an
  infeasibility damping factor; it is the project's default solver. Canonical sources are
  Thibault–Chizat–Dossal–Papadakis 2021 (the multiplicative ω-power step, global
  convergence for ω ∈ (0,2) by Lyapunov descent, adaptive ω needing no spectral
  pre-knowledge) and Lehmann–von Renesse–Sambale–Uschmajew 2022 (local rate and
  `ω_opt = 2/(1+√(1−ϑ₂))` via Young's SOR theorem), with Soma–Uschmajew 2024 extending to
  operator scaling.
- Explicit disclaimer: arXiv:2011.14312 (iEPPA) is NOT the algorithm source. Its headline
  contribution — an outer inexact entropic-proximal-point loop for a transport-cost LP —
  is not implemented; at `C = 0` that loop is inert, and only the inner Sinkhorn-style
  scaling sweep survives, which is just IPF.

## ORIS θ₂ estimator does not port to other solvers
- source: docs/methods/oris.md, docs/methods/sinkhorn.md, docs/methods/greenkhorn.md
- The free-coordinate iterate-change observable works for ORIS because ORIS's sweep is a
  STATELESS fixed-point map, so the move cleanly reflects `ρ(M_II)`. On Sinkhorn the same
  observable conflates genuine sweep slowness with Dykstra/bisection correction grind,
  inflating the estimate toward 1, prescribing a near-ceiling ω and driving a limit cycle.
  On Greenkhorn, greedy single-coordinate selection is outside the proven full-sweep
  regime; raising the exponent above 1 was measured INERT (no iteration-count improvement,
  no divergence) across K ∈ {3,6,9}, tight and loose bounds, multiple seeds. Both therefore
  ship plain `ω = 1`.

## Three IPF schemes, two orthogonal axes
- source: docs/methods/oris.md
- Axis 1 (coordinate selection): Sinkhorn normalises all rows then all columns; cyclic
  Gauss–Seidel IPF (the ORIS core) normalises margins one at a time round-robin or
  greedy-priority; Greenkhorn normalises only the single worst-violated margin. Axis 2
  (step length): over-relaxation raises the ratio to power ω, independent of which
  coordinate is updated. ORIS is the only one that can combine greedy selection (Axis 1)
  with over-relaxation (Axis 2) — and that combination sits OUTSIDE the proven regime, with
  the per-margin ω adaptation and stall-revert guards as the engineering safety net rather
  than a theorem (stated confidence 85).

## Known IPF failure modes (literature)
- source: docs/methods/oris.md, docs/methods/raking.md
- Structural/sampling zeros: IPF preserves the seed's zero pattern by induction, so a
  positive target over an empty sample cell cannot converge. Oscillation between two
  accumulation points when margins are structurally incompatible (Pukelsheim). Bounds
  interference: if the trimmed feasible set does not intersect the marginal-constraint set,
  the bounded algorithm cannot converge. Margin inconsistency across sources: IPF has no
  mechanism to detect control totals that disagree and will cycle or fail silently.
  Category sparseness below ~1–2% creates near-zero denominators. Margin ORDERING affects
  transient behaviour; no universal ordering rule exists. Variance estimation after raking
  invalidates naive standard errors — replication methods are the practical fallback, and a
  replicate that fails to converge biases the variance estimator upward.

## How leafblower deviates from standard IPF
- source: docs/methods/oris.md
- Cell-level rather than observation-level iteration (10–100× faster at high compression,
  at the cost that all obs in a cell get an identical relative adjustment). SOR plus
  infeasibility damping instead of a plain multiplicative ratio. Bounds deferred to
  finalize rather than folded into each iteration. SRAA-m acceleration, absent from CRAN
  packages. Explicit `RK_ERR_INFEAS` with (k,j) enumeration instead of NaN propagation.
  Metric-improvement plus plateau detection instead of a bare max-error threshold. No
  variance estimation — out of scope by design.

## Raking identity
- source: docs/methods/raking.md
- Classical cyclic IPF with hard per-cell box constraints enforced INSIDE each margin step
  by `water_fill_cat`, a single-scaling projection solving
  `Σ_c clamp(X_orig[c]·m, L_c, U_c) = T_kj` by fixing bound-hitting cells and rescaling
  the free remainder until the active set stabilises. The per-margin sweep `F_eval` is
  STATELESS (no correction vectors carried), which is exactly what permits the SRAA-m
  wrapper. `inner_max_iter` is the budget; `outer_max_iter` is unused.
- Deviation from `autumn`: autumn's redistribution can push weights slightly ABOVE the cap
  (post-redistribution overflow); leafblower's free-pool rescale only touches unclamped
  cells, so clamped cells stay exactly at the bound.

## Sinkhorn identity
- source: docs/methods/sinkhorn.md
- Alternating Sinkhorn marginal rescales plus Dykstra's algorithm for the per-cell capacity
  box, with the capacity step solved by bisection on a single log-scale multiplier μ. The
  Dykstra correction `a[c]` is iterate history, which makes the solver stateful and
  therefore NOT SRAA-able. Standard alternating KL projections converge to the intersection
  only for AFFINE constraints; capacity boxes are inequalities and require the Dykstra
  correction — omitting it yields incorrect projections (Benamou et al. 2015).
- Not a general entropic-OT solver: there is no cost matrix, no Gibbs kernel and no ε, so
  it does not suffer entropic blur; the KL geometry is implicit in the marginal sweeps.

## Greenkhorn identity and its objective mismatch
- source: docs/methods/greenkhorn.md
- Greedy coordinate descent on the dual: update only `argmax_k errRp[k]`, then maintain
  `W` and `S_flat` incrementally. Deterministic greedy ⇒ stationary map ⇒ SRAA-able.
  Crucially it does NOT minimise KL directly — `best_objective` stays ∞ — so its natural
  convergence criterion is max-residual and reporting KL/χ² requires a separate tracker.
  Known pathology: ping-pong between two near-tied margins (stall_kind = 2), partially
  mitigated by the SRAA outer-stall revert. GPU-hostile: the argmax scan is sequential
  with very low arithmetic intensity; POT recommends it only for CPU and n < 10³.
- Survey calibration and statistical raking are NOT represented in the Greenkhorn
  literature as application domains; the complexity bound is inherited, not re-proved for
  the K-margin bounded-weight setting, and this is flagged as a known open gap.

## Chebyshev identity and its practical limit
- source: docs/methods/chebyshev.md
- The only interior-point/LP solver in the suite; minimises the worst-case (L∞) marginal
  error via Mehrotra predictor–corrector with a second-order centring correction, using
  reference-category elimination to break the normalisation degeneracy. Structural exploits:
  the 0/1 cell-margin incidence matrix gives a very sparse `AΘAᵀ`, the δ variable enters by
  a rank-1 Sherman–Morrison update, and no crossover is needed because calibration weights
  are continuous masses.
- Documented external-validity limit: the source itself records non-convergence on dense
  K ≥ 9 overlapping-margin systems, hard-capped at 500 iterations. Engineering departures
  from production LP solvers: no Ruiz/MC64 pre-scaling, a hard iteration cap instead of a
  strategy switch, best-iterate by calibration error rather than the LP objective, and no
  homogeneous self-dual embedding — so infeasible inputs diverge rather than being
  certified infeasible.

## GREG identity
- source: docs/methods/greg.md
- The Deville–Särndal linear (χ²) calibration estimator solved by active-set Newton, giving
  `X_new[c] = X_init[c]·(1 + Σ_k λ[k, g_k(c)])`. Refactor only on active-set change (R1
  cache); convergence is "no active-set change in one sweep" = KKT optimum. Optional
  Tikhonov ridge `st.ridge_lambda`.
- χ² weights can go negative or explode without bounds, so the box plus active set are
  essential; tight boxes cause many refactors. No anti-cycling guard — cycling is
  theoretically possible for degenerate QPs. Where active-set loses to standard NR on the
  dual: NR is O(K³) per iteration independent of M, which is cheaper when M ≫ K.

## Logit identity
- source: docs/methods/logit.md
- Bounds satisfied by construction through the logistic link, with Newton on the normal
  equations plus Armijo. Guard machinery is substantial and load-bearing: `D_eff` floor at
  `1e-6·range` (the inverse-scaling alternative was tested and REJECTED), an Armijo norm
  guard capping any `z` shift at ~2, a Layer-2 warm start rejected if any component exceeds
  10, and a 50%-saturation early exit at `|z| > 650`.
- Deviation from CALMAR-style implementations: leafblower applies the logistic link to the
  ABSOLUTE cell weight, not to the g-weight ratio `w_k/d_k`. The two parametrisations agree
  when `d_k = 1` or when design-weight information is absorbed into `L_c`/`U_c`, but differ
  numerically when per-cell bounds encode heterogeneous design weights. This is deliberate:
  it makes the bound exact in the weight scale.

## Newton-KL identity
- source: docs/methods/newton_kl.md
- Second-order attack on the smooth dual `g(λ) = log Z(λ) − T·λ` with log-sum-exp
  stabilisation, TSVD regularisation of the eigendecomposed Hessian, Levenberg–Marquardt
  damping and a Steihaug–CG trust region. Best-iterate is tracked in DUAL space
  (`lam_best`, `best_gap`) separately from the weight-space tracker. Confined by design to
  the zero-compression, non-saturating regime; hard bounds are left to finalize.
- Stated as having NO published precedent in the survey-calibration literature: no prior
  trust-region Newton dual implementation exists in the field. It applies a standard
  optimisation technique to a domain where the well-conditioned case does not need it but
  the ill-conditioned case benefits — framed explicitly as an engineering deviation, not a
  methodological one. Known TSVD caveat: hard spectral thresholding discards whole
  components and introduces bias in small-but-nonzero singular directions.

## Bounded-raking survey (research report)
- source: docs/raking.md
- Wide-scope research report on CPU-bound algorithms for bounded multi-marginal raking and
  their optimal-transport equivalences: the Deville–Särndal distance-function framework,
  the logit link `F(u) = [L(U−1) + U(1−L)exp(Au)] / [(U−1) + (1−L)exp(Au)]` with
  `A = (U−L)/((U−1)(1−L))`, the heteroskedasticity parameter `q_k` (setting it very large
  pins a unit's weight), blockwise and sequential-blockwise IPF for memory-bound settings,
  and a solver comparison table (Newton-Raphson O(p³); L-BFGS-B O(mp); projected/semismooth
  Newton; dual coordinate descent O(nnz) with `O(log(1/ε))` global linear convergence and
  good cache locality).
- CAUTION: this document's §8.2 attributes the Chu et al. inexact-entropic-proximal-point
  algorithm — outer proximal step, dual BCD subsolver, inexact stopping criteria — to
  "ORIS", which directly contradicts `docs/methods/oris.md` and the rename spec. See
  `../INGEST-CONFLICTS.md`.

## Raking hybrid critique (superseded analysis)
- source: docs/raking_bounds.md
- Prose analysis of the pre-Bregman implementation, identifying the core defect: IPF steps
  minimise KL (multiplicative) while the Dykstra box and hyperplane steps minimise
  Euclidean distance (additive), so the objective changes mid-loop and no theorem prevents
  cycling. Prescribes exactly the fix later specified as
  `2026-04-27-raking-bregman-dykstra-design.md` — initialise corrections to 1.0, replace
  `y = w + q` with `y = w·q` and `q = y − w_new` with `q = y/w_new`, and make the
  hyperplane projection a scale with `q_hyp = 1/scale`.
- Two structural notes on the analysis: it is written against an OBSERVATION-level
  implementation (`std::vector<double> q(st.n)`, `w[i]`), whereas the shipped raking solver
  is cell-level; and it also notes the absence of design-effect/variance output, which
  `docs/methods/raking.md` confirms is out of scope by design.
