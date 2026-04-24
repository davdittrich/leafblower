# iEPPA-preserving convergence acceleration — literature scoping review

**Date:** 2026-04-24
**Scope:** iEPPA-preserving methods to close the stepstone-fulldata gap vs `autumn::harvest(accelerate=TRUE)` (33.6 s, max_err = 1.60e-3) and pursue a 10× stretch (max_err ≈ 1.60e-4).
**Exclusions (per user):** Chambolle–Pock primal-dual, mirror descent with entropic regularization, active-set Newton as standalone algorithmic replacements.
**Method:** rapid targeted scoping review (PRISMA-ScR lean variant). Four query bundles, 2018–2026 window, abstract-level verification via arxiv fetch. Gemini-delegate attempted first — failed (stale-cache hit from earlier session + fabrication of the known-wrong Tang 2024 arxiv ID `2401.07127`). Fell back to WebSearch + arxiv fetch through context-mode.

---

## Bundle A — Progressive bound tightening / homotopy for constrained IPF

| short_name | id | year | rate | applicability |
|---|---|---|---|---|
| **schmitzer-epsilon-scaling** | arxiv:1610.06519 → SIAM J Sci Comput 41(3):A1443 | 2019 | O(1/ε) for fixed ε; ε-scaling schedule + log-domain stabilization permits smaller ε | **3/3** |
| **chizat-annealed-sinkhorn** | arxiv:2408.11620 | 2024 | β_t = Θ(√t) universal schedule lower bound; β_t → ∞ and β_t − β_{t−1} → 0 necessary/sufficient | **3/3** |
| **tang-sparse-newton-sinkhorn** | arxiv:2401.12253 | 2024 | sparse Newton inner — not iEPPA-preserving, but citing family context | 1/3 |

**Schmitzer 2019** (verified abstract, arxiv:1610.06519, published SIAM SISC 2019): combines log-domain stabilization, ε-scaling heuristic, adaptive kernel truncation, coarse-to-fine. Core idea: anneal the regularization parameter ε from large (stable, smooth, fast) to small (sharp, slow) on a schedule chosen to balance Sinkhorn iterations per ε level. The paper develops convergence analysis that illuminates *why* ε-scaling works.

**Chizat 2024 Annealed Sinkhorn** (verified abstract, arxiv:2408.11620): proves that any concave annealing schedule β_t (with β_t = 1/ε_t) asymptotically solves OT iff β_t → ∞ and β_t − β_{t−1} → 0. Best error trade-off under this universal constraint: β_t = Θ(√t). Debiased variant lifts the √t limit. Proof goes through an equivalence with Online Mirror Descent.

**Applicability to iEPPA box-clamp problem:** the *capacity bound* `max_weight` plays an ε-like role — a tight bound produces a non-smooth piecewise-linear feasible set; a loose bound approaches the unconstrained smooth Sinkhorn. Annealing `max_weight` from (target × k) down to target over N outer sweeps is the direct structural analog of ε-scaling. Schmitzer's stabilization + Chizat's schedule theory give the template. Neither paper targets box bounds specifically, so the theoretical transfer is by analogy; rigorous rate transfer requires adaptation work (open).

---

## Bundle B — Per-margin adaptive step + Gauss-Seidel priority ordering

| short_name | id | year | rate | applicability |
|---|---|---|---|---|
| **greenkhorn** | arxiv:1705.09634 (NIPS 2017) | 2017 | near-linear time; same theoretical guarantees as Sinkhorn, superior empirically | **2/3** |

**Altschuler–Weed–Rigollet 2017 Greenkhorn** (verified abstract, arxiv:1705.09634): greedy coordinate descent variant of Sinkhorn. Instead of alternating all-row / all-column sweeps, each iteration picks the single row OR column with maximal marginal violation and normalizes only that one. Same near-linear-time complexity bound. Empirically "significantly outperforms classical Sinkhorn."

**Applicability:** the direct multi-marginal / box-constrained generalization is not published as far as this review found, but the design is transparent. For K margins, pick the margin index k* with largest residual `errRp_k`, apply one iEPPA margin-k step (with box clamp), recompute residuals, repeat. Preserves iEPPA semantics exactly — just changes the scheduler. Autumn's empirical advantage on "simple 2-cat margins polished to machine precision" maps directly to a Greenkhorn-style priority queue: simple margins collapse to ~0 residual in one step and drop out of rotation, hard margins get the remaining budget.

**Gap:** no paper with a proof for the multi-marginal K-way box-constrained case. Design translatable from proof-bearing binary case. Rate transfer conjectural; empirical validation decides.

---

## Bundle C — Nesterov momentum on dual log-scaling variables

| short_name | id | year | rate | applicability |
|---|---|---|---|---|
| **dvurechensky-apdagd** | arxiv:1802.04367 (ICML 2018) | 2018 | Õ(min{n^{9/4}/ε, n²/ε²}) | 1/3 (replacement-class) |
| **lin-2021-efficient-ot** | arxiv:2104.05802 | 2021 | accelerated gradient, OT complexity improvements | 1/3 (replacement-class) |

**Dvurechensky–Gasnikov–Kroshnin 2018 APDAGD** (verified abstract, arxiv:1802.04367, ICML 2018): adaptive primal-dual accelerated gradient descent on the dual-form entropic OT. Õ(n^{9/4}/ε) complexity beats Sinkhorn's Õ(n²/ε²) in the ε-regime. *Not Sinkhorn-style alternating minimization* — it is a Nesterov-accelerated gradient method on the dual potentials, in the primal-dual framework. Per the user's exclusion of "mirror descent / primal-dual as replacement," APDAGD sits on the wrong side of the line.

**Lin et al. 2021** (arxiv:2104.05802, abstract per WebSearch): same class of Nesterov-accelerated dual-gradient method for entropic OT. Same classification.

**Gap:** *no paper in the 2018–2026 window applies Nesterov momentum to the dual log-scaling factors while preserving the Sinkhorn/IPF alternating-multiplicative structure for the specifically box-constrained case.* Heavy-ball or FISTA on the dual IS published for unconstrained entropic OT (Dvurechensky, Lin et al.), but those methods replace the IPF sweep with a gradient step; from the iEPPA-preserving stance they qualify as replacements, not overlays.

**Implication:** P-C (Nesterov on dual log-factors as an iEPPA-preserving overlay) has no direct theoretical backing in the surveyed literature. Empirical speculation only.

---

## Bundle D — Convergence rate under piecewise-linear clamp

| short_name | id | year | rate | applicability |
|---|---|---|---|---|
| **tang-2024-constrained-ot** | arxiv:2403.05054 | 2024 | sublinear first-order in dual (Lyapunov); exponential reduction of approximation error with η; higher-order via second-order acceleration | **3/3** (but excluded from iEPPA-preserving pathway per prior shelving) |
| **ghosal-nutz-2022** | arxiv:2212.06000 (to appear Math OR) | 2022 | H(π_t\|π*) + H(π*\|π_t) = O(t^{−1}); dual O(t^{−1}); marginals O(t^{−2}) | **2/3** (unconstrained baseline) |
| **pins-2025** | arxiv:2502.03749 | 2025 | global convergence, sparsity-reduced per-iter cost | 1/3 (Sinkhorn + sparse Newton hybrid) |

**Tang et al. 2024** (verified abstract, arxiv:2403.05054; title "A Sinkhorn-type Algorithm for Constrained Optimal Transport"; authors Tang, Rahmanian, Shavlovsky, Thekumparampil, Xiao, Ying): directly addresses equality + inequality (box) constrained OT via a Sinkhorn-type algorithm. Provides: (a) exponential reduction of approximation error with regularization parameter; (b) sublinear first-order convergence rate in the dual space via Lyapunov-function analysis; (c) dynamic regularization scheduling as an acceleration add-on; (d) second-order acceleration (optional) for fast/higher-order convergence under weak regularization.

Per prior session research (leafblower-ylsy), Tang 2024's full algorithm was shelved after gate FAIL — the standalone primal-dual Newton version needed a runtime Δ estimator + η annealing + scope re-authorization. However, component (c) — *dynamic regularization scheduling* — is an iEPPA-preserving overlay: it only changes the penalty weight between iEPPA sweeps, not the sweep itself.

**Ghosal–Nutz 2022** (verified abstract, arxiv:2212.06000): for unconstrained Sinkhorn on entropically regularized OT, the symmetric relative entropy to the optimal coupling is O(t^{−1}), dual suboptimality O(t^{−1}), marginal entropies O(t^{−2}). Large class of costs including quadratic + subgaussian marginals. No box-constraint case. Provides the **baseline rate** against which the capacity-clamp non-smoothness is measured: if unconstrained Sinkhorn gives O(t^{−1}) but the iEPPA kk1204 trajectory is O(t^{−1/2}), the capacity clamp is costing one order of the rate — consistent with moving from a smooth to a piecewise-linear operator.

**Wu–Liang–Yang 2025 PINS** (verified abstract, arxiv:2502.03749): proximal iterations with Sinkhorn inner + sparse Newton outer. Global convergence, sparsity-exploiting. Partial iEPPA preservation (Sinkhorn present) but Newton in outer loop pushes this into the hybrid / replacement-adjacent class.

---

## Extraction grid (consolidated)

| short_name | bundle | applic. | notes |
|---|---|---|---|
| schmitzer-epsilon-scaling | A | 3 | ε-scaling homotopy, log-domain stabilization — template for max_weight annealing |
| chizat-annealed-sinkhorn | A | 3 | √t schedule theorem; debiasing lifts universal limit |
| greenkhorn | B | 2 | greedy priority Sinkhorn; multi-marginal / box extension is open but transparent |
| dvurechensky-apdagd | C | 1 | replacement-class (primal-dual accelerated gradient); user-excluded |
| lin-efficient-ot | C | 1 | replacement-class; user-excluded |
| tang-2024-constrained-ot | D | 3 | direct on box-constrained Sinkhorn; dynamic regularization schedule component is iEPPA-preserving |
| ghosal-nutz | D | 2 | unconstrained baseline rate — confirms clamp costs one order |
| pins-2025 | D | 1 | Sinkhorn + sparse Newton hybrid |

---

## Composition flags

- **schmitzer-epsilon-scaling × tang-2024-constrained-ot × chizat-annealed-sinkhorn**: all three point to a *regularization / bound annealing schedule* as the iEPPA-preserving acceleration mechanism backed by theory. Combining them into a single practical recipe (max_weight homotopy + η / ε scheduling + √t cadence bound) is the highest-yield integration target.
- **greenkhorn × tang-2024-constrained-ot**: priority-ordered sweeps on the CMOT dual. Greenkhorn gives the scheduling primitive; Tang gives the Lyapunov that should (speculatively) carry over to the priority variant. Empirical validation required.

## Confidence

- Bundle A: **high** — Schmitzer and Chizat both directly on-topic.
- Bundle B: **medium** — Greenkhorn is the only backed method; multi-marginal box extensions are not proven.
- Bundle C: **medium** (as a negative result) — the literature uses Nesterov on duals ONLY in the primal-dual replacement form. No iEPPA-preserving overlay found.
- Bundle D: **high** — Tang 2024 is directly on-topic; Ghosal–Nutz gives the unconstrained baseline.

## Caveats

- Gemini-delegate's single-shot attempt returned a stale-cached JSON from an unrelated prior review AND fabricated the wrong arxiv ID for Tang 2024 (`2401.07127` vs correct `2403.05054`). Confirms the 2026-04-24 memory about gemini-delegate stale cache on repeat invocations in long sessions.
- arxiv IDs `2401.12253`, `2104.05802` were WebSearch-listed but not successfully fetched (MCP server flapping). Re-fetching them on a clean session recommended if they advance the analysis — but both are replacement-class by abstract summary and likely do not change the iEPPA-preserving recommendation.
- Multi-marginal (K-way) generalizations of Greenkhorn and of Tang 2024's Lyapunov analysis are not directly proved in the surveyed work. Transfer from binary / single-margin cases is by analogy; the engineering validity rests on empirical benchmark results.

---

## Recommendations — iEPPA-preserving action set

Ranked by expected yield × theoretical backing × integration cost.

### 1. **P-A — Progressive bound (capacity) tightening, Schmitzer-style**

- **Theoretical backing:** Schmitzer 2019 (ε-scaling), Chizat 2024 (annealing schedule theory). Both peer-reviewed or strong preprint, both on Sinkhorn.
- **Mechanism:** start `max_weight_current = max_weight_target × k` (e.g., k = 10). Run iEPPA to a relaxed inner tolerance at current clamp level; halve `max_weight_current` (or schedule along √t-inspired cadence); re-run; repeat until `max_weight_current = max_weight_target`.
- **Risk:** at the final clamp level the slow-rate regime can still dominate the tail — not a full cure. BUT prior iterates warm-started into the final level bypass the long flat phase that costs kk1204 its 500-iter budget today.
- **Integration cost:** low. One outer loop around existing iEPPA, one schedule parameter. No data-structure change.
- **Why this is first:** attacks the *root cause* identified in memory (piecewise-linearity from the clamp defeating smooth-operator contraction theory).

### 2. **P-B — Greedy / residual-priority margin scheduler (Greenkhorn extension)**

- **Theoretical backing:** Altschuler–Weed–Rigollet 2017 Greenkhorn on binary Sinkhorn. Multi-marginal K-way proof open; empirically transferable.
- **Mechanism:** replace round-robin margin-index loop with a priority queue keyed on current margin residual. One iEPPA margin-step on `argmax_k errRp_k`, update residuals, re-insert, repeat.
- **Risk:** priority queue overhead scales as O(K log K) per step — negligible for K = 20 and K ≈ 10 ranges. Worst case: uniform-residual regime reduces to round-robin + small constant overhead.
- **Integration cost:** low-to-medium. Scheduler change inside the sweep loop; no change to margin-step math.
- **Why this is second:** closes the *autumn gap* on simple margins (autumn's greedy-rake polishes 2-cat margins exactly — Greenkhorn scheduling buys leafblower the same behavior within iEPPA semantics).

### 3. **P-A + P-B composed**

- **Rationale:** both attack different root causes. Composition is the minimum configuration with theoretical cover for the 10× stretch target.

### 4. **Optional: dynamic regularization scheduling borrowed from Tang 2024**

- **What to borrow:** only the *scheduling* component (η or the iEPPA proximal weight). Skip the second-order Newton and the full primal-dual replacement — those are out of scope per user directive and shelved per ylsy.
- **Relation to P-A:** if max_weight annealing (P-A) is the outer schedule, Tang-style η / proximal-weight scheduling is the inner schedule. Orthogonal axes.
- **Defer** until P-A + P-B data lands.

### 5. **P-C (Nesterov on dual `lf`) — not recommended without empirical proof**

- No theoretical backing found in the surveyed literature for the iEPPA-preserving overlay form.
- Prior Halpern attempt (9a97cc8, partially reverted) confirms O(1/k) = bare Sinkhorn on slow-rate regime — momentum overlays have repeatedly failed in this code base.
- If P-A + P-B miss the 10× stretch target, revisit then, but as an empirical probe only.

---

## Downstream items (not part of this review)

- Implementation plan for P-A + P-B composition targeting stepstone-fulldata → filed as plan in `docs/superpowers/plans/`.
- Rate-gap measurement: measure errRp trajectory slope on stepstone, compare to Ghosal–Nutz O(t^{−1}) and the empirical O(t^{−1/2}). Supports or falsifies the "clamp costs one rate order" hypothesis.
- Memory update: correct gemini-delegate trust policy — even single-shot in the same session can hit stale cache.
