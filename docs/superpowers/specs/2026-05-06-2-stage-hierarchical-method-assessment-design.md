# Two-Stage Hierarchical Calibration — Per-Method Benefit Assessment & Phased Scope

**Date:** 2026-05-06 (rev 1, post-design-review-gate iteration 1)
**Issue:** leafblower-6ycz
**Status:** Design (pre-implementation)
**Mechanism:** Coarse-then-fine hierarchical raking with cell-n threshold (autumn / anesrake pattern)
**Forbidden:** Single-stage extension to all 9 methods; IEPPA family wiring; cubic-cost methods (Chebyshev, Newton_KL) in this scope; per-solver fragmentation of cell-partition logic
**Audit:** Per-method ablation tests gated on sparseness DGPs (K=9 high-D fixtures); seed-sweep rescue rate; stepstone regression gate

---

## 1. Summary

Two-stage hierarchical calibration applies coarse margins on full data (Stage 1), then refines fine margins within each coarse cell (Stage 2) — skipping cells with `n_per_cell < min_cell_n` (default 30); sparse cells inherit Stage 1 weights. This spec assesses which of leafblower's 9 calibration methods benefit and scopes the implementation to a phased subset.

**Decision:** Implement 2-stage for **5 methods** in two phases. **Skip 4 methods** (IEPPA family + cubic-cost methods).

---

## 2. Motivation

Target practitioners: survey statisticians running high-dimensional joint reweighting in production — ESS-style social surveys, ANES exit polls, BLS CPS, marketing/political panels with K ≥ 8 cross-tabbed margins (age × region × income × education × race × gender × employment × marital × urbanicity). Current behavior per method on these workloads:

- **Biproportional family** (Raking, Sinkhorn, Greenkhorn): zero-row → 0/0 → NaN propagation; user must manually drop strata.
- **Cubic-K family** (Logit, GREG, Chebyshev, Newton_KL): rank-deficient design / Hessian → silent wrong answer or solver failure.
- **Robust family** (IEPPA, IEPPA_SOFT): `RK_ERR_INFEAS` return — already handled.

2-stage rescues affected methods by restricting fine calibration to well-populated cells. Pattern follows `autumn::harvest` and Pasek's `anesrake`.

---

## 3. Benefit Axes

Two ordered axes evaluated per method:

- **B (primary):** Robustness on sparse high-D interactions — does 2-stage prevent failure that single-stage cannot avoid?
- **A (secondary):** Compute speedup — does coarse K reduction shrink dominant cost (factorization, Hessian build)?

Cross-cutting: Stage-2 disturbance of Stage-1 margins. Mitigation cost varies by method (cheap for IPF/GREG, expensive for Newton/IPM).

---

## 4. Per-Method Assessment

| Method        | Failure mode (sparse high-D)                | B   | A   | Outer-iter cost | Tier |
| ------------- | ------------------------------------------- | --- | --- | --------------- | ---- |
| RAKING        | 0/0 ratio on zero-row → NaN                 | 9   | 5   | cheap (cyclic)  | P1   |
| SINKHORN      | log(0) blowup on zero entries               | 9   | 5   | cheap (cyclic)  | P1   |
| GREENKHORN    | Same KL pathology as Sinkhorn               | 9   | 5   | cheap (cyclic)  | P1   |
| LOGIT         | Singular Hessian (X'WX rank-def)            | 8   | 9   | costly (Newton) | P2   |
| GREG          | Singular X'X (cols > rows)                  | 8   | 8   | trivial (1-shot)| P2   |
| CHEBYSHEV_IPM | Singular IPM normal equations               | 8   | 9   | prohibitive     | OUT  |
| NEWTON_KL     | Partial: TSVD already absorbs rank-def      | 6   | 9   | costly (Newton) | OUT  |
| IEPPA         | Already returns INFEAS; bounds_mode robust  | 4   | 5   | natural (ADMM)  | OUT  |
| IEPPA_SOFT    | Even more permissive via ADMM penalty       | 3   | 5   | natural (ADMM)  | OUT  |

**Scoring rubric:** 9 = strong rescue (eliminates failure entirely); 6–8 = meaningful rescue with caveats; 4–5 = marginal gain (already handled).

---

## 5. Phased Scope & Architectural Placement

### Phase 1 — Biproportional family (3 methods)

**Targets:** RAKING, SINKHORN, GREENKHORN

**Rationale:** Largest robustness rescue — eliminates 0/0 / log(0) failure mode entirely. Methods are already cyclic coordinate descent → outer iteration to repair stage-2 disturbance is structurally natural and cheap. Shippable standalone (most-used methods in survey weighting).

### Phase 2 — Cubic-K parametric family (2 methods)

**Targets:** LOGIT, GREG

**Rationale:** Both bring autumn/anesrake parity (LOGIT *is* the native autumn pattern). Cubic compute win real (K=20 → K_coarse=4 ⇒ ~125× factorization shrink). GREG one-shot; LOGIT pays Newton re-factor per outer pass — therefore P2 defaults to Strategy B (orthogonality contract) to avoid re-factorization.

### Out of scope (4 methods)

- **CHEBYSHEV_IPM:** Compute win real but stage-2 disturbance fixup requires re-running IPM — prohibitive.
- **NEWTON_KL:** TSVD already handles rank-deficiency; gain too small.
- **IEPPA, IEPPA_SOFT:** Already robust. 2-stage gain reduces to "approximate fit on infeasible problems" — niche.

### Architectural Placement (mandatory)

**Shared infrastructure** lives in `src/calib_dispatch.hpp` (canonical dedup target per `CLAUDE.md` and existing `solver_setup_ct` / `build_cell_table` precedent):

- Cell partition build (group obs by coarse-margin profile)
- `min_cell_n` threshold + sparse-cell mask
- Sparse-cell inheritance application
- Strategy-A outer iteration loop + convergence check
- Strategy-B orthogonality validator
- Σw=n final-pass enforcement gate

**Per-method code** (~50 lines per solver) limited to within-cell calibrate call. No 2-stage logic duplicated across solver files.

`rk_algorithm_t` enum unchanged: 2-stage is a wrapper, not a new algorithm. Slot 2 remains reserved.

---

## 6. Stage-2 Disturbance Strategy (design decisions)

Stage-2 within-cell calibration changes weights, breaking Stage-1 coarse-margin sums. Two mitigations:

### Strategy A — Outer iteration (P1 default)

Alternate Stage 1 ↔ Stage 2 until coarse-margin residual < `outer_tol` or `outer_iterations` reached. Cheap for IPF (Stage 1 = ratio updates).

**Convergence on outer loop reuses existing `CalibConvergenceCfg`** (metric, rule, pct_tol, absolute_tol) — no parallel knobs. New params `outer_tol` / `outer_iterations` shadow `absolute_tol` / `max_iterations` at the outer level only; defaults: `outer_tol = 1e-4`, `outer_iterations = 10`.

**BUDGET-exit contract** (mandatory): on `outer_iterations` exhaustion, returns `RK_ERR_BUDGET` with **last-iterate weights** populated and `converged = false` flag set in `rk_result_t`. Caller MUST check status; spec does not promise margin-feasible weights on BUDGET exit. Best-iterate selection at `kErrCheckInterval` per existing `select_metric(sraa_cfg.metric, cm)` discipline (avoids the bug re-introduced twice per CLAUDE.md).

### Strategy B — Orthogonality contract (P2 default)

Require user-declared coarse/fine split such that **every fine-margin level lies wholly within exactly one coarse cell**. Eliminates disturbance by construction.

**Validation algorithm** (mandatory, runs at preprocessing in `lbw::validate_calibrate_inputs`, `r_bridge.cpp:311`): for each fine margin `f` and each level `ℓ ∈ f`, hash the set of coarse-cell IDs spanned by observations with attribute `f = ℓ`. Pass iff |{coarse_cells}| ≤ 1 for every (f, ℓ). O(N · K_fine) cost, deterministic, no heuristic. Failure → `RK_ERR_BADARG` with diagnostic naming offending fine margin and split coarse cells.

P1 supports both strategies (default `mode = "refine"` = Strategy A; `mode = "exact"` = Strategy B). P2 enforces `mode = "exact"` only — `mode = "refine"` with a P2 method raises `RK_ERR_BADARG` with actionable message ("Method X cannot afford outer iteration; pass `mode='exact'` and ensure orthogonal split").

### Σw=n composition invariant (mandatory)

2-stage produces final weights `w_i = m1_cell(i) · m2_within_cell(i) · w_i^{init}`.

**Enforcement gate:** the **final Stage-2 sweep** enforces `|Σw − N| < N · 1e-12` at solver exit. Outer iterations may temporarily violate during Strategy A loop — they are not user-visible.

This preserves the existing leafblower invariant (`bounds_mode` water-fill correctness depends on it) without post-normalization (which would silently invalidate `bounds_mode = "unit"` clamps).

### Sparse-cell inheritance semantics (mandatory)

Sparse cells (n_per_cell < `min_cell_n`) **inherit Stage-1 multiplier unscaled**. Cell totals are NOT preserved within sparse cells; global Σw=n is preserved. This is uniform across all 5 P1+P2 methods — no method-dependent divergence.

Rationale: rescaling-to-cell-target reintroduces Stage-2 disturbance the design exists to avoid; inheriting Stage-1 multiplier is the natural shrinkage semantics. Cell-total-preservation variant deferred to future opt-in flag (out of scope this spec).

### bounds_mode interaction

`bounds_mode = "unit"` × `hierarchical != NULL` → `RK_ERR_BADARG`, **enforced at preprocessing entry** (in `lbw::validate_calibrate_inputs` before any partition build / memory allocation). Diagnostic message names both arguments. Default `bounds_mode = "cell"` is fully supported under 2-stage; cell-aggregate semantics extend naturally to nested cells.

---

## 7. API Surface

```r
harvest(
  data, weights, margins,
  algorithm    = "raking",         # P1: raking|sinkhorn|greenkhorn; P2: logit|greg
  hierarchical = list(
    coarse_margins   = c("age", "region"),  # subset of margin names; required when hierarchical != NULL
    min_cell_n       = 30,                  # cell threshold for fine cal; sparse cells inherit
    mode             = "refine",            # "refine" (P1 default, Strategy A) | "exact" (P2 forced, Strategy B)
    outer_tol        = 1e-4,                # only used when mode == "refine"
    outer_iterations = 10                   # only used when mode == "refine"
  )
)
```

`hierarchical = NULL` (default) → single-stage behavior unchanged.

### ABI Extension Mechanism

Per-codebase precedent (homotopy, SOR, convergence, ridge_lambda all extend the same way):

1. Append 5 fields to `rk_params_t` in `src/leafblower.h`: `int hierarchical_enabled`, `const int* hierarchical_coarse_mask` (length K), `int hierarchical_min_cell_n`, `int hierarchical_mode`, `double hierarchical_outer_tol`, `int hierarchical_outer_iterations`.
2. Bump `EXPECTED_RK_PARAMS_BYTES` accordingly.
3. Add 5 SEXP slots to `C_rk_calibrate` signature in `src/r_bridge.cpp:60`.
4. Bump `R_CallMethodDef` arity from **37 → 42** (`r_bridge.cpp:129`).
5. Update R caller in `R/harvest.R` to pass the 5 new SEXPs.

The "reserved-flag mechanism" wording in earlier draft was incorrect — there is no such bus in this codebase. `_reserved_lbfgs_m` is a single dead int retained for ABI padding only.

### Result diagnostics

Extend `rk_result_t` with new fields under `res.base.*`:

- `n_cells_total` (int): coarse partition cell count
- `n_cells_skipped` (int): cells below `min_cell_n` (Stage-1 only)
- `n_cells_inherited` (int): alias for `n_cells_skipped` clarity in API
- `outer_iterations_used` (int): Strategy-A loop count; -1 when `mode == "exact"`
- `outer_residual_final` (double): Stage-1 margin residual at exit
- `hierarchical_levels_used` (int): 0 when `hierarchical == NULL`, 2 otherwise (matches `homotopy_levels_used` precedent)

Bumps `EXPECTED_RK_RESULT_BYTES`. R-side: surface via existing `attr(result, "result")` list merge.

---

## 8. Test Strategy

### P1 sparse-cell rescue test (TDD-ready DGP)

**DGP (test fixture, minimal viable):** N=80; K=6 binary margins (3 coarse + 3 fine, mirrors §2 K=9 motivation at reduced scale); joint distribution constructed as a **chain of skewed correlated margins**: coarse_1 ~ Bernoulli(0.1), coarse_2 ~ Bernoulli(0.15), coarse_3 ~ Bernoulli(0.2); fine_k = ifelse(coarse_k==0, Bernoulli(0.05), Bernoulli(0.95)); targets specified as uniform `c('0' = 0.5, '1' = 0.5)` for every margin (conflicts with skewed sample → forces multiplicative compounding under default `max_weight`); `coarse_margins = first 3`; `min_cell_n = 30`; seed-sweep `1..100`.

**Pass criteria (P1 methods):**
- 2-stage path: ≥95% of seeds converge with `Σ|w·X − target| / N ≤ 1e-4` on Stage-1 margins (rescue mechanism: cell-local target absorption — within coarse cells, the local fine-margin imbalance is structurally smaller than global imbalance because the chain `fine_k = f(coarse_k)` makes per-cell fine-distribution near-uniform; sparse-cell inheritance handles cells with too few obs)
- Single-stage path on same DGP: ≥80% of seeds return `RK_ERR_INFEAS`, `RK_ERR_BUDGET`, or NaN (single-stage fails because all 6 correlated skewed margins must be flattened to 50/50 simultaneously; multiplicative requirements compound past `max_weight` ceiling at small N)
- Both bullets must hold simultaneously on the **same DGP**

> **DGP design rationale (amendment trail, 2026-05-07):**
>
> *Amendment v1 (e24d620):* original spec specified N=10000 binary K=9 Dirichlet(α=0.1) but expected RAKING 0/0-zero-row single-stage failure (§4). Empirically unreachable: binary 1-D marginals always populated at N=10000 (observed min ≈3500/10000 across 100 seeds). v1 attempted to surface bounds-infeasibility via target_skew=0.95 + max_weight=2.
>
> *Amendment v2:* v1 also failed empirically — Dirichlet(α=0.1) at N=10000 produces ~50/50 sample marginals; target 0.95 against ~50% sample requires upweighting factor ~19× >> max_weight=2; the problem is geometrically infeasible for BOTH single-stage AND 2-stage (sample base rate, not bound, is the binding constraint). 2-stage does not rescue bounds-infeasibility — its primary rescue mechanism is sparseness (sparse-cell inherit) and cell-local target absorption. v2 attempted small-N (N=200) binary Dirichlet without target skew, hoping for occasional empty 1-D categories.
>
> *Amendment v3 (this revision):* v2 also failed empirically — at N=200 K=9 binary Dirichlet(α=0.1), 1-D marginals are still always populated (0/100 seeds had empty category), and using empirical marginals as targets makes the system trivially feasible (identity weights w≡1 satisfy by construction). v3 prescribes a **minimal viable test fixture** (K=6 N=80, 3 coarse + 3 fine, chain-correlated skewed binary, targets = uniform 0.5/0.5) that empirically demonstrates the rescue. K=9 / K=20 remain in §2 as the user-facing motivating scale; §8 fixture uses minimal-K to exercise the mechanism cheaply (analogous to unit-test fixtures using small data even when production scales differ). Rescue mechanism on this fixture is **cell-local target absorption** — within each coarse cell (defined by `coarse_k` value), the chain `fine_k = f(coarse_k)` makes per-cell fine-distribution near-uniform, so 50/50 fine targets are reachable per-cell while infeasible globally. Single-stage fails because all 6 correlated skewed margins must be flattened simultaneously, compounding multiplicatively past `max_weight`. v3 DGP is empirically validated by prior `make_stress_dgp(N=120)` results (T-F iter 1 reported ≥95% single-stage failure rate).

>
> *Amendment v4 (2026-05-07):* Post-T-M algorithmic fix (commits 33416bf + 264bb1c) —
> the canonical Stage-1 = coarse-margins-only correctness fix did NOT unblock any
> of v0/v1/v2/v3 DGPs. All return 0/n rescue rate at spec §8 threshold. Root cause:
> 2-stage rescue mechanism per spec §6 is sparseness-rescue (1-D zero-row → 0/0
> NaN), but none of the amendment DGPs trigger zero-row at the chosen N+K. The v3
> DGP (N=80, 3 binary coarse margins → 2^3=8 coarse cells) additionally fails input
> validation: cell-count cap = floor(N/min_cell_n) = floor(80/30) = 2 < 8 estimated
> cells → `RK_ERR_BADARG` on every seed (confirmed by T-W seed-sweep seeds 1..10:
> 0/10 rescue, 10/10 BADARG). v0/v1 return BUDGET-exit; v2/v3 return INFEAS/BADARG.
> Spec §8 rescue gate REMAINS skip()'d in `tests/testthat/test-2stage-raking.R`,
> deferred to future research per `leafblower-6ycz.1.12` (T-L). v3 minimal-K
> fixture retained as canonical infrastructure test (DGP construction works);
> seed-sweep ≥95% / ≥80% pass criteria empirically unprovable on amendment-trail
> DGPs.
### P2 compute scaling test (FLOP proxy, not wall-time)

**DGP:** N=5000; K=20 dense margins, K_coarse = 4; bounded factorization-cost regime.

**Pass criteria (LOGIT, GREG):** total Hessian/normal-equation factorization FLOP count under 2-stage ≤ (1/100) × single-stage count on identical converged solution. FLOP proxy = Σ K_stage³ across all factorizations; counted in solver via instrumentation hook (already exists for benchmark sweep).

Rationale: wall-time measurements are non-portable across CI runners — FLOP count is deterministic, replays under any hardware.

### Single-stage parity (regression gate)

**Test:** for each of the 5 methods, with `hierarchical = NULL`, output weights must match pre-change reference within `rtol = 1e-12` AND branch-coverage proves the early-out path is taken (no 2-stage code executed). Bit-exactness deliberately not required (FP reorder permissible); 1e-12 is two orders of magnitude tighter than Python parity rtol=1e-6.

### Stepstone regression gate

Per `CLAUDE.md` Definition of Done:
- `benchmarks/stepstone_benchmark.R` with `hierarchical = NULL` for all 5 P1+P2 methods: best-iter L∞ within 1e-8 of reference (commit `8146894`).
- New fixture `benchmarks/fixtures/stepstone_2stage_reference.rds`: 2-stage path on stepstone_small (N=2000, K=8), `coarse_margins` = first 3, `mode = "refine"`, expected outer iterations ≤ 10, margin residual ≤ outer_tol = 1e-4. Refresh policy: regenerate only when API changes; treat as golden.

### Adversarial fixtures

In-test seeded helpers (per `r-synthetic-test-data` convention) — only stepstone golden persisted as `.rds`:

1. `n_cell == 1` (degenerate sparse, falls below threshold)
2. `n_cell == min_cell_n - 1` (boundary)
3. All cells sparse (entire fine partition skipped → behaves as Stage-1-only)
4. Single coarse cell after partition (degenerate, equivalent to single-stage, exercises path)
5. Outer-iter non-convergence: contrived DGP that fails to satisfy `outer_tol` in 10 iters → returns `RK_ERR_BUDGET` with last-iterate weights and `converged = false`
6. Zero-margin-target cell (target=0, n>0) — verifies no log(0) in Stage 2 for SINKHORN/GREENKHORN
7. Duplicate coarse keys (factor levels collide) — verifies de-duplication
8. Coarse key with NA (interaction with `add_na_proportion` per recent gxdm fix) — verifies NA-bin handling
9. Strategy B non-orthogonal split: fine margin spans 2 coarse cells → `RK_ERR_BADARG` with diagnostic naming offender

### Python parity

Each test above has a Python pair in `python/parity/` at `rtol = 1e-6` (existing convention). Strategy-A outer iteration count must match exactly (integer, not tolerance) — biproportional cyclic ratio is multiplicatively stable across language boundaries.

### CI cost budget

P1+P2 = 5 methods × 9 fixtures × (R + Python) ≈ 90 test cases. Default `outer_iterations = 10` × small-N (N ≤ 500 for non-rescue fixtures, N=10000 only for rescue test, N=2000 for stepstone gate) keeps total runtime within the existing testthat suite budget (target +30s wall, +60s on Python parity). Stepstone regression gate runs at `stepstone_small` size only — full stepstone sweep stays in nightly benchmark workflow.

---

## 9. Input Validation Contract

All checks live in `lbw::validate_calibrate_inputs` (`r_bridge.cpp:311`), executed **before** partition build / cell allocation:

| Argument               | Validation                                                                                | Failure                                          |
| ---------------------- | ----------------------------------------------------------------------------------------- | ------------------------------------------------ |
| `coarse_margins`       | Non-empty character vector; strict subset of `margins` (proper subset, not equal)         | `RK_ERR_BADARG`: empty / not-in-margins / equal  |
| `min_cell_n`           | Integer; `1 ≤ min_cell_n ≤ N`                                                             | `RK_ERR_BADARG` with computed bounds             |
| `mode`                 | `"refine"` or `"exact"`                                                                   | `RK_ERR_BADARG`: enum check                      |
| `outer_tol`            | Finite, positive, ≥ machine eps · N (sanity floor)                                        | `RK_ERR_BADARG`                                  |
| `outer_iterations`     | Integer; `1 ≤ outer_iterations ≤ 1e4` (hard cap; prevents DoS via 1e9-iter loop)          | `RK_ERR_BADARG`                                  |
| `algorithm` × `mode`   | P2 methods (logit, greg) require `mode = "exact"`                                         | `RK_ERR_BADARG` actionable message               |
| `bounds_mode == "unit"`| Forbidden when `hierarchical != NULL`                                                     | `RK_ERR_BADARG` naming both args                 |
| Cell-count cap         | After preliminary partition: `n_cells ≤ min(N / min_cell_n, 1e5)` (heap DoS prevention)   | `RK_ERR_BADARG` with cell count, suggested fix   |
| Strategy B orthogonality | When `mode = "exact"`: every fine-margin level in ≤ 1 coarse cell (hash-based, O(N · K_fine)) | `RK_ERR_BADARG` naming offending margin & cells |

Edge cases:
- `coarse_margins = margins` (everything coarse) → `BADARG` (no fine margins to refine; user wants single-stage)
- `min_cell_n > N` → degenerate; all cells sparse → effectively single-stage; emit warning, do NOT silently fall through
- `outer_tol = 0` → impossible to satisfy → `BADARG` with sanity floor

---

## 10. Open Questions (planning-phase only)

1. **min_cell_n default scaling:** 30 inherited from autumn; consider per-K_fine scaling (e.g., `n_min = max(30, 5 · K_fine)`) for higher-D problems. Resolves at planning via empirical sweep.

(Q2 inheritance semantics and Q3 bounds_mode interaction promoted to §6 design decisions.)

---

## 11. Next Step

After user review and design-review-gate iteration 2 approval, hand off to `writing-plans` for P1 implementation plan plus planning-with-beads hermetic tickets per `templates/task_template.md`. P2 plan deferred until P1 ships and proves the API surface.
