# Two-Stage Hierarchical Calibration — Per-Method Benefit Assessment & Phased Scope

**Date:** 2026-05-06
**Issue:** leafblower-6ycz
**Status:** Design (pre-implementation)
**Mechanism:** Coarse-then-fine hierarchical raking with cell-n threshold (autumn / anesrake pattern)
**Forbidden:** Single-stage extension to all 9 methods; IEPPA family wiring; cubic-cost methods (Chebyshev, Newton_KL) in this scope
**Audit:** Per-method ablation tests gated on sparseness DGPs (K=9 high-D fixtures)

---

## 1. Summary

Two-stage hierarchical calibration applies coarse margins on full data (Stage 1), then refines fine margins within each coarse cell (Stage 2) — skipping cells with `n_per_cell < n_min` (default 30); sparse cells inherit Stage 1 weights. This spec assesses which of leafblower's 9 calibration methods benefit from the 2-stage extension and scopes implementation to a phased subset.

**Decision:** Implement 2-stage for **5 methods** in two phases. **Skip 4 methods** (IEPPA family + cubic-cost methods).

---

## 2. Motivation

High-dimensional joint calibration (K ≥ 8 dense margins) produces inevitable cell sparseness in cross-product strata. Current behavior per method:

- **Biproportional family** (Raking, Sinkhorn, Greenkhorn): zero-row → 0/0 → NaN propagation
- **Cubic-K family** (Logit, GREG, Chebyshev, Newton_KL): rank-deficient design / Hessian → silent wrong answer or solver failure
- **Robust family** (IEPPA, IEPPA_SOFT): `RK_ERR_INFEAS` return — already handled

2-stage rescues affected methods by restricting fine calibration to well-populated cells.

---

## 3. Benefit Axes

Two ordered axes evaluated per method:

- **B (primary):** Robustness on sparse high-D interactions — does 2-stage prevent failure that single-stage cannot avoid?
- **A (secondary):** Compute speedup — does coarse K reduction shrink dominant cost (factorization, Hessian build)?

Cross-cutting concern: **Stage-2 disturbance** of Stage-1 margins. Mitigation cost varies by method (cheap for IPF/GREG, expensive for Newton/IPM).

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

## 5. Phased Scope

### Phase 1 — Biproportional family (3 methods)

**Targets:** RAKING, SINKHORN, GREENKHORN

**Rationale:** Largest robustness rescue. 2-stage eliminates the 0/0 / log(0) failure mode entirely. Methods are already cyclic coordinate descent → outer iteration to repair stage-2 disturbance is structurally natural and cheap. Simplest code path: margin loop already exists; 2-stage = outer wrapper with cell partition.

**Deliverable:** Coarse/fine split argument; cell-aggregation; n_min threshold; sparse-cell inheritance; per-cell residual fine cal.

### Phase 2 — Cubic-K parametric family (2 methods)

**Targets:** LOGIT, GREG

**Rationale:** Both bring autumn/anesrake parity (LOGIT *is* the native autumn pattern; GREG mirrors anesrake closed-form). Cubic compute win real (K=20 → K_coarse=4 ⇒ ~125× factorization shrink). GREG is one-shot (no outer-iter cost); LOGIT pays Newton re-factor per outer pass — accept as cost of correctness.

**Deliverable:** Same surface as P1 + per-method Hessian/design-matrix construction within cells.

### Out of scope (4 methods)

- **CHEBYSHEV_IPM:** Compute win exists but stage-2 disturbance fixup requires re-running IPM — prohibitively expensive. Defer until either (a) disturbance-free margin contract enforced or (b) explicit user demand.
- **NEWTON_KL:** TSVD already handles rank-deficiency; robustness gain too small. Compute win does not justify wiring cost given P1/P2 cover the use case.
- **IEPPA, IEPPA_SOFT:** Already robust. 2-stage gain reduces to "approximate fit on infeasible problems" — niche workflow. User can request opt-in extension if needed.

---

## 6. Stage-2 Disturbance Strategy

Stage-2 within-cell calibration changes weights, breaking Stage-1 coarse-margin sums. Two mitigations available; this spec selects **Strategy A** for P1 and **Strategy B** for P2.

### Strategy A — Outer iteration (P1 default)

Alternate Stage 1 ↔ Stage 2 until coarse-margin residual < `outer_tol` or `outer_max_iter` reached. Cheap for IPF (Stage 1 = ratio updates); composes with existing convergence framework.

### Strategy B — Orthogonality contract (P2 default)

Require user-declared coarse/fine split such that fine margins are within-cell-orthogonal to coarse partition (i.e., fine margins are exact splits of coarse cells). User contract — validated at preprocessing; non-orthogonal splits raise `RK_ERR_BADARG`. Eliminates disturbance by construction. autumn enforces this implicitly.

P2 methods cannot afford repeated Newton re-factorization, so contract is mandatory.

P1 supports both: Strategy A by default, Strategy B available via flag (skips outer loop).

---

## 7. API Surface (preview)

```r
calibrate(
  data, weights, margins,
  algorithm     = "raking",        # P1: raking|sinkhorn|greenkhorn; P2: logit|greg
  hierarchical  = list(
    coarse        = c("age", "region"),  # subset of margin names
    n_min         = 30,                  # cell threshold
    disturbance   = "iterate",           # "iterate" (P1) | "orthogonal" (P2)
    outer_tol     = 1e-4,                # only used when disturbance = "iterate"
    outer_max_iter = 10
  )
)
```

`hierarchical = NULL` (default) → single-stage behavior unchanged. ABI: extends `rk_params_t` via reserved-flag mechanism (not by adding fixed-arity SEXP — see `c-rk-calibrate-arity-is-fixed-at-36` memory).

---

## 8. Test Strategy

### Public-boundary tests (per phase)

- **Sparse-cell rescue:** K=9 DGP with empty cross-product cells → P1 methods converge with valid weights vs single-stage NaN/INFEAS
- **Compute scaling:** K=20 dense margins, K_coarse=4 → P2 methods show ≥ 50× wall-time win (factorization-bounded regime)
- **Stage-1 margin preservation:** disturbance-iterate path → coarse-margin residual < `outer_tol` at exit
- **Orthogonality contract:** non-orthogonal fine margins under `disturbance = "orthogonal"` → `RK_ERR_BADARG`
- **Single-stage parity:** `hierarchical = NULL` → bit-exact match to current solver output

### Adversarial fixtures

- Cell with n=1 (degenerate sparse)
- Cell with n_min - 1 (boundary)
- Coarse partition that produces all sparse cells (degenerate; should fall back to Stage 1 only)
- Disturbance-iterate that fails to converge in `outer_max_iter` → `RK_ERR_BUDGET` with diagnostic

---

## 9. Out of Scope (this design)

- IEPPA / IEPPA_SOFT 2-stage extension
- CHEBYSHEV / NEWTON_KL 2-stage extension
- Auto-detection of optimal coarse/fine split (user-specified only)
- Three-or-more-stage extensions
- Adaptive `n_min` per cell
- Bounds_mode = "unit" interaction with cell partition (treat as default-unsupported; raise `RK_ERR_BADARG` if combined)

---

## 10. Open Questions

1. **n_min default:** 30 inherited from autumn; validate against leafblower's typical N=10⁴–10⁶ regimes — may want per-K_fine scaling (e.g., n_min = max(30, 5·K_fine))
2. **Inheritance semantics:** sparse cells inherit Stage-1 weights *exactly*, or inherit-and-rescale-to-cell-target? Latter preserves cell totals but reintroduces disturbance.
3. **Bounds_mode interaction:** default `"cell"` already operates at cell aggregate — does 2-stage need cell-of-cells reinterpretation, or is the existing aggregation sufficient?

These resolve at planning phase, not blocking design approval.

---

## 11. Next Step

After user review and design-review-gate (5-agent: PM, Architect, Designer, Security, CTO) approval, hand off to `writing-plans` for P1 implementation plan. P2 plan deferred until P1 ships and proves the API surface.
