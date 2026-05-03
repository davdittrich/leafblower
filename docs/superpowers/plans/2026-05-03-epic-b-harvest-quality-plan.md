# Epic B: harvest.R quality + performance — Implementation Plan

**Date:** 2026-05-03
**Beads epic:** `leafblower-fmot`
**Beads tasks:** `fmot.2` (B1), `fmot.1` (B2)
**Motivation:** Code review 2026-05-03 — SIMPLIFY + SUGGESTIONS in R/harvest.R.

## Mechanism

Two sequential harvest.R refactors. No algorithm changes. Metric values numerically identical before/after.

**Mechanism per WU:**
- B1: Extract 28-line metric block into `compute_quality_metrics(weights, target, data)` helper (in `# --- Helpers ---` section). Table-drive `convergence_reason` 8-branch if/else via named-vector lookup keyed on `(status, accelerate_bool)`.
- B2: Replace K separate `tapply` calls in margin_kl with single-pass aggregation. Dispatch between single-pass and K-pass based on K and NA-presence.

**Forbidden:** Metric value changes (must be identical to within 1e-10 relative). Algorithm changes. Any C++ changes.

**Audit:**
- B1: `all.equal(compute_quality_metrics_pre_extract(w, tgt, df), compute_quality_metrics(w, tgt, df))` → TRUE on stepstone.
- B2: `all.equal(margin_kl_old, margin_kl_new)` → TRUE; `microbenchmark` shows improvement on stepstone K=9 (target: ≥ 2× speedup on K ≥ 5).

## Work Units

| WU | Bead | Title | Dep | Model | Wall |
|---|---|---|---|---|---|
| B1 | `leafblower-fmot.2` | Extract compute_quality_metrics + table-drive convergence_reason | — | Haiku | ~45min |
| B2 | `leafblower-fmot.1` | Single-pass margin_kl aggregation | B1 | Gemini | ~1h |

**Dep chain:** B1 → B2 (linear).

## Decision Rule

Epic B closes PASS iff:
- All 2 WU tickets closed.
- `devtools::test()` FAIL=2 only (pre-existing T2 basin floor).
- Metric values identical (all.equal TRUE).
- B2: documented speedup ≥ 2× on K ≥ 5 fixture.

## Risks

| # | Risk | Mitigation |
|---|---|---|
| R1 | B2 single-pass changes NA handling subtly | B2 step 3: test on fixture with NA-coded observations; verify Inf/NA behavior preserved |
| R2 | compute_quality_metrics scoping: `data` variable collides with R base `data()` | Use `df` parameter name inside helper; avoid `data` as formal arg |
| R3 | B2 single-pass slower on K=1 or K=2 | Dispatch: K < 3 → K-pass (existing); K ≥ 3 → single-pass |
