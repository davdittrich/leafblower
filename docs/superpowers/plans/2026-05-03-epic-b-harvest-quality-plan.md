# Epic B: harvest.R quality + performance — Implementation Plan (rev 2)

**Date:** 2026-05-03
**Beads epic:** `leafblower-fmot`
**Beads tasks:** `fmot.2` (B1), `fmot.1` (B2)
**Motivation:** Code review 2026-05-03 — SIMPLIFY + SUGGESTIONS in R/harvest.R.
**Revision:** rev 2 closes 3 plan-review-gate iter 1 Feasibility blockers (B2 NA semantics non-equivalence; speedup unfounded; rowsum NA warning).

## Mechanism

### B1 — extract compute_quality_metrics + table-drive convergence_reason

**Extraction:** Move R/harvest.R:536-575 (fixed metric block from 51d3ed3) to a new helper `compute_quality_metrics <- function(weights, target_list, df)` in the `# --- Helpers ---` section. Parameters named `target_list` (not `target` — avoids shadowing) and `df` (not `data` — avoids collision with base `data()`). Caller in harvest.R body: `qm <- compute_quality_metrics(weights, target, data); calib_result[names(qm)] <- qm`.

**convergence_reason:** Current code (harvest.R:438-447) has a 6-branch if/else + else = 7 outcomes on `(s, accelerate_bool)`. Table-drive via a 2-column lookup (status int → reason string), with the `s==5L && accelerate_bool` distinction handled by a pre-check before the lookup. Simpler refactor: use `switch(as.character(s), "0"="criterion", "2"="infeasible", "3"="error", "4"="budget", "5"=..., "1"="legacy", "legacy")` with the accelerate_bool modifier applied after. Block shrinks from 12 lines to ~6.

**Audit:** Capture `res_before <- compute_quality_metrics(w, tgt, df)` with old inline code; compare `all.equal(res_before, compute_quality_metrics(w, tgt, df))` after extraction → TRUE.

### B2 — single-pass margin_kl (CONDITIONAL, dispatch-guarded)

**Root state:** B1's `compute_quality_metrics` calls per-margin `tapply` with NA-exclusion per margin. `K` separate hash-builds each O(n).

**B2 scope revision (rev 2):**

The `interaction()`-based approach is BLOCKED due to two issues:
1. **NA non-equivalence**: `interaction(col_a, col_b, ...)` where row i has NA in col_a produces `<NA>` level — row i is excluded from ALL K margins, not just margin a. This differs from per-margin NA exclusion.
2. **Combinatorial explosion**: for K=9 margins with L levels each, `interaction()` produces `prod(L_k)` factor levels. For stepstone K=9 × ~10 levels = 10^9 potential cells — O(n log n) sort vs K×O(n) hash is not obviously faster.

**Revised B2 mechanism (conservative approach):**
- Profile first: `microbenchmark::microbenchmark(compute_quality_metrics(w, tgt, df_sub), times=20)` on stepstone.
- If any column has NAs (`anyNA(df[names(target)])`): use per-margin K-pass (existing) — semantics correct.
- If no NAs AND K ≥ 5: try `rowsum(weights, do.call(paste, c(df[names(target)], sep="\x01")))` approach — single string-concatenation hash instead of K factor hashes. Reduces K hash builds to 1 string pass. Profile: speedup must exceed 1.5×.
- Suppress `rowsum` NA-group warning via `withCallingHandlers(..., warning = function(w) if (grepl("missing values", conditionMessage(w))) invokeRestart("muffleWarning") else w)` (not blanket `suppressWarnings`).
- **STOP rule**: if measured speedup < 1.5× on stepstone, revert B2 entirely and close with DEFERRED note. Do NOT force.

**Audit:** `all.equal(margin_kl_old, margin_kl_new, tolerance = 1e-10)` on stepstone (no NAs) AND on a synthetic fixture with 5% NA observations (must fall back to K-pass, yielding identical output to K-pass baseline).

## Work Units

| WU | Bead | Title | Dep | Model | Wall |
|---|---|---|---|---|---|
| B1 | `leafblower-fmot.2` | Extract compute_quality_metrics helper + table-drive convergence_reason | — | Haiku | ~45min |
| B2 | `leafblower-fmot.1` | Single-pass margin_kl (conditional, with STOP rule) | B1 | Gemini | ~1h |

**Dep chain:** B1 → B2 (linear).

## Decision Rule

Epic B closes PASS iff:
- Both WU tickets closed (or B2 closed DEFERRED if speedup < 1.5×).
- `devtools::test()` FAIL=2 only (pre-existing T2 basin floor).
- Metric values identical on stepstone AND NA fixture.
- B1: helper present in Helpers section; convergence_reason simplified.

## Risks

| # | Risk | Mitigation |
|---|---|---|
| R1 | B2 NA non-equivalence | dispatch: `if (anyNA(df[names(target)])) k_pass() else single_pass()`; NA fixture test required in DoD |
| R2 | B2 speedup not realized | STOP rule: profile first, revert if <1.5× |
| R3 | B2 rowsum NA warning | Targeted muffleWarning (not suppressWarnings) on "missing values" pattern |
| R4 | B1 parameter name collision (`data` / `target` shadow base namespace) | Use `df` and `target_list` as formal params; verified in code review |
| R5 | convergence_reason switch misses future status codes | Document that switch has an else clause that falls back to "legacy" |
