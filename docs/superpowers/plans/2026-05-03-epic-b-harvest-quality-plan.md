# Epic B: harvest.R quality + performance — Implementation Plan (rev 2)

**Date:** 2026-05-03
**Beads epic:** `leafblower-fmot`
**Beads tasks:** `fmot.2` (B1), `fmot.1` (B2)
**Motivation:** Code review 2026-05-03 — SIMPLIFY + SUGGESTIONS in R/harvest.R.
**Revision:** rev 2 closes 3 plan-review-gate iter 1 Feasibility blockers (B2 NA semantics non-equivalence; speedup unfounded; rowsum NA warning).

## Mechanism

### B1 — extract compute_quality_metrics + table-drive convergence_reason

**Extraction:** Move R/harvest.R:536-575 (fixed metric block from 51d3ed3) to a new helper `compute_quality_metrics <- function(weights, target_list, df)` in the `# --- Helpers ---` section.

**Signature locked:** formal parameters MUST be named `target_list` (not `target` — avoids shadowing the parent-scope `target` variable) and `df` (not `data` — avoids collision with base R `data()` function). This is a DoD requirement.

**Pre-extraction snapshot (required before B1 commit):** Extend existing `tests/testthat/test-quality-metrics.R` (already exists — DO NOT create new file) with a new `test_that("compute_quality_metrics extraction: values identical to inline", { ... })` block that captures expected metric values from `attr(r, "result")` on the T1 fixture. This test must PASS before the extraction commit AND after. Commit the test extension in the SAME commit as the extraction.

Caller in harvest.R body: `qm <- compute_quality_metrics(weights, target, data); calib_result[names(qm)] <- qm`.

**convergence_reason:** Current code (harvest.R:438-447) has a 6-branch if/else + else = 7 outcomes on `(s, accelerate_bool)`. Table-drive via a 2-column lookup (status int → reason string), with the `s==5L && accelerate_bool` distinction handled by a pre-check before the lookup. Simpler refactor: use `switch(as.character(s), "0"="criterion", "2"="infeasible", "3"="error", "4"="budget", "5"=..., "1"="legacy", "legacy")` with the accelerate_bool modifier applied after. Block shrinks from 12 lines to ~6.

**Audit:** Capture `res_before <- compute_quality_metrics(w, tgt, df)` with old inline code; compare `all.equal(res_before, compute_quality_metrics(w, tgt, df))` after extraction → TRUE.

### B2 — single-pass margin_kl (CONDITIONAL, dispatch-guarded)

**Root state:** B1's `compute_quality_metrics` calls per-margin `tapply` with NA-exclusion per margin. `K` separate hash-builds each O(n).

**B2 scope revision (rev 2):**

The `interaction()`-based approach is BLOCKED due to two issues:
1. **NA non-equivalence**: `interaction(col_a, col_b, ...)` where row i has NA in col_a produces `<NA>` level — row i is excluded from ALL K margins, not just margin a. This differs from per-margin NA exclusion.
2. **Combinatorial explosion**: for K=9 margins with L levels each, `interaction()` produces `prod(L_k)` factor levels. For stepstone K=9 × ~10 levels = 10^9 potential cells — O(n log n) sort vs K×O(n) hash is not obviously faster.

**Revised B2 mechanism (conservative approach):**
**Steps (locked contract):**

1. **Profile baseline first** (no code changes): `Rscript -e 'library(microbenchmark); library(leafblower); source("benchmarks/profile_margin_kl.R"); print(mb)' 2>&1 | tee benchmarks/results/margin_kl_profile_baseline.txt`. Commit `benchmarks/profile_margin_kl.R` script and `benchmarks/results/margin_kl_profile_baseline.txt`. This output is the STOP RULE gate.

2. **Dispatch guard**: `anyNA(df[names(target)])` checks ALL target columns for any NA. When TRUE → K-pass (semantically correct, per per-margin NA exclusion). When FALSE AND K ≥ 5 → single-pass.

3. **Single-pass via paste+rowsum with per-margin projection**: 
   ```r
   cell_key <- do.call(paste, c(df[names(target)], list(sep = "\x01")))
   W_cell   <- rowsum(weights, cell_key, na.action = NULL)  # named vector: cell → w_sum
   Z        <- sum(weights)
   # Project back to per-margin W_kj: for each margin k, sum W_cell entries
   # whose cell_key has obs_k == level_j. Efficient via split on margin k alone.
   margin_kl_single <- sum(sapply(names(target), function(k) {
     T_k   <- target[[k]]
     # Aggregate: for each level j of margin k, sum W_cell over all cells with obs[k]==j
     key_k_only <- df[[k]]   # factor; maps obs → level of margin k
     W_k   <- tapply(as.numeric(W_cell[cell_key]), key_k_only, sum, na.rm = TRUE) / Z
     # (tapply on cell_key subset — still 1 hash over cells not obs)
     # NB: this is NOT 1 hash per margin; it's 1 hash per margin-level group.
     # TRUE single-pass would use pre-aggregated W_cell to avoid per-margin tapply.
     # Implement as: W_k <- vapply(split(W_cell, tapply(seq_along(cell_key), key_k_only, identity)), sum, 0)
     common <- intersect(names(T_k), names(W_k))
     if (any(T_k[setdiff(names(T_k), names(W_k))] > 0)) return(Inf)
     pos <- T_k[common] > 0
     sum(T_k[common][pos] * log(T_k[common][pos] / pmax(W_k[common][pos], 1e-15)))
   }))
   ```
   The full projection algorithm is: (1) cell-aggregate weights via paste-rowsum once; (2) per margin k, group cell-level weights by margin-k level to get W_kj. Step 2 is K tapply-over-M_cell (not K tapply-over-n) — the speedup is O(M_cell × K) vs O(n × K) where M_cell << n (but M_cell ≈ n on zero-compression fixtures). Speedup only materializes when M_cell/n < 0.5; profile gate below. NA keys treated as their own level — safe because dispatch ensures no NAs reach this path. Suppress the rowsum "missing values for group" warning via targeted handler: `withCallingHandlers(rowsum(...), warning = function(w) { if (grepl("missing values", conditionMessage(w))) invokeRestart("muffleWarning"); w })`.

4. **STOP RULE**: After single-pass implementation, profile again: `benchmarks/results/margin_kl_profile_single_pass.txt`. If median speedup on stepstone < 1.5×, REVERT implementation, close B2 DEFERRED with profile evidence, do NOT commit single-pass code.

5. **Dispatch test (DoD requirement)**: `test_that("B2 dispatch: NA data uses K-pass", { ... stopifnot(!environment_was_single_pass) ... })` — add to `tests/testthat/test-metrics.R`. Verify single-pass branch NOT taken when data has NAs (use mockery or spy on environment variable set in single-pass branch).

6. **muffleWarning test (DoD requirement)**: `test_that("B2 muffleWarning: suppresses rowsum warning; passes other warnings", ...)` — verifies targeted handler fires only on "missing values" text.

**NA fixture spec (DoD requirement):** Synthetic fixture in B2 test — `set.seed(42L); n <- 500L; K <- 3L; df_na <- data.frame(a=factor(sample(c("x","y","z"), n, TRUE)), b=factor(sample(c("M","F"), n, TRUE)), c=factor(sample(c("1","2","3"), n, TRUE))); df_na$a[sample(n, 25L)] <- NA; tgt_na <- list(a=c(x=0.5,y=0.3,z=0.2), b=c(M=0.6,F=0.4), c=c("1"=0.4,"2"=0.4,"3"=0.2))`. Verify K-pass output matches dispatch output byte-for-byte.

**Audit:** `all.equal(margin_kl_old, margin_kl_new, tolerance = 1e-10)` on stepstone (no NAs) + NA fixture (must fall back → identical to K-pass baseline).

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
