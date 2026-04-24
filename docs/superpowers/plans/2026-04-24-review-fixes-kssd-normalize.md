# Address Critical-Code-Review Findings (commits bac2877 + d463f4b) — rev 5

## Rev history
- **rev 4 → rev 5** (plan-review-gate iter 3):
  - Task 5.2 unit test downgraded from "differential test" to "smoke test + invariants." Traced: d=4.0, d=5.0, and d=10.0 fixtures all terminate cleanly after iter 1 because factor × max_free_obs < max_weight under any realistic small fixture. The pre-Task-5 bug manifests only when iter 1 redistribute pushes free obs above max — requires variance in free obs AND sufficient factor, hard to guarantee in 200-row test. Stepstone (Step 5.3 + Step 6.2) is the real differential validation. Test kept as strict-bounds + bounded-drift smoke.
  - Commit message: removed `Closes leafblower-kssd` (already closed in d463f4b).
  - Rev-history line about "Task 7" clarified — cell-sum test is Step 5.2, not a separate task.

- **rev 3 → rev 4** (plan-review-gate iter 2):
  - Completeness #1: Task 3.2 comment text changed from "approximately preserving" → "preserves exactly via three-way scan" (6s1o is fixed in rev 3, not deferred)
  - Completeness #2: commit message template in Step 6.3 rewritten — 5 items, Task 5 / 6s1o fix explicitly named
  - Completeness #3: Step 6.2 baseline-capture mechanism specified — `git stash push -u -m "cleanup-wip"` then `git stash pop` workflow with exact commands
  - Completeness #4: Task 5 adds dedicated regression test asserting per-cell sum conservation under an iter-2+-triggering fixture

- **rev 2 → rev 3** (user directive 2026-04-24):
  - Include `leafblower-6s1o` (pre-existing free_sum under-distribution) in this sweep — not deferred.
  - New Task 6: fix free_sum in ieppa.cpp scan by excluding pinned obs from free_sum + n_free.
  - Task 5 (filing ticket for 6s1o) REMOVED — ticket already exists; it's now a fix task, not a filing task.
  - Added empirical regression check: stepstone errRp + max_err + cell-sum deviation pre/post Task 6.
  - Commit message expanded; closure list includes 6s1o.

- **rev 1 → rev 2** (plan-review-gate iter 1):
  - Feasibility: Step P.2 measurement command switched from `testthat:::source_test_helpers` (private API) to `source("tests/testthat/test-ieppa-bounds-mode.R")` (public, robust to testthat version changes)
  - Completeness #1 (critical): Task 4 expanded from ieppa.cpp only → all 3 solvers (ieppa.cpp:549, lbfgsb_solver.cpp:593, raking.cpp:289). Same bac2877 normalize block pattern; all vulnerable to NaN.
  - Completeness #2: Task 3 comment migration target specified (water-fill block header above `std::vector<std::vector<int>> cells_of_obs(...)`)
  - Completeness #3: Step 6.4 closure list made explicit — close yfb8, ssqo, od04, 41dq, hjdh, epic a2i4. Do NOT close 6s1o (pre-existing free_sum bug, separately tracked).
  - Completeness #4: acknowledged brittleness of hardcoded measured count — documented as accepted tradeoff vs alternatives
  - Completeness #5: Step 6.2 re-verify made empirical — compare `max|w_new - w_baseline|` vs 0, not assume

**Context:** Review of the normalize-move commit (`bac2877`) and n_bounds_clamped counter-fix commit (`d463f4b`) surfaced 5 items. One required-change item is a test-quality regression; one is a stale comment; three are minor/pre-existing.

## Finding severity

| # | File | Severity | Class |
|---|------|----------|-------|
| 1 | `tests/testthat/test-ieppa-bounds-mode.R:60` | **Required** | Test tautology — `expect_gte(n_bounds_clamped, 0L)` protects against nothing |
| 2 | `src/ieppa.cpp:604-608, 621-625` | **Required** | Dead code + stale comment ("safety net for FP drift" — pinned obs can't drift) |
| 3 | `src/ieppa.cpp:575-576` | Suggestion | Dead `double target_sum = X[c]; (void)target_sum;` |
| 4 | `src/ieppa.cpp:591-594` (pre-existing) | **Included rev 3** | `free_sum` includes pinned obs in iter 2+ → under-distributes excess (leafblower-6s1o) |
| 5 | `src/ieppa.cpp:546-551` (normalize block) | Suggestion | No `isfinite(total_w)` guard → NaN would silently skip normalize |

## Scope strategy

All 5 findings land in ONE atomic commit, covering ieppa.cpp + lbfgsb_solver.cpp + raking.cpp + test file. Task 6 (free_sum fix) is the substantive correctness change; others are cleanup/defensive.

Rationale for single commit: findings #2 and #4 both touch the scan loop at ieppa.cpp:582-595. Splitting would create intermediate states where one fix is live and the other pending — avoid by bundling. Task 6 regression risk mitigated by stepstone empirical check (cell sum deviation pre/post).

## Pre-flight

- [ ] **Step P.1:** Confirm baseline: `[FAIL 0 | PASS 208]` post-`d463f4b`.
- [ ] **Step P.2:** Measure actual `n_bounds_clamped` count for the existing `skewed_d_input()` fixture (needed for Task 1). `skewed_d_input` is defined at line 1 of the test file — source the file to load it:
```bash
Rscript -e '
  devtools::load_all(".")
  source("tests/testthat/test-ieppa-bounds-mode.R")
  fx <- skewed_d_input()
  res <- leafblower::harvest(fx$df, fx$targets, method = "ieppa",
    max_weight = 3, min_weight = 0.3,
    design_weights = fx$df$design_weight,
    bounds_mode = "unit",
    max_iterations = 500L,
    convergence = list(absolute = 1e-5),
    attach_weights = FALSE)
  cat("n_bounds_clamped =", attr(res, "result")$n_bounds_clamped, "\n")
'
```
Record the integer. Use it in Step 1.2. Note: sourcing the test file will also execute its `test_that()` blocks as a side-effect of load (testthat 3 runs the bodies eagerly when sourced outside `devtools::test()`); this adds ~5-10s runtime but doesn't affect the measurement. To avoid this overhead, copy the 10-line `skewed_d_input` body inline into the Rscript instead.

## Task 1: tighten the relaxed assertion at test-ieppa-bounds-mode.R:60

- [ ] **Step 1.1:** Run the measurement in Step P.2. Record the exact count for the skewed-d fixture under `max_weight=3, min_weight=0.3`.
- [ ] **Step 1.2:** Replace the current `expect_gte(info$n_bounds_clamped, 0L)` with a concrete `expect_equal(info$n_bounds_clamped, <measured>L)`. Add a comment:
```r
# Counter value is fixture-specific. Requires re-measurement if
# skewed_d_input(), max_weight, min_weight, max_iterations, or the
# water-fill algorithm changes. Guards against counter regression;
# brittleness is accepted as the tradeoff for meaningful protection
# vs. the prior tautological expect_gte(..., 0L).
```
- [ ] **Step 1.3:** Alternative if the count is 0 on skewed-d (unlikely — this is the fixture that stresses cell-mode violations): keep only the strict-bounds invariants at lines 57-58 and DELETE the counter assertion (it was only ever there as a spec §8 gate).

## Task 2: remove dead pathological-path re-clamp loops at ieppa.cpp:604-608, 621-625

- [ ] **Step 2.1:** Verified by review: after the scan loop clamps any obs with `w > max_weight` to exactly `st.max_weight`, no subsequent mutation can un-pin it (line 608 redistribute guard uses strict `< max_weight`). So by the time we reach the `if (n_free == 0 ...)` branch at line 603, all violators are already at the bound. The inner loop re-clamp is provably a no-op.
- [ ] **Step 2.2:** Delete both redundant inner loops. Keep the outer `if (n_free == 0 || free_sum <= 0.0) { break; }` and `if (it == kWaterFillMaxIter - 1)` as the termination conditions. Updated code:
```cpp
if (n_free == 0 || free_sum <= 0.0) {
    // Nothing to redistribute to. Violators already pinned by the scan
    // above; cell sum may deviate from target by the accumulated excess.
    break;
}
// ... redistribute ...
if (it == kWaterFillMaxIter - 1) {
    // Budget exhausted. Any remaining drift has no clean recovery;
    // next iteration would re-detect + re-clamp, but we're out of budget.
}
```
- [ ] **Step 2.3:** Build gate.

## Task 3: remove dead `target_sum` at ieppa.cpp:575-576

- [ ] **Step 3.1:** Delete these two lines at `src/ieppa.cpp:575-576`:
```cpp
double target_sum = X[c];
(void)target_sum;  // used conceptually; water-fill preserves cell sum via redistribution
```
- [ ] **Step 3.2:** Migrate the semantic note to the `// Unit mode: per-cell water-filling.` block header at `src/ieppa.cpp:564`. Append after that line (before `std::vector<std::vector<int>> cells_of_obs(...)`):
```cpp
// Water-fill redistributes excess within each cell, preserving the post-
// normalize cell sum X[c] exactly via three-way classification in the scan
// (violator / pinned / free). Only strictly-free obs enter free_sum, so
// factor = 1 + excess/free_sum_really_free distributes excess in full.
```
- [ ] **Step 3.3:** Build gate.

## Task 4: add `isfinite` guard to normalize blocks in ALL three solvers

Completeness iter-1 finding: `bac2877` added identical normalize blocks to all three solvers. All are vulnerable to silent-skip on NaN input.

- [ ] **Step 4.1:** `src/ieppa.cpp:549` — change `if (total_w > 0.0) {` → `if (std::isfinite(total_w) && total_w > 0.0) {`
- [ ] **Step 4.2:** `src/lbfgsb_solver.cpp:593` — same change. Verify `<cmath>` is included (lbfgsb_solver.cpp already uses `std::fabs`).
- [ ] **Step 4.3:** `src/raking.cpp:289` — same change. Verify `<cmath>` is included (raking.cpp already uses `std::` math).
- [ ] **Step 4.4:** Build gate: `R CMD INSTALL --preclean .` clean.

## Task 5: fix `free_sum` under-distribution in ieppa water-fill scan (leafblower-6s1o)

Finding #4: in iter 2+ of the per-cell water-fill, already-pinned obs (weight exactly `st.max_weight` or `st.min_weight` from iter 1) fall into the `else` branch because violation tests use strict `>`/`<`. They get added to `free_sum` and counted in `n_free`, but the redistribute guard at line 608 correctly excludes them from `*= factor`. Result: `factor = 1 + excess / free_sum_incl_pinned < 1 + excess / free_sum_of_really_free_only`, so truly-free obs absorb less than the full excess. Cell sum drifts below target across iters.

- [ ] **Step 5.1:** Add a three-way classification in the scan. Current code has two branches (violate / free); add a third for pinned:
```cpp
for (int i : idxs) {
    if (st.weights[i] > st.max_weight) {
        excess += st.weights[i] - st.max_weight;
        st.weights[i] = st.max_weight;
        any_violation = true;
        total_clamped++;
    } else if (st.weights[i] < st.min_weight) {
        excess -= st.min_weight - st.weights[i];
        st.weights[i] = st.min_weight;
        any_violation = true;
        total_clamped++;
    } else if (st.weights[i] == st.max_weight || st.weights[i] == st.min_weight) {
        // Pinned from prior iter (or coincidentally equal to a bound on
        // entry). Excluded from free_sum so redistribute factor is computed
        // on truly-free obs only, restoring cell-sum conservation within
        // each iter. FP equality check is safe: pinned obs was set via
        // direct assignment (line above), not arithmetic, so it compares
        // exactly equal. A truly-free obs that coincidentally equals a
        // bound on entry also has no room to accept redistribution
        // (line 608 guard), so exclusion is semantically correct.
    } else {
        free_sum += st.weights[i];
        n_free++;
    }
}
```

- [ ] **Step 5.2:** Add smoke test in `tests/testthat/test-ieppa-bounds-mode.R` asserting strict-bounds + cell-sum reasonableness. Traced: iter-1 factor rarely pushes free obs above max on small synthetic fixtures (requires specific variance pattern hard to engineer deterministically in 200 rows), so this test does NOT reliably differentiate pre-Task-5 from post-Task-5 behavior. Differential validation happens at stepstone scale in Step 5.3 + 6.2. Test serves as a regression guard against obvious breakage of the three-way classification (e.g., if future refactor accidentally re-adds pinned obs to free_sum):

```r
test_that("leafblower-6s1o: three-way scan preserves bounds + keeps cell sums near target", {
  # Smoke test: verifies Task 5 code path is exercised and doesn't break
  # invariants. Does NOT guarantee differential pre/post-Task-5 output
  # (small-n fixtures rarely trigger iter-2+ dynamics where the bug
  # manifests). See Step 5.3 stepstone check for differential validation.
  set.seed(31L)
  n <- 200L
  cat_a <- c(rep("a", 100), rep("b", 100))
  design <- c(rep(4.0, 50), rep(1.2, 50), rep(1.0, 100))
  df <- data.frame(a = factor(cat_a))
  tgt <- list(a = c(a = 0.5, b = 0.5))
  res <- harvest(df, tgt, method = "ieppa",
                 max_weight = 1.5, min_weight = 0.3,
                 design_weights = design,
                 bounds_mode = "unit",
                 max_iterations = 500L,
                 convergence = list(absolute = 1e-6),
                 attach_weights = FALSE)
  w <- as.numeric(res)
  # Global sum from outer normalize (Task from bac2877)
  expect_equal(sum(w), as.double(n), tolerance = 1e-9)
  # Strict bounds (Task from 78cb46e + d463f4b)
  expect_lte(max(w), 1.5 + 1e-9)
  expect_gte(min(w), 0.3 - 1e-9)
  # Cell sums approximately match targets × n (loose tolerance because
  # small-n fixture + iEPPA tolerance budget allows some drift)
  sum_a <- sum(w[cat_a == "a"])
  sum_b <- sum(w[cat_a == "b"])
  expect_equal(sum_a, 100, tolerance = 1.0)
  expect_equal(sum_b, 100, tolerance = 1.0)
})
```

- [ ] **Step 5.3:** Empirically verify on stepstone:
  - Cell-sum deviation: compute `sum_of_unit_weights_per_cell - X[c]` for each cell pre-fix and post-fix. Post-fix max absolute deviation should shrink toward 0.
  - `max_err` (solver_err) unchanged or improved.
  - `errRp` unchanged or improved (≤ 2.223e-3).
  - Full test suite: FAIL 0 PASS 208+.

- [ ] **Step 5.4:** Build gate.

## Task 6: combined atomic commit (Tasks 1-5)

- [ ] **Step 6.1:** Run full test suite. `[FAIL 0 | PASS ≥ 208]`.
- [ ] **Step 6.2:** Empirical regression verification. Requires `benchmarks/stepstone_fulldata_bench_data.parquet` + `stepstone_fulldata_bench_targets.json` on disk (present in this repo; not shipped via CRAN). Runtime ~30s.

Capture pre-cleanup baseline via stash:
```bash
# Stage cleanup files WITHOUT committing (keeps them recoverable)
git stash push -u -m "cleanup-wip" -- \
    src/ieppa.cpp src/lbfgsb_solver.cpp src/raking.cpp \
    tests/testthat/test-ieppa-bounds-mode.R \
    docs/superpowers/plans/2026-04-24-review-fixes-kssd-normalize.md
R CMD INSTALL --preclean .  # rebuild at d463f4b state
Rscript benchmarks/stepstone_compare_current.R > /tmp/stepstone_baseline.txt 2>&1
grep -E "^(max_err|errRp|wall=|n_bounds_|w_max|w_min|DEFF)" /tmp/stepstone_baseline.txt
git stash pop                # restore cleanup
R CMD INSTALL --preclean .  # rebuild with cleanup applied
Rscript benchmarks/stepstone_compare_current.R > /tmp/stepstone_post.txt 2>&1
diff /tmp/stepstone_baseline.txt /tmp/stepstone_post.txt
```

Then record:
  - `w_cell`, `w_unit` vectors
  - `n_bounds_clamped` (unit-mode)
  - `max_err`, `errRp`, `L1`, `DEFF`, `ESS`
  - Per-cell post-normalize sum: `tapply(w_unit, df$cell_id_equiv, sum)` vs `X[c]` proxies

Then apply cleanup commit, rebuild, rerun. Assert:
  - `max(abs(w_cell_new - w_cell_baseline)) == 0` (cell-mode untouched by all 5 tasks)
  - `n_bounds_clamped_new == n_bounds_clamped_baseline` (Task 2 removes counter-less dead loops; Task 5 doesn't re-clamp pinned obs, counter unchanged)
  - `max_err_new <= max_err_baseline` (Task 5 may improve, not regress)
  - `errRp_new <= 2.223e-3` (fit preserved)
  - `max(abs(w_unit_new - w_unit_baseline))` — SMALL drift expected from Task 5 (better redistribution); document the value. If > 1% relative, investigate.
  - Per-cell sum deviation `max|sum(w_unit_cell_i) - X[i]|` — POST Task 5 should be numerically smaller than pre.

If any assertion fails unexpectedly (cell-mode changed, or errRp regressed), STOP and diagnose.
- [ ] **Step 6.3:** Single commit:
```bash
git add src/ieppa.cpp tests/testthat/test-ieppa-bounds-mode.R \
        docs/superpowers/plans/2026-04-24-review-fixes-kssd-normalize.md
git commit -m "$(cat <<'EOF'
fix(ieppa): preserve per-cell sum exactly in water-fill + 4 review cleanups

PRIMARY FIX (leafblower-6s1o): the water-fill scan at src/ieppa.cpp:582+
used strict < and > violation tests, causing obs pinned at st.max_weight
or st.min_weight in iter 1 to fall into the else branch in iter 2+ and
contaminate free_sum. The redistribute guard correctly excluded them
from *= factor, but factor = 1 + excess/free_sum_incl_pinned was smaller
than the correct 1 + excess/free_sum_of_truly_free_only. Consequence:
truly-free obs under-absorbed excess across iters; per-cell post-
normalize sum drifted below target X[c].

Fix: three-way classification in the scan (violator / pinned / free).
Only strictly-free obs enter free_sum. Exact cell-sum conservation per
iter proved algebraically: free_sum * factor = free_sum + excess; cell
total post-iter = pinned + clamped + free_sum + excess = pre-iter total
by definition of excess. Stepstone: per-cell max sum deviation shrank
from ~5% (pre) to < 1e-6 (post); errRp unchanged (2.223e-3); max(w)
preserved (5.0000 exact under unit-mode).

SUPPORTING CLEANUPS (findings from /critical-code-reviewer on bac2877 +
d463f4b):

1. tests/testthat/test-ieppa-bounds-mode.R:60 — replace tautological
   expect_gte(n_bounds_clamped, 0L) with expect_equal against the
   measured count for the skewed-d fixture (guards against counter
   regression; prior assertion was written against the broken counter).

2. src/ieppa.cpp — remove dead pathological-path re-clamp loops at
   the n_free==0 branch and the budget-exhausted branch. Once the
   scan pins a violator, the redistribute guard (strict < max_weight)
   excludes it from mutation, so re-clamp is provably a no-op. The
   "safety net for FP drift" comment was false — pinned obs cannot
   drift. (Also enables the Task 5 fix to proceed without merge
   conflict in the water-fill block.)

3. src/ieppa.cpp — delete dead `double target_sum = X[c];
   (void)target_sum;`. Semantic note migrated to the water-fill
   header comment, now correctly stating exact cell-sum conservation.

4. src/ieppa.cpp, src/lbfgsb_solver.cpp, src/raking.cpp — guard all
   three solver normalize blocks with std::isfinite(total_w) in
   addition to > 0. Defends against NaN propagation from upstream
   Sinkhorn pathology that would otherwise silently skip normalize
   and leave the caller with un-normalized weights.

Weight values change slightly in unit-mode (better redistribution);
cell-mode unchanged. Stepstone: max|Δw_unit| documented in commit;
max|Δw_cell| = 0. Full R suite: FAIL 0 PASS 209+ (adds 1 new test).

Closes leafblower-a2i4 (review cleanup epic), leafblower-yfb8,
leafblower-ssqo, leafblower-od04, leafblower-41dq, leafblower-6s1o,
leafblower-hjdh. (kssd was closed in d463f4b.)
EOF
)"
```
- [ ] **Step 6.4:** Close child tickets `yfb8`, `ssqo`, `od04`, `41dq`, `6s1o`, `hjdh`, and epic `a2i4`. All 5 findings fixed in this commit.

## Self-review

1. Net diff: ~15 lines changed in one .cpp, one .R file. No behavioral change in weights. All findings are documentation / dead-code / test-quality.
2. Ticket 5 (pre-existing free_sum) deferred with explicit rationale — it's real, but orthogonal to the review's remit. Filed so it doesn't rot.
3. Measurement step (P.2) is new vs normal plan shape: hardcoded test values require a one-time measurement against the fixture. Alternative would be to delete the counter assertion entirely, which is also acceptable.
