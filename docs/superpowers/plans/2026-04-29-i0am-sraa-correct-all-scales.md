# SRAA-m Correct Acceleration for All K Scales

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix SRAA-m giving 35% worse quality than plain on K>=3 overlapping-margin problems by making f_eval_sraa's sort order adaptive (Track 2), with a plateau-gated outer revert as a safety net (Track 1) if needed.

**Architecture:** Two tracks. Track 2: 5-line change in greenkhorn.cpp f_eval_sraa making it use adaptive sort (same as plain greenkhorn) — gives AA a unique max_err-optimal fixed point. Track 1: plateau-gated AA activation + outer revert in both greenkhorn.cpp and raking.cpp outer loops, with aa_unlocked lifecycle gate in SRAAState.

**Tech Stack:** C++17 (greenkhorn.cpp, raking.cpp, sraa.hpp), R testthat3.

---

## Spec Reference

`docs/superpowers/specs/2026-04-29-i0am-sraa-correct-all-scales.md`

## Mechanism / Forbidden / Audit

- **Mechanism:** Adaptive max_err-greedy step ordering inside `f_eval_sraa` (Track 2); plateau-gated `aa_unlocked` lifecycle + outer revert to `X_best`/`W_best` (Track 1).
- **Forbidden:** Surface patches to AA acceptance criteria (kSRAARelaxFactor tweaks); changing the global SRAA descent monitor; shrinking K-step inner loop; introducing new heuristics outside the spec; tuning constants to make tests pass.
- **Audit:** RED tests first; benchmark `max_error` on the StepStone K=9 fixture as ground truth; iteration-count proxy for the plateau gate; `bench::mark` with interleaved before/after for any timing claim.

---

## Epic A — Track 2: Adaptive Sort (implement first)

Estimated effort: ~30 minutes. Strict TDD.

### Task A1 — Pre-implementation baseline

- [ ] **Goal:** Capture the current pre-fix `max_error` numbers so the regression and the fix are both falsifiable. No code change.
- [ ] **Files touched:** none (read-only).
- [ ] **Bead:** `bd create -t task -T "A1: SRAA Track 2 baseline"` with body
  `Task adaptive-sort-baseline ! tuning-constants. Capture greenkhorn+plain vs greenkhorn+SRAA max_err on stepstone K=9 prior to fix.`
- [ ] **Commands & expected outputs:**
  - Build current state once:
    `R CMD INSTALL --preclean . 2>&1 | tail -3` -> last line `* DONE (leafblower)`.
  - Capture baseline:
    ```bash
    Rscript -e '
    suppressWarnings(suppressMessages({library(arrow); library(jsonlite); library(leafblower)}))
    df  <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet"); df$uuid <- NULL
    tgt <- lapply(jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json"),
                  function(t){t<-unlist(t); t/sum(t)})
    for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
    r_aa    <- suppressWarnings(harvest(df,tgt,method="greenkhorn",accelerate=TRUE,
                  max_weight=5,min_weight=0,max_iterations=5000L,attach_weights=FALSE,verbose=0))
    r_plain <- suppressWarnings(harvest(df,tgt,method="greenkhorn",accelerate=FALSE,
                  max_weight=5,min_weight=0,max_iterations=5000L,attach_weights=FALSE,verbose=0))
    cat(sprintf("BASELINE me_aa=%.4e me_plain=%.4e ratio=%.3f\n",
        attr(r_aa,"result")$max_error, attr(r_plain,"result")$max_error,
        attr(r_aa,"result")$max_error / attr(r_plain,"result")$max_error))
    ' 2>&1 | tail -1
    ```
    Expected: line beginning `BASELINE me_aa=` with `ratio` >> 1.0 (~1.30 confirming the regression). Any ratio <= 1.001 means the bug is not reproducible — HALT and re-confirm spec.
- [ ] **Compile gate:** `R CMD INSTALL --preclean . 2>&1 | tail -3` -> `* DONE (leafblower)`.
- [ ] **Test gate (no new tests yet):** `Rscript -e "devtools::test()" 2>&1 | tail -3` -> existing FAIL count recorded.
- [ ] **Self-review checklist:**
  - [ ] Baseline numbers logged into the bead comment.
  - [ ] No edits to source.
- [ ] **Commit:** none — baseline is read-only. Append baseline to the bead via `bd comment <id> "BASELINE me_aa=... me_plain=... ratio=..."`.

---

### Task A2 — RED tests

- [ ] **Goal:** Write two failing tests that pin down the K=9 quality regression for greenkhorn+SRAA and raking+SRAA on the StepStone fixture.
- [ ] **Files touched:** `tests/testthat/test-sraa-correct-all-scales.R` (new).
- [ ] **Bead:** `bd create -t task -T "A2: SRAA Track 2 RED tests"` body
  `Task RED-tests ! tuning-tolerance. Two stepstone K=9 tests assert AA<=plain*1.001+1e-10.`
- [ ] **Exact file contents** (`tests/testthat/test-sraa-correct-all-scales.R`):

```r
test_that("T_sraa_adaptive_K9: greenkhorn+AA K=9 max_err <= plain (stepstone)", {
  skip_if_not_installed("arrow"); skip_if_not_installed("jsonlite")
  skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
          "stepstone benchmark data not available")
  df  <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
  df$uuid <- NULL
  tgt <- lapply(jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json"),
                function(t) { t <- unlist(t); t / sum(t) })
  for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
  r_aa    <- suppressWarnings(harvest(df, tgt, method="greenkhorn", accelerate=TRUE,
                                       max_weight=5, min_weight=0, max_iterations=5000L,
                                       attach_weights=FALSE, verbose=0))
  r_plain <- suppressWarnings(harvest(df, tgt, method="greenkhorn", accelerate=FALSE,
                                       max_weight=5, min_weight=0, max_iterations=5000L,
                                       attach_weights=FALSE, verbose=0))
  me_aa    <- attr(r_aa,    "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  expect_lte(me_aa, me_plain * 1.001 + 1e-10,
    label=sprintf("K=9 SRAA (%.4e) must not exceed plain (%.4e)", me_aa, me_plain))
})

test_that("T_sraa_raking_K9: raking+AA K=9 max_err <= plain (stepstone)", {
  skip_if_not_installed("arrow"); skip_if_not_installed("jsonlite")
  skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
          "stepstone benchmark data not available")
  df  <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
  df$uuid <- NULL
  tgt <- lapply(jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json"),
                function(t) { t <- unlist(t); t / sum(t) })
  for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
  r_aa    <- suppressWarnings(harvest(df, tgt, method="raking", accelerate=TRUE,
                                       max_weight=5, min_weight=0, max_iterations=5000L,
                                       attach_weights=FALSE, verbose=0))
  r_plain <- suppressWarnings(harvest(df, tgt, method="raking", accelerate=FALSE,
                                       max_weight=5, min_weight=0, max_iterations=5000L,
                                       attach_weights=FALSE, verbose=0))
  me_aa    <- attr(r_aa,    "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  expect_lte(me_aa, me_plain * 1.001 + 1e-10,
    label=sprintf("K=9 raking+AA (%.4e) must not exceed plain (%.4e)", me_aa, me_plain))
})
```

- [ ] **Compile gate:** `R CMD INSTALL --preclean . 2>&1 | tail -3` -> `* DONE (leafblower)`.
- [ ] **Test gate:** `Rscript -e "devtools::test()" 2>&1 | tail -3` -> `[ FAIL 3 | ...]` (the two new tests RED + any prior failure).
  - If FAIL is fewer than 2 from these new tests, HALT — the bug is not exercised.
- [ ] **Self-review checklist:**
  - [ ] Test names match: `T_sraa_adaptive_K9`, `T_sraa_raking_K9`.
  - [ ] Tolerance is exactly `* 1.001 + 1e-10` (no inflation).
  - [ ] No skip on stepstone fixture aside from existence guard.
- [ ] **Commit:**
  ```bash
  git add tests/testthat/test-sraa-correct-all-scales.R
  git commit -m "$(cat <<'EOF'
  test(sraa): RED tests for K=9 SRAA correctness on stepstone

  Pin down greenkhorn+SRAA and raking+SRAA max_err regression vs plain at K=9
  on the StepStone full-data fixture. Tests assert AA <= plain * 1.001 + 1e-10.
  Currently RED — Track 2 adaptive-sort fix in A3 makes them GREEN.
  EOF
  )"
  ```

---

### Task A3 — Implement adaptive sort in f_eval_sraa

- [ ] **Goal:** Replace the fixed K-step pre-sort in `f_eval_sraa` with the same adaptive `max_element` greedy step used by plain greenkhorn, so SRAA's fixed point coincides with plain's max_err-optimal fixed point.
- [ ] **Files touched:** `src/greenkhorn.cpp`.
- [ ] **Bead:** `bd create -t task -T "A3: SRAA Track 2 adaptive sort in f_eval_sraa"` body
  `Task adaptive-sort ! pre-sort-by-errRp. Remove order_sraa; loop K times calling argmax(errRp) -> greenkhorn_step. Identical step semantics to plain greenkhorn outer loop.`
- [ ] **Pre-edit verification:**
  - Read `src/greenkhorn.cpp` lines 95–195 to confirm landmarks (greenkhorn_step lambda at 98, order_sraa decl at 127, init at 131, sort at 146–148, plain step at 186).
  - Confirm `<algorithm>` is already included for `std::max_element`.
- [ ] **Edits (use Edit tool, not Write):**

  **Edit 1 — remove `order_sraa` declaration (around line 127):**
  Delete the line:
  ```cpp
      std::vector<int>    order_sraa;
  ```

  **Edit 2 — remove `order_sraa.assign(K, 0);` initialization (around line 131):**
  Delete the line:
  ```cpp
          order_sraa.assign(K, 0);
  ```

  **Edit 3 — replace iota+stable_sort+ki-loop body (lines ~146–149) with adaptive loop:**
  Replace:
  ```cpp
          std::iota(order_sraa.begin(), order_sraa.end(), 0);
          std::stable_sort(order_sraa.begin(), order_sraa.end(),
              [&](int a, int b){ return errRp[a] > errRp[b]; });
          for (int ki = 0; ki < K; ki++) greenkhorn_step(order_sraa[ki]);
  ```
  With:
  ```cpp
          for (int ki = 0; ki < K; ki++) {
              int k_star = static_cast<int>(
                  std::max_element(errRp.begin(), errRp.end()) - errRp.begin());
              greenkhorn_step(k_star);
          }
  ```

- [ ] **Post-edit verification:**
  - `grep -n "order_sraa\|stable_sort" src/greenkhorn.cpp` -> no matches.
  - `grep -n "std::max_element(errRp" src/greenkhorn.cpp` -> exactly two hits (existing plain-loop call + new f_eval_sraa call).
- [ ] **Compile gate:** `R CMD INSTALL --preclean . 2>&1 | tail -3` -> `* DONE (leafblower)`.
- [ ] **Test gate:** `Rscript -e "devtools::test()" 2>&1 | tail -3` -> `[ FAIL 0 | ...]` (or pre-existing unrelated failures only — count must not increase).
  - Specifically the two new tests must pass:
    `Rscript -e 'testthat::test_file("tests/testthat/test-sraa-correct-all-scales.R")' 2>&1 | tail -5` -> 2 PASS, 0 FAIL.
- [ ] **Self-review checklist:**
  - [ ] No lingering `order_sraa` references.
  - [ ] greenkhorn_step semantics unchanged.
  - [ ] No iota / stable_sort residue.
  - [ ] errRp updates inside greenkhorn_step still drive the next argmax (verify by inspecting the greenkhorn_step lambda body for the in-place errRp[k_step] write).
  - [ ] No K-related allocations leak from removed `order_sraa.assign`.
- [ ] **Commit:**
  ```bash
  git add src/greenkhorn.cpp
  git commit -m "$(cat <<'EOF'
  fix(sraa): adaptive max_err-greedy step ordering in f_eval_sraa

  Replace fixed K-step pre-sort with the same adaptive argmax(errRp) loop used by
  plain greenkhorn. f_eval_sraa's fixed point now coincides with plain's
  max_err-optimal fixed point, restoring AA quality on K>=3 overlapping-margin
  problems (StepStone K=9). RED tests T_sraa_adaptive_K9 / T_sraa_raking_K9 GREEN.
  EOF
  )"
  ```

---

### Task A4 — Verify Track 2 and decide on Track 1

- [ ] **Goal:** Re-run StepStone K=9 measurement and full test suite; decide whether Epic B is needed.
- [ ] **Files touched:** none.
- [ ] **Bead:** `bd create -t task -T "A4: SRAA Track 2 verification"` body
  `Task verification ! pivot-without-data. Re-measure K=9 ratio. Decide Epic B necessity.`
- [ ] **Commands & expected outputs:**
  - Re-run the A1 baseline snippet (same script) -> `BASELINE me_aa=... me_plain=... ratio=R`.
  - Pass criterion: `R <= 1.001`.
  - Full suite:
    `Rscript -e "devtools::test()" 2>&1 | tail -3` -> FAIL count not greater than the A1-recorded baseline (excluding the 2 new tests, which are now PASS).
- [ ] **Decision:**
  - If `ratio <= 1.001` AND no other regressions: STOP. Epic B is not needed — close A4 and skip Epic B.
  - If `ratio > 1.001` OR a different SRAA test fails: HALT, log finding, proceed to Epic B (do NOT pivot to ad-hoc constant tuning).
- [ ] **Compile gate:** `R CMD INSTALL --preclean . 2>&1 | tail -3` -> `* DONE (leafblower)`.
- [ ] **Test gate:** `Rscript -e "devtools::test()" 2>&1 | tail -3` -> FAIL=0 (or only pre-existing unrelated failures).
- [ ] **Self-review checklist:**
  - [ ] Decision recorded in bead comment with both ratios.
  - [ ] No silent constant tuning.
- [ ] **Commit:** none (read-only verification). Use `bd comment` to log the decision.

---

## Epic B — Track 1: Plateau-Gated AA + Outer Revert

Run only if Epic A's A4 verification fails. Strict TDD.

### Task B1 — Add `aa_unlocked` lifecycle gate + constants to sraa.hpp

- [ ] **Goal:** Introduce a lifecycle bit on `SRAAState` so the outer loops can defer AA activation until plateau is detected, and centralize the four Track-1 constants.
- [ ] **Files touched:** `src/sraa.hpp`.
- [ ] **Bead:** `bd create -t task -T "B1: aa_unlocked + plateau/revert constants"` body
  `Task aa-unlocked-gate ! reset-on-clear. clear() must NOT touch aa_unlocked; outer loop owns lifecycle via enable_aa()/disable_aa().`
- [ ] **Pre-edit verification:**
  - `grep -n "aa_accepted_count\|kSRAARestartGamma\|sraa_step\|kSRAAMinCount" src/sraa.hpp` -> note exact line numbers for each landmark.
  - Read the SRAAState struct, the constants block, and the `sraa_step` "Step 5" guard.
- [ ] **Edits:**

  **Edit 1 — add `aa_unlocked` field + enable/disable methods (insert immediately after `aa_accepted_count`):**

  ```cpp
  // Lifecycle gate owned by the outer loop. NOT reset by clear() — see
  // enable_aa()/disable_aa() below. Outer loop sets it true after a plateau
  // is detected and false after an outer revert, so AA cannot fire on the
  // initial transient or after the iterate is restored to a known basin.
  bool aa_unlocked = false;
  void enable_aa()  { aa_unlocked = true;  }
  void disable_aa() { aa_unlocked = false; }
  ```

  **Edit 2 — add four constants immediately after `kSRAARestartGamma`:**

  ```cpp
  // Track 1: plateau-gated AA + outer revert.
  static constexpr double kSRAAPlateauEps      = 1e-3;  // 0.1%/iter improvement => plateau
  static constexpr int    kSRAAPlateauWindow   = 4;     // 4 consecutive plateau iters => AA enabled
  static constexpr double kSRAAOuterSlack      = 0.10;  // 10% above best_errRp => regression
  static constexpr int    kSRAAOuterStallWindow = 5;    // 5 regressed outer iters => revert
  ```

  **Edit 3 — extend Step 5 guard in `sraa_step` to also short-circuit when `!state.aa_unlocked`:**
  Replace:
  ```cpp
      if (state.count < kSRAAMinCount) {
          std::swap(X, state.F_cur);
          return {false, 1, err_plain};
      }
  ```
  With:
  ```cpp
      // Not enough history OR outer loop has not unlocked AA -> plain step
      if (state.count < kSRAAMinCount || !state.aa_unlocked) {
          std::swap(X, state.F_cur);
          return {false, 1, err_plain};
      }
  ```

  **Edit 4 — confirm `clear()` does NOT touch `aa_unlocked`:**
  - Read the existing `clear()` body. If it touches `aa_unlocked`, HALT — that is the spec contract.
  - Add a comment immediately above `clear()`:
    ```cpp
    // Resets AA history (count, F/G/E buffers, accepted/rejected counters).
    // Does NOT reset aa_unlocked — lifecycle is owned by the outer loop.
    ```

- [ ] **Compile gate:** `R CMD INSTALL --preclean . 2>&1 | tail -3` -> `* DONE (leafblower)`.
- [ ] **Test gate (no behaviour change yet because nobody calls enable_aa):**
  - With Track 2 already merged, all stepstone K=9 tests still PASS (B1 is purely adding dormant infrastructure).
  - `Rscript -e "devtools::test()" 2>&1 | tail -3` -> FAIL not greater than post-A4.
- [ ] **Self-review checklist:**
  - [ ] `aa_unlocked` default = `false`.
  - [ ] No outer-loop call to `enable_aa()` yet (callers come in B2/B3).
  - [ ] Constants in the same `static constexpr` block as `kSRAARestartGamma`.
  - [ ] Step 5 guard uses logical OR with `!state.aa_unlocked`.
  - [ ] No edits to other AA acceptance criteria.
- [ ] **Commit:**
  ```bash
  git add src/sraa.hpp
  git commit -m "$(cat <<'EOF'
  feat(sraa): aa_unlocked lifecycle gate + plateau/revert constants

  Add an outer-loop-owned aa_unlocked bit on SRAAState plus four Track-1
  constants (kSRAAPlateauEps, kSRAAPlateauWindow, kSRAAOuterSlack,
  kSRAAOuterStallWindow). sraa_step now short-circuits to a plain step when
  the outer loop has not yet unlocked AA. clear() deliberately does NOT
  touch aa_unlocked — outer loop owns lifecycle.
  EOF
  )"
  ```

---

### Task B2 — Plateau detection + outer revert in greenkhorn.cpp

- [ ] **Goal:** Wire the outer loop to detect plateau and trigger AA, and to revert to `X_best` when AA escapes the basin. **MUST MERGE** with the existing `curr_max` / `best_errRp` block in greenkhorn.cpp (~lines 192–196), not duplicate it.
- [ ] **Files touched:** `src/greenkhorn.cpp`.
- [ ] **Bead:** `bd create -t task -T "B2: greenkhorn outer plateau gate + revert"` body
  `Task outer-revert-greenkhorn ! duplicate-curr_max. MERGE plateau/revert into existing best_errRp update at lines ~192-196. Use kSRAA* constants from sraa.hpp.`
- [ ] **Pre-edit verification:**
  - Read `src/greenkhorn.cpp` lines 60–210 — locate the existing block:
    ```cpp
    double curr_max = *std::max_element(errRp.begin(), errRp.end());
    if (curr_max < best_errRp) {
        best_errRp = curr_max; res.best_error = best_errRp;
        res.best_iter = res.iterations; X_best = X;
    }
    ```
  - Locate where outer-loop locals are declared (`best_errRp`, `X_best`, `S_flat`, `W`, `S_stride`, `M`, `K`, `errRp`, `compute_errRp_k`, `ct.g_per_cell`, `st.cat_counts`, `st.accelerate`, `grk_sraa`).
  - Confirm `<limits>` is included (for `std::numeric_limits<double>::infinity()`).
- [ ] **Edits:**

  **Edit 1 — add three locals immediately above the outer `for` loop:**

  ```cpp
  double prev_outer_quality = std::numeric_limits<double>::infinity();
  int    plateau_count      = 0;
  int    outer_stall_count  = 0;
  ```

  **Edit 2 — REPLACE the existing curr_max/best_errRp block (lines ~192–196) with the merged plateau-gate + revert block:**

  Old (delete exactly):
  ```cpp
          double curr_max = *std::max_element(errRp.begin(), errRp.end());
          if (curr_max < best_errRp) {
              best_errRp = curr_max; res.best_error = best_errRp;
              res.best_iter = res.iterations; X_best = X;
          }
  ```

  New (insert in its place):
  ```cpp
          // Track outer quality (single curr_max — do not duplicate above).
          double curr_max = *std::max_element(errRp.begin(), errRp.end());
          if (curr_max < best_errRp) {
              best_errRp     = curr_max;
              res.best_error = best_errRp;
              res.best_iter  = res.iterations;
              X_best         = X;
          }

          if (st.accelerate && K > 0) {
              // Plateau detection: enable AA when outer-quality improvement saturates.
              if (std::isfinite(prev_outer_quality)) {
                  double impr = (prev_outer_quality - curr_max)
                              / std::max(prev_outer_quality, 1e-15);
                  plateau_count = (impr < kSRAAPlateauEps) ? plateau_count + 1 : 0;
                  if (plateau_count >= kSRAAPlateauWindow) grk_sraa.enable_aa();
              }
              prev_outer_quality = curr_max;

              // Outer revert: AA escaped the correct basin -> snap back to X_best.
              if (grk_sraa.aa_unlocked
                  && curr_max > best_errRp * (1.0 + kSRAAOuterSlack)) {
                  if (++outer_stall_count >= kSRAAOuterStallWindow) {
                      X = X_best;
                      grk_sraa.clear();        // history reset; aa_unlocked untouched
                      grk_sraa.disable_aa();   // re-require plateau before AA fires again
                      plateau_count     = 0;
                      outer_stall_count = 0;

                      // Rebuild W, S_flat, errRp from X = X_best.
                      W = 0.0;
                      std::fill(S_flat.begin(), S_flat.end(), 0.0);
                      for (int c = 0; c < M; c++) {
                          W += X[c];
                          for (int k = 0; k < K; k++) {
                              int g = ct.g_per_cell[k][c];
                              if (g >= 0 && g < st.cat_counts[k])
                                  S_flat[k * S_stride + g] += X[c];
                          }
                      }
                      for (int k = 0; k < K; k++) errRp[k] = compute_errRp_k(k);
                      best_errRp         = *std::max_element(errRp.begin(), errRp.end());
                      prev_outer_quality = best_errRp;
                  }
              } else {
                  outer_stall_count = 0;
              }
          }
  ```

- [ ] **Post-edit verification:**
  - `grep -nc "double curr_max = \*std::max_element(errRp" src/greenkhorn.cpp` -> `1` (no duplicate).
  - `grep -n "grk_sraa.enable_aa\|grk_sraa.disable_aa\|grk_sraa.aa_unlocked" src/greenkhorn.cpp` -> three matches.
  - `grep -n "kSRAAPlateauEps\|kSRAAPlateauWindow\|kSRAAOuterSlack\|kSRAAOuterStallWindow" src/greenkhorn.cpp` -> four matches.
- [ ] **Compile gate:** `R CMD INSTALL --preclean . 2>&1 | tail -3` -> `* DONE (leafblower)`.
- [ ] **Test gate:** `Rscript -e "devtools::test()" 2>&1 | tail -3` -> FAIL count <= post-A4 (no regression introduced; B4 will add the dedicated B-track tests).
- [ ] **Self-review checklist:**
  - [ ] Exactly one `curr_max` computation per outer iteration.
  - [ ] `clear()` is called before `disable_aa()` to keep aa_unlocked semantics explicit.
  - [ ] Revert rebuilds W, S_flat, errRp from `X = X_best` and reseeds `best_errRp`/`prev_outer_quality`.
  - [ ] No edits outside the outer loop.
  - [ ] No silent constant tweaks.
- [ ] **Commit:**
  ```bash
  git add src/greenkhorn.cpp
  git commit -m "$(cat <<'EOF'
  feat(sraa): plateau-gated AA + outer revert in greenkhorn outer loop

  Merge plateau detection and X_best revert into the existing curr_max /
  best_errRp block. Defer AA activation until outer-quality improvement
  saturates (kSRAAPlateauEps over kSRAAPlateauWindow iters). On a sustained
  regression > kSRAAOuterSlack for kSRAAOuterStallWindow iters, snap back to
  X_best, rebuild W/S_flat/errRp, clear AA history and re-lock aa_unlocked.
  EOF
  )"
  ```

---

### Task B3 — Same pattern in raking.cpp (W_best / best_metric_seen / rk_sraa)

- [ ] **Goal:** Mirror B2 in `src/raking.cpp` using raking's existing names.
- [ ] **Files touched:** `src/raking.cpp`.
- [ ] **Bead:** `bd create -t task -T "B3: raking outer plateau gate + revert"` body
  `Task outer-revert-raking ! diverge-from-greenkhorn. Same pattern as B2 with W_best, best_metric_seen, rk_sraa. Revert: X = W_best, rebuild W/S_flat/errRp from X.`
- [ ] **Pre-edit verification:**
  - `grep -n "W_best\|best_metric_seen\|rk_sraa\|std::max_element(errRp" src/raking.cpp` — confirm exact landmark lines.
  - Read the analogous outer block; identify the single existing curr_max/best update.
- [ ] **Edits — apply the same pattern as B2 with substitutions:**

  - `X_best`         -> `W_best`
  - `best_errRp`     -> `best_metric_seen`
  - `grk_sraa`       -> `rk_sraa`
  - Revert assignment is `X = W_best;` (raking's iterate name).
  - Rebuild `W`, `S_flat`, `errRp` from the (now-restored) `X = W_best` using the same per-cell loop pattern raking already uses elsewhere — confirm by reading the existing init block once and matching it byte-for-byte (no reformatting).
  - Insert the same three locals (`prev_outer_quality`, `plateau_count`, `outer_stall_count`) immediately above raking's outer `for`.
  - MERGE with the existing single curr_max/best_metric_seen block — do NOT duplicate.

- [ ] **Post-edit verification:**
  - `grep -nc "double curr_max = \*std::max_element(errRp" src/raking.cpp` -> `1`.
  - `grep -n "rk_sraa.enable_aa\|rk_sraa.disable_aa\|rk_sraa.aa_unlocked" src/raking.cpp` -> three matches.
  - `grep -n "W_best" src/raking.cpp` -> existing references plus the one new `X = W_best;` line.
- [ ] **Compile gate:** `R CMD INSTALL --preclean . 2>&1 | tail -3` -> `* DONE (leafblower)`.
- [ ] **Test gate:** `Rscript -e "devtools::test()" 2>&1 | tail -3` -> FAIL count not greater than post-B2.
- [ ] **Self-review checklist:**
  - [ ] Names substituted consistently (no stray `X_best`/`grk_sraa`/`best_errRp` in raking.cpp).
  - [ ] Single curr_max in raking.cpp.
  - [ ] Rebuild block matches raking's existing init style.
  - [ ] No edits outside the raking outer loop.
- [ ] **Commit:**
  ```bash
  git add src/raking.cpp
  git commit -m "$(cat <<'EOF'
  feat(sraa): plateau-gated AA + outer revert in raking outer loop

  Mirror the greenkhorn Track-1 changes using raking's existing names
  (W_best, best_metric_seen, rk_sraa). Single merged curr_max block; revert
  snaps X back to W_best and rebuilds W/S_flat/errRp before re-locking
  aa_unlocked.
  EOF
  )"
  ```

---

### Task B4 — Plateau-gate + outer-revert tests + verify

- [ ] **Goal:** Add two RED tests that pin the plateau gate (proxy: iteration count) and the K=6 outer-revert quality, then confirm all four SRAA correctness tests are GREEN.
- [ ] **Files touched:** `tests/testthat/test-sraa-correct-all-scales.R` (append), no other code edits.
- [ ] **Bead:** `bd create -t task -T "B4: plateau-gate + outer-revert tests"` body
  `Task track1-tests ! tighter-tolerance. Iteration-count proxy for plateau gate; K=6 stepstone-style synthetic for outer revert.`
- [ ] **Append exactly** to `tests/testthat/test-sraa-correct-all-scales.R`:

```r
test_that("T_sraa_plateau_gate: greenkhorn+AA early iters run plain-only before plateau", {
  # With kSRAAPlateauWindow=4, AA cannot fire in first 4 outer iterations.
  # Proxy: max_iterations=4 means 4 outer iters; iters_aa must equal iters_plain
  # (all plain: K*1 each, not K*2 for AA accepted steps).
  set.seed(99); n <- 2000L
  df  <- data.frame(x=factor(sample(letters[1:3],n,TRUE)), y=factor(sample(c("M","F"),n,TRUE)))
  tgt <- list(x=c(a=0.3,b=0.4,c=0.3), y=c(M=0.5,F=0.5))
  r_aa    <- suppressWarnings(harvest(df,tgt,method="greenkhorn",accelerate=TRUE,
                                       max_iterations=4L,attach_weights=FALSE))
  r_plain <- suppressWarnings(harvest(df,tgt,method="greenkhorn",accelerate=FALSE,
                                       max_iterations=4L,attach_weights=FALSE))
  expect_equal(attr(r_aa,"result")$iterations, attr(r_plain,"result")$iterations,
    label="Before plateau, AA must not fire — iterations must match plain (K*1 per outer iter)")
})

test_that("T_sraa_outer_revert: greenkhorn+AA K=6 quality recovers to <= plain with Track 1", {
  set.seed(42); n <- 8000L
  df <- data.frame(gender=factor(sample(c("M","F"),n,TRUE)),
                   time=factor(sample(1:3,n,TRUE)), age=factor(sample(1:4,n,TRUE)))
  df$gt <- factor(paste0(df$gender,df$time))
  df$ga <- factor(paste0(df$gender,df$age))
  df$ta <- factor(paste0(df$time,df$age))
  gt_t <- table(df$gt)/n; ga_t <- table(df$ga)/n; ta_t <- table(df$ta)/n
  tgt <- list(gender=c(M=0.48,F=0.52), time=setNames(c(0.4,0.35,0.25),1:3),
    age=setNames(c(0.3,0.25,0.25,0.2),1:4),
    gt={t<-setNames(as.numeric(gt_t)*c(0.95,1.02,0.98,1.03,0.97,1.05),names(gt_t));t/sum(t)},
    ga={t<-setNames(as.numeric(ga_t),names(ga_t));t/sum(t)},
    ta={t<-setNames(as.numeric(ta_t),names(ta_t));t/sum(t)})
  tgt <- lapply(tgt, function(t) t/sum(t))
  r_aa    <- suppressWarnings(harvest(df,tgt,method="greenkhorn",accelerate=TRUE,
                                       max_iterations=200L,attach_weights=FALSE))
  r_plain <- suppressWarnings(harvest(df,tgt,method="greenkhorn",accelerate=FALSE,
                                       max_iterations=200L,attach_weights=FALSE))
  me_aa    <- attr(r_aa,   "result")$max_error
  me_plain <- attr(r_plain,"result")$max_error
  expect_lte(me_aa, me_plain * 1.001 + 1e-10,
    label=sprintf("K=6 Track1: AA (%.2e) must not exceed plain (%.2e)", me_aa, me_plain))
})
```

- [ ] **Compile gate:** `R CMD INSTALL --preclean . 2>&1 | tail -3` -> `* DONE (leafblower)`.
- [ ] **Test gate:** `Rscript -e "devtools::test()" 2>&1 | tail -3` -> `[ FAIL 0 | ...]` (or only pre-existing unrelated failures).
  - Targeted: `Rscript -e 'testthat::test_file("tests/testthat/test-sraa-correct-all-scales.R")' 2>&1 | tail -5` -> 4 PASS, 0 FAIL.
- [ ] **Self-review checklist:**
  - [ ] All four tests live in one file.
  - [ ] Plateau-gate test asserts equality of iteration counts (no fragile timing).
  - [ ] K=6 test seed-deterministic.
  - [ ] No constant tuning to make tests pass — if a test fails, fix the algorithm or HALT.
  - [ ] All four SRAA correctness tests still PASS together (Track-2 tests must not regress).
- [ ] **Commit:**
  ```bash
  git add tests/testthat/test-sraa-correct-all-scales.R
  git commit -m "$(cat <<'EOF'
  test(sraa): plateau-gate iteration proxy + K=6 outer-revert correctness

  Add T_sraa_plateau_gate (asserts AA cannot fire before kSRAAPlateauWindow,
  via iteration-count equality with plain) and T_sraa_outer_revert (K=6
  overlapping-margin synthetic where Track 1 must keep AA<=plain*1.001).
  Together with T_sraa_adaptive_K9 / T_sraa_raking_K9 the suite pins
  correct-all-scales behaviour.
  EOF
  )"
  ```

---

## Done Definition (entire effort)

- [ ] Track 2 (Epic A) merged. Stepstone K=9 ratio AA/plain <= 1.001 on both greenkhorn and raking.
- [ ] If A4 fails, Track 1 (Epic B) merged. All four tests in `tests/testthat/test-sraa-correct-all-scales.R` GREEN.
- [ ] No SRAA constants tuned to "make tests pass" — every constant traces to spec.
- [ ] No surface patches to AA acceptance criteria.
- [ ] `R CMD INSTALL --preclean . 2>&1 | tail -3` -> `* DONE (leafblower)`.
- [ ] `Rscript -e "devtools::test()" 2>&1 | tail -3` shows no new FAIL beyond pre-existing unrelated failures.
- [ ] Each commit message follows Conventional Commits and explains "why", not "what".
- [ ] No `--no-verify`, no `--amend`, no `git push --force`.
- [ ] After all commits, push: `git pull --rebase && git push`.

## Halt Rules (`SPEC_FAILURE`)

- If A3's adaptive sort change does not bring K=9 ratio to <= 1.001 AND the K=9 raking test still fails after Epic B: emit `SPEC_FAILURE`. Do NOT pivot to constant tuning, do NOT relax tolerances, do NOT change SRAA acceptance criteria.
- If `clear()` is found to reset `aa_unlocked` (contract violation): emit `SPEC_FAILURE`.
- If a duplicate `curr_max` computation appears in either greenkhorn.cpp or raking.cpp after B2/B3: emit `SPEC_FAILURE`.
