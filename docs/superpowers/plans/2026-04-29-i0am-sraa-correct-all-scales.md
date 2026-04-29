# SRAA-m Correct Acceleration for All K Scales — Combined Fix

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix SRAA-m 35% quality regression on K>=3 overlapping-margin problems by combining adaptive-sort f_eval (eliminates basin multiplicity) with an outer-loop revert mechanism (catches transient AA overshoots).

**Architecture:** Two changes applied simultaneously: (1) f_eval_sraa in greenkhorn.cpp switches from fixed sort (K margins sorted once at round entry) to adaptive sort (re-sort errRp after each greedy step = plain greenkhorn). (2) Both greenkhorn.cpp and raking.cpp add a 15-line outer revert: if outer quality regresses >10% above best_errRp for 5 consecutive outer iterations, X reverts to X_best and AA history clears. No plateau gating, no aa_unlocked field — plateau gating was needed only to protect best_errRp from wrong-basin contamination, which adaptive sort eliminates.

**Tech Stack:** C++17 (sraa.hpp, greenkhorn.cpp, raking.cpp), R testthat3, stepstone benchmark fixture.

**Spec:** `docs/superpowers/specs/2026-04-29-i0am-sraa-correct-all-scales.md`

---

## Task 1: Pre-baseline + RED Tests

**Mechanism:** Capture pre-fix quality baseline on stepstone fixture, then add 3 failing testthat3 tests that encode the combined-fix acceptance criteria. Tests must FAIL on `master` (RED) before any C++ changes land in Task 2.

**Forbidden:** No C++ edits in this task. No skipping the baseline capture. No `expect_true(TRUE)` placeholders. Do not gate on parquet absence with `skip()` silently — emit `skip_if_not(file.exists(...))` with explicit reason.

**Audit:** Tests must directly drive `leafblower::harvest()` (or the SRAA-m–exposing R wrapper) on the stepstone parquet fixture; assertions must compare per-K final `max(errRp)` against numerically explicit thresholds derived from the spec, not against placeholder constants.

**Beads ticket:** `bd create --title "SRAA-m T1: pre-baseline + RED tests" --description "Task [stepstone baseline + 3 RED tests in test-sraa-global.R] ! [skipping baseline; expect_true(TRUE); silent skip]"` — capture the issued ID into `${T1_BEAD_ID}` for the close at the end of the task.

### Steps

- [ ] **1.1 Read inputs (no edits)**
  - [ ] Read `docs/superpowers/specs/2026-04-29-i0am-sraa-correct-all-scales.md` end-to-end. Extract: (a) per-K acceptance thresholds for `T_sraa_adaptive_K9`, `T_sraa_raking_K9`, `T_sraa_outer_revert`; (b) stepstone fixture path; (c) which R entry point exposes SRAA-m (`lba(..., method = "...")` vs internal helper).
  - [ ] Read `tests/testthat/test-sraa-global.R` fully (current header, helpers, existing test names) to align style and avoid name collisions.
  - [ ] Read `src/greenkhorn.cpp` lines 1–230 to confirm the line-number anchors used by Task 2 (`order_sraa` decl ~127, init ~131, iota+sort ~146–149, best_errRp block ~192–196, secondary `X_best=X` ~206). If anchors drift, record actual line numbers in `.wolf/memory.md` and proceed — Task 2 patches are described by surrounding code, not line numbers alone.
  - [ ] Read `src/raking.cpp` SRAA path to locate the `if (r.err_result < best_metric_seen) { ... W_best = X; }` block referenced by Task 3.
  - [ ] Read `src/sraa.hpp` to locate `kSRAARestartGamma` (insertion anchor for the two new constants).

- [ ] **1.2 Capture pre-fix baseline**
  - [ ] From repo root: `R CMD INSTALL --preclean .`. Compile must succeed on master before measuring baseline.
  - [ ] Run the stepstone benchmark exactly as the spec prescribes; redirect stdout+stderr to `benchmarks/sraa-m-baseline-pre.log`:
    ```bash
    Rscript benchmarks/stepstone_all_methods.R 2>&1 | tee benchmarks/sraa-m-baseline-pre.log
    ```
    If `benchmarks/stepstone_all_methods.R` does not exist, halt and grep `bench/` and `benchmarks/` for the script name actually referenced by the spec; do not invent one.
  - [ ] Extract per-K final `max(errRp)` for SRAA-m at K ∈ {3, 5, 9} and the corresponding plain-greenkhorn / raking baselines. Persist to `benchmarks/sraa-m-baseline-pre.csv` with columns `method,K,max_errRp,iters,wallclock_s`. Use `Rscript -e` writing via `data.table::fwrite()`; do not hand-craft CSV.
  - [ ] `git add benchmarks/sraa-m-baseline-pre.{log,csv}` (do not commit yet — the RED-test commit owns these alongside the tests).

- [ ] **1.3 Append RED tests**
  - [ ] Append to `tests/testthat/test-sraa-global.R`:
    - `test_that("T_sraa_adaptive_K9: SRAA-m matches plain greenkhorn within 2% on K=9 overlapping-margin stepstone", { ... })` — drives `harvest()` with SRAA-m and with plain greenkhorn on the K=9 stepstone fixture; asserts `expect_lt(max_errRp_sraa_m, max_errRp_greenkhorn * 1.02)`.
    - `test_that("T_sraa_raking_K9: SRAA-raking matches plain raking within 2% on K=9 fixture", { ... })` — analogous, comparing SRAA-accelerated raking vs plain raking.
    - `test_that("T_sraa_outer_revert: outer quality never exceeds best-seen by more than 10% for 5 consecutive iters", { ... })` — uses K=6 cross-margin DGP (set.seed(42), n=8000); asserts me_aa <= me_plain * 1.001 + 1e-10. harvest() has no per-iteration trace — use final attr(r,"result")$max_error.
  - [ ] Each test must `skip_if_not(file.exists(stepstone_parquet_path()), "stepstone parquet fixture not present")` with `stepstone_parquet_path()` defined once at the top of the file (or already present — reuse it).
  - [ ] No `expect_true(TRUE)` filler. Every assertion must reference a real numeric output of `harvest()`.

- [ ] **1.4 Verify RED locally**
  - [ ] `R CMD INSTALL --preclean .` (already current; cheap reinstall to prove the build still works).
  - [ ] `Rscript -e 'testthat::test_file("tests/testthat/test-sraa-global.R", reporter = "summary")' 2>&1 | tee benchmarks/sraa-m-red-evidence.log`.
  - [ ] Confirm the three new tests are present and **FAIL** (or `skip` only if parquet truly absent — record which in the log). All three must be RED on master, not skipped, when the parquet exists. If any of them passes accidentally, the assertion threshold is wrong — tighten before proceeding.

- [ ] **1.5 Commit RED state**
  - [ ] Files staged: `tests/testthat/test-sraa-global.R`, `benchmarks/sraa-m-baseline-pre.log`, `benchmarks/sraa-m-baseline-pre.csv`, `benchmarks/sraa-m-red-evidence.log`.
  - [ ] Commit:
    ```bash
    git commit -m "$(cat <<'EOF'
    test(sraa): T1 — RED tests + pre-fix stepstone baseline

    Adds T_sraa_adaptive_K9, T_sraa_raking_K9, T_sraa_outer_revert to
    tests/testthat/test-sraa-global.R. Captures pre-fix per-K max(errRp)
    on stepstone K∈{3,5,9} for SRAA-m vs plain greenkhorn/raking
    (benchmarks/sraa-m-baseline-pre.{log,csv}). Tests fail on master,
    establishing the RED baseline for the combined fix in T2/T3.
    EOF
    )"
    ```
  - [ ] `bd update ${T1_BEAD_ID} --status closed --comment "RED tests appended; baseline captured; failing as expected."`

### Compile gate

`R CMD INSTALL --preclean .` succeeds.

### Test gate

`Rscript -e 'testthat::test_file("tests/testthat/test-sraa-global.R", reporter = "summary")'` shows all three new tests FAIL (parquet present) or SKIP with explicit reason (parquet absent). No new test passes.

---

## Task 2: Combined Fix in greenkhorn.cpp + sraa.hpp

**Mechanism:** Apply (A) two new constants in `sraa.hpp`, (B) adaptive-sort replacement in `greenkhorn.cpp::f_eval_sraa`, and (C) outer-loop revert merged into the existing best-tracking block. All three sub-changes ship in one commit because they are interdependent: the revert depends on the constants, and removing the fixed sort without the revert leaves no fallback for AA overshoots.

**Forbidden:** No plateau gating. No `aa_unlocked` field. No new helper functions. No reformatting of unrelated lines. Do not modify the secondary `X_best = X` assignment at ~line 206 inside the convergence check. No tweaking of `kSRAARestartGamma` or any other existing constant. No partial commits — A, B, C land together.

**Audit:** A unit-level spy is impractical for inline C++ logic; instead, T1's `T_sraa_adaptive_K9` and `T_sraa_outer_revert` tests directly observe the post-conditions (per-step argmax dispatch via the quality match; revert via the run-length assertion). The compile gate plus the full testthat suite is the authoritative audit.

**Beads ticket:** `bd create --title "SRAA-m T2: adaptive-sort + outer revert in greenkhorn.cpp/sraa.hpp" --description "Task [adaptive argmax f_eval + 15-line outer revert merged with best-tracking; 2 constants in sraa.hpp] ! [plateau gating; aa_unlocked field; new helpers; partial commits; touching X_best=X at ~line 206]"` — capture as `${T2_BEAD_ID}`.

### Steps

- [ ] **2.1 Re-read anchors**
  - [ ] Re-read `src/sraa.hpp` around `kSRAARestartGamma` to confirm exact insertion site and surrounding `static constexpr` style (semicolons, comment style, indentation).
  - [ ] Re-read `src/greenkhorn.cpp` around: `f_eval_sraa` (lines ~120–155 for the order_sraa removal); the outer for-loop and the existing `if (curr_max < best_errRp) { ... }` block (~lines 188–200); and the secondary `X_best = X` at ~line 206 (verify it is inside the convergence/return path and confirm the "harmless no-op after revert" claim in the spec by tracing control flow).
  - [ ] Confirm `S_stride`, `S_flat`, `W`, `ct.g_per_cell`, `st.cat_counts`, `compute_errRp_k`, and `errRp` are all in scope where the revert block lands. If any name differs (e.g., `S_stride` is locally `K_stride` or computed inline), record the actual identifier in `.wolf/memory.md` and use the actual name in the patch — do not invent.

- [ ] **2.2 sraa.hpp — add 2 constants**
  - [ ] Edit `src/sraa.hpp`. Immediately after the line declaring `kSRAARestartGamma`, insert:
    ```cpp
    static constexpr double kSRAAOuterSlack = 0.10;
    static constexpr int    kSRAAOuterStallWindow = 5;
    ```
  - [ ] Match surrounding indentation exactly. No new include. No comment unless the file's local style demands one — then mirror the form used near `kSRAARestartGamma`.

- [ ] **2.3 greenkhorn.cpp — adaptive sort in f_eval_sraa**
  - [ ] In `src/greenkhorn.cpp::f_eval_sraa`:
    - [ ] **REMOVE** the declaration `std::vector<int> order_sraa;` (~line 127).
    - [ ] **REMOVE** the init line `order_sraa.assign(K, 0);` (~line 131).
    - [ ] **REMOVE** the fixed-sort block (~lines 146–149): the `std::iota(order_sraa.begin(), order_sraa.end(), 0);` followed by `std::stable_sort(order_sraa.begin(), order_sraa.end(), [&](int a, int b){ return errRp[a] > errRp[b]; });` and the dependent `for (int idx = 0; idx < K; ++idx) { int k = order_sraa[idx]; greenkhorn_step(k); }`.
    - [ ] **ADD** in place of the removed loop:
      ```cpp
      for (int ki = 0; ki < K; ki++) {
          int k_star = (int)(std::max_element(errRp.begin(), errRp.end()) - errRp.begin());
          greenkhorn_step(k_star);
      }
      ```
  - [ ] Confirm `<algorithm>` is already included (it is required for `std::max_element`); if not, add it once at the top of `greenkhorn.cpp`. Do not add `<numeric>` — `std::iota` is no longer needed in this function but may be used elsewhere; leave the include alone.

- [ ] **2.4 greenkhorn.cpp — outer revert merged with best-tracking**
  - [ ] Immediately before the outer `for` loop that bumps the SRAA outer iteration counter, insert: `int outer_stall_count = 0;`.
  - [ ] Locate the existing block (~lines 192–196):
    ```cpp
    if (curr_max < best_errRp) {
        best_errRp = curr_max;
        X_best = X;
    }
    ```
  - [ ] Extend it in place to read exactly:
    ```cpp
    if (curr_max < best_errRp) {
        best_errRp = curr_max;
        X_best = X;
        outer_stall_count = 0;
    } else if (curr_max > best_errRp * (1.0 + kSRAAOuterSlack)) {
        if (++outer_stall_count >= kSRAAOuterStallWindow) {
            X = X_best;
            grk_sraa.clear();
            outer_stall_count = 0;
            W = 0.0;
            std::fill(S_flat.begin(), S_flat.end(), 0.0);
            for (int c = 0; c < M; c++) {
                W += X[c];
                for (int k = 0; k < K; k++) {
                    int g = ct.g_per_cell[k][c];
                    if (g >= 0 && g < st.cat_counts[k]) S_flat[k * S_stride + g] += X[c];
                }
            }
            for (int k = 0; k < K; k++) errRp[k] = compute_errRp_k(k);
            best_errRp = *std::max_element(errRp.begin(), errRp.end());
        }
    } else {
        outer_stall_count = 0;
    }
    ```
  - [ ] Do **not** modify the secondary `X_best = X;` at ~line 206. The spec confirms it is a harmless no-op after revert (X already equals X_best at that point in the convergence path).

- [ ] **2.5 Compile gate**
  - [ ] `R CMD INSTALL --preclean .` from repo root. Must succeed with zero warnings introduced by this change. If a warning fires (e.g., unused-variable on a leftover `order_sraa`-related symbol), halt and re-audit the removals — do not paper over with `(void)`.

- [ ] **2.6 Run targeted tests**
  - [ ] `Rscript -e 'testthat::test_file("tests/testthat/test-sraa-global.R", reporter = "summary")' 2>&1 | tee benchmarks/sraa-m-t2-tests.log`.
  - [ ] Expectation: `T_sraa_adaptive_K9` and `T_sraa_outer_revert` now PASS. `T_sraa_raking_K9` may still fail (it depends on Task 3). If `T_sraa_adaptive_K9` or `T_sraa_outer_revert` still fail, **do not pivot**: halt, output `SPEC_FAILURE`, log root cause to `.wolf/buglog.json`, and stop.
  - [ ] Run the full package suite: `Rscript -e 'devtools::test()' 2>&1 | tee benchmarks/sraa-m-t2-fulltest.log`. Zero new failures outside the still-RED `T_sraa_raking_K9`.

- [ ] **2.7 Commit**
  - [ ] Stage exactly: `src/sraa.hpp`, `src/greenkhorn.cpp`, `benchmarks/sraa-m-t2-tests.log`, `benchmarks/sraa-m-t2-fulltest.log`.
  - [ ] Commit:
    ```bash
    git commit -m "$(cat <<'EOF'
    fix(sraa): T2 — adaptive-sort f_eval + outer revert in greenkhorn

    Replaces fixed K-margin sort in f_eval_sraa with per-step argmax
    dispatch (eliminates basin multiplicity that drove the K>=3 quality
    regression). Adds a 15-line outer-loop revert merged with the
    existing best_errRp tracking block: if curr_max exceeds best_errRp
    by >10% for 5 consecutive outer iters, X reverts to X_best, AA
    history clears, and S_flat/W/errRp are rebuilt from scratch. Adds
    kSRAAOuterSlack and kSRAAOuterStallWindow to sraa.hpp.

    Greens T_sraa_adaptive_K9 and T_sraa_outer_revert; T_sraa_raking_K9
    remains red until T3.
    EOF
    )"
    ```
  - [ ] `bd update ${T2_BEAD_ID} --status closed --comment "Adaptive sort + outer revert landed; 2/3 SRAA tests green."`

### Compile gate

`R CMD INSTALL --preclean .` succeeds, no new warnings.

### Test gate

`T_sraa_adaptive_K9` PASS; `T_sraa_outer_revert` PASS; full suite shows no regressions outside `T_sraa_raking_K9` (expected red until T3).

---

## Task 3: Outer Revert in raking.cpp + Verify

**Mechanism:** Mirror the outer-revert pattern in `raking.cpp`'s SRAA path. Raking's quality metric is the scalar `r.err_result`, not a vector argmax — so no `S_flat` rebuild is needed (the next `F_eval` recomputes derived state). Then run the full quality-gate stack: compile, full test suite, and the post-fix benchmark to compare against the T1 baseline.

**Forbidden:** No `S_flat` rebuild in raking (raking's `F_eval` reconstructs derived state on the next call — adding a manual rebuild duplicates work and risks divergence). No introduction of vector `errRp` semantics into raking. No changes to `greenkhorn.cpp` or `sraa.hpp` in this task. No skipping the full benchmark re-run.

**Audit:** `T_sraa_raking_K9` directly observes the post-condition (raking SRAA matches plain raking within 2% on K=9). The benchmark CSV diff (`benchmarks/sraa-m-baseline-pre.csv` vs `benchmarks/sraa-m-baseline-post.csv`) audits across all K.

**Beads ticket:** `bd create --title "SRAA-m T3: outer revert in raking.cpp + verify" --description "Task [scalar-metric outer revert in raking SRAA path; full suite + post-fix benchmark] ! [S_flat rebuild in raking; touching greenkhorn/sraa.hpp; skipping benchmark]"` — capture as `${T3_BEAD_ID}`.

### Steps

- [ ] **3.1 Re-read anchor**
  - [ ] Re-read `src/raking.cpp` SRAA block. Locate:
    - the `int rk_outer` (or equivalent) loop entry,
    - the `if (r.err_result < best_metric_seen) { ... W_best = X; }` block (this is the insertion target),
    - the `rk_sraa` history container (used for `.clear()`).
    Record actual identifier names in `.wolf/memory.md` if any differ from the spec's shorthand; use actual names in the patch.

- [ ] **3.2 raking.cpp — scalar outer revert**
  - [ ] At SRAA block entry (just before the outer SRAA loop), add: `int rk_outer_stall_count = 0;`.
  - [ ] Immediately after the existing `if (r.err_result < best_metric_seen) { ... W_best = X; }` block, append:
    ```cpp
    {
        double curr_quality_rk = r.err_result;
        if (curr_quality_rk > best_metric_seen * (1.0 + kSRAAOuterSlack)) {
            if (++rk_outer_stall_count >= kSRAAOuterStallWindow) {
                X = W_best;          // F_eval rebuilds derived state on next call
                rk_sraa.clear();
                rk_outer_stall_count = 0;
            }
        } else {
            rk_outer_stall_count = 0;
        }
    }
    ```
  - [ ] Also reset `rk_outer_stall_count = 0;` inside the existing improvement branch (to mirror the greenkhorn semantics: a successful improvement clears the stall counter). If the existing block is `if (r.err_result < best_metric_seen) { best_metric_seen = r.err_result; W_best = X; }`, extend it to `{ best_metric_seen = r.err_result; W_best = X; rk_outer_stall_count = 0; }`.
  - [ ] Confirm `kSRAAOuterSlack` and `kSRAAOuterStallWindow` are visible via the existing `sraa.hpp` include. If `raking.cpp` does not already include `sraa.hpp`, add `#include "sraa.hpp"` at the top in alphabetical order with the other project includes.

- [ ] **3.3 Compile gate**
  - [ ] `R CMD INSTALL --preclean .`. Zero new warnings.

- [ ] **3.4 Full test suite**
  - [ ] `Rscript -e 'devtools::test()' 2>&1 | tee benchmarks/sraa-m-t3-fulltest.log`.
  - [ ] All three SRAA tests must be GREEN. No regression elsewhere. If `T_sraa_raking_K9` still fails, halt with `SPEC_FAILURE`, log to `.wolf/buglog.json`, do **not** pivot.

- [ ] **3.5 Post-fix benchmark + comparison**
  - [ ] `Rscript benchmarks/stepstone_all_methods.R 2>&1 | tee benchmarks/sraa-m-baseline-post.log` (same script used in T1).
  - [ ] Persist post-fix per-K results to `benchmarks/sraa-m-baseline-post.csv` with the same schema as `benchmarks/sraa-m-baseline-pre.csv`.
  - [ ] Generate `benchmarks/sraa-m-baseline-diff.md` via:
    ```bash
    Rscript -e 'pre <- data.table::fread("benchmarks/sraa-m-baseline-pre.csv"); post <- data.table::fread("benchmarks/sraa-m-baseline-post.csv"); m <- merge(pre, post, by = c("method","K"), suffixes = c("_pre","_post")); m[, delta_pct := 100 * (max_errRp_post - max_errRp_pre) / max_errRp_pre]; data.table::fwrite(m, "benchmarks/sraa-m-baseline-diff.csv"); cat("# SRAA-m baseline diff\n\n", file = "benchmarks/sraa-m-baseline-diff.md"); knitr::kable(m) |> as.character() |> cat(sep = "\n", file = "benchmarks/sraa-m-baseline-diff.md", append = TRUE)'
    ```
  - [ ] Inspect the diff. Required outcomes (else `SPEC_FAILURE`, do not pivot):
    - SRAA-m at K=9: `max_errRp_post` within 2% of plain greenkhorn at K=9.
    - SRAA-raking at K=9: `max_errRp_post` within 2% of plain raking at K=9.
    - No K (3, 5, 9) regresses by more than the spec-stated tolerance vs the pre-fix value of plain greenkhorn / raking.
    - Wallclock for SRAA-m at K=9 is not >2x plain greenkhorn (sanity guardrail; the revert should be rare).
  - [ ] Trace every code path of the diff in writing in the commit body: outer-improve branch (counter resets), small-stall branch (counter increments but no revert), revert-fires branch, raking scalar-metric branch. Explicitly state outcomes for unoptimized (already-converged) paths.

- [ ] **3.6 Commit**
  - [ ] Stage exactly: `src/raking.cpp`, `benchmarks/sraa-m-t3-fulltest.log`, `benchmarks/sraa-m-baseline-post.log`, `benchmarks/sraa-m-baseline-post.csv`, `benchmarks/sraa-m-baseline-diff.csv`, `benchmarks/sraa-m-baseline-diff.md`.
  - [ ] Commit:
    ```bash
    git commit -m "$(cat <<'EOF'
    fix(sraa): T3 — outer revert in raking SRAA path + verify

    Mirrors the greenkhorn outer-revert pattern in raking.cpp using the
    scalar r.err_result metric. No S_flat rebuild — F_eval reconstructs
    derived state on the next call. Closes T_sraa_raking_K9.

    Post-fix stepstone benchmark (benchmarks/sraa-m-baseline-{post,diff}.{csv,md})
    confirms SRAA-m and SRAA-raking at K=9 land within 2% of their plain
    counterparts, eliminating the 35% K>=3 regression.

    Trace per branch: improve resets counter; small stall increments
    only; revert fires after 5 consecutive >10% excursions and clears
    AA history; raking scalar path uses W_best (F_eval rebuilds state).
    EOF
    )"
    ```
  - [ ] `bd update ${T3_BEAD_ID} --status closed --comment "raking outer revert landed; full suite green; baseline diff verifies K∈{3,5,9}."`

### Compile gate

`R CMD INSTALL --preclean .` succeeds, zero new warnings.

### Test gate

`devtools::test()` fully green; all three new SRAA tests pass; no regression elsewhere.

### Benchmark gate

`benchmarks/sraa-m-baseline-diff.md` shows: SRAA-m K=9 within 2% of plain greenkhorn K=9; SRAA-raking K=9 within 2% of plain raking K=9; no K regresses beyond spec tolerance; wallclock not >2x plain.

---

## Session Close (after Task 3)

- [ ] `git pull --rebase`
- [ ] `bd dolt push`
- [ ] `git push`
- [ ] `git status` — must show `up to date with origin`.
- [ ] If push fails, resolve and retry. Do not stop with work stranded locally.

PLAN_COMPLETE
