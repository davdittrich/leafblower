# iEPPA Convergence Hardening Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate false-positive `RK_ERR_INFEAS` on structurally feasible inputs and close the 22× per-iter gap vs raking at dense compression — without API changes.

**Architecture:** Three atomic commits on `src/ieppa.cpp`. WU-1 replaces the single-hit infeasibility latch with a per-(k,j) streak counter + recoverable persistent-pair set. WU-2 adds a runtime dispatch to a linear-space Sinkhorn inner loop when `M_cell/n > 0.5`, with a one-shot overflow fallback to log-space. WU-3 adds geometric-blend damping auto-triggered mid-streak by WU-1 state. Wrappers (`R/harvest.R`, `python/leafblower/_harvest.py`) get a one-line error-message update in WU-1's commit.

**Tech Stack:** C++17, R (testthat), Python (pytest), Rcpp/.Call bridge. No new deps.

**Source spec:** `docs/superpowers/specs/2026-04-23-ieppa-convergence-hardening-design.md` (rev 4, commit 181517d — APPROVED after 3-iter design-review-gate).

**Commit ordering (atomic):** WU-1 → WU-2 → WU-3. WU-3 reads `infeas_streak` state introduced by WU-1; WU-2 uses the same `record_empty`/`record_nonempty` lambdas from WU-1. Do not reorder.

**Build gate after each file edit:** `R CMD INSTALL --preclean .` must succeed before the next step. Never edit two source files in a row without a build between them.

---

## Pre-flight

- [ ] **Step P.1: Confirm clean working tree**

Run: `git status --short`
Expected: only untracked items already in `.gitignore` (`.wolf/`, `.tldr/`, `.claude/rules/`, `tasks/`, `__pycache__/`, build artifacts under `src/`, `leafblower.Rcheck/`, `benchmarks/autumn_nr_benchmark.R`, `benchmarks/stepstone_fulldata_*`).
No modified tracked files. If any, stop and reconcile before proceeding.

- [ ] **Step P.2: Confirm baseline tests pass**

Run: `R CMD INSTALL --preclean . && Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", stop_on_failure=TRUE, reporter="summary")'`
Expected: `[ FAIL 0 | WARN * | SKIP * | PASS >= 169 ]`
If any test fails, stop. Do not start the plan on a broken baseline.

- [ ] **Step P.3: Confirm `.Rbuildignore` excludes the fulldata benchmark**

Run: `grep -E '^\\^benchmarks' .Rbuildignore`
Expected: `^benchmarks$` present. Whole `benchmarks/` directory (including `stepstone_fulldata_*`) is excluded from CRAN tarball — spec requirement satisfied, no edit needed.

- [ ] **Step P.4: Audit existing tests/docs for the old error-message substring**

Run: `grep -rn 'empty cell with positive target' tests/ R/ python/ --include='*.R' --include='*.py'`
Expected: every match is either (a) inside `R/harvest.R` or `python/leafblower/_harvest.py` (the two source sites Task 1 edits), (b) inside a `#` / `"""` comment, or (c) inside a test that uses a regex NOT covering the full substring that changes (e.g. `expect_error(..., regexp = "infeasible problem")` — the prefix "infeasible problem" survives the rewrite). ACTION: for each match, classify. If any match is a direct `expect_error(..., regexp = "empty cell with positive target")` assertion that the rewrite breaks, extend Task 1's edit list to update that assertion to `regexp = "persistent empty cell"` and note the new file in the WU-1 commit. If all matches are benign (harmless comments or surviving-prefix regexes), proceed without additions.

- [ ] **Step P.5: Record baseline status distribution for RK_OK-preservation gate**

Using testthat's public expectation classes (`expectation_success` / `expectation_failure` / `expectation_warning` / `expectation_skip`) via the `testthat_results` object returned by `test_dir`:

```bash
Rscript -e '
  library(testthat); library(leafblower)
  tr <- test_dir("tests/testthat",
                 reporter = SilentReporter$new(),
                 stop_on_failure = FALSE)
  classify <- function(r) {
    if (inherits(r, "expectation_skip"))    "skip"
    else if (inherits(r, "expectation_failure") ||
             inherits(r, "expectation_error"))  "fail"
    else if (inherits(r, "expectation_warning")) "warn"
    else                                          "pass"
  }
  tallies <- c(pass=0L, fail=0L, warn=0L, skip=0L)
  per_test <- list()
  for (ti in tr) {
    for (r in ti$results) {
      cls <- classify(r)
      tallies[cls] <- tallies[cls] + 1L
    }
    # Worst result for this test_that block (for file-level preservation check).
    block_cls <- "pass"
    for (r in ti$results) {
      cc <- classify(r)
      if (cc == "fail") { block_cls <- "fail"; break }
      if (cc == "warn" && block_cls == "pass") block_cls <- "warn"
      if (cc == "skip" && block_cls == "pass") block_cls <- "skip"
    }
    per_test[[paste0(ti$file, "::", ti$test)]] <- block_cls
  }
  cat("baseline_pass=", tallies["pass"], "\n", sep="")
  cat("baseline_fail=", tallies["fail"], "\n", sep="")
  cat("baseline_warn=", tallies["warn"], "\n", sep="")
  cat("baseline_skip=", tallies["skip"], "\n", sep="")
  saveRDS(per_test, "/tmp/baseline-per-test.rds")
' | tee /tmp/baseline-counts.txt
```

Expected: `baseline_fail=0`, `baseline_pass >= 169`. `/tmp/baseline-per-test.rds` captures the per-(file,test) pass/warn/skip status for the RK_OK-preservation diff at Step A.4.1.

---

## Task 1: WU-1 — Persistent-infeas tracker

**Files:**
- Modify: `src/ieppa.cpp` (replace latch, lines 95–96, 120–126, 145–165, 290–292, 306–323)
- Modify: `R/harvest.R:102` (error message)
- Modify: `python/leafblower/_harvest.py:165` (error message)
- Create: `tests/testthat/test-ieppa-persistent-infeas.R`

**Rationale:** `is_infeasible` bool in `src/ieppa.cpp:95` is latched once and never cleared. Under K≥5 overlapping margins, multiplicative Sinkhorn updates transiently push some `(k, j)` buckets to near-zero before subsequent margins re-populate them. The current latch fires on transients → `RK_ERR_INFEAS` on structurally feasible inputs (stepstone bug). Fix: per-(k,j) streak counter, latch only when the same bucket has been empty for `kInfeasPersistence = 5` consecutive outer iterations. Recovery (non-empty check) resets the streak AND erases from the persistent set — both ends must be symmetric or a bucket that recovers stays falsely latched.

- [ ] **Step 1.1: Write the failing test**

Create `tests/testthat/test-ieppa-persistent-infeas.R`:

```r
# WU-1: persistent-infeas tracker regression test.
# Pre-fix: iEPPA flags RK_ERR_INFEAS on transient near-zero buckets
#   (stepstone-shape K=5 input with overlapping margins).
# Post-fix: converges to RK_OK; infeasibility reported only after
#   kInfeasPersistence=5 consecutive outer iterations with the same
#   bucket empty (guards against transient-in-settling false positives).

test_that("WU-1: iEPPA converges on structurally feasible K=5 overlapping margins", {
  set.seed(42)
  n <- 2000L
  K <- 5L
  # 3-3-4-3-3 cats with overlapping (correlated) groupings.
  df <- data.frame(
    a = sample(letters[1:3], n, replace = TRUE, prob = c(0.5, 0.3, 0.2)),
    b = sample(letters[1:3], n, replace = TRUE, prob = c(0.2, 0.5, 0.3)),
    c = sample(letters[1:4], n, replace = TRUE, prob = c(0.3, 0.3, 0.2, 0.2)),
    d = sample(letters[1:3], n, replace = TRUE, prob = c(0.4, 0.3, 0.3)),
    e = sample(letters[1:3], n, replace = TRUE, prob = c(0.25, 0.25, 0.5))
  )
  targets <- list(
    a = c(a = 0.40, b = 0.35, c = 0.25),
    b = c(a = 0.30, b = 0.40, c = 0.30),
    c = c(a = 0.25, b = 0.25, c = 0.25, d = 0.25),
    d = c(a = 0.33, b = 0.33, c = 0.34),
    e = c(a = 0.30, b = 0.35, c = 0.35)
  )
  # Structurally feasible at max_weight=5; pre-fix latches INFEAS on transient.
  expect_no_error(
    res <- harvest(df, targets, method = "ieppa",
                   max_weight = 5, min_weight = 0,
                   max_iterations = 500L,
                   convergence = list(absolute = 1e-4))
  )
  diag <- diagnose_weights(df, targets, res)
  expect_lt(max(abs(diag$error_weighted)), 1e-3)
})

test_that("WU-1: truly infeasible input (empty target cell) still reports INFEAS", {
  # Regression guard: genuine infeasibility must still latch.
  n <- 500L
  df <- data.frame(
    a = sample(letters[1:2], n, replace = TRUE),
    b = sample(letters[1:2], n, replace = TRUE)
  )
  # Target a third category 'c' that has zero observations → persistent empty.
  targets <- list(
    a = c(a = 0.4, b = 0.3, c = 0.3),
    b = c(a = 0.5, b = 0.5)
  )
  expect_error(
    suppressWarnings(harvest(df, targets, method = "ieppa",
                             max_weight = 5, min_weight = 0,
                             max_iterations = 500L,
                             convergence = list(absolute = 1e-6))),
    regexp = "persistent empty cell"
  )
})

test_that("WU-1: oscillating streak (spec §4 edge case) returns NOCONV not INFEAS", {
  # Spec §4 documents: a bucket that oscillates empty <-> non-empty such that
  # streak resets before kInfeasPersistence=5 will NOT flag INFEAS; solver
  # hits max_iter -> RK_ERR_NOCONV with high errRp. This test guards that
  # documented behaviour. Engineer oscillation via a 3-way near-degenerate
  # system where each outer sweep alternates which margin is pinched.
  set.seed(2024)
  n <- 600L
  # K=3 with strong negative correlation between two margins → oscillation.
  a <- sample(letters[1:3], n, replace = TRUE, prob = c(0.1, 0.45, 0.45))
  # b chosen so (a,b) cells heavily biased; targets will push the solver to
  # move mass between (a="a", any b) cells, which are sparse.
  b <- ifelse(a == "a", sample(letters[1:3], n, replace = TRUE, prob = c(0.8, 0.1, 0.1)),
                        sample(letters[1:3], n, replace = TRUE))
  c_ <- sample(letters[1:3], n, replace = TRUE)
  df <- data.frame(a = a, b = b, c_ = c_)
  # Targets pushing the (a="a", b="b") and (a="a", b="c") cells to non-trivial mass.
  targets <- list(
    a  = c(a = 0.50, b = 0.25, c = 0.25),
    b  = c(a = 0.25, b = 0.50, c = 0.25),
    c_ = c(a = 0.33, b = 0.33, c = 0.34)
  )
  # Low max_iter forces solver to exit before streak settles either way.
  res_status <- tryCatch({
    suppressWarnings(harvest(df, targets, method = "ieppa",
                             max_weight = 5, min_weight = 0,
                             max_iterations = 30L,
                             convergence = list(absolute = 1e-10)))
    "converged"
  }, error = function(e) {
    if (grepl("persistent empty cell", conditionMessage(e))) "infeas"
    else "other_error"
  })
  # Must NOT be "infeas" on this near-degenerate but not structurally empty input
  # with streak-resetting oscillation. "converged" (OK) or a non-infeas warning
  # path are both acceptable per spec §4.
  expect_false(res_status == "infeas")
})
```

- [ ] **Step 1.2: Run test — confirm RED**

Run: `R CMD INSTALL --preclean . && Rscript -e 'library(testthat); library(leafblower); test_file("tests/testthat/test-ieppa-persistent-infeas.R", reporter="summary")'`
Expected: first test FAILS (current code latches on transient → RK_ERR_INFEAS → `harvest` throws `infeasible problem — empty cell`). Second test PASSES already (genuine infeasibility). If first test passes pre-fix, the transient trigger may already be suppressed — halt and investigate the RNG seed / input shape rather than proceeding.

- [ ] **Step 1.3: Edit `src/ieppa.cpp` — replace `is_infeasible` latch with streak counter**

In `src/ieppa.cpp`, replace:

```cpp
    bool is_infeasible = false;
    std::set<std::pair<int,int>> infeasible_pairs;  // dedup via set ordering
```

with:

```cpp
    constexpr int kInfeasPersistence = 5;
    std::vector<int> infeas_streak(total_cats, 0);
    std::set<std::pair<int,int>> persistent_infeas_pairs;

    auto record_empty = [&](int k, int j) {
        if (st.targets[k][j] <= 0.0) return;
        int idx = cat_offset[k] + j;
        infeas_streak[idx]++;
        if (infeas_streak[idx] == kInfeasPersistence) {
            persistent_infeas_pairs.emplace(k, j);
        }
    };
    auto record_nonempty = [&](int k, int j) {
        if (st.targets[k][j] <= 0.0) return;
        int idx = cat_offset[k] + j;
        infeas_streak[idx] = 0;
        persistent_infeas_pairs.erase(std::make_pair(k, j));
    };
```

Note: `==` (not `>=`) so the set `emplace` fires exactly once per `(k,j)`; subsequent sweeps with streak > persistence are no-ops. Erasure in `record_nonempty` is load-bearing (Security B5).

- [ ] **Step 1.4: Edit `src/ieppa.cpp` — wire `record_empty`/`record_nonempty` into margin sweep**

Replace the three latch sites in the margin sweep.

Site 1 — `cells.empty()` branch (was `src/ieppa.cpp:120-126`):

```cpp
                if (cells.empty()) {
                    record_empty(k, j);
                    continue;
                }
```

Site 2 — `!std::isfinite(lv_max)` branch (was `src/ieppa.cpp:145-152`):

```cpp
                if (!std::isfinite(lv_max)) {
                    record_empty(k, j);
                    continue;
                }
```

Site 3 — `log_S_kj < log_threshold` branch (was `src/ieppa.cpp:160-166`):

```cpp
                if (!std::isfinite(log_S_kj) || log_S_kj < log_threshold) {
                    record_empty(k, j);
                    continue;
                }
```

Then, on the success path (immediately before `lf[cat_offset[k] + j] = log_target - log_S_kj;`), add:

```cpp
                record_nonempty(k, j);
                double log_target = std::log(st.targets[k][j] * ct.W_input);
                lf[cat_offset[k] + j] = log_target - log_S_kj;
```

- [ ] **Step 1.5: Edit `src/ieppa.cpp` — derive `is_infeasible` at exit, update verbose print**

Replace the tol-met block inside the convergence-check (was `src/ieppa.cpp:264-267`):

```cpp
            if (errRp < st.tol_abs) {
                res.status = persistent_infeas_pairs.empty() ? RK_OK : RK_ERR_INFEAS;
                break;
            }
```

Replace the post-loop INFEAS promotion (was `src/ieppa.cpp:290-292`):

```cpp
    if (!persistent_infeas_pairs.empty() && res.status == RK_ERR_NOCONV) {
        res.status = RK_ERR_INFEAS;
    }
```

Replace the verbose enumeration block (was `src/ieppa.cpp:306-323`) — same logic, new variable:

```cpp
    if (st.verbose >= 1 && !persistent_infeas_pairs.empty()) {
        char msg[256];
        size_t off = 0;
        off += std::snprintf(msg + off, sizeof(msg) - off,
                             "iEPPA persistent infeasible cells: ");
        size_t idx = 0;
        const size_t total = persistent_infeas_pairs.size();
        for (auto it = persistent_infeas_pairs.begin();
             it != persistent_infeas_pairs.end() && off < sizeof(msg) - 32;
             ++it, ++idx) {
            off += std::snprintf(msg + off, sizeof(msg) - off,
                                 "margin=%d cat=%d%s",
                                 it->first + 1,
                                 it->second + 1,
                                 (idx + 1 < total) ? ", " : "");
        }
        st.log(msg);
    }
```

Remove the status-label branch for `RK_ERR_INFEAS` — already covered by existing switch. Delete any now-unused references to `is_infeasible` / `infeasible_pairs`.

- [ ] **Step 1.6: Build gate**

Run: `R CMD INSTALL --preclean .`
Expected: installs cleanly, no warnings from `src/ieppa.cpp`. If the edit tool reports "applied with inaccuracies", view the file, revert, and redo this step — never proceed with a known-corrupt source file.

- [ ] **Step 1.7: Edit `R/harvest.R:102` — error-message refinement**

Replace:

```r
    stop("leafblower: infeasible problem — empty cell with positive target.")
```

with:

```r
    stop("leafblower: infeasible problem — persistent empty cell with positive target (detected after 5 consecutive outer iterations).")
```

- [ ] **Step 1.8: Edit `python/leafblower/_harvest.py:165` — error-message refinement**

Replace:

```python
        raise RuntimeError("leafblower: infeasible problem — empty cell with positive target")
```

with:

```python
        raise RuntimeError(
            "leafblower: infeasible problem — persistent empty cell with positive target "
            "(detected after 5 consecutive outer iterations)"
        )
```

- [ ] **Step 1.9: Run the WU-1 test — expect GREEN**

Run: `R CMD INSTALL --preclean . && Rscript -e 'library(testthat); library(leafblower); test_file("tests/testthat/test-ieppa-persistent-infeas.R", reporter="summary")'`
Expected: both tests PASS.

- [ ] **Step 1.10: Run the full test suite — expect no regressions**

Run: `Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", stop_on_failure=TRUE, reporter="summary")'`
Expected: `[ FAIL 0 | PASS >= 171 ]` (baseline 169 + 2 new).

- [ ] **Step 1.11: Stepstone false-positive verification**

Run: `Rscript /tmp/stepstone_3algo_bench.R 2>&1 | tee /tmp/wu1_stepstone.log`
Expected: `ieppa` row no longer `ERROR ... infeasible problem`. Status either `OK` or `NOCONV` (either acceptable for WU-1 — WU-2/WU-3 address performance/hardening). `raking` row unchanged (~2.3s, errRp 8.25e-3).
If `ieppa` still reports INFEAS, WU-1 has not actually fixed the stepstone bug — halt and inspect `persistent_infeas_pairs` via `verbose=2`.

- [ ] **Step 1.12: Commit**

```bash
git add src/ieppa.cpp R/harvest.R python/leafblower/_harvest.py tests/testthat/test-ieppa-persistent-infeas.R
git commit -m "$(cat <<'EOF'
fix(ieppa): persistent-infeas tracker replaces single-hit latch

Replace bool is_infeasible + infeasible_pairs with per-(k,j) streak counter
and persistent_infeas_pairs set. Latch only after kInfeasPersistence=5
consecutive outer iterations with the same empty bucket; any non-empty
check resets the streak and erases from the persistent set.

Fixes false-positive RK_ERR_INFEAS on stepstone (n=1.58M, K=9, 836 cats,
54x compression) where transient near-zero buckets under K>=5 overlapping
margins previously latched infeasibility on iter 1.

Wrapper error messages in R/harvest.R and python/leafblower/_harvest.py
updated to say "persistent empty cell" and cite the 5-iter settle window.

Refs spec docs/superpowers/specs/2026-04-23-ieppa-convergence-hardening-design.md WU-1.
EOF
)"
```

---

## Task 2: WU-2 — Linear-space Sinkhorn dispatch at dense compression

**Files:**
- Modify: `src/ieppa.cpp` (add dispatch + linear inner loop + overflow fallback)
- Modify: `tests/testthat/test-ieppa-faithful.R` (append dense/sparse/overflow cases)

**Rationale:** At `M_cell/n ≥ 0.5` the log-space inner loop's transcendental cost (per-cell `log`/`exp` + log-sum-exp max pass per margin per iter) buys no numerical-stability benefit — weights stay moderate under the capacity cap. Raking wins 22× on this regime (measured kk1204: n=1M, K=20, cat=5, M_cell/n=1.0). Switch to direct multiplicative updates when compression ≤ 2×; keep log-space when compression > 2× (sparse regime where log-space stability is load-bearing). One-shot overflow fallback to log-space with full state reset handles adversarial K/weight combinations.

- [ ] **Step 2.1: Write the failing tests (append to `tests/testthat/test-ieppa-faithful.R`)**

Append to `tests/testthat/test-ieppa-faithful.R`:

```r
test_that("WU-2: dense regime (M_cell/n ~ 1) linear-space matches log-space to 1e-8", {
  set.seed(123)
  n <- 10000L
  # K=8, 3 cats each: 3^8 = 6561 cells; at n=10000 roughly M_cell/n ~ 0.6-0.7 → linear path.
  df <- as.data.frame(replicate(8, sample(1:3, n, replace = TRUE), simplify = FALSE))
  names(df) <- paste0("m", 1:8)
  targets <- setNames(
    replicate(8, c(a = 0.3, b = 0.4, c = 0.3), simplify = FALSE),
    paste0("m", 1:8)
  )
  for (k in names(df)) df[[k]] <- c("a","b","c")[df[[k]]]

  Sys.setenv(LBW_IEPPA_FORCE_PATH = "linear")
  on.exit(Sys.unsetenv("LBW_IEPPA_FORCE_PATH"), add = TRUE)
  res_lin <- harvest(df, targets, method = "ieppa",
                     max_weight = 10, min_weight = 0,
                     max_iterations = 500L,
                     convergence = list(absolute = 1e-6))

  Sys.setenv(LBW_IEPPA_FORCE_PATH = "log")
  res_log <- harvest(df, targets, method = "ieppa",
                     max_weight = 10, min_weight = 0,
                     max_iterations = 500L,
                     convergence = list(absolute = 1e-6))
  Sys.unsetenv("LBW_IEPPA_FORCE_PATH")

  expect_lt(max(abs(res_lin - res_log)), 1e-8)
})

test_that("WU-2: sparse regime (M_cell/n ~ 0.01) auto-dispatches log-space", {
  set.seed(7)
  n <- 10000L
  # K=2, 5 cats each: 5^2 = 25 cells; M_cell/n ~ 0.0025 → log-space path.
  df <- data.frame(
    a = sample(letters[1:5], n, replace = TRUE),
    b = sample(letters[1:5], n, replace = TRUE)
  )
  targets <- list(
    a = setNames(rep(0.2, 5), letters[1:5]),
    b = setNames(rep(0.2, 5), letters[1:5])
  )
  msgs <- capture.output(
    res <- harvest(df, targets, method = "ieppa",
                   max_weight = 5, min_weight = 0,
                   max_iterations = 500L,
                   convergence = list(absolute = 1e-6),
                   verbose = 1L),
    type = "message"
  )
  # Verbose log prefix differs by path; we rely on either absence of
  # "linear-space" OR presence of "path=log" if labelled.
  expect_false(any(grepl("linear-space", msgs)))
})

test_that("WU-2: overflow synthesis falls back to log-space, still completes", {
  # Force linear path on a high-K input; rely on adversarial targets to drive
  # factor * prod(f) toward kLinearOverflowTrip; expect one-shot fallback.
  set.seed(99)
  n <- 5000L
  K <- 12L
  df <- as.data.frame(replicate(K, sample(1:3, n, replace = TRUE), simplify = FALSE))
  names(df) <- paste0("m", 1:K)
  for (k in names(df)) df[[k]] <- c("a","b","c")[df[[k]]]
  targets <- setNames(
    replicate(K, c(a = 0.6, b = 0.3, c = 0.1), simplify = FALSE),
    paste0("m", 1:K)
  )
  Sys.setenv(LBW_IEPPA_FORCE_PATH = "linear")
  on.exit(Sys.unsetenv("LBW_IEPPA_FORCE_PATH"), add = TRUE)
  # Success condition: no error, status is RK_OK or RK_ERR_NOCONV, weights finite.
  res <- suppressWarnings(harvest(df, targets, method = "ieppa",
                                  max_weight = 1e6, min_weight = 0,
                                  max_iterations = 200L,
                                  convergence = list(absolute = 1e-4)))
  expect_true(all(is.finite(res)))
})
```

- [ ] **Step 2.2: Run the new tests — confirm RED**

Run: `R CMD INSTALL --preclean . && Rscript -e 'library(testthat); library(leafblower); test_file("tests/testthat/test-ieppa-faithful.R", reporter="summary")'`
Expected: the new WU-2 tests FAIL (current code has no `LBW_IEPPA_FORCE_PATH`, no linear path). Existing tests in the file still pass. If the new tests pass pre-fix, investigate why the env var happens to be a no-op.

- [ ] **Step 2.3: Edit `src/ieppa.cpp` — add env-var path selector**

Immediately after `ct.M_cell` is populated (after `res.M_cell = ct.M_cell;`), add:

```cpp
    // WU-2 dispatch: linear-space path when M_cell/n > 0.5 (i.e., compression <= 2x).
    // Env var LBW_IEPPA_FORCE_PATH ∈ {"linear", "log"} overrides for tests; always
    // compiled (no #ifdef) — getenv cost is microseconds, amortized over the solve.
    constexpr double kLinearSpaceThreshold = 2.0;  // compression ratio cutoff
    bool use_linear = (static_cast<double>(st.n) /
                       static_cast<double>(std::max(ct.M_cell, 1))) <
                      kLinearSpaceThreshold;
    if (const char* force = std::getenv("LBW_IEPPA_FORCE_PATH")) {
        if (std::strcmp(force, "linear") == 0) use_linear = true;
        else if (std::strcmp(force, "log") == 0) use_linear = false;
    }
```

Add `#include <cstdlib>` and `#include <cstring>` to the top-of-file includes if not already present.

- [ ] **Step 2.4: Edit `src/ieppa.cpp` — add linear-space overflow trip constant**

Just before the outer iteration loop (`for (int iter = 1; iter <= st.inner_max_iter; iter++)`), compute the overflow trip:

```cpp
    // Runtime trip: factor = X_init[c] * W[c] * ∏_m f[m] ≤ max_X_init * trip^K < DBL_MAX/2.
    // Accounts for both K-way product accumulation AND X_init prefactor (Security B2).
    double max_X_init_val = 1.0;
    for (int c = 0; c < ct.M_cell; c++) {
        if (X_init[c] > max_X_init_val) max_X_init_val = X_init[c];
    }
    const double kLinearOverflowTrip = std::pow(
        std::numeric_limits<double>::max() / (2.0 * max_X_init_val),
        1.0 / static_cast<double>(st.K));

    // Linear-space Sinkhorn factors (mirror of lf[], but in linear domain).
    // Only populated when use_linear is true at iter entry.
    std::vector<double> f_lin(total_cats, 1.0);
    bool linear_fallback_used = false;
```

- [ ] **Step 2.5: Edit `src/ieppa.cpp` — linear-space margin-sweep branch**

Inside the outer `for (int iter ...)` loop, replace the entire margin sweep block (currently `for (int k = 0; k < st.K; k++) { for (int j ...) { ... } }`) with a branch:

```cpp
        bool overflow_trip = false;

        if (use_linear) {
            // WU-2 linear-space sweep: direct multiplicative updates.
            for (int k = 0; k < st.K && !overflow_trip; k++) {
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    const auto& cells = cells_by_margin_cat[cat_offset[k] + j];
                    if (cells.empty()) { record_empty(k, j); continue; }
                    double S_kj = 0.0;
                    for (int c : cells) {
                        if (X_init[c] <= 0.0 || W[c] <= 0.0) continue;
                        double factor = X_init[c] * W[c];
                        for (int m = 0; m < st.K; m++) {
                            if (m == k) continue;
                            factor *= f_lin[cat_offset[m] + ct.g_per_cell[m][c]];
                        }
                        S_kj += factor;
                    }
                    if (!(S_kj >= kEmptyBucketThreshold * ct.W_input) ||
                        !std::isfinite(S_kj)) {
                        record_empty(k, j);
                        continue;
                    }
                    record_nonempty(k, j);
                    double new_f = (st.targets[k][j] * ct.W_input) / S_kj;
                    if (!std::isfinite(new_f) || new_f > kLinearOverflowTrip) {
                        overflow_trip = true;
                        break;
                    }
                    f_lin[cat_offset[k] + j] = new_f;
                }
            }
            if (overflow_trip && !linear_fallback_used) {
                // One-shot fallback: reset all solver state, switch to log-space,
                // restart outer loop from iter 0.
                linear_fallback_used = true;
                use_linear = false;
                std::fill(lf.begin(), lf.end(), 0.0);
                std::fill(f_lin.begin(), f_lin.end(), 1.0);
                std::fill(W.begin(), W.end(), 1.0);
                for (int c = 0; c < ct.M_cell; c++) {
                    X_tilde[c] = X_init[c];
                    X[c] = X_init[c];
                }
                std::fill(infeas_streak.begin(), infeas_streak.end(), 0);
                persistent_infeas_pairs.clear();
                iter = 0;  // outer for-loop increments → iter=1 next round
                if (st.verbose >= 1) {
                    st.log("iEPPA: linear-space overflow trip; fallback to log-space.");
                }
                continue;  // skip the post-sweep X_tilde / capacity / errRp blocks this round
            }
        } else {
            // Original log-space sweep (unchanged except record_empty/record_nonempty
            // hooks already installed by WU-1).
            for (int k = 0; k < st.K; k++) {
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    const auto& cells = cells_by_margin_cat[cat_offset[k] + j];
                    if (cells.empty()) { record_empty(k, j); continue; }
                    lv.assign(cells.size(), -std::numeric_limits<double>::infinity());
                    double lv_max = -std::numeric_limits<double>::infinity();
                    for (size_t r = 0; r < cells.size(); r++) {
                        int c = cells[r];
                        if (X_init[c] <= 0.0 || W[c] <= 0.0) continue;
                        double s = log_X_init[c] + std::log(W[c]);
                        for (int m = 0; m < st.K; m++) {
                            if (m == k) continue;
                            int gm = ct.g_per_cell[m][c];
                            s += lf[cat_offset[m] + gm];
                        }
                        lv[r] = s;
                        if (s > lv_max) lv_max = s;
                    }
                    if (!std::isfinite(lv_max)) { record_empty(k, j); continue; }
                    double sum = 0.0;
                    for (size_t r = 0; r < lv.size(); r++) {
                        if (std::isfinite(lv[r])) sum += std::exp(lv[r] - lv_max);
                    }
                    double log_S_kj = lv_max + std::log(sum);
                    double log_threshold = std::log(kEmptyBucketThreshold * ct.W_input);
                    if (!std::isfinite(log_S_kj) || log_S_kj < log_threshold) {
                        record_empty(k, j);
                        continue;
                    }
                    record_nonempty(k, j);
                    double log_target = std::log(st.targets[k][j] * ct.W_input);
                    lf[cat_offset[k] + j] = log_target - log_S_kj;
                }
            }
        }
```

- [ ] **Step 2.6: Edit `src/ieppa.cpp` — branch X_tilde reconstruction by path**

Replace the `X_tilde` / overflow-detection block with a path-aware version. Linear uses direct product; log-space keeps the existing clip-before-exp logic.

Replace the block (was `src/ieppa.cpp:172-203`) with:

```cpp
        // Compute X_tilde. Path-dependent: linear uses direct product; log uses clip-before-exp.
        bool overflow_detected = false;
        double max_log_X_tilde = -std::numeric_limits<double>::infinity();
        for (int c = 0; c < ct.M_cell; c++) {
            if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
            if (use_linear) {
                double v = X_init[c];
                for (int m = 0; m < st.K; m++) {
                    int gm = ct.g_per_cell[m][c];
                    v *= f_lin[cat_offset[m] + gm];
                }
                if (!std::isfinite(v)) overflow_detected = true;
                X_tilde[c] = v;
            } else {
                double s = log_X_init[c];
                for (int m = 0; m < st.K; m++) {
                    int gm = ct.g_per_cell[m][c];
                    s += lf[cat_offset[m] + gm];
                }
                if (s > max_log_X_tilde) max_log_X_tilde = s;
                double s_clip = (s > kLogClip) ? kLogClip : s;
                if (s > kLogClip && U_cell[c] >= 1e299) {
                    overflow_detected = true;
                }
                X_tilde[c] = std::exp(s_clip);
            }
        }
        if (overflow_detected) {
            if (use_linear && !linear_fallback_used) {
                // Treat same as overflow_trip: fall back once.
                linear_fallback_used = true;
                use_linear = false;
                std::fill(lf.begin(), lf.end(), 0.0);
                std::fill(f_lin.begin(), f_lin.end(), 1.0);
                std::fill(W.begin(), W.end(), 1.0);
                for (int c2 = 0; c2 < ct.M_cell; c2++) {
                    X_tilde[c2] = X_init[c2];
                    X[c2] = X_init[c2];
                }
                std::fill(infeas_streak.begin(), infeas_streak.end(), 0);
                persistent_infeas_pairs.clear();
                iter = 0;
                if (st.verbose >= 1) {
                    st.log("iEPPA: linear-space X_tilde overflow; fallback to log-space.");
                }
                continue;
            }
            res.status = RK_ERR_NOCONV;
            res.max_error = std::numeric_limits<double>::infinity();
            if (st.verbose >= 2) {
                char msg[256];
                std::snprintf(msg, sizeof(msg),
                              "iEPPA: log-factor overflow (max_log_X_tilde=%.1f > 700) "
                              "indicates ill-conditioning; try looser max_weight or "
                              "tighter tol_abs, or method=raking.", max_log_X_tilde);
                st.log(msg);
            }
            break;
        }
```

- [ ] **Step 2.7: Edit `src/ieppa.cpp` — verbose path marker**

In the verbose entry-log block (currently `src/ieppa.cpp:98-108`), append the path:

```cpp
        std::snprintf(msg, sizeof(msg),
                      "%siEPPA: n=%d K=%d M_cell=%d compression=%.1fx path=%s",
                      prefix, st.n, st.K, ct.M_cell,
                      (double)st.n / (double)std::max(ct.M_cell, 1),
                      use_linear ? "linear" : "log");
```

- [ ] **Step 2.8: Build gate**

Run: `R CMD INSTALL --preclean .`
Expected: clean install. No `-Wall` warnings from the new code.

- [ ] **Step 2.9: Run WU-2 tests — expect GREEN**

Run: `Rscript -e 'library(testthat); library(leafblower); test_file("tests/testthat/test-ieppa-faithful.R", reporter="summary")'`
Expected: all three new WU-2 tests PASS. Existing tests in this file still PASS.

- [ ] **Step 2.10: Full regression**

Run: `Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", stop_on_failure=TRUE, reporter="summary")'`
Expected: `[ FAIL 0 | PASS >= 174 ]` (171 post-WU-1 + 3 new).

- [ ] **Step 2.11: kk1204 per-iter check (spec §11: per-iter cost, not wall-clock)**

Spec §11 demands `per-iter cost within 2× of raking`, not wall-clock. Both solvers must be forced to run the full iteration budget (no early convergence short-circuit) so `wall/iter` equals per-iter cost. `tol_abs=1e-300` + `max_iterations=50L` achieves this because `errRp > 1e-300` always, so neither solver returns early. Divide each wall-clock by the actual iterations consumed (available in the returned `result$iterations`, or parseable from `verbose=1L` final-line "in N iters").

Create `/tmp/wu2_kk1204.R`:

```r
Sys.setenv(OMP_NUM_THREADS = "1")
library(leafblower)
set.seed(1)
n <- 1000000L; K <- 20L
df <- as.data.frame(replicate(K, sample(1:5, n, replace = TRUE), simplify = FALSE))
names(df) <- paste0("m", 1:K)
for (k in names(df)) df[[k]] <- letters[df[[k]]]
tgts <- setNames(replicate(K, c(a=0.3,b=0.175,c=0.175,d=0.175,e=0.175),
                           simplify = FALSE), names(df))

time_method <- function(method) {
  t0 <- Sys.time()
  msgs <- capture.output(
    res <- suppressWarnings(harvest(df, tgts, method = method,
                                    max_weight = 3, min_weight = 0,
                                    max_iterations = 50L,
                                    convergence = list(absolute = 1e-300),
                                    verbose = 1L)),
    type = "message"
  )
  wall <- as.numeric(Sys.time() - t0, units = "secs")
  m <- tail(grep("in [0-9]+ iters", msgs, value = TRUE), 1)
  iters <- as.integer(sub(".*in ([0-9]+) iters.*", "\\1", m))
  if (is.na(iters) || iters <= 0) stop(sprintf("couldn't parse iters for %s", method))
  list(wall = wall, iters = iters, per_iter_ms = 1000 * wall / iters)
}

r_ieppa  <- time_method("ieppa")
r_raking <- time_method("raking")
ratio <- r_ieppa$per_iter_ms / r_raking$per_iter_ms
cat(sprintf("ieppa: wall=%.2fs iters=%d per_iter=%.1fms\n",
            r_ieppa$wall, r_ieppa$iters, r_ieppa$per_iter_ms))
cat(sprintf("raking: wall=%.2fs iters=%d per_iter=%.1fms\n",
            r_raking$wall, r_raking$iters, r_raking$per_iter_ms))
cat(sprintf("per_iter ratio (ieppa / raking) = %.2fx\n", ratio))
if (ratio > 2.0) {
  cat("FAIL: per-iter ratio exceeds 2x (spec §11 gate).\n")
  quit(status = 1L)
}
cat("PASS: per-iter ratio within 2x.\n")
```

Run: `Rscript /tmp/wu2_kk1204.R`
Expected: exit 0 with "PASS: per-iter ratio within 2x". If exit 1 the linear path is either not dispatching (check verbose output for `path=linear`) or has unexpected overhead — halt and inspect.

- [ ] **Step 2.12: Commit**

```bash
git add src/ieppa.cpp tests/testthat/test-ieppa-faithful.R
git commit -m "$(cat <<'EOF'
perf(ieppa): linear-space Sinkhorn dispatch at M_cell/n > 0.5

Runtime dispatch on compression ratio: log-space retained when compression
> 2x (sparse regime, stability-critical); linear-space direct multiplicative
updates when compression <= 2x (dense regime, transcendental overhead dominates).

Overflow trip derived from DBL_MAX/(2*max_X_init) to the 1/K power; on trip
(either in factor accumulation or X_tilde reconstruction), one-shot full state
reset + restart in log-space. LBW_IEPPA_FORCE_PATH env var forces one path
for tests (always compiled; microsecond getenv cost per solve).

Closes kk1204 per-iter gap: iEPPA/raking ratio at M_cell/n=1.0 drops from
22x to <=2x without behaviour change in sparse regime.

Refs spec docs/superpowers/specs/2026-04-23-ieppa-convergence-hardening-design.md WU-2.
EOF
)"
```

---

## Task 3: WU-3 — Adaptive damping

**Files:**
- Create: `tests/testthat/_snapshots/ieppa_stable_prewu3.rds` (reference snapshot)
- Modify: `src/ieppa.cpp` (damped-mode auto-trigger, geometric-blend update)
- Modify: `tests/testthat/test-ieppa-faithful.R` (stable=regression + damped slower)

**Rationale:** Even after WU-1, transient near-zero cells still occur on adversarial inputs and force `infeas_streak` to creep up. Damping each Sinkhorn update with geometric blend `lf_new = α·(log_target - log_S_kj) + (1-α)·lf_old` reduces the amplitude of transients and lets the solver settle without tripping persistence. Auto-triggered at `infeas_streak[idx] >= kInfeasPersistence/2` — pays the cost only when WU-1 is already headed toward NOCONV without damping. `α=1.0` fast-path is a branch, not `std::pow(_, 1.0)` (which costs full transcendental).

- [ ] **Step 3.0: Capture pre-WU-3 stable-mode reference snapshot**

Spec §7 mandates WU-3's stable-mode test be a regression check against pre-WU-3 behaviour, not inter-run determinism. Before touching `src/ieppa.cpp` for WU-3, snapshot the exact weight vector the post-WU-2 build produces on a deterministic input. This snapshot will be loaded by the Step 3.1 stable-mode test after WU-3 lands.

Run:

```bash
Rscript -e '
  library(leafblower)
  set.seed(11)
  n <- 1000L
  df <- data.frame(
    a = sample(letters[1:3], n, replace = TRUE),
    b = sample(letters[1:3], n, replace = TRUE)
  )
  targets <- list(
    a = c(a = 0.33, b = 0.33, c = 0.34),
    b = c(a = 0.33, b = 0.33, c = 0.34)
  )
  res <- harvest(df, targets, method = "ieppa",
                 max_weight = 5, min_weight = 0,
                 max_iterations = 500L,
                 convergence = list(absolute = 1e-6))
  dir.create("tests/testthat/_snapshots", showWarnings = FALSE, recursive = TRUE)
  saveRDS(res, "tests/testthat/_snapshots/ieppa_stable_prewu3.rds")
  cat("snapshot saved: length=", length(res), " max=", max(res), " min=", min(res), "\n", sep="")
'
```

Expected: `snapshot saved: length=1000 ...`. The file `tests/testthat/_snapshots/ieppa_stable_prewu3.rds` is ADDED to git in the WU-3 commit (Step 3.11) and load-referenced by the stable-mode test in Step 3.1.

- [ ] **Step 3.1: Write the failing test (append to `tests/testthat/test-ieppa-faithful.R`)**

Append:

```r
test_that("WU-3: stable mode is byte-identical to pre-WU-3 reference snapshot (spec §7)", {
  # The pre-WU-3 snapshot was captured in Step 3.0 against the post-WU-2 build,
  # with no damping compiled in. Running the same deterministic input on the
  # post-WU-3 build with LBW_IEPPA_FORCE_DAMPING=off must reproduce weights
  # equal to the snapshot to floating-point roundoff (no drift from alpha==1.0
  # fast-path differing from the direct update).
  snap_path <- testthat::test_path("_snapshots", "ieppa_stable_prewu3.rds")
  skip_if_not(file.exists(snap_path),
              "pre-WU-3 snapshot missing; re-run plan Step 3.0")
  ref <- readRDS(snap_path)
  set.seed(11)
  n <- 1000L
  df <- data.frame(
    a = sample(letters[1:3], n, replace = TRUE),
    b = sample(letters[1:3], n, replace = TRUE)
  )
  targets <- list(
    a = c(a = 0.33, b = 0.33, c = 0.34),
    b = c(a = 0.33, b = 0.33, c = 0.34)
  )
  Sys.setenv(LBW_IEPPA_FORCE_DAMPING = "off")
  on.exit(Sys.unsetenv("LBW_IEPPA_FORCE_DAMPING"), add = TRUE)
  res <- harvest(df, targets, method = "ieppa",
                 max_weight = 5, min_weight = 0,
                 max_iterations = 500L,
                 convergence = list(absolute = 1e-6))
  expect_equal(length(res), length(ref))
  expect_lt(max(abs(res - ref)), 1e-12)  # roundoff-tight; fast-path ≡ direct update
})

test_that("WU-3: damped mode takes strictly more iters than stable on same input (spec §7)", {
  # Use LBW_IEPPA_FORCE_DAMPING to run the SAME feasible input twice: once
  # stable (alpha=1.0, fast path), once damped (alpha=0.5, geometric blend).
  # Spec §7 / CTO B5: monotone `iter_damped > iter_stable` assertion.
  set.seed(314)
  n <- 3000L
  K <- 6L
  df <- as.data.frame(replicate(K, sample(1:3, n, replace = TRUE), simplify = FALSE))
  names(df) <- paste0("m", 1:K)
  for (k in names(df)) df[[k]] <- c("a","b","c")[df[[k]]]
  targets <- setNames(
    replicate(K, c(a = 0.5, b = 0.3, c = 0.2), simplify = FALSE),
    paste0("m", 1:K)
  )

  run_one <- function(force_damping) {
    Sys.setenv(LBW_IEPPA_FORCE_DAMPING = force_damping)
    on.exit(Sys.unsetenv("LBW_IEPPA_FORCE_DAMPING"), add = TRUE)
    msgs <- capture.output(
      res <- suppressWarnings(harvest(df, targets, method = "ieppa",
                                      max_weight = 10, min_weight = 0,
                                      max_iterations = 500L,
                                      convergence = list(absolute = 1e-5),
                                      verbose = 1L)),
      type = "message"
    )
    # Final verbose line: "iEPPA <status> in <N> iters, errRp=..."
    m <- tail(grep("in [0-9]+ iters", msgs, value = TRUE), 1)
    iter <- as.integer(sub(".*in ([0-9]+) iters.*", "\\1", m))
    list(res = res, iter = iter)
  }

  r_stable <- run_one("off")
  r_damped <- run_one("on")
  expect_true(all(is.finite(r_stable$res)))
  expect_true(all(is.finite(r_damped$res)))
  expect_gt(r_damped$iter, r_stable$iter)  # monotone; spec §7
})
```

- [ ] **Step 3.2: Run the new tests — confirm RED**

Run: `R CMD INSTALL --preclean . && Rscript -e 'library(testthat); library(leafblower); test_file("tests/testthat/test-ieppa-faithful.R", reporter="summary")'`
Expected: stable-mode test (identical result across two runs) PASSES pre-WU-3 already — deterministic input. Damped-mode test FAILS pre-WU-3 because `LBW_IEPPA_FORCE_DAMPING` has no effect (damping code path does not exist), so `r_stable$iter` and `r_damped$iter` are equal — `expect_gt(r_damped$iter, r_stable$iter)` fails. If the damped-mode test passes pre-fix, halt and diagnose (either damping already got sneaked in elsewhere, or the input isn't stressing enough — spec §7 requires a falsifiable monotone assertion).

- [ ] **Step 3.3: Edit `src/ieppa.cpp` — add alpha state, auto-trigger, and test-only force override**

Immediately after the `linear_fallback_used` declaration (from Step 2.4), add:

```cpp
    // WU-3: adaptive damping. Default alpha=1.0 (stable mode, byte-identical to pre-WU-3).
    // Auto-switch to alpha=0.5 (damped mode) when any bucket's streak reaches
    // kInfeasPersistence/2 = 2. Latched per-solve; does not revert.
    double alpha = 1.0;
    bool damped_latched = false;
    // Test-only override (parallel to LBW_IEPPA_FORCE_PATH): "on"|"off"|unset.
    // Always compiled; microsecond getenv cost. Enables falsifiable
    // iter_damped > iter_stable assertion (spec §7, CTO B5).
    const char* force_damp = std::getenv("LBW_IEPPA_FORCE_DAMPING");
    bool force_damping_on  = (force_damp != nullptr && std::strcmp(force_damp, "on")  == 0);
    bool force_damping_off = (force_damp != nullptr && std::strcmp(force_damp, "off") == 0);
    if (force_damping_on) { alpha = 0.5; damped_latched = true; }
    auto maybe_engage_damping = [&]() {
        if (damped_latched || force_damping_off) return;
        for (int idx = 0; idx < total_cats; idx++) {
            if (infeas_streak[idx] >= kInfeasPersistence / 2) {
                alpha = 0.5;
                damped_latched = true;
                if (st.verbose >= 1) {
                    st.log("iEPPA: mid-streak detected; damping engaged (alpha=0.5).");
                }
                return;
            }
        }
    };
```

- [ ] **Step 3.4: Edit `src/ieppa.cpp` — damp the log-space Sinkhorn update**

Inside the log-space sweep (else branch from Step 2.5), replace the success-path update:

```cpp
                    record_nonempty(k, j);
                    double log_target = std::log(st.targets[k][j] * ct.W_input);
                    lf[cat_offset[k] + j] = log_target - log_S_kj;
```

with:

```cpp
                    record_nonempty(k, j);
                    double log_target = std::log(st.targets[k][j] * ct.W_input);
                    if (alpha == 1.0) {
                        lf[cat_offset[k] + j] = log_target - log_S_kj;
                    } else {
                        double lf_old = lf[cat_offset[k] + j];
                        lf[cat_offset[k] + j] =
                            (1.0 - alpha) * lf_old
                            + alpha * (log_target - log_S_kj);
                    }
```

- [ ] **Step 3.5: Edit `src/ieppa.cpp` — damp the linear-space Sinkhorn update**

Inside the linear-space sweep (from Step 2.5), replace:

```cpp
                    record_nonempty(k, j);
                    double new_f = (st.targets[k][j] * ct.W_input) / S_kj;
                    if (!std::isfinite(new_f) || new_f > kLinearOverflowTrip) {
                        overflow_trip = true;
                        break;
                    }
                    f_lin[cat_offset[k] + j] = new_f;
```

with:

```cpp
                    record_nonempty(k, j);
                    double naive = (st.targets[k][j] * ct.W_input) / S_kj;
                    if (!std::isfinite(naive) || naive > kLinearOverflowTrip) {
                        overflow_trip = true;
                        break;
                    }
                    double new_f;
                    if (alpha == 1.0) {
                        new_f = naive;
                    } else {
                        double f_old = f_lin[cat_offset[k] + j];
                        new_f = std::pow(f_old, 1.0 - alpha)
                              * std::pow(naive, alpha);
                    }
                    if (!std::isfinite(new_f) || new_f > kLinearOverflowTrip) {
                        overflow_trip = true;
                        break;
                    }
                    f_lin[cat_offset[k] + j] = new_f;
```

- [ ] **Step 3.6: Edit `src/ieppa.cpp` — call `maybe_engage_damping` once per outer iter**

At the top of the outer iter loop (just after `res.iterations = iter;`), add:

```cpp
        maybe_engage_damping();
```

- [ ] **Step 3.6.5: Patch WU-2 fallback blocks to reset `alpha`/`damped_latched` (spec §5 state-clean checklist)**

WU-2 landed before WU-3 and therefore could not reset WU-3-only state. Two fallback sites added in Step 2.5 and Step 2.6 need amending.

In Step 2.5's overflow-trip fallback (just after `persistent_infeas_pairs.clear();`), insert:

```cpp
                alpha = 1.0;
                damped_latched = false;
                if (force_damping_on) { alpha = 0.5; damped_latched = true; }
```

Apply the identical insertion to Step 2.6's X_tilde-overflow fallback (same location: just after `persistent_infeas_pairs.clear();`).

This honours spec §5 exactly ("`alpha` reset to 1.0 (stable mode; WU-3 re-triggers if streak re-builds)"). The `force_damping_on` re-application preserves the test-only forced mode across fallback — consistent with `LBW_IEPPA_FORCE_PATH` semantics (env var dominates auto-dispatch).

- [ ] **Step 3.7: Build gate**

Run: `R CMD INSTALL --preclean .`
Expected: clean install, no warnings.

- [ ] **Step 3.8: Run WU-3 tests + full regression — expect GREEN**

Run: `Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", stop_on_failure=TRUE, reporter="summary")'`
Expected: `[ FAIL 0 | PASS >= 176 ]` (174 post-WU-2 + 2 new).

- [ ] **Step 3.9: Stepstone merge-gate verification**

Run: `Rscript /tmp/stepstone_3algo_bench.R 2>&1 | tee /tmp/wu3_stepstone.log`
Expected acceptance criteria (all must hold):
- `ieppa` row: status NOT `ERROR ... infeasible`; status is `OK` or `NOCONV`.
- `ieppa` max_err ≤ `raking` max_err on same input (raking ~8.25e-3).
- `ieppa` wall-clock ≤ 2× raking's (raking ~2.3s → iEPPA ≤ ~5s target).
If any fails, halt and report which gate broke; do not commit until the acceptance block holds.

- [ ] **Step 3.10: kk1204 per-iter re-check**

Run: `Rscript /tmp/wu2_kk1204.R` (same script from Step 2.11)
Expected: `ratio` still ≤ 2.0× (damping must not regress per-iter ratio).

- [ ] **Step 3.11: Commit**

```bash
git add src/ieppa.cpp tests/testthat/test-ieppa-faithful.R \
        tests/testthat/_snapshots/ieppa_stable_prewu3.rds
git commit -m "$(cat <<'EOF'
feat(ieppa): adaptive damping auto-triggered on persistence stress

Geometric-blend damping applied to Sinkhorn update (both log-space and
linear-space paths). Default alpha=1.0 fast-path via explicit branch
(std::pow(_, 1.0) is NOT a no-op — costs transcendental overhead).

Auto-engage at infeas_streak[idx] >= kInfeasPersistence/2; latched
per-solve, resets between solves. Shares infeas_streak state introduced
by WU-1; no additional public API, no rk_params_t change.

Convergence rate slows by ~2x in damped mode (Peyré-Cuturi 2019 §4.4 eq
4.52 "softened updates"); damped mode engages only when stable mode is
already headed toward NOCONV.

Refs spec docs/superpowers/specs/2026-04-23-ieppa-convergence-hardening-design.md WU-3.
EOF
)"
```

---

## Post-implementation acceptance

- [ ] **Step A.1: Final regression sweep**

Run: `Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", stop_on_failure=TRUE, reporter="summary")'`
Expected: `[ FAIL 0 | PASS >= 176 ]`.

- [ ] **Step A.2: R CMD check**

Run: `R CMD build . && R CMD check leafblower_*.tar.gz --no-manual --as-cran`
Expected: 0 ERRORs, 0 WARNINGs. NOTES unchanged from baseline (if any). `benchmarks/` remains excluded from the tarball per `.Rbuildignore`.

- [ ] **Step A.3: Python build sanity (no wheel release; just sanity)**

Run: `pip install -e python/ --quiet && python -c "import leafblower; print(leafblower.__version__)"`
Expected: installs cleanly, prints version. If install fails, confirm `python/CMakeLists.txt` still lists `cell_table.cpp`, `raking.cpp`, `ieppa.cpp` and has `-mavx2`; this plan does not add C++ files, so no CMakeLists change is needed.

- [ ] **Step A.4: Acceptance matrix (from spec §11)**

Capture final numbers:
- Stepstone full data (`Rscript /tmp/stepstone_3algo_bench.R`): iEPPA returns RK_OK or RK_ERR_NOCONV (NOT RK_ERR_INFEAS); errRp ≤ raking errRp; wall-clock ≤ 2× raking.
- kk1204 regime (`Rscript /tmp/wu2_kk1204.R`): iEPPA/raking per-iter ratio ≤ 2×.
- Full suite: FAIL=0, PASS ≥ post-baseline + 6 (2 WU-1 + 3 WU-2 + 2 WU-3 − 1 smoke-only reshape; baseline from Step P.5).

Record the numbers in the merge PR body.

- [ ] **Step A.4.1: RK_OK-preservation regression gate (spec §11 last bullet)**

Spec §11 explicitly requires "any test previously returning RK_OK must still return RK_OK (no new RK_ERR_INFEAS from tightened persistence logic)". FAIL=0 alone does NOT prove this — a flipped status could still pass if a test only checks `expect_no_error` or `expect_true(is.finite(...))`. Compare per-(file,test) status against the baseline captured in Step P.5 using the `/tmp/baseline-per-test.rds` snapshot:

```bash
Rscript -e '
  library(testthat); library(leafblower)
  tr <- test_dir("tests/testthat",
                 reporter = SilentReporter$new(),
                 stop_on_failure = FALSE)
  classify <- function(r) {
    if (inherits(r, "expectation_skip"))    "skip"
    else if (inherits(r, "expectation_failure") ||
             inherits(r, "expectation_error"))  "fail"
    else if (inherits(r, "expectation_warning")) "warn"
    else                                          "pass"
  }
  post <- list()
  for (ti in tr) {
    block_cls <- "pass"
    for (r in ti$results) {
      cc <- classify(r)
      if (cc == "fail") { block_cls <- "fail"; break }
      if (cc == "warn" && block_cls == "pass") block_cls <- "warn"
      if (cc == "skip" && block_cls == "pass") block_cls <- "skip"
    }
    post[[paste0(ti$file, "::", ti$test)]] <- block_cls
  }
  baseline <- readRDS("/tmp/baseline-per-test.rds")
  regressions <- list()
  for (k in names(baseline)) {
    b <- baseline[[k]]
    p <- post[[k]]
    if (is.null(p)) { regressions[[k]] <- paste0(b, "->MISSING"); next }
    if (p == b) next  # unchanged: OK
    # PASS -> anything other than PASS is a regression per spec §11.
    if (b == "pass" && p != "pass") regressions[[k]] <- paste0(b, "->", p)
  }
  if (length(regressions) > 0) {
    cat("RK_OK PRESERVATION REGRESSIONS (spec §11 last bullet):\n")
    for (k in names(regressions)) cat("  ", k, ": ", regressions[[k]], "\n", sep="")
    quit(status = 1L)
  }
  cat("RK_OK preservation: ", length(baseline), " pre-existing test blocks all preserved.\n", sep="")
'
```

Expected: exit 0 with "`all preserved`" line. Exit 1 with a list of regressing blocks means a pre-existing test flipped from `pass` to something worse — HALT and investigate. This is the falsifiable gate for spec §11's last bullet.

- [ ] **Step A.5: Update beads**

```bash
bd close leafblower-bel --reason="WU-1 persistent-infeas tracker lands; stepstone no longer false-positive INFEAS"
bd close leafblower-g8f --reason="WU-2 linear-space dispatch lands; kk1204 per-iter ratio <=2x"
```

(Do not close `leafblower-3c2` — AUTO routing fallback is a separate ticket.)

---

## Self-review checklist

Run through these once after finishing the plan, fix inline:

1. **Spec coverage.** WU-1 (persistent-infeas + wrapper error msg) → Task 1 ✓. WU-2 (linear-space dispatch + overflow fallback + env var) → Task 2 ✓. WU-3 (geometric-blend damping + alpha=1 fast path + auto-trigger) → Task 3 ✓. Merge-gate acceptance (stepstone no INFEAS, errRp ≤ raking, wall-clock ≤ 2×, 169+ tests pass, kk1204 ≤ 2×) → Step A.4 ✓.

2. **No placeholders.** Every code step shows actual code. Commands are exact. No "TBD" or "similar to Task N".

3. **Type consistency.** `infeas_streak` (vector<int>, size=total_cats), `persistent_infeas_pairs` (set<pair<int,int>>), `record_empty`/`record_nonempty` (lambdas), `kInfeasPersistence` (constexpr int = 5), `kLinearSpaceThreshold` (constexpr double = 2.0), `kLinearOverflowTrip` (runtime double), `f_lin` (vector<double>, size=total_cats), `alpha` (double), `damped_latched` (bool), `force_damping_on`/`force_damping_off` (bool), `linear_fallback_used` (bool), `use_linear` (bool), `LBW_IEPPA_FORCE_PATH` (env var: "linear"|"log"), `LBW_IEPPA_FORCE_DAMPING` (env var: "on"|"off") — all used consistently across Tasks 1–3, with WU-2 fallback patched by Step 3.6.5 to include WU-3 state.

4. **Atomic ordering respected.** WU-1 commit lands first (primary correctness). WU-2 references WU-1's lambdas. WU-3 reads WU-1's `infeas_streak` and branches on `alpha`. Tasks cannot be reordered.
