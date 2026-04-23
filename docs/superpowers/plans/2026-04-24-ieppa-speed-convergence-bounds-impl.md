# iEPPA Speed / Convergence / Bounds Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close kk1204 per-iter gap, reach tol_abs=1e-3 on kk1204, and add strict per-observation bounds option — without regressing stepstone or breaking existing API callers.

**Architecture:** Four atomic commits on `src/ieppa.cpp` (+ wrapper files for P3.1). P1.1 fuses three post-sweep O(M_cell) passes into one. P2.1 replaces hard-0.5 damping latch with smooth schedule `alpha = 1/(1+β·stress)`. P2.2 adds Anderson(m=5) acceleration on `lf` via LAPACK `dgels` with triple-layer guard (INFO, isfinite, γ-norm). P3.1 adds `bounds_mode` param with intra-cell water-filling for `bounds_mode="unit"`.

**Tech Stack:** C++17, R (testthat), Python (pytest). LAPACK via `R_ext/Lapack.h`. No new deps.

**Source spec:** `docs/superpowers/specs/2026-04-24-ieppa-speed-convergence-bounds-design.md` (rev 5, commit 32987ab). Post-review-gate at 4 iterations with user auth to exceed the 3-iter cap.

**Commit ordering (atomic):** P1.1 → P2.1 → P2.2 → P3.1. P2.2 depends on P2.1's `alpha` state. P3.1 is file-scope-disjoint from P2.x. P1.1 must land before P2.2 because Anderson reads `n_cap_active` from the fused capacity block.

**Build gate after each source edit:** `R CMD INSTALL --preclean .` must succeed before the next step. Never edit a second source file before confirming the first compiles.

**Baseline:** commit `1b29df1` (WU-4 split structural/transient infeas landed). Post-WU-4: 181 tests green, stepstone errRp=2.21e-3, kk1204 per-iter 2.17× raking, kk1204 NOCONV at 500 iter.

---

## Pre-flight

- [ ] **Step P.1: Confirm clean working tree**

Run: `git status --short`
Expected: only untracked items in `.gitignore` (build artifacts, `.wolf/`, etc.). No modified tracked files.

- [ ] **Step P.2: Baseline test suite**

Run: `R CMD INSTALL --preclean . && Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", stop_on_failure=TRUE, reporter="summary")'`
Expected: `[ FAIL 0 | WARN 2 | SKIP 0 | PASS 181 ]`.

- [ ] **Step P.3: Capture baseline metrics**

Run:
```bash
Rscript /tmp/wu2_kk1204.R | tee /tmp/baseline_kk1204.log
Rscript /tmp/stepstone_2algo.R | tee /tmp/baseline_stepstone.log
```
Record values in `/tmp/baseline-summary.txt`:
- kk1204 ieppa per_iter_ms
- kk1204 raking per_iter_ms
- kk1204 ratio (should be ~2.17×)
- stepstone ieppa errRp (should be ~2.21e-3)
- stepstone wall-clock (should be ~3.5s)

These anchor per-pass regression checks.

---

## Task 1: P1.1 — Fuse post-sweep X_tilde + capacity + X_cur rebuild

**Files:**
- Modify: `src/ieppa.cpp` (linear path post-sweep; roughly lines 362–447)
- Modify: `src/ieppa.hpp` (add `IEPPAResult::n_xcur_writes_per_iter_linear`)
- Modify: `tests/testthat/test-ieppa-faithful.R` (new RED test on counter)

**Rationale:** Current linear path runs three O(M_cell) passes after the sweep: X_tilde compute, capacity block, X_cur rebuild. Each loads/stores the same cells. Spec §4.1: collapse to one pass using `X_cur[c] / W[c]` as inline X_tilde — the sweep already maintains `X_cur = X_init · W · ∏ f_lin`. Memory traffic 3× → 1× for the fused section.

**Clean-code discipline:** Single-responsibility function; names reveal intent (`X_tilde_c` local, not `tmp`); explicit overflow branch; no slop comments (spec rationale lives in the spec).

### Step 1.1: Add counter field to IEPPAResult

Modify `src/ieppa.hpp` — append to the struct:

```cpp
int n_xcur_writes_per_iter_linear;  // 0 outside linear path; counter for P1.1 RED test
```

Zero-initialize in `src/ieppa.cpp::ieppa_solve` near `res.n_cap_active = 0;`:

```cpp
res.n_xcur_writes_per_iter_linear = 0;
```

### Step 1.2: Write the failing RED test

Append to `tests/testthat/test-ieppa-faithful.R`:

```r
test_that("P1.1: linear path writes X_cur exactly M_cell times per iter (fused block)", {
  # Dense-compression input routes to linear path. Counter asserts the fused
  # block touches each X_cur[c] exactly once per outer iter (not 2x or 3x).
  Sys.setenv(LBW_IEPPA_FORCE_PATH = "linear")
  on.exit(Sys.unsetenv("LBW_IEPPA_FORCE_PATH"), add = TRUE)
  set.seed(991)
  n <- 5000L
  K <- 4L
  df <- as.data.frame(replicate(K, sample(1:3, n, replace = TRUE), simplify = FALSE))
  names(df) <- paste0("m", 1:K)
  for (k in names(df)) df[[k]] <- c("a","b","c")[df[[k]]]
  targets <- setNames(
    replicate(K, c(a = 0.4, b = 0.35, c = 0.25), simplify = FALSE),
    paste0("m", 1:K)
  )
  # attach_weights=FALSE → res is a numeric vector; but we need the result struct.
  # Access via diagnose-return harness: harvest returns data frame; fetch diagnostic
  # counter from the attribute set by the C bridge (add this via r_bridge.cpp plumb).
  res <- harvest(df, targets, method = "ieppa",
                 max_weight = 5, min_weight = 0,
                 max_iterations = 20L,
                 convergence = list(absolute = 1e-300),
                 attach_weights = FALSE)
  # The C result plumbs n_xcur_writes_per_iter_linear onto attr(res, "result")$n_xcur_writes_per_iter_linear
  result_info <- attr(res, "result")
  expect_true(!is.null(result_info$n_xcur_writes_per_iter_linear))
  # M_cell per iter × iter count. M_cell for this input: at most 3^4 = 81 cells;
  # actual is obs count, measured via C_leafblower_cell_table_probe in other tests.
  # We assert the ratio (writes / iter) equals the number of cells written per iter,
  # which for the fused block is exactly M_cell. Pre-P1.1 it was 2 or 3× M_cell.
  stopifnot(result_info$iterations > 0)
  writes_per_iter <- result_info$n_xcur_writes_per_iter_linear / result_info$iterations
  # M_cell bounded above by n; lower bound is 1. Assert ratio ≤ M_cell (not 2×M_cell).
  # Fetch M_cell via the probe C call.
  gid_list <- lapply(names(targets), function(nm) {
    lv <- names(targets[[nm]])
    idx <- match(as.character(df[[nm]]), lv) - 1L
    idx[is.na(idx)] <- -1L
    as.integer(idx)
  })
  probe <- .Call("C_leafblower_cell_table_probe", gid_list, n, PACKAGE = "leafblower")
  expect_equal(writes_per_iter, probe$M_cell)
})
```

### Step 1.3: Run test — expect RED

Run: `R CMD INSTALL --preclean . && Rscript -e 'library(testthat); library(leafblower); test_file("tests/testthat/test-ieppa-faithful.R", reporter="summary")'`
Expected: the new test FAILS. Pre-fix, the counter is zero (no writes tracked) OR it equals 2×M_cell or 3×M_cell (post-sweep loops count). The `expect_equal` fails.

If the counter is plumbed correctly but the test passes vacuously, halt and check: the test-only counter field must be readable by the R bridge.

### Step 1.4: Plumb counter to R via r_bridge.cpp

Check `src/r_bridge.cpp` — where `IEPPAResult` is consumed. Append `n_xcur_writes_per_iter_linear` to the R list returned alongside `iterations`, `max_error`, etc.

Find the block that builds the R result list; add:

```cpp
SET_VECTOR_ELT(result, <next_idx>, PROTECT(Rf_ScalarInteger(ieppa_result.n_xcur_writes_per_iter_linear)));
SET_STRING_ELT(result_names, <next_idx>, Rf_mkChar("n_xcur_writes_per_iter_linear"));
```

Increment the list size and nprotect count accordingly.

### Step 1.5: Edit src/ieppa.cpp — replace three post-sweep passes with fused block

Locate the current linear-path post-sweep code (roughly `ieppa.cpp:362–447`):
1. X_tilde compute loop
2. Capacity block (sets W, X)
3. X_cur rebuild loop

Replace all three with:

```cpp
        if (use_linear) {
            // P1.1 fused block: X_tilde derived inline as X_cur / W; capacity + X_cur rebuild fused.
            bool overflow_detected = false;
            int n_cap = 0;
            for (int c = 0; c < ct.M_cell; c++) {
                if (X_init[c] <= 0.0 || W[c] <= 0.0) {
                    X[c] = 0.0;
                    X_cur[c] = 0.0;
                    W[c] = 1.0;
                    continue;
                }
                double X_tilde_c = X_cur[c] / W[c];
                if (!std::isfinite(X_tilde_c) || X_tilde_c > kLinearOverflowTrip) {
                    overflow_detected = true;
                    break;
                }
                if (X_tilde_c <= 0.0) {
                    X[c] = 0.0;
                    X_cur[c] = 0.0;
                    W[c] = 1.0;
                    continue;
                }
                double xc = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);
                X[c] = xc;
                W[c] = xc / X_tilde_c;
                X_cur[c] = xc;
                res.n_xcur_writes_per_iter_linear++;
                if (xc != X_tilde_c) n_cap++;
            }
            res.n_cap_active = n_cap;
            if (overflow_detected) {
                // Full state reset on mid-loop break; partial writes to W/X/X_cur undone.
                std::fill(X_cur.begin(),   X_cur.end(),   0.0);
                std::fill(W.begin(),       W.end(),       1.0);
                std::fill(X.begin(),       X.end(),       0.0);
                std::fill(X_tilde.begin(), X_tilde.end(), 0.0);
                std::fill(lf.begin(),      lf.end(),      0.0);
                std::fill(f_lin.begin(),   f_lin.end(),   1.0);
                std::fill(infeas_streak.begin(), infeas_streak.end(), 0);
                res.n_xcur_writes_per_iter_linear = 0;
                use_linear = false;
                linear_fallback_used = true;
                if (st.verbose >= 1) st.log("iEPPA: linear-space overflow trip; fallback to log-space.");
                continue;
            }
        } else {
            // Log-path: X_tilde + capacity + X_cur unchanged from current implementation.
            // ... (leave existing log-path code intact) ...
        }
```

Keep the existing log-path code inside the `else` branch untouched. The existing WU-2 fallback block (the other overflow detector that fires during sweep) also remains unchanged.

### Step 1.6: Build gate

Run: `R CMD INSTALL --preclean .`
Expected: clean install, no warnings on `src/ieppa.cpp` or `src/r_bridge.cpp`. Halt on any warning on the edited code paths.

### Step 1.7: Run P1.1 test — expect GREEN

Run: `Rscript -e 'library(testthat); library(leafblower); test_file("tests/testthat/test-ieppa-faithful.R", reporter="summary")'`
Expected: the new test PASSES. `writes_per_iter == M_cell`.

### Step 1.8: Full regression

Run: `Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", stop_on_failure=TRUE, reporter="summary")'`
Expected: `[ FAIL 0 | PASS ≥ 182 ]` (baseline 181 + 1 new).

### Step 1.9: Speed + stepstone regression

Run:
```bash
Rscript /tmp/wu2_kk1204.R
Rscript /tmp/stepstone_2algo.R
```
Expected:
- kk1204 per-iter ratio ≤ 1.5× raking (spec §8 gate; down from 2.17×).
- stepstone ieppa errRp ≤ 2.21e-3 + 1e-4 (no regression).
- stepstone wall ≤ 4.0s (no regression vs 3.5s baseline).

Halt if any gate fails.

### Step 1.10: Commit

```bash
git add src/ieppa.cpp src/ieppa.hpp src/r_bridge.cpp tests/testthat/test-ieppa-faithful.R
git commit -m "$(cat <<'EOF'
perf(ieppa): fuse post-sweep X_tilde + capacity + X_cur rebuild (P1.1)

Three O(M_cell) passes collapse to one by deriving X_tilde inline as
X_cur[c] / W[c] — the sweep already maintains X_cur = X_init · W · ∏ f_lin
invariant per WU-2. Memory traffic drops from 3× to 1× on linear path.

Linear-overflow detection integrated inline with full state reset
(X_cur, W, X, X_tilde, lf, f_lin, infeas_streak) before falling back to
log-space. Ensures mid-loop break does not carry partial cell writes into
the log-space restart.

Adds IEPPAResult::n_xcur_writes_per_iter_linear counter for the RED test
asserting the fused block touches each X_cur[c] exactly once per outer iter.

kk1204 per-iter ratio 2.17× → ≤ 1.5× raking (spec §8 gate).

Refs spec docs/superpowers/specs/2026-04-24-ieppa-speed-convergence-bounds-design.md P1.1.
EOF
)"
```

---

## Task 2: P2.1 — Adaptive damping schedule

**Files:**
- Modify: `src/ieppa.cpp` (remove `damped_latched`; replace `alpha = 0.5` latch with schedule)
- Modify: `src/ieppa.hpp` (add `double min_alpha_seen`, `double final_alpha`)
- Modify: `tests/testthat/test-ieppa-faithful.R` (append RED tests)

**Rationale:** WU-3 latched `alpha = 0.5` permanently on first stress. Stepstone trajectory probe (2026-04-24 session) showed damping engaged at iter 260, then stayed damped for 240 more iters without recovery — halving rate unnecessarily. Smooth schedule `alpha = 1/(1+β·stress)` recovers as streaks reset.

**Clean-code discipline:** Remove `damped_latched` (dead state). Single `alpha_for_iter(stress, beta)` helper if helpful; otherwise inline the formula. `force_damping_on/off` semantics documented in one place.

### Step 2.1: Add alpha state fields to IEPPAResult

Modify `src/ieppa.hpp` — append:

```cpp
double min_alpha_seen;   // min alpha over all sweeps; 1.0 if damping never engaged
double final_alpha;      // alpha at solver exit (after last sweep)
```

Initialize in `ieppa_solve`:

```cpp
res.min_alpha_seen = 1.0;
res.final_alpha = 1.0;
```

### Step 2.2: Write failing RED tests

Append to `tests/testthat/test-ieppa-faithful.R`:

```r
test_that("P2.1: benign input keeps alpha == 1.0 (fast path)", {
  set.seed(55)
  n <- 500L
  df <- data.frame(a = sample(letters[1:3], n, TRUE),
                   b = sample(letters[1:3], n, TRUE))
  targets <- list(a = c(a=0.33,b=0.33,c=0.34), b = c(a=0.33,b=0.33,c=0.34))
  res <- harvest(df, targets, method = "ieppa",
                 max_weight = 5, min_weight = 0,
                 max_iterations = 500L,
                 convergence = list(absolute = 1e-6),
                 attach_weights = FALSE)
  info <- attr(res, "result")
  expect_equal(info$min_alpha_seen, 1.0)  # no stress → no damping
  expect_equal(info$final_alpha, 1.0)
})

test_that("P2.1: stress input engages damping (alpha < 1.0) with smooth schedule", {
  set.seed(314)
  n <- 3000L
  K <- 6L
  df <- as.data.frame(replicate(K, sample(1:3, n, TRUE), simplify=FALSE))
  names(df) <- paste0("m", 1:K)
  for (k in names(df)) df[[k]] <- c("a","b","c")[df[[k]]]
  targets <- setNames(
    replicate(K, c(a=0.7, b=0.2, c=0.1), simplify=FALSE),
    paste0("m", 1:K)
  )
  res <- suppressWarnings(harvest(df, targets, method = "ieppa",
                                  max_weight = 3, min_weight = 0,
                                  max_iterations = 500L,
                                  convergence = list(absolute = 1e-4),
                                  attach_weights = FALSE))
  info <- attr(res, "result")
  expect_lt(info$min_alpha_seen, 1.0)   # stress engaged damping at some point
  # Unlatched schedule: if streaks subside before exit, alpha recovers.
  # Final alpha may be 1.0 (full recovery) or intermediate. Not asserted strict.
  expect_true(info$min_alpha_seen > 0.0)  # sanity: formula is bounded below
})

test_that("P2.1: LBW_IEPPA_FORCE_DAMPING=on forces alpha = 0.5", {
  Sys.setenv(LBW_IEPPA_FORCE_DAMPING = "on")
  on.exit(Sys.unsetenv("LBW_IEPPA_FORCE_DAMPING"), add = TRUE)
  set.seed(55)
  n <- 500L
  df <- data.frame(a = sample(letters[1:3], n, TRUE),
                   b = sample(letters[1:3], n, TRUE))
  targets <- list(a = c(a=0.33,b=0.33,c=0.34), b = c(a=0.33,b=0.33,c=0.34))
  res <- harvest(df, targets, method = "ieppa",
                 max_weight = 5, min_weight = 0,
                 max_iterations = 500L,
                 convergence = list(absolute = 1e-6),
                 attach_weights = FALSE)
  info <- attr(res, "result")
  expect_equal(info$min_alpha_seen, 0.5)
  expect_equal(info$final_alpha, 0.5)
})
```

### Step 2.3: Run test — expect RED

Run: `R CMD INSTALL --preclean . && Rscript -e 'library(testthat); library(leafblower); test_file("tests/testthat/test-ieppa-faithful.R", reporter="summary")'`
Expected: three new tests FAIL — `min_alpha_seen` / `final_alpha` fields don't exist or aren't populated correctly (old code latches at 0.5 permanently, failing test 2's recovery semantics and test 3's force-damping=0.5 specific).

### Step 2.4: Edit src/ieppa.cpp — replace damping state with schedule

Locate the existing damping block (search for `damped_latched` and `alpha = 0.5`). Replace with:

```cpp
    // P2.1 adaptive damping schedule. alpha = 1 / (1 + β · stress); stress = max streak.
    // β = 0.5: at stress=2, alpha=0.5; at stress=10, alpha≈0.17. Unlatched — recovers as
    // streaks reset. Preserves Peyré-Cuturi §4.4 convergence (alpha ∈ (0, 1]).
    constexpr double kBeta = 0.5;
    double alpha = 1.0;

    const char* force_damp = std::getenv("LBW_IEPPA_FORCE_DAMPING");
    bool force_damping_on  = (force_damp != nullptr && std::strcmp(force_damp, "on")  == 0);
    bool force_damping_off = (force_damp != nullptr && std::strcmp(force_damp, "off") == 0);

    auto compute_alpha = [&]() -> double {
        if (force_damping_on)  return 0.5;
        if (force_damping_off) return 1.0;
        int stress = 0;
        for (int idx = 0; idx < total_cats; idx++) {
            if (infeas_streak[idx] > stress) stress = infeas_streak[idx];
        }
        if (stress == 0) return 1.0;
        return 1.0 / (1.0 + kBeta * static_cast<double>(stress));
    };
```

At the top of the outer iter loop (before the sweep):

```cpp
        alpha = compute_alpha();
        if (alpha < res.min_alpha_seen) res.min_alpha_seen = alpha;
```

At the bottom of the outer iter loop (just before the convergence check or loop end):

```cpp
        res.final_alpha = alpha;
```

Remove all references to `damped_latched` (field, updates, checks).

In both the log-space and linear-space Sinkhorn update sites, the existing `if (alpha == 1.0) { ... naive ... } else { ... damped ... }` branch already works with the new `alpha`. No change needed there.

### Step 2.5: Plumb min_alpha_seen / final_alpha to R

Modify `src/r_bridge.cpp` — add fields to the returned R list:

```cpp
SET_VECTOR_ELT(result, <idx>, PROTECT(Rf_ScalarReal(ieppa_result.min_alpha_seen)));
SET_STRING_ELT(result_names, <idx>, Rf_mkChar("min_alpha_seen"));
SET_VECTOR_ELT(result, <idx+1>, PROTECT(Rf_ScalarReal(ieppa_result.final_alpha)));
SET_STRING_ELT(result_names, <idx+1>, Rf_mkChar("final_alpha"));
```

### Step 2.6: Build gate + tests

Run: `R CMD INSTALL --preclean . && Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", stop_on_failure=TRUE, reporter="summary")'`
Expected: `[ FAIL 0 | PASS ≥ 185 ]`.

### Step 2.7: P2.1 interim assertion

kk1204 with default (auto) damping, 500 iter, `ACCEL_ANDERSON=off`: errRp @ 500 iter ≤ 2.5e-3. Spec §9 interim-state guard.

Run:
```bash
LBW_IEPPA_ACCEL_ANDERSON=off Rscript /tmp/wu2_kk1204.R
```
Wait — that script doesn't assert errRp, only per-iter. Add a one-shot check:

```bash
Rscript -e '
  Sys.setenv(OMP_NUM_THREADS = "1", LBW_IEPPA_ACCEL_ANDERSON = "off")
  library(leafblower)
  set.seed(1); n <- 1000000L; K <- 20L
  df <- as.data.frame(replicate(K, sample(1:5, n, replace=TRUE), simplify=FALSE))
  names(df) <- paste0("m", 1:K)
  for (k in names(df)) df[[k]] <- letters[df[[k]]]
  tgts <- setNames(replicate(K, c(a=0.3,b=0.175,c=0.175,d=0.175,e=0.175), simplify=FALSE),
                   names(df))
  res <- suppressWarnings(harvest(df, tgts, method="ieppa",
    max_weight=3, min_weight=0, max_iterations=500L,
    convergence=list(absolute=1e-6), attach_weights=FALSE))
  info <- attr(res, "result")
  cat(sprintf("max_err=%.3e iters=%d min_alpha=%.3f\n",
              info$max_error, info$iterations, info$min_alpha_seen))
  stopifnot(info$max_error <= 2.5e-3)
'
```
Expected: `max_err ≤ 2.5e-3`. If higher, the adaptive schedule is regressing convergence — halt.

### Step 2.8: Stepstone regression

Run: `Rscript /tmp/stepstone_2algo.R`
Expected: ieppa errRp ≤ 2.21e-3 + 1e-4. No regression.

### Step 2.9: Commit

```bash
git add src/ieppa.cpp src/ieppa.hpp src/r_bridge.cpp tests/testthat/test-ieppa-faithful.R
git commit -m "$(cat <<'EOF'
feat(ieppa): adaptive damping schedule α = 1/(1 + β · stress) (P2.1)

Replace WU-3's hard-0.5 alpha latch with a smooth, unlatched schedule
where β = 0.5 and stress = max_{k,j} infeas_streak[cat_offset[k]+j].
At stress=0 alpha=1.0 (fast path), at stress=2 alpha=0.5, at stress=10
alpha≈0.17. Recovers monotonically as streaks reset — solver resumes
full Sinkhorn rate once transients subside.

Removes damped_latched state variable. LBW_IEPPA_FORCE_DAMPING retains
its test-only override: "on" → alpha=0.5 constant, "off" → alpha=1.0
constant regardless of streak.

Adds IEPPAResult::min_alpha_seen and final_alpha for struct-based
regression tests. No log-string parsing required.

Preserves Peyré-Cuturi 2019 §4.4 convergence (α ∈ (0, 1]).

Refs spec §5.1 P2.1.
EOF
)"
```

---

## Task 3: P2.2 — Anderson(m=5) acceleration

**Files:**
- Modify: `src/ieppa.cpp` (add Anderson block after capacity; shadow-lf for linear path)
- Modify: `src/ieppa.hpp` (add counters: `n_anderson_iters_engaged`, `n_anderson_nan_fallbacks`)
- Create: none (Anderson state is local to `ieppa_solve`)
- Modify: `tests/testthat/test-ieppa-faithful.R` (append RED tests)
- Modify: `src/r_bridge.cpp` (plumb two counters)

**Rationale:** Classical Sinkhorn is sublinear. Anderson mixing on `lf` accelerates via O(m) LS solve per iter. kk1204 at K=20 stays NOCONV at 500 iter baseline; Anderson engaged on uncapacitated iters targets 2× fewer iters to RK_OK.

**Clean-code discipline:** One clearly-named function `apply_anderson(lf, history, m_active)` if the code grows beyond ~40 lines inline. Constants like `kAndersonWarmup = 5`, `kGammaNormMax = 1e4` named in `constexpr`. No cargo-cult; every LAPACK flag explained at its call site.

### Step 3.1: Add counter fields

Modify `src/ieppa.hpp`:

```cpp
int n_anderson_iters_engaged;   // iters where Anderson actually fired (not warmup, not skip)
int n_anderson_nan_fallbacks;   // iters Anderson skipped due to INFO/isfinite/γ-norm guard
```

Initialize to zero in `ieppa_solve`.

### Step 3.2: Write failing RED tests

```r
test_that("P2.2: ACCEL_ANDERSON=off → zero engaged iters", {
  Sys.setenv(LBW_IEPPA_ACCEL_ANDERSON = "off")
  on.exit(Sys.unsetenv("LBW_IEPPA_ACCEL_ANDERSON"), add = TRUE)
  set.seed(77)
  n <- 1000L
  df <- data.frame(a = sample(letters[1:3], n, TRUE),
                   b = sample(letters[1:3], n, TRUE))
  targets <- list(a = c(a=0.4,b=0.35,c=0.25), b = c(a=0.4,b=0.35,c=0.25))
  res <- harvest(df, targets, method = "ieppa",
                 max_weight = 5, min_weight = 0,
                 max_iterations = 50L,
                 convergence = list(absolute = 1e-300),
                 attach_weights = FALSE)
  info <- attr(res, "result")
  expect_equal(info$n_anderson_iters_engaged, 0)
})

test_that("P2.2: ACCEL_ANDERSON=on engages Anderson post-warmup on uncapacitated iters", {
  Sys.setenv(LBW_IEPPA_ACCEL_ANDERSON = "on")
  on.exit(Sys.unsetenv("LBW_IEPPA_ACCEL_ANDERSON"), add = TRUE)
  set.seed(77)
  n <- 1000L
  df <- data.frame(a = sample(letters[1:3], n, TRUE),
                   b = sample(letters[1:3], n, TRUE))
  targets <- list(a = c(a=0.4,b=0.35,c=0.25), b = c(a=0.4,b=0.35,c=0.25))
  res <- harvest(df, targets, method = "ieppa",
                 max_weight = 5, min_weight = 0,
                 max_iterations = 50L,
                 convergence = list(absolute = 1e-300),
                 attach_weights = FALSE)
  info <- attr(res, "result")
  expect_gt(info$n_anderson_iters_engaged, 0)  # some iters fired
  # Warmup skips first 5 iters. Plus any cap-active iters. Net: at most 50-5 = 45.
  expect_lte(info$n_anderson_iters_engaged, info$iterations - 5)
})

test_that("P2.2: kk1204 with Anderson converges in ≤ 400 iter AND ≥ 2× fewer than off", {
  skip_if_not(file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
              "kk1204 probe data not available")
  set.seed(1)
  n <- 200000L  # scaled-down kk1204-shape for CI speed
  K <- 20L
  df <- as.data.frame(replicate(K, sample(1:5, n, replace=TRUE), simplify=FALSE))
  names(df) <- paste0("m", 1:K)
  for (k in names(df)) df[[k]] <- letters[df[[k]]]
  tgts <- setNames(replicate(K, c(a=0.3,b=0.175,c=0.175,d=0.175,e=0.175), simplify=FALSE),
                   names(df))

  Sys.setenv(LBW_IEPPA_ACCEL_ANDERSON = "off")
  res_off <- suppressWarnings(harvest(df, tgts, method="ieppa",
    max_weight=3, min_weight=0, max_iterations=500L,
    convergence=list(absolute=1e-3), attach_weights=FALSE))
  info_off <- attr(res_off, "result")

  Sys.setenv(LBW_IEPPA_ACCEL_ANDERSON = "on")
  res_on <- suppressWarnings(harvest(df, tgts, method="ieppa",
    max_weight=3, min_weight=0, max_iterations=500L,
    convergence=list(absolute=1e-3), attach_weights=FALSE))
  info_on <- attr(res_on, "result")
  Sys.unsetenv("LBW_IEPPA_ACCEL_ANDERSON")

  expect_equal(info_on$status_code, 0L)  # RK_OK = 0
  expect_lte(info_on$iterations, 400L)
  # At least 2× fewer: iter_on ≤ iter_off / 2 (both must be finite positive)
  expect_lte(info_on$iterations, info_off$iterations / 2L)
})
```

### Step 3.3: Run — expect RED

Run: `R CMD INSTALL --preclean . && Rscript -e 'library(testthat); library(leafblower); test_file("tests/testthat/test-ieppa-faithful.R", reporter="summary")'`
Expected: new tests FAIL. No Anderson code exists; fields are zero; convergence gate not reached.

### Step 3.4: Implement Anderson in src/ieppa.cpp

At the top of `ieppa_solve`, after existing WU-3 state:

```cpp
    // P2.2 Anderson(m=5) acceleration on lf via LAPACK dgels.
    constexpr int    kAndersonM        = 5;
    constexpr int    kAndersonWarmup   = 5;
    constexpr double kGammaNormMax     = 1e4;

    const char* force_accel = std::getenv("LBW_IEPPA_ACCEL_ANDERSON");
    const bool  anderson_off = (force_accel != nullptr && std::strcmp(force_accel, "off") == 0);
    const bool  anderson_enabled = !anderson_off;  // default on

    // History buffers: column-major, total_cats rows × kAndersonM cols.
    // F[:,j] = r_{t-j} (residual difference); X_hist[:,j] = lf_{t-j} (iterate difference).
    std::vector<double> lf_prev(total_cats, 0.0);    // lf at start of previous sweep
    std::vector<double> r_prev(total_cats, 0.0);     // previous residual
    std::vector<double> F_hist(total_cats * kAndersonM, 0.0);
    std::vector<double> X_hist(total_cats * kAndersonM, 0.0);
    std::vector<double> r_curr(total_cats, 0.0);
    std::vector<double> gamma(kAndersonM, 0.0);
    std::vector<double> lapack_work;  // sized by workspace query on first call
    int m_active = 0;                  // number of history columns currently filled
```

Inside the outer iter loop, AFTER the capacity block (fused in P1.1), BEFORE the errRp check:

```cpp
        const bool can_anderson = anderson_enabled
                                  && iter > kAndersonWarmup
                                  && res.n_cap_active == 0
                                  && total_cats > 1;
        if (!can_anderson) {
            // Still-capacitated or warmup: clear history (residuals mixed different W no good).
            m_active = 0;
        } else {
            // Materialize lf for the linear path via shadow-log.
            if (use_linear) {
                for (int kj = 0; kj < total_cats; kj++) {
                    lf[kj] = (f_lin[kj] > 0.0) ? std::log(f_lin[kj]) : -kLogClip;
                }
            }
            // Current residual r = G(lf) - lf = lf_this_iter - lf_prev.
            for (int kj = 0; kj < total_cats; kj++) {
                r_curr[kj] = lf[kj] - lf_prev[kj];
            }
            // Shift history right by one column; append r_curr - r_prev (residual diff) and lf - lf_prev.
            if (m_active > 0) {
                // Shift columns [0..m_active-1] → [1..m_active]; new col at 0.
                int shift_cols = std::min(m_active, kAndersonM - 1);
                for (int j = shift_cols; j > 0; j--) {
                    std::memcpy(&F_hist[j * total_cats], &F_hist[(j - 1) * total_cats],
                                total_cats * sizeof(double));
                    std::memcpy(&X_hist[j * total_cats], &X_hist[(j - 1) * total_cats],
                                total_cats * sizeof(double));
                }
            }
            for (int kj = 0; kj < total_cats; kj++) {
                F_hist[kj] = r_curr[kj] - r_prev[kj];
                X_hist[kj] = lf[kj]      - lf_prev[kj];
            }
            if (m_active < kAndersonM) m_active++;

            // Cap m_active for the shape guard (dgels requires m ≥ n for overdetermined LS).
            int m_solve = std::min(m_active, total_cats - 1);
            if (m_solve <= 0) {
                lf_prev = lf;
                r_prev  = r_curr;
            } else {
                // Copy r_curr into B buffer (dgels overwrites B with solution's first m_solve entries).
                std::vector<double> B_buf = r_curr;  // length total_cats
                std::vector<double> A_buf(F_hist.begin(), F_hist.begin() + total_cats * m_solve);

                // Workspace query.
                int n_rows = total_cats, n_cols = m_solve, nrhs = 1, lda = total_cats, ldb = total_cats, info = 0;
                double wkopt = 0.0; int lwork = -1;
                F77_CALL(dgels)("N", &n_rows, &n_cols, &nrhs,
                                A_buf.data(), &lda, B_buf.data(), &ldb,
                                &wkopt, &lwork, &info FCONE);
                if (info == 0) {
                    lwork = static_cast<int>(wkopt);
                    if (static_cast<int>(lapack_work.size()) < lwork) lapack_work.resize(lwork);
                    F77_CALL(dgels)("N", &n_rows, &n_cols, &nrhs,
                                    A_buf.data(), &lda, B_buf.data(), &ldb,
                                    lapack_work.data(), &lwork, &info FCONE);
                }

                // Three-layer guard.
                bool ok = (info == 0);
                double gamma_absmax = 0.0;
                if (ok) {
                    for (int j = 0; j < m_solve; j++) {
                        double g = B_buf[j];
                        if (!std::isfinite(g)) { ok = false; break; }
                        if (std::fabs(g) > gamma_absmax) gamma_absmax = std::fabs(g);
                    }
                    if (ok && gamma_absmax > kGammaNormMax) ok = false;
                }

                if (!ok) {
                    m_active = 0;   // clear history; use plain iterate
                    res.n_anderson_nan_fallbacks++;
                } else {
                    // Apply Anderson mixing: lf_new = lf_current - X_hist · γ + small correction via F_hist.
                    // Standard AA: lf_{t+1}_AA = lf + r - (X_hist + F_hist) · γ.
                    for (int kj = 0; kj < total_cats; kj++) {
                        double delta = 0.0;
                        for (int j = 0; j < m_solve; j++) {
                            delta += (X_hist[j * total_cats + kj] + F_hist[j * total_cats + kj]) * B_buf[j];
                        }
                        lf[kj] = lf[kj] + r_curr[kj] - delta;
                    }
                    res.n_anderson_iters_engaged++;
                }

                lf_prev = lf;     // after Anderson update
                r_prev  = r_curr;
            }

            // Rematerialize f_lin on linear path after lf mutation.
            if (use_linear) {
                for (int kj = 0; kj < total_cats; kj++) {
                    f_lin[kj] = std::exp(lf[kj]);
                    // Rebuild X_cur consistently: X_cur[c] = X_init[c] · W[c] · ∏_m f_lin[m][g_m(c)]
                    // is not needed here — next sweep recomputes from scratch.
                }
            }
        }
```

Header add: `#include <cstring>` (for memcpy) if not already present. `FCONE` macro comes from `R_ext/RS.h` (usually via `Rcpp.h` or `R_ext/Lapack.h`); if not defined, `#define FCONE FCLEN` may be needed — check R version compatibility.

### Step 3.5: Plumb counters + status_code to R bridge

Add to `r_bridge.cpp`:

```cpp
SET_VECTOR_ELT(result, <idx>,     PROTECT(Rf_ScalarInteger(ieppa_result.n_anderson_iters_engaged)));
SET_STRING_ELT(result_names, <idx>, Rf_mkChar("n_anderson_iters_engaged"));
SET_VECTOR_ELT(result, <idx+1>,   PROTECT(Rf_ScalarInteger(ieppa_result.n_anderson_nan_fallbacks)));
SET_STRING_ELT(result_names, <idx+1>, Rf_mkChar("n_anderson_nan_fallbacks"));
SET_VECTOR_ELT(result, <idx+2>,   PROTECT(Rf_ScalarInteger(ieppa_result.status)));
SET_STRING_ELT(result_names, <idx+2>, Rf_mkChar("status_code"));
```

### Step 3.6: Build gate + tests

Run: `R CMD INSTALL --preclean . && Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", stop_on_failure=TRUE, reporter="summary")'`
Expected: `[ FAIL 0 | PASS ≥ 188 ]`. All 3 new Anderson tests pass.

### Step 3.7: kk1204 convergence gate

Run the kk1204 convergence check:

```bash
Rscript -e '
  Sys.setenv(OMP_NUM_THREADS = "1", LBW_IEPPA_ACCEL_ANDERSON = "on")
  library(leafblower)
  set.seed(1); n <- 1000000L; K <- 20L
  df <- as.data.frame(replicate(K, sample(1:5, n, replace=TRUE), simplify=FALSE))
  names(df) <- paste0("m", 1:K)
  for (k in names(df)) df[[k]] <- letters[df[[k]]]
  tgts <- setNames(replicate(K, c(a=0.3,b=0.175,c=0.175,d=0.175,e=0.175), simplify=FALSE),
                   names(df))
  res <- suppressWarnings(harvest(df, tgts, method="ieppa",
    max_weight=3, min_weight=0, max_iterations=500L,
    convergence=list(absolute=1e-3), attach_weights=FALSE))
  info <- attr(res, "result")
  cat(sprintf("status=%d iters=%d max_err=%.3e anderson_fired=%d fallbacks=%d\n",
              info$status_code, info$iterations, info$max_error,
              info$n_anderson_iters_engaged, info$n_anderson_nan_fallbacks))
  stopifnot(info$status_code == 0L, info$iterations <= 400L)
'
```
Expected: `status=0 iters≤400`. If NOCONV or iters > 400: halt and inspect. Check `n_anderson_iters_engaged > 0` (Anderson actually engaged at all).

### Step 3.8: Also verify Anderson OFF baseline

Same script with `LBW_IEPPA_ACCEL_ANDERSON=off` — record `iters_off`. Spec §8: require `iters_on ≤ iters_off / 2`.

### Step 3.9: Stepstone regression

Run: `Rscript /tmp/stepstone_2algo.R`
Expected: errRp ≤ 2.21e-3 + 1e-4.

### Step 3.10: Commit

```bash
git add src/ieppa.cpp src/ieppa.hpp src/r_bridge.cpp tests/testthat/test-ieppa-faithful.R
git commit -m "$(cat <<'EOF'
perf(ieppa): Anderson(m=5) acceleration on lf via LAPACK dgels (P2.2)

Adds Anderson type-II acceleration to the Sinkhorn iterate after each
outer iter's capacity block. Operates on lf (log-space factors) for both
paths: linear-path maintains shadow-lf = log(f_lin) with f_lin > 0 guard.

Engagement gated on n_cap_active == 0 (capacity BCD changes the effective
operator; mixing residuals across different W breaks contraction). History
cleared on any capacity-active iter; warmup 5 iters.

Three-layer rank-deficiency / near-singular guard after dgels: INFO == 0
AND all(isfinite(γ)) AND max|γ| ≤ 1e4. Any failure → clear history + plain
iterate + increment n_anderson_nan_fallbacks.

Shape guard: m_active = min(m, total_cats - 1) prevents dgels entering the
underdetermined branch (min-norm solutions ≠ LS intent).

LBW_IEPPA_ACCEL_ANDERSON ∈ {on, off, unset} env var. Default on.

Targets kk1204 iter-to-RK_OK ≤ 400 AND ≥ 2× fewer than off-baseline
(spec §8 gate).

Refs spec §5.2 P2.2.
EOF
)"
```

---

## Task 4: P3.1 — `bounds_mode` param + intra-cell water-filling

**Files:**
- Modify: `src/leafblower.h` (enum + struct field + static_assert)
- Modify: `src/c_api.cpp` (memset + bounds_mode plumbing)
- Modify: `src/ieppa.cpp` (water-filling expansion for unit mode; diagnostic scan for cell mode)
- Modify: `src/ieppa.hpp` (add `n_bounds_violated`, `n_bounds_clamped`)
- Modify: `src/r_bridge.cpp` (accept bounds_mode arg; return bounds fields)
- Modify: `R/harvest.R` (parse_bounds_mode helper + harvest() param)
- Modify: `python/leafblower/_harvest.py` (bounds_mode kwarg)
- Create: `tests/testthat/test-ieppa-bounds-mode.R`

**Rationale:** Current cell-aggregate bounding can produce individual `w_i ∉ [min_weight, max_weight]` when base weights `d_i` are skewed within a cell. Opt-in `bounds_mode="unit"` enforces strict per-obs bounds via intra-cell water-filling (spec assessment Solution 2).

**Clean-code discipline:** Enum values named for statistical meaning (`RK_BOUNDS_CELL`/`RK_BOUNDS_UNIT`), not numbers. `water_fill_cell(cell_idx, obs_indices, min_w, max_w, target_sum)` extracted as a function with a single responsibility if the loop grows. R helper returns a validated character (matches `map_method` pattern).

### Step 4.1: Extend C API header

Modify `src/leafblower.h`:

```c
/* Bounds enforcement mode (appended field in rk_params_t; default RK_BOUNDS_CELL) */
typedef enum {
    RK_BOUNDS_CELL = 0,  /* cell-aggregate bounds (current default; may violate per-obs) */
    RK_BOUNDS_UNIT = 1   /* per-observation strict bounds via intra-cell water-filling */
} rk_bounds_mode_t;
```

Append to `rk_params_t` (after `log_ctx`):

```c
rk_bounds_mode_t bounds_mode;  /* default RK_BOUNDS_CELL */
```

Extend `rk_result_t`:

```c
int n_bounds_violated;  /* cell-mode diagnostic: count of w_i outside bounds (no action) */
int n_bounds_clamped;   /* unit-mode action: count of w_i clamped after water-fill exhausted */
```

At end of header, before `#endif`:

```c
/* ABI tripwires. If EXPECTED_RK_PARAMS_BYTES fails, a new field was added —
 * update this value after auditing ABI consumers. */
#ifdef __cplusplus
static_assert(RK_ALG_AUTO == 0, "memset(0) default must equal RK_ALG_AUTO");
/* Compute sizeof(rk_params_t) on the target platform at implementation time
 * and hard-code it here. Record the value in a comment. Example:
 *   Linux x86_64 GCC 13, verified 2026-04-24: 72 bytes (6 doubles + 3 ints +
 *   1 enum + 1 fn-ptr + 1 void* + 1 enum + padding).
 * After measuring, replace the placeholder below with the actual value. */
#define EXPECTED_RK_PARAMS_BYTES 0 /* TODO: replace with measured sizeof after first build */
/* static_assert(sizeof(rk_params_t) == EXPECTED_RK_PARAMS_BYTES,
 *               "rk_params_t size changed; check ABI consumers"); */
#endif
```

Uncomment the static_assert in a second pass after measuring the actual size.

### Step 4.2: Patch rk_params_init

Modify `src/c_api.cpp`:

Prepend `memset(p, 0, sizeof(*p));` as the first statement of `rk_params_init`:

```cpp
void rk_params_init(rk_params_t* p) {
    if (p == nullptr) return;
    memset(p, 0, sizeof(*p));
    /* ... existing field assignments ... */
    p->bounds_mode = RK_BOUNDS_CELL;  // explicit for clarity; memset=0 already gives this
    /* ... rest unchanged ... */
}
```

Ensure `#include <cstring>` is present.

### Step 4.3: Add IEPPAResult fields + zero-init

Modify `src/ieppa.hpp`:

```cpp
int n_bounds_violated;
int n_bounds_clamped;
```

Zero-init in `ieppa_solve`.

### Step 4.4: Write failing RED tests

Create `tests/testthat/test-ieppa-bounds-mode.R`:

```r
skewed_d_input <- function(n = 500L, seed = 404) {
  set.seed(seed)
  d <- c(runif(n - 5, 0.5, 1.5), rexp(5, rate = 0.2))  # 5 heavy-tailed weights
  df <- data.frame(
    a = sample(letters[1:3], n, TRUE, prob = c(0.5, 0.3, 0.2)),
    b = sample(letters[1:3], n, TRUE, prob = c(0.3, 0.4, 0.3))
  )
  df$design_weight <- d
  list(df = df, targets = list(a = c(a=0.4,b=0.35,c=0.25),
                               b = c(a=0.33,b=0.33,c=0.34)))
}

test_that("P3.1: default bounds_mode='cell' preserves current behaviour", {
  fx <- skewed_d_input()
  res_default <- harvest(fx$df, fx$targets, method = "ieppa",
                         max_weight = 3, min_weight = 0.3,
                         design_weights = fx$df$design_weight,
                         max_iterations = 500L,
                         convergence = list(absolute = 1e-5),
                         attach_weights = FALSE)
  res_explicit <- harvest(fx$df, fx$targets, method = "ieppa",
                          max_weight = 3, min_weight = 0.3,
                          design_weights = fx$df$design_weight,
                          bounds_mode = "cell",
                          max_iterations = 500L,
                          convergence = list(absolute = 1e-5),
                          attach_weights = FALSE)
  expect_lt(max(abs(as.numeric(res_default) - as.numeric(res_explicit))), 1e-12)
})

test_that("P3.1: cell-mode emits warning + n_bounds_violated > 0 on skewed-d", {
  fx <- skewed_d_input()
  expect_warning(
    res <- harvest(fx$df, fx$targets, method = "ieppa",
                   max_weight = 3, min_weight = 0.3,
                   design_weights = fx$df$design_weight,
                   bounds_mode = "cell",
                   max_iterations = 500L,
                   convergence = list(absolute = 1e-5),
                   attach_weights = FALSE),
    regexp = "cell-mode bounds"
  )
  info <- attr(res, "result")
  expect_gt(info$n_bounds_violated, 0)
  expect_equal(info$n_bounds_clamped, 0)  # no clamping in cell mode
})

test_that("P3.1: unit-mode produces strict per-obs bounds", {
  fx <- skewed_d_input()
  res <- harvest(fx$df, fx$targets, method = "ieppa",
                 max_weight = 3, min_weight = 0.3,
                 design_weights = fx$df$design_weight,
                 bounds_mode = "unit",
                 max_iterations = 500L,
                 convergence = list(absolute = 1e-5),
                 attach_weights = FALSE)
  expect_lte(max(as.numeric(res)), 3 + 1e-9)
  expect_gte(min(as.numeric(res)), 0.3 - 1e-9)
  info <- attr(res, "result")
  # Benign-ish skew: degenerate single-obs cells rare → clamped count small.
  expect_lt(info$n_bounds_clamped, 0.001 * nrow(fx$df))
})

test_that("P3.1: invalid bounds_mode raises clear error", {
  fx <- skewed_d_input()
  expect_error(
    harvest(fx$df, fx$targets, method = "ieppa", bounds_mode = "invalid"),
    regexp = "should be one of"  # match.arg message
  )
})

test_that("P3.1: C API integer mapping (raw .Call with bounds_mode = 1L)", {
  fx <- skewed_d_input()
  # Directly invoke the C bridge with raw integer 1 → unit mode.
  # This validates the C enum integer mapping independently of the R helper.
  # ... (call via .Call wrapper; exact form depends on r_bridge.cpp signature) ...
  # For now: assert that harvest(..., bounds_mode="unit") and a direct int-based call
  # produce identical results. Skip if bridge doesn't expose int-based entry.
  skip("requires direct .Call bridge for bounds_mode; implement if supported")
})
```

### Step 4.5: Run — expect RED

Run: `R CMD INSTALL --preclean . && Rscript -e 'library(testthat); library(leafblower); test_file("tests/testthat/test-ieppa-bounds-mode.R", reporter="summary")'`
Expected: new tests FAIL — `bounds_mode` argument doesn't exist yet, `n_bounds_violated`/`n_bounds_clamped` zero.

### Step 4.6: Implement water-filling expansion in src/ieppa.cpp

Replace the existing final expansion block (roughly:
```cpp
for (int c = 0; c < ct.M_cell; c++) mult[c] = ...;
for (int i = 0; i < st.n; i++) st.weights[i] = st.weights[i] * mult[ct.cell_of[i]];
```
) with path-branched logic:

```cpp
    // Expansion to observation weights.
    std::vector<double> mult(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        mult[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
    }
    for (int i = 0; i < st.n; i++) {
        st.weights[i] = st.weights[i] * mult[ct.cell_of[i]];
    }

    if (st.bounds_mode == RK_BOUNDS_CELL) {
        // Diagnostic scan: count violations without action.
        int violations = 0;
        for (int i = 0; i < st.n; i++) {
            if (st.weights[i] > st.max_weight || st.weights[i] < st.min_weight) {
                violations++;
            }
        }
        res.n_bounds_violated = violations;
    } else {
        // Unit mode: per-cell water-filling.
        // Build cells_of_obs (list of obs indices per cell) in one pass.
        std::vector<std::vector<int>> cells_of_obs(ct.M_cell);
        for (int i = 0; i < st.n; i++) cells_of_obs[ct.cell_of[i]].push_back(i);

        constexpr int kWaterFillMaxIter = 50;
        int total_clamped = 0;
        for (int c = 0; c < ct.M_cell; c++) {
            const auto& idxs = cells_of_obs[c];
            if (idxs.empty()) continue;

            double target_sum = X[c];
            for (int it = 0; it < kWaterFillMaxIter; it++) {
                double excess = 0.0;
                double free_sum = 0.0;
                int    n_free = 0;
                bool   any_violation = false;
                for (int i : idxs) {
                    if (st.weights[i] > st.max_weight) {
                        excess += st.weights[i] - st.max_weight;
                        st.weights[i] = st.max_weight;
                        any_violation = true;
                    } else if (st.weights[i] < st.min_weight) {
                        excess -= st.min_weight - st.weights[i];
                        st.weights[i] = st.min_weight;
                        any_violation = true;
                    } else {
                        free_sum += st.weights[i];
                        n_free++;
                    }
                }
                if (!any_violation) break;
                if (n_free == 0 || free_sum <= 0.0) {
                    // Pathological: no room to redistribute. Last-resort clamp survives.
                    for (int i : idxs) {
                        if (st.weights[i] > st.max_weight) { st.weights[i] = st.max_weight; total_clamped++; }
                        else if (st.weights[i] < st.min_weight) { st.weights[i] = st.min_weight; total_clamped++; }
                    }
                    break;
                }
                // Redistribute excess proportionally over free observations.
                double factor = 1.0 + excess / free_sum;
                for (int i : idxs) {
                    if (st.weights[i] > st.min_weight && st.weights[i] < st.max_weight) {
                        st.weights[i] *= factor;
                    }
                }
                if (it == kWaterFillMaxIter - 1) {
                    // Inner budget exhausted with violations remaining — final clamp.
                    for (int i : idxs) {
                        if (st.weights[i] > st.max_weight) { st.weights[i] = st.max_weight; total_clamped++; }
                        else if (st.weights[i] < st.min_weight) { st.weights[i] = st.min_weight; total_clamped++; }
                    }
                }
            }
        }
        res.n_bounds_clamped = total_clamped;
    }
```

### Step 4.7: Wire bounds_mode through c_api.cpp

Verify `st.bounds_mode` is populated from `rk_params_t::bounds_mode` in the adapter path. If `CalibState` does not already include `bounds_mode`, add it and copy from params.

### Step 4.8: R wrapper — helper + harvest() param

Modify `R/harvest.R`:

Add helper (near `map_method`):

```r
parse_bounds_mode <- function(x = c("cell", "unit")) {
  match.arg(x)
}
```

Add `bounds_mode = "cell"` to the `harvest(...)` signature.

In the body, before `.Call`:

```r
bounds_mode_char <- parse_bounds_mode(bounds_mode)
bounds_mode_int  <- match(bounds_mode_char, c("cell", "unit")) - 1L
```

Pass `bounds_mode_int` into the `.Call` alongside other params. Add to the R bridge `.Call` signature.

Add roxygen `@param bounds_mode`:

```r
#' @param bounds_mode
#'   Bounds enforcement mode. "cell" (default): per-cell aggregate bounds —
#'   sum of weights within each cross-classified demographic cell is bounded,
#'   but individual weights may fall outside [min_weight, max_weight] if
#'   initial design weights are skewed within a cell. "unit": per-observation
#'   bounds via intra-cell water-filling redistribution. In degenerate single-
#'   obs cells, strict bounds cannot always be guaranteed;
#'   `result$n_bounds_clamped` reports residual clamp count.
```

Emit warnings at harvest() exit:

```r
if (!is.null(calib_result$n_bounds_violated) && calib_result$n_bounds_violated > 0) {
  warning(sprintf(
    "cell-mode bounds: %d weights fell outside [%.3f, %.3f] due to skewed base weights within cells. Consider bounds_mode='unit' for strict per-observation bounds.",
    calib_result$n_bounds_violated, min_weight, max_weight))
}
if (!is.null(calib_result$n_bounds_clamped) && calib_result$n_bounds_clamped > 0) {
  warning(sprintf(
    "unit-mode bounds: %d weights clamped to nearest bound after water-fill exhausted; degenerate cells could not be fully redistributed.",
    calib_result$n_bounds_clamped))
}
```

### Step 4.9: Python wrapper

Modify `python/leafblower/_harvest.py`:

Add `bounds_mode: str = "cell"` to `harvest()` signature. Near top:

```python
if bounds_mode not in ("cell", "unit"):
    raise ValueError(f"bounds_mode must be 'cell' or 'unit', got {bounds_mode!r}")
_bounds_mode_int = {"cell": 0, "unit": 1}[bounds_mode]
```

Pass `_bounds_mode_int` into the C API call. Docstring mirrors R roxygen.

Post-result warning:

```python
if result_dict.get("n_bounds_violated", 0) > 0:
    warnings.warn(
        f"cell-mode bounds: {result_dict['n_bounds_violated']} weights fell outside "
        f"[{min_weight:.3f}, {max_weight:.3f}]. Consider bounds_mode='unit'.",
        UserWarning, stacklevel=2)
if result_dict.get("n_bounds_clamped", 0) > 0:
    warnings.warn(
        f"unit-mode bounds: {result_dict['n_bounds_clamped']} weights clamped after water-fill exhausted.",
        UserWarning, stacklevel=2)
```

### Step 4.10: Plumb fields to R bridge

Modify `src/r_bridge.cpp` — accept `bounds_mode_int` as a new arg; write `n_bounds_violated`, `n_bounds_clamped` into the returned list.

### Step 4.11: Measure sizeof(rk_params_t) and re-enable static_assert

After a successful build of Steps 4.1–4.10:

```bash
Rscript -e '.Call("C_leafblower_params_sizeof", PACKAGE="leafblower")'
```

(Add this probe if it doesn't exist — one-line helper in `c_api.cpp`.) Alternatively compile a tiny probe:

```bash
g++ -E -c -x c -I./src - <<'EOF'
#include "src/leafblower.h"
#include <cstddef>
_Static_assert(sizeof(rk_params_t) == 72, "measure this");
EOF
```

Record the measured value. Replace `#define EXPECTED_RK_PARAMS_BYTES 0` with the actual value in `src/leafblower.h`. Uncomment the `static_assert`. Rebuild.

### Step 4.12: Build + tests

Run:
```bash
R CMD INSTALL --preclean .
Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", stop_on_failure=TRUE, reporter="summary")'
```
Expected: `[ FAIL 0 | PASS ≥ 192 ]`.

### Step 4.13: Stepstone + kk1204 regression

Run:
```bash
Rscript /tmp/wu2_kk1204.R
Rscript /tmp/stepstone_2algo.R
```
Expected: no regression on any metric.

### Step 4.14: Commit

```bash
git add src/leafblower.h src/c_api.cpp src/ieppa.hpp src/ieppa.cpp src/r_bridge.cpp \
        R/harvest.R python/leafblower/_harvest.py \
        tests/testthat/test-ieppa-bounds-mode.R
git commit -m "$(cat <<'EOF'
feat(ieppa): bounds_mode parameter + intra-cell water-filling (P3.1)

Adds rk_bounds_mode_t {RK_BOUNDS_CELL=0, RK_BOUNDS_UNIT=1} to rk_params_t.
Default RK_BOUNDS_CELL preserves current behaviour byte-for-byte.

RK_BOUNDS_UNIT engages post-expansion intra-cell water-filling: per cell,
clamp out-of-bound weights and redistribute excess over in-bound observations
until no violations remain (50-iter fallback). Degenerate single-obs cells
may require final clamp; count reported via result.n_bounds_clamped.

Cell-mode diagnostic: scan expanded weights for violations, emit count via
result.n_bounds_violated (no clamping). R/Python wrappers warn when either
counter > 0 so users can see when skewed base weights produce out-of-bound
individual weights under the default mode.

rk_params_init now memsets the struct before field assignments — prevents
uninitialized bounds_mode in stack-allocated callers. Paired static_assert
on RK_ALG_AUTO == 0 and sizeof(rk_params_t) tripwires.

parse_bounds_mode R helper returns character ("cell" | "unit"); conversion
to C enum integer at .Call site. Python wrapper uses str type.

ABI note: rk_params_t grows by sizeof(int) + padding. Raw-struct callers
not using rk_params_init() must recompile.

Refs spec §6 P3.1.
EOF
)"
```

---

## Post-implementation acceptance

- [ ] **Step A.1: Final regression sweep**

Run: `Rscript -e 'library(testthat); library(leafblower); test_dir("tests/testthat", stop_on_failure=TRUE, reporter="summary")'`
Expected: `[ FAIL 0 | PASS ≥ 192 ]`.

- [ ] **Step A.2: R CMD check**

Run: `R CMD build . && R CMD check leafblower_*.tar.gz --no-manual --as-cran`
Expected: 0 ERRORs, 0 WARNINGs.

- [ ] **Step A.3: Python sanity**

Run: `pip install -e python/ --quiet && python -c "import leafblower; help(leafblower.harvest)"`
Expected: docstring includes `bounds_mode` param.

- [ ] **Step A.4: Merge-gate matrix**

Capture final numbers from a single run:

```bash
Rscript /tmp/wu2_kk1204.R      # per-iter ratio ≤ 1.5×
Rscript /tmp/stepstone_2algo.R # errRp ≤ 2.21e-3 + 1e-4
# Plus the kk1204 convergence gate from Step 3.7
```

Record in PR body:
- kk1204 per-iter ratio (P1.1 gate)
- kk1204 iter-to-RK_OK with Anderson on/off (P2.2 gate)
- Stepstone errRp + wall (no regression)
- n_bounds_clamped on unit-mode skewed input (< 0.1% of n)

- [ ] **Step A.5: Env var interaction matrix (spec §7.4)**

Run the 7 required combos (T1–T7) from a parameterized script. Assert the output weights for each combo are finite; no crashes; T6/T7 are the bare-baseline regression guards.

- [ ] **Step A.6: Update beads**

Close the tickets:

```bash
bd close <ticket-id-for-P1.1> --reason="P1.1 landed <sha>; kk1204 per-iter <ratio>"
bd close <ticket-id-for-P2.1> --reason="P2.1 landed <sha>; min_alpha_seen observability"
bd close <ticket-id-for-P2.2> --reason="P2.2 landed <sha>; kk1204 RK_OK in <iters> iter"
bd close <ticket-id-for-P3.1> --reason="P3.1 landed <sha>; bounds_mode unit/cell ship"
```

---

## Self-review checklist

1. **Spec coverage.** P1.1 fused block → Task 1. P2.1 schedule → Task 2. P2.2 Anderson → Task 3. P3.1 bounds + water-fill → Task 4. Merge-gate items §8 → Steps 1.9, 2.7, 3.7, 4.13, A.4. Env-var matrix §7.4 → Step A.5. ✓

2. **Placeholder scan.** No "TBD", no "Add X" without code, no "similar to Task N" without inline code. EXPECTED_RK_PARAMS_BYTES is an explicit two-step instruction (measure + replace in Step 4.11), not a placeholder. ✓

3. **Type consistency.** `n_xcur_writes_per_iter_linear` (int), `min_alpha_seen`/`final_alpha` (double), `n_anderson_iters_engaged`/`n_anderson_nan_fallbacks` (int), `n_bounds_violated`/`n_bounds_clamped` (int), `alpha` (double), `kBeta`/`kGammaNormMax` (double), `kAndersonM`/`kAndersonWarmup`/`kWaterFillMaxIter` (int), `rk_bounds_mode_t` (enum), `bounds_mode` (C enum, R char, Python str), `LBW_IEPPA_ACCEL_ANDERSON` (env var) — all consistent across all tasks. ✓

4. **Atomic ordering.** P1.1 → P2.1 → P2.2 → P3.1. P2.2 uses `n_cap_active` from P1.1's fused block. P3.1 is file-scope-disjoint from P1/P2 but depends on none of them algorithmically. ✓

5. **Clean-code discipline (/clean-code invoked).**
   - Meaningful names: `X_tilde_c`, `kAndersonWarmup`, `parse_bounds_mode`, `n_anderson_nan_fallbacks` — all intention-revealing.
   - Small functions: `compute_alpha` lambda, `water_fill` inlined only because it's used once. Extract if grows.
   - No slop comments: comments explain the WHY (stress schedule rationale, capacity-BCD interaction), not the WHAT.
   - No comments restating the variable name.
   - Law of Demeter: `res.n_anderson_iters_engaged++` — one hop, fine.
   - Error handling: no bare exceptions; early returns on null ptr; explicit fallbacks on LAPACK failure.
   - Three-law TDD: every task writes failing test before implementation. ✓
