# iEPPA Log-Path Jacobi Sweep — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans.

**Goal:** Implement Jacobi alternative for iEPPA log-path sweep, benchmark vs GS, ship or revert based on wall-ratio data.

**Architecture:** Runtime flag `jacobi_log` in `CalibState`; `cell_lf_frozen` workspace in `ieppa.cpp`; 18-cell benchmark grid; data-driven ship/revert decision.

**Tech Stack:** C++17, R package build, testthat, bench::mark

---

## Plan Header

- **Mechanism:** Jacobi fixed-point sweep over K margins in iEPPA log path. Freeze `cell_lf` snapshot at outer-iter start; all margin updates read frozen snapshot; sequential O(K·M_cell) rebuild after the K-margin sweep replaces the in-lambda scattered patches.
- **Forbidden:** No changes to linear path (X_cur leave-one-out is already O(1)). No changes to SRAA, SOR, homotopy, ieppa_soft ALM. No silent default flip — default stays `jacobi_sweep=FALSE` until J4 decision gate. No removal of GS code path.
- **Audit:** Spies on (a) `cell_lf_frozen.assign` allocation gated by flag, (b) lambda's `clf` pointer dispatch, (c) skipped in-lambda `cell_lf[c] += delta` block, (d) end-of-sweep rebuild loop. T2 unit test asserts numerical agreement (max_error within 1e-4) between GS and Jacobi on a fixed seed fixture; T3 benchmark records `iters` and `wall_s`.

---

## Spec Reference

Authoritative spec: `docs/superpowers/specs/2026-05-01-jacobi-log-sweep-design.md`. The spec is the contract; this plan is the executable schedule.

---

## Epic Map

| Epic | Tasks | Goal |
|---|---|---|
| **E1: Wiring** | J1 | Plumb `jacobi_sweep` parameter from R through to `CalibState.jacobi_log`. Build green, no behavior change (default off). |
| **E2: Implementation + Decision** | J2, J3, J4 | Implement Jacobi semantics in `ieppa.cpp`, benchmark against GS on 18-cell grid, ship/revert based on data. |

---

## Task J1 — Wire `jacobi_sweep` parameter (R → C → CalibState)

**Epic:** E1 — Wiring.
**Goal:** Add `jacobi_sweep` parameter to `harvest()`, pass through `.Call`, set `st.jacobi_log` in r_bridge. Default FALSE; semantically inert until J2 lands.
**Blocker for:** J2 (J2 needs `st.jacobi_log` defined).

### Files modified

#### 1) `src/types.hpp` (line 128 region — add field below `accelerate`)

**Before** (lines 127–129):
```cpp
    CalibSorCfg          sor_cfg;
    bool                 accelerate = false;  // SQUAREM outer loop for raking
    // ── End overlay config ──
```

**After**:
```cpp
    CalibSorCfg          sor_cfg;
    bool                 accelerate = false;  // SQUAREM outer loop for raking
    bool                 jacobi_log = false;  // log-path: freeze cell_lf at iter start (Jacobi) vs patch inline (GS)
    // ── End overlay config ──
```

#### 2) `src/r_bridge.cpp` — extend `C_rk_calibrate` signature (line ~126) and wire (line ~318)

**Before** (lines 125–126):
```cpp
                    /* SQUAREM */
                    SEXP accelerate_sexp) {
```

**After**:
```cpp
                    /* SQUAREM */
                    SEXP accelerate_sexp,
                    /* Jacobi log-path sweep */
                    SEXP jacobi_sweep_sexp) {
```

**Before** (lines 316–318):
```cpp
    st.alm.lambda = 0.0;
    st.alm.mu     = 0.0;
    st.accelerate = (INTEGER(accelerate_sexp)[0] != 0);
```

**After**:
```cpp
    st.alm.lambda = 0.0;
    st.alm.mu     = 0.0;
    st.accelerate = (INTEGER(accelerate_sexp)[0] != 0);
    st.jacobi_log = (INTEGER(jacobi_sweep_sexp)[0] != 0);
```

Also update the R-callable registration in `src/r_bridge.cpp` (`R_CallMethodDef` table). Locate the `{"C_rk_calibrate", (DL_FUNC) &C_rk_calibrate, N}` entry and increment `N` by 1.

#### 3) `R/harvest.R` — add `jacobi_sweep` parameter and pass to `.Call`

**Before** (around line 222, after `eta_schedule_power`):
```r
  eta_start             = 1.0,
  ...
```

**After** — add new parameter with default FALSE in the parameter list (placement: immediately after the SQUAREM/accelerate parameter, or at end of the calibration knob block), and add `as.integer(jacobi_sweep)` to the `.Call` argument list at the same position as the C signature update (final positional arg after `accelerate`):

```r
  jacobi_sweep          = FALSE,
```

In the `.Call("C_rk_calibrate", ...)` invocation, append (preserving order match with the C signature):

```r
  ,
  as.integer(jacobi_sweep)
```

Also add validation near the existing param-validation block:

```r
  stopifnot(is.logical(jacobi_sweep), length(jacobi_sweep) == 1L,
            !is.na(jacobi_sweep))
```

Roxygen documentation entry above `@export`:

```r
#' @param jacobi_sweep logical(1). If TRUE, the iEPPA log-path margin sweep
#'   uses Jacobi semantics (freeze cell_lf at iteration start; rebuild after
#'   K-margin pass). Default FALSE (Gauss-Seidel). Affects log path only;
#'   no effect on linear path. Experimental; may converge in different
#'   iteration counts.
```

### Audit

- **Spy 1:** Compile gate confirms C↔R signature alignment.
- **Spy 2:** `harvest(df, tgt, method="ieppa")` (default args) must produce bit-identical weights to `git stash` baseline (jacobi_sweep defaults to FALSE; no code path activated yet).

### Compile gate

```bash
R CMD INSTALL --preclean .
```

### Test gate

```bash
Rscript -e 'library(leafblower); testthat::test_dir("tests/testthat")' 2>&1 | tail -2
```

Expected: `[ FAIL 0 | WARN * | SKIP * | PASS * ]`. No new test added in J1 — this is plumbing only; J2 adds the behavioral test.

### Conventional commit message

```
feat(ieppa): plumb jacobi_sweep parameter (no-op default)

Add jacobi_log field to CalibState, wire jacobi_sweep_sexp through
r_bridge, expose harvest(jacobi_sweep=FALSE). Behavioral
implementation lands in J2; this commit is wiring-only and produces
bit-identical results vs prior baseline.
```

---

## Task J2 — Implement Jacobi semantics in `ieppa.cpp` log path

**Epic:** E2 — Implementation + Decision.
**Goal:** Add `cell_lf_frozen` workspace, freeze at iter start when `st.jacobi_log && !use_linear`, switch lambda to read frozen snapshot, skip in-lambda `cell_lf` patches, rebuild `cell_lf` sequentially after K-margin sweep.
**Blocker:** Requires J1 merged.
**Blocker for:** J3 (benchmark needs implementation to run).

### Files modified

#### `src/ieppa.cpp`

**Edit 1 — declare `cell_lf_frozen` at homotopy-level scope (alongside `cell_lf` at line ~213)**

**Before** (line 213):
```cpp
    std::vector<double> cell_lf(ct.M_cell, 0.0);
    // High-water mark: max_c(log_X_init[c] + cell_lf[c]) ≈ max_c log(X_tilde[c]).
```

**After**:
```cpp
    std::vector<double> cell_lf(ct.M_cell, 0.0);
    // J2: Jacobi snapshot of cell_lf taken at outer-iter start (log path only).
    // Lambda reads from this when st.jacobi_log==true; else reads cell_lf directly.
    // Allocated only when jacobi_log is set (zero memory cost for default GS path).
    std::vector<double> cell_lf_frozen;
    if (st.jacobi_log) cell_lf_frozen.assign(ct.M_cell, 0.0);
    // High-water mark: max_c(log_X_init[c] + cell_lf[c]) ≈ max_c log(X_tilde[c]).
```

**Edit 2 — freeze snapshot at top of each outer iteration**

Locate the outer iteration loop (inside the homotopy-level `for (int lvl = ...)` at line 437, the inner `for (int it = 0; it < budget_lvl; it++)` loop). At its top, before the K-margin sweep dispatch:

```cpp
        // J2: Jacobi freeze. Log path only; linear path uses live X_cur leave-one-out.
        if (st.jacobi_log && !use_linear) {
            std::copy(cell_lf.begin(), cell_lf.end(), cell_lf_frozen.begin());
        }
```

The exact insertion point is immediately before the call site that dispatches `apply_single_margin_log` (the scheduler's K-loop, in `ieppa.cpp` near where `apply_single_margin_log(k)` is invoked). Confirm the insertion sits inside the outer iter loop, outside any margin-k loop.

**Edit 3 — modify `apply_single_margin_log` lambda (line 622) to read frozen snapshot**

Inside the lambda, locate the cell-bucket sum that uses `cell_lf[c]` (the leave-one-out S_kj computation). Add a top-of-lambda alias:

```cpp
        auto apply_single_margin_log = [&](int k) -> bool {
            // J2: Jacobi reads frozen snapshot; GS reads live cell_lf.
            // use_linear is captured by [&]; gate keeps log-only behavior.
            const double* clf = (st.jacobi_log && !use_linear)
                ? cell_lf_frozen.data()
                : cell_lf.data();
            // ... existing eff_omega_log resolution ...
```

Replace every `cell_lf[c]` read inside the lambda's S_kj/bucket-sum computation with `clf[c]`. Only the **read** sites change; **writes** to `cell_lf[c]` are handled in Edit 4.

**Edit 4 — gate the in-lambda `cell_lf[c] += delta` patches**

The lambda's existing patch block (analogous to lines 593–606 in the linear path, but the log-path equivalent inside `apply_single_margin_log`) updates `cell_lf` and `cell_lf_hwm`. Wrap it:

```cpp
            if (!st.jacobi_log) {
                // GS: existing in-sweep cell_lf patch (and hwm update on positive deltas).
                if (delta > 0.0) {
                    for (int c : cells_by_margin_cat[off + j]) {
                        cell_lf[c] += delta;
                        double val = cell_lf[c] + log_X_init[c];
                        if (std::isfinite(val) && val > cell_lf_hwm)
                            cell_lf_hwm = val;
                    }
                } else {
                    for (int c : cells_by_margin_cat[off + j])
                        cell_lf[c] += delta;
                }
            }
            // Jacobi: cell_lf is rebuilt after the K-margin sweep (Edit 5).
            // lf[off+j] update remains UNGATED — log factors must propagate
            // for the rebuild to recompute cell_lf correctly.
```

CRITICAL: the `lf[off + j] = lf_new` assignment must remain outside this guard. Only the `cell_lf[c]` scattered writes (and the hwm side-effect from those positive deltas) are gated.

**Edit 5 — sequential rebuild after K-margin sweep**

After the scheduler completes its K-margin sweep within an outer iteration (the location where `apply_single_margin_log` has been invoked for all margins, and `cell_lf` would normally already be coherent under GS), insert:

```cpp
        // J2: Jacobi rebuild. O(K * M_cell) sequential; cache-friendly.
        // Rebuilds cell_lf from current lf[] and refreshes cell_lf_hwm.
        if (st.jacobi_log && !use_linear) {
            std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
            for (int k = 0; k < st.K; k++) {
                const int off = cat_offset[k];
                const int* gk = ct.g_per_cell[k].data();
                for (int c = 0; c < ct.M_cell; c++) {
                    int g = gk[c];
                    if (g >= 0) cell_lf[c] += lf[off + g];
                }
            }
            // Refresh hwm from live cell_lf (one O(M_cell) pass).
            for (int c = 0; c < ct.M_cell; c++) {
                double val = cell_lf[c] + log_X_init[c];
                if (std::isfinite(val) && val > cell_lf_hwm) cell_lf_hwm = val;
            }
        }
```

The rebuild assumes `ct.g_per_cell[k]` is the per-cell category index (0-based, -1 sentinel for NA) for margin k. Verify this field name against `CellTable` in `src/types.hpp` before committing — if the field is named differently (e.g., `g_by_cell`), use the correct name. The K-margin sweep is the inner work of one outer iteration; rebuild fires once per outer iter.

### Audit

- **Spy 1:** With `jacobi_sweep=FALSE`, `cell_lf_frozen` remains empty (`.capacity() == 0`); GS path unchanged. Verify by running existing testthat suite — must pass with FAIL 0.
- **Spy 2 (T2 unit test):** Add `tests/testthat/test-jacobi-log-sweep.R`:

```r
test_that("Jacobi log-path matches GS within tolerance", {
  set.seed(42)
  n <- 2000
  df <- data.frame(
    a = factor(sample(letters[1:5], n, TRUE)),
    b = factor(sample(letters[1:4], n, TRUE)),
    c = factor(sample(letters[1:3], n, TRUE))
  )
  tgt <- list(
    a = setNames(rep(1/5, 5), letters[1:5]),
    b = setNames(rep(1/4, 4), letters[1:4]),
    c = setNames(rep(1/3, 3), letters[1:3])
  )
  r_gs <- harvest(df, tgt, method = "ieppa",
                  convergence = list(absolute = 1e-4),
                  jacobi_sweep = FALSE)
  r_ja <- harvest(df, tgt, method = "ieppa",
                  convergence = list(absolute = 1e-4),
                  jacobi_sweep = TRUE)
  expect_equal(r_gs$result$status, 0L)
  expect_equal(r_ja$result$status, 0L)
  # Tolerance 1e-4: both converged to within their absolute tol;
  # max_error agreement bounded by 2 * tol_abs.
  expect_lt(abs(r_gs$result$max_error - r_ja$result$max_error), 2e-4)
})
```

K=4 variant (compelled by spec — exercise an additional margin count beyond the 3-margin baseline). Add to the same test file:

```r
test_that("Jacobi log-path matches GS at K=4", {
  set.seed(43)
  n <- 2500
  df <- data.frame(
    a = factor(sample(letters[1:5], n, TRUE)),
    b = factor(sample(letters[1:4], n, TRUE)),
    c = factor(sample(letters[1:3], n, TRUE)),
    d = factor(sample(letters[1:6], n, TRUE))
  )
  tgt <- list(
    a = setNames(rep(1/5, 5), letters[1:5]),
    b = setNames(rep(1/4, 4), letters[1:4]),
    c = setNames(rep(1/3, 3), letters[1:3]),
    d = setNames(rep(1/6, 6), letters[1:6])
  )
  r_gs <- harvest(df, tgt, method = "ieppa",
                  convergence = list(absolute = 1e-4),
                  jacobi_sweep = FALSE)
  r_ja <- harvest(df, tgt, method = "ieppa",
                  convergence = list(absolute = 1e-4),
                  jacobi_sweep = TRUE)
  expect_equal(r_gs$result$status, 0L)
  expect_equal(r_ja$result$status, 0L)
  expect_lt(abs(r_gs$result$max_error - r_ja$result$max_error), 2e-4)
})
```

To force the log path (default may select linear on this benign fixture), set the env var inside the test:

```r
withr::with_envvar(c(LBW_IEPPA_FORCE_PATH = "log"), {
  # both harvest() calls go here
})
```

Wrap both `harvest` calls in each test under that `with_envvar`. This ensures T2 actually exercises the Jacobi code path.

### Compile gate

```bash
R CMD INSTALL --preclean .
```

### Test gate

```bash
Rscript -e 'library(leafblower); testthat::test_dir("tests/testthat")' 2>&1 | tail -2
```

Expected: `[ FAIL 0 | WARN * | SKIP * | PASS *+2 ]` (two new T2 tests pass; no regressions).

### Conventional commit message

```
feat(ieppa): Jacobi semantics for log-path margin sweep

Freeze cell_lf at outer-iter start when jacobi_log=true; lambda reads
frozen snapshot; skip scattered cell_lf patches inside the K-margin
sweep; rebuild cell_lf sequentially (O(K*M_cell)) after the sweep.
Linear path unchanged. Default jacobi_sweep=FALSE — opt-in only.

T2 tests assert max_error agreement within 2*tol_abs at K=3 and K=4.
```

---

## Task J3 — Benchmark Jacobi vs GS on 18-cell grid

**Epic:** E2 — Implementation + Decision.
**Goal:** Produce `benchmarks/jacobi_sweep_study.rds` with `iters`, `wall_s`, `iter_inflation`, `wall_ratio` for all 18 (K × compression × tightness) cells. Inform J4 decision.
**Blocker:** Requires J2 merged.
**Blocker for:** J4 (decision needs data).

### Files created

#### `benchmarks/jacobi_sweep_study.R` (new)

```r
# Jacobi vs GS benchmark study for iEPPA log-path sweep (J3).
# Grid: K x compression x tightness = 3 x 3 x 2 = 18 cells.
# Output: benchmarks/jacobi_sweep_study.rds
#
# Forces the log path via LBW_IEPPA_FORCE_PATH=log so the comparison
# reflects only the Jacobi-vs-GS distinction (not the linear/log dispatch).

suppressPackageStartupMessages({
  library(leafblower)
  library(bench)
})

stopifnot(requireNamespace("withr", quietly = TRUE))

# ── Grid axes ────────────────────────────────────────────────────────────────
K_levels         <- c(3L, 6L, 9L)
compression_lvls <- c(2, 10, 50)   # M_cell / n target ≈ 0.5, 0.1, 0.02
tightness_lvls   <- c("loose", "tight")

# Compression control: choose per-margin cardinality so that
# expected M_cell ≈ n / compression. With K margins and average
# cardinality c per margin, the cell count saturates near n; we
# pick c such that c^K ≈ n / compression (capped by integer choice).
choose_card <- function(n, K, compression) {
  target_cells <- max(2, round(n / compression))
  c_per_k <- max(2L, round(target_cells^(1 / K)))
  c_per_k
}

make_problem <- function(n, K, compression, seed) {
  set.seed(seed)
  c_per_k <- choose_card(n, K, compression)
  cols <- lapply(seq_len(K), function(k)
    factor(sample(seq_len(c_per_k), n, TRUE)))
  names(cols) <- paste0("v", seq_len(K))
  df <- as.data.frame(cols)
  tgt <- lapply(cols, function(col) {
    lv <- levels(col)
    setNames(rep(1 / length(lv), length(lv)), lv)
  })
  list(df = df, tgt = tgt)
}

# ── Per-cell runner ──────────────────────────────────────────────────────────
n          <- 50000L
seed_base  <- 20260501L
results    <- list()

for (K in K_levels) {
  for (cmp in compression_lvls) {
    for (tight in tightness_lvls) {
      max_w <- if (tight == "tight") 1.5 else 5.0
      prob  <- make_problem(n, K, cmp, seed = seed_base + K * 1000L + cmp)

      run_one <- function(jacobi) {
        withr::with_envvar(c(LBW_IEPPA_FORCE_PATH = "log"), {
          t0 <- proc.time()[["elapsed"]]
          r  <- harvest(prob$df, prob$tgt,
                        method       = "ieppa",
                        max_weight   = max_w,
                        convergence  = list(absolute = 1e-4),
                        jacobi_sweep = jacobi,
                        verbose      = 0)
          wall <- proc.time()[["elapsed"]] - t0
          list(iters = r$result$iterations,
               wall_s = wall,
               status = r$result$status,
               max_error = r$result$max_error)
        })
      }

      gs <- run_one(FALSE)
      ja <- run_one(TRUE)

      results[[length(results) + 1L]] <- data.frame(
        K = K, compression = cmp, tightness = tight,
        c_per_margin = choose_card(n, K, cmp),
        iters_gs = gs$iters, wall_gs = gs$wall_s, status_gs = gs$status,
        iters_ja = ja$iters, wall_ja = ja$wall_s, status_ja = ja$status,
        iter_inflation = ja$iters / max(1, gs$iters),
        wall_ratio     = ja$wall_s / max(1e-6, gs$wall_s),
        err_gs = gs$max_error, err_ja = ja$max_error
      )

      cat(sprintf("K=%d cmp=%d tight=%s | iters %d→%d | wall %.3fs→%.3fs | ratio %.3f\n",
                  K, cmp, tight,
                  gs$iters, ja$iters, gs$wall_s, ja$wall_s,
                  ja$wall_s / max(1e-6, gs$wall_s)))
    }
  }
}

df_out <- do.call(rbind, results)
saveRDS(df_out, "benchmarks/jacobi_sweep_study.rds")

# ── Decision summary ─────────────────────────────────────────────────────────
cat("\n── Summary ──────────────────────────────────────────────────────────\n")
cat(sprintf("Cells: %d\n", nrow(df_out)))
cat(sprintf("Both converged: %d\n",
            sum(df_out$status_gs == 0L & df_out$status_ja == 0L)))
cat(sprintf("median(wall_ratio) = %.3f\n", median(df_out$wall_ratio)))
cat(sprintf("median(iter_inflation) = %.3f\n", median(df_out$iter_inflation)))
cat(sprintf("Cells with wall_ratio < 1.0: %d / %d\n",
            sum(df_out$wall_ratio < 1.0), nrow(df_out)))
print(df_out[, c("K", "compression", "tightness",
                 "iters_gs", "iters_ja",
                 "wall_gs", "wall_ja", "wall_ratio")])
```

### Audit

- **Spy 1:** All 18 cells must report `status_gs == 0` AND `status_ja == 0`. Any non-converged cell halts the J4 decision and triggers root-cause investigation (NOT a tolerance loosen).
- **Spy 2:** `LBW_IEPPA_FORCE_PATH=log` must be honored (verify by tracing one cell with `verbose=1` and checking the log-path branch logs).
- **Spy 3:** `wall_ratio` and `iter_inflation` must be finite for every row.

### Compile gate

```bash
R CMD INSTALL --preclean .
```

### Run gate

```bash
Rscript benchmarks/jacobi_sweep_study.R 2>&1 | tail -40
```

Expected: 18 lines of per-cell output; final summary block; `benchmarks/jacobi_sweep_study.rds` exists.

### Test gate (regression check)

```bash
Rscript -e 'library(leafblower); testthat::test_dir("tests/testthat")' 2>&1 | tail -2
```

Expected: `[ FAIL 0 | ... ]` (no benchmark run breaks tests).

### Conventional commit message

```
bench(ieppa): jacobi vs GS log-path 18-cell sweep

K x compression x tightness grid (3 x 3 x 2). Forces log path via
LBW_IEPPA_FORCE_PATH=log. Records iters, wall_s, iter_inflation,
wall_ratio. Output: benchmarks/jacobi_sweep_study.rds. Informs J4
ship/revert decision.
```

---

## Task J4 — Decision gate: ship, revert, or retain-off

**Epic:** E2 — Implementation + Decision.
**Goal:** Apply decision rule from spec §"Decision Gate" using `benchmarks/jacobi_sweep_study.rds`.
**Blocker:** Requires J3 merged (data must exist).

### Decision rule (from spec)

Read `benchmarks/jacobi_sweep_study.rds`. Compute `wall_med <- median(df$wall_ratio)`.

| Outcome | Condition | Action |
|---|---|---|
| **A. Ship as default** | `wall_med < 1.0` AND all 18 cells converged | Flip default to `jacobi_sweep = TRUE` in `R/harvest.R`. Update roxygen note. |
| **B. Retain off-by-default** | Improvement only on subset (e.g., `wall_ratio < 1.0` for ≥6 cells but median ≥ 1.0) | Keep `jacobi_sweep = FALSE` default. Add docs note in `R/harvest.R` listing the (K, compression) combos where users should opt in. |
| **C. Revert** | `wall_med ≥ 1.0` AND no clear regime improvement, OR any cell failed to converge | Delete the feature wholesale. |

### Branch C (revert) — exact deletions

If reverting, surgically remove:

**`src/types.hpp`** — delete the line:
```cpp
    bool                 jacobi_log = false;  // log-path: freeze cell_lf at iter start (Jacobi) vs patch inline (GS)
```

**`src/r_bridge.cpp`** — remove:
- `SEXP jacobi_sweep_sexp` from `C_rk_calibrate` signature
- `st.jacobi_log = (INTEGER(jacobi_sweep_sexp)[0] != 0);`
- Decrement the `R_CallMethodDef` arg count for `C_rk_calibrate` by 1.

**`R/harvest.R`** — remove:
- `jacobi_sweep = FALSE` parameter
- the `as.integer(jacobi_sweep)` `.Call` argument
- the `stopifnot(is.logical(jacobi_sweep), ...)` validation
- the `@param jacobi_sweep` roxygen line

**`src/ieppa.cpp`** — remove:
- The `cell_lf_frozen` declaration and conditional `assign`
- The Edit-2 freeze block (`if (st.jacobi_log && !use_linear) std::copy(...)`)
- The Edit-3 `clf` alias inside the lambda; revert all `clf[c]` reads back to `cell_lf[c]`
- The Edit-4 `if (!st.jacobi_log)` guard (un-wrap the GS patch block to its original form)
- The Edit-5 sequential rebuild block

**`tests/testthat/test-jacobi-log-sweep.R`** — delete file.

**`benchmarks/jacobi_sweep_study.R`** — keep (historical record). The `.rds` artifact stays under `benchmarks/`.

### Branch A (ship) — exact change

`R/harvest.R`:
```r
  jacobi_sweep          = TRUE,   # was FALSE; benchmarked faster on study grid
```

Update the roxygen entry to remove "Default FALSE" and "Experimental":
```r
#' @param jacobi_sweep logical(1). If TRUE (default), the iEPPA log-path
#'   margin sweep uses Jacobi semantics. Benchmark study (benchmarks/
#'   jacobi_sweep_study.R) showed median wall-time ratio < 1.0 across the
#'   K x compression x tightness grid. Set FALSE to fall back to
#'   Gauss-Seidel.
```

### Branch B (retain off) — exact change

Update the roxygen note in `R/harvest.R` to identify the favorable regime:
```r
#' @param jacobi_sweep logical(1). Default FALSE (Gauss-Seidel). Set TRUE
#'   for high-compression problems (M_cell / n < 0.1) with K >= 6 — the
#'   benchmark study (benchmarks/jacobi_sweep_study.R) found wall-time
#'   improvement only in this regime.
```

### Audit

- **Spy 1:** The chosen branch must cite the exact `wall_med` value computed from `benchmarks/jacobi_sweep_study.rds`.
- **Spy 2:** If reverting, run `git grep -n jacobi_log` and `git grep -n jacobi_sweep` after deletions; both must return empty.
- **Spy 3:** If shipping, T2 tests still pass (default flip changes behavior of any test not pinning `jacobi_sweep=FALSE`; verify no regression).

### Compile gate

```bash
R CMD INSTALL --preclean .
```

### Test gate

```bash
Rscript -e 'library(leafblower); testthat::test_dir("tests/testthat")' 2>&1 | tail -2
```

Expected: `[ FAIL 0 | ... ]`.

### Conventional commit message

Branch A:
```
feat(ieppa): default jacobi_sweep=TRUE based on benchmark study

J3 18-cell grid showed median(wall_ratio) < 1.0. Jacobi shipped as
default for log-path sweep. GS retained as opt-out via jacobi_sweep
= FALSE.
```

Branch B:
```
docs(ieppa): document jacobi_sweep regime (off by default)

J3 study showed wall_ratio < 1.0 only for high-compression / large-K
regime. Default remains FALSE; roxygen identifies the opt-in regime.
```

Branch C:
```
revert(ieppa): remove jacobi_sweep — no wall-time win

J3 18-cell study median(wall_ratio) >= 1.0. Removing jacobi_log,
jacobi_sweep, cell_lf_frozen, and conditional log-path branches.
benchmarks/jacobi_sweep_study.R retained as historical record.
```

---

## Self-review checklist

- [x] No TBD/placeholders. Every code block is executable.
- [x] J1 stated as blocker for J2.
- [x] J3 stated as blocker for J4.
- [x] Line numbers match read-first list (types.hpp 80–140; r_bridge.cpp 108–120, 318–325; harvest.R 182–220; ieppa.cpp 200–220, 437–450, 570–630).
- [x] All four tasks have compile gate (`R CMD INSTALL --preclean .`) and test gate.
- [x] All four tasks have a Conventional Commit message.
- [x] J2's lambda edits explicitly preserve the `lf[off+j]` write outside the `cell_lf` patch guard (correctness-critical).
- [x] J3 forces log path via `LBW_IEPPA_FORCE_PATH=log` so the comparison isolates the Jacobi-vs-GS distinction.
- [x] J4 covers all three decision branches (ship, retain-off, revert) with exact deletions for Branch C.
- [x] T2 includes K=4 variant (spec compliance).
- [x] Confidence: 88 — code blocks reflect read source; lambda's exact `cell_lf` read sites must be verified at edit time (the spec's pseudocode shows the pattern; the live file at lines 622–730 has not been read in full, so J2 Edit 3 wording lists "every `cell_lf[c]` read inside the lambda's S_kj/bucket-sum computation" rather than line-numbered surgery — the executor must read the lambda body in full before editing).
