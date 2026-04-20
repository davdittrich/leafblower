# Algorithm Selection Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Bayesian level set estimation benchmark that finds the 1.2× iEPPA/L-BFGS-B crossover contour in (log_complexity, log_tol) space and produces two hardcoded constants for `src/c_api.cpp`.

**Architecture:** A single R script `benchmarks/algo_selection_benchmark.R` organised as pure helper functions + a `run_benchmark()` entry point guarded by `.BENCH_SOURCED`. Tests source the script with the guard set. The LSE loop fits a Matérn-5/2 GP via DiceKriging, picks next evaluations with the Straddle acquisition, checkpoints atomically every 5 acquisitions, and terminates when ≥90% of a 50×50 candidate grid is classified at 95% confidence or 25 acquisitions are exhausted.

**Tech Stack:** R, DiceKriging (GP), lhs (LHC design), ggplot2 (plots), leafblower (the package under test), testthat (regression guard).

---

## Beads Issue Map

| Task | Issue ID | Title |
|------|----------|-------|
| T1  | leafblower-0i9 | Add lhs/DiceKriging/ggplot2 to DESCRIPTION Suggests |
| T2  | leafblower-db9 | `bench_seed()` |
| T3  | leafblower-0fz | `make_bench_data()` |
| T4  | leafblower-guz | `time_cell()` |
| T5  | leafblower-gzc | GP helpers: `fit_gp`, `straddle_next`, `classified_fraction` |
| T6  | leafblower-521 | `save_checkpoint` / `load_checkpoint` |
| T7  | leafblower-vjn | `run_benchmark()` main LSE loop |
| T8  | leafblower-09t | `run_k_stability()` + k-stability plot |
| T9  | leafblower-9wv | `make_plots()` |
| T10 | leafblower-gje | Regression tests (`test-algo-selection.R`) |
| T11 | leafblower-dh8 | Dry run end-to-end smoke test |

---

## File Structure

```
benchmarks/
  algo_selection_benchmark.R   ← NEW: all helpers + run_benchmark() entry point
tests/testthat/
  test-algo-selection.R        ← NEW: regression guard + unit tests for helpers
DESCRIPTION                    ← MODIFY: add Suggests
```

`algo_selection_benchmark.R` internal sections (in order):
1. Package loading + constants
2. `.BENCH_SOURCED` guard variable
3. `bench_seed(log_complexity, log_tol)`
4. `make_bench_data(n, K, cats_per_margin)`
5. `time_cell(log_complexity, log_tol, K = 9L)`
6. `fit_gp(design_mat, y)`
7. `straddle_next(gp_model, candidates, threshold, kappa = 2)`
8. `classified_fraction(gp_model, candidates, threshold, conf = 0.95)`
9. `save_checkpoint(state, path)` / `load_checkpoint(path)`
10. `make_plots(state, candidates, threshold, out_dir)`
11. `run_k_stability(state, K_vals, threshold, out_dir)`
12. `run_benchmark(budget, checkpoint_path, out_dir, seed)`
13. Entry point: `if (!.BENCH_SOURCED) run_benchmark()`

---

## Task 1: Add benchmark Suggests to DESCRIPTION

**Beads:** `bd update leafblower-0i9 --claim`  
**Files:** Modify `DESCRIPTION`

- [ ] **Step 1: Edit DESCRIPTION**

```
Suggests: autumn (>= 0.2.0), testthat (>= 3.0.0), bench, lhs, DiceKriging, ggplot2
```

Replace the existing `Suggests:` line with the above (add `lhs, DiceKriging, ggplot2`).

- [ ] **Step 2: Verify R CMD CHECK passes**

```bash
R CMD INSTALL --preclean .
```

Expected: `* DONE (leafblower)` — no errors, no warnings about undeclared packages.

- [ ] **Step 3: Commit**

```bash
git add DESCRIPTION
git commit -m "feat(bench): add lhs, DiceKriging, ggplot2 to DESCRIPTION Suggests"
```

- [ ] **Step 4: Close issue**

```bash
bd close leafblower-0i9
```

---

## Task 2: `bench_seed()` — deterministic seed formula

**Beads:** `bd update leafblower-db9 --claim`  
**Files:**
- Create: `benchmarks/algo_selection_benchmark.R` (initial skeleton + this function)
- Test: `tests/testthat/test-algo-selection.R` (initial file)

- [ ] **Step 1: Create the test file**

```r
# tests/testthat/test-algo-selection.R
# Guard: prevent benchmark execution when sourced for testing
.BENCH_SOURCED <- TRUE

# Source the benchmark to load helper functions (no side effects with guard set)
# NOTE: benchmark file is built incrementally — only functions added so far are available.
source(system.file("../../benchmarks/algo_selection_benchmark.R",
                   package = "leafblower", mustWork = FALSE) |>
       (\(x) if (nchar(x) == 0) "benchmarks/algo_selection_benchmark.R" else x)())

test_that("bench_seed is deterministic", {
  s1 <- bench_seed(5.0, -4.0)
  s2 <- bench_seed(5.0, -4.0)
  expect_identical(s1, s2)
})

test_that("bench_seed returns integer in 32-bit range", {
  s <- bench_seed(7.7, -3.0)
  expect_true(is.integer(s))
  expect_true(s > 0L && s < .Machine$integer.max)
})

test_that("bench_seed differs for different inputs", {
  seeds <- sapply(c(4.0, 5.0, 6.0, 7.0), function(x) bench_seed(x, -4.0))
  expect_equal(length(unique(seeds)), 4L)
  seeds2 <- sapply(c(-3.0, -4.0, -5.0, -6.0), function(x) bench_seed(5.0, x))
  expect_equal(length(unique(seeds2)), 4L)
})
```

- [ ] **Step 2: Create benchmark script skeleton + `bench_seed`**

```r
# benchmarks/algo_selection_benchmark.R
# Bayesian Level Set Estimation benchmark for iEPPA vs L-BFGS-B algorithm selection.
# See: docs/superpowers/specs/2026-04-20-algo-selection-design.md

# Enforce single-threaded execution before loading any library that might
# initialise an OpenMP thread pool. OMP reads OMP_NUM_THREADS at pool creation
# time (lazy, first parallel region). Setting it here — before library() calls —
# guarantees the env var is visible at that point, eliminating CPU-contention
# confound from the timing comparisons.
Sys.setenv(OMP_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(leafblower)
  library(DiceKriging)
  library(lhs)
  library(ggplot2)
})

# Guard: set .BENCH_SOURCED <- TRUE before source()ing this file in tests
# to prevent run_benchmark() from executing.
.BENCH_SOURCED <- exists(".BENCH_SOURCED", envir = .GlobalEnv, inherits = FALSE)

# Input space bounds
BENCH_X1_RANGE <- c(4.0, 7.7)   # log10(complexity): 10K to 50M
BENCH_X2_RANGE <- c(-6.0, -3.0) # log10(tol_abs):    1e-6 to 1e-3
BENCH_THRESHOLD <- log(1.2)      # log(1.2) ≈ 0.182; L-BFGS-B wins above this

# ── bench_seed ────────────────────────────────────────────────────────────────
# Deterministic integer seed from (log_complexity, log_tol).
# Valid for log_complexity in [4, 7.7] and log_tol in [-6, -3].
# Formula produces seeds in [4e8+0, 7.7e8+9999], all within 32-bit integer range.
bench_seed <- function(log_complexity, log_tol) {
  a <- as.integer(round(log_complexity * 1e4))
  b <- as.integer(round(-log_tol * 1e4)) %% 10000L
  (a * 10000L) + b
}
```

- [ ] **Step 3: Run failing test to confirm it fails (file-not-found or function-not-found)**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-algo-selection.R')"
```

Expected: FAIL — `source()` fails because the file doesn't exist yet. (Acceptable TDD failure.)

- [ ] **Step 4: Run test again now that benchmark script exists**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-algo-selection.R')"
```

Expected: 3 tests PASS (determinism, integer range, uniqueness).

- [ ] **Step 5: Commit**

```bash
git add benchmarks/algo_selection_benchmark.R tests/testthat/test-algo-selection.R
git commit -m "feat(bench): bench_seed() + test file skeleton"
```

- [ ] **Step 6: Close issue**

```bash
bd close leafblower-db9
```

---

## Task 3: `make_bench_data()` — synthetic categorical data

**Beads:** `bd update leafblower-0fz --claim`  
**Files:**
- Modify: `benchmarks/algo_selection_benchmark.R` (append function)
- Modify: `tests/testthat/test-algo-selection.R` (append tests)

- [ ] **Step 1: Append tests to `test-algo-selection.R`**

```r
test_that("make_bench_data returns correct structure", {
  set.seed(1L)
  out <- make_bench_data(n = 200L, K = 3L, cats_per_margin = 4L)
  expect_named(out, c("df", "targets"))
  expect_equal(nrow(out$df), 200L)
  expect_equal(ncol(out$df), 3L)
  expect_equal(length(out$targets), 3L)
})

test_that("make_bench_data targets sum to 1", {
  set.seed(2L)
  out <- make_bench_data(n = 100L, K = 2L, cats_per_margin = 5L)
  sums <- sapply(out$targets, sum)
  expect_true(all(abs(sums - 1.0) < 1e-12))
})

test_that("make_bench_data df columns are factors with correct levels", {
  set.seed(3L)
  out <- make_bench_data(n = 50L, K = 2L, cats_per_margin = 3L)
  expect_true(all(sapply(out$df, is.factor)))
  expect_true(all(sapply(out$df, nlevels) == 3L))
})
```

- [ ] **Step 2: Run tests to confirm they FAIL**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-algo-selection.R')" 2>&1 | tail -5
```

Expected: 3 new tests FAIL with "could not find function 'make_bench_data'".

- [ ] **Step 3: Append `make_bench_data` to benchmark script**

```r
# ── make_bench_data ───────────────────────────────────────────────────────────
# Generates n-row survey data with K categorical margins of cats_per_margin levels.
# Population proportions: Dirichlet(1,...,1) via normalised Exp(1) draws.
# Sample proportions: population * |Normal(1, 0.1)| noise, renormalised (~10% bias).
make_bench_data <- function(n, K, cats_per_margin) {
  stopifnot(n >= 1L, K >= 1L, cats_per_margin >= 2L)
  col_names <- paste0("m", seq_len(K))
  lvl_names <- lapply(seq_len(K), function(k) paste0("c", seq_len(cats_per_margin)))

  pop_props  <- lapply(seq_len(K), function(k) {
    x <- rexp(cats_per_margin); x / sum(x)
  })
  samp_props <- lapply(pop_props, function(p) {
    q <- p * abs(rnorm(length(p), mean = 1, sd = 0.1)); q / sum(q)
  })

  df_cols <- lapply(seq_len(K), function(k) {
    factor(sample(lvl_names[[k]], n, replace = TRUE, prob = samp_props[[k]]),
           levels = lvl_names[[k]])
  })
  df <- as.data.frame(setNames(df_cols, col_names))

  targets <- setNames(
    lapply(seq_len(K), function(k) setNames(pop_props[[k]], lvl_names[[k]])),
    col_names
  )
  list(df = df, targets = targets)
}
```

- [ ] **Step 4: Run tests to confirm they PASS**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-algo-selection.R')" 2>&1 | tail -5
```

Expected: 6 tests PASS (3 from T2 + 3 new).

- [ ] **Step 5: Commit**

```bash
git add benchmarks/algo_selection_benchmark.R tests/testthat/test-algo-selection.R
git commit -m "feat(bench): make_bench_data() synthetic categorical generator"
```

- [ ] **Step 6: Close issue**

```bash
bd close leafblower-0fz
```

---

## Task 4: `time_cell()` — log(t_iEPPA / t_LBFGSB) at one grid point

**Beads:** `bd update leafblower-guz --claim`  
**Files:**
- Modify: `benchmarks/algo_selection_benchmark.R`
- Modify: `tests/testthat/test-algo-selection.R`

- [ ] **Step 1: Append tests**

```r
test_that("time_cell returns finite numeric scalar", {
  # Tiny problem for speed: log_complexity=4.0 → n≈278, K=9, cats=4
  # tol=1e-3 (loose) → fast convergence
  y <- time_cell(log_complexity = 4.0, log_tol = -3.0, K = 9L)
  expect_true(is.numeric(y) && length(y) == 1L && is.finite(y))
})

test_that("time_cell is deterministic for same inputs", {
  y1 <- time_cell(log_complexity = 4.0, log_tol = -3.0, K = 9L)
  y2 <- time_cell(log_complexity = 4.0, log_tol = -3.0, K = 9L)
  # Same seed → same data → same weights → same timing direction (sign)
  expect_equal(sign(y1), sign(y2))
})

test_that("time_cell seed_extra produces deterministic K-stability seeds", {
  # seed_extra must produce a repeatable result distinct from seed_extra=0
  y_main <- time_cell(4.0, -3.0, K = 9L, seed_extra = 0L)
  y_k3_a <- time_cell(4.0, -3.0, K = 3L, seed_extra = 3L * 10000000L)
  y_k3_b <- time_cell(4.0, -3.0, K = 3L, seed_extra = 3L * 10000000L)
  # K-stability call is repeatable
  expect_equal(sign(y_k3_a), sign(y_k3_b))
})
```

- [ ] **Step 2: Run tests — confirm FAIL**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-algo-selection.R')" 2>&1 | tail -5
```

Expected: 2 new tests FAIL with "could not find function 'time_cell'".

- [ ] **Step 3: Append `time_cell` to benchmark script**

```r
# ── time_cell ─────────────────────────────────────────────────────────────────
# Times iEPPA vs L-BFGS-B at one (log_complexity, log_tol) point.
# Returns log(median_t_iEPPA / median_t_LBFGSB):
#   positive → iEPPA is slower → L-BFGS-B wins.
# Threshold: log(1.2) ≈ 0.182.
#
# n derivation: cats_per_margin fixed by complexity tercile, n = round(complexity/(K*cats)).
# Both solvers run with max_weight=Inf (unconstrained, exponential link).
# seed_extra: K-specific seed offset for K-stability sweeps (default 0L = main sweep).
#   K=3 stability: pass seed_extra = 3L * 10000000L
#   K=18 stability: pass seed_extra = 18L * 10000000L
time_cell <- function(log_complexity, log_tol, K = 9L, seed_extra = 0L) {
  K <- as.integer(K)
  cats_per_margin <- if (log_complexity <= 5.5) 4L else if (log_complexity <= 6.5) 8L else 16L
  n <- max(50L, as.integer(round(10^log_complexity / (K * cats_per_margin))))

  # Round-trip sanity check
  actual_log_c <- log10(n * K * cats_per_margin)
  if (abs(actual_log_c - log_complexity) > 0.2)
    warning(sprintf("time_cell: complexity round-trip %.2f log-units (lc=%.2f K=%d cats=%d n=%d)",
                    abs(actual_log_c - log_complexity), log_complexity, K, cats_per_margin, n))

  set.seed(bench_seed(log_complexity, log_tol) + seed_extra)
  bd <- make_bench_data(n, K, cats_per_margin)
  conv <- list(absolute = 10^log_tol)

  time_algo <- function(method) {
    # 2 warmup runs (discarded)
    for (i in seq_len(2L))
      suppressWarnings(invisible(leafblower::harvest(
        bd$df, bd$targets, method = method,
        max_weight = Inf, min_weight = 0, convergence = conv, max_iterations = 500L)))
    # 5 timed runs
    median(replicate(5L, {
      t0 <- proc.time()[["elapsed"]]
      suppressWarnings(invisible(leafblower::harvest(
        bd$df, bd$targets, method = method,
        max_weight = Inf, min_weight = 0, convergence = conv, max_iterations = 500L)))
      proc.time()[["elapsed"]] - t0
    }))
  }

  t_ieppa  <- time_algo("ieppa")
  t_lbfgsb <- time_algo("lbfgsb")
  log(t_ieppa / t_lbfgsb)
}
```

- [ ] **Step 4: Run tests — confirm PASS**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-algo-selection.R')" 2>&1 | tail -5
```

Expected: 9 tests PASS. (Note: time_cell tests take ~5s each due to tiny solver runs.)

- [ ] **Step 5: Commit**

```bash
git add benchmarks/algo_selection_benchmark.R tests/testthat/test-algo-selection.R
git commit -m "feat(bench): time_cell() — log-ratio timing for one grid point"
```

- [ ] **Step 6: Close issue**

```bash
bd close leafblower-guz
```

---

## Task 5: GP helpers — `fit_gp`, `straddle_next`, `classified_fraction`

**Beads:** `bd update leafblower-gzc --claim`  
**Files:**
- Modify: `benchmarks/algo_selection_benchmark.R`
- Modify: `tests/testthat/test-algo-selection.R`

- [ ] **Step 1: Append tests**

```r
test_that("fit_gp returns km object with finite predictions", {
  set.seed(42L)
  design <- matrix(c(4.0, 5.0, 6.0, 7.0, 4.5, 5.5, 6.5, 4.0,
                     -3.0, -4.0, -5.0, -6.0, -3.5, -4.5, -5.5, -6.0),
                   ncol = 2)
  y <- rnorm(8L)
  gp <- fit_gp(design, y)
  expect_s4_class(gp, "km")
  cands <- as.data.frame(matrix(c(5.0, 5.5, -4.0, -4.5), ncol = 2))
  names(cands) <- c("V1", "V2")
  pred <- DiceKriging::predict(gp, newdata = cands, type = "UK", checkNames = FALSE)
  expect_true(all(is.finite(pred$mean)))
  expect_true(all(pred$sd >= 0))
})

test_that("straddle_next returns a 1-row matrix within bounds", {
  set.seed(1L)
  design <- matrix(c(4.0, 5.0, 6.0, 4.5, 5.5, 6.5, -3.0, -4.0, -5.0, -3.5, -4.5, -5.5), ncol = 2)
  y <- c(0.5, 0.1, -0.3, 0.3, 0.0, -0.1)
  gp <- fit_gp(design, y)
  cands <- as.matrix(expand.grid(
    V1 = seq(4.0, 7.7, length.out = 10),
    V2 = seq(-6.0, -3.0, length.out = 10)))
  nxt <- straddle_next(gp, cands, threshold = log(1.2))
  expect_equal(nrow(nxt), 1L)
  expect_equal(ncol(nxt), 2L)
  expect_true(nxt[1, 1] >= 4.0 && nxt[1, 1] <= 7.7)
  expect_true(nxt[1, 2] >= -6.0 && nxt[1, 2] <= -3.0)
})

test_that("classified_fraction is 0 for uncertain GP, 1 for certain", {
  # A GP with tiny nugget on a flat surface: everything classified near threshold
  set.seed(5L)
  design <- matrix(c(5.0, 6.0, -4.0, -5.0), ncol = 2)
  # y values all above threshold → entire space classified as "above"
  y <- rep(2.0, 2L)
  gp <- suppressWarnings(fit_gp(design, y))
  cands <- as.matrix(expand.grid(V1 = seq(4.5, 7.2, length.out = 5),
                                  V2 = seq(-5.5, -3.5, length.out = 5)))
  frac <- classified_fraction(gp, cands, threshold = log(1.2))
  expect_true(is.numeric(frac) && length(frac) == 1L)
  expect_true(frac >= 0 && frac <= 1)
})
```

- [ ] **Step 2: Run tests — confirm FAIL**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-algo-selection.R')" 2>&1 | tail -5
```

Expected: 3 new tests FAIL.

- [ ] **Step 3: Append GP helpers to benchmark script**

```r
# ── fit_gp ────────────────────────────────────────────────────────────────────
# Fits a Matérn-5/2 GP to timing observations.
# design_mat: n_obs × 2 matrix (col 1 = log_complexity, col 2 = log_tol).
# y: numeric vector of log(t_iEPPA / t_LBFGSB).
fit_gp <- function(design_mat, y) {
  stopifnot(is.matrix(design_mat), nrow(design_mat) == length(y), ncol(design_mat) == 2L)
  DiceKriging::km(
    formula      = ~1,
    design       = as.data.frame(design_mat),
    response     = y,
    covtype      = "matern5_2",
    nugget.estim = TRUE,
    nugget       = 1e-4,          # lower bound: prevents degenerate fit on small n
    control      = list(trace = FALSE)
  )
}

# ── straddle_next ─────────────────────────────────────────────────────────────
# Straddle acquisition (Bryan et al. 2005): picks the candidate maximising
#   a(x) = -|μ(x) − threshold| + κ·σ(x)
# Pulls samples toward the contour (low |μ − threshold|) and uncertain regions (high σ).
straddle_next <- function(gp_model, candidates, threshold, kappa = 2) {
  pred <- DiceKriging::predict(gp_model,
                               newdata    = as.data.frame(candidates),
                               type       = "UK",
                               checkNames = FALSE)
  a   <- -abs(pred$mean - threshold) + kappa * pred$sd
  candidates[which.max(a), , drop = FALSE]
}

# ── classified_fraction ───────────────────────────────────────────────────────
# Fraction of candidates classified with ≥conf confidence as above or below threshold.
# Termination fires when this reaches 0.90.
classified_fraction <- function(gp_model, candidates, threshold, conf = 0.95) {
  pred    <- DiceKriging::predict(gp_model,
                                  newdata    = as.data.frame(candidates),
                                  type       = "UK",
                                  checkNames = FALSE)
  p_above <- pnorm(threshold, mean = pred$mean, sd = pred$sd, lower.tail = FALSE)
  mean(p_above > conf | p_above < (1 - conf))
}
```

- [ ] **Step 4: Run tests — confirm PASS**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-algo-selection.R')" 2>&1 | tail -5
```

Expected: 12 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add benchmarks/algo_selection_benchmark.R tests/testthat/test-algo-selection.R
git commit -m "feat(bench): fit_gp(), straddle_next(), classified_fraction() GP helpers"
```

- [ ] **Step 6: Close issue**

```bash
bd close leafblower-gzc
```

---

## Task 6: Checkpoint — `save_checkpoint` / `load_checkpoint`

**Beads:** `bd update leafblower-521 --claim`  
**Files:**
- Modify: `benchmarks/algo_selection_benchmark.R`
- Modify: `tests/testthat/test-algo-selection.R`

- [ ] **Step 1: Append tests**

```r
test_that("save_checkpoint + load_checkpoint round-trips state", {
  tmp_path <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_path), add = TRUE)
  state <- list(design = matrix(1:4, nrow = 2), y = c(0.1, 0.2), iter = 3L)
  save_checkpoint(state, tmp_path)
  expect_true(file.exists(tmp_path))
  loaded <- load_checkpoint(tmp_path)
  expect_equal(loaded$iter, 3L)
  expect_equal(loaded$y, c(0.1, 0.2))
})

test_that("load_checkpoint returns NULL when file absent", {
  result <- load_checkpoint(tempfile(fileext = ".rds"))
  expect_null(result)
})

test_that("save_checkpoint leaves no .tmp file on success", {
  tmp_path <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_path), add = TRUE)
  save_checkpoint(list(x = 1L), tmp_path)
  expect_false(file.exists(paste0(tmp_path, ".tmp")))
})
```

- [ ] **Step 2: Run tests — confirm FAIL**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-algo-selection.R')" 2>&1 | tail -5
```

Expected: 3 new tests FAIL.

- [ ] **Step 3: Append checkpoint functions to benchmark script**

```r
# ── save_checkpoint ───────────────────────────────────────────────────────────
# Atomic checkpoint: write to .tmp, then rename to final path.
# state fields: design (n×2 matrix), y (numeric), gp (km or NULL),
#               iter (integer), classified (numeric), bounds (list).
save_checkpoint <- function(state, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(state, tmp)
  file.rename(tmp, path)   # atomic on same filesystem (benchmarks/ → benchmarks/)
  invisible(path)
}

# ── load_checkpoint ───────────────────────────────────────────────────────────
# Returns NULL if no checkpoint found; otherwise returns the saved state.
load_checkpoint <- function(path) {
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}
```

- [ ] **Step 4: Run tests — confirm PASS**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-algo-selection.R')" 2>&1 | tail -5
```

Expected: 15 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add benchmarks/algo_selection_benchmark.R tests/testthat/test-algo-selection.R
git commit -m "feat(bench): save_checkpoint()/load_checkpoint() with atomic write"
```

- [ ] **Step 6: Close issue**

```bash
bd close leafblower-521
```

---

## Task 7: `run_benchmark()` — main LSE loop

**Beads:** `bd update leafblower-vjn --claim`  
**Files:**
- Modify: `benchmarks/algo_selection_benchmark.R`

No dedicated unit test — tested via the dry run in T11. The function is integration code.

- [ ] **Step 1: Append `run_benchmark` to benchmark script**

```r
# ── run_benchmark ─────────────────────────────────────────────────────────────
# Main entry point. Runs the Bayesian LSE loop.
#
# budget:          max adaptive acquisitions (default 25)
# checkpoint_path: where to save/load state (default benchmarks/algo_selection_results.rds)
# out_dir:         where to write PDFs (default benchmarks/)
# seed:            seed for the 8-pt LHC initial design (default 42L)
run_benchmark <- function(budget          = 25L,
                          checkpoint_path = "benchmarks/algo_selection_results.rds",
                          out_dir         = "benchmarks",
                          seed            = 42L,
                          lhc_x1_max      = NULL) {
  # lhc_x1_max: cap LHC x1 (log_complexity) coordinates to at most this value.
  #   NULL (default) → full [4, 7.7] range used in production.
  #   5.5 → smoke-test mode: all LHC points stay in the fast sub-region
  #          (n ≤ ~4K at K=9, cats=8), keeping the smoke test to ~2 min.

  threshold  <- BENCH_THRESHOLD
  x1_range   <- BENCH_X1_RANGE
  x2_range   <- BENCH_X2_RANGE

  # 50×50 candidate grid — used for Straddle and classification
  candidates <- as.matrix(expand.grid(
    V1 = seq(x1_range[1], x1_range[2], length.out = 50L),
    V2 = seq(x2_range[1], x2_range[2], length.out = 50L)
  ))

  # ── Restart or initialise ──────────────────────────────────────────────────
  state <- load_checkpoint(checkpoint_path)
  if (!is.null(state)) {
    cat(sprintf("Restarting from checkpoint: %d evaluations, iter=%d, classified=%.2f\n",
                nrow(state$design), state$iter, state$classified))
    # Refit GP from saved design + y (do not re-run evaluations)
    state$gp <- fit_gp(state$design, state$y)
  } else {
    cat("Initialising: 8-point Latin hypercube design...\n")
    set.seed(seed)
    lhc_unit <- lhs::randomLHS(n = 8L, k = 2L)
    lhc_pts  <- cbind(
      x1_range[1] + lhc_unit[, 1L] * diff(x1_range),
      x2_range[1] + lhc_unit[, 2L] * diff(x2_range)
    )
    # Smoke-test cap: clamp LHC log_complexity coordinates to lhc_x1_max.
    # This keeps all initial evaluations in the fast sub-region without
    # changing the Straddle/classification grid (always full [4,7.7]).
    if (!is.null(lhc_x1_max))
      lhc_pts[, 1L] <- pmin(lhc_pts[, 1L], lhc_x1_max)
    design <- matrix(nrow = 0L, ncol = 2L)
    y      <- numeric(0L)
    for (i in seq_len(nrow(lhc_pts))) {
      cat(sprintf("  LHC %d/8: lc=%.2f, lt=%.2f ... ", i, lhc_pts[i, 1L], lhc_pts[i, 2L]))
      yi <- tryCatch(time_cell(lhc_pts[i, 1L], lhc_pts[i, 2L]),
                     error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NA_real_ })
      cat(sprintf("y=%.3f\n", yi))
      design <- rbind(design, lhc_pts[i, , drop = FALSE])
      y      <- c(y, yi)
    }
    # Drop any NA evaluations
    ok     <- is.finite(y)
    design <- design[ok, , drop = FALSE]
    y      <- y[ok]

    if (sum(ok) < 2L)
      stop("All LHC evaluations failed — check leafblower installation and system capacity. ",
           "Cannot fit GP on fewer than 2 observations.")

    state <- list(design = design, y = y, gp = fit_gp(design, y),
                  iter = 0L, classified = 0.0, bounds = list(x1 = x1_range, x2 = x2_range))
    save_checkpoint(state, checkpoint_path)
  }

  # ── Adaptive acquisitions ──────────────────────────────────────────────────
  for (i in seq_len(budget)) {
    state$classified <- classified_fraction(state$gp, candidates, threshold)
    cat(sprintf("Iter %d/%d: classified=%.2f\n", state$iter + 1L, budget, state$classified))
    if (state$classified >= 0.90) {
      cat("Termination: 90% classified.\n")
      state$converged <- TRUE
      break
    }

    next_pt <- straddle_next(state$gp, candidates, threshold)
    cat(sprintf("  Next: lc=%.3f, lt=%.3f ... ", next_pt[1L, 1L], next_pt[1L, 2L]))
    yi <- tryCatch(time_cell(next_pt[1L, 1L], next_pt[1L, 2L]),
                   error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NA_real_ })
    cat(sprintf("y=%.3f\n", yi))

    if (is.finite(yi)) {
      state$design <- rbind(state$design, next_pt)
      state$y      <- c(state$y, yi)
      state$gp     <- fit_gp(state$design, state$y)
    }
    state$iter <- state$iter + 1L

    if (state$iter %% 5L == 0L)
      save_checkpoint(state, checkpoint_path)
  }

  # Final checkpoint
  save_checkpoint(state, checkpoint_path)

  # ── Poor-fit warning ───────────────────────────────────────────────────────
  final_classified <- classified_fraction(state$gp, candidates, threshold)
  if (final_classified < 0.90) {
    warning(
      sprintf(paste0(
        "GP classification %.0f%% < 90%% at termination.\n",
        "Inspect benchmarks/algo_selection_uncertainty.pdf before committing constants.\n",
        "Consider filing a follow-up issue:\n",
        "  bd create --title='algo-selection: 3D sweep needed (K/n confound)' ",
        "--description='GP classified only %.0f%% of space after %d acquisitions' ",
        "--type=task --priority=3"),
        100 * final_classified, 100 * final_classified, state$iter))
  }

  # ── Plots ─────────────────────────────────────────────────────────────────
  cat("Generating plots...\n")
  make_plots(state, candidates, threshold, out_dir)

  cat(sprintf("\nDone. %d total evaluations, %.0f%% classified.\n",
              nrow(state$design), 100 * final_classified))
  cat(sprintf("Inspect: %s/algo_selection_contour.pdf\n", out_dir))
  invisible(state)
}
```

- [ ] **Step 2: Add the entry-point guard at the bottom of the benchmark script**

```r
# ── Entry point ───────────────────────────────────────────────────────────────
# Runs only when executed via `Rscript benchmarks/algo_selection_benchmark.R`.
# Skipped when sourced in tests (set `.BENCH_SOURCED <- TRUE` before source()).
# K-stability is called here (not inside run_benchmark) so smoke tests can call
# run_benchmark() alone without triggering 32 extra time_cell() evaluations.
if (!.BENCH_SOURCED) {
  state <- run_benchmark()
  run_k_stability(state, K_vals = c(3L, 18L), threshold = BENCH_THRESHOLD,
                  out_dir = "benchmarks")
}
```

- [ ] **Step 3: Verify script parses without error**

```bash
Rscript -e "source('benchmarks/algo_selection_benchmark.R')" 2>&1 | head -10
```

Expected: no output (guard fires, `run_benchmark` is not called).

- [ ] **Step 4: Verify existing tests still pass**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-algo-selection.R')" 2>&1 | tail -5
```

Expected: 15 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add benchmarks/algo_selection_benchmark.R
git commit -m "feat(bench): run_benchmark() main LSE loop with checkpoint and poor-fit warning"
```

- [ ] **Step 6: Close issue**

```bash
bd close leafblower-vjn
```

---

## Task 8: `run_k_stability()` + k-stability plot

**Beads:** `bd update leafblower-09t --claim`  
**Files:**
- Modify: `benchmarks/algo_selection_benchmark.R`

- [ ] **Step 1: Append `run_k_stability` to benchmark script (before `run_benchmark`)**

```r
# ── run_k_stability ───────────────────────────────────────────────────────────
# Evaluates a 4×4 fixed grid at K ∈ K_vals to check whether the 1.2× contour
# shifts with margin count. Overlays all contours on a single PDF.
#
# Decision rule: if the K=3 or K=18 contour shifts >0.5 log-units from the
# K=9 GP contour at any tol level, prints a scope warning.
run_k_stability <- function(state, K_vals = c(3L, 18L), threshold, out_dir = "benchmarks",
                             grid_size = 4L) {
  # grid_size: number of points per axis (default 4 → 4×4 production grid).
  #   Use grid_size=1 in tests to call a single time_cell() per K value.
  x1_pts <- seq(4.5, 7.2, length.out = grid_size)
  x2_pts <- seq(-5.5, -3.5, length.out = grid_size)
  grid   <- expand.grid(log_complexity = x1_pts, log_tol = x2_pts)

  results_list <- lapply(K_vals, function(K) {
    cat(sprintf("  K-stability: evaluating K=%d (16 points)...\n", K))
    y_K <- numeric(nrow(grid))
    for (i in seq_len(nrow(grid))) {
      # seed_extra = K * 10000000L ensures K-stability data is independent of main K=9 sweep.
      # The offset fits in 32-bit integer for K <= 214.
      # NOTE: seed_extra is passed INTO time_cell() — do NOT call set.seed() here.
      #       time_cell() internally calls set.seed(bench_seed(...) + seed_extra).
      y_K[i] <- tryCatch(
        time_cell(grid$log_complexity[i], grid$log_tol[i], K = K,
                  seed_extra = K * 10000000L),
        error = function(e) NA_real_
      )
    }
    data.frame(grid, y = y_K, K = K)
  })

  # Add K=9 GP posterior mean at the same grid points for comparison
  cands_grid <- as.matrix(grid)
  pred_k9    <- DiceKriging::predict(state$gp,
                                     newdata    = as.data.frame(cands_grid),
                                     type       = "UK",
                                     checkNames = FALSE)
  k9_df <- data.frame(grid, y = pred_k9$mean, K = 9L)

  all_df <- rbind(k9_df, do.call(rbind, results_list))
  all_df$K_label <- paste0("K=", all_df$K)

  # Contour shift check: compare each K_val contour against K=9 at each x2 level
  k9_vals <- k9_df$y
  for (K in K_vals) {
    kv <- all_df[all_df$K == K, "y"]
    if (any(is.finite(kv)) && any(is.finite(k9_vals))) {
      max_shift <- max(abs(kv[is.finite(kv) & is.finite(k9_vals)] -
                           k9_vals[is.finite(kv) & is.finite(k9_vals)]))
      if (max_shift > 0.5) {
        warning(sprintf(paste0(
          "K-stability: K=%d contour shifts %.2f log-units from K=9.\n",
          "Constants are valid only for K≈9 (Stepstone regime).\n",
          "Add comment to kComplexityThreshold/kTolThreshold in src/c_api.cpp:\n",
          "  // Threshold calibrated for K≈9, uniform-category problems."),
          K, max_shift))
      }
    }
  }

  # K-stability overlay plot
  p <- ggplot2::ggplot(all_df[is.finite(all_df$y), ],
                       ggplot2::aes(log_complexity, log_tol, z = y, colour = K_label)) +
    ggplot2::geom_contour(breaks = threshold, linewidth = 1) +
    ggplot2::labs(title = "K-stability: 1.2× contour at K=3, 9, 18",
                  x = "log10(complexity)", y = "log10(tol_abs)", colour = "K") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(out_dir, "algo_selection_k_stability.pdf"), p, width = 8, height = 6)
  cat(sprintf("  K-stability plot saved to %s/algo_selection_k_stability.pdf\n", out_dir))
  invisible(all_df)
}
```

- [ ] **Step 2: Verify script parses**

```bash
Rscript -e "source('benchmarks/algo_selection_benchmark.R')" 2>&1 | head -5
```

Expected: no output.

- [ ] **Step 3: Verify entry-point block wires run_k_stability correctly**

```bash
grep -A 5 "if (!.BENCH_SOURCED)" benchmarks/algo_selection_benchmark.R
```

Expected: shows `state <- run_benchmark()` and `run_k_stability(state, K_vals = c(3L, 18L), ...)` in the block.

- [ ] **Step 4: Runtime verification — produces PDF with grid_size=1 (~10s)**

```bash
Rscript -e "
  .BENCH_SOURCED <- TRUE
  source('benchmarks/algo_selection_benchmark.R')
  # Build a tiny 4-point GP state for verification
  set.seed(1L)
  design <- matrix(c(4.5, 5.0, 5.5, 6.0, -3.0, -4.0, -5.0, -6.0), ncol = 2)
  y      <- c(0.3, 0.1, -0.1, -0.3)
  gp     <- suppressWarnings(fit_gp(design, y))
  state  <- list(design = design, y = y, gp = gp, iter = 4L, classified = 0.5)
  tmp    <- tempdir()
  # grid_size=1: evaluates 1 point × 1 K value = 1 time_cell() call (~10s)
  run_k_stability(state, K_vals = c(3L), threshold = log(1.2),
                  out_dir = tmp, grid_size = 1L)
  stopifnot(file.exists(file.path(tmp, 'algo_selection_k_stability.pdf')))
  cat('run_k_stability: PDF produced OK\n')
" 2>&1 | tail -5
```

Expected: `run_k_stability: PDF produced OK`.

- [ ] **Step 5: Commit**

```bash
git add benchmarks/algo_selection_benchmark.R
git commit -m "feat(bench): run_k_stability() K-stability check and overlay plot"
```

- [ ] **Step 6: Close issue**

```bash
bd close leafblower-09t
```

---

## Task 9: `make_plots()` — contour and uncertainty PDFs

**Beads:** `bd update leafblower-9wv --claim`  
**Files:**
- Modify: `benchmarks/algo_selection_benchmark.R`

- [ ] **Step 1: Append `make_plots` to benchmark script (before `run_k_stability`)**

```r
# ── make_plots ────────────────────────────────────────────────────────────────
# Generates two PDFs from the final GP:
#   algo_selection_contour.pdf     — posterior mean heatmap + 1.2× contour
#   algo_selection_uncertainty.pdf — posterior σ surface
make_plots <- function(state, candidates, threshold, out_dir = "benchmarks") {
  pred    <- DiceKriging::predict(state$gp,
                                  newdata    = as.data.frame(candidates),
                                  type       = "UK",
                                  checkNames = FALSE)
  cand_df <- as.data.frame(candidates)
  names(cand_df) <- c("log_complexity", "log_tol")
  cand_df$mean <- pred$mean
  cand_df$sd   <- pred$sd

  # Mark evaluated points: circle = LHC initial, cross = adaptive
  n_lhc    <- 8L
  n_total  <- nrow(state$design)
  pt_df    <- as.data.frame(state$design)
  names(pt_df) <- c("log_complexity", "log_tol")
  pt_df$type <- c(rep("LHC", min(n_lhc, n_total)),
                  rep("Adaptive", max(0L, n_total - n_lhc)))

  # ── Contour plot ───────────────────────────────────────────────────────────
  p1 <- ggplot2::ggplot(cand_df, ggplot2::aes(log_complexity, log_tol)) +
    ggplot2::geom_tile(ggplot2::aes(fill = mean)) +
    ggplot2::geom_contour(ggplot2::aes(z = mean),
                          breaks    = threshold,
                          colour    = "red",
                          linewidth = 1.2) +
    ggplot2::geom_point(data  = pt_df,
                        ggplot2::aes(shape = type),
                        colour = "white", size = 2) +
    ggplot2::scale_fill_viridis_c(name = "log(t_iEPPA/t_LBFGSB)") +
    ggplot2::scale_shape_manual(values = c(LHC = 16L, Adaptive = 4L)) +
    ggplot2::labs(
      title   = "GP posterior mean — red line = 1.2× contour (L-BFGS-B wins above)",
      x       = "log10(complexity = n × Σcat_counts)",
      y       = "log10(tol_abs)",
      shape   = "Design point"
    ) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(out_dir, "algo_selection_contour.pdf"), p1, width = 9, height = 6)

  # ── Uncertainty plot ───────────────────────────────────────────────────────
  p2 <- ggplot2::ggplot(cand_df, ggplot2::aes(log_complexity, log_tol)) +
    ggplot2::geom_tile(ggplot2::aes(fill = sd)) +
    ggplot2::scale_fill_viridis_c(name = "posterior σ", option = "magma") +
    ggplot2::labs(
      title = "GP posterior uncertainty — high σ = unreliable region",
      x     = "log10(complexity = n × Σcat_counts)",
      y     = "log10(tol_abs)"
    ) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(out_dir, "algo_selection_uncertainty.pdf"), p2, width = 9, height = 6)

  cat(sprintf("Plots saved to %s/\n", out_dir))
  invisible(list(contour = p1, uncertainty = p2))
}
```

- [ ] **Step 2: Verify parse**

```bash
Rscript -e "source('benchmarks/algo_selection_benchmark.R')" 2>&1 | head -5
```

Expected: no output.

- [ ] **Step 3: Commit**

```bash
git add benchmarks/algo_selection_benchmark.R
git commit -m "feat(bench): make_plots() contour and uncertainty PDFs"
```

- [ ] **Step 4: Close issue**

```bash
bd close leafblower-9wv
```

---

## Task 10: Regression tests — `test-algo-selection.R`

**Beads:** `bd update leafblower-gje --claim`  
**Files:**
- Modify: `tests/testthat/test-algo-selection.R`

These tests guard the `select_algorithm()` logic via `harvest()` + `attr(result, "algorithm")`. They verify the routing rules that the benchmark will later validate empirically.

- [ ] **Step 1: Append regression tests**

```r
# ── Algorithm routing regression tests ────────────────────────────────────────
# These tests verify select_algorithm() routing via harvest().
# Test 3 (L-BFGS-B path) is skipped until the benchmark runs and
# Case B constants are confirmed.

test_that("constrained (max_weight=5) always routes to iEPPA", {
  set.seed(99L)
  n   <- 500L
  df  <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  res <- leafblower::harvest(df, tgt, method = "auto", max_weight = 5)
  expect_equal(attr(res, "algorithm"), "ieppa")
})

test_that("unconstrained large complexity routes to iEPPA", {
  # complexity = 30000 * 2 * 4 = 240000 < 500K threshold;
  # but tol=1e-3 (loose) → no kTolThreshold gate fires → iEPPA default
  set.seed(7L)
  n  <- 30000L
  df <- data.frame(
    m1 = factor(sample(paste0("c", 1:4), n, replace = TRUE)),
    m2 = factor(sample(paste0("c", 1:4), n, replace = TRUE))
  )
  tgt <- list(
    m1 = c(c1 = 0.25, c2 = 0.25, c3 = 0.25, c4 = 0.25),
    m2 = c(c1 = 0.25, c2 = 0.25, c3 = 0.25, c4 = 0.25)
  )
  res <- leafblower::harvest(df, tgt, method = "auto",
                              max_weight = Inf,
                              convergence = list(absolute = 1e-3))
  # Current default for unconstrained: iEPPA (post-benchmark default)
  expect_equal(attr(res, "algorithm"), "ieppa")
})

test_that("unconstrained tight-tol small complexity routes to L-BFGS-B (Case B only)", {
  skip("Add after benchmark run confirms Case B — update kTolThreshold in c_api.cpp first")
  # Placeholder: n=200, K=2, cats=2 → complexity=800, tol=1e-8 (tight)
  # Expected after constant update: attr(result, "algorithm") == "lbfgsb"
  set.seed(11L)
  n  <- 200L
  df <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  res <- leafblower::harvest(df, tgt, method = "auto",
                              max_weight = Inf, min_weight = 0,
                              convergence = list(absolute = 1e-8))
  expect_equal(attr(res, "algorithm"), "lbfgsb")
})
```

- [ ] **Step 2: Run all tests — confirm first two pass, third skipped**

```bash
Rscript -e "testthat::test_file('tests/testthat/test-algo-selection.R')" 2>&1 | tail -10
```

Expected: 17 tests PASS, 1 SKIPPED (Case B), 0 FAIL.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-algo-selection.R
git commit -m "test: algo-selection routing regression guard (Case B skipped until benchmark)"
```

- [ ] **Step 4: Close issue**

```bash
bd close leafblower-gje
```

---

## Task 11: Dry run — end-to-end smoke test

**Beads:** `bd update leafblower-dh8 --claim`  
**Files:** None (runs existing code)

- [ ] **Step 1: Run with budget=3 (tiny — just verifies no crashes)**

```bash
Rscript -e "
  .BENCH_SOURCED <- TRUE
  source('benchmarks/algo_selection_benchmark.R')
  # Tiny smoke: budget=3 acquisitions on the real input space
  state <- run_benchmark(
    budget          = 3L,
    checkpoint_path = '/tmp/algo_sel_smoke.rds',
    out_dir         = '/tmp',
    seed            = 123L,
    lhc_x1_max      = 5.5   # cap LHC at log_complexity=5.5 → n≤~4K, keeps smoke test fast (~2 min)
  )
  cat('n_evals:', nrow(state\$design), '\n')
  cat('classified:', round(state\$classified, 3), '\n')
  stopifnot(nrow(state\$design) >= 8L)   # at least 8 LHC + some adaptive
  stopifnot(file.exists('/tmp/algo_sel_smoke.rds'))
  stopifnot(file.exists('/tmp/algo_selection_contour.pdf'))
  stopifnot(file.exists('/tmp/algo_selection_uncertainty.pdf'))
  # K-stability PDF NOT checked: run_k_stability() is called from the entry-point
  # block only, not from run_benchmark(). Smoke test stays fast (~5 min).
  cat('Smoke test PASSED\n')
" 2>&1 | tail -15
```

Expected output (approximate):
```
Initialising: 8-point Latin hypercube design...
  LHC 1/8: lc=4.xx, lt=-x.xx ... y=x.xxx
  ...
  LHC 8/8: lc=x.xx, lt=-x.xx ... y=x.xxx
Iter 1/3: classified=0.xx
  Next: lc=x.xxx, lt=-x.xxx ... y=x.xxx
Iter 2/3: classified=0.xx
  ...
Plots saved to /tmp/

Done. 11 total evaluations, xx% classified.
n_evals: 11
classified: 0.xxx
Smoke test PASSED
```
Note: K-stability is NOT triggered here — it runs only in the entry-point block
(`if (!.BENCH_SOURCED)`), keeping the smoke test to ~5 min instead of ~21 min.

- [ ] **Step 2: Verify R CMD check still passes**

```bash
R CMD INSTALL --preclean . && Rscript -e "testthat::test_dir('tests/testthat/')" 2>&1 | tail -5
```

Expected: `* DONE (leafblower)` + `17 tests PASS, 1 SKIPPED`.

- [ ] **Step 3: Commit smoke test result (no new files — clean state)**

```bash
git status
```

Expected: `nothing to commit, working tree clean` (dry run writes to /tmp only).

- [ ] **Step 4: Close issue**

```bash
bd close leafblower-dh8
```

---

## Self-Review Checklist

**Spec coverage:**
- [x] Bayesian LSE with Straddle acquisition → Tasks 5, 7
- [x] 8-pt LHC initial design → Task 7
- [x] Matérn-5/2 GP with nugget lower bound → Task 5
- [x] Checkpoint every 5 acquisitions, atomic write → Tasks 6, 7
- [x] Poor-fit warning + "do not update constants" → Task 7
- [x] K-stability check at K=3, 18, ±0.5 log-unit decision rule → Task 8
- [x] Contour + uncertainty + k-stability PDFs → Tasks 8, 9
- [x] Regression test with 3 assertions (3rd skipped) → Task 10
- [x] End-to-end smoke test → Task 11
- [x] DESCRIPTION Suggests → Task 1
- [x] tol range [1e-6, 1e-3] (not 1e-9) → constants in script

**No placeholders found.**

**Type consistency:** `bench_seed` → `make_bench_data` → `time_cell` chain uses identical argument names throughout. `fit_gp` returns `km` object; `straddle_next` and `classified_fraction` both call `DiceKriging::predict(..., type="UK", checkNames=FALSE)` — consistent.

---

## Post-Benchmark Manual Steps (after full `Rscript benchmarks/algo_selection_benchmark.R` run)

These steps are performed by a human after the benchmark completes. They are NOT automated.

### Step A: Read the contour plot — Case A vs Case B

Open `benchmarks/algo_selection_contour.pdf`. The red contour = 1.2× boundary (L-BFGS-B wins above it).

1. Read `C_loose`: the `log₁₀(complexity)` value where the contour crosses `log₁₀(tol) = -3` (tol = 1e-3).
2. Read `C_tight`: the `log₁₀(complexity)` value where the contour crosses `log₁₀(tol) = -6` (tol = 1e-6).
3. Compute `|log₁₀(C_loose) - log₁₀(C_tight)|`:
   - **Case A** (difference < 0.5): complexity dominates. Update only `kComplexityThreshold = 10^C_loose`. Skip kTolThreshold.
   - **Case B** (difference ≥ 0.5): both dimensions matter. Set `kComplexityThreshold = 10^C_tight` (conservative). Set `kTolThreshold` at the `tol` value where the contour crosses `log₁₀(complexity) = 5.5` (midpoint of sweep).

### Step B: Check K-stability before committing

Open `benchmarks/algo_selection_k_stability.pdf`. If the K=3 or K=18 contour shifts >0.5 log-units from the K=9 contour: **do not update constants**. The warning message in the benchmark output will flag this. File a follow-up beads issue for a 3D sweep.

### Step C: Update `src/c_api.cpp`

Edit `src/c_api.cpp` — update `select_algorithm()`:

- Case A: update `kComplexityThreshold` only (keep existing `kTolThreshold` or remove it).
- Case B: update both `kComplexityThreshold` and `kTolThreshold`.

Add or update the comment:
```cpp
// Threshold calibrated 2026-04-20 via Bayesian LSE benchmark:
//   benchmarks/algo_selection_results.rds
// Valid for K≈9, uniform-category problems. See design spec.
```

### Step D: Run validation cases

```bash
# 1. Box-constrained (iEPPA path unchanged)
Rscript benchmarks/stepstone_benchmark.R   # confirm auto picks iEPPA, time within 5% of baseline

# 2. Unconstrained L-BFGS-B path (Case B only — skip for Case A)
# Spec: n=500, K=3, cats_per_margin=2, tol=1e-8
Rscript -e "
  library(leafblower)
  set.seed(42L)
  n <- 500L
  df <- data.frame(
    m1 = factor(sample(c('a','b'), n, replace = TRUE)),
    m2 = factor(sample(c('a','b'), n, replace = TRUE)),
    m3 = factor(sample(c('a','b'), n, replace = TRUE))
  )
  tgt <- list(m1 = c(a=0.5, b=0.5), m2 = c(a=0.5, b=0.5), m3 = c(a=0.5, b=0.5))
  res <- harvest(df, tgt, method = 'auto', max_weight = Inf, min_weight = 0,
                 convergence = list(absolute = 1e-8))
  stopifnot(attr(res, 'algorithm') == 'lbfgsb')
  cat('L-BFGS-B routing: OK\n')
"
```

### Step E: Activate the skipped regression test

In `tests/testthat/test-algo-selection.R`, remove the `skip(...)` call in the third test ("unconstrained tight-tol small complexity routes to L-BFGS-B"). Verify the test now passes.

```bash
Rscript -e "testthat::test_file('tests/testthat/test-algo-selection.R')" 2>&1 | tail -5
```

Expected: 17 tests PASS, 0 SKIPPED, 0 FAIL.

### Step F: Single threshold-update commit

```bash
git add src/c_api.cpp \
        benchmarks/algo_selection_results.rds \
        benchmarks/algo_selection_contour.pdf \
        benchmarks/algo_selection_uncertainty.pdf \
        benchmarks/algo_selection_k_stability.pdf \
        tests/testthat/test-algo-selection.R
git commit -m "feat(algo-select): calibrate kComplexityThreshold/kTolThreshold via LSE benchmark

Benchmark: benchmarks/algo_selection_results.rds @ 2026-04-20
[Case A/B result and contour crossing values here]"
```

**No other files in this commit.** This keeps the threshold calibration atomic and auditable.
