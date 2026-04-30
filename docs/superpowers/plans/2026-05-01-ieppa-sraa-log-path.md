# iEPPA SRAA-m Log-Path Extension — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans.

**Goal:** Extend SRAA-m Anderson Acceleration to ieppa's log path (high-compression surveys where n/M_cell >= 2), removing the `&& use_linear` restriction.

**Architecture:** `f_eval_lf` dispatches based on the current `use_linear` flag; `apply_single_margin_log` is hoisted to homotopy scope; `unpack_lf` extended (via separate post-step in `f_eval_lf`) to rebuild `X_tilde`/`log_W` for log-path consistency, and `cell_lf_hwm` reset is folded into `unpack_lf` itself.

**Tech Stack:** C++17, R package build (`R CMD INSTALL --preclean .`), testthat 3.

**Mechanism:** Branch `f_eval_lf` on the `use_linear` capture; reuse existing `apply_single_margin_log`, `bulk_log`, `bulk_exp_clipped`, and `cells_by_margin_cat` infrastructure. No new algorithmic component.

**Forbidden:**
- Do NOT duplicate the SRAA while-loop into a separate log-path block. Single dispatched `f_eval_lf` is the contract.
- Do NOT introduce a parallel `apply_single_margin_log_for_sraa` lambda. The hoisted lambda is the single source of truth for the log sweep.
- Do NOT hard-code `use_linear` capture by value at lambda creation. SRAA fallback may flip it; capture by reference (`[&]` already in use).
- Do NOT remove or weaken the `cell_lf_hwm` invariant when the log path becomes SRAA-driven.

**Audit:**
- Spy: `aa_count > 0` (via `IEPPAResult::aa_accepted_count`) on a forced-log-path problem proves SRAA actually ran the log iterate.
- Spy: `LBW_IEPPA_FORCE_PATH=log` env var with `accelerate=TRUE` exercises the new branch deterministically.
- Spy: residual `max_error` and weight-sum invariants (`sum(w_calibrated) == n`) prove correctness of the dispatched sweep.

---

## Task LL1: Hoist `apply_single_margin_log` to homotopy scope

**Files modified:**
- `src/ieppa.cpp` (move lambda from ~line 774 to ~line 600, immediately after `apply_single_margin_linear`)

**Captured-variable audit (verified by grep):**
- `lv` — declared at line 283 (homotopy scope, BEFORE the `for(iter_in_lvl)` loop at line 762). Already accessible.
- `log_empty_threshold` — declared at line 419 (homotopy scope). Already accessible.
- `sor_active`, `sor_auto_v`, `omega_fixed_v`, `sor_omega` — all homotopy-scope.
- `cells_by_margin_cat`, `cat_offset`, `st.cat_counts`, `st.K`, `st.targets`, `ct.g_per_cell` — all outer-scope.
- `X_init`, `W`, `log_W`, `log_X_init`, `lf`, `cell_lf` — all homotopy-scope.
- `alpha` — declared at line 319, homotopy-scope. **Note:** mutated by `compute_alpha()` inside the for-loop at line 765. Hoisting changes nothing for the SRAA call site (SRAA uses `alpha` from the previous step's update; same as linear path already does).
- `record_empty`, `record_nonempty` — closure helpers at homotopy scope.
- **No reference** to `iter_in_lvl`, `overflow_trip`, or any for-loop-local variable. **Safe to hoist.**

**Before** (line ~774, inside `if (!sraa_active_lvl) { for (...) {`):

```cpp
        // Do NOT read st.max_weight here — homotopy levels pass
        // current_max_weight indirectly via the shared U_cell already built.

        auto apply_single_margin_log = [&](int k) -> bool {
            // ... 60-line body ...
            return false;  // log path does not trip overflow mid-sweep
        };
```

**After** (move verbatim to ~line 600, immediately after the closing brace of `apply_single_margin_linear` and BEFORE the comment block introducing `f_eval_lf`):

```cpp
        // ────────────────────────────────────────────────────────────────────
        // SRAA-m fixed-point map for the log path of this homotopy level.
        // Hoisted out of the inner for-loop so the SRAA while-loop can invoke
        // it via f_eval_lf when use_linear=false. Captures lf, cell_lf, log_W,
        // X_init, log_X_init, lv, cells_by_margin_cat, alpha, sor_*, by [&].
        // Returns false unconditionally — log path does not trip overflow
        // mid-sweep (LSE stabilization absorbs all magnitude).
        // ────────────────────────────────────────────────────────────────────
        auto apply_single_margin_log = [&](int k) -> bool {
            // ... identical body ...
            return false;
        };
```

Then at the original site (~line 774), delete the lambda definition; the call sites inside the for-loop continue to bind to the hoisted symbol via name lookup.

**Compile gate:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -20
```
Must produce `* DONE (leafblower)` with zero new warnings.

**Test gate:**
```bash
Rscript -e 'testthat::test_dir("tests/testthat", filter = "ieppa")' 2>&1 | tail -2
```
Existing ieppa tests must still pass — this task is a pure refactor.

**Commit:**
```
refactor(ieppa): hoist apply_single_margin_log to homotopy scope

Move the log-path margin lambda from inside the non-SRAA for-loop to the
homotopy level scope, immediately after apply_single_margin_linear. This
prepares LL2 to dispatch f_eval_lf on use_linear; no behavior change.
Verified no for-loop-local captures (lv, log_empty_threshold, alpha all
already at homotopy scope).
```

---

## Task LL2: Extend `f_eval_lf` to dispatch on `use_linear`

**Files modified:**
- `src/ieppa.cpp` (`unpack_lf` helper at lines 87-115, `f_eval_lf` lambda at lines 614-646)

### LL2.a — Extend `unpack_lf` to reset `cell_lf_hwm` when needed

The hwm is read by `apply_single_margin_linear` (line ~583) but the log path does not currently maintain it. After SRAA extrapolation overwrites `lf`, the hwm can be stale-low (extrapolated `cell_lf` values may exceed the previously cached max). Linear path is unaffected because the hwm is monotone-nondecreasing under in-place margin updates AND the linear path does not depend on hwm precision (only used for overflow forecasting). **Decision:** recompute hwm inside `unpack_lf` so both paths are SRAA-safe.

**Before** (lines 87-115):

```cpp
static inline void unpack_lf(const std::vector<double>& src,
                             std::vector<double>& lf,
                             std::vector<double>& f_lin,
                             std::vector<double>& cell_lf,
                             std::vector<double>& X_cur,
                             const lbw::CellTable& ct,
                             const std::vector<double>& X_init,
                             int K,
                             const std::vector<int>& cat_offset) {
    const int total_cats = cat_offset[K];
    for (int i = 0; i < total_cats; i++) {
        lf[i]    = src[i];
        f_lin[i] = std::exp(src[i]);
    }
    const int M = ct.M_cell;
    std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
    for (int k = 0; k < K; k++) {
        const int* gk  = ct.g_per_cell[k].data();
        const int  off = cat_offset[k];
        for (int c = 0; c < M; c++) {
            int g = gk[c];
            if (g < 0) continue;
            cell_lf[c] += src[off + g];
        }
    }
    for (int c = 0; c < M; c++) {
        X_cur[c] = (X_init[c] > 0.0) ? X_init[c] * std::exp(cell_lf[c]) : 0.0;
    }
}
```

**After:** add `log_X_init` and `cell_lf_hwm` outparams; recompute hwm in the same pass that fills `cell_lf`.

```cpp
// unpack_lf: rebuild derived state for either path. Sets lf, f_lin, cell_lf,
// X_cur, and refreshes cell_lf_hwm = max_c(log_X_init[c] + cell_lf[c]).
// hwm refresh is mandatory because SRAA extrapolation can produce new lf
// values whose cell_lf sums exceed the previously cached monotone high-water
// mark; without refresh, the linear-path overflow guard reads a stale-low
// bound and fails to trip when it should.
static inline void unpack_lf(const std::vector<double>& src,
                             std::vector<double>& lf,
                             std::vector<double>& f_lin,
                             std::vector<double>& cell_lf,
                             std::vector<double>& X_cur,
                             const lbw::CellTable& ct,
                             const std::vector<double>& X_init,
                             const std::vector<double>& log_X_init,
                             int K,
                             const std::vector<int>& cat_offset,
                             double& cell_lf_hwm) {
    const int total_cats = cat_offset[K];
    for (int i = 0; i < total_cats; i++) {
        lf[i]    = src[i];
        f_lin[i] = std::exp(src[i]);
    }
    const int M = ct.M_cell;
    std::fill(cell_lf.begin(), cell_lf.end(), 0.0);
    for (int k = 0; k < K; k++) {
        const int* gk  = ct.g_per_cell[k].data();
        const int  off = cat_offset[k];
        for (int c = 0; c < M; c++) {
            int g = gk[c];
            if (g < 0) continue;
            cell_lf[c] += src[off + g];
        }
    }
    double hwm = std::numeric_limits<double>::lowest();
    for (int c = 0; c < M; c++) {
        if (X_init[c] > 0.0) {
            X_cur[c] = X_init[c] * std::exp(cell_lf[c]);
            double v = log_X_init[c] + cell_lf[c];
            if (std::isfinite(v) && v > hwm) hwm = v;
        } else {
            X_cur[c] = 0.0;
        }
    }
    cell_lf_hwm = hwm;
}
```

**Update call sites** — three of them in `ieppa.cpp`:
1. Inside `f_eval_lf` (line 615): add `log_X_init` and `cell_lf_hwm` args.
2. Inside the SRAA outer-stall revert path (line 706): add same args.
3. Any other callers — grep `unpack_lf(` to confirm only two sites.

```bash
grep -n "unpack_lf(" src/ieppa.cpp
```

### LL2.b — Dispatch `f_eval_lf` on `use_linear`

**Before** (lines 614-646):

```cpp
auto f_eval_lf = [&](std::vector<double>& flat) -> double {
    unpack_lf(flat, lf, f_lin, cell_lf, X_cur, ct, X_init,
              st.K, cat_offset);
    bool overflow = false;
    for (int k = 0; k < st.K && !overflow; k++) {
        if (apply_single_margin_linear(k)) overflow = true;
    }
    pack_lf(lf, flat);
    if (overflow) return std::numeric_limits<double>::infinity();

    double W_total = 0.0;
    for (int c = 0; c < ct.M_cell; c++) W_total += X_cur[c];
    if (!(W_total > 0.0)) return std::numeric_limits<double>::infinity();
    double errRp = 0.0;
    for (int k = 0; k < st.K; k++) {
        const int nj = st.cat_counts[k];
        std::fill(S_lin.begin(), S_lin.begin() + nj, 0.0);
        const int* gk = ct.g_per_cell[k].data();
        for (int c = 0; c < ct.M_cell; c++) {
            int j = gk[c];
            if (j >= 0 && j < nj) S_lin[j] += X_cur[c];
        }
        for (int j = 0; j < nj; j++) {
            double e = std::fabs(S_lin[j] / W_total - st.targets[k][j]);
            if (e > errRp) errRp = e;
        }
    }
    return errRp;
};
```

**After:** branch on `use_linear`. Log path needs `X_tilde`/`log_W` rebuilt before the sweep (because `apply_single_margin_log` reads `log_W[c]` from the previous P1.1 update — but during SRAA, that update was skipped, so we must seed `log_W = log(W)` or `log_W = 0` if no capacity correction has run yet). Allocate `X_tilde` lazily on the first log-path SRAA evaluation.

```cpp
auto f_eval_lf = [&](std::vector<double>& flat) -> double {
    unpack_lf(flat, lf, f_lin, cell_lf, X_cur, ct, X_init, log_X_init,
              st.K, cat_offset, cell_lf_hwm);

    bool overflow = false;
    if (use_linear) {
        for (int k = 0; k < st.K && !overflow; k++) {
            if (apply_single_margin_linear(k)) overflow = true;
        }
    } else {
        // Log-path setup: X_tilde = X_init * exp(cell_lf) (== X_cur from
        // unpack_lf); log_W must reflect the last capacity correction. If
        // no correction has run this level (W still 1), log_W is 0.
        // bulk_log is idempotent and cheap (O(M_cell)).
        if (X_tilde.size() != static_cast<size_t>(ct.M_cell)) {
            X_tilde.assign(ct.M_cell, 0.0);
        }
        std::copy(X_cur.begin(), X_cur.end(), X_tilde.begin());
        lbw::bulk_log(W.data(), log_W.data(), ct.M_cell);

        for (int k = 0; k < st.K && !overflow; k++) {
            if (apply_single_margin_log(k)) overflow = true;  // never true
        }
        // After the log sweep, refresh X_tilde from updated cell_lf so the
        // errRp computation below sees the post-sweep cell masses.
        for (int c = 0; c < ct.M_cell; c++) {
            X_tilde[c] = (X_init[c] > 0.0)
                ? X_init[c] * std::exp(cell_lf[c]) : 0.0;
        }
    }
    pack_lf(lf, flat);
    if (overflow) return std::numeric_limits<double>::infinity();

    // errRp: read from X_cur (linear) or X_tilde (log).
    const std::vector<double>& X_eval = use_linear ? X_cur : X_tilde;
    double W_total = 0.0;
    for (int c = 0; c < ct.M_cell; c++) W_total += X_eval[c];
    if (!(W_total > 0.0)) return std::numeric_limits<double>::infinity();
    double errRp = 0.0;
    for (int k = 0; k < st.K; k++) {
        const int nj = st.cat_counts[k];
        std::fill(S_lin.begin(), S_lin.begin() + nj, 0.0);
        const int* gk = ct.g_per_cell[k].data();
        for (int c = 0; c < ct.M_cell; c++) {
            int j = gk[c];
            if (j >= 0 && j < nj) S_lin[j] += X_eval[c];
        }
        for (int j = 0; j < nj; j++) {
            double e = std::fabs(S_lin[j] / W_total - st.targets[k][j]);
            if (e > errRp) errRp = e;
        }
    }
    return errRp;
};
```

**Compile gate:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -20
```

**Test gate:**
```bash
Rscript -e 'testthat::test_dir("tests/testthat", filter = "ieppa")' 2>&1 | tail -2
```

**Commit:**
```
feat(ieppa): dispatch f_eval_lf on use_linear, hwm-refresh in unpack_lf

unpack_lf now refreshes cell_lf_hwm in its existing O(M_cell) pass.
f_eval_lf branches on use_linear: linear path unchanged, log path runs
apply_single_margin_log against rebuilt X_tilde and log_W. errRp reads
from the path-correct cell-mass vector.
```

---

## Task LL3: Activate SRAA on the log path; fix X←X_tilde sync

**Files modified:**
- `src/ieppa.cpp` line 653 (`sraa_active_lvl` declaration)
- `src/ieppa.cpp` line 758 (post-SRAA X sync)

**Before** (line 653):

```cpp
const bool sraa_active_lvl = st.accelerate && use_linear;
```

**After:**

```cpp
// SRAA-m is path-agnostic: f_eval_lf dispatches on use_linear internally,
// so accelerate=TRUE activates SRAA on both linear and log paths.
const bool sraa_active_lvl = st.accelerate;
```

**Before** (line 755-758):

```cpp
            // Sync X_cur → X so the post-loop expansion at line ~1624 uses the
            // final linear-path cell masses. The non-SRAA path updates X[c] inside
            // the P1.1 fused block; SRAA skips that block.
            std::copy(X_cur.begin(), X_cur.end(), X.begin());
```

**After:**

```cpp
            // Sync working cell masses → X for the post-loop expansion at
            // ~line 1624. Linear path's X_cur and log path's X_tilde both
            // hold X_init * exp(cell_lf) post-sweep; pick the one consistent
            // with the path that just ran.
            if (use_linear) {
                std::copy(X_cur.begin(), X_cur.end(), X.begin());
            } else {
                std::copy(X_tilde.begin(), X_tilde.end(), X.begin());
            }
```

**Path-flip safety (case F from the spec):** If `use_linear` is flipped from `true` to `false` mid-solve by a linear-overflow fallback, the existing fallback site already calls `ieppa_sraa.clear()`. After the fallback, `f_eval_lf` re-evaluates with `use_linear=false`, which triggers the log-path setup branch (lazy `X_tilde` allocation, `log_W` refresh from `W`). No additional code needed — the dispatch is dynamic by construction (capture-by-reference via `[&]`).

**Compile gate:**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -20
```

**Test gate:**
```bash
Rscript -e 'testthat::test_dir("tests/testthat", filter = "ieppa")' 2>&1 | tail -2
```

**Manual smoke:**
```bash
LBW_IEPPA_FORCE_PATH=log Rscript -e '
suppressMessages(library(leafblower))
set.seed(1); n <- 5000
df <- data.frame(a = sample(1:4, n, TRUE), b = sample(1:5, n, TRUE))
res <- rake(df, list(a = c(0.25,0.25,0.25,0.25),
                     b = c(0.2,0.2,0.2,0.2,0.2)),
            method = "ieppa", accelerate = TRUE, verbose = 1L)
cat("aa_count =", res$aa_accepted_count,
    " max_error =", res$max_error, "\n")
'
```
`aa_count > 0` proves SRAA actually ran on the log path.

**Commit:**
```
feat(ieppa): enable SRAA-m on log path, fix X sync for both paths

Drop the use_linear gate from sraa_active_lvl. Post-SRAA sync copies the
path-correct cell-mass vector (X_cur for linear, X_tilde for log) into X
so the post-loop expansion sees consistent state. Path flip during the
solve is handled by capture-by-reference; the existing fallback clear()
already resets SRAA history.
```

---

## Task LL4: Tests + benchmark

**Files modified/created:**
- `tests/testthat/test-ieppa-sraa-log-path.R` (new)
- `bench/ieppa-sraa-log-path.R` (new, optional bench run)

### LL4.a — Regression test

```r
# tests/testthat/test-ieppa-sraa-log-path.R
test_that("SRAA-m activates on the log path (forced)", {
  withr::local_envvar(LBW_IEPPA_FORCE_PATH = "log")
  set.seed(2026)
  n <- 50000  # high compression: n/M_cell will exceed 2
  df <- data.frame(
    a = sample(1:4, n, TRUE),
    b = sample(1:5, n, TRUE),
    c = sample(1:3, n, TRUE)
  )
  targets <- list(
    a = rep(0.25, 4),
    b = rep(0.2, 5),
    c = rep(1/3, 3)
  )
  res <- leafblower::rake(
    df, targets, method = "ieppa",
    accelerate = TRUE, verbose = 0L
  )
  expect_true(is.finite(res$max_error))
  expect_lt(res$max_error, 1e-6)
  expect_equal(sum(res$weights), n, tolerance = 1e-6)
  expect_gt(res$aa_accepted_count, 0L,
            label = "SRAA must accept at least one log-path step")
})

test_that("SRAA-m log path matches non-accelerated log path", {
  withr::local_envvar(LBW_IEPPA_FORCE_PATH = "log")
  set.seed(2026)
  n <- 20000
  df <- data.frame(a = sample(1:5, n, TRUE), b = sample(1:6, n, TRUE))
  tg <- list(a = rep(0.2, 5), b = rep(1/6, 6))
  r_off <- leafblower::rake(df, tg, method = "ieppa", accelerate = FALSE)
  r_on  <- leafblower::rake(df, tg, method = "ieppa", accelerate = TRUE)
  expect_equal(r_on$weights, r_off$weights, tolerance = 1e-4)
})
```

### LL4.b — Benchmark scaffold (optional, run manually)

```r
# bench/ieppa-sraa-log-path.R
set.seed(2026)
n <- 200000
df <- data.frame(
  a = sample(1:8, n, TRUE),
  b = sample(1:10, n, TRUE),
  c = sample(1:6, n, TRUE),
  d = sample(1:5, n, TRUE)
)
tg <- list(a = rep(1/8,8), b = rep(0.1,10), c = rep(1/6,6), d = rep(0.2,5))
Sys.setenv(LBW_IEPPA_FORCE_PATH = "log")
res <- bench::mark(
  off = leafblower::rake(df, tg, method = "ieppa", accelerate = FALSE),
  on  = leafblower::rake(df, tg, method = "ieppa", accelerate = TRUE),
  iterations = 5, check = FALSE
)
print(res[, c("expression", "median", "n_itr")])
```
Interleaved before/after via the single `bench::mark` call.

**Test gate:**
```bash
Rscript -e 'testthat::test_dir("tests/testthat", filter = "ieppa-sraa-log-path")' 2>&1 | tail -10
Rscript -e 'testthat::test_dir("tests/testthat")' 2>&1 | tail -2
```
Full suite must remain green; the new file must show 2 passing tests with `aa_accepted_count > 0`.

**Stepstone smoke:**
```bash
Rscript -e 'source("bench/stepstone_small.R")' 2>&1 | tail -5
```
Existing stepstone benchmark must still complete without overflow or convergence regression.

**Commit:**
```
test(ieppa): regression guards for SRAA-m log-path acceleration

Two tests under LBW_IEPPA_FORCE_PATH=log: (1) accelerate=TRUE on a
high-compression problem produces aa_accepted_count > 0 and converges,
(2) accelerated and non-accelerated log paths agree to 1e-4 on weights.
Adds bench/ieppa-sraa-log-path.R for interleaved before/after timing.
```

---

## Self-review

**Spec coverage A–G:**
- **A (hoist apply_single_margin_log):** LL1.
- **B (rebuild X_tilde/log_W in f_eval_lf):** LL2.b — done explicitly in `f_eval_lf`, not folded into `unpack_lf` (per spec design decision).
- **C (dispatch f_eval_lf on use_linear):** LL2.b.
- **D (drop `&& use_linear` from sraa_active_lvl):** LL3.
- **E (X←X_tilde sync at SRAA exit):** LL3.
- **F (path-flip mid-solve):** LL3 — handled by capture-by-reference; verified no extra code needed.
- **G (cell_lf_hwm reset in unpack_lf):** LL2.a.

**Placeholders:** none. Every code block is concrete C++ or R.

**Variable consistency:** `use_linear`, `cell_lf_hwm`, `X_tilde`, `log_W`, `bulk_log`, `apply_single_margin_log`, `f_eval_lf` match existing names in `src/ieppa.cpp`.

**Risk flags:**
- LL2.a changes `unpack_lf`'s signature — call-site grep is mandatory before edit. Failure mode is a clean compile error, not silent breakage.
- LL2.b's `bulk_log(W, log_W, M_cell)` per SRAA step is O(M_cell). For large problems with many SRAA steps this adds up; if profiling shows it dominates, hoist to once-per-`f_eval_lf` only when `W` has changed (a `W_dirty` flag set by the P1.1 block). Out of scope for this plan.
- The `X_tilde` post-sweep refresh (second loop in LL2.b log branch) is redundant with the values that `apply_single_margin_log` implicitly maintains via `cell_lf`, but materializing them keeps the errRp loop path-uniform. Cost: O(M_cell exp). Acceptable.

**Confidence:** 88. Architecture is mechanical given the existing linear-path scaffolding. Lowered from 95 because (a) the `log_W` semantics during SRAA without an intervening capacity correction may need empirical verification — the assumption "log_W = log(W); W=1 initially → log_W=0" is correct by construction but the first SRAA step on the log path is the first place this is exercised end-to-end; (b) `cell_lf_hwm` recompute cost is amortized but not measured.
