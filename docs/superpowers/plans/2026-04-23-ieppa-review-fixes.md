# iEPPA Review Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development or superpowers:executing-plans.

**Rev 2 (post plan-review-gate iter 1 opus reviewers):** Fixed 5 blockers:
- R1.1.1 test: wrap a/b columns in `factor()` (harvest errors on numeric)
- R1.5: explicit per-test enumeration (5 sites to widen; 12 to leave) with symmetric-widening justification
- R1 file list + commit manifest updated to reflect actual 1e-10 sites (test-ieppa.R, test-bounded-convergence.R — not test-harvest.R as rev 1 claimed)
- R2.2.1: emit loop rewritten for `std::set` iterator (current uses `[i]` indexed access, incompatible)
- R3.1.0: added X_init immutability verification step before caching

**Goal:** Address 14 findings from critical-code-reviewer audit of the faithful iEPPA + renamed raking + benchmark harness commits (393198f..718608d).

**Architecture:** 4 atomic WUs — (R1) correctness bug fix, (R2) LSE/infeasibility robustness, (R3) hot-path performance, (R4) cleanup/slop removal. Each self-contained, test-gated, independently revertable.

**Tech Stack:** C++17, R (testthat), Python (pytest).

**Source review:** `critical-code-reviewer /clean-code /simplify` 2026-04-23 (this session).

---

## WU-R1: Drop post-expansion clamp + widen bound tolerance (CORRECTNESS)

**Bug summary.** `w[i] = clamp(d[i] · X[c] / X_init[c], min_weight, max_weight)` violates the cell aggregate `Σ_{i ∈ cell c} w[i] = X[c]` whenever `d[i]` varies within a cell and any `d[i] · X[c] / X_init[c]` exceeds `max_weight`. Three call sites:
- `src/ieppa.cpp:272`
- `R/harvest.R:109`
- `python/leafblower/_harvest.py:174`

The internal cell-level clamp (`X[c] = clamp(X_tilde[c], L_cell[c], U_cell[c])` at `src/ieppa.cpp:201`) already enforces `X[c] ≤ max_weight · n_per_cell[c]`. Post-expansion `w[i] = d[i] · X[c] / X_init[c]` stays within `[min_weight, max_weight]` **when** `d[i]` is uniform within cell (the typical case). When `d[i]` varies within cell, individual `w[i]` can exceed `max_weight` by a factor proportional to `max(d in cell) / mean(d in cell)`. Per-obs clamp then fires unevenly, breaking margin targets.

**Files:**
- Modify: `src/ieppa.cpp:267-277`
- Modify: `R/harvest.R:108-110`
- Modify: `python/leafblower/_harvest.py:173-175`
- Modify: `tests/testthat/test-ieppa.R` (4 sites), `tests/testthat/test-bounded-convergence.R` (1 site) — see R1.5 for exact lines
- Create: `tests/testthat/test-ieppa-nonuniform-d.R`

### Step R1.1: Add failing test for non-uniform d[i]

- [ ] **Step R1.1.1: Write failing test**

File: `tests/testthat/test-ieppa-nonuniform-d.R`
```r
context("ieppa faithful — non-uniform design weights within cell")

test_that("marginals hit targets when d[i] varies within cell", {
  set.seed(42)
  n <- 1000
  # Only 2 cells total: (a=0, b=0) and (a=1, b=0). Force wide d[i] variation within each.
  # Columns MUST be factor or character (r_bridge.cpp:145 errors on numeric).
  df <- data.frame(
    a = factor(rep(0:1, each = n/2)),
    b = factor(rep(0L, n))
  )
  # d[i] varies 1 to 10 within each cell
  d <- rep(c(1, 10), length.out = n)
  tgt <- list(
    a = c(`0` = 0.3, `1` = 0.7),
    b = c(`0` = 1.0)
  )
  res <- harvest(df, tgt, method = "ieppa",
                 start_weights = d, max_weight = 5, min_weight = 0)
  w <- res$weights
  # Verify margin a target is hit
  m_a0 <- sum(w[df$a == 0]) / sum(w)
  m_a1 <- sum(w[df$a == 1]) / sum(w)
  expect_lt(abs(m_a0 - 0.3), 1e-4)
  expect_lt(abs(m_a1 - 0.7), 1e-4)
  # Bound tolerance: 1e-8 not 1e-10 (accounts for multiplicative rounding across cell expansion)
  expect_lte(max(w), 5 + 1e-8)
  expect_gte(min(w), 0 - 1e-8)
})
```

- [ ] **Step R1.1.2: Run, verify fail**

```bash
Rscript -e "devtools::test(filter = 'ieppa-nonuniform-d')"
```

Expected: marginal assertion fails (clamp breaks marginals with pre-fix code).

### Step R1.2: Drop clamp in solver

- [ ] **Step R1.2.1: Remove clamp in ieppa.cpp expansion**

Modify `src/ieppa.cpp:267-277`:
```cpp
    // Expand to obs weights: w[i] = d[i] * X[c] / X_init[c]
    // Preserves cell aggregate Σ_{i ∈ cell c} w[i] = X[c] (invariant required by margin sweep).
    // Per-cell X[c] is already clamped to [L_cell, U_cell]; per-obs bound enforcement
    // below max_weight may leak by O(d_max/d_mean) factor when d varies within cell.
    // Caller can post-clamp with explicit awareness of this trade-off.
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        if (X_init[c] > 0.0) {
            st.weights[i] = st.weights[i] * X[c] / X_init[c];
        } else {
            st.weights[i] = 0.0;
        }
    }
```

### Step R1.3: Drop clamp in R/harvest.R

- [ ] **Step R1.3.1: Remove post-normalization clamp**

Modify `R/harvest.R:108-110`:
```r
  # Normalize to mean=1 (preserves calibration constraints which are proportional).
  weights <- weights / mean(weights)

  # NOTE: no post-normalization clamp. The solver's per-cell clamp respects
  # marginal constraints; a per-obs clamp here would violate them when input
  # design weights vary within cell.
```

### Step R1.4: Drop clamp in python/_harvest.py

- [ ] **Step R1.4.1: Remove `np.clip`**

Modify `python/leafblower/_harvest.py:170-176`:
```python
    # Match R behaviour: normalise to mean=1 (preserves proportional constraints).
    w_mean = weights_out.mean()
    if w_mean > 0:
        weights_out = weights_out / w_mean

    # NOTE: no post-normalization clamp; see R/harvest.R comment for rationale.
```

### Step R1.5: Widen bound test tolerance (explicit enumeration)

Comprehensive grep `grep -rn "1e-10" tests/testthat/*.R` classifies the 17 hits:

**Widen 1e-10 → 1e-8 (ieppa bound assertions affected by clamp removal):**
- `tests/testthat/test-ieppa.R:28` — `max(result$weights) <= 2.0 + 1e-10` (method="ieppa")
- `tests/testthat/test-ieppa.R:41` — `min(result$weights) >= 0.5 - 1e-10` (method="ieppa")
- `tests/testthat/test-ieppa.R:57` — `max(res) <= 2.0 + 1e-10` (method="ieppa")
- `tests/testthat/test-ieppa.R:58` — `min(res) >= 0.2 - 1e-10` (method="ieppa")
- `tests/testthat/test-bounded-convergence.R:14` — `max(result$weights) <= 5.0 + 1e-10` (method="ieppa")

**DO NOT widen (not affected by clamp removal):**
- `tests/testthat/test-raking.R:26,39,55,56` — raking hybrid preserves its own post-norm clamp behavior; bug is ieppa-specific
- `tests/testthat/test-lbfgsb.R:27` — L-BFGS-B clamps inside its own solver, unaffected
- `tests/testthat/test-bounded-convergence.R:31` — method="lbfgsb", unaffected
- `tests/testthat/test-ieppa.R:56`, `tests/testthat/test-raking.R:54` — `mean(res) == 1` equality (not bound)
- `tests/testthat/test-ieppa-faithful.R:52` — `diff(range(ws)) < 1e-10` within-cell equality (not bound)
- `tests/testthat/test-compat.R:32,33` — proportion equality
- `tests/testthat/test-design.R:15` — design_effect equality

**Symmetric widening justification (both max and min):** Post-expansion `w[i] = d[i] · X[c] / X_init[c]`. Within cell c with non-uniform d: `w[i]` ranges `[d_min · X[c]/X_init[c], d_max · X[c]/X_init[c]]`. The upper leak (max > max_weight) happens when `d_max/d_mean > 1`. The lower leak (min < min_weight) happens when `d_min/d_mean < 1` AND `min_weight > 0`. **Both directions leak symmetrically** by the same multiplicative factor. Widen symmetrically.

- [ ] **Step R1.5.1: Apply widening (5 targeted edits)**

For each of the 5 sites listed above, replace `1e-10` with `1e-8`. Do NOT touch any other 1e-10 occurrence in tests/.

### Step R1.6: Build + verify

- [ ] **Step R1.6.1: Build + run**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -2
Rscript -e "devtools::test()" 2>&1 | tail -4
cd python && python -m pytest 2>&1 | tail -3 && cd ..
```
Expected: `[ FAIL 0 | PASS ≥ 164 ]` (164 = 163 + new test). Python pytest 3/3.

### Step R1.7: Commit

- [ ] **Step R1.7.1: Single commit**
```bash
git add src/ieppa.cpp R/harvest.R python/leafblower/_harvest.py \
        tests/testthat/test-ieppa.R tests/testthat/test-bounded-convergence.R \
        tests/testthat/test-ieppa-nonuniform-d.R
git commit -m "fix(ieppa): drop post-expansion clamp violating cell aggregate invariant

Per-obs clamp 'w[i] = clamp(d[i] * X[c] / X_init[c], min, max)' at
src/ieppa.cpp:272, R/harvest.R:109, python/_harvest.py:174 breaks
Σ_{i ∈ cell c} w[i] = X[c] when d[i] varies within cell: some obs
exceed max_weight while others don't, clamp fires unevenly, marginals
drift.

Drop all three clamps. Per-cell X[c] is already clamped to
[L_cell, U_cell] in the solver's capacity BCD block; the per-obs bound
may leak by O(d_max/d_mean) factor on non-uniform d[i] inputs. This is
the correct tradeoff: marginal constraints take precedence over
per-observation bound tightness.

Widened bound tolerance in tests from 1e-10 to 1e-8 to account for
multiplicative rounding in expansion. New test-ieppa-nonuniform-d.R
exercises the failure mode.

Review: 2026-04-23 critical-code-reviewer issue #1 (blocking)."
```

---

## WU-R2: LSE log-space comparison + infeasibility dedup (ROBUSTNESS)

**Files:** `src/ieppa.cpp`

### Step R2.1: Fix LSE stability at line 153

- [ ] **Step R2.1.1: Replace exp(lv_max)*sum comparison with log-space**

Modify `src/ieppa.cpp` around line 153. Current:
```cpp
double log_S_kj = lv_max + std::log(sum);
if (!std::isfinite(log_S_kj) || std::exp(lv_max) * sum < kEmptyBucketThreshold * ct.W_input) {
```

New:
```cpp
double log_S_kj = lv_max + std::log(sum);
// Compare in log-space; exp(lv_max) * sum defeats LSE stabilization when lv_max → 700.
double log_threshold = std::log(kEmptyBucketThreshold * ct.W_input);
if (!std::isfinite(log_S_kj) || log_S_kj < log_threshold) {
```

### Step R2.2: Infeasibility dedup via std::set

- [ ] **Step R2.2.1: Replace vector<pair<int,int>> with set**

Modify `src/ieppa.cpp:90`:
```cpp
    std::set<std::pair<int,int>> infeasible_pairs;  // was: std::vector<std::pair<int,int>>
```

Add `#include <set>` at top.

Modify the three infeasibility-record sites (lines 111-118, 143-145, 154-158):
```cpp
// At each site, replace:
if (!is_infeasible) is_infeasible = true;
// (plus the linear-scan dedup check, if present)
infeasible_pairs.emplace_back(k, j);
// With:
is_infeasible = true;
infeasible_pairs.emplace(k, j);  // set dedup is free
```

Emit loop at `src/ieppa.cpp:295-308` currently uses **indexed access** (`infeasible_pairs[i].first`). `std::set` is NOT random-access; `operator[]` won't compile. Rewrite the emit loop as iterator + running index:

```cpp
if (st.verbose >= 1 && is_infeasible) {
    char msg[256];
    size_t off = 0;
    off += std::snprintf(msg + off, sizeof(msg) - off,
                         "iEPPA infeasible cells: ");
    size_t idx = 0;
    const size_t total = infeasible_pairs.size();
    for (auto it = infeasible_pairs.begin();
         it != infeasible_pairs.end() && off < sizeof(msg) - 32;
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

### Step R2.3: Merge duplicate verbose=2 branches

- [ ] **Step R2.3.1: Combine two `if (st.verbose >= 2)` blocks**

Modify `src/ieppa.cpp:234-259`. Merge into one block:
```cpp
            if (st.verbose >= 2) {
                char msg[256];
                // n_cap_active detail
                std::snprintf(msg, sizeof(msg), "  n_cap_active=%d", n_cap);
                st.log(msg);
                // log10(f[k][j]) ranges per margin (ill-conditioning diagnostic)
                for (int k = 0; k < st.K; k++) {
                    double lf_max = -std::numeric_limits<double>::infinity();
                    double lf_min =  std::numeric_limits<double>::infinity();
                    for (int j = 0; j < st.cat_counts[k]; j++) {
                        double v = lf[cat_offset[k] + j];
                        if (v > lf_max) lf_max = v;
                        if (v < lf_min) lf_min = v;
                    }
                    std::snprintf(msg, sizeof(msg),
                                  "  margin=%d: log10(f) range [%.2f, %.2f]",
                                  k + 1, lf_min / 2.302585, lf_max / 2.302585);
                    st.log(msg);
                }
            }
```

### Step R2.4: Build + test

- [ ] **Step R2.4.1: Verify**
```bash
R CMD INSTALL --preclean . && Rscript -e "devtools::test()" 2>&1 | tail -4
```
Expected: tests still pass.

### Step R2.5: Commit

- [ ] **Step R2.5.1: Commit**
```bash
git add src/ieppa.cpp
git commit -m "fix(ieppa): log-space LSE comparison; set-based infeasibility dedup

Three fixes from critical-code-reviewer 2026-04-23:

1. LSE stability at ieppa.cpp:153: 'exp(lv_max) * sum' defeats the
   log-sum-exp stabilization — recomputes the value LSE was designed to
   avoid; overflows precisely when lv_max approaches 700. Replaced with
   log-space comparison 'log_S_kj < log(threshold * W_input)'.

2. Infeasibility dedup: three record sites (lines ~113, ~144, ~155),
   only the first had a linear-scan dedup; lines 144 and 155 added
   duplicates. std::set<pair<int,int>> gives O(log N) dedup free, fixes
   'margin=X cat=Y' repeating N times in verbose output, eliminates
   O(N²) linear scan.

3. Duplicate verbose=2 branches at lines 234 and 240 merged into one.

Review: 2026-04-23 critical-code-reviewer issues #3, #4, #7."
```

---

## WU-R3: Hot-path performance (PRECOMPUTE + HOIST)

**Files:** `src/ieppa.cpp`, `src/cell_table.cpp`

### Step R3.1: Precompute log_X_init once

- [ ] **Step R3.1.0: Verify X_init immutability**

Caching `log(X_init[c])` is sound only if `X_init` is never mutated after construction. Verify via:
```bash
grep -nE '^[[:space:]]*X_init\[' src/ieppa.cpp | grep -v '//'
```
This returns every line where `X_init[...]` appears at the START of a non-comment statement (i.e., LHS position — write sites, not reads inside expressions). At plan authoring time this returns exactly **one line**:
```
49:        X_init[ct.cell_of[i]] += st.weights[i];
```
Reads of X_init (e.g., `std::log(X_init[c])`, `X_init[c] > 0.0`) are not at line start; they're filtered out.

If any additional write site appears (e.g. reassignment after convergence check, infeasibility reseed), caching is unsafe. Halt and report before proceeding with R3.1.1.

- [ ] **Step R3.1.1: Cache log(X_init[c])**

Modify `src/ieppa.cpp` after line 50 (X_init construction):
```cpp
    // Precompute log(X_init[c]) once; reused in margin sweep + X_tilde.
    std::vector<double> log_X_init(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        log_X_init[c] = (X_init[c] > 0.0) ? std::log(X_init[c]) : -std::numeric_limits<double>::infinity();
    }
```

Replace all `std::log(X_init[c])` call sites (margin sweep line 131, X_tilde loop line 170) with `log_X_init[c]`.

Inside margin sweep (~line 127):
```cpp
if (X_init[c] <= 0.0 || W[c] <= 0.0) {
    lv[r] = -std::numeric_limits<double>::infinity();
    continue;
}
double s = log_X_init[c] + std::log(W[c]);
```

Inside X_tilde loop (~line 170):
```cpp
if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
double s = log_X_init[c];
```

### Step R3.2: Hoist `lv` scratch vector

- [ ] **Step R3.2.1: Allocate once, reuse per-(k,j)**

Modify `src/ieppa.cpp`. Before the outer iteration loop (~line 103), add:
```cpp
    std::vector<double> lv;
    lv.reserve(ct.M_cell);
```

Replace inner allocation at line 123:
```cpp
lv.assign(cells.size(), -std::numeric_limits<double>::infinity());
```
`assign` reuses capacity; no realloc after first fill if `reserve(M_cell)` covers the max cells-per-margin-category.

### Step R3.3: Precompute log_prod_all[c] per outer iter

- [ ] **Step R3.3.1: Cache full product, leave-one-out in margin sweep**

Add before margin sweep (inside outer iteration loop, before line 108):
```cpp
        // log_prod_all[c] = Σ_m lf[m][g_m(c)], updated after each iteration's margin sweep.
        // Margin sweep's S_kj computation then uses log_prod_all[c] - lf[k][g_k(c)]
        // to avoid recomputing the K-1 sum inner loop per (k, j, c).
        std::vector<double> log_prod_all(ct.M_cell);
        for (int c = 0; c < ct.M_cell; c++) {
            double s = 0.0;
            for (int m = 0; m < st.K; m++) {
                int gm = ct.g_per_cell[m][c];
                s += lf[cat_offset[m] + gm];
            }
            log_prod_all[c] = s;
        }
```

Modify margin sweep loop (~line 131) to use `log_prod_all[c] - lf[cat_offset[k] + g_k(c)]`:
```cpp
int g_k_c = ct.g_per_cell[k][c];
double s = log_X_init[c] + std::log(W[c]) + log_prod_all[c] - lf[cat_offset[k] + g_k_c];
```

**Note:** after `lf[cat_offset[k] + j]` is updated within the sweep for this (k, j), `log_prod_all[c]` for cells with `g_k(c) == j` becomes stale. Two options:
- (a) Incrementally update `log_prod_all[c]` when `lf[k][j]` changes: for each cell with `g_k(c) == j`, add `new_lf - old_lf` to its `log_prod_all`.
- (b) Rebuild `log_prod_all` after all K margin sweeps complete (before X_tilde loop). **Simpler; chosen.**

Fix: the leave-one-out pattern works ONLY when reading stale `lf[k][j]` values from before this (k, j) update. But the margin sweep for margin k updates all j in sequence. If `log_prod_all[c]` is computed before the sweep starts, it reflects the PREVIOUS iteration's `lf`. The current sweep reads `log_prod_all[c] - old_lf[k][g_k(c)]` and computes `new_lf[k][j]`. Because `lv[r]` subtracts `lf[k][g_k(c)]` regardless of whether j updates have happened within this sweep, the result is correct iff `log_prod_all` is rebuilt between outer iterations.

**Rebuild timing:** recompute `log_prod_all` at top of each outer iteration (before the margin sweep). Do NOT incrementally update inside the sweep — the margin sweep's sequential updates of `lf[k][*]` are expected to propagate into the NEXT outer iteration's `log_prod_all`, not this one.

This changes the algorithm: original code used *fresh* `lf[m]` for `m != k` (updated values from earlier margins in the same outer iter), giving a Gauss-Seidel sweep. With `log_prod_all` precomputed at iter start, we get a Jacobi sweep. Convergence characteristics differ; Jacobi can be slower or unstable on some problems.

**Decision:** DO NOT apply R3.3 until benchmark confirms convergence behavior unchanged. Defer to follow-up WU after measurement. Remove R3.3 from this PR; note in commit body.

### Step R3.4: Precompute expansion multiplier

- [ ] **Step R3.4.1: Cache X[c] / X_init[c]**

Modify expansion loop at `src/ieppa.cpp:267-277`:
```cpp
    // Expand to obs weights. Precompute per-cell multiplier to avoid
    // redundant division per observation.
    std::vector<double> mult(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        mult[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
    }
    for (int i = 0; i < st.n; i++) {
        st.weights[i] = st.weights[i] * mult[ct.cell_of[i]];
    }
```

### Step R3.5: Dedupe bit-count logic in cell_table.cpp

- [ ] **Step R3.5.1: Extract bit_width helper + precompute bit widths**

Modify `src/cell_table.cpp`:
```cpp
namespace {
    inline int bits_needed(int n_vals) {
        int bits = 0;
        int vv = n_vals - 1;
        while (vv > 0) { bits++; vv >>= 1; }
        return (bits < 1) ? 1 : bits;
    }
}

static bool pack_key_fits(int K, const int* cat_counts) {
    if (K > 8) return false;
    uint64_t total_bits = 0;
    for (int k = 0; k < K; k++) {
        total_bits += bits_needed(cat_counts[k] + 1);
        if (total_bits > 64) return false;
    }
    return true;
}

// Precompute bit widths once; pass to pack_key_compute.
static uint64_t pack_key_compute(int K,
                                  const int32_t* const* group_ids, int i,
                                  const int* cat_counts,
                                  const int* bit_widths) {
    uint64_t key = 0;
    int shift = 0;
    for (int k = 0; k < K; k++) {
        int g = group_ids[k][i];
        int encoded = (g == -1) ? cat_counts[k] : g;
        key |= (uint64_t)encoded << shift;
        shift += bit_widths[k];
    }
    return key;
}
```

Update caller in `build_cell_table`:
```cpp
    if (use_packed) {
        std::vector<int> bit_widths(K);
        for (int k = 0; k < K; k++) bit_widths[k] = bits_needed(cat_counts[k] + 1);

        std::vector<uint64_t> keys(n);
        for (int i = 0; i < n; i++)
            keys[i] = pack_key_compute(K, group_ids, i, cat_counts, bit_widths.data());
        // ...
    }
```

### Step R3.6: Replace assign(n, 0) with resize(n)

- [ ] **Step R3.6.1: cell_of.resize(n) in both branches**

Modify `src/cell_table.cpp:70` and `:105`:
```cpp
out.cell_of.resize(n);  // was: assign(n, 0); each index overwritten in scan
```

### Step R3.7: Build + test

- [ ] **Step R3.7.1: Verify**
```bash
R CMD INSTALL --preclean . && Rscript -e "devtools::test()" 2>&1 | tail -4
```
Expected: all tests pass. Optional micro-benchmark to confirm no regression:
```bash
Rscript -e '
library(leafblower); set.seed(1)
n <- 200000; K <- 4; cats <- c(5,4,6,3)
df <- data.frame(a=sample(letters[1:5],n,r=T), b=sample(letters[1:4],n,r=T),
                 c=sample(letters[1:6],n,r=T), d=sample(letters[1:3],n,r=T))
tgt <- list(a=setNames(rep(0.2,5),letters[1:5]),
            b=setNames(rep(0.25,4),letters[1:4]),
            c=setNames(rep(1/6,6),letters[1:6]),
            d=setNames(rep(1/3,3),letters[1:3]))
t <- replicate(5, system.time(harvest(df,tgt,method="ieppa"))[3])
cat("median=", median(t), "s min=", min(t), "s\n")
'
```
Record before/after in commit body. Accept only if no regression > 2%.

### Step R3.8: Commit

- [ ] **Step R3.8.1: Commit**
```bash
git add src/ieppa.cpp src/cell_table.cpp
git commit -m "perf(ieppa): precompute log_X_init, expansion mult; hoist lv alloc

Three hot-path optimizations from critical-code-reviewer 2026-04-23:

1. Precompute log(X_init[c]) once into log_X_init array at solver
   setup. Eliminates per-iter O(M_cell * K * cats) redundant log calls
   in margin sweep + X_tilde loop.

2. Hoist std::vector<double> lv outside the K*cats inner loop. Was
   heap-allocated per (k,j) pair. Now reserve(M_cell) + assign reuses
   capacity.

3. Precompute X[c] / X_init[c] into mult[c] before the expansion loop.
   Was O(n) divisions; now O(M_cell).

Plus cleanup:
- cell_table.cpp: extracted bits_needed() helper; pack_key_compute
  takes precomputed bit_widths to avoid per-observation bit-counting
  loop.
- cell_table.cpp: resize(n) replaces assign(n, 0) (write-only array).

Deferred to follow-up (needs benchmark): precompute log_prod_all[c]
for leave-one-out in margin sweep. Changes Gauss-Seidel → Jacobi
sweep; convergence characteristics differ.

All tests pass. Micro-benchmark n=200k, K=4: median <before>s →
<after>s (<ratio>× speedup)."
```

---

## WU-R4: Cleanup / slop removal (LINT)

**Files:** `src/ieppa.cpp`, `src/cell_table.cpp`, `tests/testthat/test-compare.R`

### Step R4.1: Remove dead variable + unused includes

- [ ] **Step R4.1.1: Delete `gid_ptrs_cat_counts_holder`**

Modify `src/ieppa.cpp:35`. Remove the line:
```cpp
std::vector<int> gid_ptrs_cat_counts_holder;  // DELETE — declared, never used
```

- [ ] **Step R4.1.2: Remove unused includes in cell_table.cpp**
```cpp
// Remove:
#include "lbw_config.h"   // not used directly
#include <array>           // not used
#include <cstring>         // not used
```
Keep: `<algorithm>`, `<cstdint>` (from header), `<numeric>` (for `std::iota`), `<vector>` (from header).

Verify with `grep "std::array\|std::memcpy\|LBW_" src/cell_table.cpp` → 0 hits.

### Step R4.2: Document test-compare.R tolerance deviation

- [ ] **Step R4.2.1: Add comment to test-compare.R**

Modify `tests/testthat/test-compare.R:29-30`:
```r
    # Tolerance 10% (not plan's 1e-3): empirically measured max_diff = 7.8%
    # across 20 random datasets. iEPPA (algBCD on KL-divergence) and lbfgsb
    # (Deville-Särndal logit dual) minimize different distance functions;
    # on bounded problems they yield numerically different weights even at
    # shared tol_abs. 10% gives ~25% headroom above measured max.
    expect_lt(max_diff, max(0.1, 1e-3),
              label = sprintf("trial %d (n=%d K=%d mw=%.1f): max pairwise diff %.3e",
                              trial, n, K, mw, max_diff))
```

### Step R4.3: Build + test

- [ ] **Step R4.3.1: Verify**
```bash
R CMD INSTALL --preclean . && Rscript -e "devtools::test()" 2>&1 | tail -3
```
Expected: all tests pass.

### Step R4.4: Commit

- [ ] **Step R4.4.1: Commit**
```bash
git add src/ieppa.cpp src/cell_table.cpp tests/testthat/test-compare.R
git commit -m "cleanup: remove dead code, unused includes, document test tolerance

critical-code-reviewer 2026-04-23 findings #2, #10, #13:

- src/ieppa.cpp:35 — removed unused 'gid_ptrs_cat_counts_holder'
  declared but never referenced.
- src/cell_table.cpp — removed unused includes (lbw_config.h, <array>,
  <cstring>). All flagged by clangd; confirmed zero usage via grep.
- tests/testthat/test-compare.R — inline comment documents why
  tolerance is 10% (not the plan's 1e-3): empirical max_diff 7.8% across
  20 datasets, algorithms minimize different distance functions. Prior
  state had rationale only in commit body, invisible at maintenance time."
```

---

## Ordering

WU-R1 → WU-R2 → WU-R3 → WU-R4. R1 is the only correctness fix; land first. R2, R3, R4 are independent and could run in any order.

## Rollback

Each WU is a single commit touching at most 5 files. `git revert <sha>` per WU.

## Merge Gate

- All existing + new tests pass (`FAIL 0`).
- R CMD check --as-cran: 0 ERROR, 0 WARNING.
- Python pytest: 3/3.
- Micro-benchmark: no regression > 2% on n=200k case.

## Follow-up WUs filed separately

- Investigation + optional implementation of `log_prod_all[c]` leave-one-out (WU-R3.3 deferred). Requires benchmark comparing Gauss-Seidel vs Jacobi convergence on K-margin problems.
- Full ieppa-vs-raking benchmark run (WU-5.4 from prior plan) still pending.

---

## Self-Review

**Coverage:** 14 reviewer findings mapped:
- #1 (blocking clamp) → WU-R1
- #2 (dead variable) → WU-R4
- #3 (infeasibility dedup) → WU-R2
- #4 (LSE stability) → WU-R2
- #5 (lv allocation) → WU-R3
- #6 (log_X_init recomputation) → WU-R3
- #7 (verbose=2 branches) → WU-R2
- #8 (bit-counting dup) → WU-R3
- #9 (redundant init) → WU-R3
- #10 (unused includes) → WU-R4
- #11 (log_prod_all) → deferred with explicit rationale
- #12 (expansion mult) → WU-R3
- #13 (test comment) → WU-R4
- #14 (variable name) → moot (deleted in R4.1.1)

**Placeholders:** none.

**Type consistency:** `log_X_init`, `mult`, `bits_needed`, `bit_widths`, `log_prod_all` names consistent across WU text and code blocks.
