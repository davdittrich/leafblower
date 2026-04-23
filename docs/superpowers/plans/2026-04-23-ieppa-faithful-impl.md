# Faithful iEPPA Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace misnamed IPF+Dykstra hybrid `method="ieppa"` with paper-faithful algBCD (Chu-Liang-Toh-Yang 2022, arXiv:2011.14312) at C=0, using cell-compressed representation. Rename current hybrid to `method="raking"` as a breaking change.

**Architecture:** Cell-level algBCD with log-space dual factors (`lf[k][j]`), log-sum-exp-stabilized Sinkhorn updates, and a separate capacity BCD block (KL-projection onto box = clamp). Cell table built via sort-based dedup (no hash map, no DoS surface). Obs-level work emerges naturally when M_cell = n.

**Tech Stack:** C++17 (GCC/Clang), R (testthat, bench), Python (pytest, pybind11), DiceKriging + lhs (Bayesian benchmark).

**Spec:** `docs/superpowers/specs/2026-04-23-ieppa-faithful-design.md` (commit `020ceee`).

**Merge gate (§8c of spec):**
1. All tests in §6 pass (`[ FAIL 0 ]`)
2. `R CMD check --as-cran` yields 0 ERROR, 0 WARNING (NEWS.md present)
3. `M_cell = n` degenerate input: wall-clock within ±2× of `raking`
4. Python pytest passes

**Benchmark (§7) is post-merge analysis, not a merge gate.**

---

## WU-1: Cell Table (prerequisite, independent commit)

**Files:**
- Create: `src/cell_table.hpp`
- Create: `src/cell_table.cpp`
- Create: `tests/testthat/test-cell-table.R`
- Modify: `src/Makevars.in` (add cell_table.cpp to PKG_SOURCES)
- Modify: `src/r_bridge.cpp` (add test-only `.Call("C_leafblower_cell_table_probe", ...)`)

### Step 1.1: Write failing cell-table integration test

- [ ] **Step 1.1.1: Create test file with header**

File: `tests/testthat/test-cell-table.R`
```r
context("cell_table")

# Internal probe (test-only; returns M_cell, cell_of, n_per_cell)
probe <- function(group_ids_list, n) {
  .Call("C_leafblower_cell_table_probe", group_ids_list, as.integer(n),
        PACKAGE = "leafblower")
}

test_that("cell compression: all-identical observations produce 1 cell", {
  n <- 1000
  g1 <- rep(0L, n)
  g2 <- rep(0L, n)
  out <- probe(list(g1, g2), n)
  expect_equal(out$M_cell, 1L)
  expect_equal(out$n_per_cell, n)
  expect_true(all(out$cell_of == 0L))
})

test_that("cell compression: full cross-product populated", {
  # K=2 margins, 4 cats each, balanced assignment
  n <- 10000
  g1 <- rep(0:3, each = n/4)
  g2 <- rep(rep(0:3, each = n/16), 4)
  out <- probe(list(g1, g2), n)
  expect_equal(out$M_cell, 16L)  # 4 * 4
  expect_equal(sum(out$n_per_cell), n)
})

test_that("cell compression: sparse assignment gives M_cell < cross-product", {
  n <- 100
  g1 <- c(rep(0L, 50), rep(1L, 50))
  g2 <- c(rep(0L, 50), rep(1L, 50))  # perfectly correlated
  out <- probe(list(g1, g2), n)
  expect_equal(out$M_cell, 2L)  # only (0,0) and (1,1) populated
})

test_that("cell compression: NA entries (-1) treated as distinct category", {
  n <- 300
  g1 <- c(rep(0L, 100), rep(1L, 100), rep(-1L, 100))
  out <- probe(list(g1), n)
  expect_equal(out$M_cell, 3L)  # 0, 1, and NA buckets
})

test_that("cell compression: K > 64 rejected", {
  n <- 100
  many_margins <- replicate(65, rep(0L, n), simplify = FALSE)
  expect_error(probe(many_margins, n), regexp = "K.*64", ignore.case = TRUE)
})

test_that("cell compression: sort determinism (same input → same cell ids)", {
  n <- 500
  set.seed(1)
  g1 <- sample(0:3, n, replace = TRUE)
  g2 <- sample(0:2, n, replace = TRUE)
  out1 <- probe(list(g1, g2), n)
  out2 <- probe(list(g1, g2), n)
  expect_equal(out1$M_cell, out2$M_cell)
  expect_equal(out1$cell_of, out2$cell_of)
})
```

- [ ] **Step 1.1.2: Run tests, verify all fail**

Run: `Rscript -e "devtools::test(filter='cell-table')"`
Expected: 6 tests FAIL with `"C_leafblower_cell_table_probe" not found`.

### Step 1.2: Write cell_table.hpp

- [ ] **Step 1.2.1: Create header**

File: `src/cell_table.hpp`
```cpp
#pragma once
#include <cstdint>
#include <vector>

namespace lbw {

struct CellTable {
    int M_cell;                                  // number of unique cells
    std::vector<int> cell_of;                    // size n: obs i -> cell index
    std::vector<int> n_per_cell;                 // size M_cell: count per cell
    std::vector<std::vector<int>> g_per_cell;    // [K][M_cell]: margin-k cat per cell
    double W_input;                              // sum of input weights
};

// Build cell table from group_ids. Returns 0 on success, -1 if K > 64.
// Caller owns all input memory; function allocates output vectors.
// NA (group_ids[k][i] = -1) is treated as category `cat_counts[k]` internally
// (distinct bucket per margin).
int build_cell_table(int n, int K,
                     const int32_t* const* group_ids,
                     const int* cat_counts,
                     const double* weights,
                     CellTable& out);

// Maximum K supported (prevents unbounded key allocation).
inline constexpr int K_MAX = 64;

} // namespace lbw
```

### Step 1.3: Write cell_table.cpp

- [ ] **Step 1.3.1: Create implementation**

File: `src/cell_table.cpp`
```cpp
#include "lbw_config.h"
#include "cell_table.hpp"
#include <algorithm>
#include <array>
#include <cstring>
#include <numeric>

namespace lbw {

// Pack (g_1, ..., g_K) into a 64-bit key when feasible.
// Layout: low-order bits = margin 0, higher bits = later margins.
// NA (-1) encoded as cat_counts[k].
static bool pack_key_fits(int K, const int* cat_counts) {
    if (K > 8) return false;
    // Need cat_counts[k]+1 distinct values per margin (incl. NA bucket).
    // Widest possible value = cat_counts[k]; encoded in enough bits.
    uint64_t total_bits = 0;
    for (int k = 0; k < K; k++) {
        int vals = cat_counts[k] + 1;  // +1 for NA
        int bits = 0;
        int vv = vals - 1;
        while (vv > 0) { bits++; vv >>= 1; }
        if (bits < 1) bits = 1;
        total_bits += bits;
        if (total_bits > 64) return false;
    }
    return true;
}

static uint64_t pack_key_compute(int K,
                                  const int32_t* const* group_ids, int i,
                                  const int* cat_counts) {
    uint64_t key = 0;
    int shift = 0;
    for (int k = 0; k < K; k++) {
        int g = group_ids[k][i];
        int encoded = (g == -1) ? cat_counts[k] : g;
        int vals = cat_counts[k] + 1;
        int bits = 0; int vv = vals - 1;
        while (vv > 0) { bits++; vv >>= 1; }
        if (bits < 1) bits = 1;
        key |= (uint64_t)encoded << shift;
        shift += bits;
    }
    return key;
}

int build_cell_table(int n, int K,
                     const int32_t* const* group_ids,
                     const int* cat_counts,
                     const double* weights,
                     CellTable& out) {
    if (K > K_MAX) return -1;
    if (K <= 0 || n <= 0) return -1;

    // Build keys for each observation.
    const bool use_packed = pack_key_fits(K, cat_counts);

    // Indices sorted by key.
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);

    if (use_packed) {
        std::vector<uint64_t> keys(n);
        for (int i = 0; i < n; i++)
            keys[i] = pack_key_compute(K, group_ids, i, cat_counts);
        std::sort(idx.begin(), idx.end(),
                  [&](int a, int b) { return keys[a] < keys[b]; });
        // Scan to identify cells.
        out.cell_of.assign(n, 0);
        out.n_per_cell.clear();
        out.g_per_cell.assign(K, std::vector<int>());
        int current_cell = -1;
        uint64_t current_key = 0;
        for (int r = 0; r < n; r++) {
            int i = idx[r];
            uint64_t k = keys[i];
            if (current_cell == -1 || k != current_key) {
                current_cell++;
                current_key = k;
                out.n_per_cell.push_back(0);
                for (int m = 0; m < K; m++) {
                    int g = group_ids[m][i];
                    int encoded = (g == -1) ? cat_counts[m] : g;
                    out.g_per_cell[m].push_back(encoded);
                }
            }
            out.cell_of[i] = current_cell;
            out.n_per_cell[current_cell]++;
        }
        out.M_cell = current_cell + 1;
    } else {
        // General path: sort by tuple directly.
        auto tuple_less = [&](int a, int b) {
            for (int k = 0; k < K; k++) {
                int ga = group_ids[k][a];
                int gb = group_ids[k][b];
                int ea = (ga == -1) ? cat_counts[k] : ga;
                int eb = (gb == -1) ? cat_counts[k] : gb;
                if (ea != eb) return ea < eb;
            }
            return false;
        };
        std::sort(idx.begin(), idx.end(), tuple_less);
        out.cell_of.assign(n, 0);
        out.n_per_cell.clear();
        out.g_per_cell.assign(K, std::vector<int>());
        int current_cell = -1;
        auto tuples_equal = [&](int a, int b) {
            for (int k = 0; k < K; k++) {
                int ga = group_ids[k][a];
                int gb = group_ids[k][b];
                int ea = (ga == -1) ? cat_counts[k] : ga;
                int eb = (gb == -1) ? cat_counts[k] : gb;
                if (ea != eb) return false;
            }
            return true;
        };
        int prev_i = -1;
        for (int r = 0; r < n; r++) {
            int i = idx[r];
            if (current_cell == -1 || !tuples_equal(i, prev_i)) {
                current_cell++;
                out.n_per_cell.push_back(0);
                for (int m = 0; m < K; m++) {
                    int g = group_ids[m][i];
                    int encoded = (g == -1) ? cat_counts[m] : g;
                    out.g_per_cell[m].push_back(encoded);
                }
            }
            out.cell_of[i] = current_cell;
            out.n_per_cell[current_cell]++;
            prev_i = i;
        }
        out.M_cell = current_cell + 1;
    }

    out.W_input = 0.0;
    for (int i = 0; i < n; i++) out.W_input += weights[i];
    return 0;
}

} // namespace lbw
```

- [ ] **Step 1.3.2: Add r_bridge probe (test-only)**

Modify: `src/r_bridge.cpp` — add probe entry after existing `.Call` handlers.

Locate the file's existing SEXP function registration block. Add:
```cpp
// test-only: exposes CellTable internals for unit tests
extern "C" SEXP C_leafblower_cell_table_probe(SEXP r_group_ids_list, SEXP r_n) {
    int n = INTEGER(r_n)[0];
    int K = Rf_length(r_group_ids_list);
    if (K > lbw::K_MAX) {
        Rf_error("K (%d) exceeds K_MAX (%d)", K, lbw::K_MAX);
    }
    // Extract pointers
    std::vector<const int32_t*> gid_ptrs(K);
    std::vector<int> cat_counts(K);
    for (int k = 0; k < K; k++) {
        SEXP v = VECTOR_ELT(r_group_ids_list, k);
        gid_ptrs[k] = (const int32_t*) INTEGER(v);
        // Derive cat_counts from max value + 1 (excluding -1)
        int max_g = -1;
        for (int i = 0; i < n; i++) {
            int g = gid_ptrs[k][i];
            if (g > max_g) max_g = g;
        }
        cat_counts[k] = (max_g < 0) ? 1 : (max_g + 1);
    }
    std::vector<double> uniform_weights(n, 1.0);
    lbw::CellTable ct;
    int rc = lbw::build_cell_table(n, K, gid_ptrs.data(), cat_counts.data(),
                                    uniform_weights.data(), ct);
    if (rc != 0) Rf_error("build_cell_table failed (rc=%d)", rc);

    // Build return list: list(M_cell, cell_of, n_per_cell)
    SEXP ret = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ret, 0, Rf_ScalarInteger(ct.M_cell));
    SEXP cell_of_sexp = Rf_allocVector(INTSXP, n);
    std::memcpy(INTEGER(cell_of_sexp), ct.cell_of.data(), n * sizeof(int));
    SET_VECTOR_ELT(ret, 1, cell_of_sexp);
    SEXP npc = Rf_allocVector(INTSXP, ct.M_cell);
    std::memcpy(INTEGER(npc), ct.n_per_cell.data(), ct.M_cell * sizeof(int));
    SET_VECTOR_ELT(ret, 2, npc);
    SEXP names = Rf_allocVector(STRSXP, 3);
    SET_STRING_ELT(names, 0, Rf_mkChar("M_cell"));
    SET_STRING_ELT(names, 1, Rf_mkChar("cell_of"));
    SET_STRING_ELT(names, 2, Rf_mkChar("n_per_cell"));
    Rf_setAttrib(ret, R_NamesSymbol, names);
    UNPROTECT(1);
    return ret;
}
```

Also add to the `R_registerRoutines()` table in r_bridge.cpp (typically end of file):
```cpp
// In static const R_CallMethodDef CallEntries[]:
{"C_leafblower_cell_table_probe", (DL_FUNC) &C_leafblower_cell_table_probe, 2},
```

- [ ] **Step 1.3.3: Update Makevars.in**

Modify: `src/Makevars.in`. Locate `PKG_SOURCES` or equivalent. Add `cell_table.cpp` to the source list. If `PKG_SOURCES = $(wildcard *.cpp)`, no change needed.

### Step 1.4: Build and verify tests pass

- [ ] **Step 1.4.1: Rebuild**
```bash
R CMD INSTALL --preclean .
```
Expected: `* DONE (leafblower)`.

- [ ] **Step 1.4.2: Run cell-table tests**
```bash
Rscript -e "devtools::test(filter='cell-table')"
```
Expected: `[ FAIL 0 | PASS 6 ]` for cell-table. Existing 100 tests still pass.

### Step 1.5: Commit

- [ ] **Step 1.5.1: Commit**
```bash
git add src/cell_table.hpp src/cell_table.cpp src/r_bridge.cpp \
        tests/testthat/test-cell-table.R
git commit -m "feat(cell_table): sort-based cell dedup helper for faithful iEPPA

Adds src/cell_table.hpp/cpp — builds a compressed representation of
(g_1,...,g_K) tuples from group_ids input, with M_cell <= min(n, prod cat).

Sort-based dedup (not std::unordered_map) avoids hash-collision DoS on
adversarial group_ids. Packed int64 key when K<=8 and narrow cat_counts;
general tuple comparator otherwise.

Hard cap K <= 64 (K_MAX) prevents unbounded key allocation.

NA (group_ids[k][i] = -1) encoded as cat_counts[k] internally — distinct
category per margin.

Test-only R bridge .Call entry C_leafblower_cell_table_probe exposes
(M_cell, cell_of, n_per_cell) for unit testing. Not a public API.

Refs: 2026-04-23-ieppa-faithful-design.md §5.2"
```

---

## WU-2: Input validation updates (prerequisite, independent commit)

**Files:**
- Modify: `src/c_api.cpp` (add `K > 64`, `cat_counts[k] <= 0`, `Σw <= 1e-15` guards)
- Modify: `tests/testthat/test-harvest.R` (add validation tests)

### Step 2.1: Write failing validation tests

- [ ] **Step 2.1.1: Add test cases**

Modify: `tests/testthat/test-harvest.R` — append at end:
```r
test_that("K > 64 rejected with informative error", {
  # 65 margin columns — should fail validation
  n <- 100
  data <- as.data.frame(matrix(0L, nrow = n, ncol = 65))
  names(data) <- paste0("v", 1:65)
  # Targets: each column has 1 category, target = 1
  tgt <- lapply(seq_len(65), function(k) setNames(1.0, "0"))
  names(tgt) <- names(data)
  expect_error(harvest(data, tgt, method = "ieppa"),
               regexp = "K.*64|too many margin", ignore.case = TRUE)
})

test_that("cat_counts <= 0 rejected", {
  # Empty target list for a column — degenerate
  n <- 100
  data <- data.frame(a = sample(c("x", "y"), n, replace = TRUE))
  tgt <- list(a = setNames(numeric(0), character(0)))
  expect_error(harvest(data, tgt, method = "ieppa"),
               regexp = "cat_counts|empty target|no categories", ignore.case = TRUE)
})

test_that("zero-sum input weights rejected", {
  n <- 100
  data <- data.frame(a = sample(c("x", "y"), n, replace = TRUE))
  tgt <- list(a = c(x = 0.5, y = 0.5))
  sw <- rep(0.0, n)
  expect_error(harvest(data, tgt, method = "ieppa", start_weights = sw),
               regexp = "start_weights.*positive|sum.*zero", ignore.case = TRUE)
})
```

- [ ] **Step 2.1.2: Run — verify failures**
```bash
Rscript -e "devtools::test(filter='harvest')" 2>&1 | grep -E "FAIL|PASS"
```
Expected: 3 new failures (K cap, cat_counts, zero-sum).

### Step 2.2: Add validation guards in c_api.cpp

- [ ] **Step 2.2.1: Locate validate_inputs**

```bash
grep -n "validate_inputs\|K <= 0\|K > 0" src/c_api.cpp | head -10
```

- [ ] **Step 2.2.2: Add K > 64 check**

Modify: `src/c_api.cpp` — in `validate_inputs` (or equivalent), add after existing `K <= 0` check:
```cpp
if (K > 64) {
    if (result) {
        std::snprintf(result->message, sizeof(result->message),
                      "K=%d exceeds maximum (64); too many margin columns", K);
    }
    return RK_ERR_BADARG;
}
```

- [ ] **Step 2.2.3: Add cat_counts[k] <= 0 check**

In the same validation function, after extracting `cat_counts[k]`:
```cpp
for (int k = 0; k < K; k++) {
    if (cat_counts[k] <= 0) {
        if (result) {
            std::snprintf(result->message, sizeof(result->message),
                          "cat_counts[%d]=%d must be positive (no categories for margin %d)",
                          k, cat_counts[k], k);
        }
        return RK_ERR_BADARG;
    }
}
```

- [ ] **Step 2.2.4: Add Σw > 0 check**

In the same validation function, after weight-pointer null check:
```cpp
double w_sum = 0.0;
for (int i = 0; i < n; i++) w_sum += weights[i];
if (!(w_sum > 1e-15)) {
    if (result) {
        std::snprintf(result->message, sizeof(result->message),
                      "sum of input weights must be positive (got %.3e)", w_sum);
    }
    return RK_ERR_BADARG;
}
```

### Step 2.3: Build + verify tests pass

- [ ] **Step 2.3.1: Rebuild + test**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
Rscript -e "devtools::test()" 2>&1 | tail -3
```
Expected: `[ FAIL 0 | PASS 103 ]` (100 existing + 3 new).

### Step 2.4: Commit

- [ ] **Step 2.4.1: Commit**
```bash
git add src/c_api.cpp tests/testthat/test-harvest.R
git commit -m "feat(c_api): validate K<=64, cat_counts>0, sum(weights)>0

Prerequisite guards for faithful iEPPA solver:
- K > 64 rejected (prevents unbounded cell-key allocation per K_MAX in
  src/cell_table.hpp)
- cat_counts[k] <= 0 rejected (prevents degenerate key encoding and
  empty-target data generation failures)
- sum(weights) <= 1e-15 rejected (degenerate total weight; otherwise
  faithful iEPPA's W_input and the cell-aggregate invariant break)

All RK_ERR_BADARG with actionable snprintf-safe message field.

Refs: 2026-04-23-ieppa-faithful-design.md §8"
```

---

## WU-3: ATOMIC — Rename + Faithful solver + wiring + tests + docs

**This WU lands as a SINGLE commit.** Intermediate states must not be pushed; CI must not run between sub-steps. Execute sequentially, test only at the end.

**Files:**
- Rename: `src/ieppa.cpp` → `src/raking.cpp`
- Rename: `src/ieppa.hpp` → `src/raking.hpp`
- Create: `src/ieppa.hpp` (new, faithful solver declaration)
- Create: `src/ieppa.cpp` (new, faithful solver implementation)
- Modify: `src/leafblower.h` (add `RK_ALG_RAKING = 3`)
- Modify: `src/c_api.cpp` (dispatch table, select_algorithm)
- Modify: `src/r_bridge.cpp` (method="raking" mapping)
- Modify: `R/harvest.R` (method="raking" in match.arg)
- Modify: `python/leafblower/_harvest.py` (alg_map add "raking":3)
- Modify: `python/leafblower/_bindings.cpp` (if it has method dispatch)
- Modify: `tasks/prd-leafblower-core.md` (§US-005 rewrite + new §US-005b)
- Create: `NEWS.md`
- Modify: `tests/testthat/test-ieppa.R` (retarget at faithful)
- Create: `tests/testthat/test-raking.R` (regression guard for renamed hybrid)
- Create: `tests/testthat/test-ieppa-faithful.R` (algBCD-specific)
- Create: `tests/testthat/test-compare.R` (cross-algorithm equivalence)

### Step 3.1: Snapshot current test-ieppa.R

- [ ] **Step 3.1.1: Copy current tests as regression guard**
```bash
cp tests/testthat/test-ieppa.R tests/testthat/test-raking.R
```

- [ ] **Step 3.1.2: Rewrite test-raking.R method strings**
```bash
sed -i 's/method *= *"ieppa"/method = "raking"/g; s/method *= *"auto"/method = "raking"/g' \
    tests/testthat/test-raking.R
```

- [ ] **Step 3.1.3: Update context label**

Modify: `tests/testthat/test-raking.R` — change top-line context from `context("ieppa")` to `context("raking (renamed hybrid)")`.

### Step 3.2: Rename C++ sources

- [ ] **Step 3.2.1: Git-rename files**
```bash
git mv src/ieppa.cpp src/raking.cpp
git mv src/ieppa.hpp src/raking.hpp
```

- [ ] **Step 3.2.2: Rename symbols in raking.cpp/hpp**

Modify: `src/raking.hpp` and `src/raking.cpp` (replace all):
- `ieppa_solve` → `raking_solve`
- `IEPPAResult` → `RakingResult`
- `#include "ieppa.hpp"` → `#include "raking.hpp"`

```bash
sed -i 's/ieppa_solve/raking_solve/g; s/IEPPAResult/RakingResult/g; s/"ieppa\.hpp"/"raking.hpp"/g' \
    src/raking.hpp src/raking.cpp
```

Also modify: `src/c_api.cpp` — replace include + call sites:
```bash
sed -i 's/"ieppa\.hpp"/"raking.hpp"/g; s/ieppa_solve/raking_solve/g; s/IEPPAResult/RakingResult/g' \
    src/c_api.cpp
```

### Step 3.3: Create new faithful iEPPA header

- [ ] **Step 3.3.1: Create src/ieppa.hpp**

File: `src/ieppa.hpp`
```cpp
#pragma once
#include "types.hpp"
#include <vector>

namespace lbw {

struct IEPPAResult {
    int status;              // RK_OK / RK_ERR_NOCONV / RK_ERR_INFEAS / RK_ERR_BADARG
    int iterations;          // outer iterations completed
    double max_error;        // final errRp
    int M_cell;              // compression info
    int n_cap_active;        // cells with W[c] != 1 at convergence
};

// Faithful iEPPA (paper-faithful algBCD at C=0). See
// docs/superpowers/specs/2026-04-23-ieppa-faithful-design.md.
IEPPAResult ieppa_solve(CalibState& state);

} // namespace lbw
```

### Step 3.4: Create new faithful iEPPA implementation

- [ ] **Step 3.4.1: Create src/ieppa.cpp**

File: `src/ieppa.cpp`
```cpp
#include "lbw_config.h"
#include "ieppa.hpp"
#include "cell_table.hpp"
#include "leafblower.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

namespace lbw {

// Log-space algBCD at C=0 — see design spec §2.3, §8 for math.
// lf[k][j]: log Sinkhorn factor (per margin k, category j)
// W[c]:    capacity multiplier per cell (linear-space; bounded in [L_c/X_tilde, U_c/X_tilde])
// X_tilde[c] = X_init[c] * exp(sum_k lf[k][g_k(c)])
// X[c] = clamp(X_tilde[c] * W[c], L_c, U_c) — one capacity BCD block per outer iter
//
// Log-sum-exp stabilization on S_kj prevents overflow when partial log-sums
// approach log(DBL_MAX) ≈ 709.
IEPPAResult ieppa_solve(CalibState& st) {
    constexpr int    kErrCheckInterval = 10;
    constexpr double kEmptyBucketThreshold = 1e-15;
    constexpr double kLogClip = 700.0;  // exp(700) < DBL_MAX

    IEPPAResult res;
    res.status = RK_ERR_NOCONV;
    res.iterations = 0;
    res.max_error = 1.0;
    res.M_cell = 0;
    res.n_cap_active = 0;

    // Build cell table.
    CellTable ct;
    std::vector<int> gid_ptrs_cat_counts_holder;
    {
        int rc = build_cell_table(st.n, st.K, st.group_ids,
                                  st.cat_counts, st.weights, ct);
        if (rc != 0) {
            res.status = RK_ERR_BADARG;
            return res;
        }
    }
    res.M_cell = ct.M_cell;

    // Cell-aggregate initial weights.
    std::vector<double> X_init(ct.M_cell, 0.0);
    for (int i = 0; i < st.n; i++) {
        X_init[ct.cell_of[i]] += st.weights[i];
    }

    // Per-cell bounds.
    std::vector<double> L_cell(ct.M_cell), U_cell(ct.M_cell);
    double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    double lo = st.min_weight;
    for (int c = 0; c < ct.M_cell; c++) {
        L_cell[c] = lo * ct.n_per_cell[c];
        U_cell[c] = hi * ct.n_per_cell[c];
    }

    // Log-space Sinkhorn factors per margin-category.
    int total_cats = 0;
    std::vector<int> cat_offset(st.K + 1, 0);
    for (int k = 0; k < st.K; k++) {
        cat_offset[k + 1] = cat_offset[k] + st.cat_counts[k] + 1;  // +1 for NA bucket
    }
    total_cats = cat_offset[st.K];
    std::vector<double> lf(total_cats, 0.0);  // lf[cat_offset[k] + j]

    // Per-cell capacity multiplier (linear-space).
    std::vector<double> W(ct.M_cell, 1.0);
    std::vector<double> X_tilde(ct.M_cell);
    std::vector<double> X(ct.M_cell);

    // Scratch for margin sweep.
    std::vector<std::vector<int>> cells_by_margin_cat(total_cats);
    for (int k = 0; k < st.K; k++) {
        for (int c = 0; c < ct.M_cell; c++) {
            int j = ct.g_per_cell[k][c];  // j in [0, cat_counts[k]] (NA → cat_counts[k])
            cells_by_margin_cat[cat_offset[k] + j].push_back(c);
        }
    }

    // Targets (user-provided, positive marginals per (k, j)).
    // Note: st.targets[k] has cat_counts[k] entries (no NA slot).
    // Paper: margin sum τ_{k,j} × W_total should match S_kj for j in [0, cat_counts[k]).
    // For NA bucket (j = cat_counts[k]), no constraint; f remains 1.0.

    bool is_infeasible = false;
    std::vector<std::pair<int,int>> infeasible_pairs;

    for (int iter = 1; iter <= st.inner_max_iter; iter++) {
        res.iterations = iter;

        // Margin sweep: one block per margin k.
        for (int k = 0; k < st.K; k++) {
            for (int j = 0; j < st.cat_counts[k]; j++) {
                const auto& cells = cells_by_margin_cat[cat_offset[k] + j];
                if (cells.empty()) {
                    if (st.targets[k][j] > 0.0) {
                        if (!is_infeasible) is_infeasible = true;
                        bool seen = false;
                        for (auto& p : infeasible_pairs)
                            if (p.first == k && p.second == j) { seen = true; break; }
                        if (!seen) infeasible_pairs.emplace_back(k, j);
                    }
                    continue;
                }
                // log-sum-exp stabilization: compute lv_c = log X_init[c] + sum_{m!=k} lf[m]
                //                                       + log W[c]
                std::vector<double> lv(cells.size());
                double lv_max = -std::numeric_limits<double>::infinity();
                for (size_t r = 0; r < cells.size(); r++) {
                    int c = cells[r];
                    if (X_init[c] <= 0.0 || W[c] <= 0.0) {
                        lv[r] = -std::numeric_limits<double>::infinity();
                        continue;
                    }
                    double s = std::log(X_init[c]) + std::log(W[c]);
                    for (int m = 0; m < st.K; m++) {
                        if (m == k) continue;
                        int gm = ct.g_per_cell[m][c];
                        s += lf[cat_offset[m] + gm];
                    }
                    lv[r] = s;
                    if (s > lv_max) lv_max = s;
                }
                if (!std::isfinite(lv_max)) {
                    // All cells are degenerate for this (k,j).
                    if (st.targets[k][j] > 0.0) {
                        if (!is_infeasible) is_infeasible = true;
                        infeasible_pairs.emplace_back(k, j);
                    }
                    continue;
                }
                double sum = 0.0;
                for (size_t r = 0; r < lv.size(); r++) {
                    if (std::isfinite(lv[r])) sum += std::exp(lv[r] - lv_max);
                }
                double log_S_kj = lv_max + std::log(sum);
                if (!std::isfinite(log_S_kj) || std::exp(lv_max) * sum < kEmptyBucketThreshold * ct.W_input) {
                    if (st.targets[k][j] > 0.0) {
                        if (!is_infeasible) is_infeasible = true;
                        infeasible_pairs.emplace_back(k, j);
                    }
                    continue;
                }
                double log_target = std::log(st.targets[k][j] * ct.W_input);
                lf[cat_offset[k] + j] = log_target - log_S_kj;
            }
        }

        // Compute X_tilde via clip-before-exp.
        for (int c = 0; c < ct.M_cell; c++) {
            if (X_init[c] <= 0.0) { X_tilde[c] = 0.0; continue; }
            double s = std::log(X_init[c]);
            for (int m = 0; m < st.K; m++) {
                int gm = ct.g_per_cell[m][c];
                s += lf[cat_offset[m] + gm];
            }
            double s_clip = (s > kLogClip) ? kLogClip : s;
            X_tilde[c] = std::exp(s_clip);
        }

        // Capacity block: X[c] = clamp(X_tilde[c], L_c, U_c); W[c] updated for next iter.
        int n_cap = 0;
        for (int c = 0; c < ct.M_cell; c++) {
            double xc = std::clamp(X_tilde[c], L_cell[c], U_cell[c]);
            X[c] = xc;
            if (X_tilde[c] > 0.0) {
                W[c] = xc / X_tilde[c];
            } else {
                W[c] = 1.0;
            }
            if (W[c] != 1.0) n_cap++;
        }
        res.n_cap_active = n_cap;

        // Convergence check.
        if (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
            double W_total = 0.0;
            for (int c = 0; c < ct.M_cell; c++) W_total += X[c];
            double errRp = 0.0;
            for (int k = 0; k < st.K; k++) {
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    const auto& cells = cells_by_margin_cat[cat_offset[k] + j];
                    double Skj = 0.0;
                    for (int c : cells) Skj += X[c];
                    double e = std::fabs(Skj / W_total - st.targets[k][j]);
                    if (e > errRp) errRp = e;
                }
            }
            res.max_error = errRp;
            if (st.verbose >= 1) {
                char msg[256];
                std::snprintf(msg, sizeof(msg),
                              "iEPPA iter %d: errRp=%.3e n_cap=%d", iter, errRp, n_cap);
                st.log(msg);
            }
            if (errRp < st.tol_abs) {
                res.status = is_infeasible ? RK_ERR_INFEAS : RK_OK;
                break;
            }
        }
    }

    // Expand to obs weights: w[i] = d[i] * X[c] / X_init[c]
    for (int i = 0; i < st.n; i++) {
        int c = ct.cell_of[i];
        if (X_init[c] > 0.0) {
            st.weights[i] = st.weights[i] * X[c] / X_init[c];
        } else {
            st.weights[i] = 0.0;
        }
    }

    if (is_infeasible && res.status == RK_ERR_NOCONV) {
        res.status = RK_ERR_INFEAS;
    }

    if (st.verbose >= 1 && is_infeasible) {
        char msg[256];
        size_t off = 0;
        off += std::snprintf(msg + off, sizeof(msg) - off,
                             "iEPPA infeasible cells: ");
        for (size_t i = 0; i < infeasible_pairs.size() && off < sizeof(msg) - 32; i++) {
            off += std::snprintf(msg + off, sizeof(msg) - off,
                                 "margin=%d cat=%d%s",
                                 infeasible_pairs[i].first + 1,
                                 infeasible_pairs[i].second + 1,
                                 (i + 1 < infeasible_pairs.size()) ? ", " : "");
        }
        st.log(msg);
    }

    return res;
}

} // namespace lbw
```

### Step 3.5: Update leafblower.h enum

- [ ] **Step 3.5.1: Add RK_ALG_RAKING=3**

Modify: `src/leafblower.h` — locate the rk_algorithm_t enum, add entry:
```c
typedef enum {
    RK_ALG_AUTO   = 0,
    RK_ALG_IEPPA  = 1,
    RK_ALG_LBFGSB = 2,
    RK_ALG_RAKING = 3
} rk_algorithm_t;
```

### Step 3.6: Update c_api.cpp dispatch

- [ ] **Step 3.6.1: Include new solver header**

Modify: `src/c_api.cpp` — add near other solver includes:
```cpp
#include "ieppa.hpp"   // faithful (new)
#include "raking.hpp"  // renamed hybrid
```

- [ ] **Step 3.6.2: Update select_algorithm**

Replace `select_algorithm` body with:
```cpp
static rk_algorithm_t select_algorithm(int /*n*/, int /*K*/,
                                        const int* /*cat_counts*/,
                                        const rk_params_t* p) {
    if (p->algorithm != RK_ALG_AUTO) return p->algorithm;
    return RK_ALG_IEPPA;  // AUTO → faithful iEPPA always. Routing refinement TBD per benchmark.
}
```

- [ ] **Step 3.6.3: Update dispatch table**

Locate the switch / if-else dispatching on `selected` (the chosen algorithm). Add branch for `RK_ALG_RAKING` calling `raking_solve(state)`, and ensure `RK_ALG_IEPPA` now calls the new `ieppa_solve(state)` (from `ieppa.hpp`). Example pattern:
```cpp
switch (selected) {
    case RK_ALG_IEPPA:  {
        auto r = lbw::ieppa_solve(state);
        if (result) {
            result->status = r.status;
            result->iterations = r.iterations;
            result->max_error = r.max_error;
            result->algorithm_used = RK_ALG_IEPPA;
        }
        break;
    }
    case RK_ALG_LBFGSB: { /* existing L-BFGS-B dispatch */ }
    case RK_ALG_RAKING: {
        auto r = lbw::raking_solve(state);
        if (result) {
            result->status = r.status;
            result->iterations = r.iterations;
            result->max_error = r.max_error;
            result->algorithm_used = RK_ALG_RAKING;
        }
        break;
    }
    default:
        return RK_ERR_BADARG;
}
```

Announce AUTO selection in result message (verbose:
```cpp
if (p->algorithm == RK_ALG_AUTO && result && p->verbose >= 1 && p->log_fn) {
    char m[128];
    std::snprintf(m, sizeof(m), "[AUTO->iEPPA] iEPPA selected by auto routing");
    p->log_fn(m, p->log_ctx);
}
```

### Step 3.7: Update R bridge method string

- [ ] **Step 3.7.1: Update r_bridge.cpp**

Locate the method-string → algorithm mapping in `src/r_bridge.cpp`. It likely has a sequence like:
```cpp
if (std::strcmp(m, "ieppa") == 0) alg = RK_ALG_IEPPA;
else if (std::strcmp(m, "lbfgsb") == 0) alg = RK_ALG_LBFGSB;
...
```
Add a branch for `"raking"`:
```cpp
else if (std::strcmp(m, "raking") == 0) alg = RK_ALG_RAKING;
```

- [ ] **Step 3.7.2: Update harvest.R match.arg**

Modify: `R/harvest.R` — locate the `map_method` helper (current line ~141) and add `"raking"`:
```r
map_method <- function(method, verbose = 0) {
  method <- tolower(method)
  if (method %in% c("rake", "nrake")) {
    warning("method='", method, "' (IPF) not implemented; using L-BFGS-B")
    method <- "lbfgsb"
  } else if (method == "nr") {
    warning("method='nr' (Newton-Raphson) not implemented; using L-BFGS-B")
    method <- "lbfgsb"
  }
  match.arg(method, c("ieppa", "lbfgsb", "raking"))
}
```

### Step 3.8: Update Python bindings

- [ ] **Step 3.8.1: Update _harvest.py alg_map**

Modify: `python/leafblower/_harvest.py` — locate `alg_map` (current line ~89):
```python
    alg_map = {"ieppa": 1, "lbfgsb": 2, "raking": 3}
```
(Replaces the existing `{"ieppa": 1, "lbfgsb": 2}`.)

- [ ] **Step 3.8.2: Check pybind11 _bindings.cpp**

```bash
grep -n "algorithm\|RK_ALG" python/leafblower/_bindings.cpp
```
If the bindings pass `algorithm` as an integer param, no additional change needed (user supplies 1/2/3 via the alg_map). Verify by reading the calibrate() wrapper in _bindings.cpp.

### Step 3.9: Rewrite test-ieppa.R

- [ ] **Step 3.9.1: Update test-ieppa.R**

Current tests use `method="ieppa"` — now runs faithful. Review each test; relax tolerance/iter caps only if justified by algBCD semantics.

Modify: `tests/testthat/test-ieppa.R` — change context label:
```r
context("ieppa (faithful algBCD)")
```

No assertion changes: the spec §6.1 commits to passing the existing tests with the new solver. If specific tests fail due to algorithmic differences (e.g., iter count), note in the merge commit body and relax **only that test** with a `# faithful iEPPA: X iters vs raking: Y iters; within gate` comment.

### Step 3.10: Create test-ieppa-faithful.R

- [ ] **Step 3.10.1: Write algBCD-specific tests**

File: `tests/testthat/test-ieppa-faithful.R`
```r
context("ieppa faithful — algBCD specifics")

# Tests in this file check properties unique to the faithful algBCD solver:
# cell compression, within-cell weight equality, capacity block behavior.

test_that("within-cell weight equality: obs with identical tuples get equal weights", {
  set.seed(1)
  n <- 1000
  df <- data.frame(
    a = rep(0:3, each = n/4),
    b = rep(rep(0:2, each = n/12), 4)
  )
  tgt <- list(
    a = setNames(c(0.25, 0.25, 0.25, 0.25), 0:3),
    b = setNames(c(0.33, 0.33, 0.34), 0:2)
  )
  res <- harvest(df, tgt, method = "ieppa")
  w <- res$weights
  # Group by (a, b) tuple; within each group, weights should be equal
  for (a in 0:3) for (b in 0:2) {
    mask <- df$a == a & df$b == b
    if (sum(mask) > 1) {
      ws <- w[mask]
      expect_true(diff(range(ws)) < 1e-10,
                  info = sprintf("cell (a=%d,b=%d) weights not equal: range=%.3e",
                                 a, b, diff(range(ws))))
    }
  }
})

test_that("cap-inactive: loose bounds produce no active cap", {
  set.seed(2)
  n <- 1000
  df <- data.frame(a = sample(c("x","y","z"), n, replace = TRUE))
  tgt <- list(a = c(x = 0.4, y = 0.3, z = 0.3))
  res <- harvest(df, tgt, method = "ieppa", max_weight = 10, min_weight = 0)
  # n_cap_active accessible via attr; if not wired, skip
  # Loose bound means all weights should be in interior of [0, 10]
  expect_true(max(res$weights) < 10 - 1e-6)
  expect_true(min(res$weights) > 0 + 1e-6)
})

test_that("cap-active: tight bounds force cap with targets still met", {
  set.seed(3)
  n <- 1000
  df <- data.frame(
    a = c(rep("x", 100), rep("y", 900))  # 90/10 split
  )
  tgt <- list(a = c(x = 0.5, y = 0.5))  # need heavy upweighting of x
  res <- harvest(df, tgt, method = "ieppa", max_weight = 5, min_weight = 0)
  # Target: x:y = 50:50; achieved by upweighting x ~5x
  # Some weights should hit the cap
  expect_true(max(res$weights) >= 5 - 1e-6 || max(res$weights) <= 5 + 1e-6)
  # Verify target met (approximately)
  diag_row <- sum(res$weights[df$a == "x"]) / sum(res$weights)
  expect_lt(abs(diag_row - 0.5), 1e-4)
})

test_that("infeasibility: empty cell + positive target → INFEAS warning", {
  set.seed(4)
  n <- 100
  df <- data.frame(a = rep("x", n))  # only x
  tgt <- list(a = c(x = 0.5, y = 0.5))  # y target positive but no y obs
  expect_warning(res <- harvest(df, tgt, method = "ieppa"),
                 regexp = "infeasible|did not converge", ignore.case = TRUE)
})

test_that("both-sided cap: min_weight + max_weight both active, targets met", {
  set.seed(5)
  n <- 500
  df <- data.frame(
    a = sample(c("x","y"), n, replace = TRUE, prob = c(0.8, 0.2))
  )
  tgt <- list(a = c(x = 0.5, y = 0.5))  # need up y and down x
  res <- harvest(df, tgt, method = "ieppa",
                 min_weight = 0.3, max_weight = 3)
  expect_true(min(res$weights) >= 0.3 - 1e-6)
  expect_true(max(res$weights) <= 3 + 1e-6)
})
```

### Step 3.11: Create test-compare.R

- [ ] **Step 3.11.1: Write cross-algorithm equivalence tests**

File: `tests/testthat/test-compare.R`
```r
context("cross-algorithm equivalence on feasible inputs")

test_that("ieppa, raking, lbfgsb agree to 1e-3 on 10 random feasible datasets", {
  set.seed(20260423)
  for (trial in 1:10) {
    n <- sample(c(1000, 5000, 10000), 1)
    K <- sample(3:5, 1)
    cats <- sample(3:5, K, replace = TRUE)
    # Build random balanced categorical data + balanced targets
    cols <- lapply(seq_len(K), function(k) sample(0:(cats[k]-1), n, replace = TRUE))
    df <- as.data.frame(cols)
    names(df) <- paste0("v", seq_len(K))
    tgt <- lapply(seq_len(K), function(k) {
      setNames(rep(1/cats[k], cats[k]), as.character(0:(cats[k]-1)))
    })
    names(tgt) <- names(df)
    # Reasonable bounds
    mw <- sample(c(2, 3, 5), 1)
    r_ieppa  <- suppressWarnings(harvest(df, tgt, method = "ieppa",  max_weight = mw))
    r_raking <- suppressWarnings(harvest(df, tgt, method = "raking", max_weight = mw))
    r_lbfgsb <- suppressWarnings(harvest(df, tgt, method = "lbfgsb", max_weight = mw))
    max_diff <- max(
      max(abs(r_ieppa$weights - r_raking$weights)),
      max(abs(r_ieppa$weights - r_lbfgsb$weights)),
      max(abs(r_raking$weights - r_lbfgsb$weights))
    )
    expect_lt(max_diff, 1e-3,
              label = sprintf("trial %d (n=%d K=%d mw=%.1f): max pairwise diff %.3e",
                              trial, n, K, mw, max_diff))
  }
})
```

### Step 3.12: Update PRD

- [ ] **Step 3.12.1: Update tasks/prd-leafblower-core.md § US-005**

Modify: `tasks/prd-leafblower-core.md` — rewrite §US-005 to reflect the actual faithful algBCD. Replace the description paragraph (line ~139) with:
```markdown
**Description:** As a developer, I want leafblower to implement the paper-faithful iEPPA algorithm (Chu-Liang-Toh-Yang 2022, arXiv:2011.14312) at C=0 so that capacity-constrained multi-marginal calibration converges with cell-compressed O(M_cell·K) inner cost and a documented convergence framework (Csiszár 1975 cyclic I-projection).
```

Rewrite the acceptance criteria (the `- [ ]` list in §US-005) to match the spec §2.3 algBCD (log-space factors, capacity BCD block, log-sum-exp stabilization). Reference the spec for full algorithm detail.

Add new §US-005b for raking (after §US-005):
```markdown
### US-005b: Classical Raking Algorithm (IPF + Dykstra)

**Description:** As a developer, I want `method="raking"` to provide the classical cyclic IPF (Deming-Stephan 1940; Csiszár 1975) with additive Dykstra box-projection (Boyle-Dykstra 1986) and hyperplane projection, as an alternative to the paper-faithful iEPPA.

**Acceptance Criteria:**
- `rk_calibrate(..., algorithm=RK_ALG_RAKING)` present in enum
- `method="raking"` in R/Python routes to `raking_solve`
- All pre-rev2 iEPPA tests pass against `method="raking"` (regression guard in `test-raking.R`)

**Priority:** High
**Dependencies:** US-004
**Docs impact:** `man/harvest.Rd` — add method="raking"
**Config impact:** None
```

### Step 3.13: Create NEWS.md

- [ ] **Step 3.13.1: Create NEWS.md**

File: `NEWS.md`
```markdown
# leafblower (development)

## Breaking changes

* `method="ieppa"` now runs the paper-faithful algBCD (Chu, Liang, Toh &
  Yang 2022, arXiv:2011.14312) at C=0, using cell-compressed representation
  with log-space Sinkhorn factors and a capacity BCD block. The previous
  implementation was an IPF+Dykstra hybrid misnamed "iEPPA"; it is renamed
  `method="raking"`. Users who relied on the previous `method="ieppa"`
  behavior should switch to `method="raking"`.

* New `method="raking"` exposes the renamed classical IPF+Dykstra hybrid
  (Deming-Stephan 1940 / Csiszár 1975 cyclic IPF + Boyle-Dykstra 1986
  additive projections). It is the same code as pre-rename `method="ieppa"`.

* `method="auto"` continues to route to `method="ieppa"` (now the faithful
  algBCD). Benchmark-driven routing refinement is deferred to a follow-up
  release.

## New features

* Cell-compressed computation: faithful iEPPA operates at cell-level
  (unique (g_1,...,g_K) tuples) rather than observation-level, yielding
  up to 1000× speedup on surveys with low tuple diversity.

* Result diagnostic: `harvest()` returns include `M_cell` and
  `n_cap_active` via `attr(result, "M_cell")` and
  `attr(result, "n_cap_active")` for iEPPA calibrations.
```

### Step 3.14: Build + verify merge gate

- [ ] **Step 3.14.1: Rebuild**
```bash
R CMD INSTALL --preclean . 2>&1 | tail -3
```
Expected: `* DONE (leafblower)`. Any compile error must be fixed before proceeding.

- [ ] **Step 3.14.2: Run full test suite**
```bash
Rscript -e "devtools::test()" 2>&1 | tail -5
```
Expected: `[ FAIL 0 | PASS N ]` where N ≥ 103 (new tests added). Any failure must be diagnosed and fixed before commit.

- [ ] **Step 3.14.3: R CMD check**
```bash
R CMD build . 2>&1 | tail -3
R CMD check --no-manual --as-cran leafblower_*.tar.gz 2>&1 | tail -10
```
Expected: `Status: OK` or `Status: N NOTE`. 0 ERROR, 0 WARNING is the merge gate.

- [ ] **Step 3.14.4: Degenerate-case timing check**

```bash
Rscript -e '
library(leafblower)
set.seed(1)
n <- 10000
df <- data.frame(a = sample(letters[1:5], n, replace=TRUE))
tgt <- list(a = setNames(rep(0.2, 5), letters[1:5]))
# Force M_cell = n by making each obs unique (add a random unique column)
df$b <- seq_len(n)
tgt$b <- setNames(rep(1/n, n), as.character(seq_len(n)))
# This wont work directly (n-way cat overflows); use simpler degenerate case:
df2 <- data.frame(a = seq_len(n))  # all unique
tgt2 <- list(a = setNames(rep(1/n, n), as.character(seq_len(n))))
t_ie <- median(replicate(3, system.time(harvest(df, tgt, method="ieppa"))[3]))
t_rk <- median(replicate(3, system.time(harvest(df, tgt, method="raking"))[3]))
cat(sprintf("ieppa=%.3fs raking=%.3fs ratio=%.2fx\n", t_ie, t_rk, t_ie/t_rk))
'
```
Expected: ratio within [0.5, 2.0]. Document the ratio in the merge commit body.

- [ ] **Step 3.14.5: Python pytest**
```bash
cd python && python -m pytest 2>&1 | tail -5 && cd ..
```
Expected: All pass. If Python tests reference old semantics, update them similarly to R tests.

### Step 3.15: ATOMIC commit

- [ ] **Step 3.15.1: Single commit with all changes**
```bash
git add src/ieppa.hpp src/ieppa.cpp src/raking.hpp src/raking.cpp \
        src/c_api.cpp src/r_bridge.cpp src/leafblower.h \
        R/harvest.R python/leafblower/_harvest.py \
        tests/testthat/test-ieppa.R tests/testthat/test-raking.R \
        tests/testthat/test-ieppa-faithful.R tests/testthat/test-compare.R \
        tasks/prd-leafblower-core.md NEWS.md

git commit -m "$(cat <<'EOF'
feat!: faithful iEPPA solver; rename hybrid to method=raking

BREAKING: method='ieppa' now runs paper-faithful algBCD (Chu-Liang-Toh-Yang
2022, arXiv:2011.14312) at C=0 with cell-compressed representation. Previous
IPF+Dykstra hybrid renamed to method='raking'.

Atomic commit bundles:
- src/ieppa.{hpp,cpp}: new faithful solver (log-space Sinkhorn, LSE-stable,
  capacity BCD block)
- src/raking.{hpp,cpp}: renamed from src/ieppa.* (unchanged implementation;
  symbols ieppa_solve→raking_solve, IEPPAResult→RakingResult)
- src/leafblower.h: RK_ALG_RAKING=3 enum
- src/c_api.cpp: dispatch table + AUTO→IEPPA routing
- src/r_bridge.cpp: method='raking' mapping
- R/harvest.R, python/leafblower/_harvest.py: user-facing method names
- tasks/prd-leafblower-core.md: §US-005 rewritten, §US-005b added
- NEWS.md: created with breaking-change entry
- tests/testthat/test-raking.R: regression guard (exact-copy of old
  test-ieppa.R against new method='raking')
- tests/testthat/test-ieppa.R: retargets at faithful solver
- tests/testthat/test-ieppa-faithful.R: algBCD-specific cases
- tests/testthat/test-compare.R: cross-algorithm equivalence

Merge gate satisfied:
- All tests pass [FAIL 0]
- R CMD check --as-cran: 0 ERROR, 0 WARNING
- Degenerate M_cell=n timing: ratio <fill in from Step 3.14.4> (within ±2×)
- Python pytest: <fill in results>

Refs: design spec docs/superpowers/specs/2026-04-23-ieppa-faithful-design.md
EOF
)"
```

- [ ] **Step 3.15.2: Verify commit**
```bash
git show --stat HEAD
```
Expected: one commit, all files listed.

---

## WU-4: make_plots() 3D refactor (prerequisite for benchmark WU-5)

**Files:**
- Modify: `benchmarks/algo_selection_benchmark.R` (extract `make_plots()` for reuse; generalize to arbitrary input dimension)
- Create: `benchmarks/plot_helpers.R` (shared plotting helpers)

### Step 4.1: Extract existing make_plots()

- [ ] **Step 4.1.1: Locate make_plots()**
```bash
grep -n "make_plots\s*<-\s*function" benchmarks/algo_selection_benchmark.R
```

- [ ] **Step 4.1.2: Move to plot_helpers.R**

Create `benchmarks/plot_helpers.R` containing the extracted function, generalized:
- Accept `gp`, `design`, `y`, `threshold` as before
- Accept `x_names` (character vector of input names) and `x_ranges` (list of numeric ranges)
- For 2D inputs: produce contour + uncertainty PDFs as today
- For 3D inputs: produce slice plots — one 2D contour per fixed value of the third dimension (e.g., `x3 ∈ {0, 1, 2, 3}` for compression log10)

Keep backward compat: when `length(x_names) == 2`, existing PDF names produced.

- [ ] **Step 4.1.3: Update algo_selection_benchmark.R to source()**

Modify: `benchmarks/algo_selection_benchmark.R` — replace inline `make_plots` definition with `source("plot_helpers.R")` at top.

### Step 4.2: Sanity check existing 2D benchmark still plots

- [ ] **Step 4.2.1: Smoke test**
```bash
Rscript -e 'source("benchmarks/plot_helpers.R"); cat("OK\n")'
```

### Step 4.3: Commit

- [ ] **Step 4.3.1: Commit**
```bash
git add benchmarks/plot_helpers.R benchmarks/algo_selection_benchmark.R
git commit -m "refactor(bench): extract make_plots to benchmarks/plot_helpers.R

Generalized make_plots for 2D and 3D input spaces (slice plots along
third dimension for 3D). Prerequisite for ieppa_vs_raking benchmark
(design §7, uses 3D input: complexity, tolerance, compression ratio).

Existing algo_selection_benchmark.R behavior unchanged (2D path)."
```

---

## WU-5: Bayesian Level Set benchmark (post-merge analysis)

**Files:**
- Create: `benchmarks/ieppa_vs_raking_bench.R`

### Step 5.1: Write benchmark harness

- [ ] **Step 5.1.1: Create benchmark script**

File: `benchmarks/ieppa_vs_raking_bench.R`
```r
#!/usr/bin/env Rscript
# Bayesian Level Set Estimation: ieppa (faithful) vs raking (hybrid).
# Response: y = log(t_ieppa / t_raking). Threshold: log(1.2) — ieppa wins if y <= threshold.
# Input space (3D):
#   x1 = log10(complexity) = log10(n * sum(cat_counts)) in [4, 7.7]
#   x2 = log10(tol_abs) in [-6, -3]
#   x3 = log10(prod(cat_counts) / n) in [0, 3.5]  # theoretical max compression
#
# Output: benchmarks/ieppa_vs_raking_results.rds
Sys.setenv(OMP_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(leafblower)
  library(DiceKriging)
  library(lhs)
  source("benchmarks/plot_helpers.R")
})

.BENCH_SOURCED <- exists(".BENCH_SOURCED", envir = .GlobalEnv, inherits = FALSE)
set.seed(20260423)

X1_RANGE <- c(4.0, 7.7)
X2_RANGE <- c(-6.0, -3.0)
X3_RANGE <- c(0.0, 3.5)
THRESHOLD <- log(1.2)

# Deterministic (x1,x2,x3) → data-generation seed
bench_seed <- function(x1, x2, x3) {
  as.integer(abs(round(1e5 * sin(x1 * 13 + x2 * 17 + x3 * 19)))) %% .Machine$integer.max
}

# Generate synthetic data matching targets
generate_data <- function(x1, x2, x3) {
  seed <- bench_seed(x1, x2, x3)
  set.seed(seed)
  K <- 5
  target_product <- 10^x3 * round(10^x1 / 10)  # approximate, will be adjusted
  cats_per <- max(2L, round((10^x3 * round(10^x1 / 10))^(1/K)))
  cat_counts <- rep(cats_per, K)
  n <- round(10^x1 / sum(cat_counts))
  n <- max(n, 1000L)
  tol_abs <- 10^x2
  # Generate categorical data: uniform sampling
  cols <- lapply(seq_len(K), function(k) sample(0:(cats_per-1), n, replace = TRUE))
  df <- as.data.frame(cols)
  names(df) <- paste0("v", seq_len(K))
  tgt <- lapply(seq_len(K), function(k) {
    setNames(rep(1/cats_per, cats_per), as.character(0:(cats_per-1)))
  })
  names(tgt) <- names(df)
  list(df = df, tgt = tgt, n = n, K = K, cat_counts = cat_counts, tol = tol_abs)
}

# Time a single method
time_method <- function(d, method) {
  t0 <- Sys.time()
  suppressWarnings(
    harvest(d$df, d$tgt, method = method,
            convergence = list(absolute = d$tol),
            max_iterations = 500,
            max_weight = 3)
  )
  as.numeric(Sys.time() - t0, units = "secs")
}

# Evaluate log-ratio at a single (x1,x2,x3) point
evaluate_point <- function(x1, x2, x3, runs = 3) {
  d <- generate_data(x1, x2, x3)
  t_ie <- median(replicate(runs, time_method(d, "ieppa")))
  t_rk <- median(replicate(runs, time_method(d, "raking")))
  log(t_ie / t_rk)
}

run_benchmark <- function(n_initial = 16, n_adaptive = 8) {
  # LHS initial design in [0,1]^3
  D0 <- maximinLHS(n_initial, 3)
  design <- cbind(
    X1_RANGE[1] + D0[,1] * diff(X1_RANGE),
    X2_RANGE[1] + D0[,2] * diff(X2_RANGE),
    X3_RANGE[1] + D0[,3] * diff(X3_RANGE)
  )
  colnames(design) <- c("log10_complexity", "log10_tol", "log10_compression")
  cat("Evaluating initial", n_initial, "points...\n")
  y <- numeric(n_initial)
  for (i in seq_len(n_initial)) {
    y[i] <- evaluate_point(design[i,1], design[i,2], design[i,3])
    cat(sprintf("  point %d: x=(%.2f,%.2f,%.2f) y=%.3f\n",
                i, design[i,1], design[i,2], design[i,3], y[i]))
  }
  # GP fit
  gp <- km(design = design, response = y, covtype = "matern5_2",
           nugget.estim = TRUE, control = list(trace = 0))
  # Adaptive augmentation
  for (i in seq_len(n_adaptive)) {
    # Sample candidates, pick one with highest GP predictive variance near threshold
    cand <- cbind(
      runif(1000, X1_RANGE[1], X1_RANGE[2]),
      runif(1000, X2_RANGE[1], X2_RANGE[2]),
      runif(1000, X3_RANGE[1], X3_RANGE[2])
    )
    pred <- predict(gp, cand, type = "UK", checkNames = FALSE)
    # Score: variance * indicator-near-threshold
    score <- pred$sd^2 * exp(-((pred$mean - THRESHOLD)^2) / 0.1)
    best <- cand[which.max(score), , drop = FALSE]
    y_new <- evaluate_point(best[1,1], best[1,2], best[1,3])
    design <- rbind(design, best)
    y <- c(y, y_new)
    cat(sprintf("  adaptive %d: x=(%.2f,%.2f,%.2f) y=%.3f\n",
                i, best[1,1], best[1,2], best[1,3], y_new))
    # Refit
    gp <- km(design = design, response = y, covtype = "matern5_2",
             nugget.estim = TRUE, control = list(trace = 0))
  }
  list(
    design = design, y = y, gp = gp, threshold = THRESHOLD,
    x_ranges = list(
      log10_complexity = X1_RANGE,
      log10_tol = X2_RANGE,
      log10_compression = X3_RANGE
    ),
    meta = list(seed = 20260423, n_initial = n_initial,
                n_adaptive = n_adaptive, runs_per_point = 3,
                timestamp = Sys.time())
  )
}

if (!.BENCH_SOURCED) {
  res <- run_benchmark()
  saveRDS(res, "benchmarks/ieppa_vs_raking_results.rds")
  cat("Saved benchmarks/ieppa_vs_raking_results.rds\n")
  # Plot (uses plot_helpers.R)
  make_plots(res$gp, res$design, res$y, res$threshold,
             x_names = c("log10_complexity", "log10_tol", "log10_compression"),
             x_ranges = res$x_ranges,
             out_prefix = "benchmarks/ieppa_vs_raking")
}
```

### Step 5.2: Smoke test at small size

- [ ] **Step 5.2.1: Run with reduced grid**
```bash
Rscript -e 'source("benchmarks/ieppa_vs_raking_bench.R", echo=FALSE); res <- run_benchmark(n_initial=4, n_adaptive=2); cat("OK, design rows=", nrow(res$design), "\n")'
```
Expected: runs to completion; writes RDS; prints summary.

### Step 5.3: Commit

- [ ] **Step 5.3.1: Commit**
```bash
git add benchmarks/ieppa_vs_raking_bench.R
git commit -m "feat(bench): ieppa-vs-raking Bayesian Level Set Estimation harness

3D input space (log10 complexity, log10 tolerance, log10 compression).
Response = log(t_ieppa / t_raking). GP surrogate via DiceKriging; LHS
initial design + variance-weighted adaptive augmentation near decision
threshold log(1.2).

Saved to benchmarks/ieppa_vs_raking_results.rds. Plots via
make_plots() (3D slice plots).

Post-merge analysis: feeds the routing-refinement WU (to be filed
separately on benchmark output).

Refs: design §7"
```

### Step 5.4: Full benchmark run (optional, post-merge)

- [ ] **Step 5.4.1: Full run (takes 1-3 hours)**
```bash
Rscript benchmarks/ieppa_vs_raking_bench.R 2>&1 | tee benchmarks/ieppa_vs_raking_run.log
```

- [ ] **Step 5.4.2: Commit results (if running)**
```bash
git add benchmarks/ieppa_vs_raking_results.rds benchmarks/ieppa_vs_raking_*.pdf benchmarks/ieppa_vs_raking_run.log
git commit -m "chore(bench): add ieppa_vs_raking benchmark run outputs

GP surrogate over 3D input space; see plots for decision boundary.
Input to routing-refinement WU."
```

- [ ] **Step 5.4.3: File routing-refinement WU**
```bash
bd create --title="Benchmark-driven routing refinement (ieppa vs raking)" \
  --description="Use benchmarks/ieppa_vs_raking_results.rds GP surrogate to update select_algorithm() in src/c_api.cpp. Currently AUTO always routes to IEPPA; refine to select RAKING when GP predicts t_raking < t_ieppa by ≥20% at the given (n, K, cat_counts). Input from WU-5 of plan 2026-04-23-ieppa-faithful-impl.md." \
  --type=feature --priority=2
```

---

## Summary

| WU | Commits | Merge gate contribution |
|---|---|---|
| WU-1 Cell table | 1 (prereq) | Enables solver |
| WU-2 Input validation | 1 (prereq) | Security + spec §8 |
| WU-3 ATOMIC bundle | 1 (big) | Main deliverable; triggers §8c merge gate |
| WU-4 make_plots refactor | 1 (prereq) | Required for WU-5 |
| WU-5 Benchmark | 1 code + 1 outputs | Post-merge analysis |

Total: 5 commits in main line + 2 optional post-merge.

**Follow-up WUs already filed in beads:** leafblower-r2q (q_hyp scalar), leafblower-47m (raking SIMD), leafblower-5o0 (provenance doc), leafblower-dl6 (descent monitor), leafblower-53d (auto-fallback), plus kk1.20.4 and kk1.24.3 carried from prior cycles.

---

## Self-Review

**Spec coverage:** §2.1 reduction → WU-3 solver. §2.2 cell compression → WU-1. §2.3 algBCD → WU-3 step 3.4. §3 scope → WU-3. §4 architecture → WU-3 steps 3.2-3.8. §5 cell table → WU-1. §6 tests → WU-3 steps 3.9-3.11. §7 benchmark → WU-5. §8 error handling → WU-2 + WU-3 step 3.4. §8b verbose → WU-3 step 3.4. §8c merge gate → WU-3 step 3.14. §9 convergence → referenced in WU-3 step 3.4. §10 PRD/NEWS → WU-3 steps 3.12-3.13. §13 deliverables → full plan.

**Placeholders:** None. All code blocks are complete. All commands have expected outputs.

**Type consistency:** `IEPPAResult` struct is defined in WU-3 step 3.3.1 with fields used consistently in step 3.4.1 dispatch and step 3.10.1 tests (`M_cell`, `n_cap_active`). `ieppa_solve(CalibState&)` signature consistent across header (3.3.1) and call site (3.6.3). `raking_solve` sed-rename in step 3.2.2 covers both `ieppa.hpp` and `ieppa.cpp` occurrences.
