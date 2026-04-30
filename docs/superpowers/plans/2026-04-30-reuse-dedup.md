# Code Reuse & Deduplication — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extract 7 repeated C++ patterns from 6 calibration solvers into shared helpers in `calib_dispatch.hpp` and `cell_table.hpp`.

**Architecture:** All helpers go into `calib_dispatch.hpp` (free functions in `lbw` namespace) or `cell_table.hpp` (CellTable member/free function). No new files. Each task is one helper extraction + call-site replacement + compile gate.

**Tech Stack:** C++17, R package build (`R CMD INSTALL --preclean .`), testthat via `Rscript`

**Mechanism:** Header-only free functions in `lbw` namespace — inline, zero ABI impact.
**Forbidden:** New `.cpp` files; separate translation units for helpers; touching convergence logic or algorithm correctness.
**Audit:** R-level regression tests per solver (each task confirms same numeric output before/after).

---

## Task 1 — fcbo.1: Extract `solver_setup_ct()` helper

**Ticket:** fcbo.1

**Pattern to extract:** The build-cell-table preamble common to every solver entry point:
1. `build_cell_table(...)` call and error-wrap
2. `X_init` aggregation from `st.weights`
3. `resolve_hi` (see Task 2, but inline the scalar here first)
4. `L_cell` / `U_cell` allocation and fill
5. `cat_offset` build + `n_cats_total`
6. `calib_validate_preentry(...)` call and error-wrap

**Sites:**
| File | Lines (approx) |
|------|---------------|
| `src/chebyshev.cpp` | 38–116 (build_cell_table → validate; see note below) |
| `src/sinkhorn.cpp` | 96–117 |
| `src/greg.cpp` | 18–49 |
| `src/greenkhorn.cpp` | 22–62 |
| `src/logit_calib.cpp` | 35–76 |

**Note — chebyshev.cpp preamble boundary:** Lines 35–36 are a solver-specific metric override (`if (st.convergence_cfg.metric == L1_WEIGHT) st.convergence_cfg.metric = MAX_ERR`) that must stay inline before the helper call. The extractable preamble is lines 38–116 only (build_cell_table → L_cell/U_cell → cat_offset → calib_validate_preentry). Lines 118–124 are a warm-start override that aggregates ieppa obs-level weights into X_init; this must also remain inline after the helper call. The X_init construction in the helper (fresh aggregate from st.weights) is correct and matches all 5 solvers: sinkhorn (lines 104–106) builds X fresh from st.weights then copies to X_init — identical semantics to the helper skeleton.

**Note:** `raking.cpp` splits the preamble and omits `cat_offset` (no linear system) — it shares steps 1–2 and 4 only. Use `solver_setup_ct_base` (see second overload below) for raking; leave cat_offset inline in raking.

**Actual signatures (from reading source):**

```cpp
// calib_dispatch.hpp — namespace lbw
//
// solver_setup_ct: performs the standard solver entry preamble.
//
// Out-params: ct (populated), X_init (M_cell), L_cell, U_cell, cat_offset, n_cats_total.
// hi_eff: resolved hi bound = isfinite(st.max_weight) ? st.max_weight : 1e300.
//
// Returns RK_OK (0) on success; on error fills err_status + err_message and returns != RK_OK.
// The caller must immediately propagate err_status + err_message to its own result struct
// and return when solver_setup_ct returns != RK_OK.
//
// Template param ResT: any result struct that has .status (int) and .message (char[]).
template <typename ResT>
inline int solver_setup_ct(
    CalibState&              st,
    CellTable&               ct,
    std::vector<double>&     X_init,     // out: size ct.M_cell after call
    double&                  hi_eff,     // out: resolved hi bound
    std::vector<double>&     L_cell,     // out: size ct.M_cell
    std::vector<double>&     U_cell,     // out: size ct.M_cell
    std::vector<int>&        cat_offset, // out: size st.K
    int&                     n_cats_total,
    ResT&                    res) noexcept;

// solver_setup_ct_base: narrower overload for solvers that do NOT use a linear system
// (no cat_offset, no n_cats_total, no calib_validate_preentry call).
// Use this for raking.cpp and any future solver without a constraint matrix.
template <typename ResT>
inline int solver_setup_ct_base(
    CalibState&              st,
    CellTable&               ct,
    std::vector<double>&     X_init,     // out: size ct.M_cell after call
    double&                  hi_eff,     // out: resolved hi bound
    std::vector<double>&     L_cell,     // out: size ct.M_cell
    std::vector<double>&     U_cell,     // out: size ct.M_cell
    ResT&                    res) noexcept;
```

**Implementation skeleton:**

```cpp
template <typename ResT>
inline int solver_setup_ct(
    CalibState& st, CellTable& ct,
    std::vector<double>& X_init, double& hi_eff,
    std::vector<double>& L_cell, std::vector<double>& U_cell,
    std::vector<int>& cat_offset, int& n_cats_total,
    ResT& res) noexcept
{
    if (build_cell_table(st.n, st.K, st.group_ids, st.cat_counts, st.weights, ct) != 0) {
        res.status = RK_ERR_BADARG;
        std::strncpy(res.message, "build_cell_table failed", sizeof(res.message) - 1);
        return RK_ERR_BADARG;
    }
    X_init.assign(ct.M_cell, 0.0);
    for (int i = 0; i < st.n; i++) X_init[ct.cell_of[i]] += st.weights[i];
    hi_eff = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    L_cell.resize(ct.M_cell); U_cell.resize(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        L_cell[c] = st.min_weight * static_cast<double>(ct.n_per_cell[c]);
        U_cell[c] = hi_eff      * static_cast<double>(ct.n_per_cell[c]);
    }
    cat_offset.resize(st.K); n_cats_total = 0;
    for (int k = 0; k < st.K; k++) { cat_offset[k] = n_cats_total; n_cats_total += st.cat_counts[k]; }
    rk_result_t tmp = {};
    if (calib_validate_preentry(ct, st, &tmp, X_init.data(), n_cats_total) != RK_OK) {
        res.status = tmp.status;
        std::strncpy(res.message, tmp.message, sizeof(res.message) - 1);
        return tmp.status;
    }
    return RK_OK;
}
```

**`solver_setup_ct_base` implementation skeleton (raking / no-linear-system variant):**

```cpp
template <typename ResT>
inline int solver_setup_ct_base(
    CalibState& st, CellTable& ct,
    std::vector<double>& X_init, double& hi_eff,
    std::vector<double>& L_cell, std::vector<double>& U_cell,
    ResT& res) noexcept
{
    if (build_cell_table(st.n, st.K, st.group_ids, st.cat_counts, st.weights, ct) != 0) {
        res.status = RK_ERR_BADARG;
        std::strncpy(res.message, "build_cell_table failed", sizeof(res.message) - 1);
        return RK_ERR_BADARG;
    }
    X_init.assign(ct.M_cell, 0.0);
    for (int i = 0; i < st.n; i++) X_init[ct.cell_of[i]] += st.weights[i];
    hi_eff = std::isfinite(st.max_weight) ? st.max_weight : 1e300;
    L_cell.resize(ct.M_cell); U_cell.resize(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        L_cell[c] = st.min_weight * static_cast<double>(ct.n_per_cell[c]);
        U_cell[c] = hi_eff      * static_cast<double>(ct.n_per_cell[c]);
    }
    // NOTE: no cat_offset, no calib_validate_preentry — raking does not use a linear system.
    return RK_OK;
}
```

**Required includes in calib_dispatch.hpp:** `"calib_validate.hpp"` (already transitively included via `leafblower.h`? — verify before adding).

**Steps:**
- [ ] Read `src/calib_validate.hpp` to confirm `calib_validate_preentry` signature.
- [ ] Add `solver_setup_ct` template to `src/calib_dispatch.hpp` (after existing helpers, before closing `}`).
- [ ] Add `solver_setup_ct_base` template to `src/calib_dispatch.hpp` (immediately after `solver_setup_ct`).
- [ ] Replace preamble in `src/sinkhorn.cpp` with `solver_setup_ct` call; delete displaced lines.
- [ ] `R CMD INSTALL --preclean .` — must pass.
- [ ] Replace preamble in `src/greg.cpp`; compile.
- [ ] Replace preamble in `src/greenkhorn.cpp`; compile.
- [ ] Replace preamble in `src/logit_calib.cpp`; compile.
- [ ] Replace preamble in `src/chebyshev.cpp`: keep lines 35–36 (metric override) before the helper call; replace lines 38–116 with `solver_setup_ct` call; keep lines 118–124 (warm-start override) after the helper call; compile.
- [ ] Read `src/raking.cpp` lines 73–115 to identify exact preamble extent.
- [ ] Replace `raking.cpp` preamble (build_cell_table + X_init + L_cell/U_cell, approx lines 73–92) with `solver_setup_ct_base` call; SKIP cat_offset step — it stays inline in raking; compile.
- [ ] `Rscript -e 'testthat::test_dir("tests/testthat")'` — all tests green.
- [ ] `git commit -m "refactor(solvers): extract solver_setup_ct() preamble helper — fcbo.1"`

---

## Task 2 — fcbo.2: Extract `lbw::resolve_hi()` + `lbw::compute_cell_bounds()`

**Ticket:** fcbo.2

**Pattern:** Two 1–3 line idioms repeated in every solver:

```cpp
// resolve_hi — 6 sites
const double hi = std::isfinite(st.max_weight) ? st.max_weight : 1e300;

// compute_cell_bounds — 6 sites (loop)
for (int c = 0; c < ct.M_cell; c++) {
    L_cell[c] = lo * static_cast<double>(ct.n_per_cell[c]);
    U_cell[c] = hi * static_cast<double>(ct.n_per_cell[c]);
}
```

**Note:** If Task 1 is done first, `solver_setup_ct` already inlines both. These helpers are still worth extracting for callers that only need the bounds (e.g., post-loop obs-expansion in `chebyshev.cpp:837`). Task 2 is independent of Task 1 and can be done before or after.

**Signatures:**

```cpp
namespace lbw {

/// Returns st.max_weight if finite, else 1e300.
inline double resolve_hi(const CalibState& st) noexcept {
    return std::isfinite(st.max_weight) ? st.max_weight : 1e300;
}

/// Fills L[c] = lo*n_per_cell[c], U[c] = hi*n_per_cell[c] for c in [0, M_cell).
inline void compute_cell_bounds(
    const CellTable& ct,
    double lo, double hi,
    std::vector<double>& L,
    std::vector<double>& U) noexcept
{
    L.resize(ct.M_cell); U.resize(ct.M_cell);
    for (int c = 0; c < ct.M_cell; c++) {
        L[c] = lo * static_cast<double>(ct.n_per_cell[c]);
        U[c] = hi * static_cast<double>(ct.n_per_cell[c]);
    }
}

} // namespace lbw
```

**Files modified:**
- `src/calib_dispatch.hpp` — add two helpers (after `apply_obs_expansion`).
- `src/chebyshev.cpp` — replace lines 46–52 and line 837 `hi_obs` computation with `resolve_hi(st)`.
- `src/sinkhorn.cpp` — replace lines 119–125.
- `src/greg.cpp` — replace lines 29–35.
- `src/greenkhorn.cpp` — replace lines 65–68.
- `src/logit_calib.cpp` — replace lines 72–76.
- `src/raking.cpp` — replace lines 109–115.

**Steps:**
- [ ] Add `resolve_hi` and `compute_cell_bounds` to `src/calib_dispatch.hpp`.
- [ ] Replace each site one file at a time; `R CMD INSTALL --preclean .` after each file.
- [ ] `Rscript -e 'testthat::test_dir("tests/testthat")'` — green.
- [ ] `git commit -m "refactor(solvers): extract resolve_hi() + compute_cell_bounds() — fcbo.2"`

---

## Task 3 — fcbo.3: Extract `lbw::build_cat_offset()` + `lbw::max_cats_count()`

**Ticket:** fcbo.3

**Pattern:**

```cpp
// build_cat_offset — 4 sites: chebyshev, greg, logit_calib, plus inline in sinkhorn
std::vector<int> cat_offset(st.K);
int n_cats_total = 0;
for (int k = 0; k < st.K; k++) { cat_offset[k] = n_cats_total; n_cats_total += st.cat_counts[k]; }

// max_cats_count — 5 sites (each solver computes max_element over cat_counts)
int max_cats = *std::max_element(st.cat_counts, st.cat_counts + st.K);
```

**Signatures:**

```cpp
namespace lbw {

/// Builds prefix-sum cat_offset[k] and returns n_cats_total.
/// cat_offset[k] = sum of cat_counts[0..k-1].
inline int build_cat_offset(int K, const int* cat_counts,
                             std::vector<int>& cat_offset) noexcept
{
    cat_offset.resize(K);
    int nct = 0;
    for (int k = 0; k < K; k++) { cat_offset[k] = nct; nct += cat_counts[k]; }
    return nct;
}

/// Returns max(cat_counts[0..K-1]), or 0 when K==0.
inline int max_cats_count(int K, const int* cat_counts) noexcept {
    if (K == 0) return 0;
    return *std::max_element(cat_counts, cat_counts + K);
}

} // namespace lbw
```

**Files modified:**
- `src/calib_dispatch.hpp` — add helpers.
- `src/chebyshev.cpp` — replace lines 54–56 and line 199.
- `src/sinkhorn.cpp` — replace line 132 (`max_cats`).
- `src/greg.cpp` — replace lines 38–40 and line 54.
- `src/greenkhorn.cpp` — replace lines 35–36.
- `src/logit_calib.cpp` — replace lines 51–52 and lines 63–65.
- `src/raking.cpp` — replace lines 126.

**Steps:**
- [ ] Add `build_cat_offset` and `max_cats_count` to `src/calib_dispatch.hpp`.
- [ ] Replace each site; `R CMD INSTALL --preclean .` after each file.
- [ ] `Rscript -e 'testthat::test_dir("tests/testthat")'` — green.
- [ ] `git commit -m "refactor(solvers): extract build_cat_offset() + max_cats_count() — fcbo.3"`

---

## Task 4 — fcbo.4: Extract `lbw::build_cells_per_cat()` to `cell_table.hpp`

**Ticket:** fcbo.4

**Pattern:** Three solvers build an identical `cells_per_cat[k][j]` index structure:

```cpp
// greenkhorn.cpp:39-46, raking.cpp:90-98, logit_calib.cpp:53-60
std::vector<std::vector<std::vector<int>>> cells_per_cat(K);
for (int k = 0; k < K; k++) {
    cells_per_cat[k].assign(st.cat_counts[k], {});
    for (int c = 0; c < M; c++) {
        int g = ct.g_per_cell[k][c];
        if (g >= 0 && g < st.cat_counts[k]) cells_per_cat[k][g].push_back(c);
    }
}
```

**Target file:** `src/cell_table.hpp` — this is a structural property of `CellTable`, not a calibration dispatch concern.

**Signature:**

```cpp
// cell_table.hpp — namespace lbw (after CellTable struct)

/// Build cells_per_cat[k][j] = sorted list of cell indices c where g_per_cell[k][c] == j.
/// K and cat_counts must match the CellTable that was built from the same inputs.
inline std::vector<std::vector<std::vector<int>>>
build_cells_per_cat(const CellTable& ct, int K, const int* cat_counts)
{
    std::vector<std::vector<std::vector<int>>> cpc(K);
    for (int k = 0; k < K; k++) {
        cpc[k].assign(cat_counts[k], {});
        for (int c = 0; c < ct.M_cell; c++) {
            int g = ct.g_per_cell[k][c];
            if (g >= 0 && g < cat_counts[k]) cpc[k][g].push_back(c);
        }
    }
    return cpc;
}
```

**Files modified:**
- `src/cell_table.hpp` — add `build_cells_per_cat` free function after the `estimate_M_cell` declaration.
- `src/greenkhorn.cpp` — replace lines 39–46 with `auto cells_per_cat = lbw::build_cells_per_cat(ct, K, st.cat_counts);`
- `src/raking.cpp` — replace lines 90–98 with `auto cells_per_cat = lbw::build_cells_per_cat(ct, st.K, st.cat_counts);`
- `src/logit_calib.cpp` — replace lines 53–60 with `auto cells_per_cat = lbw::build_cells_per_cat(ct, K, st.cat_counts);`

**Steps:**
- [ ] Verify callers already include `cell_table.hpp`: `grep -n '#include.*cell_table' src/greenkhorn.cpp src/raking.cpp src/logit_calib.cpp` — add `#include "cell_table.hpp"` to any file missing it before proceeding.
- [ ] Add `build_cells_per_cat` to `src/cell_table.hpp`.
- [ ] Replace `greenkhorn.cpp`; `R CMD INSTALL --preclean .`.
- [ ] Replace `raking.cpp`; compile.
- [ ] Replace `logit_calib.cpp`; compile.
- [ ] `Rscript -e 'testthat::test_dir("tests/testthat")'` — green.
- [ ] `git commit -m "refactor(cell_table): extract build_cells_per_cat() helper — fcbo.4"`

---

## Task 5 — fcbo.5: Hoist `compute_weight_kl` to `calib_dispatch.hpp`

**Ticket:** fcbo.5

**Three sites:** Each solver defines a local lambda with the same semantic (weight-space KL divergence Σ X[c]·log(X[c]/X_init[c])/n), but with two implementations:

*Simple variant* — `ieppa.cpp:337–345`, `sinkhorn.cpp:141–149`:
```cpp
auto compute_weight_kl = [&]() -> double {
    double wkl = 0.0;
    const double inv_n = 1.0 / static_cast<double>(st.n);
    for (int c = 0; c < ct.M_cell; c++) {
        if (X_init[c] > 0.0 && X[c] > 0.0)
            wkl += X[c] * std::log(X[c] / X_init[c]) * inv_n;
    }
    return std::isfinite(wkl) ? wkl : 0.0;
};
```

*Vectorized (bulk_log) variant* — `raking.cpp:175–190`:
```cpp
auto compute_weight_kl = [&]() -> double {
    const double inv_n = 1.0 / static_cast<double>(st.n);
    int valid_count = 0;
    for (int c = 0; c < ct.M_cell; c++) {
        if (X_init[c] > 0.0 && X[c] > 0.0) {
            kl_ratio_scratch[valid_count] = X[c] / X_init[c];
            kl_weight_scratch[valid_count] = X[c];
            valid_count++;
        }
    }
    lbw::bulk_log(kl_ratio_scratch.data(), kl_ratio_scratch.data(), valid_count);
    double wkl = 0.0;
    for (int i = 0; i < valid_count; i++)
        wkl += kl_weight_scratch[i] * kl_ratio_scratch[i] * inv_n;
    return std::isfinite(wkl) ? wkl : 0.0;
};
```

**Canonical:** Use the raking `bulk_log` version — it is vectorization-friendly and more numerically stable for large M_cell. The scratch vectors (`ratio_buf`, `weight_buf`) are passed by the caller to avoid heap allocation per call.

**Signature:**

```cpp
// calib_dispatch.hpp — add after compute_cell_metrics
// Requires: #include "lbw_math.hpp" (for bulk_log)

/// Weight-space KL: Σ_c X[c]*log(X[c]/X_init[c]) / n.
/// ratio_buf and weight_buf are caller-owned scratch of size >= ct.M_cell.
/// Returns 0.0 when result is non-finite.
inline double compute_weight_kl(
    const std::vector<double>& X,
    const std::vector<double>& X_init,
    int M_cell, int n,
    double* ratio_buf,
    double* weight_buf) noexcept
{
    const double inv_n = 1.0 / static_cast<double>(n);
    int valid = 0;
    for (int c = 0; c < M_cell; c++) {
        if (X_init[c] > 0.0 && X[c] > 0.0) {
            ratio_buf[valid]  = X[c] / X_init[c];
            weight_buf[valid] = X[c];
            valid++;
        }
    }
    lbw::bulk_log(ratio_buf, ratio_buf, valid);
    double wkl = 0.0;
    for (int i = 0; i < valid; i++) wkl += weight_buf[i] * ratio_buf[i] * inv_n;
    return std::isfinite(wkl) ? wkl : 0.0;
}
```

**Files modified:**
- `src/calib_dispatch.hpp` — add `compute_weight_kl` free function; add `#include "lbw_math.hpp"` if not already present.
- `src/raking.cpp` — remove local lambda `compute_weight_kl` (lines 175–190); replace calls with `lbw::compute_weight_kl(X, X_init, ct.M_cell, st.n, kl_ratio_scratch.data(), kl_weight_scratch.data())`. The scratch vectors `kl_ratio_scratch` / `kl_weight_scratch` are already declared; keep them.
- `src/sinkhorn.cpp` — remove local lambda (lines 141–149); add two scratch vectors `std::vector<double> kl_ratio_buf(ct.M_cell), kl_weight_buf(ct.M_cell);` before line 151 (the existing `bucket` declaration); replace lambda calls.
- `src/ieppa.cpp` — remove local lambda (lines 337–345); add two scratch vectors `std::vector<double> kl_ratio_buf(ct.M_cell), kl_weight_buf(ct.M_cell);` before line 347 (the homotopy outer-driver comment); replace lambda calls.

**Steps:**
- [ ] Verify `lbw_math.hpp` is already `#include`d in `calib_dispatch.hpp` or add it.
- [ ] Add `compute_weight_kl` to `src/calib_dispatch.hpp`.
- [ ] Update `src/raking.cpp`; `R CMD INSTALL --preclean .`.
- [ ] Update `src/sinkhorn.cpp`; compile.
- [ ] Update `src/ieppa.cpp`; compile.
- [ ] `Rscript -e 'testthat::test_dir("tests/testthat")'` — green.
- [ ] `git commit -m "refactor(solvers): hoist compute_weight_kl() to calib_dispatch.hpp — fcbo.5"`

---

## Task 6 — fcbo.6: Extract `lbw::aggregate_to_margin()` helper

**Ticket:** fcbo.6

**Pattern:** The "accumulate cell masses into a margin bucket" loop appears 7+ times across all solvers:

```cpp
// Generic form — seen in compute_cell_metrics, sinkhorn, raking, greg, chebyshev, greenkhorn, logit_calib
std::fill(bucket.begin(), bucket.begin() + nj, 0.0);
for (int c = 0; c < ct.M_cell; c++) {
    int g = ct.g_per_cell[k][c];
    if (g >= 0 && g < nj) bucket[g] += X[c];
}
```

**Signature:**

```cpp
namespace lbw {

/// Aggregate cell masses X[] into margin k bucket[0..cat_counts[k]).
/// bucket must be pre-allocated to at least cat_counts[k] elements.
/// Zeroes bucket[0..nj) before accumulating.
inline void aggregate_to_margin(
    const CellTable& ct,
    const std::vector<double>& X,
    int k, int nj,
    double* bucket) noexcept
{
    std::fill(bucket, bucket + nj, 0.0);
    for (int c = 0; c < ct.M_cell; c++) {
        int g = ct.g_per_cell[k][c];
        if (g >= 0 && g < nj) bucket[g] += X[c];
    }
}

} // namespace lbw
```

**Files modified (replace sites):**
- `src/calib_dispatch.hpp` — add `aggregate_to_margin` (and update `compute_cell_metrics` internally to call it).
- `src/sinkhorn.cpp` — replace the inner-loop bucket fill at lines 163–168 for each margin k.
- `src/raking.cpp` — the `compute_errRp_ct` static helper at lines 42–46; and the convergence-check bucket loop.
- `src/greenkhorn.cpp` — bucket sums at lines 80–81 and recompute sites.
- `src/logit_calib.cpp` — margin defect accumulation loop.
- `src/greg.cpp` — lines 90–95 (`b[]` build).

**Steps:**
- [ ] Add `aggregate_to_margin` to `src/calib_dispatch.hpp`.
- [ ] Update `compute_cell_metrics` in the same file to call `aggregate_to_margin`.
- [ ] Replace sites in `src/sinkhorn.cpp`; `R CMD INSTALL --preclean .`.
- [ ] Replace sites in `src/raking.cpp`; compile.
- [ ] Replace sites in `src/greg.cpp`; compile.
- [ ] Replace sites in `src/greenkhorn.cpp`; compile.
- [ ] Replace sites in `src/logit_calib.cpp`; compile.
- [ ] `Rscript -e 'testthat::test_dir("tests/testthat")'` — green.
- [ ] `git commit -m "refactor(solvers): extract aggregate_to_margin() helper — fcbo.6"`

---

## Task 7 — fcbo.7: Extract `lbw::mark_converged()` helper

**Ticket:** fcbo.7

**Pattern:** A 4–6 line block that populates the convergence metadata fields of any result struct. Repeated verbatim in every solver that uses `check_convergence`:

```cpp
// sinkhorn.cpp:245-251 (representative)
res.status             = RK_OK;
res.convergence_metric = static_cast<int>(st.convergence_cfg.metric);
res.convergence_rule   = static_cast<int>(st.convergence_cfg.rule);
res.convergence_tol    = st.convergence_cfg.pct_tol;
res.convergence_iter   = iter;
```

**Sites:**
| File | Lines (approx) |
|------|---------------|
| `src/sinkhorn.cpp` | 245–251 |
| `src/greenkhorn.cpp` | 234–239 (sets `res.status = RK_OK; res.convergence_iter`) |
| `src/raking.cpp` | 489–494 |
| `src/chebyshev.cpp` | 284–290 |
| `src/logit_calib.cpp` | 275–278 (partial — sets `res.status = RK_OK; res.convergence_iter`) |

**Note:** `greenkhorn` and `logit_calib` currently write only a subset of the 5 fields, but both `GreenkornResult` and `LogitCalibResult` have been confirmed to carry all 5 fields (`convergence_metric`, `convergence_rule`, `convergence_tol`, `convergence_iter`, `status`). Apply `mark_converged` unconditionally to all 5 solvers — no conditional check or "leave as-is" fallback needed.

**Signature:**

```cpp
namespace lbw {

/// Populate standard convergence metadata on a result struct.
/// ResT must have: .status (int), .convergence_metric (int), .convergence_rule (int),
///                 .convergence_tol (double), .convergence_iter (int).
template <typename ResT>
inline void mark_converged(ResT& res, const CalibConvergenceCfg& cfg, int iter) noexcept {
    res.status             = RK_OK;
    res.convergence_metric = static_cast<int>(cfg.metric);
    res.convergence_rule   = static_cast<int>(cfg.rule);
    res.convergence_tol    = cfg.pct_tol;
    res.convergence_iter   = iter;
}

} // namespace lbw
```

**Files modified:**
- `src/calib_dispatch.hpp` — add `mark_converged` template.
- `src/sinkhorn.cpp` — replace convergence assignment block.
- `src/raking.cpp` — replace convergence assignment block.
- `src/chebyshev.cpp` — replace convergence assignment block (note: chebyshev also sets `res.convergence_rule` and `res.convergence_tol` — they are covered by the template).
- `src/greenkhorn.cpp` — replace lines 234–239 with `lbw::mark_converged(res, st.convergence_cfg, iter)` (GreenkornResult confirmed to have all 5 fields).
- `src/logit_calib.cpp` — replace lines 275–278 with `lbw::mark_converged(res, st.convergence_cfg, iter)` (LogitCalibResult confirmed to have all 5 fields).

**Steps:**
- [ ] Add `mark_converged` to `src/calib_dispatch.hpp`.
- [ ] Replace block in `src/sinkhorn.cpp`; `R CMD INSTALL --preclean .`.
- [ ] Replace block in `src/raking.cpp`; compile.
- [ ] Replace block in `src/chebyshev.cpp`; compile.
- [ ] Replace block in `src/greenkhorn.cpp`; compile.
- [ ] Replace block in `src/logit_calib.cpp`; compile.
- [ ] `Rscript -e 'testthat::test_dir("tests/testthat")'` — green.
- [ ] `git commit -m "refactor(solvers): extract mark_converged() template helper — fcbo.7"`

---

## Cross-cutting: TDD protocol for all tasks

Because these are internal C++ helpers with no R-visible API surface, the "test" is the existing R-level testthat suite. For each task:

1. **Baseline:** Before any edits, confirm all tests pass:
   ```bash
   Rscript -e 'testthat::test_dir("tests/testthat")'
   ```
2. **Implement** the helper in the header.
3. **Replace one call site** and compile:
   ```bash
   R CMD INSTALL --preclean .
   ```
4. **Verify** R output is bit-identical for the affected solver using the relevant test file:
   - fcbo.1–4: `test-calibration-solvers.R`, `test-raking.R`, `test-logit.R`
   - fcbo.5: `test-raking.R`, `test-ieppa.R`, `test-calibration-solvers.R`
   - fcbo.6: all of the above
   - fcbo.7: `test-calibration-solvers.R`, `test-raking.R`
5. **Repeat** for remaining call sites.
6. **Full suite** must be green before commit.

**Compile command:** `cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean .`

**Test command:** `Rscript -e 'testthat::test_dir("tests/testthat")'`

---

## Execution order (recommended)

Tasks 2–4 are independent and can be executed in any order before or after Task 1. Task 5 depends on `lbw_math.hpp` being present (it is — already used in raking). Tasks 6 and 7 are best done last as they touch the most sites.

Recommended: 2 → 3 → 4 → 1 → 5 → 6 → 7

- Tasks 2–3 are the smallest (1–3 line patterns) — lowest risk, fastest validation.
- Task 4 (`build_cells_per_cat`) adds to `cell_table.hpp` with no deps on dispatch.
- Task 1 supersedes parts of 2–3 in the preamble; do it after to avoid double-replacing lines.
- Tasks 5–7 modify algorithm loop bodies — do last with fresh eyes.

---

## Deferred Findings

The following findings from the 2026-04-30 Opus review are explicitly out-of-scope for this reuse epic and deferred to separate tickets. They must not be addressed incidentally during Tasks 1–7.

### F11 — Duplicate target-sum validator in validation.hpp + calib_validate.cpp

`validation.hpp:74-86` and `calib_validate.cpp:78-88` both perform a target-sum-to-1 check. The check runs twice per R-bridge call with no guard. This is a correctness/maintenance risk but not a regression introduced by this epic.

**Ticket:** `leafblower-c8dg`

### F14 — `compute_errRp_ct` in raking.cpp potentially dead code

`raking.cpp:34-53` defines `compute_errRp_ct` (errRp only, O(K×M_cell)). The Opus review found it may have no callers in the main raking body (possibly replaced inline during an earlier refactor). Requires a `grep` audit before deletion.

**Ticket:** `leafblower-9jmj`
