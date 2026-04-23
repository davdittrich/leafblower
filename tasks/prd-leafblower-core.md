# Leafblower Core — Product Requirements Document

**Status:** Draft v3 (post design-review-gate round 2)  
**Author:** Dennis Alexis Valin Dittrich  
**Date:** 2026-04-18  
**Epic:** leafblower-kk1  
**Changes v3:** Resolved all 16 round-2 blockers + 3 round-3 blockers: error-code dedup, struct rename (inner_max_iter), normU defined (=max_weight, ‖U‖_∞ derivation), complexity formula (int64_t, init=0), KKT residuals simplified (errRp only), S_j threshold (1e-15×W), logit singularity guard ("either...or", independent checks), alloc overflow + zero-weight guard, H(u) corrected (L·u + (U-L)/A·ln(…)), anesrake full mapping table, Python copy contract + GIL trampoline spec, Windows Python wheel non-goal, TDD RED phase all 3 phases, thread safety, US-003 wall-clock, autumn (>=0.2.0) pinned.

---

## 1. Executive Summary

### Problem

Survey calibration (raking) in R is dominated by `autumn`, which uses iterative proportional fitting (IPF) — an O(n·K) per-iteration algorithm that degrades to non-convergence on large datasets with tight weight bounds, many margins, or sparse cells. No production-ready Python equivalent exists. Neither package exposes a lower-bound (`min_weight`) constraint.

### Solution

**Leafblower** is a high-performance survey calibration library with:
- A **C++17 core** exposing a **C API** (`leafblower.h`), shared by both language bindings
- An **R package** that is a drop-in for `autumn::harvest()` (same API shape, adds `min_weight`)
- A **Python package** with a pandas-first interface over the same C core
- Two algorithms — **iEPPA** and **L-BFGS-B** — selected automatically based on problem size and bound tightness

### Success Criteria

| Criterion | Target |
|-----------|--------|
| R API compatibility | All `autumn` exported functions present with same signatures |
| Performance — medium | 100K rows, 5 margins (3–5 cats each) < 1 s |
| Performance — large | 1M rows, 20 margins (5 cats each) < 30 s |
| Weight bound enforcement | `max(w) ≤ max_weight` and `min(w) ≥ min_weight` always, within 1e-10 |
| Convergence | Calibration error `max_k max_j |Σw·1[g_k=j]/Σw − τ_j^(k)| < 1e-6` at reported convergence |
| Distribution | CRAN submission-ready R package; PyPI-publishable Python wheel |

### Timeline

Phased delivery — see §11.

---

## 2. Goals

- **G-1:** Replace `autumn::harvest()` in existing R workflows with a zero-code-change drop-in (excluding `min_weight` addition)
- **G-2:** Handle 1M-row, 20-margin datasets within 30 s on a modern 8-core CPU
- **G-3:** Enforce hard bounds `[min_weight, max_weight]` on calibrated weights in both algorithms
- **G-4:** Provide a Python API for the same calibration workflow, accepting pandas DataFrames
- **G-5:** Share a single compiled C++ core (no algorithm code duplicated between R and Python)
- **G-6:** Ship CRAN-ready R package and PyPI-ready Python wheel in v1

---

## 3. User Stories

### US-001: R Drop-In Replacement

**Description:** As an R survey analyst using `autumn`, I want to replace `library(autumn)` with `library(leafblower)` and call `harvest(data, target)` identically so that I get calibrated weights without changing my script.

**Acceptance Criteria:**
- [ ] `harvest(data, target)` returns a data frame with a `weights` column using the same default parameters as `autumn::harvest()`
- [ ] All `autumn`-exported functions present: `harvest`, `anesrake`, `diagnose_weights`, `design_effect`, `effective_sample_size`, `get_current_miss`, `weighted_pct`
- [ ] The following `autumn` parameters are **accepted and honoured**: `attach_weights`, `weight_column`, `max_weight`, `max_iterations`, `verbose`, `start_weights`, `target_map`
- [ ] The following `autumn` parameters are **accepted and silently ignored** (with `verbose ≥ 2` note): `select_params`, `select_function`, `error_function`, `adaptive_order`, `enforce_mean`, `accelerate` — iEPPA/L-BFGS-B calibrate all variables simultaneously and cannot honour per-variable selection logic
- [ ] The following `autumn` parameters are **accepted with a deprecation warning**: `convergence["pct"]` (use `convergence["absolute"]` instead; `pct` criterion not applicable to iEPPA/L-BFGS-B); `auto_collapse`, `collapse_vars`, `add_na_proportion` — not in v1, raises informative error
- [ ] `method = "rake"` maps to `"lbfgsb"` with warning: `"method='rake' (IPF) not implemented; using L-BFGS-B (method='lbfgsb')"`
- [ ] `method = "nr"` maps to `"lbfgsb"` with warning: `"method='nr' (Newton-Raphson) not implemented; using L-BFGS-B (method='lbfgsb')"`
- [ ] Calibrated weights satisfy `max(w) ≤ max_weight` within 1e-10
- [ ] Weights normalized to mean ≈ 1 on output
- [ ] `R CMD check` passes with 0 errors, 0 warnings

**Priority:** High  
**Dependencies:** US-004 (C API), US-005 (iEPPA), US-006 (L-BFGS-B)  
**Docs impact:** `README.md` must document drop-in status and parameter differences  
**Config impact:** `DESCRIPTION Suggests: autumn` (needed for integration tests)

---

### US-002: Lower Bound on Weights

**Description:** As a survey methodologist, I want to specify `min_weight` so that calibrated weights never fall below a floor (preventing near-zero influential observations from being effectively excluded).

**Acceptance Criteria:**
- [ ] `harvest(data, target, min_weight = 0.5)` produces weights where `min(w) ≥ 0.5` within 1e-10
- [ ] Default `min_weight = 0` means no lower bound (only `max_weight` enforced); internally treated as clamping to max(0, w_i) — no special-casing needed since weights start positive
- [ ] `min_weight` accepted by `rk_calibrate()` C API (as `rk_params_t.min_weight`) and by Python `harvest()`
- [ ] **`min_weight ≥ max_weight` returns `RK_ERR_BADARG`** immediately with message: `"min_weight must be strictly less than max_weight"`
- [ ] `diagnose_weights()` output unchanged by `min_weight`
- [ ] Routing: `min_weight > 0` forces iEPPA (see US-007)

**Priority:** High  
**Dependencies:** US-004 (C API)  
**Docs impact:** `man/harvest.Rd` and Python docstring  
**Config impact:** None

---

### US-003: Large-Scale Calibration

**Description:** As a census microsimulation researcher, I want to calibrate 1M+ observations across 20+ margins in under 30 seconds so that I can iterate quickly on synthetic population models.

**Acceptance Criteria:**
- [ ] `harvest(data, target)` with n=1,000,000, K=20, 5 categories each completes in < 30 s wall-clock time (single-threaded)
- [ ] Algorithm automatically selects iEPPA or L-BFGS-B per US-007 routing rules
- [ ] `verbose = 1` prints selected algorithm name and routing reason
- [ ] Convergence criterion met (`max calibration error < 1e-6`) or informative warning issued

**Priority:** High  
**Dependencies:** US-005 (iEPPA), US-006 (L-BFGS-B), US-007 (routing)  
**Docs impact:** `vignette/performance.Rmd` (new) documents routing thresholds  
**Config impact:** None

---

### US-004: C API

**Description:** As a library integrator, I want a stable C API header (`leafblower.h`) so that I can call the calibration engine from any language with C FFI support without duplicating algorithm code.

**Acceptance Criteria:**
- [ ] `leafblower.h` is a valid C99 header (compiles with `gcc -std=c99 -Werror`)
- [ ] `rk_calibrate()` modifies weight array in-place and returns a status code
- [ ] `rk_params_init()` fills a `rk_params_t` struct with safe defaults
- [ ] `rk_result_t` includes `algorithm_used` field so callers know which algorithm ran
- [ ] Return codes: `RK_OK=0`, `RK_ERR_NOCONV=1`, `RK_ERR_INFEAS=2`, `RK_ERR_BADARG=3`
- [ ] `RK_ERR_INFEAS` triggered when: feasible set is provably empty (e.g., `min_weight ≥ max_weight`, or all observations in a category are absent and target > 0)
- [ ] `RK_ERR_BADARG` triggered when: NULL pointer for required params, `n ≤ 0`, `K ≤ 0`, any `cat_counts[k] ≤ 0`, NaN/Inf in `targets[]` or initial `weights[]`, `targets[k]` does not sum to 1 ± 1e-8
- [ ] Each `group_ids[k]` is a pointer to `n` contiguous `int32` elements; each `targets[k]` is a pointer to `cat_counts[k]` contiguous `double` elements; each `group_ids[k][i]` ∈ `{-1, 0, ..., cat_counts[k]-1}` where -1 = NA/OOV; values outside this range return `RK_ERR_BADARG`
- [ ] No C++ symbols leak into the header (all wrapped in `extern "C"`)
- [ ] R `.Call()` bridge and Python pybind11 both call the same `rk_calibrate()` symbol

**Priority:** High  
**Dependencies:** None  
**Docs impact:** `src/leafblower.h` is self-documenting via doxygen comments  
**Config impact:** None

---

### US-005: Faithful iEPPA Algorithm (Paper-Faithful algBCD at C=0)

**Description:** As a developer, I want leafblower to implement the paper-faithful iEPPA algorithm (Chu-Liang-Toh-Yang 2022, arXiv:2011.14312) at C=0 so that capacity-constrained multi-marginal calibration converges with cell-compressed O(M_cell·K) inner cost and a documented convergence framework (Csiszár 1975 cyclic I-projection).

**Acceptance Criteria:**
- [ ] `rk_calibrate(..., algorithm=RK_ALG_IEPPA)` converges on: 1M rows, 20 margins, max_weight=3 within 30 s
- [ ] Cell compression: unique (g_1, ..., g_K) tuples deduplicated via sort-based cell table; M_cell ≤ min(n, ∏ cat_counts); single cell-level weight expansion per outer iteration (O(M_cell·K) core, O(n) expansion)
- [ ] Log-space Sinkhorn factors lf[k][j] per margin-category; X_tilde[c] = X_init[c] · exp(Σ_k lf[k][g_k(c)])
- [ ] Log-sum-exp stabilization: partial sums clipped to [−∞, 700] to avoid overflow on exp() (IEEE double overflow at exp(709))
- [ ] Capacity BCD block: W[c] per cell, updated each inner iteration; X[c] = clamp(X_tilde[c] · W[c], L_c, U_c) where L_c = min_weight · n_per_cell[c], U_c = max_weight · n_per_cell[c]
- [ ] Convergence check: errRp = max_k max_j |S_kj / W_total − τ_kj| < tol_abs (default 1e-6)
- [ ] Verbose logging: verbose=1 entry (compression ratio), per-iter errRp, exit status; verbose=2 adds n_cap_active, log10(f[k][j]) range per margin for ill-conditioning debug; [AUTO→iEPPA] prefix when routed via AUTO
- [ ] Overflow detection: if log-factor drift max_log_X_tilde > 700 AND uncapped cell detected, return RK_ERR_NOCONV with actionable message (direct to raking or looser bounds)
- [ ] Returns `RK_ERR_NOCONV` (not crash) if `inner_max_iter` exhausted; output weights reflect last iterate; RK_ERR_INFEAS if empty cell with positive target
- [ ] All weights satisfy `[min_weight, max_weight]` at every inner iteration, not only at convergence

**Priority:** High  
**Dependencies:** US-004 (C API)  
**Docs impact:** `man/harvest.Rd` — add method="ieppa" (faithful algBCD); new §3 of `docs/iEPPA/` README  
**Config impact:** None

---

### US-005b: Classical Raking Algorithm (IPF + Dykstra)

**Description:** As a developer, I want `method="raking"` to provide the classical cyclic IPF (Deming-Stephan 1940; Csiszár 1975) with additive Dykstra box-projection (Boyle-Dykstra 1986) and hyperplane projection, as an alternative to the paper-faithful iEPPA.

**Acceptance Criteria:**
- [ ] `rk_calibrate(..., algorithm=RK_ALG_RAKING)` present in enum
- [ ] `method="raking"` in R/Python routes to `raking_solve`
- [ ] All pre-rev2 iEPPA tests pass against `method="raking"` (regression guard in `test-raking.R`)

**Priority:** High  
**Dependencies:** US-004 (C API)  
**Docs impact:** `man/harvest.Rd` — add method="raking"  
**Config impact:** None

---

### US-006: L-BFGS-B Algorithm

**Description:** As a developer, I want leafblower to implement L-BFGS-B over the Deville-Sarndal logit dual so that standard-scale problems converge faster with provably bounded weights.

**Acceptance Criteria:**
- [ ] `rk_calibrate(..., algorithm=RK_ALG_LBFGSB)` converges on: 100K rows, 5 margins within 1 s
- [ ] For `min_weight = 0` (default): use exponential link `F(u) = exp(u)` (logit limit as L→0); this avoids the undefined `A` parameter; weights naturally bounded above by `max_weight` via gradient projection, no lower bound enforced
- [ ] For `min_weight > 0` and `max_weight < Inf`: use logit link `F(u) = [L(U-1) + U(1-L)·exp(A·u)] / [(U-1) + (1-L)·exp(A·u)]` with `A = (U-L)/((U-1)(1-L))`, `L = min_weight`, `U = max_weight` (weights normalized to mean 1 before entry, so L and U are in same units)
- [ ] For `max_weight = Inf` (unbounded): use exponential link `F(u) = exp(u)`; this special-case must be checked before computing `A` to avoid `∞` in the formula
- [ ] `exp(A·u)` is clamped to `exp(700)` before use in all formulas (overflow guard; IEEE 754 double overflows at exp(709))
- [ ] Dual objective `φ(λ) = Σ_{k,j} T_kj · λ_kj − Σ_i d_i · H(u_i(λ))`, where `H(u) = ∫F(u)du`:
  - Exponential link: `H(u) = exp(u)` (antiderivative of exp is exp)
  - Logit link: `H(u) = L·u + (U-L)/A · ln(((U-1) + (1-L)·exp(Au)) / (U-L))` (analytically derived: H'(u) = F(u) ✓; H(0) = (U-L)/A · ln((U-L)/(U-L)) = 0 ✓)
- [ ] Gradient `∂φ/∂λ_kj = T_kj − Σ_{i: g_k(i)=j} w_i(λ)` where `w_i = d_i · F(Σ_k λ_{k,g_k(i)})`
- [ ] L-BFGS 2-loop recursion: history of last `m = lbfgs_m` (default 10) `(s,y)` pairs; Hessian scaling by `(s·y)/(y·y)` 
- [ ] Wolfe line search: sufficient decrease `c1=1e-4`, curvature `c2=0.9`
- [ ] Final weights: `w_i = clamp(d_i · F(u_i), min_weight, max_weight)` — the clamp handles any minor numerical violations from the gradient-projection dual solution
- [ ] Convergence: gradient norm `‖∇φ‖_∞ < tol_abs` or `max_k max_j |T_kj/Σw − τ_j^(k)| < tol_abs` — these are equivalent at convergence since `∂φ/∂λ_kj = T_kj − S_kj`, so `‖∇φ‖_∞ = max calibration error`; L∞ is canonical throughout

**Priority:** High  
**Dependencies:** US-004  
**Docs impact:** None  
**Config impact:** None

---

### US-007: Automatic Algorithm Selection

**Description:** As a user, I want `harvest()` to automatically pick the right algorithm so that I don't need to know the difference between iEPPA and L-BFGS-B.

**Acceptance Criteria:**
- [ ] Default `method = "auto"` in `harvest()` and `algorithm = RK_ALG_AUTO` in C API
- [ ] Routing rule (implemented in `c_api.cpp`, not R/Python layer):
  ```
  complexity = n × sum(cat_counts[k] for k in 1..K)
  use_ieppa = (complexity > 500000) OR (max_weight < 3.0) OR (min_weight > 0.0)
  ```
  Note: the 500K threshold is empirically chosen; `verbose ≥ 1` must report the actual `complexity` value alongside the routing decision so users can understand and override it.
- [ ] `method = "ieppa"` and `method = "lbfgsb"` force specific algorithm; `method = "rake"` and `method = "nr"` map to `"lbfgsb"` with warning (see US-001)
- [ ] `rk_result_t.algorithm_used` is set to the algorithm that actually ran (`RK_ALG_IEPPA` or `RK_ALG_LBFGSB`), not `RK_ALG_AUTO`
- [ ] R bridge reads `result.algorithm_used` and includes it in `attr(data, "algorithm")` on the returned data frame
- [ ] `verbose ≥ 1` prints: `"Auto-selected iEPPA: complexity=1500000 > 500000 threshold"` or `"Auto-selected L-BFGS-B: complexity=200000 ≤ threshold, max_weight=5 ≥ 3, min_weight=0"`

**Priority:** Medium  
**Dependencies:** US-005, US-006  
**Docs impact:** `man/harvest.Rd` documents `method` param and routing rules  
**Config impact:** None

---

### US-008: Python Package — pandas Interface

**Description:** As a Python survey analyst, I want to call `leafblower.harvest(df, targets)` with a pandas DataFrame so that I can calibrate weights in a Pythonic workflow.

**Acceptance Criteria:**
- [ ] `harvest(data, targets)` accepts `pd.DataFrame` for `data` (auto-detected via `isinstance`) and a dict-of-dicts for `targets` (e.g., `{"age": {"18-34": 0.3, "35-54": 0.45, "55+": 0.25}}`)
- [ ] Also accepts `data` as `{col: list}` dict (converted to DataFrame internally)
- [ ] `targets` sum-to-1 validation performed in Python before calling C API (raises `ValueError` with column name if violated)
- [ ] `group_ids` length validation: `len(group_ids[k]) == n` for all k before calling `rk_calibrate()`
- [ ] Returns DataFrame with `weights` column appended (or `numpy.ndarray` dtype `float64` if `attach_weights=False`); `weights` array is always a **copy** — never a view into the input buffer; caller may mutate it without affecting internal state
- [ ] `min_weight`, `max_weight`, `method`, `verbose`, `max_iterations`, `convergence` parameters present; `convergence` is a dict accepting key `"absolute"` (maps to `tol_abs`); key `"pct"` raises `DeprecationWarning` and is ignored; any other unknown key raises `ValueError: unknown convergence key '{k}'`
- [ ] `diagnose_weights(data, targets, weights)` returns `pd.DataFrame` with same columns as R equivalent
- [ ] `pip install leafblower` installs a self-contained wheel (no separate C library install)
- [ ] Python ≥ 3.9 supported; CI tests 3.9, 3.10, 3.11, 3.12, 3.13
- [ ] pybind11 module named `_leafblower`; import path: `from leafblower._leafblower import calibrate`; `__init__.py` re-exports as public API `leafblower.harvest`

**Priority:** High  
**Dependencies:** US-004, US-005, US-006  
**Docs impact:** `python/README.md` (new)  
**Config impact:** None

---

### US-009: Diagnostic Functions

**Description:** As an R or Python user, I want `diagnose_weights()` and `design_effect()` to work exactly as in `autumn` so that I can evaluate calibration quality without changing my analysis code.

**Acceptance Criteria:**
- [ ] `diagnose_weights(data, target, weights)` returns data frame with columns: `variable`, `level`, `prop_original`, `prop_weighted`, `target`, `error_original`, `error_weighted`
- [ ] `design_effect(weights)` returns Kish (1992) estimator: `n × Σw² / (Σw)²`
- [ ] `design_effect(weights, outcome, data, target)` returns Henry & Valliant (2015) estimator
- [ ] `effective_sample_size(weights)` = `length(weights) / design_effect(weights)`
- [ ] `get_current_miss(data, target, weights)` exported with this exact name (matching autumn's export)
- [ ] `weighted_pct(x, weights)` returns named numeric vector of weighted proportions
- [ ] R implementations are pure R (no C calls required for diagnostics)

**Priority:** Medium  
**Dependencies:** US-001  
**Docs impact:** `man/diagnose_weights.Rd` etc.  
**Config impact:** `DESCRIPTION Suggests: autumn` (integration tests compare output)

---

### US-010: CRAN + PyPI Distribution

**Description:** As a package maintainer, I want leafblower to be distributable on CRAN and PyPI so that users can install it with `install.packages()` and `pip install`.

**Acceptance Criteria:**
- [ ] `R CMD check --as-cran` passes with 0 errors, 0 warnings, ≤ 1 note (new submission note acceptable)
- [ ] `configure` script detects `CXX17` support; generates `src/Makevars`; falls back to `CXX14` with a note if C++17 unavailable (no compilation failure)
- [ ] No vendored dependencies > 5 MB
- [ ] `.Rbuildignore` excludes `.beads/`, `.claude/`, `.wolf/`, `tasks/`, `python/`, `docs/iEPPA/`
- [ ] `pyproject.toml` with scikit-build-core produces a valid wheel via `pip wheel .`
- [ ] Wheel bundles compiled extension; `python -c "import leafblower; leafblower.harvest"` succeeds after `pip install`
- [ ] Python package version matches R package version (`DESCRIPTION` Version field)
- [ ] `cran-comments.md` documents any remaining notes for CRAN reviewers

**Priority:** Medium  
**Dependencies:** US-001, US-008  
**Docs impact:** `NEWS.md`, `cran-comments.md`  
**Config impact:** `DESCRIPTION` Imports/Suggests: `autumn (>= 0.2.0)` (pinned to audited API version); `pyproject.toml` build-system; `.Rbuildignore`

---

## 4. Functional Requirements

### Core Engine

| ID | Requirement |
|----|-------------|
| FR-1 | `rk_calibrate(n, K, weights, group_ids, cat_counts, targets, params, result)` — calibrates `weights[n]` in place; returns status code |
| FR-2 | `rk_params_init(p)` fills defaults: `min_weight=0.0`, `max_weight=5.0`, `algorithm=RK_ALG_AUTO`, `inner_max_iter=500`, `outer_max_iter=50`, `tol_abs=1e-6`, `epsilon=0.05`, `lbfgs_m=10`, `verbose=0`, `log_fn=NULL`. (R bridge maps autumn's `max_iterations` → `inner_max_iter`) |
| FR-3 | `group_ids[k]` = pointer to `n` contiguous `int32`; `group_ids[k][i]` ∈ `{-1, 0, …, cat_counts[k]-1}`; -1 = NA/OOV (observation skipped for margin k) |
| FR-4 | Input validation before any computation (all checks return `RK_ERR_BADARG`): NULL pointer for `weights`, `group_ids`, `cat_counts`, or `targets`; `n ≤ 0` or `K ≤ 0`; any `cat_counts[k] ≤ 0`; any `cat_counts[k] > n` (more categories than observations); any `group_ids[k][i] ≥ cat_counts[k]`; any `group_ids[k][i] < -1` (only -1 is valid NA); NaN or Inf in `targets[]` or initial `weights[]`; any `targets[k]` not sum to 1 ± 1e-8; `min_weight ≥ max_weight`; **logit singularity guard**: when `min_weight > 0 AND isfinite(max_weight)`, either `min_weight == 1.0` (makes `(1-L)=0` in A denominator) or `max_weight == 1.0` (makes `(U-1)=0` in A denominator) must be rejected with `RK_ERR_BADARG: "logit link undefined when min_weight=1 or max_weight=1"` — reject each independently, not together; `group_ids` range validation is a full O(n×K) pass completed **before** any weight modification; **allocation overflow guard**: compute `total_cats = Σ cat_counts[k]` as `size_t`; if `(size_t)n * total_cats > SIZE_MAX / 2` return `RK_ERR_BADARG` with message `"problem too large for platform size_t"`; **zero-weight guard**: `Σ weights[i] < 1e-15` (all-zero or driven-to-zero input) → `RK_ERR_BADARG: "total weight is zero or negative"` |
| FR-5 | `min_weight = 0` means no lower bound; clamp applied as `max(min_weight, w_i)` — when `min_weight = 0`, this is a no-op (weights start positive) |
| FR-6 | On `RK_ERR_NOCONV`, `weights[]` reflects last iterate; `rk_result_t.max_error` reflects last computed calibration error |
| FR-7 | `rk_result_t.algorithm_used` set to actual algorithm run (`RK_ALG_IEPPA` or `RK_ALG_LBFGSB`), never `RK_ALG_AUTO` |
| FR-8 | `rk_result_t.message[256]` always null-terminated even on truncation (use `snprintf`, never `sprintf`) |
| FR-9 | Verbose output via `log_fn` callback if non-NULL; falls back to `fprintf(stderr, …)` when `log_fn = NULL` and `verbose > 0`. R bridge passes R's `Rprintf`-backed callback. **Python bridge GIL contract**: `rk_calibrate()` is called from Python with the GIL held (no `py::gil_scoped_release`); `_bindings.cpp` registers a static C trampoline `void py_log_trampoline(const char* msg, void* ctx)` that casts `ctx` to `PyObject*` (a borrowed reference to a Python callable held alive for the duration of the call) and invokes it; because the GIL is held throughout, no `PyGILState_Ensure/Release` is needed; the Python `harvest()` wrapper passes `print` as the callable when `verbose > 0` |
| FR-10 | `RK_ERR_INFEAS` triggered by: a category has `S_j < 1e-15 × W` with `τ_j > 0` at the end of an inner BCD sweep (empty cell with positive target — calibration geometrically impossible). Note: `min_weight ≥ max_weight` is `RK_ERR_BADARG` (FR-4), not infeasibility |

### iEPPA Algorithm (FR-11 to FR-19)

Source: Chu, Liang, Toh & Yang (2022), arXiv:2011.14312. Reference implementation: `docs/iEPPA/code/ExpSynthesisData/solvers/ieppa_2d.m` and `bcd_2d.m`.

| ID | Requirement |
|----|-------------|
| FR-11 | **Initialization**: `w_i = d_i` (initial weights, normalized to mean 1); `ε = rk_params_t.epsilon` (fixed throughout, default 0.05); outer proximal center `w^0 = w` |
| FR-12 | **Outer loop**: runs for at most `outer_max_iter` (default 50) iterations; each outer iteration calls the inner BCD to approximately solve the entropic proximal subproblem |
| FR-13 | **Inner BCD — one sweep over K margins**: for each margin k = 1..K, for each category j = 1..c_k: (a) `S_j = Σ_{i: g_k(i)=j} w_i`; (b) `T_j = τ_j^(k) × W` where `W = Σ_i w_i`; (c) if `S_j < 1e-15 × W` and `T_j > 0`: set infeasibility flag; skip; do not use exact `== 0` (denormal sums are non-zero but near ε_machine); (d) `scale = T_j / S_j`; (e) `w_i ← w_i × scale` for all `i` where `g_k(i) = j`; (f) `w_i ← clamp(w_i, min_weight, max_weight)` |
| FR-14 | **Inner stopping criterion** (from paper §4, `bcd_2d.m` lines 32–47): `errRp = max_k max_j |S_j/W − τ_j^(k)|`; stop inner when `errRp < tolRp` AND Bregman distance `D(w̃, w^k) / (1 + normU) < tolRb`; where `D(w̃, w^k) = Σ_i [w̃_i log(w̃_i/w_i^k) − w̃_i + w_i^k]`; `normU = max_weight` — derivation: in `bcd_2d.m` the paper computes `normU = norm(U(:), inf)` where U is the capacity upper-bound matrix; for survey raking every weight has the same upper bound `max_weight`, so U is a scalar (or constant vector) and `‖U‖_∞ = max_weight`; this correctly scales the Bregman tolerance relative to the problem's weight magnitude; `tolRp` is adaptive: `1.0` on the first outer iteration, then `max(1e-6, errRp_prev / 1.5)` on subsequent iterations; `tolRb = 1 / outer_iter^1.1` |
| FR-15 | **Outer stopping criterion**: For categorical survey raking, `errRd`, `errRu`, `errRc` (dual infeasibility, capacity complementarity, KL complementarity) are identically 0 — no continuous transport plan, no vector capacity constraint, and KL complementarity is implicit in the Sinkhorn update. Outer stopping criterion simplifies to: `errRp < tol_abs` where `errRp = max_k max_j |S_j/W − τ_j^(k)|` |
| FR-16 | **NA/OOV rows**: `group_ids[k][i] = -1`; skip in `S_j` computation; do not update weight; weight remains unchanged for that margin pass |
| FR-17 | If outer loop exhausted without convergence: return `RK_ERR_NOCONV`; set `rk_result_t.message` to `"iEPPA: outer loop exhausted after N iterations, max_error=E"` |
| FR-18 | `rk_params_t.inner_max_iter` (default 500; R bridge maps autumn's `max_iterations` here) caps BCD sweeps per outer iteration regardless of stopping criterion |
| FR-19 | All weights satisfy `[min_weight, max_weight]` after every inner BCD step (invariant maintained throughout, not just at convergence) |

### L-BFGS-B Algorithm (FR-20 to FR-28)

| ID | Requirement |
|----|-------------|
| FR-20 | **Link function selection**: `min_weight = 0` OR `max_weight = Inf` → exponential link `F(u) = exp(u)` (check `!isfinite(max_weight)` first); otherwise logit link; this choice is made once at algorithm entry |
| FR-21 | **Exponential link**: `F(u) = exp(u)`, `F'(u) = exp(u)`, `H(u) = exp(u)` (antiderivative = itself); `exp(u)` clamped to `exp(700)` before use |
| FR-22 | **Logit link** (L > 0, U < ∞, and L ≠ 1 and U ≠ 1 — enforced by FR-4): `A = (U-L)/((U-1)(1-L))`; `F(u) = [L(U-1) + U(1-L)·exp(Au)] / [(U-1) + (1-L)·exp(Au)]`; `F'(u) = A·(F(u)-L)·(U-F(u))/(U-L)`; `H(u) = L·u + (U-L)/A · ln(((U-1) + (1-L)·exp(Au)) / (U-L))` — derived by integrating F(u); verified: H'(u) = L + (U-L)·(1-L)·exp(Au)/((U-1)+(1-L)·exp(Au)) = F(u) ✓; H(0) = 0 ✓ (since (U-1+1-L) = (U-L)); `exp(Au)` clamped to `exp(700)` before use |
| FR-23 | **Dual variables**: `λ ∈ R^{Σ c_k}` (one scalar per category per margin); `u_i = Σ_k λ_{k, g_k(i)}` (dual aggregate for obs i); `w_i = d_i · F(u_i)` |
| FR-24 | **Objective**: `φ(λ) = Σ_{k,j} T_kj · λ_kj − Σ_i d_i · H(u_i(λ))`; maximized (or equivalently, minimize `-φ`) |
| FR-25 | **Gradient**: `∂φ/∂λ_kj = T_kj − Σ_{i: g_k(i)=j} w_i(λ)` |
| FR-26 | **L-BFGS 2-loop recursion**: history of last `lbfgs_m` `(s,y)` pairs where `s^t = λ^{t+1} - λ^t`, `y^t = ∇φ(λ^{t+1}) - ∇φ(λ^t)`; initial Hessian scale `γ = (s·y)/(y·y)` |
| FR-27 | **Wolfe line search**: backtracking with sufficient decrease `c1=1e-4` and curvature `c2=0.9`; max 20 step halvings |
| FR-28 | **Final weights**: `w_i = clamp(d_i · F(u_i), min_weight, max_weight)`; convergence declared when `‖∇φ‖_∞ < tol_abs` or `max calibration error < tol_abs` |

### R Package (FR-29 to FR-35)

| ID | Requirement |
|----|-------------|
| FR-29 | `harvest()` parameters (beyond autumn compat): `min_weight=0`, `method="auto"` (also accepts `"ieppa"`, `"lbfgsb"`, `"rake"`, `"nr"`) |
| FR-30 | `harvest()` normalizes `start_weights` to mean 1 before passing to C API |
| FR-31 | Factor/character columns in `data` encoded to 0-indexed `int32` arrays by R bridge; NA → -1; OOV (level not in target) → -1 |
| FR-32 | Return value structure identical to `autumn::harvest()` output; additionally sets `attr(data, "algorithm")` to the algorithm name string |
| FR-33 | `useDynLib(leafblower, .registration=TRUE)` in NAMESPACE; routines registered via `R_registerRoutines()` in `zzz.R` |
| FR-34 | `anesrake()` is a wrapper around `harvest()` mapping anesrake parameter names. Full mapping: `inputter` → `data`; `weightvec` → `start_weights` (NULL → uniform); `caseid` → ignored (row label, not used by calibration); `targets` → `target`; `pctlim` → `convergence["pct"]` (deprecated warning issued); `cap` → `max_weight`; `choosemethod` ("rake"/"nrake") → `method` (both map to `"lbfgsb"` with warning); `type` → ignored; `nlim` → `max_iterations`; `iterate` → ignored (always iterates); `threads` → ignored (single-threaded v1). Unknown anesrake params → warning, ignored |
| FR-35 | `get_current_miss()` exported exactly (matching autumn's export name) |

### Python Package (FR-36 to FR-40)

| ID | Requirement |
|----|-------------|
| FR-36 | `harvest(data, targets, ...)` auto-detects `pd.DataFrame`; also accepts `{col: list}` dict (converted to DataFrame before processing) |
| FR-37 | Categorical columns encoded as 0-indexed `numpy.int32` arrays; NA → -1 |
| FR-38 | `python/CMakeLists.txt` explicitly excludes `r_bridge.cpp` from compiled sources (that file includes `Rinternals.h` which is unavailable in Python build environment) |
| FR-39 | pybind11 module `_leafblower` exposes `calibrate(n, K, weights_np, group_ids_list, cat_counts_list, targets_list, params_dict)` returning `(status, weights_out, result_dict)`; `weights_out` is a **new** `numpy.ndarray` (dtype=float64, C-contiguous copy) — never a view into `weights_np`; the binding enforces `weights_np.dtype == float64` and C-contiguous layout before passing the raw pointer to `rk_calibrate()`, raising `ValueError` otherwise |
| FR-40 | `diagnose_weights()` returns `pd.DataFrame`; `design_effect()` and `effective_sample_size()` implemented in pure Python |

---

## 5. Non-Goals (v1)

- **GPU acceleration** — CPU-only
- **Bounded IPF (water-filling) as a named `method`** — not exposed; `method="rake"` maps to L-BFGS-B with warning
- **SQUAREM acceleration** — not ported from autumn
- **`auto_collapse`** — sparse level merging not in v1 (raises informative error if requested)
- **`add_na_proportion`** — NA-as-category not in v1 (raises informative error if requested)
- **Variance estimation** — jackknife/bootstrap not in v1
- **Multi-threading** — single-threaded in v1
- **Windows CRAN binary** — source-only CRAN submission acceptable for v1; `Makevars.win` deferred
- **Windows Python wheel** — PyPI wheel for Windows deferred to v2; Windows users must build from source; CI does not include Windows targets in v1; this is an explicit non-goal (not an oversight)

---

## 6. Design Considerations

### C API Contract (`leafblower.h`)

```c
// Return codes
#define RK_OK         0  // Success
#define RK_ERR_NOCONV 1  // Did not converge within outer_max_iter
#define RK_ERR_INFEAS 2  // Infeasible: empty cell (S_j < 1e-15*W) with positive target
#define RK_ERR_BADARG 3  // Invalid argument: NULL pointer, bad dimensions, NaN/Inf, targets not sum to 1

typedef enum {
    RK_ALG_AUTO   = 0,
    RK_ALG_IEPPA  = 1,
    RK_ALG_LBFGSB = 2
} rk_algorithm_t;

typedef struct {
    double          min_weight;      // default 0.0 (no lower bound)
    double          max_weight;      // default 5.0
    int             inner_max_iter;  // inner BCD max sweeps per outer iter, default 500 (R bridge maps max_iterations here)
    int             outer_max_iter;  // outer EPP max iters (iEPPA) / L-BFGS max iters, default 50
    double          tol_abs;         // convergence tolerance, default 1e-6
    rk_algorithm_t  algorithm;       // default RK_ALG_AUTO
    int             verbose;         // 0=silent, 1=progress, 2=debug
    // iEPPA-specific
    double          epsilon;         // entropic proximal parameter, default 0.05
    // L-BFGS-B-specific
    int             lbfgs_m;         // history size, default 10
    // Logging callback (NULL = use fprintf(stderr) when verbose > 0)
    void            (*log_fn)(const char* msg, void* ctx);
    void*           log_ctx;         // passed as second arg to log_fn
} rk_params_t;

typedef struct {
    int             status;           // return code (RK_OK etc.)
    int             iterations;       // outer iterations completed
    double          max_error;        // max calibration error at last iterate
    rk_algorithm_t  algorithm_used;   // actual algorithm that ran (never RK_ALG_AUTO)
    char            message[256];     // null-terminated status message (snprintf-safe)
} rk_result_t;

// Fill *p with safe defaults
void rk_params_init(rk_params_t* p);

// Calibrate survey weights in-place
// weights[n]: on input = initial (design) weights; on output = calibrated weights
// group_ids[k]: pointer to n contiguous int32; value -1 = NA/OOV
// cat_counts[k]: number of categories in margin k
// targets[k]: pointer to cat_counts[k] contiguous doubles; must sum to 1 ± 1e-8
int rk_calibrate(
    int n, int K,
    double* weights,
    const int** group_ids,
    const int* cat_counts,
    const double** targets,
    const rk_params_t* params,   // NULL = use defaults
    rk_result_t* result          // NULL = ignore
);
```

### Algorithm Routing (Auto Mode, in `c_api.cpp`)

```c
static rk_algorithm_t select_algorithm(int n, int K,
                                        const int* cat_counts,
                                        const rk_params_t* p) {
    int64_t complexity = 0;  // int64_t not long — long is 32-bit on Windows
    for (int k = 0; k < K; k++) complexity += (int64_t)n * cat_counts[k];
    // complexity = n × Σ cat_counts[k]  (total category-observation pairs)
    // Use iEPPA when: large problem, tight upper bound, or lower bound present
    if (complexity > 500000L || p->max_weight < 3.0 || p->min_weight > 0.0)
        return RK_ALG_IEPPA;
    return RK_ALG_LBFGSB;
}
```

### Link Function Specification (L-BFGS-B only)

```
Case 1: min_weight == 0 OR !isfinite(max_weight)
  → Exponential link: F(u) = min(exp(u), exp(700))
    H(u) = exp(u)  [antiderivative = itself]

Case 2: min_weight > 0 AND isfinite(max_weight)  [L > 0, U < ∞]
  L = min_weight, U = max_weight
  A = (U - L) / ((U - 1) * (1 - L))
  e = min(exp(A*u), exp(700))   [clamped]
  F(u) = (L*(U-1) + U*(1-L)*e) / ((U-1) + (1-L)*e)
  F'(u) = A * (F(u) - L) * (U - F(u)) / (U - L)
  H(u) = L*u + (U-L)/A * ln(((U-1) + (1-L)*e) / (U-L))
          [H(0) = 0 by construction; H'(u) = F(u) verified analytically]
```

**Verification (case 2):** At `u=0`: `F(0) = (L(U-1) + U(1-L)) / ((U-1) + (1-L)) = (LU-L+U-UL) / (U-L) = (U-L)/(U-L) = 1` ✓ (weight = 1 when dual = 0, consistent with normalization). At `u→+∞`: `F(u) → U` ✓. At `u→−∞`: `F(u) → L` ✓.

### R Package File Map

```
leafblower/               (R package root)
├── DESCRIPTION           (NeedsCompilation: yes; SystemRequirements: C++17; Suggests: autumn)
├── NAMESPACE
├── configure             → generates src/Makevars from src/Makevars.in
├── R/
│   ├── harvest.R         → harvest() (all autumn params + min_weight + method)
│   ├── anesrake.R        → anesrake() compat wrapper
│   ├── diagnose_weights.R
│   ├── design_effect.R   → design_effect() + effective_sample_size()
│   ├── weighted_pct.R
│   ├── current_miss.R    → get_current_miss() (exported, matching autumn name)
│   └── zzz.R             → R_registerRoutines() + useDynLib
├── src/
│   ├── Makevars.in       → R build flags template
│   ├── leafblower.h      → PUBLIC C API (also included by Python pybind11)
│   ├── types.hpp         → internal structs (CalibState, MarginCache)
│   ├── margin.hpp/cpp    → encode categorical columns → group_id arrays + category sums
│   ├── logit.hpp/cpp     → link functions F(u), F'(u), H(u) for L-BFGS-B
│   ├── ieppa.hpp/cpp     → iEPPA outer EPP + inner BCD (Sinkhorn-style)
│   ├── lbfgsb_solver.hpp/cpp → L-BFGS-B 2-loop recursion + Wolfe line search
│   ├── c_api.cpp         → rk_calibrate(), rk_params_init(), select_algorithm()
│   └── r_bridge.cpp      → .Call() entry points (includes Rinternals.h; excluded from Python build)
├── tests/testthat/
│   ├── helper.R          → shared helpers (load autumn data if available)
│   ├── test-logit.R      → logit link: F(0)=1, F(u)∈[L,U], H'(u)=F(u)
│   ├── test-ieppa.R      → iEPPA convergence + bound enforcement
│   ├── test-lbfgsb.R     → L-BFGS-B convergence + bound enforcement
│   ├── test-harvest.R    → harvest() API + min_weight + method warnings
│   └── test-design.R     → design_effect() numerical verification
└── python/
    ├── CMakeLists.txt    → pybind11; sources: all src/*.cpp EXCEPT r_bridge.cpp
    ├── pyproject.toml    → scikit-build-core; version from DESCRIPTION
    └── leafblower/
        ├── __init__.py   → re-exports harvest(), diagnose_weights()
        ├── _harvest.py   → harvest() implementation wrapping _leafblower.calibrate
        └── _bindings.cpp → pybind11 bindings over C API
```

---

## 7. Technical Specifications

### Build System

**C++17 features used** (inventory for C++14 fallback validation):
- `if constexpr` — link function dispatch in `logit.hpp` template; C++14 fallback: `std::enable_if` or runtime branch
- Structured bindings (`auto [s, y] = ...`) — L-BFGS pair iteration; C++14 fallback: explicit `.first`/`.second`
- `std::optional<double>` — Wolfe search return value; C++14 fallback: `std::pair<bool, double>`
- `[[nodiscard]]` attribute — `rk_calibrate()` return code; C++14 fallback: omit attribute (no semantic change)
- All uses are in `logit.hpp`, `lbfgsb_solver.hpp`, `c_api.cpp`; `r_bridge.cpp` and `ieppa.cpp` use C++11-compatible code only
- C++14 fallback: `configure` sets `-std=c++14` flag; the four above features have `#if __cplusplus >= 201703L` guards with explicit C++14 equivalents in each file

**R Package:**
- `configure`: detects `CXX17` via `$CXX --std=c++17 -x c++ /dev/null -o /dev/null 2>/dev/null`; falls back to `CXX14` with printed note; generates `src/Makevars`
- `src/Makevars.in`: `PKG_CXXFLAGS = @CXXFLAGS_STD@ -I. -O3 -DSTRICT_R_HEADERS`; `PKG_SOURCES` lists all `.cpp` files except `r_bridge.cpp` treated specially (it includes R headers)
- `DESCRIPTION`: `Imports: Matrix`; `Suggests: autumn, testthat, bench`

**Python Package:**
- scikit-build-core + CMake
- `python/CMakeLists.txt`: adds `../src/*.cpp` EXCLUDING `../src/r_bridge.cpp`; finds pybind11 via `find_package(pybind11)`; pybind11 version ≥ 2.11 required
- Wheel bundles compiled extension; no shared library separate install needed
- `python/pyproject.toml`: `build-backend = "scikit_build_core.build"`; version read from `../DESCRIPTION`

### Memory Contract

- All arrays passed by pointer; caller owns all memory
- `rk_calibrate()` may allocate internal workspace (vectors of size n and Σc_k) using `new`/`delete`; all freed before return
- No heap allocation survives the call
- `rk_result_t` is caller-allocated stack struct; `rk_calibrate()` writes into it if non-NULL

### Thread Safety

`rk_calibrate()` is **reentrant** — it uses no global or static mutable state; all workspace is allocated on the heap within the call and freed before return. Concurrent calls from separate threads on separate data are safe. R packages using `parallel` or `future` may call `rk_calibrate()` concurrently; this is explicitly supported. The `log_fn` callback must remain valid for the entire duration of its `rk_calibrate()` call; R and Python bridges must not pass pointers to stack-allocated objects or GC-eligible objects that could be collected mid-call.

### Convergence Reporting

`rk_result_t.max_error` = `max_k max_j |Σ_{i:g_k(i)=j} w_i / Σ w_i − τ_j^(k)|` computed on final iterate.

---

## 8. Testing & Validation

### Unit Tests

| File | Test | Input | Pass Criterion |
|------|------|-------|----------------|
| `test-logit.R` | F(0) = 1 | L=0.5, U=5 | `abs(F(0) - 1) < 1e-12` |
| `test-logit.R` | F(u) ∈ [L,U] | 1000 random u | All values in [L, U] |
| `test-logit.R` | H'(u) = F(u) | numerical diff | `abs(H'(u) - F(u)) < 1e-8` |
| `test-logit.R` | exp link: F = exp | L=0, U=Inf | `F(u) == exp(u)` for u=0,1,-1 |
| `test-ieppa.R` | 1 margin, 2 cats, no bounds | n=100 | `max_error < 1e-6` |
| `test-ieppa.R` | 3 margins, tight bounds (max=2) | n=10000 | `max_error < 1e-6`, `max(w) ≤ 2` |
| `test-ieppa.R` | min_weight=0.5 | n=10000, 5 margins | `min(w) ≥ 0.5 - 1e-10` |
| `test-ieppa.R` | Large: 1M rows, 20 margins | max_weight=3 | `max_error < 1e-6`, time < 30s |
| `test-lbfgsb.R` | 3 margins, no bounds | n=50000 | `max_error < 1e-6`, time < 2s |
| `test-lbfgsb.R` | max_weight=5 (default) | n=100000 | `max(w) ≤ 5 + 1e-10` |
| `test-harvest.R` | Drop-in: `harvest(data, target)` | autumn's `ns_target` if available | Returns data frame with `weights` col, `mean(w) ≈ 1` |
| `test-harvest.R` | `min_weight=0.5` | any data | `min(result$weights) ≥ 0.5 - 1e-10` |
| `test-harvest.R` | `min_weight ≥ max_weight` | min=5, max=5 | Error with informative message |
| `test-harvest.R` | `method="rake"` warning | any | Warning mentioning "L-BFGS-B" |
| `test-harvest.R` | `method="nr"` warning | any | Warning mentioning "L-BFGS-B" |
| `test-harvest.R` | Bound property (random) | 50 random datasets | `max(w) ≤ max_weight + 1e-10` always |
| `test-design.R` | Kish estimator | `w = c(1,2,3,4)` | `n*sum(w^2)/sum(w)^2 ± 1e-12` |

### Integration Tests (require `autumn` in `Suggests`)

- Load `autumn::ns_target`; run `leafblower::harvest(ns_data, ns_target)` and `autumn::harvest(ns_data, ns_target)`; compare `diagnose_weights()` output — weighted errors must agree within 1e-4 (algorithms differ, so exact agreement not required)
- Python: same data as pandas DataFrame; `leafblower.harvest()` weights compared to R output — must agree within 1e-4

### Performance Benchmarks (CI artifacts, not gates)

```r
bench::mark(
  lbfgsb_100k = harvest(data_100k, target_5, method="lbfgsb"),
  ieppa_1m    = harvest(data_1m, target_20, method="ieppa"),
  iterations  = 3,
  check       = FALSE
)
```

Phase 2 gate: `ieppa_1m` median < 30 s.

---

## 9. Risks & Mitigation

| Risk | Impact | Probability | Mitigation |
|------|--------|-------------|------------|
| iEPPA non-convergence on tight bounds with many NA rows | High | Medium | Return `RK_ERR_NOCONV` with last iterate; `verbose=1` reports iteration count and error |
| `exp(A·u)` overflow in logit link (L-BFGS-B) | High | Medium | Clamp `exp(A·u)` to `exp(700)` before any division; verified in `test-logit.R` |
| CRAN rejection (C++17 not available on older builders) | Medium | Low | `configure` gracefully falls back to C++14 with a note |
| `snprintf` into `message[256]` truncation | Low | Low | All format strings verified to produce < 200 characters; `snprintf` truncates safely |
| Python wheel missing `r_bridge.cpp` exclusion | High | High (known) | FR-38 explicitly requires exclusion; Python CI fails fast if this is missing |
| `ρ` underflow to zero (N/A) | N/A | N/A | iEPPA does not use ρ — uses fixed ε; no underflow risk from parameter schedule |

---

## 10. Phased Rollout Plan

### Phase 1 — Core Engine + L-BFGS-B (Gate: L-BFGS-B works for standard cases)

**TDD RED phase** — write and verify these tests FAIL before writing implementation:
- `test-logit.R::test_F_at_zero` — fails: function not defined
- `test-logit.R::test_H_prime_equals_F` — fails: function not defined
- `test-lbfgsb.R::test_3margin_no_bounds` — fails: rk_calibrate not compiled
- `test-harvest.R::test_min_weight_badarg` — fails: harvest not exported
- `test-harvest.R::test_method_rake_warning` — fails: harvest not exported

**Implementation:**
- C API header + input validation + stub implementations
- Logit link functions (`logit.hpp/cpp`) with test-logit.R passing
- L-BFGS-B solver + `harvest()` with `method="lbfgsb"`
- `diagnose_weights()`, `design_effect()`, all autumn-compat functions
- R unit tests passing
- **Gate:** `harvest(ns_data, ns_target, method="lbfgsb")` produces `max_error < 1e-6`; all `test-lbfgsb.R` and `test-harvest.R` pass; `R CMD check` clean

### Phase 2 — iEPPA + Auto-Routing + min_weight (Gate: Large-scale benchmark passes)

**TDD RED phase** — write and verify these tests FAIL before writing iEPPA:
- `test-ieppa.R::test_1margin_2cat_no_bounds` — fails: RK_ALG_IEPPA not implemented
- `test-ieppa.R::test_3margin_tight_bounds` — fails: same
- `test-ieppa.R::test_min_weight_0.5` — fails: same
- `test-harvest.R::test_routing_auto_large` — fails: routing not wired

**Implementation:**
- iEPPA outer EPP + inner BCD (`ieppa.hpp/cpp`)
- Auto-routing logic in `c_api.cpp`
- `min_weight` parameter wired end-to-end
- Large-scale test (1M rows) passes
- **Gate:** `harvest(data_1m, target_20, method="ieppa")` median < 30 s; `min(w) ≥ min_weight − 1e-10` on all tests; `test-ieppa.R` all pass

### Phase 3 — Python Package (Gate: pip-installable wheel)

**TDD RED phase** — write and verify these tests FAIL before writing Python bindings:
- `test_python.py::test_harvest_returns_copy` — fails: `_leafblower` not importable
- `test_python.py::test_convergence_unknown_key_raises` — fails: `_leafblower` not importable
- `test_python.py::test_min_weight_badarg_python` — fails: `_leafblower` not importable

**Implementation:**
- pybind11 bindings (`_bindings.cpp`) — explicitly excludes `r_bridge.cpp`; GIL held during call; `py_log_trampoline` per FR-9
- `_harvest.py`, `diagnose_weights()` in Python
- scikit-build-core wheel; Python tests passing
- **Gate:** `pip install . && python -c "from leafblower import harvest; print('OK')"` succeeds on Linux and macOS CI

### Phase 4 — CRAN + PyPI (Gate: Submitted)

- `R CMD check --as-cran` clean; `cran-comments.md`
- `.Rbuildignore` excludes dev artifacts
- `NEWS.md`, vignette (`vignette/performance.Rmd`)
- `pyproject.toml` classifiers and metadata complete
- **Gate:** CRAN submission made (not necessarily accepted — "submitted and tracking" is the v1 gate per §11 KPI table)

---

## 11. Success Metrics & KPIs

| Metric | Target | Measurement |
|--------|--------|-------------|
| R API compat | All 7 autumn functions present | `R CMD check` + `test-harvest.R` |
| L-BFGS-B convergence | 100K rows, 5 margins < 2 s, `max_error < 1e-6` | `test-lbfgsb.R` Phase 1 gate |
| iEPPA convergence | 1M rows, 20 margins, max_weight=3 < 30 s | `test-ieppa.R` Phase 2 gate |
| Weight bound enforcement | `max(w) ≤ max_weight`, `min(w) ≥ min_weight` — 50 random datasets | property-based test in `test-harvest.R` |
| CRAN check | 0 errors, 0 warnings | `R CMD check --as-cran` |
| Python wheel | Installs on Linux/macOS, Python 3.9–3.13 | CI matrix |

---

## 12. Open Questions

| # | Question | Decision |
|---|----------|----------|
| OQ-1 | `method="rake"` / `"nr"` — silent map or warn? | **Resolved:** warn with informative message (see US-001 AC) |
| OQ-2 | `min_weight=0` logit degenerate case | **Resolved:** use exponential link when `min_weight=0` or `max_weight=Inf` (see FR-20, FR-21) |
| OQ-3 | iEPPA ε schedule — fixed or adaptive? | **Resolved:** fixed at `epsilon=0.05` per paper; outer stopping is adaptive via `tolRp` and `tolRb` |
| OQ-4 | verbose output: callback or `message()`/`print()`? | **Resolved:** `log_fn` callback in `rk_params_t` (NULL → `fprintf(stderr)`); R bridge passes R `message()` equivalent; Python bridge passes `print` (see FR-9) |
| OQ-5 | Windows CRAN binary | **Resolved as Non-Goal:** source-only for v1; `Makevars.win` deferred |
