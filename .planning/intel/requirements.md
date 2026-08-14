# Requirements (PRD-derived)

Source: **one** PRD-classified document — `tasks/prd-leafblower-core.md`
("Leafblower Core — Product Requirements Document", Status: Draft v3, dated 2026-04-18,
Epic `leafblower-kk1`). It is the ONLY PRD in the 44-doc ingest set and the OLDEST
document in it.

**Read this before using any entry below.** The PRD predates every SPEC in the corpus
(SPECs run 2026-04-18 … 2026-08-14) and describes a **two-algorithm** package; the
package now ships **eight solvers**. Precedence (`SPEC > PRD`) and document date agree,
so SPEC content wins wherever the two differ. Superseded PRD content is NOT extracted as
a requirement — it is listed in § Superseded / withdrawn at the foot of this file and
recorded as INFO (not as a conflict) in `../INGEST-CONFLICTS.md`.

Naming: entries keep the PRD's own `US-xxx` identifier inside the `REQ-` slug, because
the SPEC corpus cites those identifiers verbatim (`2026-04-23-ieppa-faithful-design.md`
§3.1/§4.3/§10 cite "PRD § US-001 AC", "§ US-005", "§ US-005b", "FR entries").

---

## REQ-us001-autumn-drop-in
- source: tasks/prd-leafblower-core.md § US-001 (also § G-1, FR-29 … FR-35)
- description: An R survey analyst replaces `library(autumn)` with
  `library(leafblower)` and calls `harvest(data, target)` identically, with no other
  script change (the sole addition being `min_weight`).
- acceptance:
  - `harvest(data, target)` returns a data frame with a `weights` column using the same
    default parameters as `autumn::harvest()`.
  - All 7 `autumn`-exported functions present: `harvest`, `anesrake`, `diagnose_weights`,
    `design_effect`, `effective_sample_size`, `get_current_miss`, `weighted_pct`.
  - Accepted **and honoured**: `attach_weights`, `weight_column`, `max_weight`,
    `max_iterations`, `verbose`, `start_weights`, `target_map`.
  - Accepted **and silently ignored** (noted at `verbose >= 2`): `select_params`,
    `select_function`, `error_function`, `adaptive_order`, `enforce_mean`, `accelerate` —
    rationale: the solvers calibrate all variables simultaneously and cannot honour
    per-variable selection logic.
  - Accepted **with a deprecation warning / informative error**: `convergence["pct"]`;
    `auto_collapse`, `collapse_vars`, `add_na_proportion`.
  - Calibrated weights satisfy `max(w) <= max_weight` within 1e-10.
  - Weights normalized to mean ~= 1 on output.
  - `R CMD check` passes with 0 errors, 0 warnings.
  - FR-30: `harvest()` normalizes `start_weights` to mean 1 before passing to the C API.
  - FR-31: factor/character columns encoded to 0-indexed `int32` by the R bridge;
    NA -> -1; OOV (level absent from target) -> -1.
  - FR-32: return structure identical to `autumn::harvest()`, plus
    `attr(data, "algorithm")` set to the algorithm name string.
  - FR-33: `useDynLib(leafblower, .registration=TRUE)`; routines registered via
    `R_registerRoutines()` in `zzz.R`.
  - FR-34: `anesrake()` wraps `harvest()` with the mapping `inputter`->`data`,
    `weightvec`->`start_weights` (NULL -> uniform), `caseid` ignored, `targets`->`target`,
    `pctlim`->`convergence["pct"]` (deprecated), `cap`->`max_weight`, `type` ignored,
    `nlim`->`max_iterations`, `iterate` ignored, `threads` ignored; unknown anesrake
    params warn and are ignored.
  - FR-35: `get_current_miss()` exported under exactly that name.
- scope: R package surface, autumn compatibility, `R/harvest.R`, `R/anesrake.R`, NAMESPACE
- supersession: the `method = "rake"` / `method = "nr"` acceptance criteria (map to
  `"lbfgsb"` with a named warning) are WITHDRAWN — see § Superseded / withdrawn. The
  `choosemethod` leg of the FR-34 mapping is withdrawn for the same reason. Everything
  else in this entry stands.

## REQ-us002-min-weight-lower-bound
- source: tasks/prd-leafblower-core.md § US-002 (also § G-3, FR-5)
- description: A survey methodologist specifies `min_weight` so that calibrated weights
  never fall below a floor, preventing near-zero observations from being effectively
  excluded. This is the package's headline addition over `autumn` (neither `autumn` nor
  `anesrake` exposes a lower bound).
- acceptance:
  - `harvest(data, target, min_weight = 0.5)` produces `min(w) >= 0.5` within 1e-10.
  - Default `min_weight = 0` means no lower bound; the clamp `max(min_weight, w_i)` is a
    no-op at 0 because weights start positive (FR-5).
  - `min_weight` accepted by `rk_calibrate()` as `rk_params_t.min_weight` and by Python
    `harvest()`.
  - `min_weight >= max_weight` returns `RK_ERR_BADARG` immediately with the message
    "min_weight must be strictly less than max_weight".
  - `diagnose_weights()` output is unchanged by `min_weight`.
- scope: bounds enforcement, C API params, R + Python surfaces
- supersession: the AC "Routing: `min_weight > 0` forces iEPPA (see US-007)" is
  superseded — AUTO routing now spans eight solvers per
  `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md`. The bound
  requirement itself is unaffected.

## REQ-us003-large-scale-performance
- source: tasks/prd-leafblower-core.md § US-003 (also § G-2, § 1 Success Criteria, § 11)
- description: A census microsimulation researcher calibrates 1M+ observations across
  20+ margins in under 30 seconds, to allow fast iteration on synthetic population models.
- acceptance:
  - `harvest(data, target)` with n = 1,000,000, K = 20, 5 categories each completes in
    < 30 s wall-clock, single-threaded, on a modern 8-core CPU.
  - Medium case: 100K rows, 5 margins (3–5 categories each) < 1 s (§1 Success Criteria;
    § 11 states < 2 s for the same shape — see note below).
  - Convergence criterion met (`max calibration error < 1e-6`) or an informative warning
    is issued.
  - `verbose = 1` prints the selected algorithm name and the routing reason.
- scope: performance targets, benchmark gates, AUTO routing observability
- note: the PRD states the medium-scale target twice with different numbers — §1
  Success Criteria says "100K rows, 5 margins < 1 s", §11 KPI says "100K rows, 5 margins
  < 2 s". Both were written against the removed `lbfgsb` solver. Recorded as INFO, not a
  competing variant: the target must be re-benchmarked against a live solver before it
  can be used as a gate.

## REQ-us004-c-api-contract
- source: tasks/prd-leafblower-core.md § US-004 and FR-1 … FR-10, § 6 Design
  Considerations, § 7 Memory Contract / Thread Safety / Convergence Reporting
- description: A stable C API header (`leafblower.h`) lets any language with C FFI call
  the calibration engine without duplicating algorithm code. Both the R `.Call()` bridge
  and the Python pybind11 module call the same `rk_calibrate()` symbol.
- acceptance:
  - `leafblower.h` is a valid C99 header (compiles under `gcc -std=c99 -Werror`); no C++
    symbols leak (everything inside `extern "C"`).
  - FR-1: `rk_calibrate(n, K, weights, group_ids, cat_counts, targets, params, result)`
    calibrates `weights[n]` **in place** and returns a status code.
  - FR-2: `rk_params_init(p)` fills defaults `min_weight=0.0`, `max_weight=5.0`,
    `algorithm=RK_ALG_AUTO`, `inner_max_iter=500`, `outer_max_iter=50`, `tol_abs=1e-6`,
    `epsilon=0.05`, `verbose=0`, `log_fn=NULL`. The R bridge maps autumn's
    `max_iterations` -> `inner_max_iter`.
  - FR-3: `group_ids[k]` is a pointer to `n` contiguous `int32`;
    `group_ids[k][i]` in `{-1, 0, …, cat_counts[k]-1}`; `-1` = NA/OOV (observation
    skipped for margin k). `targets[k]` is `cat_counts[k]` contiguous doubles.
  - FR-4 (validation — **all** return `RK_ERR_BADARG`, all checked **before** any weight
    is modified): NULL pointer for `weights`/`group_ids`/`cat_counts`/`targets`;
    `n <= 0` or `K <= 0`; any `cat_counts[k] <= 0`; any `cat_counts[k] > n`; any
    `group_ids[k][i] >= cat_counts[k]`; any `group_ids[k][i] < -1`; NaN or Inf in
    `targets[]` or initial `weights[]`; any `targets[k]` not summing to 1 +/- 1e-8;
    `min_weight >= max_weight`; allocation-overflow guard —
    `total_cats = sum(cat_counts)` as `size_t`, reject if
    `(size_t)n * total_cats > SIZE_MAX / 2` with "problem too large for platform size_t";
    zero-weight guard — `sum(weights) < 1e-15` rejected with "total weight is zero or
    negative". The `group_ids` range check is a full O(n*K) pass.
  - FR-6: on `RK_ERR_NOCONV`, `weights[]` reflects the last iterate and
    `rk_result_t.max_error` the last computed calibration error.
  - FR-8: `rk_result_t.message[256]` is always null-terminated even on truncation — use
    `snprintf`, never `sprintf`.
  - FR-9: verbose output goes through the `log_fn` callback when non-NULL, else
    `fprintf(stderr, …)` when `verbose > 0`. The R bridge passes an `Rprintf`-backed
    callback. **Python GIL contract:** `rk_calibrate()` is called with the GIL held (no
    `py::gil_scoped_release`); `_bindings.cpp` registers a static C trampoline
    `void py_log_trampoline(const char* msg, void* ctx)` that casts `ctx` to a borrowed
    `PyObject*` callable kept alive for the call; because the GIL is held throughout, no
    `PyGILState_Ensure/Release` is needed.
  - FR-10 / `RK_ERR_INFEAS`: triggered when a category has `S_j < 1e-15 * W` with
    `tau_j > 0` at the end of an inner sweep (empty cell with positive target —
    geometrically impossible). `min_weight >= max_weight` is `RK_ERR_BADARG`, NOT
    infeasibility.
  - Memory contract: caller owns all arrays; `rk_calibrate()` may allocate internal
    workspace of size n and sum(c_k) and frees all of it before return; no heap
    allocation survives the call; `rk_result_t` is caller-allocated.
  - Thread safety: `rk_calibrate()` is **reentrant** — no global or static mutable state;
    concurrent calls from separate threads on separate data are explicitly supported
    (R `parallel` / `future` is a supported caller). `log_fn` must stay valid for the
    whole call; bridges must not pass stack-allocated or GC-eligible pointers.
  - Convergence reporting: `rk_result_t.max_error =
    max_k max_j |sum_{i: g_k(i)=j} w_i / sum w_i - tau_j^(k)|` on the final iterate.
- scope: `src/leafblower.h`, `src/c_api.cpp`, `src/r_bridge.cpp`,
  `python/leafblower/_bindings.cpp`
- supersession: the base error codes `RK_OK=0`, `RK_ERR_NOCONV=1`, `RK_ERR_INFEAS=2`,
  `RK_ERR_BADARG=3` stand and are EXTENDED (not replaced) by `RK_ERR_BUDGET=4` /
  `RK_ERR_STALL=5` from `2026-04-28-convergence-status-design.md`. The `rk_algorithm_t`
  enum given in § 6, the `lbfgs_m` default in FR-2, the FR-7 statement that
  `algorithm_used` is "`RK_ALG_IEPPA` or `RK_ALG_LBFGSB`", and the FR-4 logit-singularity
  guard (`min_weight == 1` / `max_weight == 1`, which exists only for the L-BFGS-B logit
  link) are all superseded — see § Superseded / withdrawn. `rk_params_t` and
  `rk_result_t` have since grown many fields; the PRD layout in § 6 is the v1 seed, not
  the current ABI.

## REQ-us005-oris-capacity-constrained-solver
- source: tasks/prd-leafblower-core.md § US-005 (FR-11 … FR-19)
- description: A capacity-constrained multi-marginal calibration solver with
  cell-compressed O(M_cell * K) inner cost. Named `iEPPA` in the PRD; this is the solver
  now named **ORIS** (renamed per
  `docs/superpowers/specs/2026-05-30-oris-rename-design.md`; enum value 1 unchanged).
  Same solver, not a competing variant.
- acceptance:
  - Cell compression: unique `(g_1, …, g_K)` tuples deduplicated via a sort-based cell
    table; `M_cell <= min(n, prod cat_counts)`; one cell-level weight expansion per outer
    iteration (O(M_cell*K) core, O(n) expansion).
  - Log-space Sinkhorn factors `lf[k][j]` per margin-category;
    `X_tilde[c] = X_init[c] * exp(sum_k lf[k][g_k(c)])`.
  - Log-sum-exp stabilization: partial sums clipped to `[-Inf, 700]` (IEEE double
    overflows at `exp(709)`).
  - Capacity block: `X[c] = clamp(X_tilde[c] * W[c], L_c, U_c)` with
    `L_c = min_weight * n_per_cell[c]`, `U_c = max_weight * n_per_cell[c]`.
  - Convergence: `errRp = max_k max_j |S_kj / W_total - tau_kj| < tol_abs` (default 1e-6).
  - FR-13(c): empty-category test uses `S_j < 1e-15 * W`, never exact `== 0` (denormal
    sums are non-zero but near machine epsilon).
  - FR-16: NA/OOV rows (`group_ids[k][i] = -1`) are skipped in `S_j` and their weight is
    not updated for that margin pass.
  - FR-19 / US-005 last AC: all weights satisfy `[min_weight, max_weight]` after **every**
    inner step — an invariant, not just an exit condition.
  - FR-17/FR-18: outer loop exhausted -> `RK_ERR_NOCONV` with
    "outer loop exhausted after N iterations, max_error=E"; `inner_max_iter` (default 500)
    caps sweeps per outer iteration regardless of the stopping criterion.
  - Overflow detection: if log-factor drift `max_log_X_tilde > 700` AND an uncapped cell
    is detected, return `RK_ERR_NOCONV` with an actionable message (direct the user to
    `raking` or looser bounds).
  - Verbose: `verbose=1` prints the compression ratio, per-iteration `errRp` and exit
    status; `verbose=2` adds `n_cap_active` and the `log10(f[k][j])` range per margin.
- scope: `src/oris.cpp` (was `src/ieppa.cpp`), cell compression, bounds invariant
- supersession: **partial.** The PRD's framing — "paper-faithful iEPPA (Chu-Liang-Toh-Yang
  2022, arXiv:2011.14312) at C=0" with an outer entropic-proximal-point loop (FR-11,
  FR-12, FR-14 `tolRb`/`normU`/Bregman-distance inner stop) — is explicitly repudiated by
  `docs/methods/oris.md` and by the rename spec: that outer loop is NOT implemented, is
  mathematically inert at `C = 0`, and citing the paper over-claims. Treat FR-11/FR-12/
  FR-14 as historical. The acceptance criteria listed above are the durable part and are
  corroborated by `docs/methods/oris.md` and `docs/methods/00-overview.md`.

## REQ-us005b-classical-raking
- source: tasks/prd-leafblower-core.md § US-005b
- description: `method="raking"` provides classical cyclic IPF (Deming-Stephan 1940;
  Csiszar 1975) with box projection and hyperplane projection, as an alternative to the
  ORIS solver. This is the requirement the SPEC corpus cites as "§ US-005b".
- acceptance:
  - `rk_calibrate(..., algorithm=RK_ALG_RAKING)` present in the enum.
  - `method="raking"` in R and Python routes to `raking_solve`.
  - All pre-rev2 iEPPA tests pass against `method="raking"` (regression guard in
    `test-raking.R`).
- scope: `src/raking.cpp`, `method="raking"`, R + Python routing
- supersession: the projection GEOMETRY named here (additive Dykstra box projection,
  Boyle-Dykstra 1986) is superseded by
  `docs/superpowers/specs/2026-04-27-raking-bregman-dykstra-design.md`, which replaces the
  Euclidean corrections with multiplicative KL-Bregman projections at cell level. The
  existence and routing requirement is unaffected.

## REQ-us008-python-pandas-interface
- source: tasks/prd-leafblower-core.md § US-008 and FR-36 … FR-40 (also § G-4)
- description: A Python survey analyst calls `leafblower.harvest(df, targets)` with a
  pandas DataFrame and gets calibrated weights in a Pythonic workflow, over the same
  compiled C core as R (§ G-5: no algorithm code duplicated between the bindings).
- acceptance:
  - FR-36: `harvest(data, targets)` auto-detects `pd.DataFrame` via `isinstance`; also
    accepts `{col: list}` (converted to a DataFrame internally). `targets` is a
    dict-of-dicts, e.g. `{"age": {"18-34": 0.3, "35-54": 0.45, "55+": 0.25}}`.
  - `targets` sum-to-1 validation happens in Python before the C call and raises
    `ValueError` naming the offending column.
  - `len(group_ids[k]) == n` validated for all k before calling `rk_calibrate()`.
  - FR-37: categorical columns encoded as 0-indexed `numpy.int32`; NA -> -1.
  - Returns a DataFrame with a `weights` column appended, or a `numpy.ndarray`
    (`float64`) when `attach_weights=False`. The `weights` array is **always a copy**,
    never a view into the input buffer; the caller may mutate it freely.
  - Parameters present: `min_weight`, `max_weight`, `method`, `verbose`,
    `max_iterations`, `convergence`. `convergence` is a dict accepting `"absolute"`
    (-> `tol_abs`); `"pct"` raises `DeprecationWarning` and is ignored; any other key
    raises `ValueError: unknown convergence key '{k}'`.
  - FR-39: pybind11 module `_leafblower` exposes
    `calibrate(n, K, weights_np, group_ids_list, cat_counts_list, targets_list,
    params_dict)` returning `(status, weights_out, result_dict)`; `weights_out` is a NEW
    C-contiguous `float64` ndarray, never a view; the binding enforces `float64` +
    C-contiguity before taking the raw pointer, raising `ValueError` otherwise.
  - FR-38: `python/CMakeLists.txt` explicitly EXCLUDES `r_bridge.cpp` (it includes
    `Rinternals.h`, unavailable in the Python build).
  - FR-40: `diagnose_weights()` returns a `pd.DataFrame` with the same columns as the R
    equivalent; `design_effect()` and `effective_sample_size()` are pure Python.
  - `pip install leafblower` installs a self-contained wheel (no separate C library).
  - Python >= 3.9; CI tests 3.9, 3.10, 3.11, 3.12, 3.13.
  - Import path `from leafblower._leafblower import calibrate`; `__init__.py` re-exports
    the public `leafblower.harvest`.
- scope: `python/leafblower/_harvest.py`, `python/leafblower/_bindings.cpp`,
  `python/CMakeLists.txt`, wheel packaging
- note: FR-38's exclusion rule is the inverse of the live convention recorded in
  `CLAUDE.md` — the Python build does NOT glob; `python/CMakeLists.txt` carries an
  explicit `CORE_SOURCES` list, so a new `src/*.cpp` must be ADDED there or the pybind11
  link fails. Same outcome (no `r_bridge.cpp` in the Python build), different mechanism.

## REQ-us009-diagnostic-functions
- source: tasks/prd-leafblower-core.md § US-009 (also FR-35, FR-40)
- description: `diagnose_weights()` and `design_effect()` behave exactly as in `autumn`
  so that calibration quality can be evaluated without changing analysis code.
- acceptance:
  - `diagnose_weights(data, target, weights)` returns a data frame with columns
    `variable`, `level`, `prop_original`, `prop_weighted`, `target`, `error_original`,
    `error_weighted`.
  - `design_effect(weights)` returns the Kish (1992) estimator
    `n * sum(w^2) / (sum(w))^2`.
  - `design_effect(weights, outcome, data, target)` returns the Henry & Valliant (2015)
    estimator.
  - `effective_sample_size(weights) = length(weights) / design_effect(weights)`.
  - `get_current_miss(data, target, weights)` exported under exactly that name.
  - `weighted_pct(x, weights)` returns a named numeric vector of weighted proportions.
  - R implementations are pure R — no C call required for diagnostics.
- scope: `R/diagnose_weights.R`, `R/design_effect.R`, `R/weighted_pct.R`,
  `R/current_miss.R`, Python equivalents

## REQ-us010-cran-pypi-distribution
- source: tasks/prd-leafblower-core.md § US-010 (also § G-6, § 10 Phase 4)
- description: leafblower is distributable on CRAN and PyPI so users install it with
  `install.packages()` and `pip install`.
- acceptance:
  - `R CMD check --as-cran` passes with 0 errors, 0 warnings, <= 1 note (new-submission
    note acceptable).
  - `configure` detects C++17 support and generates `src/Makevars`; falls back to C++14
    with a printed note rather than failing to compile.
  - No vendored dependency > 5 MB.
  - `.Rbuildignore` excludes `.beads/`, `.claude/`, `.wolf/`, `tasks/`, `python/`.
  - `pyproject.toml` with scikit-build-core produces a valid wheel via `pip wheel .`;
    the wheel bundles the compiled extension and
    `python -c "import leafblower; leafblower.harvest"` succeeds after install.
  - Python package version matches the R `DESCRIPTION` Version field.
  - `cran-comments.md` documents any remaining notes for CRAN reviewers.
- scope: `configure`, `DESCRIPTION`, `.Rbuildignore`, `python/pyproject.toml`, `NEWS.md`
- supersession: the § 7 build-flag statement `PKG_CXXFLAGS = ... -O3` is superseded — the
  R build now deliberately sets NO `-O` level (CRAN's `tools:::.check_make_vars` rejects
  `-O` flags in `PKG_CXXFLAGS`; R supplies the user/site level via `$(CXXFLAGS)`). Only
  `python/CMakeLists.txt` sets `-O3`. The version-sync AC stands and is manual — there is
  no automation (`DESCRIPTION` and `python/pyproject.toml` are bumped by hand).
- note: the `.Rbuildignore` list also names `docs/iEPPA/`, a directory the rename removed.

## REQ-kpi-success-metrics
- source: tasks/prd-leafblower-core.md § 11 Success Metrics & KPIs (and § 1 Success
  Criteria)
- description: The project-level acceptance table. Every row names both a target and the
  artefact that measures it.
- acceptance:
  - R API compat — all 7 autumn functions present — measured by `R CMD check` +
    `test-harvest.R`.
  - Weight bound enforcement — `max(w) <= max_weight` and `min(w) >= min_weight`,
    within 1e-10, over 50 random datasets — measured by a property-based test in
    `test-harvest.R`.
  - Convergence — `max_k max_j |sum w * 1[g_k=j] / sum w - tau_j^(k)| < 1e-6` at reported
    convergence.
  - Large-scale — 1M rows, 20 margins, `max_weight = 3`, < 30 s — measured by the
    Phase-2 gate.
  - CRAN check — 0 errors, 0 warnings — measured by `R CMD check --as-cran`.
  - Python wheel — installs on Linux/macOS, Python 3.9–3.13 — measured by the CI matrix.
  - Performance benchmarks are CI **artifacts, not gates** (§ 8); only the phase gates
    are binding.
- scope: project-level definition of done
- supersession: the "L-BFGS-B convergence — 100K rows, 5 margins < 2 s" row is WITHDRAWN
  with the solver. Its performance target survives only as the unattributed medium-scale
  target in `REQ-us003-large-scale-performance`, where it must be re-benchmarked against a
  live solver before use as a gate.

---

## Superseded / withdrawn PRD content (recorded, not extracted)

None of the following is a conflict requiring resolution. The PRD is the oldest document
in the corpus and precedence already places SPEC above PRD; date agrees. Each item is
recorded as INFO in `../INGEST-CONFLICTS.md`.

1. **US-006 (whole user story) and FR-20 … FR-28 — L-BFGS-B. WITHDRAWN.** The solver was
   removed; `RK_ALG_LBFGSB = 2` is a permanently reserved hole
   (`docs/superpowers/specs/2026-08-14-removed-solver-slots-supersession.md`,
   `src/leafblower.h:44`). Everything specific to it goes with it: the Deville-Sarndal
   logit/exponential link selection (FR-20 … FR-22), the dual `phi(lambda)` formulation
   and gradient (FR-23 … FR-25), the L-BFGS 2-loop recursion and `lbfgs_m` (FR-26), the
   Wolfe line search `c1=1e-4`/`c2=0.9` (FR-27), and the FR-28 final clamp. Note the
   Deville-Sarndal logit DISTANCE lives on in the separate `logit` solver
   (`RK_ALG_LOGIT = 10`, `src/logit_calib.cpp`) — a different solver, not a survival of
   US-006.
2. **§ 6 enum `RK_ALG_AUTO=0, RK_ALG_IEPPA=1, RK_ALG_LBFGSB=2`. SUPERSEDED.** This PRD is
   the ORIGIN of the slot-2 hole. The live enum is AUTO 0, ORIS 1, RAKING 3, SINKHORN 4,
   CHEBYSHEV 5, GREG 6, ORIS_SOFT 8, GREENKHORN 9, LOGIT 10, NEWTON_KL 11 — eight solvers,
   values frozen, slots 2 and 7 permanently reserved.
3. **US-007 auto-routing (select between `RK_ALG_IEPPA` and `RK_ALG_LBFGSB` on
   `complexity = n * sum(cat_counts)` vs 500000, `max_weight < 3.0`, `min_weight > 0`).
   SUPERSEDED.** Routing now spans eight solvers; the current rule is the three-way test
   in `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` (K >= 5,
   `M_cell/n >= 0.9`, `target_skew` vs 5.0). FR-7's "never `RK_ALG_AUTO`" invariant
   survives; its enumeration of possible values does not.
4. **US-001 `method="rake"` / `method="nr"` -> `"lbfgsb"` with warning. SUPERSEDED.**
   `lbfgsb` does not exist. The autumn synonym handling still needs a destination; the
   PRD no longer supplies one.
5. **`iEPPA` throughout = the solver now named `oris`.** Renamed per
   `2026-05-30-oris-rename-design.md`; enum value 1 unchanged. Terminology, not a
   contradiction.
6. **§ 5 Non-Goals — "SQUAREM acceleration (not ported)" and "Bounded IPF (water-filling)
   as a named method (not exposed)". BOTH SHIPPED SINCE.** SQUAREM landed via
   `2026-04-28-raking-squarem-design.md` (and has since itself been superseded by SRAA-m);
   water-filling is the ORIS `bounds_mode="unit"` path. The non-goals are superseded
   scope, NOT violated requirements.
7. **§ 5 Non-Goals — `auto_collapse` and `add_na_proportion` "not in v1".** Later work
   exists on both. Superseded scope, not a contradiction. The US-001 AC that they "raise
   an informative error" is therefore also stale.
8. **§ 7 Build System — `PKG_CXXFLAGS = ... -O3`. SUPERSEDED.** The R build sets no `-O`
   level of its own by design (CRAN portability check); only `python/CMakeLists.txt` sets
   `-O3`. The `PKG_SOURCES` list in `Makevars.in` is likewise decorative — R auto-globs
   `src/*.cpp`.
9. **§ 5 Non-Goal "Multi-threading — single-threaded in v1"** is consistent with the
   thread-safety contract in REQ-us004 (reentrancy is a caller-side guarantee, not
   internal parallelism) and with the single-thread-BLAS requirement for deterministic
   R<->Python parity. Recorded for the roadmapper: not superseded, still in force.
10. **§ 9 Risks and § 10 Phased Rollout (Phases 1–4).** Phase 1 is L-BFGS-B and is void;
    Phases 2–4 describe work long since delivered. The TDD RED-phase discipline they
    encode (write and verify the test FAILS before implementing) is the project's standing
    methodology and outlives the phase plan.
11. **§ 12 Open Questions OQ-1 … OQ-5.** OQ-1 (rake/nr mapping) and OQ-2 (logit degenerate
    case) resolved into withdrawn L-BFGS-B content. OQ-3 (fixed `epsilon = 0.05`), OQ-4
    (`log_fn` callback) and OQ-5 (Windows source-only for v1) stand.
