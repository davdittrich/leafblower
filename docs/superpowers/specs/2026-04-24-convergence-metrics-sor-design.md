# Convergence Criterion Reform + SOR Stabilization + Best-Iterate Design

**Date:** 2026-04-24
**Status:** Draft — revision 3, pending design-review-gate re-run
**Applies to:** all three solvers (iEPPA, raking, lbfgsb) unless noted

## Problem

iEPPA solver on stepstone-fulldata (1.58M rows, K=9, max_weight=5) diverges:
errRp increases from ~2.22e-3 at iter 50 to ~6.5e-3 at iter 3000. Root causes:

1. **Wrong stopping criterion:** `errRp` (L∞ marginal error) oscillates in the
   capacity-clamp non-smooth regime. Solver declares NOCONV but reports the
   diverged final iterate rather than the best-seen iterate.
2. **No oscillation damping:** Fixed multiplicative step size ω=1.0 in iEPPA.
3. **Incomplete quality reporting:** Only `max_error` returned across all solvers.

## Goals

1. **pct = 0.001 as default convergence criterion** across all three solvers —
   oscillation-resistant, industry-standard (GREG, survey::calibrate, anesrake).
2. **User-selectable criteria:** `pct`, `max_err`, `mean_err`, `kl`, `chi2`.
3. **All five quality metrics** always reported at exit across all three solvers.
4. **SOR adaptive under-relaxation** (iEPPA only) to dampen oscillation.
5. **Best-iterate tracking** (W at minimum observed errRp) across all solvers.

## Non-goals

- Do not change iEPPA/raking/lbfgsb core algorithm structures.
- No SOR for raking (convergence guaranteed) or lbfgsb (has Wolfe line search).
- No Polyak-Ruppert averaging.

## Migration

**Naming collision:** `harvest.R:244` currently warns on `convergence["pct"]` as
deprecated and discards it. That deprecation warning is **removed** and
`parse_convergence()` is **rewritten** to treat `pct` as a primary key.

**Silent semantic reversal risk:** users who previously set `convergence =
list(pct = X)` received a deprecation warning and had X silently discarded (only
`absolute` was honored). After this change, `pct = X` will actively fire. Add
to NEWS.md: *"Users who suppressed the convergence['pct'] deprecation warning:
your pct value will now be honored as a convergence threshold."*

**Default criterion change:** `convergence = list()` (or no `convergence` arg)
previously activated `max_err < 1e-6`. Now activates `pct < 0.001`. All callers
not passing an explicit `convergence` argument will switch criterion.

**Affected test files (6):** `test-ieppa.R`, `test-compare.R`,
`test-config-defaults.R`, `test-ieppa-nonuniform-d.R`, `test-compat.R`,
`test-bounded-convergence.R`. Implementer MUST audit: if any assertion checks
`calib_result$max_error` or iteration count without pinning convergence, add
`convergence = list(absolute = 1e-6)` to that call. A1 verified by
`devtools::test()` FAIL 0 after audit.

**Python:** mirror the R/CRAN migration note in the Python package changelog.
Remove any equivalent `pct`-deprecated handling in `python/leafblower/_harvest.py`.

**CRAN NEWS.md template:**
> Convergence default changed from `absolute = 1e-6` (max marginal error) to
> `pct = 0.001` (max proportional weight change). Old behavior:
> `convergence = list(absolute = 1e-6)`. New criterion string: `"kl"` (was
> unpublished internal). New result fields: `mean_error`, `kl`, `chi2`,
> `pct_change`, `best_error`, `best_iter`, `best_weights`.

---

## Design

### §1 Convergence Criterion System (all three solvers)

**Default:** `convergence = list(pct = 0.001)` — stop when
`max_i |W_i[t] - W_i[t-1]| / W_i[t-1] < 0.001`.

**Pluggable criterion:**

```r
convergence = list(pct = 0.001)                         # new default
convergence = list(absolute = 1e-6)                     # backward compat (max_err)
convergence = list(absolute = 1e-3, criterion = "chi2")
convergence = list(absolute = 1e-3, criterion = "mean_err")
convergence = list(absolute = 1e-3, criterion = "kl")
convergence = list(pct = 0.001, absolute = 1e-6, stop_when = "any")
```

**Criterion strings** (lowercase, `match.arg()`-validated):
`"pct"` / `"max_err"` / `"mean_err"` / `"kl"` / `"chi2"`.

**`stop_when`:** `"any"` (stop when EITHER threshold fires, default) or `"all"`
(stop only when BOTH fire). Silently ignored when only one threshold is non-zero.

**New `CalibConvergenceCfg` struct (added to `CalibState`):**

| field | type | default |
|---|---|---|
| `pct_tol` | double | 0.001 |
| `absolute_tol` | double | 0.0 |
| `criterion` | enum | PCT |
| `stop_when` | enum | ANY |

**`parse_convergence()` contract** (R/harvest.R rewrite — replaces deprecated stub):

```r
parse_convergence <- function(convergence) {
  explicit_pct <- !is.null(convergence[["pct"]])
  explicit_abs <- !is.null(convergence[["absolute"]])
  # pct default 0.001 only when neither pct nor absolute is specified (bare list()).
  # When user specifies absolute without pct, pct is OFF (backward compat).
  pct_tol      <- if (explicit_pct) convergence[["pct"]]
                  else if (!explicit_abs) 0.001
                  else 0.0
  absolute_tol <- convergence[["absolute"]] %||% 0.0
  criterion    <- match.arg(
    convergence[["criterion"]] %||% (if (explicit_pct || !explicit_abs) "pct" else "max_err"),
    c("pct","max_err","mean_err","kl","chi2"))
  stop_when    <- match.arg(convergence[["stop_when"]] %||% "any", c("any","all"))
  list(pct_tol=pct_tol, absolute_tol=absolute_tol, criterion=criterion, stop_when=stop_when)
}
```

Edge case verification:
- `list()` → pct_tol=0.001, abs=0.0, criterion="pct" (new default)
- `list(absolute=1e-6)` → pct_tol=0.0, abs=1e-6, criterion="max_err" (backward compat ✅)
- `list(pct=0.001)` → pct_tol=0.001, abs=0.0, criterion="pct" ✅
- `list(pct=P, absolute=A)` → pct_tol=P, abs=A, criterion="pct", stop_when="any" ✅

**Alternative criteria formulas** (computed at `kErrCheckInterval` cadence):

- **max_err:** `max_k max_j |Ŝ_kj/W_total - T_kj|`
- **mean_err:** `mean_k(max_j |Ŝ_kj/W_total - T_kj|)`
- **kl:** `max_k Σ_j T_kj * log((T_kj + ε) / (Ŝ_kj/W_total + ε))`, ε = 1e-10;
  skip j where T_kj = 0 (contribution = 0)
- **chi2:** `Σ_k Σ_j (Ŝ_kj - T_kj * W_total)² / (T_kj * W_total + 1)`;
  floor of 1 in denominator. Note: scales with W_total ≈ n — user must supply
  n-scaled `absolute` threshold. chi2 can return `Inf` on degenerate cells;
  this is a diagnostic value (non-crashing in R). Document in roxygen.

The **active criterion** (per `criterion` field) is compared to the threshold for
stopping. All criteria are computed regardless and returned in §2 metrics.

**pct computation** — runs at `kErrCheckInterval`, after capacity block. See
solver-local storage in §8 for `W_prev` details.

```cpp
double pct_change = 0.0;
for each weight unit (see §8 for per-solver indexing):
    double rel = fabs(W_current - W_prev) / max(W_prev, 1e-12);
    if (rel > pct_change) pct_change = rel;
W_prev = W_current;   // update AFTER computing pct_change
```

---

### §2 Quality Metric Reporting (all three solvers)

Computed at every solver's exit from the **final iterate** weights. Always in
`calib_result`:

| R field | metric | formula |
|---|---|---|
| `max_error` | L∞ marginal | unchanged `errRp` |
| `mean_error` | L1-over-margins | `mean_k(max_j |Ŝ_kj/W_total - T_kj|)` |
| `kl` | max KL | `max_k Σ_j T_kj log((T_kj+ε)/(Ŝ_kj/W_total+ε))` |
| `chi2` | total χ² | `Σ_k Σ_j (Ŝ_kj-T_kj·W_total)²/(T_kj·W_total+1)` |
| `pct_change` | weight pct | last iteration's `max_i |ΔW_i/W_i|` |

**Naming:** field `kl` matches criterion string `"kl"` exactly. Field `chi2`
matches criterion string `"chi2"`. No mismatch. Old field names (`kl_divergence`,
`chi_square`) do not exist and are not introduced.

**All three solvers** (iEPPA, raking, lbfgsb) compute all five metrics at exit.
For raking and lbfgsb, margins Ŝ_kj are accumulated from obs-level weights
(no cell compression); the formulas are identical since they only reference Ŝ_kj
and T_kj proportions.

---

### §3 SOR Adaptive Under-Relaxation (iEPPA only)

Replace iEPPA margin step:

```cpp
// Before:
f_k_new[j] = f_k_old[j] * (T_kj / S_kj);

// After (with guard and performance branch):
if (S_kj < kEmptyBucketThreshold * ct.W_input) {
    // skip — empty bucket already guarded by existing linear-path logic
    // kEmptyBucketThreshold = 1e-15 (ieppa.cpp:65)
} else {
    const double ratio = T_kj / S_kj;
    const double r = (omega[k] == 1.0) ? ratio : std::pow(ratio, omega[k]);
    f_k_new[j] = f_k_old[j] * r;
}
```

The `omega[k] == 1.0` branch prevents `std::pow()` on the non-damped path
(not all R-target compilers elide `pow(x,1.0)` at -O2; this guarantee is unsafe
cross-platform).

**Default: `sor = list(auto = TRUE, omega_min = 0.3)`**

Auto mode (local to `ieppa_solve()`, K-element arrays `omega[K]` and
`prev_errRp_k[K]` — see §8):

- `omega[K]` initialized to `omega_init = 1.0`
- `prev_errRp_k[K]` initialized to ∞
- Burn-in: first `burnin = 20` outer iterations — ω unchanged
- After burn-in, at each `kErrCheckInterval`:
  - Compute per-margin errRp_k = `max_j |Ŝ_kj/W_total - T_kj|`
  - Sign flip (errRp_k[t] > errRp_k[t-1] AND prev was decreasing):
    `omega[k] = max(omega_min, omega[k] * 0.7)`
  - Monotone decrease: `omega[k] = min(1.0, omega[k] * 1.05)`
  - `prev_errRp_k[k] = errRp_k[t]`

**omega_min = 0.3 rationale:** Conservative empirical lower bound from pilot
runs on stepstone-fulldata: ω < 0.3 over-dampens non-oscillating margins,
increasing wall time >2×; ω = 0.3 stabilizes oscillation within ~100 iters.
No analytic SOR theorem applies to non-linear clamped iEPPA; 0.3 is the
empirical safe floor. Reproducible via `benchmarks/stepstone_fulldata_homotopy.R`
with `sor = list(omega = X)` sweep.

**R API:**

```r
sor = list(auto = TRUE, omega_min = 0.3)    # default
sor = list(omega = 0.5)                     # fixed uniform ω, no per-margin adaptation
sor = NULL                                  # disable (sole off state)
```

**`CalibSorCfg` struct (added to `CalibState`):**

| field | type | default |
|---|---|---|
| `enabled` | bool | true |
| `auto_adapt` | bool | true |
| `omega_init` | double | 1.0 |
| `omega_min` | double | 0.3 |
| `omega_fixed` | double | -1.0 (sentinel: use auto) |
| `burnin` | int | 20 |

**SOR diagnostics in `calib_result$sor`** (nested list, mirrors input structure):

```r
calib_result$sor$min_omega   # lowest ω across any margin × any iter (1.0 if SOR off)
calib_result$sor$n_damped    # total (margin, iter) pairs where ω < 1 (0 if SOR off)
```

**C ABI additions to `rk_result_t`:**

```c
double sor_min_omega;   /* 1.0 for non-iEPPA; zero-initialized by c_api.cpp */
int    sor_n_damped;    /* 0 for non-iEPPA; zero-initialized by c_api.cpp */
```

Raking and lbfgsb callers never set these fields; `c_api.cpp` zero-initializes
`rk_result_t` so sentinel values (1.0, 0) are automatic.

---

### §4 Best-Iterate Tracking (all three solvers)

Always enabled. Storage: local vectors, O(size) per-solver copy-on-improvement.

**Implementation (solver-local, NOT CalibState — see §8 for size):**

```cpp
// Locals:
double best_errRp = std::numeric_limits<double>::infinity();
int    best_iter  = 0;
std::vector<double> W_best(SIZE, 0.0);  // SIZE per solver: see §8

// After errRp computed at kErrCheckInterval:
if (errRp < best_errRp) {
    best_errRp = errRp;
    best_iter  = iter;
    W_best     = W_current;
}
```

**W_best expansion at exit — iEPPA:**

```cpp
// Expand cell-level W_best to obs-level best_weights:
std::vector<double> best_weights(n);
for (int i = 0; i < n; i++)
    best_weights[i] = st.weights[i] * W_best[ct.cell_of[i]];  // scalar mult
// Normalize sum = n (same normalization as W_final):
double s = 0.0;
for (double w : best_weights) s += w;
for (double &w : best_weights) w *= (n / s);
// DO NOT apply water-fill or bounds clamping to W_best.
// W_best is a mid-loop snapshot; re-clamping changes its meaning.
```

**W_best expansion at exit — raking and lbfgsb:** W_best is already obs-level
(`n`-element). Apply only the sum-normalization (`n / sum(W_best)`). No other
post-processing.

**Return in `calib_result`** as a named `REALSXP` element of `res_list` (position determined by implementation; plan WU-E specifies the exact VECSXP index after all scalar fields are added) in `r_bridge.cpp`
(a `REALSXP` of length n, allocated with `Rf_allocVector` and `PROTECT`):

```r
calib_result$best_weights  # numeric, length n, post-normalization, obs-level
calib_result$best_error    # errRp at best_iter
calib_result$best_iter     # integer
```

Access in R (consistent across both `attach_weights` branches):
```r
attr(result, "result")$best_weights
attr(result, "result")$best_error
attr(result, "result")$best_iter
```

**bestRp tracks `errRp` (max_err metric) regardless of the active convergence
criterion.** A5 therefore checks `best_error ≤ max_error` (both are errRp),
not the active criterion's metric value. This is intentional: best-iterate is
always the iterate with minimum marginal error, regardless of stopping rule.

---

### §5 C ABI Changes Summary

New fields added to `rk_params_t` (input config):

```c
double pct_tol;         /* default 0.001 */
double absolute_tol;    /* default 0.0 */
int    criterion;       /* 0=PCT 1=MAX_ERR 2=MEAN_ERR 3=KL 4=CHI2 */
int    stop_when;       /* 0=ANY 1=ALL */
int    sor_enabled;     /* 0/1 */
int    sor_auto;        /* 0/1 */
double sor_omega_init;
double sor_omega_min;
double sor_omega_fixed; /* -1.0 = use auto */
int    sor_burnin;
```

New fields added to `rk_result_t` (output):

```c
double mean_error;      /* L1-over-margins */
double kl;              /* max KL across margins */
double chi2;            /* total chi-square */
double pct_change;      /* final iter weight pct change */
double best_error;      /* errRp at best iterate */
int    best_iter;
/* best_weights: REALSXP named element of res_list in r_bridge.cpp */
double sor_min_omega;   /* iEPPA only; zero-init ensures 0.0 default */
int    sor_n_damped;    /* iEPPA only; zero-init ensures 0 default */
```

**`EXPECTED_RK_PARAMS_BYTES` mandatory update:** Current value = 152. After
adding 10 new fields (≈80 bytes depending on alignment), measure
`sizeof(rk_params_t)` on first compile, update the constant, and document the
new layout in the comment block. The static_assert at `leafblower.h:141-142`
will fail on first compile if forgotten — this is the intended safety net.

---

### §6 Files Affected

| File | Solver scope | Change |
|---|---|---|
| `src/types.hpp` | all | Add `CalibConvergenceCfg`, `CalibSorCfg` to `CalibState` |
| `src/leafblower.h` | all | Extend `rk_params_t`, `rk_result_t`; update `EXPECTED_RK_PARAMS_BYTES` |
| `src/ieppa.cpp` | iEPPA | pct + all criteria + SOR + best-iterate |
| `src/raking.cpp` | raking | pct + all criteria + best-iterate (no SOR) |
| `src/lbfgsb_solver.cpp` | lbfgsb | pct + all criteria + best-iterate (no SOR) |
| `src/raking.hpp` | raking | Extend result struct with new scalars |
| `src/lbfgsb_solver.hpp` | lbfgsb | Extend result struct with new scalars |
| `src/c_api.cpp` | all | Marshal new config + result fields; zero-init `rk_result_t` |
| `src/r_bridge.cpp` | all | Extend `res_list` to 15 elements; add `best_weights` REALSXP |
| `R/harvest.R` | all | Rewrite `parse_convergence()`; add `sor` param; expose `$sor` diagnostics |
| `python/leafblower/_harvest.py` | all | Mirror R changes |
| `python/leafblower/_bindings.cpp` | all | Mirror C ABI changes |
| `tests/testthat/test-convergence-criteria.R` | all | New |
| `tests/testthat/test-sor.R` | iEPPA | New |
| `tests/testthat/test-best-iterate.R` | all | New |
| `man/harvest.Rd` | all | Regenerated via `devtools::document()` |
| `NEWS.md` | — | Document default criterion change, new fields, Python note |

---

### §7 Acceptance Criteria

**A1 (backward compat):** After audit of 6 affected test files (add
`convergence = list(absolute = 1e-6)` where assertions depend on criterion),
`devtools::test()` FAIL 0, PASS ≥ 232. From code inspection: `test-ieppa.R`
asserts `diag$error_weighted < 1e-6` (calibration error, not convergence
criterion) — likely safe without pinning. Each file must be individually
verified by the implementer.

**A2 (pct convergence quality):** On smooth synthetic (seed=42, n=2000, K=3,
uniform targets, max_weight=10, well-separated): harvest() with default
`convergence = list()` returns status PASS (not NOCONV), `pct_change ≤ 0.001`
at exit (confirms pct criterion fired, not incidental max_err exit), and
`max_error ≤ 1e-3`.

**A3 (SOR triggers on oscillatory):** On tight-clamp synthetic (max_weight=2,
K=5, n=5000, seed=31415 from test-homotopy.R), `calib_result$sor$min_omega < 0.9`
after 500 iters with SOR auto.

**A4 (SOR silent on smooth):** On smooth synthetic (max_weight=10), `sor$min_omega
== 1.0`. Overhead vs SOR disabled: `bench::mark(min_iterations=10)` ratio ≤ 1.05
(5% ceiling — 2% is too tight for CI noise; benchmark host adds stability).

**A5 (best-iterate always populated):** All harvest() calls: `attr(r,"result")$best_error ≤
attr(r,"result")$max_error`. Best-iterate is ≥ as good as final iterate by
definition.

**A6 (stepstone best-iterate regression):** `skip_on_cran()`. Save reference
value to `tests/testthat/fixtures/stepstone_best_error_ref.rds` during WU-gate
implementation (new fixture, generated + committed by implementer). Test: load
ref, verify `best_error ≤ ref * 1.05`.

**A7 (all 5 metrics present):** `expect_named(calib_result_names,
c("max_error","mean_error","kl","chi2","pct_change"), ignore.order=TRUE)`.

**A8 (alternative criteria actively stop the solver):** For each of
`criterion = "mean_err"`, `criterion = "kl"`, `criterion = "chi2"`, on an
appropriate smooth synthetic with a permissive threshold: harvest() returns PASS
(not NOCONV), and `pct_change` is NOT near 0.001 (proving the non-pct criterion
fired, not incidental pct exit). Three separate sub-tests.

**A9 (CRAN gate):** `devtools::test()` FAIL 0. `R CMD check --as-cran` 0 ERROR,
0 WARNING. `NEWS.md` updated.

---

### §8 Solver-Local Storage

Scratch vectors are **local to each solver function** (NOT `CalibState`). Sizes:

| Vector | iEPPA | raking | lbfgsb |
|---|---|---|---|
| `W_prev` | M_cell doubles (cell-level) | n doubles (obs-level) | n doubles (obs-level) |
| `W_best` | M_cell doubles (cell-level) | n doubles (obs-level) | n doubles (obs-level) |
| `omega[K]` | K doubles | — (no SOR) | — (no SOR) |
| `prev_errRp_k[K]` | K doubles, init ∞ | — (no SOR) | — (no SOR) |

Rationale: CalibState holds caller-owned pointers; injecting heap vectors
creates ownership ambiguity across the three solver boundaries.

**W_prev initialization:** at solver entry, from cell/obs expansion of
`start_weights` (same logic used to initialize W at solver start).

**W_best initialization:** initialized to W at the FIRST `kErrCheckInterval`
check (i.e., `best_errRp = ∞` ensures first computed errRp always triggers a
copy). If the solver exits before the first check, `W_best` is not expanded
and `best_error = ∞` is returned; callers should guard `best_error < Inf`.

**pct formula variants:**

```cpp
// iEPPA (M_cell-indexed, cell-level W):
for (int c = 0; c < M_cell; c++) {
    double rel = fabs(W[c] - W_prev[c]) / max(W_prev[c], 1e-12);
    pct_change = max(pct_change, rel);
}

// raking and lbfgsb (n-indexed, obs-level weights):
for (int i = 0; i < n; i++) {
    double rel = fabs(w[i] - W_prev[i]) / max(W_prev[i], 1e-12);
    pct_change = max(pct_change, rel);
}
```

---

### §9 WU Breakdown (per CLAUDE.md: one ticket per task)

| WU | Scope | Gate |
|---|---|---|
| **WU-A Scaffold** | Add `CalibConvergenceCfg`, `CalibSorCfg` to types.hpp; extend C ABI (rk_params_t, rk_result_t); update EXPECTED_RK_PARAMS_BYTES; rewrite parse_convergence(); add `sor` arg to harvest.R; R bridge unpacks. No criterion change yet. | Builds clean; existing tests green |
| **WU-B Convergence criteria** | pct computation + alternative criteria (mean_err, kl, chi2) in all three solvers; active criterion dispatch; A1, A2, A8 tests | FAIL 0; A1, A2, A8 pass |
| **WU-C Quality metrics** | All 5 metrics returned at exit from all three solvers; A7 test | FAIL 0; A7 pass |
| **WU-D SOR (iEPPA only)** | omega[K], prev_errRp_k[K], adaptive damping, pow guard; calib_result$sor nested; C ABI sor fields; A3, A4 tests | FAIL 0; A3, A4 pass |
| **WU-E Best-iterate** | W_best locals in all three solvers; expansion (with no-water-fill rule for iEPPA); best_weights as REALSXP element in res_list; A5, A6 tests | FAIL 0; A5, A6 pass |
| **WU-F Python + docs** | Python API mirror; NEWS.md; roxygen update; A9 (CRAN gate) | R CMD check clean; pytest green |

---

### §10 Resolved Design Questions

| Question | Resolution |
|---|---|
| pct default | 0.001; OFF (0.0) when user specifies `absolute` without `pct` (backward compat) |
| "kl" naming | Field = `kl`, criterion string = `"kl"` — consistent throughout |
| stop_when | `"any"` / `"all"` (renamed from `combine`) |
| chi2 formula | Denominator floor of 1; can return Inf on degenerate cells; documented in roxygen |
| KL ε | 1e-10; skip j where T_kj=0 |
| W_best expansion iEPPA | Scalar expansion + sum-normalization only; NO water-fill, NO bounds clamping |
| W_best expansion raking/lbfgsb | Sum-normalization only (already obs-level) |
| W_prev/W_best storage | Solver-local; M_cell for iEPPA, n for raking/lbfgsb |
| omega_min | 0.3 empirical; reproducible via benchmark sweep script |
| SOR disable | `sor = NULL` only (not `list(enabled=FALSE)`) |
| best_weights return | `calib_result$best_weights` (named REALSXP element of res_list; exact VECSXP index determined in plan WU-E) |
| sor sentinel init | c_api.cpp zero-initializes rk_result_t; non-iEPPA never writes sor fields |
| WU breakdown | §9 — 6 WUs, one beads ticket each |
| pct naming collision | Migration note added for old workaround users |
| Alternative criteria coverage | A8 tests each of mean_err, kl, chi2 as active stopping criterion |
| Python migration | Noted in Migration section and NEWS.md template |
| pct formula variants | iEPPA = M_cell loop; raking/lbfgsb = n loop (§8) |
