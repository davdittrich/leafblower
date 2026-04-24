# Convergence Criterion Reform + SOR Stabilization + Best-Iterate Design

**Date:** 2026-04-24
**Status:** Draft — revision 2, pending design-review-gate re-run
**Applies to:** all three solvers (iEPPA, raking, lbfgsb) unless noted otherwise

## Problem

iEPPA solver on stepstone-fulldata (1.58M rows, K=9, max_weight=5) diverges:
errRp increases from ~2.22e-3 at iter 50 to ~6.5e-3 at iter 3000. Root causes:

1. **Wrong stopping criterion:** `errRp` (L∞ marginal error) oscillates in the
   capacity-clamp non-smooth regime. Solver declares NOCONV but reports the
   diverged final iterate rather than the best-seen iterate.
2. **No oscillation damping:** Fixed multiplicative step size ω=1.0 in iEPPA
   amplifies oscillation.
3. **Incomplete quality reporting:** Only `max_error` returned; practitioners
   need `mean_error`, `kl_div`, `chi2` for cross-solver comparison and
   survey-statistics compliance.
4. **Criterion not portable across solvers:** raking and lbfgsb also only
   expose `max_error`; users cannot compare solvers on a common metric.

## Goals

1. **Replace `max_err` as default convergence criterion across all three solvers**
   with weight-change `pct = 0.001` — oscillation-resistant, industry-standard
   (GREG, survey::calibrate, anesrake, autumn).
2. **Allow user-selectable convergence criteria:** `pct`, `max_err`, `mean_err`,
   `kl`, `chi2`.
3. **Always report all five quality metrics** at exit across all three solvers.
4. **Auto-enable SOR adaptive under-relaxation** (iEPPA only — Thibault et al.
   2021, Soma & Uschmajew 2024) to dampen oscillation.
5. **Always track best-iterate** (W at minimum observed errRp) across all three
   solvers and expose it in `calib_result`.

## Non-goals

- Do not change the multiplicative IPF / iEPPA core algorithm structure.
- Do not add SOR to raking (convergence is theoretically guaranteed; smooth
  operators do not need it). SOR is iEPPA-only. Follow-on ticket if needed.
- Do not add SOR to lbfgsb (has its own Wolfe-condition line search).
- Do not implement Polyak-Ruppert iterate averaging — best-iterate is simpler
  and sufficient.

## Migration

`convergence['pct']` is currently a **deprecated key** in `harvest.R:244`
that warns and redirects to `absolute`. That deprecation warning must be
**removed entirely** and `parse_convergence()` rewritten to treat `pct` as the
primary key.

**Tests affected by default change:** six files call `harvest()` without an
explicit `convergence` argument — `test-ieppa.R`, `test-compare.R`,
`test-config-defaults.R`, `test-ieppa-nonuniform-d.R`, `test-compat.R`,
`test-bounded-convergence.R`. These will switch from `max_err < 1e-6` to
`pct < 0.001`. Implementer MUST audit each file: if any assertion checks
`calib_result$max_error` or iteration count without pinning convergence, add
`convergence = list(absolute = 1e-6)` to that test call to preserve
prior behavior. A1 is verified by `devtools::test()` FAIL 0 after this audit.

**CRAN:** Default-criterion change is a behavioral breaking change. Add
`NEWS.md` entry under next version bump:
> Convergence default changed from `absolute = 1e-6` (max marginal error) to
> `pct = 0.001` (max proportional weight change). To preserve prior behavior:
> `convergence = list(absolute = 1e-6)`.

---

## Design

### §1 Convergence Criterion System (all three solvers)

**New default:** `convergence = list(pct = 0.001)` — stop when
`max_i |W_i[t] - W_i[t-1]| / W_i[t-1] < 0.001` (0.1% max proportional weight
change per iteration).

**Pluggable criterion:**

```r
convergence = list(pct = 0.001)                         # new default
convergence = list(absolute = 1e-6)                     # backward compat: max_err
convergence = list(absolute = 1e-3, criterion = "chi2")
convergence = list(absolute = 1e-3, criterion = "mean_err")
convergence = list(absolute = 1e-3, criterion = "kl")
convergence = list(pct = 0.001, absolute = 1e-6, stop_when = "any")  # OR logic
```

**`stop_when`** (renamed from `combine`): relevant only when BOTH `pct` and
`absolute` are non-zero. Values: `"any"` (stop when EITHER fires, default) or
`"all"` (stop only when BOTH fire). Ignored when only one threshold is set.

**Criterion strings** (lowercase, `match.arg()`-validated in R):
`"pct"` / `"max_err"` / `"mean_err"` / `"kl"` / `"chi2"`.

**New `CalibConvergenceCfg` struct (added to `CalibState`):**

| field | type | default | meaning |
|---|---|---|---|
| `pct_tol` | double | 0.001 | weight-change threshold (0 = off) |
| `absolute_tol` | double | 0.0 | marginal error threshold (0 = off) |
| `criterion` | enum | PCT | PCT / MAX_ERR / MEAN_ERR / KL / CHI2 |
| `stop_when` | enum | ANY | ANY / ALL |

**pct computation** — runs inside inner loop, after capacity block, before
convergence check. `W_prev` is a local `std::vector<double>` in each solver
function (NOT in CalibState — see §8).

```cpp
// After capacity block, at errRp-check iterations:
double pct_change = 0.0;
for (int c = 0; c < M_cell; c++) {
    double rel = std::fabs(W[c] - W_prev[c]) / std::max(W_prev[c], 1e-12);
    if (rel > pct_change) pct_change = rel;
}
W_prev = W;   // update AFTER computing pct_change
```

**Alternative criteria formulas:**

- **max_err:** `max_k max_j |Ŝ_kj/W_total - T_kj|` — unchanged `errRp`
- **mean_err:** `mean_k(max_j |Ŝ_kj/W_total - T_kj|)` — L1-over-margins,
  L∞-within-margin
- **kl:** `max_k Σ_j T_kj * log((T_kj + ε) / (Ŝ_kj/W_total + ε))`, ε = 1e-10;
  undefined when T_kj = 0 for some j — skip those j in the sum (contribution = 0)
- **chi2:** `Σ_k Σ_j (Ŝ_kj - T_kj * W_total)² / (T_kj * W_total + 1)`;
  denominator floor of 1 prevents divide-by-zero on empty target cells.
  Note: chi2 scales with `W_total ≈ n`; users must supply n-scaled `absolute`
  threshold. Document in `@param convergence` roxygen.

All alternative criteria computed at the same `kErrCheckInterval` cadence as
`errRp`. The active criterion (per `CalibConvergenceCfg.criterion`) is compared
to the threshold; the other criteria are still computed for reporting (§2).

**Backward compatibility:** `convergence = list(absolute = X)` without
`criterion` field → `{criterion=MAX_ERR, absolute_tol=X, pct_tol=0}` — exact
current behavior. `pct_tol=0` means pct criterion is disabled; `stop_when` is
ignored.

**New `parse_convergence()` contract** (R/harvest.R rewrite):

```r
parse_convergence <- function(convergence) {
  # Returns list: pct_tol, absolute_tol, criterion, stop_when
  pct_tol      <- convergence[["pct"]]      %||% 0.001
  absolute_tol <- convergence[["absolute"]] %||% 0.0
  criterion    <- match.arg(
    convergence[["criterion"]] %||% (if (!is.null(convergence[["pct"]])) "pct" else "max_err"),
    c("pct","max_err","mean_err","kl","chi2"))
  stop_when    <- match.arg(convergence[["stop_when"]] %||% "any", c("any","all"))
  list(pct_tol=pct_tol, absolute_tol=absolute_tol, criterion=criterion, stop_when=stop_when)
}
```

The old `pct`-deprecated path is removed. `convergence[["pct"]]` is now the
primary key with no warning.

---

### §2 Quality Metric Reporting (all three solvers)

Computed at every solver's exit from the **final iterate's** W/weights. Always
present in `calib_result`:

| R field | metric | formula |
|---|---|---|
| `max_error` | L∞ marginal | unchanged `errRp` |
| `mean_error` | L1-over-margins | `mean_k(max_j |Ŝ_kj/W_total - T_kj|)` |
| `kl_div` | max KL | `max_k Σ_j T_kj log((T_kj+ε)/(Ŝ_kj/W_total+ε))` |
| `chi2` | total χ² | `Σ_k Σ_j (Ŝ_kj - T_kj·W_total)²/(T_kj·W_total + 1)` |
| `pct_change` | weight pct | last iteration's `max_i |ΔW_i/W_i|` |

**Field naming change:** `chi_square` → `chi2`, `kl_divergence` → `kl_div`
in the C ABI `rk_result_t` and R bridge. Matches criterion string names.

`DEFF` and `ESS` remain R-side helpers requiring the returned weights.

---

### §3 SOR Adaptive Under-Relaxation (iEPPA only)

**Theory:** Thibault–Chizat–Dossal–Papadakis (HAL 2017, Algorithms 14(5):143
2021); Soma–Uschmajew arxiv:2410.14104 (2024). Replace iEPPA margin step:

```cpp
// Before:
f_k_new[j] = f_k_old[j] * (T_kj / S_kj);

// After:
const double ratio = T_kj / S_kj;
// Guard: empty bucket makes S_kj → 0; clamp ratio using existing threshold
if (S_kj < kEmptyBucketThreshold * ct.W_input) {
    // skip update — bucket already guarded by existing empty-bucket logic
} else {
    const double r = (omega[k] == 1.0) ? ratio
                                       : std::pow(ratio, omega[k]);
    f_k_new[j] = f_k_old[j] * r;
}
```

The `kEmptyBucketThreshold` check is already present in the linear-path sweep;
the SOR guard reuses it. Using `omega[k] == 1.0` branch avoids `std::pow()`
on the non-oscillating common path — not all R-target compilers (LLVM/Windows)
elide `pow(x, 1.0)` at `-O2`.

**Default: `sor = list(auto = TRUE, omega_min = 0.3)`**

Auto mode per margin k (local to `ieppa_solve()`, K-element arrays):

- `omega[K]` initialized to `omega_init = 1.0`
- `prev_errRp_k[K]` initialized to ∞ (K-element local array, required for sign-flip detection)
- Burn-in: first `burnin = 20` outer iterations, all ω unchanged
- After burn-in, at each `kErrCheckInterval` check:
  - Compute per-margin errRp_k = `max_j |Ŝ_kj/W_total - T_kj|`
  - Sign flip (errRp_k[t] > errRp_k[t-1]): `omega[k] = max(omega_min, omega[k] * 0.7)`
  - Monotone decrease: `omega[k] = min(1.0, omega[k] * 1.05)` (slow recovery)
  - Store `prev_errRp_k[k] = errRp_k`

**omega_min = 0.3 rationale:** Conservative empirical lower bound. Pilot runs
on stepstone-fulldata showed ω < 0.3 over-dampens non-oscillating margins,
slowing convergence by >2×. ω = 0.3 stabilizes oscillation within ~100 iters
on the benchmark dataset. SOR theory (Young 1954, Thibault 2021) gives optimal
ω for smooth linear problems; for non-linear clamped iEPPA no analytic formula
exists — 0.3 is the safe empirical floor from pilot testing.

**R API:**

```r
sor = list(auto = TRUE, omega_min = 0.3)    # default
sor = list(omega = 0.5)                     # fixed uniform ω, no per-margin adaptation
sor = NULL                                  # disable (canonical off state)
```

`list(enabled = FALSE)` is NOT supported — use `sor = NULL`.

**`CalibSorCfg` struct (added to `CalibState`):**

| field | type | default |
|---|---|---|
| `enabled` | bool | true |
| `auto_adapt` | bool | true |
| `omega_init` | double | 1.0 |
| `omega_min` | double | 0.3 |
| `omega_fixed` | double | -1 (sentinel: use auto) |
| `burnin` | int | 20 |

**SOR diagnostics in `calib_result$sor`** (nested list):

```r
calib_result$sor$min_omega   # lowest ω across any margin × any iter
calib_result$sor$n_damped    # count of (margin, iter) pairs where ω < 1
```

Nested under `$sor` to mirror the `sor = list(...)` input structure and avoid
polluting the top-level result namespace.

**C ABI additions** to `rk_result_t`:

```c
double sor_min_omega;   /* iEPPA only; 1.0 if SOR inactive */
int    sor_n_damped;    /* iEPPA only; 0 if SOR inactive */
```

---

### §4 Best-Iterate Tracking (all three solvers)

Always enabled. Zero toggle overhead. Storage: O(M_cell) doubles copied only
when errRp improves (monotone copy-on-improvement).

**Implementation in each solver** (locals, NOT CalibState):

```cpp
// Local to solver function:
double best_errRp = std::numeric_limits<double>::infinity();
int    best_iter  = 0;
std::vector<double> W_best(M_cell, 0.0);  // initialized at first check

// After errRp computed:
if (errRp < best_errRp) {
    best_errRp = errRp;
    best_iter  = iter;
    W_best     = W;           // copy M_cell doubles (cell-level weights, PRE-normalization)
}

// At exit:
// Expand W_best → obs-level using same formula as W_final expansion.
// Apply the same post-loop normalization (sum(w) = n) to W_best.
// best_weights (obs-level) is therefore POST-normalization, consistent with w_final.
```

**W_best storage clarification:** stored as cell-level W[] (capacity multipliers,
pre-loop normalization). Expansion + normalization at exit makes `best_weights`
(obs-level) post-normalization — identical treatment as `weights` (final iterate).

**Return via `calib_result`** — NOT via `attr()`. Promoted to a list element
inside the named list returned by `.Call`:

```r
calib_result$best_weights  # numeric vector, length n, obs-level, post-normalization
calib_result$best_error    # errRp at best_iter
calib_result$best_iter     # integer
```

Both `attach_weights = TRUE` and `attach_weights = FALSE` branches access
`best_weights` via `attr(result, "result")$best_weights` (consistent with
existing `attr(result, "result")$max_error` pattern). No separate `attr()`
entry for `best_weights`.

**Memory:** M_cell << n for compressed cell-table inputs. For stepstone-fulldata:
M_cell is the count of unique demographic cell combinations; much less than
1.58M. The spec mandates that implementer confirms M_cell size from
`build_cell_table()` output before finalizing. If M_cell > 500k, escalate
to a separate HUMAN ticket.

---

### §5 C ABI Changes Summary

New fields added to `rk_params_t` (input):

```c
/* Convergence config */
double pct_tol;         /* default 0.001 */
double absolute_tol;    /* default 0.0 */
int    criterion;       /* 0=PCT 1=MAX_ERR 2=MEAN_ERR 3=KL 4=CHI2 */
int    stop_when;       /* 0=ANY 1=ALL */

/* SOR config (iEPPA only; ignored by raking/lbfgsb) */
int    sor_enabled;
int    sor_auto;
double sor_omega_init;
double sor_omega_min;
double sor_omega_fixed; /* -1.0 = use auto */
int    sor_burnin;
```

New fields added to `rk_result_t` (output):

```c
double mean_error;      /* L1-over-margins */
double kl_div;          /* max KL across margins */
double chi2;            /* total chi-square */
double pct_change;      /* final iter weight pct change */
double best_error;      /* errRp at best iterate */
int    best_iter;
/* best_weights: returned as a separate SEXP in the R bridge result list */
double sor_min_omega;   /* iEPPA only; 1.0 otherwise */
int    sor_n_damped;    /* iEPPA only; 0 otherwise */
```

**`EXPECTED_RK_PARAMS_BYTES` update (mandatory):** After adding new fields,
measure `sizeof(rk_params_t)` and update the `static_assert` constant in
`src/leafblower.h`. Document new field layout in the layout comment block.

---

### §6 Files Affected

| File | Change |
|---|---|
| `src/types.hpp` | Add `CalibConvergenceCfg`, `CalibSorCfg` to `CalibState` |
| `src/leafblower.h` | Extend `rk_params_t`, `rk_result_t`; update `EXPECTED_RK_PARAMS_BYTES` |
| `src/ieppa.cpp` | pct + alternative criteria + SOR + best-iterate |
| `src/raking.cpp` | pct + alternative criteria + best-iterate (NO SOR) |
| `src/lbfgsb_solver.cpp` | pct + alternative criteria + best-iterate (NO SOR) |
| `src/raking.hpp` | Extend result struct with new scalar fields |
| `src/lbfgsb_solver.hpp` | Extend result struct with new scalar fields |
| `src/c_api.cpp` | Marshal all new config + result fields |
| `src/r_bridge.cpp` | Unpack new args; build extended result list incl. `best_weights` SEXP |
| `R/harvest.R` | Rewrite `parse_convergence()`; add `sor` arg; expose `$sor` diagnostics |
| `python/leafblower/_harvest.py` | Mirror R API changes |
| `python/leafblower/_bindings.cpp` | Mirror C ABI changes |
| `tests/testthat/test-convergence-criteria.R` | New tests for pct, alternative criteria, backward compat |
| `tests/testthat/test-sor.R` | New tests for SOR auto trigger / no-trigger |
| `tests/testthat/test-best-iterate.R` | New tests for best_weights, best_error, best_iter |
| `man/harvest.Rd` | Regenerated via `devtools::document()` |
| `NEWS.md` | Document default criterion change |

---

### §7 Acceptance Criteria

**A1 (backward compat):** `convergence = list(absolute = X)` (no criterion
field) = exact current behavior. All 6 test files that omit convergence are
audited and pinned; `devtools::test()` FAIL 0 at 232+ PASS after audit.

**A2 (pct convergence quality):** On a smooth synthetic (uniform targets
achievable within max_weight, no binding capacity constraints, seed=42,
n=2000, K=3), `pct = 0.001` converges (reports PASS not NOCONV) and the
final `max_error ≤ 1e-3`. No iteration-count parity claim (pct and max_err
measure different quantities; no analytic equivalence exists).

**A3 (SOR triggers on oscillatory):** On the tight-clamp synthetic
(max_weight=2, K=5, n=5000, seed=31415), `calib_result$sor$min_omega < 0.9`
after 500 iterations with SOR auto.

**A4 (SOR silent on smooth):** On smooth synthetic (max_weight=10,
well-separated targets), `calib_result$sor$min_omega == 1.0` (no damping
triggered). `devtools::bench::mark()` shows ≤ 2% overhead vs SOR disabled.

**A5 (best-iterate always populated):** In all harvest() calls across all
three methods, `attr(result,"result")$best_error ≤ attr(result,"result")$max_error`.

**A6 (stepstone best-iterate regression):** `skip_on_cran()`. Load
`tests/testthat/fixtures/stepstone_best_error_ref.rds` (pre-saved reference:
2.22e-3); verify `best_error ≤ ref * 1.05` (within 5% of reference).

**A7 (all 5 metrics present):** `expect_named(names(attr(result,"result")),
c("max_error","mean_error","kl_div","chi2","pct_change"), ignore.order=TRUE)`.

**A8 (CRAN gate):** `devtools::test()` FAIL 0. `R CMD check --as-cran` 0 ERROR,
0 WARNING. `NEWS.md` updated.

---

### §8 Solver-Local Storage (NOT CalibState)

The following scratch vectors are **local to each solver function**, NOT members
of CalibState:

- `W_prev` (M_cell doubles): initialized to start_weights cell expansion at
  solver entry. Updated after each pct check.
- `W_best` (M_cell doubles): initialized to W at first errRp check. Updated
  on improvement.
- `omega[K]` (K doubles, iEPPA only): initialized to `sor.omega_init`.
- `prev_errRp_k[K]` (K doubles, iEPPA SOR only): initialized to ∞. Tracks
  per-margin errRp at last check for sign-flip detection.

Rationale: CalibState holds caller-owned pointers with caller-controlled lifetimes.
Heap vectors with solver-local lifetimes belong in solver-local scope. raking and
lbfgsb allocate their own W_prev / W_best without touching CalibState.

---

### §9 Resolved Design Questions

| Question | Resolution |
|---|---|
| pct default value | 0.001 (0.1% max weight change). Tighter than anesrake (0.01); matches precision-calibration practice |
| chi2 denominator | T_kj·W_total + 1 (floor of 1; unusual choice documented in roxygen — comparable to Yates-corrected Pearson χ²) |
| KL ε | 1e-10; skips j where T_kj = 0 (contribution = 0 by convention) |
| W_best storage | Cell-level (M_cell); expanded + normalized at exit; post-normalization obs-level result |
| omega_min | 0.3 — empirical conservative floor; no analytic SOR theorem applies to non-linear clamped operator |
| SOR on raking/lbfgsb | Not in this spec. Follow-on ticket if needed. |
| best_weights return path | `calib_result$best_weights` (list element) not `attr()` — consistent across both attach_weights branches |
| criterion string case | Lowercase `"kl"` / `"chi2"` / `"max_err"` / `"mean_err"` / `"pct"` — R convention, `match.arg()` validated |
| Metric field names | `kl_div` / `chi2` — match criterion string names; replace prior `kl_divergence` / `chi_square` |
| stop_when | Renamed from `combine`; values `"any"` / `"all"` |
