# Convergence Criteria Redesign — Spec

**Date:** 2026-04-25
**Status:** Draft v2 — post design-review-gate revision (all 5 reviewer blockers addressed)

## Problem

Current convergence design conflates two orthogonal concerns:

1. **What to measure** (loss function / metric)
2. **How to stop** (stopping rule)

The current `criterion` enum mixes both: `PCT` = L∞ weight change + threshold rule;
`MAX_ERR` = L∞ marginal error + threshold rule; etc. This prevents users from choosing,
e.g., "stop when KL-divergence stops improving" or "stop when L1 weight change plateaus."

Additionally, the current `pct` definition (max relative weight change per obs,
`max_i |Δw_i/w_i|`) does not match any industry standard. The standard (autumn,
anesrake) is L1 total weight change sum `Σ|Δw_i|`.

## Research Grounding

### Literature
- **Csiszár (1975)** / **Ireland-Kullback (1968)**: KL-divergence is monotone decreasing
  under IPF for feasible problems. Theoretical convergence criterion = `D(p*||p^t) → 0`.
  No explicit threshold given.
- **Deville-Särndal (1992)**: Calibration convergence measured by constraint residual
  `||Σ w_k x_k - X||`. Newton steps converge in 2–4 iterations for well-posed problems.
- **Gietl-Fröhlich (2013)**: When exact convergence is impossible (bounded weights,
  structural zeros), IPF accumulates at TWO limit points — oscillation is structural.
  Neither limit point satisfies all margins. This mandates best-iterate tracking.

### Software implementations
| Package | Metric | Rule | Default tol |
|---|---|---|---|
| survey::grake | `|misfit_k|/(1+|pop_k|)` | threshold | 1e-7 |
| ipfp | `‖Ax−y‖₂` | threshold | 1.49e-8 |
| mipfp | `max|x^t − x^{t-1}|` | threshold | 1e-10 |
| ipfn | `max_k |ΣS_k/T_k − 1|` | improvement+threshold | 1e-5 + 1e-8 |
| icarus | `‖Ax−b‖` | threshold | 1e-6 |
| anesrake | `Σ|Δw|` | improvement | pct=0.01 |
| autumn | `Σ|Δw|` | both | user-set |

**Key gap:** No package uses "rate of change of loss function" as stopping criterion —
but `ipfn`'s `rate_tolerance` is the closest analog.

## Design

### §1 Orthogonal Decomposition

Replace single `criterion` enum with two independent fields:

```c
/* Metric: what quantity is evaluated */
int metric;  /* see CalibMetric enum */

/* Rule: how to stop */
int rule;    /* see CalibRule enum */
```

#### §1.1 Metrics (CalibMetric)

| Value | Name | Formula | Grounding |
|---|---|---|---|
| 0 | `max_err` | `max_k max_j |Ŝ_kj/W - T_kj|` (errRp) | Deville-Sarndal constraint residual |
| 1 | `mean_err` | `mean_k(max_j |Ŝ_kj/W - T_kj|)` | Mean constraint satisfaction |
| 2 | `kl` | `max_k Σ_j T_kj log((T_kj+ε)/(Ŝ_kj/W+ε))` | Csiszár I-divergence |
| 3 | `chi2` | `Σ_k Σ_j (Ŝ_kj - T_kj·W)²/(T_kj·W+1)` | Chi-squared misfit |
| 4 | `grake_norm` | `max_k |misfit_k|/(1+|pop_k|)` | survey::grake normalized residual |
| 5 | `l1_weight` | `Σ_i |w_i[t] - w_i[t-1]| / n` | autumn/anesrake industry standard |

**`l1_weight`:** normalized by n to make scale-free (same interpretation regardless of
sample size). Equivalent to autumn's `Σ|Δw|` divided by n.

#### §1.2 Stopping Rules (CalibRule)

| Value | Name | Condition | Grounding |
|---|---|---|---|
| 0 | `threshold` | `metric_t < tol` | Absolute threshold on loss |
| 1 | `improvement` | `|metric_t - metric_{t-1}| / metric_{t-1} < tol` | Rate of change of loss; ipfn `rate_tolerance` |
| 2 | `plateau` | `metric_t > metric_{t-1} * (1 - tol)` | Autumn/anesrake pct criterion |

**`improvement`:** Absolute relative change. Fires when loss barely changes in EITHER
direction (including oscillation once amplitude < tol). Conservative — does not stop
immediately on divergence; waits for oscillation to stabilize.

**`plateau`:** Fires when metric stops IMPROVING (autumn semantics). Fires immediately
if metric increases. More aggressive than `improvement`.

**Shorthand keys** (software not yet released — clean break, no backward compat shims needed):

| Shorthand key | Equivalent |
|---|---|
| `list(pct = X)` | `metric="l1_weight", rule="plateau", tol=X` |
| `list(absolute = X)` | `metric="max_err", rule="threshold", tol=X` |
| `list(improvement = X)` | `metric="max_err", rule="improvement", tol=X` |

### §2 Default Convergence

```r
convergence = list(metric = "max_err", rule = "improvement", tol = 0.001)
```

**Rationale:**
- `max_err` = constraint residual = direct measure of calibration quality
  (Deville-Sarndal). Not a proxy.
- `improvement` = absolute relative change of constraint residual between checks.
  Fires when residual stops changing meaningfully (plateau OR oscillation stabilization).
  Conservative: oscillating problems run to max_iter, but best_iterate tracking returns
  the best-seen weights.
- `tol = 0.001` (0.1% relative change): for errRp ≈ 1e-3, 0.1% corresponds to
  ~1e-6 absolute — comparable to survey::grake's default `epsilon=1e-7`.

**Active criterion in result:** `attr(r,"result")$convergence_used` exposes what fired:
```r
list(metric="max_err", rule="improvement", tol=0.001, fired_at_iter=47)
```
Users can inspect this to understand why a run converged or ran to max_iter.

### §3 R API

```r
# New unified long-form API
convergence = list(
  metric = "max_err",      # "max_err"|"mean_err"|"kl"|"chi2"|"grake_norm"|"l1_weight"
  rule   = "improvement",  # "threshold"|"improvement"|"plateau"
  tol    = 0.001
)

# Shorthand keys (parsed to long-form; take precedence over metric/rule if both present)
convergence = list(improvement = 0.001)    # → max_err+improvement, tol=0.001 [NEW DEFAULT]
convergence = list(pct = 0.001)            # → l1_weight+plateau, tol=0.001 [RE-DEFINED; warns]
convergence = list(absolute = 1e-6)        # → max_err+threshold, tol=1e-6 [UNCHANGED]

# stop_when: combines two independent criteria (each with own metric/tol)
# When stop_when is present, a second threshold is always max_err+threshold:
# convergence = list(improvement = 0.001, absolute = 1e-3, stop_when = "any")
# → fires when EITHER (max_err+improvement < 0.001) OR (max_err+threshold < 1e-3)
```

**`stop_when` semantics with new API:**
`stop_when` is retained for the specific case of combining an improvement/plateau rule
with a hard absolute fallback on `max_err`. When `stop_when` is present:
- Primary criterion: the shorthand or long-form `metric+rule+tol`
- Secondary criterion (always): `metric=max_err, rule=threshold, tol=absolute`
- `stop_when="any"` (default): fire when EITHER criterion is met
- `stop_when="all"`: fire when BOTH criteria are met

This is the only supported multi-criterion combination. Specifying `stop_when` without
`absolute` is an error.

**`parse_convergence()` contract (complete):**

Valid keys: `c("metric", "rule", "tol", "pct", "absolute", "improvement", "stop_when")`

Unknown keys → `stop()` with "Unknown convergence key(s): X. Valid keys: ..."

`metric` validation: `match.arg(metric, c("max_err","mean_err","kl","chi2","grake_norm","l1_weight"))`
`rule` validation: `match.arg(rule, c("threshold","improvement","plateau"))`

**`tol` interpretation by rule:**
- `threshold`: `metric_t < tol`
- `improvement`: `|metric_t - metric_{t-1}| / metric_{t-1} < tol` (absolute relative change)
- `plateau`: `metric_t > metric_{t-1} * (1 - tol)` (fire when improvement stalls)

**`tol` range constraint:** `tol` must be in `(0, 1)` for `plateau` rule (tol ≥ 1 makes
`1 - tol ≤ 0` and condition degenerates). Validated in `parse_convergence()`.

**grake_norm user explanation (for roxygen):** "survey::grake-compatible normalized margin
residual: `max_k |misfit_k| / (1 + |population_k|)`. Use `tol=1e-7` for convergence
equivalent to `survey::calibrate(epsilon=1e-7)`."

### §4 C ABI Changes

**Struct change:** `criterion` (int) removed; `metric` (int) + `rule` (int) added.
Net change: remove 1 int, add 2 ints = **+4 bytes**.

```c
/* BEFORE: criterion (4B) at offset ~208 */
/* AFTER: metric (4B) + rule (4B) at same region */
int metric;      /* CalibMetric: 0=max_err 1=mean_err 2=kl 3=chi2 4=grake_norm 5=l1_weight */
int rule;        /* CalibRule: 0=threshold 1=improvement 2=plateau */
/* tol_abs / pct_tol fields retained; stop_when retained */
```

**`EXPECTED_RK_PARAMS_BYTES` = 216 + 4 = 220** (verified on Linux x86_64; static_assert
in `leafblower.h` enforces this at compile time).

**Validation guards (explicit):**
```c
if (p->metric < 0 || p->metric > 5)
    return err("metric out of range [0,5]: 0=max_err 1=mean_err 2=kl 3=chi2 4=grake_norm 5=l1_weight");
if (p->rule < 0 || p->rule > 2)
    return err("rule out of range [0,2]: 0=threshold 1=improvement 2=plateau");
```

**Backward compat helper** (declared in `leafblower.h`, implemented in `c_api.cpp`):
```c
/* Map old criterion int [0,4] to new metric+rule. Returns RK_OK or RK_ERR_BADARG. */
LBW_API int rk_params_set_legacy_convergence(rk_params_t* p, int old_criterion);
```
Mapping: 0=PCT→metric=5,rule=2; 1=MAX_ERR→metric=0,rule=0; 2=MEAN_ERR→metric=1,rule=0;
3=KL→metric=2,rule=0; 4=CHI2→metric=3,rule=0. Old values 5+ → `RK_ERR_BADARG`.

**rk_result_t field rename** (clean break — no alias needed):
```c
double l1_weight_change;  /* renamed from pct_change; L1 normalized Σ|Δw|/n */
```
`pct_change` removed entirely from `rk_result_t` and from R result list.
`attr(r,"result")$l1_weight_change` replaces `attr(r,"result")$pct_change`.

**rk_result_t new output field for `convergence_used`:**
```c
int    convergence_metric;  /* CalibMetric used at convergence */
int    convergence_rule;    /* CalibRule used at convergence */
double convergence_tol;     /* threshold that fired */
int    convergence_iter;    /* iteration at which convergence fired (-1 if max_iter) */
```
These populate `attr(r,"result")$convergence_used` in R.

### §5 Solver Implementation

All three solvers (iEPPA, raking, lbfgsb) compute all 6 metrics at every
`kErrCheckInterval`. The active stopping condition uses `metric` + `rule` from
`CalibState.convergence_cfg`.

#### §5.1 `l1_weight` metric (renamed from pct_change)

**iEPPA:** cell-level mass change, normalized by `ct.W_input` (= Σ X_init[c] = n for
unit-sum weights):
```cpp
double l1w = 0.0;
for (int c = 0; c < ct.M_cell; c++) l1w += std::fabs(X[c] - X_prev[c]);
if (ct.W_input > 0.0) l1w /= ct.W_input;  // normalize; ct.W_input > 0 guaranteed by validation
```

**raking/lbfgsb:** obs-level, normalized by `st.n`:
```cpp
double l1w = 0.0;
for (int i = 0; i < st.n; i++) l1w += std::fabs(w[i] - w_prev[i]);
l1w /= static_cast<double>(st.n);
```

**Existing `X_prev`/`w_prev` tracking reused** (already declared in WU-B infrastructure).

#### §5.2 `grake_norm` metric

Per-margin iteration over categories; requires `W_total` (already computed at this block):
```cpp
double gn = 0.0;
for (int k = 0; k < st.K; k++) {
    for (int j = 0; j < st.cat_counts[k]; j++) {
        double pop_kj = st.targets[k][j] * W_total;
        double misfit = S_kj[j] - pop_kj;  // S_kj already in margin accumulation buffer
        double norm_misfit = std::fabs(misfit) / (1.0 + std::fabs(pop_kj));
        if (norm_misfit > gn) gn = norm_misfit;
    }
}
```

`user-benefit:` Enables direct comparison with `survey::calibrate(epsilon=X)` — same
stopping criterion, same iteration count. Validated by A6.

#### §5.3 `improvement` rule (new tracking)

`prev_metric` is solver-local, declared BEFORE the outer (homotopy) loop, initialized
to `+∞`. Reset at each homotopy level boundary (same as `X_prev` reset pattern):
```cpp
double prev_metric_for_improvement = std::numeric_limits<double>::infinity();
// ... at homotopy level boundary:
prev_metric_for_improvement = std::numeric_limits<double>::infinity();
```

Inside `kErrCheckInterval` block:
```cpp
double rel_change = 1.0;  // default: no convergence on first check (prev=∞)
if (std::isfinite(prev_metric_for_improvement) &&
    prev_metric_for_improvement > 1e-15) {  // guard: avoid explosion near zero
    rel_change = std::fabs(curr_metric - prev_metric_for_improvement)
                 / prev_metric_for_improvement;
}
// Check isfinite(curr_metric) before convergence:
if (std::isfinite(curr_metric)) {
    prev_metric_for_improvement = curr_metric;
    converged = (rel_change < cfg.tol);
}
```

#### §5.4 `plateau` rule (autumn-equivalent)

```cpp
bool improved_enough = (curr_metric < prev_metric * (1.0 - cfg.tol));
// prev_metric initialized to +∞; first check: ∞*(1-tol) = ∞ → improved_enough=true
// → converged=false. Correct.
converged = !improved_enough;
prev_metric = curr_metric;
```

`cfg.tol` validated in `(0,1)` by R-side `parse_convergence()` before reaching C.

Note: on first check (prev_metric = ∞), `curr_metric < ∞ * (1 - tol)` is always true →
improved_enough=true → converged=false. Correct: never converge on first check.

### §6 Breaking Changes (pre-release — no user migration needed)

Software has not been released. All breaking changes are clean cuts:

1. `criterion` C ABI field removed; replaced by `metric` + `rule`.
2. `pct_change` result field renamed to `l1_weight_change` everywhere (C, R, Python).
3. Default convergence: `list()` → `max_err+improvement+0.001` (was max_err threshold).
4. `pct` key re-defined: now `l1_weight+plateau` (was max-per-obs relative change).
5. All existing tests referencing `pct_change` or old `criterion` strings updated in WU-F.

No shims, no aliases, no deprecation warnings required.

### §7 Acceptance Criteria

**A1 (default converges smooth):** `harvest()` with no `convergence` arg on well-posed
synthetic (K=3, n=2000, max_weight=10): status=OK, `max_error < 1e-3`,
`attr(r,"result")$convergence_used$rule == "improvement"`.

**A2 (default handles oscillation):** On tight-clamp synthetic (max_weight=2, K=5,
n=2000, seed=31415): `attr(r,"result")$best_error ≤ attr(r,"result")$max_error`.
If NOCONV: `best_error < 0.9 * max_error` (best iterate meaningfully better than final).

**A3 (pct = autumn L1 improvement semantics):** `list(pct=0.001)` on stepstone-small:
verify `attr(r,"result")$l1_weight_change < 0.001` at exit (pct criterion fired).

**A4 (all 6 metrics present):** `attr(r,"result")` contains `max_error`, `mean_error`,
`kl`, `chi2`, `grake_norm`, `l1_weight_change` (renamed from `pct_change`).

**A5 (absolute shorthand):** `list(absolute=1e-6)` → `max_err+threshold`; status=OK
on smooth input; `max_error < 1e-6` at exit.

**A6 (grake_norm matches survey::grake):** `skip_if_not_installed("survey")`. On smooth
synthetic (K=2, n=1000): leafblower `metric="grake_norm", rule="threshold", tol=1e-7`
converges within ±2 iterations of `survey::calibrate(epsilon=1e-7)`.

**A7 (CRAN gate):** `devtools::test()` FAIL 0; `R CMD check --as-cran` 0 ERROR, 0 WARNING.

### §8 WU Breakdown

**Pre-step (before WU-A):** Pin all existing `harvest()` test calls without explicit
convergence to `convergence=list(absolute=1e-6)` so the test suite stays green during
the scaffold window (WU-A removes `criterion` field before WU-C implements new dispatch).

| WU | Scope | Files | Gate |
|---|---|---|---|
| **WU-A** | Scaffold: CalibMetric + CalibRule enums in types.hpp; CalibConvergenceCfg update; leafblower.h: replace `criterion` with `metric`+`rule`; `EXPECTED_RK_PARAMS_BYTES = 220`; `rk_params_set_legacy_convergence()` declared in leafblower.h + implemented in c_api.cpp; **r_bridge.cpp: replace `p.criterion` with `p.metric`+`p.rule`**; c_api.cpp: update `convergence_cfg.criterion` → `.metric`/`.rule` assignments | types.hpp, leafblower.h, c_api.cpp, r_bridge.cpp | builds; static_assert passes at 220; existing pinned tests green |
| **WU-B** | c_api.cpp: `rk_params_init` defaults (`metric=0, rule=1, tol=0.001`); validation.hpp: update guards to `metric∈[0,5]` and `rule∈[0,2]`; `rk_result_t`: rename `pct_change` → `l1_weight_change`; add `convergence_metric/rule/tol/iter` fields; c_api.cpp: `rk_result_init` update | c_api.cpp, validation.hpp, leafblower.h | builds |
| **WU-C** | ieppa.cpp: compute all 6 metrics at kErrCheckInterval; `improvement` + `plateau` dispatch; `prev_metric` solver-local with homotopy reset | ieppa.cpp | A1, A2 pass |
| **WU-D** | raking.cpp + lbfgsb_solver.cpp: mirror WU-C | raking.cpp, lbfgsb_solver.cpp | A3 pass |
| **WU-E** | r_bridge.cpp: unpack new result fields (`l1_weight_change`, `grake_norm`, `convergence_used`); R/harvest.R: `parse_convergence()` rewrite with `metric+rule+tol` API; remove `criterion` key; `pct` now → l1_weight+plateau; `criterion_int` → `metric_int`/`rule_int` maps | r_bridge.cpp, R/harvest.R | A1–A5 pass |
| **WU-F** | Tests: one ticket for new A1–A6 tests; one ticket for test migration (update pinned tests to use explicit `convergence=list(metric=..., rule=..., tol=...)` or remove pins where they can run with new default); add `survey` to DESCRIPTION Suggests; wrap A6 in `skip_if_not_installed("survey")` | tests/testthat/*.R, DESCRIPTION | FAIL 0; A6 pass |
| **WU-G** | python/_harvest.py: update ctypes Structure (`criterion` → `metric`+`rule` fields); update `_parse_convergence()`; update result field mapping (`pct_change` → `l1_weight_change`); NEWS.md; roxygen regen | python/leafblower/_harvest.py, NEWS.md, man/ | pytest green; A7 pass |

### §9 Resolved Design Questions

| Question | Resolution |
|---|---|
| pct re-definition | `Σ|Δw|/n` + plateau rule = autumn/anesrake L1 standard (normalized for scale-free) |
| improvement rule | Absolute relative change `|metric_t - metric_{t-1}| / metric_{t-1} < tol` |
| default | `metric=max_err, rule=improvement, tol=0.001` — theoretically grounded, handles oscillation via best_iterate |
| grake_norm | Metric 4; enables survey::grake-equivalent stopping; validated by A6 |
| l1_weight user benefit | Auth for `pct` key; matches industry-standard raking package behavior |
| C ABI compat | `rk_params_set_legacy_convergence()` for any future C API users; not required for current R/Python paths |
| pct_change rename | Clean break: `pct_change` → `l1_weight_change` everywhere; no alias (unreleased software) |
| oscillation | Improvement rule + best_iterate: conservative (runs to max_iter on oscillation, returns best seen) |
| broken build window | WU-A updates r_bridge.cpp + c_api.cpp criterion refs atomically; pre-step pins tests |
| survey in DESCRIPTION | Added to Suggests in WU-F; A6 wrapped in skip_if_not_installed |
| Python ctypes | WU-G explicitly updates ctypes Structure to add metric+rule, remove criterion |
| stop_when semantics | Combines primary criterion with secondary `max_err+threshold+absolute`; only valid when `absolute` also set |
| tol validation | `(0,1)` enforced by parse_convergence() for `plateau` rule; no C-side enforcement needed |
