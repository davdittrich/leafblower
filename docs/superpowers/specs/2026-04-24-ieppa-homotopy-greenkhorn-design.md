# Design Spec: iEPPA Homotopy + Greenkhorn Priority Scheduler + Tang Dynamic-η

**Date:** 2026-04-24  
**Revision:** 1  
**Status:** Approved (plan rev 2)  
**Author:** WU-1 scaffolding  
**Literature backing:** `docs/investigations/2026-04-24-ieppa-accel-research.md`

---

## §1 Motivation

The stepstone-fulldata benchmark shows leafblower iEPPA at max_err ≈ 1.60e-3 after 500 iterations
vs `autumn::harvest(accelerate=TRUE)` reaching the same error in ~34 s wall time. Three
iEPPA-preserving overlays compose to close this gap without altering the core algBCD invariants:

- **P-A Homotopy:** progressively tighten `max_weight` from `start_factor` to `end_factor` across
  `n_levels` levels, interpolating geometrically. Each level warm-starts from the previous solution.
  Inspired by Chizat et al. 2018 (√t-annealing for unbalanced OT) and Section 3 of Tang 2024
  (arXiv:2403.05054) which shows constraint relaxation + tightening yields faster convergence for
  constrained Sinkhorn.
- **P-B Greenkhorn Priority Scheduler:** replace the fixed round-robin margin sweep with an argmax
  selector that sweeps the margin with the highest residual first, capped at K steps per pass.
  Greedy scheduling was shown to accelerate Sinkhorn by up to 4× on structured matrices (Altschuler
  et al. 2017; Peyré & Cuturi 2019 textbook §4.4). Internal `residual_recheck_fraction = 0.1`
  controls how often the priority queue is refreshed — this is NOT a user-visible parameter.
- **Tang-η Dynamic Schedule:** between homotopy levels, multiply the ALM penalty `alm_mu` by a
  factor computed as `eta_start * (l / (n_levels - 1))^schedule_power` at level `l`. When
  `n_levels = 1` (homotopy disabled) the schedule degenerates to the constant `eta_start = 1.0`
  (identity). Tang 2024 §3.2 demonstrates that a decaying penalty schedule speeds dual convergence
  for constrained optimal transport.

All three overlays are default-off (identity). Existing behaviour is preserved exactly.

**ylsy-memory citation:** decision to compose three overlays (not replace the solver) recorded in
session memory 2026-04-24: "iEPPA-preserving overlays only; no algorithmic replacement."

---

## §2 Config Contract

### C++ (internal, `src/types.hpp`)

```cpp
struct HomotopyConfigLbw {
    int    n_levels        = 1;     // 1 = disabled (single level = current behaviour)
    double start_factor    = 1.0;   // starting max_weight multiplier
    double end_factor      = 1.0;   // ending max_weight multiplier
    double budget_split_p  = 0.5;   // Chizat-inspired budget split (0.5 = equal)
    bool   enabled         = false; // master toggle; set true when n_levels > 1
};

enum class SchedulerMode   : int { ROUND_ROBIN = 0, GREEDY = 1 };
enum class EtaScheduleMode : int { FIXED = 0, TANG_DYNAMIC = 1 };

struct SchedulerConfigLbw {
    SchedulerMode mode                    = SchedulerMode::ROUND_ROBIN;
    double        residual_recheck_fraction = 0.1;  // INTERNAL; not exposed to ABI or R
};

struct EtaScheduleConfigLbw {
    EtaScheduleMode mode          = EtaScheduleMode::FIXED;
    double          eta_start     = 1.0;
    double          eta_end       = 1.0;
    double          schedule_power = 0.5;
};
```

These three structs are members of `CalibState` (not exposed through the public C ABI directly —
they are filled by `c_api.cpp` from the new `rk_params_t` fields).

### C ABI (`src/leafblower.h`)

```c
typedef enum { RK_SCHED_ROUND_ROBIN = 0, RK_SCHED_GREEDY = 1 } rk_scheduler_t;
typedef enum { RK_ETA_FIXED = 0, RK_ETA_TANG_DYNAMIC = 1 } rk_eta_mode_t;

typedef struct {
    int    n_levels;       /* default 1 */
    double start_factor;   /* default 1.0 */
    double end_factor;     /* default 1.0 */
    double budget_split_p; /* default 0.5 */
    int    enabled;        /* 0/1; default 0 */
} rk_homotopy_cfg_t;
```

`rk_params_t` gains: `rk_homotopy_cfg_t homotopy; rk_scheduler_t scheduler; rk_eta_mode_t eta_mode;
double eta_start; double eta_end; double eta_schedule_power;`

`residual_recheck_fraction` is INTERNAL to `CalibState` only — never on the ABI, never in R wrapper.

### R API (`R/harvest.R`)

New named arguments (all default off / identity):

| Arg | Default | Maps to |
|-----|---------|---------|
| `homotopy_levels` | `1` | `homotopy.n_levels` |
| `homotopy_start_factor` | `1.0` | `homotopy.start_factor` |
| `homotopy_end_factor` | `1.0` | `homotopy.end_factor` |
| `homotopy_budget_p` | `0.5` | `homotopy.budget_split_p` |
| `scheduler` | `"round_robin"` | `scheduler.mode` |
| `eta_schedule` | `"fixed"` | `eta_schedule.mode` |
| `eta_start` | `1.0` | `eta_schedule.eta_start` |
| `eta_end` | `1.0` | `eta_schedule.eta_end` |
| `eta_schedule_power` | `0.5` | `eta_schedule.schedule_power` |

`scheduler` and `eta_schedule` use `match.arg()`.

---

## §3 Homotopy Semantics

1. **Level sequence:** generate `n_levels` values of `max_weight_l` by geometric interpolation:
   `max_weight_l = start_factor * (end_factor / start_factor)^(l / (n_levels - 1))` for
   `l = 0 .. n_levels - 1`. At `l = n_levels - 1`: `max_weight_l = end_factor * original_max_weight`.
2. **Budget split:** allocate total iteration budget `T` across levels using a Chizat-√t-inspired
   split: level `l` gets `floor(budget_split_p * T * (l+1) / sum(1..n_levels))` iterations (exact
   formula TBD in WU-3). The final level gets all remaining budget.
3. **Warm start:** after each level, pass the converged weight vector as the starting point for the
   next level without re-normalization (preserves sum = n invariant).
4. **Fallback:** if `n_levels == 1` (default), skip the outer homotopy loop entirely and call
   `inner_solve_one_level` once with the user-supplied `max_weight`. Net effect: zero change.

---

## §4 Priority Scheduler (Greenkhorn)

At each outer pass, instead of iterating margins `k = 0, 1, ..., K-1` in fixed order:

1. Compute `errRp[k]` (max absolute calibration error for margin k) for all K margins.
2. Select the argmax: `k* = argmax_k errRp[k]`.
3. Apply the BCD update for margin `k*` only.
4. Refresh the priority queue every `ceil(residual_recheck_fraction * K)` updates (internal
   constant = 0.1; gives full refresh at K ≥ 10). For K < 10, recompute every step.
5. Cap at K BCD steps per outer iteration to preserve the per-iteration cost bound.

When `scheduler = "round_robin"` (default), behaviour is identical to current iEPPA: margins
visited in order 0..K-1. Net effect at default: zero change.

---

## §5 Tang-η Dynamic Schedule

Between homotopy levels `l-1` and `l`, update `alm_mu`:

```
alm_mu_l = alm_mu_0 * eta_start * ((l) / (n_levels - 1))^schedule_power
```

When `eta_schedule = "fixed"` (default) and `eta_start = eta_end = 1.0`, the multiplier is always
1.0, so `alm_mu` is unchanged from its initial value. Net effect at default: zero change.

When `n_levels = 1` (no homotopy), the schedule produces a single multiplier applied at level 0
(no-op at defaults).

---

## §6 KKT Invariants at Each Level

After each homotopy level `l` converges, the following invariants must hold:

1. **Sum invariant:** `|sum(weights) - n| < 1e-10`.
2. **Calibration residual:** `max_k max_j |S_kj / W - tau_kj| <= tol_abs * sqrt(n_levels)` (the
   tolerance is relaxed at intermediate levels; the final level uses the user-supplied `tol_abs`).
3. **Bounds:** each calibrated weight satisfies `min_weight <= w_i <= max_weight_l` at level `l`.
4. **Non-negativity:** `min(weights) >= 0`.

---

## §7 Acceptance Criteria (A1–A6)

- **A1 (net-zero):** `harvest(data, target, ..., homotopy_levels=1, scheduler="round_robin",
  eta_schedule="fixed")` produces weights identical to the no-overlay call within 1e-12.
- **A2 (pairwise toggleability):** each overlay can be enabled independently without error and
  produces weights within `max_weight` bounds.
- **A3 (regression):** all 163 prior tests remain PASS; no new FAIL.
- **A4 (ABI size):** `sizeof(rk_params_t)` matches the updated `EXPECTED_RK_PARAMS_BYTES` constant.
- **A5 (struct defaults):** `rk_params_init()` initializes all new fields to their identity
  defaults (verifiable by memset-then-init pattern).
- **A6 (no C++ ABI break):** zero changes to existing struct member offsets; new fields appended
  only at the end of `rk_params_t`.

---

## §8 Per-Overlay Toggle Semantics

Each overlay is independently controlled by its own field(s):

| Overlay | Off (default) | On |
|---------|---------------|----|
| Homotopy (P-A) | `homotopy_levels = 1` (implies `enabled = false`) | `homotopy_levels > 1` |
| Scheduler (P-B) | `scheduler = "round_robin"` | `scheduler = "greedy"` |
| Tang-η (P-C) | `eta_schedule = "fixed"` | `eta_schedule = "tang_dynamic"` |

When all three are at defaults, behaviour is strictly identical to the current iEPPA implementation.
The overlays may be combined freely (each is a separate code path; no interaction assumptions are
encoded at the config layer — interaction correctness is validated in WU-6).

---

## §9 Rollback

WU-1 is pure scaffolding (no behavioural change). Rollback = `git revert <WU-1 commit SHA>`. No
data migration or ABI consumer update required because the new fields are appended at the end of
`rk_params_t` and `rk_params_init()` sets them to identity defaults — zero-initialised existing
callers continue to work.
