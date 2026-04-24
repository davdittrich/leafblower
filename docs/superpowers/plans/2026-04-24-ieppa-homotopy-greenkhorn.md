# iEPPA Homotopy + Greenkhorn + Tang Dynamic-η Implementation Plan (rev 2)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the stepstone-fulldata convergence gap vs `autumn::harvest(accelerate=TRUE)` — at least match its `max_err = 1.60e-3` at `max_iterations = 3000`, stretch 10× to `1.60e-4` — without replacing the iEPPA multiplicative-scaling solver core.

**Architecture:** Compose three iEPPA-preserving overlays backed by Schmitzer 2019 (arxiv:1610.06519), Chizat 2024 (arxiv:2408.11620), Altschuler–Weed–Rigollet 2017 (arxiv:1705.09634), and Tang et al. 2024 (arxiv:2403.05054):

1. **P-A (outer homotopy):** progressive `max_weight` tightening across `N_levels` outer stages, warm-started across levels, attacking the piecewise-linear capacity-clamp non-smoothness. Budget split per level follows a Chizat-√t-inspired schedule.
2. **P-B (Greenkhorn priority scheduler):** replace the round-robin margin loop inside one inner iEPPA pass with a max-residual priority-ordered single-margin update. Extends binary-Sinkhorn Greenkhorn to the K-way box-constrained case.
3. **Tang dynamic-η:** borrow only the dynamic regularization / proximal-weight schedule component (NOT Newton, NOT primal-dual replacement). Annealed proximal-weight schedule across outer homotopy levels.

The three overlays are orthogonal — homotopy is outer, priority is scheduling inside one inner pass, dynamic-η is the per-level proximal-weight schedule. All three are **default-off** so the existing 163-test baseline is preserved.

**Tech Stack:** C++17 (`src/ieppa.cpp`, `src/ieppa.hpp`, `src/types.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`, `src/leafblower.h`), R wrapper (`R/harvest.R`), Python wrapper (`python/leafblower/_harvest.py`, `python/leafblower/_bindings.cpp`), testthat 3, pytest, bench::mark, stepstone-fulldata fixtures at `benchmarks/stepstone_fulldata_bench_*`.

**Research basis:** `docs/investigations/2026-04-24-ieppa-accel-research.md`.

**Ground-truth R API (from `R/harvest.R` lines 37–62):** argument names are `target` (singular), `method`, `max_iterations`, `convergence = list(absolute = ..., pct = ...)`, `bounds_mode`, `attach_weights`, `weight_column`, `verbose`. All new overlay knobs added by this plan use these exact names and route through the `...` pass-through plus explicit new parameters where needed.

**Merge gate (must pass before the bundled PR merges):**

1. `devtools::test()` green.
2. `pytest` green.
3. `R CMD check --as-cran`: 0 ERROR, 0 WARNING, NOTEs ≤ baseline.
4. **Benchmark floor:** stepstone-fulldata at `max_iterations=3000, max_weight=5, convergence=list(absolute=1e-10)`: leafblower AB config (P-A + P-B) produces `max_err ≤ 1.60e-3`.
5. **kk1204 non-regression:** `max_err ≤ 1.322e-3 at max_iterations=500` (current baseline from leafblower-ylsy memory).
6. **Reference-commit agreement:** Pearson r ≥ 0.99 between new AB weights and `tests/testthat/fixtures/stepstone_reference.rds` (commit 8146894).

**Reported-but-not-gating metrics (for regression awareness, not merge-blocking):**
- Wall time of AB vs autumn (information only; the user said wall-time is non-comparable while iEPPA is NOCONV, so the ceiling is a monitoring metric, not a gate).
- Stretch: AB+η `max_err ≤ 1.60e-4`.
- Rate slope on stepstone-small trajectory (linear fit of `log(errRp)` vs `log(iter)`).
- Degenerate `M_cell = n` micro-bench wall-time ratio vs raking.

**Rule reminders (global CLAUDE.md + memory):**

- One beads ticket per task (plan task = one bead ticket; no bundling).
- Atomic per-WU commits. Subagent commit-merging defeats selective revert.
- `R CMD INSTALL --preclean .` gates every C++ edit before advancing.
- Full-budget trajectory probe required before asserting any rate claim.
- `harvest()` returns a data.frame (`attach_weights=TRUE`) by default — tests use `res$weights` or pass `attach_weights=FALSE` (memory).
- `R capture.output(type='output')` for `Rprintf` verbose logs (memory).
- Plan-review-gate approval already obtained before implementation begins.

---

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `src/types.hpp` | modify | Add `HomotopyConfigLbw`, `SchedulerConfigLbw`, `EtaScheduleConfigLbw` inside `CalibState`. |
| `src/ieppa.hpp` | modify | Extend `IEPPAResult` with `homotopy_levels_used`, `homotopy_final_factor`, `greedy_sweeps_taken`, `eta_final`. Trajectory logging stays internal (no new public field). |
| `src/ieppa.cpp` | modify | Refactor `ieppa_solve` into (a) `inner_solve_one_level` helper + (b) outer homotopy loop; add priority-scheduler code path; add dynamic-η schedule between levels; internal env-var-driven trajectory probe. |
| `src/c_api.cpp` | modify | Extend C ABI to carry the three config structs + new diagnostics fields. |
| `src/leafblower.h` | modify | Add `rk_scheduler_t`, `rk_eta_mode_t` enums and the new config-struct fields on `rk_calib_config_t` (or whatever the existing ABI struct is named). |
| `src/r_bridge.cpp` | modify | Marshal R-side arguments into the extended `CalibState`. |
| `R/harvest.R` | modify | Accept overlay knobs ONLY where they belong in user API: `homotopy_levels`, `homotopy_start_factor`, `homotopy_end_factor`, `homotopy_budget_p`, `scheduler`, `eta_schedule`, `eta_start`, `eta_end`, `eta_schedule_power`. Internal knobs (`residual_recheck_frac`, trajectory probe iterations) stay out of the signature. |
| `python/leafblower/_bindings.cpp` | modify | Mirror the new C ABI fields. |
| `python/leafblower/_harvest.py` | modify | Python wrapper accepts same overlay knobs (snake_case). |
| `tests/testthat/test-homotopy.R` | **create** | WU-3 behavioural tests. |
| `tests/testthat/test-priority-sweep.R` | **create** | WU-4 behavioural tests. |
| `tests/testthat/test-eta-schedule.R` | **create** | WU-5 behavioural tests. |
| `tests/testthat/test-convergence-trajectory.R` | **create** | WU-2 env-var-driven trajectory probe test. |
| `tests/testthat/test-config-defaults.R` | **create** | WU-1 net-zero + pairwise toggleability tests. |
| `tests/testthat/test-bench-gate.R` | **create** | WU-6 merge-gate assertions. |
| `python/leafblower/test_python.py` | modify | Parity tests mirroring the R ones. |
| `benchmarks/make_stepstone_small_fixture.R` | **create** | One-time helper to produce 10k-row stepstone-small fixture + targets. |
| `benchmarks/probe_baseline.R` | **create** | Trajectory probe runner. |
| `benchmarks/baseline_tuning_sweep.R` | **create** | WU-2.5 parameter-only-baseline falsification on stepstone-fulldata. |
| `benchmarks/stepstone_fulldata_homotopy.R` | **create** | WU-6 benchmark script (AB + autumn + baseline). |
| `docs/superpowers/specs/2026-04-24-ieppa-homotopy-greenkhorn-design.md` | **create** | Design spec, written in WU-1. |

---

## Work Units

Seven WUs. One beads ticket per WU. Atomic commit per WU. Do **not** merge WUs.

---

### WU-1: Design spec + CalibState config scaffolding (no behavioural change)

**Beads:** `bd create --title "WU-1: iEPPA homotopy/priority/eta config scaffolding" --type task --priority 2`

**Files:** create design spec, modify `src/types.hpp`, `src/leafblower.h`, `src/ieppa.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`, `R/harvest.R`. Create `tests/testthat/test-config-defaults.R`.

- [ ] **Step 1: Write the design spec.**

Create `docs/superpowers/specs/2026-04-24-ieppa-homotopy-greenkhorn-design.md` with these required sections (`§1` motivation + ylsy-memory citation, `§2` config contract with exact struct fields, `§3` homotopy semantics incl. geometric `max_weight` interpolation + Chizat-√t-inspired budget split + warm start, `§4` priority scheduler with argmax/K-step cap/`residual_recheck_frac` internal default 0.1, `§5` Tang-η schedule as `alm_mu` multiplier, `§6` KKT invariants at each level, `§7` acceptance criteria A1–A6, `§8` per-overlay toggle semantics). Reference `docs/investigations/2026-04-24-ieppa-accel-research.md` for the literature backing.

- [ ] **Step 2: Write failing net-zero + pairwise toggleability test.**

`tests/testthat/test-config-defaults.R`:

```r
test_that("all overlays default-off: baseline identical to defaulted", {
  set.seed(1)
  n <- 2000
  data <- data.frame(
    a = sample(1:3, n, replace = TRUE),
    b = sample(1:2, n, replace = TRUE)
  )
  target <- list(a = c(0.3, 0.5, 0.2), b = c(0.6, 0.4))
  baseline <- leafblower::harvest(
    data, target, max_weight = 3, method = "ieppa",
    attach_weights = FALSE
  )
  defaulted <- leafblower::harvest(
    data, target, max_weight = 3, method = "ieppa",
    homotopy_levels = 1,
    scheduler      = "round_robin",
    eta_schedule   = "fixed",
    attach_weights = FALSE
  )
  expect_equal(as.numeric(baseline), as.numeric(defaulted), tolerance = 1e-12)
})

test_that("pairwise toggleability: each overlay can be on/off independently", {
  set.seed(2)
  n <- 2000
  data <- data.frame(
    a = sample(1:3, n, replace = TRUE),
    b = sample(1:2, n, replace = TRUE)
  )
  target <- list(a = c(0.3, 0.5, 0.2), b = c(0.6, 0.4))
  common <- list(data = data, target = target, max_weight = 3,
                 method = "ieppa", attach_weights = FALSE)
  # Only A on; only B on; only eta on — all must run without error and
  # stay within user-supplied max_weight.
  a_only <- do.call(leafblower::harvest,
    c(common, list(homotopy_levels = 3, homotopy_start_factor = 5,
                   homotopy_end_factor = 1)))
  b_only <- do.call(leafblower::harvest,
    c(common, list(scheduler = "greedy")))
  e_only <- do.call(leafblower::harvest,
    c(common, list(homotopy_levels = 3, homotopy_start_factor = 5,
                   homotopy_end_factor = 1,
                   eta_schedule = "tang_dynamic",
                   eta_start = 5, eta_end = 1)))
  expect_true(max(a_only) <= 3 + 1e-10)
  expect_true(max(b_only) <= 3 + 1e-10)
  expect_true(max(e_only) <= 3 + 1e-10)
})
```

- [ ] **Step 3: Run tests to confirm failure (new args unused).**

Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-config-defaults.R")'`
Expected: FAIL with `unused argument (homotopy_levels = 1)` or similar; tests cannot exercise new args yet.

- [ ] **Step 4: Extend `src/types.hpp`.**

Append inside `namespace lbw`:

```cpp
struct HomotopyConfigLbw {
    int    n_levels        = 1;
    double start_factor    = 1.0;
    double end_factor      = 1.0;
    double budget_split_p  = 0.5;
    bool   enabled         = false;
};
enum class SchedulerMode   : int { ROUND_ROBIN = 0, GREEDY = 1 };
enum class EtaScheduleMode : int { FIXED = 0, TANG_DYNAMIC = 1 };
struct SchedulerConfigLbw {
    SchedulerMode mode = SchedulerMode::ROUND_ROBIN;
    double        residual_recheck_fraction = 0.1;  // internal
};
struct EtaScheduleConfigLbw {
    EtaScheduleMode mode = EtaScheduleMode::FIXED;
    double eta_start      = 1.0;
    double eta_end        = 1.0;
    double schedule_power = 0.5;
};
```

Attach to `CalibState` as members:

```cpp
HomotopyConfigLbw    homotopy;
SchedulerConfigLbw   scheduler;
EtaScheduleConfigLbw eta_schedule;
```

- [ ] **Step 5: Extend `src/leafblower.h`.**

Add enums + fields to the existing C-ABI calibration-config struct (the one consumed by `c_api.cpp`):

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

Extend the existing config struct with: `rk_homotopy_cfg_t homotopy; rk_scheduler_t scheduler; rk_eta_mode_t eta_mode; double eta_start; double eta_end; double eta_schedule_power;`. Internal `residual_recheck_fraction` is NOT on the ABI — it is a CalibState default hidden from callers.

- [ ] **Step 6: Extend `IEPPAResult` in `src/ieppa.hpp`.**

```cpp
struct IEPPAResult {
    // ... existing fields unchanged ...
    int    homotopy_levels_used  = 0;   // 0 iff homotopy disabled
    double homotopy_final_factor = 1.0;
    int    greedy_sweeps_taken   = 0;   // per last inner pass
    double eta_final             = 0.0;
};
```

Trajectory logging is internal (WU-2); it does not extend `IEPPAResult`.

- [ ] **Step 7: Extend `src/c_api.cpp` + `src/r_bridge.cpp`.**

Unpack the new fields from the R-side list / argument vector and fill `CalibState`. Keep defaults so existing callers see no change.

- [ ] **Step 8: Extend `R/harvest.R` signature.**

Insert new named arguments before the existing `...`:

```r
harvest <- function(
  data, target,
  min_weight       = 0,
  max_weight       = 5,
  method           = "ieppa",
  verbose          = 0,
  max_iterations   = 500,
  start_weights    = NULL,
  attach_weights   = TRUE,
  weight_column    = "weights",
  convergence      = list(),
  bounds_mode      = "cell",
  # --- new overlay knobs (all default off / identity) ---
  homotopy_levels       = 1,
  homotopy_start_factor = 1.0,
  homotopy_end_factor   = 1.0,
  homotopy_budget_p     = 0.5,
  scheduler             = c("round_robin", "greedy"),
  eta_schedule          = c("fixed", "tang_dynamic"),
  eta_start             = 1.0,
  eta_end               = 1.0,
  eta_schedule_power    = 0.5,
  # --- end new ---
  select_params = NULL, select_function = NULL, error_function = NULL,
  adaptive_order = NULL, enforce_mean = TRUE, accelerate = FALSE,
  add_na_proportion = FALSE, auto_collapse = FALSE, collapse_vars = NULL,
  target_map = NULL, design_weights = NULL,
  ...
) {
  scheduler    <- match.arg(scheduler)
  eta_schedule <- match.arg(eta_schedule)
  # ... forward through .Call's argument list ...
}
```

Update the `.Call("C_rk_calibrate", ...)` argument list in `R/harvest.R` to include the new fields; match positional ordering on the `r_bridge.cpp` side.

- [ ] **Step 9: Build + run the config-default tests.**

Run: `R CMD INSTALL --preclean .`
Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-config-defaults.R")'`
Expected: both tests PASS.

- [ ] **Step 10: Full regression suite — confirm baseline preserved.**

Run: `Rscript -e 'devtools::test()'`
Expected: no new FAILs relative to `FAIL 0 | PASS 163` baseline; PASS count ≥ 165 (prior 163 + two new).

- [ ] **Step 11: Commit WU-1 atomically.**

```bash
git add src/types.hpp src/leafblower.h src/ieppa.hpp src/c_api.cpp \
        src/r_bridge.cpp R/harvest.R \
        tests/testthat/test-config-defaults.R \
        docs/superpowers/specs/2026-04-24-ieppa-homotopy-greenkhorn-design.md
git commit -m "$(cat <<'EOF'
feat(ieppa): config scaffolding for homotopy/priority/eta overlays

Adds HomotopyConfigLbw, SchedulerConfigLbw, EtaScheduleConfigLbw structs to
CalibState and the corresponding C ABI / R wrapper. All overlays default off;
existing 163 tests remain green plus two new net-zero / pairwise-toggle tests
(165). Design spec:
docs/superpowers/specs/2026-04-24-ieppa-homotopy-greenkhorn-design.md.
EOF
)"
```

- [ ] **Step 12: Close beads ticket.**

```bash
bd close <WU-1 ticket-id>
```

---

### WU-2: Internal `errRp` trajectory probe (env-var driven)

**Beads:** `bd create --title "WU-2: iEPPA trajectory probe via env var" --type task --priority 2`

**Rationale:** memory mandates full-budget trajectory probe before any rate claim. The probe is an internal validation hook — NOT a user-facing argument. Activated by environment variable `LBW_TRAJECTORY_AT` (comma-separated iteration numbers) and written to `LBW_TRAJECTORY_OUT` (file path). Default off.

**Files:** modify `src/ieppa.cpp`. Create `tests/testthat/test-convergence-trajectory.R`, `benchmarks/make_stepstone_small_fixture.R`, `benchmarks/probe_baseline.R`.

- [ ] **Step 1: Write failing trajectory-probe test.**

`tests/testthat/test-convergence-trajectory.R`:

```r
test_that("LBW_TRAJECTORY_AT env var produces probe file", {
  skip_if_not_installed("arrow")
  fx <- test_path("fixtures/stepstone_small.parquet")
  tg <- test_path("fixtures/stepstone_small_targets.rds")
  skip_if(!file.exists(fx) || !file.exists(tg))
  data   <- arrow::read_parquet(fx)
  target <- readRDS(tg)
  out <- tempfile(fileext = ".csv")
  withr::with_envvar(
    c(LBW_TRAJECTORY_AT  = "1,10,50,100,200,500",
      LBW_TRAJECTORY_OUT = out),
    {
      leafblower::harvest(
        data, target, max_weight = 5,
        method = "ieppa",
        max_iterations = 500,
        convergence = list(absolute = 1e-12),   # unreachable → full budget
        attach_weights = FALSE
      )
    }
  )
  expect_true(file.exists(out))
  probes <- utils::read.csv(out)
  expect_named(probes, c("iter", "errRp"))
  expect_equal(probes$iter, c(1, 10, 50, 100, 200, 500))
  expect_length(probes$errRp, 6)
  expect_true(all(diff(probes$errRp) <= 1e-12))   # monotone non-increasing
})
```

- [ ] **Step 2: Run test to confirm failure (env var is a no-op currently).**

Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-trajectory.R")'`
Expected: FAIL at `file.exists(out)`.

- [ ] **Step 3: Create small fixture.**

`benchmarks/make_stepstone_small_fixture.R`:

```r
src <- "benchmarks/stepstone_fulldata_bench_data.parquet"
stopifnot(file.exists(src))
data <- arrow::read_parquet(src)
set.seed(42)
small <- data[sample.int(nrow(data), 10000), , drop = FALSE]
arrow::write_parquet(small, "tests/testthat/fixtures/stepstone_small.parquet")
target <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
saveRDS(target, "tests/testthat/fixtures/stepstone_small_targets.rds")
```

Run once: `Rscript benchmarks/make_stepstone_small_fixture.R`.

- [ ] **Step 4: Implement probe in `src/ieppa.cpp`.**

Near `ieppa_solve` entry, parse env vars once:

```cpp
static std::vector<int> parse_trajectory_iters() {
    const char* s = std::getenv("LBW_TRAJECTORY_AT");
    if (!s || !*s) return {};
    std::vector<int> out;
    std::string buf;
    for (const char* p = s;; ++p) {
        if (*p == ',' || *p == '\0') {
            if (!buf.empty()) out.push_back(std::stoi(buf));
            buf.clear();
            if (*p == '\0') break;
        } else buf.push_back(*p);
    }
    std::sort(out.begin(), out.end());
    return out;
}
```

At the end of `ieppa_solve`, if probes were populated, write `LBW_TRAJECTORY_OUT` as a CSV with `iter,errRp` header. Inside the main iteration loop, after the `errRp` computation at each check interval, if the current iter matches the next pending probe iter, append `{iter, errRp}` to an in-memory `std::vector<std::pair<int,double>>`.

Critical: the probe must record `errRp` *as computed at that iter* — must respect `kErrCheckInterval`. If a requested probe iter is not on a check interval, compute `errRp` inline at that iter (one extra pass).

- [ ] **Step 5: Build + run test.**

Run: `R CMD INSTALL --preclean .`
Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-convergence-trajectory.R")'`
Expected: PASS.

- [ ] **Step 6: Full regression.**

Run: `Rscript -e 'devtools::test()'`
Expected: PASS count ≥ 166 (prior 165 + one new).

- [ ] **Step 7: Record stepstone-small baseline trajectory.**

`benchmarks/probe_baseline.R`:

```r
data   <- arrow::read_parquet("tests/testthat/fixtures/stepstone_small.parquet")
target <- readRDS("tests/testthat/fixtures/stepstone_small_targets.rds")
out <- "benchmarks/baseline_trajectory_stepstone_small.csv"
Sys.setenv(LBW_TRAJECTORY_AT  = "1,10,50,100,200,500,1000,2000,3000",
           LBW_TRAJECTORY_OUT = out)
leafblower::harvest(
  data, target, max_weight = 5,
  method = "ieppa",
  max_iterations = 3000,
  convergence = list(absolute = 1e-12),
  attach_weights = FALSE
)
cat("probe written to ", out, "\n")
```

Run once: `Rscript benchmarks/probe_baseline.R`.

- [ ] **Step 8: Commit WU-2.**

```bash
git add src/ieppa.cpp tests/testthat/test-convergence-trajectory.R \
        tests/testthat/fixtures/stepstone_small.parquet \
        tests/testthat/fixtures/stepstone_small_targets.rds \
        benchmarks/make_stepstone_small_fixture.R \
        benchmarks/probe_baseline.R \
        benchmarks/baseline_trajectory_stepstone_small.csv
git commit -m "$(cat <<'EOF'
feat(ieppa): internal errRp trajectory probe via env vars

Activated by LBW_TRAJECTORY_AT (csv iters) + LBW_TRAJECTORY_OUT (file path).
Writes iter,errRp CSV. Internal-only — no new public API. Baseline probe on
stepstone-small captured for rate comparison.
EOF
)"
```

- [ ] **Step 9: Close beads ticket.**

---

### WU-2.5: Parameter-only-baseline falsification on stepstone-fulldata

**Beads:** `bd create --title "WU-2.5: falsify simpler alternative — existing-knob tuning on stepstone" --type task --priority 2`

**Rationale:** before investing in three overlays, empirically confirm that tuning existing knobs (`bounds_mode`, damping, `max_iterations` budget) cannot reach the 1.60e-3 floor on stepstone-fulldata. If they can, the overlays are unjustified; halt.

**Files:** create `benchmarks/baseline_tuning_sweep.R`. No source edits.

- [ ] **Step 1: Write tuning sweep.**

`benchmarks/baseline_tuning_sweep.R`:

```r
library(leafblower)
data   <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
target <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")

max_err_of <- function(w, data, target) {
  errs <- numeric(length(target))
  for (k in seq_along(target)) {
    tab <- tapply(w, data[[names(target)[k]]], sum) / sum(w)
    errs[k] <- max(abs(tab - target[[k]]))
  }
  max(errs)
}

grid <- expand.grid(
  bounds_mode    = c("cell", "unit"),
  max_iterations = c(3000, 6000, 10000),
  stringsAsFactors = FALSE
)

results <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, , drop = FALSE]
  t0 <- Sys.time()
  w <- leafblower::harvest(
    data, target, max_weight = 5,
    method         = "ieppa",
    bounds_mode    = g$bounds_mode,
    max_iterations = g$max_iterations,
    convergence    = list(absolute = 1e-10),
    attach_weights = FALSE
  )
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  data.frame(g, wall_s = wall, max_err = max_err_of(w, data, target))
}))

print(results, digits = 4)
saveRDS(results, "benchmarks/baseline_tuning_sweep.rds")

hit <- any(results$max_err <= 1.60e-3)
if (hit) {
  stop("FALSIFIED: existing knobs reach 1.60e-3 on stepstone-fulldata. ",
       "Overlays may be unnecessary. Halt and revisit plan.")
} else {
  cat("Confirmed: no existing-knob configuration reaches 1.60e-3. ",
      "Overlays justified.\n", sep = "")
}
```

- [ ] **Step 2: Run sweep.**

Run: `Rscript benchmarks/baseline_tuning_sweep.R`
Expected: all configurations miss 1.60e-3 (i.e., `max_err > 1.60e-3`); script prints the confirmation message. If the sweep falsifies the hypothesis (any config hits the floor), STOP, file a beads HUMAN ticket, and revisit the plan — the overlays may be unjustified.

- [ ] **Step 3: Commit WU-2.5.**

```bash
git add benchmarks/baseline_tuning_sweep.R benchmarks/baseline_tuning_sweep.rds
git commit -m "$(cat <<'EOF'
bench(ieppa): falsify simpler-alternative hypothesis on stepstone-fulldata

Sweep of bounds_mode x max_iterations on stepstone-fulldata confirms no
existing-knob configuration reaches 1.60e-3 max_err. Justifies proceeding
with P-A / P-B / Tang-eta overlays.
EOF
)"
```

- [ ] **Step 4: Close beads ticket.**

---

### WU-3: P-A homotopy outer loop

**Beads:** `bd create --title "WU-3: iEPPA P-A progressive max_weight tightening" --type task --priority 2`

**Files:** modify `src/ieppa.cpp`. Create `tests/testthat/test-homotopy.R`.

- [ ] **Step 1: Refactor inner iteration into a helper.**

Hoist the body of `for (int iter = 1; iter <= st.inner_max_iter; iter++)` (inside `ieppa_solve`) into a private static function:

```cpp
namespace {
struct InnerResult { int iters_used; double errRp; int status; bool early_converged; };
static InnerResult inner_solve_one_level(
    CalibState& st,
    CellTable& ct,
    std::vector<double>& W,
    std::vector<double>& X_cur,
    std::vector<double>& lf,
    int budget,
    double current_max_weight,
    double tol_for_this_level,
    IEPPAResult& res,
    std::vector<std::pair<int,double>>& probe_samples,
    const std::vector<int>& probe_iters);
}
```

Caller-owned `W`, `X_cur`, `lf` enable warm-start across levels. `current_max_weight` replaces the literal `st.max_weight` in the capacity-clamp block. `tol_for_this_level` enables early exit. `early_converged=true` if inner reaches `tol_for_this_level` before exhausting `budget`. Probe-writing hook passes through so the trajectory CSV still receives writes across all levels.

This refactor is net-zero behaviourally when called with `budget = st.inner_max_iter`, `current_max_weight = st.max_weight`, `tol_for_this_level = st.tol_abs`.

- [ ] **Step 2: Run existing suite to confirm refactor is inert.**

Run: `R CMD INSTALL --preclean .`
Run: `Rscript -e 'devtools::test()'`
Expected: PASS count unchanged (≥ 166).

- [ ] **Step 3: Add outer homotopy driver inside `ieppa_solve`.**

```cpp
IEPPAResult ieppa_solve(CalibState& st) {
    // ... existing preamble: CellTable, W, lf initialization ...

    IEPPAResult res{};
    std::vector<std::pair<int,double>> probe_samples;
    std::vector<int> probe_iters = parse_trajectory_iters();

    const int N = (st.homotopy.enabled && st.homotopy.n_levels > 1)
                  ? st.homotopy.n_levels : 1;
    const double k_start = st.homotopy.start_factor;
    const double k_end   = st.homotopy.end_factor;
    const double p       = st.homotopy.budget_split_p;

    double weight_sum = 0.0;
    for (int i = 0; i < N; i++) weight_sum += std::pow(i + 1, p);

    double alm_mu_base = st.alm_mu;

    int iters_absorbed = 0;
    for (int lvl = 0; lvl < N; lvl++) {
        const double frac = (N == 1) ? 0.0 :
            static_cast<double>(lvl) / (N - 1);
        const double factor = (N == 1) ? 1.0 :
            k_start * std::pow(k_end / k_start, frac);
        const double cur_maxw = st.max_weight * factor;
        const int budget_i = std::max(1,
            static_cast<int>(std::round(
                st.inner_max_iter * std::pow(lvl + 1, p) / weight_sum)));
        const double tol_i = (lvl == N - 1) ? st.tol_abs
                                            : std::max(st.tol_abs, 1e-5);

        // Tang dynamic-eta (WU-5 writes this branch; safe no-op until then).
        if (st.eta_schedule.mode == EtaScheduleMode::TANG_DYNAMIC && N > 1) {
            const double eta_i = st.eta_schedule.eta_start * std::pow(
                st.eta_schedule.eta_end / st.eta_schedule.eta_start, frac);
            st.alm_mu = eta_i * alm_mu_base;
            res.eta_final = eta_i;
        }

        InnerResult ir = inner_solve_one_level(
            st, ct, W, X_cur, lf, budget_i, cur_maxw, tol_i, res,
            probe_samples, probe_iters);

        iters_absorbed += ir.iters_used;
        res.iterations = iters_absorbed;
        res.homotopy_levels_used = lvl + 1;
        res.homotopy_final_factor = factor;
        res.max_error = ir.errRp;

        if (ir.status == RK_ERR_INFEAS || ir.status == RK_ERR_BADARG) {
            res.status = ir.status;
            break;
        }
        if (ir.early_converged && lvl < N - 1) {
            // Level converged to intermediate tol_i early; proceed to next level
            // warm-started. Skip the rest of budget_i.
            continue;
        }
        if (lvl == N - 1 && ir.errRp < st.tol_abs) {
            res.status = RK_OK;
        }
    }

    write_trajectory_csv_if_enabled(probe_samples);
    // ... existing post-inner expansion + normalization ...
    return res;
}
```

In `inner_solve_one_level`, replace the literal `st.max_weight` everywhere inside the capacity-clamp block with `current_max_weight`. Add the early-exit condition: after the `errRp < tol_for_this_level` branch, return `InnerResult{...,.early_converged=true}`.

- [ ] **Step 4: Write failing homotopy behaviour test.**

`tests/testthat/test-homotopy.R`:

```r
test_that("P-A homotopy reaches target max_weight and improves errRp on stepstone-small", {
  fx <- test_path("fixtures/stepstone_small.parquet")
  tg <- test_path("fixtures/stepstone_small_targets.rds")
  skip_if(!file.exists(fx) || !file.exists(tg))
  data   <- arrow::read_parquet(fx)
  target <- readRDS(tg)

  tmp_base <- tempfile(fileext = ".csv")
  withr::with_envvar(
    c(LBW_TRAJECTORY_AT  = "1,10,50,100,200,500",
      LBW_TRAJECTORY_OUT = tmp_base),
    leafblower::harvest(data, target, max_weight = 5,
                        method = "ieppa",
                        max_iterations = 500,
                        convergence = list(absolute = 1e-12),
                        attach_weights = FALSE))
  base <- utils::read.csv(tmp_base)

  tmp_homo <- tempfile(fileext = ".csv")
  withr::with_envvar(
    c(LBW_TRAJECTORY_AT  = "1,10,50,100,200,500",
      LBW_TRAJECTORY_OUT = tmp_homo),
    homo <- leafblower::harvest(data, target, max_weight = 5,
                                method = "ieppa",
                                max_iterations = 500,
                                convergence = list(absolute = 1e-12),
                                homotopy_levels = 5,
                                homotopy_start_factor = 10,
                                homotopy_end_factor = 1,
                                homotopy_budget_p = 0.5,
                                attach_weights = FALSE))
  expect_true(max(homo) <= 5 + 1e-10)
  probe_homo <- utils::read.csv(tmp_homo)
  base_err <- tail(base$errRp, 1)
  homo_err <- tail(probe_homo$errRp, 1)
  expect_lt(homo_err, 0.7 * base_err)   # ≥30% reduction (A2)
})
```

- [ ] **Step 5: Run failing test.**

Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-homotopy.R")'`
Expected: FAIL until Step 3 is wired end-to-end.

- [ ] **Step 6: Iterate implementation → PASS.**

Build + run. If the 30% gate is consistently missed, the hypothesis is falsified — halt and file a beads HUMAN ticket with the trajectory CSVs attached.

- [ ] **Step 7: Full regression.**

All tests green; PASS count ≥ 167.

- [ ] **Step 8: Commit WU-3.**

```bash
git add src/ieppa.cpp tests/testthat/test-homotopy.R
git commit -m "$(cat <<'EOF'
feat(ieppa): P-A progressive max_weight homotopy outer loop

Outer driver invokes inner_solve_one_level across N levels with warm-started
W/X_cur/lf. max_weight tightens geometrically from start_factor*target to
end_factor*target. Per-level budget follows Chizat-sqrt(t) inspired schedule
(budget_split_p). Intermediate-level tolerance permits early exit and
warm-jump to the next level. stepstone-small errRp reduced by >=30% vs
baseline within same iter budget.
EOF
)"
```

- [ ] **Step 9: Close beads ticket.**

---

### WU-4: P-B Greenkhorn priority scheduler

**Beads:** `bd create --title "WU-4: iEPPA P-B greedy residual-priority margin scheduler" --type task --priority 2`

**Files:** modify `src/ieppa.cpp`, `src/ieppa.hpp`, `R/harvest.R`, `src/r_bridge.cpp` (to surface `attr(result, "iterations")`). Create `tests/testthat/test-priority-sweep.R`.

- [ ] **Step 1: Hoist the single-margin sweep body into a helper.**

Extract the existing linear-path sweep inner body (src/ieppa.cpp approx lines 237–297) into:

```cpp
static bool apply_single_margin_linear(
    CalibState& st, CellTable& ct,
    int k,
    std::vector<double>& X_cur,
    std::vector<double>& W,
    std::vector<double>& lf,
    std::vector<double>& S_lin,
    std::vector<double>& inv_f_old_lin,
    const std::vector<int>& cat_offset,
    double current_max_weight
);
```

Returns `true` on overflow trip (existing semantics). Do the same for the log-space path. The round-robin sweep becomes `for k in 0..K-1 overflow_trip |= apply_single_margin_*(...)`.

- [ ] **Step 2: Run suite → refactor inert.**

Run: `R CMD INSTALL --preclean .`; `Rscript -e 'devtools::test()'`. All green.

- [ ] **Step 3: Add per-margin residual helper.**

```cpp
static double compute_margin_errRp(
    const CalibState& st, const CellTable& ct,
    int k, const std::vector<double>& X, double W_total,
    std::vector<double>& S_lin,
    const std::vector<int>& cat_offset);
```

Body: accumulate `S_lin[j]` via sequential `X` scan as in the existing errRp block, return `max_j |S_lin[j]/W_total - st.targets[k][j]|`.

- [ ] **Step 4: Implement greedy scheduler inside `inner_solve_one_level`.**

```cpp
if (st.scheduler.mode == SchedulerMode::GREEDY) {
    std::vector<double> per_margin_err(st.K);
    double W_total = 0.0;
    for (int c = 0; c < ct.M_cell; c++) W_total += X[c];
    for (int k = 0; k < st.K; k++) {
        per_margin_err[k] = compute_margin_errRp(
            st, ct, k, X, W_total, S_lin, cat_offset);
    }
    const double initial_sum = std::accumulate(
        per_margin_err.begin(), per_margin_err.end(), 0.0);
    const double stop_threshold =
        st.scheduler.residual_recheck_fraction * initial_sum;
    int greedy_steps = 0;
    while (greedy_steps < st.K) {
        int k_star = static_cast<int>(std::distance(
            per_margin_err.begin(),
            std::max_element(per_margin_err.begin(), per_margin_err.end())));
        bool trip = use_linear
            ? apply_single_margin_linear(
                st, ct, k_star, X_cur, W, lf, S_lin, inv_f_old_lin,
                cat_offset, current_max_weight)
            : apply_single_margin_log(/* analogous args */);
        res.greedy_sweeps_taken++;
        if (trip) { overflow_trip = true; break; }
        greedy_steps++;
        W_total = 0.0;
        for (int c = 0; c < ct.M_cell; c++) W_total += X[c];
        per_margin_err[k_star] = compute_margin_errRp(
            st, ct, k_star, X, W_total, S_lin, cat_offset);
        const double cur_sum = std::accumulate(
            per_margin_err.begin(), per_margin_err.end(), 0.0);
        if (cur_sum < stop_threshold) break;
    }
} else {
    for (int k = 0; k < st.K && !overflow_trip; k++) {
        overflow_trip = apply_single_margin_linear(/* ... */);
    }
}
```

- [ ] **Step 5: Surface `iterations` attribute.**

In `R/harvest.R`, after the `.Call`, attach `attr(data, "iterations") <- calib_result$iterations`. In the `attach_weights = FALSE` branch, attach to the weights vector.

- [ ] **Step 6: Write failing priority-scheduler test.**

`tests/testthat/test-priority-sweep.R`:

```r
test_that("greedy preserves feasibility + Pearson agreement vs round-robin", {
  set.seed(7)
  n <- 5000
  data <- data.frame(
    a = sample(1:2, n, replace = TRUE),
    b = sample(1:5, n, replace = TRUE),
    c = sample(1:4, n, replace = TRUE)
  )
  target <- list(a = c(0.6, 0.4), b = rep(0.2, 5),
                 c = c(0.1, 0.2, 0.3, 0.4))
  rr  <- leafblower::harvest(data, target, max_weight = 4,
                             method = "ieppa",
                             max_iterations = 200,
                             convergence = list(absolute = 1e-5),
                             scheduler = "round_robin",
                             attach_weights = FALSE)
  grd <- leafblower::harvest(data, target, max_weight = 4,
                             method = "ieppa",
                             max_iterations = 200,
                             convergence = list(absolute = 1e-5),
                             scheduler = "greedy",
                             attach_weights = FALSE)
  expect_true(max(grd) <= 4 + 1e-10)
  expect_gt(cor(as.numeric(rr), as.numeric(grd)), 0.995)
})

test_that("greedy reaches tolerance in fewer inner iterations on 2-cat-heavy problem", {
  set.seed(11)
  n <- 4000
  data <- data.frame(
    easy  = sample(1:2, n, replace = TRUE, prob = c(0.51, 0.49)),
    hardA = sample(1:10, n, replace = TRUE),
    hardB = sample(1:10, n, replace = TRUE)
  )
  target <- list(easy = c(0.5, 0.5),
                 hardA = rep(0.1, 10),
                 hardB = rep(0.1, 10))
  rr  <- leafblower::harvest(data, target, max_weight = 3,
                             method = "ieppa",
                             max_iterations = 300,
                             convergence = list(absolute = 1e-4),
                             scheduler = "round_robin",
                             attach_weights = FALSE)
  grd <- leafblower::harvest(data, target, max_weight = 3,
                             method = "ieppa",
                             max_iterations = 300,
                             convergence = list(absolute = 1e-4),
                             scheduler = "greedy",
                             attach_weights = FALSE)
  expect_lt(attr(grd, "iterations"), attr(rr, "iterations"))
})
```

- [ ] **Step 7: Run failing tests → iterate → PASS.**

- [ ] **Step 8: Full regression.**

- [ ] **Step 9: Trajectory probe with greedy alone on stepstone-small.**

```bash
LBW_TRAJECTORY_AT=1,10,50,100,200,500 \
LBW_TRAJECTORY_OUT=benchmarks/greedy_trajectory_stepstone_small.csv \
Rscript -e 'data <- arrow::read_parquet("tests/testthat/fixtures/stepstone_small.parquet"); \
target <- readRDS("tests/testthat/fixtures/stepstone_small_targets.rds"); \
leafblower::harvest(data, target, max_weight = 5, method = "ieppa", \
  max_iterations = 500, convergence = list(absolute = 1e-12), \
  scheduler = "greedy", attach_weights = FALSE)'
```

Confirm terminal errRp at iter=500 is lower than baseline.

- [ ] **Step 10: Commit WU-4.**

```bash
git add src/ieppa.cpp src/ieppa.hpp src/r_bridge.cpp R/harvest.R \
        tests/testthat/test-priority-sweep.R \
        benchmarks/greedy_trajectory_stepstone_small.csv
git commit -m "$(cat <<'EOF'
feat(ieppa): P-B Greenkhorn priority-ordered margin scheduler

Replaces round-robin with argmax-residual margin selection inside the inner
iEPPA pass. K-way extension of Altschuler-Weed-Rigollet Greenkhorn; stops at
residual_recheck_fraction of initial sum (internal 0.1 default) or after K
steps. Default scheduler remains round_robin — no behavioural change unless
opted in. Exposes attr(result, "iterations").
EOF
)"
```

- [ ] **Step 11: Close beads ticket.**

---

### WU-5: Tang dynamic-η inner schedule

**Beads:** `bd create --title "WU-5: iEPPA Tang 2024 dynamic-eta schedule (schedule only)" --type task --priority 2`

**Scope note:** borrow ONLY the scheduling of the proximal weight η from Tang 2024 §4. No Newton. No primal-dual. η multiplies existing `alm_mu`. Default `eta_schedule = "fixed"` is net-zero.

**Files:** modify `src/ieppa.cpp` (remove the WU-3 stub Tang-η block if any remains; wire end-to-end). Create `tests/testthat/test-eta-schedule.R`.

- [ ] **Step 1: Confirm the Tang-η branch inserted in WU-3 Step 3 works end-to-end.**

Already wired in WU-3 Step 3 as a gated branch on `st.eta_schedule.mode`. Confirm `alm_mu_base = st.alm_mu` is captured at function entry before the outer loop.

- [ ] **Step 2: Write failing eta-schedule test.**

`tests/testthat/test-eta-schedule.R`:

```r
test_that("Tang dynamic-eta reduces errRp at matched iter budget on stepstone-small", {
  fx <- test_path("fixtures/stepstone_small.parquet")
  tg <- test_path("fixtures/stepstone_small_targets.rds")
  skip_if(!file.exists(fx) || !file.exists(tg))
  data   <- arrow::read_parquet(fx)
  target <- readRDS(tg)
  common <- list(data = data, target = target,
                 max_weight = 5, method = "ieppa",
                 max_iterations = 500,
                 convergence = list(absolute = 1e-12),
                 homotopy_levels = 5,
                 homotopy_start_factor = 10,
                 homotopy_end_factor = 1,
                 attach_weights = FALSE)
  tmp_fixed <- tempfile(fileext = ".csv")
  tmp_dyn   <- tempfile(fileext = ".csv")
  withr::with_envvar(
    c(LBW_TRAJECTORY_AT = "100,200,500", LBW_TRAJECTORY_OUT = tmp_fixed),
    do.call(leafblower::harvest,
            c(common, list(eta_schedule = "fixed"))))
  withr::with_envvar(
    c(LBW_TRAJECTORY_AT = "100,200,500", LBW_TRAJECTORY_OUT = tmp_dyn),
    do.call(leafblower::harvest,
            c(common, list(eta_schedule = "tang_dynamic",
                           eta_start = 10, eta_end = 1,
                           eta_schedule_power = 0.5))))
  fixed_err <- tail(utils::read.csv(tmp_fixed)$errRp, 1)
  dyn_err   <- tail(utils::read.csv(tmp_dyn)$errRp, 1)
  expect_lt(dyn_err, fixed_err)
})
```

- [ ] **Step 3: Run failing → iterate → PASS.**

Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-eta-schedule.R")'`
Expected: PASS after wiring.

- [ ] **Step 4: Full regression + `R CMD check --as-cran` NOTE count unchanged.**

Run: `Rscript -e 'devtools::test()'`
Run: `R CMD build .`
Run: `R CMD check --as-cran leafblower_*.tar.gz`
Expected: 0 ERROR, 0 WARNING, NOTEs ≤ baseline.

- [ ] **Step 5: Commit WU-5.**

```bash
git add src/ieppa.cpp tests/testthat/test-eta-schedule.R
git commit -m "$(cat <<'EOF'
feat(ieppa): Tang 2024 dynamic-eta schedule (schedule-only borrow)

Adds eta_schedule="tang_dynamic" that anneals the proximal ALM weight across
homotopy levels per Tang 2024 sec 4 scheduling. Newton and primal-dual
components of Tang 2024 are intentionally NOT borrowed (excluded by user
directive, shelved per leafblower-ylsy). Default eta_schedule="fixed" keeps
prior behaviour.
EOF
)"
```

- [ ] **Step 6: Close beads ticket.**

---

### WU-6: Integrated stepstone-fulldata benchmark + merge gate

**Beads:** `bd create --title "WU-6: stepstone-fulldata merge gate + kk1204 non-regression" --type task --priority 1`

**Files:** create `benchmarks/stepstone_fulldata_homotopy.R`, `tests/testthat/test-bench-gate.R`.

- [ ] **Step 1: Write benchmark script.**

`benchmarks/stepstone_fulldata_homotopy.R`:

```r
library(leafblower)
stopifnot(requireNamespace("autumn", quietly = TRUE))
data   <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
target <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")

common_lb <- list(data = data, target = target, max_weight = 5,
                  method = "ieppa",
                  max_iterations = 3000,
                  convergence = list(absolute = 1e-10),
                  attach_weights = FALSE)

max_err_of <- function(w, data, target) {
  errs <- numeric(length(target))
  for (k in seq_along(target)) {
    tab <- tapply(w, data[[names(target)[k]]], sum) / sum(w)
    errs[k] <- max(abs(tab - target[[k]]))
  }
  max(errs)
}

time_it <- function(expr) {
  t0 <- Sys.time()
  val <- force(expr)
  list(value = val, wall_s = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# autumn reference
autumn_run <- time_it(
  autumn::harvest(data, target = target, max_weight = 5,
                  convergence = list(absolute = 1e-10),
                  max_iterations = 3000,
                  accelerate = TRUE,
                  attach_weights = FALSE))
# baseline leafblower — all overlays off
base_run <- time_it(do.call(leafblower::harvest, common_lb))
# P-A alone
A_run   <- time_it(do.call(leafblower::harvest,
  c(common_lb, list(homotopy_levels = 5, homotopy_start_factor = 10,
                    homotopy_end_factor = 1, homotopy_budget_p = 0.5))))
# P-A + P-B
AB_run  <- time_it(do.call(leafblower::harvest,
  c(common_lb, list(homotopy_levels = 5, homotopy_start_factor = 10,
                    homotopy_end_factor = 1, homotopy_budget_p = 0.5,
                    scheduler = "greedy"))))
# P-A + P-B + Tang-eta
ABE_run <- time_it(do.call(leafblower::harvest,
  c(common_lb, list(homotopy_levels = 5, homotopy_start_factor = 10,
                    homotopy_end_factor = 1, homotopy_budget_p = 0.5,
                    scheduler = "greedy",
                    eta_schedule = "tang_dynamic",
                    eta_start = 10, eta_end = 1, eta_schedule_power = 0.5))))

report <- data.frame(
  config  = c("autumn_accel","leafblower_base","leafblower_A",
              "leafblower_AB","leafblower_ABE"),
  wall_s  = c(autumn_run$wall_s, base_run$wall_s, A_run$wall_s,
              AB_run$wall_s, ABE_run$wall_s),
  max_err = c(max_err_of(autumn_run$value, data, target),
              max_err_of(base_run$value,   data, target),
              max_err_of(A_run$value,      data, target),
              max_err_of(AB_run$value,     data, target),
              max_err_of(ABE_run$value,    data, target)))
# Pearson r vs commit-8146894 reference
ref <- readRDS("tests/testthat/fixtures/stepstone_reference.rds")
pearson_r <- sapply(
  list(base = base_run$value, A = A_run$value,
       AB = AB_run$value, ABE = ABE_run$value),
  function(w) cor(as.numeric(w), as.numeric(ref)))
attr(report, "pearson_r") <- pearson_r
print(report, digits = 4)
print(pearson_r, digits = 4)
saveRDS(report, "benchmarks/stepstone_fulldata_homotopy_report.rds")

# Merge gate assertions (fatal; matches test-bench-gate.R).
stopifnot(report$max_err[report$config == "leafblower_AB"] <= 1.60e-3)
stopifnot(pearson_r[["AB"]] >= 0.99)
cat("Merge gate: PASS (AB max_err floor + Pearson-ref).\n")
if (report$max_err[report$config == "leafblower_ABE"] <= 1.60e-4)
  cat("Stretch: A4 PASS.\n")
else cat(sprintf("Stretch: A4 NOT met (%.3e); reported, not required.\n",
                 report$max_err[report$config == "leafblower_ABE"]))
# Reported-but-not-gating wall-time ratio
cat(sprintf("Wall-time ratio AB/autumn = %.2fx (reported only)\n",
            AB_run$wall_s / autumn_run$wall_s))
```

- [ ] **Step 2: Run benchmark.**

Run: `Rscript benchmarks/stepstone_fulldata_homotopy.R`
Expected: AB `max_err ≤ 1.60e-3` AND `pearson_r[["AB"]] ≥ 0.99`. If both pass, the merge floor is met.

If A3 fails, iterate on `homotopy_levels` (3–10), `homotopy_budget_p` (0.3–1.0), `residual_recheck_fraction` (internal; 0.05–0.5). Do NOT relax the gate. If no in-range setting passes, hypothesis is falsified — file beads HUMAN ticket and halt.

- [ ] **Step 3: Write bench-gate test.**

`tests/testthat/test-bench-gate.R`:

```r
test_that("stepstone-fulldata AB config meets merge floor + Pearson agreement", {
  skip_on_cran()
  rpt_path <- "benchmarks/stepstone_fulldata_homotopy_report.rds"
  skip_if(!file.exists(rpt_path))
  report <- readRDS(rpt_path)
  ab_row <- report[report$config == "leafblower_AB", ]
  pearson <- attr(report, "pearson_r")
  expect_lte(ab_row$max_err, 1.60e-3)
  expect_gte(pearson[["AB"]], 0.99)
})

test_that("kk1204 non-regression: max_err <= 1.322e-3 at max_iterations=500", {
  skip_on_cran()
  skip_if_not_installed("arrow")
  # Synthetic kk1204 regeneration (K=20, n=1M, max_weight=3). Matches the
  # memory-recorded config. Uses fixed seed for reproducibility.
  set.seed(1204)
  n <- 1000000
  K <- 20
  cats <- 5
  data <- as.data.frame(
    replicate(K, sample(seq_len(cats), n, replace = TRUE), simplify = FALSE))
  names(data) <- paste0("m", seq_len(K))
  target <- stats::setNames(
    lapply(seq_len(K), function(k) rep(1/cats, cats)),
    names(data))
  t0 <- Sys.time()
  w <- leafblower::harvest(
    data, target, max_weight = 3,
    method = "ieppa",
    max_iterations = 500,
    convergence = list(absolute = 1e-10),
    # Overlays default OFF — this test validates the baseline path unchanged.
    attach_weights = FALSE)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # max_err computation
  errs <- sapply(seq_along(target), function(k) {
    tab <- tapply(w, data[[names(target)[k]]], sum) / sum(w)
    max(abs(tab - target[[k]]))
  })
  cat(sprintf("kk1204 non-regression: max_err=%.3e, wall=%.1fs\n",
              max(errs), wall))
  expect_lte(max(errs), 1.322e-3)   # A5 non-regression
})
```

Run: `Rscript -e 'devtools::test_active_file("tests/testthat/test-bench-gate.R")'`
Expected: both PASS.

- [ ] **Step 4: Rate-slope numeric gate (non-blocking but logged).**

After WU-6 Step 2, compute the effective convergence slope from the final-config trajectory:

```r
probe <- utils::read.csv("benchmarks/ABE_trajectory_stepstone_small.csv")
fit   <- stats::lm(log(errRp) ~ log(iter), data = probe)
slope <- stats::coef(fit)["log(iter)"]
cat(sprintf("Measured slope (log errRp vs log iter) = %.3f\n", slope))
stopifnot(slope <= -0.75)   # must beat the O(t^-0.5) baseline
saveRDS(list(slope = slope, fit = fit),
        "benchmarks/final_rate_slope.rds")
```

The `-0.75` threshold is the midpoint between O(t^-0.5) (baseline, clamp-penalised) and O(t^-1) (Ghosal–Nutz unconstrained). Any slope steeper than `-0.75` is evidence that the overlays have partially restored the smooth-operator rate. This IS a gate (`stopifnot`).

Generate `benchmarks/ABE_trajectory_stepstone_small.csv` alongside the AB benchmark by running the same config once under env-var probe on stepstone-small.

- [ ] **Step 5: Record degenerate `M_cell = n` regression metric (informational).**

Run the existing degenerate `M_cell = n` micro-benchmark (from earlier plan §8c infrastructure) once; record `ratio = wall_time_ieppa / wall_time_raking`. Save to `benchmarks/degenerate_ratio.rds`. No assertion — regression awareness only.

- [ ] **Step 6: Commit WU-6.**

```bash
git add benchmarks/stepstone_fulldata_homotopy.R \
        benchmarks/stepstone_fulldata_homotopy_report.rds \
        benchmarks/ABE_trajectory_stepstone_small.csv \
        benchmarks/final_rate_slope.rds \
        benchmarks/degenerate_ratio.rds \
        tests/testthat/test-bench-gate.R
git commit -m "$(cat <<'EOF'
bench(ieppa): stepstone-fulldata merge gate + kk1204 non-regression

AB config meets max_err 1.60e-3 floor and Pearson r>=0.99 vs commit-8146894
reference. kk1204 non-regression test confirms max_err<=1.322e-3 at
max_iterations=500 with overlays off. Rate-slope gate (<=-0.75) confirms
overlays partially restore smooth-operator rate. Wall-time ratio and
degenerate-M_cell ratio reported-only.
EOF
)"
```

- [ ] **Step 7: Close beads ticket.**

---

### WU-7: Python parity + documentation

**Beads:** `bd create --title "WU-7: Python parity and docs update for overlays" --type task --priority 2`

**Files:** modify `python/leafblower/_bindings.cpp`, `python/leafblower/_harvest.py`, `python/leafblower/test_python.py`, `R/harvest.R` (roxygen block), `README.md` (or docs index).

- [ ] **Step 1: Extend pybind bindings with new C ABI fields.**

Mirror the changes in `src/leafblower.h` inside `python/leafblower/_bindings.cpp`. Expose enums `rk_scheduler_t` and `rk_eta_mode_t`. Accept homotopy / scheduler / eta-schedule kwargs from Python.

- [ ] **Step 2: Rebuild Python extension.**

Run: `pip install -e .` from the repo root.
Expected: clean build; `python -c "import leafblower; print(leafblower.__version__)"` succeeds.

- [ ] **Step 3: Extend `python/leafblower/_harvest.py`.**

Add new kwargs: `homotopy_levels=1, homotopy_start_factor=1.0, homotopy_end_factor=1.0, homotopy_budget_p=0.5, scheduler="round_robin", eta_schedule="fixed", eta_start=1.0, eta_end=1.0, eta_schedule_power=0.5`. Forward to the extension.

- [ ] **Step 4: Write pytest parity tests.**

`python/leafblower/test_python.py` additions:

```python
import os, tempfile
import numpy as np
import leafblower

def _make_fixture():
    rng = np.random.default_rng(1)
    n = 2000
    import pandas as pd
    return pd.DataFrame({
        "a": rng.integers(1, 4, n),
        "b": rng.integers(1, 3, n),
    })

def _target_fixture():
    return {"a": [0.3, 0.5, 0.2], "b": [0.6, 0.4]}

def test_homotopy_default_off_is_identical():
    data = _make_fixture()
    tg   = _target_fixture()
    base = leafblower.harvest(data, tg, max_weight=3, method="ieppa",
                              attach_weights=False)
    defaulted = leafblower.harvest(data, tg, max_weight=3, method="ieppa",
                                   homotopy_levels=1,
                                   scheduler="round_robin",
                                   eta_schedule="fixed",
                                   attach_weights=False)
    np.testing.assert_allclose(np.asarray(base), np.asarray(defaulted),
                               atol=1e-12)

def test_homotopy_reduces_errRp_stepstone_small():
    fx = "tests/testthat/fixtures/stepstone_small.parquet"
    tg_path = "tests/testthat/fixtures/stepstone_small_targets.rds"
    if not (os.path.exists(fx) and os.path.exists(tg_path)):
        import pytest
        pytest.skip("stepstone-small fixture missing")
    import pyarrow.parquet as pq
    import pyreadr
    data   = pq.read_table(fx).to_pandas()
    target = pyreadr.read_r(tg_path)[None]
    with tempfile.TemporaryDirectory() as d:
        base_csv = os.path.join(d, "base.csv")
        homo_csv = os.path.join(d, "homo.csv")
        os.environ["LBW_TRAJECTORY_AT"]  = "100,200,500"
        os.environ["LBW_TRAJECTORY_OUT"] = base_csv
        leafblower.harvest(data, target, max_weight=5, method="ieppa",
                           max_iterations=500,
                           convergence={"absolute": 1e-12},
                           attach_weights=False)
        os.environ["LBW_TRAJECTORY_OUT"] = homo_csv
        leafblower.harvest(data, target, max_weight=5, method="ieppa",
                           max_iterations=500,
                           convergence={"absolute": 1e-12},
                           homotopy_levels=5,
                           homotopy_start_factor=10,
                           homotopy_end_factor=1,
                           homotopy_budget_p=0.5,
                           attach_weights=False)
        import pandas as pd
        base_err = pd.read_csv(base_csv)["errRp"].iloc[-1]
        homo_err = pd.read_csv(homo_csv)["errRp"].iloc[-1]
    assert homo_err < 0.7 * base_err
```

- [ ] **Step 5: Run pytest.**

Run: `pytest python/leafblower/ -v`
Expected: all pass (including the two new parity tests).

- [ ] **Step 6: Update R roxygen for `harvest()`.**

Add `@param` lines for every new argument, with example values. Add `@references` citing Schmitzer 2019 (arxiv:1610.06519), Chizat 2024 (arxiv:2408.11620), Altschuler–Weed–Rigollet 2017 (arxiv:1705.09634), Tang 2024 (arxiv:2403.05054).

Run: `Rscript -e 'devtools::document()'`
Expected: `man/harvest.Rd` regenerated.

- [ ] **Step 7: Update `README.md` with a usage recipe.**

Add a "Convergence acceleration" section showing the recommended recipe:

```r
res <- leafblower::harvest(
  data, target, max_weight = 5,
  method = "ieppa",
  homotopy_levels = 5, homotopy_start_factor = 10,
  homotopy_end_factor = 1,
  scheduler = "greedy"
)
```

- [ ] **Step 8: Final CRAN check.**

Run: `R CMD build .`
Run: `R CMD check --as-cran leafblower_*.tar.gz`
Expected: 0 ERROR, 0 WARNING, NOTEs ≤ baseline.

- [ ] **Step 9: Commit WU-7.**

```bash
git add python/leafblower/_bindings.cpp python/leafblower/_harvest.py \
        python/leafblower/test_python.py R/harvest.R man/harvest.Rd \
        README.md
git commit -m "$(cat <<'EOF'
docs+bindings: Python parity + roxygen for overlay knobs

Python harvest() accepts the same homotopy_*/scheduler/eta_* kwargs as R;
parity tests mirror test-homotopy.R and test-priority-sweep.R. Roxygen cites
Schmitzer 1610.06519, Chizat 2408.11620, Altschuler-Weed-Rigollet 1705.09634,
Tang 2403.05054. README gains a convergence-acceleration recipe.
EOF
)"
```

- [ ] **Step 10: Close beads ticket.**

---

## Final verification

- [ ] `devtools::test()` green
- [ ] `pytest` green
- [ ] `R CMD check --as-cran` — 0 ERROR, 0 WARNING, NOTEs ≤ baseline
- [ ] `benchmarks/stepstone_fulldata_homotopy_report.rds` committed; AB max_err ≤ 1.60e-3 AND pearson r ≥ 0.99
- [ ] `benchmarks/final_rate_slope.rds` committed; slope ≤ -0.75
- [ ] `tests/testthat/test-bench-gate.R` kk1204 non-regression PASS
- [ ] All seven beads tickets closed

## Self-review

- **Spec coverage:** P-A WU-3; P-B WU-4; Tang-η WU-5; scaffolding WU-1; internal probe WU-2; falsification WU-2.5; benchmark + kk1204 WU-6; Python parity + docs WU-7. Merge gate items 1–6 mapped: (1) WU-1 Step 10; (2) WU-7 Step 5; (3) WU-5 Step 4 + WU-7 Step 8; (4) WU-6 Step 2–3; (5) WU-6 Step 3; (6) WU-6 Step 2 (Pearson check).
- **Placeholders:** none — every step has concrete code/commands. "TBD", "similar to Task N", "add appropriate" absent.
- **Type consistency:** `HomotopyConfigLbw`, `SchedulerConfigLbw`, `EtaScheduleConfigLbw` used consistently in C++; `homotopy_levels`, `homotopy_start_factor`, `homotopy_end_factor`, `homotopy_budget_p`, `scheduler`, `eta_schedule`, `eta_start`, `eta_end`, `eta_schedule_power` consistent across R and Python wrappers. `inner_solve_one_level`, `apply_single_margin_linear`, `apply_single_margin_log`, `compute_margin_errRp` named consistently across WU-3 and WU-4.
- **API names ground-truthed:** `target` (singular), `method`, `max_iterations`, `convergence = list(absolute = ...)`, `bounds_mode`, `attach_weights`, `weight_column`, `verbose` — verified against `R/harvest.R:37–62`. No fabricated arguments.
- **Internal vs public knobs:** `residual_recheck_fraction` internal-only; trajectory probe via env vars (`LBW_TRAJECTORY_AT`, `LBW_TRAJECTORY_OUT`) — no public API expansion beyond the three overlay knob groups.
- **Reported vs gating metrics:** wall-time ratio and degenerate-M_cell ratio reported-only; `max_err`, Pearson-r, kk1204 non-regression, rate-slope are gates.
