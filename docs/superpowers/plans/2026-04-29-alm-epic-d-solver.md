# ALM ieppa_soft — Epic D: Solver Implementation (T7)

> **For agentic workers:** Use `superpowers:subagent-driven-development` to implement this plan.

**Goal:** Implement the ALM soft capacity block in `src/ieppa.cpp` so that `ieppa_soft` strictly outperforms `ieppa` on tight-bounds problems (T3 GREEN).
**Mechanism:** Linearized Newton ALM primal step + per-cell dual ascent + adaptive μ growth + final hard-clamp projection
**Forbidden:** Modifying `src/ieppa.cpp` before Epic A (fixture) and Epics B/C (plumbing — `src/types.hpp`, `src/leafblower.h`, `src/cell_table.*`, `src/r_bridge.cpp`, `src/c_api.cpp`, `R/harvest.R`) are committed. Do not reuse `alm_lambda`/`alm_mu` CalibState fields (those belong to lbfgsb).
**Audit:** T3 must be observed FAILING before any implementation line is written. T3 must be GREEN after Step 8. All other ALM tests (T4–T7, T9, T10) must be GREEN by end of task.

**Prerequisite chain:**

| Epic | Ticket | Status required |
|------|--------|-----------------|
| A — Fixture capture | leafblower-epic-a | CLOSED (committed) |
| B — Types/ABI (`src/types.hpp`, `src/leafblower.h`, `src/cell_table.*`) | leafblower-epic-b | CLOSED |
| C — Bridge/API (`src/r_bridge.cpp`, `src/c_api.cpp`, `R/harvest.R`) | leafblower-epic-c | CLOSED |
| **D — Solver** (`src/ieppa.cpp`) | **leafblower-epic-d / T7** | **THIS TASK** |

---

## RED Gate (must run BEFORE touching `src/ieppa.cpp`)

```bash
# Confirm pre-requisites landed:
grep -c 'capacity_mu' src/types.hpp      # must be >= 1
grep -c 'capacity_penalty' src/leafblower.h  # must be >= 1
grep -c 'capacity_mu_auto' src/cell_table.hpp  # must be >= 1

# Confirm ALM solver block is NOT yet in ieppa.cpp:
grep -c 'lambda_cell\|alm_active\|capacity_mu_adaptive' src/ieppa.cpp
# Expected: 0 (clean pre-D state)
```

Run the RED test and confirm it FAILS before touching any solver code:

```bash
Rscript -e "
  library(leafblower)
  set.seed(3); n <- 5000L
  df  <- data.frame(v1=factor(sample(5, n, TRUE)))
  tgt <- list(v1=setNames(c(0.4,0.3,0.15,0.1,0.05), as.character(1:5)))
  r_hard <- harvest(df, tgt, method='ieppa',
    max_weight=1.8, min_weight=0, max_iterations=500, attach_weights=FALSE)
  r_soft <- harvest(df, tgt, method='ieppa_soft',
    max_weight=1.8, min_weight=0, max_iterations=500, attach_weights=FALSE)
  me_hard <- attr(r_hard,'result')\$max_error
  me_soft <- attr(r_soft,'result')\$max_error
  cat('hard:', me_hard, '\nsoft:', me_soft, '\nidentical:', me_hard == me_soft, '\n')
"
# Expected: soft == hard == 0.0273... (identical — T3 fails)
```

In `tests/testthat/test-calibration-solvers.R`, write the T3 test block (or verify it already exists) and run it:

```r
# In tests/testthat/test-calibration-solvers.R
test_that("T3: ieppa_soft strictly better than ieppa on tight-bounds", {
  set.seed(3); n <- 5000L
  df <- data.frame(v1=factor(sample(5, n, TRUE)))
  tgt <- list(v1=setNames(c(0.4,0.3,0.15,0.1,0.05), as.character(1:5)))
  r_hard <- harvest(df, tgt, method="ieppa",
    max_weight=1.8, min_weight=0, max_iterations=500, attach_weights=FALSE)
  r_soft <- harvest(df, tgt, method="ieppa_soft",
    max_weight=1.8, min_weight=0, max_iterations=500, attach_weights=FALSE)
  me_hard <- attr(r_hard,"result")$max_error
  me_soft <- attr(r_soft,"result")$max_error
  expect_lt(me_soft, me_hard - 1e-6)
})
```

```bash
Rscript -e "devtools::test(filter='calibration-solvers')" 2>&1 | grep -E 'T3|FAIL|PASS|Error'
# Expected: T3 FAILS (me_soft == me_hard — ALM not yet active)
```

**Do not proceed to Step 1 until T3 FAILURE is confirmed.**

---

## Task 7: ALM Solver Implementation (single ticket, 9 sub-steps)

**Branch:** `fix/correctness-performance-2026-04-28`
**File:** `src/ieppa.cpp`
**Compile gate:** `R CMD INSTALL --preclean . 2>&1 | tail -5` — must show `* DONE (leafblower)` after each sub-step that introduces a structural change.

---

### Step 1 — Declare persistent ALM state before homotopy loop

**Location:** After the existing vector declarations block (around line 140 in `src/ieppa.cpp`, after `X_cur`, `W`, `lambda` arrays are declared, before the homotopy loop begins).

**Mechanism:** Single `const bool alm_active` guard — all ALM state is conditionally initialized. `lambda_cell` is a `std::vector<double>` sized `ct.M_cell`, zeroed at construction, lives on the heap (no stack-overflow risk even for M_cell = 1M).

```cpp
// ALM persistent state (ieppa_soft only)
const bool alm_active = st.use_admm_capacity && st.capacity_mu > 0.0;
const double capacity_mu_base   = st.capacity_mu;
double capacity_mu_adaptive     = capacity_mu_base;
double eta_i_current            = 1.0;
int    alm_violation_streak     = 0;
const bool tang_active = (st.eta_schedule.mode == EtaScheduleMode::TANG_DYNAMIC &&
                          N_levels > 1);
std::vector<double> lambda_cell;
if (alm_active) lambda_cell.assign(ct.M_cell, 0.0);
```

**Compile gate:**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
# Must show: * DONE (leafblower)
```

---

### Step 2 — Homotopy level entry: update `eta_i_current` + reset `lambda_cell`

**Location:** Replace the existing `alm_mu` scaling block at approximately lines 381–391 in `src/ieppa.cpp` (inside the homotopy level loop, at the block where `st.capacity_mu` is currently adjusted for the level).

**Exact match target — find the existing block:**

```bash
grep -n 'alm_mu\|capacity_mu\|alm_lambda' src/ieppa.cpp | head -20
# Identify the current per-level scaling block. It will be inside the for/while
# homotopy loop, before the inner Sinkhorn loop starts.
```

**Replacement:**

```cpp
if (alm_active) {
    if (tang_active) {
        const double scaled_frac = std::pow(frac, st.eta_schedule.schedule_power);
        eta_i_current = st.eta_schedule.eta_start *
            std::pow(st.eta_schedule.eta_end / st.eta_schedule.eta_start, scaled_frac);
        res.eta_final  = eta_i_current;
        st.capacity_mu = eta_i_current * capacity_mu_adaptive;
    } else {
        eta_i_current  = 1.0;
        st.capacity_mu = capacity_mu_adaptive;
    }
    alm_violation_streak = 0;
    std::fill(lambda_cell.begin(), lambda_cell.end(), 0.0);
}
```

**Unchanged components:** The existing non-ALM homotopy level logic (homotopy factor computation, sinkhorn sweep setup, convergence criterion setup) is untouched.

**Compile gate:**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

---

### Step 3 — Replace P1.1 capacity block with ALM primal+dual update

**Location:** Lines ~804–813 in `src/ieppa.cpp` — the existing ADMM conditional block that reads `X_tilde_c` and writes `X[c]`, `W[c]`, `X_cur[c]`. This is inside the per-cell capacity enforcement loop.

**Exact match target — verify the block before replacing:**

```bash
sed -n '800,820p' src/ieppa.cpp
# Confirm it contains the full ADMM conditional:
#   if (st.use_admm_capacity) {
#       double z = std::clamp(X_tilde_c + u[c], L_cell[c], U_cell[c]);
#       u[c] += X_tilde_c - z;
#       X[c] = z; W[c] = z / X_tilde_c; X_cur[c] = z;
#   } else {
#       double xc = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);
#       X[c] = xc; W[c] = xc / X_tilde_c; X_cur[c] = xc;
#   }
```

The current code is an ADMM conditional (NOT a bare hard-clamp). The replacement transforms the `if (st.use_admm_capacity)` branch into the ALM primal+dual update, gated on `alm_active` (= `use_admm_capacity && capacity_mu > 0`). The `else` branch (bare hard-clamp) is kept as the fallback for `alm_active=false` and `X_tilde_c <= 0`.

**Replacement — entire capacity block (replaces both ADMM branches):**

```cpp
if (alm_active && X_tilde_c > 0.0) {
    // ALM: linearized Newton step from s=1 (single-iter approximation, see spec §Algorithm)
    const double z   = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);
    const double rho = st.capacity_mu * X_tilde_c;  // dimensionless: KL-Hess / penalty-Hess ratio
    double X_alm = X_tilde_c * (1.0 - lambda_cell[c] + st.capacity_mu * z) / (1.0 + rho);
    // NaN/Inf/negative guard: revert to hard clamp on numerical failure
    if (!std::isfinite(X_alm) || X_alm <= 0.0) X_alm = z;
    X[c]     = X_alm;
    W[c]     = X_alm / X_tilde_c;
    X_cur[c] = X_alm;
    // Dual ascent: lambda += mu * (X - z); standard ALM update
    const double lambda_cap = 10.0 * capacity_mu_base * st.max_weight;
    lambda_cell[c] += st.capacity_mu * (X_alm - z);
    lambda_cell[c]  = std::clamp(lambda_cell[c], -lambda_cap, lambda_cap);
} else {
    // Hard clamp: ieppa default path (use_admm_capacity=false OR capacity_mu=0)
    // AND fallback when X_tilde_c <= 0
    const double xc = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);
    X[c] = xc; W[c] = xc / X_tilde_c; X_cur[c] = xc;
}
```

**Note on u[c]:** The ADMM `u[c]` dual vector from the original `if (st.use_admm_capacity)` branch is entirely replaced. The new `lambda_cell[c]` is the ALM dual and must NOT coexist with `u[c]` in the same block. Confirm:

```bash
grep -n 'u\[c\]' src/ieppa.cpp | head -20
# lambda_cell and u[c] must not both appear in the same capacity-update block
```

**Compile gate:**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

**Smoke test after this step (first meaningful behavioral check):**

```bash
Rscript -e "
  library(leafblower)
  set.seed(3); n <- 5000L
  df  <- data.frame(v1=factor(sample(5, n, TRUE)))
  tgt <- list(v1=setNames(c(0.4,0.3,0.15,0.1,0.05), as.character(1:5)))
  r_soft <- harvest(df, tgt, method='ieppa_soft',
    max_weight=1.8, min_weight=0, max_iterations=500, attach_weights=FALSE)
  cat('soft max_err:', attr(r_soft,'result')\$max_error, '\n')
  cat('ieppa_soft ran without crash\n')
"
```

---

### Step 3b — Update log-path capacity block with ALM logic

**Location:** Lines ~871–886 in `src/ieppa.cpp` — the log-path capacity block that reads `X_tilde[c]` (note: array, not scalar `X_tilde_c`) and writes `X[c]`, `W[c]`. This is a SECOND capacity enforcement site distinct from the linear-path block at Step 3.

**Exact match target — verify the block before replacing:**

```bash
sed -n '863,895p' src/ieppa.cpp
# Confirm it contains the log-path hard-clamp:
#   double xc = std::clamp(X_tilde[c], L_cell[c], U_cell[c]);
#   X[c] = xc;
#   if (X_tilde[c] > 0.0) {
#       W[c] = xc / X_tilde[c];
#   } else {
#       W[c] = 1.0;
#   }
```

**Replacement — apply identical ALM logic as Step 3:**

```cpp
if (alm_active && X_tilde[c] > 0.0) {
    // ALM: same linearized Newton step as linear-path (Step 3)
    const double z   = std::clamp(X_tilde[c], L_cell[c], U_cell[c]);
    const double rho = st.capacity_mu * X_tilde[c];
    double X_alm = X_tilde[c] * (1.0 - lambda_cell[c] + st.capacity_mu * z) / (1.0 + rho);
    // NaN/Inf/negative guard: revert to hard clamp on numerical failure
    if (!std::isfinite(X_alm) || X_alm <= 0.0) X_alm = z;
    X[c] = X_alm;
    W[c] = X_alm / X_tilde[c];
    // Dual ascent: same formula as Step 3
    const double lambda_cap = 10.0 * capacity_mu_base * st.max_weight;
    lambda_cell[c] += st.capacity_mu * (X_alm - z);
    lambda_cell[c]  = std::clamp(lambda_cell[c], -lambda_cap, lambda_cap);
} else {
    // Hard clamp fallback: log-path default AND X_tilde[c] <= 0
    const double xc = std::clamp(X_tilde[c], L_cell[c], U_cell[c]);
    X[c] = xc;
    W[c] = (X_tilde[c] > 0.0) ? xc / X_tilde[c] : 1.0;
}
```

**Note:** The log-path uses array `X_tilde[c]` (not the scalar `X_tilde_c` used in the linear-path). The `lambda_cell[c]` dual is shared across both paths — this is correct: dual state accumulated during log-path iterations informs the next linear-path entry. The Step 4 reset (on linear→log fallback) remains the authoritative reset point.

**Compile gate:**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
# Must show: * DONE (leafblower)
```

---

### Step 4 — Reset `lambda_cell` on linear→log fallback

**Location:** Find the two overflow fallback sites in `src/ieppa.cpp` that reset `u[]` (the existing ADMM dual). These are at approximately lines 692 and 708 — search:

```bash
grep -n 'overflow\|fallback\|u\[' src/ieppa.cpp | grep -i 'fill\|reset\|fallback' | head -10
```

At each such site (where `std::fill(u.begin(), u.end(), 0.0)` or equivalent appears), add immediately after:

```cpp
if (alm_active) std::fill(lambda_cell.begin(), lambda_cell.end(), 0.0);
```

**Rationale (spec §u[c] Reset Policy rev 2):** Linear-path X̃ scale ≠ log-path X̃ scale. Carrying linear-scale dual into log-path produces NaN or wrong-magnitude penalty. Reset accepts loss of dual signal at fallback; the fallback is itself a numerical reset event.

**Unchanged:** The existing `u[]` reset lines — keep them exactly as-is. Only ADD the `lambda_cell` reset immediately after.

**Compile gate:**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

---

### Step 5 — Adaptive μ growth block (after per-cell capacity loop)

**Location:** After the per-cell capacity loop closes (the `for (int c = 0; c < ct.M_cell; c++)` loop that contains Step 3's ALM block), inside the outer Sinkhorn iteration loop, before the convergence check.

```cpp
if (alm_active) {
    constexpr int    kAlmPersistenceThreshold = 5;
    constexpr double kAlmGrowthFactor         = 2.0;
    constexpr double kAlmMaxScale             = 1000.0;
    const double     mean_n_per_cell =
        static_cast<double>(st.n) / std::max(1, ct.M_cell);
    const double     tol_primal =
        0.01 * st.max_weight * mean_n_per_cell;

    double max_violation = 0.0;
    for (int c = 0; c < ct.M_cell; c++) {
        const double v = std::max(X[c] - U_cell[c], L_cell[c] - X[c]);
        if (std::isfinite(v)) max_violation = std::max(max_violation, v);
    }

    if (max_violation > tol_primal) {
        if (++alm_violation_streak >= kAlmPersistenceThreshold &&
            capacity_mu_adaptive < capacity_mu_base * kAlmMaxScale) {
            capacity_mu_adaptive *= kAlmGrowthFactor;
            st.capacity_mu = eta_i_current * capacity_mu_adaptive;
            res.alm_n_growth_events++;
            alm_violation_streak = 0;
            if (st.verbose >= 2) {
                char msg[128];
                std::snprintf(msg, sizeof(msg),
                    "[ieppa_soft] mu growth: adaptive=%.4e",
                    capacity_mu_adaptive);
                st.log(msg);
            }
        }
    } else {
        alm_violation_streak = 0;
    }
}
```

**Compile gate:**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

---

### Step 6 — Final projection (after while-loop, before result population)

**Location:** After the Sinkhorn while-loop exits, before `res.max_error` and other result fields are populated. This is the last ALM-specific block.

```cpp
if (alm_active) {
    // Hard-clamp all cells: guarantees IEEE-exact bound satisfaction at output
    for (int c = 0; c < ct.M_cell; c++)
        X[c] = std::clamp(X[c], L_cell[c], U_cell[c]);

    // Bounded clamp+rescale loop to restore sum-to-n (up to kRescaleTol * n drift)
    constexpr int    kMaxRescaleIters = 3;
    constexpr double kRescaleTol      = 1e-12;
    double prev_total = 0.0;
    for (int iter = 0; iter < kMaxRescaleIters; iter++) {
        double total = 0.0;
        for (int c = 0; c < ct.M_cell; c++) total += X[c];
        if (std::abs(total - prev_total) < kRescaleTol * st.n || total <= 0.0) break;
        prev_total = total;
        const double rescale = static_cast<double>(st.n) / total;
        for (int c = 0; c < ct.M_cell; c++) {
            X[c] = std::clamp(X[c] * rescale, L_cell[c], U_cell[c]);
        }
    }
    // Compute and store sum drift
    double final_total = 0.0;
    for (int c = 0; c < ct.M_cell; c++) final_total += X[c];
    res.alm_sum_drift = std::abs(final_total - static_cast<double>(st.n));
    if (res.alm_sum_drift > 1e-6 * st.n && st.verbose >= 1) {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
            "[ieppa_soft] final projection sum drift = %.2e", res.alm_sum_drift);
        st.log(msg);
    }
}
```

**Post-condition:** After this block, `L_cell[c] <= X[c] <= U_cell[c]` holds exactly for all c. `sum(X)` may differ from n by at most `kRescaleTol * n * kMaxRescaleIters + single-clamp-drift ≈ 3e-12 * n`.

**Compile gate:**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

---

### Step 7 — Populate `alm_*` result fields (at existing result population site)

**Location:** At the existing result population block (where `res.iterations`, `res.max_error`, etc. are set), add after those assignments:

```cpp
if (alm_active) {
    res.alm_capacity_mu_final = capacity_mu_adaptive;
    double max_dual = 0.0;
    for (int c = 0; c < ct.M_cell; c++)
        max_dual = std::max(max_dual, std::abs(lambda_cell[c]));
    res.alm_max_dual_norm = max_dual;
    // alm_n_growth_events: incremented in Step 5 inline — no assignment needed here
    // alm_sum_drift: set in Step 6 — no assignment needed here
}
```

**Note:** `res.alm_n_growth_events` and `res.alm_sum_drift` are set in Steps 5 and 6 respectively. For `method="ieppa"` (alm_active=false), these fields default to 0/0.0 from `rk_result_t` zero-initialization — verify this:

```bash
grep -n 'alm_capacity_mu_final\|alm_n_growth_events\|alm_max_dual_norm\|alm_sum_drift' src/leafblower.h
# Must show field declarations; confirm struct has zero-init (= 0 or = 0.0 initializers)
```

**Compile gate:**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -5
```

---

### Step 8 — Run T3 and confirm GREEN

```bash
# Install clean:
R CMD INSTALL --preclean . 2>&1 | tail -5
# Must show: * DONE (leafblower)

# Run T3:
Rscript -e "devtools::test(filter='calibration-solvers')" 2>&1 | tail -20
```

**Expected T3 output:**

```
[ PASS ] T3: ieppa_soft strictly better than ieppa on tight-bounds
```

**If T3 is still FAILING, diagnose in order:**

1. `grep -c 'lambda_cell' src/ieppa.cpp` — must be >= 3 (declare, use in primal, use in dual)
2. `grep -c 'alm_active' src/ieppa.cpp` — must be >= 4 (one per step)
3. Check that `alm_active` condition is true at runtime:
   ```bash
   Rscript -e "
     library(leafblower)
     set.seed(3); n <- 5000L
     df  <- data.frame(v1=factor(sample(5, n, TRUE)))
     tgt <- list(v1=setNames(c(0.4,0.3,0.15,0.1,0.05), as.character(1:5)))
     r   <- harvest(df, tgt, method='ieppa_soft',
       max_weight=1.8, min_weight=0, max_iterations=500,
       attach_weights=FALSE, verbose=2)
     res <- attr(r,'result')
     cat('capacity_mu_final:', res\$alm_capacity_mu_final, '\n')
     cat('n_growth_events:',   res\$alm_n_growth_events, '\n')
     cat('max_dual_norm:',     res\$alm_max_dual_norm, '\n')
   "
   ```
   If `alm_capacity_mu_final == 0`: the `use_admm_capacity` flag or `capacity_mu` field was not set by Epic C. Verify Epics B/C are correctly committed.
4. If `capacity_mu_final > 0` but T3 still fails: the initial `capacity_mu_base` (= `M_cell/n = 5/5000 = 0.001`) may be too small to escape the fixed point in 500 iterations. Try:
   ```bash
   Rscript -e "
     library(leafblower)
     set.seed(3); n <- 5000L
     df  <- data.frame(v1=factor(sample(5, n, TRUE)))
     tgt <- list(v1=setNames(c(0.4,0.3,0.15,0.1,0.05), as.character(1:5)))
     r   <- harvest(df, tgt, method='ieppa_soft',
       capacity_penalty=10.0,  # manual override, 10000x auto
       max_weight=1.8, min_weight=0, max_iterations=500, attach_weights=FALSE)
     me  <- attr(r,'result')\$max_error
     cat('soft max_err at manual penalty=10:', me, '\n')
   "
   ```
   If manual `capacity_penalty=10.0` passes T3, the adaptive growth in Step 5 has a bug (threshold too high, or streak counter not incrementing). Add `verbose=2` and check for growth event log lines.

**Do not proceed to Step 9 until T3 is GREEN.**

---

### Step 9 — Full test suite + commit

```bash
# Full test suite:
Rscript -e "devtools::test()" 2>&1 | tail -5
# FAIL count must equal pre-merge baseline (3 exactly — see spec §Acceptance Criteria #9)
# If FAIL > 3: investigate regressions before committing.
```

**Run the full ALM acceptance matrix:**

```bash
Rscript -e "devtools::test(filter='calibration-solvers')" 2>&1 | grep -E 'PASS|FAIL|SKIP'
# Expected PASS: T3, T4 (both sub-tests), T5 (both sub-tests), T7, T9, T10
# Expected SKIP: T6, T8 (require stepstone fixture, skip_on_cran)
# Expected FAIL: 0 new failures
```

**Run backward compat test explicitly:**

```bash
Rscript -e "
  library(leafblower)
  fixture_path <- 'tests/testthat/fixtures/ieppa_pre_alm_ref.rds'
  if (!file.exists(fixture_path)) stop('fixture missing — Epic A not complete')
  ref <- readRDS(fixture_path)
  r   <- harvest(ref\$df, ref\$tgt, method='ieppa',
                 max_weight=ref\$max_weight, min_weight=ref\$min_weight,
                 max_iterations=ref\$max_iterations,
                 convergence=ref\$convergence, attach_weights=FALSE)
  w   <- as.numeric(r)
  stopifnot(isTRUE(all.equal(w, ref\$weights, tolerance=1e-12)))
  cat('T9 backward compat: PASS\n')
"
```

**Only commit after all of the above pass:**

```bash
git add src/ieppa.cpp
git commit -m "$(cat <<'EOF'
feat(ieppa): ALM soft capacity — linearized Newton step, dual ascent, adaptive mu, final projection

Replace hard-clamp at P1.1 with per-cell ALM primal update (linearized Newton
from s=1, closed-form: X_alm = X_tilde*(1 - lambda + mu*z)/(1+rho)) plus raw
dual ascent (lambda += mu*(X-z), capped at 10*mu_base*max_weight). Declare
persistent capacity_mu_adaptive and lambda_cell vectors before homotopy loop;
reset lambda_cell on level transitions and linear->log fallback (rev-2 decision
5b'). Add adaptive mu growth (doubles every 5 consecutive violation iterations,
cap 1000x base). Add final clamp+rescale projection (3 iterations) to guarantee
IEEE-exact bounds at output with bounded sum-drift. Populate alm_* result fields.

T3 now passes: ieppa_soft max_err < ieppa max_err by >=1e-6 on tight-bounds
5-category problem (max_weight=1.8). Method='ieppa' produces bit-identical
weights to pre-ALM fixture (T9 backward compat confirmed).
EOF
)"
```

---

## Acceptance Criteria (Epic D / T7)

| # | Criterion | Command |
|---|-----------|---------|
| 1 | T3 GREEN | `devtools::test(filter='calibration-solvers')` — T3 PASS |
| 2 | T4 GREEN | capacity_penalty=NULL routes to auto; invalid inputs error |
| 3 | T5 GREEN | Final weights in [L,U] exactly; sum drift < 1e-6*n; differ from ieppa |
| 4 | T7 GREEN | Adaptive growth fires on adversarial start (capacity_penalty=1e-6) |
| 5 | T9 GREEN | ieppa produces bit-identical weights to pre-ALM fixture |
| 6 | T10 GREEN | capacity_penalty + non-ieppa_soft method emits warning |
| 7 | No regressions | `devtools::test()` FAIL == 3 (pre-merge baseline, not ≤) |
| 8 | Compile clean | `R CMD INSTALL --preclean .` exits 0, no warnings in ALM block |
| 9 | ALM diagnostics visible | `attr(r,"result")$alm_capacity_mu_final > 0` for ieppa_soft; 0 for ieppa |
| 10 | No lbfgsb collision | `grep 'alm_lambda\|alm_mu' src/ieppa.cpp` returns 0 (using capacity_mu only) |

---

## Risk Register

| Risk | Mitigation |
|------|-----------|
| `alm_active=false` at runtime (use_admm_capacity not set by Epic C) | Verify Epic C committed; check `capacity_mu_final==0` diagnostic; add verbose=2 logging |
| Auto μ too small to escape fixed point in 500 iters (T3 slow convergence) | Adaptive growth catches persistent violations; if not, T3 diagnosis protocol (Step 8 troubleshoot) identifies and caps at manual penalty |
| `lambda_cell` size mismatch if M_cell changes mid-solve | `lambda_cell.assign(ct.M_cell, 0.0)` at solver entry; M_cell is const for the solve lifetime |
| ABI tripwire fires (EXPECTED_RK_RESULT_BYTES mismatch) | Epic B must update to 480; this epic only reads fields it does not declare |
| Convergence degradation: final rescale moves weights | Rescale loop bounded to 3 iters; T5 sum-drift test catches any drift > 1e-6*n |
| ieppa backward compat broken by stray side-effect | T9 fixture test; `alm_active=false` for ieppa (use_admm_capacity=false) means all new blocks are no-ops |
