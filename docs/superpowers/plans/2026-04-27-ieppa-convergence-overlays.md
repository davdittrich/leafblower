# iEPPA Convergence Overlays Evaluation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Evaluate four convergence overlays against current iEPPA baseline. Primary metric: marginal KL. Secondary: max_err, weight KL, wall.

**Key finding during planning:** `scheduler="greedy"` (P-B), `eta_schedule="tang_dynamic"` (Tang-η), and `homotopy_levels` (P-A) are ALL already implemented in C++ (ieppa.cpp) and exposed in harvest.R. Only Augmented Lagrangian is genuinely new C++ work.

**Architecture:**
- Task 0: Add `marginal_kl` metric + fix best_iter selection — prerequisite baseline change.
- Task 1: Benchmark existing overlays (P-A via `homotopy_levels`, P-B via `scheduler="greedy"`, Tang-η via `eta_schedule="tang_dynamic"`) — R-only parameter grid, no C++ changes.
- Task 2: Implement Augmented Lagrangian ADMM capacity block — only new C++ overlay.
- Task 3: Full benchmark comparison of all combinations.

**Metric definitions:**
- **Marginal KL** (primary): `Σ_k Σ_j t_kj × log(t_kj / (S_kj/W))` — calibration quality.
- **Weight KL** (secondary): from `result$convergence_used$objective` — sample-design preservation.
- **max_err** (secondary): `max_k,j |S_kj/W - t_kj|`.

**Baseline (current ieppa, default params, stepstone-fulldata):** marginal KL = 5.004e-3, max_err = 2.74e-3, wall = 2.3s.
**ieppa weight KL at best_iter (A1 fixture):** 3.008e-3.
**Autumn reference:** max_err = 1.60e-3, wall = 30.4s.
**Target:** marginal KL < 5.004e-3 and/or max_err ≤ 1.60e-3.

**Tech Stack:** C++17 (`src/ieppa.cpp`, `src/ieppa.hpp`, `src/types.hpp`), R (`R/harvest.R`, `benchmarks/stepstone_all_methods.R`). Build: `R CMD INSTALL --preclean . 2>&1 | tail -3`.

**Existing overlay APIs (already wired):**
- `homotopy_levels=3, homotopy_start_factor=5.0, homotopy_end_factor=1.0` — P-A
- `scheduler="greedy"` — P-B (max-residual priority BCD)
- `eta_schedule="tang_dynamic", eta_schedule_power=0.5` — Tang-η (requires homotopy_levels > 1)

---

## File Map

| File | Task | Change |
|------|------|--------|
| `src/types.hpp` | 0 | Add `MARGINAL_KL = 7` to CalibMetric enum |
| `src/ieppa.hpp` | 0 | Add `marginal_kl_at_iter` field to IEPPAResult |
| `src/ieppa.cpp` lines ~850-885 | 0 | Compute marginal KL; add MARGINAL_KL best-iter block |
| `R/harvest.R` | 0 | Change ieppa default metric from `"kl"` to `"marginal_kl"` |
| `benchmarks/stepstone_all_methods.R` | 0, 1, 3 | Add marginal KL computation + overlay grid |
| `src/ieppa.cpp` lines ~137, ~759-775, ~663, ~773 | 2 | Augmented Lagrangian: add u[c] + ADMM P1.1 |

---

## Task 0: marginal_kl convergence metric + best_iter fix (leafblower-djvp)

**Why:** Best_iter currently selects minimum `errRp` (max marginal error = L∞). Marginal KL measures all margins jointly (L2-style). In the infeasible case, marginal KL identifies the solution with best overall calibration quality, not just best worst-margin. Default stopping criterion `"kl"` currently tracks WEIGHT KL — change to `"marginal_kl"` for calibration-quality-driven stopping.

### Step 0.1: Add MARGINAL_KL to CalibMetric enum

Read `src/types.hpp`. Find `enum class CalibMetric`. Confirm existing entries end at value ≤ 6. Add:
```cpp
MARGINAL_KL = 7,   // Σ_k Σ_j t_kj log(t_kj / achieved_kj) — calibration quality
```

### Step 0.2: Add marginal_kl_at_iter to IEPPAResult

Read `src/ieppa.hpp`. Find `struct IEPPAResult`. Add:
```cpp
double marginal_kl_at_iter = 0.0;  // marginal KL at current iteration (updated each outer iter)
```

### Step 0.3: Add marginal KL computation in convergence block

Read `src/ieppa.cpp` lines 840-885. Find the log-path errRp block (the `} else {` branch at ~line 863 that uses `cells_by_margin_cat` to compute `Skj += X[c]`):

```cpp
        } else {
            for (int k = 0; k < st.K; k++) {
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    const auto& cells = cells_by_margin_cat[cat_offset[k] + j];
                    double Skj = 0.0;
                    for (int c : cells) Skj += X[c];
                    double e = std::fabs(Skj / W_total - st.targets[k][j]);
                    if (e > errRp) errRp = e;
                }
            }
        }
```

Replace the inner body to accumulate marginal KL alongside errRp:
```cpp
        } else {
            double marg_kl = 0.0;
            for (int k = 0; k < st.K; k++) {
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    const auto& cells = cells_by_margin_cat[cat_offset[k] + j];
                    double Skj = 0.0;
                    for (int c : cells) Skj += X[c];
                    double e = std::fabs(Skj / W_total - st.targets[k][j]);
                    if (e > errRp) errRp = e;
                    double tkj = st.targets[k][j];
                    double skj = (W_total > 0.0) ? Skj / W_total : 0.0;
                    if (tkj > 1e-300 && skj > 1e-300)
                        marg_kl += tkj * std::log(tkj / skj);
                }
            }
            res.marginal_kl_at_iter = marg_kl;
        }
```

Also declare `double marg_kl = 0.0;` before the `if/else` block and add `res.marginal_kl_at_iter = marg_kl;` in the linear-path branch too (search for the other errRp block that uses `S_lin` or `g_per_cell` — add the same KL accumulation using the same per-cell X[c] data).

### Step 0.4: Add MARGINAL_KL best-iterate tracking

After the existing MAX_ERR best-iterate block (lines ~876-885):
```cpp
            // MARGINAL_KL best-iterate: minimum marginal calibration loss.
            if (st.convergence_cfg.metric == lbw::CalibMetric::MARGINAL_KL) {
                if (res.marginal_kl_at_iter < best_metric_seen) {
                    best_metric_seen     = res.marginal_kl_at_iter;
                    best_iter_val        = iter;
                    for (int c = 0; c < ct.M_cell; c++)
                        W_best[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
                }
            }
```

Also add MARGINAL_KL to the convergence-firing dispatch (the section that computes `conv_metric` and calls `check_convergence`): when metric == MARGINAL_KL, use `res.marginal_kl_at_iter` as the conv_metric value.

### Step 0.5: Update harvest.R ieppa default

Find `R/harvest.R` lines ~163-166 where ieppa's default metric is set. Current code:
```r
if (is.null(convergence$metric) && method == "ieppa") convergence$metric <- "kl"
```
Change to:
```r
if (is.null(convergence$metric) && method == "ieppa") convergence$metric <- "marginal_kl"
```

Also add `"marginal_kl"` to the valid metric names vector in harvest.R (search for where other metric names like `"max_err"`, `"kl"` are listed for validation).

### Step 0.6: Compile gate
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
```

### Step 0.7: Full test suite
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -4
```
Expected: FAIL 0, PASS ≥ 382.

### Step 0.8: Update benchmark script to output marginal KL separately

The `fit_metrics()` function in `benchmarks/stepstone_all_methods.R` already computes the KL field (which IS marginal KL from `tk*log(tk/Sr)`). Rename it to `marg_kl` for clarity. Also add weight KL from `result$convergence_used$objective`:

```r
fit_metrics <- function(w, df, tgt, res=NULL) {
  # ... existing computation ...
  list(max_err=max_err, marg_kl=KL, chi2=chi2, L1=L1,
       weight_kl = if (!is.null(res) && !is.null(res$convergence_used$objective))
                     res$convergence_used$objective else NA,
       DEFF=DEFF, ESS=ESS, wmin=wmin, wmed=wmed, wmax=wmax)
}
```

Update `run()` to pass `res` to `fit_metrics`.

### Step 0.9: Commit
```bash
git add src/types.hpp src/ieppa.hpp src/ieppa.cpp R/harvest.R benchmarks/stepstone_all_methods.R
git commit -m "feat(ieppa): marginal_kl convergence metric + best_iter selection

Add MARGINAL_KL=7 to CalibMetric. Compute Σ_k Σ_j t_kj log(t_kj/achieved_kj)
in convergence block. Best_iter now selects minimum marginal KL (calibration
quality) not max_err (L-inf). Default metric changed from weight KL to
marginal_kl+improvement — consistent with calibration objective as loss.

Marginal KL: calibration quality loss (all margins jointly, L2-style).
Weight KL: sample-design preservation (solver's internal objective).

Closes: leafblower-djvp"
```

---

## Task 1: Benchmark existing overlays — P-A, P-B, Tang-η (no new C++)

**Key insight:** All three overlays already exist:
- `homotopy_levels=3, homotopy_start_factor=5.0` → P-A (progressive max_weight)
- `scheduler="greedy"` → P-B (max-residual priority BCD, Greenkhorn-style)
- `eta_schedule="tang_dynamic"` → Tang-η (dynamic ALM schedule, requires homotopy)

**Files:** `benchmarks/stepstone_all_methods.R` only.

### Step 1.1: Add overlay methods to benchmark script

Extend `benchmarks/stepstone_all_methods.R` to include overlay combinations:

```r
cat("\n=== ieppa overlay combinations ===\n")

# P-B: greedy scheduler only
run("ieppa",
    scheduler = "greedy",
    max_iterations = ITERS)

# P-A: progressive bounds (3 levels, k=5)
run("ieppa",
    homotopy_levels = 3L,
    homotopy_start_factor = 5.0,
    homotopy_end_factor = 1.0,
    max_iterations = ITERS)

# P-A + Tang-eta
run("ieppa",
    homotopy_levels = 3L,
    homotopy_start_factor = 5.0,
    homotopy_end_factor = 1.0,
    eta_schedule = "tang_dynamic",
    eta_schedule_power = 0.5,
    max_iterations = ITERS)

# P-A + P-B + Tang-eta (combined)
run("ieppa",
    scheduler = "greedy",
    homotopy_levels = 3L,
    homotopy_start_factor = 5.0,
    homotopy_end_factor = 1.0,
    eta_schedule = "tang_dynamic",
    max_iterations = ITERS)
```

### Step 1.2: Run and record results
```bash
OMP_NUM_THREADS=1 Rscript benchmarks/stepstone_all_methods.R > /tmp/overlay_bench.log 2>&1 &
tail -f /tmp/overlay_bench.log | grep -E "ieppa|sinkhorn|autumn" --line-buffered
```

### Step 1.3: Commit
```bash
git add benchmarks/stepstone_all_methods.R
git commit -m "bench(ieppa): add P-A, P-B, Tang-eta overlay combinations to benchmark

Uses existing harvest.R APIs (scheduler, homotopy_levels, eta_schedule).
No new C++ — all three overlays already implemented."
```

---

## Task 2: Augmented Lagrangian ADMM capacity block (leafblower-rn0c)

**Why:** P1.1 Euclidean hard clamp is correct at each cell individually but disrupts KL monotone descent globally across BCD sweeps. ADMM soft enforcement maintains a dual variable `u[c]` that accumulates constraint violations, allowing temporary soft-boundary excursions that reduce interference with KL descent. Fixed point is identical to the hard clamp (u[c]→0 at convergence).

**Mathematical form (Douglas-Rachford splitting):**
```
X̃[c]  = X_cur[c] / W[c]              x-update (trivial: prox of KL = identity)
z[c]  = clamp(X̃[c] + u[c], L, U)   z-update  (project onto capacity box)
u[c] += X̃[c] - z[c]                dual-update (accumulate violation)
X[c]  = z[c],  W[c] = z[c]/X̃[c]    consensus
```

At convergence: u[c]→0 → z[c] = clamp(X̃[c], L, U). Identical fixed point to current, better convergence path.

### Step 2.1: Add dual variable u

After line ~137 (near W declaration), add:
```cpp
    // ADMM dual variable for capacity constraint.
    // u[c] accumulates X_tilde - z violations; converges to 0 at fixed point.
    std::vector<double> u(ct.M_cell, 0.0);
```

### Step 2.2: Replace P1.1 clamp with ADMM update

Read lines ~759-775. Find:
```cpp
                double xc = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);
                X[c] = xc;
                W[c] = xc / X_tilde_c;
                X_cur[c] = xc;
```

Replace:
```cpp
                double z = std::clamp(X_tilde_c + u[c], L_cell[c], U_cell[c]);
                u[c] += X_tilde_c - z;   // dual update: accumulate violation
                X[c] = z;
                W[c] = z / X_tilde_c;
                X_cur[c] = z;
```

### Step 2.3: Reset u in both fallback blocks

In fallback block 1 (~line 663, after cell_lf reset):
```cpp
                std::fill(u.begin(), u.end(), 0.0);
```

In fallback block 2 (~line 773, after cell_lf reset):
```cpp
                std::fill(u.begin(), u.end(), 0.0);
```

Both resets must appear BEFORE the W reset in their respective blocks (u is coupled to W).

### Step 2.4: Compile gate
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
```

### Step 2.5: Full test suite
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -4
```
Expected: FAIL 0, PASS ≥ 382.

### Step 2.6: Add ADMM variants to benchmark

In `benchmarks/stepstone_all_methods.R`:
```r
cat("\n=== ieppa ADMM capacity block ===\n")
run("ieppa", max_iterations = ITERS)                          # baseline (now ADMM)
run("ieppa", scheduler="greedy",
    homotopy_levels=3L, homotopy_start_factor=5.0,
    eta_schedule="tang_dynamic", max_iterations=ITERS)        # all combined + ADMM
```

Note: ADMM is now the default for all ieppa runs (always-on after this task). The "baseline" after Task 2 is ieppa+ADMM.

### Step 2.7: Commit
```bash
git add src/ieppa.cpp
git commit -m "feat(ieppa): Aug. Lagrangian ADMM capacity block

Replace hard P1.1 clamp with ADMM Douglas-Rachford: z=clamp(X_tilde+u,L,U),
u+=X_tilde-z. u[c] accumulates constraint violations across iterations,
allowing KL descent to temporarily explore near-bound region. Fixed point
identical to hard clamp (u→0). 8MB extra storage (u vector = same as W).
Both fallback blocks reset u before W.

Closes: leafblower-rn0c"
```

---

## Task 3: Benchmark comparison (leafblower-fvjy)

### Step 3.1: Run full comparison (after Tasks 0-2)
```bash
OMP_NUM_THREADS=1 Rscript benchmarks/stepstone_all_methods.R > /tmp/final_bench.log 2>&1
cat /tmp/final_bench.log | grep -E "^(ieppa|sinkhorn|autumn|===)" | head -30
```

### Step 3.2: Expected results table

Collect:

| Method | Wall | max_err | marginal KL | weight KL | DEFF |
|--------|------|---------|-------------|-----------|------|
| ieppa (baseline, Task 0) | 2.3s | 2.74e-3 | 5.00e-3 | 3.01e-3 | 1.949 |
| ieppa + greedy (P-B) | ? | ? | ? | ? | ? |
| ieppa + homotopy (P-A) | ? | ? | ? | ? | ? |
| ieppa + P-A + Tang-η | ? | ? | ? | ? | ? |
| ieppa + P-A+P-B+Tang-η | ? | ? | ? | ? | ? |
| ieppa + ADMM (Task 2) | ? | ? | ? | ? | ? |
| ieppa + ADMM + P-A+P-B+Tang-η | ? | ? | ? | ? | ? |
| sinkhorn | 1.7s | 6.01e-3 | 6.73e-3 | ? | 1.934 |
| autumn | 30.4s | 1.60e-3 | ? | — | 1.995 |

### Step 3.3: Close tickets
```bash
bd close leafblower-bjvc leafblower-axll leafblower-fvjy --reason="Evaluated via benchmark run. Results: [fill actual numbers]"
```

---

## Self-Review

**Spec coverage:**

| Item | Task |
|------|------|
| `MARGINAL_KL = 7` in CalibMetric | 0.1 |
| `marginal_kl_at_iter` in IEPPAResult | 0.2 |
| Marginal KL in convergence block (both paths) | 0.3 |
| MARGINAL_KL best-iterate block | 0.4 |
| ieppa default `"kl"` → `"marginal_kl"` | 0.5 |
| P-A via existing `homotopy_levels` API | 1 |
| P-B via existing `scheduler="greedy"` API | 1 |
| Tang-η via existing `eta_schedule="tang_dynamic"` | 1 |
| u[c] ADMM dual variable declaration | 2.1 |
| ADMM P1.1 replacement | 2.2 |
| u reset before W in both fallback blocks | 2.3 |
| Benchmark with marginal KL primary + weight KL | 3 |

**Placeholder scan:** None.

**Type consistency:** `marginal_kl_at_iter` as `double`. `u` as `std::vector<double>`, same size as W and cell_lf.

**Dependency order:** Task 0 before all others. Tasks 1 and 2 independent after Task 0. Task 3 after Tasks 0-2.
