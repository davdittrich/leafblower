# iEPPA Convergence Overlays Evaluation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Evaluate four convergence overlays (P-A progressive bounds, P-B Greenkhorn, Tang-η dynamic schedule, Augmented Lagrangian ADMM) against current iEPPA baseline. Primary metric: marginal KL. Secondary: max_err, weight KL, wall.

**Architecture:**
- Task 0: Add `marginal_kl` metric to iEPPA (fix best_iter + stopping criterion). This is the prerequisite baseline.
- Task 1 (P-A): Pure R overlay — progressive max_weight annealing, no C++.
- Task 2 (P-B): C++ Greenkhorn — sort K-margin BCD loop by max-residual priority.
- Task 3 (Tang-η): Combined with P-A in R — tighten inner tol each level.
- Task 4 (Aug.Lagrangian): C++ ADMM — replace P1.1 hard clamp with soft enforcement via dual variable u[c].
- Task 5: Benchmark all on stepstone-fulldata; compare marginal KL, max_err.

**Metric definitions:**
- **Marginal KL** (primary): `Σ_k Σ_j t_kj × log(t_kj / (S_kj/W))` — calibration quality.
- **Weight KL** (secondary): from `result$convergence_used$objective` — sample-design preservation.
- **max_err** (secondary): max over k,j of `|S_kj/W - t_kj|`.

**Baseline (current ieppa on stepstone-fulldata):** marginal KL = 5.004e-3, max_err = 2.74e-3, wall = 2.3s.
**Target:** beat autumn max_err = 1.60e-3 and/or marginal KL < 5.004e-3.

**ieppa KL reference:** `tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds` → weight KL at best_iter = 3.008e-3.

**Tech Stack:** C++17 (`src/ieppa.cpp`, `src/types.hpp`), R (`R/harvest.R`, `benchmarks/stepstone_all_methods.R`). Build: `R CMD INSTALL --preclean . 2>&1 | tail -3`.

**Research basis:** `docs/investigations/2026-04-24-ieppa-accel-research.md`, `docs/superpowers/specs/2026-04-24-ieppa-speed-convergence-bounds-design.md`.

---

## File Map

| File | Task | Change |
|------|------|--------|
| `src/types.hpp` | 0 | Add `MARGINAL_KL` to CalibMetric enum |
| `src/ieppa.cpp` lines ~850-885 | 0 | Add marginal KL computation; update best_iter |
| `R/harvest.R` | 0 | Update ieppa default metric to `marginal_kl` |
| `benchmarks/stepstone_all_methods.R` | 0, 5 | Add marginal KL computation + overlay methods |
| `src/ieppa.cpp` lines ~620-630 | 2 | P-B: priority sort of BCD loop |
| `src/ieppa.cpp` lines ~740-775 | 4 | Aug.Lagrangian: replace clamp with ADMM |

---

## Task 0: marginal_kl metric + best_iter fix (leafblower-djvp)

**Why:** `best_iter` currently selects weights at minimum `errRp` (max marginal error). Marginal KL is the complete calibration loss — it measures ALL margins jointly, not just the worst. Stopping at min(marginal KL) gives the solution closest to target distributions.

Each method is already optimized wrt its own internal objective (weight KL for ieppa/sinkhorn, chi2 for greg, max_err for chebyshev). But best_iter selection and external stopping should use marginal KL as the shared calibration quality standard.

### Step 0.1: Add MARGINAL_KL to CalibMetric enum

Find `src/types.hpp`. Locate `enum class CalibMetric`. Add:
```cpp
MARGINAL_KL = 7,   // Σ_k Σ_j t_kj log(t_kj / achieved_kj) — calibration quality
```
after the existing entries. Verify with `grep -n "CalibMetric\|MARGINAL" src/types.hpp`.

### Step 0.2: Add marginal KL computation in convergence block

Read `src/ieppa.cpp` lines 840-885 to confirm the errRp computation block and BLOCK 1. The block computes per-category `S_kj = Σ_{c:g_k(c)=j} X[c]` and `errRp = max |S_kj/W - t_kj|`.

In the errRp loop (both linear and log path variants at lines ~851-873), add marginal KL accumulation. Find:
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
        res.max_error = errRp;
```

Replace the inner body with marginal KL accumulation alongside errRp:
```cpp
        } else {
            double marginal_kl = 0.0;
            for (int k = 0; k < st.K; k++) {
                for (int j = 0; j < st.cat_counts[k]; j++) {
                    const auto& cells = cells_by_margin_cat[cat_offset[k] + j];
                    double Skj = 0.0;
                    for (int c : cells) Skj += X[c];
                    double e = std::fabs(Skj / W_total - st.targets[k][j]);
                    if (e > errRp) errRp = e;
                    // Marginal KL: t log(t/achieved) — calibration quality
                    double tkj = st.targets[k][j];
                    double achieved = (W_total > 0.0) ? Skj / W_total : 0.0;
                    if (tkj > 1e-300 && achieved > 1e-300)
                        marginal_kl += tkj * std::log(tkj / achieved);
                }
            }
            res.marginal_kl_at_iter = marginal_kl;  // NEW field — added below
        }
        res.max_error = errRp;
```

Also add `marginal_kl = 0.0` declaration before the `if/else` block (same scope as `errRp`).

Apply identically to the OTHER errRp block (the linear path variant at ~lines 851-862) — find it by searching for the `use_linear` branch containing `S_lin`.

### Step 0.3: Add marginal_kl to IEPPAResult struct

Find `src/ieppa.hpp`. Add to `IEPPAResult`:
```cpp
double marginal_kl_at_iter = 0.0;  // latest iteration marginal KL
```

### Step 0.4: Update BLOCK 1 — add MARGINAL_KL best-iterate selection

Find lines ~876-885 (BLOCK 1: MAX_ERR best-iterate). After the existing MAX_ERR block, add:
```cpp
            // MARGINAL_KL best-iterate: select weights at minimum marginal calibration loss.
            // Preferred over MAX_ERR: measures all margins jointly, not just worst.
            if (st.convergence_cfg.metric == lbw::CalibMetric::MARGINAL_KL) {
                if (res.marginal_kl_at_iter < best_metric_seen) {
                    best_metric_seen = res.marginal_kl_at_iter;
                    best_iter_val    = iter;
                    for (int c = 0; c < ct.M_cell; c++)
                        W_best[c] = (X_init[c] > 0.0) ? X[c] / X_init[c] : 0.0;
                }
            }
```

Also add `MARGINAL_KL` to the convergence dispatch (find the metric comparison section that fires when `conv_metric < best_metric_seen` for non-MAX_ERR metrics). Add MARGINAL_KL handling alongside KL: use `res.marginal_kl_at_iter` as the convergence value.

### Step 0.5: Compile gate
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
```

### Step 0.6: Update harvest.R ieppa default to marginal_kl

Find `R/harvest.R` — the section that sets ieppa's default convergence metric. Change:
```r
# Before:
if (is.null(convergence$metric) && method == "ieppa") convergence$metric <- "max_err"
# After:
if (is.null(convergence$metric) && method == "ieppa") convergence$metric <- "marginal_kl"
```
(Adapt to the actual harvest.R code structure — search for `ieppa.*metric` or `max_err.*ieppa`.)

### Step 0.7: Add marginal KL to benchmark script

Update `benchmarks/stepstone_all_methods.R` to report marginal KL alongside other metrics. The existing `fit_metrics` function already computes KL (marginal). Rename the field to `marg_kl` to be unambiguous. Update the output format string:
```r
cat(sprintf("%-12s  wall=%6.1fs  iters=%4s  status=%d  max_err=%.4e  marg_kl=%.4e  DEFF=%.4f\n",
    method, wall, ..., m$max_err, m$KL, m$DEFF))
```

### Step 0.8: Run full test suite
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -4
```
Expected: FAIL 0 (or same pre-existing baseline), PASS ≥ 382.

### Step 0.9: Commit
```bash
git add src/types.hpp src/ieppa.hpp src/ieppa.cpp R/harvest.R benchmarks/stepstone_all_methods.R
git commit -m "feat(ieppa): marginal_kl convergence metric + best_iter selection

Add MARGINAL_KL to CalibMetric enum. Compute Σ_k Σ_j t_kj log(t_kj/achieved_kj)
in convergence block alongside errRp. Update best_iter to select minimum marginal
KL (calibration quality) rather than minimum max_err. Change ieppa default metric
to marginal_kl+improvement — consistent with calibration objective.

Marginal KL measures all margins jointly; max_err only worst. In the infeasible
case (tight bounds), marginal KL gives better solution selection.

Closes: leafblower-djvp"
```

---

## Task 1: P-A progressive bound tightening (leafblower-bjvc)

**Why:** Smooth objective near loose bounds (large max_weight) → fast convergence to a good solution. Progressively tighten bounds while warm-starting. Schmitzer 2019 epsilon-scaling + Chizat 2024 sqrt-t schedule theory.

**Files:** `R/harvest.R` only (pure R wrapper — no C++ changes).

### Step 1.1: Add `homotopy` parameter to harvest()

Find `R/harvest.R` function signature. Add:
```r
homotopy = NULL,   # list(n_levels=3, k_init=5) or NULL to disable
```

### Step 1.2: Implement P-A logic in harvest()

After argument parsing and before the `.Call()`, add:
```r
  # P-A: progressive bound tightening (Schmitzer/Chizat)
  if (!is.null(homotopy) && method == "ieppa") {
    n_levels <- homotopy$n_levels %||% 3L
    k_init   <- homotopy$k_init   %||% 5.0
    budget   <- max_iterations %/% n_levels
    w_warm   <- NULL

    for (l in seq_len(n_levels)) {
      # Chizat sqrt-t schedule: level l, factor = k_init^(1 - (l-1)/(n_levels-1))
      # Level 1: max_weight * k_init (loosest), Level n: max_weight * 1 (target)
      factor_l <- k_init ^ (1.0 - (l - 1L) / (n_levels - 1L))
      mw_l <- max_weight * factor_l
      r_l  <- harvest(df = df, target = target, method = method,
                       max_weight = mw_l, min_weight = min_weight,
                       max_iterations = budget,
                       start_weights = w_warm,
                       convergence = convergence,
                       bounds_mode = bounds_mode,
                       attach_weights = FALSE, verbose = 0L)
      w_warm <- as.numeric(r_l)
    }
    # Final level: run at target max_weight from warm start
    return(harvest(df = df, target = target, method = method,
                    max_weight = max_weight, min_weight = min_weight,
                    max_iterations = budget,
                    start_weights = w_warm,
                    convergence = convergence,
                    bounds_mode = bounds_mode,
                    attach_weights = attach_weights,
                    verbose = verbose, ...))
  }
```

### Step 1.3: Add P-A to benchmark script

In `benchmarks/stepstone_all_methods.R`, add a `run_pa()` function:
```r
run_pa <- function(n_levels=3L, k_init=5.0) {
  cat(sprintf("%-22s ...", paste0("ieppa-PA(k=",k_init,",L=",n_levels,")"))); flush.console()
  t0 <- proc.time()["elapsed"]
  r <- suppressWarnings(
    leafblower::harvest(df, tgt, method="ieppa",
      max_weight=5, min_weight=0, max_iterations=ITERS,
      homotopy=list(n_levels=n_levels, k_init=k_init),
      attach_weights=FALSE, verbose=0))
  wall <- proc.time()["elapsed"] - t0
  w <- as.numeric(r); m <- fit_metrics(w, df, tgt)
  res <- attr(r, "result")
  cat(sprintf("  wall=%6.1fs  iters=%4s  status=%d  max_err=%.4e  marg_kl=%.4e  DEFF=%.4f\n",
    wall, res$iterations, res$status, m$max_err, m$KL, m$DEFF))
  invisible(list(w=w, wall=wall, m=m))
}
```

### Step 1.4: Run benchmark
```bash
OMP_NUM_THREADS=1 Rscript -e '
source("benchmarks/stepstone_all_methods.R")
run("ieppa", max_iterations=5000L)      # baseline
run_pa(n_levels=3L, k_init=5.0)         # P-A (3 levels, k=5)
run_pa(n_levels=5L, k_init=10.0)        # P-A (5 levels, k=10)
' 2>&1 | grep -E "ieppa|PA" | head -5
```

### Step 1.5: Commit
```bash
git add R/harvest.R benchmarks/stepstone_all_methods.R
git commit -m "eval(ieppa): P-A progressive bound tightening overlay

R-only wrapper: N_levels outer stages with max_weight × k_init^(1-(l-1)/(L-1))
schedule (Chizat sqrt-t). Warm-start across levels. homotopy=list(n_levels=3,
k_init=5) activates P-A for method='ieppa'.

Closes: leafblower-bjvc"
```

---

## Task 2: P-B Greenkhorn priority scheduler (leafblower-axll)

**Why:** Round-robin K-margin BCD updates all margins equally. Greenkhorn orders by max-residual (highest errRp first) — concentrates compute on the most violated margin. O(K log K) overhead per sweep (negligible for K ≤ 20).

**Files:** `src/ieppa.cpp` only.

### Step 2.1: Add per_margin_err_prev storage

Near the `W_best` declaration (~line 280), add:
```cpp
    // P-B: per-margin errRp from previous sweep for priority ordering.
    // Initialized to 1/K (uniform priority) so first sweep is round-robin.
    std::vector<double> per_margin_err_prev(st.K, 1.0 / static_cast<double>(st.K));
```

### Step 2.2: Replace BCD loop with priority-sorted loop

Find the linear-path BCD loop (~line 622):
```cpp
            } else {
                for (int k = 0; k < st.K && !overflow_trip; k++) {
                    if (apply_single_margin_linear(k)) overflow_trip = true;
                }
            }
```

Replace with priority-sorted version:
```cpp
            } else {
                // P-B: Greenkhorn — update margins in decreasing errRp order.
                // Uses per_margin_err_prev from previous iteration's convergence block.
                std::vector<int> margin_order(st.K);
                std::iota(margin_order.begin(), margin_order.end(), 0);
                std::sort(margin_order.begin(), margin_order.end(),
                    [&per_margin_err_prev](int a, int b) {
                        return per_margin_err_prev[a] > per_margin_err_prev[b];
                    });
                for (int ki = 0; ki < st.K && !overflow_trip; ki++) {
                    if (apply_single_margin_linear(margin_order[ki])) overflow_trip = true;
                }
            }
```

Apply identical change to the log-path BCD loop (search for `apply_single_margin_log` in the round-robin loop).

### Step 2.3: Update per_margin_err_prev after convergence block

After the errRp+marginal_kl computation block (after `res.max_error = errRp` at line ~874), add:
```cpp
            // P-B: update per-margin errRp for next sweep's priority ordering.
            // Reuse cells_by_margin_cat (no extra pass needed beyond errRp loop).
            if (W_total > 0.0) {
                for (int k = 0; k < st.K; k++) {
                    double ek = 0.0;
                    for (int j = 0; j < st.cat_counts[k]; j++) {
                        const auto& cells = cells_by_margin_cat[cat_offset[k] + j];
                        double Skj = 0.0;
                        for (int c : cells) Skj += X[c];
                        double e = std::fabs(Skj / W_total - st.targets[k][j]);
                        if (e > ek) ek = e;
                    }
                    per_margin_err_prev[k] = ek;
                }
            }
```

Note: this adds one O(K × M_cell / K) = O(M_cell) pass per iteration. Same order as errRp computation already done. Acceptable.

### Step 2.4: Compile gate
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
```

### Step 2.5: Add to benchmark, run comparison
```bash
OMP_NUM_THREADS=1 Rscript -e '
source("benchmarks/stepstone_all_methods.R")
run("ieppa", max_iterations=5000L)   # baseline
# P-B is always-on — just run ieppa again and compare results
# (P-B changes the convergence path but not the API)
' 2>&1 | grep "ieppa" | head -3
```

Note: P-B is embedded in ieppa — no API change needed. Run the existing A1 test to confirm sinkhorn KL acceptance:
```bash
Rscript -e 'devtools::test_active_file("tests/testthat/test-calibration-solvers.R")' 2>&1 | grep -E "A1|FAIL|PASS" | head -5
```

### Step 2.6: Full test suite + commit
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -4
git add src/ieppa.cpp
git commit -m "eval(ieppa): P-B Greenkhorn priority scheduler in BCD loop

Replace round-robin K-margin loop with max-residual priority ordering.
per_margin_err_prev[k] tracks previous iteration's errRp_k; margins sorted
descending. First sweep is uniform (1/K init). O(K log K) overhead per
sweep. Combined with P-A and marginal_kl best_iter.

Closes: leafblower-axll"
```

---

## Task 3: Tang-η dynamic schedule (combined with P-A in R)

**No separate ticket** — Tang-η is a parameter of the P-A wrapper. Extend Task 1's `homotopy` argument with `eta_schedule`:

```r
# In the P-A loop, add per-level inner tolerance (Tang-η):
# Tighten convergence tolerance with level: tol_l = base_tol / l
tol_l <- (convergence$tol %||% 0.001) / l   # tighter each level
conv_l <- modifyList(convergence, list(tol = tol_l))

r_l <- harvest(df = df, target = target, method = method,
                max_weight = mw_l, min_weight = min_weight,
                max_iterations = budget,
                start_weights = w_warm,
                convergence = conv_l,   # ← per-level tol
                ...)
```

Add `tang_eta = TRUE` flag to `homotopy` list. When TRUE, apply per-level tol schedule. Update `run_pa()` benchmark function to accept `tang_eta` parameter.

---

## Task 4: Augmented Lagrangian ADMM capacity block (leafblower-rn0c)

**Why:** P1.1 Euclidean hard clamp disrupts KL monotone descent — it correctly projects each cell individually but breaks the global KL descent across BCD sweeps. ADMM soft enforcement accumulates constraint violations in a dual variable `u[c]`, allowing the KL descent to temporarily explore the super-bound region. At convergence u[c]→0 and the solution satisfies hard bounds exactly.

**Mathematical formulation (Douglas-Rachford splitting):**
```
x-update:   X̃[c]  = X_cur[c] / W[c]           (unconstrained update, trivial)
z-update:   z[c]  = clamp(X̃[c] + u[c], L, U)  (project onto capacity box)
dual-update: u[c] += X̃[c] - z[c]              (accumulate constraint violation)
consensus:  X[c] = z[c], W[c] = z[c]/X̃[c], X_cur[c] = z[c]
```

At convergence: `u[c] → 0`, so `z[c] = clamp(X̃[c], L, U)` — same fixed point as current hard clamp but better convergence path.

**Files:** `src/ieppa.cpp` only.

### Step 4.1: Add dual variable u

Near the `W` declaration (~line 137), add:
```cpp
    // Aug.Lagrangian ADMM: dual variable for capacity constraint.
    // u[c] accumulates X_tilde - z violations; converges to 0 at fixed point.
    // Reset at fallback (same as W reset).
    std::vector<double> u(ct.M_cell, 0.0);
```

### Step 4.2: Replace P1.1 clamp with ADMM update

Find the P1.1 capacity block (~lines 759-775):
```cpp
                double X_tilde_c = X_cur[c] / W[c];
                if (!std::isfinite(X_tilde_c) || X_tilde_c > kLinearOverflowTrip) { ... }
                if (X_tilde_c <= 0.0) { ... }
                double xc = std::clamp(X_tilde_c, L_cell[c], U_cell[c]);
                X[c] = xc;
                W[c] = xc / X_tilde_c;
                X_cur[c] = xc;
```

Replace the clamp line with ADMM:
```cpp
                double X_tilde_c = X_cur[c] / W[c];
                if (!std::isfinite(X_tilde_c) || X_tilde_c > kLinearOverflowTrip) { ... }
                if (X_tilde_c <= 0.0) { ... }
                // ADMM: z-update (project adjusted X_tilde onto capacity box)
                double z = std::clamp(X_tilde_c + u[c], L_cell[c], U_cell[c]);
                u[c] += X_tilde_c - z;   // dual update: accumulate violation
                double xc = z;
                X[c] = xc;
                W[c] = xc / X_tilde_c;
                X_cur[c] = xc;
```

### Step 4.3: Reset u at both fallback blocks

In fallback block 1 (after `std::fill(cell_lf...`, ~line 663):
```cpp
                std::fill(u.begin(), u.end(), 0.0);
```
In fallback block 2 (~line 773):
```cpp
                std::fill(u.begin(), u.end(), 0.0);
```

### Step 4.4: Compile gate
```bash
cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
```

### Step 4.5: Full test suite
```bash
Rscript -e 'devtools::test()' 2>&1 | tail -4
```
Expected: FAIL 0, PASS ≥ 382.

### Step 4.6: Commit
```bash
git add src/ieppa.cpp
git commit -m "eval(ieppa): Aug. Lagrangian ADMM capacity block

Replace hard P1.1 clamp with ADMM Douglas-Rachford: z=clamp(X_tilde+u,L,U),
u+=X_tilde-z. Dual variable u[c] accumulates constraint violations, allowing
temporary soft-boundary excursions during KL descent. Fixed point identical
to hard clamp (u→0 at convergence). Reduces KL disruption from capacity step.

Closes: leafblower-rn0c"
```

---

## Task 5: Benchmark comparison (leafblower-fvjy)

### Step 5.1: Run full comparison
```bash
OMP_NUM_THREADS=1 Rscript benchmarks/stepstone_all_methods.R > /tmp/overlay_bench.log 2>&1 &
tail -f /tmp/overlay_bench.log | grep -E "ieppa|sinkhorn|autumn|PA|PB|ADMM|status" --line-buffered
```

Methods to compare: ieppa (baseline), ieppa-PA, ieppa-PA+η, ieppa (with P-B embedded), ieppa-ADMM, ieppa-PA+ADMM, sinkhorn, autumn.

### Step 5.2: Collect results table

Expected output format:
```
method                  wall    iters  status  max_err   marg_kl   weight_kl  DEFF
ieppa (baseline)         2.3s    190      0    2.74e-3   5.00e-3   ?          1.949
ieppa-PA(k=5,L=3)       6.9s    570      0    ?         ?         ?          ?
ieppa-PA+eta             ?       ?        0    ?         ?         ?          ?
ieppa-ADMM               ?       ?        0    ?         ?         ?          ?
ieppa-PA+ADMM+PB+eta     ?       ?        0    ?         ?         ?          ?
sinkhorn                 1.7s     60      0    6.01e-3   6.73e-3   ?          1.934
autumn                  30.4s    —        0    1.60e-3   ?         —          1.995
```

### Step 5.3: Record weight KL for each method

Add weight KL measurement to benchmark:
```r
# Weight KL (what solvers minimize): only meaningful for ieppa/sinkhorn
# result$convergence_used$objective gives this when available
wkl <- if (!is.null(res$convergence_used) && !is.null(res$convergence_used$objective))
         res$convergence_used$objective else NA
```

### Step 5.4: Close ticket
```bash
bd close leafblower-fvjy --reason="Benchmark complete. Results: [fill in actual numbers]"
```

---

## Self-Review

**Spec coverage:**
| Spec item | Task |
|-----------|------|
| marginal_kl CalibMetric enum entry | Task 0.1 |
| marginal KL computation in convergence block | Task 0.2 |
| IEPPAResult.marginal_kl_at_iter field | Task 0.3 |
| best_iter selection at min(marginal_kl) | Task 0.4 |
| ieppa default metric → marginal_kl | Task 0.6 |
| P-A R wrapper with Chizat schedule | Task 1 |
| P-B Greenkhorn BCD priority sort | Task 2 |
| per_margin_err_prev update per iteration | Task 2.3 |
| Tang-η per-level tol tightening | Task 3 |
| u[c] dual variable declaration | Task 4.1 |
| ADMM P1.1 replacement | Task 4.2 |
| u reset at both fallback blocks | Task 4.3 |
| Benchmark with marginal KL primary | Task 5 |

**Placeholder scan:** None.

**Type consistency:** `marginal_kl_at_iter` as `double` everywhere. `per_margin_err_prev` as `std::vector<double>`. `u` as `std::vector<double>` with same size as W.

**Dependency order:** Task 0 must precede all others (baseline change). Tasks 1-4 are independent. Task 5 requires Tasks 0-4.
