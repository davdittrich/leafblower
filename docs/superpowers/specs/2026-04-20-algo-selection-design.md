# Algorithm Selection Benchmark Design

**Date**: 2026-04-20  
**Status**: Revised (gate iteration 1 → 2)  
**Scope**: Unconstrained calibration only (`max_weight=Inf`, `min_weight=0`)

---

## Goal

Find the 1.2× iso-contour in `(log_complexity, log_tol_abs)` space that separates "L-BFGS-B is ≥20% faster than iEPPA" from "iEPPA wins". Use the contour to hand-tune two constants in `src/c_api.cpp`.

**Win criterion**: L-BFGS-B wins at a grid point if `t_iEPPA / t_LBFGSB ≥ 1.2`.

---

## Scope

### In scope
- Unconstrained problems: `max_weight=Inf`, `min_weight=0`
- Primary GP axes: `log₁₀(complexity = n × Σcat_counts)` and `log₁₀(tol_abs)`
- Secondary K-stability check: run two additional fixed-grid sweeps at K=3 and K=18 to validate that the contour does not shift significantly with margin count (see Known Limitation below)
- Algorithm: Bayesian level set estimation (Straddle acquisition) with a GP surrogate

### Out of scope
- Box-constrained problems — iEPPA is unconditionally correct there; no benchmark needed
- Automating the constant update — human reads the plot, updates `c_api.cpp`

### Known limitation: complexity axis conflation

`complexity = n × Σcat_counts` collapses two distinct algorithmic cost drivers:

- **iEPPA** inner sweep is O(n × K): dominated by `n` (number of observations per IPF sweep). `Σcat_counts` contributes only negligibly (bucket accumulation O(cat_counts[k]) per margin).
- **L-BFGS-B** gradient pass is O(n × K) — same as iEPPA — but the L-BFGS 2-loop adds O(m × Σcat_counts) per outer step (m=10 history pairs, Σcat_counts = dual dimension). For large Σcat_counts (e.g., 835 in Stepstone), this term becomes significant relative to the gradient pass.

At equal `complexity`, a problem with large n and small Σcat_counts penalises both solvers equally on the gradient pass, but barely penalises L-BFGS-B on the 2-loop. A problem with small n and large Σcat_counts is cheap on the gradient pass but increasingly expensive for L-BFGS-B's 2-loop. The GP will therefore learn a boundary specific to the n/Σcat_counts regime of the data generated in the benchmark.

**Mitigation**: the main 2D GP sweep fixes K=9 with uniform `cats_per_margin` (Stepstone regime). The K-stability check (K ∈ {3, 18}) tests whether the contour shifts with margin count at fixed complexity. If the contour is stable across K, the constants are general; if it shifts, the constants are labelled as Stepstone-calibrated and the comment in `c_api.cpp` documents the scope.

**If the K-stability check shows contour shift > 0.5 log-units**: do not update constants; file a follow-up issue for a 3D sweep over `(log_complexity, log_tol, log(Σcat_counts/n))`.

---

## Deliverables

| File | Purpose |
|---|---|
| `benchmarks/algo_selection_benchmark.R` | LSE benchmark script + K-stability check |
| `benchmarks/algo_selection_results.rds` | Raw timing data + GP object (reproducibility) |
| `benchmarks/algo_selection_contour.pdf` | GP posterior mean + 1.2× contour |
| `benchmarks/algo_selection_uncertainty.pdf` | GP posterior σ (shows fit quality) |
| `benchmarks/algo_selection_k_stability.pdf` | Contour at K=3 vs K=9 vs K=18 overlay |
| `tests/testthat/test-algo-selection.R` | Regression guard for the updated constants |
| `src/c_api.cpp` | Updated constants (manual commit after plot review) |

---

## Benchmark Harness

### Input space (2D GP)
- `x₁ = log₁₀(complexity)` ∈ [4, 7.7] — covers complexity 10K to 50M
- `x₂ = log₁₀(tol_abs)` ∈ [−6, −3] — covers tol 1e-6 to 1e-3

**Lower bound rationale**: the existing benchmark (`stepstone_benchmark.R`, line 18–19) documents that iEPPA stalls beyond tol ≈ 1e-4 on the 9-margin overlapping system due to cyclic IPF oscillation. Below ~1e-5, `t_iEPPA` reflects a convergence failure, not a genuine speed comparison. The lower bound is set to 1e-6 (one decade of safety margin) to avoid contaminating the GP with failure-mode observations. If users request tol < 1e-6, they are already outside normal survey-calibration usage.

### Synthetic data generation

For a given `(log_complexity, log_tol)` evaluation point:

1. Derive `n` by holding `K=9` and `cats_per_margin = round(sqrt(10^log_complexity / (9 * n)))`. In practice: fix `cats_per_margin ∈ {4, 8, 16}` based on the complexity tercile, then solve for `n = round(10^log_complexity / (9 * cats_per_margin))`. Round-trip check: `abs(log10(n * 9 * cats_per_margin) - log_complexity) < 0.15`; if violated, skip point and log warning.

2. Set seed deterministically: `set.seed(as.integer(round(log_complexity * 1e4)) * 10000L + as.integer(round(-log_tol * 1e4)) %% 10000L)`. No external hash function required.

3. Generate independent categorical margins (same structure as Stepstone): for each margin, sample population proportions from `rdirichlet(1, rep(1, cats_per_margin))`, then generate survey sample with ~10% relative bias.

4. Compute targets from population proportions (sum to 1 exactly via `proportions()` normalization).

5. Pass `max_weight = Inf`, `min_weight = 0`, `convergence = list(absolute = 10^log_tol)` to `harvest()` — explicit, not via defaults.

6. Two warmup runs (discarded) per algo to absorb JIT and cache effects.

### Timing function `time_cell(log_complexity, log_tol, K = 9)`

```r
time_cell <- function(log_complexity, log_tol, K = 9) {
  # 1. Derive and validate n, cats_per_margin
  # 2. set.seed deterministically (formula above)
  # 3. Generate data
  # 4. 2 warmup runs per algo (discarded)
  # 5. 5 timed runs per algo, return median
  # Returns: log(median_t_iEPPA / median_t_LBFGSB)
  #   positive = iEPPA slower = L-BFGS-B wins
}
```

- Single-threaded (no CPU contention confound)
- `method="ieppa"` and `method="lbfgsb"` both called with `max_weight=Inf`
- Returns `y = log(t_iEPPA / t_LBFGSB)` — contour of interest is `y = log(1.2) ≈ 0.182`

---

## Bayesian Level Set Estimation Loop

### Initial design
8-point Latin hypercube over `[4, 7.7] × [−6, −3]` via `lhs::randomLHS(n=8, k=2)`, scaled to the input bounds. Evaluated sequentially.

### GP model
- Package: `DiceKriging::km()`
- Kernel: Matérn-5/2
- Nugget: estimated (`nugget.estim = TRUE`) with lower bound `nugget = 1e-4` to prevent degenerate fits on small-n datasets
- Refit from scratch each iteration on the full accumulated dataset (cheap at ≤33 points)

### Straddle acquisition
Evaluated on a 50×50 candidate grid over `[4, 7.7] × [−6, −3]`:

```r
threshold <- log(1.2)
pred <- predict(gp_model, newdata = candidates, type = "UK")
a    <- -abs(pred$mean - threshold) + 2 * pred$sd  # κ=2 (Bryan et al. 2005)
next_point <- candidates[which.max(a), ]
```

### Termination
After each acquisition, classify the 50×50 grid using the GP posterior CDF:

```r
p_above    <- pnorm(threshold, mean = pred$mean, sd = pred$sd, lower.tail = FALSE)
classified <- mean(p_above > 0.95 | p_above < 0.05)
```

Stop when `classified ≥ 0.90` **or** 25 acquisitions exhausted (whichever first).

**Poor-fit protocol**: if the 25-acquisition cap is reached with `classified < 0.90`:
1. Print a prominent warning: "GP classification < 90% at termination — inspect uncertainty plot before committing constants."
2. Do **not** update `src/c_api.cpp` constants.
3. File a follow-up issue for the 3D sweep.

### Checkpointing
Every 5 acquisitions: write to a temp file then rename atomically:
```r
tmp <- paste0(checkpoint_path, ".tmp")
saveRDS(state, tmp)
file.rename(tmp, checkpoint_path)
```
`state` contains: design matrix, observed y values, serialized GP object, iteration count, classified fraction. On restart: load checkpoint, refit GP from saved design+y (do not re-run evaluations), resume acquisition loop.

### Total budget
8 initial + ≤25 adaptive = ≤33 evaluations × ~30s each ≈ 15–20 min wall time.

---

## K-Stability Check

After the main 2D GP sweep, run a fixed 4×4 grid of `(log_complexity, log_tol)` points at K ∈ {3, 18} (holding the same seeds as the corresponding K=9 points via `set.seed(... + K * 1e7L)`). Plot the three 1.2× contours overlaid in `algo_selection_k_stability.pdf`.

**Decision rule**:
- If the K=3 and K=18 contours lie within ±0.5 log-units of the K=9 contour at every tol level → constants are general, no scope caveat needed.
- If they diverge beyond 0.5 log-units → add a comment to `select_algorithm()`: `// Threshold calibrated for K≈9, uniform-category problems (Stepstone regime). Re-derive for K<<3 or K>>18.`

---

## Output & Threshold Update

### Plots
1. **`algo_selection_contour.pdf`**: GP posterior mean heatmap with `log(1.2)` contour overlaid; design points marked (circle = initial LHC, cross = adaptive). Axes: `log₁₀(complexity)` × `log₁₀(tol_abs)`.
2. **`algo_selection_uncertainty.pdf`**: posterior σ surface.
3. **`algo_selection_k_stability.pdf`**: three contours (K=3, 9, 18) overlaid.

### Contour → rectangular gate approximation

The GP contour is a curve. The production gate is a rectangle (`&&` of two thresholds). The approximation protocol:

1. Read the contour curve from the plot.
2. Find the complexity value at `tol = 1e-3` (loose) and `tol = 1e-6` (tight). Call them `C_loose` and `C_tight`.
3. If `|log₁₀(C_loose) - log₁₀(C_tight)| < 0.5` → Case A (complexity dominates): update `kComplexityThreshold = C_loose` only. No `kTolThreshold`.
4. If the difference exceeds 0.5 log-units → Case B: set `kComplexityThreshold = C_tight` (tighter, conservative) and `kTolThreshold` at the tol value where the contour crosses `log₁₀(complexity) = 5.5` (the midpoint of the sweep). The AND-gate is deliberately conservative: it routes to L-BFGS-B only when both conditions are clearly satisfied, and defaults to iEPPA otherwise.

### Full proposed `select_algorithm()` body (Case B)

```cpp
static rk_algorithm_t select_algorithm(int n, int K,
                                        const int* cat_counts,
                                        const rk_params_t* p,
                                        int64_t& complexity_out) {
    complexity_out = INT64_C(0);
    for (int k = 0; k < K; k++) complexity_out += (int64_t)n * cat_counts[k];

    if (p->algorithm != RK_ALG_AUTO) return p->algorithm;

    // Box-constrained: iEPPA always (Dykstra exact, no logit singularity).
    if (std::isfinite(p->max_weight) || p->min_weight > 0.0)
        return RK_ALG_IEPPA;

    // Unconstrained: L-BFGS-B only when superlinear convergence pays off.
    // Threshold calibrated 2026-04-20 via Bayesian LSE benchmark:
    //   benchmarks/algo_selection_results.rds
    // Valid for K≈9, uniform-category problems. See design spec.
    // NOTE: default for large/loose unconstrained problems is now iEPPA,
    // replacing the prior default of L-BFGS-B. The benchmark confirms
    // iEPPA is competitive or faster in that regime.
    if (p->tol_abs < kTolThreshold && complexity_out < kComplexityThreshold)
        return RK_ALG_LBFGSB;

    return RK_ALG_IEPPA;
}
```

**Semantic change from current code**: the prior default for unconstrained + small complexity was `RK_ALG_LBFGSB`. The updated code defaults to `RK_ALG_IEPPA` for unconstrained + large complexity or loose tolerance, and to `RK_ALG_LBFGSB` only for the tight-tolerance + small-complexity intersection. This is the correct outcome if the benchmark confirms iEPPA is competitive in the loose-tolerance regime.

### Validation after constant update

Two validation cases:

1. **Box-constrained (iEPPA path unchanged)**: run `benchmarks/stepstone_benchmark.R` (n=200K, tol=1e-3, max_weight=5). Confirm auto-selection picks iEPPA; median harvest time within 5% of pre-update baseline.

2. **Unconstrained L-BFGS-B path**: run `harvest(data_small, targets, max_weight=Inf, method="auto", convergence=list(absolute=1e-8))` on a synthetic n=500, K=3, cats_per_margin=2 problem. Confirm `attr(result, "algorithm") == "lbfgsb"`. This exercises the new `kTolThreshold` branch.

### Regression test (`tests/testthat/test-algo-selection.R`)

Added to deliverables. Contents: call `select_algorithm()` (or a thin R wrapper) at three known grid points and assert correct algorithm:

| Point | Expected |
|---|---|
| Constrained (max_weight=5, any tol) | iEPPA |
| Unconstrained, tol=1e-3, complexity=10M | iEPPA |
| Unconstrained, tol=1e-8, complexity=50K | L-BFGS-B (Case B only) |

The third assertion is only added after the benchmark confirms L-BFGS-B wins there. If Case A applies (no kTolThreshold), only the first two assertions are added.

### Threshold update commit

Single commit containing:
- Updated constants in `src/c_api.cpp` with comment `// Benchmark: benchmarks/algo_selection_results.rds @ 2026-04-20`
- `benchmarks/algo_selection_results.rds`
- All three PDF plots
- `tests/testthat/test-algo-selection.R`
- No other changes

---

## Dependencies

| Package | Use |
|---|---|
| `lhs` | Latin hypercube initial design |
| `DiceKriging` | GP model fitting and prediction |
| `ggplot2` | Contour and uncertainty plots |
| `leafblower` | The package under test (installed build) |

All are CRAN packages. Add to `Suggests:` in `DESCRIPTION` (benchmark-only, not runtime deps).
