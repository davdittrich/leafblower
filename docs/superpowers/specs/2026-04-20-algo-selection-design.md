# Algorithm Selection Benchmark Design

**Date**: 2026-04-20  
**Status**: Approved  
**Scope**: Unconstrained calibration only (`max_weight=Inf`, `min_weight=0`)

---

## Goal

Find the 1.2× iso-contour in `(log_complexity, log_tol_abs)` space that separates "L-BFGS-B is ≥20% faster than iEPPA" from "iEPPA wins". Use the contour to hand-tune two constants in `src/c_api.cpp`.

**Win criterion**: L-BFGS-B wins at a grid point if `t_iEPPA / t_LBFGSB ≥ 1.2`.

---

## Scope

### In scope
- Unconstrained problems: `max_weight=Inf`, `min_weight=0`
- Input axes: `log₁₀(complexity = n × Σcat_counts)` and `log₁₀(tol_abs)`
- Algorithm: Bayesian level set estimation (Straddle acquisition) with a GP surrogate

### Out of scope
- Box-constrained problems — iEPPA is unconditionally correct there; no benchmark needed
- Automating the constant update — human reads the plot, updates `c_api.cpp`
- Varying margin structure beyond the `complexity` collapse (known limitation, see below)

### Known limitation
`complexity = n × Σcat_counts` collapses two distinct algorithmic factors: `n` (cost per dual-space operation) and `Σcat_counts` (dual dimension for L-BFGS-B memory). The benchmark fixes `K=9` margins with `cats_per_margin` varied to hit target complexity, which mirrors the Stepstone use case. If the GP posterior shows systematic scatter (poor fit), a 3D sweep adding `Σcat_counts / n` as a separate axis is the follow-up.

---

## Deliverables

| File | Purpose |
|---|---|
| `benchmarks/algo_selection_benchmark.R` | LSE benchmark script |
| `benchmarks/algo_selection_results.rds` | Raw timing data + GP object (reproducibility) |
| `benchmarks/algo_selection_contour.pdf` | GP posterior mean + 1.2× contour |
| `benchmarks/algo_selection_uncertainty.pdf` | GP posterior σ (shows fit quality) |
| `src/c_api.cpp` | Updated `kComplexityThreshold` + new `kTolThreshold` (manual commit after review) |

---

## Benchmark Harness

### Input space
- `x₁ = log₁₀(complexity)` ∈ [4, 7.7] — covers complexity 10K to 50M
- `x₂ = log₁₀(tol_abs)` ∈ [−9, −3] — covers tol 1e-9 to 1e-3

### Synthetic data generation
For a given `(log_complexity, log_tol)` evaluation point:
1. Derive `n` and `cats_per_margin` from `complexity = n × K × cats_per_margin` with `K=9` fixed
2. Generate a synthetic population and biased survey sample with `set.seed(hash(log_complexity, log_tol))` — deterministic per grid point
3. Compute targets from population proportions (always sum to 1)
4. Margins are independent categoricals — same structure as Stepstone

### Timing function `time_cell(log_complexity, log_tol)`
```r
time_cell <- function(log_complexity, log_tol) {
  # 1. Derive problem dimensions
  # 2. Generate data
  # 3. Warmup run (discarded)
  # 4. 5 timed runs per algo, take median
  # Returns: log(t_iEPPA / t_LBFGSB)
}
```

- Single-threaded (no CPU contention confound)
- `method="ieppa"` and `method="lbfgsb"` called back-to-back
- Returns `y = log(t_iEPPA / t_LBFGSB)` — positive means iEPPA is slower (L-BFGS-B wins)
- Contour of interest: `y = log(1.2) ≈ 0.182`

---

## Bayesian Level Set Estimation Loop

### Initial design
8-point Latin hypercube over the input space via `lhs::randomLHS()`. Evaluated sequentially.

### GP model
- Package: `DiceKriging::km()`
- Kernel: Matérn-5/2
- Nugget: estimated (timing noise is real and heteroskedastic)
- Refit from scratch each iteration on the full accumulated dataset (cheap at ≤33 points)

### Straddle acquisition
Evaluated on a 50×50 candidate grid over `[4, 7.7] × [−9, −3]`:

```r
threshold <- log(1.2)
pred <- predict(gp_model, newdata = candidates, type = "UK")
a    <- -abs(pred$mean - threshold) + 2 * pred$sd  # κ=2
next_point <- candidates[which.max(a), ]
```

Pulls evaluations toward the contour (low `|μ − threshold|`) while maintaining exploration (high `σ`). `κ=2` is standard for Straddle.

### Termination
After each acquisition, classify each of the 50×50 candidate grid points using the GP posterior CDF:

```r
p_above <- pnorm(threshold, mean = pred$mean, sd = pred$sd, lower.tail = FALSE)
classified <- mean(p_above > 0.95 | p_above < 0.05)
```

Stop when `classified ≥ 0.90` **or** 25 acquisitions exhausted (whichever first). The cap prevents runaway on poor GP fits.

### Checkpointing
Every 5 acquisitions: write `benchmarks/algo_selection_results.rds` containing the full design matrix and serialized GP object. Script is restartable from checkpoint.

### Total budget
8 initial + ≤25 adaptive = ≤33 evaluations × ~30s each ≈ 15–20 min wall time.

---

## Output & Threshold Update

### Plots
1. **`algo_selection_contour.pdf`**: GP posterior mean heatmap with `log(1.2)` contour overlaid; design points marked (symbol = initial vs adaptive). Axes: `log₁₀(complexity)` × `log₁₀(tol_abs)`.
2. **`algo_selection_uncertainty.pdf`**: posterior σ surface — shows where the GP is well-fit vs uncertain. Guides judgment on whether the contour estimate is trustworthy.

### Reading the constants
The contour is a curve in 2D. Approximate it as a step function:

**Case A — contour approximately vertical** (complexity dominates, tol_abs irrelevant):
- Update `kComplexityThreshold` only
- No `kTolThreshold` needed

**Case B — contour has clear tol dependence** (expected from theory):
```cpp
static constexpr int64_t kComplexityThreshold = /* read from plot */;
static constexpr double   kTolThreshold        = /* read from plot */;

// In select_algorithm():
if (p->tol_abs < kTolThreshold && complexity < kComplexityThreshold)
    return RK_ALG_LBFGSB;
return RK_ALG_IEPPA;
```

The existing box-constraint guard remains unchanged and fires first:
```cpp
if (std::isfinite(p->max_weight) || p->min_weight > 0.0)
    return RK_ALG_IEPPA;
```

### Validation after constant update
Run `benchmarks/stepstone_benchmark.R` (n=200K, tol=1e-3, constrained — the default use case). Confirm:
1. Auto-selection picks iEPPA (box-constrained → unchanged)
2. Median harvest time within 5% of pre-update baseline

### Threshold update commit
Single commit containing:
- Updated constants in `src/c_api.cpp` with comment citing benchmark date and plot
- `benchmarks/algo_selection_results.rds`
- Both PDF plots
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
