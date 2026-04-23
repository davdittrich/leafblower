# Faithful iEPPA Solver Design

**Status:** Draft rev 3 (post design-review-gate iter 2)
**Author:** Dennis Alexis Valin Dittrich
**Date:** 2026-04-23
**Supersedes:** None (new solver)
**Related:** The current `src/ieppa.cpp` is a misnamed IPF+Dykstra hybrid
(commit `b13fda4`). This spec introduces the paper-faithful solver and renames
the existing one to `raking`.
**Rev 3 changes (iter 2 gate):** (a) X̃ overflow — replaced incorrect
global-scalar renormalization with log-space factors + clip-before-exp
(architect critical finding); (b) `"raking"` name clarified in §4.3 as "IPF
+ Dykstra" hybrid with explicit docstring note; (c) §1b tightened: cell-path
degenerates gracefully for n < 10k, not "no benefit"; (d) verbose
infeasibility enumerates all (k,j) pairs (designer); (e) §12 aligned with §10
(NEWS.md required, version bump skipped); (f) NaN check on S_kj and on
group_ids entries (security); (g) merge gate explicit (CTO).

**Rev 2 changes:** Fixed 13 design-review-gate iter-1 findings.

---

## 1. Problem Statement

The current solver labeled `ieppa` implements a hybrid of cyclic IPF +
Euclidean Dykstra box projection + Dykstra hyperplane projection (Deming-Stephan
1940 + Boyle-Dykstra 1986). It shares no ingredients with the iEPPA paper
(Chu-Liang-Toh-Yang 2022, arXiv:2011.14312 — "inexact Entropic Proximal Point
Algorithm") and has no published convergence proof.

Goal: implement a **paper-faithful** iEPPA solver that (a) empirically
outperforms the current hybrid on compression-friendly workloads, (b) carries
a documented convergence guarantee, and (c) becomes the default `ieppa` method.
The current hybrid is renamed `raking` and kept as an explicit alternative.

Back-compat is explicitly dropped. Users of `method="ieppa"` will get the new
algorithm. No version bump, no migration shim. A `NEWS.md` entry documents
the breaking change (CRAN policy compliance).

---

## 1b. User-facing Use Cases

**UC-1 (census microsimulation):** A researcher running `harvest(data, targets,
method="ieppa")` on n=1M, K=5 margins × 4 categories today waits ~60s.
With faithful iEPPA + cell-compression (M_cell ≤ 1024 vs n=1M = 1000×
compression), expected wall-clock drops to ~5-10s. **Visible benefit:**
interactive iteration on population synthesis workflows becomes feasible.

**UC-2 (tight-bound calibration in panel surveys):** A survey methodologist
calibrating a rotating panel with `max_weight=1.3` currently hits cycling or
non-convergence on the hybrid (see leafblower-370 probe evidence: 3/4 hard
scenarios hit `max_iter`). Faithful iEPPA's BCD-block capacity structure is
documented in the paper as convergent on these regimes (Csiszár 1975 cyclic
I-projection on affine + log-convex box). **Visible benefit:** fewer NOCONV
warnings on production calibration of bounded weights.

**UC-3 (large-K studies):** Multi-dimensional calibration (K=10-20 demographic
margins) on n=500k-1M with moderate cat_counts — current hybrid is O(n·K) per
sweep; faithful iEPPA is O(M_cell·K) with `M_cell << n`. **Visible benefit:**
enables more margins per calibration without linear cost growth.

**Graceful degeneracy (documented):** for n < 10k or `M_cell ≈ n`, faithful
iEPPA's cell-path collapses to obs-level work (M_cell = n → no compression
speedup, but asymptotic cost identical to raking). Measured wall-clock is
comparable (± small constant), not a regression. AUTO routing returns IEPPA
in this regime; it is **not worse**, just not measurably faster. Routing
refinement is a follow-up WU gated on benchmark.

---

## 2. Algorithmic Foundation

### 2.1 Reduction of paper iEPPA to calibration

The paper solves capacity-constrained multi-marginal optimal transport (CMOT):

```
min ⟨C, X⟩   s.t. A^(k)(X) = b^(k),  0 ≤ X ≤ U
```

Setting `C = 0` degenerates this to I-projection onto the intersection of the
affine constraint set and the box — the natural formulation of bounded raking.
At `C = 0`, the paper's outer proximal-point loop is mathematically inert
(the regularization term `ε_k D_KL(X, X^k)` scales the objective but does not
move the argmin when the LP cost vanishes). The inner subproblem,
**Algorithm 2 (algBCD)**, directly solves the calibration problem.

### 2.2 Cell compression

Observations with identical `(g_1, ..., g_K)` tuples receive identical weights
under any multiplicative method. Reformulating at cell level:

- `M_cell` = number of unique tuples = `|{(g_1,...,g_K)}|` ≤ `min(n, ∏ cat_counts)`
- Per-cell state is `O(M_cell · K)` instead of `O(n · K)`
- Setup cost `O(n · K)` (one pass to build the cell table)
- Expansion cost `O(n)` (one pass to distribute cell weights to observations)

Worst case `M_cell = n` (all observations unique): algorithm degenerates
to obs-level, identical per-iteration cost to current hybrid. No fallback
branch needed.

### 2.3 algBCD at C = 0

**State:**
- `X_init[c]` — initial cell weight (sum of obs initial weights in cell c)
- `f[k][j]` — Sinkhorn factor per (margin k, category j), initialized to 1.0
- `W[c]` — capacity multiplier per cell, initialized to 1.0
- Current cell weight: `X[c] = X_init[c] · W[c] · ∏_k f[k][g_k(c)]`

**Per outer iteration:**

1. **Margin sweep (K blocks):** for each margin `k`, for each category `j`:
   ```
   S_kj = Σ_{c: g_k(c)==j} X_init[c] · W[c] · ∏_{m≠k} f[m][g_m(c)]
   if S_kj < 1e-15 · W_total:
       if τ[k][j] > 0: is_infeasible = true
       continue
   f[k][j] = (τ[k][j] · W_total) / S_kj
   ```

2. **Uncapped X̃ computation:**
   ```
   X̃[c] = X_init[c] · ∏_k f[k][g_k(c)]
   ```
   Full recompute per outer iteration (O(M_cell · K)). Incremental per-block
   update is not used — equal asymptotic cost, much higher bug risk. Cache
   `prod_all[c]` is rebuilt after each capacity-block application.

3. **Capacity block (KL-projection onto box = clamp):**
   ```
   X[c] = clamp(X̃[c], L_cell[c], U_cell[c])
   W[c] = X[c] / X̃[c]   # multiplier, fed back into next iteration's S_kj
   ```
   `L_cell[c] = min_weight · n_per_cell[c]`
   `U_cell[c] = max_weight · n_per_cell[c]`

4. **Convergence check** (iter 1 and every 10 iters thereafter):
   ```
   errRp = max_{k,j} | (Σ_{c: g_k(c)==j} X[c]) / (Σ X[c]) − τ[k][j] |
   if errRp < tol_abs: break
   ```

### 2.4 Why this avoids the `b13fda4` cycling bug

The predecessor of the current hybrid clamped weights **inside** the margin
update, which violated the Sinkhorn invariant. Paper's capacity block is a
**separate BCD block** whose output `W[c]` feeds back into the next margin
update via the `W[c]` factor in `S_kj`. The margin sweep re-balances non-capped
cells to hit targets, accounting for the cap. Formally: cyclic KL-projection
onto affine sets (marginal constraints) plus a log-convex box, with
non-empty intersection. Csiszár (1975) + Csiszár-Tusnády (1984) give
convergence; paper's Prop 3 (inner subproblem) formalizes it for this BCD
structure.

---

## 3. Scope

### 3.1 In scope

- New `src/ieppa.cpp` + `src/ieppa.hpp` (faithful algBCD)
- New `src/cell_table.hpp/cpp` (cell compression helper)
- Rename existing `src/ieppa.cpp` → `src/raking.cpp` (unchanged implementation)
- Rename `ieppa_solve` → `raking_solve`, `IEPPAResult` → `RakingResult`
- New `RK_ALG_RAKING = 3` enum entry; existing `RK_ALG_IEPPA = 1` retargets to
  faithful solver
- Update `c_api.cpp::select_algorithm` — `AUTO` returns `RK_ALG_IEPPA` always
  (no heuristic exceptions until benchmark justifies them)
- Update R/Python user-facing method names: `"ieppa"` → faithful, `"raking"` → hybrid
- PRD § US-005 rewrite; new § for `raking`
- Test migration + new tests (§6)
- Bayesian Level Set benchmark (§7)

### 3.2 Out of scope

- Benchmark-driven routing refinement (filed as follow-up WU gated on
  benchmark output)
- Cell-compression aware size diagnostic (follow-up)
- SIMD hints on new solver's hot loops (follow-up after correctness proven)
- Back-compat shims, migration tooling, version bump, NEWS.md entry
- Outer proximal-point loop (dropped; mathematically inert at C=0)
- GPU, multi-threading, vectorized Python path, anything not in current
  PRD scope

---

## 4. Architecture

### 4.1 File layout

```
src/
├── ieppa.cpp          NEW: paper-faithful algBCD at C=0
├── ieppa.hpp          NEW: ieppa_solve declaration, IEPPAResult struct
├── raking.cpp         RENAMED from src/ieppa.cpp (current hybrid, unchanged)
├── raking.hpp         RENAMED
├── cell_table.cpp     NEW: cell compression helper
├── cell_table.hpp     NEW
├── c_api.cpp          UPDATED: select_algorithm routes to faithful iEPPA
├── r_bridge.cpp       UPDATED: method string mapping
└── leafblower.h       UPDATED: enum adds RK_ALG_RAKING
```

### 4.2 C API

```c
typedef enum {
    RK_ALG_AUTO   = 0,
    RK_ALG_IEPPA  = 1,   // paper-faithful algBCD (semantics changed)
    RK_ALG_LBFGSB = 2,
    RK_ALG_RAKING = 3    // NEW: renamed current hybrid
} rk_algorithm_t;
```

Select logic (c_api.cpp):
```c
static rk_algorithm_t select_algorithm(...) {
    if (params->algorithm != RK_ALG_AUTO) return params->algorithm;
    return RK_ALG_IEPPA;  // default. No heuristic exceptions.
}
```

### 4.3 R/Python user-facing

- `method = "ieppa"` → `RK_ALG_IEPPA` (paper-faithful algBCD at C=0)
- `method = "lbfgsb"` → `RK_ALG_LBFGSB` (Deville-Särndal dual; unchanged)
- `method = "raking"` → `RK_ALG_RAKING` (IPF + Dykstra box + Dykstra
  hyperplane; classical raking family — autumn compatible. The name covers
  multiplicative marginal projection + additive box/hyperplane corrections;
  docstring in `man/harvest.Rd` states this explicitly)
- `method = "auto"` → `RK_ALG_AUTO` (= IEPPA unconditionally until
  benchmark-driven routing refinement)

Existing autumn synonyms (`"rake"`, `"nr"`, `"nrake"`) continue to warn
+ map to `"lbfgsb"` per existing PRD § US-001 AC.

---

## 5. Cell Table

### 5.1 Structure

```cpp
struct CellTable {
    int M_cell;                                 // unique cells
    std::vector<int> cell_of;                   // size n: obs -> cell index
    std::vector<int> n_per_cell;                // size M_cell: obs count per cell
    std::vector<std::vector<int>> g_per_cell;   // [K][M_cell]: margin-k category per cell
    double W_input;                             // Σ input weights (preserved sum)
};
```

### 5.2 Construction (sort-based, not hash map)

**Sort-based dedup** (replaces hash map per security review — avoids
`std::unordered_map` collision DoS on adversarial group_ids):

Algorithm:
1. Build `keys[n]` array: row-major encoding of `(g_1, ..., g_K)`. Encoding:
   packed int64 when `K ≤ 8` and all `cat_counts[k] ≤ 256`; else fixed-stride
   int32 tuples stored in a `std::vector<std::array<int32_t, K_MAX>>`.
2. Sort observations by `keys[i]` using `std::sort` (O(n log n · K)
   comparator cost).
3. Scan sorted sequence, increment cell index when key changes, populate
   `cell_of[i]` via reverse-mapping back to original obs indices, accumulate
   `n_per_cell[c]`.
4. No hash map used; no DoS surface.

**Complexity:** O(n · K · log n) construction, O(n) storage for sort index.
At n=1M, K=5: ~20M comparisons ~= single-digit seconds one-time. Amortizes
over BCD iterations.

**Safety:**
- `K > K_MAX` (= 64) rejected at `validate_inputs` (§8) with `RK_ERR_BADARG`.
  No unbounded allocation possible.
- NA (`group_ids[k][i] = -1`) treated as a distinct category per margin
  (encoded as `cat_counts[k]` in the key, one past the valid range).
- Deterministic: same input → identical cell_table output.

### 5.3 Expansion (handles non-uniform initial weights within cell)

Correct multiplicative recovery — obs i in cell c gets `w[i] = d[i] · m_c`
where `m_c = ∏_k f[k][g_k(c)] · W[c]` is the cell multiplier. Using the
cell aggregate identity `X[c] = m_c · X_init[c]`:

```cpp
for (int i = 0; i < n; i++) {
    int c = cell_of[i];
    if (X_init[c] > 0) {
        weights[i] = input_weight[i] * X[c] / X_init[c];  // d[i] · m_c
    } else {
        weights[i] = 0.0;  // zero-weight cell: no contribution
    }
}
```

**This formula is correct regardless of initial-weight uniformity within
cell.** Previous `w[i] = X[c] / n_per_cell[c]` was wrong when `d[i]` varied
within cell — fixed per architect review.

Obs i, j in same cell with identical `d[i] = d[j]` receive identical output
weights (the special case). Obs with different `d` within cell receive
proportionally different output weights.

---

## 6. Testing

### 6.1 Migrated test suite

**`tests/testthat/test-ieppa.R`** retargets faithful iEPPA. Existing cases:

| Existing test | Expected outcome | Notes |
|---|---|---|
| 1 margin, 2 cats, no bounds | pass, `max_err < 1e-6` | pure Sinkhorn; behaviorally identical to raking |
| 3 margins, tight bounds (max=2) | pass, `max_err < 1e-6` | cap active; faithful iterates differ from raking but converge |
| `min_weight=0.5` | pass, `min(w) ≥ 0.5 - 1e-10` | both-sided cap |
| 1M rows, 20 margins, max_weight=3 | pass, time < 30s | cell-compression should yield speedup vs raking |

Tolerance / iteration caps reviewed per test; relaxed only when
justified by algBCD semantics.

**`tests/testthat/test-raking.R`** (NEW, regression guard): exact copy of
old `test-ieppa.R` assertions against `method="raking"`. Locks in current
hybrid behavior.

### 6.2 New: algBCD-specific tests

**`tests/testthat/test-ieppa-faithful.R`:**

| Test | Assertion |
|---|---|
| Cell compression correctness | `M_cell` on synthetic (n=1000, K=3, 4 cats, identical obs) = 64 |
| Cell compression extreme | `M_cell` on all-unique synthetic ≈ n (bounded by `∏ cat_counts`) |
| Within-cell weight equality | obs in same cell have identical weights to 1e-12 |
| Cap-inactive scenario | loose bounds: `W[c] == 1` for all c at convergence |
| Cap-active scenario | tight bounds: `Σ 1[W[c] ≠ 1] > 0`; targets still met |
| Infeasibility detection | empty cell + positive target → `RK_ERR_INFEAS`, not cycle |
| Both-sided cap | `min_weight > 0` ∧ `max_weight < Inf`: some cells at L, some at U; targets met |

**`tests/testthat/test-compare.R`** (NEW, cross-algorithm functional
equivalence):

On 20 random feasible datasets with **`set.seed(20260423)`** at test entry
(fixed for reproducibility per CTO C4) — varied n ∈ {1000, 10000, 100000},
K ∈ {3, 5, 10}, cat_counts ∈ {3, 5, 8}, bound tightness `max_weight` ∈
{1.5, 3, 5}. Run `ieppa`, `raking`, `lbfgsb` on each. All three must produce
weights with `max_err < 1e-6`; pairwise max absolute weight difference must
be < 1e-3 (algorithms differ, but on feasible problems must agree to three
sig figs).

---

## 7. Benchmark (Bayesian Level Set Estimation)

Reuses framework from `benchmarks/algo_selection_benchmark.R`:
`lhs` for space-filling initial design, `DiceKriging::km` GP surrogate,
adaptive follow-up sampling around the decision boundary.

### 7.1 Input space (3D)

- `x1 = log10(complexity) = log10(n · Σ cat_counts)` ∈ [4, 7.7] (10K to 50M)
- `x2 = log10(tol_abs)` ∈ [-6, -3]
- `x3 = log10(n / M_cell)` ∈ [0, 3.5] (1× to ~3000× compression)

### 7.2 Response

```
y = log(t_ieppa / t_raking)
```

Negative = faithful wins; positive = raking wins.

### 7.3 Decision threshold

`log(1.2)` (20% speedup margin to trigger switch) — same as existing
`algo_selection_benchmark.R`.

### 7.4 Design protocol

- Initial: Latin Hypercube Sample of 16 points
- Adaptive augmentation: ~8 follow-up points via maximin LHS around
  predicted decision boundary (mirrors `run_adaptive()` from existing script)
- K-stability: re-run at K ∈ {3, 10, 20} to verify GP surrogate insensitive
  to parametrization (mirrors `run_k_stability()`)

### 7.5 Data generation per design point (deterministic, fixed seed)

**x3 redefined** (per architect review A3): `x3 = log10(∏ cat_counts / n)`
— theoretical max compression. Actual `M_cell` depends on cell occupancy
and is reported separately in the output RDS.

Generation procedure:
1. Fix `K = 5` for primary benchmark.
2. Derive `cat_counts[k]` so that `∏ cat_counts = round(n · 10^x3)`; default
   all equal: `cat_counts[k] = round((n · 10^x3)^(1/K))`.
3. `n = round(10^x1 / Σ cat_counts)`.
4. `tol_abs = 10^x2`.
5. `max_weight = 3` fixed (bound-tightness a follow-up WU).
6. **Data sampling:** if `∏ cat_counts ≤ n/2` (high expected occupancy):
   generate full-coverage via Latin Hypercube over the cell grid then pad
   with uniform-random extras (ensures all cells populated). Else: uniform
   random sampling (actual `M_cell ≈ min(n, ∏ cat_counts)`).
7. `set.seed(42)` at script entry; per-point seed derived from
   `(x1, x2, x3)` via deterministic hash (reproducibility across K-stability
   re-runs).
8. Record actual `M_cell` in the output for each design point.

### 7.6 Output artifacts

```r
list(
  design = matrix,          # LHS input points (16 + adaptive)
  y = numeric,              # log-ratio responses
  gp = km object,           # DiceKriging GP fit
  threshold = log(1.2),
  x_ranges = list(x1=..., x2=..., x3=...),
  meta = list(seed, runs_per_point, n_adaptive_iters, ...)
)
```

Saved to `benchmarks/ieppa_vs_raking_results.rds`.

Contour + uncertainty PDFs via `make_plots()` helper (refactored out of
`algo_selection_benchmark.R` for reuse).

### 7.7 Report

Markdown summary of fitted GP: decision boundary as a function of
`(complexity, tol, compression)`. This report is the **input** to the
benchmark-driven routing refinement WU (filed separately, gated on
benchmark output).

---

## 8. Error Handling

| Condition | Code | Detection |
|---|---|---|
| NULL, bad dims, NaN/Inf, targets !sum to 1 | `RK_ERR_BADARG` | Shared input validation |
| `min_weight ≥ max_weight` | `RK_ERR_BADARG` | Input validation |
| **`K > 64`** (new) | `RK_ERR_BADARG` | Input validation. Prevents unbounded cell key allocation |
| **`Σ weights ≤ 1e-15`** (new) | `RK_ERR_BADARG` | Input validation. Degenerate total weight |
| Empty cell (`X_init[c] = 0`) | skipped | Cell omitted from margin sweep; `weights[i] = 0` on expansion per §5.3 |
| Empty cell with positive target | `RK_ERR_INFEAS` | Margin sweep; latched flag, reported after loop |
| `inner_max_iter` exhausted, `errRp ≥ tol_abs` | `RK_ERR_NOCONV` | Loop termination |
| `M_cell · K · sizeof(double) > SIZE_MAX/2` | `RK_ERR_BADARG` | `CellTable` construction |
| **`X̃[c]` overflow** (new handling) | `RK_ERR_NOCONV` with message | See below |

**X̃[c] overflow handling** (rewritten rev 3 — log-space factors, no
renormalization. Architect iter-2 finding: global-scalar renormalization
scaled X̃ by S^K, not 1, altering the primal iterate).

Work with log-space dual factors. Maintain `lf[k][j] = log f[k][j]` rather
than `f[k][j]`. Key quantities:

- Margin sweep update: `lf[k][j] = log(τ[k][j] · W_total) − log(S_kj)` where
  `S_kj = Σ_{c: g_k(c)=j} X_init[c] · exp(Σ_{m≠k} lf[m][g_m(c)]) · W[c]`
- `log X̃[c] = log X_init[c] + Σ_k lf[k][g_k(c)]`
- `X̃[c] = exp(min(log X̃[c], 700.0))` — clip-before-exp guarantees no IEEE
  754 overflow. `exp(700) ≈ 1.01e304`, safely below `DBL_MAX`.
- `X[c] = clamp(X̃[c], L_cell[c], U_cell[c])` — when clipping fires, X̃ is
  at 1e304, upper bound is ≤ `max_weight · n_per_cell[c]` (much smaller for
  realistic inputs), so `X[c] = U_cell[c]`. Correct clamp behavior.
- `W[c] = X[c] / X̃[c]` — safe when `X̃[c] > 0`. If `X̃[c] = 0` (empty
  cell): set `W[c] = 1.0` (multiplier is immaterial; cell aggregate is 0).
  If `X̃[c] = NaN` (shouldn't occur given log-space): exit with
  `RK_ERR_NOCONV` and message `"iEPPA: NaN detected in X̃; numerical
  breakdown. Try looser tol_abs or method=raking."`

**Overflow detection:** if `max_c log X̃[c] > 700` AND the corresponding
cell is NOT capped (i.e., U_cell[c] would not fire within reasonable
magnitude), problem is ill-conditioned beyond double precision. Exit with
`RK_ERR_NOCONV`, message `"iEPPA: log-factor overflow indicates
ill-conditioning; try looser max_weight or tighter tol_abs."` In practice
this is rare because the capacity cap fires first on tight-bound problems.

No renormalization of `lf[k][j]` is performed. The log-space formulation
avoids the multiplicative drift that required renormalization in rev 2.

**NaN safety** (security review rev 2):
- Before division in Sinkhorn update: check `!std::isnan(S_kj) && !std::isinf(S_kj)`.
  On failure: set `is_infeasible = true`, skip update. Treated as degenerate
  cell same as `S_kj < 1e-15·W_total` path.
- `validate_inputs` rejects NaN/Inf in `group_ids[k][i]` (any such value
  encoded as `int32` is undefined; we treat any `group_ids[k][i] < -1` or
  `>= cat_counts[k]` as `RK_ERR_BADARG` already per existing FR-4).

**Infeasibility flag semantics:**
`is_infeasible` is **latched** (set-once, never cleared). Triggering condition
`S_kj < 1e-15·W_total ∧ τ_kj > 0` evaluated every margin sweep; once set,
overrides `RK_OK`/`RK_ERR_NOCONV` with `RK_ERR_INFEAS` on return. All
infeasible (k,j) pairs are recorded in a list; `verbose≥1` enumerates them
all at exit (designer review iter 2).

---

## 8b. Verbose Output Contract

Defined at design time per designer review D2.

**`verbose = 0` (default):** silent. Errors via return code only.

**`verbose = 1` (progress):**
- On entry: `"iEPPA: n=N K=K M_cell=M compression=R×"` (R = n/M_cell)
- Every 10 outer iters + final: `"iEPPA iter T: errRp=E"` (matches raking
  solver's existing log cadence for consistency)
- On exit: `"iEPPA converged in T iters, errRp=E"` OR `"iEPPA: max_iter
  exhausted, errRp=E (status=NOCONV)"`
- On infeasibility: enumerate **all** infeasible (k,j) pairs — e.g.
  `"iEPPA: infeasible cells detected: margin=1 cat=3, margin=4 cat=2, ... (3 total, status=INFEAS)"` — not just the first one, so users can diagnose all constraint violations in a single run.

**`verbose = 2` (debug):** all of the above plus per-iter `f[k]` max/min
summary, per-iter count of cells with `W[c] ≠ 1` (capacity-active), and
renormalization event log.

Emitted via `st.log()` (types.hpp:30), verbose-gated; consistent with existing
iEPPA/raking/L-BFGS-B verbose patterns.

---

## 8c. Merge Gate (explicit per CTO review iter 2)

Faithful iEPPA may merge to main when:

1. All tests in §6 pass (`test-ieppa.R`, `test-raking.R`, `test-ieppa-faithful.R`,
   `test-compare.R`): `FAIL 0 | PASS N`
2. `R CMD check --as-cran` yields 0 ERROR, 0 WARNING (NEWS.md entry present)
3. On `M_cell = n` degenerate input: iteration count and wall-clock
   comparable to `raking` (within ±2×). Documented in merge commit body.
4. Python pytest suite passes (where applicable, after bindings wiring).

Benchmark (§7) is **post-merge analysis** feeding the routing refinement WU,
not a merge gate on this design.

---

## 9. Convergence Guarantees

- Inner algBCD converges to the I-projection of `X_init` onto the intersection
  of affine marginal constraints and the log-convex box (Csiszár 1975,
  Csiszár-Tusnády 1984). Paper's Prop 3 (inner subproblem) provides the
  formal statement in the BCD variable setting.
- Rate: linear under Slater (strictly interior feasible point exists).
  Degrades toward sub-linear as feasibility boundary tightens.
- Outer proximal-point loop from paper is **not** implemented. At `C=0` it is
  mathematically inert. Documented in solver's leading comment and in PRD.

---

## 10. PRD + Documentation Updates (same commit as rename)

- § US-005 rewritten to describe actual faithful algBCD (no misstatement
  about the current hybrid)
- New § US-005b: `method="raking"` for the renamed IPF+Dykstra hybrid
- FR entries updated to match
- **`NEWS.md` entry** (required per CRAN policy; one line + explanation):
  ```
  # leafblower (development)
  * BREAKING: method="ieppa" now runs the paper-faithful algBCD (Chu et al.
    2022, arXiv:2011.14312). The previous implementation was an IPF+Dykstra
    hybrid misnamed "iEPPA"; it is renamed method="raking". Users relying on
    the previous "ieppa" behavior should switch to method="raking".
  ```
- `man/harvest.Rd`: updated `method` parameter description covers both
  `"ieppa"` (faithful) and `"raking"` (renamed hybrid), explicit note that
  `"auto"` currently returns `"ieppa"` unconditionally until benchmark-driven
  refinement.

---

## 11. Non-Goals

As in current PRD plus:

- Outer proximal-point loop
- Variable-ε schedule
- Entropic regularization tuning
- Benchmark-driven routing (separate follow-up WU)
- Migration tooling or back-compat shims
- Bound-tightness sensitivity benchmark (follow-up)

---

## 12. Open Questions Resolved in Brainstorm

| Question | Decision |
|---|---|
| Cost matrix for calibration | `C = 0`; outer PP loop dropped (inert at C=0) |
| Cell-level vs obs-level | Cell-level always; degenerates gracefully when `M_cell = n` |
| Back-compat for `method="ieppa"` | None; breaking change is intentional |
| AUTO routing heuristic | None until benchmark justifies exceptions |
| Version bump | Skipped; user directive |
| NEWS.md entry | **Required** (CRAN policy compliance) — text in §10 |
| Runtime deprecation warning on method="ieppa" | Skipped; user directive. Communication via NEWS.md only |
| Migration shim | Skipped; user directive (breaking change intentional) |

---

## 13. Deliverables

Planner (writing-plans skill) decomposes into WUs. Expected atomic units:

**Prerequisite (small, early):**
- `cell_table.hpp/cpp` implementation + tests (sort-based dedup)
- Input validation update: K ≤ 64 + Σweights > 0 guards

**Core solver:**
- `ieppa_faithful.hpp/cpp` implementation + unit tests (§6)
- Renormalization + overflow detection per §8
- Verbose contract per §8b

**Atomicity constraint (hard, per CTO C3):**
- **Single atomic commit** bundles: rename `src/ieppa.cpp` → `src/raking.cpp`,
  symbol rename, new `src/ieppa.cpp`, enum addition `RK_ALG_RAKING=3`,
  `c_api.cpp::select_algorithm` update, R/Python method-name wiring, PRD
  update, NEWS.md entry. Intermediate states must not break the build; CI
  must not run between sub-steps.

**Testing:**
- `test-raking.R` (regression guard) + `test-ieppa.R` retargeting +
  `test-ieppa-faithful.R` + `test-compare.R` (§6)

**Benchmark:**
- `benchmarks/ieppa_vs_raking_bench.R` (§7)
- **Prerequisite sub-WU:** refactor `make_plots()` from
  `benchmarks/algo_selection_benchmark.R` to handle 3D input space (slice
  plots over compression dimension). Scoped as its own WU per CTO C2.

**Follow-up WUs filed separately (not this design):**
- Benchmark-driven routing refinement (gated on §7 output)
- SIMD hints on faithful solver hot loops (post-merge, after correctness)
- Bound-tightness sensitivity benchmark (vary `max_weight` dimension)
- Compression-ratio diagnostic in verbose output + return attribute
  `attr(result, "compression_ratio")` per Designer review
- Log-space formulation of f/X (alternative to renormalization, if overflow
  triggers frequently in practice)

Implementation plan and WU issue creation follow the writing-plans skill.
