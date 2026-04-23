# Faithful iEPPA Solver Design

**Status:** Draft (post-brainstorm, pre-plan-review-gate)
**Author:** Dennis Alexis Valin Dittrich
**Date:** 2026-04-23
**Supersedes:** None (new solver)
**Related:** The current `src/ieppa.cpp` is a misnamed IPF+Dykstra hybrid
(commit `b13fda4`). This spec introduces the paper-faithful solver and renames
the existing one to `raking`.

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
algorithm. No version bump, no migration note — change is intentional.

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
   Efficient via cached product `prod_all[c]`, incrementally updated during
   margin sweep.

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

- `method = "ieppa"` → `RK_ALG_IEPPA` (faithful)
- `method = "lbfgsb"` → `RK_ALG_LBFGSB` (unchanged)
- `method = "raking"` → `RK_ALG_RAKING` (hybrid)
- `method = "auto"` → `RK_ALG_AUTO` (= IEPPA)

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

### 5.2 Construction

Single O(n · K) pass over observations. Key format:
- Packed 64-bit integer when `K ≤ 8` and all `cat_counts[k] ≤ 256`
- `std::string` of concatenated int32 otherwise
- NA (`group_ids[k][i] = -1`) treated as a distinct category per margin

Hash map: `std::unordered_map<Key, int>`. Insertion-order determines
cell indexing; deterministic given input.

### 5.3 Expansion

```cpp
for (int i = 0; i < n; i++) {
    weights[i] = X[cell_of[i]] / n_per_cell[cell_of[i]];
}
```

Output obs weights equal within any cell — a provable property of any
multiplicative calibration method, not a new restriction.

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

On 20 random feasible datasets (varied n, K, cat_counts, bound tightness),
run `ieppa`, `raking`, `lbfgsb`. All three must produce weights with
`max_err < 1e-6`; pairwise max absolute weight difference must be < 1e-3
(not zero — algorithms differ, but on feasible problems must agree to
three sig figs).

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

### 7.5 Data generation per design point

Deterministic `(x1, x2, x3) → (n, K, cat_counts, tol_abs, target_M_cell)`:
- `n = round(10^x1 / Σ cat_counts)` with `K=5`, cat_counts chosen so that
  `∏ cat_counts / n = 10^-x3` (i.e., matches target compression)
- `tol_abs = 10^x2`
- `max_weight = 3` baseline (bound-tightness sensitivity is a follow-up WU,
  not part of primary benchmark)

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
| Empty cell with positive target | `RK_ERR_INFEAS` | Margin sweep; reported after loop |
| `inner_max_iter` exhausted, `errRp ≥ tol_abs` | `RK_ERR_NOCONV` | Loop termination |
| `M_cell · K · sizeof(double) > SIZE_MAX/2` | `RK_ERR_BADARG` | `CellTable` construction |
| `X̃[c]` overflow to `+Inf` | non-fatal | Clamp to 1e300; continue |

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

## 10. PRD Updates (same commit as rename)

- § US-005 rewritten to describe actual faithful algBCD (no misstatement
  about the current hybrid)
- New § US-005b: `method="raking"` for the renamed IPF+Dykstra hybrid
- FR entries updated to match

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
| Version bump / NEWS | Skipped; change is intentional |

---

## 13. Deliverables

Planner (writing-plans skill) will decompose into WUs. Expected atomic units
include:
- Cell table implementation + tests
- Faithful solver implementation + unit tests
- Rename refactor (file + symbol + PRD) as one atomic commit
- Enum addition + routing update
- R/Python method-name wiring
- Test migration + new tests (§6)
- Bayesian benchmark harness (§7)
- Follow-up WU filings: routing refinement, SIMD hints, bound-tightness
  sensitivity

Implementation plan and WU issue creation follow the writing-plans skill.
