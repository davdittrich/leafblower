# ylsy: CP + IPM Research Spike — Design (rev 2)

**Date:** 2026-05-02
**Beads:** `leafblower-ylsy` (kk1204 convergence research)
**Status:** Design — pre-plan, gate iteration 2
**Predecessor:** Epic-H WH-g (AUTO routing safety; routing-only fix; kk1204 still gap=5.01e-2 under ieppa+sraa)
**Revision history:** rev 1 → rev 2 addresses 13 blockers from design-review-gate iteration 1 (R5 isolation broken via package-source globbing; CP/IPM objective divergence; shim signature mismatch; status codes; numerical safety; sanity-test policy; atomic WU decomposition).

## 1. Problem & Goals

### Problem

The kk1204 fixture (n=1M, K=20 dense margins, 5 categories per margin, target=(0.6, 0.2, 0.1, 0.07, 0.03), max_weight=3, target_skew=20) defeats every solver in `harvest()`:

| Solver (master, post-Epic-H) | kk1204 max_err | status | wall |
|---|---|---|---|
| `newton_kl` | 6.24e-2 | 1 (NOCONV) | — |
| `ieppa` + `accelerate=TRUE` (auto via WH-g) | 5.01e-2 | 0 (converged to fixed point) | 11.4s |
| `lbfgsb` (auto-fallback) | 0.31 | 4 | 372s |

The PARTIAL gate target (max_err < 1e-3) is unreached. Underlying cause: kk1204 is PURE slow-rate Sinkhorn (errRp ≈ iter^{-0.5}, monotone, no oscillation; per `bd memory` `fixed-point-accelerators-apva-anderson-joint-anderson-on`). Smooth fixed-point accelerators (APVA, Joint Anderson, Halpern) failed because piecewise-linear bound projections defeat contraction theory.

### Analyst-facing impact

Severe-skew K≥5 calibrations currently return weights with 5% marginal error — unusable for downstream weighting and inadequate for survey publication. PARTIAL gate (max_err < 1e-3) brings kk1204 into the same accuracy class as stepstone (1.13e-4); PASS gate restores publication parity. Without the spike, AUTO routing routes kk1204-class problems to ieppa+sraa as a degraded-quality safe fallback.

### Goal

Add (do not replace) one or two new calibration methods that break the kk1204 ceiling. Two candidate algorithms with orthogonal paradigms:

- **Chambolle-Pock primal-dual (PDHG)** — first-order saddle-point method designed for non-smooth (LP/QP-class) problems. Bypasses the contraction-theory failure mode.
- **Interior-Point Newton (IPM)** — second-order central-path method. Log-barrier smooths the bounds; Newton handles strong convexity of KL.

### Decision Rule

The spike runs CP first; gates on CP verdict; then runs IPM iff CP < PASS.

| Outcome on kk1204 | Stepstone | Verdict | Action |
|---|---|---|---|
| max_err < 1e-3 AND walltime ≤ 30s AND rate exponent β ≤ −0.8 (last-iterate) or ≤ −1.0 (ergodic) AND R² ≥ 0.9 | within 1.5× of ieppa+sraa (≤ 1.7e-4) | **PASS** | File Epic-I (productionize) |
| max_err < 1e-3 AND walltime ≤ 30s AND rate exponent / R² as above | regression beyond 1.5× ieppa+sraa | **PASS-kk1204-specialist** | Productionize ONLY for severe-skew K≥5 routing; AUTO map carve-out documented |
| max_err in [1e-3, 1e-2) AND walltime ≤ 60s | any | **PARTIAL** | One follow-up tuning spike permitted (chain depth cap = 1) |
| max_err ≥ 1e-2 OR walltime > 60s | any | **FAIL** | Investigation report; ylsy stays open or closes BLOCKED |

**Stepstone reference baseline** (Epic-Dβ NEWS): ieppa+sraa achieves max_err ≈ 1.13e-4 on stepstone K=9 fulldata. Parity tolerance × 1.5 → cutoff ≈ 1.7e-4.

**Rate exponent**: slope $\beta$ of $\log(\mathrm{err}) = \alpha + \beta \log(\mathrm{iter})$. Sinkhorn slow-rate gives $\beta \approx -0.5$. PDHG ergodic gives $\beta \approx -1$; PDHG last-iterate is typically $\beta \in (-1, -0.5)$ — relaxed last-iterate threshold $\beta \le -0.8$ accounts for this. Trace must record both x̄ (ergodic) and x (last-iterate); compare beta from both. Fit on iters $\in (50, 0.9 \times n_{\text{iter}})$ requiring $n_{\text{fit\_points}} \ge 30$ and $R^2 \ge 0.9$ for the fit to qualify as PASS-supporting evidence.

**R3 falsification criterion** (joint 5e-2 fixed point): max_err agreement to within 1% across CP, IPM, ieppa+sraa, AND newton_kl on kk1204 → R3 confirmed (basin floor fundamental). Otherwise R3 is ruled out and FAIL is solver-specific.

## 2. Architecture & Isolation Contract

### File layout (package-root `research/`, NOT `src/research/`)

```
research/                              (package root, EXCLUDED via .Rbuildignore)
├── cp_calib.hpp / cp_calib.cpp
├── ipm_calib.hpp / ipm_calib.cpp
├── research_bridge.cpp                (.Call shim: cp_solve_R, ipm_solve_R)
└── Makefile                           (standalone: R CMD SHLIB)

benchmarks/research/
├── ylsy_cp_bench.R
├── ylsy_ipm_bench.R                   (only if CP < PASS)
├── ylsy_compare.R
├── sanity_t1_recovery.R               (CTO blocker: w=1 ground-truth test)
└── results/                           (CSVs)

tools/
└── check_research_isolation.R         (CI symbol-level audit)

.Rbuildignore                          (append: ^research$, ^tools/check_research_isolation\.R$)

docs/investigations/
└── 2026-05-02-ylsy-cp-ipm-spike-result.md
```

### Why package-root `research/` (not `src/research/`)

R's `R CMD INSTALL` globs `src/*.cpp` automatically into the package OBJECTS — `PKG_SOURCES` in `src/Makevars.in` is a custom variable that R does not honor. Verified by Architect reviewer: `src/raking.cpp` is absent from `PKG_SOURCES` yet still compiled and linked into `leafblower.so`. Therefore any `src/research/*.cpp` file would also be globbed in.

`.Rbuildignore` excludes paths from `R CMD build` (source tarball assembly), but `R CMD INSTALL` on the working tree still picks up everything under `src/`. The only safe relocation is OUTSIDE `src/` — package-root `research/` is not on R's compilation path at all.

### Isolation enforcement (mechanical CI gate)

- **Verbatim `.Rbuildignore` lines** (added by WU-1):
  ```
  ^research$
  ^tools/check_research_isolation\.R$
  ```
- **CI script** `tools/check_research_isolation.R`:
  ```r
  # Pre-commit + pre-push gate. Exits 1 if package .so contains research symbols.
  syms <- system2("nm", c("-D", "src/leafblower.so"), stdout = TRUE)
  forbidden <- c("cp_solve_R", "ipm_solve_R", "cp_calibrate", "ipm_calibrate")
  hits <- sapply(forbidden, function(s) any(grepl(s, syms, fixed = TRUE)))
  if (any(hits)) {
    stop("Research symbols leaked into package .so: ",
         paste(forbidden[hits], collapse = ", "))
  }
  cat("research isolation OK\n")
  ```
- **Pre-commit hook** wires `Rscript tools/check_research_isolation.R` after `R CMD INSTALL --preclean .`. WU-1 includes a green-bar test of this gate.
- **Smoke test in WU-1**: install package, dump `.so` symbols, assert no `cp_*` or `ipm_*` exports.

### R/C boundary contract (applies to `research_bridge.cpp` and any C++ TU exposing SEXP)

All `.Call` shim translation units MUST:

1. `#define STRICT_R_HEADERS 1` and `#define R_NO_REMAP 1` at top of file (before any R header).
2. Include `<R.h>`, `<Rinternals.h>`, `<Rdefines.h>` only.
3. Heap allocations: `std::vector<double>` or `Eigen::VectorXd` only — raw `new` / `delete` forbidden.
4. SEXP allocations: every `Rf_alloc*` / `PROTECT` balanced with `UNPROTECT(n)` before return; reviewer audits the count manually before commit.
5. Input dimension validation: check `Rf_length(df) == n_expected` and `Rf_isReal(tgt)` etc. before any pointer dereference. Use `Rf_error("...")` on mismatch.
6. No `Rcpp::` types (avoids LinkingTo dependency drift). Use plain `Rinternals.h` SEXP only.

### `.Call` shim signature (locked contract)

**Input:**
- `A_csr` — sparse design matrix in CSR triplet form `list(p, j, n_row, n_col)` where `p[i:i+1]` gives row offsets and `j[p[i]:p[i+1]]` gives column indices for row `i`. Each row has exactly `K` ones (block-incidence). Pre-built by R caller (`build_design_matrix(df)` in `benchmarks/research/utils.R`), decoupling research code from `harvest()`'s preprocessing.
- `b` — target marginal mass vector, `length = n_col = ΣJ_k`, type `numeric`.
- `d` — base weights (length n, `numeric`); use `rep(1, n)` for unweighted.
- `lo`, `hi` — box bounds, length-n `numeric` (or scalar broadcast).
- `max_iterations` — integer (matches `harvest()` naming).
- `capture_trace` — logical scalar.
- `seed` — integer for any internal RNG (none expected; reserved).

**Output**: named R list:
- `weights` — length-n `numeric`
- `status_code` — integer (0 = converged, 1 = max_iterations exhausted, 2 = NaN/Inf detected, 3 = infeasible bounds, 4 = power-iter divergence)
- `status_msg` — character scalar
- `iterations` — integer (total iterations performed)
- `wall_time_ms` — double
- `trace` — matrix with columns documented per algorithm (Sec 3); zero-row matrix iff `capture_trace=FALSE`. Per-iter scalar aggregates ONLY (no per-element residual vectors); ≤ 10 doubles per iter; max 5000 rows.

## 3. Algorithm Specifications

### Common objective (CP and IPM share)

$$\min_w \; \sum_i d_i \, \mathrm{KL}(w_i / d_i \| 1) \quad \text{s.t.} \quad A^\top w = b, \; \ell \le w \le u$$

Expansion: $d_i \, \mathrm{KL}(w_i/d_i \| 1) = w_i \log(w_i / d_i) - w_i + d_i$.

For the unweighted base case ($d_i = 1$), this collapses to $w_i \log w_i - w_i + 1$. CP and IPM both consume `d` as input; if `d = rep(1, n)` the design-weighted KL reduces to the textbook KL — both forms covered by one implementation.

### Initialization (both algorithms)

$$w^0_i = \mathrm{clamp}(d_i, \, \ell_i + \delta, \, u_i - \delta), \quad \delta = 10^{-8}$$

If $d_i$ falls inside $[\ell_i + \delta, u_i - \delta]$, $w^0_i = d_i$. If $d_i < \ell_i$ (over-bounded), $w^0_i = \ell_i + \delta$. Strict interior start guarantees $\log w^0_i$ and $1/(w^0_i - \ell_i)$ both finite.

Pre-flight check (WU-1, in shim): if $\ell_i \ge u_i - 2\delta$ for any $i$, return `status_code = 3` (infeasible bounds).

### 3.1 Chambolle-Pock (PDHG)

**Reference:** Chambolle & Pock 2011 (J. Math. Imaging Vision 40:120–145), Algorithm 1.

**Saddle-point form:**
$$\min_{w \in [\ell, u]} \max_y \; f(w) + \langle A^\top w - b, \, y \rangle$$
with $f(w) = \sum_i d_i \mathrm{KL}(w_i/d_i \| 1) + I_{[\ell, u]}(w)$. The dual variable $y$ is unconstrained (multiplier on equality constraint $A^\top w = b$); $g^*(y) = b^\top y$.

**Iteration:**

```
||A|| computed once via 50-iter power method on A^T A applied to a random unit vector.
Stopping criterion for power iter: |λ_{k+1} - λ_k| / λ_k < 1e-6 OR k=50.
If not converged at 50 iters AND relative delta > 1e-3: return status_code=4 (R9).
Apply 1.05 safety factor: ||A||_used = 1.05 * ||A||_estimated.
σ = τ = 1 / ||A||_used    (so σ τ ||A||² ≈ 1/1.05² < 1, strictly inside theory).

Init: w⁰ as above, w̄⁰ = w⁰, y⁰ = 0.
For k = 0 .. max_iterations:
  y^{k+1} = y^k + σ (A^T w̄^k − b)
  w^{k+1} = prox_{τ f}(w^k − τ A y^{k+1})
  w̄^{k+1} = 2 w^{k+1} − w^k
  if k % trace_stride == 0 AND capture_trace:
    record (k, time_ms, max_err(w^{k+1}), max_err(ergodic_avg(w)), primal_resid, dual_resid)
  if ‖A^T w^{k+1} - b‖_∞ / max(Z, 1) < 1e-7: status_code = 0; break
  if any NaN/Inf: status_code = 2; break
```

`Z = sum(d)` (precondition: $Z > 0$, asserted by `Rf_error` in shim).

**`prox_{τ f}` per component (direct Newton, NO Lambert-W):**

Solve scalar:
$$\arg\min_{w \in [\ell_i, u_i]} \; \tau d_i \, \mathrm{KL}(w/d_i \| 1) + 0.5 (w - z_i)^2$$

First-order condition: $\tau d_i \log(w/d_i) + (w - z_i) = 0$.

Newton iteration (≤5 iters typical):
- Init: $w_0 = \mathrm{clamp}(d_i, \ell_i + \delta, u_i - \delta)$ (interior).
- Update: $w_{n+1} = w_n - \frac{\tau d_i \log(w_n/d_i) + (w_n - z_i)}{\tau d_i / w_n + 1}$.
- Overflow guard: if $|z_i / (\tau d_i) - 1| > 700$, switch to asymptote $w \approx z_i$ then clamp.
- Stop when $|f'(w_n)| < 10^{-12}$.
- Final clamp to $[\ell_i + \delta, u_i - \delta]$.

**Per-iter cost:** 2 sparse matvecs $A w̄$, $A^\top y$ ($O(nK)$ each on block-incidence A) + n scalar prox ops ($O(n)$ each). Total $O(nK)$ per outer iteration.

**One-time setup cost:** 50 sparse matvecs for power iter ≈ 1s on n=1M, K=20.

**Diagnostic outputs (per traced iter):**

| col | meaning |
|---|---|
| `iter` | iteration index k |
| `time_ms` | cumulative wall time (ms) |
| `max_err_last` | $\max_{k,j} \lvert \hat T_{kj}(w^{k+1}) - T_{kj} \rvert$ |
| `max_err_ergodic` | $\max_{k,j} \lvert \hat T_{kj}(\bar w^k) - T_{kj} \rvert$ where $\bar w^k = \frac{1}{k+1} \sum_{i=0}^k w^i$ |
| `primal_resid` | $\| A^\top w^{k+1} - b \|_\infty$ |
| `dual_resid` | $\| w^{k+1} - w^k \|_\infty / \tau$ (proxy for primal stationarity) |

`trace_stride = max(1, max_iterations / 1000)` to cap trace size at ~1000 rows.

### 3.2 Interior-Point Newton (IPM)

**Same objective + constraints** as CP (Sec 3 common form).

**Log-barrier:** $\Phi(w) = -\mu \sum_i \log(w_i - \ell_i) - \mu \sum_i \log(u_i - w_i)$.

**KKT residual** (with margin multiplier $\lambda \in \mathbb{R}^{\sum J_k}$):

$$r_w = \log(w / d) + \mu \left(\frac{1}{u - w} - \frac{1}{w - \ell}\right) + A \lambda$$
$$r_\lambda = A^\top w - b$$

(The $\log(w/d)$ term comes from $\partial_w (w \log(w/d) - w) = \log(w/d)$; consistent with the design-weighted KL objective.)

**Newton step (Schur on diagonal Hessian):**

```
H_w = diag( 1/w + μ/(w−ℓ)² + μ/(u−w)² )         (diagonal, n entries)
S   = A^T H_w⁻¹ A                                (m × m, m = ΣJ_k = K · J_k for kk1204 = 20·5 = 100)
S   = TSVD(S, ratio=1e-8)                        (handle rank deficiency, mirror Epic-Dβ WL-1)
Δλ  = S⁻¹ (A^T H_w⁻¹ r_w − r_λ)
Δw  = −H_w⁻¹ (r_w + A Δλ)
α_max = max α ∈ (0, 1] s.t. ℓ + δ ≤ w + α Δw ≤ u − δ        (fraction-to-boundary, δ=1e-8)
α   = 0.99 · α_max                                          (single rule, no second min())
w   ← w + α Δw,  λ ← λ + α Δλ
```

**Outer loop:** $\mu_{k+1} = 0.5 \mu_k$, $\mu_0 = 1$, until $\mu < 10^{-9}$.
**Inner loop:** 1–3 Newton steps per $\mu$ until $\|r_w\|_\infty + \|r_\lambda\|_\infty < 10 \mu$.
**Total iters cap:** outer ≤ 50, inner ≤ 5 (hard cap; status_code=1 if exceeded).

**Stop:** $\|r_w\|_\infty + \|r_\lambda\|_\infty < 10^{-7}$ AND $\mu < 10^{-9}$ → status_code=0.
Any NaN/Inf in $w$, $\lambda$, $r_w$, $r_\lambda$ → status_code=2.

**Per-iter cost:** 1 dense $m \times m$ TSVD via LAPACK `dsyevd` ($O(m^3) = 10^6$ FLOPs) + 2 sparse matvecs ($O(nK) = 2 \times 10^7$ FLOPs each). Dominant cost: matvecs. Per-Newton: ~50 ms. Total 30–90 Newton solves × 50 ms = 1.5–4.5 s. Comfortable under 30 s budget.

**Dimensional analysis of $\mu_0 = 1$:** objective magnitude on kk1204 is $\sum_i d_i \mathrm{KL}(w_i/d_i \| 1) \approx n \cdot O(1) = 10^6$. Initial barrier weight $\mu = 1$ contributes $\approx n \cdot \mu = 10^6$ — same order as objective. Adequate centering. Decay 0.5 takes barrier to $10^{-9}$ in 30 outer iterations.

**Diagnostic outputs (per outer iter):**

| col | meaning |
|---|---|
| `iter_outer`, `iter_inner` | (k_outer, k_inner) |
| `time_ms` | cumulative wall time (ms) |
| `mu` | barrier parameter |
| `max_err` | margin error |
| `kkt_resid` | $\|r_w\|_\infty + \|r_\lambda\|_\infty$ |
| `n_projected_dims` | TSVD truncation count (logged each outer iter even when zero — silent rank-deficiency promotion is observable) |
| `alpha_step` | fraction-to-boundary step length |

### 3.3 Faithful-Textbook Constraint

No hyperparameter tuning, no momentum / restart heuristics, no warm-start, no accelerated variants in the spike. Default config only:

- CP: $\sigma = \tau = 1/(\|A\| \cdot 1.05)$, no acceleration.
- IPM: $\mu_0 = 1$, decay $0.5$, vanilla central path (no Mehrotra predictor-corrector).
- All state in `double` (R10: float-precision avoidance — extrapolation $\bar w = 2 w^{k+1} - w^k$ is catastrophic-cancellation-prone in single).

Tuning lands in productionization phase if spike PARTIAL.

## 4. Bench Harness

### Fixtures

| Fixture | Source | Role |
|---|---|---|
| `t1_small` | inline (n=1000, K=3, ncat=3, mild target) | Sanity — any solver hits machine epsilon |
| `stepstone_K9` | `benchmarks/stepstone_bench_data.parquet` + targets JSON | Parity guard against ieppa+sraa baseline |
| `kk1204_K20` | `benchmarks/kk1204_make.R::make_kk1204(n=1e6, K=20, nj=5)` (existing fixture generator at `benchmarks/tsvd_kk1204_K20.R`) | Primary gate |

Fixture generator pinned by name + commit hash in bench script header (R6 reproducibility).

### Sanity-test policy (CTO blocker)

`benchmarks/research/sanity_t1_recovery.R` — runs CP and IPM on a fixture where the analytical solution is $w = d$ (i.e., observed marginals already match targets exactly). Asserts `max_err < 1e-8` and `status_code = 0`. WU-3 / WU-6 must run sanity first; bench gate halts if sanity fails.

### CSV schema (locked, used by `ylsy_compare.R`)

`results/<solver>_<fixture>_trace.csv`:
- algorithm-specific columns per Sec 3 diagnostic tables.

`results/<solver>_<fixture>_summary.csv` (one row per (solver × fixture × repetition)):

| column | meaning |
|---|---|
| `solver` | "cp" / "ipm" / "ieppa+sraa" / "newton_kl" / "lbfgsb" |
| `fixture` | "t1_small" / "stepstone_K9" / "kk1204_K20" |
| `rep` | repetition index (1..3 for kk1204, 1 otherwise) |
| `status_code` | 0 / 1 / 2 / 3 / 4 |
| `status_msg` | character |
| `converged` | logical (status_code == 0) |
| `max_err` | margin max error (last iterate) |
| `max_err_ergodic` | CP only; NA otherwise |
| `wall_s` | elapsed seconds |
| `n_iter` | iteration count |
| `final_primal_resid` | CP only |
| `final_dual_resid_or_kkt` | unified residual field |
| `n_projected_dims_max` | IPM only; NA otherwise |
| `rate_exponent` | $\beta$ from log-log fit |
| `rate_R2` | fit $R^2$ |
| `n_fit_points` | iters used in fit (NA if < 30) |
| `stepstone_parity_ratio` | only for stepstone_K9 row: max_err / 1.13e-4 |
| `seed` | R seed |
| `omp_threads` | OMP_NUM_THREADS |
| `blas_threads` | OPENBLAS / MKL / RhpcBLASctl::blas_get_num_procs() |
| `git_sha` | `git rev-parse --short HEAD` at run time |

`ylsy_compare.R` is a trivial `rbind` over per-solver CSVs (R8: shared schema with NA for non-applicable columns).

### Repetitions and walltime

kk1204 measured 3× with median walltime (R7: thermal noise on 30s gate). `Sys.setenv(OMP_NUM_THREADS = "1")` and `RhpcBLASctl::blas_set_num_threads(1L)` set in bench preamble (R10: BLAS thread contention). Hardware metadata logged: `Sys.info()`, CPU model from `/proc/cpuinfo`, frequency mode if scriptable.

### Pre-flight memory check

`benchmarks/research/utils.R::check_memory(n, K, max_iter)`:
```r
needed <- n * K * 8 + n * 8 * 5 + max_iter * 10 * 8   # sparse A + state + trace
avail  <- as.numeric(system2("free", "-b", stdout=TRUE)[2] |> strsplit("\\s+")[[1]][7])
if (needed > 0.5 * avail) stop("insufficient RAM: need ", needed/1e9, " GB; have ", avail/1e9)
```

R6 codified as gate.

## 5. Out of Scope

| Item | Reason |
|---|---|
| `harvest()` / `r_bridge.cpp` modifications | Productionization-only |
| `leafblower.h` enum (`RK_ALG_CP`, `RK_ALG_IPM`) | Same |
| AUTO routing changes | Same |
| Test additions in `tests/testthat/` | Spike validates algorithm; production tests later. Sanity tests live in `benchmarks/research/` instead. |
| `NEWS.md` / spec changes (other than this design + investigation report) | Investigation report is the artefact |
| Python bindings | Productionization-only |
| Multi-threading / SIMD optimisation | Faithful textbook only; tuning later |
| Hyperparameter sweeps | Default config only |
| Warm-start from ieppa / newton_kl | Cold start only |
| Accelerated CP variant (Algorithm 2 of Chambolle-Pock) | Standard PDHG only |
| Predictor-corrector IPM (Mehrotra) | Vanilla central path only |
| Cell compression in CP / IPM | kk1204 has M_cell/n=1.0; compression buys nothing |
| Convergence improvements to existing methods | New methods only, in `research/` |
| Comparison vs external solvers (CVXPY, IPOPT, MOSEK) | Out of spike; could be follow-up if PASS triggers productionization |

## 6. Risks & Mitigations

| # | Risk | L | I | Mitigation |
|---|---|---|---|---|
| R1 | CP step sizes σ τ ‖A‖² < 1 too conservative → slow per-iter on kk1204 | H | M | Tight ‖A‖ via 50-iter power method; safety factor 1.05; allow σ/τ override via .Call args; report ergodic vs last-iterate max_err separately |
| R2 | IPM Schur $S = A^\top H_w^{-1} A$ rank-deficient on overlapping margins | M | M | TSVD ratio 1e-8 (mirror WL-1); diagnostic `n_projected_dims_max`; logged each outer iter |
| R3 | CP and IPM both converge to same intrinsic 5e-2 fixed point (basin floor fundamental) | M | H | Acceptable outcome — closes ylsy as "kk1204 genuinely hard". Falsification criterion: 1% agreement across all 4 solvers (CP, IPM, ieppa+sraa, newton_kl). |
| R4 | Faithful implementation underperforms; production-grade variants needed | L-M | M | Decision rule allows PARTIAL → triggers tuning follow-up (chain depth cap = 1) |
| R5 | `research/` not excluded → `R CMD INSTALL` picks up files, breaks build | L (post-relocation) | H | Mechanical CI gate: `tools/check_research_isolation.R` runs after `R CMD INSTALL --preclean .`, fails if `cp_*` or `ipm_*` symbols appear in `src/leafblower.so`. Pre-commit hook wires this. WU-1 includes green-bar test. |
| R6 | n=1M kk1204 OOMs; trace storage explodes | L | M | Pre-flight memory check; trace = scalar aggregates only (≤ 10 doubles/iter, capped at 1000 rows via stride); sparse A in CSR, ~240 MB |
| R7 | Walltime noise on 30s gate (thermal throttling, CPU governor) | M | L | 3× repetition with median for kk1204; log hardware metadata; OMP_NUM_THREADS=1 |
| R8 | Eigen / LAPACK linkage in standalone Makefile diverges from `Makevars.in` | L | L | Mirror flags exactly; lock via comment header in research/Makefile; pin Eigen via header path or vendored snapshot |
| R9 | Power iteration on $\|A\|$ fails to converge → CP step sizes wrong → divergence | L | H | Stop on $\|\lambda_{k+1} - \lambda_k\| / \lambda_k < 10^{-6}$ OR k=50; if not converged AND relative delta > $10^{-3}$ at k=50, return `status_code=4` (no CP run) |
| R10 | BLAS thread contention distorts walltime; CP extrapolation in single-precision = catastrophic cancellation | M | L | OMP_NUM_THREADS=1; RhpcBLASctl::blas_set_num_threads(1); all state in `double` (no float anywhere); spec asserts in WU-2/WU-5 implementer prompt |
| R11 | PROTECT/UNPROTECT imbalance corrupts R GC, crashes session | L | H | R/C boundary contract Sec 2; reviewer audits PROTECT count manually before commit; ASan+UBSan target in research/Makefile for one-time pre-commit smoke run |

**Discontinuation triggers:**
- R5 fires (build break) → halt, fix isolation, do not proceed.
- R9 fires (power iter divergence) → CP unimplementable; skip to IPM.
- CP and IPM both hit R3 (joint 5e-2 fixed point) → write report, close ylsy.

## 7. Final Deliverable

`docs/investigations/2026-05-02-ylsy-cp-ipm-spike-result.md` containing:

1. **Verdict per prototype**: PASS / PASS-kk1204-specialist / PARTIAL / FAIL against the decision rule.
2. **Quantitative table** across (CP, IPM, ieppa+sraa, newton_kl, lbfgsb) × (t1_small, stepstone_K9, kk1204_K20): max_err (last + ergodic for CP), walltime, n_iter, rate exponent, rate $R^2$, parity ratio, status code.
3. **Trajectory plots**: max_err vs iter, log-log for rate fitting; both last-iterate and ergodic for CP. Trace CSVs included as artefacts so reader can re-fit.
4. **Hypothesis postmortem if FAIL**: which assumption broke (basin floor, conditioning, stepsize, bound activity, rank deficiency, init point)?
5. **Recommendation**: productionize CP / productionize IPM / productionize for kk1204-specialist routing only / tune-and-retry / alternative algorithm / close ylsy BLOCKED.

### Reversibility / FAIL-artefact policy

On FAIL or PARTIAL-without-followup verdict:
- `research/` code STAYS on master under `.Rbuildignore` for traceability and future-spike baselining (no orphan deletion).
- Investigation report explicitly marks each prototype as ABANDONED with `git_sha` of last commit.
- A `bd` ticket records the closure (status FAIL → ylsy stays open with comment; status BLOCKED → ylsy closes BLOCKED).
- Future spikes can resume from the captured trace CSVs without re-implementation — that is the value of committing FAIL artefacts.

On PASS verdict:
- File Epic-I (productionization plan: method enum, full integration, AUTO routing, tests, docs).
- `research/` code remains as authoritative reference until productionization fully replaces it; only then archived.

## 8. Implementation Phases (Atomic Work Units)

One bd ticket per WU. No bundling.

| WU | Title | Deps | Scope | Decision Gate |
|---|---|---|---|---|
| **WU-1** | Skeleton + isolation enforcement | — | Create `research/` directory, `research/Makefile`, `tools/check_research_isolation.R`, `.Rbuildignore` lines. Stub `cp_solve_R` returning `list(status_code=99)`. R CMD INSTALL --preclean . → run isolation gate → assert no leak. | Build clean; gate green. |
| **WU-2** | CP implementation | WU-1 | `research/cp_calib.{hpp,cpp}` + `research_bridge.cpp` `cp_solve_R`. Power iter, prox Newton, full algorithm per Sec 3.1. Compiles standalone. | Compiles + sanity recovers w=1 on t1_small to max_err < 1e-8. |
| **WU-3** | CP bench | WU-2 | `benchmarks/research/ylsy_cp_bench.R` + `sanity_t1_recovery.R`. Run on t1_small, stepstone_K9, kk1204_K20 (3×). Produce summary CSV + trace CSVs. | All three fixtures complete; CSVs written. |
| **WU-4** | CP gate | WU-3 | Read CP CSVs; compute verdict via decision rule. Output `research/cp_verdict.txt` (PASS / PASS-specialist / PARTIAL / FAIL). | Verdict set. If PASS → skip WU-5/WU-6. |
| **WU-5** | IPM implementation | WU-1, WU-4 ≠ PASS | `research/ipm_calib.{hpp,cpp}` + shim. Schur+TSVD per Sec 3.2. Sanity recovers w=1 to max_err < 1e-8. | Compiles + sanity passes. |
| **WU-6** | IPM bench | WU-5 | `benchmarks/research/ylsy_ipm_bench.R`. Same fixtures, same repetitions. CSVs. | All three fixtures complete. |
| **WU-7** | Investigation report | WU-3 OR WU-6 (whichever final) | `ylsy_compare.R` produces unified table. Write `docs/investigations/2026-05-02-ylsy-cp-ipm-spike-result.md` with verdict, table, plots, postmortem, recommendation. Update bd ylsy ticket. | Report committed; ticket updated; epic closed (PASS / FAIL / PARTIAL). |

### Subagent routing

- WU-1 (mechanical scaffolding, isolation tests): **Haiku**.
- WU-2, WU-3, WU-5, WU-6 (implementation, bench): **Gemini** (per memory rule — Tier 2 delegated work to Gemini, not Sonnet).
- WU-4 (gate logic, verdict computation): **Haiku** (mechanical: read CSV, apply rule).
- WU-7 (analysis, postmortem): **Opus** (Tier 3 — judgment, recommendations).
- This design (rev 2) + plan-review-gate + adversarial review of investigation report: **Opus**.

### Sequential learnings

If WU-4 verdict ≠ PASS, IPM design (Sec 3.2) is FROZEN at this spec — no implicit redesign during WU-5. If CP exposes a pathology that invalidates the IPM math (e.g., overlapping-margin case showing different rank structure), IPM design changes get a separate gate (rev 3 of this spec + new plan-review-gate run). No silent hyperparameter sweeps.

## 9. Predecessor & Future Context

- **Predecessor**: Epic-H WH-g landed on master at `3265a53` (AUTO routes severe-skew K≥5 to ieppa+sraa). Verified live via `git log --oneline 3265a53^..3265a53` before merging this spec rev. Empirical kk1204 verify: `bd memory kk1204-k-20-wh-g-empirical-auto-correctly`.
- **Reference fixed-point analysis**: `bd memory fixed-point-accelerators-apva-anderson-joint-anderson-on` documents why APVA / Joint Anderson / Halpern fail on piecewise-linear bound problems — the analysis CP and IPM are designed to bypass.
- **Stepstone K=9 basin floor**: `2.61e-4` (Epic-Dβ NEWS PARTIAL verdict). CP/IPM may also help stepstone, but stepstone is NOT this spike's primary gate; productionization is.

## 10. Q&A History (resolves CTO Q5 reference)

The design was brainstormed via `superpowers:brainstorming` with 6 multiple-choice questions:

- **Q1**: Scope — single method / multi / spike. Picked **C (spike)**.
- **Q2**: Candidate pair — picked **A (CP + IPM)**.
- **Q3**: Implementation language — picked **B (C++ standalone via Rcpp/cpp11)**. Note: rev 2 relocates from `src/research/` to package-root `research/` per Architect blocker.
- **Q4**: Success criteria — picked **A (kk1204 < 1e-3 in 30s + parity + rate ≥ 1/k)**. Refined in rev 2 with R² floor + last-iterate β.
- **Q5**: Code location + commit policy — picked **A (commit to master under research/, including FAIL artefacts)**.
- **Q6**: Implementation order — picked **A (CP first → IPM if needed)**.

Both Architect blocker findings (R5 broken; objective mismatch) and Designer/Security/CTO blockers were folded into this rev 2.
