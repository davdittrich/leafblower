# ylsy: CP + IPM Research Spike — Design

**Date:** 2026-05-02
**Beads:** `leafblower-ylsy` (kk1204 convergence research)
**Status:** Design — pre-plan
**Predecessor:** Epic-H WH-g (AUTO routing safety; routing-only fix; kk1204 still gap=5.01e-2 under ieppa+sraa)

## 1. Problem & Goals

### Problem

The kk1204 fixture (n=1M, K=20 dense margins, 5 categories per margin, target=(0.6, 0.2, 0.1, 0.07, 0.03), max_weight=3, target_skew=20) defeats every solver in `harvest()`:

| Solver (master, post-Epic-H) | kk1204 max_err | status | wall |
|---|---|---|---|
| `newton_kl` | 6.24e-2 | 1 (NOCONV) | — |
| `ieppa` + `accelerate=TRUE` (auto via WH-g) | 5.01e-2 | 0 (converged) | 11.4s |
| `lbfgsb` (auto-fallback) | 0.31 | 4 | 372s |

The PARTIAL gate target (max_err < 1e-3) is unreached.

Underlying cause: kk1204 is PURE slow-rate Sinkhorn (errRp ≈ iter^{-0.5}, monotone, no oscillation; per `bd memory` `fixed-point-accelerators-apva-anderson-joint-anderson-on`). All explored fixed-point accelerators (APVA, Joint Anderson, Halpern) failed because piecewise-linear bound projections defeat smooth contraction theory.

### Goal

Add (do not replace) one or two new calibration methods that break the kk1204 ceiling. Two candidate algorithms with orthogonal paradigms:

- **Chambolle-Pock primal-dual (PDHG)** — first-order saddle-point method designed for non-smooth (LP/QP-class) problems. Avoids the contraction-theory failure mode.
- **Interior-Point Newton (IPM)** — second-order central-path method. Log-barrier smooths the bounds; Newton handles strong convexity of KL.

### Decision Rule

| Outcome on kk1204 | Verdict | Action |
|---|---|---|
| max_err < 1e-3 AND walltime ≤ 30s AND stepstone max_err ≤ ieppa+sraa × 1.5 AND rate exponent ≤ −1.0 (i.e. err decays at least as 1/k) | **PASS** | File Epic-I (productionize: method enum + tests + AUTO + docs) |
| max_err in [1e-3, 1e-2) AND walltime ≤ 60s | **PARTIAL** | Document. Possibly tune hyperparams in follow-up spike. |
| max_err ≥ 1e-2 OR walltime > 60s OR stepstone max_err > ieppa+sraa × 1.5 | **FAIL** | Investigation report. ylsy stays open or closes BLOCKED. |

**Stepstone reference baseline** (from `bench memory` and Epic-Dβ NEWS): ieppa+sraa achieves max_err ≈ 1.13e-4 on stepstone K=9 fulldata. Parity tolerance × 1.5 → cutoff ≈ 1.7e-4.

**Rate exponent**: slope $\beta$ of $\log(\mathrm{err}) = \alpha + \beta \log(\mathrm{iter})$. $\beta = -0.5$ corresponds to sinkhorn slow-rate (the problem). $\beta \le -1$ corresponds to 1/k or faster (PDHG ergodic). Fit on iters ∈ (50, 0.9 × n_iter) with reported $R^2$.

CP runs first (cheaper to implement). IPM runs only if CP < PASS (sequential per Q6). If CP PASS: skip IPM in this spike.

## 2. Architecture & Isolation Contract

```
src/research/                          (NEW — excluded from PKG_SOURCES)
├── cp_calib.hpp / cp_calib.cpp
├── ipm_calib.hpp / ipm_calib.cpp
├── research_bridge.cpp                (.Call shim: cp_solve_R, ipm_solve_R)
└── Makefile                           (standalone: R CMD SHLIB)

benchmarks/research/
├── ylsy_cp_bench.R
├── ylsy_ipm_bench.R                   (only if CP < PASS)
├── ylsy_compare.R
└── results/                           (CSVs)

docs/investigations/
└── 2026-05-02-ylsy-cp-ipm-spike-result.md
```

**Isolation invariants** (verified before commit):

- `src/research/Makefile` standalone. Mirrors `src/Makevars.in` flags (Eigen path, LAPACK linkage) exactly.
- `R CMD INSTALL --preclean .` does **not** pick up `src/research/*`. Either: not listed in `PKG_SOURCES`, or under `.Rbuildignore`. Verify by running install after committing skeleton.
- Bench scripts call `dyn.load("src/research/leafblower_research.so")` + `.Call`.
- Zero coupling to `leafblower.h` enum, `r_bridge.cpp`, `harvest.R`, `tests/testthat/`.
- I/O surface: `(df, tgt, max_weight, min_weight, max_iter, capture_trace)` → `(weights, trace_matrix)`.

## 3. Algorithm Specifications

### 3.1 Chambolle-Pock (PDHG)

**Reference:** Chambolle & Pock 2011 (J. Math. Imaging Vision 40:120–145), Algorithm 1.

**Problem (saddle-point form):**

$$\min_{w \in [\ell, u]} \; f(w) + \langle A^\top w - b, y \rangle + I_{*}(y)$$

with $f(w) = \sum_i d_i \, \mathrm{KL}(w_i / d_i \| 1)$ on the box $[\ell, u]$, and $A \in \mathbb{R}^{n \times \sum_k J_k}$ the block-incidence design matrix (one column per (margin $k$, level $j$)). $b = T \cdot Z$ is the target marginal mass vector, with $Z = $ design-weighted sample size.

**Iteration:**

```
Init: x⁰ = 1ₙ, x̄⁰ = 1ₙ, y⁰ = 0
||A|| computed once via 50-iter power method on A^T A.
Step sizes: σ = τ = 1 / ||A||, ensuring σ τ ||A||² = 1 (boundary of theory).
For k = 0 .. max_iter:
  y^{k+1} = y^k + σ (A^T x̄^k − b)             # g* linear, no prox
  x^{k+1} = prox_{τ f}(x^k − τ A y^{k+1})      # KL+box prox per element
  x̄^{k+1} = 2 x^{k+1} − x^k                    # extrapolation
```

**`prox_{τ f}` per component:** scalar Newton (≤5 iters) on
$$\arg\min_{w \in [\ell_i, u_i]} \tau d_i (w \log(w/d_i) - w + d_i) + 0.5(w - z_i)^2$$
followed by clamp to $[\ell_i, u_i]$. Closed form initial guess: $w_0 = \mathrm{lambertW}(\tau d_i \exp(z_i / (\tau d_i) - 1))$ for unbounded case; one Newton step typically suffices.

**Convergence:** O(1/k) ergodic on primal-dual gap. (Accelerated $O(1/k^2)$ variant deferred to productionization.)

**Per-iter cost:** 2 sparse matvecs ($A x̄$, $A^\top y$) + n scalar prox ops = $O(n K)$.

**Stop:** $\|A^\top x - b\|_\infty / Z < 10^{-7}$, OR $\|x^{k+1} - x^k\|_\infty < 10^{-12}$, OR `max_iter`.

**Diagnostic outputs:** per-iter `(time_ms, max_err, primal_resid = ||A^T x − b||_∞, dual_resid = ||x − prox|| / τ)`.

### 3.2 Interior-Point Newton (IPM)

**Problem:**
$$\min_w \; \sum_i (w_i \log w_i - w_i) \quad \text{s.t.} \quad A^\top w = b, \; \ell \le w \le u$$

**Log-barrier:** $\Phi(w) = -\mu \sum_i [\log(w_i - \ell_i) + \log(u_i - w_i)]$.

**KKT residual** with margin multiplier $\lambda$:

$$r_w(w, \lambda) = \log w + \mu \left(\frac{1}{u - w} - \frac{1}{w - \ell}\right) + A \lambda$$
$$r_\lambda(w) = A^\top w - b$$

**Newton step (Schur on diagonal Hessian):**

```
H_w = diag(1/w + μ/(w−ℓ)² + μ/(u−w)²)        (n diag entries)
S   = A^T H_w⁻¹ A                            (m × m, m = ΣJ_k ≈ 100)
S   = TSVD(S, ratio=1e-8)                    (handles rank deficiency, mirror Epic-Dβ WL-1)
Δλ  = S⁻¹ (A^T H_w⁻¹ r_w − r_λ)
Δw  = −H_w⁻¹ (r_w + A Δλ)
α   = max α ∈ (0, 1] s.t. ℓ + ε ≤ w + α Δw ≤ u − ε   (fraction-to-boundary; ε = 1e-12)
α   ← min(α, 0.99 × α_max)                            (strict feasibility)
w   ← w + α Δw,  λ ← λ + α Δλ
```

**Outer loop:** $\mu_{k+1} = 0.5 \mu_k$, $\mu_0 = 1$, until $\mu < 10^{-9}$.
**Inner loop:** 1–3 Newton steps per $\mu$ until $\|r_w\|_\infty + \|r_\lambda\|_\infty < 10 \mu$.

**Stop:** $\|r_w\|_\infty + \|r_\lambda\|_\infty < 10^{-7}$ AND $\mu < 10^{-9}$, OR outer iter $\ge 50$.

**Per-iter cost:** $O(n K)$ for matvecs + $O(m^3)$ for dense $m \times m$ system, $m \approx 100$ — Schur cheap.

**Diagnostic outputs:** per-(outer,inner) `(time_ms, mu, max_err, kkt_resid, n_projected_dims, alpha_step)`.

### 3.3 Faithful-Textbook Constraint

No hyperparameter tuning, no momentum / restart heuristics, no warm-start, no accelerated variants in the spike. Default config only:

- CP: $\sigma = \tau = 1/\|A\|$, no acceleration.
- IPM: $\mu_0 = 1$, decay $0.5$, vanilla central path (no Mehrotra predictor-corrector).

Tuning lands in productionization phase if spike PARTIAL.

## 4. Bench Harness

### Fixtures

| Fixture | Source | Role |
|---|---|---|
| `t1_small` | inline (n=1000, K=3) | Sanity — any solver hits machine epsilon |
| `stepstone_K9` | `benchmarks/stepstone_bench_data.parquet` + targets JSON | Parity guard against ieppa+sraa baseline |
| `kk1204_K20` | inline (n=1M, K=20, target=(0.6, 0.2, 0.1, 0.07, 0.03)) | Primary gate |

### Per-fixture captured fields

| Field | Description |
|---|---|
| `iter` | iteration index |
| `time_ms` | cumulative wall time (ms) |
| `max_err` | $\max_{k,j} \lvert \hat{T}_{kj} - T_{kj} \rvert / \max(T_{kj}, 1\mathrm{e}-12)$ |
| `primal_resid`, `dual_resid` | CP only |
| `kkt_resid`, `mu`, `alpha_step` | IPM only |
| `n_projected_dims` | IPM TSVD truncation count |

### Combined comparison

`benchmarks/research/ylsy_compare.R` produces a single CSV `ylsy_comparison.csv`:

| solver | fixture | max_err | wall_s | n_iter | rate_exponent | rate_R² |
|---|---|---|---|---|---|---|

`rate_exponent` = slope of `lm(log(max_err) ~ log(iter), data = trace[trace$iter > 50, ])`. Drop final 10% of trace to avoid boundary artefacts.

Optional plots saved to `benchmarks/research/results/plots/` (log-log trajectory per fixture, per solver).

### Reproducibility

- All fixtures use `set.seed(1)` matching existing benches (`benchmarks/tsvd_kk1204_K20.R`, `benchmarks/stepstone_bench_data.parquet`).
- `Sys.setenv(OMP_NUM_THREADS = "1")` for fair walltime.

## 5. Out of Scope

| Item | Reason |
|---|---|
| `harvest()` / `r_bridge.cpp` modifications | Productionization-only |
| `leafblower.h` enum (`RK_ALG_CP`, `RK_ALG_IPM`) | Same |
| AUTO routing changes | Same |
| Test additions in `tests/testthat/` | Spike validates algorithm; production tests later |
| `NEWS.md` / spec changes | Investigation report is the artefact |
| Python bindings | Productionization-only |
| Multi-threading / SIMD optimisation | Faithful textbook only; tuning later |
| Hyperparameter sweeps | Default config only |
| Warm-start from ieppa / newton_kl | Cold start only |
| Accelerated CP variant (Algorithm 2 of Chambolle-Pock) | Standard PDHG only |
| Predictor-corrector IPM (Mehrotra) | Vanilla central path only |
| Cell compression in CP / IPM | kk1204 has M_cell/n=1.0; compression buys nothing |
| Convergence improvements to existing methods | New methods only, in `src/research/` |

## 6. Risks & Mitigations

| # | Risk | L | I | Mitigation |
|---|---|---|---|---|
| R1 | CP step sizes σ τ ‖A‖² < 1 too conservative → slow per-iter on kk1204 | H | M | Tight ‖A‖ via 50-iter power method; allow σ/τ override via .Call args; report ergodic vs last-iterate err separately |
| R2 | IPM Schur $S = A^\top H_w^{-1} A$ rank-deficient on overlapping margins | M | M | TSVD ratio 1e-8 (mirror WL-1); diagnostic `ipm_n_projected_dims` |
| R3 | CP and IPM both converge to same intrinsic 5e-2 fixed point (basin floor fundamental) | M | H | Acceptable outcome — closes ylsy as "kk1204 is genuinely hard". Joint barrier documented. |
| R4 | Faithful implementation underperforms; production-grade variants (acc-PDHG, Mehrotra) needed | L-M | M | Decision rule allows PARTIAL → triggers tuning follow-up |
| R5 | `src/research/` not excluded → `R CMD INSTALL` picks up files, breaks build | L | H | Explicit `R CMD INSTALL --preclean .` smoke test after committing skeleton; verify also via `.Rbuildignore` if needed |
| R6 | n=1M kk1204 OOMs without cell compression | L | M | Pre-flight: n × K × 8 = 160 MB scratch — fits comfortably on 16 GB |
| R7 | Rate fitting noise-dominated for low iter counts | L | L | Fit only iters > 50, drop final 10%; report R² |
| R8 | Eigen / LAPACK linkage in standalone `Makefile` diverges from `Makevars.in` | L | L | Mirror flags exactly; lock via comment header in Makefile |

**Discontinuation triggers:**
- R5 fires (build break) → halt, fix isolation, do not proceed.
- CP and IPM both hit R3 (joint 5e-2 fixed point) → write report, close ylsy.

## 7. Final Deliverable

`docs/investigations/2026-05-02-ylsy-cp-ipm-spike-result.md` containing:

1. **Verdict per prototype**: PASS / PARTIAL / FAIL against the decision rule.
2. **Quantitative table** across (CP, IPM, ieppa+sraa, newton_kl, lbfgsb) × (t1_small, stepstone_K9, kk1204_K20): max_err, walltime, n_iter, rate exponent, rate $R^2$.
3. **Trajectory plots**: max_err vs iter, log-log for rate fitting.
4. **Hypothesis postmortem if FAIL**: which assumption broke (basin floor, conditioning, stepsize, bound activity, rank deficiency)?
5. **Recommendation**: productionize CP / productionize IPM / tune-and-retry / alternative algorithm / close ylsy BLOCKED.

## 8. Implementation Order (Subagent Strategy)

CP first → measure → if not PASS, IPM second → measure → write report.

Both prototypes implemented by **Gemini** (per global routing rule until further notice). **Opus** for design (this doc) + plan-review-gate + adversarial review of spike result. Haiku for mechanical bench-harness scaffolding if needed.

No parallel execution. Sequential learnings: IPM design adjusts on CP findings if CP fails (e.g., if CP exposes overlapping-margin pathology, IPM TSVD ratio may need adjustment).

## 9. Predecessor & Future Context

- **Predecessor**: Epic-H WH-g landed on master at `3265a53` (AUTO routes severe-skew K≥5 to ieppa+sraa). Empirical kk1204 verify shipped as `bd memory` `kk1204-k-20-wh-g-empirical-auto-correctly`.
- **Future**: Stepstone K=9 basin floor (`2.61e-4`) is a separate investigation — Epic-E candidate. CP/IPM may also help stepstone, but stepstone is NOT this spike's primary gate.
- **Reference fixed-point analysis**: `bd memory fixed-point-accelerators-apva-anderson-joint-anderson-on` documents why APVA / Joint Anderson / Halpern fail on piecewise-linear bound problems — the analysis CP and IPM are designed to bypass.
