# leafblower (development version)

## Hierarchical calibration

* `harvest()` gains `hierarchical = list(...)` for two-stage coarse-then-fine
  calibration on RAKING. Coarse margins are calibrated globally (Stage 1); fine
  margins are then calibrated within each coarse cell (Stage 2). Use
  `mode = 0L` (refine, default) for iterative outer convergence, or
  `mode = 1L` (exact) for a single-pass orthogonal split. Sparse cells
  (fewer than `min_cell_n` observations) inherit Stage-1 multipliers.

## Breaking changes

* `harvest()` default changed: `sor` is now `NULL` (disabled) instead of
  `list(auto = TRUE, omega_min = 0.3)`. SOR caused 2–3× slowdowns at loose
  bounds (`max_weight ≥ 5`, `K ≤ 3`) — the common calibration scenario.
  To opt in: `harvest(..., sor = list(auto = TRUE))`.

* `anesrake()` now defaults to `choosemethod = "ieppa"` (was `"rake"` which triggered a deprecation warning).

## Newton-KL calibration

* Epic-Dβ verdict: **PARTIAL** — Newton-KL pipeline shipped with LM scale-invariant
  damping (Marquardt gain ratio), truncated-SVD pseudoinverse (LAPACK `dsyevd`),
  and Steihaug-CG trust-region in retained subspace. LM damping adapts via gain-ratio
  ρ (ρ>0.75 ⇒ Δ ×2; ρ<0.25 ⇒ Δ ÷4), with LSE stabilization in dual evaluation
  and weight recovery (fixes NaN weights on K=20 severe-skew). Best-iterate fallback
  rolls back rank-deficient drift; `LEAFBLOWER_NEWTON_TRACE` audit shows fires
  across all scenarios (25 T1/T4/T5, 3 T7-K4, 9 T2-stepstone, 37 regression).
  Diagnostic fields: `n_projected_dims` (truncation ratio `1e-8 × λ_max`), `lm_mu_final`.
  Stepstone K=9 basin floor at ~2.6e-4 appears intrinsic to dual landscape
  (from 2.79e-4 master, 6.5% improvement). Full regression FAIL=0 outside
  documented basin; T7 K=4 over-projection PASS (2.84e-7); kk1204 K=20 severe-skew
  converges status=0 at gap=6.24e-2 (vs master diverged). Closure to <1e-4
  deferred to Epic-E (continuation methods, multi-start, alternative algorithm).

## Breaking changes

* `harvest()` parameter `jacobi_sweep` removed. The parameter has been
  C++-inert since its introduction (the underlying sweep order is always
  Gauss-Seidel). Callers passing `jacobi_sweep=` will have it silently
  absorbed by `...` (no behavior change vs prior silent no-op).

* `harvest(method="greenkhorn", accelerate=TRUE)` and `harvest(method="raking", accelerate=TRUE)`:
  calibrated weights will differ from previous versions. The prior SQUAREM/CBB
  acceleration overshot the bounded optimum (greenkhorn: 35% worse quality than plain).
  Replaced by Safeguarded Regularized Anderson Acceleration (SRAA-m, m=5), which
  guarantees quality >= plain per super-step. Reproducible pipelines using
  `set.seed() + accelerate=TRUE` will produce different (more accurate) results.

* **AUTO routing for K≥5 severe-skew problems** (`max_T/min_T > 5`): `harvest(method="auto", ...)`
  now selects `method="ieppa"` with `accelerate=TRUE` instead of `method="newton_kl"`.
  Motivated by Epic-Dβ kk1204 K=20 evidence (Newton-KL converges to high-error fixed point
  at gap=6.24e-2 on severe-skew). Moderate-skew K≥5 (max_T/min_T ≤ 5) still routes to
  Newton-KL. Affected users: AUTO callers with K≥5 categorical strata and order-of-magnitude
  target imbalance. Migration: pin `method="newton_kl"` explicitly to retain prior behavior.

## Acceleration

* Replaced SQUAREM/CBB with SRAA-m (Type II Anderson Acceleration, window m=5) for
  `method="greenkhorn"` and `method="raking"`. SRAA-m fixes the 35%/9% quality
  regression on bounded calibration problems. Uses 2 F-evals/accepted step vs
  SQUAREM's 3. Automatic restart on residual divergence (||R_k||^2 > 4*||R_{k-1}||^2)
  and Tikhonov regularization on the m*m Gram system.

---

# leafblower (development)

## Breaking changes

* **Convergence API redesigned** — `criterion` key replaced by `metric` + `rule`:
  - **Default changed** to `metric = "max_err"`, `rule = "improvement"`, `tol = 0.001`
    (was: `pct = 0.001`, i.e., l1_weight plateau). To restore the old `pct` default:
    `convergence = list(pct = 0.001)` (still accepted as a shorthand).
  - `convergence[["criterion"]]` is removed. Replace with `metric`.
  - `convergence[["pct"]]` is now a shorthand for `metric="l1_weight"` +
    `rule="plateau"`. Previously it controlled an L1 weight-change threshold.
  - Valid `metric` values: `"max_err"`, `"mean_err"`, `"kl"`, `"chi2"`,
    `"grake_norm"`, `"l1_weight"`.
  - Valid `rule` values: `"improvement"`, `"threshold"`, `"plateau"`.
  - `tol` replaces the implicit tolerance in `pct`/`absolute` shortcuts.
  - Python `harvest()` convergence API updated in lockstep.

* **`pct_change` removed from result** — renamed to `l1_weight_change` to
  match the field name in the C struct (`rk_result_t::l1_weight_change`).
  Update any code reading `attr(result, "result")$pct_change` or
  `df.attrs["result"]["pct_change"]` to use `l1_weight_change`.

* `method="ieppa"` now runs the paper-faithful algBCD (Chu, Liang, Toh &
  Yang 2022, arXiv:2011.14312) at C=0, using cell-compressed representation
  with log-space Sinkhorn factors and a capacity BCD block. The previous
  implementation was an IPF+Dykstra hybrid misnamed "iEPPA"; it is renamed
  `method="raking"`.

* New `method="raking"` exposes the renamed classical IPF+Dykstra hybrid
  (Deming-Stephan 1940 / Csiszár 1975 cyclic IPF + Boyle-Dykstra 1986
  additive projections). Same implementation as pre-rename `method="ieppa"`.

* `method="auto"` continues to route to `method="ieppa"` (now the faithful
  algBCD).

## New features

* Pluggable convergence via `convergence = list(metric = "kl", rule = "threshold", tol = 1e-4)`:
  six metrics (`max_err`, `mean_err`, `kl`, `chi2`, `grake_norm`, `l1_weight`),
  three rules (`improvement`, `threshold`, `plateau`). `stop_when` controls
  whether any or all conditions must fire.

* `convergence_used` nested list in `attr(result, "result")$convergence_used`
  (R) and `df.attrs["result"]["convergence_used"]` (Python) documents the
  metric, rule, tol, and `fired_at_iter` that caused the solver to stop.

* Six quality metrics always returned in `attr(result, "result")`:
  `max_error`, `mean_error`, `kl`, `chi2`, `l1_weight_change`, `grake_norm`.

* SOR adaptive under-relaxation for iEPPA via `sor` argument (default:
  auto-enabled). Diagnostics in `attr(result, "result")$sor`.

* Best-iterate tracking: `attr(result, "result")$best_weights` populated
  with weights at minimum observed marginal error.

* Python `harvest()` updated with matching `convergence` and `sor` arguments;
  new result fields exposed via `df.attrs["result"]`.

* Cell-compressed computation: faithful iEPPA operates at cell-level
  (unique (g_1,...,g_K) tuples) rather than observation-level, yielding
  up to 1000× speedup on surveys with low tuple diversity.
