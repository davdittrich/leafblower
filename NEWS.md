# leafblower (development version)

## Newton-KL calibration

* `method="newton_kl"` adds Steihaug-CG trust-region step in retained subspace
  when the truncated-SVD pseudoinverse exceeds the trust radius. Trust radius
  adapts via Marquardt gain ratio (ρ>0.75 ⇒ Δ doubles; ρ<0.25 ⇒ Δ /= 4).
* `method="newton_kl"` now applies truncated-SVD pseudoinverse to the LM-damped
  Hessian to handle rank-deficient Hessian directions on overlapping-margin
  fixtures. Eigendecomposition via LAPACK's `dsyevd`. Default truncation ratio
  `1e-8 × λ_max`. New diagnostic `n_projected_dims` on result.
* `method="newton_kl"` now uses Levenberg-Marquardt scale-invariant damping
  with Marquardt gain-ratio adaptive μ (commits 88033d1 + e4276cc). Replaces
  the prior trust-region clip + pivot-bump bandaid. Adds LSE stabilization in
  the dual evaluation and weight recovery (fixes NaN weights on K=20 severe-skew
  inputs that previously short-circuited convergence via false gap=0). Adds
  best-iterate fallback so rank-deficient drift past the optimum is auto-rolled
  back to the lowest-gap λ. New diagnostic field `lm_mu_final` on the C++
  `NewtonCalibResult` (R-side surfacing pending in a follow-up).
* Stepstone K=9 max_err: 0.988 (regression from rev-1 LM) → 2.8e-4 (final).
  Still above the 1e-4 gate; closure depends on Epic-B (target homotopy) —
  see `docs/investigations/2026-05-01-newton-kl-lm-result.md`.

## Breaking changes

* `harvest(method="greenkhorn", accelerate=TRUE)` and `harvest(method="raking", accelerate=TRUE)`:
  calibrated weights will differ from previous versions. The prior SQUAREM/CBB
  acceleration overshot the bounded optimum (greenkhorn: 35% worse quality than plain).
  Replaced by Safeguarded Regularized Anderson Acceleration (SRAA-m, m=5), which
  guarantees quality >= plain per super-step. Reproducible pipelines using
  `set.seed() + accelerate=TRUE` will produce different (more accurate) results.

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
