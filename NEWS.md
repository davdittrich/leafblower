# leafblower 0.2.1 (2026-05-07)

Remediation wave for P1 epic `leafblower-6ycz.1`. Bug fixes only; no new public
API. Covers 11 tickets (T-M through T-W) identified during T-K code-review pass.

## Bug fixes

* Stage-1 hierarchical calibration now fits coarse margins only (was: all
  margins). Affected RAKING, SINKHORN, and GREENKHORN (`leafblower-6ycz.1.13`,
  T-M).
* Cell-aggregate Stage-1 multiplier replaces last-observation overwrite for
  sparse-cell weight inheritance (`leafblower-6ycz.1.13`).
* `outer_residual_final` now reports spec §8 `Σ|w·X − target|/N` L1 metric
  (was: cell-balance proxy). BUDGET-exit recomputes residual on returned weights
  for caller-visible consistency (`leafblower-6ycz.1.15`, T-O).
* `hier_levels_used` returns 2 when hierarchical calibration is enabled (was: 1)
  (`leafblower-6ycz.1.14`, T-N).
* `outer_iterations_used` returns -1 in `mode = "exact"` branches (was: 1)
  (`leafblower-6ycz.1.17`, T-Q).
* Σw=N invariant violations are now surfaced via `result$message` and an
  `Rprintf` warning; new `enforce_sigmaw_eq_n_diag` overload returns
  `(passed, dev)` (`leafblower-6ycz.1.16`, T-P).
* Python `HierarchicalConfig` rejects non-finite `outer_tol` via
  `math.isfinite` check (`leafblower-6ycz.1.21`, T-U).

## Internal

* Drop duplicate `LBW_MAX_HIER_CELLS` define from `validation.hpp`; canonical
  site is `calib_dispatch.hpp` (`leafblower-6ycz.1.18`, T-R).
* Right-size `fine_gids_buf` / `fine_weights_buf` to `max_n_per_cell`
  (~196 KB saved per call) (`leafblower-6ycz.1.19`, T-S).
* `CalibState` sub-state uses explicit field-by-field init (was: `sub = st;`
  raw pointer aliasing) (`leafblower-6ycz.1.20`, T-T).
* `lbw::CalibResult` gained `char message[256]` field (`leafblower-6ycz.1.16`).

## Test / docs

* Spec §8 amended to v4 (v1–v3 history preserved) post-T-M empirical re-test
  (`leafblower-6ycz.1.23`, T-W); v5 post-T-V stepstone-gate decision
  (`leafblower-6ycz.1.22`).
* Stepstone regression-gate fixture regenerated at T-M + T-O baseline;
  `RESID_SLACK = 0.05` retained — structural floor on 9-margin overlapping
  geometry (`leafblower-6ycz.1.22`).
* Vignette Stage-1 description corrected to coarse-only semantics
  (`leafblower-6ycz.1.13`).

## Known limitations

* Spec §8 P1 sparse-cell rescue seed-sweep gate remains `skip()`'d in
  `tests/testthat/test-2stage-raking.R`. DGP discovery deferred to follow-up
  `leafblower-6ycz.1.12` (T-L; closed with disposition note: 2-stage rescues
  sparseness, not bounds-infeasibility — none of v0/v1/v2/v3 DGPs achieve rescue
  at spec §8 thresholds under the corrected algorithm).

# leafblower 0.2.0 (2026-05-07)

P1 epic `leafblower-6ycz.1` — hierarchical two-stage calibration for RAKING,
SINKHORN, and GREENKHORN. Backward-compatible: `hierarchical = NULL` (default)
preserves single-stage behavior identical to 0.1.x. Issues: `leafblower-6ycz`.

## New features

* `harvest()` gains `hierarchical = list(coarse_margins, fine_margins, ...)` for
  two-stage coarse-then-fine calibration. Coarse margins are calibrated globally
  (Stage 1); fine margins are calibrated within each coarse cell (Stage 2). Sparse
  cells (fewer than `min_cell_n` observations) inherit Stage-1 multipliers unchanged.
  Supported for `method = "raking"`, `"sinkhorn"`, and `"greenkhorn"`.

* `harvest()` gains `mode` argument controlling Stage-2 solver strategy:
  `mode = "refine"` (default, `0L`) — iterative outer convergence;
  `mode = "exact"` (`1L`) — single-pass orthogonal split.

* Six convergence/quality diagnostics always returned in `attr(result, "result")`:
  `max_error`, `mean_error`, `kl`, `chi2`, `l1_weight_change`, `grake_norm`.
  Two hierarchical diagnostics added: `hier_outer_iters`, `hier_outer_residual`.

* Newton-KL calibration (`method = "newton_kl"`) ships with LM scale-invariant
  damping (Marquardt gain ratio), truncated-SVD pseudoinverse (LAPACK `dsyevd`),
  and Steihaug-CG trust-region in retained subspace. Fixes NaN weights on K=20
  severe-skew. Diagnostic fields: `n_projected_dims`, `lm_mu_final`.

* Replaced SQUAREM/CBB with SRAA-m (Safeguarded Regularized Anderson Acceleration,
  window m=5) for `method = "greenkhorn"` and `method = "raking"`. SRAA-m guarantees
  quality >= plain per super-step (fixes prior 35%/9% quality regression on bounded
  problems). Uses 2 F-evals/accepted step vs SQUAREM's 3.

## API additions

* `hierarchical` argument: `list(coarse_margins, fine_margins, min_cell_n, mode)`.
  Python: `HierarchicalConfig` dataclass with identical fields.

* `mode` argument: `"refine"` | `"exact"` (also accepted as integer `0L` / `1L`).

* `bounds_mode = "unit"` combined with `hierarchical` raises `BADARG` immediately —
  per-observation water-fill and hierarchical cell dispatch are mutually exclusive.

* Cell count cap: `LBW_MAX_HIER_CELLS = 100000`; exceeded at dispatch raises `BADARG`.

* `anesrake()` now defaults to `choosemethod = "ieppa"` (was `"rake"`, which
  triggered a deprecation warning).

* **AUTO routing for K≥5 severe-skew** (`max_T/min_T > 5`): `harvest(method="auto")`
  now selects `method = "ieppa"` with `accelerate = TRUE`. Moderate-skew K≥5
  (ratio ≤ 5) still routes to Newton-KL. Migration: pin `method = "newton_kl"`
  to retain prior behavior.

## Breaking changes

* `harvest()` default changed: `sor` is now `NULL` (disabled). SOR caused 2–3×
  slowdowns at loose bounds (`max_weight ≥ 5`, `K ≤ 3`). To opt in:
  `harvest(..., sor = list(auto = TRUE))`.

* `harvest()` parameter `jacobi_sweep` removed (was C++-inert since introduction).
  Callers passing `jacobi_sweep =` will have it silently absorbed by `...`.

* `harvest(method="greenkhorn", accelerate=TRUE)` and
  `harvest(method="raking", accelerate=TRUE)`: calibrated weights will differ from
  0.1.x due to SRAA replacing SQUAREM/CBB. Reproducible pipelines using
  `set.seed() + accelerate=TRUE` will produce different (more accurate) results.

## Tests / CI

* Stepstone regression gate added (`tests/testthat/test-stepstone.R`): 9 adversarial
  fixtures per method (RAKING, SINKHORN, GREENKHORN) validated against stored RDS
  baselines. Gate uses a 5% slack regression band; `outer_residual_final` does not
  enforce 1e-4 absolute on stepstone (intrinsic basin floor ~2.6e-4 on K=9).

* Python parity suite (`python/parity/`) passes 50 tests (3 skipped, 0 failed)
  at `rtol = 1e-6` vs R reference output. Outer iteration counts are integer-exact.

## Known limitations

* **Sparse-cell rescue / seed-sweep gate deferred** to `leafblower-6ycz.1.12` (T-L).
  When a Stage-2 cell has fewer than `min_cell_n` observations and the inherited
  Stage-1 multiplier produces a poor seed for the within-cell solver, recovery is
  currently limited to inheriting Stage-1 weights unchanged. A DGP-discovery sweep
  to characterize failure modes is pending.

* Newton-KL K=9 stepstone basin floor (~2.6e-4) is intrinsic to the dual landscape.
  Closure to `< 1e-4` deferred to Epic-E (continuation methods, multi-start,
  alternative algorithm).

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
