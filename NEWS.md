# leafblower (development)

## Breaking changes

* **Default convergence criterion changed** from `absolute = 1e-6` (max marginal
  error) to `pct = 0.001` (max proportional weight change). To restore prior
  behaviour: `convergence = list(absolute = 1e-6)`. The Python `harvest()`
  default is updated in lockstep.

* `convergence[["pct"]]` was previously deprecated and discarded. It is now
  the primary convergence threshold — any code that suppressed the deprecation
  warning will observe the pct value being honoured.

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

* Pluggable convergence criteria via `convergence = list(criterion = "kl")`:
  accepts `"pct"` (default), `"max_err"`, `"mean_err"`, `"kl"`, `"chi2"`.
  `"stop_when"` controls whether any or all criteria must fire.

* Five quality metrics always returned in `attr(result, "result")`:
  `max_error`, `mean_error`, `kl`, `chi2`, `pct_change`.

* SOR adaptive under-relaxation for iEPPA via `sor` argument (default:
  auto-enabled). Diagnostics in `attr(result, "result")$sor`.

* Best-iterate tracking: `attr(result, "result")$best_weights` populated
  with weights at minimum observed marginal error.

* Python `harvest()` updated with matching `convergence` and `sor` arguments;
  new result fields exposed via `df.attrs["result"]`.

* Cell-compressed computation: faithful iEPPA operates at cell-level
  (unique (g_1,...,g_K) tuples) rather than observation-level, yielding
  up to 1000× speedup on surveys with low tuple diversity.
