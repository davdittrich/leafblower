# Constraints (SPEC-derived)

Extracted from the 32 `SPEC`-classified documents (2026-04-18 … 2026-08-14). Entries are
grouped by subsystem for navigation only; every entry is self-contained and carries its
own `source:`. Where a later spec supersedes an earlier one on the same subject, the
supersession is stated inside `content:` — the superseded entry is retained (specs
append, never rewrite history) and marked.

---

## Weight finalization order: normalize then bounds
- source: docs/superpowers/specs/2026-04-24-ieppa-speed-convergence-bounds-design.md
- type: protocol
- content: `bounds_mode` enum `RK_BOUNDS_CELL = 0` (default) / `RK_BOUNDS_UNIT = 1`
  appended to `rk_params_t`. Cell mode enforces the cell-aggregate contract
  `sum_{i in c} w_i in [L*|c|, U*|c|]` and only *counts* per-obs violations into
  `n_bounds_violated`; unit mode enforces strict per-observation
  `w_i in [min_weight, max_weight]` via intra-cell water-fill. Water-fill inner loop is
  bounded at 50 iterations; on exhaustion residuals are clamped to the nearest bound and
  `n_bounds_clamped` is incremented. The two counters are deliberately distinct fields
  (§6.2): `n_bounds_violated` = diagnosis (cell mode, no action), `n_bounds_clamped` =
  action taken (unit mode). Single-observation cells cannot redistribute — clamp + warn.

## Unit-mode water-fill invariant
- source: docs/superpowers/specs/2026-04-24-ieppa-speed-convergence-bounds-design.md
- type: nfr
- content: After water-filling, `sum_{i in c} w_i = X[c]` is preserved (redistribution
  conserves mass) and `w_i in [min_weight, max_weight]` for all i, up to the documented
  clamp failure in pathological cells. Merge gate: `n_bounds_clamped == 0` on benign
  unit-mode input (uniform-d, dense cells); `n_bounds_clamped < 0.001 * n` on skewed-d
  stress input. Default `bounds_mode="cell"` must preserve pre-existing weight vectors to
  1e-8 (not 1e-12 — that tolerance was declared untestable against
  `test-ieppa-faithful.R:109-136`).

## Cell-aggregate expansion formula (non-uniform design weights)
- source: docs/superpowers/specs/2026-04-23-ieppa-faithful-design.md
- type: protocol
- content: Obs `i` in cell `c` receives `w[i] = input_weight[i] * X[c] / X_init[c]`, and
  `w[i] = 0` when `X_init[c] == 0`. The earlier `w[i] = X[c] / n_per_cell[c]` form is
  WRONG whenever `d_i` varies within a cell (fixed per architect review). Re-asserted at
  `2026-04-25-calibration-solvers-design.md` §1 as the shared cell-table exit expansion
  `w_i = d_i * M[cell_of[i]]`, then `w_i *= n / sum(w)`.

## Cell table construction: sort-based, no hash map
- source: docs/superpowers/specs/2026-04-23-ieppa-faithful-design.md
- type: schema
- content: `CellTable { int M_cell; vector<int> cell_of; vector<int> n_per_cell;
  vector<vector<int>> g_per_cell; double W_input; }`. Construction is sort-based dedup
  (`std::sort` on a packed int64 key when `K <= 8` and all `cat_counts[k] <= 256`, else
  fixed-stride int32 tuples) — explicitly NOT `std::unordered_map`, to remove the
  collision-DoS surface on adversarial `group_ids`. Complexity O(n·K·log n) once.
  Deterministic: same input yields identical output. NA (`group_ids[k][i] = -1`) is a
  distinct category per margin, encoded as `cat_counts[k]` (one past valid range).

## Input validation floor
- source: docs/superpowers/specs/2026-04-23-ieppa-faithful-design.md
- type: api-contract
- content: `RK_ERR_BADARG` for: NULL/bad dims/NaN/Inf, targets not summing to 1,
  `min_weight >= max_weight`, `K > 64` (K_MAX — prevents unbounded cell-key allocation),
  `sum(weights) <= 1e-15`, `cat_counts[k] <= 0`, and
  `M_cell * K * sizeof(double) > SIZE_MAX/2`. Empty cell with a positive target is
  `RK_ERR_INFEAS`, not BADARG.

## Cell-table capacity auto-penalty field
- source: docs/superpowers/specs/2026-04-29-ieppa-alm-soft-capacity.md
- type: schema
- content: `CellTable` gains `double capacity_mu_auto`, computed in `build_cell_table`
  as `M_cell / n` when `n > 0 && M_cell > 0`, else 1.0. This is the auto value for the
  `oris_soft` ALM penalty and is chosen so the penalty Hessian matches the KL Hessian
  (`rho ≈ 1`, balanced blend) — justified against Boyd et al. 2011 §3.4.

---

## Convergence: metric and rule are orthogonal
- source: docs/superpowers/specs/2026-04-25-convergence-redesign.md
- type: api-contract
- content: The single `criterion` enum is REMOVED and replaced by two independent int
  fields, `metric` (CalibMetric 0=max_err 1=mean_err 2=kl 3=chi2 4=grake_norm
  5=l1_weight) and `rule` (CalibRule 0=threshold 1=improvement 2=plateau). Validation
  guards: `metric in [0,5]`, `rule in [0,2]`. `tol` must be in `(0,1)` for the `plateau`
  rule. Legacy mapping helper `rk_params_set_legacy_convergence(p, old_criterion)` maps
  0=PCT→(5,2), 1=MAX_ERR→(0,0), 2=MEAN_ERR→(1,0), 3=KL→(2,0), 4=CHI2→(3,0); values ≥5
  return `RK_ERR_BADARG`. SUPERSEDES the five-value `criterion` enum of
  `2026-04-24-convergence-metrics-sor-design.md` §1.

## Default convergence criterion
- source: docs/superpowers/specs/2026-04-25-convergence-redesign.md
- type: api-contract
- content: Default is `metric="max_err", rule="improvement", tol=0.001`. Shorthand keys:
  `list(improvement=X)` → max_err+improvement; `list(absolute=X)` → max_err+threshold;
  `list(pct=X)` → l1_weight+plateau (REDEFINED — `pct` previously meant max per-obs
  relative weight change). SUPERSEDES the `pct = 0.001` default of
  `2026-04-24-convergence-metrics-sor-design.md`, which itself superseded
  `absolute = 1e-6`. `stop_when` ("any"/"all") is retained only for combining the primary
  criterion with a secondary `max_err+threshold+absolute`; `stop_when` without `absolute`
  is an error.

## Metric formulas (all six, computed every check interval)
- source: docs/superpowers/specs/2026-04-25-convergence-redesign.md
- type: protocol
- content: `max_err = max_k max_j |Ŝ_kj/W − T_kj|`;
  `mean_err = mean_k(max_j |Ŝ_kj/W − T_kj|)`;
  `kl = max_k Σ_j T_kj log((T_kj+ε)/(Ŝ_kj/W+ε))` with ε=1e-10, skipping `T_kj = 0`;
  `chi2 = Σ_k Σ_j (Ŝ_kj − T_kj·W)²/(T_kj·W + 1)` (denominator floor 1; may return Inf on
  degenerate cells — documented, non-crashing);
  `grake_norm = max_k |misfit_k|/(1 + |pop_k|)` (survey::grake-compatible);
  `l1_weight = Σ_i |w_i[t] − w_i[t−1]| / n` (cell-indexed and normalised by `ct.W_input`
  for cell-level solvers, obs-indexed and normalised by `n` otherwise).

## Stopping-rule semantics
- source: docs/superpowers/specs/2026-04-25-convergence-redesign.md
- type: protocol
- content: `threshold`: `metric_t < tol`. `improvement`:
  `|metric_t − metric_{t−1}| / metric_{t−1} < tol`, with `prev_metric` initialised to +∞,
  reset at every homotopy level boundary, and a `prev > 1e-15` guard against explosion
  near zero. `plateau`: `converged = !(metric_t < prev_metric*(1−tol))`; on the first
  check `prev = ∞` makes `improved_enough` true, so the solver can never converge on the
  first check.

## Result field rename: pct_change → l1_weight_change
- source: docs/superpowers/specs/2026-04-25-convergence-redesign.md
- type: api-contract
- content: `pct_change` is REMOVED from `rk_result_t` and from the R result list; the
  replacement is `l1_weight_change`. Clean break, no alias, no deprecation — declared
  acceptable because the package was unreleased at spec time. `rk_result_t` also gains
  `convergence_metric`, `convergence_rule`, `convergence_tol`, `convergence_iter`
  (−1 when max_iter was hit), surfaced in R as `result$convergence_used`.

## Solver objective is decoupled from the stopping metric
- source: docs/superpowers/specs/2026-04-27-sinkhorn-a1-fix-ieppa-admm-method.md
- type: api-contract
- content: `convergence_objective` is renamed `convergence_solver_objective` (rename
  forces compile-time discovery of callers) and now reports the solver's MATHEMATICAL
  objective, not the stopping-criterion value. Mapping via
  `select_solver_objective(alg_id, m)` in `calib_dispatch.hpp`: ORIS/RAKING/SINKHORN →
  `m.kl` (weight KL); GREG → `m.chi2`; GRAKE → `m.grake_norm`; CHEBYSHEV → `m.errRp`;
  default → `m.errRp`. The three assignments `best_metric_seen`, `best_iter_val`,
  `best_objective_seen` MUST stay co-located in one named block. Non-finite objective is
  coerced to 0.0. SUPERSEDED IN PART by `docs/superpowers/specs/2026-08-14-removed-solver-slots-supersession.md`:
  the `case RK_ALG_GRAKE: return m.grake_norm;` arm is WITHDRAWN (slot 7 removed
  pre-release, commit `9a67891`). The `grake_norm` METRIC is unaffected and remains live
  (`src/greg.cpp:153`, `src/logit_calib.cpp:558`, packed at `src/c_api.cpp:74/:110/:426`) —
  do not confuse the live metric with the removed solver. Every other arm stands.

## Best-iterate tracking is always on
- source: docs/superpowers/specs/2026-04-24-convergence-metrics-sor-design.md
- type: protocol
- content: Every solver keeps `best_errRp`, `best_iter`, `W_best` (solver-local vectors,
  NOT `CalibState`, to avoid ownership ambiguity across the solver boundary; sized
  `M_cell` for cell-level solvers, `n` for obs-level). Best-iterate tracks `errRp`
  regardless of the active convergence criterion — so `best_error <= max_error` always
  holds, both being errRp. `W_best` expansion applies scalar expansion plus
  sum-normalisation ONLY: **no water-fill, no bounds clamping** — `W_best` is a mid-loop
  snapshot and re-clamping changes its meaning.

## Best-iterate metric must be the configured metric, not the errRp proxy
- source: docs/methods/greenkhorn.md
- type: protocol
- content: Two best-iterate trackers are kept deliberately distinct: `sraa_best` on the
  cheap `errRp` scale (drives the SRAA outer-stall revert) and the reported `best` on the
  configured metric. Conflating them is the regression re-introduced twice per
  `CLAUDE.md`; the correct call is `select_metric(sraa_cfg.metric, cm)` evaluated at
  `kErrCheckInterval` boundaries.

## Status codes: budget and stall are distinct from NOCONV
- source: docs/superpowers/specs/2026-04-28-convergence-status-design.md
- type: api-contract
- content: Appended (no renumbering of 0–3): `RK_ERR_BUDGET = 4` (budget exhausted while
  loss still decreasing — increase `max_iterations`) and `RK_ERR_STALL = 5` (loss
  plateau at the constrained optimum — weights are valid, accept them).
  `RK_ERR_NOCONV = 1` is retained for ABI compatibility but is NOT emitted by updated
  solvers. Priority: criterion (0) > stall (breaks inside the loop) > budget (set after
  the loop only if status is still the initial NOCONV). Breaking change for callers
  testing `status == 1`; callers testing `status > 0` are unaffected.

## Stall criterion: weight KL, not errRp (flat loop)
- source: docs/superpowers/specs/2026-04-28-convergence-status-design.md
- type: protocol
- content: The flat raking loop replaces the `min_errRp_window` stall detector with a
  weight-KL stall, because errRp is not monotone near the constrained KL optimum and
  fires false stalls while weight KL is still strictly decreasing. Mandatory ordering:
  the convergence criterion is evaluated BEFORE the stall detector on every check
  interval, and `wkl <= tol_abs` exits as `RK_OK` — so perfect calibration cannot be
  misreported as a stall. SQUAREM/accelerated paths keep an errRp-based stall (CBB
  extrapolation is not KL-monotone); `convergence_reason` distinguishes `"stall_kl"`
  from `"stall_errRp"`.

## convergence_reason string contract
- source: docs/superpowers/specs/2026-04-28-convergence-status-design.md
- type: api-contract
- content: `result$convergence_used$convergence_reason` is a non-NA character for every
  exit path: `"criterion"` (0), `"budget"` (4), `"stall_kl"` (5, flat), `"stall_errRp"`
  (5, accelerated), `"infeasible"` (2), `"error"` (3), `"legacy"` (1). R behaviour:
  status 2 and 3 `stop()`; status 4 and 5 `warning()` and return valid weights.

## Weight-change stall for the accelerated path
- source: docs/superpowers/specs/2026-04-28-squarem-geometry-fix.md
- type: protocol
- content: SUPERSEDED IN PART by the SRAA migration. For the accelerated branch the
  errRp stall is replaced by an obs-level L1 weight change
  `wchange = Σ_c |X[c] − X_prev[c]| / n_per_cell[c] / n`, fed through the same
  relative-improvement window as the flat KL stall (`min_loss_window` initialised +∞,
  `kMaxNoImprove` consecutive non-improvements → `RK_ERR_STALL`). All changes are
  confined to the accelerated branch; the flat loop's `compute_weight_kl()` stall is
  explicitly unchanged.

---

## ORIS core update rule
- source: docs/methods/oris.md
- type: protocol
- content: Per (margin k, category j) the naive Sinkhorn ratio is
  `naive = T_kj·W_input / S_kj` with `S_kj = Σ_{c∈(k,j)} X_cur[c] / f_old`. The applied
  step is `f_new = f_old^(1−α·ω) · naive^(α·ω)` (linear) or
  `lf_new = (1−α·ω)·lf_old + α·ω·(log T_kj·W − log S_kj)` (log). `ω` is the SOR factor,
  `α = 1/(1 + β·stress)` the infeasibility damping with `β = kAlphaBeta = 0.5`. The NET
  exponent `α·ω` must stay inside `(0, 2)` — the window with a global-convergence proof
  (Thibault et al. 2021). `ω = 1, α = 1` is exactly Sinkhorn–Knopp. Margins are swept
  Gauss–Seidel.

## Linear/log space dispatch threshold
- source: docs/superpowers/specs/2026-04-23-ieppa-convergence-hardening-design.md
- type: protocol
- content: `kLinearSpaceThreshold = 2.0`: run the linear-space inner loop when
  `M_cell/n > 0.5` (compression ≤ 2×), otherwise stay in log space. Below 2× compression
  log-space transcendental overhead dominates; above it, per-iteration cost is dominated
  by the O(K) inner loop instead. The linear inner loop MUST use prefactored products —
  maintain `X_cur[c] = X_init[c]·W[c]·Π_m f_lin[m][g_m(c)]`, compute
  `S_kj = Σ_{c∈bucket} X_cur[c] / f_lin[k][j]` (one division per cell), then rescale only
  that bucket by `new_f/old_f`. This is O(K·M_cell) per iter. The naive per-cell K-way
  product is O(K²·M_cell) and is explicitly marked DO NOT IMPLEMENT (measured 14.45× on
  kk1204, failing the ≤2× gate at K=20).

## Linear-space overflow detection and correction
- source: docs/superpowers/specs/2026-04-27-ieppa-linear-overflow-fix.md
- type: protocol
- content: Overflow is a PRODUCT phenomenon, not a per-factor one — individual `f_lin`
  values sit at 5–19 while `Π_k f_lin` reaches ~2.5e15. The per-margin geometric-mean
  trigger (T1.A) cannot detect it and MUST BE REMOVED BEFORE the replacement is inserted
  (the two mechanisms are not composable — T1.A scales per margin, T1.B applies a global
  scale). Replacement (T1.B): maintain `cell_lf[c] = Σ_k lf[k][g_k(c)]` and a high-water
  mark `cell_lf_hwm = max_c(log_X_init[c] + cell_lf[c])` updated only on positive deltas;
  when `cell_lf_hwm >= log(kLinearOverflowThreshold)` apply `shift`, distributing
  `−shift/K` across all K margins (INCLUDING each margin's NA bucket, else the invariant
  breaks for NA cells) and scaling `X_cur *= exp(−shift)`. `W` is intentionally NOT
  updated — the invariant `X_cur = X_init·W·Π f_lin` is preserved by the paired update.

## Linear→log fallback is one-shot and state-clean
- source: docs/superpowers/specs/2026-04-23-ieppa-convergence-hardening-design.md
- type: protocol
- content: On overflow trip the solver aborts the linear path, resets EVERY piece of
  linear-path state — `lf[]` zeroed, `f_lin[]` to 1.0, `X_cur[]` to 0.0, `W[c]` to 1.0,
  `X_tilde[c]`/`X[c]` to `X_init[c]`, `infeas_streak[]` to 0,
  `persistent_infeas_pairs` cleared, iteration counter to 0, `alpha` to 1.0 — and
  restarts in log space. A `linear_fallback_used` boolean prevents re-entry, guaranteeing
  termination. `2026-04-29-ieppa-alm-soft-capacity.md` adds: `lambda_cell[]` is also reset
  on fallback (the linear-scale dual is invalid in log-path X̃ and would produce silent
  NaN or wrong-magnitude penalty).

## Persistent-infeasibility tracker
- source: docs/superpowers/specs/2026-04-23-ieppa-convergence-hardening-design.md
- type: protocol
- content: The single latched `is_infeasible` bool is replaced by per-(k,j) persistence:
  `kInfeasPersistence = 5` consecutive empty observations promote a pair into
  `persistent_infeas_pairs`; insertion uses `==` not `>=` so it fires exactly once per
  pair. `record_nonempty` MUST reset the streak AND erase the pair from the set —
  otherwise a recovered bucket stays latched and a max-iter exit falsely reports
  `RK_ERR_INFEAS`, defeating the whole work unit. Categories with `target <= 0` are
  skipped entirely. Known accepted gap: a bucket oscillating faster than the persistence
  window is never flagged; the solver returns NOCONV with high errRp instead.

## Adaptive damping is unlatched
- source: docs/superpowers/specs/2026-04-24-ieppa-speed-convergence-bounds-design.md
- type: protocol
- content: `alpha = 1.0 / (1.0 + beta * stress)` with `beta = 0.5` and
  `stress = min(max_{k,j} infeas_streak, kInfeasPersistence)` over buckets with positive
  target. `alpha` is recomputed each sweep; `damped_latched` is DELETED. The `min()` cap
  is load-bearing: uncapped, stepstone drove stress to 490+, `alpha ≈ 0.004`, and errRp
  regressed 2.21e-3 → 2.28e-3; with the cap it improves to 2.14e-3. Fast path
  `if (stress == 0)` must branch — `std::pow(x, 1.0)` is not free. SUPERSEDES the hard
  `alpha = 0.5` latch of `2026-04-23-ieppa-convergence-hardening-design.md` §6.

## Fused post-sweep block (2 passes → 0)
- source: docs/superpowers/specs/2026-04-24-ieppa-speed-convergence-bounds-design.md
- type: nfr
- content: The separate X̃ compute and X_cur rebuild passes are eliminated; capacity is
  computed inline from the sweep's final `X_cur` via `X_tilde_c = X_cur[c] / W[c]`.
  Degenerate cells (`X_init <= 0`, `W <= 0`, or `X_tilde_c <= 0`) reset `W[c] = 1.0`
  explicitly to prevent stale-W propagation. Framing correction: the sweep itself is NOT
  fused and remains the dominant pass — the claim is "2 post-sweep passes → 0", not
  "3 passes → 1". Target: kk1204 per-iter wall-clock ratio ≤ 1.5× raking (from 2.17×);
  P1.1 alone yields only ~5%, the rest must come from iteration reduction.

## SOR over-relaxation window and mode table
- source: docs/methods/oris.md
- type: api-contract
- content: Over-relaxation is OPT-IN via `sor = list(auto = TRUE)`; a bare
  `harvest(method="oris")` runs plain IPF at `ω = 1`. Three strategies selected by
  `omega_mode_id`: `0` heuristic per-margin nudge (×1.05 on errRp decrease, ×0.7 on
  sign-flip, ceiling `omega_max` 1.5); `1` fixed per-margin jump to `omega_max` while
  errRp decreases; `2` iterate-change single global ω with ceiling
  `kSorProdCeiling = 1.8` — documented as the SHIPPED DEFAULT. All modes observe a
  burn-in (default 20) and damp toward `omega_min` (default 0.3) on oscillation. NOTE:
  this reverses the direction of `2026-04-24-convergence-metrics-sor-design.md` §3, which
  specified SOR as *under*-relaxation (`omega_init = 1.0`, adapt downward only,
  `omega_min = 0.3`) to damp oscillation.

## Iterate-change θ₂ estimator (mode 2)
- source: docs/superpowers/specs/2026-06-02-oris-iterate-change-omega-design.md
- type: protocol
- content: The observable is the FREE-COORDINATE ITERATE CHANGE, not the marginal
  residual: `S_dX = Σ_{c: !is_pinned[c]} (X[c] − X_snapshot[c])²`, measured every
  `I = kErrCheckInterval` sweeps, `ratio = S_dX(m)/S_dX(m−1) → ρ^(2I)`, recovered by the
  block root `theta2 = ratio^(1/I)`, then `ω = 2/(1 + √(1 − clamp(theta2, 0, 1−1e-9)))`
  capped at 1.8. ONE global ω is written to every margin. Root cause of the predecessor's
  failure, verified three ways: `R2_free` with the clamped mass folded into the free
  target is algebraically IDENTICAL to the global marginal residual, which plateaus at a
  nonzero floor on an infeasible-after-clamp margin (ratio → 1 → ω → 2 → stall), while
  `‖ΔX_free‖²` decays at the true free-block rate and is feasibility-agnostic. Also fixes
  a latent cadence bug: treating a lag-`I` ratio as ρ² systematically under-reads ρ.
  SUPERSEDES the marginal-residual estimator of
  `2026-05-31-oris-free-subspace-omega-design.md`.

## ω guard envelope (10 ordered gates)
- source: docs/superpowers/specs/2026-06-02-oris-iterate-change-omega-design.md
- type: protocol
- content: In order, before writing ω: (1) `alm_active` → block not reached (oris_soft
  out of scope); (2) `n_free == 0` → ω=1, reset EMA; (3) warm-up / `S_dX_prev` unset →
  ω=1; (4) `!isfinite(S_dX)` or `S_dX_prev < kResidFloor (1e-12)` or free mass
  `< kMinSafeTotalWeight (1e-100)` → ω=1, reset EMA; (6) `ratio >= 1` → HARD reset
  `theta2_ema = 0`, ω=1, `consec_up++` (a growing increment can never produce elevated
  ω, and no stale EMA survives); (7) EMA with `kSorEmaAlpha = 0.2`; (8) formula, ceiling
  `kSorProdCeiling = 1.8` — explicitly NOT `kSorSpectralCeiling = 1.99`, which IS the
  known oscillation value; (9) oscillation damp ×0.7 on trend sign-flip; (10)
  `consec_up >= 3` → ω=1 plus `kSorCooldown = 5`, and `cooldown_trips >= 3` →
  `sor_latched = true`, permanent ω=1 for the solve. The latch makes a limit cycle of any
  period structurally impossible rather than merely asserted away.

## Active-set pinning definition
- source: docs/superpowers/specs/2026-05-31-oris-free-subspace-omega-design.md
- type: protocol
- content: `is_pinned[c] = (X[c] >= U_cell[c]*(1−kPinTol)) || (X[c] <= L_cell[c]*(1+kPinTol))`
  with `kPinTol = 1e-9` relative, written in the flat capacity-clamp loop and READ in the
  convergence sweep — two different loops separated by the ALM/metrics blocks; no fused
  single loop may be assumed. Solve-local `std::vector<char>(M_cell)`. The theoretical
  basis: with per-cell clamping the active set is a hard wall, error propagates only
  through the free block, and the active-set-dependent optimum is
  `ω_opt(I) = 2/(1+√(1−ρ(M_II)²))` on the principal submatrix `M_II`. No closed-form
  box-constrained `ω_opt` exists in the literature (Lehmann 2022 and Soma–Uschmajew 2024
  do not treat inequality constraints).

## SOR observability fields
- source: docs/superpowers/specs/2026-05-31-oris-free-subspace-omega-design.md
- type: api-contract
- content: OUTPUT-ABI additive: `sor_omega_mean` plus one counter per distinct fallback
  cause so every silent all-fallback mode is individually diagnosable —
  `sor_n_pinned_fallback`, `sor_n_warmup_fallback`, `sor_n_converged_fallback`,
  `sor_n_resid_grew` (the mj1p.2 canary), `sor_n_monotone_cooldown`, plus a latch flag.
  Requires updating `EXPECTED_RK_RESULT_BYTES` and wiring through both `r_bridge.cpp` and
  `c_api.cpp`/pybind11 in the same commit. No new INPUT ABI —
  `sor_omega_mode_id` already exists in `rk_params_t`.

## omega_mode_id must be selectable and defaulted identically in R, C and Python
- source: docs/superpowers/specs/2026-06-02-oris-iterate-change-omega-design.md
- type: api-contract
- content: Documented live divergence: the C default is 2 (`src/types.hpp:74`,
  `src/c_api.cpp:114`), R `parse_sor` defaults to 1L, and the Python `_parse_sor` returns
  a 6-tuple that DROPS `omega_mode_id` entirely, so `_bindings.cpp` never forwards it —
  a Python call silently runs the C default while the identical R call runs mode 1. Two
  required corrections: (a) `_parse_sor` 6-tuple → 7-tuple and `_bindings.cpp` forwards
  `p.sor_omega_mode_id`; (b) the ship decision flips the default in ALL THREE sites
  together (`R/harvest.R parse_sor`, `src/types.hpp:74`, `src/c_api.cpp:114`) — SHIP → 2
  everywhere, NO-SHIP → 1 everywhere, no site left behind.

## Jacobi log sweep (optional, benchmark-gated)
- source: docs/superpowers/specs/2026-05-01-jacobi-log-sweep-design.md
- type: protocol
- content: `CalibState.jacobi_log` (R: `jacobi_sweep = FALSE`) selects Jacobi semantics
  for the LOG path only — freeze `cell_lf` into `cell_lf_frozen` at iteration start,
  sweep all K margins against the snapshot with no mid-sweep patches, then rebuild
  `cell_lf` once in a sequential O(K·M_cell) pass. The linear path is out of scope (it
  already has O(1) leave-one-out via `X_cur[c]/f_lin[k][j]`). Jacobi is a DIFFERENT
  fixed-point iteration and may land on a slightly different fixed point; acceptance is
  `max_error` within 1e-4 of Gauss–Seidel on the same fixture. Ship-as-default only if
  `median(wall_ratio) < 1.0` across the 18-cell grid; otherwise flag-off or removed.

---

## SRAA-m replaces SQUAREM
- source: docs/superpowers/specs/2026-04-29-i0am-sraa-acceleration.md
- type: protocol
- content: SQUAREM is removed from greenkhorn and raking and replaced by Safeguarded
  Regularized Anderson Acceleration. Root cause: greenkhorn SQUAREM had no backtracking
  and accepted anything within `err_star <= err_w2 * 1.01`, accumulating degradation over
  540 super-steps (max_err 2.12e-3 vs plain 1.57e-3, +35% WORSE); raking's step-halving
  only partly mitigated it (+9%). The CBB scalar α is calibrated for smooth operators and
  overestimates step length at the bounds kink. Cost: AA path is exactly 2 F-evals, plain
  path exactly 1 (SQUAREM needed 3). Constants: `kSRAAm = 5`, `kSRAAMaxM = 10`,
  `kSRAAMinCount = 2`, `kSRAAdeltaReg = 1e-10`, `kSRAARestartGamma = 2.0`. Memory at
  stepstone ≈ 176 MB (+63 MB vs SQUAREM). BREAKING: existing `accelerate=TRUE` weight
  vectors change.

## SRAA state and step invariants
- source: docs/superpowers/specs/2026-04-29-i0am-sraa-acceleration.md
- type: schema
- content: `SRAAState` is a pre-allocated slot-contiguous circular buffer
  (`buf[slot*M + cell]`, so dot products iterate sequentially and vectorise) with zero
  per-step heap allocation; the m×m Gram/rhs/gamma arrays are stack-resident. Invariants:
  ACCEPT_PLAIN always uses `std::swap`, never `operator=` (O(1), not O(M)); `f_eval` is
  NEVER called with the outer `X` — always `state.F_cur` or `state.scratch`, which must be
  distinct objects; the restart check is guarded by `has_prev` so it can never fire on the
  first call; history is cleared on restart, safeguard rejection, LDLT failure and
  non-finite `err_AA`.

## SRAA global safeguard and revert-to-best
- source: docs/superpowers/specs/2026-04-29-i0am-sraa-global-safeguard.md
- type: protocol
- content: The per-step local safeguard `err_AA <= err_plain` is replaced by a GLOBAL
  quality floor `err_AA <= best_err_seen * (1 + kSRAAGlobalEps)` with
  `kSRAAGlobalEps = 1e-3`. `best_err_seen` and `X_best` are monotone across the whole
  solve and are deliberately NOT reset by `clear()` — after any restart the plain steps
  must beat the floor again before AA can re-fire. `stall_count >= kSRAAStallWindow (15)`
  triggers `X = X_best` by COPY (not swap, so `X_best` stays valid across repeated
  reverts) plus `clear()`. Rationale recorded: the sibling project's Aitken guard
  over-rejected valid steps; this guard was under-rejecting invalid ones — both from
  using a purely LOCAL criterion.

## SRAA F_eval must use adaptive sort
- source: docs/superpowers/specs/2026-04-29-i0am-sraa-correct-all-scales.md
- type: protocol
- content: `f_eval_sraa` in greenkhorn must re-sort `errRp` after EACH greedy step
  (identical to plain greenkhorn), not sort once at round entry. Fixed-sort and
  adaptive-sort operators have DIFFERENT fixed points — AA extrapolation converges into
  the KL-optimal basin (2.12e-3) instead of the max_err-optimal basin (1.57e-3), and no
  per-step safeguard can distinguish the two basins in early iterations. Raking's
  `F_eval` already used adaptive order and is unchanged. Paired with an outer revert:
  quality above `best * (1 + kSRAAOuterSlack = 0.10)` for `kSRAAOuterStallWindow = 5`
  consecutive outer iterations → revert to `X_best`, clear AA history, rebuild `W`,
  `S_flat` and `errRp` from the reverted `X`.

## SRAA on ORIS operates in log-factor space
- source: docs/superpowers/specs/2026-04-30-ieppa-sraa-acceleration.md
- type: protocol
- content: The SRAA iterate for ORIS/oris_soft is `lf` (dimension `n_cats_total_with_na`
  ≈ 50–500), NOT `X_cur` (dimension `M_cell` ≈ 1e6) — history is ~13 KB instead of ~64 MB
  at stepstone, and `lf` uniquely parameterises the active state. `sraa_step` gains
  `bool apply_clamp = true`; when false, `L_cell`/`U_cell` are never dereferenced so
  empty vectors are safe. `pack_lf`/`unpack_lf` rebuild `f_lin`, `cell_lf` and `X_cur`
  from the flat iterate. On F_eval overflow the partial `lf` is packed and `+inf` returned
  so the safeguard rejects and reverts. `lf_best` is a SEPARATE saved copy — it cannot be
  derived from `W_best`, since `lf` is not recoverable from `X_cur` without solving a
  system.

## Accelerate implies SOR off and greedy downgraded (ORIS)
- source: docs/superpowers/specs/2026-04-30-ieppa-sraa-acceleration.md
- type: api-contract
- content: With `accelerate=TRUE` on ORIS/oris_soft: `sor_auto = st.sor_cfg.auto_adapt && !st.accelerate`
  so ω stays at `omega_init`; a `greedy` scheduler is silently downgraded to round-robin
  with a `verbose>=1` log (AA requires a deterministic sweep map);
  `res$iterations` reports F-EVALS consumed (AA step = 2, plain = 1), not outer BCD
  iterations; `res$aa_accepted_count` reports accepted AA steps; SRAA history is reset at
  every homotopy level boundary and on linear→log fallback. Confirmed by
  `docs/methods/00-overview.md`: enabling Anderson disables SOR — the two accelerators
  target the same fixed-point map and do not compound.

## Only stateless solvers may be SRAA-wrapped
- source: docs/methods/00-overview.md
- type: protocol
- content: SRAA-m applies only where the sweep is a pure fixed-point map with no carried
  correction vectors: ORIS, Raking, Greenkhorn. Sinkhorn (Dykstra correction state) and
  the Newton/IPM family (iterate-dependent factorisations) cannot be wrapped.
  Over-relaxation is likewise ORIS-only: Sinkhorn's stateful Dykstra+bisection breaks the
  power-law ≡ log-SOR equivalence (the rate estimator misreads correction grind as
  slowness and over-shoots), and Greenkhorn's greedy single-coordinate selection falls
  outside the proven full-sweep regime with no published theorem and no measured speed-up.

---

## Raking uses Bregman (multiplicative) Dykstra, not Euclidean
- source: docs/superpowers/specs/2026-04-27-raking-bregman-dykstra-design.md
- type: protocol
- content: Box correction `p[c]` is initialised to 1.0 and updated multiplicatively:
  `yc = X[c]*p[c]; Xc = clamp(yc, L_cell[c], U_cell[c]); p[c] = (Xc > 0) ? yc/Xc : 1.0`
  (the guard prevents a silent 0/0 NaN when `min_weight = 0` and a cell is a structural
  zero). Hyperplane correction `q_hyp` is initialised to 1.0 and becomes a scale:
  `X[c] *= q_hyp; scale = n/s; X[c] *= scale; q_hyp = (scale > 0) ? 1/scale : 1.0`. The
  same multiplicative pattern applies in the post-loop finaliser. Rationale: the additive
  hyperplane correction normalises by SUBTRACTION, which is not the KL projection onto
  `{sum = n}` — mixing Euclidean Dykstra with multiplicative IPF leaves the algorithm with
  no unified convergence theorem (Bauschke & Lewis 2000 supplies one for the Bregman form).

## SOR and greedy scheduler wired into raking
- source: docs/superpowers/specs/2026-04-27-raking-bregman-dykstra-design.md
- type: api-contract
- content: `CalibState::sor_cfg` is wired into raking's IPF marginal step as
  `X[c] *= pow(T/S, eff_omega)` with `eff_omega <= 1`, active only when bounds are
  present (`st.min_weight > 0.0 || hi < 1e300`) — without bounds no oscillation occurs and
  `eff_omega = 1.0`. Guard: `if (S <= 0.0) continue;` before the `pow` for structurally
  empty groups. `CalibState::scheduler` sorts the K-margin sweep by descending per-margin
  `errRp_k`, computed DURING the sweep from the `bucket[j]` sums already in hand (no extra
  cell pass). The `@param sor` docstring must stop saying "iEPPA only".

## Raking per-obs bounds via post-exit clamp
- source: docs/superpowers/specs/2026-04-25-calibration-solvers-design.md
- type: protocol
- content: After the `w_i = d_i × M[c]` expansion, a per-obs clamp
  `w_i = clamp(d_i × M[c], min_weight, max_weight)` guarantees strict per-obs bounds. The
  documented trade-off: for non-uniform `d_i` the clamp slightly distorts margins in cells
  where some obs were clamped, bounded by within-cell `d_i` variance × |M[c] − clamp|, and
  reported in `max_error`.

## SQUAREM SqS3 (historical, raking)
- source: docs/superpowers/specs/2026-04-28-raking-squarem-design.md
- type: protocol
- content: SUPERSEDED by SRAA-m (`2026-04-29-i0am-sraa-acceleration.md`). Recorded for
  provenance: `F` = one full inner iteration (K-margin Bregman IPF sweep + box + hyperplane
  correction); CBB step `α = −‖r‖₂/‖v‖₂` capped to `[−1000, −1]`; the fixed-point guard is
  RELATIVE (`‖v‖₂/(‖w2‖₂ + 1e-300) < 1e-10`); step-halving uses `α ← (α−1)/2`, a
  contraction with fixed point −1, up to 16 times, comparing against the already-computed
  `‖r‖₂` so no 4th F-eval is needed; snapshots of `X`, `p`, `q_hyp` are taken AFTER the
  ‖v‖ guard and BEFORE extrapolation. `accelerate=TRUE` on a non-raking method warns and
  falls back. `rk_params_t` is NOT modified; only the `.Call` arity changes.

## SQUAREM α must be computed in observation geometry
- source: docs/superpowers/specs/2026-04-28-squarem-geometry-fix.md
- type: protocol
- content: Two separate norms are required. α uses OBS-level weighted norms
  `Σ_c r_c²/n_per_cell[c]` (because all obs in a cell share the multiplier, so
  `r_obs[i] = r_cell[c]/n_c`); step-halving uses CELL-level unweighted norms so both sides
  of the comparison live in the same space. Cell-level α is systematically too large in
  magnitude, over-extrapolates, triggers halving, collapses α to ≈ −1 and yields no
  acceleration — measured 113 F-evals for a WORSE answer (1.71e-3) than the flat loop's
  50 F-evals (1.60e-3). `n_per_cell[c] >= 1` is guaranteed by `build_cell_table`.

---

## Sinkhorn capacity projection is KL-Dykstra with bisection
- source: docs/superpowers/specs/2026-04-25-calibration-solvers-design.md
- type: protocol
- content: The Euclidean water-fill is replaced by a log-domain Dykstra projection:
  bisect `μ` such that `Σ_c clamp(X[c]·exp(a[c]+μ), L_c, U_c) = target_mass`, where
  `target_mass` is the mass BEFORE projection (mass-preserving), then accumulate
  `a[c] += log(X[c]) − log(X_proj[c])`. Bisection bounds
  `μ ∈ [log(L_min/X_max), log(U_max/X_min)]`, ~40 iterations to 1e-12, O(M_cell × 40) per
  capacity step. Rationale: in KL geometry the correct projection onto the box is
  multiplicative; Euclidean Dykstra on a KL problem breaks monotone descent and produces
  the two-accumulation-point oscillation (Gietl–Fröhlich 2013). Infeasibility pre-check:
  `Σ U_c < target_mass || Σ L_c > target_mass → RK_ERR_INFEAS`.

## Sinkhorn preserves Σ X = n by construction
- source: docs/methods/sinkhorn.md
- type: nfr
- content: The bisection targets `n`, so `Σ X[c] = n` holds every iteration and NO final
  renormalisation is applied (source comment). `f(μ)` is strictly increasing and the
  bracket is derived from per-cell clamp thresholds so `f(lo) <= 0 <= f(hi)` is
  guaranteed — bracketing cannot fail. Cells permanently at their lower bound are frozen
  (`at_lower`) and excluded from further Dykstra corrections to stop `a[c]` accumulating
  without bound. `convergence_rule = 0` — sinkhorn uses no improvement rule.

## Sinkhorn default stopping metric is kl
- source: docs/superpowers/specs/2026-04-27-sinkhorn-a1-fix-ieppa-admm-method.md
- type: api-contract
- content: `R/harvest.R` sets the sinkhorn default stopping metric to `"kl"` rather than
  `"max_err"`, because sinkhorn's KL is monotone-decreasing and `kl+improvement` is the
  criterion matching its mathematical objective.

---

## Shared normal-equations kernel
- source: docs/superpowers/specs/2026-04-25-calibration-solvers-design.md
- type: api-contract
- content: `src/calib_linalg.hpp` provides
  `compute_normal_equations(ct, D, N, cat_offset, n_cats_total)`,
  `ldlt_factor_inplace(N, n, eps_perturb)` and `ldlt_solve(L, d, rhs, n)`. LDLT with
  Gill–Murray diagonal perturbation (NOT plain Cholesky) is required to survive degenerate
  margin groups (zero-mass cells give a zero row/col in N). ALL dimension parameters are
  `size_t`, never `int`, to prevent int32 overflow above `n_cats_total = 46340`. Hard
  memory guard `kNCatsTotalMax = 2048` → `RK_ERR_BADARG` with the message pointing the
  user at `method='ieppa'` or `'raking'`. NOTE: the same spec's resolved-questions table
  states the cap as 8192, contradicting the 2048 in its own code block.

## Chebyshev reference-category elimination
- source: docs/superpowers/specs/2026-04-26-chebyshev-nu-reference-elimination.md
- type: protocol
- content: For every margin with `Σ_j T[k][j] = 1`, the normalisation row is a linear
  combination of the margin rows, so `schur_nu = D_ν − eᵀN⁻¹e ≈ 0` and the ν
  Sherman–Morrison degenerates — the Newton step stops being sum-preserving and W drifts
  to 0 after ~500 iterations on K≥4. Fix: drop the LAST category of every margin with
  `cat_counts[k] >= 2` from the Newton system, giving
  `nct_red = nct − count(k : cat_counts[k] >= 2)`. Margins with `cat_counts[k] < 2` are
  NOT reduced (dropping their only category would remove all Newton correction for that
  margin). All LP slacks and duals stay full-size; `dlambda_full[is_ref] = 0`; `dnu`
  applies to ALL cells in the `dX` reconstruction. Guard `schur_nu > 1e-8`. The renorm
  workaround branch is superseded and must be discarded.

## Chebyshev Mehrotra predictor-corrector
- source: docs/superpowers/specs/2026-04-29-chebyshev-greg-fix.md
- type: protocol
- content: The single-step barrier update with fixed `μ *= 0.1` is replaced by Mehrotra's
  two-phase step. Invariants: `N = A·D·Aᵀ` is rebuilt FRESH each iteration (so Jacobi
  scaling is never applied twice to the same matrix); `D` is held FIXED across Phase A and
  Phase B; Phase B REUSES Phase A's LDLT factorisation, only the RHS changes;
  `σ = clamp((μ_aff/μ)³, 1e-8, 1.0)`; `μ_aff` division is guarded and clamped to
  `[0, μ·100]`; `m == 0` (no complementarity pairs) skips Mehrotra entirely. Jacobi
  preconditioning `D_jac[j] = 1/sqrt(max(N[j][j], 1e-12))` is applied to N and rhs and
  unscaled after the solve. The second-order term is `−Δs_aff·Δy_aff`.

## Chebyshev corrector cross-term is correct as written
- source: docs/methods/chebyshev.md
- type: protocol
- content: `corr = −Δs_aff·Δy_aff` is the genuine Mehrotra cross-term and must NOT be
  "fixed" (also recorded in `CLAUDE.md`). A dual-explosion guard reverts to a scaled unit
  diagonal when `μ_new > 100·μ`. Hard cap `kMaxIpm = 500`; convergence at
  `kTolMu = 1e-6`. Best-iterate is tracked by calibration error `errRp`, NOT the LP
  objective `δ` — δ can bottom out while the primal is still improving.

## Chebyshev warm start from ORIS
- source: docs/superpowers/specs/2026-04-29-chebyshev-greg-fix.md
- type: protocol
- content: Before dispatching chebyshev, run a short pre-solve
  (`inner_max_iter = max(5, min(100, budget/10))`) and pass its obs-level best weights
  plus `delta_warm = ieppa_max_err * 1.5`. CRITICAL: the pre-solve MUTATES `st.weights`
  in place, so a deep copy must be handed to it and the temporary `CalibState` must not
  escape the block (dangling-pointer risk). Aggregation to cell masses happens INSIDE
  chebyshev after its own `build_cell_table`, avoiding a second cell-table build, and is
  followed by a clamp plus mass-preserving renormalisation. Empty or wrong-length warm
  weights fall back to cold start (uniform init, `delta_0 = 1.0`). The source spec writes
  "chebyshev/grake"; the grake half is WITHDRAWN per `docs/superpowers/specs/2026-08-14-removed-solver-slots-supersession.md`
  (slot 7 removed). The contract applies to chebyshev only.

## GREG quality warning
- source: docs/superpowers/specs/2026-04-29-chebyshev-greg-fix.md
- type: api-contract
- content: `harvest()` emits a `warning()` (never a fallback or a `stop()`) when
  `method == "greg"` returns a finite `max_error > 0.05`, naming the observed error, K and
  `max_weight`, and suggesting `method='raking'` or `'ieppa'`. Threshold rationale: 0.05 is
  ~50× typical raking quality (1e-3), so it flags only catastrophic failures. Root cause
  documented: for K≥5 overlapping margins with tight bounds the Gram inverse is
  near-singular and GREG's Newton converges to the wrong interior point.

## Logit link makes bounds structural
- source: docs/superpowers/specs/2026-04-29-greenkhorn-solver.md
- type: protocol
- content: `w_c = L_cell[c] + (U_cell[c] − L_cell[c])·σ(z_c)` with
  `z_c = Σ_k λ[cat_offset[k] + g_per_cell[k][c]]`, so `w_c ∈ [L,U]` for ALL λ — no
  clamping, no water-fill, no active set. Newton weight
  `D_eff[c] = (U−L)·σ(z)·(1−σ(z))` is exactly `∂w_c/∂z_c` (Deville–Särndal 1992 eq. 8),
  not an approximation. `z` is clamped to ±700 before `exp`. Convergence MUST route
  through the shared `check_convergence` — hand-rolling `max_resid < pct_tol * n`
  implements one rule and silently ignores the others. `RK_ERR_NOCONV = 1` is never
  returned; the loop exit classifies into `RK_ERR_BUDGET` (residual improved by >0.1%) or
  `RK_ERR_STALL`.

## Logit Newton stabilisation: three layers
- source: docs/superpowers/specs/2026-04-29-sv89-logit-newton-fix.md
- type: protocol
- content: (1) Armijo backtracking with `c = 0.01`, `kMaxHalvings = 10`, accepting when
  `‖b_trial‖² < ‖b‖²(1 − 0.01α)`, searching in `[0, alpha_max]` where `alpha_max` caps the
  per-cell shift at `kMaxDeltaZ = 2.0` (a huge step is damaging even at α = 2⁻¹⁰).
  (2) Design-weight initialisation: solve `(AAᵀ)λ₀ = A·z_target` with
  `z_target = logit(clip((X_init−L)/(U−L), 1e-4, 1−1e-4))`, placing λ₀ in the convergence
  basin. (3) `D_eff` floor `max(1e-6·range, range·σ(1−σ))`, active only when σ is within
  1e-6 of 0 or 1 (|z| > ~13.8), where cells are genuinely at bounds and contribute nothing.
  Armijo scratch buffers (`w_trial`, `b_trial`, `lambda_trial`) MUST be pre-allocated
  outside the Newton loop — per-halving allocation is 1.08 GB of transient churn at
  M_cell = 1.58M. Root cause fixed: λ=0 midpoint init → enormous first step → σ saturates
  → `D_eff → 0` → N near-singular → Tikhonov dominates → `Δλ ≈ 0` → false convergence at
  max_err = 1.0, DEFF = 527,577.

## Greenkhorn greedy step and incremental bookkeeping
- source: docs/superpowers/specs/2026-04-29-greenkhorn-solver.md
- type: protocol
- content: Each step updates ONLY `k* = argmax_k errRp[k]`, with incremental maintenance:
  per touched cell, `delta = X_new − X_old` is added to `W` and to `S[k2][g2(c)]` for
  every other margin, and `S[k*][j]` is recomputed from its bucket. `W <= 0` (all cells
  clamped to zero) is a structural infeasibility → `RK_ERR_INFEAS` with a message.
  `cells_per_cat[k][j]` is NOT in `CellTable` and must be built locally. `X_best` (lowest
  max-errRp seen), not `X`, is returned. `method="greenkhorn"` must warn-and-ignore
  `accelerate`, `scheduler`, `homotopy_*`, `eta_schedule`, `capacity_penalty`,
  `sor_enabled`. `min_weight >= max_weight` must error BEFORE the solver — `std::clamp`
  with `L > U` is UB in C++17.

## Greenkhorn round-level SQUAREM amendment
- source: docs/superpowers/specs/2026-04-29-greenkhorn-solver.md
- type: protocol
- content: SUPERSEDED by SRAA-m plus the adaptive-sort fix. Recorded for provenance: the
  amendment defined `accelerate=TRUE` for greenkhorn as round-level SQUAREM where one
  round = K greedy steps with the margin order sorted ONCE at round entry, arguing that
  fixing the sort at entry makes `F_eval` stationary. `2026-04-29-i0am-sraa-correct-all-scales.md`
  later showed that this very fixed-sort operator has a different fixed point from plain
  adaptive-sort greenkhorn and is the cause of the 35% quality regression.

## Newton-KL is obs-level, never cell-level
- source: docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md
- type: protocol
- content: Accumulation iterates raw observations via `group_ids[k][i]`, NOT
  `g_per_cell[k][c]`; no `cell_lf` is used. Dual `g(λ) = log Z(λ) − Σ T_kj λ_kj`,
  `Z = Σ_i d_i exp(u_i)`, `u_i = Σ_k λ_{k,j_k(i)}`. Reference elimination fixes
  `λ_{k,0} = 0` per margin, giving `n_λ = Σ_k (cat_counts[k] − 1)`. Single pass
  accumulates Z, gradient and the symmetric upper-triangle Hessian, then `G /= Z`,
  `H /= Z`, `H −= G⊗G`, `G −= T`. Costs: gradient O(n·K), Hessian O(n·K²/2).

## Newton-KL LM damping and trust region
- source: docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md
- type: protocol
- content: `H_damped[a,a] = max(H[a,a]·(1+lm_mu), lm_mu·d_floor)` with
  `d_floor = mean(diag(H))` — the multiplicative term preserves per-coordinate scale where
  H is well-conditioned, the additive floor covers rank-collapsed directions. Adaptive
  `lm_mu` by Marquardt gain ratio: `ρ > 0.75` and full Armijo accept → `lm_mu = max(lm_mu/3, 1e-12)`;
  `ρ < 0.25` or line search exhausted → `lm_mu = min(lm_mu·10, 1e12)`. Init 1.0, bounds
  `[1e-12, 1e12]`. A failed line search retries the LDLT solve with `lm_mu *= 10` in the
  SAME iteration, max 3 retries; three consecutive failed iterations return
  `RK_ERR_NOCONV`. Explicitly ruled out: multiplicative-only damping, spherical `μ·I`,
  trust-region clipping, tiny-step fallback.

## Newton-KL bounds fallback
- source: docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md
- type: api-contract
- content: Newton solves the UNCONSTRAINED KL dual; bounds are not enforced during
  optimisation. After convergence, if
  `count(w_i > max_weight or w_i < min_weight)/n > 0.05`, return `RK_ERR_NOCONV` with the
  Newton weights and a warning. The returned weights are ALWAYS the Newton weights
  (`w_i = f_i/Z × n`) — never a secondary solve. The status field, not the weights, tells
  the caller whether bounds held. Convergence: `‖∇g‖_∞ < tol_abs` primary,
  `step_norm < 1e-12` secondary, `max_iter = min(max_iterations, 50)`.

## Newton-KL warm start and target homotopy are NOT current behaviour
- source: docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md
- type: protocol
- content: Both subsections carry explicit negative results. Target homotopy
  (ε ∈ {0.5, 0.1, 0.03, 0.01, 0.003, 0}) was BLOCKED. The ORIS warm-start (8 SRAA sweeps,
  `λ[j] = lf_best[off+j] − lf_best[off+0]`) is falsified by an in-document ERRATUM dated
  2026-05-02: it regresses stepstone K=9 from 2.79e-4 (cold) to 4.39e-4 because ORIS's
  basin floor sits OUTSIDE Newton's quadratic basin. The section is retained for
  historical reference only. Any downstream plan that revives either mechanism must treat
  it as new work, not as documented behaviour.

## AUTO routing rule (three-way, both dispatch sites)
- source: docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md
- type: api-contract
- content: `K >= 5 && M_cell/n >= 0.9 && target_skew <= 5.0` → `RK_ALG_NEWTON_KL`;
  `K >= 5 && M_cell/n >= 0.9 && target_skew > 5.0` → `RK_ALG_IEPPA` with
  `st.accelerate = true`; `K < 5 && M_cell/n >= 0.9` → `RK_ALG_RAKING`; any
  `M_cell/n < 0.9` → `RK_ALG_IEPPA`. `target_skew = max_T / max(min_T, 1e-12)`; the floor
  guards division by zero and routes exact-zero-min targets to the severe-skew arm by
  construction. The IDENTICAL rule must be implemented in BOTH `src/c_api.cpp` and
  `src/r_bridge.cpp`. Severe-skew rationale: Newton-KL converges to a high-error fixed
  point at gap ≈ 6.24e-2 on severe-skew dual landscapes; ORIS+SRAA escapes via primal
  projection. SUPERSEDES `2026-04-23-ieppa-faithful-design.md` §4.2 ("AUTO returns
  IEPPA unconditionally, no heuristic exceptions") and the two-constant complexity/tol
  gate of `2026-04-20-algo-selection-design.md`.

---

## ORIS_SOFT ALM primal/dual update
- source: docs/superpowers/specs/2026-04-29-ieppa-alm-soft-capacity.md
- type: protocol
- content: `z = clamp(X̃, L_cell, U_cell)`; `rho = capacity_mu · X̃`;
  `X_alm = X̃·(1 − λ_cell + capacity_mu·z)/(1 + rho)`; non-finite or non-positive `X_alm`
  falls back to `z`. Dual ascent `λ_cell += capacity_mu·(X_alm − z)`, clamped to
  `±10·capacity_mu_base·max_weight` so a single step cannot exceed the cap. Derivation is
  a SINGLE Newton step from `s = 1` on `f(X) = KL(X|X̃) + λ(X−z) + (μ/2)(X−z)²`, exact at
  `X = X̃` and good while `|X − X̃|/X̃` is small — not the exact minimiser (that requires
  solving a transcendental equation). Limits: `rho → 0` gives `X = X̃(1−λ)`;
  `rho → ∞` recovers the hard clamp. Cross-referenced by `docs/methods/oris.md` and
  `CLAUDE.md`: do NOT "correct" this formula — it is right for the un-normalized-KL
  generator.

## ORIS_SOFT final projection
- source: docs/superpowers/specs/2026-04-29-ieppa-alm-soft-capacity.md
- type: protocol
- content: After the solver loop: clamp all cells, then iterate clamp+rescale at most
  `kMaxRescaleIters = 3` with tolerance `kRescaleTol = 1e-12`, breaking early when the
  total stabilises or is non-positive. Post-condition: bounds hold EXACTLY (the last
  operation is a clamp); `sum(X)` may differ from n by up to ~`3e-12·n`. Drift above
  `1e-6·n` is logged at `verbose >= 1` and reported in `alm_sum_drift`. "Soft" refers to
  the convergence path, NOT to soft bounds.

## ALM penalty schedule
- source: docs/superpowers/specs/2026-04-29-ieppa-alm-soft-capacity.md
- type: protocol
- content: Composed, not alternative: `capacity_mu = eta_i_current · capacity_mu_adaptive`
  where `eta_i_current` is the Tang-η level factor (1.0 when the schedule is fixed or
  `n_levels == 1`) and `capacity_mu_adaptive` persists across homotopy levels. Adaptive
  growth: `kAlmPersistenceThreshold = 5` consecutive iterations with
  `max_violation > kAlmTolPrimalRel(0.01)·max_weight·(n/M_cell)` multiply
  `capacity_mu_adaptive` by `kAlmGrowthFactor = 2.0`, capped at `kAlmMaxScale = 1000×`
  base. `lambda_cell` and the violation streak reset at every level transition.

## capacity_penalty API
- source: docs/superpowers/specs/2026-04-29-ieppa-alm-soft-capacity.md
- type: api-contract
- content: `harvest(..., capacity_penalty = NULL)` — NULL means auto (`M_cell/n`), a
  positive finite scalar overrides. Rejected with an error: non-numeric, length != 1,
  non-finite, `<= 0`, `> 1e15`; warned: `< 1e-15`. Passing it with any method other than
  `oris_soft` warns and ignores. Field naming is load-bearing: `alm_lambda`/`alm_mu` in
  `CalibState` are ACTIVELY USED by the lbfgsb sum-to-n ALM and must not be reused —
  the new field is `capacity_mu`, and the per-cell dual `lambda_cell[c]` lives in a
  solver-local vector, not `CalibState`. SUPERSEDED IN PART by
  `docs/superpowers/specs/2026-08-14-removed-solver-slots-supersession.md`: `lbfgsb` is removed (slot 2 reserved,
  annotated at `src/leafblower.h:44`), so the lbfgsb sum-to-n ALM that justified
  reserving `alm_lambda`/`alm_mu` no longer exists. Re-check those two fields against the
  current source before reuse or deletion — the protection note has no live consumer. The
  `capacity_mu` / `lambda_cell[c]` contract itself is unaffected.

## ALM diagnostics
- source: docs/superpowers/specs/2026-04-29-ieppa-alm-soft-capacity.md
- type: api-contract
- content: `rk_result_t` gains `alm_capacity_mu_final`, `alm_n_growth_events`,
  `alm_max_dual_norm`, `alm_sum_drift`, zero/NA for non-`oris_soft` methods. Interpretation
  contract published to users: ratio `≈ 1` means the auto value sufficed; ratio `>= 1000`
  means the growth ceiling was hit and `capacity_penalty` should be raised 10×;
  `alm_max_dual_norm` near `10·capacity_mu_base·max_weight` indicates very binding
  constraints or infeasibility; `alm_sum_drift > 1e-6·n` indicates infeasible bounds.

## Homotopy / scheduler / eta overlays are identity by default
- source: docs/superpowers/specs/2026-04-24-ieppa-homotopy-greenkhorn-design.md
- type: api-contract
- content: Three independent overlays, all default-off, all exactly identity at default:
  homotopy (`n_levels = 1` skips the outer loop entirely),
  scheduler (`ROUND_ROBIN = 0` visits margins 0..K−1),
  eta schedule (`FIXED = 0` with `eta_start = eta_end = 1.0`). Level sequence is geometric:
  `max_weight_l = start_factor·(end_factor/start_factor)^(l/(n_levels−1))`; warm start
  passes the converged vector forward WITHOUT renormalisation. Per-level KKT invariants:
  `|Σw − n| < 1e-10`, `max_k max_j |S_kj/W − τ_kj| <= tol_abs·√n_levels` at intermediate
  levels (final level uses the user tolerance), `min_weight <= w_i <= max_weight_l`, and
  `min(weights) >= 0`. `residual_recheck_fraction = 0.1` is INTERNAL to `CalibState` —
  never on the ABI, never in the R wrapper. New `rk_params_t` fields are APPENDED only;
  existing member offsets must not move.

---

## ABI tripwires are mandatory and must be re-measured
- source: docs/superpowers/specs/2026-04-25-calibration-solvers-design.md
- type: schema
- content: `EXPECTED_RK_RESULT_BYTES` is added mirroring the existing
  `EXPECTED_RK_PARAMS_BYTES`, both enforced by `static_assert`. Every field addition
  requires measuring `sizeof()` on first compile, hard-coding the observed value, and
  recording it in a comment adjacent to the assert. Recorded values drift across the
  corpus: params 152 → 216+4 = 220 → 224 → 232; result 448 → 480. `rk_params_init()` must
  `memset(p, 0, sizeof(*p))` FIRST, guarded by
  `static_assert(RK_ALG_AUTO == 0, "memset(0) default must equal RK_ALG_AUTO")` —
  otherwise stack-allocated callers get garbage in newly appended fields. Raw-struct
  callers that bypass `rk_params_init` must recompile; that is the accepted cost of adding
  API without a major version.

## Method enum values are frozen across the rename
- source: docs/superpowers/specs/2026-05-30-oris-rename-design.md
- type: schema
- content: `RK_ALG_IEPPA = 1 → RK_ALG_ORIS = 1` and
  `RK_ALG_IEPPA_SOFT = 8 → RK_ALG_ORIS_SOFT = 8` — identifiers change, VALUES DO NOT.
  Slot 2 stays reserved (LBFGSB). This keeps every serialized fixture that stores
  `algorithm_used` as an integer valid without regeneration; only fixtures storing the
  literal STRING `"ieppa"` must be regenerated. Enforced by
  `static_assert(RK_ALG_ORIS == 1 && RK_ALG_ORIS_SOFT == 8, "enum values frozen")` in a
  compiled TU. The rename is behaviour-neutral: the diff must be identifier-only, verified
  by R tests green, Python parity at rtol=1e-6, and no stepstone regression.
  EXTENDED by `docs/superpowers/specs/2026-08-14-removed-solver-slots-supersession.md`: slot 7 is likewise
  permanently reserved (was `RK_ALG_GRAKE`), and slot 12 was never occupied.

## Live algorithm enum and permanently reserved slots (authoritative)
- source: docs/superpowers/specs/2026-08-14-removed-solver-slots-supersession.md
- type: schema
- content: As of 2026-08-14, `src/leafblower.h:40-53` defines exactly eight solvers plus
  the AUTO dispatcher: `RK_ALG_AUTO = 0`, `RK_ALG_ORIS = 1`, `RK_ALG_RAKING = 3`,
  `RK_ALG_SINKHORN = 4`, `RK_ALG_CHEBYSHEV = 5`, `RK_ALG_GREG = 6`,
  `RK_ALG_ORIS_SOFT = 8`, `RK_ALG_GREENKHORN = 9`, `RK_ALG_LOGIT = 10`,
  `RK_ALG_NEWTON_KL = 11`. Slots 2 (was `RK_ALG_LBFGSB`) and 7 (was `RK_ALG_GRAKE`) are
  permanently reserved holes and MUST NOT be reused; slot 12 (`RK_ALG_CP`) was never
  occupied. Enum values are frozen (commit `77d0614`). `RK_ALG_RAKING` is 3 — the
  parenthetical giving it as 2 in `2026-04-27-sinkhorn-a1-fix-ieppa-admm-method.md` is a
  transcription error and is corrected here. Any R/Python `alg_names`-style positional
  lookup MUST carry holes at 2 and 7; the eleven-element table in
  `2026-04-29-greenkhorn-solver.md` (which pins `"lbfgsb"` at 2 and `"grake"` at 7) is
  WITHDRAWN.

## Removed solver dispositions (grake, lbfgsb, cp)
- source: docs/superpowers/specs/2026-08-14-removed-solver-slots-supersession.md
- type: api-contract
- content: `grake` — specified, partially built, removed before any release by commit
  `9a67891` (2026-04-30, "no release, no ABI constraint"); `src/grake.cpp` and
  `src/grake.hpp` do not exist and NO grake acceptance criterion carries forward
  (including A4 of `2026-04-25-calibration-solvers-design.md`, match
  `survey::calibrate(epsilon=1e-10)` within 1%). `lbfgsb` — removed; doc drift purged in
  commit `7fa211c`; a residual header note lives at `src/leafblower.h:65`. `cp` — landed
  in `00a3f10` (2026-05-03) and reverted by `3fac1d6` the same day; `grep -c 'RK_ALG_CP'
  src/leafblower.h` returns 0. Treat `2026-05-02-epic-k-cp-productionization-design.md`
  as a WITHDRAWN PROPOSAL, not pending work; if Epic-K is deliberately revived, slot 12
  becomes available and this record must be amended. This record supersedes every
  earlier spec's treatment of these three as live methods.

## Rename scope: live docs only, history preserved
- source: docs/superpowers/specs/2026-05-30-oris-rename-design.md
- type: protocol
- content: Historical records are NOT edited — dated plans, specs, investigations,
  `.beads/plans`, `tasks/`, `.wolf/`, ticket-prefixed benchmark subdirs, `docs/iEPPA/`,
  the `chu2022ieppa` bib key, `_92c4f45.*` snapshots and `*.Rcheck/` are all excluded from
  the grep-clean completion gate. Rewriting a dated spec to say "oris" would falsify the
  project record. Live docs ARE in scope, explicitly including `docs/raking.md`,
  `docs/ieppa_assessmen.md`, `docs/internal/r_bridge_floor.md` and top-level benchmark
  scripts. Completion gate: the two-stage filtered `grep -rIi ieppa` returns only the live
  files plus the intentional "renamed from iEPPA" notes.

## Rename rationale (naming honesty)
- source: docs/superpowers/specs/2026-05-30-oris-rename-design.md
- type: protocol
- content: The shipped solver is the Chu–Liang–Toh–Yang paper's INNER dual block-coordinate
  descent specialised to `C = 0` — i.e. iterative matrix scaling (RAS / Sinkhorn–Knopp /
  IPF) on log-factors with SOR and infeasibility damping. The paper's defining OUTER
  inexact-entropic-proximal-point loop is mathematically inert at `C = 0` and is NOT
  implemented. Carrying the `iEPPA` name over-claims a paper whose headline contribution
  the code does not contain. `ORIS = Over-Relaxed Iterative Scaling` states the lineage and
  the flagship accelerator without naming a paper. No deprecation alias and no migration
  shim, on the recorded basis that the package has no external users.

## Two build source lists must stay in sync
- source: docs/superpowers/specs/2026-05-30-oris-rename-design.md
- type: nfr
- content: R auto-globs `src/*.cpp` (so `PKG_SOURCES` in `Makevars.in` is decorative and
  silently drifts), while the Python build uses an explicit `CORE_SOURCES` list in
  `python/CMakeLists.txt`. Every renamed or added `.cpp` must be updated in `CORE_SOURCES`
  or the pybind11 link fails with undefined symbols. Both `R CMD INSTALL --preclean .` and
  the Python editable install are completion gates. Corollary recorded at
  `2026-05-02-ylsy-cp-ipm-spike-design.md` §2: because R globs `src/`, research code
  CANNOT be isolated at `src/research/` — only a package-root directory excluded via
  `.Rbuildignore` is off R's compilation path.

## Research isolation gate
- source: docs/superpowers/specs/2026-05-02-ylsy-cp-ipm-spike-design.md
- type: nfr
- content: `tools/check_research_isolation.R` dumps `nm -D src/leafblower.so` and exits 1
  if any of `cp_solve_R`, `ipm_solve_R`, `cp_calibrate`, `ipm_calibrate` appears, wired
  into the pre-commit hook after `R CMD INSTALL --preclean .`. `2026-05-02-epic-k-cp-productionization-design.md`
  R17 amends the list on productionisation: `cp_solve_R` and `cp_calibrate` are REMOVED
  (they become legitimate `src/` symbols) while `ipm_*` stay forbidden — and this edit must
  be the FIRST sub-step of the move, or the gate blocks every commit attempt.

## R/C boundary contract for any SEXP-exposing TU
- source: docs/superpowers/specs/2026-05-02-ylsy-cp-ipm-spike-design.md
- type: api-contract
- content: `#define STRICT_R_HEADERS 1` and `#define R_NO_REMAP 1` before any R header;
  include only `<R.h>`, `<Rinternals.h>`, `<Rdefines.h>`; heap allocation via
  `std::vector<double>` or Eigen only, raw `new`/`delete` forbidden; every
  `Rf_alloc*`/`PROTECT` balanced by `UNPROTECT(n)` before return with a manual reviewer
  audit; input dimensions validated before any pointer dereference; no `Rcpp::` types
  (avoids LinkingTo dependency drift). Hard input bounds enforced in the shim before any
  allocation: `n <= 1e8`, `n_col <= 1e6`, `max_iterations <= 1e5`.

## Chambolle-Pock production contract (WITHDRAWN PROPOSAL — never landed)
- source: docs/superpowers/specs/2026-05-02-epic-k-cp-productionization-design.md
- type: api-contract
- content: `RK_ALG_CP = 12` appended to `rk_algorithm_t`; `method="cp"` added to
  `kAlgMap`, to the `match.arg` whitelist (11 → 12 entries) and to the `accelerate`
  whitelists. `accelerate` defaults to TRUE for CP ONLY, implemented as
  `accelerate_explicit <- !missing(accelerate); if (method == "cp" && !accelerate_explicit) accelerate <- TRUE`
  so the global `accelerate = FALSE` signature default is unchanged for every other method.
  No `rk_params_t` field changes — `cp_safety_factor = 1.05` is hardcoded and
  `EXPECTED_RK_PARAMS_BYTES` is unchanged. The dispatch arm must be placed BEFORE the
  catch-all `else`. Result list grows from 37 to 47 elements (slots 37–46:
  `n_cells`, `algorithm_requested`, `algorithm_used`, `A_norm_estimate`, `n_power_iter`,
  `final_theta`, `final_tau`, `final_sigma`, `fell_back_to_pdhg`, `wall_time_ms`), NA/0/""
  for non-CP solvers.

## Chambolle-Pock algorithm and fallbacks (WITHDRAWN PROPOSAL — never landed)
- source: docs/superpowers/specs/2026-05-02-epic-k-cp-productionization-design.md
- type: protocol
- content: Algorithm 1 (vanilla PDHG, O(1/k)): `σ = τ = 1/(‖A‖·1.05)`,
  `y ← y + σ(Aᵀw̄ − b)`, `w ← prox_{τf}(w − τAy)`, `w̄ ← 2w − w`. Algorithm 2
  (accelerated, O(1/k²)) requires γ-strong convexity, `γ = 1/max(u_i)` taken from the
  SCALAR `st.max_weight`; `θ = 1/√(1+2γτ)`, `τ ← θτ`, `σ ← σ/θ`,
  `w̄ ← w + θ(w − w_prev)`. Two mandatory fallbacks: `max_weight = Inf` → γ = 0 → run
  Algorithm 1, set `fell_back_to_pdhg = true`; `θ_k < 1e-15` → RESET σ,τ to the fixed
  Algorithm-1 values and continue (do NOT freeze τ at the underflowed value — that stalls
  progress). `prox` is a direct scalar Newton (NO Lambert-W), ≤5 iterations, with an
  asymptote branch when `|z/(τd) − 1| > 700`. Cell compression is default; obs-level is
  forced when `M_cell·10 > n·9` or `bounds_mode = "unit"` (per-obs bounds are not
  cell-aggregable). Power iteration: 50 iterations, stop at 1e-6 relative, `status_code = 4`
  if the relative delta still exceeds 1e-3 at k = 50.

## CP is opt-in and is not wall-time competitive (WITHDRAWN PROPOSAL — never landed)
- source: docs/superpowers/specs/2026-05-02-epic-k-cp-productionization-design.md
- type: nfr
- content: AUTO routing is explicitly UNTOUCHED — `method="cp"` must be passed
  explicitly, and `test-algo-selection.R` must still pass. The published claim is
  quality-at-budget, NOT speed: on stepstone_K9, cp reaches 5.08e-5 at 5000 iterations
  versus ORIS+SRAA's 4.39e-4 bounded fixed point at 200 iterations — ~2× tighter weights —
  while taking 52s against ORIS+SRAA's 0.34s. Gate T2 is therefore explicitly labelled
  "quality-at-budget (NOT wall-fair)" and paired with T2b, a loose 90s wall-time SANITY
  ceiling (73% headroom over the measured 52s), not a fairness gate. CP is NOT recommended
  for severe-skew K≥5 (`target_skew > 5`), where the spike showed it fails to converge.

## No per-iteration trace in production
- source: docs/superpowers/specs/2026-05-02-epic-k-cp-productionization-design.md
- type: nfr
- content: `CpCalibResult` stores final-iteration scalars only; memory is bounded by
  O(n + ΣJ_k). The 1000-row spike trace (per-iter `max_err_last`, `max_err_ergodic`,
  `primal_resid`, `primal_stationarity_proxy`) is research-only. Rate exponents
  (β_last, β_ergodic, R²) are spike diagnostics and are NOT runtime-exposed — production
  users get `max_error` and final residuals.

## Faithful-textbook constraint for research spikes
- source: docs/superpowers/specs/2026-05-02-ylsy-cp-ipm-spike-design.md
- type: protocol
- content: No hyperparameter tuning, no momentum/restart heuristics, no warm start, no
  accelerated variants inside a spike: CP runs `σ = τ = 1/(‖A‖·1.05)` unaccelerated, IPM
  runs `μ₀ = 1` with 0.5 decay and a vanilla central path (no Mehrotra). ALL state in
  `double` — the extrapolation `w̄ = 2w^{k+1} − w^k` is catastrophic-cancellation-prone in
  single precision. Sanity gate before any benchmark: on a fixture whose analytical
  solution is `w = d`, both solvers must reach `max_err < 1e-8` with `status_code = 0`;
  the bench halts if sanity fails. Decision rule is pre-registered with a rate-exponent
  floor (β ≤ −0.8 last-iterate, ≤ −1.0 ergodic, R² ≥ 0.9, ≥ 30 fit points) so a PASS
  cannot be claimed on max_err alone.

## Determinism and measurement protocol
- source: docs/superpowers/specs/2026-05-02-ylsy-cp-ipm-spike-design.md
- type: nfr
- content: `OMP_NUM_THREADS = 1` and `RhpcBLASctl::blas_set_num_threads(1)` in every bench
  preamble; kk1204 measured 3× with the median taken (thermal noise on a 30s gate);
  hardware metadata (`Sys.info()`, CPU model, frequency mode) logged; pre-flight memory
  check `needed = n·K·8 + n·8·5 + max_iter·10·8` must be under half of available RAM.
  Interleaved comparison via `bench::mark(min_iterations = 3)` is required for wall-time
  claims — sequential `system.time()` is dominated by GC and scheduler variance.

## Benchmark-gate discipline
- source: docs/superpowers/specs/2026-04-20-algo-selection-design.md
- type: protocol
- content: Constants derived from a benchmark are committed ONLY when the fit is
  trustworthy: the Bayesian level-set loop stops at 90% GP classification or 25
  acquisitions, and on reaching the cap with `classified < 0.90` it prints a prominent
  warning, does NOT update `src/c_api.cpp`, and files a follow-up. K-stability is checked
  at K ∈ {3, 18} against K = 9; contour shift beyond 0.5 log-units means the constants are
  labelled Stepstone-calibrated rather than general. Checkpoints every 5 acquisitions are
  written to a temp file then renamed atomically. The lower tolerance bound is 1e-6
  precisely to keep convergence-failure observations out of the GP.

## Test-only environment overrides
- source: docs/superpowers/specs/2026-04-24-ieppa-speed-convergence-bounds-design.md
- type: protocol
- content: `LBW_ORIS_FORCE_PATH` (linear|log|unset), `LBW_ORIS_FORCE_DAMPING`
  (on|off|unset) and the acceleration toggle are read once per solver entry via
  `std::getenv` — microsecond cost, no compile-time divergence, so the CRAN-shipped binary
  and the tested binary are identical (a `#ifdef LBW_TEST_HOOKS` macro cannot cleanly
  distinguish `R CMD INSTALL` from a CRAN build). Documented as test-only in source
  comments, absent from user-facing docs. Only `FORCE_PATH × FORCE_DAMPING` is
  orthogonal-by-construction; acceleration interacts with damping (it reads `alpha` and
  maintains the shadow `lf`) and needs the 7 explicitly enumerated combination tests.

## RED-test mechanics must be struct-based
- source: docs/superpowers/specs/2026-04-24-ieppa-speed-convergence-bounds-design.md
- type: protocol
- content: Compile-time instrumentation and log-string parsing are BANNED as test
  mechanisms; assertions read integer/double counters on the result struct that are always
  computed and zero when not applicable — `n_xcur_writes_per_iter_linear`,
  `n_anderson_iters_engaged`, `n_anderson_nan_fallbacks`, `min_alpha_seen`, `final_alpha`.
  Iteration-count assertions must be monotone (`iter_damped > iter_stable`), never
  ratio-bounded, because ratios flake across OS and BLAS.

## Fixture regeneration is a gated manual step
- source: docs/superpowers/specs/2026-04-27-raking-bregman-dykstra-design.md
- type: protocol
- content: Changing the fixed point invalidates reference fixtures. The protocol is:
  verify the quality acceptance criterion FIRST, halt if it fails, and only then rerun the
  generator — never regenerate to make a test pass. `2026-04-27-sinkhorn-a1-fix-ieppa-admm-method.md`
  adds the ordering for objective-reporting fixes: implement, compile, regenerate, commit
  fixture and script together, then run tests — regenerating against old code bakes in
  wrong values. `2026-04-29-ieppa-alm-soft-capacity.md` adds the inverse case: capture the
  pre-change reference fixture and COMMIT IT BEFORE any source edit (verified by
  `git log` ordering).
