# Benchmark Study Design — leafblower vs competing R/Python calibration software

**Epic:** leafblower-2ouc · **Status:** design v5 (post external Gemini review; gate-converged v4 + methodological hardening) · **Date:** 2026-07-08
**Downstream:** scientific article (journal-agnostic harness; venue chosen post-results).
**Revision note:** v2 addresses design-review round 1 (all-5 NEEDS_REVISION). Change log at §14.

---

## 0. Alternatives Considered (build-vs-adopt)

| Option | Fit | Verdict |
|---|---|---|
| **Benchopt** (NeurIPS'22 cross-language solver bench) | Python-centric; R solvers run only via `rpy2` shims. But **most competitors are R-only** (survey, sampling, ebal, optweight, GECal, ReGenesees, icarus, laeken, jointCalib). Benchopt cannot host a first-class R arm, and its objective/solver abstraction assumes one shared objective — we deliberately span *different* objectives (KL/χ²/logit/minimax). | **Rejected** — cannot host the two-language, multi-objective design. Adopt its *methodology* (convergence-vs-time curves, reproducibility discipline), not the framework. |
| **`asv` / airspeed-velocity** | Python microbenchmark harness; no cross-language, no solver-quality axis. | Rejected — timing-only, single-language. |
| **Extend `benchmarks/stepstone_all_methods.R`** | Already has per-method `run()` with INFEAS `tryCatch`, `fit_metrics()` (margin-KL/ESS/DEFF/χ²/L∞/L1), the Pearson-r within-family agreement block, and the R→Python subprocess bridge. | **Partially adopt** — the new harness **extends/refactors** these (see §8), not parallel-stacks them. |
| **Bespoke harness (this design)** | Two arms (R, Py), adapter-normalised, reusing the proven metric code. | **Selected** — only option that hosts a first-class R arm + multi-objective + reuses trusted metrics. |

**Selected mechanism:** bespoke two-arm harness that *reuses* `fit_metrics`/`compute_metrics` and the `run()` proto-adapter, and *adopts* Benchopt-style reporting (Dolan–Moré profiles, work-precision, frozen reproducibility bundle).

---

## Mechanism / Forbidden / Audit

- **Mechanism:** language-split, adapter-normalised harness. Each competitor + leafblower method wrapped in a uniform adapter `run(problem) -> {weights_ref, iterations|NA, status, converged, error_message, wall_time, peak_rss}`. Metrics recomputed from returned weights by shared code refactored out of the existing `fit_metrics`/`compute_metrics`. **Primary comparison = work-precision (achieved accuracy vs cost)**; wall-time is secondary and always annotated with independently-recomputed accuracy.
- **Forbidden:** cross-language timing comparison; renormalising weights after a solver's own finalize; comparing weights *across objective families* as error; trusting a solver's self-reported `converged`; per-solver hand-tuned tolerances; letting competitor deps enter DESCRIPTION/`pyproject` (would break the DoD gate); running any timed measurement before the spec+matrix+metrics are git-frozen.
- **Audit:** every metric recomputed independently of the solver; R and Python metric impls validated against a **hand-computed golden** on the toy (not just mutual R↔Py agreement) then to `rtol=1e-6` (repo parity bar; 1e-9 not required — 1.58M-row summation order differs); each adapter validated on **both** the shared toy **and its own package's documented home-turf result**; environment frozen via lockfiles + git SHA.

---

## 1. Objectives & research questions

Compare leafblower (8 methods + `oris`, `oris_soft`) against direct competing OSS within each language, on (a) each competitor's canonical problem and (b) stepstone (unbounded + bounded), plus (c) an **independent non-home-turf synthetic** scaling problem.

- **RQ1 (speed):** work-precision (accuracy vs wall-time) primary; wall-time secondary; FFI overhead separated *for leafblower only* (§5). Threads ∈ {1,4}.
- **RQ2 (memory):** peak RSS.
- **RQ3 (work):** iteration count where exposed; margin-error-vs-iteration curves for instrumented solvers.
- **RQ4 (quality):** margin KL, ESS (↑), design effect (↓), margin L∞/L1.
- **RQ5 (agreement):** within an objective family on the same problem, are weights tightly correlated? **Pre-registered thresholds:** "tight" ⟺ Pearson ≥ 0.999 **and** max|Δw|/mean(w) ≤ 1e-3 **and** Spearman ≥ 0.999.
- **RQ6 (SOTA gap):** **narrative deliverable** (§9), not a measured RQ — enumerated gaps for the article.

---

## 2. Competitor selection (canonical + maintained alt; locked)

✓ maintained · ⚠ stale/dead (2026). Applicability is gated by **objective family**, **bounds capability**, and **margin count K** (critical: see OT note).

### KL / IPF family — competes with ORIS, ORIS_soft, Raking, Newton-KL (same KL optimum, any K)
| Arm | Package | Bounds | K-margin? | Note |
|---|---|---|---|---|
| R | `survey::calibrate(raking)` ✓ | yes | yes | authoritative canonical |
| R | `anesrake` ⚠2018 | cap | yes | most-cited; also **canary** (known-dominated harness check) |
| R | `ipfr` ✓ | ratio box | yes | maintained raking alt |
| Py | `ipfn` ⚠2021 | no | yes | canonical Python IPF (in venv) |
| Py | `weightipy` ✓ | no | yes | maintained RIM alt |

### Entropic OT — NO valid external competitor for K-margin calibration (redundant-IPF cross-check only)
> **Reframed per external review:** at K=2 with zero cost, balanced Sinkhorn's fixed point **is** the 2-margin IPF/raking solution, so POT/ott reproduce leafblower's raking weights *by construction* (`w_i=u_i·d_i` just re-derives IPF). They therefore test "another IPF library," not entropic OT. **Conclusion (a §9 finding, not a benchmark claim): there is no external K-margin entropic-OT competitor for leafblower's `sinkhorn`/`greenkhorn`.** POT/ott are retained only as a **redundant-IPF cross-check** on the 2-margin reduction + their OT home-turf (1D-Gauss), explicitly labelled as such — not presented as a head-to-head.
POT `ot.sinkhorn`/`ot.bregman.greenkhorn` and `ott-jax` solve **2-marginal** entropic OT with a transport **cost matrix** and return a coupling. leafblower's `sinkhorn` is a **K-margin Dykstra projection with no transport cost** — a *different object*. They are directly comparable **only at K≤2 with zero cost**. Therefore:
- POT/ott run on their **OT home-turf** (1D-Gauss) and a **2-margin calibration reduction** of stepstone (project any two margins; spec `id` distinct from the K=9 id so weight files never collide) — **not** the K=9 head-to-head.
- **Coupling→weight extraction** (adapter golden must pin this): a *balanced* Sinkhorn coupling has row-sum `Σ_j P_ij ≡ a_i` exactly, so a naive row-sum ratio degenerates to `w_i = d_i` (uncalibrated). The calibration weight must therefore come from the **left Sinkhorn scaling potential** `u` (`w_i = u_i · d_i`) or from an **unbalanced** OT (`ot.unbalanced`, free row marginal). **The OT-adapter WU pins the exact extraction (potential vs unbalanced) against a hand-solved 2-margin IPF golden**; if neither reproduces the IPF weights, the OT arm degrades to the §9-item-4 "no comparable K-margin competitor" finding rather than shipping a wrong mapping.
- On K≥3 calibration there is **no external Sinkhorn/Greenkhorn competitor** — itself a §9 finding.
- ott-jax is **excluded from ranked {1,4}-thread timing** (XLA ignores OMP/BLAS env; CPU-vs-GPU not comparable) — reported as a separate scaling vignette with `XLA_FLAGS` thread pinning disclosed.

### χ² linear — GREG
R: `survey::calibrate(linear)` ✓, `sampling::calib(linear)` ✓, `ReGenesees` (⚠ not on CRAN 2026 — install from GitHub DiegoZardetto/ReGenesees or drop). Py: `svy` ✓ (samplics successor).

### Logit (bounded Deville–Särndal)
R: `survey::calibrate(logit,bounds)` ✓ canonical, `sampling::calib(logit)` ✓, `icarus`/`laeken`/`ReGenesees` ✓. **Py: gap** — no maintained bounded DS logit (§9 finding).

### Newton-KL (trust-region dual)
R: `jointCalib` ✓ / `nonprobsvy` ✓ (double-dogleg TR dual — closest analog). Py: `scipy.optimize(trust-constr,bounds)`.

### Chebyshev (minimax L∞)
R: **`optweight(norm="linf")`** ✓ (only direct competitor), `WeightIt`. Py: `cvxpy` LP baseline.

### ORIS_soft (ADMM soft capacity)
R: `ebal` ✓, `optweight(norm="entropy")` ✓, `WeightIt`. Py: Meta `balance` ✓, POT `ot.unbalanced`.
> ebal/optweight use causal means-balancing framing; the adapter expands each categorical margin to indicator columns (mean ≡ target proportion) — **validated on the toy** (not asserted).

### Unifying baseline
R: `GECal` ✓ (generalized entropy: χ²+raking+logit).

---

## 3. Problems & spec schema

**(A) Canonical home-turf** — survey `apistrat` (200×39), sampling `belgianmunicipalities` (589×17), anesrake `anes04` (1212×5), ebal/optweight `lalonde` (614×9), ipfn 3×4 seed, POT 1D-Gauss (100-bin), balance `simulate_data()`.
**(B) stepstone** — reuse existing `benchmarks/stepstone_fulldata_bench_data.parquet` + `_targets.json` (n≈1.58M, K=9), unbounded + bounded (incl. eb79.23 INFEAS max_weight=5).
**(C) Parametric instance family (non-home-turf)** — a **suite of 30–100 synthetic instances** sweeping n ∈ {10k,100k,1.58M}, K, cell sparsity, margin skew, constraint-matrix conditioning, and infeasibility slack (fresh seeds, not derived from leafblower fixtures). **This is what makes the Dolan–Moré performance profiles (§7) statistically valid** — a profile over only ~4 problems is noise; profiles need a problem *distribution*. The conditioning/infeasibility sweeps double as the ill-conditioning/robustness axis (§9 item 6). n-scaling curves are the marginal-over-n view of this suite. **DGP frozen with the other specs and finalized BEFORE planning/WU decomposition** (not deferred); **preferably sourced from an externally published synthetic-calibration generator** (cited), or if in-house, its margin count/skew/cell-sparsity recipe designed independently of leafblower's known convergence behavior and reviewed as such — recipe published in supplementary. Stepstone conclusions scoped as **applicability demonstration**, not "general superiority."

**Spec schema (`spec/*.json`)** — inline only scalars/small arrays; **reference** bulk data via a typed `data_ref` with a **loader convention** covering all three dataset origins:
- `"file:<path.parquet>"` — bulk columnar (stepstone).
- `"pkg:survey::apistrat"` — package-bundled canonical dataset (loaded via the package API in each arm; the loader field maps to `data(apistrat, package="survey")` / Python equivalent). Snapshotted into the frozen spec bundle to guard against in-place upstream data changes.
- `"gen:<recipe_id>"` — synthetic (§3C), generated from a frozen recipe + seed.
- inline arrays — small toy only.
```json
{ "id":"stepstone_bounded", "data_ref":"file:benchmarks/stepstone_fulldata_bench_data.parquet",
  "design_weights":"uniform|<column>", "margins":["age","region",...],
  "targets_ref":"..._targets.json", "bounds":{"min":0,"max":5},
  "tol":1e-8, "objective_families":["kl","chi2","logit"], "K":9 }
```
Note: harvest's per-obs arg is **`design_weights=`** (no `weights=` — silently swallowed).

---

## 4. Run matrix — CORE vs STRETCH (feasibility-gated)

**Cost estimate:** stepstone-fulldata single passes already cost seconds–minutes per method (chebyshev ~12s, greg ~2s in `stepstone_all_methods.R`); ~20 solvers × 2 stepstone forms × threads × reps ⇒ full matrix plausibly **hours–days**. Therefore phased:

- **CORE (v1, satisfies the article):** every applicable solver on (its canonical home-turf) + (stepstone unbounded + bounded) + (**the §3C parametric instance family, 30–100 instances** — required for valid Dolan–Moré profiles) at **1 thread**, **weights once** (deterministic) + **timing reps=10 on canonical/family, reps=5 on 1.58M**. The instance family is deliberately **weighted toward small/medium n** (most instances 1k–100k; only a handful at 1.58M) so profile *resolution* comes from instance *count*, not per-instance cost.
- **Numeric CORE budget (enumerated ceiling):** ~22 applicable (solver,arm) cells. Canonical: ~22×10 reps × ≤1s ≈ minutes. stepstone×2 forms × ~15 applicable solvers × 5 reps; dominated by second-order solvers (chebyshev ~12s, greg ~2s per pass on 1.58M) → ≈ 15×2×5×~10s ≈ **~25–40 min compute**, ×2 leafblower build variants for leafblower rows only. Scaling {10k,100k,1.58M} adds ≈15 min. Instance family: 30–100 mostly-small instances × ~15 solvers × 10 reps × ≤~0.1s median ≈ **~10–30 min**. **CORE ceiling ≈ 2.5 h wall-clock**; if a projected cell set exceeds it, cap the family at 30 instances / drop STRETCH duplicates first. STRETCH (thread=4 sweep) roughly doubles stepstone time.
- **STRETCH (appendix):** thread=4 sweep; maintained-alt duplicates (ipfr, weightipy, ott-jax scaling vignette); GECal unifying baseline.

Applicability gated by a **machine-readable registry** (`registry.json`: solver → {arm, families, bounds, K_max, applicable_problems}) built against **actually-installed versions** (not the research table).

---

## 5. Measurement protocol (within-language)

- **Cell-table boundary (standardised):** the cell/indicator construction (groupby over rows) is built **once, outside the timer, for every adapter**, so "solve time" is comparable. This fixes the `python_ipf_benchmark.py:66-77` asymmetry (Polars groupby timed inside the solve). **leafblower requires a new capability:** its public `harvest()` builds the cell table *inside* the C++ solve (`build_cell_table`, `src/cell_table.hpp:19`, which already accepts pre-coded `group_ids`). The seam is separable, so **WU-0 adds a thin lower-level R/Python entry point that accepts a pre-materialised CellTable and runs solve-only** — a new public surface, so it carries its own parity test (R↔Py + vs `harvest()`) under the DoD gate. **Precise boundary** (WU-0 ticket): "cell construction" = `build_cell_table` + the O(n) `X_init` scatter over `cell_of` (`calib_dispatch.hpp`); bounds/cat_offset/validate (O(M_cell)) stay **inside** the timer as genuine solve setup — else the fix over-credits leafblower. **Fallback** if WU-0 is descoped: measure leafblower's internal cell-build separately (instrument the C++ boundary) and subtract; documented per result.
- **leafblower binary under test:** the Python arm MUST import the **venv-built** `leafblower` (`python/.venv/bin/python`, after `uv pip install -e . --reinstall-package leafblower`), never the stale `~/.local` shadow `.so` that bare `python`/`pytest` picks up (CLAUDE.md build note). Measurement-validity precondition, asserted at harness start (import path + `__file__` check).
- **Timing:** in-process high-res (`bench::mark`/`time.perf_counter`), ≥2 warmups discarded; **median + 5/95 quantiles + bootstrap CI** (never mean±sd). **Minimum-timed-duration rule**: sub-ms canonical solves are inner-loop-batched to ≥50ms total; per-result timer-resolution diagnostic recorded. FFI/serialization split is measured **only for leafblower** (Rcpp/pybind11 boundary); **N/A for pure-R/opaque-C competitors** (no separable boundary) — stated, not faked.
- **Memory:** **one fresh subprocess per `(solver,problem,thread,build)` cell** (VmHWM is a monotone process high-water mark — batching solvers in one session would inherit prior peaks and mis-attribute RSS); peak RSS via `/proc/self/status:VmHWM`; `gc()` before timer; **per-language empty-process baseline subtracted** (R and Python startup footprints differ — one baseline each, not pooled; DoD item).
- **Iterations:** solver's exposed count; `NA` (typed null in parquet, never 0/"NA") where hidden.
- **Threads:** sweep {1,4}; set `OMP_/OPENBLAS_/MKL_NUM_THREADS` + per-package controls identically; enumerate packages with independent parallelism (foreach/joblib) and pin them. **ott-jax excluded from ranked timing** (XLA pinning disclosed separately).
- **Tolerance handling:** the shared `tol` is passed to each solver's native stopping parameter **with a documented per-solver mapping table** (what quantity it bounds: margin L∞ vs duality gap vs KL). Because equal `tol` ≠ equal achieved accuracy, **work-precision is the primary axis** and every wall-time is annotated with independently-recomputed achieved margin error (**achieved-accuracy-at-nominal-tol** diagnostic reported per solver).
- **Competitor hyperparameters:** **pre-registered table** (Sinkhorn ε, max_iter, optweight tols, trust-constr settings) fixed to each package's *documented/recommended* values, committed before runs.
- **Build-flag parity (control, not disclosure):** wall-time is part of the primary work-precision axis, so the `-O3 -march=native` (leafblower `~/.R/Makevars`) vs shipped-portable-binary asymmetry is a real confound. Control: build **two leafblower variants** — (i) `native` (`-O3 -march=native`) and (ii) **`portable`** (CRAN/PyPI default `-O2`, no `-march=native`). **Headline RQ1 figures use the `portable` variant** (matched to competitors' shipped builds); `native` reported separately as a "tuned-build" delta. The exact `CXXFLAGS`/`-march`/`~/.R/Makevars` of both, **plus the linked BLAS/LAPACK backend + version + runtime CPU-dispatch mode** (OpenBLAS/MKL auto-select native SIMD kernels at runtime regardless of caller build flags — a residual native-instruction path for leafblower *and* competitors), are captured in `environment.json`. Article speed language is capped to **accuracy-normalised, build-flag-annotated** claims — no bare "N× faster".

---

## 6. Quality metrics (refactored from existing `fit_metrics`/`compute_metrics`; golden-checked)

Recomputed from returned weights, independent of the solver, cancellation-free (`p(1-p)`, zero-design rows excluded):
- **margin KL** `KL_k = Σ_j T_kj log(T_kj/p_kj)`, `p_kj=Σ_i w_i[i∈j]/Σw`; mean/max over k. **Divergent terms are NOT dropped** — when a failing solver starves a category (`p_kj→0`, `T_kj>0`) the term → +∞; report as a large finite sentinel + a `divergent` flag, never silently excluded (the existing `t>0 & S>0` gate in `fit_metrics:34` / `compute_metrics:46` **flatters catastrophic failures** and must be removed — the golden asserts the +∞ behavior on a starved-category toy). Aligns with the project no-hiding-failures ethos.
- **ESS** (Kish) `(Σw)²/Σw²` (↑). **Kish weighting DEFF / UWE** `= 1+CV²(w) = n·Σw²/(Σw)² = n/ESS` (↓) — **explicitly labelled Kish weighting-loss, NOT the true design effect** (which depends on the estimand and the w–y correlation). Since DEFF ≡ n/ESS, report **one** (ESS) as the RQ4 quality axis and DEFF only as its labelled reciprocal, not a second independent metric.
- **margin error** L∞, L1. **closeness-to-design** `Σ w log(w/d)` with the **per-problem design weights `d_i`** (NOT hardcoded `d_i=1` as in `fit_metrics:37` / `compute_metrics:48`; apistrat/belgianmunicipalities/lalonde have real `d_i≠1` — the golden pins `d_i` per problem).
- **bound violation** (bounded problems): count and max/mean magnitude of `max(0, L_c−w_i, w_i−U_c)`.
- **Weight agreement (RQ5):** within objective family — Pearson, Spearman, max|Δw|, cosine — against pre-registered thresholds (§1). **Anchored to an independent high-precision convex reference solve** (CVXR/ECOS to ~1e-12 for the family's objective) on **canonical + instance-family problems only** (a 1e-12 convex solve is infeasible at stepstone 1.58M — that scale has no external anchor, stated as a limitation); the reference weights are **produced + stored as a pseudo-solver row** so reporting consumes stored vectors, not a re-solve, **never to an in-house leafblower method** (the existing `Pearson-r-vs-ORIS` anchor at `stepstone_all_methods.R:108` measures closeness-to-leafblower, not correctness — replaced). Cross-family reported as "different optima," never scored as error.

Metric code = `common/metrics.{R,py}`. **Three existing impls must be reconciled** (they differ): `stepstone_all_methods.R:fit_metrics`, `python_ipf_benchmark.py:compute_metrics`, and `allmethod_bench.R:compute_metrics` (differently-normalised `weight_kl` and `chi2`). The **hand-computed golden on the toy is the canonical arbiter** — it selects the correct definition; all three are validated against it (not merely against each other, which two identical-but-wrong formulas would pass), then R↔Py to rtol=1e-6 (1e-9 unreachable over 1.58M-row summation order).

---

## 7. Reporting (full SOTA)

Dolan–Moré **performance profiles** (over the §3C instance family — dozens of instances, not the ~4 canonical problems — computed **within objective-family strata**, not pooled); **work-precision diagrams** plotting each family's error on **its own native divergence** (KL for the KL family, χ² for GREG, logit-distance for Logit), with **margin-L∞ used only as a neutral constraint-satisfaction axis carrying an explicit caveat that Chebyshev (and any margin-L∞-stopped solver) has a home-field trajectory advantage on that axis**; RQ1 owns the *wall-time*-axis panel, RQ3 the *iteration*-axis panel (no duplicate figures); convergence curves (instrumented solvers); quality tables (median+quantile); weight-agreement heatmaps per family. Stale (⚠) packages labelled "historical baseline." **Every speed figure caption states the leafblower build variant** (`portable` for headline; `native` deltas flagged) and that timings are accuracy-normalised.
**INFEAS scoring (pre-registered):** on the bounded problem, a solver that cannot return bounded weights is **DNF (τ=∞)** in the performance profile; leafblower's INFEAS-halt is reported **both** as a profile DNF **and** narratively (guarded refusal vs competitors' silent bound violation) — the usability downside of halting is stated, not hidden.

---

## 8. Harness architecture (reuse-first)

```
benchmarks/study/            # ISOLATED — no deps leak into DESCRIPTION/pyproject (DoD gate stays green)
  spec/                      # problem specs (JSON, data_ref to parquet)
  registry.json              # solver applicability (arm/families/bounds/K_max/problems)
  R/competitors.R            # adapters, extend existing run() proto-adapter + INFEAS tryCatch
  R/run_arm.R                # driver: solver × problem × threads × reps
  python/competitors.py, run_arm.py
  common/metrics.{R,py}      # REFACTORED from fit_metrics / compute_metrics
  common/problem_io.*        # spec loader (wraps existing parquet/targets)
  results/
    runs.parquet             # tidy: one row per (solver,problem,thread,build,rep). Columns (dtypes):
                             #   solver:str, problem:str, thread:int, build:str{portable|native|na},
                             #   rep:int, weights_ref:str, iterations:int?, status:cat, converged:bool,
                             #   error_message:str?, wall_time_s:f64, peak_rss_bytes:i64,
                             #   + denormalised quality-metric cols (marg_kl,ess,deff,linf,l1,bound_viol,...)
    weights/<solver>__<problem>__t<thread>__<build>.parquet  # length-n weights, ONE per
                             #   (solver,problem,thread,build) — deterministic across reps; build="na"
                             #   for competitors; weights_ref = this key
  report/                    # quarto: profiles, work-precision, tables, heatmaps
  env/                       # renv.lock (R), uv.lock (Py), git SHA, environment.json
                             #   (incl. CXXFLAGS/-march/~/.R/Makevars for BOTH leafblower variants)
  run_all.sh                 # orchestrator + env capture; --smoke (n=5000, reps=2) fast dry run
```

**Adapter contract:** `run(problem) -> {weights_ref, iterations|NA, status, converged, error_message, wall_time, peak_rss}`.
- `status`: **harmonized enum** (stored as string category) {converged, no_conv, infeasible, bound_violation, bad_arg, budget, stall, error} with a **per-adapter mapping table** (leafblower RK_ERR_* → enum; each competitor's native codes/exceptions → enum). `bound_violation` = returned weights breach `[L_c,U_c]` beyond tol (the unbounded-competitor-on-bounded-problem case). `error_message` captures `conditionMessage(e)`/exception text (as `allmethod_bench.R:57-68` already does).
- `converged`: **harness-recomputed** — a **uniform margin-L∞ threshold at the spec `tol`** applied to ALL solvers (so `converged=FALSE` reflects achieved accuracy, not each solver's native tol-semantics), never the solver's self-report.
- `weights_ref` = the `(solver,problem,thread)` key resolving to one `weights/` file; the ≤10 reps of a given `(solver,problem,thread)` all share it (weights deterministic across reps). `wall_time`/`peak_rss` stored as **seconds float64 / bytes int64**; `iterations` as typed-null-able int.
- Field name **`iterations`** (matches `src/types.hpp:89` `CalibResult`, `_harvest.py`), not `n_iter`.
- **OT/ebal adapter semantics** (cost-matrix mapping; indicator-column expansion) validated on the toy + each package's home-turf golden.

---

## 9. SOTA gap analysis (RQ6 narrative deliverable)

> All claims below are **draft framings pending §6/§7 results** — each is confirmed against the measured pipeline (esp. items 4–5, now backed by the bound-violation metric and K-margin applicability registry) before appearing in the article.

1. **Bounded Python gap** — no maintained Python bounded Deville–Särndal logit/GREG (GECal R-only).
2. **Minimax survey calibration — a software-availability gap, not a methodological novelty.** L∞ calibration is a plain LP (the design itself lists a cvxpy LP baseline); only `optweight` exposes it in the calibration/causal ecosystem. Framed as "no off-the-shelf survey library offers L∞," **not** as a new method.
3. **TR-Newton on dual** — dogleg exists (`jointCalib`); novelty framed as *true trust-region Newton on the classical DS dual* — **restricted to the strictly-convex smooth families (KL, χ², logit)**; the claim does **not** extend to minimax/L∞ (non-smooth, inequality-constrained dual — not twice-differentiable). Narrative claim, not a measured RQ.
4. **No external K-margin entropic-OT competitor** — POT/ott are 2-marginal balanced OT that reduce to IPF; leafblower's K-margin `sinkhorn`/`greenkhorn` have no OT counterpart to benchmark against (§2).
5. **Infeasibility robustness** — leafblower INFEAS pre-check + guards vs silent violation/divergence (both sides reported; backed by the bound-violation metric §6).
6. **Downstream variance/replicate-weight estimation** (measured, not just acknowledged) — the reason practitioners use `survey`/GREG. Measure the cost of the generalized-Jacobian / B replicate re-solves; this reweights the whole speed story (a solver cheap to fit once but expensive to re-solve B times ranks differently).
7. **Other acknowledged axes:** GPU/JAX scaling (leafblower CPU-only); streaming/out-of-core.

---

## 10. Risks & controls

| Risk | Control |
|---|---|
| Heterogeneous optima → invalid weight comparison | family-scoped agreement only (§6) |
| OT adapter mis-maps calibration semantics | toy + home-turf golden gate; POT/ott K≤2 only (§2) |
| Shared tol ≠ equal accuracy → biased speed | work-precision primary; achieved-accuracy diagnostic (§5) |
| Home-turf bias on stepstone | independent synthetic scaling problem (§3C); stepstone scoped as demonstration |
| Selective reporting | **git-freeze** spec+matrix+metrics+hyperparams before any timed run, anchored by an **external immutable timestamp** — **Zenodo/OSF DOI deposit of the frozen spec bundle is the primary, mandatory anchor** (tamper-evident); a public git-mirror tag is a best-effort supplement, not a substitute — cited in the article; an automated pre-run check asserts the working tree matches the frozen SHA before any timed run executes; **DoD commits to publish all results incl. losses** |
| Build-flag asymmetry (leafblower `-O3 -march=native` vs CRAN portable binaries) | **controlled, not just disclosed**: headline RQ1 uses the `portable` leafblower build matched to competitors' shipped flags; `native` reported as a separate tuned-build delta; both builds' flags in `environment.json` (§5) |
| Future-dated versions (survey 4.5, scipy 1.18…) | registry built against **installed** versions; lockfiles pin exact builds |
| Matrix runtime blowup | CORE/STRETCH split; reduced reps on 1.58M; `--smoke` mode |
| Competitor deps break DoD gate | strict `benchmarks/study/` isolation (§8); assert no DESCRIPTION/pyproject leakage |
| Subprocess RSS overhead | empty-process baseline subtracted |

---

## 11. Definition of Done

- [ ] **WU-0**: pre-materialised-CellTable entry point (R+Py) so cell-build is outside the timer for leafblower; parity-tested vs `harvest()` under the DoD gate.
- [ ] **Two leafblower builds** (`portable`, `native`); both flag-sets in `environment.json`; Python arm asserts venv `.so` (not `~/.local` shadow).
- [ ] Selected packages installed (R user-lib, Py venv); **registry built against installed versions**; **renv.lock + uv.lock + git SHA** captured; `environment.json` complete.
- [ ] Canonical + stepstone + synthetic-scaling problem specs emitted (data_ref, not inlined); scaling DGP frozen + certified non-home-turf (§3C).
- [ ] `common/metrics.{R,py}` reconcile **all three** existing impls to a **hand-computed golden**, then R↔Py rtol=1e-6; incl. **bound-violation** metric.
- [ ] Adapters for every competitor + leafblower method; contract incl. error_message / harmonized status (+bound_violation) / uniform-L∞ recomputed converged; validated on toy **and each package's home-turf golden** (fallback for packages w/o citable numbers: agreement vs `survey::calibrate` on the same problem within family-tol — applied **only to non-KL-family** competitors to avoid validating an adapter against a package that is itself in the KL comparison set; KL-family adapters without a citable golden fall back to the hand-computed toy golden instead).
- [ ] Hyperparameter pre-registration table + tol-mapping table committed.
- [ ] **spec + registry + metrics + hyperparams git-frozen + externally timestamp-anchored before any timed run; automated pre-run SHA check.**
- [ ] Harness runs CORE matrix; emits tidy results + per-(solver,problem,thread) weight store.
- [ ] Reports: work-precision (primary) + Dolan–Moré profiles (family-stratified, INFEAS=DNF) + quality tables + agreement heatmaps.
- [ ] Empty-process RSS baseline measured + subtracted.
- [ ] One-command replication verified from clean (lockfile restore).
- [ ] Competitor deps confirmed absent from DESCRIPTION/pyproject; DoD gate (R CMD INSTALL + testthat + pytest) still green.
- [ ] §9 SOTA gap section drafted; **full result set incl. leafblower losses committed.**

---

## 12. Flow

design-review-gate (this, round 2) → `planning-with-beads` (WU decomposition, one ticket each) → plan-review-gate → orchestrated execution (IMPLEMENT→VALIDATE→ADVERSARIAL REVIEW→COMMIT).

## 13. Open questions (deferred to planning, non-blocking)
- Container (Docker/Apptainer) in addition to lockfiles? (repro strength vs effort)
- Exact synthetic DGP for §3C scaling problem (margin structure, skew).

## 14. Change log

**v4→v5** (external Gemini 3.1 Pro review): §3C single-synthetic → **parametric instance family (30–100)** so Dolan–Moré profiles are statistically valid (§3C, §4, §7); **metric-formula bug fixes** — `marg_kl` no longer drops divergent (S→0,t>0) terms (reports +∞/sentinel, golden asserts), `weight_kl` uses per-problem `d_i` (not hardcoded 1), DEFF relabelled **Kish weighting DEFF/UWE** ≡ n/ESS (report ESS once) (§6); RQ5 agreement **anchored to an independent 1e-12 convex solve**, not in-house ORIS (§6); work-precision plotted on each family's **native divergence**, margin-L∞ only a neutral axis w/ Chebyshev-advantage caveat (§7); OT arm **relabelled redundant-IPF cross-check** + "no external K-margin OT competitor" promoted to §9 finding (§2, §9-4); §9-2 minimax reframed as **software-availability gap not novelty**; §9-3 TR-Newton novelty **restricted to smooth families** (not minimax); **downstream variance/replicate-weight cost promoted to a measured axis** (§9-6). Rejected Gemini's "survey::calibrate missing" — it is the first-listed Raking/GREG/Logit competitor (§2).

**v3→v4** (round 3 — 4/5 approved; Designer mechanical + Architect correctness): `build` variant added to weights filename + runs.parquet key (portable/native/na, no collision); runs.parquet full column+dtype list (§8); spec `data_ref` loader convention for file:/pkg:/gen:/inline (§3); OT extraction corrected — balanced row-sum degenerates, use left potential `u` or unbalanced OT, pinned against IPF golden (§2); one fresh subprocess per (solver,problem,thread,build) for correct VmHWM (§5); WU-0 precise cell-construction line boundary (§5); BLAS backend+version+dispatch captured (§5); §3C DGP externalised + finalized-before-planning; Zenodo DOI = primary freeze anchor (§10); survey golden-fallback restricted to non-KL family (§11).

**v2→v3** (round 2 blockers): WU-0 pre-materialised-CellTable entry point so leafblower's cell-build is outside the timer (§5, §11) + fallback; weights store keyed `(solver,problem,thread)` with explicit filename, `run_id`/`weights_ref` grain resolved (§8); build-flag confound *controlled* via `portable`+`native` leafblower builds, headline uses `portable`, claims capped (§5, §7, §10); bound-violation metric + `bound_violation` status state backing §9-finding-5 (§6, §8); external timestamp-anchored freeze + pre-run SHA check (§10, §11); three metric impls reconciled to golden (§6); venv-`.so` precondition + capture `CXXFLAGS`/`-march` (§5, §8); per-language RSS baseline (§5); OT coupling→weight extraction pinned (§2); numeric CORE budget ~2h + reps=5 on 1.58M (§4); uniform-L∞ recomputed converged + status dtype + wall_time unit (§8); metrics denormalised on runs.parquet (§8); home-turf golden fallback (§11); §3C DGP frozen + non-home-turf certification; §9 claims gated on results; iterations ref fixed to `types.hpp:89`; RQ1/RQ3 panel ownership (§7); `--smoke` sized (§8).

**v1→v2** (round 1 blockers):
Added §0 Alternatives; reuse-first §8 (refactor fit_metrics/compute_metrics, extend run()); OT restricted to K≤2 + ott excluded from ranked timing (§2); spec data_ref (§3); weights store (§8); adapter contract +error_message/status-enum/iterations/recomputed-converged (§8); work-precision primary + tol-mapping + hyperparameter pre-registration + min-duration rule (§5); INFEAS=DNF scoring + family-stratified profiles (§7); CORE/STRETCH matrix + cost estimate (§4); independent synthetic scaling problem + stepstone reframed (§3C); renv/uv lockfiles + git SHA (§8/§11); hand-computed golden + rtol relaxed to 1e-6 (§6); RQ6 = narrative deliverable (§1/§9); DoD-gate isolation (§8/§11); RQ5 pre-registered thresholds (§1); freeze + publish-losses (§10/§11); canary competitor (§2).
