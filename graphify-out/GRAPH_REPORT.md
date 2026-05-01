# Graph Report - leafblower  (2026-05-01)

## Corpus Check
- 86 files · ~461,708 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 459 nodes · 610 edges · 36 communities detected
- Extraction: 91% EXTRACTED · 9% INFERRED · 0% AMBIGUOUS · INFERRED: 55 edges (avg confidence: 0.8)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Community 0|Community 0]]
- [[_COMMUNITY_Community 1|Community 1]]
- [[_COMMUNITY_Community 2|Community 2]]
- [[_COMMUNITY_Community 3|Community 3]]
- [[_COMMUNITY_Community 4|Community 4]]
- [[_COMMUNITY_Community 5|Community 5]]
- [[_COMMUNITY_Community 6|Community 6]]
- [[_COMMUNITY_Community 7|Community 7]]
- [[_COMMUNITY_Community 8|Community 8]]
- [[_COMMUNITY_Community 9|Community 9]]
- [[_COMMUNITY_Community 10|Community 10]]
- [[_COMMUNITY_Community 11|Community 11]]
- [[_COMMUNITY_Community 12|Community 12]]
- [[_COMMUNITY_Community 13|Community 13]]
- [[_COMMUNITY_Community 14|Community 14]]
- [[_COMMUNITY_Community 15|Community 15]]
- [[_COMMUNITY_Community 16|Community 16]]
- [[_COMMUNITY_Community 17|Community 17]]
- [[_COMMUNITY_Community 18|Community 18]]
- [[_COMMUNITY_Community 19|Community 19]]
- [[_COMMUNITY_Community 20|Community 20]]
- [[_COMMUNITY_Community 21|Community 21]]
- [[_COMMUNITY_Community 57|Community 57]]
- [[_COMMUNITY_Community 62|Community 62]]
- [[_COMMUNITY_Community 63|Community 63]]
- [[_COMMUNITY_Community 64|Community 64]]
- [[_COMMUNITY_Community 65|Community 65]]
- [[_COMMUNITY_Community 66|Community 66]]
- [[_COMMUNITY_Community 67|Community 67]]
- [[_COMMUNITY_Community 68|Community 68]]
- [[_COMMUNITY_Community 69|Community 69]]
- [[_COMMUNITY_Community 70|Community 70]]
- [[_COMMUNITY_Community 71|Community 71]]
- [[_COMMUNITY_Community 72|Community 72]]
- [[_COMMUNITY_Community 73|Community 73]]
- [[_COMMUNITY_Community 74|Community 74]]

## God Nodes (most connected - your core abstractions)
1. `CellTable data structure (sort-based cell deduplication)` - 23 edges
2. `iEPPA Algorithm (Iterative EPPA)` - 22 edges
3. `src/ieppa.cpp — iEPPA solver implementation` - 21 edges
4. `src/lbfgsb_solver.cpp` - 18 edges
5. `build_cell_table()` - 15 edges
6. `harvest()` - 14 edges
7. `method='ieppa' — paper-faithful algBCD solver` - 14 edges
8. `convergence API (metric + rule + tol)` - 14 edges
9. `src/c_api.cpp — C API + algorithm selection` - 13 edges
10. `harvest() calibration entry point` - 12 edges

## Surprising Connections (you probably didn't know these)
- `ieppa_solve()` --calls--> `build_cell_table()`  [INFERRED]
  ieppa_92c4f45.cpp → src/cell_table.cpp
- `Raking Calibration Algorithm` --compares_ieppa_against--> `GP posterior mean log10_compression=0.0: iEPPA vs raking time ratio; raking wins (purple) at high complexity ~10^7`  [EXTRACTED]
  docs/superpowers/plans/2026-04-26-calibration-solvers-B-raking-cell-table.md → benchmarks/ieppa_vs_raking_3d_slice_0_0.pdf
- `Raking Calibration Algorithm` --compares_ieppa_against--> `GP posterior mean log10_compression=1.0: iEPPA vs raking ratio; raking advantage at high complexity persists with moderate compression`  [EXTRACTED]
  docs/superpowers/plans/2026-04-26-calibration-solvers-B-raking-cell-table.md → benchmarks/ieppa_vs_raking_3d_slice_1_0.pdf
- `Raking Calibration Algorithm` --compares_ieppa_against--> `GP posterior mean log10_compression=2.0: iEPPA vs raking ratio; crossover region shifts with higher compression`  [EXTRACTED]
  docs/superpowers/plans/2026-04-26-calibration-solvers-B-raking-cell-table.md → benchmarks/ieppa_vs_raking_3d_slice_2_0.pdf
- `Raking Calibration Algorithm` --compares_ieppa_against--> `GP posterior mean log10_compression=3.0: iEPPA vs raking ratio; at max compression raking dominates at complexity >10^6.5`  [EXTRACTED]
  docs/superpowers/plans/2026-04-26-calibration-solvers-B-raking-cell-table.md → benchmarks/ieppa_vs_raking_3d_slice_3_0.pdf

## Hyperedges (group relationships)
- **** — method_ieppa, method_raking, method_lbfgsb [INFERRED 1.00]
- **** — convergence_metric_max_err, convergence_metric_kl, convergence_metric_l1_weight, convergence_rule_improvement, convergence_rule_plateau [INFERRED 1.00]
- **** — cell_compression, prefactored_products, linear_space_dispatch [INFERRED 0.90]
- **Calibration solver suite: iEPPA + raking + sinkhorn + greg + chebyshev + grake share calib_linalg.cpp infrastructure** — algo_ieppa, algo_raking, algo_sinkhorn, algo_greg, algo_chebyshev, algo_grake, src_calib_linalg_cpp [INFERRED]
- **Convergence redesign: metric (6 values) x rule (3 values) orthogonal for all solvers** — convergence_redesign_metric_rule, convergence_metric_kl, convergence_metric_max_err, convergence_metric_chi2, convergence_criterion_improvement, src_types_hpp [INFERRED]
- **iEPPA acceleration: homotopy + greenkhorn + tang-eta compose to close 10x gap vs autumn** —  [EXTRACTED 1.00]
- **Convergence reform triad: pluggable criteria + SOR + best-iterate across all solvers** —  [EXTRACTED 1.00]
- **Three calibration solvers: iEPPA, raking, L-BFGS-B — all dispatched via harvest()** —  [EXTRACTED 1.00]
- **iEPPA+BCD framework for CMOT: iEPPA outer loop + dual BCD subsolver + CMOT problem structure together achieve provable convergence without small proximal parameter** —  [EXTRACTED 1.00]
- **GP-based algorithm selection system: GP emulator trained on iEPPA/raking/LBFGSB benchmark runs, uncertainty map, and K-stability validation jointly determine routing decision boundary** —  [INFERRED 0.80]
- **Tomography reconstruction evaluation: ground-truth images + PSNR metric + tomography projection operator together define the discrete tomography benchmark for iEPPA** —  [EXTRACTED 1.00]

## Communities

### Community 0 - "Community 0"
Cohesion: 0.05
Nodes (50): WU-3: Adaptive damping (geometric blend, auto-trigger), algBCD at C=0 (Sinkhorn Block Coordinate Descent), GP posterior mean (log10_compression=0): contour map of log(t_iEPPA/t_raking) over complexity x tol space; raking faster at high complexity, Algorithm Selection Benchmark Design (2026-04-20), K-stability plot: 1.2x contour at K=3,9,18 showing robustness of GP algorithm selection boundary across GP ensemble size, GP posterior uncertainty: high posterior sd (orange/yellow) indicates unreliable algorithm selection regions in complexity x tol space, Bayesian Level Set Estimation benchmark (GP surrogate), Best-iterate tracking (W at minimum observed errRp) (+42 more)

### Community 1 - "Community 1"
Cohesion: 0.05
Nodes (47): Algorithm Selection Benchmark (Bayesian LSE), Plan: Algorithm Selection Benchmark (Bayesian LSE, K-stability), Assumption 1: A^(i)_j binary entries, non-overlapping Hadamard patterns — enables closed-form BCD subproblem solutions, Plan: Fix best_iter to Track Active Metric, Best-Iterate Tracking (best_weights/best_error/best_iter), Plan: Bounded Convergence Fix (Dykstra + TDD), Bregman distance with Boltzmann-Shannon entropy kernel phi(X), Bug: bcd_sweep clamp cycling near constraint boundaries (~2.3e-3 residual) (+39 more)

### Community 2 - "Community 2"
Cohesion: 0.1
Nodes (29): GREG calibration (Deville-Sarnal 1992, Newton QP for chi2), Modified LDLT factorization (Gill-Murray diagonal perturbation), Normal equations (compute_normal_equations), Bug: sinkhorn Dykstra correction a[c] unbounded accumulation causing exp overflow, Plan D: calib_linalg + method=greg (Newton QP for chi2), Sinkhorn Correctness Fixes (Dykstra overflow, bisection), apply_obs_expansion(), compute_normal_equations() (+21 more)

### Community 3 - "Community 3"
Cohesion: 0.07
Nodes (35): Dykstra alternating projections, Bug: iEPPA RK_ERR_INFEAS never returned (pre-clean-code-4), Bug: W_best snapshot incorrect across homotopy levels, Bug: wrapper normalization after solver breaks bounds_mode=unit, calib_dispatch.hpp — shared metric/rule dispatch header, Clean Code Fixes Round 2 Plan (2026-04-18), Clean Code Fixes Round 3 Plan (2026-04-19), Plan: Code Review Fixes — ieppa.cpp + configure (6 tasks) (+27 more)

### Community 4 - "Community 4"
Cohesion: 0.09
Nodes (29): bench_harvest(), design_effect(), effective_sample_size(), load_data(), main(), diagnose_weights(), harvest(), _parse_convergence() (+21 more)

### Community 5 - "Community 5"
Cohesion: 0.12
Nodes (25): Bug: int overflow in M_cell*10 > n*9 routing, Bug: raking post-clamp normalization violates bounds, Clean Code Fixes Round 1 Plan (2026-04-18), Fix: Eliminate per-iteration vector copies in L-BFGS-B (std::swap), Plan: L-BFGS-B Audit Fixes (rev3, WU-1 to WU-3), Plan: L-BFGS-B Phase 2, L-BFGS-B Solver, PYBIND11_MODULE() (+17 more)

### Community 6 - "Community 6"
Cohesion: 0.09
Nodes (27): Augmented Lagrangian Method (ALM for sum(w)=n in L-BFGS-B), Cell table (sort-based dedup, cell-compressed representation), Chebyshev calibration (true L-inf minimum), GRAKE calibration (normalized Chebyshev), iEPPA (Sinkhorn-BCD / algBCD calibration), L-BFGS-B with logit/exp dual link, Logarithmic barrier method (central-path, LP for L-inf), Raking (multiplicative iterative proportional fitting) (+19 more)

### Community 7 - "Community 7"
Cohesion: 0.24
Nodes (19): leafblower C Performance Optimization (SIMD, libmvec), Post-P3 Code Review Fixes (2026-04-22), build_offsets(), compute_du(), compute_final_weights_and_error(), compute_targets_abs(), compute_u(), src/lbfgsb_solver.cpp (+11 more)

### Community 8 - "Community 8"
Cohesion: 0.12
Nodes (18): APVA (Asymmetric Partial-Variable Anderson) — failed approach, bounds_mode parameter ('cell' | 'unit'), Cell compression (unique tuple optimization), Cell-level bounding (vs obs-level bounding), Dual-domain (linear/log) execution, iEPPA code assessment: deviations from standard raking, iEPPA Speed/Convergence/Bounds Hardening Design (2026-04-24), P1.1 — Fuse post-sweep X_tilde + capacity inline (speed) (+10 more)

### Community 9 - "Community 9"
Cohesion: 0.15
Nodes (8): CellTable data structure (sort-based cell deduplication), bits_needed(), build_cell_table(), pack_key_compute(), pack_key_fits(), ieppa_solve(), src/cell_table.hpp + cell_table.cpp, chebyshev_solve()

### Community 10 - "Community 10"
Cohesion: 0.15
Nodes (15): Bug: bcd_sweep clamps during IPF violating Sinkhorn invariant, Bounded Convergence Fix Design (2026-04-18), Bregman Dykstra (multiplicative KL-space Dykstra), calib_linalg.hpp — shared normal equations kernel, Calibration Solvers Redesign spec (2026-04-25), Calibration Solvers Plan D: calib_linalg + method=greg (2026-04-26), Dykstra's alternating projections algorithm, Bug: L-BFGS-B uses exp link when max_weight finite (should use logit) (+7 more)

### Community 11 - "Community 11"
Cohesion: 0.21
Nodes (10): aggregate_to_margin(), apply_rule(), build_cat_offset(), check_convergence(), compute_cell_bounds(), compute_cell_metrics(), resolve_hi(), select_metric() (+2 more)

### Community 12 - "Community 12"
Cohesion: 0.19
Nodes (13): Anderson acceleration (APVA, m=5, via dgels LAPACK), Halpern mixing (O(1/k) fixed-point acceleration), Intra-cell water-filling (bounds_mode=unit), Benchmark: kk1204 (n=1M, K=20, cat=5, dense compression), Benchmark: stepstone fulldata (n=1.58M, K=9, 836 cats), Bug: false-positive RK_ERR_INFEAS on transient empty buckets, Parameter: bounds_mode (unit vs cell), iEPPA Convergence Hardening (WU-1/2/3) (+5 more)

### Community 13 - "Community 13"
Cohesion: 0.43
Nodes (7): bench_harvest(), design_effect(), effective_sample_size(), load_data(), main(), Time leafblower.harvest() over n_runs; return timing + result stats., Kish (1965) design effect: n * sum(w²) / sum(w)².

### Community 14 - "Community 14"
Cohesion: 0.25
Nodes (8): Discrete 2D tomography application: image reconstruction from N projections using iEPPA+BCD, Ground-truth image: flower (rose with stem and leaves, 256x256), used in iEPPA tomography demo, Ground-truth image: tree (pine/fir silhouette, 256x256), used in iEPPA tomography demo, Hcap1: Linear operator A^(i)(X) for direction (1,0) — row-sum projection of matrix X onto marginal vector, Hcap2: Linear operator A^(i)(X) for direction (2,1) — diagonal-sum projection of matrix X onto marginal vector, Figure 4: PSNR vs N projections (N=10..90) for 5 ground-truth images; PSNR increases monotonically with N for all images, PSNR (Peak Signal-to-Noise Ratio): reconstruction quality metric for tomography experiments, Figure 5: Tomography reconstructions for flower/tree/animals/brain/lung at N=20,50,80 vs ground-truth

### Community 15 - "Community 15"
Cohesion: 0.53
Nodes (5): compute_metrics(), load_data(), main(), ipfn.IPFN in DataFrame mode on compressed cell table (28,905 cells)., run_ipfn()

### Community 16 - "Community 16"
Cohesion: 0.33
Nodes (1): bulk_scaled_exp()

### Community 17 - "Community 17"
Cohesion: 0.67
Nodes (2): patch_wolfe(), patch_wolfe_line_search()

### Community 18 - "Community 18"
Cohesion: 0.5
Nodes (4): bd (beads) issue tracker, graphify knowledge graph tool, leafblower R/Python package, OpenWolf context management system

### Community 19 - "Community 19"
Cohesion: 0.5
Nodes (4): Convergence criterion: improvement (relative errRp decrease), Convergence criterion: pct (weight-change stopping), Improvement-Based Convergence Criterion, Rationale: improvement criterion measures actual calibration progress; pct (weight change) can misfire when weights stall but error remains high

### Community 20 - "Community 20"
Cohesion: 1.0
Nodes (1): Root conftest: remove python/ src tree from sys.path so the installed wheel is i

### Community 21 - "Community 21"
Cohesion: 1.0
Nodes (1): Ensure the installed wheel's leafblower package is found, not the local source t

### Community 57 - "Community 57"
Cohesion: 1.0
Nodes (1): src/leafblower.h — public C header + enum definitions

### Community 62 - "Community 62"
Cohesion: 1.0
Nodes (1): Ensure the installed wheel's leafblower package is found, not the local source t

### Community 63 - "Community 63"
Cohesion: 1.0
Nodes (1): Time leafblower.harvest() over n_runs; return timing + result stats.

### Community 64 - "Community 64"
Cohesion: 1.0
Nodes (1): Kish (1965) design effect: n * sum(w²) / sum(w)².

### Community 65 - "Community 65"
Cohesion: 1.0
Nodes (1): leafblower: high-performance survey calibration.

### Community 66 - "Community 66"
Cohesion: 1.0
Nodes (1): weights_out must be a copy, not a view into input.

### Community 67 - "Community 67"
Cohesion: 1.0
Nodes (1): Balanced 2-level fixture for parity tests.

### Community 68 - "Community 68"
Cohesion: 1.0
Nodes (1): Default convergence: metric=max_err, rule=improvement (WU-G).

### Community 69 - "Community 69"
Cohesion: 1.0
Nodes (1): l1_weight_change and grake_norm present; pct_change absent (WU-G).

### Community 70 - "Community 70"
Cohesion: 1.0
Nodes (1): All 6 quality metrics must be present in result attrs (WU-G).

### Community 71 - "Community 71"
Cohesion: 1.0
Nodes (1): best_error <= max_error always (best iterate tracked throughout).

### Community 72 - "Community 72"
Cohesion: 1.0
Nodes (1): Diagnose calibration quality (Python equivalent of R diagnose_weights()).      P

### Community 73 - "Community 73"
Cohesion: 1.0
Nodes (1): # NOTE: No post-normalization clamp to [min_weight, max_weight]. Clamping

### Community 74 - "Community 74"
Cohesion: 1.0
Nodes (1): Parameter: method (ieppa/raking/lbfgsb/sinkhorn/greg/chebyshev/grake)

## Knowledge Gaps
- **129 isolated node(s):** `Root conftest: remove python/ src tree from sys.path so the installed wheel is i`, `Time leafblower.harvest() over n_runs; return timing + result stats.`, `Kish (1965) design effect: n * sum(w²) / sum(w)².`, `ipfn.IPFN in DataFrame mode on compressed cell table (28,905 cells).`, `Ensure the installed wheel's leafblower package is found, not the local source t` (+124 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **Thin community `Community 16`** (6 nodes): `lbw_math.hpp`, `bulk_exp_clipped()`, `bulk_log()`, `bulk_scaled_exp()`, `bulk_scaled_log()`, `lbw_math.hpp`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 17`** (4 nodes): `patch_wolfe()`, `patch_wolfe_line_search()`, `patch_wolfe.py`, `patch_wolfe.py`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 20`** (2 nodes): `conftest.py`, `Root conftest: remove python/ src tree from sys.path so the installed wheel is i`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 21`** (2 nodes): `conftest.py`, `Ensure the installed wheel's leafblower package is found, not the local source t`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 57`** (1 nodes): `src/leafblower.h — public C header + enum definitions`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 62`** (1 nodes): `Ensure the installed wheel's leafblower package is found, not the local source t`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 63`** (1 nodes): `Time leafblower.harvest() over n_runs; return timing + result stats.`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 64`** (1 nodes): `Kish (1965) design effect: n * sum(w²) / sum(w)².`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 65`** (1 nodes): `leafblower: high-performance survey calibration.`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 66`** (1 nodes): `weights_out must be a copy, not a view into input.`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 67`** (1 nodes): `Balanced 2-level fixture for parity tests.`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 68`** (1 nodes): `Default convergence: metric=max_err, rule=improvement (WU-G).`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 69`** (1 nodes): `l1_weight_change and grake_norm present; pct_change absent (WU-G).`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 70`** (1 nodes): `All 6 quality metrics must be present in result attrs (WU-G).`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 71`** (1 nodes): `best_error <= max_error always (best iterate tracked throughout).`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 72`** (1 nodes): `Diagnose calibration quality (Python equivalent of R diagnose_weights()).      P`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 73`** (1 nodes): `# NOTE: No post-normalization clamp to [min_weight, max_weight]. Clamping`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 74`** (1 nodes): `Parameter: method (ieppa/raking/lbfgsb/sinkhorn/greg/chebyshev/grake)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `CellTable data structure (sort-based cell deduplication)` connect `Community 9` to `Community 0`, `Community 2`, `Community 3`, `Community 5`, `Community 8`, `Community 11`?**
  _High betweenness centrality (0.193) - this node is a cross-community bridge._
- **Why does `harvest() calibration entry point` connect `Community 0` to `Community 1`, `Community 3`, `Community 5`, `Community 6`?**
  _High betweenness centrality (0.154) - this node is a cross-community bridge._
- **Why does `method='ieppa' — paper-faithful algBCD solver` connect `Community 0` to `Community 8`, `Community 9`, `Community 3`, `Community 6`?**
  _High betweenness centrality (0.133) - this node is a cross-community bridge._
- **Are the 10 inferred relationships involving `build_cell_table()` (e.g. with `ieppa_solve()` and `ieppa_solve()`) actually correct?**
  _`build_cell_table()` has 10 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Root conftest: remove python/ src tree from sys.path so the installed wheel is i`, `Time leafblower.harvest() over n_runs; return timing + result stats.`, `Kish (1965) design effect: n * sum(w²) / sum(w)².` to the rest of the system?**
  _129 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.05 - nodes in this community are weakly interconnected._
- **Should `Community 1` be split into smaller, more focused modules?**
  _Cohesion score 0.05 - nodes in this community are weakly interconnected._