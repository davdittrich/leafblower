# Graph Report - leafblower  (2026-04-30)

## Corpus Check
- 86 files · ~413,267 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 302 nodes · 384 edges · 25 communities detected
- Extraction: 94% EXTRACTED · 6% INFERRED · 0% AMBIGUOUS · INFERRED: 22 edges (avg confidence: 0.8)
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
- [[_COMMUNITY_Community 21|Community 21]]
- [[_COMMUNITY_Community 22|Community 22]]
- [[_COMMUNITY_Community 23|Community 23]]
- [[_COMMUNITY_Community 24|Community 24]]
- [[_COMMUNITY_Community 25|Community 25]]

## God Nodes (most connected - your core abstractions)
1. `iEPPA Algorithm (Iterative EPPA)` - 22 edges
2. `src/ieppa.cpp — iEPPA solver implementation` - 15 edges
3. `harvest()` - 14 edges
4. `method='ieppa' — paper-faithful algBCD solver` - 14 edges
5. `convergence API (metric + rule + tol)` - 14 edges
6. `harvest() calibration entry point` - 12 edges
7. `Calibration Refactor (remove AUTO, iEPPA default, ALM)` - 9 edges
8. `Raking Calibration Algorithm` - 8 edges
9. `Raking documentation: survey calibration algorithms` - 7 edges
10. `Calibration Solvers Redesign spec (2026-04-25)` - 7 edges

## Surprising Connections (you probably didn't know these)
- `method='ieppa' — paper-faithful algBCD solver` --validates--> `Stepstone full dataset benchmark (n=1.58M, K=9, M_cell/n=0.018)`  [INFERRED]
  NEWS.md → docs/superpowers/specs/2026-04-23-ieppa-convergence-hardening-design.md
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
Cohesion: 0.06
Nodes (36): Assumption 1: A^(i)_j binary entries, non-overlapping Hadamard patterns — enables closed-form BCD subproblem solutions, Plan: Bounded Convergence Fix (Dykstra + TDD), Bregman distance with Boltzmann-Shannon entropy kernel phi(X), Bug: bcd_sweep clamp cycling near constraint boundaries (~2.3e-3 residual), Bug: n_bounds_clamped undercount (leafblower-kssd), Capacity Constrained Multi-Marginal Optimal Transport (CMOT), Theorem 1: iEPPA convergence — X^k and X-tilde^k converge to same optimal solution under summable parameter sequences, Theorem 2: Dual BCD R-linear convergence to optimal solution of dual problem (3.4) (+28 more)

### Community 1 - "Community 1"
Cohesion: 0.09
Nodes (32): Cell table (sort-based dedup, cell-compressed representation), Chebyshev calibration (true L-inf minimum), GRAKE calibration (normalized Chebyshev), GREG calibration (Deville-Sarnal 1992, Newton QP for chi2), iEPPA (Sinkhorn-BCD / algBCD calibration), Modified LDLT factorization (Gill-Murray diagonal perturbation), Logarithmic barrier method (central-path, LP for L-inf), Normal equations (compute_normal_equations) (+24 more)

### Community 2 - "Community 2"
Cohesion: 0.09
Nodes (30): Dykstra alternating projections, Bug: iEPPA RK_ERR_INFEAS never returned (pre-clean-code-4), Bug: W_best snapshot incorrect across homotopy levels, Bug: wrapper normalization after solver breaks bounds_mode=unit, calib_dispatch.hpp — shared metric/rule dispatch header, Clean Code Fixes Round 2 Plan (2026-04-18), Clean Code Fixes Round 3 Plan (2026-04-19), convergence metric: grake_norm (survey::calibrate-compatible) (+22 more)

### Community 3 - "Community 3"
Cohesion: 0.11
Nodes (24): diagnose_weights(), harvest(), _parse_convergence(), _parse_sor(), Calibrate survey weights. Drop-in for R leafblower::harvest().      Parameters, # NOTE: No post-normalization clamp to [min_weight, max_weight]. Clamping, Derive pct_tol, absolute_tol, metric, rule, stop_when from convergence dict., Diagnose calibration quality (Python equivalent of R diagnose_weights()).      P (+16 more)

### Community 4 - "Community 4"
Cohesion: 0.08
Nodes (26): Augmented Lagrangian Method (ALM for sum(w)=n in L-BFGS-B), L-BFGS-B with logit/exp dual link, GP posterior mean (log10_compression=0): contour map of log(t_iEPPA/t_raking) over complexity x tol space; raking faster at high complexity, K-stability plot: 1.2x contour at K=3,9,18 showing robustness of GP algorithm selection boundary across GP ensemble size, GP posterior uncertainty: high posterior sd (orange/yellow) indicates unreliable algorithm selection regions in complexity x tol space, Bug: int overflow in M_cell*10 > n*9 routing, Bug: raking post-clamp normalization violates bounds, Plan: Calibration Solvers Plan B (Raking Cell-Table Migration) (+18 more)

### Community 5 - "Community 5"
Cohesion: 0.11
Nodes (24): algBCD at C=0 (Sinkhorn Block Coordinate Descent), Algorithm Selection Benchmark Design (2026-04-20), Bayesian Level Set Estimation benchmark (GP surrogate), Best-iterate tracking (W at minimum observed errRp), Capacity-Constrained Optimal Transport (CCOT), Cell compression (unique tuple optimization), CellTable data structure (sort-based cell deduplication), kComplexityThreshold constant for auto routing (+16 more)

### Community 6 - "Community 6"
Cohesion: 0.1
Nodes (20): WU-3: Adaptive damping (geometric blend, auto-trigger), Bug: bcd_sweep clamps during IPF violating Sinkhorn invariant, Bounded Convergence Fix Design (2026-04-18), Bregman Dykstra (multiplicative KL-space Dykstra), iEPPA Convergence Hardening Design (2026-04-23), Plan: Descent Monitor + kk1.20.4 Investigation (WU-P1, WU-P2), Dykstra's alternating projections algorithm, WU-1: False-positive infeasibility detection (persistent tracker) (+12 more)

### Community 7 - "Community 7"
Cohesion: 0.12
Nodes (17): APVA (Asymmetric Partial-Variable Anderson) — failed approach, bounds_mode parameter ('cell' | 'unit'), Cell-level bounding (vs obs-level bounding), Dual-domain (linear/log) execution, iEPPA code assessment: deviations from standard raking, iEPPA Speed/Convergence/Bounds Hardening Design (2026-04-24), P1.1 — Fuse post-sweep X_tilde + capacity inline (speed), P2.1 — Adaptive damping schedule (replace hard 0.5 latch) (+9 more)

### Community 8 - "Community 8"
Cohesion: 0.14
Nodes (15): Algorithm Selection Benchmark (Bayesian LSE), Plan: Algorithm Selection Benchmark (Bayesian LSE, K-stability), Plan: Fix best_iter to Track Active Metric, Best-Iterate Tracking (best_weights/best_error/best_iter), Bug: best_iter always tracks errRp regardless of active metric, Convergence Criteria Reform (pct/max_err/mean_err/kl/chi2), Plan: Convergence Reform + SOR + Best-Iterate (6 WUs), Fix: best_errRp_seen → best_metric_seen via select_metric() (+7 more)

### Community 9 - "Community 9"
Cohesion: 0.19
Nodes (13): Anderson acceleration (APVA, m=5, via dgels LAPACK), Halpern mixing (O(1/k) fixed-point acceleration), Intra-cell water-filling (bounds_mode=unit), Benchmark: kk1204 (n=1M, K=20, cat=5, dense compression), Benchmark: stepstone fulldata (n=1.58M, K=9, 836 cats), Bug: false-positive RK_ERR_INFEAS on transient empty buckets, Parameter: bounds_mode (unit vs cell), iEPPA Convergence Hardening (WU-1/2/3) (+5 more)

### Community 10 - "Community 10"
Cohesion: 0.43
Nodes (7): bench_harvest(), design_effect(), effective_sample_size(), load_data(), main(), Time leafblower.harvest() over n_runs; return timing + result stats., Kish (1965) design effect: n * sum(w²) / sum(w)².

### Community 11 - "Community 11"
Cohesion: 0.43
Nodes (8): calib_linalg.hpp — shared normal equations kernel, Calibration Solvers Redesign spec (2026-04-25), Calibration Solvers Plan D: calib_linalg + method=greg (2026-04-26), method='chebyshev' — LP for L∞ calibration, method='grake' — LP for grake_norm (survey::calibrate equivalent), method='greg' — Newton QP for chi2 (GREG calibration), method='raking' migration to cell-table (100-300x speedup), Rationale: migrate raking to cell-table for 100-300x speedup

### Community 12 - "Community 12"
Cohesion: 0.25
Nodes (8): Discrete 2D tomography application: image reconstruction from N projections using iEPPA+BCD, Ground-truth image: flower (rose with stem and leaves, 256x256), used in iEPPA tomography demo, Ground-truth image: tree (pine/fir silhouette, 256x256), used in iEPPA tomography demo, Hcap1: Linear operator A^(i)(X) for direction (1,0) — row-sum projection of matrix X onto marginal vector, Hcap2: Linear operator A^(i)(X) for direction (2,1) — diagonal-sum projection of matrix X onto marginal vector, Figure 4: PSNR vs N projections (N=10..90) for 5 ground-truth images; PSNR increases monotonically with N for all images, PSNR (Peak Signal-to-Noise Ratio): reconstruction quality metric for tomography experiments, Figure 5: Tomography reconstructions for flower/tree/animals/brain/lung at N=20,50,80 vs ground-truth

### Community 13 - "Community 13"
Cohesion: 0.53
Nodes (5): compute_metrics(), load_data(), main(), ipfn.IPFN in DataFrame mode on compressed cell table (28,905 cells)., run_ipfn()

### Community 14 - "Community 14"
Cohesion: 0.67
Nodes (5): bench_harvest(), design_effect(), effective_sample_size(), load_data(), main()

### Community 15 - "Community 15"
Cohesion: 0.4
Nodes (5): Plan: Code Review Fixes — ieppa.cpp + configure (6 tasks), Fix: configure OMP simd probe + run-gate + flag split + Makevars.in, Fix: ieppa.cpp convergence gate + W-sum helper + std::clamp, Refactor: Move weight normalization from wrappers into solvers, Rationale: Move normalization to solvers so wrappers are thin and output is consistent

### Community 16 - "Community 16"
Cohesion: 0.67
Nodes (2): patch_wolfe(), patch_wolfe_line_search()

### Community 17 - "Community 17"
Cohesion: 0.5
Nodes (4): bd (beads) issue tracker, graphify knowledge graph tool, leafblower R/Python package, OpenWolf context management system

### Community 18 - "Community 18"
Cohesion: 0.5
Nodes (4): Convergence criterion: improvement (relative errRp decrease), Convergence criterion: pct (weight-change stopping), Improvement-Based Convergence Criterion, Rationale: improvement criterion measures actual calibration progress; pct (weight change) can misfire when weights stall but error remains high

### Community 19 - "Community 19"
Cohesion: 0.67
Nodes (1): Ensure the installed wheel's leafblower package is found, not the local source t

### Community 21 - "Community 21"
Cohesion: 1.0
Nodes (1): Diagnose calibration quality (Python equivalent of R diagnose_weights()).      P

### Community 22 - "Community 22"
Cohesion: 1.0
Nodes (1): # NOTE: No post-normalization clamp to [min_weight, max_weight]. Clamping

### Community 23 - "Community 23"
Cohesion: 1.0
Nodes (1): src/leafblower.h — public C header + enum definitions

### Community 24 - "Community 24"
Cohesion: 1.0
Nodes (1): src/r_bridge.cpp

### Community 25 - "Community 25"
Cohesion: 1.0
Nodes (1): Parameter: method (ieppa/raking/lbfgsb/sinkhorn/greg/chebyshev/grake)

## Knowledge Gaps
- **120 isolated node(s):** `Time leafblower.harvest() over n_runs; return timing + result stats.`, `Kish (1965) design effect: n * sum(w²) / sum(w)².`, `ipfn.IPFN in DataFrame mode on compressed cell table (28,905 cells).`, `leafblower: high-performance survey calibration.`, `weights_out must be a copy, not a view into input.` (+115 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **Thin community `Community 16`** (4 nodes): `patch_wolfe.py`, `patch_wolfe()`, `patch_wolfe_line_search()`, `patch_wolfe.py`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 19`** (3 nodes): `conftest.py`, `Ensure the installed wheel's leafblower package is found, not the local source t`, `conftest.py`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 21`** (1 nodes): `Diagnose calibration quality (Python equivalent of R diagnose_weights()).      P`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 22`** (1 nodes): `# NOTE: No post-normalization clamp to [min_weight, max_weight]. Clamping`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 23`** (1 nodes): `src/leafblower.h — public C header + enum definitions`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 24`** (1 nodes): `src/r_bridge.cpp`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 25`** (1 nodes): `Parameter: method (ieppa/raking/lbfgsb/sinkhorn/greg/chebyshev/grake)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `harvest() calibration entry point` connect `Community 4` to `Community 0`, `Community 1`, `Community 2`, `Community 5`, `Community 6`, `Community 8`?**
  _High betweenness centrality (0.273) - this node is a cross-community bridge._
- **Why does `iEPPA Algorithm (Iterative EPPA)` connect `Community 0` to `Community 8`, `Community 4`, `Community 12`?**
  _High betweenness centrality (0.207) - this node is a cross-community bridge._
- **Why does `method='ieppa' — paper-faithful algBCD solver` connect `Community 5` to `Community 2`, `Community 4`, `Community 6`?**
  _High betweenness centrality (0.203) - this node is a cross-community bridge._
- **Are the 10 inferred relationships involving `harvest()` (e.g. with `bench_harvest()` and `main()`) actually correct?**
  _`harvest()` has 10 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Time leafblower.harvest() over n_runs; return timing + result stats.`, `Kish (1965) design effect: n * sum(w²) / sum(w)².`, `ipfn.IPFN in DataFrame mode on compressed cell table (28,905 cells).` to the rest of the system?**
  _120 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Community 0` be split into smaller, more focused modules?**
  _Cohesion score 0.06 - nodes in this community are weakly interconnected._
- **Should `Community 1` be split into smaller, more focused modules?**
  _Cohesion score 0.09 - nodes in this community are weakly interconnected._