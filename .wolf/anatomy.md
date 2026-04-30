# anatomy.md

> Auto-maintained by OpenWolf. Last scanned: 2026-04-30T08:26:26.443Z
> Files: 511 tracked | Anatomy hits: 0 | Misses: 0

## ..Rcheck/

- `00check.log` — Declares calls (~2235 tok)
- `00install.out` (~251 tok)
- `leafblower-manual.log` (~4595 tok)
- `leafblower-manual.tex` — algorithm: proportions, proportions (~3688 tok)
- `Rdlatex.log` (~5635 tok)

## ..Rcheck/R_check_bin/

- `R` (~22 tok)
- `Rscript` (~24 tok)

## ..Rcheck/leafblower/

- `DESCRIPTION` (~247 tok)
- `INDEX` (~139 tok)
- `LICENSE` — Project license (~16 tok)
- `NAMESPACE` (~64 tok)
- `NEWS.md` — leafblower (development) (~224 tok)

## ..Rcheck/leafblower/Meta/

- `features.rds` (~32 tok)
- `hsearch.rds` (~134 tok)
- `links.rds` (~60 tok)
- `nsInfo.rds` (~98 tok)
- `package.rds` (~258 tok)
- `Rd.rds` (~132 tok)

## ..Rcheck/leafblower/R/

- `leafblower` — File share/R/nspackloader.R (~283 tok)
- `leafblower.rdb` (~3933 tok)
- `leafblower.rdx` (~117 tok)

## ..Rcheck/leafblower/help/

- `aliases.rds` (~44 tok)
- `AnIndex` (~54 tok)
- `leafblower.rdb` (~3126 tok)
- `leafblower.rdx` (~87 tok)
- `paths.rds` (~55 tok)

## ..Rcheck/leafblower/html/

- `00Index.html` — R: High-Performance Survey Calibration via iEPPA and L-BFGS-B (~535 tok)
- `R.css` — Styles: 3 rules, 1 media queries (~527 tok)

## ..Rcheck/leafblower/libs/

- `symbols.rds` (~516 tok)

## ..Rcheck/tests/

- `startup.Rs` — # A custom startup file for tests (~38 tok)
- `testthat.R` — This file is part of the R package leafblower. (~49 tok)
- `testthat.Rout` — Declares to (~259 tok)
- `testthat.Rout.fail` — Declares to (~870 tok)

## ..Rcheck/tests/testthat/

- `lbfgsb_baseline_time.rds` (~27 tok)
- `lbfgsb_ref_weights.rds` (~2943 tok)
- `task1_ref.rds` (~531 tok)
- `task2_ieppa_ref.rds` (~1526 tok)
- `test-algo-selection.R` — tests/testthat/test-algo-selection.R (~2289 tok)
- `test-bench-gate.R` (~609 tok)
- `test-best-iterate.R` (~1633 tok)
- `test-bounded-convergence.R` (~400 tok)
- `test-calib-linalg.R` — Direct tests for LDLT factorization and compute_normal_equations (~1442 tok)
- `test-calibration-solvers.R` (~8311 tok)
- `test-cell-table.R` (~497 tok)
- `test-compare.R` (~393 tok)
- `test-compat.R` (~292 tok)
- `test-config-defaults.R` (~600 tok)
- `test-convergence-criteria.R` (~4070 tok)
- `test-convergence-trajectory.R` (~253 tok)
- `test-design.R` (~223 tok)
- `test-eta-schedule.R` (~402 tok)
- `test-harvest.R` (~1083 tok)
- `test-homotopy.R` (~533 tok)
- `test-ieppa-bounds-mode.R` (~2360 tok)
- `test-ieppa-faithful.R` (~1057 tok)
- `test-ieppa-nonuniform-d.R` (~447 tok)
- `test-ieppa-persistent-infeas.R` — WU-1: persistent-infeas tracker regression test. (~1171 tok)
- `test-ieppa.R` (~611 tok)
- `test-lbfgsb.R` (~957 tok)
- `test-logit.R` (~308 tok)
- `test-priority-sweep.R` (~783 tok)
- `test-quality-metrics.R` (~1266 tok)
- `test-raking.R` (~604 tok)
- `test-rk-params-passthrough.R` — Regression test: rk_params_init not called on caller-supplied params (~300 tok)
- `test-sor.R` (~459 tok)
- `testthat-problems.rds` (~4935 tok)

## ..Rcheck/tests/testthat/_problems/

- `test-ieppa-nonuniform-d-28.R` — Extracted from test-ieppa-nonuniform-d.R:28 (~264 tok)
- `test-ieppa-nonuniform-d-29.R` — Extracted from test-ieppa-nonuniform-d.R:29 (~272 tok)
- `test-raking-89.R` — Extracted from test-raking.R:89 (~295 tok)

## ..Rcheck/tests/testthat/fixtures/

- `ieppa_kl_reference_stepstone.rds` (~47 tok)
- `raking_obs_reference_stepstone.rds` (~58 tok)
- `raking_squarem_baseline.rds` (~353 tok)
- `stepstone_best_error_ref.rds` (~15 tok)
- `stepstone_reference_diag.R` — Diagnostic report on the captured stepstone_reference.rds. (~1598 tok)
- `stepstone_reference_diag.txt` (~828 tok)
- `stepstone_reference_run.R` — Reference stepstone run matching the ORIGINAL Rmd (~2101 tok)
- `stepstone_reference_summary.rds` (~218 tok)
- `stepstone_small_targets.rds` (~1737 tok)
- `stepstone_small.parquet` (~24475 tok)
- `stepstone_verify.R` — Verification-only: report n + per-margin counts from the saved parquet. (~1151 tok)

## ./

- `.gitignore` — Git ignore rules (~113 tok)
- `.Rbuildignore` (~75 tok)
- `.tldrignore` — TLDR ignore patterns (gitignore syntax) (~280 tok)
- `AGENTS.md` — Agent Instructions (~726 tok)
- `baseline_bench.R` (~191 tok)
- `cell_table_92c4f45.cpp` — include "lbw_config.h" (~1375 tok)
- `cell_table_92c4f45.hpp` — pragma once (~285 tok)
- `CLAUDE.md` — OpenWolf (~729 tok)
- `cleanup` (~8 tok)
- `code-review-findings.md` — Code Review: C++ Source Files (cell_table, chebyshev, sinkhorn, grake, greg) (~1627 tok)
- `configure` — Detect C++17 support. (~1050 tok)
- `conftest.py` — Root conftest: remove python/ src tree from sys.path so the installed wheel (~89 tok)
- `DESCRIPTION` (~228 tok)
- `GEMINI.md` — graphify (~165 tok)
- `ieppa_92c4f45.cpp` — include "lbw_config.h" (~3859 tok)
- `LICENSE` — Project license (~16 tok)
- `NAMESPACE` (~71 tok)
- `NEWS.md` — leafblower (development version) (~1013 tok)
- `patch_raking.py` — Declares CalibState (~4573 tok)
- `patch_wolfe.py` — patch_wolfe, patch_wolfe_line_search (~903 tok)
- `REVIEW_FINDINGS.md` — Code Review: harvest.R & calib_dispatch.hpp (~775 tok)
- `security-review-2026-04-27-ieppa.json` (~3036 tok)
- `test_output.log` (~2622 tok)

## .beads/

- `.gitignore` — Git ignore rules (~444 tok)
- `.local_version` (~2 tok)
- `config.yaml` — Beads Configuration File (~597 tok)
- `export-state.json` (~38 tok)
- `interactions.jsonl` (~44963 tok)
- `issues.jsonl` (~216971 tok)
- `last-touched` (~5 tok)
- `metadata.json` (~46 tok)
- `README.md` — Project documentation (~562 tok)

## .beads/embeddeddolt/

- `.lock` (~0 tok)

## .beads/embeddeddolt/leafblower/.dolt/

- `config.json` (~1 tok)
- `repo_state.json` (~24 tok)

## .beads/embeddeddolt/leafblower/.dolt/noms/

- `LOCK` (~0 tok)
- `manifest` (~40 tok)

## .beads/embeddeddolt/leafblower/.dolt/temptf/

- `dolt_embedded_metrics` (~0 tok)

## .beads/hooks/

- `post-checkout` — --- BEGIN BEADS INTEGRATION v1.0.0 --- (~795 tok)
- `post-commit` — graphify-hook-start (~625 tok)
- `post-merge` — --- BEGIN BEADS INTEGRATION v1.0.0 --- (~210 tok)
- `pre-commit` — --- BEGIN BEADS INTEGRATION v1.0.0 --- (~210 tok)
- `pre-push` — --- BEGIN BEADS INTEGRATION v1.0.0 --- (~208 tok)
- `prepare-commit-msg` — --- BEGIN BEADS INTEGRATION v1.0.0 --- (~218 tok)

## .beads/plans/

- `active-plan.md` — iEPPA Overflow Fix (T1.B) Implementation Plan (~6198 tok)

## .claude/

- `settings.json` (~661 tok)
- `settings.local.json` (~1027 tok)

## .claude/rules/

- `openwolf.md` (~313 tok)

## .gemini-bridge/

- `feedback.log` (~756 tok)

## .gemini/

- `settings.json` (~137 tok)

## .tldr/

- `languages.json` (~35 tok)

## .tldr/cache/

- `call_graph.json` (~928 tok)

## benchmarks/

- `ABE_trajectory_stepstone_small.csv` (~40 tok)
- `algo_selection_benchmark.R` — benchmarks/algo_selection_benchmark.R (~4837 tok)
- `algo_selection_results.rds` (~533 tok)
- `autumn_nr_benchmark.R` — benchmarks/autumn_nr_benchmark.R (~1030 tok)
- `baseline_trajectory_stepstone_small.csv` (~39 tok)
- `baseline_tuning_sweep.R` (~439 tok)
- `baseline_tuning_sweep.rds` (~74 tok)
- `benchmark_run.log` (~229 tok)
- `compute_rate_slope.R` — benchmarks/compute_rate_slope.R (~516 tok)
- `final_rate_slope.rds` (~369 tok)
- `ieppa_vs_raking_bench.R` — Bayesian Level Set Estimation: ieppa (faithful) vs raking (hybrid). (~1352 tok)
- `ieppa_vs_raking_results.rds` (~6523 tok)
- `make_stepstone_small_fixture.R` — Prefer fulldata if available (main repo); fall back to standard bench data in worktree. (~490 tok)
- `plot_helpers.R` — benchmarks/plot_helpers.R (~1734 tok)
- `probe_baseline.R` (~149 tok)
- `python_ipf_benchmark.py` — load_data, compute_metrics, run_ipfn, main (~1338 tok)
- `sraa-m-baseline-pre.log` (~258 tok)
- `stepstone_all_methods.R` — All-methods comparison on stepstone fulldata (1.58M rows, 9 margins). (~1566 tok)
- `stepstone_bench_targets.json` (~10918 tok)
- `stepstone_benchmark.py` — def: load_data, bench_harvest, design_effect, effective_sample_size + 1 more (~1398 tok)
- `stepstone_benchmark.R` — benchmarks/stepstone_benchmark.R (~2830 tok)
- `stepstone_compare_current.R` — 3-way stepstone comparison: autumn (ref) vs leafblower cell-mode (ref+current) (~1751 tok)
- `stepstone_fulldata_bench_targets.json` (~10952 tok)
- `stepstone_fulldata_benchmark.py` — load_data, bench_harvest, design_effect, effective_sample_size + 1 more (~1487 tok)
- `stepstone_fulldata_benchmark.R` — benchmarks/stepstone_fulldata_benchmark.R (~4895 tok)
- `stepstone_fulldata_homotopy_report.rds` (~96 tok)
- `stepstone_fulldata_homotopy_run_v2.log` (~666 tok)
- `stepstone_fulldata_homotopy.R` — benchmarks/stepstone_fulldata_homotopy.R (~1974 tok)

## data-raw/

- `gen_ieppa_kl_ref.R` — Generates tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds (~376 tok)
- `gen_ieppa_pre_alm_ref.R` (~199 tok)
- `gen_raking_obs_ref.R` — Generates tests/testthat/fixtures/raking_obs_reference_stepstone.rds (~342 tok)

## docs/

- `ieppa_assessmen.md` — To guarantee that **individual observation weights** never exceed your `min_weight` ($L$) and `max_weight` ($U$) bounds, you must address the core ... (~2824 tok)
- `raking_bounds.md` — Declares changes (~2256 tok)
- `raking.md` — Declares is (~10956 tok)

## docs/iEPPA/arxiv.2011.14312/

- `iEPPA-LP-R2-arXiv.tex` — Declares of (~50652 tok)

## docs/iEPPA/arxiv.2011.14312/references/

- `Ref_iEPPA.bib` — Declares measures (~7454 tok)

## docs/iEPPA/code/2Dtomography/

- `master.m` — Declares of (~728 tok)

## docs/iEPPA/code/2Dtomography/main/

- `EPPA.m` (~2376 tok)
- `Gprocedure.m` — Declares goes (~502 tok)

## docs/iEPPA/code/2Dtomography/tools/

- `cellOnes1D.m` (~60 tok)
- `genCollectionA.m` — Declares goes (~1373 tok)
- `genData2D.m` (~386 tok)
- `psnr2.m` — Declares goes (~56 tok)
- `sumcell.m` — Declares goes (~40 tok)
- `tableplot.m` — Declares goes (~75 tok)

## docs/iEPPA/code/ExpSynthesisData/

- `run_cot_2d.m` (~1100 tok)
- `run_cot_3d.m` (~892 tok)

## docs/iEPPA/code/ExpSynthesisData/solvers/

- `bcd_2d.m` (~370 tok)
- `bcd_3d.m` (~530 tok)
- `bregdist_2d.m` — Declares bd (~164 tok)
- `bregdist_3d.m` — Declares bd (~176 tok)
- `dykstra.m` (~656 tok)
- `ieppa_2d.m` (~1188 tok)
- `ieppa_3d.m` (~1321 tok)
- `ot_rounding_2d.m` — Declares Q (~185 tok)
- `ot_rounding_3d.m` — Declares Q (~220 tok)
- `stabilized_dykstra.m` — Declares value (~1168 tok)

## docs/investigations/

- `2026-04-23-kk1204-convergence.md` — kk1.20.4 Convergence + Wall-Clock Investigation (~2163 tok)
- `2026-04-24-ieppa-accel-research.md` — iEPPA-preserving convergence acceleration — literature scoping review (~3805 tok)

## docs/superpowers/plans/

- `.Rhistory` (~0 tok)
- `2026-04-18-bounded-convergence-fix.md` — Bounded Convergence Fix Implementation Plan (~4809 tok)
- `2026-04-18-clean-code-fixes-2.md` — Clean Code Fixes (Round 2) Implementation Plan (~6277 tok)
- `2026-04-18-clean-code-fixes.md` — Clean Code Fixes Implementation Plan (~4487 tok)
- `2026-04-18-leafblower-core.md` — Leafblower Core — Implementation Plan (~24756 tok)
- `2026-04-19-c-perf-optimization.md` — leafblower C Performance Optimization — Implementation Plan (~11563 tok)
- `2026-04-19-clean-code-fixes-3.md` — Clean Code Fixes (Round 3) Implementation Plan (~8196 tok)
- `2026-04-19-clean-code-fixes-4.md` — Clean Code Fixes (Round 4) Implementation Plan (~2244 tok)
- `2026-04-20-algo-selection-benchmark.md` — Algorithm Selection Benchmark Implementation Plan (~11806 tok)
- `2026-04-20-code-review-fixes.md` — Code Review Fixes — ieppa.cpp + configure Implementation Plan (~4514 tok)
- `2026-04-22-calibration-refactor.md` — Calibration Refactor Implementation Plan (~9329 tok)
- `2026-04-22-review-fixes.md` — Plan: Post-P3 Code Review Fixes (~2661 tok)
- `2026-04-23-descent-monitor-and-kk1204-investigation.md` — Descent Monitor + kk1.20.4 Convergence Investigation (~5947 tok)
- `2026-04-23-ieppa-convergence-hardening-impl.md` — iEPPA Convergence Hardening Implementation Plan (~14719 tok)
- `2026-04-23-ieppa-faithful-impl.md` — Faithful iEPPA Implementation Plan (~16794 tok)
- `2026-04-23-ieppa-review-fixes.md` — iEPPA Review Fixes Implementation Plan (~6692 tok)
- `2026-04-23-lbfgsb-audit-fixes.md` — Plan: L-BFGS-B Audit Fixes (2026-04-23) — rev 3 (~3198 tok)
- `2026-04-23-lbfgsb-phase2.md` — Plan: L-BFGS-B Phase 2 — O(n) SIMD + Epic Hygiene + Stall Investigation (rev 3) (~4703 tok)
- `2026-04-24-convergence-metrics-sor-impl.md` — Convergence Reform + SOR + Best-Iterate Implementation Plan (~12656 tok)
- `2026-04-24-fix-max-err-of-rerun-benchmark.md` — Fix max_err_of() name-alignment bug + re-run WU-6 benchmark (~1502 tok)
- `2026-04-24-ieppa-halpern-impl.md` — P2.2d Halpern Mixing Implementation Plan (~5026 tok)
- `2026-04-24-ieppa-homotopy-greenkhorn.md` — iEPPA Homotopy + Greenkhorn + Tang Dynamic-η Implementation Plan (rev 2) (~14279 tok)
- `2026-04-24-ieppa-speed-convergence-bounds-impl.md` — iEPPA Speed / Convergence / Bounds Implementation Plan (~19124 tok)
- `2026-04-24-move-normalization-to-solvers.md` — Move Weight Normalization from Wrappers to Solvers (rev 3) (~3651 tok)
- `2026-04-24-n_bounds_clamped-counter-fix.md` — Fix ieppa n_bounds_clamped Counter Undercount (rev 2) (~2520 tok)
- `2026-04-24-review-fixes-kssd-normalize.md` — Address Critical-Code-Review Findings (commits bac2877 + d463f4b) — rev 5 (~4834 tok)
- `2026-04-25-best-iter-active-metric.md` — Fix best_iter to Track Active Metric (~2469 tok)
- `2026-04-25-calibration-solvers-A-infrastructure.md` — Calibration Solvers — Plan A: Infrastructure + iEPPA Default Change (~7421 tok)
- `2026-04-25-convergence-redesign-impl.md` — Convergence Criteria Redesign Implementation Plan (~14497 tok)
- `2026-04-25-convergence-reform-fixes.md` — Convergence Reform Post-Review Fixes Implementation Plan (~10027 tok)
- `2026-04-25-convergence-simplify.md` — Convergence Redesign — Simplify & Fix Plan (~4612 tok)
- `2026-04-25-improvement-criterion.md` — Improvement-Based Convergence Criterion Implementation Plan (~6449 tok)
- `2026-04-26-calibration-solvers-B-raking-cell-table.md` — Calibration Solvers — Plan B: Raking Cell-Table Migration (~6642 tok)
- `2026-04-26-calibration-solvers-C-sinkhorn.md` — Calibration Solvers — Plan C: method="sinkhorn" (KL Bregman Dykstra) (~5078 tok)
- `2026-04-26-calibration-solvers-D-greg.md` — Calibration Solvers — Plan D: calib_linalg + method="greg" (Newton QP for chi2) (~6556 tok)
- `2026-04-26-calibration-solvers-E-chebyshev-grake.md` — Calibration Solvers — Plan E: method="chebyshev" + method="grake" (Mehrotra-style LP IPM) (~7504 tok)
- `2026-04-26-chebyshev-nu-fix.md` — Fix chebyshev IPM: Add Normalization Dual ν (second Sherman-Morrison) (~2118 tok)
- `2026-04-26-chebyshev-nu-reference-elimination.md` — Chebyshev ν Fix: Reference Category Elimination — Implementation Plan (~6143 tok)
- `2026-04-26-chebyshev-simplify.md` — Chebyshev/Grake Simplification Fixes (~2182 tok)
- `2026-04-26-critical-review-fixes.md` — Critical Review Fixes (~3628 tok)
- `2026-04-26-plan-a-review-fixes.md` — Plan A Review Fixes (~2210 tok)
- `2026-04-26-sinkhorn-correctness-fixes.md` — Sinkhorn Correctness Fixes (~2423 tok)
- `2026-04-27-ieppa-cellmode-bound-fix.md` — iEPPA Cell-Mode Post-Normalization Bound Fix (~1708 tok)
- `2026-04-27-ieppa-convergence-overlays.md` — iEPPA Convergence Overlays Evaluation Plan (~4202 tok)
- `2026-04-27-ieppa-overflow-fix.md` — iEPPA Overflow Fix (T1.B) Implementation Plan (~6282 tok)
- `2026-04-27-python-ipf-benchmark.md` — Python IPF Benchmark vs Leafblower Methods (~2969 tok)
- `2026-04-27-raking-bregman-dykstra.md` — Raking: Bregman Dykstra + SOR + Greedy Implementation Plan (~5110 tok)
- `2026-04-27-sinkhorn-a1-fix-ieppa-soft.md` — Sinkhorn A1 Fix + ieppa_soft Method Implementation Plan (~6520 tok)
- `2026-04-28-chebyshev-sinkhorn-greg-correctness.md` — Chebyshev, Sinkhorn, Greg Correctness — Implementation Plan (~5804 tok)
- `2026-04-28-convergence-status.md` — Convergence Status Redesign Implementation Plan (~5347 tok)
- `2026-04-28-critical-review-fixes.md` — Critical Code Review Fixes Implementation Plan (~8526 tok)
- `2026-04-28-ieppa-raking-correctness.md` — iEPPA & Raking Algorithm Correctness — Implementation Plan (~6709 tok)
- `2026-04-28-performance-correctness-fixes.md` — Performance + Correctness Fixes Implementation Plan (~9944 tok)
- `2026-04-28-raking-squarem.md` — Raking SQUAREM Acceleration Implementation Plan (~7611 tok)
- `2026-04-28-raking-water-filling.md` — Raking Water-Filling Implementation Plan (~9501 tok)
- `2026-04-28-safety-memory-fixes.md` — Critical Safety & Memory Safety — Implementation Plan (~4836 tok)
- `2026-04-28-simd-performance.md` — SIMD Performance — libmvec Vectorization — Implementation Plan (~7968 tok)
- `2026-04-28-squarem-geometry.md` — SQUAREM Geometry Fix Implementation Plan (~3407 tok)
- `2026-04-28-status-code-propagation.md` — Status Code Propagation Fixes — Implementation Plan (~4988 tok)
- `2026-04-28-validation-error-detection.md` — Validation & Error Detection Gaps — Implementation Plan (~3746 tok)
- `2026-04-29-alm-epic-a-fixture.md` — ALM ieppa_soft — Epic A: Fixture Capture (~788 tok)
- `2026-04-29-alm-epic-b-abi.md` — ALM ieppa_soft — Epic B: ABI & Data Structures (~3417 tok)
- `2026-04-29-alm-epic-c-bridges.md` — ALM ieppa_soft — Epic C: R-Bridge + C-API Wiring (~4410 tok)
- `2026-04-29-alm-epic-d-solver.md` — ALM ieppa_soft — Epic D: Solver Implementation (T7) (~5970 tok)
- `2026-04-29-alm-epic-e-r-api.md` — ALM ieppa_soft Epic E — R API Implementation Plan (~6658 tok)
- `2026-04-29-chebyshev-greg-fix.md` — Chebyshev IPM Rewrite + GREG Quality Warning (~10603 tok)
- `2026-04-29-e0zb-greenkhorn-squarem.md` — Fix leafblower-e0zb: greenkhorn accelerate=TRUE silently ignores SQUAREM (~1934 tok)
- `2026-04-29-greenkhorn-logit-epic-a-infra.md` — Greenkhorn + Logit — Epic A: Infrastructure (~2231 tok)
- `2026-04-29-greenkhorn-logit-epic-b-greenkhorn.md` — Implementation Plan: Greenkhorn Logit — Epic B (Greenkhorn Solver) (~6154 tok)
- `2026-04-29-greenkhorn-logit-epic-c-logit.md` — Epic C — Logit Calibration Solver (logit_calib.{hpp,cpp}) (~4779 tok)
- `2026-04-29-greenkhorn-logit-epic-d-bridge.md` — Greenkhorn + Logit — Epic D: Bridge Wiring (~4056 tok)
- `2026-04-29-greenkhorn-logit-epic-e-tests.md` — Epic E: R API + Tests — Greenkhorn & Logit (~3140 tok)
- `2026-04-29-i0am-global-safeguard.md` — SRAA-m Global Safeguard Implementation Plan (~4783 tok)
- `2026-04-29-i0am-sraa-acceleration.md` — SRAA-m: Replace SQUAREM with Anderson Acceleration (~6308 tok)
- `2026-04-29-i0am-sraa-correct-all-scales.md` — SRAA-m Correct Acceleration for All K Scales — Combined Fix (~5487 tok)
- `2026-04-30-efficiency-hotpath.md` — Hot-Path Efficiency Fixes — Implementation Plan (~10383 tok)
- `2026-04-30-quality-fixes.md` — Code Quality Fixes — Implementation Plan (~10279 tok)
- `2026-04-30-reuse-dedup.md` — Code Reuse & Deduplication — Implementation Plan (~6740 tok)

## docs/superpowers/specs/

- `2026-04-18-bounded-convergence-fix-design.md` — Bounded Convergence Fix — Design Spec (~1799 tok)
- `2026-04-20-algo-selection-design.md` — Algorithm Selection Benchmark Design (~3262 tok)
- `2026-04-23-ieppa-convergence-hardening-design.md` — iEPPA Convergence Hardening Design (~6217 tok)
- `2026-04-23-ieppa-faithful-design.md` — Faithful iEPPA Solver Design (~7433 tok)
- `2026-04-24-convergence-metrics-sor-design.md` — Convergence Criterion Reform + SOR Stabilization + Best-Iterate Design (~5477 tok)
- `2026-04-24-ieppa-homotopy-greenkhorn-design.md` — Design Spec: iEPPA Homotopy + Greenkhorn Priority Scheduler + Tang Dynamic-η (~2321 tok)
- `2026-04-24-ieppa-speed-convergence-bounds-design.md` — iEPPA Speed, Convergence, and Bounds Hardening Design (~13094 tok)
- `2026-04-25-calibration-solvers-design.md` — Calibration Solvers Redesign — Spec (~5665 tok)
- `2026-04-25-convergence-redesign.md` — Convergence Criteria Redesign — Spec (~4632 tok)
- `2026-04-26-chebyshev-nu-reference-elimination.md` — Chebyshev ν Fix: Reference Category Elimination (~2670 tok)
- `2026-04-27-ieppa-linear-overflow-fix.md` — iEPPA Linear-Space Overflow Fix + Log-Path Acceleration (~3695 tok)
- `2026-04-27-raking-bregman-dykstra-design.md` — Raking: Bregman Dykstra + SOR + Greedy A/B Test (~2877 tok)
- `2026-04-27-sinkhorn-a1-fix-ieppa-admm-method.md` — Sinkhorn A1 Fix + ieppa_soft Method + Decoupled Solver Objective (~3213 tok)
- `2026-04-28-convergence-status-design.md` — Convergence Status Redesign (~2479 tok)
- `2026-04-28-raking-squarem-design.md` — Raking SQUAREM Acceleration Design (~2081 tok)
- `2026-04-28-squarem-geometry-fix.md` — SQUAREM Geometry Fix: Obs-Level α + Weight-Change Stall (~2555 tok)
- `2026-04-29-chebyshev-greg-fix.md` — Chebyshev IPM Rewrite + GREG Quality Warning (~3901 tok)
- `2026-04-29-greenkhorn-solver.md` — Greenkhorn + Logit Calibration Solvers — Design Spec (rev 5 — gate approved) (~7978 tok)
- `2026-04-29-i0am-sraa-acceleration.md` — i0am: Safeguarded Regularized Anderson Acceleration (SRAA-m) (~4818 tok)
- `2026-04-29-i0am-sraa-correct-all-scales.md` — i0am-C: SRAA-m Correct Acceleration for All K Scales (~2788 tok)
- `2026-04-29-i0am-sraa-global-safeguard.md` — i0am-B: SRAA-m Global Safeguard + Revert-to-Best (~2357 tok)
- `2026-04-29-ieppa-alm-soft-capacity.md` — iEPPA Soft: ALM Capacity Enforcement — Design Spec (rev 3) (~10424 tok)
- `2026-04-29-sv89-logit-newton-fix.md` — sv89: Logit Newton Stabilization — Design Spec (~3315 tok)

## graphify-out/

- `.graphify_python` (~4 tok)
- `.graphify_root` (~7 tok)
- `cost.json` (~59 tok)
- `GRAPH_REPORT.md` — Graph Report - leafblower  (2026-04-30) (~4494 tok)
- `graph.html` — graphify - graphify-out/graph.html (~66730 tok)
- `graph.json` (~85961 tok)
- `manifest.json` (~4366 tok)

## graphify-out/cache/

- `05af7a3aaea37e474fc2766605d77e5cc4d1df90de94fe4242cbd1289eecd605.json` (~251 tok)
- `067cd5f7a20a4f8472bbc0682732858e5bbadce017b9d5cba512ab46f0741bc2.json` (~106 tok)
- `071290fa1cbce02404c823717e24221209b8df8b9c320410e289af449b900291.json` (~248 tok)
- `1b3e9e1a82592492dc6b287dbf3ff657b882e0eb1b3e16b1bf87288779d658b4.json` (~220 tok)
- `1d715d94623974fea4df98023c6235318255c88d0640a433f30666e8a3911b50.json` (~244 tok)
- `206cc168b07cdd7ae02da9fad4162a3e9abb796fb51b9fd445ed8dd283a7c9a5.json` (~221 tok)
- `251563671063c8a0af4a0f3c012f5d2798d69e9603e5340dedb50774b9488ca0.json` (~457 tok)
- `2a114e7b14e999722092101b2ea11d3785c4c62b49f94cc4dc454733043c1003.json` (~249 tok)
- `2b91668e75d470a2334f96d9bf7a934d54ffe72935c5c46829b65aa423b30019.json` (~273 tok)
- `303eed22a8255fa872eda5884957f6ee5749656f1197310e5e941e70405e6dde.json` (~229 tok)
- `3326662435f5d2ffb3307052f2c8b7a66136d5dd2d54f6b87f0c89cbd355a090.json` (~246 tok)
- `388e4251a9bed70c8ee57bfddcc26a009b4bfb8904a7134483427fa991c386e9.json` (~195 tok)
- `39558b68d7749ed324b28379d129892fb4264f74ed1e7eeccb0aeae8c36134fd.json` (~7046 tok)
- `40da61f854d251272ade2923ab5a9100e744e5b859cd901f32448bfe7c34ce21.json` (~110 tok)
- `42a08c21d0168461659e5fe894b61ade7e74ba69c79d60e65fef265e4f7d2ff6.json` (~104 tok)
- `43271132868b745ab526038bc618ef64faa0e62ba56548cc02cf6e2980800b92.json` (~6912 tok)
- `47119dbe4a955e5fa28a905b1e85ba578a64b533c95a4ea209b4b3108432e3e7.json` (~7962 tok)
- `47b54ab49aaf46fbe2a1ce561ae2c483489b94c8b6cbd31850eafc9800a666ff.json` (~241 tok)
- `47f1dd7fb41578da933ac0fbd968e5f0603f53bc0d5886179ab47870b4712a1a.json` (~224 tok)
- `490b54a28002327d8de5db4aae04af6632f246d45aae36b7738e4c7eda26e66d.json` (~5920 tok)
- `4f6cae39d93e2edae46447110080b617b581b71c6a344f02543072dccd600dbc.json` (~460 tok)
- `593fcfbc33cd6dc7d05a774697df7c69e19613b4b8cd9d4d2f3910bc9a8878de.json` (~4805 tok)
- `63a49f2fc3a7ca5bae578b895b4bfddf67d296a76181a130889e8aeda89910a4.json` (~109 tok)
- `65afdee86b707f89e0d71524513f1f01a54b904c5112275d1359f56048233992.json` (~264 tok)
- `6630caf93978794a744dcf827cba4bf332c8ae93328202b124f5dedaaf01b558.json` (~322 tok)
- `6922c8de6abc125dff563b2c5c285781ad44db7d24935413e21f691d3de6c244.json` (~385 tok)
- `69d065ffc04054b2159bbeb53e1996dc07d13f358a69d4826d91295733b8912e.json` (~264 tok)
- `71242c807c2a0b64bf935b507686edb182b0cd482cc53df157fcb691252ed9d7.json` (~129 tok)
- `72ffe7b0f7722f9b26fb6c8fc27346ca81de6842a0ccfa5c2f90008f5be59714.json` (~285 tok)
- `7379c9885df632f63a86e695139b7f15a7843e5b31dfd23f5fad65943b4a18dd.json` (~146 tok)
- `7574ae05b2f08c595f7d618a8e96630cf06711cfcdc3d87144a87a4dbd406d6a.json` (~234 tok)
- `759349fedde8acdaf373aad643a2a08a744cc12ad475bcc283c91e66da52a896.json` (~488 tok)
- `7ffd8df39d0655b00796982cd2d2d841cbd37321a0b377010e38c081741d3507.json` (~249 tok)
- `80d276df1d21262627f98f43367b5b14f939cd58848dd7f95a9d6521a04d56df.json` (~246 tok)
- `82776448d5bfc1e1ede366018facdb7da2fe91980d74954c5f64c3272511dad5.json` (~91 tok)
- `82ea2765b2060c996d7959110f285246008076f70929f41061cc254c72290491.json` (~219 tok)
- `83867a04920f97d4d907227dec186789c9f6b69762414c62b70b396647f420e4.json` (~1164 tok)
- `85944fea414ea5915df287815d58fb9ea8a9bc70c72d61b9f5995472c7f92a85.json` (~5323 tok)
- `87ee12e470724d10c6ac1658ee24df0bff20a6daa253fcae6b030c901e9d08c2.json` (~438 tok)
- `94146a918c7988939a12299e8fd8cfe852d05f7944d159028ce07f0e2e8d5f3f.json` (~253 tok)
- `9cc5b9ffb2e7ef172d940c48c217f9a8fed8c7fc4f1bd276450007b99202b350.json` (~436 tok)
- `a435f2a247ce20fe812c72fa5b87c327a80f1892b0075df409f2cdff117d0858.json` (~3908 tok)
- `a66a725fed2f9bf2ce4a5f7a62f43021ea38a3db6a35fd8b3b8563b20ba10207.json` (~260 tok)
- `ad7f176b6a2371e4af0cb85368f40060dcc7c1d3cb3dfe14633091618ed142ec.json` (~259 tok)
- `aed46ee284bab950f15aaa2a43723f5d32f2cfef18e5ee81b3009a6a8db1bf99.json` (~5665 tok)
- `af71e51320828961f76cef6ced16ac138d28c2a4ee9fdc809f8da46c62360650.json` (~264 tok)
- `b1098bb08b286006da67cbf54337c5708163de8d1b56679fd458fa6bde978427.json` (~76 tok)
- `b2b3d802863b99a7e6854eeb337aecfc816881fa610c720f01295d0f7727bb34.json` (~337 tok)
- `b85338137e339486974a442de497e16bf31bda3577d7ec604ea89b2b0f365f97.json` (~104 tok)
- `bfa3120a76b2fb060e442cfe9336c266111360a135b5d0d817ac979b60bc0e93.json` (~159 tok)
- `c2ec4fa1ee0153bb604bbbf81491dd4955cfcae51cbe3f655a711d027b1ae400.json` (~296 tok)
- `c4b85f4e7f1d3834c24a5d2c42247a4b20c8abe2b80be681c0e122dc28305a8b.json` (~422 tok)
- `cac8835293fa982b7572f9a1921a0c97b69c2fc278da0537ac93d922e96477a6.json` (~154 tok)
- `cbd54d66af0e6ef43a82441132d01a40c1a38d7f783e5985751b882f108ef139.json` (~253 tok)
- `cc23001c2ebeb49c086d34c8e79549d5b9b88cf01be8fca6b53d1c178b119bcc.json` (~440 tok)
- `ccd663ace2f49a1ef194fbac6dbe9d45d295970018443021133565a33dd5b852.json` (~1282 tok)
- `d448c5fd9a908ff739852d41a700a9aa76cf2946b190f3f84c9b780cd3664cdc.json` (~254 tok)
- `d582fd91bc7fc0ba5cd2934e927b60430c3294e3265e041e20c895150f4fa52d.json` (~290 tok)
- `d629ad1907eec236d1daf2abb6f376de28badaec9c2606066ec6bbb82a4d215c.json` (~7046 tok)
- `d659b00f194eb6937556bcebee970b10aa028a0da06e50c93a96aff967965c6b.json` (~264 tok)
- `daec21cc4827209b3f5674844b235a9da3a38e4b400ff067303bbbd902b798bf.json` (~84 tok)
- `e109d63534e374b21bcc982159c71b32e6c6a1fddf085dcb189ac3f415f25b3b.json` (~415 tok)
- `e26fad4259c1b70b23573ce71d4d509d3e7b99b9ce7411f5901553acff54f81a.json` (~253 tok)
- `e2bae133e05964458a098ce14ff0c02ca6ddab8bf5ea5ce65bb3a54bb622a7eb.json` (~106 tok)
- `e33c2d6aeca974a95d9b9cd1dd83d8a6fefd82d76d00332ed07004e09a3fa7cf.json` (~152 tok)
- `e5b1c9fdaeb2b43a6944ad6f15c92e9f4125d971dc7dc46055b9756631c0e13a.json` (~248 tok)
- `e732e086fed8bd665b6e0270103e36246e6707131dd48d63b7cb099bc046a178.json` (~145 tok)
- `ed4eff8e622cd97a68367d48ceb6858fcb430e7164b1e9c218ac75464ae5855a.json` (~340 tok)
- `ee778ee1eda010c0228177439180e6c0c435335a35afa2e6c4814c7a115ee5ef.json` (~146 tok)
- `f0d0393d9147bcf0a2cfa649e87c27566894b9e9f175f19de0b37be829f391d2.json` (~73 tok)
- `f0d48a76004fd82176282e6bb3ee31810da32fdb4b7e8354b6919eb34bd5b9b3.json` (~378 tok)
- `f61fae7539b8e2c3c65488e27b76da70254022824981191cf2b5ffd67a2be264.json` (~7962 tok)
- `f82c016626a1b588ee0d332b0350dda2b927b504eae3c591499d57ae0dcc1126.json` (~597 tok)
- `fa41a656cfb40c4feadb89d71ab8e945428de96568cf0485adc4877b5c232966.json` (~288 tok)
- `fd57a7f239a4090c2b7f48bb929e37e505fa05761a963eb7c4e9652f7e4d6889.json` (~74 tok)

## leafblower.Rcheck/

- `00check.log` — Declares calls (~1084 tok)
- `00install.out` (~742 tok)

## leafblower.Rcheck/00_pkg_src/leafblower/

- `.tldrignore` — TLDR ignore patterns (gitignore syntax) (~280 tok)
- `cleanup` (~8 tok)
- `configure` — Detect C++17 support. (~1050 tok)
- `DESCRIPTION` (~246 tok)
- `LICENSE` — Project license (~16 tok)
- `NAMESPACE` (~64 tok)
- `NEWS.md` — leafblower (development) (~224 tok)
- `patch_wolfe.py` — patch_wolfe, patch_wolfe_line_search (~903 tok)

## leafblower.Rcheck/00_pkg_src/leafblower/.gemini-bridge/

- `feedback.log` (~80 tok)

## leafblower.Rcheck/00_pkg_src/leafblower/.pytest_cache/

- `CACHEDIR.TAG` (~51 tok)
- `README.md` — Project documentation (~76 tok)

## leafblower.Rcheck/00_pkg_src/leafblower/.pytest_cache/v/cache/

- `lastfailed` (~1 tok)
- `nodeids` (~57 tok)

## leafblower.Rcheck/00_pkg_src/leafblower/.tldr/

- `languages.json` (~35 tok)

## leafblower.Rcheck/00_pkg_src/leafblower/.tldr/cache/

- `call_graph.json` (~708 tok)

## leafblower.Rcheck/00_pkg_src/leafblower/R/

- `anesrake.R` — ' anesrake compatibility wrapper (~489 tok)
- `current_miss.R` — ' Current calibration miss (matching autumn's export name) (~190 tok)
- `design_effect.R` — ' Kish design effect (1-argument) or Henry-Valliant (4-argument) (~408 tok)
- `diagnose_weights.R` — ' Diagnose calibration quality (~466 tok)
- `harvest.R` — ' Generate calibrated weights (drop-in for autumn::harvest) (~2019 tok)
- `weighted_pct.R` — ' Weighted proportions (~145 tok)
- `zzz.R` — R_init_leafblower() in r_bridge.cpp is called automatically by R when the (~86 tok)

## leafblower.Rcheck/00_pkg_src/leafblower/man/

- `anesrake.Rd` (~312 tok)
- `design_effect.Rd` (~240 tok)
- `diagnose_weights.Rd` (~176 tok)
- `effective_sample_size.Rd` (~103 tok)
- `get_current_miss.Rd` (~141 tok)
- `harvest.Rd` (~722 tok)
- `weighted_pct.Rd` (~109 tok)

## leafblower.Rcheck/00_pkg_src/leafblower/src/

- `c_api.cpp` — include "leafblower.h" (~2373 tok)
- `c_api.o` (~2132 tok)
- `cell_table.cpp` — include "lbw_config.h" (~1397 tok)
- `cell_table.hpp` — pragma once (~285 tok)
- `cell_table.o` (~4840 tok)
- `ieppa.cpp` — include "lbw_config.h" (~3649 tok)
- `ieppa.hpp` — pragma once (~156 tok)
- `ieppa.o` (~5742 tok)
- `lbfgsb_solver.cpp` — include "lbw_config.h" (~6918 tok)
- `lbfgsb_solver.hpp` — pragma once (~74 tok)
- `lbfgsb_solver.o` (~12986 tok)
- `lbw_config.h` — ifndef LBW_CONFIG_H (~56 tok)
- `lbw_math.hpp` — pragma once (~260 tok)
- `leafblower.h` — ifndef LEAFBLOWER_H (~945 tok)
- `logit.cpp` — include "logit.hpp" (~112 tok)
- `logit.hpp` — pragma once (~921 tok)
- `logit.o` (~250 tok)
- `Makevars` (~50 tok)
- `Makevars.in` (~58 tok)
- `r_bridge.cpp` — include "leafblower.h" (~3080 tok)
- `r_bridge.o` (~10410 tok)
- `raking.cpp` — include "lbw_config.h" (~2397 tok)
- `raking.hpp` — pragma once (~50 tok)
- `raking.o` (~2859 tok)
- `symbols.rds` (~516 tok)
- `types.hpp` — pragma once (~334 tok)

## leafblower.Rcheck/00_pkg_src/leafblower/tests/

- `testthat.R` — This file is part of the R package leafblower. (~49 tok)

## leafblower.Rcheck/00_pkg_src/leafblower/tests/testthat/

- `lbfgsb_baseline_time.rds` (~27 tok)
- `lbfgsb_ref_weights.rds` (~2943 tok)
- `task1_ref.rds` (~531 tok)
- `task2_ieppa_ref.rds` (~1526 tok)
- `test-algo-selection.R` — tests/testthat/test-algo-selection.R (~2289 tok)
- `test-bounded-convergence.R` (~400 tok)
- `test-cell-table.R` (~497 tok)
- `test-compare.R` (~393 tok)
- `test-compat.R` (~292 tok)
- `test-design.R` (~223 tok)
- `test-harvest.R` (~1083 tok)
- `test-ieppa-faithful.R` (~1057 tok)
- `test-ieppa.R` (~611 tok)
- `test-lbfgsb.R` (~957 tok)
- `test-logit.R` (~308 tok)
- `test-raking.R` (~604 tok)

## leafblower.Rcheck/R_check_bin/

- `R` (~22 tok)
- `Rscript` (~24 tok)

## leafblower.Rcheck/leafblower/

- `DESCRIPTION` (~263 tok)
- `INDEX` (~139 tok)
- `LICENSE` — Project license (~16 tok)
- `NAMESPACE` (~64 tok)
- `NEWS.md` — leafblower (development) (~224 tok)

## leafblower.Rcheck/leafblower/Meta/

- `features.rds` (~32 tok)
- `hsearch.rds` (~134 tok)
- `links.rds` (~60 tok)
- `nsInfo.rds` (~98 tok)
- `package.rds` (~271 tok)
- `Rd.rds` (~132 tok)

## leafblower.Rcheck/leafblower/R/

- `leafblower` — File share/R/nspackloader.R (~283 tok)
- `leafblower.rdb` (~3935 tok)
- `leafblower.rdx` (~120 tok)

## leafblower.Rcheck/leafblower/help/

- `aliases.rds` (~44 tok)
- `AnIndex` (~54 tok)
- `leafblower.rdb` (~3176 tok)
- `leafblower.rdx` (~87 tok)
- `paths.rds` (~59 tok)

## leafblower.Rcheck/leafblower/html/

- `00Index.html` — R: High-Performance Survey Calibration via iEPPA and L-BFGS-B (~535 tok)
- `R.css` — Styles: 3 rules, 1 media queries (~527 tok)

## leafblower.Rcheck/leafblower/libs/

- `symbols.rds` (~516 tok)

## leafblower.Rcheck/tests/

- `startup.Rs` — # A custom startup file for tests (~38 tok)
- `testthat.R` — This file is part of the R package leafblower. (~49 tok)
- `testthat.Rout` — Declares to (~259 tok)

## leafblower.Rcheck/tests/testthat/

- `lbfgsb_baseline_time.rds` (~27 tok)
- `lbfgsb_ref_weights.rds` (~2943 tok)
- `task1_ref.rds` (~531 tok)
- `task2_ieppa_ref.rds` (~1526 tok)
- `test-algo-selection.R` — tests/testthat/test-algo-selection.R (~2289 tok)
- `test-bounded-convergence.R` (~400 tok)
- `test-cell-table.R` (~497 tok)
- `test-compare.R` (~393 tok)
- `test-compat.R` (~292 tok)
- `test-design.R` (~223 tok)
- `test-harvest.R` (~1083 tok)
- `test-ieppa-faithful.R` (~1057 tok)
- `test-ieppa.R` (~611 tok)
- `test-lbfgsb.R` (~957 tok)
- `test-logit.R` (~308 tok)
- `test-raking.R` (~604 tok)

## man/

- `anesrake.Rd` (~312 tok)
- `design_effect.Rd` (~240 tok)
- `diagnose_weights.Rd` (~176 tok)
- `effective_sample_size.Rd` (~103 tok)
- `get_current_miss.Rd` (~141 tok)
- `harvest.Rd` (~3256 tok)
- `weighted_pct.Rd` (~109 tok)

## python/

- `CMakeLists.txt` — CMake build configuration (~284 tok)
- `conftest.py` — Ensure the installed wheel's leafblower package is found, not the local (~74 tok)
- `pyproject.toml` — Python project configuration (~297 tok)

## python/.pytest_cache/

- `.gitignore` — Git ignore rules (~10 tok)
- `CACHEDIR.TAG` (~51 tok)
- `README.md` — Project documentation (~76 tok)

## python/.pytest_cache/v/cache/

- `lastfailed` (~11 tok)
- `nodeids` (~51 tok)

## python/leafblower/

- `__init__.py` — leafblower: high-performance survey calibration. (~42 tok)
- `_bindings.cpp` — include <pybind11/pybind11.h> (~2328 tok)

## src/

- `ieppa.cpp` — include "lbw_config.h" (~20966 tok)
- `ieppa.hpp` — pragma once (~795 tok)
- `lbfgsb_solver.cpp` — include "lbw_config.h" (~8343 tok)

## tests/testthat/

- `test-rk-params-passthrough.R` (~427 tok)
