# Graph Report - leafblower  (2026-05-05)

## Corpus Check
- 426 files · ~726,206 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 230 nodes · 273 edges · 23 communities detected
- Extraction: 85% EXTRACTED · 15% INFERRED · 0% AMBIGUOUS · INFERRED: 40 edges (avg confidence: 0.8)
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
- [[_COMMUNITY_Community 16|Community 16]]
- [[_COMMUNITY_Community 17|Community 17]]
- [[_COMMUNITY_Community 18|Community 18]]
- [[_COMMUNITY_Community 68|Community 68]]
- [[_COMMUNITY_Community 69|Community 69]]
- [[_COMMUNITY_Community 70|Community 70]]
- [[_COMMUNITY_Community 71|Community 71]]
- [[_COMMUNITY_Community 72|Community 72]]

## God Nodes (most connected - your core abstractions)
1. `harvest()` - 16 edges
2. `lbfgsb_solve_inner()` - 10 edges
3. `C_rk_calibrate()` - 9 edges
4. `build_cell_table()` - 9 edges
5. `rk_calibrate()` - 8 edges
6. `ieppa_solve()` - 7 edges
7. `solver_setup_ct()` - 7 edges
8. `main()` - 6 edges
9. `main()` - 6 edges
10. `_make_fixture()` - 6 edges

## Surprising Connections (you probably didn't know these)
- `ieppa_solve()` --calls--> `build_cell_table()`  [INFERRED]
  ieppa_92c4f45.cpp → src/cell_table.cpp
- `bench_harvest()` --calls--> `harvest()`  [INFERRED]
  benchmarks/stepstone_benchmark.py → python/leafblower/_harvest.py
- `main()` --calls--> `harvest()`  [INFERRED]
  benchmarks/stepstone_benchmark.py → python/leafblower/_harvest.py
- `bench_harvest()` --calls--> `harvest()`  [INFERRED]
  benchmarks/stepstone_fulldata_benchmark.py → python/leafblower/_harvest.py
- `main()` --calls--> `harvest()`  [INFERRED]
  benchmarks/stepstone_fulldata_benchmark.py → python/leafblower/_harvest.py

## Communities

### Community 0 - "Community 0"
Cohesion: 0.11
Nodes (24): diagnose_weights(), harvest(), _parse_convergence(), _parse_sor(), Calibrate survey weights. Drop-in for R leafblower::harvest().      Parameters, Derive pct_tol, absolute_tol, metric, rule, stop_when from convergence dict., # NOTE: No post-normalization clamp to [min_weight, max_weight]. Clamping, Diagnose calibration quality (Python equivalent of R diagnose_weights()).      P (+16 more)

### Community 1 - "Community 1"
Cohesion: 0.19
Nodes (16): PYBIND11_MODULE(), pack_solver_result(), rk_calibrate(), rk_params_init(), rk_result_init(), validate_inputs(), C_leafblower_cell_table_probe(), C_logit_F_at_zero() (+8 more)

### Community 2 - "Community 2"
Cohesion: 0.13
Nodes (12): apply_obs_expansion(), mark_converged(), cholesky_factor_inplace(), cholesky_solve(), compute_normal_equations(), chebyshev_ipm(), chebyshev_solve(), greg_solve() (+4 more)

### Community 3 - "Community 3"
Cohesion: 0.19
Nodes (11): aggregate_to_margin(), apply_rule(), build_cat_offset(), check_convergence(), compute_cell_bounds(), compute_cell_metrics(), resolve_hi(), select_metric() (+3 more)

### Community 4 - "Community 4"
Cohesion: 0.3
Nodes (14): build_offsets(), compute_du(), compute_final_weights_and_error(), compute_targets_abs(), compute_u(), dot(), lbfgs_direction(), lbfgsb_solve() (+6 more)

### Community 5 - "Community 5"
Cohesion: 0.47
Nodes (6): ieppa_solve(), bits_needed(), build_cell_table(), estimate_M_cell(), pack_key_compute(), pack_key_fits()

### Community 6 - "Community 6"
Cohesion: 0.43
Nodes (7): bench_harvest(), design_effect(), effective_sample_size(), load_data(), main(), Time leafblower.harvest() over n_runs; return timing + result stats., Kish (1965) design effect: n * sum(w²) / sum(w)².

### Community 7 - "Community 7"
Cohesion: 0.48
Nodes (5): ieppa_solve(), pack_lf(), parse_trajectory_iters(), unpack_lf(), write_trajectory_csv()

### Community 8 - "Community 8"
Cohesion: 0.53
Nodes (5): compute_metrics(), load_data(), main(), ipfn.IPFN in DataFrame mode on compressed cell table (28,905 cells)., run_ipfn()

### Community 9 - "Community 9"
Cohesion: 0.67
Nodes (5): bench_harvest(), design_effect(), effective_sample_size(), load_data(), main()

### Community 10 - "Community 10"
Cohesion: 0.33
Nodes (1): bulk_scaled_exp()

### Community 11 - "Community 11"
Cohesion: 0.8
Nodes (4): bits_needed(), build_cell_table(), pack_key_compute(), pack_key_fits()

### Community 12 - "Community 12"
Cohesion: 0.7
Nodes (3): compute_errRp(), raking_solve(), sum_weights_ilp()

### Community 13 - "Community 13"
Cohesion: 0.6
Nodes (4): _make_synthetic(), _r_weights(), Cross-language weight parity: R and Python must return identical weights for all, test_weight_parity()

### Community 14 - "Community 14"
Cohesion: 0.67
Nodes (3): post_hoc_mkl(), Compute obs-level marginal_kl from returned weights via groupby., run()

### Community 16 - "Community 16"
Cohesion: 0.67
Nodes (1): LinkFn()

### Community 17 - "Community 17"
Cohesion: 1.0
Nodes (1): Root conftest: remove python/ src tree from sys.path so the installed wheel is i

### Community 18 - "Community 18"
Cohesion: 1.0
Nodes (1): Ensure the installed wheel's leafblower package is found, not the local source t

### Community 68 - "Community 68"
Cohesion: 1.0
Nodes (1): Derive pct_tol, absolute_tol, metric, rule, stop_when from convergence dict.

### Community 69 - "Community 69"
Cohesion: 1.0
Nodes (1): Mirror R parse_sor(): returns (enabled, auto, omega_init, omega_min, omega_fixed

### Community 70 - "Community 70"
Cohesion: 1.0
Nodes (1): Calibrate survey weights. Drop-in for R leafblower::harvest().      Parameters

### Community 71 - "Community 71"
Cohesion: 1.0
Nodes (1): Diagnose calibration quality (Python equivalent of R diagnose_weights()).      P

### Community 72 - "Community 72"
Cohesion: 1.0
Nodes (1): # NOTE: No post-normalization clamp to [min_weight, max_weight]. Clamping

## Knowledge Gaps
- **24 isolated node(s):** `Root conftest: remove python/ src tree from sys.path so the installed wheel is i`, `Time leafblower.harvest() over n_runs; return timing + result stats.`, `Kish (1965) design effect: n * sum(w²) / sum(w)².`, `ipfn.IPFN in DataFrame mode on compressed cell table (28,905 cells).`, `Compute obs-level marginal_kl from returned weights via groupby.` (+19 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **Thin community `Community 10`** (6 nodes): `lbw_math.hpp`, `bulk_exp_clipped()`, `bulk_log()`, `bulk_scaled_exp()`, `bulk_scaled_log()`, `lbw_math.hpp`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 16`** (3 nodes): `logit.hpp`, `logit.hpp`, `LinkFn()`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 17`** (2 nodes): `conftest.py`, `Root conftest: remove python/ src tree from sys.path so the installed wheel is i`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 18`** (2 nodes): `conftest.py`, `Ensure the installed wheel's leafblower package is found, not the local source t`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 68`** (1 nodes): `Derive pct_tol, absolute_tol, metric, rule, stop_when from convergence dict.`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 69`** (1 nodes): `Mirror R parse_sor(): returns (enabled, auto, omega_init, omega_min, omega_fixed`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 70`** (1 nodes): `Calibrate survey weights. Drop-in for R leafblower::harvest().      Parameters`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 71`** (1 nodes): `Diagnose calibration quality (Python equivalent of R diagnose_weights()).      P`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 72`** (1 nodes): `# NOTE: No post-normalization clamp to [min_weight, max_weight]. Clamping`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `newton_calibrate()` connect `Community 2` to `Community 3`, `Community 4`?**
  _High betweenness centrality (0.048) - this node is a cross-community bridge._
- **Why does `solver_setup_ct()` connect `Community 3` to `Community 2`, `Community 5`?**
  _High betweenness centrality (0.040) - this node is a cross-community bridge._
- **Why does `harvest()` connect `Community 0` to `Community 9`, `Community 13`, `Community 6`, `Community 14`?**
  _High betweenness centrality (0.038) - this node is a cross-community bridge._
- **Are the 12 inferred relationships involving `harvest()` (e.g. with `bench_harvest()` and `main()`) actually correct?**
  _`harvest()` has 12 INFERRED edges - model-reasoned connections that need verification._
- **Are the 4 inferred relationships involving `C_rk_calibrate()` (e.g. with `rk_params_init()` and `rk_calibrate()`) actually correct?**
  _`C_rk_calibrate()` has 4 INFERRED edges - model-reasoned connections that need verification._
- **Are the 4 inferred relationships involving `build_cell_table()` (e.g. with `ieppa_solve()` and `ieppa_solve()`) actually correct?**
  _`build_cell_table()` has 4 INFERRED edges - model-reasoned connections that need verification._
- **Are the 2 inferred relationships involving `rk_calibrate()` (e.g. with `PYBIND11_MODULE()` and `C_rk_calibrate()`) actually correct?**
  _`rk_calibrate()` has 2 INFERRED edges - model-reasoned connections that need verification._