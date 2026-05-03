# Graph Report - leafblower  (2026-05-03)

## Corpus Check
- 224 files · ~584,443 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 225 nodes · 275 edges · 18 communities detected
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
- [[_COMMUNITY_Community 15|Community 15]]
- [[_COMMUNITY_Community 16|Community 16]]
- [[_COMMUNITY_Community 17|Community 17]]

## God Nodes (most connected - your core abstractions)
1. `harvest()` - 14 edges
2. `lbfgsb_solve_inner()` - 11 edges
3. `rk_calibrate()` - 9 edges
4. `build_cell_table()` - 9 edges
5. `ieppa_solve()` - 7 edges
6. `solver_setup_ct()` - 7 edges
7. `main()` - 6 edges
8. `main()` - 6 edges
9. `_make_fixture()` - 6 edges
10. `C_rk_calibrate()` - 6 edges

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
Nodes (24): diagnose_weights(), harvest(), _parse_convergence(), _parse_sor(), Calibrate survey weights. Drop-in for R leafblower::harvest().      Parameters, # NOTE: No post-normalization clamp to [min_weight, max_weight]. Clamping, Derive pct_tol, absolute_tol, metric, rule, stop_when from convergence dict., Diagnose calibration quality (Python equivalent of R diagnose_weights()).      P (+16 more)

### Community 1 - "Community 1"
Cohesion: 0.13
Nodes (12): apply_obs_expansion(), mark_converged(), compute_normal_equations(), ldlt_factor_inplace(), ldlt_solve(), chebyshev_ipm(), chebyshev_solve(), greg_solve() (+4 more)

### Community 2 - "Community 2"
Cohesion: 0.38
Nodes (14): build_offsets(), compute_du(), compute_final_weights_and_error(), compute_targets_abs(), compute_u(), dot(), lbfgs_direction(), lbfgsb_solve() (+6 more)

### Community 3 - "Community 3"
Cohesion: 0.19
Nodes (11): aggregate_to_margin(), apply_rule(), build_cat_offset(), check_convergence(), compute_cell_bounds(), compute_cell_metrics(), resolve_hi(), select_metric() (+3 more)

### Community 4 - "Community 4"
Cohesion: 0.23
Nodes (11): ieppa_solve(), bits_needed(), build_cell_table(), estimate_M_cell(), pack_key_compute(), pack_key_fits(), ieppa_solve(), pack_lf() (+3 more)

### Community 5 - "Community 5"
Cohesion: 0.33
Nodes (8): PYBIND11_MODULE(), pack_lbfgsb_result(), pack_solver_result(), rk_calibrate(), rk_params_init(), rk_result_init(), validate_inputs(), C_rk_calibrate()

### Community 6 - "Community 6"
Cohesion: 0.43
Nodes (7): bench_harvest(), design_effect(), effective_sample_size(), load_data(), main(), Time leafblower.harvest() over n_runs; return timing + result stats., Kish (1965) design effect: n * sum(w²) / sum(w)².

### Community 7 - "Community 7"
Cohesion: 0.43
Nodes (6): C_leafblower_cell_table_probe(), C_logit_F_at_zero(), C_logit_Hprime_check(), C_logit_range_check(), R_init_leafblower(), r_log_trampoline()

### Community 8 - "Community 8"
Cohesion: 0.29
Nodes (4): cp_calibrate(), ipm_calibrate(), cp_solve_R(), ipm_solve_R()

### Community 9 - "Community 9"
Cohesion: 0.53
Nodes (5): compute_metrics(), load_data(), main(), ipfn.IPFN in DataFrame mode on compressed cell table (28,905 cells)., run_ipfn()

### Community 10 - "Community 10"
Cohesion: 0.67
Nodes (5): bench_harvest(), design_effect(), effective_sample_size(), load_data(), main()

### Community 11 - "Community 11"
Cohesion: 0.33
Nodes (1): bulk_scaled_exp()

### Community 12 - "Community 12"
Cohesion: 0.8
Nodes (4): bits_needed(), build_cell_table(), pack_key_compute(), pack_key_fits()

### Community 13 - "Community 13"
Cohesion: 0.7
Nodes (3): compute_errRp(), raking_solve(), sum_weights_ilp()

### Community 14 - "Community 14"
Cohesion: 0.67
Nodes (2): patch_wolfe(), patch_wolfe_line_search()

### Community 15 - "Community 15"
Cohesion: 0.67
Nodes (1): LinkFn()

### Community 16 - "Community 16"
Cohesion: 1.0
Nodes (1): Root conftest: remove python/ src tree from sys.path so the installed wheel is i

### Community 17 - "Community 17"
Cohesion: 1.0
Nodes (1): Ensure the installed wheel's leafblower package is found, not the local source t

## Knowledge Gaps
- **17 isolated node(s):** `Root conftest: remove python/ src tree from sys.path so the installed wheel is i`, `Time leafblower.harvest() over n_runs; return timing + result stats.`, `Kish (1965) design effect: n * sum(w²) / sum(w)².`, `ipfn.IPFN in DataFrame mode on compressed cell table (28,905 cells).`, `Ensure the installed wheel's leafblower package is found, not the local source t` (+12 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **Thin community `Community 11`** (6 nodes): `lbw_math.hpp`, `bulk_exp_clipped()`, `bulk_log()`, `bulk_scaled_exp()`, `bulk_scaled_log()`, `lbw_math.hpp`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 14`** (4 nodes): `patch_wolfe()`, `patch_wolfe_line_search()`, `patch_wolfe.py`, `patch_wolfe.py`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 15`** (3 nodes): `logit.hpp`, `logit.hpp`, `LinkFn()`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 16`** (2 nodes): `conftest.py`, `Root conftest: remove python/ src tree from sys.path so the installed wheel is i`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Community 17`** (2 nodes): `conftest.py`, `Ensure the installed wheel's leafblower package is found, not the local source t`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `newton_calibrate()` connect `Community 1` to `Community 2`, `Community 3`?**
  _High betweenness centrality (0.052) - this node is a cross-community bridge._
- **Why does `solver_setup_ct()` connect `Community 3` to `Community 1`, `Community 4`?**
  _High betweenness centrality (0.043) - this node is a cross-community bridge._
- **Why does `build_cell_table()` connect `Community 4` to `Community 3`?**
  _High betweenness centrality (0.038) - this node is a cross-community bridge._
- **Are the 10 inferred relationships involving `harvest()` (e.g. with `bench_harvest()` and `main()`) actually correct?**
  _`harvest()` has 10 INFERRED edges - model-reasoned connections that need verification._
- **Are the 2 inferred relationships involving `rk_calibrate()` (e.g. with `PYBIND11_MODULE()` and `C_rk_calibrate()`) actually correct?**
  _`rk_calibrate()` has 2 INFERRED edges - model-reasoned connections that need verification._
- **Are the 4 inferred relationships involving `build_cell_table()` (e.g. with `ieppa_solve()` and `ieppa_solve()`) actually correct?**
  _`build_cell_table()` has 4 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Root conftest: remove python/ src tree from sys.path so the installed wheel is i`, `Time leafblower.harvest() over n_runs; return timing + result stats.`, `Kish (1965) design effect: n * sum(w²) / sum(w)².` to the rest of the system?**
  _17 weakly-connected nodes found - possible documentation gaps or missing edges._