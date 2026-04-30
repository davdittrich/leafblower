# Code Review: harvest.R & calib_dispatch.hpp

## Issues Summary (8 total; ranked by severity)

| Rank | Severity | Issue | File:Line | Confidence |
|------|----------|-------|-----------|------------|
| 1 | HIGH | Chebyshev objective returns errRp, should return grake_norm (L∞) | src/calib_dispatch.hpp | 80 |
| 2 | HIGH | apply_rule 1e-15 threshold universal for all metrics; chi2 threshold metric-dependent | src/calib_dispatch.hpp | 85 |
| 3 | MEDIUM | PCT stall guard missing chi2, grake_norm, l1_weight metrics | R/harvest.R:360 | 70 |
| 4 | MEDIUM | check_convergence side effect (prev_metric update) undocumented | src/calib_dispatch.hpp | 65 |
| 5 | MEDIUM | design_weights silently ignored when start_weights provided; no warning | R/harvest.R:103 | 60 |
| 6 | LOW | stall_wchange never fires for fixed-point guard; unclear intent | R/harvest.R:~420 | 50 |
| 7 | LOW | pct metric alias duplicate of l1_weight; inconsistent naming downstream | R/harvest.R:198 | 45 |
| 8 | LOW | map_method unknown fallback to match.arg; defensive but appropriate | R/harvest.R:391 | 40 |

## Details

### Issue #1: Chebyshev Objective (HIGH)
- **File:** src/calib_dispatch.hpp, select_solver_objective()
- **Problem:** `case RK_ALG_CHEBYSHEV: return m.errRp;` but chebyshev minimizes L∞, not errRp.
- **Fix:** Return `m.grake_norm` for chebyshev (L∞ norm).

### Issue #2: apply_rule Zero-Guard (HIGH)
- **File:** src/calib_dispatch.hpp, apply_rule()
- **Problem:** `if (curr <= 1e-15) converged = true` fires for all metrics. CHI2 is sample-size scaled; 1e-15 threshold inappropriate.
- **Fix:** Document threshold or move zero-guard inside metric-specific blocks.

### Issue #3: PCT Stall Guard Incomplete (MEDIUM)
- **File:** R/harvest.R:360
- **Problem:** Guard only checks max_err, mean_err. Should include chi2, grake_norm, l1_weight.
- **Fix:** Expand guard to all valid metrics.

### Issue #4: Side Effect Undocumented (MEDIUM)
- **File:** src/calib_dispatch.hpp, check_convergence()
- **Problem:** `apply_rule()` updates `prev_metric` even when convergence=false. Not documented.
- **Fix:** Add comment in function header.

### Issue #5: design_weights Silent Override (MEDIUM)
- **File:** R/harvest.R:103
- **Problem:** When both design_weights and start_weights provided, design_weights silently ignored. No warning.
- **Fix:** Add warning.

### Issue #6: stall_wchange Logic Unclear (LOW)
- **File:** R/harvest.R:~420
- **Problem:** stall_wchange set when status==5 && accelerate=TRUE. But fixed-point guard sets status=RK_OK, not 5. Reason never fires for fixed-point case.
- **Fix:** Clarify intent or add separate convergence reason for fixed-point guard.

### Issue #7: pct Alias Duplicate (LOW)
- **File:** R/harvest.R:198
- **Problem:** `pct = 5L` duplicates `l1_weight = 5L`. Creates naming inconsistency downstream.
- **Fix:** Normalize to l1_weight after parsing; document alias in comment.

### Issue #8: map_method Fallback (LOW)
- **File:** R/harvest.R:391
- **Problem:** Unknown methods fall through to match.arg(). No bug; defensive code works as intended.
- **Fix:** None needed.
