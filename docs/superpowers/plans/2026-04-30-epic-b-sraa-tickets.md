# Epic B: SRAA-m iEPPA Implementation — Ticket Summary

**Created:** 2026-04-30  
**Spec:** docs/superpowers/specs/2026-04-30-ieppa-sraa-acceleration.md  
**Epic ID:** `leafblower-9zfx`

## Epic Overview

Core implementation of SRAA-m (Safeguarded Regularized Anderson Acceleration) for iEPPA using lf-space iterates (log-factors). Replaces the `iter_in_lvl` for-loop with a SRAA while-loop when `accelerate=TRUE && use_linear`. Estimated scope: 5 sequential tasks with compile gates between each.

## Tickets

### Epic Parent
- **leafblower-9zfx** [P1 epic] — Overall coordination
  - Success: accelerate=TRUE converges on stepstone_small; res$aa_accepted_count > 0; log-path fallback works; ALM reseed works; all tests pass

### Task B1: Helpers
- **leafblower-9zfx.1** [P1 task] — `pack_lf` / `unpack_lf` helpers in ieppa.cpp
  - **Input:** src/ieppa.cpp
  - **Output:** Two static helper functions (O(n_cats_total_with_na) and O(n_cats_total + K·M_cell) respectively)
  - **Compile gate:** R CMD INSTALL --preclean .
  - **Key constraints:** preserve alpha damping; no changes to non-SRAA paths

### Task B2: F_eval Lambda
- **leafblower-9zfx.2** [P1 task] — `f_eval_lf` lambda (captures sweep state, handles overflow)
  - **Input:** src/ieppa.cpp (inside homotopy level loop)
  - **Output:** Lambda closure that evals F in lf-space, calls apply_single_margin_linear for each margin, packs result, returns infinity on overflow
  - **Compile gate:** R CMD INSTALL --preclean .
  - **Key constraints:** capture by reference; overflow → infinity; partial lf state preserved for SRAA safeguard

### Task B3: SRAA Outer Loop
- **leafblower-9zfx.3** [P1 task] — SRAA outer loop with `lf_best` tracking (replaces for-loop)
  - **Input:** src/ieppa.cpp (homotopy level, after linear phase)
  - **Output:** if/else branch: accelerate=TRUE path uses SRAA while-loop; accelerate=FALSE path uses existing for-loop
  - **Compile gate:** R CMD INSTALL --preclean .
  - **Key constraints:**
    - F_cur seeded ONCE before loop (not reset inside)
    - lf_best = lf_flat snapshot (for SRAA stall revert), separate from W_best (normalized ratios for output)
    - Outer stall guard with revert-to-best logic
    - res.base.iterations reports F-evals when accelerate=TRUE

### Task B4: ALM Reseed
- **leafblower-9zfx.4** [P1 task] — ALM ieppa_soft: clear+reseed on capacity_mu update
  - **Input:** src/ieppa.cpp (ALM violation handler)
  - **Output:** Logic block: lambda reset → SRAA clear → pack_lf reseed
  - **Compile gate:** R CMD INSTALL --preclean .
  - **Key constraints:**
    - Order matters: lambda reset before pack_lf
    - aa_accepted_count persists across clear() (cumulative across ALM levels)

### Task B5: SOR / Greedy
- **leafblower-9zfx.5** [P1 task] — SOR disable & greedy downgrade when accelerate=TRUE
  - **Input:** src/ieppa.cpp (solver initialization)
  - **Output:**
    - SOR disable: `sor_auto_v = st.sor_cfg.auto_adapt && !st.accelerate`
    - Greedy downgrade: log message at verbose≥1, use round-robin instead
  - **Compile gate:** R CMD INSTALL --preclean .
  - **Key constraints:**
    - accelerate=FALSE preserves existing SOR behavior
    - Log message only at verbose≥1

## Execution Order

```
B1 (helpers) → B2 (lambda) → B3 (outer loop) → B4 (ALM) → B5 (SOR/greedy)
```

Each task has an implicit compile gate: R CMD INSTALL --preclean . must succeed before proceeding to the next task. If any task fails to compile, the execution stops and the issue is logged.

## Key Design Points

1. **lf-space iteration:** n_cats_total_with_na ≈ 50–500 dims (vs M_cell ≈ 1M on large surveys). SRAA history: ~13KB on stepstone vs ~64MB in X-space.

2. **F_cur seeding:** Seed F_cur = lf_flat ONCE before loop. sraa_step carries F_cur forward via swap inside its state; do NOT reset F_cur inside the while-loop.

3. **Two tracking variables:**
   - `W_best[c] = X[c] / X_init[c]` (normalized ratios) — existing field, unchanged, used for res.best_weights output
   - `lf_best` (lf-space iterate) — NEW, used only for SRAA stall revert; cannot be derived from W_best alone

4. **Overflow handling:** When apply_single_margin_linear returns true, lf is partially updated. pack_lf writes this partial state. SRAA safeguard detects err_AA=inf > err_plain, reverts lf_flat to F_cur, restoring consistency.

5. **Log-path fallback:** When linear overflow → use_linear=false, SRAA history is cleared and aa_accepted_count is snapshotted. Control falls through to non-accelerated for-loop.

6. **ALM constraint:** Both lambda_cell reset AND SRAA clear must occur together in the same code block.

## Testing

New file: tests/testthat/test-ieppa-sraa.R

Key assertions:
- Convergence parity (accelerate vs plain at same max_iterations)
- Iteration reduction (≥30% fewer F-evals on tight problems)
- ieppa_soft convergence with bounds
- Greedy downgrade (no error; verbose message)
- Output correlation (w_plain vs w_sraa > 0.9999)
- aa_accepted_count ≥ 5 on tight problems
- Regression: equal budget no worse than plain

## Spec Reference

Full specification: docs/superpowers/specs/2026-04-30-ieppa-sraa-acceleration.md

Depends on: Epic A (apply_clamp parameter + aa_accepted_count field in IEPPAResult)
