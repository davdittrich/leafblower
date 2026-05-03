# Cerebrum

> OpenWolf's learning memory. Updated automatically as the AI learns from interactions.
> Do not edit manually unless correcting an error.
> Last updated: 2026-05-03

## User Preferences

<!-- How the user likes things done. Code style, tools, patterns, communication. -->

## Key Learnings

- **Refactoring Pattern (fcbo.2):** When extracting duplicated code across multiple solvers (chebyshev, sinkhorn, greg, raking), use context-mode batch_execute with grep to find all sites first, then add helpers to calib_dispatch.hpp, then replace sites one file at a time with Edit tool. Prefer using references to vectors passed by reference rather than initializing them with sizes, letting the helper function do the resizing.
- **compute_weight_kl vectorization (fcbo.5):** When hoisting lambdas that compute weight-space KL, the vectorized version uses bulk_log for performance. Simple scalar paths (sinkhorn/ieppa) can be unified with vectorized path by passing caller-owned scratch buffers. Raking already had scratch buffers; sinkhorn/ieppa needed new vector allocations. All 8 call sites across 3 solvers use identical function signature.

## Do-Not-Repeat

<!-- Mistakes made and corrected. Each entry prevents the same mistake recurring. -->
<!-- Format: [YYYY-MM-DD] Description of what went wrong and what to do instead. -->

[2026-04-27] Task 4 A1 fixture fix: gen_ieppa_kl_ref.R was using `r$best_error` instead of `r$convergence_used$solver_objective`. Also the A1 test was using `r_s$convergence_used$objective` instead of `attr(w_s,"result")$convergence_used$solver_objective`. FIXED: both references now point to correct field.

## Decision Log

<!-- Significant technical decisions with rationale. Why X was chosen over Y. -->