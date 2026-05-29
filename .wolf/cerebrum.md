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

[2026-05-29] Do NOT "fix" the chebyshev.cpp Phase-B corrector by dropping the linear `y*dS_aff` term (leafblower-xiox proposed this, conf 60, citing N-W 19.55). The full corrector is `corr = y*dS_aff + y*dS_aff^2/s = -dS_aff*dY_aff` where the affine dual step is `dY_aff = -y - y*dS_aff/s`. The ticket's derivation dropped the `-y`. Verified correct vs CLARABEL LP reference (solver hits optimum ~1e-11 on converged problems). The K≥9 stall is a separate conditioning issue, NOT the corrector. Guard comment now at the corrector. Lesson: reviewer formula claims (even Opus, even with a citation) require independent first-principles re-derivation + external numerical reference before acting — see [[feedback_subagent_numerical_verification]].

[2026-05-29] Do NOT add a "-rho" term to the ieppa.cpp ALM linearized Newton step (leafblower-7emq proposed it, conf 50). The code X=X_tilde*(1-lambda+mu*z)/(1+rho) is the correct linearization of iEPPA's UN-normalized KL / I-divergence generator (dD/dX=log(X/X_tilde), not log+1). The ticket's -rho derives from the normalized KL — wrong generator. 2-cell sim: both forms share fixed point clamp(X_tilde,L,U), rate is a wash. Guard comment at the step (commit 7c62d50). SECOND over-stated formula ticket this session (after xiox): the pattern is clear — SUSPECTED/low-confidence formula tickets need derivation-from-the-actual-generator + numerical sim BEFORE editing; the reviewer's premise (which divergence/objective) is the thing to check first.

## Decision Log

<!-- Significant technical decisions with rationale. Why X was chosen over Y. -->