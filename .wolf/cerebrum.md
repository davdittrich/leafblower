# Cerebrum

> OpenWolf's learning memory. Updated automatically as the AI learns from interactions.
> Do not edit manually unless correcting an error.
> Last updated: 2026-04-24

## User Preferences

<!-- How the user likes things done. Code style, tools, patterns, communication. -->

## Key Learnings

- **Project:** leafblower

## Do-Not-Repeat

<!-- Mistakes made and corrected. Each entry prevents the same mistake recurring. -->
<!-- Format: [YYYY-MM-DD] Description of what went wrong and what to do instead. -->

[2026-04-27] Task 4 A1 fixture fix: gen_ieppa_kl_ref.R was using `r$best_error` instead of `r$convergence_used$solver_objective`. Also the A1 test was using `r_s$convergence_used$objective` instead of `attr(w_s,"result")$convergence_used$solver_objective`. FIXED: both references now point to correct field.

## Decision Log

<!-- Significant technical decisions with rationale. Why X was chosen over Y. -->