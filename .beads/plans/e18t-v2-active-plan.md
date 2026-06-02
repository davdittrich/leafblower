# Active Plan — e18t v2 iterate-change ω (gate-approved)
<!-- approved: 2026-06-02 -->
<!-- gate-iterations: 2 -->
<!-- user-approved: design yes; plan gate-approved iter 2 -->
<!-- status: ready -->

## Root cause (3-way verified)
Spectral theory correct; mode-2 estimator used wrong observable. Marginal "free residual"
≡ global residual (clamped mass cancels) → plateaus at infeasibility floor on bounded
problems → ratio→1 → ω→1.8 → stall (mj1p.2 reborn). Verified: algebra + SageMath
(/tmp/e18t_rootcause/) + NotebookLM Lehmann-SOR notebook 1e3036a1.

## Fix
Replace mode-2 marginal θ₂ with free-coordinate iterate-change ‖ΔX_free‖² power-iteration
estimator, single global ω, with ^(1/I) cadence recovery (also fixes latent e18t.3
lag-10=ρ²⁰ bug). Spec: docs/superpowers/specs/2026-06-02-oris-iterate-change-omega-design.md (rev 2).

## Tickets (chain 7→8→10→9)
- e18t.7  Derive iterate-change θ₂ in SageMath (Phase-1 GO gate)
- e18t.8  Replace mode-2 marginal estimator with global iterate-change in oris.cpp
- e18t.10 Wire omega_mode_id through Python binding (R/Python parity; no new ABI)
- e18t.9  Re-run ship gate, decide default, flip omega_mode_id in all 3 sites
          (R parse_sor + types.hpp + c_api.cpp), finalize docs/tests/NEWS

## Gate (plan-review-gate, iter 2): Feasibility PASS / Completeness PASS / Scope PASS
