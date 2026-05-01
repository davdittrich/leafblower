# Active Plan
<!-- approved: 2026-05-01 -->
<!-- gate-iterations: 2 -->
<!-- user-approved: pending -->
<!-- status: ready_for_execution -->
<!-- epic: leafblower-usg8 -->
<!-- plan-doc: docs/superpowers/plans/2026-05-01-newton-kl-ieppa-warmstart-plan.md -->

# Newton-KL IEPPA Warm-Start (Epic-C)

**Triggered by:** Epic-B (target homotopy) BLOCKED — T2 stepstone regressed 2.8e-4 → 2.29e-3.
**Plan rev 3 approved:** plan-review-gate iter 2, all 3 reviewers PASS (Feasibility, Completeness, Scope & Alignment).

## Sequence (9 atomic tickets)

WI-0 (usg8.1) spec amend → WI-0b (usg8.8) basin-overlap kill-switch → WI-1 (usg8.2) ieppa lf-capture shim + K_warm sweep → WI-2 (usg8.3) wiring + rename + NEWS + diagnostic + #include → WI-1c (usg8.9) r_bridge SEXP surface → WI-3 (usg8.4) T6-T10 tests → WI-4 (usg8.5) verify suite → WI-5a (usg8.6) bench → WI-5b (usg8.7) verdict + epic close.

## Mechanism

Run IEPPA+SRAA inner for K_warm=8 sweeps at original T → convert lf to Newton's λ via `λ_{k,j} = lf[cat_offset_ieppa[k]+j] - lf[cat_offset_ieppa[k]+0]` for j ≥ 1 (rebased to Newton's lam_off layout) → run existing `run_newton_inner` from this λ. Mathematical no-op handoff: lf and λ produce identical weights up to LSE-absorbed constant.

Three-tier fallback: SRAA-fail → plain IEPPA → λ=0 cold start. Strictly additive — warm-start can only improve, never regress.

## Verdict gates (WI-5b)

- **GATE_MET:** stepstone <1e-4 AND kk1204 max_err<1e-4 AND kk1204 wall<3s — close epic.
- **PARTIAL:** stepstone <1e-4 AND (kk1204 max_err<1e-4 with wall ≥3s, OR max_err ∈ [1e-4, 1e-3]) — close epic, file Epic-D follow-up.
- **BLOCKED:** stepstone ≥1e-4 OR kk1204 diverges — file Epic-E (alternative path).

Full plan: `docs/superpowers/plans/2026-05-01-newton-kl-ieppa-warmstart-plan.md` (rev 3).
