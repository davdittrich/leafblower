# Active Plan
<!-- approved: 2026-05-02 -->
<!-- gate-iterations: 1 -->
<!-- user-approved: pending -->
<!-- status: ready_for_execution -->
<!-- epic: leafblower-8aex -->
<!-- plan-doc: docs/superpowers/plans/2026-05-02-newton-kl-hygiene-plan.md -->

# Newton-KL post-Epic-Dβ hygiene + AUTO safety (Epic-H)

**Triggered by:** 5 Newton-KL epics shipped (Epic-A LM, Epic-B target homotopy BLOCKED, Epic-C IEPPA warm-start BLOCKED, Epic-Dα cheap-first ESCALATE, Epic-Dβ TSVD+CG PARTIAL). Stepstone basin floor at 2.61e-4 unsolved. kk1204 K=20 severe-skew converges to high-error fixed point. Epic-H closes 4 cleanup items deferred across the chain.

**Plan rev 1 approved:** plan-review-gate iter 1, all 3 reviewers PASS.

## Sequence (4 atomic tickets, sequential)

WH-c (8aex.1) NEWS.md consolidation + IEPPA warm-start spec erratum → WH-d (8aex.2) stale field cleanup (delete `n_homotopy_levels_used`; surface `lm_mu_final` via r_bridge) → WH-e (8aex.3) `newton_tsvd_ratio` user parameter + I1 `dsy_info` diagnostic + I2 degenerate λ_max≤0 message → WH-g (8aex.4) AUTO routing target-skew gate (severe-skew K≥5 → ieppa+sraa) — BREAKING.

Sequential dep chain justified: shared NEWS.md surface + struct delete-then-add ordering on `NewtonCalibResult`.

## Verdict gates

Implicit per Epic Success Criteria checkbox. No formal GATE_MET/PARTIAL/BLOCKED tree (hygiene work — pass/fail per checkbox).

Full plan: `docs/superpowers/plans/2026-05-02-newton-kl-hygiene-plan.md`.
