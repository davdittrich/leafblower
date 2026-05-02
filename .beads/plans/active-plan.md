# Active Plan
<!-- approved: 2026-05-02 -->
<!-- gate-iterations: 2 -->
<!-- user-approved: pending -->
<!-- status: ready_for_execution -->
<!-- epic: leafblower-wkmq -->
<!-- plan-doc: docs/superpowers/plans/2026-05-02-newton-kl-tsvd-cg-plan.md -->

# Newton-KL TSVD + Steihaug-CG (Epic-Dβ)

**Triggered by:** Epic-B (target homotopy) BLOCKED, Epic-C (IEPPA warm-start) BLOCKED, Epic-Dα (cheap-first lm_mu/pivot sweep) ESCALATE per `benchmarks/results/lm_pivot_sweep.csv` showing 2.79e-04 plateau across {1e-12, 1e-10, 1e-8}. Bottleneck is structural Newton-step direction, not conditioning. User picked B3 composed pipeline; LAPACK explicitly approved.

**Plan rev 2 approved:** plan-review-gate iter 2, all 3 reviewers PASS.

## Sequence (7 atomic tickets)

WL-0 (wkmq.1) Makevars+Lapack wrappers → WL-1 (wkmq.2) eigendecomp+TSVD pinv + r_bridge SEXP-pack + n_projected_dims field → WL-2 (wkmq.3) T2 stepstone pinv-only measurement gate → WL-3 (wkmq.4) Steihaug-CG trust-region → WL-4 (wkmq.5) T7 K=4 over-projection HARD gate → WL-5 (wkmq.6) kk1204 K=20 PARTIAL gate → WL-6 (wkmq.7) verdict + best-iterate audit + cleanup ticket filing.

## Mechanism

1. LAPACK `dsyevd` eigendecomposition of `H_pre = V Λ V^T` (symmetric PSD, n_λ ≤ 80).
2. Truncated-SVD: drop directions with `λ_i < 1e-8 × λ_max`. Project `g_keep = V_keep^T G`, compute `Λ_damped = Λ_keep × (1+μ) + μ × d_floor_retained`, get `δ_keep = g_keep / Λ_damped`. Back-project `δ = V_keep · δ_keep`.
3. Steihaug-CG trust-region: when `||δ_pinv||₂ > Δ`, use CG iterates bounded by Δ. Diagonal H_proj makes math simple.
4. ρ-formula: `(g_curr - g_trial) / (-G^T·δ - 0.5·δ^T·H·δ)`. Δ adapts: ρ>0.75 ⇒ Δ ← 2Δ (no upper cap, Nocedal-Wright Alg 4.1); ρ<0.25 ⇒ Δ ← Δ/4.
5. LM damping composed in eigenbasis (preserves scale-invariance from Epic-A LM rev 2).

## Verdict gates (WL-6)

- **GATE_MET:** T2 stepstone <1e-4 AND T7 K=4 over-projection (n_projected_dims==0) AND kk1204 K=20 max_err<1e-4 (any wall) → close epic.
- **PARTIAL:** T2 <1e-4 AND T7 OK AND kk1204 ∈ [1e-4, 1e-3] → close epic, file Epic-E for kk1204 follow-up.
- **BLOCKED:** T2 ≥1e-4 even with composed pipeline → SPEC_FAILURE per discipline §9; file Epic-F (alternative path).

## Honest gates

- T2 stepstone <1e-4 HARD (close 13% gap from current 2.79e-4 master baseline).
- T7 K=4 well-conditioned `n_projected_dims == 0` HARD (over-projection regression guard).
- T8 kk1204 K=20 severe-skew <1e-3 PARTIAL (master diverges; PARTIAL = strict improvement).

Full plan: `docs/superpowers/plans/2026-05-02-newton-kl-tsvd-cg-plan.md` (rev 2).
