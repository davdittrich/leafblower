# Active Plan — ORIS/OT-family over-relaxation optimizations
<!-- approved: 2026-05-31 -->
<!-- design-review-gate: n/a (ticket-creation task, not a design doc) -->
<!-- plan-review-gate: PASSED 3/3 on iteration 2 (round-1 Completeness blocker resolved) -->
<!-- user-approved: pending -->
<!-- status: backlog -->

Impact-ranked over-relaxation optimizations, grounded in docs/methods/oris.md canonical
sources (Thibault 2021, Lehmann 2022, Soma-Uschmajew 2024). Two epics by approach; each
task carries an experiment gate (build + devtools::test + python parity rtol=1e-6 +
stepstone no-regression, plus a benchmark pass/fail or GO/NO-GO verdict).

## Epic leafblower-mj1p — ORIS-internal over-relaxation (P1)
- leafblower-mj1p.1 (P1, Rank1): unlock omega>1 — wire existing types.hpp omega_max into
  the recovery clamp (src/oris.cpp:837 SRAA, :1564 flat; currently std::min(1.0,...)).
  HIGH impact / LOW cost. The headline fix: the solver named "over-relaxed" never over-relaxes.
- leafblower-mj1p.2 (P2, Rank2, DEPENDS-ON mj1p.1): spectral optimal-omega from the
  successive-residual ratio (Lehmann omega_opt=2/(1+sqrt(1-theta2))) replacing the
  0.7/1.05 heuristic. NO-GO if it doesn't beat fixed omega_max.

## Epic leafblower-e65t — OT-family transfer (P2, research-gated)
- leafblower-e65t.1 (P2, Rank3): Anderson(SRAA) x over-relaxation interplay — 4-arm study;
  currently SOR adaptation is disabled on AA-accepted steps (oris.cpp:822). GO/NO-GO.
- leafblower-e65t.2 (P3, Rank4): over-relaxed greedy Greenkhorn — research probe; greedy
  breaks the full-sweep Lyapunov proof; NO-GO is the expected base case.
- leafblower-e65t.3 (P3, Rank5): over-relaxed Dykstra-Sinkhorn feasibility probe — theory
  gate first (Dykstra statefulness blocks direct transfer); likely documented NO-GO.

## Global guards (all tickets)
Fixed point unchanged (weights within 1e-8 of baseline); net exponent alpha*eff_omega in
(0,2); enum values/bounds semantics/public API frozen; no default behavior change without
a GO verdict + stepstone no-regression. Full hermetic bodies in beads.
