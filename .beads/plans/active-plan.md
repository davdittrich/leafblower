# Active Plan
<!-- approved: 2026-04-24T14:59:26+02:00 -->
<!-- gate-iterations: 2 (iter1 FAIL all 3, iter2 PASS all 3) -->
<!-- user-approved: true -->
<!-- status: in-progress -->
<!-- execution-mode: subagent-driven -->
<!-- previous-plan: docs/superpowers/plans/2026-04-23-ieppa-faithful-impl.md (status: complete 2026-04-23) -->

See: `docs/superpowers/plans/2026-04-24-ieppa-homotopy-greenkhorn.md` (rev 2)

Research basis: `docs/investigations/2026-04-24-ieppa-accel-research.md`

## WU Status

- [x] WU-1 (leafblower-d88d) — iEPPA homotopy/priority/eta config scaffolding — commits e4694bc + 6c6e286 + ed50889 — spec ✅ + quality ✅
- [x] WU-2 (leafblower-82oe) — errRp trajectory probe via env vars — commits cf3c36a + 0954f8b — spec ✅ + quality ✅
- [x] WU-2.5 (leafblower-xb81) — falsify simpler alternative — commit cc3c26b — max_err=0.165 at 10k iters, DIVERGES → overlays justified ✅
- [x] WU-3 (leafblower-s6hi) — P-A homotopy outer loop — commits cb51974 + 6b5676d + 0b079c4 — 56.9% errRp reduction on tight-clamp synthetic — spec ✅ + quality ✅
- [x] WU-4 (leafblower-sa62) — P-B Greenkhorn scheduler — commits ccb0488 + 9bf6b5e — 40% sweep savings on 2-cat-heavy fixture — spec ✅ + quality ✅
- [x] WU-5 (leafblower-u7u4) — Tang dynamic-eta (damping via beta) — commits b6d2c71 + cbf5d0c + 5338fd5 — 25.3% errRp reduction at eta_start=20 — spec ✅ + quality ✅
- [x] WU-6 (leafblower-aa9b) — merge gate — commit d9dc5f9 — GATE FAILS: AB errRp=6.567e-3 (threshold 1.60e-3), slope=-0.285 (threshold -0.75), Pearson=0.977 (threshold 0.99). Homotopy+greedy DEGRADES vs baseline 2.223e-3. kk1204 non-regression PASS. HUMAN ticket filed.
- [ ] WU-7 — Python parity + roxygen + README

## Merge gate

1. `devtools::test()` green
2. `pytest` green
3. `R CMD check --as-cran` — 0 ERROR, 0 WARNING, NOTEs ≤ baseline
4. stepstone-fulldata AB config `max_err ≤ 1.60e-3`
5. kk1204 non-regression `max_err ≤ 1.322e-3 @ max_iterations=500`
6. Pearson r ≥ 0.99 vs commit-8146894 reference

## Reported-only (non-gating)

- Wall time AB vs autumn
- Stretch: ABE `max_err ≤ 1.60e-4`
- Rate slope on stepstone-small trajectory (gate at ≤ -0.75 per WU-6 Step 4)
- Degenerate M_cell=n ratio vs raking
