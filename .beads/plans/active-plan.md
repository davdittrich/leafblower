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

- [ ] WU-1 — iEPPA homotopy/priority/eta config scaffolding (scaffolding; no behavioural change)
- [ ] WU-2 — errRp trajectory probe via env vars (internal diagnostic)
- [ ] WU-2.5 — falsify simpler-alternative hypothesis on stepstone-fulldata
- [ ] WU-3 — P-A progressive max_weight homotopy outer loop
- [ ] WU-4 — P-B Greenkhorn priority-ordered margin scheduler
- [ ] WU-5 — Tang 2024 dynamic-eta schedule (schedule-only borrow)
- [ ] WU-6 — stepstone-fulldata merge gate + kk1204 non-regression
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
