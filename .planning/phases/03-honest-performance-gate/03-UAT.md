---
status: complete
phase: 03-honest-performance-gate
source: [03-01-SUMMARY.md, 03-02-SUMMARY.md, 03-03-SUMMARY.md, 03-04-SUMMARY.md]
started: 2026-08-15T21:31:38Z
updated: 2026-08-15T21:45:00Z
---

## Current Test

[testing complete]

## Tests

### 1. D1 (03-01) — medium_100k_5margins headline numbers
expected: leafblower_oris_soft (canonical marginal_kl/improvement convergence) vs survey_calibrate (canonical epsilon=1e-3), identical per-observation bounds, identically-computed max_error/max_w/n_eff — measured fresh, not asserted.
result: issue
reported: "the benchmark scope was interpreted too narrow. I want to see for EACH method in leafblower how they compare to their closest competitors each."
severity: major

### 2. D2 (03-02) — known_limit_k20_uniform substitution
expected: m_cell_over_n=1.0000 measured (every observation its own cell), oris_soft and raking_accelerated both report finite max_error above 1e-6 under a bounded 500-iteration budget — confirms the retired <30s AND <1e-6 gate is structurally unachievable, backed by a number not a citation. Fixture substitution (RESEARCH.md Pitfall 2 divergence) confirmed as correct resolution.
result: pass

### 3. D1 (03-03) — published threshold decision
expected: "paired" framing selected (speed never asserted without accuracy/bound/n_eff), wall_s<=0.5s ceiling and n_eff>=60000 floor with derivations, decided before any code was touched.
result: pass

### 4. D3 (03-03) — docs/performance.md read
expected: headline claim (paired framing), per-class results tables, verbatim machine/version provenance, methodology, one-command reproduction, known-limit section cross-linked without restating, competitors section grounded in docs/methods/oris.md — every figure traceable, no detached speed number.
result: issue
reported: "the comparison of full set of methods vs their closest competitor each is missing"
severity: major

### 5. D5 (03-04) — Definition-of-Done before/after count sanity
expected: devtools::test() SKIP 12→13 (+1, kk1204 now skips by default), PASS 1839→1836 (-3, the three kk1204 assertions no longer running under a plain run) — reads as "the gate correctly stopped running by default," not a silent break. Python parity 160/0 unaffected.
result: pass

### 6. D2 (03-01) — honest-gate wrapper + regeneration
expected: benchmarks/run_honest_gate.sh regenerates comparison artifacts under single-thread BLAS.
result: pass
source: automated
coverage_id: D2

### 7. D3 (03-01) — testthat honest-gate assertion + negative control
expected: LBW_BENCH_GATE=1 assertion reads the CSV and fails on breach; unset it skips.
result: pass
source: automated
coverage_id: D3

### 8. D1 (03-02) — large_stepstone_fulldata measured fresh
expected: leafblower_oris_soft measured on the 1.58M-row real-survey fixture, competitor infeasibility recorded as computed, not omitted.
result: pass
source: automated
coverage_id: D1

### 9. D3 (03-02) — full D-07 competitor set on medium class
expected: oris_soft, survey_calibrate, icarus_calibration, ReGenesees_e_calibrate all graded by the same metrics; neither icarus nor ReGenesees added to DESCRIPTION/pyproject.toml.
result: pass
source: automated
coverage_id: D3

### 10. D4 (03-02) — Definition-of-Done unbroken
expected: R CMD INSTALL + devtools::test() 0 FAIL; plan 01's wrapper/assertion still exit 0; DESCRIPTION/pyproject.toml diff empty.
result: pass
source: automated
coverage_id: D4

### 11. D2 (03-03) — honest-gate wall_s/n_eff assertions + negative controls
expected: wall_s<=0.5s and n_eff>=60000 assertions added with derivations; negative controls confirm they actually fire; gate stays opt-in.
result: pass
source: automated
coverage_id: D2

### 12. D4 (03-03) — Definition-of-Done unbroken
expected: R CMD INSTALL + devtools::test() unchanged from 03-02 baseline; run_honest_gate.sh still exits 0.
result: pass
source: automated
coverage_id: D4

### 13. D1 (03-04) — README.md published
expected: package name, orienting sentence, paired headline claim, link to docs/performance.md, <=25 lines, no Phase 5 scope, no autumn mention.
result: pass
source: automated
coverage_id: D1

### 14. D2 (03-04) — PRD contradiction sites marked superseded
expected: three "100K rows, 5 margins" sites in tasks/prd-leafblower-core.md marked superseded in place; void test-lbfgsb.R reference annotated, not deleted.
result: pass
source: automated
coverage_id: D2

### 15. D3 (03-04) — REQUIREMENTS.md rewritten to match what Phase 3 measured
expected: US-003/KPI-04 entries name live artifacts, restate K=20 finding as documented known limit.
result: pass
source: automated
coverage_id: D3

### 16. D4 (03-04) — kk1204 block re-gated
expected: skip_if(CI) replaced with skip_if(LBW_BENCH_GATE=='') gating; description/comment rewritten to describe what it actually measures.
result: pass
source: automated
coverage_id: D4

## Summary

total: 16
passed: 13
issues: 2
pending: 0
skipped: 0
auto_passed: 11

## Gaps

- gap_id: G-03-1
  truth: "leafblower_oris_soft benchmarked against its closest generic competitors (survey_calibrate, icarus_calibration, ReGenesees_e_calibrate) on medium/large fixtures"
  status: failed
  reason: "User reported: the benchmark scope was interpreted too narrow. I want to see for EACH method in leafblower how they compare to their closest competitors each."
  severity: major
  test: 1
  artifacts: []
  missing:
    - "Per-method benchmark coverage: leafblower exposes method= {raking, greenkhorn, sinkhorn, oris, oris_soft, newton_kl, chebyshev} (R/harvest.R:273 default + docs/methods/*.md), but benchmarks/oris_soft_vs_competitors.R and its CSV only exercise oris_soft (3 fixtures) and raking_accelerated (1 fixture, known-limit only). greenkhorn, sinkhorn, plain oris, newton_kl, chebyshev have zero competitor-comparison rows."
    - "Closest-competitor mapping per method not yet established: e.g. greenkhorn's closest competitor is likely POT's ot.bregman.greenkhorn / ot.sinkhorn(method='greenkhorn') per docs/methods/greenkhorn.md's own citation (flamary2021pot), not survey::calibrate — each method doc may already name its own natural competitor and needs to be read to build the mapping before benchmarking."

- gap_id: G-03-4
  truth: "docs/performance.md's competitors section covers every leafblower method against its own closest competitor"
  status: failed
  reason: "User reported: the comparison of full set of methods vs their closest competitor each is missing"
  severity: major
  test: 4
  artifacts: []
  missing:
    - "Same root cause as G-03-1 (per-method benchmark coverage) — docs/performance.md can only publish a per-method comparison table once the underlying per-method benchmark data (G-03-1) exists. Resolve G-03-1 first; this is its documentation-side consequence."
