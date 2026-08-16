---
phase: 03-honest-performance-gate
verified: 2026-08-16T01:02:07Z
status: passed
score: 6/6 must-haves verified
behavior_unverified: 0
overrides_applied: 0
re_verification:
  previous_status: passed (2026-08-15, pre-gap-closure — superseded by UAT gaps G-03-1/G-03-4 opened after that verification ran)
  previous_score: 4/4
  gaps_closed:
    - "G-03-1: per-method benchmark coverage — every leafblower method (oris, oris_soft, raking, newton_kl, greg, logit, chebyshev, sinkhorn, greenkhorn) now has a measured row against its own doc-named closest competitor"
    - "G-03-4: docs/performance.md's Competitors section now covers all 9 methods with per-method mapping, not just the original 3-package oris_soft comparison"
  gaps_remaining: []
  regressions: []
---

# Phase 3: Honest Performance Gate Verification Report

**Phase Goal:** Publish honest, benchmark-grounded performance claims for every leafblower
method against a doc-named closest competitor, with no silently dropped scope caveats or
bound violations.
**Verified:** 2026-08-16T01:02:07Z
**Status:** passed
**Re-verification:** Yes — after gap closure (03-05..08) and a subsequent code review
(03-REVIEW.md) whose two Warning fixes (WR-01, WR-02) landed in commits `af6913e`/`c4b511a`.

This verification supersedes `03-VERIFICATION.md`'s 2026-08-15 run, which certified only the
original 4 plans (03-01..04, scoped to `oris_soft` vs 3 generic survey packages). UAT
(`03-UAT.md`) subsequently opened gaps G-03-1/G-03-4 because that scope excluded 6 of
leafblower's 9 methods. Gap-closure plans 03-05..08 and a follow-up code review with two
just-landed fixes are what this run independently re-verifies, per the request scope.

## Goal Achievement

### Observable Truths

| # | Truth | Status | Evidence |
|---|---|---|---|
| 1 | Every leafblower method (`oris`, `oris_soft`, `raking`, `newton_kl`, `greg`, `logit`, `chebyshev`, `sinkhorn`, `greenkhorn`) has a fresh, measured row against its own doc-named closest competitor | ✓ VERIFIED | `docs/performance.md`'s Competitors table (lines 259-291) lists all 9 methods with per-method mapping (`survey::calibrate`/`icarus`/`ReGenesees` for the raking family, distance-matched `survey::calibrate(calfun=)` for greg/logit, `optweight::optweight.svy(norm='linf')` for chebyshev, POT's `ot.bregman.greenkhorn`/`ot.sinkhorn` for greenkhorn/sinkhorn), each citing the source `docs/methods/*.md` "Practitioner implementations" section. Cross-checked against live-regenerated `benchmarks/results/oris_soft_vs_competitors.csv` (13 rows on `medium_100k_5margins` alone) and `benchmarks/results/greenkhorn_sinkhorn_vs_pot.csv` (4 rows) — every figure transcribed verbatim. |
| 2 | Bound violations are never silently reported `ok=TRUE` — a solver's own self-reported status is independently re-graded against `max_weight` | ✓ VERIFIED (behaviorally, live re-run) | Ran `OMP/OPENBLAS/MKL_NUM_THREADS=1 Rscript benchmarks/oris_soft_vs_competitors.R` live during this verification. Regenerated CSV: `leafblower_newton_kl` reports `ok=FALSE` with `max_w=6.2105 > max_weight=3`, `note` field states "max_w=6.2105 exceeds max_weight=3 -- newton_kl reports rather than clamping bound violations ... see leafblower-73d7". Every other bounded arm (`leafblower_oris_soft`, `survey_calibrate`, `icarus_calibration`, `ReGenesees_e_calibrate`, `leafblower_oris`, `leafblower_raking`, `leafblower_greg`, `leafblower_logit`, both `survey_calibrate_linear/logit`, `leafblower_chebyshev`) correctly reports `ok=TRUE` at or within bound. This is WR-01's fix (`af6913e`), confirmed landed in source (`oris_soft_vs_competitors.R:386-393`, `ok_lb <- isTRUE(res_lb$status %in% c(0L,5L)) && max(w_lb_n) <= 3 + 1e-6`) and confirmed correct at runtime, not just present. |
| 3 | Bound-compliance checks (`ok_*`) and reported `max_w`/`min_w` are computed on the same weight scale, not silently mismatched | ✓ VERIFIED | `grep -n "w_sv_n" benchmarks/oris_soft_vs_competitors.R` shows both `survey_calibrate` call sites (pre-existing block, line ~173; new `survey_calfun_arm_row()`, line ~440) now compute `ok_sv <- all(is.finite(w_sv_n)) && max(w_sv_n) <= max_weight + 1e-6` on the normalized `w_sv_n` vector — the same vector reported as `max_w`/`min_w` in the same row. This is WR-02's fix (`c4b511a`), confirmed landed at both call sites (previously only one site checked the raw, unnormalized `w_sv`). |
| 4 | Known scope caveats (objective mismatch, missing-bound-argument competitors, cross-language timer methodology, K=2-only equivalence) are stated on the published page, not left silently dropped | ✓ VERIFIED | `docs/performance.md` Competitors table's `optweight_linf` row states the weight-deviation-vs-margin-error objective mismatch and the missing `max.w` argument verbatim from the CSV note; the `sinkhorn`/`greenkhorn` rows state POT's lack of a bounds mechanism and the K=2-only equivalence scope, also verbatim from CSV notes. The `known_limit_k20_uniform` section is unchanged from the original verification (structural limit, not silently dropped). |
| 5 | The one-command reproduction (`benchmarks/run_honest_gate.sh`) covers every published figure, including the newly-added POT comparison | ✓ VERIFIED | `run_honest_gate.sh` (last touched `78103f8`, 03-07) exports single-thread BLAS + `LBW_BENCH_GATE`/`NOT_CRAN`, then runs the stepstone regression gate, `oris_soft_vs_competitors.R`, and `python/.venv/bin/python benchmarks/greenkhorn_sinkhorn_vs_pot.py` in sequence — `python/.venv/bin/python` confirmed to exist on disk. |
| 6 | WR-03 (stale `1e-10` tolerance comment) is genuinely cosmetic, not masking a real gap | ✓ VERIFIED | `oris_soft_vs_competitors.R:65-68`'s comment on `normalize_to_n()` still misattributes a `1e-10` tolerance to "this script itself" — confirmed unfixed, matching 03-REVIEW.md's description exactly. Independently confirmed every actual `ok_*` bound check in the file (`grep -n "+ 1e-6"` at 4 call sites: lines ~173, ~228 icarus, ~293 ReGenesees, ~386 lb_only_arm_row, ~440 survey_calfun_arm_row — all consistently `+ 1e-6`) uses the same tolerance; the `1e-10` figure genuinely belongs to a different file (`tests/testthat/test-bench-gate.R:45`), which does independently enforce its own `1e-10` gate on the `oris_soft` row only. The comment error does not change any computed `ok`/`max_w` value in this script — it is a documentation-accuracy defect in a comment, not a logic defect. |

**Score:** 6/6 truths verified (0 present, behavior-unverified)

### Required Artifacts

| Artifact | Expected | Status | Details |
|---|---|---|---|
| `benchmarks/oris_soft_vs_competitors.R` | 9-method-capable measurement script with independently re-graded `ok` flags | ✓ VERIFIED | 13 arms on `medium_100k_5margins` alone (raking family x3 + competitors x3, greg/logit x2 + distance-matched competitors x2, chebyshev + optweight, plus known-limit rows). `lb_only_arm_row()`/`survey_calfun_arm_row()` helpers confirmed present and correctly re-grading bound compliance post-fix. |
| `benchmarks/greenkhorn_sinkhorn_vs_pot.py` | New file, K=2 POT-equivalence measurement | ✓ VERIFIED | Exists (13,230 bytes), CSV output confirmed substantive (4 rows, real numeric equivalence to ~15 significant figures across leafblower/POT arm pairs). |
| `benchmarks/results/oris_soft_vs_competitors.csv` | 13+ rows spanning 9 methods | ✓ VERIFIED | Regenerated live during this verification; content matches every figure in `docs/performance.md`. |
| `benchmarks/results/greenkhorn_sinkhorn_vs_pot.csv` | 4 rows (2 leafblower arms x 2 POT arms) | ✓ VERIFIED | Present, content matches `docs/performance.md`'s `k2_margin_pot_equiv` table. |
| `docs/performance.md` | Full 9-method Competitors table + extended Results | ✓ VERIFIED | 298 lines; "per-method competitor coverage" table (lines 69-104) and Competitors section (lines 259-291) both present and populated for all 9 methods, transcribed verbatim from the CSVs. |
| `benchmarks/run_honest_gate.sh` | One command reproduces every figure | ✓ VERIFIED | Wired step for the POT comparison added (03-07); executable, all three steps present. |

### Key Link Verification

| From | To | Via | Status | Details |
|---|---|---|---|---|
| `run_honest_gate.sh` | `greenkhorn_sinkhorn_vs_pot.py` | invocation | ✓ WIRED | Confirmed in script body (`python/.venv/bin/python benchmarks/greenkhorn_sinkhorn_vs_pot.py`), venv binary confirmed to exist. |
| `docs/performance.md` Competitors table | `benchmarks/results/*.csv` note fields | verbatim transcription | ✓ WIRED | Every caveat sentence in the Competitors table (optweight objective mismatch, POT no-bounds, newton_kl bound violation) matches a live-read CSV `note` field byte-for-byte on the load-bearing clauses. |
| `lb_only_arm_row()`'s `ok_lb` | `max_weight` bound | independent re-grading | ✓ WIRED (behaviorally confirmed live) | Live-run CSV shows `ok_lb` diverges from `status`-only truth exactly where expected (`newton_kl`), proving the check is load-bearing, not a no-op. |
| `survey_calfun_arm_row()`'s `ok_sv` | normalized weight vector `w_sv_n` | scale consistency | ✓ WIRED | Both call sites read from the same `w_sv_n` used for the reported `max_w`/`min_w`. |

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|---|---|---|---|
| WR-01 fix, live | `OMP/OPENBLAS/MKL_NUM_THREADS=1 Rscript benchmarks/oris_soft_vs_competitors.R` | `leafblower_newton_kl` row: `ok=FALSE`, `max_w=6.2105 > max_weight=3`; all other bounded arms `ok=TRUE` | ✓ PASS |
| Determinism | git diff of regenerated CSV vs. committed CSV after live re-run | No diff | ✓ PASS |
| WR-02 fix, static | `grep -n "w_sv_n" benchmarks/oris_soft_vs_competitors.R` | Both `survey_calibrate` call sites compute `ok_sv` on `w_sv_n` | ✓ PASS |
| Commits touched no package source | `git show --stat af6913e c4b511a \| grep -E "^ (R/\|src/\|python/leafblower/)"` | no matches | ✓ PASS (confirms no DoD re-run needed for these two fix commits; prior verification's `R CMD INSTALL`/`devtools::test()`/pytest 0-FAIL results stand unaffected) |
| Anti-pattern scan | `grep -nE "TBD\|FIXME\|XXX\|TODO\|HACK\|PLACEHOLDER"` on all 4 phase-modified files | 0 matches | ✓ PASS |

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|---|---|---|---|---|
| US-003 | 03-01..03-08 | 1M+ obs / 20+ margins < 30s | ✓ SATISFIED (as "Partial," honestly) | `.planning/REQUIREMENTS.md:110-131` names live artefacts and states plainly no single fixture combines "1M+ rows AND 20+ margins" — this remains true after gap closure (which added per-method competitor coverage, not a new large+20-margin fixture) and was not required to change by 03-05..08's own `must_haves` (none listed REQUIREMENTS.md as an artifact). Traceability table (line 272) correctly still shows `Partial`. |
| KPI-04 | 03-01..03-08 | Large-scale 1M rows/20 margins/<30s | ✓ SATISFIED (as "Partial," honestly) | Same reasoning as US-003; `.planning/REQUIREMENTS.md:174-` and traceability line 282 (`Partial`) unaffected and still accurate. |

No orphaned requirements found for Phase 3 in `.planning/REQUIREMENTS.md`.

### Anti-Patterns Found

None blocking. No `TBD`/`FIXME`/`XXX`/`TODO`/`HACK`/`PLACEHOLDER` markers in any of the 4 files reviewed by `03-REVIEW.md` or touched by the WR-01/WR-02 fix commits. WR-03 (stale `1e-10` comment) is a documentation-accuracy defect confirmed cosmetic — see Truth 6 above — left unfixed by deliberate choice, consistent with the review's own "Fix: [optional]" framing (WR-03 has a `Fix:` suggestion but no severity that blocks the phase; it was not marked blocking in `03-REVIEW.md`'s frontmatter, which records `critical: 0`).

### Human Verification Required

None. All six truths are independently verifiable via file content, CSV/doc cross-referencing, and one live re-execution of the benchmark script under the mandated single-thread BLAS envelope, all performed during this verification.

### Gaps Summary

None. Both UAT gaps (G-03-1, G-03-4) are closed: `docs/performance.md` and the two regenerated CSVs now cover all 9 leafblower methods, each against its own doc-named closest competitor, with scope caveats (objective mismatch, missing bound arguments, K=2-only equivalence, cross-language timer methodology) carried into the published page rather than left only in CSV `note` fields. Both Warning-severity code-review findings (WR-01: self-reported-status-only `ok` flag silently misrepresenting a visible bound violation; WR-02: `ok` computed on a different weight scale than the reported `max_w`/`min_w`) are confirmed landed in source and behaviorally correct via a live re-run, not merely claimed by commit messages. WR-03 (stale tolerance comment) is confirmed genuinely cosmetic — it does not affect any computed value, and every actual `ok_*` check in the file uses a consistent, correct `1e-6` tolerance.

---

*Verified: 2026-08-16T01:02:07Z*
*Verifier: Claude (gsd-verifier)*
