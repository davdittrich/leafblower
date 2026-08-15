---
phase: 04-truthful-surface
verified: 2026-08-15T18:30:00Z
status: passed
score: 4/4 must-haves verified
behavior_unverified: 0
overrides_applied: 0
---

# Phase 4: Truthful Surface Verification Report

**Phase Goal:** Nothing documented is untrue, and nothing the API silently swallows stays
silent — a reader can trust the docs and a user cannot get a plausible wrong answer by typo.
**Verified:** 2026-08-15T18:30:00Z
**Status:** passed
**Re-verification:** No — initial verification

## Goal Achievement

### Observable Truths (ROADMAP Success Criteria)

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| SC1 | No live document attributes the unimplemented outer entropic-proximal-point loop, or any other unshipped capability, to ORIS or to a removed solver | ✓ VERIFIED | `docs/raking.md`: `grep -c "renamed from iEPPA"` = 0, `grep -ci Gurobi` = 0, `grep -c "definitive state-of-the-art"` = 0, `grep -c "8.2 The ORIS Solver"` = 0, `grep -c "12. Synthesis and Final Conclusions"` = 0. Sentinel lines confirming surrounding content untouched: `"strictly respecting the capacity constraints"` = 1, `"external hardware accelerators"` = 1 (now file's last sentence, `tail -c 300` confirms). `git show 3143594` matches the two-passage deletion described in the plan/summary, zero replacement prose inserted. |
| SC2 | A developer reading `rk_algorithm_t` sees why slot 7 is a hole, as they already do for slot 2, and cannot reuse either value by accident | ✓ VERIFIED | `src/leafblower.h:49`: `/* 7 = removed (was RK_ALG_GRAKE) */` present between `RK_ALG_GREG = 6` and `RK_ALG_ORIS_SOFT = 8`; enum values unchanged (comment-only, matches slot-2's pre-existing `/* 2 = removed (was RK_ALG_LBFGSB) */` convention). `CLAUDE.md:66`: "Algorithm slots 2 and 7 are reserved (LBFGSB and GRAKE removed, respectively)." — both slots named in one sentence. |
| SC3 | `harvest(..., weights = w)` raises an informative error naming `design_weights=` instead of silently ignoring the argument and returning a plausible unweighted result | ✓ VERIFIED | `R/harvest.R:309-312`: `if ("weights" %in% names(dots)) stop("leafblower: unrecognized argument 'weights' — did you mean 'design_weights'? ...")` executes ahead of the generic RVAL.2 `warning()` at line 313-315 (confirmed by source order and by the `RVAL.4` regression test in `tests/testthat/test-harvest-rval.R:153-162`, which asserts `expect_error(..., regexp = "design_weights")`). Re-ran (not trusted from SUMMARY): `testthat::test_dir(filter='harvest-rval')` → `FAIL 0 \| WARN 1 \| SKIP 0 \| PASS 11` (the 1 warning is the pre-existing unrelated RVAL.2 sparse-category warning at line 26, not from the new guard). `eb79.18` in `test-logit.R:195` renamed to `design_weights = base_w`; re-ran `testthat::test_dir(filter='logit')` → `FAIL 0 \| WARN 0 \| SKIP 0 \| PASS 54`. `NEWS.md:67-68` records the breaking change under a `## Breaking changes` heading with the literal `design_weights=` string. |
| SC4 | An audit of `README`, `NEWS.md`, `man/` and `docs/` finds no surviving reference to `grake`, `lbfgsb` or `cp` as available methods | ✓ VERIFIED | Re-ran (not trusted from SUMMARY) case-insensitive `grake\|lbfgsb\|l-bfgs-b` sweep across `README.md`, `NEWS.md`, `man/*.Rd`, `docs/methods/*.md`, `docs/raking.md`, `docs/performance.md`. All 17 hits classify into the six documented false-positive classes: `grake_norm` (live convergence-metric field, `NEWS.md`/`man/harvest.Rd`), `survey::grake`/`grake.R` (comparative reference to the R `survey` package's own function, `docs/methods/logit.md`), `L-BFGS-B` correctly-documented-removed (`docs/methods/00-overview.md:3,81`), 9 `L-BFGS-B` mentions confined to `docs/raking.md` §6.2/§11 (general-technique discussion, explicitly out of D-01/D-02 scope — count matches the plan's predicted 9 exactly), and `regrake`/`rswjax` (third-party package names, `docs/methods/raking.md`). Word-boundary `\bcp\b` sweep: one hit, `docs/performance.md:176`, the withdrawn `ylsy-cp-ipm-spike-result.md` research-spike report path — pre-classified false positive, correctly framed as rejected. No new true-positive hit found. README.md: zero hits (clean). |

**Score:** 4/4 truths verified (0 present, behavior-unverified)

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `R/harvest.R` | RVAL.2 gains a `stop()` guard naming `design_weights=` | ✓ VERIFIED | Guard present, wired ahead of generic warning, exercised by `RVAL.4` test |
| `tests/testthat/test-harvest-rval.R` | New `RVAL.4` regression test | ✓ VERIFIED | Present at line 153, passes |
| `tests/testthat/test-logit.R` | `eb79.18` renamed to `design_weights=` | ✓ VERIFIED | Line 195, passes |
| `NEWS.md` | Breaking-change bullet documenting the guard | ✓ VERIFIED | Line 67-68, under `## Breaking changes` |
| `docs/raking.md` | §8.2/§12 misattribution deleted | ✓ VERIFIED | All acceptance greps return 0; sentinels intact |
| `src/leafblower.h` | Slot 7 comment-annotated | ✓ VERIFIED | Comment-only insertion, enum values unchanged |
| `CLAUDE.md` | Slot 2 + slot 7 named in one sentence | ✓ VERIFIED | Line 66; unrelated deletion (WR-01) fixed by commit 5200621 |

### Key Link Verification

| From | To | Via | Status | Details |
|------|-----|-----|--------|---------|
| `R/harvest.R` RVAL.2 dots-check | `stop()` guard | Guard `if` clause precedes generic `warning()` in source order | ✓ WIRED | Lines 309-312 execute before 313-315; confirmed by both source read and RVAL.4/eb79.18 test outcomes |
| `src/leafblower.h` enum comment | R + Python build sites | Both build sites compile against this single header | ✓ WIRED | `R CMD INSTALL --preclean .` rebuilt cleanly; full R suite green; Python parity suite (which forces a rebuild via `uv pip install -e . --reinstall-package leafblower`) green |

### Behavioral Spot-Checks

| Behavior | Command | Result | Status |
|----------|---------|--------|--------|
| SC3: bare `weights=` errors naming `design_weights` | `testthat::test_dir(filter='harvest-rval')` (single named-filter run, not full suite) | `FAIL 0 \| WARN 1 \| SKIP 0 \| PASS 11` | ✓ PASS |
| SC3: `eb79.18` still converges under `design_weights=` | `testthat::test_dir(filter='logit')` | `FAIL 0 \| WARN 0 \| SKIP 0 \| PASS 54` | ✓ PASS |
| DoD gate: full R suite (run once) | `Rscript -e "devtools::test()"` | `FAIL 0 \| WARN 141 \| SKIP 13 \| PASS 1837` | ✓ PASS |
| DoD gate: R build | `R CMD INSTALL --preclean .` | `* DONE (leafblower)` | ✓ PASS |
| DoD gate: Python parity (single-thread BLAS, run once) | `.venv/bin/python -m pytest -q` | `160 passed, 0 failed` | ✓ PASS |

All five numbers match SUMMARY.md's claimed figures exactly — independently re-executed, not
taken on trust.

### Requirements Coverage

No requirement IDs mapped to Phase 4 (`.planning/ROADMAP.md:196` states "Requirements: (none —
defect-driven)"). `REQUIREMENTS.md:271` lists Phase 4's beads tickets (`leafblower-05ha`,
`leafblower-x2iq`, `weights=` guard) as the tracking mechanism, not formal requirement IDs. No
orphaned requirements found.

### Anti-Patterns Found

None. Scanned all seven files touched across both plans (`R/harvest.R`,
`tests/testthat/test-harvest-rval.R`, `tests/testthat/test-logit.R`, `NEWS.md`,
`docs/raking.md`, `src/leafblower.h`, `CLAUDE.md`) for TBD/FIXME/XXX/TODO/HACK/PLACEHOLDER
markers and stub patterns — zero hits in the diffs introduced by this phase.

### Code Review Follow-up (04-REVIEW.md)

`04-REVIEW.md` (standing review, not re-litigated here) found one WARNING (WR-01: commit
`2a06f1f` silently deleted an unrelated `bd remember`/MEMORY.md governance bullet from
`CLAUDE.md` while editing the adjacent slot-7 sentence) and one INFO (IN-01: the `weights=`
guard is exact-match only, does not cover near-typos like `weight=`/`wt=` — explicitly scoped
out of this phase per the ticket, not a gap).

**WR-01 fix confirmed:** commit `5200621` ("fix(claude-md): restore accidentally-deleted bd
remember governance line") restores the exact deleted line (`git show 5200621` — 1 file
changed, 1 insertion, 0 deletions). Post-fix state verified: `grep -n "bd remember\|MEMORY.md"
CLAUDE.md` → line 27, the restored bullet. The slot-2/7 sentence added by `2a06f1f` remains
intact at line 66 — the fix commit touched only the WR-01 line, no other regression introduced.

**IN-01:** no fix required — the review itself classifies this as out-of-scope residual
coverage, not a phase-4 gap (D-04's decision scope was the one confirmed in-repo caller using
the exact string `weights=`).

### Human Verification Required

None. All four Success Criteria are grep/test-verifiable; the guard's error-raising behavior is
exercised by an automated regression test (`RVAL.4`), not merely present-and-wired.

### Gaps Summary

None. All four ROADMAP Success Criteria verified against the live codebase (not SUMMARY.md
claims): guard behavior independently re-tested, doc-deletion boundaries independently
re-greped with sentinel checks, enum annotation independently read, SC4 audit independently
re-run across all seven in-scope file globs including a fresh word-boundary `cp` check. The
one code-review WARNING (WR-01) has a confirmed, minimal, correctly-scoped fix already
committed. Full DoD gate (R build, R test suite, Python parity suite) re-run once each in this
verification pass — 0 FAIL in R (1837 pass), 0 failed in Python (160 pass) — matching
SUMMARY.md's claimed figures exactly.

---

_Verified: 2026-08-15T18:30:00Z_
_Verifier: Claude (gsd-verifier)_
