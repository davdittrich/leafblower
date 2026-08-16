---
phase: 03-honest-performance-gate
reviewed: 2026-08-16T02:13:58+02:00
depth: standard
files_reviewed: 4
files_reviewed_list:
  - benchmarks/oris_soft_vs_competitors.R
  - benchmarks/greenkhorn_sinkhorn_vs_pot.py
  - benchmarks/run_honest_gate.sh
  - docs/performance.md
findings:
  critical: 0
  warning: 3
  info: 2
  total: 5
status: issues_found
---

# Phase 03: Honest Performance Gate — Code Review Report (plans 03-05..03-08)

**Reviewed:** 2026-08-16T02:13:58+02:00
**Depth:** standard
**Files Reviewed:** 4
**Status:** issues_found

## Summary

Reviewed the gap-closure additions from plans 03-05 through 03-08:
`benchmarks/oris_soft_vs_competitors.R` (new `lb_only_arm_row()` /
`survey_calfun_arm_row()` helpers and the chebyshev/optweight arms),
`benchmarks/greenkhorn_sinkhorn_vs_pot.py` (new file), `benchmarks/run_honest_gate.sh`
(one added step), and `docs/performance.md` (Competitors section rewrite +
two new Results subsections).

Every published number I cross-checked against the two regenerated CSVs
(`benchmarks/results/oris_soft_vs_competitors.csv`,
`benchmarks/results/greenkhorn_sinkhorn_vs_pot.csv`) transcribed exactly —
all 9 rows of the per-method table and all 4 rows of the POT table match to
the displayed precision, and the doc's verbatim quotes from
`docs/methods/{oris,newton_kl,chebyshev,greenkhorn}.md` are faithful. No
hardcoded secrets, injection, or dangerous-function patterns found. No
security-relevant surface in these four files (local benchmark scripts and a
docs page).

The real defect is in the benchmark script's own honesty machinery: the new
`lb_only_arm_row()` helper (introduced in plan 03-05, reused by every
plan-05/06 arm) computes its `ok` column from the solver's self-reported
`status` alone, never independently checking the returned weights against
`max_weight` — exactly the "don't trust a package's own convergence claim"
failure mode the file's `Methodology` section (and `docs/performance.md`'s
own prose) explicitly calls out for competitor packages, but does not apply
to leafblower's own arms. This produces a `leafblower_newton_kl` CSV row
with `ok=TRUE` even though its `max_w=6.2105` visibly exceeds the fixture's
stated `max_weight=3` bound, and the caveat that would explain this to a
reader (newton_kl reports rather than clamps bound violations — flagged
explicitly by 03-05-SUMMARY.md as "worth flagging for the docs plan") never
made it into either the CSV `note` field or `docs/performance.md`'s
published table. Two lesser findings (a scale-mismatch between an `ok` check
and its reported `max_w`, and a stale comment) round out the list.

## Warnings

### WR-01: `leafblower_newton_kl`'s `ok=TRUE` in the CSV misrepresents a known, visible bound violation — and the docs page drops the caveat that would explain it

**File:** `benchmarks/oris_soft_vs_competitors.R:366-392` (new in plan 03-05, `lb_only_arm_row()`), consumed at `docs/performance.md:83` and `:91-94`

**Issue:** `lb_only_arm_row()` sets `ok_lb <- isTRUE(res_lb$status %in% c(0L, 5L))` (line 379) — it checks only the solver's own reported status code, never `max(w_lb_n) <= max_weight + <tol>` the way every competitor arm in the same file does (`icarus_calibration`'s `ok_ic`, `ReGenesees_e_calibrate`'s `ok_rg`, `survey_calibrate`'s `ok_sv`, `optweight_linf`'s explicit `ok=NA` — all of these independently re-grade the returned weight vector against the bound). This asymmetry means leafblower's own arms get exactly the "self-reported convergence taken at face value" treatment the script's header comment and `docs/performance.md`'s own Methodology section (*"A package that reports 'converged' is not taken at its word; its returned weights are re-graded independently."*) explicitly disclaim for competitors.

The concrete manifestation is in the regenerated CSV: the `leafblower_newton_kl` row (added by this plan) has `max_w=6.21049873735179` against a fixture `max_weight=3`, yet `ok=TRUE, status=0` — a solver that is 2x over its stated per-observation bound is reported as fully converged/ok. This is a documented, intentional *solver* behavior (`newton_kl` reports `RK_ERR_NOCONV` rather than clamping — 03-05-SUMMARY.md names `leafblower-73d7`/`PROJECT.md` as the source of that contract and explicitly flags it as "worth flagging for the docs plan"), but neither the CSV `note` field for this row nor `docs/performance.md`'s "Notes carried verbatim from the CSV's `note` field" paragraph (`docs/performance.md:91-94`, which only discusses `leafblower_greg` and `optweight_linf`) carries that caveat forward. A reader of the published table (`docs/performance.md:83`, `leafblower_newton_kl | 0.0492 | 4.220e-09 | 6.2105 | 0.1460 | 66188.7`) sees a `max_w` that visibly exceeds the stated `[0,3]` bound used by every other row on the same table, with zero explanation — the exact "speed claim detached from what it cost in accuracy" failure the page's own Headline-claim section says this phase exists to end.

**Fix:** Two independent fixes, either is sufficient, both are cheap:
1. Make `ok_lb` in `lb_only_arm_row()` (and the pre-existing `run_input_class` Arm-1 block, same pattern) check bound compliance the same way every competitor arm does:
```r
ok_lb <- isTRUE(res_lb$status %in% c(0L, 5L)) &&
  max(w_lb_n) <= max_weight + 1e-6
```
2. Carry the known caveat into the row's `note` field so it survives into the CSV and can be transcribed into `docs/performance.md`:
```r
note_lb <- sprintf(
  "convergence=list() (per-method natural default per R/harvest.R:424); status=%d, iterations=%d%s",
  res_lb$status, res_lb$iterations,
  if (method_name == "newton_kl" && max(w_lb_n) > 3 + 1e-6)
    "; max_w exceeds max_weight=3 -- newton_kl reports (RK_ERR_NOCONV) rather than clamping bound violations, per leafblower-73d7 (documented solver contract, not a benchmark-script bug)"
  else "")
```
and add the equivalent sentence to `docs/performance.md`'s "Notes carried verbatim" paragraph (`docs/performance.md:91-94`).

---

### WR-02: `ok` computed on unnormalized weights while the reported `max_w`/`min_w` use the normalized scale (scale mismatch)

**File:** `benchmarks/oris_soft_vs_competitors.R:173` (pre-existing `run_input_class` survey_calibrate arm) and `:429` (new in plan 03-05, `survey_calfun_arm_row()`)

**Issue:** Both `survey_calibrate` competitor blocks compute the bound-compliance flag on the *raw* returned weight vector, but report `max_w`/`min_w` on the *normalized* (`normalize_to_n()`) vector used everywhere else in the file for cross-arm comparability:

```r
w_sv_n <- normalize_to_n(as.numeric(w_sv), n)
max_error_sv <- margin_max_error(w_sv_n, df, tgt)
ok_sv <- all(is.finite(w_sv)) && max(w_sv) <= 3 + 1e-6   # <- raw w_sv, not w_sv_n
...
row <- arm_row(..., max(w_sv_n), min(w_sv_n), ...)         # <- w_sv_n reported
```
(identical pattern at line 173 for the pre-existing `survey_calibrate` arm, and freshly copy-pasted into the new `survey_calfun_arm_row()` at line 429, which is what the greg/logit competitor rows in this plan's own new code use). Every other arm in the file (`icarus_calibration`, `ReGenesees_e_calibrate`, `optweight_linf`) computes `ok` on the *same* normalized vector it reports — this pair is the outlier. Because `survey::calibrate`'s calibrated weights already sum very close to `n` (population totals were themselves scaled to `n`), the two scales happen to coincide numerically in this run, so it produced no visible discrepancy this time — but the check and the reported metric are testing two different numbers, and a future run where `survey::calibrate` converges to a weight sum further from `n` (e.g. under `epsilon`-tolerance slack) could produce an `ok` flag that disagrees with the `max_w` column shown in the same row.

**Fix:** Use the normalized vector consistently, matching the icarus/ReGenesees pattern:
```r
ok_sv <- all(is.finite(w_sv_n)) && max(w_sv_n) <= max_weight + 1e-6
```
(and the `3 + 1e-6` variant in `survey_calfun_arm_row`).

---

### WR-03: stale comment claims a `1e-10` tolerance "this script itself asserts" — the script uses `1e-6` everywhere

**File:** `benchmarks/oris_soft_vs_competitors.R:58-69` (comment on `normalize_to_n()`)

**Issue:** The comment justifying why `normalize_to_n()` is applied only to competitor arms says: *"re-scaling leafblower's own output pushed a weight sitting exactly at the max_weight bound to 3.0000000036, failing the max_w <= max_weight + 1e-10 check this script itself asserts."* No `1e-10` tolerance exists anywhere in this script — every `ok_*` bound check in this file (lines 173, 228, 293, 429, and the hardcoded `+ 1e-6` on the chebyshev/optweight block) uses `1e-6`. The `1e-10` tolerance the comment describes is `tests/testthat/test-bench-gate.R:45` (`expect_lte(r$max_w, r$max_weight + 1e-10, ...)`) — a different file, with a different (stricter) tolerance, that this script does not itself enforce or reference. The comment's underlying point (do not renormalize leafblower's own already-`sum(w)==n` output) is correct, but attributing the `1e-10` check to "this script itself" is inaccurate and could mislead a future reader about what tolerance this file's own `ok` columns actually use.

**Fix:** Correct the comment to name the actual source of the `1e-10` check:
```r
# ... failing the max_w <= max_weight + 1e-10 check
# tests/testthat/test-bench-gate.R asserts (this script's own ok_* checks
# use a looser 1e-6 tolerance, see below).
```

## Info

### IN-01: cross-language wall-time comparison uses two different timer APIs without a stated caveat

**File:** `benchmarks/greenkhorn_sinkhorn_vs_pot.py:153-161` (`timed_median()`, `time.perf_counter()`) vs. `benchmarks/oris_soft_vs_competitors.R` (`bench::mark()` throughout)

**Issue:** The `k2_margin_pot_equiv` table directly juxtaposes `leafblower_greenkhorn`/`leafblower_sinkhorn` wall times (measured in-process via Python's `time.perf_counter()`, harness = `harvest()` called through the Python bindings) against `pot_greenkhorn`/`pot_sinkhorn` (same timer, same process — internally consistent) — that part is fine. But the same table is presented alongside the R-measured `oris_soft_vs_competitors.csv` figures elsewhere on the page, and neither `docs/performance.md` nor the script states that R's `bench::mark()` (which does GC-aware bookkeeping around each call) and Python's bare `perf_counter()` loop are not directly interchangeable measurement methodologies. This is a minor methodological caveat, not a data-correctness issue — the two tables are never compared to each other numerically on the page — but is worth a one-line note given the phase's stated bar for "every wall-time figure ... never alone."

**Fix:** Optional one-sentence caveat in `docs/performance.md`'s Methodology section noting the Python-side timer differs from R's `bench::mark()`, or leave as-is since the two measurement classes (`medium_100k_5margins` vs `k2_margin_pot_equiv`) are never compared against each other in the same row.

### IN-02: `docs/performance.md`'s Competitors table row for `oris`/`oris_soft`/`raking` repeats "Not run at `large_stepstone_fulldata` scale" identically across three rows via "Same as `oris`"

**File:** `docs/performance.md:271-274`

**Issue:** Minor duplication style note, not a defect: three consecutive rows (`oris`, `oris_soft`, `raking`) all resolve to the same competitor set and the same caveat, with two of the three rows saying "Same as `oris`." This is intentional and legible (per-method table honesty > DRY, consistent with the project's own "pragmatic DRY" convention), not something to change — noted only because a first read looks like copy-paste until the second row's "Same as `oris`" is read.

**Fix:** None needed — flagged for completeness, not actionable.

---

_Reviewed: 2026-08-16T02:13:58+02:00_
_Reviewer: Claude (gsd-code-reviewer)_
_Depth: standard_
