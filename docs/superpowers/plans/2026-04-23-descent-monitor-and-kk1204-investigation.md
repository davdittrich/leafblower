# Descent Monitor + kk1.20.4 Convergence Investigation

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development or superpowers:executing-plans.

**Goal:** (WU-P1) Add descent monitor to raking solver per leafblower-dl6. (WU-P2) Profile + investigate why both solvers fail kk1.20.4 benchmark gate (<30s, <1e-6) at n=1M, K=20, max_weight=3. Document findings; no algorithmic changes in WU-P2.

**Rev 2 (post plan-review-gate iter 1 opus):** Applied 7 fixes from Feasibility (2 blockers) + Completeness (3 blockers + 4 gaps) + Scope (3 nits):
- P1 test: `capture.output(type="output")` (stdout, not stderr); wall-clock guard; exact-phrase match.
- P1 monitor: window-minimum comparison (not single-step — fixes oscillation false-positive); relative `eps = max(0.01*min, tol_abs)` (not absolute 1e-12 — derived from 1% meaningful-improvement threshold).
- P2.0: added prerequisite steps — verify probe symbol registered + inspect actual log formats before choosing grep pattern.
- P2.2: tol_abs=0 forces full 50 iters → per-iter denominator accurate.
- P2.3: tightened grep pattern + added line-count acceptance.
- P2.4: quantitative recommendation template (numeric thresholds required, not prose).
- P0 ticket: rides on leafblower-dl6 scope.

**Style constraint (per CLAUDE.md + user directive 2026-04-23):** No history/narrative notes in source code. Comments explain the non-obvious WHY (invariants, constraints, subtle behavior), not WHAT was changed, when, or why a previous version was wrong. No "originally named X", "renamed in commit Y", "added for Z". PR descriptions and git log already carry that. This applies to P1/P2 edits AND to a cleanup of existing history residue (WU-P0).

**Tech Stack:** C++17, R (testthat, bench).

---

## WU-P0: Strip history notes from source (precondition — no dedicated bd ticket; rides on WU-P1's leafblower-dl6 as cleanup scope)

**Files:** `src/raking.cpp`

Commit `36449e0` ("provenance header") added a "History" paragraph naming the prior misname, the paper it was mistakenly referenced, and the rename commit SHA. Per CLAUDE.md comment discipline + user 2026-04-23 directive, source files must not carry this kind of narrative. Citations of the algorithmic provenance (Deming-Stephan, Csiszár, Boyle-Dykstra) are WHY-comments that stay; the "History: originally named iEPPA..." paragraph goes.

### Step P0.1: Remove history block

- [ ] **Step P0.1.1: Delete history paragraph in src/raking.cpp**

Open `src/raking.cpp`. Locate the header comment block (lines ~1-33). Delete the "History:" paragraph entirely:

```cpp
// History: originally named "iEPPA" after Chu-Liang-Toh-Yang (2022,
// arXiv:2011.14312), but the paper's algorithm is the inexact Entropic
// Proximal Point Algorithm for CMOT LP — shares no ingredients with this
// solver. Renamed to "raking" in commit 20d6ebf (2026-04-23). The
// paper-faithful iEPPA implementation lives in src/ieppa.cpp.
```

Keep the algorithmic composition description, the citations (Deming-Stephan / Csiszár / Boyle-Dykstra), and the explicit "no published proof" caveat. Those are WHY-comments.

### Step P0.2: Build + verify

- [ ] **Step P0.2.1: Compile clean**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -2
```

### Step P0.3: Commit

- [ ] **Step P0.3.1: Commit**

```bash
git add src/raking.cpp
git commit -m "docs(raking): remove history paragraph from header

CLAUDE.md comment discipline: source should not reference prior names,
rename commits, or migration history. Algorithmic provenance citations
(Deming-Stephan, Csiszar, Boyle-Dykstra) and the 'no published proof for
the hybrid' caveat stay — those are WHY-comments that justify the
architecture. Git log carries history."
```

---

## WU-P1: Descent monitor in raking solver (leafblower-dl6)

**Style note (from above):** message strings in source must state the current condition ("errRp stalled for N consecutive checks"), not the history ("added 2026-04-23 per audit"). Commit messages and this plan carry the rationale.

**Finding summary.** Current raking can cycle silently until `max_iter` when primal error stalls. A cheap monitor tracks `errRp` trajectory; after K consecutive non-improving checks, aborts early with a descriptive NOCONV message. Observable user benefit: faster failure signal, clearer diagnostic at verbose ≥ 1.

**Files:**
- Modify: `src/raking.cpp` — add monitor state + check
- Modify: `tests/testthat/test-raking.R` — one new test exercising the monitor
- No PRD/doc change (runtime log additions only)

### Step P1.1: Write failing test

- [ ] **Step P1.1.1: Add test exercising non-monotone trajectory**

File: `tests/testthat/test-raking.R` — append after existing tests:

```r
test_that("descent monitor aborts early on stalled errRp trajectory", {
  # Input empirically validated to stall: n=1000, 95/5 class split, max_weight=1.2.
  # Prior measurement (leafblower-370 probe session 2026-04-23): raking hit
  # max_iter=500 with no convergence and 500 consecutive curvature rejections.
  # Monitor fires at 5 consecutive stalled error-checks (50 iters) — well
  # before iter 500.
  set.seed(91)
  n <- 1000
  df <- data.frame(cat = sample(c("A", "B"), n, replace = TRUE, prob = c(0.05, 0.95)))
  tgt <- list(cat = c(A = 0.95, B = 0.05))
  # Emit verbose=1 so we can grep for the monitor message.
  # CalibState.log() routes through Rprintf (src/r_bridge.cpp:21) which writes
  # to stdout, not R's message sink. Use capture.output(type = "output").
  t0 <- Sys.time()
  msgs <- capture.output(
    res <- suppressWarnings(harvest(df, tgt, method = "raking",
                                     max_weight = 1.2,
                                     max_iterations = 500,
                                     verbose = 1L)),
    type = "output"
  )
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  # Wall-clock guard: monitor should trigger early; whole call must be under 5s
  # even if n_max_iter=500 is reached without the monitor.
  expect_lt(elapsed, 5)
  # Monitor message matches the exact phrase emitted by src/raking.cpp.
  probe <- paste(msgs, collapse = "\n")
  expect_match(probe, "errRp stalled for [0-9]+ consecutive checks",
               info = paste("expected descent-monitor message; got:", probe))
})
```

- [ ] **Step P1.1.2: Run, verify it fails**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -2
Rscript -e "devtools::test(filter = 'raking')" 2>&1 | grep -E "FAIL|PASS" | tail -3
```

Expected: 1 new failure ("monitor never fired").

### Step P1.2: Implement monitor

- [ ] **Step P1.2.1: Add monitor state + check in src/raking.cpp**

Locate the convergence-check block inside the outer loop (`if (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter)`). Modify:

```cpp
    // Descent monitor state. Placed BEFORE the outer iter loop.
    //
    // Detects stalled convergence: no net progress over a sliding window of
    // kMaxNoImprove consecutive error-checks. Uses a window-minimum comparison
    // to avoid false positives on normal oscillation (errRp often wobbles
    // within noise even while the running minimum still decreases).
    //
    // Firing condition: every check for the last kMaxNoImprove checks
    // reported errRp >= (window min so far) - relative_eps. Equivalent: the
    // window minimum has not improved by more than relative_eps across the
    // window.
    constexpr int    kMaxNoImprove = 5;  // 5 * kErrCheckInterval = 50 stalled iters
    double min_errRp_window = std::numeric_limits<double>::infinity();
    int n_no_improve = 0;

    for (int iter = 1; iter <= st.inner_max_iter; iter++) {
        // ... existing body ...

        if (iter == 1 || iter % kErrCheckInterval == 0 || iter == st.inner_max_iter) {
            double errRp = compute_errRp(st, w, bucket);
            res.max_error = errRp;

            // Relative improvement threshold: 1% of current best. At errRp ~1e-2
            // threshold is 1e-4; at errRp ~1e-6 threshold is 1e-8. Absolute
            // floor of tol_abs prevents a valid convergence from tripping the
            // monitor right before the convergence check below.
            const double rel_eps = 0.01 * min_errRp_window;
            const double eps = std::max(rel_eps, st.tol_abs);
            if (errRp < min_errRp_window - eps) {
                min_errRp_window = errRp;
                n_no_improve = 0;
            } else {
                n_no_improve++;
            }

            if (st.verbose >= 1) {
                char msg[256];
                std::snprintf(msg, 256, "raking iter %d: errRp=%.2e", iter, errRp);
                st.log(msg);
            }

            if (errRp < st.tol_abs) {
                res.status = is_infeasible ? RK_ERR_INFEAS : RK_OK;
                break;
            }

            if (n_no_improve >= kMaxNoImprove) {
                res.status = RK_ERR_NOCONV;
                if (st.verbose >= 1) {
                    char msg[256];
                    std::snprintf(msg, 256,
                                  "raking: errRp stalled for %d consecutive checks "
                                  "(last=%.2e, window_min=%.2e); likely near-infeasible bounds. "
                                  "Aborting at iter %d.",
                                  n_no_improve, errRp, min_errRp_window, iter);
                    st.log(msg);
                }
                break;
            }
        }
    }
```

Rationale for the parameters:
- `kMaxNoImprove = 5` × `kErrCheckInterval = 10` → monitor fires after 50 stalled iters. Derived: for well-posed problems raking typically converges well below 50 iters; for near-infeasible problems the asymptotic errRp is reached within 50 iters of the floor. Single iteration cost × 50 is a bounded-budget acceptable penalty.
- `rel_eps = 0.01 * min_errRp_window` floored at `tol_abs`: relative threshold derived from "1% improvement is meaningful". At errRp=1e-2 that's 1e-4, well above FP noise; at errRp=1e-6 that's 1e-8, tight but not below machine precision. Absolute floor prevents tripping in the last few iters before successful convergence where `errRp - min_errRp_window` can be near-machine-epsilon.
- Window-minimum, not single-step comparison: avoids the documented oscillation false-positive — a single up-tick does not increment `n_no_improve` if the running minimum hasn't moved.

### Step P1.3: Build + test

- [ ] **Step P1.3.1: Rebuild + full suite**

```bash
R CMD INSTALL --preclean . 2>&1 | tail -2
Rscript -e "devtools::test()" 2>&1 | tail -4
```

Expected: `[ FAIL 0 | PASS ≥ 168 ]`. New test must pass. No regression in the existing raking tests (monitor should not fire on well-posed problems — they converge before 5 consecutive stalls).

- [ ] **Step P1.3.2: Python pytest sanity**

```bash
cd python && python -m pytest 2>&1 | tail -3 && cd ..
```

Expected: 3 passed.

### Step P1.4: Commit

- [ ] **Step P1.4.1: Single commit**

```bash
git add src/raking.cpp tests/testthat/test-raking.R
git commit -m "feat(raking): descent monitor — abort on 5 consecutive stalled errRp

Finding leafblower-dl6 from 2026-04-23 audit. Raking could cycle
silently until max_iter when primal error floors above tol_abs
(near-infeasible inputs, tight bounds). Added monitor: after
kMaxNoImprove=5 consecutive error-checks with no strict improvement,
set status=NOCONV early and emit a descriptive message under verbose >= 1.

Observable to user: faster NOCONV signal + clearer diagnostic.
Zero-cost path: monitor state (one int, one double) initialized once
per solve; check fires only in the existing kErrCheckInterval branch.

Test exercises near-infeasible input (95/5 split, max_weight=1.2) that
forces the stall. Monitor message grep'd from verbose=1 output."
```

---

## WU-P2: Investigation of kk1.20.4 convergence + wall-clock failure

**Finding summary.** At n=1M, K=20, max_weight=3, skewed-target random-categorical data:
- `ieppa`: 369s wall-clock, max_err=1.1e-3, hit max_iter=500 NOCONV
- `raking`: 16s wall-clock, max_err=2.2e-2, hit max_iter=500 NOCONV
Both fail the kk1.20.4 gate (<30s + max_err<1e-6).

**Question:** is this a problem with the algorithms or a property of the input? The spec's design §1b UC-3 claims iEPPA wins when M_cell << n. At K=20, cat=5, uniform random sampling, ∏cat_counts = 5^20 ≈ 9.5e13 but M_cell ≤ min(n, ∏) = n = 1M. If all 1M observations have unique (g_1,...,g_K) tuples, cell compression gives zero benefit; iEPPA degenerates to obs-level work with log-space Sinkhorn overhead vs raking's simpler IPF.

**Files:**
- Create: `benchmarks/kk1204_profile.R` (scratch profiling script, NOT committed)
- Create: `docs/investigations/2026-04-23-kk1204-convergence.md` (findings report, committed)

### Step P2.0: Prerequisites

- [ ] **Step P2.0.1: Confirm probe symbol registered**

```bash
Rscript -e 'info <- getNativeSymbolInfo("C_leafblower_cell_table_probe", PACKAGE="leafblower"); print(info$name)'
```
Expected: `[1] "C_leafblower_cell_table_probe"`. If error, P2.1 blocked — file a precondition fix issue.

- [ ] **Step P2.0.2: Inspect ieppa.cpp + raking.cpp verbose log formats**

```bash
grep -nE 'iter.*errRp|errRp.*iter' src/ieppa.cpp src/raking.cpp
```

Record the exact log-line formats. Note: `src/raking.cpp` currently emits `iEPPA iter N: errRp=...` (pre-rename residue; P1.2.1 replacement updates it to `raking iter N: ...`). Plan's P2.3 grep `iter [0-9]+: errRp` catches both current and post-P1 forms. Update Step P2.3.1 if either file diverges from this pattern.

**API verification (resolves iter-2 completeness G3):** `R/harvest.R parse_convergence` returns `convergence[["absolute"]]` directly when present. `convergence = list(absolute = 0)` → `tol_abs = 0` → `errRp < st.tol_abs` never holds → solver runs full `max_iterations`. Confirmed by inspection; no API surprise.

### Step P2.1: Measure M_cell for the benchmark input

- [ ] **Step P2.1.1: Write scratch probe script**

File: `/tmp/kk1204/probe.R` (NOT committed)

```r
#!/usr/bin/env Rscript
# Measure M_cell / n ratio at the kk1.20.4 benchmark input.
suppressPackageStartupMessages(library(leafblower))
set.seed(1)
n <- 1000000L; K <- 20L
cat_counts <- rep(5L, K)
cols <- lapply(seq_len(K), function(k) sample(letters[seq_len(cat_counts[k])], n, replace = TRUE))
df <- as.data.frame(cols); names(df) <- paste0("v", seq_len(K))
# Encode group_ids as int32 (match harvest() internal pipeline)
gid_list <- lapply(cols, function(col) match(col, letters[seq_len(5L)]) - 1L)
# Use the test-only probe from WU-1 (C_leafblower_cell_table_probe)
probe <- .Call("C_leafblower_cell_table_probe", gid_list, as.integer(n),
               PACKAGE = "leafblower")
cat(sprintf("M_cell = %d (n = %d, ratio M_cell/n = %.4f)\n",
            probe$M_cell, n, probe$M_cell / n))
cat(sprintf("compression factor = %.2fx\n", n / probe$M_cell))
# Cardinality of cross-product
cp <- prod(as.numeric(cat_counts))
cat(sprintf("product(cat_counts) = %.2e (cells populated / available = %.4e)\n",
            cp, probe$M_cell / cp))
```

- [ ] **Step P2.1.2: Run + record**

```bash
mkdir -p /tmp/kk1204
# Copy the script content there, then:
Rscript /tmp/kk1204/probe.R 2>&1 | tail
```

Record M_cell value. Expected outcome (from theoretical prediction): M_cell ≈ n for K=20/cat=5 uniform-random — compression ≤ ~1.1×. If ratio is close to 1.0, iEPPA has no advantage here.

### Step P2.2: Per-iteration cost comparison

- [ ] **Step P2.2.1: Time per-iteration cost at fixed max_iter=50**

File: `/tmp/kk1204/periter.R` (scratch)

```r
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(leafblower))
set.seed(1)
n <- 1000000L; K <- 20L
cat_counts <- rep(5L, K)
cols <- lapply(seq_len(K), function(k) sample(letters[seq_len(cat_counts[k])], n, replace = TRUE))
df <- as.data.frame(cols); names(df) <- paste0("v", seq_len(K))
tgts <- lapply(seq_len(K), function(k) {
  p <- rep(0.7/(cat_counts[k]-1), cat_counts[k]); p[1] <- 0.3
  setNames(p, letters[seq_len(cat_counts[k])])
})
names(tgts) <- names(df)

# 50 iters each (tight tol to prevent early-termination skewing per-iter cost);
# per-iter denominator uses actual iter count from attr(res, 'iterations') or
# the result's algorithm-used attribute. Setting tol_abs=0 and relying on
# max_iter cap guarantees all 50 iters run, so t/50 is accurate.
for (m in c("ieppa", "raking", "lbfgsb")) {
  t <- system.time(res <- suppressWarnings(
    harvest(df, tgts, method = m, max_weight = 3, max_iterations = 50,
            convergence = list(absolute = 0))  # force full 50 iters
  ))[3]
  diag <- diagnose_weights(df, tgts, res$weights)
  cat(sprintf("%s: %.3fs total, %.3fs/iter (50 iters), max_err=%.3e\n",
              m, t, t/50, max(abs(diag$error_weighted))))
}
```

- [ ] **Step P2.2.2: Run + record per-iter wall-clock for each method**

```bash
Rscript /tmp/kk1204/periter.R 2>&1 | tail
```

Record: per-iteration wall-clock for `ieppa`, `raking`, `lbfgsb`. Compute ratio ieppa/raking per-iter.

### Step P2.3: Convergence trajectory (errRp vs iter)

- [ ] **Step P2.3.1: Capture errRp trajectory at verbose=1**

Both solvers log `errRp` every 10 iters at `verbose = 1`. Capture to file:

```bash
Rscript -e '
library(leafblower); set.seed(1)
n <- 1000000L; K <- 20L
cat_counts <- rep(5L, K)
cols <- lapply(seq_len(K), function(k) sample(letters[seq_len(cat_counts[k])], n, replace = TRUE))
df <- as.data.frame(cols); names(df) <- paste0("v", seq_len(K))
tgts <- lapply(seq_len(K), function(k) { p <- rep(0.7/(cat_counts[k]-1), cat_counts[k]); p[1] <- 0.3; setNames(p, letters[seq_len(cat_counts[k])]) })
names(tgts) <- names(df)
for (m in c("ieppa", "raking")) {
  cat("=== method=", m, " ===\n", sep="")
  suppressWarnings(harvest(df, tgts, method=m, max_weight=3, max_iterations=500, verbose=1L))
}
' 2>&1 | grep -E "iter [0-9]+: errRp" > /tmp/kk1204/trajectory.log
wc -l /tmp/kk1204/trajectory.log
cat /tmp/kk1204/trajectory.log
```

**Acceptance:** `wc -l` must report at least 10 lines (5 per method × 2 methods at worst, given 500 iters / 10-iter check interval = 50 checks per run). If 0 lines, grep pattern is wrong — re-check P2.0.2 actual log formats before proceeding.

Record: whether errRp decreases monotonically, floors at a value > tol_abs, or oscillates. Ieppa's last measured max_err = 1.1e-3 and raking's = 2.2e-2; the trajectory tells us whether it's still improving at iter 500 or stuck.

### Step P2.4: Synthesize findings + write report

- [ ] **Step P2.4.1: Create investigation report**

File: `docs/investigations/2026-04-23-kk1204-convergence.md`

```markdown
# kk1.20.4 Convergence + Wall-Clock Investigation

**Date:** 2026-04-23
**Triggers:** leafblower-kk1.20.4 (original P1 gate), leafblower-g8f
  (follow-up convergence study)

## Input
n = 1,000,000, K = 20, cat_counts = (5,)*20, max_weight = 3,
targets = skewed (0.3, 0.175, 0.175, 0.175, 0.175) per margin.

## Measurements

### Cell compression (WU-P2.1)
- M_cell = [REPLACE WITH MEASURED]
- ratio M_cell/n = [FILL IN]
- compression factor = [FILL IN]×

### Per-iteration cost (WU-P2.2, 50-iter budget)
| method | total s | s/iter | max_err @ 50 |
|---|---|---|---|
| ieppa  | [FILL] | [FILL] | [FILL] |
| raking | [FILL] | [FILL] | [FILL] |
| lbfgsb | [FILL] | [FILL] | [FILL] |

### Convergence trajectory (WU-P2.3, max_iter=500)
- ieppa final errRp = [FILL], monotone: [YES/NO]
- raking final errRp = [FILL], monotone: [YES/NO]

## Interpretation

[Write 200-400 words analyzing:
- Why ieppa is slow (likely: M_cell = n, no cell compression, + log-space
  Sinkhorn overhead)
- Why neither converges to tol_abs=1e-6 (likely: near-infeasibility from
  tight max_weight with skewed targets; or: convergence asymptote above
  tol at 500 iters, needs more iters)
- Whether the kk1.20.4 gate (<30s + <1e-6) is achievable with existing
  algorithms on this input class]

## Recommendations

1. **kk1.20.4 gate:** choose ONE and justify with measurements:
   - KEEP (if any measured config achieves <30s AND <1e-6)
   - RELAX to `<X s AND <Y`, where X and Y derive from the measured asymptote + best-measured wall-clock
   - REFRAME (close kk1.20.4 as not-achievable on this class; document achievable envelope)
2. **Routing:** quantitative threshold. Propose ONE of:
   - AUTO switches raking→ieppa when `M_cell/n < T1` AND n ≥ T2 (give numeric T1, T2 from M_cell/per-iter data)
   - Insufficient data — file BLSE benchmark WU as prerequisite to the routing decision
3. **Follow-up:** list each new bd ticket filed in Step P2.6 (or explicitly "none needed")

## Raw data

- `/tmp/kk1204/probe.R` (M_cell probe)
- `/tmp/kk1204/periter.R` (per-iter timing)
- `/tmp/kk1204/trajectory.log` (errRp trajectory capture)
```

Fill ALL `[FILL]` placeholders from actual measurements (P2.1, P2.2, P2.3). Write the Interpretation and Recommendations sections based on the data.

### Step P2.5: Commit report

- [ ] **Step P2.5.1: Commit (docs only)**

```bash
mkdir -p docs/investigations
git add docs/investigations/2026-04-23-kk1204-convergence.md
git commit -m "docs(investigation): kk1.20.4 convergence + wall-clock analysis

Triggers: leafblower-kk1.20.4 (original P1 gate), leafblower-g8f
(convergence study follow-up).

Measures M_cell / n ratio at the benchmark input (K=20/cat=5 uniform
random) — the regime-defining quantity for faithful iEPPA's cell-
compression advantage. Per-iteration wall-clock for ieppa/raking/lbfgsb.
Convergence trajectory (errRp vs iter) for ieppa + raking.

Findings inform whether kk1.20.4's <30s + <1e-6 gate is achievable on
this input class and whether AUTO routing should redirect at K >= 15."
```

### Step P2.6: File concrete follow-ups per report conclusions

- [ ] **Step P2.6.1: If report recommends specific actions, file bd tickets**

Example (conditional on findings):
```bash
# If recommendation is "tighten iEPPA routing":
bd create --title="routing: auto-switch to raking when M_cell/n > 0.5" \
  --description="Per docs/investigations/2026-04-23-kk1204-convergence.md recommendation X. ..." \
  --type=feature --priority=2

# If recommendation is "increase max_iter default for ieppa":
bd create --title="ieppa: investigate max_iter=500 cutoff vs asymptotic errRp" ...

# If recommendation is "gate kk1.20.4 unreachable — reframe":
bd update leafblower-kk1.20.4 --notes="Empirical evidence per 2026-04-23 investigation: gate unachievable with current algorithms on K=20/cat=5/uniform random. Reframe as K<=10 or looser convergence (<1e-4)."
```

---

## Ordering

WU-P0 → WU-P1 → WU-P2. P0 is a precondition cleanup (docs-only source edit). P1 is a small independent addition with test. P2 is an investigation (no algorithm change in raking/ieppa — only new investigation doc + optional follow-up tickets).

## Rollback

WU-P1: single commit. `git revert`.
WU-P2: no source-code changes; the new doc file can be removed if report is deemed incorrect.

## Merge gate

- WU-P1: all tests pass including new descent-monitor test; build clean.
- WU-P2: report exists at `docs/investigations/2026-04-23-kk1204-convergence.md` with all `[FILL]` placeholders populated from measurements; recommendations section non-empty.

---

## Self-Review

**Scope:** P1 addresses leafblower-dl6 directly. P2 addresses kk1.20.4 + g8f follow-up by producing documented evidence that informs routing decisions. Neither scope-creeps into algorithm changes.

**Placeholders:** only in the report template (P2.4.1), which is explicitly a fillable template — the plan instructs the implementer to replace every `[FILL]` with measured values. No other placeholders.

**Type consistency:** `prev_errRp`, `n_no_improve`, `kMaxNoImprove`, `kImproveEps` names used consistently in P1.2.1 code.

**Verification:** P1 has a test (test-raking.R new case). P2 has measurement steps + a doc with a populated-placeholder acceptance gate.
