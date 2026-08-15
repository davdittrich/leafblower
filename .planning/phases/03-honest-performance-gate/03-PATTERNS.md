# Phase 3: Honest Performance Gate - Pattern Map

**Mapped:** 2026-08-15
**Files analyzed:** 6
**Analogs found:** 6 / 6

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|----------------|---------------|
| `benchmarks/oris_soft_vs_competitors.R` | benchmark script | batch (cross-solver timing loop) | `benchmarks/newton_kl_bench.R` | exact (same repo, same `bench::mark()` cross-method loop shape) |
| `benchmarks/run_honest_gate.sh` | wrapper/config | batch (shell orchestration) | `benchmarks/run_allmethod.sh` | exact (same 2-step `set -euo pipefail` + `Rscript` wrapper shape) |
| `tests/testthat/test-bench-gate.R` (extend, or adjacent new test) | test (opt-in perf gate) | request-response (assert on solver output) | `tests/testthat/test-bench-gate.R` itself (existing `LBW_BENCH_GATE=1` test in the same file) | exact — modify in place, reuse existing skip idiom |
| `README.md` (new) | doc/config | — | none in-repo (confirmed never existed via `git log --all --full-history -- README.md`) | no analog — use RESEARCH.md D-13 shape directly |
| `docs/performance.md` (new, exact name Claude's discretion) | doc | — | `docs/methods/oris.md` (competitor table + citation style) | role-match (doc structure/citation convention) |
| `.planning/REQUIREMENTS.md` KPI-04 row (edit) | doc/config | — | same file, US-003 row (lines 100-112, already rewritten this session by researcher) | exact — same file, same table conventions |

## Pattern Assignments

### `benchmarks/oris_soft_vs_competitors.R` (benchmark script, batch)

**Analog:** `benchmarks/newton_kl_bench.R` (full file read, 58 lines)

**Header/env pattern** (lines 1-15):
```r
Sys.setenv(OMP_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(leafblower)
  library(bench)
})
```
Note: CLAUDE.md's DoD command additionally sets `OPENBLAS_NUM_THREADS`/`MKL_NUM_THREADS`
at the shell level (`benchmarks/stepstone_all_methods.R:3` comment: "set by caller") — the
new script should rely on the shell wrapper (`run_honest_gate.sh`) exporting all three, not
just `OMP_NUM_THREADS` in-script, to match the stricter CLAUDE.md envelope.

**Core cross-solver timing loop** (lines 17-48, this is the pattern to copy verbatim in shape):
```r
run_one <- function(skew_name, p_skew, n = 1e6L, K = 20L, nj = 5L,
                    methods = c("newton_kl", "oris")) {
  set.seed(1)
  df <- as.data.frame(lapply(seq_len(K), function(k)
    factor(sample(letters[seq_len(nj)], n, TRUE))))
  names(df) <- paste0("m", seq_len(K))
  ...
  for (m in methods) {
    res <- bench::mark(
      run = harvest(df, tgt, method = m, max_weight = 3,
                    max_iterations = 50, accelerate = (m == "oris")),
      iterations = 2, check = FALSE, memory = FALSE, filter_gc = FALSE
    )
    r <- harvest(df, tgt, method = m, max_weight = 3,
                 max_iterations = 50, accelerate = (m == "oris"))
    R <- attr(r, "result")
    cat(sprintf("  %-12s wall=%6.2fs status=%d max_err=%.3e iters=%d\n",
                m, as.numeric(res$median), R$status, R$max_error, R$iterations))
    out[[m]] <- data.frame(fixture = skew_name, method = m,
      wall_s = as.numeric(res$median), status = R$status,
      max_error = R$max_error, iterations = R$iterations, n = n, K = K, nj = nj)
  }
  do.call(rbind, out)
}
```
For the new script, `methods` becomes the doc-named competitors (`oris_soft`, plus one call
per competitor package: `survey::calibrate()`, `icarus::calibration()`,
`ReGenesees::e.calibrate()` — these are NOT `leafblower::harvest(method=...)` calls, they are
separate package calls, so the loop body needs a `switch`/dispatch per package rather than a
uniform `harvest(method=m)` call). Per RESEARCH.md Pitfall 3, do NOT reuse this file's own
K=20-uniform-random fixture generator — use a `stepstone_benchmark.R`-style realistic fixture
instead (see below).

**Output/write pattern** (lines 50-57):
```r
dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
write.csv(results, "benchmarks/results/newton_kl_kk1204.csv", row.names = FALSE)
cat("\nWrote benchmarks/results/newton_kl_kk1204.csv\n")
```
Copy this shape, writing to a new `benchmarks/results/oris_soft_vs_competitors.csv`.

**Fixture-generation analog (for a realistic, non-degenerate input class per Pitfall 3):**
`benchmarks/stepstone_benchmark.R` (read lines 1-229) — reuse its `arrow::write_parquet`,
`jsonlite::write_json` fixture I/O pattern (lines 158-162) and its `time_one()`/`report()`
helpers (lines 171-185) if the new script needs a purpose-built large fixture rather than the
existing `benchmarks/stepstone_bench_data.parquet`. Its competitor-timing loop compares
against `autumn` — **D-05 forbids citing autumn**; do not copy the `autumn::harvest` calls
(lines 187-201) or the `speedup <- r_autumn$ms / r_lb$ms` framing (line 223) — only the
env/fixture/report scaffolding, not the comparison target.

**Package guard pattern (new — no exact analog, D-09 requirement):** the script must guard
`icarus`/`ReGenesees` calls with `requireNamespace()` + informative skip, since RESEARCH.md's
recommendation is NOT to add them to package `DESCRIPTION` `Suggests:`. No existing file in
this repo does this guard pattern for a benchmark competitor (all existing benchmarks assume
their deps installed) — this is genuinely new code, not a copy.

---

### `benchmarks/run_honest_gate.sh` (wrapper, batch)

**Analog:** `benchmarks/run_allmethod.sh` (full file, 6 lines)
```bash
#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
echo '=== R benchmark ===' && Rscript benchmarks/allmethod_bench.R
echo '=== Python benchmark ===' && python benchmarks/allmethod_bench.py
```
Copy this shape exactly, substituting the two steps per D-12 (SC4's "one command"):
```bash
#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
echo '=== Stepstone regression gate ===' && \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  LBW_BENCH_GATE=1 Rscript -e "testthat::test_dir('tests/testthat', filter='bench-gate')"
echo '=== oris_soft vs competitors ===' && \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  Rscript benchmarks/oris_soft_vs_competitors.R
```
This is the single-thread-BLAS envelope (CLAUDE.md build/test section, RESEARCH.md Pattern 2)
applied at the shell level, matching `run_allmethod.sh`'s two-step `echo && Rscript` idiom.

---

### `tests/testthat/test-bench-gate.R` (test, opt-in perf gate — extend in place)

**Analog:** the file's own existing two tests (full file, 54 lines read).

**Existing opt-in gate pattern to copy for the new oris_soft headline assertion** (lines 1-17):
```r
test_that("stepstone-fulldata AB config meets merge floor (errRp) + Pearson agreement", {
  skip_on_cran()
  skip_if(Sys.getenv("LBW_BENCH_GATE") == "")
  rpt_path <- "benchmarks/stepstone_fulldata_homotopy_report.rds"
  skip_if(!file.exists(rpt_path))
  ...
})
```

**Existing kk1204 test to reconcile, NOT delete blindly** (lines 28-54):
```r
test_that("kk1204 gate: n=500k K=20 converges in <30s with best_error<1e-3", {
  skip_on_cran()
  skip_if(Sys.getenv("CI") != "")
  ...
  expect_lte(elapsed, 30, label = "elapsed: speed gate <30s on n=500k K=20")
  expect_lte(r$best_error, 1e-3, label = "best_error: quality gate below 1e-3")
})
```
Per D-01/D-14 and RESEARCH.md Pitfall 2, this test's fate (keep as internal regression floor
with corrected label / relabel / retire) is a planner decision, not a pattern-mapping one —
flagged here so the plan references this exact code, not a re-derived version.

**New assertion to add** (mirrors the `LBW_BENCH_GATE=1` skip idiom above, NOT the `CI`-skip
idiom): asserts the `oris_soft` headline number (wall time and/or DEFF/ESS per D-06) against
a fixed threshold, reading from `benchmarks/results/oris_soft_vs_competitors.csv` or invoking
`harvest()` directly — follow the first test's `skip_if(Sys.getenv("LBW_BENCH_GATE") == "")`
pattern exactly (D-10/D-11: this new assertion must be opt-in, not always-on).

---

### `README.md` (new file, doc)

**No analog** — confirmed via `git log --all --full-history -- README.md` (empty, cited in
RESEARCH.md Sources). Use RESEARCH.md D-13's shape directly: one headline performance line +
link to `docs/performance.md`. Do not invent structure from a different project; this is a
first-of-its-kind file in this repo, scoped narrowly to the SC1/SC2 headline per RESEARCH.md
Finding 1 (do not let it balloon into a full package README — that's Phase 5's job).

---

### `docs/performance.md` (new, doc)

**Analog:** `docs/methods/oris.md` competitor table and citation style (lines 216-226,
253-267 read this session).

**Competitor table pattern to reference/link, not duplicate** (lines 216-226):
```markdown
| Package | Language | Function | Algorithm | Bounds support | Citation |
|---------|----------|----------|-----------|----------------|---------|
| `survey` | R | `rake()` / `calibrate()` | Raking (multiplicative IPF) + GREG + logit | Yes — `bounds=` arg; logit mandatory | [lumley2010survey] |
| `icarus` | R | `calibration()` | GREG + raking + logit (Calmar-inspired) | Yes — simplex / bisection | [rebecq2017icarus] |
| `ReGenesees` | R | `e.calibrate()` | GREG + raking + logit via Newton–Raphson | Yes — mandatory for logit | [zardetto2015regenesees] |
```

**Bounds-handling deviation row to cite (grounds the D-07 competitor selection)** (line 261):
```
| **Bounds handling** | ... bounds folded into every iteration via bisection or logit
transform (`survey`, `icarus`, `ReGenesees`) | Bounds deferred to `finalize_weights` after
convergence (ORIS_SOFT uses ALM/ADMM to enforce bounds inside the loop) | ... |
```

Per RESEARCH.md Open Question 2: link to `docs/methods/oris.md`'s existing citations
(`[zardetto2015regenesees]`, `[rebecq2017icarus]`, `[lumley2010survey]` via
`docs/methods/references.bib`) rather than restating the bibliography in the new page.

---

### `.planning/REQUIREMENTS.md` KPI-04 row (edit)

**Analog:** the same file's US-003 row, already rewritten with the correct "no live artefact"
framing this research session (lines 100-112, quoted above) — copy its style: bold
**Partial**/status marker, cite the specific investigation doc and commit hash, name the void
artefact explicitly (`test-lbfgsb.R` void → parallel for the new artefact name). The KPI-04
row should follow the same convention: name the NEW live artefact
(`tests/testthat/test-bench-gate.R`'s new assertion + `benchmarks/oris_soft_vs_competitors.R`)
in place of the removed `lbfgsb`-era one.

## Shared Patterns

### Single-thread BLAS determinism envelope
**Source:** CLAUDE.md build/test section; `benchmarks/stepstone_all_methods.R:3` comment
**Apply to:** `benchmarks/oris_soft_vs_competitors.R`, `benchmarks/run_honest_gate.sh`, and
the new `test-bench-gate.R` assertion (any invocation route)
```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript ...
```

### `bench::mark()` cross-solver timing loop (not interleaved before/after)
**Source:** `benchmarks/newton_kl_bench.R:28-33`
**Apply to:** `benchmarks/oris_soft_vs_competitors.R`
```r
bench::mark(run = harvest(...), iterations = 2, check = FALSE, memory = FALSE, filter_gc = FALSE)
```

### `LBW_BENCH_GATE=1` opt-in skip idiom (do not invent a second gate convention)
**Source:** `tests/testthat/test-bench-gate.R:1-17`
**Apply to:** any new assertion in `tests/testthat/test-bench-gate.R`
```r
skip_on_cran()
skip_if(Sys.getenv("LBW_BENCH_GATE") == "")
```

### Shell wrapper: two-step `echo && Rscript`, `set -euo pipefail`
**Source:** `benchmarks/run_allmethod.sh` (full file)
**Apply to:** `benchmarks/run_honest_gate.sh`

## No Analog Found

| File | Role | Data Flow | Reason |
|------|------|-----------|--------|
| `README.md` | doc | — | Has never existed in this repo (`git log` empty); build from RESEARCH.md D-13 shape, not a copied analog. |
| `requireNamespace()` guard for `icarus`/`ReGenesees` inside a benchmark script | utility/guard | — | No existing benchmark script in this repo guards an unregistered (non-`Suggests:`) competitor package; this is new code per D-09. |

## Metadata

**Analog search scope:** `benchmarks/`, `tests/testthat/`, `docs/methods/`, `.planning/`
**Files scanned:** `benchmarks/newton_kl_bench.R`, `benchmarks/run_allmethod.sh`,
`benchmarks/stepstone_benchmark.R`, `tests/testthat/test-bench-gate.R`,
`docs/methods/oris.md`, `.planning/REQUIREMENTS.md`
**Pattern extraction date:** 2026-08-15
