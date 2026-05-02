# ylsy CP+IPM Research Spike — Implementation Plan (rev 2)

**Date:** 2026-05-02
**Beads epic:** `leafblower-y2ls`
**Beads tasks:** `leafblower-y2ls.1` through `.7`
**Spec:** `docs/superpowers/specs/2026-05-02-ylsy-cp-ipm-spike-design.md` (rev 2; design-review-gate APPROVED iter 2 — all 5 reviewers PASS)
**Predecessor:** Epic-H WH-g (master `3265a53`); kk1204 empirical state in `bd memory kk1204-k-20-wh-g-empirical-auto-correctly`
**Revision history:** rev 1 → rev 2 closes 7 blockers from plan-review-gate iter 1 (Eigen include path; R package preinstall; ergodic-rate column; SKIPPED protocol; hermetic baselines; parquet pre-flight; WU-7 dep edge).

## Mechanism

**Target pattern:** package-root `research/` directory containing two standalone C++ solvers (Chambolle-Pock primal-dual; Interior-Point Newton) compiled into `research/leafblower_research.so` via standalone `Makefile`. R-callable through `dyn.load` + `.Call("cp_solve_R", ...)` and `.Call("ipm_solve_R", ...)`.

**Library dependencies:**

| Component | Source | Resolution |
|---|---|---|
| Eigen (header-only) | `RcppEigen` R package | `EIGEN_INCLUDE := $(shell Rscript -e "cat(system.file('include', package='RcppEigen'))")` injected into `research/Makefile`; `PKG_CXXFLAGS += -I$(EIGEN_INCLUDE)` |
| LAPACK | system | `R_ext/Lapack.h` + `-llapack -lblas` (matches `src/Makevars.in` `PKG_LIBS`) |
| BLAS thread control | `RhpcBLASctl` R package | preinstalled by WU-1 |
| Parquet I/O | `arrow` R package | already available |
| JSON I/O | `jsonlite` R package | already available |

**WU-1 prerequisite install step** (mandatory, NOT mirroring Makevars.in which has no Eigen):
```r
needed <- c("RcppEigen", "RhpcBLASctl")
missing <- needed[!needed %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
stopifnot(all(needed %in% rownames(installed.packages())))
```

**Audit strategy:** mechanical CI gate `tools/check_research_isolation.R` runs `nm -D src/leafblower.so` and asserts `cp_solve_R`, `ipm_solve_R`, `cp_calibrate`, `ipm_calibrate` symbols are NOT present (research must NEVER leak into the package `.so`). WU-1 also appends invocation to `.git/hooks/pre-commit`.

## Forbidden

- **No edits to** `src/leafblower.h`, `src/Makevars.in`, `src/r_bridge.cpp`, `R/harvest.R`, or any file under `tests/testthat/`, `NEWS.md` — the spike does NOT productionize.
- **No method enum entries** (`RK_ALG_CP`, `RK_ALG_IPM` are not added).
- **No AUTO routing changes**.
- **No accelerated PDHG variant** — vanilla Algorithm 1 only.
- **No Mehrotra predictor-corrector IPM** — vanilla central-path only.
- **No warm-start** between solvers.
- **No hyperparameter tuning** during the spike.
- **No raw `new`/`delete` in any C++ file** — `std::vector` / `Eigen` only.
- **No Lambert-W in CP prox** — direct Newton with overflow-asymptote guard.
- **No `--no-verify` on commits** — pre-commit hook runs the isolation gate.
- **No bundling of WUs** — one bd ticket per atomic deliverable.

## Audit (Spy/Mock Strategy)

- **Sanity recovery test** (WU-2 + WU-5): `benchmarks/research/sanity_t1_recovery.R` constructs a fixture where targets equal observed marginals → analytical solution `w = d`. Both CP and IPM must converge to `max_err < 1e-8`. Hard halt before kk1204 bench if sanity fails.
- **Isolation gate** (WU-1, re-run by WU-2/WU-5): `tools/check_research_isolation.R`.
- **Trace capture** (WU-3, WU-6): per-iter scalar aggregates ≤10 doubles, hard cap 1000 rows via `trace_stride`. **Two rate fits per fixture**: `lm(log(max_err_last) ~ log(iter))` AND `lm(log(max_err_ergodic) ~ log(iter))` (CP only); IPM has only `max_err_last`. Both with $R^2 \ge 0.9$ and $n_{\text{fit}} \ge 30$ floor.
- **R3 falsification check** (WU-7): do CP, IPM (if ran), ieppa+sraa, newton_kl on kk1204 agree within 1%? If yes → joint 5e-2 fixed point CONFIRMED.
- **Verdict reproducibility** (WU-7): `ylsy_compare.R` produces unified comparison CSV; trajectory plots saved to `benchmarks/research/results/plots/`.
- **Parquet pre-flight** (WU-3, WU-6): existence check on `benchmarks/stepstone_bench_data.parquet` and `benchmarks/stepstone_bench_targets.json` before running stepstone fixture; halt with explicit error if missing.

## Work Units (one bd ticket each)

| WU | Bead | Title | Hard deps | Model | Wall |
|---|---|---|---|---|---|
| WU-1 | `leafblower-y2ls.1` | Skeleton + isolation enforcement + R package install | — | Haiku | ~45 min |
| WU-2 | `leafblower-y2ls.2` | CP implementation | WU-1 | Gemini | ~3 h |
| WU-3 | `leafblower-y2ls.3` | CP bench (3 fixtures, 3× kk1204) | WU-2 | Gemini | ~1 h |
| WU-4 | `leafblower-y2ls.4` | CP verdict + WU-5/6 SKIP closure | WU-3 | Haiku | ~15 min |
| WU-5 (cond) | `leafblower-y2ls.5` | IPM implementation | WU-1, WU-4 | Gemini | ~4 h |
| WU-6 (cond) | `leafblower-y2ls.6` | IPM bench | WU-5 | Gemini | ~1 h |
| WU-7 | `leafblower-y2ls.7` | Investigation report (incl. inline baselines) | WU-4, WU-6 | Opus | ~3 h |

**Conditional skip path:** if WU-4 verdict ∈ {PASS, PASS-kk1204-specialist}, WU-4 itself runs `bd close leafblower-y2ls.5 leafblower-y2ls.6 --reason "SKIPPED: WU-4 verdict <V> per spec Sec 5"` to unblock WU-7. WU-7 reads `research/cp_verdict.txt` line 2 (`wu5_skip=true|false`) to decide whether to ingest IPM data.

**`bd dep` graph (post rev 2):** WU-2 → WU-1; WU-3 → WU-2; WU-4 → WU-3; WU-5 → {WU-1, WU-4}; WU-6 → WU-5; WU-7 → {WU-4, WU-6}. Closure of WU-5/WU-6 as SKIPPED satisfies WU-7's hard dep without running them.

**Total wall:** ~5–6 h if CP PASS at WU-4 (skip IPM); ~10–12 h if both ran.

## Decision Rule (verbatim from spec Sec 1)

| Outcome on kk1204 | Stepstone | Verdict | Action |
|---|---|---|---|
| max_err < 1e-3 AND walltime ≤ 30s AND (β_last ≤ -0.8 OR β_ergodic ≤ -1.0) AND R² ≥ 0.9 | within 1.5× of ieppa+sraa (≤ 1.7e-4) | **PASS** | File Epic-K (productionize) |
| max_err < 1e-3 AND walltime ≤ 30s AND rate as above | regression beyond 1.5× | **PASS-kk1204-specialist** | Productionize for severe-skew K≥5 only (carve-out scope-deferred to Epic-K) |
| max_err in [1e-3, 1e-2) AND walltime ≤ 60s | any | **PARTIAL** | One follow-up tuning spike permitted (chain depth = 1) |
| max_err ≥ 1e-2 OR walltime > 60s | any | **FAIL** | Investigation report; ylsy stays open or closes BLOCKED |

## SKIPPED Protocol (rev 2 fix)

Single source of truth: `research/cp_verdict.txt` line 2 = `wu5_skip=true` iff WU-4 verdict ∈ {PASS, PASS-kk1204-specialist}; else `wu5_skip=false`.

WU-4 (after writing the verdict file) executes:
```bash
if grep -q '^wu5_skip=true$' research/cp_verdict.txt; then
  bd close leafblower-y2ls.5 leafblower-y2ls.6 \
    --reason "SKIPPED: WU-4 verdict $(head -1 research/cp_verdict.txt) per spec Sec 5"
fi
```

WU-7 reads `wu5_skip` to decide IPM data ingestion; if true, ipm_summary.csv is absent and the report's quantitative table omits IPM rows.

## Hermetic Baselines for WU-7 (rev 2 fix)

WU-7 produces the comparison table by running baselines INLINE in `benchmarks/research/ylsy_compare.R`:

```r
# In-process baselines via package's harvest() — guarantees reproducibility on the bench host
for (fixture in c("t1_small", "stepstone_K9", "kk1204_K20")) {
  for (method in c("ieppa", "newton_kl", "lbfgsb")) {
    r <- suppressWarnings(harvest(fx[[fixture]]$df, fx[[fixture]]$tgt,
      method = method,
      max_weight = fx[[fixture]]$max_weight,
      attach_weights = FALSE,
      max_iterations = if (method == "ieppa") 50L else 200L,
      verbose = 0L,
      accelerate = (method == "ieppa")))
    res <- attr(r, "result")
    write_baseline_row(method, fixture, res)
  }
}
```

This eliminates the rev 1 hermetic gap ("manual entries from prior bench results") — every baseline is regenerated against `git rev-parse HEAD` at WU-7 run time.

## Reversibility / FAIL-artefact policy (spec Sec 7)

On FAIL or PARTIAL-without-followup verdict (closure):
- `research/` code STAYS on master under `.Rbuildignore` for traceability.
- WU-7 investigation report explicitly marks each prototype ABANDONED with `git_sha` of last `research/` commit.
- WU-7 updates bd `leafblower-ylsy` ticket with verdict + report path; orchestrator closes Epic-J + ylsy per verdict.

On PASS verdict:
- File Epic-K (productionization plan: method enum, AUTO carve-out, tests, docs).

## Risks & Mitigations (spec Sec 6, R1–R11; reproduced)

| # | Risk | Mitigation |
|---|---|---|
| R1 | CP step sizes too conservative | Power-iter ‖A‖ + 1.05 safety; report ergodic AND last-iterate err |
| R2 | IPM Schur rank-deficient | TSVD ratio 1e-8 (mirror Epic-Dβ WL-1 in `src/newton_calib.cpp:351`) |
| R3 | Both CP+IPM hit joint 5e-2 fixed point | Falsification: 1% agreement across 4 solvers |
| R4 | Faithful impl underperforms | PARTIAL → tune-and-retry follow-up (chain depth = 1) |
| R5 | research/ leaks into src/leafblower.so | Mechanical CI gate `tools/check_research_isolation.R` (WU-1) |
| R6 | OOM on kk1204 trace capture | Pre-flight memory check; trace ≤ 10 doubles/iter, capped 1000 rows |
| R7 | Walltime noise on 30s gate | 3× repetition median for kk1204; OMP=1, BLAS=1 |
| R8 | Eigen/LAPACK linkage drift | RcppEigen include path injected via `system.file()` (rev 2 fix); LAPACK matches Makevars.in `-llapack -lblas` |
| R9 | Power-iter divergence | Stop on rel-delta < 1e-6 OR k=50 with rel-delta > 1e-3 → status_code=4 |
| R10 | BLAS thread contention; float-precision cancellation | OMP=1, BLAS=1; all state in `double` |
| R11 | PROTECT/UNPROTECT imbalance crashes R | Manual reviewer audit; ASan+UBSan smoke target in research/Makefile |

## Success Criteria (Epic-J close)

- [ ] All 7 WU tickets closed (or WU-5/WU-6 closed SKIPPED via WU-4 step)
- [ ] `docs/investigations/2026-05-02-ylsy-cp-ipm-spike-result.md` committed
- [ ] `bd update leafblower-ylsy --notes "..."` records verdict + report path
- [ ] `bd close leafblower-y2ls --reason="<verdict>: <one-sentence summary>"` closes Epic-J
- [ ] PASS / PASS-specialist: file Epic-K productionization ticket
- [ ] FAIL: optionally `bd close leafblower-ylsy --reason BLOCKED`

## Post-merge action

Plan-review-gate iter 2 must PASS (all 3 adversarial reviewers) before any WU starts.
