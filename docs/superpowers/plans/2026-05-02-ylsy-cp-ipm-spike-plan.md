# ylsy CP+IPM Research Spike — Implementation Plan

**Date:** 2026-05-02
**Beads epic:** `leafblower-y2ls`
**Beads tasks:** `leafblower-y2ls.1` through `.7`
**Spec:** `docs/superpowers/specs/2026-05-02-ylsy-cp-ipm-spike-design.md` (rev 2; design-review-gate APPROVED iter 2 — all 5 reviewers PASS)
**Predecessor:** Epic-H WH-g (master `3265a53`); kk1204 empirical state captured in `bd memory kk1204-k-20-wh-g-empirical-auto-correctly`

## Mechanism

**Target pattern:** package-root `research/` directory containing two standalone C++ solvers (Chambolle-Pock primal-dual; Interior-Point Newton) compiled into `research/leafblower_research.so` via standalone `Makefile`. R-callable through `dyn.load` + `.Call("cp_solve_R", ...)` and `.Call("ipm_solve_R", ...)`.

**Library dependencies:** Eigen (header-only, mirror existing `Makevars.in` linkage); LAPACK via `R_ext/Lapack.h` (`dsyevd` for IPM Schur TSVD).

**Audit strategy:** mechanical CI gate `tools/check_research_isolation.R` runs `nm -D src/leafblower.so` and asserts `cp_solve_R`, `ipm_solve_R`, `cp_calibrate`, `ipm_calibrate` symbols are NOT present (research must NEVER leak into the package `.so`). Wired into pre-commit hook by WU-1.

## Forbidden

- **No edits to** `src/leafblower.h`, `src/Makevars.in`, `src/r_bridge.cpp`, `R/harvest.R`, or any file under `tests/testthat/` — the spike does NOT productionize.
- **No method enum entries** (`RK_ALG_CP`, `RK_ALG_IPM` are not added).
- **No AUTO routing changes**.
- **No accelerated PDHG variant** (Algorithm 2 of Chambolle-Pock 2011) — vanilla Algorithm 1 only.
- **No Mehrotra predictor-corrector IPM** — vanilla central-path only.
- **No warm-start** between solvers (cold start only).
- **No hyperparameter tuning** during the spike — productionization phase only.
- **No raw `new`/`delete` in any C++ file** — `std::vector` / `Eigen` only.
- **No Lambert-W in CP prox** — direct Newton with overflow-asymptote guard.
- **No `--no-verify` on commits** — pre-commit hook runs the isolation gate; if it fails, fix root cause.
- **No bundling of WUs** — one bd ticket per atomic deliverable.

## Audit (Spy/Mock Strategy)

- **Sanity recovery test** (WU-2 + WU-5): `benchmarks/research/sanity_t1_recovery.R` constructs a fixture where targets equal observed marginals → analytical solution is `w = d`. Both CP and IPM must converge to `max_err < 1e-8`. Hard halt before kk1204 bench if sanity fails.
- **Isolation gate** (WU-1, re-run by WU-2/WU-5): `tools/check_research_isolation.R`. Fails if research symbols leak into `src/leafblower.so`.
- **Trace capture** (WU-3, WU-6): per-iter scalar aggregates ≤10 doubles, hard cap 1000 rows via `trace_stride`. Rate-fitting via `lm(log(max_err) ~ log(iter), data = subset)` with $R^2 \ge 0.9$ + $n_{\text{fit}} \ge 30$ floor.
- **R3 falsification check** (WU-7): do CP, IPM, ieppa+sraa, newton_kl on kk1204 agree within 1%? If yes → joint 5e-2 fixed point CONFIRMED (basin floor fundamental). If no → R3 ruled out, FAIL is solver-specific.
- **Verdict reproducibility** (WU-7): `ylsy_compare.R` produces unified comparison CSV; trajectory plots saved to `benchmarks/research/results/plots/`. Reader can re-fit rate exponents from raw trace CSVs.

## Work Units (one bd ticket each — see Epic `leafblower-y2ls`)

| WU | Bead | Title | Deps | Model | Wall budget |
|---|---|---|---|---|---|
| WU-1 | `leafblower-y2ls.1` | Skeleton + isolation enforcement | — | Haiku | ~30 min |
| WU-2 | `leafblower-y2ls.2` | CP implementation | WU-1 | Gemini | ~3 h |
| WU-3 | `leafblower-y2ls.3` | CP bench (3 fixtures, 3× kk1204) | WU-2 | Gemini | ~1 h (≤30s × kk1204 × 3 + setup) |
| WU-4 | `leafblower-y2ls.4` | CP verdict | WU-3 | Haiku | ~10 min |
| WU-5 (cond) | `leafblower-y2ls.5` | IPM implementation | WU-1, WU-4 ∈ {PARTIAL, FAIL} | Gemini | ~4 h |
| WU-6 (cond) | `leafblower-y2ls.6` | IPM bench | WU-5 | Gemini | ~1 h |
| WU-7 | `leafblower-y2ls.7` | Investigation report | WU-3 (CP) ∨ WU-6 (IPM) | Opus | ~2 h |

**Conditional skip path:** if WU-4 verdict ∈ {PASS, PASS-kk1204-specialist}, WU-5 + WU-6 are SKIPPED; WU-7 runs immediately on CP-only data.

**Total wall:** ~5–6 h if CP PASS at WU-4 (skip IPM); ~10–12 h if both ran.

## Decision Rule (spec Sec 1, copied for hermetic reference)

| Outcome on kk1204 | Stepstone | Verdict | Action |
|---|---|---|---|
| max_err < 1e-3 AND walltime ≤ 30s AND rate exp β ≤ −0.8 (last) or ≤ −1.0 (ergodic) AND R² ≥ 0.9 | within 1.5× of ieppa+sraa (≤ 1.7e-4) | **PASS** | File Epic-K (productionize: method enum + tests + AUTO + docs) |
| max_err < 1e-3 AND walltime ≤ 30s AND rate exp / R² as above | regression beyond 1.5× ieppa+sraa | **PASS-kk1204-specialist** | Productionize ONLY for severe-skew K≥5 routing (AUTO carve-out, scope-deferred to Epic-K) |
| max_err in [1e-3, 1e-2) AND walltime ≤ 60s | any | **PARTIAL** | One follow-up tuning spike permitted (chain depth cap = 1) |
| max_err ≥ 1e-2 OR walltime > 60s | any | **FAIL** | Investigation report; ylsy stays open or closes BLOCKED |

## Reversibility / FAIL-artefact policy (spec Sec 7)

On FAIL or PARTIAL-without-followup verdict (closure):
- `research/` code STAYS on master under `.Rbuildignore` for traceability and future-spike baselining (no orphan deletion).
- WU-7 investigation report explicitly marks each prototype ABANDONED with `git_sha` of last `research/` commit.
- WU-7 updates bd `leafblower-ylsy` ticket: `bd update leafblower-ylsy --notes "Epic-J spike verdict: <V>; report: ..."`. Status FAIL → ylsy stays open with comment; status BLOCKED → ylsy closes BLOCKED.
- Future spikes can resume from captured trace CSVs without re-implementation.

On PASS verdict:
- File Epic-K (productionization plan: method enum, full integration, AUTO routing carve-out, tests, docs).
- `research/` code remains as authoritative reference until productionization fully replaces it; only then archived.

## Risks & Mitigations (spec Sec 6, R1–R11; reproduced for plan completeness)

| # | Risk | Mitigation |
|---|---|---|
| R1 | CP step sizes too conservative | Power-iter ‖A‖ + 1.05 safety; report ergodic AND last-iterate err |
| R2 | IPM Schur rank-deficient | TSVD ratio 1e-8 (mirror Epic-Dβ WL-1) |
| R3 | Both CP+IPM hit joint 5e-2 fixed point | Falsification: 1% agreement across 4 solvers; ylsy closes BLOCKED |
| R4 | Faithful impl underperforms; needs production-grade variants | PARTIAL → tune-and-retry follow-up (chain depth = 1) |
| R5 | research/ leaks into src/leafblower.so | Mechanical CI gate `tools/check_research_isolation.R` (WU-1) |
| R6 | OOM on kk1204 trace capture | Pre-flight memory check; trace ≤ 10 doubles/iter, capped 1000 rows |
| R7 | Walltime noise on 30s gate | 3× repetition median for kk1204; OMP=1, BLAS=1 |
| R8 | Eigen/LAPACK linkage drift | Mirror `Makevars.in` flags exactly; lock via Makefile comment header |
| R9 | Power-iter divergence | Stop on rel-delta < 1e-6 OR k=50 with rel-delta > 1e-3 → status_code=4 |
| R10 | BLAS thread contention; float-precision cancellation | OMP=1, BLAS=1; all state in `double` |
| R11 | PROTECT/UNPROTECT imbalance crashes R | Manual reviewer audit; ASan+UBSan smoke target in Makefile |

## Success Criteria (Epic-J close)

- [ ] All 7 WU tickets closed (or WU-5/WU-6 closed SKIPPED with note)
- [ ] `docs/investigations/2026-05-02-ylsy-cp-ipm-spike-result.md` committed with verdict label
- [ ] `bd update leafblower-ylsy --notes "..."` records the verdict + report path
- [ ] `bd close leafblower-y2ls --reason="<verdict>: <one-sentence summary>"` closes Epic-J
- [ ] If verdict = PASS or PASS-kk1204-specialist: file Epic-K (productionization) bd ticket as next step
- [ ] If verdict = FAIL: `bd close leafblower-ylsy --reason="..."` closes ylsy BLOCKED (or stays open with summary if user wants follow-up)

## Post-merge action

Run `/plan-review-gate` on this plan before any WU starts. All 3 adversarial reviewers (Feasibility, Completeness, Scope & Alignment) must PASS.
