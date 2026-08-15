The project's philosophy is SOTA + cutting edge + absolute statistical correctness, proven if possible, being the most efficient and fastest.
Do: numeric stability.
Stop: cancelations.

# Project Instructions for AI Agents

This file provides instructions and context for AI coding agents working on this project.

<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:ca08a54f -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol

## Session Completion

**When ending a work session**, complete the steps below. This project has **NO git remote (local-only)** — work is complete when committed locally and all quality gates pass. Do NOT `git push`/`bd dolt push`; there is nothing to push to.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **Commit locally** - Atomic per-ticket commits with explicit pathspec (`git commit <path1> <path2> -m ...`) to keep `.beads/issues.jsonl` and stray build artifacts out of the commit. NO push (no remote).
5. **Clean up** - Clear stashes; remove stray build artifacts (e.g. a regenerated `man/dot-*.Rd`)
6. **Verify** - All changes committed locally; working tree clean
7. **Hand off** - Provide context for next session

**CRITICAL RULES:**
- Local-only: complete = committed locally + quality gates green. There is NO remote to push to.
- Commit with explicit pathspec, NEVER `git add -A` (a bd/graphify hook re-stages `.beads/issues.jsonl` into the index; pathspec commits only the named paths).
- When two tickets touch one file, split-commit (reverse one ticket's edits to a backup, commit the other, restore, commit the first).
<!-- END BEADS INTEGRATION -->


## Build & Test

```bash
R CMD INSTALL --preclean .                     # R build gate — use this, NOT devtools::install
Rscript -e "devtools::test()"                  # R tests (testthat v3)
cd python && uv pip install -e . --reinstall-package leafblower   # venv is uv-managed — NO pip
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/python -m pytest  # parity (rtol=1e-6); bare python/pytest loads a STALE ~/.local shadow .so
```

**Definition of Done:** R tests pass + Python parity tests pass + stepstone benchmark shows no regression.

## Architecture

- Single C++17 core (`src/`). C API in `leafblower.h`. R bridge in `r_bridge.cpp`. Python via scikit-build + CMake.
- `calib_dispatch.hpp` = canonical home for shared solver helpers. Do NOT add shared logic to individual solver files. CellTable-specific helpers → `cell_table.hpp`. Both use `lbw` namespace.
- `CalibResult` fields live at `res.base.*` (not `res.*`) since ztid.4 — direct field access breaks silently.
- Algorithm slots 2 and 7 are reserved (LBFGSB and GRAKE removed, respectively). Do not reuse in `rk_algorithm_t` enum.
- **R build sets no `-O` optimization level** of its own (`configure` and `Makevars.in` both note this intentionally): R supplies the user/site `-O` via `$(CXXFLAGS)` in `$(ALL_CXXFLAGS)`. CRAN portability check (`tools:::.check_make_vars`) rejects `-O` flags in `PKG_CXXFLAGS`, so a user wanting `-O3` sets it in `~/.R/Makevars`. The Python build has no such constraint and hard-sets `-O3` unconditionally (`python/CMakeLists.txt:99`) — this asymmetry is deliberate and documented, not equalized (phase-02 SC2, `leafblower-qzto`); R↔Python parity tests (`python/leafblower/test_solver_parity.py`) treat their tolerances as the bound on how much it may move a result.
- **Two build sites for `src/*.cpp`:** R auto-globs `src/*.cpp` (the `PKG_SOURCES` list in `Makevars.in` is decorative — R ignores it). The Python build does NOT glob: `python/CMakeLists.txt` has an explicit `CORE_SOURCES` list. A new `src/*.cpp` MUST be added to `CORE_SOURCES` or the pybind11 link fails with undefined symbols.
- **No LTO** (`-flto` absent from `configure`/`Makevars`). A TU split only stays perf-neutral if HOT code (per-iteration loops, `static inline` kernels) stays co-located with its caller — cross-TU calls don't inline. Move only COLD (once-per-solve) code to new TUs. `oris.cpp` is split this way (`oris_finalize.cpp`, `oris_trajectory.cpp`, `oris_internal.hpp`); the hot `oris_solve` stays in `oris.cpp`.

## Conventions & Footguns

**Adding a new solver — all 8 steps required:**
1. `src/<name>.cpp` + `src/<name>.hpp` in `lbw` namespace — add the `.cpp` to `python/CMakeLists.txt` `CORE_SOURCES` (R auto-globs; Python does NOT)
2. Add enum value to `rk_algorithm_t` in `leafblower.h`
3. Wire shared helpers through `calib_dispatch.hpp`
4. R wrapper function + roxygen2 docs
5. Python binding
6. R test fixture (`.rds`)
7. Python parity test in `python/leafblower/test_*_parity.py` (e.g. `test_solver_parity.py`)
8. Benchmark fixture for stepstone regression gate

- **Version sync:** bump `DESCRIPTION` AND `python/pyproject.toml` manually — no automation.
- **Output weights:** Σw=n enforced at exit via `lbw::finalize_weights[_buf]` (calib_dispatch.hpp): a single pre-bounds scale to n, THEN `bounds_mode` dispatch (cell = count-only, unit = per-cell water-fill). The sanctioned order is normalize→bounds. FORBIDDEN: renormalizing AFTER water-fill — that silently breaks the `bounds_mode="unit"` clamps. Degenerate `total_w < kMinSafeTotalWeight` (1e-100) is left unscaled (avoids subnormal-overflow).
- **`bounds_mode`:** `"cell"` = cell-aggregate (default). `"unit"` = strict per-obs via ORIS water-fill.
- `homotopy_levels_used` returns `1` for `n_levels=1` (single-pass), not `0`. Struct comment is wrong.
- SRAA best-iterate: use `select_metric(sraa_cfg.metric, cm)` at `kErrCheckInterval` — NOT `errRp` fast proxy. Bug has been re-introduced twice.
- Do NOT "fix" two solver formulas (both verified correct, guarded by code comments): the chebyshev Mehrotra corrector's linear `y·Δs_aff` term (it's the `−Δs_aff·Δy_aff` cross-term, not a stray residual), and the oris ALM Newton step `X̃(1−λ+μz)/(1+ρ)` (correct for the un-normalized-KL generator — no missing `−ρ`). Reviewer claims to add/drop terms came from the wrong divergence/derivation; verify against the actual generator before touching solver math.
- Lambda `[&]` in `raking.cpp`: declare bool guards BEFORE `auto F_eval = [&]` definition — `[&]` captures only vars in scope at definition site, not at call site.
- **`design_weights=`** is the `harvest()` argument for per-observation design weights (the `d_i` in `Z=Σ d_i exp(u_i)`); there is NO `weights=` argument — it silently lands in `...` and is ignored.
- **Deterministic tests/parity require single-thread BLAS:** export `OMP_NUM_THREADS`/`OPENBLAS_NUM_THREADS`/`MKL_NUM_THREADS`=1 TOGETHER (in Python, before `import numpy`), else R↔Python parity and benchmarks drift.
- **No cancellations (project philosophy):** compute variances/covariances cancellation-free — e.g. the diagonal `p(1-p)` as `p*(1-p)`, NOT `p - p*p`. Exclude zero-design-weight rows from LSE sums (a `0*inf` on a divergent shift).
- **Module-side file writes MUST use C stdio** (`fopen`/`fprintf`/`fclose`), NOT `std::ofstream`/`iostream`: the Python `_leafblower.so` static-links libstdc++, so any C++ stream I/O SIGSEGVs in `std::codecvt do_unshift` when driven from the module (`oris_trajectory.cpp write_trajectory_csv`, leafblower-9nuo).

# Metaswarm

Configured via `/metaswarm:setup` (profile: `.metaswarm/project-profile.json`). Minimal scaffolding — the repo's existing `CLAUDE.md`, `.beads`, `.wolf`, and `model-routing` conventions remain authoritative.

- **Quality gate is BEHAVIORAL, not line-coverage.** The orchestrated-execution VALIDATE phase reads `.coverage-thresholds.json` → `enforcement.command`, which runs the real Definition of Done: `R CMD INSTALL --preclean .` + R testthat (0 FAIL) + Python `pytest` (0 FAIL), single-thread BLAS. Stepstone no-regression is opt-in via `LBW_BENCH_GATE=1` (heavy, local-only). Do NOT add a `--cov-fail-under` gate — no `covr`/`pytest-cov` is wired and it would only cover the Python layer.
- **External tools:** `.metaswarm/external-tools.yaml` — `gemini` enabled (CLI present), `codex` disabled (not installed). Delegated implementer = Gemini, consistent with the `model-routing` skill (Tier-2 → Gemini, Tier-3 → Opus; route web research to `agy-delegate`).
- **Orchestrated execution:** the 4-phase loop (implement → orchestrator independently validates → FRESH adversarial reviewer → commit) is the standing methodology; pair with `planning-with-beads` for 3+-file work. Commit with explicit pathspec to avoid `.beads/issues.jsonl` hook-leak; no git remote (local-only).
