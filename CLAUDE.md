# OpenWolf

@.wolf/OPENWOLF.md

This project uses OpenWolf for context management. Read and follow .wolf/OPENWOLF.md every session. Check .wolf/cerebrum.md before generating code. Check .wolf/anatomy.md before reading files.


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
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

## Session Completion

**When ending a work session**, you MUST complete ALL steps below. Work is NOT complete until `git push` succeeds.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **PUSH TO REMOTE** - This is MANDATORY:
   ```bash
   git pull --rebase
   bd dolt push
   git push
   git status  # MUST show "up to date with origin"
   ```
5. **Clean up** - Clear stashes, prune remote branches
6. **Verify** - All changes committed AND pushed
7. **Hand off** - Provide context for next session

**CRITICAL RULES:**
- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds
<!-- END BEADS INTEGRATION -->


## Build & Test

```bash
R CMD INSTALL --preclean .                     # R build gate — use this, NOT devtools::install
Rscript -e "devtools::test()"                  # R tests (testthat v3)
cd python && pip install -e . && pytest        # Python parity tests (rtol=1e-6 vs R output)
```

**Definition of Done:** R tests pass + Python parity tests pass + stepstone benchmark shows no regression.

## Architecture

- Single C++17 core (`src/`). C API in `leafblower.h`. R bridge in `r_bridge.cpp`. Python via scikit-build + CMake.
- `calib_dispatch.hpp` = canonical home for shared solver helpers. Do NOT add shared logic to individual solver files. CellTable-specific helpers → `cell_table.hpp`. Both use `lbw` namespace.
- `CalibResult` fields live at `res.base.*` (not `res.*`) since ztid.4 — direct field access breaks silently.
- Algorithm slot 2 is reserved (LBFGSB removed). Do not reuse in `rk_algorithm_t` enum.
- The package sets **no `-O` optimization level** of its own (`configure` and `Makevars.in` both note this intentionally): R supplies the user/site `-O` via `$(CXXFLAGS)` in `$(ALL_CXXFLAGS)`. CRAN portability check (`tools:::.check_make_vars`) rejects `-O` flags in `PKG_CXXFLAGS`, so a user wanting `-O3` sets it in `~/.R/Makevars`.
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

<!-- code-review-graph MCP tools -->
## MCP Tools: code-review-graph

**IMPORTANT: This project has a knowledge graph. ALWAYS use the
code-review-graph MCP tools BEFORE using Grep/Glob/Read to explore
the codebase.** The graph is faster, cheaper (fewer tokens), and gives
you structural context (callers, dependents, test coverage) that file
scanning cannot.

### When to use graph tools FIRST

- **Exploring code**: `semantic_search_nodes` or `query_graph` instead of Grep
- **Understanding impact**: `get_impact_radius` instead of manually tracing imports
- **Code review**: `detect_changes` + `get_review_context` instead of reading entire files
- **Finding relationships**: `query_graph` with callers_of/callees_of/imports_of/tests_for
- **Architecture questions**: `get_architecture_overview` + `list_communities`

Fall back to Grep/Glob/Read **only** when the graph doesn't cover what you need.

### Key Tools

| Tool | Use when |
|------|----------|
| `detect_changes` | Reviewing code changes — gives risk-scored analysis |
| `get_review_context` | Need source snippets for review — token-efficient |
| `get_impact_radius` | Understanding blast radius of a change |
| `get_affected_flows` | Finding which execution paths are impacted |
| `query_graph` | Tracing callers, callees, imports, tests, dependencies |
| `semantic_search_nodes` | Finding functions/classes by name or keyword |
| `get_architecture_overview` | Understanding high-level codebase structure |
| `refactor_tool` | Planning renames, finding dead code |

### Workflow

1. The graph auto-updates on file changes (via hooks).
2. Use `detect_changes` for code review.
3. Use `get_affected_flows` to understand impact.
4. Use `query_graph` pattern="tests_for" to check coverage.

## graphify

This project has a graphify knowledge graph at graphify-out/.

Rules:
- Before answering architecture or codebase questions, read graphify-out/GRAPH_REPORT.md for god nodes and community structure
- If graphify-out/wiki/index.md exists, navigate it instead of reading raw files
- For cross-module "how does X relate to Y" questions, prefer `graphify query "<question>"`, `graphify path "<A>" "<B>"`, or `graphify explain "<concept>"` over grep — these traverse the graph's EXTRACTED + INFERRED edges instead of scanning files
- After modifying code files in this session, run `graphify update .` to keep the graph current (AST-only, no API cost)
