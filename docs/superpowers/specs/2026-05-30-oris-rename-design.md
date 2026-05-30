# Rename `ieppa` / `ieppa_soft` → `oris` / `oris_soft` (ORIS)

**Status:** Draft (Design Review Gate pending)
**Date:** 2026-05-30
**Author:** Dennis Alexis Valin Dittrich
**Type:** Mechanical rename (no behavior change)
**Supersedes:** the `iEPPA` naming established in `docs/superpowers/specs/2026-04-23-ieppa-faithful-design.md`

---

## 1. Context & Motivation

Deep literature research (see `docs/methods/ieppa.md` "Relationship to the source paper" + `references.bib`) established that the solver currently named **iEPPA** is **not** the algorithm of Chu, Liang, Toh & Yang (arXiv:2011.14312). The shipped solver is the paper's *inner* dual block-coordinate-descent specialized to `C = 0` — i.e. **iterative matrix scaling (RAS / Sinkhorn–Knopp / IPF) on log-factors**, with **successive over-relaxation (SOR)** and an infeasibility-damping factor, bounds deferred to finalize. The paper's defining *outer* inexact-entropic-proximal-point loop is mathematically inert at `C = 0` and is not implemented. Carrying the `iEPPA` name over-claims a paper whose headline contribution this code does not contain.

**Decision:** rename to **ORIS — Over-Relaxed Iterative Scaling**. The name states the always-relevant lineage (iterative scaling = RAS/Sinkhorn–Knopp/IPF) and the flagship accelerator (over-relaxation) without naming a paper. The soft variant keeps the parallel suffix: **`oris_soft`** (ADMM/ALM soft-capacity variant).

No deprecation alias, no migration shim: the package has **no external users yet** (user-confirmed), so this is not a breaking change in practice.

## 2. Naming Decision

| Old | New | Notes |
|-----|-----|-------|
| concept "IEPPA — Iterative Entropy-Penalized Proportional Adjustment" | "ORIS — Over-Relaxed Iterative Scaling" | prose/title |
| `RK_ALG_IEPPA = 1` | `RK_ALG_ORIS = 1` | **enum value unchanged** |
| `RK_ALG_IEPPA_SOFT = 8` | `RK_ALG_ORIS_SOFT = 8` | **enum value unchanged** |
| `method = "ieppa"` | `method = "oris"` | R + Python user-facing string |
| `method = "ieppa_soft"` | `method = "oris_soft"` | R + Python |
| `ieppa_solve` | `oris_solve` | C++ symbol |
| `IEPPAResult` | `ORISResult` | C++ struct |

**Invariant: enum *values* (1, 8) are preserved.** Only identifiers/strings change. Slot 2 stays reserved (LBFGSB). This keeps serialized fixtures (`.rds`, parity CSVs) that store `algorithm_used` as an integer valid without regeneration.

## 3. Rename Surface (file-by-file)

### 3.1 C++ core (`src/`)
- `leafblower.h` — enum identifiers `RK_ALG_IEPPA`→`RK_ALG_ORIS`, `RK_ALG_IEPPA_SOFT`→`RK_ALG_ORIS_SOFT`; struct comments mentioning iEPPA.
- Rename files: `ieppa.cpp`→`oris.cpp`, `ieppa.hpp`→`oris.hpp`, `ieppa_internal.hpp`→`oris_internal.hpp`, `ieppa_finalize.cpp`→`oris_finalize.cpp`, `ieppa_trajectory.cpp`→`oris_trajectory.cpp`.
- Symbols: `ieppa_solve`→`oris_solve`, `IEPPAResult`→`ORISResult`, internal `ieppa_*` helpers/guards/include-guards → `oris_*` / `ORIS_*`.
- `#include "ieppa*.hpp"` updated at every include site.
- `calib_dispatch.hpp` — `case RK_ALG_IEPPA*` labels, dispatch calls, `select_solver_objective`/AUTO comments.
- `c_api.cpp` — `case RK_ALG_IEPPA*`, `ieppa_solve(...)` calls, `ieppa_auto_selected` field/flag, the algorithm-name string table entry `"ieppa"`→`"oris"` (and `"ieppa_soft"`), AUTO-fallback comments.

### 3.2 Build
- `python/CMakeLists.txt` — `CORE_SOURCES` list entries for the renamed `.cpp` files (R auto-globs `src/*.cpp`, so no R Makevars change required; verify `Makevars.in` PKG_SOURCES is decorative as documented).

### 3.3 R (`R/`)
- `harvest.R` — `match.arg(method, c("auto","ieppa","ieppa_soft", ...))` → `...,"oris","oris_soft",...`; the metric-routing map (`"auto" = "marginal_kl"`, any `"ieppa" = ...`, `"ieppa_soft" = ...` entries); any `method == "ieppa"` / `"ieppa_soft"` conditionals; the string→enum mapping passed to the C API; roxygen `@param method` docs + examples.
- Any other R files referencing the method strings (`current_miss.R`, `diagnose_weights.R`, etc.) — grep-driven.

### 3.4 Python (`python/`)
- Pybind bindings (`_bindings.cpp` / `python/leafblower/`) — method-string acceptance + any enum exposure; docstrings.

### 3.5 Tests / benchmarks / fixtures
- `tests/testthat/test-ieppa*.R` → `test-oris*.R`; assertions and `method=`/`algorithm=` literals.
- `test-stall-kind-*.R`, `test-rk-params-passthrough.R`, parity tests, and any test referencing `"ieppa"`/`RK_ALG_IEPPA` — grep-driven.
- `benchmarks/` fixtures and result keys referencing `ieppa` (stepstone regression gate).
- Python parity tests under `python/` referencing the method string.

### 3.6 Docs
- `docs/methods/ieppa.md` → `docs/methods/oris.md`; retitle "ORIS — Over-Relaxed Iterative Scaling"; update the `> Enum:`/`> Source:` header; keep the "Relationship to the source paper" section, adding one line that the name was changed *away from* iEPPA because the code implements the paper's inner scaling step, not its outer PPA contribution.
- `docs/methods/00-overview.md` — comparison table row, selection-guidance rows, enum-reference block, per-method-documents table link, any prose naming IEPPA.
- `docs/methods/references.bib` — unchanged (paper entry `chu2022ieppa` stays; it is still cited).
- `CLAUDE.md` — all `ieppa`/`IEPPA` mentions (TU-split note, SRAA best-iterate, the "two solver formulas" ALM-Newton note, etc.).
- `NEWS.md` — add a development entry recording the rename.
- `DESCRIPTION` / `python/pyproject.toml` — only if they name the method (likely not; verify).

### 3.7 Auto-maintained (update after file moves)
- `.wolf/anatomy.md`, `.wolf/memory.md` per OpenWolf protocol; graphify/`graphify update .`; `.beads` plan/issue references if any name the solver.

## 4. Out of Scope
- Any change to algorithm behavior, numeric output, enum values, or convergence.
- The `raking`, `sinkhorn`, `chebyshev`, `greg`, `greenkhorn`, `logit`, `newton_kl`, `auto` names.
- The historical design spec `2026-04-23-ieppa-faithful-design.md` — left intact as a record (spec amendments append, never delete history); a one-line header note may point to this rename.

## 5. Risks & Mitigations
- **R1 — missed reference breaks build/link.** Two build sites (R auto-glob; Python explicit `CORE_SOURCES`). *Mitigation:* exhaustive `grep -ri 'ieppa'` audit across the repo as the completion gate; `R CMD INSTALL --preclean .` AND `pip install -e .` both must succeed.
- **R2 — silent behavior change.** *Mitigation:* enum values frozen; diff must be identifier-only. Verify: R `devtools::test()` green, Python parity (rtol=1e-6) green, stepstone benchmark shows no regression (Definition of Done).
- **R3 — stale fixtures.** Parity CSVs/`.rds` store `algorithm_used` as int → unaffected by identifier rename (R2 invariant). Any fixture storing the *string* `"ieppa"` must be regenerated — grep fixtures for the literal.
- **R4 — docs/code drift.** *Mitigation:* update docs in the SAME commit as code (CLAUDE.md §4).

## 6. Verification / Definition of Done
1. `grep -ri 'ieppa' --include='*.{cpp,hpp,h,R,py,md,txt,in}'` returns only intentional historical references (the `chu2022ieppa` bib key, the dated 2026-04-23 spec, and the "renamed from iEPPA" notes).
2. `R CMD INSTALL --preclean .` succeeds.
3. `Rscript -e "devtools::test()"` — all green.
4. `cd python && pip install -e . && pytest` — parity green (rtol=1e-6).
5. Stepstone benchmark: no regression vs pre-rename baseline.
6. `RK_ALG_ORIS == 1` and `RK_ALG_ORIS_SOFT == 8` asserted (static_assert or test).

## 7. Rollout
Single atomic change set (rename + all references + docs + NEWS.md) so no intermediate state breaks the build. No alias, no deprecation cycle (no users). One ticket per work unit per `planning-with-beads`; this spec → `writing-plans` → plan-review-gate → orchestrated-execution.
