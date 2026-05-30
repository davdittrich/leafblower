# ORIS Rename Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rename the `ieppa` / `ieppa_soft` solver to `oris` / `oris_soft` (ORIS — Over-Relaxed Iterative Scaling) across code, API, tests, docs, and benchmarks, with zero behavior change.

**Architecture:** Pure identifier/string/file rename. Enum *values* (`RK_ALG_ORIS = 1`, `RK_ALG_ORIS_SOFT = 8`) are frozen so serialized fixtures stay valid. Single atomic change set (a partial C++ rename will not compile); work-unit phases below culminate in one verification gate (WU5). No deprecation alias (no users).

**Tech Stack:** C++17 core (R auto-globs `src/*.cpp`; Python build uses explicit `CORE_SOURCES`), R (roxygen2/testthat), Python (pybind11/pytest).

**Plan rev 2** (plan-review-gate round 1): enumerated the full test-suite touch set incl. `attr(,"algorithm")` assertions + `test-simd-math.R`/`task2_ieppa_ref.rds`; added `R/anesrake.R` roxygen; added `benchmarks/plot_helpers.R`; excluded generated `graphify-out/` from the grep gate and run `graphify update .` first.
**Plan rev 3** (plan-review-gate round 2): corrected `task2_ieppa_ref.rds` path (`tests/testthat/`, not `fixtures/`; opaque blob, no generator → move as-is); fixed dispatch-test accessor to `attr(r,"algorithm")` (was `algorithm_used`); added `src/raking.cpp` to WU1 comment-update list; scoped WU3's test grep to `tests/ python/leafblower/` and named `python/leafblower/test_python.py`.

**Authoritative scope:** `docs/superpowers/specs/2026-05-30-oris-rename-design.md` §3 + §8. Historical records are NOT edited (§8.10). This plan implements that spec; the binding completion check is the §8.10 grep-clean gate.

**Execution note:** Work on branch `oris-rename`. Per-WU commits are fine *on the branch*; the branch **tip** must be green (WU5). Because the rename is atomic at the C++ level, WU1 will not compile until complete — do the whole WU before building. Recommended: implement WU1→WU4 as edits, then run WU5 once; squash on merge if desired.

---

## File Structure (decomposition)

| Work unit | Responsibility | Build/verify at end |
|-----------|----------------|---------------------|
| WU1 | C++ core: enum, file renames, symbols, dispatch, field rename, strings/comments, build lists | `R CMD INSTALL --preclean .` compiles |
| WU2 | R + Python API surface: method strings, defaults, warnings, docstrings, pyproject | R loads; `match.arg` accepts `oris`/`oris_soft` |
| WU3 | Tests + fixtures + generators + new dispatch tests + enum assert | `devtools::test()` green; pytest parity green |
| WU4 | Docs + live benchmarks + man regeneration | `devtools::document()`; grep-clean of live docs |
| WU5 | Verification gate + atomic commit + repo-map refresh | full DoD (§6 + §8.8) passes |

---

## Task WU1 — C++ core rename

**Files (rename via `git mv`, then edit contents):**
- `git mv src/ieppa.cpp src/oris.cpp`; `ieppa.hpp→oris.hpp`; `ieppa_internal.hpp→oris_internal.hpp`; `ieppa_finalize.cpp→oris_finalize.cpp`; `ieppa_trajectory.cpp→oris_trajectory.cpp`
- Modify: `src/leafblower.h`, `src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`, `src/types.hpp`, `src/calib_validate.cpp`, `src/chebyshev.cpp`, `src/chebyshev.hpp`, `src/greg.cpp`, `src/logit_calib.cpp`, `src/newton_calib.cpp`, `src/newton_calib.hpp`, `src/sraa.hpp`, `src/cell_table.hpp`, `src/calib_validate.hpp`, `python/CMakeLists.txt`, `src/Makevars.in`

- [ ] **Step 1: Rename the 5 TU files** with `git mv` (preserves history).
```bash
cd src
for s in ieppa.cpp:oris.cpp ieppa.hpp:oris.hpp ieppa_internal.hpp:oris_internal.hpp ieppa_finalize.cpp:oris_finalize.cpp ieppa_trajectory.cpp:oris_trajectory.cpp; do git mv "${s%%:*}" "${s##*:}"; done
```

- [ ] **Step 2: Enum identifiers (values frozen)** in `src/leafblower.h`:
  - `RK_ALG_IEPPA = 1,` → `RK_ALG_ORIS = 1,`
  - `RK_ALG_IEPPA_SOFT = 8,` → `RK_ALG_ORIS_SOFT = 8,`
  - update the `/* 2 = removed (was RK_ALG_LBFGSB) */` neighbourhood comments mentioning iEPPA only if present; keep slot-2 reserved note.

- [ ] **Step 3: Symbols** (repo-wide within `src/`):
  - `ieppa_solve` → `oris_solve`; `IEPPAResult` → `ORISResult`; include guards `LBW_IEPPA_INTERNAL_HPP` → `LBW_ORIS_INTERNAL_HPP`; any `ieppa_*` internal helper names → `oris_*`.
  - `#include "ieppa*.hpp"` → `#include "oris*.hpp"` at every site (c_api.cpp, r_bridge.cpp, calib_dispatch.hpp, oris_finalize.cpp, oris_trajectory.cpp, oris_internal.hpp, newton_calib.cpp warm-start include if any).

- [ ] **Step 4: Dispatch + field**:
  - `src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`: `case RK_ALG_IEPPA*` → `RK_ALG_ORIS*`; `lbw::ieppa_solve(...)` → `lbw::oris_solve(...)`; alg-name string-table entries `"ieppa"`→`"oris"`, `"ieppa_soft"`→`"oris_soft"`.
  - `src/types.hpp`: rename field `ieppa_auto_selected` → `oris_auto_selected`; update all readers in `c_api.cpp`/`r_bridge.cpp`/`oris.cpp`. Fix `use_admm_capacity`/`capacity_mu` comments naming ieppa_soft.

- [ ] **Step 5: User-facing C string + comments**:
  - `src/calib_validate.cpp`: error string `"use method='ieppa' or 'raking'"` → `"use method='oris' or 'raking'"`.
  - Comment-only `ieppa`→`oris`/`ORIS` in: `raking.cpp` (`// … same API as ieppa`), `chebyshev.cpp/.hpp`, `greg.cpp`, `logit_calib.cpp`, `newton_calib.cpp/.hpp`, `sraa.hpp`, `cell_table.hpp`, `calib_validate.hpp`. (Run `grep -rIl -i ieppa src/ | grep -vE 'oris'` after WU1 to confirm none remain.)

- [ ] **Step 6: Build lists**:
  - `python/CMakeLists.txt` `CORE_SOURCES`: `ieppa.cpp`,`ieppa_trajectory.cpp`,`ieppa_finalize.cpp` → `oris.cpp`,`oris_trajectory.cpp`,`oris_finalize.cpp`.
  - `src/Makevars.in` `PKG_SOURCES`: same rename (decorative, for grep-cleanliness).

- [ ] **Step 7: Add enum-freeze guard.** In a compiled TU (e.g. top of `src/c_api.cpp` after the include of leafblower.h):
```cpp
static_assert(RK_ALG_ORIS == 1 && RK_ALG_ORIS_SOFT == 8, "ORIS enum values frozen for fixture compatibility");
```

- [ ] **Step 8: Compile gate.**
Run: `R CMD INSTALL --preclean .`
Expected: builds with no undefined-symbol/link errors. (If link error mentions `ieppa_solve` → a call site was missed; grep `src/` for `ieppa`.)

- [ ] **Step 9: Commit.**
```bash
git add -A src/ python/CMakeLists.txt
git commit -m "refactor(oris)!: rename ieppa C++ core → oris (enum values frozen)"
```

---

## Task WU2 — R + Python API surface

**Files:** Modify `R/harvest.R`, `R/anesrake.R`, `python/leafblower/_harvest.py`, `python/pyproject.toml`.

- [ ] **Step 1: `R/harvest.R` method enum + default**:
  - `match.arg(method, c("auto","ieppa","ieppa_soft", ...))` → `c("auto","oris","oris_soft", ...)` with `"oris"` in position 2 (no other token starts with `o` — verify).
  - default arg `method = "ieppa"` → `method = "oris"`.
  - string→enum mapping passed to the C API: `"ieppa"`→`RK_ALG_ORIS` path, `"ieppa_soft"`→`RK_ALG_ORIS_SOFT`.
  - metric-route map entries keyed on `"ieppa"`/`"ieppa_soft"` → `"oris"`/`"oris_soft"`.

- [ ] **Step 2: `R/harvest.R` user-visible warning strings**:
  - `method='ieppa+accel'` → `method='oris+accel'`
  - `method='ieppa_soft'` → `method='oris_soft'`
  - `"method='raking' or 'ieppa'"` → `"method='raking' or 'oris'"`

- [ ] **Step 3: `R/harvest.R` roxygen**: `@param method` — introduce "(ORIS: Over-Relaxed Iterative Scaling)" at first mention; `@param sor` doc "iEPPA only" → "ORIS (and raking)"; examples using `method="ieppa"` → `"oris"`.

- [ ] **Step 4: `R/anesrake.R`**: legacy remap `choosemethod <- "ieppa"` → `"oris"` (code, ~line 39); **also update the roxygen block** that documents the remap (the `\code{"ieppa"}` "recommended method" line, ~line 11) → `"oris"`, so `man/anesrake.Rd` regenerates clean in WU4.

- [ ] **Step 5: Python `_harvest.py`**: mirror the three warning strings (Step 2) and docstrings ("iEPPA only" → "ORIS"); method-string acceptance/mapping `oris`/`oris_soft`.

- [ ] **Step 6: `python/pyproject.toml`**: `description` "iEPPA" → "ORIS".

- [ ] **Step 7: Load gate.**
Run: `Rscript -e 'devtools::load_all("."); match.arg("oris", c("auto","oris","oris_soft","raking","sinkhorn","chebyshev","greg","greenkhorn","logit","newton_kl"))'`
Expected: returns `"oris"` (no error).

- [ ] **Step 8: Commit.**
```bash
git add R/harvest.R R/anesrake.R python/leafblower/_harvest.py python/pyproject.toml
git commit -m "refactor(oris)!: rename method strings + warnings in R/Python API"
```

---

## Task WU3 — Tests, fixtures, generators

**Files:** rename 8 `test-ieppa*.R`; `tests/parity/run_ieppa_soft_r.R`; fixtures `tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds`, `ieppa_pre_alm_ref.rds`, **`task2_ieppa_ref.rds`**; generators `data-raw/gen_ieppa_kl_ref.R`, `gen_ieppa_pre_alm_ref.R`; **plus content-only edits to every other test that references `ieppa`** (enumerated in Step 2b — these are NOT renamed, only edited).

- [ ] **Step 1: Rename test + parity + generator files** with `git mv`:
```bash
cd tests/testthat; for f in test-ieppa.R test-ieppa-bounds-mode.R test-ieppa-faithful.R test-ieppa-nonuniform-d.R test-ieppa-persistent-infeas.R test-ieppa-sraa.R test-ieppa-sraa-log-path.R test-ieppa-sraa-sor.R; do git mv "$f" "${f/ieppa/oris}"; done; cd -
git mv tests/parity/run_ieppa_soft_r.R tests/parity/run_oris_soft_r.R
git mv data-raw/gen_ieppa_kl_ref.R data-raw/gen_oris_kl_ref.R
git mv data-raw/gen_ieppa_pre_alm_ref.R data-raw/gen_oris_pre_alm_ref.R
```

- [ ] **Step 2: Content edits in the renamed files** — in the 8 renamed `test-oris*.R`, `run_oris_soft_r.R`, and `gen_oris_*.R`: replace `"ieppa"`/`"ieppa_soft"` → `"oris"`/`"oris_soft"`, `RK_ALG_IEPPA*`→`RK_ALG_ORIS*`, `test_that` descriptions/prose, and `saveRDS` paths → `oris_*.rds`.

- [ ] **Step 2b: Enumerate-and-edit ALL other tests referencing `ieppa`** (content only — these files are NOT renamed). First list them:
```bash
grep -rIl -i ieppa tests/ python/leafblower/ | grep -vE 'test-oris|run_oris'
```
Edit each. As of writing, the set is (verify with the grep above — which now ALSO scopes `python/leafblower/` — and edit whatever it returns):
`test-algo-selection.R`, `test-auto-routing-severe-skew.R`, `test-alm-config-grouping.R`, `test-best-iterate.R`, `test-bench-gate.R`, `test-bounded-convergence.R`, `test-calib-linalg.R`, `test-calib-result-consolidation.R`, `test-calibration-solvers.R`, `test-compare.R`, `test-config-defaults.R`, `test-convergence-criteria.R`, `test-convergence-trajectory.R`, `test-eta-schedule.R`, `test-harvest.R`, `test-homotopy.R`, `test-homotopy-enabled-field.R`, `test-method-dispatch.R`, `test-priority-sweep.R`, `test-quality-metrics.R`, `test-rk-params-passthrough.R`, `test-safety.R`, `test-select-metric-struct.R`, `test-simd-math.R`, `test-sor.R`, `tests/test_parity_weights.py`, **`python/leafblower/test_python.py`** (live `method="ieppa"` at lines ~48/59/70/80 + `("ieppa", …)`/`("ieppa_soft", …)` tuples ~106-107 + `method in ("ieppa", …)` ~158), `tests/testthat/fixtures/stepstone_reference_run.R`, plus `tests/testthat/_problems/*` if present.
  In each, replace: (a) `method = "ieppa"` / `"ieppa_soft"` → `"oris"` / `"oris_soft"`; (b) **result assertions** `attr(res, "algorithm") == "ieppa"` / `algorithm_used == "ieppa"` / `expect_equal(..., "ieppa")` → `"oris"` (these FAIL at runtime if missed); (c) `c("ieppa","ieppa_soft", ...)` lists; (d) variable names like `r_ieppa` and comments. **CRITICAL**: the `algorithm`/`algorithm_used` string assertions are the high-risk class — search each file for `algorithm` near `ieppa`.

- [ ] **Step 3: Rename + regenerate fixtures** (THREE fixtures, incl. `task2_ieppa_ref.rds`):
```bash
git mv tests/testthat/fixtures/ieppa_kl_reference_stepstone.rds tests/testthat/fixtures/oris_kl_reference_stepstone.rds
git mv tests/testthat/fixtures/ieppa_pre_alm_ref.rds tests/testthat/fixtures/oris_pre_alm_ref.rds
git mv tests/testthat/task2_ieppa_ref.rds tests/testthat/task2_oris_ref.rds   # NOTE: lives at testthat root, NOT fixtures/; opaque blob, no generator; used by test-simd-math.R
Rscript data-raw/gen_oris_kl_ref.R && Rscript data-raw/gen_oris_pre_alm_ref.R
```
Update the `test_path(...)` references: `test_path("fixtures/ieppa_*.rds")` → `fixtures/oris_*.rds` (`test-calibration-solvers.R`); `test_path("task2_ieppa_ref.rds")` → `test_path("task2_oris_ref.rds")` (`test-simd-math.R`). The `kl`/`pre_alm` fixtures are regenerated by their generators; `task2_oris_ref.rds` has **no generator** — it is moved as-is (it stores numeric SIMD-reference output; enum values are frozen so the blob stays valid).

- [ ] **Step 4: Add dispatch round-trip test** (new file `tests/testthat/test-oris-dispatch.R`):
```r
test_that("oris and oris_soft dispatch and report their names", {
  d <- data.frame(a = factor(c("x","y","x","y")), w = c(1,1,1,1))
  tgt <- list(a = c(x = 0.5, y = 0.5))
  r1 <- harvest(d, targets = tgt, weight = "w", method = "oris")
  r2 <- harvest(d, targets = tgt, weight = "w", method = "oris_soft")
  expect_match(attr(r1, "algorithm"), "oris")
  expect_match(attr(r2, "algorithm"), "oris_soft")
})
```
(Accessor verified: `harvest()` sets `attr(weights, "algorithm")` at `R/harvest.R:749` — NOT `algorithm_used` (that is the raw C list element, not exposed on the returned object). The assertion target: the returned `attr(,"algorithm")` string round-trips to `oris`/`oris_soft`. Adapt the `harvest()` arg names to the actual signature via `?harvest`.)

- [ ] **Step 5: Test gate.**
Run: `Rscript -e 'devtools::test()'`
Expected: all green (no failures referencing missing `ieppa` fixtures/strings).
Run: `cd python && pip install -e . && pytest -q`
Expected: parity green (rtol=1e-6).

- [ ] **Step 6: Commit.**
```bash
git add -A tests/ data-raw/
git commit -m "test(oris): rename ieppa tests/fixtures/generators → oris; add dispatch round-trip test"
```

---

## Task WU4 — Docs + live benchmarks + man

**Files:** `git mv docs/methods/ieppa.md docs/methods/oris.md`; modify `docs/methods/00-overview.md`, `docs/methods/{raking,sinkhorn,newton_kl,chebyshev,greg}.md`, `CLAUDE.md`, `AGENTS.md`, `NEWS.md`, `docs/raking.md`, `docs/ieppa_assessmen.md`, `docs/internal/r_bridge_floor.md`; live `benchmarks/*.{R,py}` + `git mv benchmarks/ieppa_vs_raking_bench.R benchmarks/oris_vs_raking_bench.R`.

- [ ] **Step 1: Rename + retitle the method doc.**
```bash
git mv docs/methods/ieppa.md docs/methods/oris.md
```
In `oris.md`: title → "ORIS — Over-Relaxed Iterative Scaling"; `> Enum:` line → `RK_ALG_ORIS = 1` (+ `RK_ALG_ORIS_SOFT = 8`); `> Source:` → `src/oris.cpp`, …; keep the "Relationship to the source paper" section and add one sentence: the name was changed *from* iEPPA because the code implements the paper's inner C=0 scaling step, not its outer-PPA contribution. Update inline `IEPPA`→`ORIS` (keep "renamed from iEPPA" notes).

- [ ] **Step 2: Overview + sibling docs.** `docs/methods/00-overview.md`: table row, selection rows, enum block, the per-method-documents link `ieppa.md`→`oris.md`, prose `IEPPA`→`ORIS`. In `sinkhorn.md`, `newton_kl.md`, `raking.md`, `chebyshev.md`: update "vs IEPPA" prose `IEPPA`→`ORIS` and retarget Markdown links `ieppa.md`→`oris.md` (incl. the anchor `#ieppa-vs-sinkhorn-dykstra` → `#oris-vs-sinkhorn-dykstra` if the heading changed).

- [ ] **Step 3: Agent/project docs.** `CLAUDE.md` + `AGENTS.md`: update every `ieppa`/`IEPPA` mention (TU-split note → `oris.cpp`/`oris_finalize.cpp`/`oris_trajectory.cpp`/`oris_internal.hpp`; SRAA best-iterate note; the ALM-Newton "do not fix" formula note; CORE_SOURCES note). Keep meaning; only rename.

- [ ] **Step 4: Live root docs.** `docs/raking.md`, `docs/internal/r_bridge_floor.md`: `ieppa`→`oris`/`ORIS`. `docs/ieppa_assessmen.md`: update references + add header "> Superseded by `docs/methods/oris.md`."

- [ ] **Step 5: NEWS.md** — add development entry:
```markdown
# leafblower (development)
* BREAKING: solver `method="ieppa"`/`"ieppa_soft"` renamed to `"oris"`/`"oris_soft"`
  (ORIS — Over-Relaxed Iterative Scaling). No alias; update calls. Algorithm,
  numeric output, and enum values are unchanged.
```

- [ ] **Step 6: Live benchmarks.** `git mv benchmarks/ieppa_vs_raking_bench.R benchmarks/oris_vs_raking_bench.R`; enumerate live bench scripts first: `grep -rIl -i ieppa benchmarks/ | grep -vE 'benchmarks/(2apm|yh0l|results)/'`, then in each update `"ieppa"`→`"oris"` method strings/labels. Known set incl.: `allmethod_bench.{R,py}`, `stepstone_*.{R,py}`, `parity_bench.{R,py}`, `algo_selection_benchmark.R`, `newton_kl_bench.R`, `sor_sweep.R`, and **`benchmarks/plot_helpers.R`** (the output-filename literal `sprintf("ieppa_vs_raking_3d_slice_%s.pdf", …)` → `oris_vs_raking_3d_slice_…`). Do NOT touch `benchmarks/2apm/`, `benchmarks/yh0l/`, `benchmarks/results/` (historical, §8.10).

- [ ] **Step 7: Regenerate man.**
Run: `Rscript -e 'devtools::document()'`
Expected: `man/*.Rd` regenerated; `grep -ri ieppa man/` returns nothing.

- [ ] **Step 8: Commit.**
```bash
git add -A docs/ CLAUDE.md AGENTS.md NEWS.md benchmarks/ man/
git commit -m "docs(oris): rename ieppa→oris in live docs, benchmarks, man; keep history untouched"
```

---

## Task WU5 — Verification gate + repo-map refresh

- [ ] **Step 1: Build + tests + parity.**
```bash
R CMD INSTALL --preclean .
Rscript -e 'devtools::test()'
( cd python && pip install -e . && pytest -q )
```
Expected: build clean; R tests green; Python parity green (rtol=1e-6).

- [ ] **Step 2: Stepstone no-regression** (Definition of Done):
Run the stepstone benchmark per repo convention; compare to pre-rename baseline.
Expected: no regression (rename is behavior-neutral → identical results).

- [ ] **Step 3: Grep-clean gate** (self-contained: refresh the generated graph first, then audit):
```bash
graphify update . 2>/dev/null || true   # regenerate graphify-out from renamed sources before auditing
grep -rIi ieppa . \
  --include='*.cpp' --include='*.hpp' --include='*.h' --include='*.R' --include='*.Rd' \
  --include='*.py' --include='*.toml' --include='*.in' --include='*.md' --include='*.txt' \
| grep -vE '(^|/)(\.git|\.beads|\.wolf)/' \
| grep -vE '(^|/)(docs/superpowers/(plans|specs)|docs/investigations|docs/iEPPA|tasks)/' \
| grep -vE '(^|/)benchmarks/(2apm|yh0l|results)/' \
| grep -vE '(^|/)graphify-out/' \
| grep -vE '_92c4f45\.|references\.bib|\.Rcheck/'
```
Expected: only intentional "renamed from iEPPA" notes in `oris.md`/`00-overview.md`/`CLAUDE.md`/`AGENTS.md`/`docs/ieppa_assessmen.md` and the `chu2022ieppa` citation. Any other hit = a missed live site → fix. (`graphify-out/` is generated and regenerated by Step 5's `graphify update .`; it is excluded from the gate, not hand-edited.)

- [ ] **Step 4: Enum-freeze confirmed** — the WU1 `static_assert` compiled (covered by Step 1).

- [ ] **Step 5: Repo map refresh.**
```bash
graphify update . 2>/dev/null || true
```
Update `.wolf/anatomy.md` for the file renames per OpenWolf protocol.

- [ ] **Step 6: Final commit (or squash WU1–WU5).**
```bash
git add -A
git commit -m "chore(oris): verification + repo-map refresh for ieppa→oris rename"
```

---

## Self-Review

- **Spec coverage:** WU1 ⊇ §3.1/§3.2/§8.1/§8.5; WU2 ⊇ §3.3/§3.4/§8.2/§8.3; WU3 ⊇ §3.5/§8.4; WU4 ⊇ §3.6/§3.7/§8.9/§8.10 live docs; WU5 ⊇ §6/§8.8 DoD. Historical exclusions (§8.6/§8.10) deliberately untouched.
- **Placeholders:** the only adapt-on-the-spot item is WU3 Step 4's `harvest()` result accessor — flagged because the exact algorithm-name accessor must be read from the API; the *assertion intent* is fixed.
- **Type consistency:** `RK_ALG_ORIS`/`RK_ALG_ORIS_SOFT`, `oris_solve`, `ORISResult`, `oris_auto_selected`, method strings `"oris"`/`"oris_soft"` used consistently across WUs.
- **Atomicity:** WU1 will not compile until complete (C++ rename) — build only at WU1 Step 8; branch tip green at WU5.
