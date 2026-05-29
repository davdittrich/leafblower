# Code-Review Remediation Batch — Active Plan
**Date:** 2026-05-29
**Status:** `AWAITING_PLAN_REVIEW_GATE`
**Source:** Whole-repo code review 2026-05-29 (6 parallel reviewers, all findings independently verified by orchestrator).
**Scope:** 7 epics / 25 atomic tickets. All findings VERIFIED against raw source; false positives discarded (see §7).

---

## PLAN HEADER (mandatory)

- **Mechanism:** Surgical, behavior-preserving remediation of verified review findings. Each ticket = one root cause, one file-group, one diff. Build gate after every ticket.
- **Forbidden:**
  - Touching verified-correct solver math — chebyshev Mehrotra Phase-B corrector, ieppa ALM Newton step `X̃(1−λ+μz)/(1+ρ)`, SRAA `select_metric` best-iterate, `res.base.*` access pattern. (Documented in CLAUDE.md + bd memories; re-flagged falsely ≥2× historically.)
  - `git rm` of working files (use `--cached` only). Untracking `.beads/ .wolf/ .claude/` (intentional tracked infra).
  - Dropping `-O3`/`-mavx2` (perf/intrinsics) — relocate, do not delete.
  - `-flto` anywhere. `--no-verify` on commits.
  - `git add -A`/`git add .` in any subagent — explicit per-file staging only.
- **Audit:** Every numerical/parity claim re-verified by orchestrator (independent test rerun), never self-certified by implementer. R gate `R CMD INSTALL --preclean .` + `devtools::test()`; Python gate `pip install -e . --break-system-packages && pytest`. Parity = both-converged precheck THEN `rtol=1e-6`. Stepstone benchmark no-regression for any solver-touching change.

---

## 1. Objective

Clear the 25 verified defects so the package reaches a clean, CRAN/PyPI-submittable, parity-consistent state. Directly unblocks `leafblower-kk1.24.3` (CRAN submission package) — the HYG epic is its hard prerequisite (epic→epic dep unsupported by bd; tracked here in prose).

## 2. Epics & Priority

| Epic | ID | Pri | Tickets | Theme |
| :--- | :--- | :-- | :--- | :--- |
| CR-HYG | `uj6x` | P1 | .1–.4 | Build-artifact + CRAN-flag hygiene |
| CR-PAR | `6uhm` | P1 | .1 | R↔Python AUTO-dispatch parity |
| CR-NABIN | `4ihf` | P1 | .1–.3 | NA-bin diagnostics (post-tu15 consistency) |
| CR-RVAL | `qjj4` | P2 | .1–.5 | R input validation / perf |
| CR-CXX | `5fm8` | P2 | .1–.4 | C++ robustness |
| CR-TEST | `wnn5` | P2 | .1–.5 | Test hardening + parity coverage |
| CR-META | `ccor` | P3 | .1–.3 | Metadata / minor cleanups |

## 3. Execution Order & Rationale

1. **HYG.2, HYG.3 first** (build-system). Correct compiler-flag plumbing before any code change so all later work builds on the fixed Makevars/CMake. Verify `R CMD INSTALL --preclean .` + `pip install -e .` stay green and `-O3 -mavx2` still applied.
2. **HYG.1, HYG.4** (untrack artifacts / gitignore). Independent, trivial, no build impact.
3. **PAR.1** (single comparator fix) — high-value parity, isolated.
4. **NABIN.1 + .2 + .3 as a coordinated unit** — R×2 + Python must land together; parity test gates the trio.
5. **RVAL.1–.5, CXX.1–.4** — P2 correctness/robustness, mutually independent; any order.
6. **CR-TEST**: `TEST.4` (conftest BLAS env) BEFORE `TEST.3` (parity tests) — encoded as bd dep `wnn5.3 → wnn5.4`. `TEST.1/.2/.5` independent.
7. **CR-META** last (P3, cosmetic).

## 4. Per-Epic Mechanism / Forbidden / Audit

### CR-HYG (`uj6x`)
- **Mechanism:** `git rm --cached` the 13 artifacts; relocate `@OPT_FLAGS@`/`@MAVX2_FLAG@` from `PKG_CXXFLAGS` to a non-PKG var / `CXXFLAGS +=`; arch-gate `-mavx2` in CMake.
- **Forbidden:** untracking intentional infra; deleting working artifacts; removing flags.
- **Audit:** generated `src/Makevars` `PKG_CXXFLAGS` line grep-free of `-O`/`-m`; compile log still shows `-O3 -mavx2`; `git ls-files | grep -E '\.(o|so)$'` empty; build green.

### CR-PAR (`6uhm`)
- **Mechanism:** flip r_bridge comparator `>` → `>=` (match c_api:190 `M_cell*10 >= n*9`).
- **Forbidden:** altering the 10/9 ratio or `estimate_M_cell`.
- **Audit:** boundary fixture (`M_cell/n == 0.9` exactly) → identical `algorithm_used` R vs Python.

### CR-NABIN (`4ihf`)
- **Mechanism:** `col_char <- as.character(col); col_char[is.na(col)] <- "NA"` before the level loop in current_miss.R / diagnose_weights.R; `pd.isna(series)` mask for the `'NA'` level in _harvest.py. Consistent all-obs denominator.
- **Forbidden:** changing non-NA bins; diverging R vs Python.
- **Audit:** `add_na_proportion=TRUE` fixture → NA-bin observed share ≈ `na_frac`; R==Python `rtol=1e-6`.

### CR-RVAL (`qjj4`)
- **Mechanism:** target-sum guard before na rescale; `...` unknown-arg warning; OOV warning; strided `data_codes` assignment; 1-arg weight validation.
- **Forbidden:** changing `-1` C-layer sentinel semantics; altering data_codes interleave layout; breaking legitimate `...` forwarding (default to warn, not stop).
- **Audit:** RVAL.4 — assert `data_codes` byte-identical to old loop on a small fixture before claiming speedup; deff unchanged.

### CR-CXX (`5fm8`)
- **Mechanism:** r_bridge length checks (`pre_error` pattern); greenkhorn best-iterate via `select_metric(cfg.metric,…)`; newton_calib `predicted<=0 → rho=0`; c_api `sizeof` + explicit NUL.
- **Forbidden:** touching greenkhorn greedy argmax / step; Armijo acceptance; verified-correct math.
- **Audit:** mismatched-length `.Call` returns error (no crash); newton_kl + greenkhorn fixtures show no convergence regression.

### CR-TEST (`wnn5`)
- **Mechanism:** `type="message"→"output"`; replace `expect_true(TRUE)` with bound/convergence asserts; 5 new Python parity modules; `conftest.py` BLAS env; sinkhorn quality assert.
- **Forbidden:** xfail/skip-masking a real failure; arbitrary (non-derived) thresholds.
- **Audit:** if a de-vacuumed assert now fails, REPORT as a real finding (new ticket) — do not re-mask. Parity tests: both-converged precheck + `rtol=1e-6, atol=0`.

### CR-META (`ccor`)
- **Mechanism:** description text; resolve `pct` alias + drop dup `5L`; lazy pandas import.
- **Forbidden:** changing version/deps; behavior change for `metric='pct'`.
- **Audit:** `import leafblower` succeeds without pandas; `metric='pct'` still maps to l1_weight.

## 5. Alternatives Considered

- **One bundled "cleanup" ticket** vs **25 atomic tickets** → atomic chosen (CLAUDE.md one-ticket-per-task; enables selective revert + independent review; bd memory: bundled commits defeat `git revert`).
- **HYG.2 Makevars: append to `PKG_CXXFLAGS` differently** vs **separate var + `CXXFLAGS +=`** → separate var: CRAN's check greps the literal `PKG_CXXFLAGS` assignment; only relocation (not reformatting) clears it. User-space `CXXFLAGS` is permitted to carry `-O`.
- **NABIN: rename NA bin to `"__NA__"`** (collision-safe) vs **keep `"NA"`** → keep: tu15 deliberately chose `"NA"` (documented decision); changing it is out of scope and breaks the just-landed encoding. Only the diagnostic readers are fixed.
- **CXX.1: validate in R wrapper only** vs **validate in C bridge** → C bridge: r_bridge is reachable via direct `.Call` (documented: R bypasses c_api); defense belongs at the trust boundary.
- **TEST.3: mock-based vs subprocess-Rscript parity** → subprocess-Rscript (matches existing `test_harvest_na_parity.py` harness; real cross-language parity, not a mock).

## 6. Risk Register

| Risk | Mitigation |
| :--- | :--- |
| HYG.2 breaks build (flags lost) | Verify compile log shows `-O3 -mavx2`; bench no-regression before close. |
| NABIN denominator choice shifts non-NA shares | Pick all-obs denominator once, apply to all 3 impls; parity test catches drift. |
| De-vacuumed test (TEST.1/.2) exposes latent solver bug | Treat as discovery → new ticket, not a mask. |
| CXX.2 best-iterate change alters greenkhorn outputs | Regression fixture; convergence-metric semantics documented vs sinkhorn. |
| Subagent fabricates bench/parity numbers | Orchestrator re-runs all numeric claims (bd memory mandate). |

## 7. Discarded False Positives (do NOT re-open)

- `python/CMakeLists.txt` "missing `design_effect.cpp`/`raking.cpp`" — both present (L28/L26). CORE_SOURCES complete.
- `configure:3` "stray dead line" — it is a comment; reviewer read lean-ctx-compressed output.
- "tracked `.pyc`" — `git ls-files | grep .pyc` empty (start-of-session status was stale).
- `test-priority-sweep.R:55` "wrong attr key" — harvest sets `attr(,"iterations")` (L735). Valid.
- `test-calib-linalg.R:153` "type=message" — actually `type="output"`. Valid.
- `test-ieppa-faithful.R:82/209` — :82 has no capture.output; :209 uses `"output"`. Only :153 is real.
- Untracking `.beads/ .wolf/ .claude/` — intentional tracked infra.

## 8. Definition of Done (batch)

- [ ] 25 tickets closed, each with passing build + tests.
- [ ] `R CMD INSTALL --preclean .` clean; `devtools::test()` green.
- [ ] `pytest` green incl. 5 new parity tests + conftest.
- [ ] `git ls-files | grep -E '\.(o|so)$'` empty.
- [ ] Generated `src/Makevars` `PKG_CXXFLAGS` free of `-O`/`-m`.
- [ ] Stepstone benchmark: no regression.
- [ ] HYG epic complete → notify for `leafblower-kk1.24.3` (CRAN) unblock.
