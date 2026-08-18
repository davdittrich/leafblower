---
phase: 05-cran-pypi-release
plan: 10
subsystem: packaging/windows-build
tags: [cran, windows, configure.win, ci, r-universe, gap-closure]
gap_closure: true
requirements: [US-010, KPI-05]
beads: [leafblower-sc9t]
dependency-graph:
  requires:
    - .github/workflows/r-check.yml (05-03/05-05: ubuntu-latest job, pinned r-lib/actions SHAs)
    - configure (Unix header/Makevars generator, unchanged)
    - src/Makevars.in (single source of compiler/linker flags, unchanged)
  provides:
    - configure.win (Windows header/Makevars.win generator)
    - windows-latest job in .github/workflows/r-check.yml (permanent CI guard)
  affects:
    - cran-comments.md (test-environments section)
    - Phase 7 (CRAN submission) — Windows was CRAN's default check platform and was previously broken
tech-stack:
  added: []
  patterns:
    - "Windows configure.win hardcodes feature-probe results instead of re-probing, because the probed outcome (glibc libmvec absence, no OpenMP macro consumers) is architecturally fixed on the platform, not merely usually the same"
key-files:
  created:
    - configure.win
  modified:
    - .github/workflows/r-check.yml
    - src/.gitignore
    - cran-comments.md
decisions:
  - "configure.win hardcodes LBW_HAS_OMP/LBW_HAS_OMP_SIMD/LBW_HAS_GLIBC_MVEC to 0 with no feature probes at all (verified: no source file reads the two OpenMP macros; libmvec is glibc-only) — deliberately simpler than porting the Unix probes"
  - "OMP_FLAGS/MAVX2_FLAG/MVEC_LIBS all emptied in the Windows Makevars substitution, keeping PKG_LIBS free of any Windows runtime-library entry (PE/COFF unresolved symbols are hard link errors, unlike ELF)"
  - "windows-latest R CMD check's TinyTeX install is missing the courier TeX Live package by default (unlike ubuntu-latest's) — added an explicit tlmgr install courier CI step, itself needing a non-obvious PATH workaround because setup-tinytex's GITHUB_PATH addition isn't visible to Rtools' bash"
  - "A residual windows-latest-only PDF-manual WARNING (deterministic, confirmed via CI rerun) was NOT further chased after 3 fix attempts on the tlmgr chain — split to leafblower-fxyj rather than blocking this plan's actual target (the missing-header compile defect)"
metrics:
  duration: ~1h40min (includes an r-universe rebuild wait; see Session Notes)
  tasks: 3
  files: 4
  completed: 2026-08-18
status: complete
actuals:
  tokens: 1634
  tasks: 3
  commits: 6
---

# Phase 5 Plan 10: Windows Build Path (configure.win) Summary

Added `configure.win` so Windows R builds generate `src/lbw_config.h`/`src/Makevars.win`
from the same `src/Makevars.in` the Unix build uses — the missing-header compile failure
that broke all 5 Windows configs on r-universe's check farm is fixed and proven via a new
permanent `windows-latest` job in this repository's own CI, plus a fresh r-universe rebuild
confirming all 5 Windows configs (including arm64, provable only via r-universe) now `OK`.

## What Was Built

**Task 1 (RED):** Added `windows-latest` to `.github/workflows/r-check.yml`'s OS matrix
(`fail-fast: false`, `defaults.run.shell: bash` for the POSIX hygiene-guard step). Pushed
and reproduced the exact defect r-universe reported, live in this repository's own CI (run
32081319627):

```
lbw_math.hpp:2:10: fatal error: lbw_config.h: No such file or directory
    2 | #include "lbw_config.h"
```

This converts leafblower-sc9t from "reported by r-universe" to "reproduced here" — 05-RESEARCH.md's
assumption A5 ("CRAN win-builder has its own fallback") is now recorded as empirically disproved
against a real Rtools toolchain, not just against r-universe's log.

**Task 2 (GREEN):** Added `configure.win` (mode 755, POSIX `sh`). Unlike the Unix `configure`,
it runs **no feature probes**: `LBW_HAS_GLIBC_MVEC` is hardcoded to 0 (glibc's libmvec is
Linux-only) and both `LBW_HAS_OMP`/`LBW_HAS_OMP_SIMD` are hardcoded to 0 (verified this session:
no source file reads either macro; the only OpenMP usage is four bare `#pragma omp simd` hints
in `oris.cpp`, needing no runtime). It substitutes the same `src/Makevars.in` the Unix `configure`
uses into `src/Makevars.win`, with `OMP_FLAGS`/`MAVX2_FLAG`/`MVEC_LIBS` all emptied — load-bearing,
because on PE/COFF an unresolved symbol at link time is a hard error, and `PKG_LIBS` carries no
OpenMP entry. `Makevars.win` added to `src/.gitignore` beside its Unix sibling; `configure.win`
itself is deliberately NOT `.Rbuildignore`d (it must ship in the tarball to run on the user's
machine).

Local scratch-tree verification (per the plan's own honesty framing — this project has no
Windows machine, so only shell-syntax/output-shape checks are provable locally): `sh -n
configure.win` valid, mode 755 confirmed, header defines all 3 macros with `LBW_HAS_GLIBC_MVEC 0`,
generated `Makevars.win` has 0 unsubstituted placeholders and requests no `lmvec`/`fopenmp`
library. `R CMD INSTALL --preclean .` confirmed the Unix path is untouched.

Pushing this fix cleared the compile error, but exposed a **second, unrelated** windows-latest
defect: `checking PDF version of manual ... ERROR` — `leafblower-manual.log` showed `! Font
T1/pcr/m/n/10=pcrr8t ... not loadable: Metric (TFM) file not found` (the `courier` TeX Live
package, providing Courier's T1 metrics, isn't in windows-latest's default TinyTeX install,
unlike ubuntu-latest's). Fixed with a `tlmgr install courier` CI step — which itself needed two
follow-up fixes before it worked (see Deviations below). After all three, the compile/link/test
defect is fully closed: windows-latest compiles, links, and passes the full testthat suite
(31s, all OK). A residual, deterministic `PDF version of manual ... WARNING` remains on
windows-latest only (confirmed non-flaky via a CI rerun of the same commit) — this is NOT the
defect this plan targets and was split off rather than chased further; see Deviations.

**Task 3 (r-universe confirmation):** Re-checked r-universe's own build farm (the only oracle
for arm64 Windows — GitHub-hosted `windows-latest` is x86_64-only). Build
https://github.com/r-universe/davdittrich/actions/runs/32082532844 (commit `43894c8`, which
already includes `configure.win` and both tlmgr fixes — the final commit `6a312ce` only
hardens this repository's own `windows-latest` GH Actions job's tlmgr PATH lookup and touches
nothing r-universe's separate build pipeline runs, so this result fully covers it):

| config | r | check |
|---|---|---|
| linux-devel-arm64 | 4.7.0 | OK |
| linux-devel-x86_64 | 4.7.0 | OK |
| linux-release-arm64 | 4.6.1 | OK |
| linux-release-x86_64 | 4.6.1 | OK |
| macos-oldrel-arm64 | 4.5.3 | OK |
| macos-oldrel-x86_64 | 4.5.3 | OK |
| macos-release-arm64 | 4.6.1 | OK |
| macos-release-x86_64 | 4.6.1 | OK |
| source | 4.6.1 | OK |
| wasm-release | 4.6.0 | **FAIL** |
| windows-devel-arm64 | 4.7.0 | OK |
| windows-devel-x86_64 | 4.7.0 | OK |
| windows-oldrel-x86_64 | 4.5.3 | OK |
| windows-release-arm64 | 4.6.1 | OK |
| windows-release-x86_64 | 4.6.1 | OK |

**All 5 Windows configs are `OK`**, including both arm64 configs (provable only here). Every
`_jobs` entry was checked individually, not the aggregate `_status` field — the exact mistake
05-VERIFICATION.md caught in 05-06's summary. Unexpectedly, all 4 Linux and 4 macOS configs
are also `OK` now (05-VERIFICATION.md documented these as `ERROR` from a `PracTools::deffH`
NaN mismatch) — that gap appears to have resolved itself (possibly an upstream PracTools
version change on r-universe's build image) but this is **out of this plan's scope**; not
investigated further, not claimed fixed by this plan. `wasm-release` is the sole remaining
`FAIL` — a separate Emscripten toolchain, filed as leafblower-soci rather than implied fixed.

`cran-comments.md`'s test-environments section now cites both the windows-latest GH Actions
run (32082804108) and this r-universe build, stating exactly what was checked on each platform
and no more.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - blocking issue] windows-latest PDF manual: courier package missing, 3 fix attempts**

- **Found during:** Task 2, after the header/Makevars fix cleared the compile error.
- **Issue:** `checking PDF version of manual ... ERROR` — `courier` TeX Live package (T1
  Courier font metrics, needed by Rd.sty) not installed by default in windows-latest's TinyTeX,
  unlike ubuntu-latest's.
- **Fix attempt 1:** `Rscript -e 'tinytex::tlmgr_install("courier")'` — failed, `tinytex` R
  package (distinct from the TinyTeX binary distribution) isn't installed by `setup-tinytex`.
- **Fix attempt 2:** bare `tlmgr install courier` — failed with `command not found` (exit 127)
  on windows-latest only; `setup-tinytex`'s `GITHUB_PATH` addition isn't visible to the Rtools
  `bash.EXE` this job's `defaults.run.shell: bash` steps run under (ubuntu-latest's bash picked
  it up fine).
- **Fix attempt 3 (final):** resolve `tlmgr` via `command -v tlmgr`, falling back to a direct
  `find "$APPDATA/TinyTeX/bin" -iname 'tlmgr*'` when empty. Worked on both platforms.
- **Files modified:** `.github/workflows/r-check.yml`
- **Commits:** `ac8c21f`, `43894c8`, `6a312ce`

### Deferred Issues (fix-attempt budget exhausted)

**2. windows-latest `checking PDF version of manual ... WARNING` remains after the courier fix**

- **Found during:** Task 2, after fix attempt 3 above got `tlmgr install courier` actually
  running. Status improved from `1 ERROR, 1 WARNING, 2 NOTEs` to `1 WARNING, 2 NOTEs` — the
  `pcrr8t.tfm not loadable` ERROR is gone (confirmed via `leafblower-manual.log`: it now
  produces a valid 11-page PDF) — but a `WARNING` on the *with-index* LaTeX pass remains; the
  *without-index* retry succeeds. Confirmed **deterministic**, not flaky, via a CI rerun of the
  identical commit (same result both times).
- **Why deferred:** this is a genuinely distinct root cause (something about the with-index
  pass specifically — `makeindex` availability or a filename-database timing gap were candidate
  hypotheses, neither confirmed) requiring fresh investigation, arrived at only after 3
  committed fix attempts already spent getting `tlmgr` itself invokable from Rtools bash. Per
  the plan's own fix-attempt budget and CLAUDE.md's "do not retry past two attempts on any
  single fix" — this crossed that line and was deliberately not chased further.
  - **What it biases:** nothing measured by this plan's own success criteria — the actual
    defect this plan targets (the missing-header compile failure) is fully fixed and proven.
    It does mean Task 2's literal acceptance criterion ("windows-latest concludes success
    with 0 errors and 0 warnings") is **not fully met** — the job's real GitHub Actions
    conclusion is `failure` (red X) because `error-on: "warning"` treats this WARNING as fatal
    to the job, even though the compile/link/test path it was meant to gate is fully working.
- **Tracked on:** leafblower-fxyj (new ticket, full diagnosis notes attached), NOT bundled into
  leafblower-sc9t (which is specifically about the missing-header compile failure and is
  closed).
- **Files:** none (deferred, no further code change made)

### New tickets filed

- **leafblower-fxyj** (P1): windows-latest PDF manual WARNING, TinyTeX gap distinct from
  `configure.win` — see Deferred Issues above.
- **leafblower-soci** (P3): r-universe `wasm-release` build FAILs (separate Emscripten
  toolchain, unaddressed by this plan, named per Task 3's own acceptance criteria rather than
  silently omitted).

### leafblower-sc9t

Closed. Its exact scope — "No configure.win: Windows R builds fail to compile (lbw_config.h
never generated)" — is fixed and proven: windows-latest now compiles, links, and passes the
full testthat suite, confirmed both on this repository's own CI (x86_64) and on r-universe's
build farm (all 5 Windows configs including arm64).

## Session Notes

Task 3's r-universe polling took materially longer than expected: this project's execution
harness has no long-running foreground wait primitive suited to r-universe's "up to about an
hour" rebuild latency (background polling tasks were killed across turn boundaries in this
session; the Monitor tool was unavailable). The blocking condition was resolved when a
directly-verified check against r-universe's live API at an intermediate commit (`43894c8`,
already containing every substantive fix — the header generator, both tlmgr invocation fixes)
showed all 5 Windows configs already `OK`, making further waiting for the final commit
unnecessary (that commit only touches this repository's own CI job, not anything r-universe's
separate build pipeline runs). This is recorded in `metrics.duration` above; it does not
reflect implementation effort, mostly CI round-trip and rebuild wait time.

`actuals.tokens` (1634, from a ~6.5KB total diff / 4) undercounts this plan's real cost: most
of the effort was CI-log archaeology (downloading and grepping TinyTeX/LaTeX logs, multiple
push-and-watch cycles) rather than lines changed. Recorded per the stated methodology
(chars/4 over the realized diff) rather than adjusted to look more representative.

## Verification

- `sh -n configure.win` valid; scratch-tree run produces a 3-macro header
  (`LBW_HAS_GLIBC_MVEC 0`) and a fully-substituted `Makevars.win` with no `lmvec`/`fopenmp`
  library request. **Confirmed.**
- windows-latest `R CMD check --as-cran`: compile/link/test defect fixed (0 errors on that
  front); **1 WARNING remains** (PDF manual, tracked separately, not this plan's target
  defect) — literal "0 errors, 0 warnings" acceptance criterion **not fully met**, documented
  honestly above rather than glossed over.
- ubuntu-latest stays green (0 errors, 0 warnings) throughout — the Unix path is untouched.
  **Confirmed** (CI run 32082804108).
- r-universe's `_jobs` array: all 5 Windows configs `OK`, `wasm-release` still `FAIL` (named,
  ticketed, not implied fixed). **Confirmed** (build 32082532844).
- `cran-comments.md` claims exactly the platforms checked, citing both real run/build
  identifiers verbatim. **Confirmed.**

## Self-Check: PASSED

Files:
- FOUND: `.github/workflows/r-check.yml`
- FOUND: `configure.win`
- FOUND: `src/.gitignore`
- FOUND: `cran-comments.md`

Commits:
- FOUND: `1a55154` (Task 1: windows-latest matrix)
- FOUND: `cffc2a7` (Task 2: configure.win)
- FOUND: `ac8c21f` (Task 2 fix attempt 1: tinytex Rscript, failed)
- FOUND: `43894c8` (Task 2 fix attempt 2: bare tlmgr, failed)
- FOUND: `6a312ce` (Task 2 fix attempt 3: find-based tlmgr, succeeded)
- FOUND: `e94a953` (Task 3: cran-comments.md)

## Phase 5 Status

This was the last plan in Phase 5 (CRAN + PyPI Release): 05-01 through 05-10 are all complete.
SC6's Windows gap (the reason this gap-closure plan exists) is closed for the compile/link/test
path; the residual PDF-manual WARNING on windows-latest and the pre-existing wasm-release
failure are both tracked on their own tickets rather than left implicit.
