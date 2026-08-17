---
phase: 05-cran-pypi-release
plan: 06
subsystem: distribution/r-universe
tags: [r-universe, cran-adjacent, distribution, SC6]
dependency-graph:
  requires: ["05-05"]
  provides: ["SC6: leafblower live and installable from davdittrich.r-universe.dev"]
  affects: []
tech-stack:
  added: []
  patterns: ["r-universe registration via packages.json append (D-07)", "bounded background poll for async external build (D-08)"]
key-files:
  created: []
  modified: []
external-artifacts:
  - "davdittrich/davdittrich.r-universe.dev :: packages.json — commit 8bfdbf0 on main, one entry appended"
decisions:
  - "Registered via the existing r-universe registry repo (no new repo, no App reinstall) per D-07; appended leafblower entry leaving catboostr/robscale byte-identical"
  - "Waited on r-universe's async build via a bounded (90 min budget) background poll rather than asserting immediately post-push, per D-08's explicit anti-flake warning"
metrics:
  duration: "~55min (dominated by the r-universe build wait, ~44min from push to success)"
  completed: 2026-08-17
status: complete
actuals:
  tokens: 30
  tasks: 3
  commits: 1
---

# Phase 5 Plan 06: r-universe registration (SC6) Summary

leafblower is now registered and building green on `davdittrich.r-universe.dev`, with a
clean-R-session `install.packages()` proof executed against the live feed — closing SC6.
This plan modified no file inside the `leafblower` repository itself; all durable output is
external (the `davdittrich.r-universe.dev` registry repo) plus this SUMMARY.

## What Was Built

**Task 1 — Rd2pdf pre-flight (Pitfall A6 check).** Ran the exact `R CMD Rd2pdf` step
r-universe's build pipeline runs (which this repo's own `r-check.yml` does not exercise)
locally under a 300s hard timeout: `rm -f /tmp/lbw-manual.pdf && timeout 300 R CMD Rd2pdf
--no-preview --force --output=/tmp/lbw-manual.pdf .`. It completed in well under the timeout —
exit 0, 13-page PDF, 154526 bytes. No Rd file needed an escape fix; the man/harvest.Rd
apostrophes at lines 122/302 sit in ordinary prose outside any markup span (confirmed
non-defects per the plan's explicit instruction not to preemptively rewrite them). Nothing in
`man/*.Rd` was touched.

**Task 2 — Registry entry + bounded build wait.** Cloned
`davdittrich/davdittrich.r-universe.dev` to a scratch directory outside this repository,
appended `{"package": "leafblower", "url": "https://github.com/davdittrich/leafblower"}` to
`packages.json` (D-07), leaving the `catboostr` and `robscale` entries byte-identical, and
pushed to `main` (commit `8bfdbf0`). Ran the D-08 bounded wait as a single self-contained
background command polling `https://davdittrich.r-universe.dev/api/packages/leafblower` every
60s, budgeted to 90 minutes, and blocked on it (rather than foreground-sleeping or
hand-stepping the poll) until it resolved. The build reached success in ~44 minutes:

```
"_buildurl": "https://github.com/r-universe/davdittrich/actions/runs/32069801748"
"_status": "success"
"_registered": true
```

`packages.json` on `main` verified to contain exactly 3 entries, with `catboostr` and
`robscale` byte-identical to their pre-registration values and the new `leafblower` entry
exactly `{"package": "leafblower", "url": "https://github.com/davdittrich/leafblower"}`.

**Task 3 — Clean-session install proof.** From `/tmp` (outside the repo), with
`.libPaths()` set to a fresh, empty `/tmp/lbw-runiv-lib` (asserted to contain no path
matching both `leafblower` and `projects`, i.e. structurally unreachable from this working
checkout) and single-thread BLAS env vars exported first:

```r
.libPaths("/tmp/lbw-runiv-lib")
install.packages("leafblower", repos="https://davdittrich.r-universe.dev", lib="/tmp/lbw-runiv-lib")
library(leafblower, lib.loc="/tmp/lbw-runiv-lib")
stopifnot(is.function(leafblower::harvest))
```

r-universe served a **source tarball** (`leafblower_0.1.0.tar.gz`, 646 KB) — no prebuilt
binary for this platform/R combination (R 4.6.1, x86_64-pc-linux-gnu) is in the feed, so the
install compiled the C++17 core locally (a user on this platform needs a working
C++ toolchain; macOS/Windows binaries were not checked by this task). Compilation succeeded,
`library(leafblower)` loaded from the isolated library, and `leafblower::harvest` resolved as
a function — proving the shared object actually loaded, not merely that the package
unpacked. SC6's literal wording ("clean R session, no local source checkout") is proven with
executed evidence, not a weaker paraphrase.

## Deviations from Plan

None — plan executed exactly as written. One clarifying note: this plan's `files_modified: []`
meant no per-task commit inside this repository was expected or made for Tasks 1-3 (the only
durable change is the external registry repo's commit `8bfdbf0`); the final `docs(05-06):`
commit below is the only commit landing in this repository, matching the plan's own frontmatter.

## Auth Gates

None — `gh auth status` confirmed an existing authenticated session with push access to
`davdittrich/davdittrich.r-universe.dev` (`repo` scope), so no interactive auth step was
needed.

## Known Stubs

None.

## Threat Flags

None beyond what the plan's own threat model already covers (T-05-09/T-05-10/T-05-11/T-05-12,
all addressed by this plan's Task 2 acceptance checks and Task 1's pre-flight).

## Self-Check: PASSED

- `curl -fsS https://davdittrich.r-universe.dev/api/packages/leafblower` → `_status: "success"`,
  `_registered: true`, `_buildurl` present and resolves to a real r-universe-side Actions run
  (`https://github.com/r-universe/davdittrich/actions/runs/32069801748`) — FOUND, live-verified
  at write time.
- `davdittrich/davdittrich.r-universe.dev` commit `8bfdbf0` — FOUND (`git log` on the cloned
  scratch copy shows it as HEAD of `main`, pushed successfully).
- `/tmp/lbw-runiv-lib/leafblower` — FOUND (installed during Task 3's run; a throwaway temp
  library, safe to discard).

## SC6 Status

**Met.** leafblower is registered on `davdittrich.r-universe.dev`, its r-universe-side build
reports `_status: "success"`, and `install.packages("leafblower", repos =
"https://davdittrich.r-universe.dev")` is proven, with executed evidence, to succeed in a
clean R session with no local source checkout on the library path. leafblower-bl7g (and its
three children bl7g.1/.2/.3) closed.

What this plan does not prove (explicitly out of scope, matching the plan's own
`<verification>` section): that r-universe will keep building green on future commits (that's
r-universe's own periodic-rebuild mechanism), and nothing about PyPI (SC7 — 05-07).
