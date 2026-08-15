---
schema_version: 1
open_count: 0
waived_count: 0
fixed_count: 1
total_count: 1
last_updated: 2026-08-15T21:10:56.391Z
---

# Broken Windows Ledger

> Cross-phase defect register. With `workflow.windows_enforce` enabled, `/gsd-ship` blocks while `open_count > 0`.
> Waive with `gsd-tools windows waive <id> "<reason>"` (reason required).
> Mark fixed with `gsd-tools windows fixed <id>`.

| id | phase | kind | file | line | description | status | reason | recorded_at | resolved_at |
|----|-------|------|------|------|-------------|--------|--------|-------------|-------------|
| 1 | 05 | unmet-truth | .planning/phases/05-cran-pypi-release/05-01-SUMMARY.md |  | Task 1's R CMD check --as-cran result is 1 WARNING + 3 NOTEs (not 0 warnings/<=1 NOTE): checkbashisms/tidy/V8 absent locally (sandbox blocks pacman install), plus a local-R-Makeconf compilation-flags NOTE. Expected to close on 05-03/05-04's CI matrix. | fixed |  | 2026-08-15T17:35:12.698Z | 2026-08-15T21:10:56.391Z |

````json
[
  {
    "id": 1,
    "kind": "unmet-truth",
    "phase": "05",
    "file": ".planning/phases/05-cran-pypi-release/05-01-SUMMARY.md",
    "line": null,
    "description": "Task 1's R CMD check --as-cran result is 1 WARNING + 3 NOTEs (not 0 warnings/<=1 NOTE): checkbashisms/tidy/V8 absent locally (sandbox blocks pacman install), plus a local-R-Makeconf compilation-flags NOTE. Expected to close on 05-03/05-04's CI matrix.",
    "status": "fixed",
    "reason": "",
    "recorded_at": "2026-08-15T17:35:12.698Z",
    "resolved_at": "2026-08-15T21:10:56.391Z"
  }
]
````
