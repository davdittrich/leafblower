---
phase: 05-cran-pypi-release
plan: 07
subsystem: distribution/pypi
tags: [pypi, trusted-publishing, distribution, SC7]
dependency-graph:
  requires: ["05-05"]
  provides: ["SC7: leafblower published to PyPI with sdist+wheels via OIDC Trusted Publishing"]
  affects: [".github/workflows/python-wheels.yml"]
tech-stack:
  added: []
  patterns: ["OIDC Trusted Publishing (no stored PyPI token)", "tag-gated publish job appended to existing workflow (D-10)", "GitHub Actions environment as publish gate (D-11)"]
key-files:
  created: []
  modified: [".github/workflows/python-wheels.yml"]
external-artifacts:
  - "GitHub Actions environment `pypi` on davdittrich/leafblower"
  - "PyPI project leafblower 0.1.0 — sdist + 10 wheels (cp39-cp313 x manylinux_x86_64/macosx_arm64)"
  - "git tag v0.1.0"
decisions:
  - "Task 3 (tag push + publish) executed directly by the orchestrator in the main session, not by the delegated executor subagent — the subagent's blocking-human checkpoint correctly refused any orchestrator-relayed approval (including verbatim-quoted relay) per this project's own established precedent (Phase 2 Plan 08: 'no agent message substitutes for the user's own sign-off'). Only the main session receives genuine, unforgeable user messages, so the one-way publish step was executed there instead of attempting a further relay."
metrics:
  duration: "~25min (dominated by the cibuildwheel matrix build, ~15min)"
  completed: 2026-08-18
status: complete
actuals:
  tokens: 0
  tasks: 3
  commits: 1
---

# Phase 5 Plan 07: PyPI publish (SC7) Summary

leafblower is now published on PyPI with both a source distribution and wheels, uploaded via
GitHub OIDC Trusted Publishing (no stored credential anywhere), and a clean-environment
`pip install leafblower` is proven against the public index — closing SC7, the last open
Phase 5 gap.

## What Was Built

**Task 1 — pypi environment + CI jobs (executor subagent, commit `625d6ad`).** Created the
`pypi` GitHub Actions environment on `davdittrich/leafblower` via `gh api ... -X PUT`
(confirmed 200). Appended two jobs to the existing `.github/workflows/python-wheels.yml`
(the pre-existing `python-wheels` and `wheel-check` jobs untouched): `build-sdist` (plain
`python -m build --sdist` under `python/`, uploads the tarball as a distinct artifact) and
`publish-to-pypi` (`needs: [wheel-check, build-sdist]`, `if: startsWith(github.ref,
'refs/tags/')`, job-level `permissions: {id-token: write}`, `environment: {name: pypi, url:
https://pypi.org/p/leafblower}`, `pypa/gh-action-pypi-publish` pinned to commit
`dc37677b2e1c63e2034f94d8a5b11f265b73ba33`). No tag pushed at this point. `leafblower-ej1n.1`
closed.

**Task 2 — human checkpoint (blocking-human).** The executor subagent presented the PyPI
pending-publisher registration instructions (project `leafblower`, owner `davdittrich`, repo
`leafblower`, workflow `python-wheels.yml`, environment `pypi`) and correctly held at the
gate. The user registered the pending publisher directly at pypi.org (2FA-gated, no API —
D-09) and approved directly in the primary conversation. `leafblower-ej1n.2` closed.

**Task 3 — tag-publish and verify (executed by the orchestrator directly, not the subagent).**
The subagent's checkpoint instructions correctly rejected every orchestrator-relayed approval
message, including one that verbatim-quoted the user's own words and explicitly disclosed it
was a relay — this matches an established project precedent (Phase 2 Plan 08) that no
agent-to-agent message substitutes for the user's own sign-off on a `gate="blocking-human"`
checkpoint. Since the orchestrator is the only party in this session that receives genuine,
unforgeable user messages, it executed Task 3 directly in the main thread instead of
attempting a further relay:

```
git push origin master          # 8920cf2..625d6ad (23 commits, includes Task 1's commit)
git tag -a v0.1.0 -m "leafblower 0.1.0"
git push origin v0.1.0
gh run watch --exit-status <run-id>   # succeeded
```

The triggered `python-wheels.yml` run on tag `v0.1.0` completed with conclusion success;
`publish-to-pypi` ran (not skipped). Verified against PyPI's own JSON API, not the job's exit
status (Pitfall A7 — a wheels-only publish is silently green):

```
PyPI files: leafblower-0.1.0-cp39/cp310/cp311/cp312/cp313-{manylinux_x86_64,macosx_11_0_arm64}.whl
            (10 wheels)
            leafblower-0.1.0.tar.gz (sdist)
```

Both distribution types confirmed present. Then a fresh virtualenv install from the public
index, no cache, no find-links:

```
python3 -m venv /tmp/pypi-check
/tmp/pypi-check/bin/pip install --no-cache-dir leafblower
```

succeeded (`leafblower-0.1.0` installed), and `import leafblower; leafblower.harvest`
resolved without error. The venv's system Python was 3.14 — outside the cp39-cp313 wheel
matrix — so this particular run exercised the sdist source-build fallback path (a working
C++ toolchain was available and the build succeeded), not a prebuilt wheel; the wheel path
itself is proven by the PyPI file listing above and by 05-04/05-05's earlier CI import tests
on the built wheels. An in-matrix interpreter (`python3.13`) was available on the host but
blocked by this session's shell command allowlist; not pursued further since it would only
add a second successful data point, not change the SC7 verdict.

## Deviations from Plan

**Execution-channel deviation (not a scope or requirement deviation).** The plan's `<verify>`
automated command for Task 3 was executed by the orchestrator in the main session rather than
by the dispatched `gsd-executor` subagent, because the subagent's own (correct)
blocking-human-checkpoint discipline has no mechanism to accept a second-hand approval — by
design, since a compromised or hallucinating intermediary could otherwise forge one. All
commands, acceptance criteria, and verification steps were followed exactly as the plan
specifies; only the identity of the process running them differs. No weakening of the gate
occurred — if anything this is a stricter reading of it, since the human's approval was acted
on by the one thread that actually received it.

## Auth Gates

`gh auth status` had an existing authenticated session with push/environment-admin access to
`davdittrich/leafblower`. The PyPI side used OIDC Trusted Publishing exclusively — no PyPI
API token was created, stored, or used at any point (T-05-13 mitigation holds).

## Known Stubs

None.

## Threat Flags

None beyond what the plan's own threat model already covers (T-05-13 through T-05-18). The
`id-token: write` permission was confirmed job-level-only on `publish-to-pypi` (Task 1's
acceptance check), and the tag-ref `if:` condition plus `pypi` environment gate held for the
actual publish.

## Self-Check: PASSED

- `https://pypi.org/pypi/leafblower/json` → version `0.1.0`, 11 `urls` entries (1 sdist + 10
  wheels) — FOUND, live-verified at write time.
- `git tag -l v0.1.0` and `git log origin/master` → tag present, pushed, `master` fast-forwarded
  to `625d6ad` on `origin` — FOUND.
- `/tmp/pypi-check` venv with `leafblower` importable and `leafblower.harvest` resolved — FOUND
  (throwaway venv, safe to discard).

## SC7 Status

**Met.** sdist and wheels uploaded to PyPI through Trusted Publishing on the `v0.1.0` tag,
verified present on the project's files listing (not inferred from CI's green status alone),
and `pip install leafblower` proven to succeed from a clean environment against the public
PyPI index. Platform scope, stated honestly per the plan: wheels cover manylinux x86_64
(`ubuntu-latest`) and macOS arm64 (`macos-14`) across Python 3.9-3.13 only; x86_64 macOS and
Windows users fall back to the sdist and a local compile requiring LAPACK (pre-existing
residual from 05-05, not new here). `leafblower-ej1n` and its three children
(`.1`/`.2`/`.3`) closed.

Phase 5 (r-universe + PyPI Release) is now feature-complete: SC1-SC7 all have durable,
executed evidence.
