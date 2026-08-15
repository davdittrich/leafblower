# Phase 4: Truthful Surface - Discussion Log

> **Audit trail only.** Do not use as input to planning, research, or execution agents.
> Decisions are captured in CONTEXT.md — this log preserves the alternatives considered.

**Date:** 2026-08-15
**Phase:** 4-Truthful Surface
**Areas discussed:** docs/raking.md §8.2/§12 fix, doc quality flag, grake/lbfgsb/cp audit scope, weights= guard behavior

---

## docs/raking.md §8.2/§12 fix

| Option | Description | Selected |
|--------|-------------|----------|
| Rewrite under paper's name | Relabel as Chu-Liang-Toh-Yang 2022 iEPPA (arXiv:2011.14312), keep content | |
| Delete the passage | Remove §8.2/§12 entirely rather than relabel | ✓ |

**User's choice:** Delete the passage.
**Notes:** Simpler than maintaining a second description of a paper `docs/methods/oris.md`
already covers correctly.

---

## Doc quality (rest of docs/raking.md)

| Option | Description | Selected |
|--------|-------------|----------|
| File a follow-up ticket | Note as deferred idea / new ticket — out of Phase 4 scope | ✓ |
| Ignore for now | Not this phase's problem | |

**User's choice:** File a follow-up ticket.
**Notes:** §1-7, 9-11 have zero citations (vs. 59 in `docs/methods/oris.md`) and include
unrelated content (ROAM medical-imaging architecture, generic Rust-vs-Python benchmark);
§12 calls L-BFGS-B "definitive state-of-the-art" despite removal from the package. Flagged,
not acted on.

---

## grake/lbfgsb/cp audit scope

| Option | Description | Selected |
|--------|-------------|----------|
| Shipped docs only | README, NEWS.md, man/, docs/methods/, docs/raking.md, docs/performance.md, vignettes | ✓ |
| Include docs/superpowers/specs/ | Audit internal design docs too | |

**User's choice:** Shipped docs only.
**Notes:** Internal design specs are historical planning records expected to reference
retired/future method names.

---

## weights= guard behavior

| Option | Description | Selected |
|--------|-------------|----------|
| Hard error (stop()) | Blocks the call, forces the fix | ✓ |
| Specific warning, non-fatal | Keep non-blocking but name design_weights= specifically | |

**User's choice:** Hard error (stop()).
**Notes:** A swapped `weights=` silently returning a plausible-but-wrong unweighted result is
exactly the footgun this phase exists to close.

---

## Claude's Discretion

- Exact `stop()` message wording for `weights=` (beyond naming `design_weights=` explicitly).
- How to document the slot-7 hole in `rk_algorithm_t` once its history is researched.

## Deferred Ideas

- docs/raking.md broader quality review (§1-7, 9-11) — file as a new ticket, not Phase 4.
