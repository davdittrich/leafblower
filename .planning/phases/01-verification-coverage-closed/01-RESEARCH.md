# Phase 1: Verification Coverage Closed - Research

**Researched:** 2026-08-15
**Domain:** Test-infrastructure closure for an R/C++/Python statistical calibration package —
parity-test coverage, property-based bound testing, testthat edition migration.
**Confidence:** HIGH

## Summary

This phase adds no new library, framework, or external dependency. Every success criterion
is satisfied by extending or relocating existing R (`testthat`) and Python (`pytest`) test
files that already establish the exact patterns to copy. The work is entirely mechanical
extension of proven templates — `_assert_parity()` in `test_solver_parity.py`, the
parametrize list in `test_parity_weights.py`, the `expect_true(all(w <= max))` pattern in
`test-clamp-contract.R` — plus one `DESCRIPTION` field flip gated behind a fix-first
migration protocol.

Two findings materially change what "correct" means for this phase and must reach the
planner as load-bearing, not incidental:

1. **The KPI-02 property test is only a valid, non-vacuous assertion under
   `bounds_mode="unit"`.** The package's default (`bounds_mode="cell"`) *deliberately does
   not* enforce per-observation bounds on the returned weight vector — it only counts
   violations — and `R/harvest.R` emits a warning saying so verbatim: `"cell-mode bounds: %d
   weights fell outside [%.3f, %.3f] due to skewed base weights within cells. Consider
   bounds_mode='unit' for strict per-observation bounds."` (`R/harvest.R:833-835`,
   `[VERIFIED: R/harvest.R:833-835]`). A property test written against the default mode
   would be flaky by design, not a bug detector. It must explicitly pass
   `bounds_mode = "unit"`.
2. **"Eight shipped solvers, `oris_soft` included" undercounts the enum by one.** The live
   `rk_algorithm_t` enum has nine non-`AUTO` solver values (`RK_ALG_ORIS`, `RAKING`,
   `SINKHORN`, `CHEBYSHEV`, `GREG`, `ORIS_SOFT`, `GREENKHORN`, `LOGIT`, `NEWTON_KL` —
   `[VERIFIED: src/leafblower.h:42-52]`, quoted below). SC1's literal target — every one of
   these compared R-vs-Python in the parametrized weight-parity matrix — is nine
   parametrize entries, not eight. The gap between the current 6-entry parametrize list and
   the full set is three methods: `chebyshev`, `greg`, `oris_soft` — not one.

**Primary recommendation:** Sequence exactly as CONTEXT.md's D-09 mandates —
`leafblower-x7n8` (relocate `test_parity_weights.py` so the gate collects it) lands first
and alone; every other change (matrix extension, new convergence-rule tests, tolerance
fix, property test, edition flip) lands only after that relocation is proven green, because
adding coverage to a file the gate does not run is a no-op.

## Architectural Responsibility Map

| Capability | Primary Tier | Secondary Tier | Rationale |
|------------|-------------|----------------|-----------|
| Weight-vector R↔Python parity (SC1) | Test / Verification layer | C++ core (`lbw::finalize_weights_buf`, `src/calib_dispatch.hpp`) | Test-layer change only; the shared core being compared is untouched this phase |
| Convergence-rule parity for raking/sinkhorn (SC2) | Test / Verification layer | — | New pytest cases in `python/leafblower/test_solver_parity.py`, no production code |
| logit tolerance justification (SC3) | Test / Verification layer | C++ core (`src/logit_calib.cpp` Newton+Armijo) investigated but not modified | Diagnosis touches the core only to *read* it; the fix (if any) is a documented tolerance or a filed ticket, never a solver edit this phase |
| Bound property test (SC4) | Test / Verification layer | C++ core (`lbw::finalize_weights_buf` unit-mode water-fill) | Exercises the shared core through R only (D-01); no new production code |
| testthat edition (SC5) | Build/Package config (`DESCRIPTION`) | Test / Verification layer (94 files under `tests/testthat/`) | One DCF field plus test-file deprecation fixes; no `src/`, `R/`, or `python/` runtime code |
| Gate collection gap (`leafblower-x7n8`, P0, sequenced first) | Test / Verification layer (`conftest.py` semantics, file location) | — | Pure relocation; zero assertion changes permitted |

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|---------------|
| testthat | >= 3.0.0 (installed: 3.3.2, `[VERIFIED: local R env]`) | R test framework, edition 3 target | Already the project's `Suggests:` dependency (`DESCRIPTION:19`); no version bump needed |
| pytest | project-pinned via `python/pyproject.toml` | Python test framework | Already wired into the blocking gate command |
| R `stats` (base) | ships with R 4.6.1 | `rlnorm`/`rcauchy`/`rt` for the mixture-distribution generator (D-03) | Already `Imports: stats` in `DESCRIPTION`; `rlnorm` is already used in `tests/testthat/test-unit-bounds-status-consistency.R:20` for lognormal design weights — an established in-repo pattern, no new dependency |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| — | — | — | No supporting libraries required — this phase adds zero new dependencies to either `DESCRIPTION` or `python/pyproject.toml` |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| R-only property test (D-01, locked) | `hypothesis` (Python) generative testing | Rejected by CONTEXT.md D-01: bounds are enforced in the shared C++ core, so testing from one binding is sufficient; `hypothesis` is not a declared Python dependency (verified absent from `python/pyproject.toml`) and would be a net-new dependency for zero additional coverage of the code path in question |
| Fixed 50-dataset property test (D-02, locked) | Randomized-seed CI property test (classic QuickCheck-style) | Rejected by CONTEXT.md D-02: determinism is a hard project constraint (same rationale as pinned single-thread BLAS); a bound failure must be replayable without hunting a seed |
| Hand-rolled mixture RNG (D-03) | A CRAN heavy-tail package (e.g. `actuar`, `VGAM`) | Not needed — `rlnorm` (bulk) + `rcauchy` or `rt(df=1)` (contaminant) from base `stats` already produce the required skew/fat-tail shape and match the project's existing test idiom; adding a package for two RNG calls violates the "no new dependency" framing of D-01/D-03 |

**Installation:** None required — no new packages for either R or Python.

**Version verification:** `Rscript -e 'cat(as.character(packageVersion("testthat")))'` →
`3.3.2` (`[VERIFIED: local R env, this session]`), which satisfies `DESCRIPTION`'s existing
`testthat (>= 3.0.0)` bound. No `Suggests:` version change is needed for the edition-3 flip.

## Package Legitimacy Audit

**Not applicable — this phase installs zero external packages.** All test infrastructure
(`testthat`, `pytest`, R `stats`) is already a declared, installed, in-use dependency.
`[VERIFIED: DESCRIPTION:19, python/pyproject.toml, local R env]`

## Architecture Patterns

### System Architecture Diagram

```
                         Definition-of-Done gate
                    (.coverage-thresholds.json enforcement.command)
                                    |
        +---------------------------+---------------------------+
        |                                                       |
   R CMD INSTALL --preclean .                    cd python && uv pip install -e .
        |                                                       |
   Rscript test_dir("tests/testthat")             .venv/bin/python -m pytest -q
        |                                                       |
   94 test-*.R files                    +------> rootdir = python/  (pytest.ini/pyproject scope)
        |                                |                |
   [SC4] NEW: property-based             |         collects python/leafblower/test_*.py ONLY
   bound test (bounds_mode="unit",       |                |
   50 fixed datasets, R-only,            |         [x7n8, FIRST] test_parity_weights.py
   D-01..D-05)                           |         relocated FROM <root>/tests/ INTO
        |                                |         python/leafblower/ so THIS pytest run
   [SC5] DESCRIPTION                     |         collects it (previously: 8 tests silently
   Config/testthat/edition: 3            |         un-run, gate reports green regardless)
   (fix-first-then-flip, D-10..D-12)     |                |
                                          |         [SC1] parametrize list extended:
                                          |         {greenkhorn,logit,raking,oris,sinkhorn,
                                          |          newton_kl} -> +{chebyshev,greg,oris_soft}
                                          |         = full 9-entry rk_algorithm_t solver set
                                          |                |
                                          |         [SC2] test_solver_parity.py: +raking,
                                          |         +sinkhorn convergence-rule/max_error cases
                                          |         (mirrors test_logit_default_rule_parity)
                                          |                |
                                          +---------[SC3] logit tolerance: diagnose
                                                    src/logit_calib.cpp Newton+Armijo
                                                    conditioning BEFORE touching
                                                    tol=1e-6-vs-1e-10 ternary
                                                             |
                                                    both R and Python subprocess/in-proc
                                                    call the SAME shared C++17 core
                                                    (lbw namespace, calib_dispatch.hpp)
                                                    -- untouched this phase
```

### Recommended Project Structure

No new directories. File-level changes only, inside existing structure:

```
DESCRIPTION                          # +Config/testthat/edition: 3 (SC5, LAST commit of that thread)
tests/testthat/                      # existing 94 files; NEW: one property-based bound test file
tests/parity/                        # existing R-side helpers; untouched (relocation keeps resolving them)
python/leafblower/
├── test_parity_weights.py           # git mv target from <root>/tests/test_parity_weights.py (x7n8, FIRST)
└── test_solver_parity.py            # extended: +test_raking_*, +test_sinkhorn_* (SC2)
conftest.py                          # untouched — verify relocation preserves its sys.path strip behavior
```

### Pattern 1: The four-step parity protocol (reuse verbatim for SC2)

**What:** `python/leafblower/test_solver_parity.py` documents and implements a fixed
four-step protocol for every R↔Python weight/convergence comparison.
**When to use:** Any new parity test in this phase (`raking`, `sinkhorn` convergence-rule
checks) must follow it — CONTEXT.md's canonical_refs section makes this a hard reference,
not a suggestion.
**Example (existing, to be mirrored):**
```python
# Source: python/leafblower/test_solver_parity.py:171-186 (read this session)
def _assert_parity(method: str, conv_py: dict = _CONV_PY, conv_r: str = _CONV_R):
    w_py, ri_py = _run_py(method, conv_py)
    assert ri_py.get("max_error", 1.0) < _CONV_TOL, (
        f"Python {method} did not converge: max_error={ri_py.get('max_error')}"
    )
    r_out = _run_r(method, conv_r)
    assert r_out["max_error"] < _CONV_TOL, (
        f"R {method} did not converge: max_error={r_out['max_error']}"
    )
    w_r = np.array(r_out["weights"])
    assert len(w_py) == len(w_r), (
        f"{method}: length mismatch Python={len(w_py)} R={len(w_r)}"
    )
    assert np.allclose(w_r, w_py, rtol=1e-6, atol=0.0), (
        f"{method} R↔Python mismatch: max|Δw|={np.max(np.abs(w_r - w_py)):.3e}"
    )
```
The "precheck both sides converged" step (the two `assert ... < _CONV_TOL` lines) is the
step CONTEXT.md's Notes-for-planning explicitly calls out — a non-converged run must fail
loudly, not silently pass a comparison between two garbage iterates.

For the `raking`/`sinkhorn` default-rule cases specifically, mirror
`test_logit_default_rule_parity` (`python/leafblower/test_solver_parity.py:212-242`,
`[VERIFIED: python/leafblower/test_solver_parity.py:212-242]`), which additionally asserts
the *resolved rule string* on both bindings, not just weight parity — because a fixture that
converges to machine precision cannot distinguish `rule="threshold"` from
`rule="improvement"` by weights alone.

### Pattern 2: Unit-mode is the only mode with a hard per-observation bound guarantee (SC4)

**What:** `lbw::finalize_weights_buf` (`src/calib_dispatch.hpp`) dispatches on
`bounds_mode`. `RK_BOUNDS_CELL` (the default) only counts violations; `RK_BOUNDS_UNIT`
performs an iterative per-cell water-fill that clamps into `[min_weight, max_weight]`.
**When to use:** Any test asserting `min(w) >= min_weight && max(w) <= max_weight`
unconditionally MUST pass `bounds_mode = "unit"`.
**Example (existing bound-assertion idiom to copy):**
```r
# Source: tests/testthat/test-clamp-contract.R:51-58 (read this session)
for (m in c("raking", "sinkhorn", "greg", "oris", "chebyshev")) {
  w <- as.numeric(harvest(df, target, method = m, design_weights = dw,
                          max_weight = 2.0, min_weight = 0.3, bounds_mode = "unit",
                          max_iterations = 1000L, attach_weights = FALSE))
  expect_true(all(w <= 2.0 + 1e-9),  info = sprintf("%s: max=%.6f > 2.0", m, max(w)))
  expect_true(all(w >= 0.3 - 1e-9),  info = sprintf("%s: min=%.6f < 0.3", m, min(w)))
}
```
The `+1e-9` slack in this existing test is looser than KPI-02's stated `1e-10` tolerance —
the new property test must use `1e-10` per the requirement's own wording, which is one order
tighter than the existing idiom it otherwise copies. Flag this delta explicitly rather than
copy-pasting the `1e-9` constant.

**Verbatim confirmation this is not incidental** — `R/harvest.R` warns exactly this
distinction into existence:
```r
# Source: R/harvest.R:833-835 (read this session)
warning(sprintf(
  "cell-mode bounds: %d weights fell outside [%.3f, %.3f] due to skewed base weights within cells. Consider bounds_mode='unit' for strict per-observation bounds.",
  calib_result$n_bounds_violated, min_weight, max_weight))
```

### Pattern 3: Mixture-distribution generator for skewed/fat-tailed design weights (D-03)

**What:** CONTEXT.md D-03 requires design weights drawn from "lognormal bulk plus a small
heavy-tailed contaminant fraction," never Gaussian. The existing test suite already
establishes the lognormal half of this pattern.
**Example (existing lognormal idiom, extend with a contaminant draw):**
```r
# Source: tests/testthat/test-unit-bounds-status-consistency.R:20 (read this session)
d[x == "B"] <- rlnorm(sum(x == "B"), meanlog = 0, sdlog = 1.6)
```
```r
# Source: tests/testthat/test-clamp-contract.R:48-49 (read this session)
dw <- 2 ^ (rnorm(n))            # log-normal design weights: heavy skew per cell
dw[sample(n, 60L)] <- 50.0      # extreme outliers force cascading clamps
```
Both existing patterns already combine a lognormal-shaped bulk with an explicit outlier
overlay (`dw[sample(...)] <- <large constant>` or `<- rlnorm(..., sdlog=1.6)` on a subset) —
this is precisely D-03's "mixture" in miniature, already proven to trigger clamp activity in
this codebase. The new 50-dataset generator should generalize this exact idiom
(lognormal bulk `rlnorm(n, meanlog, sdlog)`, contaminant fraction replaced via
`rcauchy()`/`rt(df=1)` or a large fixed multiplier on a `sample()`-selected subset) rather
than inventing a new distributional shape.

### Anti-Patterns to Avoid

- **Testing bound enforcement under `bounds_mode="cell"` (the default):** produces a test
  that is *correctly* sometimes-failing by design, not a defect signal — see Pattern 2.
- **Renormalizing weights after the water-fill pass:** forbidden project-wide
  (`CLAUDE.md`, `src/calib_dispatch.hpp:359-390` comment) — not something this phase's test
  code does, but if a diagnostic script for SC3 momentarily reproduces solver internals to
  measure the logit conditioning, it must not violate this ordering.
- **Loosening a newly-failing 3e assertion to make it pass (D-11):** CONTEXT.md is explicit
  — every newly-failing waldo assertion under the edition flip is a finding to investigate,
  never blanket-triaged as "3e artifact."
- **Flipping `Config/testthat/edition: 3` before migrating deprecated constructs (D-12
  violation):** breaks the "every commit stays green" invariant this project's Definition of
  Done depends on (local-only repo, "complete = committed locally + gates green").

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| R↔Python subprocess comparison harness | A new comparison framework | The existing `_run_r`/`_run_py`/`_assert_parity` triplet in `test_solver_parity.py`, or the `_r_weights`/subprocess pattern in `test_parity_weights.py` | Both patterns already handle single-thread BLAS env propagation, `Rscript` invocation, JSON/CSV marshaling, and the convergence precheck — reinventing either duplicates ~80 lines of subprocess plumbing for zero benefit |
| Skewed/fat-tailed random data generation | A custom rejection sampler or copied-in distribution code | Base R `stats::rlnorm` + `stats::rcauchy`/`stats::rt` | Already imported (`DESCRIPTION` `Imports: stats`), already used in-repo for the same purpose (Pattern 3) |
| Cross-language float comparison | A hand-rolled epsilon-diff helper | `np.allclose(..., rtol=1e-6, atol=0.0)` (Python side) / `expect_lt(max(abs(...)), tol)` idiom (R side) | Both are the project's own established idiom throughout `tests/testthat/` and `python/leafblower/test_*.py` |

**Key insight:** Every capability this phase needs already has a proven, in-repo
implementation one file away. The work is disciplined copying with parameter changes
(method name, tolerance, bounds_mode), not new engineering — consistent with the phase
being explicitly "test-layer only."

## Common Pitfalls

### Pitfall 1: Extending the wrong parity matrix, or extending it before the relocation

**What goes wrong:** Adding `oris_soft` (or any method) to
`tests/test_parity_weights.py:73`'s parametrize list *before* `leafblower-x7n8` lands is a
no-op — that file is not collected by the blocking gate's `pytest` invocation (gate runs
`cd python && ... pytest -q`, whose rootdir is `python/`, which never walks up to
`<root>/tests/`).
**Why it happens:** The file passes today; a developer extending it sees local green and
assumes the gate covers it.
**How to avoid:** Sequence D-09 literally — relocate first (`git mv` to
`python/leafblower/test_parity_weights.py`, fix `REPO_ROOT`), confirm the gate's own command
now collects 149 tests (141 + 8, `[CITED: leafblower-x7n8 ticket body]`), THEN extend the
parametrize list.
**Warning signs:** `cd python && .venv/bin/python -m pytest --collect-only -q | tail -1`
still reports 141 (not 149) after the relocation commit — the move did not actually land
inside pytest's collected rootdir.

### Pitfall 2: `REPO_ROOT` path arithmetic breaks silently after relocation (D-08's "costly reversibility")

**What goes wrong:** `tests/test_parity_weights.py` currently computes
`REPO_ROOT = Path(__file__).resolve().parent.parent` (`[VERIFIED:
tests/test_parity_weights.py:31, read this session]`) — correct when the file lives at
`<root>/tests/test_parity_weights.py` (`.parent` = `tests/`, `.parent.parent` = repo root).
After the D-08 move to `python/leafblower/test_parity_weights.py`, that same expression
resolves to `.../python`, not the repo root, and `R_HELPER = REPO_ROOT / "tests" /
"parity" / "run_parity_r.R"` becomes `.../python/tests/parity/run_parity_r.R`, which does
not exist. The three R-helper paths in this file (`R_HELPER`, `_ORIS_SOFT_R_HELPER`,
`_CHEBYSHEV_R_HELPER`) ALL derive from `REPO_ROOT` and all three need the fix.
**Why it happens:** A `.parent.parent` relative-path idiom silently degrades to a wrong-but-
valid `Path` object; it does not raise ImportError.
**How to avoid:** Change to `Path(__file__).resolve().parents[2]` (three levels up:
`python/leafblower/<file>` → `leafblower/` → `python/` → repo root), and assert the helper
path exists at *module import time* (not lazily inside a test) so a wrong path fails
collection loudly.
**Warning signs:** Per CONTEXT.md D-08 itself — "getting it wrong makes the tests skip
rather than fail" — because `RSCRIPT_AVAILABLE` / helper-missing paths in this file are
wired to `pytest.skip`, a broken `REPO_ROOT` produces silent skips, not failures, which is
exactly the invisible-gap failure mode this phase exists to close.

### Pitfall 3: Treating the logit tolerance as a magic number to tune until green (D-06/D-07)

**What goes wrong:** `tol = 1e-6 if method == "logit" else 1e-10`
(`[VERIFIED: tests/test_parity_weights.py:93]`) has no comment explaining the three-order-
of-magnitude gap. The tempting fast path is bisecting the tolerance until the test passes.
**Why it happens:** No conditioning analysis exists yet; the number "already works" at
1e-6, so narrowing it to find the real floor feels like unnecessary work.
**How to avoid:** `src/logit_calib.cpp` runs a Newton solve with Armijo line search
(`kMaxHalvings=10`, `kArmijoC=0.01`, confirmed at `src/logit_calib.cpp:152-155,
[VERIFIED: src/logit_calib.cpp:152-155]`) and clamps `z` to `±700` before `exp()`
(`src/logit_calib.cpp:241-250`). Two independently plausible, testable mechanisms exist and
should be checked in order before touching the ternary: (a) the R-vs-Python **build
asymmetry already documented in `CONCERNS.md`** — Python compiles with `-O3`
(`python/CMakeLists.txt:99`) while R gets no explicit `-O` level (`configure`/`Makevars.in`
by design, `[VERIFIED: python/CMakeLists.txt:99]` for the Python side) — different
optimization levels can produce different FMA-contraction and vectorization decisions around
the Newton/Armijo dot products, producing systematic ULP-scale drift that an iterative
solver can amplify over its Newton steps; (b) genuine algorithmic sensitivity of the
logit-link Newton iteration itself (its D_eff floor, line-search halving count, or the `±700`
clamp) to input perturbation at a scale other solvers' fixed-point/IPM iterations don't
share. Measure which one it is (e.g., compare R-vs-R across `-O2`/`-O3` locally on the same
fixture using the AVX2/no-AVX2 boundary, or Python-vs-Python across build flags) before
deciding whether to tighten, document, or ticket.
**Warning signs:** A tolerance change with no comment naming a mechanism; per D-07, any
measured divergence above `1e-10` gets ticketed with the delta and suspected mechanism
recorded, not silently absorbed into a wider tolerance.

### Pitfall 4: The testthat edition flip breaking on `context()` or bare list equality, not just numeric drift

**What goes wrong:** 11 files use the deprecated `context()` (soft-deprecated under 3e,
`[VERIFIED: grep across tests/testthat/*.R, this session — files: test-alm-config-grouping.R,
test-cell-table.R, test-compare.R, test-oris-b12-fallback-best-reset.R,
test-oris-b13-best-error-honesty.R, test-oris-dispatch.R, test-oris-faithful.R,
test-oris-nonuniform-d.R, test-oris.R, test-oris-sraa-log-path.R, test-oris-sraa.R]`), and
276 `expect_equal()` call sites across 60 files (`[VERIFIED: grep -c across
tests/testthat/*.R, this session]`) switch from `all.equal()`-tolerance semantics to
waldo's stricter structural comparison. `expect_equivalent`/`expect_that`/`expect_is`
appear in **zero** files (confirmed by the same grep pass), so that specific migration
burden does not exist — but the `context()` and `expect_equal` surfaces do.
**Why it happens:** Edition 3 is silently opt-in per `Config/testthat/edition` in
`DESCRIPTION` (absent today, confirmed: `grep -n 'Config/' DESCRIPTION` returns nothing) —
nothing forces awareness of the 11 `context()` files or the 60 `expect_equal` files until
the flip actually happens.
**How to avoid:** Follow D-12's fix-first-then-flip order exactly: migrate `context()` calls
(3e replaces them with file-based grouping — deleting the call is usually sufficient since
testthat auto-derives context from the filename) and audit `expect_equal` call sites for
tolerance-dependent numeric comparisons WHILE STILL on edition 2 and green, commit that,
THEN add `Config/testthat/edition: 3` as its own isolated commit so the waldo fallout is
attributable to that commit alone.
**Warning signs:** A single commit that both migrates deprecated APIs and flips the edition
field — CONTEXT.md D-12 requires these as separate commits precisely so a bisect can
distinguish "mechanical churn" from "real waldo-exposed bug."

## Code Examples

### Extending the weight-parity matrix to the full 9-solver set (SC1)

```python
# Source: tests/test_parity_weights.py:73-75 (current state, read this session)
@pytest.mark.parametrize("method", [
    "greenkhorn", "logit", "raking", "oris", "sinkhorn", "newton_kl",
])
def test_weight_parity(method, tmp_path):
    ...
```
Target state — add the three missing enum-backed method strings
(`[VERIFIED: src/leafblower.h:42-52]`):
```python
@pytest.mark.parametrize("method", [
    "greenkhorn", "logit", "raking", "oris", "sinkhorn", "newton_kl",
    "chebyshev", "greg", "oris_soft",
])
```
Note `oris_soft` already has a *separate*, more elaborate parity test
(`test_oris_soft_default_tol_parity`) in the same file that checks `alm_capacity_mu_final`
in addition to weights — adding it to the simple parametrized list is additive coverage
(the plain weight-vector comparison at the default fixture), not a duplicate of that
existing targeted test.

### The full solver enum (ground truth for "how many solvers")

```c
// Source: src/leafblower.h:42-52 (read this session)
RK_ALG_AUTO   = 0,
RK_ALG_ORIS  = 1,
/* 2 = removed (was RK_ALG_LBFGSB) */
RK_ALG_RAKING    = 3,
RK_ALG_SINKHORN  = 4,
RK_ALG_CHEBYSHEV = 5,
RK_ALG_GREG      = 6,
/* 7 = reserved, undocumented hole per CONCERNS.md */
RK_ALG_ORIS_SOFT = 8,   /* oris + ADMM soft capacity enforcement */
RK_ALG_GREENKHORN = 9,   /* greedy coordinate-descent IPF (autumn::harvest style) */
RK_ALG_LOGIT      = 10,  /* Deville-Sarndal 1992 logit Newton calibration (autumn::calibrate style) */
RK_ALG_NEWTON_KL  = 11   /* Newton-KL smooth dual calibration (zero-compression regime) */
```
Nine non-`AUTO` values: `ORIS, RAKING, SINKHORN, CHEBYSHEV, GREG, ORIS_SOFT, GREENKHORN,
LOGIT, NEWTON_KL`.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|---------------|--------|
| testthat edition 2 (`all.equal` tolerance semantics) | testthat edition 3 (`waldo` structural diff) | This phase (SC5) | Stricter numeric/structural comparison surfaces previously-masked drift; requires the fix-first-then-flip protocol (D-12) |
| Additive Boyle-Dykstra raking projection | Multiplicative KL-Bregman Dykstra at cell level | Pre-dates this phase (`2026-04-27-raking-bregman-dykstra-design.md`) | Not this phase's concern, but explains why `raking`'s convergence-rule test (SC2) exercises a different geometry than the PRD originally described — irrelevant to correctness, relevant to not being surprised by unfamiliar internals while diagnosing SC2/SC3 |

**Deprecated/outdated:**
- `context()` in testthat: soft-deprecated since testthat 3e-adjacent releases; file-based
  grouping replaces it. 11 files in this repo still use it (see Pitfall 4).
- The PRD's iEPPA/FR-11–FR-19 "outer entropic-proximal-point loop" framing for `oris`:
  historical, mathematically inert at `C=0` per `REQUIREMENTS.md`'s US-005 status note —
  cited here only because FR-19's literal text ("weights satisfy `[min_weight, max_weight]`
  ... not just at convergence") is the phase's nominal requirement source, but the
  *implemented* mechanism that actually delivers that invariant on the returned vector is
  `bounds_mode="unit"` finalization (Pattern 2), not the historical BCD per-step clamp FR-19
  describes.

## Assumptions Log

| # | Claim | Section | Risk if Wrong |
|---|-------|---------|---------------|
| A1 | The `-O3`-vs-no-`-O` build asymmetry (and any resulting FMA-contraction difference) is the most likely mechanism behind logit's looser parity tolerance, ahead of algorithmic sensitivity in the Newton/Armijo iteration itself | Pitfall 3 / SC3 | If wrong, the diagnosis step wastes time comparing build flags instead of instrumenting the Newton iteration directly; D-06/D-07 already require measuring before deciding, so this ordering is a starting hypothesis, not a conclusion, and the cost of being wrong is bounded to one extra investigation branch |
| A2 | Specific lognormal `sdlog`/contaminant-fraction/contaminant-distribution values are NOT prescribed here — CONTEXT.md leaves these to Claude's Discretion at plan/execute time, derived from first principles against the default `max_weight=5` | Pattern 3 / Standard Stack | If the eventual parameters don't actually push weights against the clamps on enough of the 50 datasets, SC4's property test degenerates into a vacuous pass (never a real bound-violation stress test) — the planner must require empirical verification (e.g., print/assert that at least N of 50 fixtures produce a pre-water-fill clamp event) rather than trusting the derivation alone |
| A3 | `rcauchy()` or `rt(df=1)` is an adequate "heavy-tailed contaminant" generator alongside the lognormal bulk (D-03) | Standard Stack / Pattern 3 | Low risk — both are standard fat-tail generators in base R `stats`; if a different shape is wanted, swapping the contaminant draw is a one-line change with no structural impact |

**If this table is empty:** N/A — three assumptions logged above; none are structural risks
to the phase's success criteria.

## Open Questions

1. **Does SC1's "eight shipped solvers" literally mean 8, or is it a stale echo of
   `DESCRIPTION`'s 8-name prose list that predates `oris_soft` being split out as its own
   enum value?**
   - What we know: the live enum has 9 non-`AUTO` values (`[VERIFIED: src/leafblower.h:42-
     52]`); `DESCRIPTION`'s free-text `Description:` field lists exactly 8 solver names and
     does not mention `oris_soft` separately (folding it under "ORIS ... block-coordinate
     descent" prose); the currently-missing-from-parity set is 3 methods
     (`chebyshev`, `greg`, `oris_soft`), and 6+3=9.
   - What's unclear: whether SC1 intends the parametrize list to reach exactly 9 entries
     (full enum coverage) or exactly 8 (matching `DESCRIPTION`'s prose count, which would
     require deciding which 9th method to exclude — an odd outcome inconsistent with "every
     one of the eight ... `oris_soft` included").
   - Recommendation: plan for 9 — full enum coverage is the only reading consistent with
     "`oris_soft` included" as an *addition* to a base list rather than a member already
     counted in it, and it is strictly safer (more coverage, not less) if the planner's
     interpretation differs from the roadmap author's original count.

2. **What KPI-02-mandated location is "the property-based test in `test-harvest.R`"
   literally, given D-01's "sits beside the six targeted bound tests"?**
   - What we know: `REQUIREMENTS.md`'s KPI-02 row cites `test-harvest.R` by name (echoing
     the PRD's original §11 KPI table); `tests/testthat/` follows a strict one-file-per-
     behaviour convention (`TESTING.md`), and the six existing targeted bound tests are
     spread across `test-clamp-contract.R`, `test-harvest-bounds-mode.R`,
     `test-oris-bounds-mode.R`, `test-cr-d16-nbounds.R`,
     `test-unit-bounds-status-consistency.R`, `test-newton-bounds-write-guard.R` — none of
     which is `test-harvest.R` itself (a 173-line file already covering routing/API-compat
     behavior, `[VERIFIED: tests/testthat/test-harvest.R, wc -l]`).
   - What's unclear: whether the planner should add the new property test as a `test_that()`
     block appended to the existing `test-harvest.R`, or as a new file (matching the
     project's naming convention, e.g. `test-bound-property.R` or similar) placed "beside"
     the six targeted tests as CONTEXT.md phrases it.
   - Recommendation: a new dedicated file is more consistent with the project's own stated
     convention ("one file per behaviour," `TESTING.md`) and with D-01's "sit beside" framing
     than appending to the unrelated `test-harvest.R`; the planner should choose the new-file
     path and treat the PRD's literal filename as historical/non-binding, same treatment
     already given to other stale PRD specifics in `REQUIREMENTS.md`.

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-------------------|
| KPI-02 | Weight bound enforcement — `max(w) <= max_weight` and `min(w) >= min_weight` within `1e-10`, over 50 random datasets, by a property-based test. Currently "Not located" per `REQUIREMENTS.md`. | Pattern 2 (unit-mode is the only mode that can satisfy this unconditionally — critical, verified finding) + Pattern 3 (mixture-distribution generator, reusing the in-repo `rlnorm` idiom) + Pitfall 2's `1e-9`-vs-`1e-10` tolerance delta from the nearest existing test + Open Question 2 (file placement) |
</phase_requirements>

## Environment Availability

| Dependency | Required By | Available | Version | Fallback |
|------------|-------------|-----------|---------|----------|
| Rscript / R | R testthat gate, R-side parity subprocess calls | ✓ | R 4.6.1 (2026-06-24), `[VERIFIED: R --version, this session]` | — |
| testthat | Edition-3 target (SC5) | ✓ | 3.3.2 installed, satisfies `DESCRIPTION`'s `>= 3.0.0` bound, `[VERIFIED: packageVersion, this session]` | — |
| uv | Python venv-managed install per `CLAUDE.md` | ✓ | 0.12.3, `[VERIFIED: uv --version, this session]` | — |
| `python/.venv/bin/python` | Blocking gate's Python pytest step | ✓ | present, `[VERIFIED: ls, this session]` | — |
| git | `git mv` for the x7n8 relocation, `git diff -M` verification | ✓ | 2.55.0, `[VERIFIED: git --version, this session]` | — |

**Missing dependencies with no fallback:** None.
**Missing dependencies with fallback:** None — every tool this phase needs is already
installed and working in this environment.

## Validation Architecture

### Test Framework

| Property | Value |
|----------|-------|
| Framework (R) | testthat >= 3.0.0 (installed 3.3.2); edition flip is itself SC5's deliverable |
| Framework (Python) | pytest, invoked via the venv at `python/.venv/bin/python -m pytest` |
| Config file | `DESCRIPTION` (R — `Config/testthat/edition` field, currently absent); `python/pyproject.toml` (Python) |
| Quick run command | `Rscript -e 'testthat::test_file("tests/testthat/<new-file>.R")'` (single R file); `.venv/bin/python -m pytest python/leafblower/test_solver_parity.py -q` (single Python file) |
| Full suite command | `.coverage-thresholds.json`'s `enforcement.command` verbatim (quoted in Standard Stack / Summary above) |

### Phase Requirements → Test Map

| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|--------------------|-------------|
| KPI-02 | `min(w) >= min_weight` and `max(w) <= max_weight` within `1e-10` over 50 fixed, stratified, mixture-distributed datasets, `bounds_mode="unit"` | property-based (fixed-seed, 50-case sweep) | `Rscript -e 'testthat::test_file("tests/testthat/<new-file>.R")'` | ❌ Wave 0 — new file |
| SC1 (weight-parity matrix, all 9 solvers) | R↔Python weight-vector `rtol=1e-6` parity per method | integration (subprocess) | `.venv/bin/python -m pytest python/leafblower/test_parity_weights.py -q` (post-relocation path) | ❌ Wave 0 — extend after relocation |
| SC2 (raking/sinkhorn convergence-rule parity) | per-method default-rule resolution matches across bindings | integration (subprocess) | `.venv/bin/python -m pytest python/leafblower/test_solver_parity.py -q` | ❌ Wave 0 — new test functions in existing file |
| SC3 (logit tolerance) | documented or tightened tolerance, no unexplained ternary | integration (subprocess) + investigation | same as SC1's file, tolerance line only | ✅ existing assertion, comment/value changes |
| SC5 (testthat edition) | full 94-file suite green under 3e | full-suite regression | `.coverage-thresholds.json`'s R half | ✅ existing files, edition-only config change |
| `leafblower-x7n8` (P0, sequenced first) | `test_parity_weights.py`'s 8 tests collected by the blocking gate | collection-count assertion | `cd python && .venv/bin/python -m pytest --collect-only -q \| tail -1` (expect 149 post-move) | ❌ Wave 0 — relocation itself is the fix |

### Sampling Rate

- **Per task commit:** run the single new/changed file's quick command (fast, isolates the
  change under test).
- **Per wave merge:** full `.coverage-thresholds.json` `enforcement.command` (R build +
  testthat + Python install + pytest, single-thread BLAS).
- **Phase gate:** full suite green — collected test count for the Python half must read
  **149** (141 pre-existing + 8 relocated), not 141, before `/gsd-verify-work`.

### Wave 0 Gaps

- [ ] New R file under `tests/testthat/` for the KPI-02 property test (name TBD per Open
      Question 2 — planner's call, e.g. `test-bound-property.R`).
- [ ] `git mv tests/test_parity_weights.py python/leafblower/test_parity_weights.py` with the
      `REPO_ROOT` fix (Pitfall 2) — this IS `leafblower-x7n8` and must land first per D-09.
- [ ] Three new `@pytest.mark.parametrize` entries in the relocated
      `test_parity_weights.py` (`chebyshev`, `greg`, `oris_soft`).
- [ ] Two new test functions in `python/leafblower/test_solver_parity.py`
      (`test_raking_*`, `test_sinkhorn_*`), mirroring `test_logit_default_rule_parity`'s
      resolved-rule-plus-weight-parity shape.
- [ ] No framework install needed — testthat and pytest are both already present at the
      required versions.

## Security Domain

**`security_enforcement` assumed enabled (absent from `.planning/config.json`, which does
not exist in this repo — treated as default-enabled per the researcher protocol).** This
phase, however, touches zero attack surface: it adds test files and one DCF config field.
No new input parsing, no new network/credential/file-write code, no change to the C API's
untrusted-`double*`+`int n` boundary that `CONCERNS.md` already documents as the project's
actual attack surface (`[CITED: .planning/codebase/CONCERNS.md §Security Considerations]`).

### Applicable ASVS Categories

| ASVS Category | Applies | Standard Control |
|----------------|---------|--------------------|
| V2 Authentication | No | N/A — no auth surface in this package |
| V3 Session Management | No | N/A |
| V4 Access Control | No | N/A |
| V5 Input Validation | No (this phase) | Existing: `src/calib_validate.cpp` (`RK_ERR_BADARG`/`RK_ERR_INFEAS`), untouched this phase |
| V6 Cryptography | No | N/A — no crypto in this package |

### Known Threat Patterns for {stack}

None applicable — this is a local, single-machine statistical library with no network,
credential, or multi-tenant surface (`[CITED: .planning/codebase/CONCERNS.md]`). The only
project-relevant "threat" is numerical: untrusted-length raw-pointer input into the C API,
which is out of this phase's scope (test-layer only) and already has documented mitigation
+ regression coverage (`tests/testthat/test-bridge-length-checks.R`,
`test-input-validation.R`, `test-cr-d7-nobs-guard.R`, `test-diagnostics-guards.R`,
`test-cr-wave10-call-guards.R`).

## Project Constraints (from CLAUDE.md)

Extracted directives with direct bearing on this phase, from both the global and
project-level `CLAUDE.md`:

- **Build gate:** `R CMD INSTALL --preclean .` — NOT `devtools::install` — is the R build
  gate. This phase's changes (test files, `DESCRIPTION` field) must still pass this exact
  command.
- **Python venv:** `cd python && uv pip install -e . --reinstall-package leafblower` — uv-
  managed, no bare `pip`. No package installs are needed this phase, but any verification
  run must use this venv's `python`/`pytest`, never a bare `python`/`pytest` (a stale
  `~/.local` shadow `.so` gets imported otherwise).
- **Single-thread BLAS:** `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` must
  be set together, and in Python before `import numpy` — required for every parity/property
  test run in this phase to stay deterministic.
- **No cancellations:** project-wide numeric-stability rule. Not directly invoked by writing
  test files, but relevant if SC3's diagnosis touches computed quantities (e.g., diffing
  `w_r - w_py`) — prefer the existing `np.max(np.abs(w_r - w_py))` idiom already used
  throughout, which does not introduce a new cancellation-prone computation.
- **Commit discipline:** explicit pathspec commits only, never `git add -A` (a bd/graphify
  hook re-stages `.beads/issues.jsonl` into the index). Local-only repo — no `git push`.
  D-12's "fix first, then enable" edition-3 protocol naturally produces at least two
  separate commits (deprecation fixes, then the edition flip) — keep them that way rather
  than squashing.
- **Beads-only task tracking:** no `TodoWrite`/markdown TODO lists. `leafblower-og7d` (P1,
  testthat edition) and `leafblower-x7n8` (P0, gate-collection gap) already exist as beads
  tickets with detailed executable specs (read this session, quoted in full above via `bd
  show`); the planner should reference/claim these rather than re-describing their content,
  and must file new tickets for the SC1/SC2/SC3/SC4 gaps per the roadmap's "file tickets for
  the four parity/property gaps at plan time" instruction.
- **Alternatives-considered mandate (global CLAUDE.md):** satisfied above in Standard Stack
  — every candidate addition (hypothesis, a CRAN heavy-tail package, a custom comparison
  framework) was compared against zero-new-dependency reuse of existing in-repo patterns and
  rejected on that basis, consistent with CONTEXT.md's own "no new dependency" framing.

## Sources

### Primary (HIGH confidence — direct file reads this session)

- `.planning/phases/01-verification-coverage-closed/01-CONTEXT.md` — locked decisions D-01
  through D-12, canonical references, code context.
- `.planning/REQUIREMENTS.md` — KPI-02 status, solver/enum traceability.
- `.planning/codebase/CONCERNS.md`, `.planning/codebase/TESTING.md`, `.planning/ROADMAP.md`
  §Phase 1 — coverage-gap sourcing, test-framework conventions, phase success criteria.
- `.coverage-thresholds.json` — verbatim blocking gate command.
- `DESCRIPTION`, `R/harvest.R` (lines 260-290, 825-860), `src/leafblower.h` (enum + error
  codes), `src/calib_validate.hpp`, `src/calib_dispatch.hpp` (finalize_weights_buf),
  `src/logit_calib.cpp` (Newton/Armijo constants), `python/CMakeLists.txt` (`-O3`, AVX2
  flags), `configure`/`src/Makevars.in` (AVX2, no explicit `-O`).
- `tests/test_parity_weights.py`, `python/leafblower/test_solver_parity.py`,
  `tests/testthat/test-clamp-contract.R`,
  `tests/testthat/test-unit-bounds-status-consistency.R`, `conftest.py`,
  `python/conftest.py`, `python/pyproject.toml`.
- `tasks/prd-leafblower-core.md` §FR-19 and surrounding FR-11–FR-28 (read to resolve the
  KPI-02 mechanism question).
- `bd show leafblower-og7d`, `bd show leafblower-x7n8` — full existing ticket specs (this
  session).
- Local environment probes this session: `R --version`, `packageVersion("testthat")`,
  `uv --version`, `git --version`, venv presence.

### Secondary (MEDIUM confidence)

- None used — no web/external documentation was needed; this phase's domain is entirely
  in-repo conventions and existing test patterns.

### Tertiary (LOW confidence)

- The specific causal mechanism behind the logit tolerance gap (Pitfall 3 / Assumption A1)
  is a reasoned hypothesis from documented build-asymmetry facts, not a measured conclusion
  — flagged in the Assumptions Log.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH — zero new dependencies; every tool verified present and at a
  sufficient version in this session.
- Architecture: HIGH — every pattern cited is a direct, line-numbered read of code already
  in this repository, not an external convention.
- Pitfalls: HIGH for Pitfalls 1, 2, 4 (directly verified via file reads and grep counts);
  MEDIUM for Pitfall 3's causal mechanism (investigation guidance, not a proven conclusion —
  matches D-06/D-07's own "diagnose before touching" requirement).

**Research date:** 2026-08-15
**Valid until:** 30 days (stable, internal-conventions-only domain; no external library
version drift risk since zero new dependencies are introduced)
