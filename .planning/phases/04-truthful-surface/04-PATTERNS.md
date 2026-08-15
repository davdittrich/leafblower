# Phase 4: Truthful Surface - Pattern Map

**Mapped:** 2026-08-15
**Files analyzed:** 6 (4 to modify, 1 test to add-a-block-to, 1 test to fix)
**Analogs found:** 6 / 6 (all patterns exist in-repo already — this phase edits existing files,
it creates no new files, so every "analog" is the file's own surrounding pattern)

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|--------------------|------|-----------|-----------------|---------------|
| `docs/raking.md` (delete §8.2/§12 passage) | doc | transform (prose edit) | `docs/methods/oris.md` (ground-truth correction text) | exact — direct content source |
| `src/leafblower.h` (annotate enum slot 7) | config/header | N/A (comment-only) | itself, slot-2 comment 5 lines above (`src/leafblower.h:44`) | exact — same enum, established convention |
| `CLAUDE.md` (extend slot-2 line to also name slot 7) | config/doc | N/A | itself, existing "Algorithm slot 2 is reserved..." line | exact |
| `R/harvest.R` (add `weights=` `stop()` guard) | service/validation (R API surface) | request-response (arg validation at function entry) | itself — RVAL.2 block (`R/harvest.R:303-307`) and `design_weights` warning idiom (`R/harvest.R:510-515`) | exact — same function, adjacent validation block |
| `tests/testthat/test-harvest-rval.R` (add RVAL.4 block) | test | request-response (expect_error on bad arg) | itself — RVAL.2 block (`test-harvest-rval.R:35-46`) | exact — same file, same numbering convention |
| `tests/testthat/test-logit.R:195` (fix `weights=` → `design_weights=`) | test | request-response | itself — the call site being edited | exact — no analog needed, single-token arg rename |

## Pattern Assignments

### `docs/raking.md` (doc, prose delete+replace)

**Analog:** `docs/methods/oris.md:127` (ground-truth ORIS description, already correct)

**Correction text to source the replacement from** (`docs/methods/oris.md:127`, verbatim):
```
> **Not the iEPPA paper.** The solver was originally (mis)named "iEPPA" after Chu, Liang, Toh &
> Yang, *"An efficient implementable inexact entropic proximal point algorithm…"* (arXiv:2011.14312)
> [chu2022ieppa]. That paper's **headline contribution — an outer inexact entropic-proximal-point
> loop with a re-centered Bregman term to solve a transport-cost LP (`min⟨C,X⟩`) without `ε → 0`
> blow-up — is NOT implemented here.** Survey calibration has no transport cost (`C = 0`), at
> which the outer proximal loop is mathematically inert ... So citing arXiv:2011.14312 as the
> basis over-claims a contribution the code does not contain.
```

**Rule:** per D-01, delete — do NOT copy/paraphrase this passage into `docs/raking.md`.
`docs/methods/oris.md` remains the single authoritative description; `docs/raking.md` should
either drop the §8.2/§12 blocks entirely or (per the "12. Synthesis" narrowing in RESEARCH.md)
retain only the parts that don't restate the ORIS outer-loop claim or recommend L-BFGS-B.

**Structural note (from RESEARCH.md Pitfall 2):** the file is 10 physical lines, zero markdown
headings — section markers ("8.2", "12") are inline prose, not `#`-headers. Locate and replace
by exact quoted substring, not by line/byte offset. Verification greps:
`grep -o 'L-BFGS-B' docs/raking.md | wc -l` (expect 0 after edit) and
`grep -ci Gurobi docs/raking.md` (expect 0 after edit).

---

### `src/leafblower.h` (header, enum comment)

**Analog:** itself — the slot-2 comment 5 lines above the edit point, same enum block.

**Exact current enum** (`src/leafblower.h:40-53`, read in full):
```c
typedef enum {
    RK_ALG_AUTO   = 0,
    RK_ALG_ORIS  = 1,
    /* 2 = removed (was RK_ALG_LBFGSB) */
    RK_ALG_RAKING    = 3,
    RK_ALG_SINKHORN  = 4,
    RK_ALG_CHEBYSHEV = 5,
    RK_ALG_GREG      = 6,
    RK_ALG_ORIS_SOFT = 8,   /* oris + ADMM soft capacity enforcement */
    RK_ALG_GREENKHORN = 9,   /* greedy coordinate-descent IPF (autumn::harvest style) */
    RK_ALG_LOGIT      = 10,  /* Deville-Sarndal 1992 logit Newton calibration (autumn::calibrate style) */
    RK_ALG_NEWTON_KL  = 11   /* Newton-KL smooth dual calibration (zero-compression regime) */
} rk_algorithm_t;
```

**Pattern to copy** — mirror the slot-2 comment form exactly, insert between `RK_ALG_GREG = 6,`
and `RK_ALG_ORIS_SOFT = 8,`:
```c
    RK_ALG_GREG      = 6,
    /* 7 = removed (was RK_ALG_GRAKE) */
    RK_ALG_ORIS_SOFT = 8,   /* oris + ADMM soft capacity enforcement */
```
No enum values change — comment-only insertion. Editing this header forces a full rebuild of
both R (`R CMD INSTALL --preclean .`) and Python (`uv pip install -e . --reinstall-package
leafblower`) build sites (both `#include` it).

---

### `CLAUDE.md` (project doc, one-line extension)

**Analog:** itself — existing sentence under Architecture:
> "Algorithm slot 2 is reserved (LBFGSB removed). Do not reuse in `rk_algorithm_t` enum."

**Pattern:** extend the same sentence to also name slot 7/GRAKE, matching tense and phrasing,
e.g. "Algorithm slots 2 and 7 are reserved (LBFGSB and GRAKE removed, respectively). Do not
reuse in `rk_algorithm_t` enum." Do not split into a second bullet — one line, same convention.

---

### `R/harvest.R` (service/validation, request-response)

**Analog:** itself — the RVAL.2 block immediately following the function signature, and the
`design_weights` warning 200+ lines later (both already in this file, same validation style).

**Imports/signature context** (`R/harvest.R:266-303`, read in full — no `weights=` formal
exists; a bare `weights=w` call falls into `...`):
```r
harvest <- function(
  data,
  target,
  ...
  design_weights   = NULL,
  newton_tsvd_ratio = 1e-8,
  ridge_lambda = 0.0,
  ...
) {
  # RVAL.2: warn on unknown ... args (typos / removed params)
  dots <- list(...)
  if (length(dots) > 0L)
    warning("harvest: unknown argument(s) ignored: ",
            paste(names(dots), collapse = ", "), call. = FALSE)
```

**Core validation pattern to insert (before the generic RVAL.2 warning, same block)** — matches
the file's base-R `stop("leafblower: ...", call. = FALSE)` convention used by every other
`stop()` in this file (e.g. `R/harvest.R:508`: `stop("leafblower: 'data' must be a non-empty
data.frame", call. = FALSE)`):
```r
dots <- list(...)
if ("weights" %in% names(dots))
  stop("leafblower: unrecognized argument 'weights' — did you mean 'design_weights'? ",
       "harvest() takes per-observation design weights via design_weights=, not weights=.",
       call. = FALSE)
if (length(dots) > 0L)
  warning("harvest: unknown argument(s) ignored: ",
          paste(names(dots), collapse = ", "), call. = FALSE)
```

**Sibling informative-message idiom to model tone on** (`R/harvest.R:510-515`, verbatim):
```r
  # design_weights: used as start_weights when supplied (normalized to mean=1 by normalize_start_weights)
  if (!is.null(design_weights)) {
    if (!is.null(start_weights))
      warning("leafblower: both design_weights and start_weights supplied; design_weights ignored")
    else
      start_weights <- design_weights
  }
```

**Error handling pattern (package-wide):** every `stop()` in this file uses base R
`stop("leafblower: <msg>", call. = FALSE)` — no `cli` package (confirmed absent from
`DESCRIPTION` Imports/Suggests). Do not introduce `cli::cli_abort()`.

**No new validation dependency needed** — `Don't Hand-Roll` guidance from RESEARCH.md: a
hard-coded `if ("weights" %in% names(dots))` check is correct scope; a general
fuzzy/did-you-mean argument suggester is out of scope.

---

### `tests/testthat/test-harvest-rval.R` (test, new RVAL.4 block)

**Analog:** itself — the RVAL.2 block, same file, same `test_that("RVAL.N: ...")` numbering
convention.

**Exact style to copy** (`tests/testthat/test-harvest-rval.R:35-46`, RVAL.2 block, verbatim
pattern — read via grep-confirmed test_that headers, full body at these line numbers):
```r
# ---------------------------------------------------------------------------
# RVAL.2: warn on unknown harvest() ... args
# ---------------------------------------------------------------------------
test_that("RVAL.2: unknown ... arg emits warning listing the arg name", {
  set.seed(1)
  df  <- data.frame(x = factor(sample(c("a", "b"), 100, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  expect_warning(
    harvest(df, tgt, convergence = list(absolute = 1e-4),
            ...)
  )
})
```

**New block to add, RVAL.4** (following the same section-comment + test_that naming
convention; insert after the existing RVAL.3/META.2 blocks, matching file's existing numbering
order):
```r
# ---------------------------------------------------------------------------
# RVAL.4: bare weights= errors naming design_weights
# ---------------------------------------------------------------------------
test_that("RVAL.4: bare weights= arg errors naming design_weights", {
  set.seed(1)
  df  <- data.frame(x = factor(sample(c("a", "b"), 100, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  expect_error(
    harvest(df, tgt, weights = runif(100), convergence = list(absolute = 1e-4)),
    regexp = "design_weights"
  )
})
```
Full list of existing `test_that` blocks in this file (for insertion-point/ordering reference):
RVAL.1 ×2 (lines 6, 20), RVAL.2 ×2 (lines 35, 47), RVAL.3 ×2 (lines 67, 92), META.2 ×2 (lines
113, 124), dtkn.6 ×2 (lines 134, 142).

---

### `tests/testthat/test-logit.R:195` (test, one-arg rename)

**Exact call site to fix** (`tests/testthat/test-logit.R:182-202`, full `test_that` block read):
```r
test_that("eb79.18: consistent collinear with HETEROGENEOUS design weights still reaches feasibility", {
  ...
  base_w <- runif(n, 0.5, 2.0)
  w <- suppressWarnings(harvest(data, target, method = "logit", weights = base_w,
               min_weight = 0.05, max_weight = 20, max_iterations = 500L,
               attach_weights = FALSE))
  r <- attr(w, "result")
  expect_equal(r$status, 0L, label = "heterogeneous consistent collinear: converges (eb79.18)")
  expect_lt(r$max_error, 1e-4)
})
```
**Fix:** `weights = base_w` → `design_weights = base_w`. Per RESEARCH.md Pitfall 1, this MUST
land in the SAME commit as the `R/harvest.R` guard (guard alone breaks this test today via
`suppressWarnings()` currently masking the RVAL.2 warning it trips). After the rename,
empirically re-run the test to decide whether `suppressWarnings()` is still needed — the
sibling INCONSISTENT-margins test (`test-logit.R:164-179`) also wraps in `suppressWarnings()`
with no `weights=` argument at all, so a second independent warning may be legitimate; do not
strip `suppressWarnings()` by inspection alone.

---

## Shared Patterns

### R error-message convention
**Source:** `R/harvest.R` (every existing `stop()` call, e.g. lines 447, 451, 455, 508)
**Apply to:** the new `weights=` guard
```r
stop("leafblower: <message>", call. = FALSE)
```
Base R only — no `cli` package dependency exists in this project.

### C header enum-removal comment convention
**Source:** `src/leafblower.h:44` (`/* 2 = removed (was RK_ALG_LBFGSB) */`)
**Apply to:** slot 7 annotation
```c
/* N = removed (was RK_ALG_<NAME>) */
```

### NEWS.md breaking-change entry style (discretionary — only if checkpoint:decision opts in)
**Source:** `NEWS.md:22-26, 47-51, 80-84` (three existing `## Breaking changes` sections)
```
## Breaking changes

* `harvest()` default changed: `sor` is now `NULL` (disabled) instead of
  `list(auto = TRUE, omega_min = 0.3)`. ...
```
If the D-04 checkpoint decides to log the `weights=` guard as a breaking change, add a new
bullet under a `## Breaking changes` heading following this exact style (top-of-file
dev-version section, bullet starts with the function/behavior changed, states old vs new).

### testthat `test_that("<ID>: <description>", { ... })` numbering convention
**Source:** `tests/testthat/test-harvest-rval.R` (RVAL.1/RVAL.2/RVAL.3/META.2 blocks) and
`tests/testthat/test-logit.R` (`eb79.NN` ticket-ID-prefixed blocks)
**Apply to:** the new RVAL.4 block — reuse the file's own RVAL.N id sequence, not a new scheme.

## No Analog Found

None — every file in this phase's scope is an edit to an existing file, and every edit has a
same-file or adjacent-file precedent to copy (comment style, `stop()` idiom, or `test_that`
numbering). No new file, module, or abstraction is created by this phase.

## Metadata

**Analog search scope:** `docs/raking.md`, `docs/methods/oris.md`, `src/leafblower.h`,
`CLAUDE.md`, `R/harvest.R`, `tests/testthat/test-harvest-rval.R`, `tests/testthat/test-logit.R`,
`NEWS.md`.
**Files scanned:** 8 (all read directly this session per RESEARCH.md's own verified-source list;
no new reads outside that set were needed since RESEARCH.md already file:line-verified every
finding).
**Pattern extraction date:** 2026-08-15
