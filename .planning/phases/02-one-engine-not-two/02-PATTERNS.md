# Phase 2: One Engine, Not Two - Pattern Map

**Mapped:** 2026-08-15
**Files analyzed:** 5 (3 modified C++ core files + 2 new test files)
**Analogs found:** 5 / 5

This phase is a pure internal C++17 refactor (no new files in the "new module" sense
except two new regression-prevention test files). Every modified file already contains
its own best analog — the target pattern (`finalize_weights`/`finalize_weights_buf` in
`calib_dispatch.hpp`) is proven, existing, in-repo code, not an external reference.

## File Classification

| New/Modified File | Role | Data Flow | Closest Analog | Match Quality |
|-------------------|------|-----------|-----------------|---------------|
| `src/calib_dispatch.hpp` (ADD: dispatch table + neutral result struct + extraction fn) | utility/service (shared C++ header helper) | transform (solver selection + result marshal) | `src/calib_dispatch.hpp:359-433` (`finalize_weights`/`finalize_weights_buf`, same file, existing section) | exact — same file, same "single shared helper both FFI callers invoke" shape |
| `src/c_api.cpp` (SLIM: `rk_calibrate()` dispatch body → call shared table) | controller (C-ABI entry point) | request-response | `python/leafblower/_bindings.cpp:194` (already calls `rk_calibrate()`, the target shape) — and structurally, `src/c_api.cpp`'s own current `pack_solver_result`/`pack_oris_result_c`/`pack_newton_result_c` templates (lines 57-173) as the packing-pattern to extend, not replace | role-match (bindings.cpp is the target consumer pattern); exact (self, for packing-template shape) |
| `src/r_bridge.cpp` (SLIM: `C_rk_calibrate()` dispatch body → call shared table) | controller (R `.Call()` entry point) | request-response | `src/c_api.cpp:414-440` (`RK_ALG_RAKING` enum-dispatch branch — the target per-solver shape R's `strcmp` branches converge to) | role-match — same role (dispatch+pack), different current data flow (string vs enum), converging to the same target |
| `python/leafblower/test_core_sources_sync.py` (NEW — SC4 build-list sync test) | test (regression-prevention, file-list comparison) | batch (static text comparison, no compiled artifact) | `python/leafblower/test_returned_weights_invariant.py` (small, single-assertion Python test file already in the target directory) | role-match — same directory, same "narrow single-purpose pytest module" shape; RESEARCH.md already supplies the exact target content (Code Examples section) |
| SC3 verification test (grep/assertion that all solver `.cpp` files call `finalize_weights`/`finalize_weights_buf`) | test (static/structural verification, not runtime) | batch (source-text scan) | `python/leafblower/test_core_sources_sync.py`'s own pattern (regex-over-source-tree) is the closest shape once written — no pre-existing test of this exact kind exists yet | role-match — same "assert on source text, not runtime behavior" family as the SC4 test |

## Pattern Assignments

### `src/calib_dispatch.hpp` (utility header, transform)

**Analog:** same file, `finalize_weights`/`finalize_weights_buf` (lines 359-433) — this is
the ONLY existing precedent in the codebase for "one shared C++ function, called from
both `r_bridge.cpp` and every solver `.cpp`, instead of a duplicated per-caller copy."
The new dispatch table + neutral result struct must live in this same header (per
`CLAUDE.md`'s "calib_dispatch.hpp = canonical home for shared solver helpers" rule) and
follow the same `inline` function shape — no new `.cpp` translation unit (see
RESEARCH.md's Alternatives Considered: a new TU would need its own `CORE_SOURCES` entry
and buys nothing since dispatch is cold/once-per-solve).

**Core pattern to copy (lines 359-433, verbatim structure)**:
```cpp
inline void finalize_weights_buf(double* w, int n, const CalibState& st,
                                 const CellTable& ct,
                                 int& n_bounds_violated, int& n_bounds_clamped) {
    // ... single implementation body ...
}

/// Convenience overload operating on st.weights (the final iterate).
inline void finalize_weights(CalibState& st, const CellTable& ct,
                             int& n_bounds_violated, int& n_bounds_clamped) {
    finalize_weights_buf(st.weights, st.n, st, ct, n_bounds_violated, n_bounds_clamped);
}
```
The new dispatch table should mirror this two-tier shape: one core `inline` function
that does the real work (dispatch-by-enum → solver call → extract into neutral struct),
plus thin convenience wrappers as needed by each caller — never a duplicated body.

**Existing per-solver call sites to preserve exactly** (the pattern the new dispatch
table must NOT disturb): `oris_finalize.cpp:150-159` ("Delegate to the single source"
comment), plus one call each in `raking.cpp`, `chebyshev.cpp`, `greenkhorn.cpp`,
`greg.cpp`, `logit_calib.cpp`, `sinkhorn.cpp`. `newton_calib.cpp` has zero calls
(Pitfall 2 — flag, do not silently "fix" by adding a call as part of this task unless
scoped explicitly).

---

### `src/c_api.cpp` (controller, request-response)

**Analog:** `python/leafblower/_bindings.cpp:194` (target shape: caller invokes
`rk_calibrate()`/shared dispatch, does not reimplement solver selection) — and, for the
packing-side shape to extend rather than discard, `c_api.cpp`'s own existing templates.

**Existing per-solver dispatch branch to generalize (lines 414-440, verbatim)**:
```cpp
if (alg == RK_ALG_RAKING) {
    // Classical raking: IPF + Dykstra box + Dykstra hyperplane (renamed from ORIS)
    auto res = lbw::raking_solve(st);
    status = res.base.status;
    iterations = res.base.iterations;
    max_error = res.base.max_error;
    used = RK_ALG_RAKING;
    if (result) {
        result->mean_error          = res.base.mean_error;
        result->kl                  = res.base.kl;
        result->chi2                = res.base.chi2;
        result->l1_weight_change    = res.base.l1_weight_change;
        result->grake_norm          = res.base.grake_norm;
        result->convergence_metric  = res.base.convergence_metric;
        result->convergence_rule    = res.base.convergence_rule;
        result->convergence_tol     = res.base.convergence_tol;
        // ... (continues; res.base.* -> rk_result_t field-by-field)
```

**Existing shared-pack templates to reuse, not reinvent (lines 57-173)**:
```cpp
template <typename R>
static void pack_solver_result(rk_result_t* dst, const R& src, rk_algorithm_t alg) noexcept {
    if (!dst) return;
    dst->status                       = src.base.status;
    // ... base fields ...
    if constexpr (has_n_bounds_c<R>::value) {
        dst->n_bounds_violated = src.n_bounds_violated;
        dst->n_bounds_clamped  = src.n_bounds_clamped;
    }
}
```
`pack_oris_result_c` (95-135) and `pack_newton_result_c` (146-173) are the same pattern
specialized per solver family — the neutral-struct extraction helper in
`calib_dispatch.hpp` should absorb this template family's logic (read `res.base.*` +
solver-specific fields into the neutral struct once), after which `c_api.cpp` narrows
the neutral struct into `rk_result_t` via a thinner version of these same templates.

**AUTO-routing block to consolidate (lines 280-334, verbatim)** — the target single copy
of AUTO routing logic that `r_bridge.cpp`'s duplicate (see below) converges toward:
```cpp
switch (p->algorithm) {
    case RK_ALG_RAKING:   alg = RK_ALG_RAKING; break;
    // ... one case per enum value ...
    case RK_ALG_AUTO:
    default: {
        int M_cell_est = lbw::estimate_M_cell(n, K, group_ids, cat_counts);
        if (static_cast<int64_t>(M_cell_est) * 10 >= static_cast<int64_t>(n) * 9) {
            // ... target_skew computation + severe-skew branch ...
        } else {
            alg = RK_ALG_ORIS;
        }
        auto_selected = true;
        break;
    }
}
```

**kAlgNames table (lines 30-47, verbatim)** — must be unified with `r_bridge.cpp`'s
divergent copy (different display strings for the same enum), not kept as two:
```cpp
static constexpr const char* kAlgNames[] = {
    "auto", "ORIS", "(reserved)", "raking", "sinkhorn", "chebyshev", "greg",
    "(gap)", "ORIS-soft", "greenkhorn", "logit", "newton_kl"
};
static_assert(sizeof(kAlgNames)/sizeof(kAlgNames[0]) == 12,
    "kAlgNames must cover all 12 enum slots 0..11");
```

---

### `src/r_bridge.cpp` (controller, request-response)

**Analog:** `src/c_api.cpp:414-440`'s per-solver enum-dispatch branch (above) — the shape
each `strcmp` branch converges to once R routes through the shared table.

**Current pattern being replaced (lines 654-663, verbatim — string dispatch, bypasses
`rk_calibrate()` entirely)**:
```cpp
if (strcmp(method_str, "raking") == 0) {
    auto res = lbw::raking_solve(st);
    res_status     = res.base.status;
    res_iterations = res.base.iterations;
    res_max_error  = res.base.max_error;
    res_alg_used   = (int)RK_ALG_RAKING;
    pack_solver_result(res);
    res_sraa_demoted = res.sraa_demoted ? 1 : 0;
    res_best_weights = std::move(res.base.best_weights);
} else if (strcmp(method_str, "auto") == 0) {
    // ... 124 more lines of AUTO routing duplicated from c_api.cpp:292-397
```

**Existing string→enum map to reuse as-is (lines 94-105, verbatim — already correct, do
not duplicate into the new shared table, call it)**:
```cpp
const std::unordered_map<std::string_view, rk_algorithm_t> kAlgMap = {
    {"oris",       RK_ALG_ORIS},
    {"oris_soft",  RK_ALG_ORIS_SOFT},
    {"raking",     RK_ALG_RAKING},
    {"greg",       RK_ALG_GREG},
    {"chebyshev",  RK_ALG_CHEBYSHEV},
    {"sinkhorn",   RK_ALG_SINKHORN},
    {"auto",       RK_ALG_AUTO},
    {"greenkhorn", RK_ALG_GREENKHORN},
    {"logit",      RK_ALG_LOGIT},
    {"newton_kl",  RK_ALG_NEWTON_KL},
};
```

**kAlgNames divergence to resolve (lines 943-957, verbatim)** — different strings for
same enum than `c_api.cpp`'s copy; unify per Pattern 2 in RESEARCH.md:
```cpp
static const char* kAlgNames[] = {
    "",           // 0 = RK_ALG_AUTO
    "oris",       // 1 = RK_ALG_ORIS
    "",           // 2 = (removed lbfgsb slot)
    "raking",     // 3 = RK_ALG_RAKING
    // ... 12 entries total, same reserved-slot holes at 2 and 7
};
static const int kAlgNamesLen = 12;
static_assert(RK_ALG_NEWTON_KL == 11, "kAlgNames table needs update on enum change");
```

**RAII-unwind discipline that MUST survive the refactor (error-handling pattern)** — any
new heap-backed member added to the neutral result struct must be included in BOTH swap
blocks below (`Rf_error` performs `longjmp`, skipping C++ destructors — R-exts §5.5):
- Error path: `r_bridge.cpp:914-938` — explicit `std::vector<T>().swap(...)` /
  `std::string().swap(...)` for every heap-backed local before `Rf_error()`.
- Success path: `r_bridge.cpp:983-994` — same swap-before-allocation discipline, since
  `Rf_allocVector`/`Rf_mkChar` can themselves trigger `Rf_error` (OOM).

This is the single highest-risk pattern to get wrong in this file — re-verify against
every `Rf_error()` call site after the neutral struct gains fields, per RESEARCH.md
Pitfall 4.

---

### `python/leafblower/test_core_sources_sync.py` (test, batch/regression-prevention)

**Analog:** `python/leafblower/test_returned_weights_invariant.py` (small, single-file,
single-purpose pytest module already in the same directory) for file/module shape; the
exact target content is already fully specified in RESEARCH.md's Code Examples section
(D-05) and requires no further pattern search — it is a pure regex-over-source-tree
comparison with zero runtime/compiled dependency:

```python
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

def test_core_sources_matches_src_glob():
    src_files = {p.name for p in (REPO_ROOT / "src").glob("*.cpp")} - {"r_bridge.cpp"}
    cmake_text = (REPO_ROOT / "python" / "CMakeLists.txt").read_text()
    listed = set(re.findall(r"\.\./src/(\w+\.cpp)", cmake_text))
    assert src_files == listed, (
        f"src/*.cpp vs CORE_SOURCES drift: "
        f"missing from CORE_SOURCES={src_files - listed}, "
        f"extra in CORE_SOURCES={listed - src_files}"
    )
```
Per D-05/A2, R-side placement (`testthat`) is an equally valid discretionary choice — if
chosen, mirror `tests/testthat/test-newton-kl.R`'s `test_that(...)` wrapper shape
(header comment naming the file + purpose, one `test_that` block, `expect_*` assertions)
instead of the pytest `def test_*` shape above.

---

### SC3 verification test (grep/assertion — no filename locked by CONTEXT.md)

**Analog:** same shape as `test_core_sources_sync.py` above — a static source-text scan,
not a runtime numerical test. Suggested placement: `python/leafblower/test_finalize_weights_sync.py`
(D-05's sibling, same directory, same regex-over-`src/*.cpp` technique) or an R-side
equivalent, per Claude's Discretion in CONTEXT.md.

**Pattern to copy**: reuse `test_core_sources_sync.py`'s `Path(...).glob("*.cpp")` +
`re.findall` approach, but scan each solver `.cpp` file's own text for a
`finalize_weights` / `finalize_weights_buf` call instead of comparing two file lists:
```python
_SOLVERS_EXPECTED_TO_CALL_FINALIZE = {
    "oris_finalize.cpp", "raking.cpp", "chebyshev.cpp", "greenkhorn.cpp",
    "greg.cpp", "logit_calib.cpp", "sinkhorn.cpp",
}  # newton_calib.cpp deliberately excluded — see RESEARCH.md Pitfall 2
for fname in _SOLVERS_EXPECTED_TO_CALL_FINALIZE:
    text = (REPO_ROOT / "src" / fname).read_text()
    assert "finalize_weights" in text, f"{fname} no longer calls the shared finalize helper"
```
Per D-04, if the planner finds this list needs a 9th entry (`newton_calib.cpp`) added
(Pitfall 2 scope decision (a)/(b)), that is a locked-scope decision for the planner/user,
not an implementation detail — the test's solver set must match whatever the plan
decides, and Pitfall 2's exclusion must be documented in the test itself if excluded
(inline comment, as shown above), not silently omitted.

---

## Shared Patterns

### Shared-helper-in-header (the single architectural pattern this whole phase applies)
**Source:** `src/calib_dispatch.hpp:359-433` (`finalize_weights`/`finalize_weights_buf`)
**Apply to:** the new dispatch table, the new neutral result-extraction struct/function —
every one of them belongs in this same header as an `inline` function, called identically
from `c_api.cpp` and `r_bridge.cpp`, exactly as `finalize_weights` already is from 7
solver `.cpp` files. This is the ONE pattern to imitate for the entire phase; do not
invent a second shared-code mechanism (no new `.cpp` TU, no virtual dispatch, no macro
generation).

### RAII-unwind discipline (r_bridge.cpp only)
**Source:** `src/r_bridge.cpp:914-938` (error path), `:983-994` (success path)
**Apply to:** any new heap-backed field added anywhere in the R-side result marshaling
path as part of routing through the shared dispatch table.

### `kAlgMap` string→enum (already correct — reuse, do not duplicate)
**Source:** `src/r_bridge.cpp:94-105`
**Apply to:** the shared dispatch table's string-resolution step, if the table needs to
accept a string (R side) as well as an enum (both sides) — reuse this map rather than
re-deriving it inside `calib_dispatch.hpp`.

## No Analog Found

None. Every file this phase touches already has a directly-applicable in-repo analog —
either itself (for the packing-template shape to extend) or a sibling file already
implementing the target pattern (`_bindings.cpp` for "call the shared dispatch, don't
reimplement it"; `finalize_weights` for "shared header helper, single source of truth").

## Metadata

**Analog search scope:** `src/calib_dispatch.hpp`, `src/c_api.cpp`, `src/r_bridge.cpp`,
`src/leafblower.h`, `python/leafblower/_bindings.cpp`, `python/leafblower/test_*.py`,
`tests/testthat/test-newton-kl.R`
**Files scanned:** 8 (all direct-read, verbatim-verified against `master` @ `2b94e8e`,
2026-08-15)
**Pattern extraction date:** 2026-08-15
