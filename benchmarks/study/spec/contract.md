# Adapter Contract & Result/Registry Schemas

**schema_version:** `1.0.0`
**Status:** FROZEN — this document + the three sibling JSON schemas
(`status_enum.json`, `runs_schema.json`, `registry_schema.json`) are the
interface every adapter, driver, and reporting WU under `benchmarks/study/`
must conform to. Reference: `docs/benchmark/DESIGN.md` §8 (harness
architecture), §5 (measurement protocol, two co-headline timing axes).

Boundary: this WU freezes the spec only. No adapter/driver implementation
lives here.

---

## 1. Adapter call contract

Every adapter (leafblower and every competitor, both R and Python arms)
exposes one entry point with an identical signature and return shape:

```
run(problem) -> {
  weights_ref:      str,        # see §2.1
  iterations:        int | NA,   # see §2.2
  status:            <enum>,     # see status_enum.json, §2.3
  converged:         bool,       # see §2.4 — NEVER the solver's self-report
  error_message:     str | null, # see §2.5
  wall_time_e2e_s:   f64,        # see §2.6
  wall_time_solve_s: f64 | null, # see §2.6
  peak_rss_bytes:    i64,        # see §2.7
}
```

`problem` is a loaded problem spec (per DESIGN.md §3 `spec/*.json` +
resolved `data_ref`). The adapter is responsible for translating the
problem's margins/bounds/tol into its wrapped package's native call and
for catching every exception the package can raise — `run()` itself must
never throw; all failure modes are reported through `status` +
`error_message`.

### 1.1 Two co-headline timing axes (DESIGN.md §5, Blocker A)

Both axes are recorded on every row; a solver that has no separable
aggregation seam reports `wall_time_solve_s = null`, never `0` or `NA`
string.

- **`wall_time_e2e_s`** — end-to-end, raw-rows-in to weights-out. Always
  measured, for every arm, with the groupby/contingency-table/model-matrix
  construction step INSIDE the timer. This is the fairness-neutral axis:
  it does not privilege solvers (like the IPF family) that natively accept
  pre-aggregated tables over whole-unit solvers (survey/sampling/ebal) that
  have no aggregation seam to exclude.
- **`wall_time_solve_s`** — solve-only, pre-aggregated table in to weights
  out. Measured ONLY for solvers that accept a pre-materialised table:
  leafblower (via the WU-0 CellTable entry point) and IPF-family
  competitors (ipfn/ipfr/anesrake/weightipy). `null` for whole-unit
  solvers (survey, sampling, ebal, balance, optweight, sbw, raw OT
  formulations) that have no such seam. Isolates solver-core speed from
  data-marshalling cost.

Both axes are co-equal primary figures; neither is a derived/optional
column.

---

## 2. Field semantics

### 2.1 `weights_ref: str`

Reference (relative path) to the parquet file holding this run's
length-`n` weight vector. Filename convention (DESIGN.md §8):

```
weights/<solver>__<problem>__t<thread>__<build>.parquet
```

`weights_ref` is **never null**. The `(solver,problem,thread,build)` cell
is solved once and the resulting weight file is shared by every `rep` of
that cell (weights are deterministic across reps — only timing varies);
`weights_ref` is therefore identical across all rows sharing a cell.
On a hard failure where no weights were ever computed (e.g.
`status == "bad_arg"` before any iterate exists), the adapter MUST still
write a length-`n` all-`NaN` sentinel vector to the conventional path
rather than omitting the file — downstream consumers gate on `status` /
`converged` before trusting the numeric content, but the file must exist
so `runs.parquet` never carries a dangling reference.

### 2.2 `iterations: int | NA`

The solver's own exposed iteration count (`RK_CalibResult.iterations` for
leafblower; package-native equivalent for competitors). Typed null (NA) in
`runs.parquet` — never the literal `0` or the string `"NA"` — where the
wrapped package does not expose an iteration count at all (e.g. most
whole-unit `survey`/`sampling` calibration calls).

### 2.3 `status: <harmonized enum>`

The adapter's own outcome classification, mapped from the wrapped
package's native exit code / exception class into the 8-value harmonized
enum defined in `status_enum.json`. This is what the *solver* reported.
It is deliberately kept separate from `converged` (§2.4): a package can
exit with a "converged" status and still fail the harness's own
independently-recomputed accuracy check, and conversely a package that
exits `stall`/`no_conv` may still have produced weights that meet the
harness's tolerance.

### 2.4 `converged: bool`

**Always harness-recomputed, never the solver's self-report.**
The harness independently loads `weights_ref` after the run completes and
recomputes the achieved margin error against the problem's margins, using
a single **uniform L∞ metric** (max absolute margin discrepancy across all
`k` targets), compared against the problem's `tol`. `converged = TRUE` iff
that independently-recomputed L∞ error is `<= tol`; `FALSE` otherwise —
including when `status` itself is `"error"`/`"bad_arg"`/`"infeasible"` (no
usable weights ⇒ `converged = FALSE` by construction). This uniform-axis
check is deliberately blind to each solver's own native stopping
quantity (duality gap, KL, etc.); the per-solver tol→native-quantity
mapping is recorded separately (WU-9's `tol_mapping.json`) as a diagnostic,
not as an input to this field. Rationale: never trust solver self-reported
convergence (project-wide, no-hiding-failures policy).

### 2.5 `error_message: str | null`

Diagnostic text captured on any non-clean exit: `conditionMessage(e)` in R
adapters, `str(exception)` in Python adapters, or the leafblower C API's
`message[128]` buffer content. `null` iff the run produced no diagnostic
text at all (clean `status == "converged"` exit with no warnings).
Presence of a message does not by itself imply `status != "converged"`
(a solver may emit a non-fatal warning on an otherwise-clean exit).

### 2.6 `wall_time_e2e_s: f64`, `wall_time_solve_s: f64 | null`

Seconds, `float64`, measured in-process with a high-resolution timer
(`bench::mark`-class / `time.perf_counter`), per DESIGN.md §5. Never
includes interpreter/process startup — only the timed region itself
(groupby+solve for e2e; solve-only for solve, per §1.1). See §1.1 for the
nullability rule on `wall_time_solve_s`.

### 2.7 `peak_rss_bytes: i64`

High-water-mark resident set size (`/proc/self/status:VmHWM` or
platform equivalent) of the subprocess that performed this run, sampled
at process exit, in bytes. This column stores the **raw** VmHWM of the
timed subprocess; the "data-loaded" baseline (package-import +
spec-load, exits before solving) used to disambiguate solver working-set
from library heft is captured and subtracted separately, per
`(language, package)`, as a downstream reporting artifact (DESIGN.md §5,
§11 DoD) — it does not mutate this raw column, keeping the frozen
per-run schema unambiguous.

---

## 3. `runs.parquet` — tidy one-row-per-(solver,problem,thread,build,rep)

Full column list + dtypes: `runs_schema.json`. Core identity + contract
columns are frozen (§1); denormalised quality-metric columns
(`marg_kl_mean`, `ess`, `deff`, `linf`, `bound_viol_mean`, ... —
DESIGN.md §6, owned by `common/metrics.{R,py}`, a separate WU) are an
open extension point (`additionalProperties: true` in the schema) so the
metrics WU can add columns without re-freezing this contract.

`build` is one of `{portable, native, na}` — `na` for competitor rows
(the build-variant axis only applies to leafblower).

---

## 4. `registry.json` — solver applicability

Full schema: `registry_schema.json`. One entry per solver id:

```json
{
  "<solver_id>": {
    "arm": "R" | "python",
    "families": ["kl", "chi2", "logit", "minimax", "ot", ...],
    "bounds": "unbounded" | "bounded" | "both",
    "K_max": <int> | null,
    "builds": ["portable", "native"] | ["na"],
    "applicable_problems": ["<problem_id>", ...]
  }
}
```

Built against **installed** package versions (WU-10), not aspirational
ones; gates the WU-8 run matrix so a solver is never scheduled against a
problem/family/bound combination it cannot handle. `K_max` is the largest
margin count `K` the solver is validated for (`null` = unbounded); for
`ott-jax`/2-marginal-OT-only competitors this caps applicability to `K<=2`
problems.

---

## 5. Status enum reference

See `status_enum.json` for the enum values, per-value semantics, and the
per-adapter mapping template (with leafblower's `RK_ERR_*` → enum mapping
filled in as the worked example).

---

## 6. Versioning

`schema_version` (semver) is carried in:

- this file's header (source of truth for the doc bundle version),
- a top-level `schema_version` key inside `status_enum.json`,
  `runs_schema.json`, and `registry_schema.json`,
- `registry.json`'s own top-level `schema_version` key (emitted by the WU
  that generates it, copied from this bundle's version at generation
  time),
- `runs.parquet` file-level (Arrow) key-value metadata — not a per-row
  column, to avoid redundant storage across millions of rows.

A breaking change to any field's name, type, or nullability bumps the
MAJOR version and requires every downstream adapter/driver WU to re-sync
before its next timed run (DESIGN.md §11: spec+registry+metrics+
hyperparams+tol-mapping are git-frozen together at WU-9T, before the
WU-REH rehearsal; this contract's version is part of that frozen bundle).
