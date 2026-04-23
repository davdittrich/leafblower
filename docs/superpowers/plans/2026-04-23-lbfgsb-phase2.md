# Plan: L-BFGS-B Phase 2 — O(n) SIMD + Epic Hygiene + Stall Investigation (rev 3)

## Context

Three follow-up threads. Rev 3 fixes iter-2 critical bug: WU-A3's
`reduction(-:x)` with in-loop `-=` yields `orig + Σ` (not `orig - Σ`)
under OpenMP combiner semantics `a = a - private`. Rewritten with
separate accumulators using `reduction(+:x)`. Also addresses: (a)
WU-A0 line numbers corrected (358→391); (b) kk1.21 children explicitly
accounted for in WU-B1; (c) scenario-4 rationale added; (d) grep widened
to `src/*.h`; (e) R CMD check failure protocol added; (f) WU-C3
thresholds reframed as qualitative smoke-test (literature gives no
hard numbers).

Processing order per user directive ("first 2, 3, then 1"):
thread-numbering in this plan: **A=SIMD (user #2), B=hygiene (user #3),
C=stall (user #1)**. Execution sequence: **A → B → C**.

## Scope

| Thread | Files | Risk |
|---|---|---|
| A: SIMD | `src/lbfgsb_solver.cpp` only | Reversible per-hint; benchmark-gated |
| B: Epic hygiene | `bd` metadata only | No code |
| C: Stall | scratch + optionally `src/lbfgsb_solver.cpp` | Low |

No API/ABI changes. No new dependencies. No new public tests.

## Non-scope

- No change to iEPPA solver logic, R/Python bindings, link-function math,
  ALM activation, **Wolfe acceptance/bracketing/zoom algorithm**.
  WU-A3 reorganises the Wolfe-evaluator **loop body** for SIMD but does not
  alter step acceptance, bracket updates, or zoom iteration logic.
- No change to harvest.R post-normalization contract.

## Verification gate (all WUs)

- `R CMD INSTALL --preclean .` DONE after each WU.
- `devtools::test()` ends `[ FAIL 0 | WARN 2 | SKIP 0 | PASS 100 ]`.
- Per-SIMD: benchmark delta vs post-`ded10c6` baseline. Combined SIMD
  change accepted only if ≥ 3 % speedup on **n ≥ 500k, K ≥ 10** AND no
  cell regresses by > 1 % (including small-n cells n ∈ {200k, 500k}).

---

# Thread A — O(n) SIMD

## WU-A0: Aliasing proof (one-time, prerequisite for `__restrict__`)

**Purpose.** Before adding `__restrict__` qualifiers, prove `u_base`, `du`,
`u_work` never alias at any callsite. Pure reasoning, no code change.

**Allocation sites** (`src/lbfgsb_solver.cpp:391-394`, current line numbers post-audit-fix commits):
```cpp
std::vector<double> u(st.n);       // u at current lam
std::vector<double> du(st.n);      // per-obs directional derivative
std::vector<double> u_work(st.n);  // scratch for Wolfe trial
```
Three distinct `std::vector<double>` objects. Each owns an independent
heap buffer via default `std::allocator`.

**Callsites** (`wolfe_line_search` at `:393-394`, `wolfe_zoom` via
`:323-325`, `:337-339`):
- `wolfe_line_search(..., u, du, u_work, ...)` → parameters bind to
  `u_base`, `du`, `u_work` at `:270-272`.
- `wolfe_zoom(..., u_base=u, du, u_work, ...)` → same distinct objects.

No call passes the same `std::vector` for two of these parameters. `du`
is written only by `compute_du` at `:389` from `dir` (separate buffer).
`u` is written only by `compute_u` at `:369, :419` and never aliased to
the scratch `u_work`.

**Conclusion.** `__restrict__` on all three pointers is safe. No UB risk.
Record this proof in the WU-A2 commit body.

**Acceptance.** No code change. Proof text lives in WU-A2 commit.

---

## WU-A1: Benchmark harness + baseline

**Target.** Scratch R script at `/tmp/lbfgsb-simd-bench/bench.R` (not
committed). Runs `harvest(method="lbfgsb")` across a grid:

| n | K | cats per margin |
|---|---|---|
| 200k | 4 | (5, 4, 6, 3) — D=18 |
| 500k | 10 | 5 each — D=50 |
| 500k | 20 | 5 each — D=100 |
| 1e6  | 10 | 5 each — D=50 |
| 1e6  | 20 | 5 each — D=100 |

Plus **stall-sentinel** rows (recorded but not gated): n=1000 near-sat,
n=10000 tight-tol. Purpose: detect iteration-count regression in small
cells.

**Metrics per cell.** Median and min of 5 runs. Iteration count. Saved
to `/tmp/lbfgsb-simd-bench/baseline.rds`.

**Skip rule.** If 1e6×K=20 > 60 s wall-clock, drop that cell. If
additionally 500k×K=20 > 60 s, drop it too. WU-A4 gate then evaluates
on remaining n ≥ 500k, K ≥ 10 cells. If ONLY 500k×K=10 survives,
require ≥ 3 % speedup on that sole cell. If no cell ≥ (500k, K≥10)
survives: revert WU-A2/A3 regardless — the claim cannot be validated.

**Acceptance.** 5+ cells timed, `baseline.rds` written. Iteration counts
recorded per cell (required for small-n regression check).

---

## WU-A2: SIMD u-update axpy (4 sites)

**Target.** `src/lbfgsb_solver.cpp` lines 245, 288, 323, 374. Pattern:
```cpp
        const double a = alpha;  // or alpha_accepted
        const double* __restrict__ ub = u_base.data();
        const double* __restrict__ dv = du.data();
        double*       __restrict__ uw = u_work.data();
#if LBW_HAS_OMP_SIMD
#pragma omp simd
#endif
        for (int i = 0; i < st.n; i++) uw[i] = ub[i] + a * dv[i];
```

All 4 sites are structurally identical. Aliasing safety per WU-A0.

**Acceptance.**
- `grep -c "#pragma omp simd" src/lbfgsb_solver.cpp` ≥ 4.
- Build clean, 100 tests pass.
- Benchmark gate deferred to WU-A4.

---

## WU-A3: SIMD phi/slope evaluator loops (branch-hoisted)

**Target.** `src/lbfgsb_solver.cpp` loops at `:252-263` (wolfe_zoom) and
`:330-341` (wolfe_line_search).

**Scope note.** This reorganises the Wolfe-evaluator **loop body only**.
Step acceptance, bracket updates, zoom iteration remain byte-identical.
Algorithmic Wolfe logic is NOT changed.

**Current structure** (inside each loop):
```cpp
for (int i = 0; i < st.n; i++) {
    double Fi, Hi;
    if (fn.exponential) { auto fh = fn.FH(u_work[i]); Fi = fh.F; Hi = fh.H; }
    else                { Fi = fn.F_from_e(e_vec[i]); Hi = fn.H_from_e(e_vec[i], u_work[i]); }
    phi_trial -= d[i] * Hi;
    slope     -= d[i] * Fi * du[i];
    if (st.alm_mu > 0.0) { sum_w += d[i]*Fi; sum_dw += d[i]*fn.dF(u_work[i])*du[i]; }
}
```

**Reorganisation.** Hoist `fn.exponential` branch (loop-invariant); keep
the ALM branch scalar (dead path, `alm_mu == 0.0` forced). Use
**separate accumulators with `reduction(+:...)`** to avoid the OpenMP
`reduction(-:x)` combiner bug: for `-` reduction the private copies
initialize to 0 and the final combiner is `x = x - private`, which paired
with in-loop `x -= …` gives `x = x_orig + Σ`, not `x_orig - Σ`. The fix:
accumulate contributions positively into fresh locals under `reduction(+:...)`
(associative, identity=0), then subtract once after the loop.

```cpp
double phi_acc = 0.0;       // accumulates Σ d[i]*Hi
double slope_acc = 0.0;     // accumulates Σ d[i]*Fi*du[i]
if (fn.exponential) {
#if LBW_HAS_OMP_SIMD
#pragma omp simd reduction(+:phi_acc) reduction(+:slope_acc)
#endif
    for (int i = 0; i < st.n; i++) {
        double e = std::exp(std::min(u_work[i], 700.0));
        phi_acc   += d[i] * e;
        slope_acc += d[i] * e * du[i];
    }
} else {
#if LBW_HAS_OMP_SIMD
#pragma omp simd reduction(+:phi_acc) reduction(+:slope_acc)
#endif
    for (int i = 0; i < st.n; i++) {
        double ei = e_vec[i];
        double Fi = (fn.L*(fn.U-1.0) + fn.U*(1.0-fn.L)*ei) / ((fn.U-1.0) + (1.0-fn.L)*ei);
        double Hi = fn.L*u_work[i] + (fn.U-fn.L)/fn.logit_scale
                    * std::log(((fn.U-1.0) + (1.0-fn.L)*ei) / (fn.U-fn.L));
        phi_acc   += d[i] * Hi;
        slope_acc += d[i] * Fi * du[i];
    }
}
phi_trial -= phi_acc;
slope     -= slope_acc;
// ALM path: dead at runtime (alm_mu=0.0 forced) but kept for future
// reactivation. Runs its OWN scalar loop with separate accumulators
// (sum_w, sum_dw); does NOT depend on phi_acc/slope_acc scope.
if (st.alm_mu > 0.0) {
    double sum_w = 0.0, sum_dw = 0.0;
    for (int i = 0; i < st.n; i++) {
        double Fi = fn.exponential ? fn.FH(u_work[i]).F : fn.F_from_e(e_vec[i]);
        sum_w  += d[i] * Fi;
        sum_dw += d[i] * fn.dF(u_work[i]) * du[i];
    }
    double residual = sum_w - static_cast<double>(st.n);
    double alm_scale = st.alm_lambda + st.alm_mu * residual;
    phi_trial += st.alm_lambda * residual + (st.alm_mu / 2.0) * residual * residual;
    slope     += alm_scale * sum_dw;
}
```

**e_vec precondition.** Both loop sites (`:252`, `:330`) are preceded by
`if (!fn.exponential) lbw::bulk_scaled_exp(fn.logit_scale, u_work.data(),
e_vec.data(), st.n)` at `:243` and `:319` respectively (verify with
`grep -n "bulk_scaled_exp" src/lbfgsb_solver.cpp` → expect 2 hits in
Wolfe paths). `e_vec[i]` is safe to read in the logit branch.

Inlining `F_from_e`/`H_from_e` bodies gives the compiler direct visibility
for SIMD. `std::log` inside the else branch depends on `LBW_HAS_GLIBC_MVEC`
(confirmed in lbw_config.h:6) linking `_ZGVdN4v_log`; if auto-vectorization
fails, the branch-hoist alone still removes the per-iter conditional — a
win even without SIMD.

**FP non-associativity.** `reduction(+:phi_acc)` and `reduction(+:slope_acc)`
reorder additions across SIMD lanes. Results differ from scalar in the
final ULP (Kahan summation is not used here). The Wolfe test compares
`phi_trial` against `phi_0 + c1·α·slope_0` with `c1=1e-4`; ULP-level
perturbation cannot cross this threshold. 100-test regression is the gate.

**Acceptance.**
- Tests pass.
- `grep -c "pragma omp simd" src/lbfgsb_solver.cpp` ≥ 8 (WU-A2 adds 4
  axpy pragmas; WU-A3 adds 2 sites × 2 branch variants = 4 pragmas;
  total 8).
- Benchmark gate deferred to WU-A4.

---

## WU-A4: Combined SIMD benchmark + decision

**Gate (from rev 1, unchanged).** Accept WU-A2+A3 only if:
- n ≥ 500k, K ≥ 10 cell shows ≥ 3 % wall-clock speedup, **AND**
- no cell (including 200k, stall-sentinels) regresses by > 1 % wall-clock
  OR by > 5 % iteration count.

**Bisect rule.** If combined gain ≥ 3 % but only one of A2/A3 contributes
(determined by reverting each separately and re-measuring): keep the
contributing WU, revert the other.

**Acceptance.**
- Delta table (all cells) recorded in merge commit body.
- If reverted: outcome logged in `leafblower-wje` (existing SIMD issue).

---

# Thread B — Epic hygiene (per-child verification, rev 2)

**Premise correction (iter-1 finding).** The `kk1.2X` phase trackers have
**open children** that are **duplicates** of already-closed `kk1.X`
atomic tickets (e.g., `kk1.20.1 T10 iEPPA` ≡ closed `kk1.11 T10 iEPPA`).
Each open child requires individual verification before closing.

## WU-B1: Verify + close duplicate kk1.2X children

**Per-child verification protocol.**

For each open child, execute:
1. Read the child's title/description.
2. Find the claimed deliverable in the codebase (code file, test, config,
   or artifact). Use `grep`/`ls`/`bd show` for the matching closed atomic.
3. **Close only if** the deliverable exists on master AND a closed
   atomic covers the same work.
4. **Leave open** if deliverable is missing, partial, or unclear.

**Per-child checklist:**

| Child | Claim | Verification command | Action if ✓ |
|---|---|---|---|
| `kk1.20.1` (T10 iEPPA) | iEPPA solver exists | `test -f src/ieppa.cpp && grep -q 'ieppa_solve' src/ieppa.hpp` | close, supersede=`kk1.11` |
| `kk1.20.2` (T11 auto-routing) | Auto-routing wired | `grep -q 'RK_ALG_AUTO' src/c_api.cpp` | close, supersede=`kk1.10`; note: `method="auto"` removed from public R/Python APIs per `leafblower-pvn` |
| `kk1.20.3` (T12 min_weight test) | min_weight end-to-end | `grep -rq 'min_weight' tests/testthat/` | close, supersede=`kk1.12` |
| `kk1.20.4` (T13 1M-benchmark) | 1M-row bench | Check for artifact in `benchmarks/` | close if artifact present, else leave open |
| `kk1.22.3` (T9 R CMD check) | Phase 1 gate | Run `R CMD check --as-cran .` locally (one-time) | close if passes, else leave open with evidence |
| `kk1.23.1` (T14 Python bindings) | pybind11 bindings | `test -f python/leafblower/_bindings.cpp` | close, supersede=`kk1.14` |
| `kk1.23.2` (T15 Python harvest) | Python harvest + diagnose_weights | `grep -q 'def harvest' python/leafblower/_harvest.py && grep -q 'def diagnose_weights' python/leafblower/_harvest.py` | close, supersede=`kk1.15` |
| `kk1.23.3` (T16 wheel build + tests) | Python tests pass | `cd python && python -m pytest 2>&1` | close if passes, else leave open |
| `kk1.24.1` (T17 configure) | configure script | `test -f configure || test -f src/Makevars.in` | close, supersede=`kk1.17` |
| `kk1.24.2` (T18 CRAN packaging) | `R CMD check --as-cran` clean | Run check, verify 0 ERROR/WARNING | close if clean, else leave open |
| `kk1.24.3` (T19 Final gate) | CRAN + PyPI artifacts | Check `leafblower_*.tar.gz`, wheel artifacts | close if artifacts present, else leave open |

**Close command template:**
```bash
bd close <child-id> --reason="Deliverable verified on master: <specific evidence>. Duplicate of closed atomic <kk1.X>; supersede recorded."
```

## WU-B2: Close phase trackers (kk1.2X) after all children resolved

**kk1.21 pre-condition (explicit — iter-2 blocker fix).** All 6 children
(`kk1.21.1` through `kk1.21.6`) are already ✓ per `bd show` at plan
authoring time. These do NOT appear in the WU-B1 per-child table because
the table only covers OPEN children requiring verification. Closed
children need no action. Close `kk1.21` unconditionally with reason
"all 6 atomic children complete: T1-T6".

For each remaining phase tracker (`kk1.20`, `kk1.22`, `kk1.23`, `kk1.24`):
**binary rule** — if ALL children (including those closed earlier in
WU-B1) are ✓, close the tracker. If ANY child remains open, leave the
tracker open and log which child(ren) block. No partial-complete close.

**R CMD check / pytest failure protocol** (iter-2 scope fix). When WU-B1
runs `R CMD check --as-cran` or `python -m pytest` as a close gate,
**unrelated failures** (e.g., NOTEs on `.wolf/` dev artifacts, non-standard
root files) are logged as NEW bd issues (priority P3, title
"R CMD check NOTE: &lt;summary&gt;") and do NOT block closing the child.
**Related failures** (e.g., test failure on iEPPA when closing `kk1.20.3`
"min_weight end-to-end") block closure and leave the child open.

## WU-B3: Close `leafblower-pvn` with code-presence check

**Target.** `leafblower-pvn`: "remove AUTO routing + sum(w)=n hard
constraint". All 3 phase deps ✓.

**Verification** (commands, read-only):
- Auto removed from public APIs:
  - `grep -n '"auto"' R/harvest.R` → 0 hits in user-facing paths
  - `grep -n '"auto"' python/leafblower/_harvest.py` → 0 hits
- sum(w)=n:
  - iEPPA Dykstra hyperplane: `grep -q 'q_hyp' src/ieppa.cpp` → 1+ hits
  - L-BFGS-B ALM scaffolding: `grep -q 'alm_mu' src/lbfgsb_solver.cpp` → 1+ hits (ALM inactive, verified by `alm_mu = 0.0;` forced)

**Action.** If all checks pass: `bd close leafblower-pvn --reason="..."`.
Else: leave open with specific unmet check.

---

# Thread C — Stall-risk investigation

## WU-C1: Probe instrumentation (temporary, gated)

**Target.** Add a temporary counter and `REprintf` inside
`src/lbfgsb_solver.cpp`, guarded by `#ifdef LBW_STALL_PROBE`.

```cpp
#ifdef LBW_STALL_PROBE
        static thread_local int n_rejects = 0;
        if (curv_ok) n_rejects = 0;
        else { n_rejects++; if (n_rejects % 10 == 0) REprintf("stall: %d consec rejects\n", n_rejects); }
#endif
```

Build with `-DLBW_STALL_PROBE` via a scratch `src/Makevars.local` (not
committed). Ephemeral.

**Primary metric.** Iteration count at termination (SIMD-independent; not
affected by Thread A). Wall-clock is secondary.

**Baseline.** Run the current `src/lbfgsb_solver.cpp` (post-Thread A if
A is merged; pre-Thread A if A is reverted). Either way, iteration counts
reflect the curvature-gate behaviour only.

**Acceptance.**
- Probe built, stall counts recorded for the 3 scenarios (WU-C2).
- No probe-related change committed to git.

## WU-C2: Adversarial scenarios + measurement

Run 3 scenarios each against two gate variants: (a) current HEAD with
`kCurvRel=1e-8` and (b) temporarily applied `kCurvMin=1e-20` (one-file
edit, build, measure, revert):

1. **Near-saturation**: n=1000, 1 margin, target (0.95, 0.05), max_weight=1.2.
2. **High-collinearity**: n=10000, 4 margins with overlapping design.
3. **Tight-tolerance**: n=10000, `convergence=list(absolute=1e-12)`.
4. **Sparse-margin**: n=10000, 1 margin with 20 categories, skewed
   target (geometric decay p_j = 0.5^j normalised). **Distinct from
   collinearity (scenario 2)**: scenario 2 stresses the Hessian
   approximation across multiple correlated margins; scenario 4 stresses
   it within a single margin where many categories have near-zero
   targets, producing near-flat curvature along those lambda
   dimensions. Different failure modes of `s·y` smallness.

Record per (scenario, variant): iteration count, `max_error`, wall-clock,
stall-reject count from probe.

**Results saved** to `/tmp/lbfgsb-stall/results.rds` AND the issue body:
`bd update leafblower-370 --notes="<table>"`.

## WU-C3: Decision + pre-commit cleanup

**Decision framing** (iter-2 blocker fix: thresholds reframed qualitatively).
Neither Liu-Nocedal (1989) nor Nocedal-Wright (2e §7.2) provides
hard numerical thresholds for "excessive curvature rejection"; the
literature treats repeated rejection as a diagnostic signal without
fixed bounds. Decision rule is therefore a **binary smoke test**:

- **Close leafblower-370** if NO scenario shows any of: (a) iteration
  count ≥ 2× the `kCurvMin=1e-20` baseline, (b) stall probe reports
  any run of ≥10 consecutive rejects, (c) max_error regression > 10×
  vs baseline at termination. Rationale: L-BFGS-B is typically 5-20×
  faster than pure steepest ascent on concave problems; a 2× iteration
  slowdown still retains most of the L-BFGS advantage.
- **Escalate (file new WU)** if ANY scenario triggers (a), (b), or (c).
  New WU chooses between: relax to `kCurvRel=1e-6` (Liu-Nocedal's
  looser bound), add consec-reject counter + steepest-ascent fallback,
  or tighten line search. Decision deferred to that WU with evidence.

**Pre-commit guard.** Before any commit in Thread C (iter-2 fix: widen
glob to include `.h` files):
```bash
grep -q 'LBW_STALL_PROBE' src/*.cpp src/*.hpp src/*.h && { echo "FAIL: probe not removed"; exit 1; }
```
Run as a scripted gate; block commit if instrumentation remains.

**Acceptance.** `leafblower-370` updated with evidence and decision. No
`LBW_STALL_PROBE` strings in committed source.

---

## Ordering

1. **Thread A (SIMD)**: WU-A0 (aliasing proof) → WU-A1 (baseline) → WU-A2
   (u-update) → WU-A3 (phi/slope) → WU-A4 (gate). Single commit or
   revert per gate outcome.
2. **Thread B (hygiene, metadata only)**: WU-B1 (per-child) → WU-B2
   (trackers) → WU-B3 (pvn). No git commits.
3. **Thread C (stall)**: WU-C1 (probe) → WU-C2 (scenarios) → WU-C3
   (decision + cleanup). Code commit only if C3 selects "tighten" or
   "guard". Otherwise doc-only update to `leafblower-370`.

Threads can pause between WUs; user interrupt-safe.

## Rollback

Thread A: per-WU `git revert`. Thread B: `bd update <id> --status=open`.
Thread C: either no commit or single revert + re-run pre-commit grep.
