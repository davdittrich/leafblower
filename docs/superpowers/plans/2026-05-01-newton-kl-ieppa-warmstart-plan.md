# Newton-KL IEPPA Warm-Start — Plan (rev 1)

**Epic:** `leafblower-usg8` (Epic-C; follow-up to Epic-B BLOCKED)
**Spec reference:** `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` (a new "IEPPA warm-start" subsection lands in WI-0)
**Predecessors:** epic `leafblower-5k08` (LM rev 2; landed); `leafblower-91u7` (Epic-B target homotopy; BLOCKED — see `docs/investigations/2026-05-01-newton-kl-homotopy-result.md`)
**Tasks (sequential):** WI-0 (spec) → WI-1 (ieppa-inner extraction) → WI-2 (warm-start wiring) → WI-3 (T7 + T8 regression tests) → WI-4 (verify) → WI-5a (bench) → WI-5b (verdict + investigation doc)
**Date:** 2026-05-01

---

## Mechanism Header

- **Mechanism:** Run IEPPA's coordinate-descent inner sweep for a small fixed budget (default `K_warm = 8` outer sweeps) on the *original* targets `T`, convert IEPPA's log-factor vector `lf` to Newton's reference-eliminated `λ` parameterisation, then enter the existing LM-Newton inner with that warm `λ`. Newton polishes from inside the convergence basin.
- **Forbidden:** Target homotopy of any kind (Epic-B BLOCKED). Touching the IEPPA coordinate-descent loop body. Loosening any test gate. Disabling SRAA on the IEPPA warm-start (default ON; rationale in §Mechanism). Touching `harvest()` public signature. Per-method internal homotopy-bounds reuse.
- **Audit:** Two new `NewtonCalibResult` fields — `n_warmstart_iters_used` (count of IEPPA outer sweeps actually run, including any early-exit-on-converged) and `warmstart_max_err_at_handoff` (the IEPPA-side `max_err` measured at handoff). T2 and T8 read both via testthat snapshots; T7 asserts `n_warmstart_iters_used >= 1 && status == 0` (loose, not brittle).

---

## Context

### What landed before this epic

- **Epic 5k08 (LM rev 2; `5355f69`):** Newton-KL solver with LM damping (gain-ratio adaptive `lm_mu`, scale-invariant Hessian damping with additive floor, Armijo line search, best-iterate restoration). Converges on stepstone K=9 to a *natural floor* `gap ≈ 1.4e-4` → `max_err ≈ 2.8e-4` (vs `<1e-4` gate). The floor is the rank-deficient-Hessian drift wall: Newton steps past optimum into unbounded-dual regions where best-iterate rolls back.
- **Epic-B (`91u7`; `5355f69` TH-1a refactor only):** target homotopy regressed T2 from 2.8e-4 to 2.29e-3 (8× worse). Mechanism failure: per-margin-uniform `T_eps` for the stepstone overlap structure (margins `rk_age10_gender`, `rk_gender_time`, `rk_i_loc_time_gender`) lands `λ*(T_eps→0+)` in a *structurally distinct basin* from `λ*(T_0)`. Warm-starting Newton from there is worse than starting from `λ=0`. Implementation discarded; master at TH-1a state.
- **TH-1a refactor (kept):** `src/newton_calib.cpp:139` defines `auto run_newton_inner = [&](const std::vector<double>& T_eps, int max_iter_inner) -> bool` capturing all outer-scope state by reference (`lam`, `lm_mu`, `Z_curr`, `u_max_curr`, `res`, `H`, `G`, `delta`, etc.). Single top-level call at line 371: `run_newton_inner(T, max_iter);`. `NewtonCalibResult` carries `n_homotopy_levels_used = 0` (stale) plus `lm_mu_final`, `dual_gap`, `step_norm`, `line_alpha`. This refactor is the affordance Epic-C consumes — no further changes needed in the Newton inner.

### Why IEPPA warm-start is the right next step

IEPPA (faithful coordinate descent in log-factor space; `src/ieppa.cpp:127` `ieppa_solve`) is *empirically robust on overlapping margins*. On the same stepstone fixture where target homotopy regressed by 8×, ieppa+sraa already converges to `max_err = 1.13e-4` (better than the LM baseline 2.8e-4 and approaching the `<1e-4` gate by itself). On kk1204 severe-skew (n=1M, K=20, nj=5, max_weight=3, OMP=1), ieppa+sraa is the published top performer at 3.7s / max_err 2.4e-14 (spec table). IEPPA's per-margin-per-category coordinate updates cannot overshoot across margin boundaries — a sweep is a sequence of one-dimensional moves, each on the linearly-separable problem `α* = log(T_kj · Z_total / S_kj)`. The dual basin is reached by construction, not by Hessian luck.

The *combined* picture: IEPPA brings us inside Newton's quadratic basin in O(K) outer sweeps (linear convergence rate); LM-Newton then polishes from inside the basin at the rate that originally motivated the Newton-KL solver. This is the textbook coordinate-descent → second-order combination for overlapping-margin KL minimisation.

### Why this is mathematically a no-op handoff (correctness argument)

- IEPPA carries `lf[cat_offset[k] + j]` for all `j = 0..cat_counts[k]-1` (no reference elimination). Each obs has `u_i_ieppa = Σ_k lf[cat_offset[k] + j_k(i)]` (positive `j` only; `j<0` is the NA bucket and is skipped — same convention as the Newton inner at `src/newton_calib.cpp:92-93`).
- Newton carries `λ[lam_off[k] + j - 1]` for `j = 1..cat_counts[k]-1` (reference category 0 fixed at `λ_{k,0} = 0`). Each obs has `u_i_newton = Σ_{k: j_k(i)>0} λ[lam_off[k] + j_k(i) - 1]`.
- **Conversion:** for each margin `k` and each `j ≥ 1`, set
  ```
  λ[lam_off[k] + j - 1] = lf[cat_offset[k] + j] - lf[cat_offset[k] + 0]
  ```
- Then `u_i_newton(λ) = u_i_ieppa(lf) - Σ_k lf[cat_offset[k] + 0] = u_i_ieppa(lf) - C` where `C` is a `λ`-independent constant.
- Newton's stable LSE form (`src/newton_calib.cpp:106-118`) computes `Z_stable = Σ d_i exp(u_i - u_max)`, `g(λ) = u_max + log(Z_stable) - T·λ`. Both `Z_stable` and `u_max` shift by the same constant under `u → u - C`; their cancellation in `log(Z_stable) - max` is exact. The recovered weights `w_i ∝ d_i · exp(u_i)` are scale-equivalent under any constant shift of `u`. Hence the converted `λ` produces *bit-identical primal weights* to IEPPA's `lf`.
- **Important caveat — NA buckets.** IEPPA stores `lf` for `j = cat_counts[k]` (the NA bucket at `cat_offset[k] + cat_counts[k]`). Newton has no NA-bucket dual. The Newton solver's `compute_u` already gates on `j > 0` (skips both `j == 0` reference *and* `j < 0` NA). Conversion uses only positive `j ≤ cat_counts[k] - 1`. The NA `lf` slot is silently dropped; this is correct because no obs with `j_k(i) = NA` enters Newton's dual sums. Confirmed against `src/ieppa.cpp:649` (`s += lf[cat_offset[m] + gm]` only fires when `gm >= 0`) and `src/newton_calib.cpp:91-93`.

### What the IEPPA solver gives us, structurally

- `ieppa_solve(CalibState&)` (`src/ieppa.cpp:127`, 1924 lines) is a single monolithic function. Its inner data — `lf` (`src/ieppa.cpp:202`, size `cat_offset[st.K]`), `cell_lf` (incremental), `X_cur`, SRAA-m state, SOR-omega state, bounds machinery — lives on its stack frame. The faithful inner loop body is the per-`(k, j)` sweep starting around `src/ieppa.cpp:637` (the `for (k) { for (j) { compute α* and update lf } }` block).
- IEPPA already exposes a *sweep budget* via `st.outer_max_iter`. We will **not** carve a new inner-loop helper out of `ieppa.cpp`. Instead, WI-1 introduces a thin orchestration shim that calls `ieppa_solve` against a duplicated `CalibState` configured with the warm-start budget, then extracts `lf` via a single new accessor. Rationale: extracting the per-sweep loop body cleanly would require touching ~50 captures and risk behaviour drift on every existing IEPPA caller. The shim is an additive change.

---

## Mechanism

### Conversion lf → λ (exact)

```cpp
// Inputs:
//   lf       : size cat_offset[K], lf[cat_offset[k] + j] for j=0..cat_counts[k]
//              (j == cat_counts[k] is the NA bucket; ignored)
//   cat_offset_ie : ieppa-side cat_offset (size K+1; +1 per margin for NA)
//   lam_off  : newton-side lam_off (size K+1; cat_counts[k]-1 free per margin)
// Output:
//   lam      : size lam_off[K], lam[lam_off[k] + j - 1] for j=1..cat_counts[k]-1
//
// Math: λ[k,j-1] = lf[k,j] - lf[k,0]. Constant column shift cancels in LSE.
for (int k = 0; k < K; ++k) {
    const double lf0 = lf[cat_offset_ie[k] + 0];
    for (int j = 1; j < cat_counts[k]; ++j)
        lam[lam_off[k] + j - 1] = lf[cat_offset_ie[k] + j] - lf0;
}
```

Verification before commit: a unit test asserts that for a tiny fixture (n=8, K=2, nj=3), running 1 IEPPA sweep + conversion produces `λ` whose `compute_u(λ, i)` differs from IEPPA's `cell_lf[c_of(i)]` by a single per-iter constant — i.e., `var(u_newton - u_ieppa) < 1e-12`. This is the structural correctness gate (T6, see §Plan Steps).

### IEPPA-iter-budget reasoning

IEPPA on stepstone-class problems converges at a *linear rate* with contraction factor empirically `≈ 0.5–0.7` per outer sweep (see `docs/investigations/2026-04-23-kk1204-convergence.md` figures). Starting from `lf = 0` (uniform), after `K_warm` sweeps the residual is `≈ ρ^K_warm × max_err_initial`. Stepstone initial `max_err` at `lf=0` is `O(1)`; for `ρ = 0.6` and target handoff residual `1e-2` (well inside Newton's basin given the LM solver's basin radius observed empirically as `gap < 1e-1`), `K_warm = log(1e-2) / log(0.6) ≈ 9`. Round to **8** as default; it provides a 5-10× margin within Newton's basin while costing one sub-second IEPPA pass on K=9 / one ~3s pass on K=20.

**Decision:** `K_warm = 8` is hard-coded at top of `newton_calibrate`, with the comment citing the contraction-rate calculation. Not exposed via `CalibState` in this epic — keep the API surface flat. The constant is `constexpr int kNewtonKLWarmstartIters = 8;`. Tunability is a follow-up if benchmarks show wall-time savings from `K_warm = 4` (early handoff) or quality gains from `K_warm = 16` (deeper warm-up).

**Adaptive early-exit (kept simple):** if IEPPA's internal convergence (`st.tol_abs`) trips inside the budget, IEPPA returns early. Newton sees `n_warmstart_iters_used < kNewtonKLWarmstartIters` and proceeds. No special-cased "skip Newton on early IEPPA convergence" path — see Q4 below.

### Q1: IEPPA refactor — minimal extraction shim

Instead of refactoring `ieppa_solve` into `ieppa_inner(CalibState&, int max_iters) -> lf`, introduce a thin **orchestration shim** in a new translation unit `src/newton_calib_warmstart.cpp` (or inline at top of `newton_calib.cpp` — implementer chooses):

```cpp
// Run a warm-start IEPPA pass and return the converted λ.
// Configures a budget-limited IEPPA call against a copied CalibState
// (so we don't mutate caller's max_iter / accelerate / convergence_cfg).
// Returns:
//   lam_out  : size lam_off[K], converted from IEPPA's final lf
//   n_iters  : number of IEPPA outer sweeps actually run
//   max_err  : IEPPA-side max_err at handoff (for diagnostics)
struct WarmstartResult {
    std::vector<double> lam;
    int n_iters;
    double max_err;
    int status;  // RK_OK or first IEPPA error code
};
WarmstartResult newton_kl_ieppa_warmstart(const CalibState& st, int K_warm);
```

The shim:
1. Copies `st` into a local `CalibState st_warm` (struct-copy is cheap; pointer fields are caller-owned and shared safely).
2. Sets `st_warm.outer_max_iter = K_warm` and `st_warm.accelerate = true` (SRAA on; see Q2).
3. Calls `ieppa_solve(st_warm)` and captures the returned `IEPPAResult.base.max_error`, `iterations`, and `status`.
4. Reaches into a *new* small accessor: an overload `ieppa_solve_capture_lf(CalibState&, std::vector<double>& lf_out)` that runs the same body but on its last `unpack_lf` / final lf state writes into the caller-supplied vector. **Or simpler**: introduce a single output parameter on `ieppa_solve` via overload — `ieppa_solve(CalibState&, std::vector<double>* lf_capture = nullptr)`. The default-nullptr overload preserves existing callers bit-identical; the warmstart shim passes a non-null pointer.
5. Performs the conversion above into `lam_out`.

**This is the minimal refactor.** It touches one new function declaration in `src/ieppa.hpp`, one new code path in `src/ieppa.cpp` (the lf-capture write — single line near each function exit), and adds the shim file. Existing callers (`harvest(method="ieppa")`, AUTO routing, `r_bridge.cpp`) are bit-identical. WI-1 ticket implements exactly this and only this.

### Q2: SRAA on or off for warm-start? — **ON**

`accelerate = true` (SRAA-m) on the IEPPA warm-start. Justification:
- SRAA-m converges 3-5× faster on stepstone (per `docs/superpowers/specs/2026-04-30-ieppa-sraa-acceleration.md`). The whole point of an 8-iter budget is to land deep in Newton's basin cheaply. Plain IEPPA at 8 iters gives roughly a 10× residual contraction; SRAA-m gives roughly 30-100×. Newton sees a tighter handoff for the same wall cost.
- SRAA-m's failure mode is super-step rejection (falls back to plain IEPPA on that sweep), not divergence. There is no scenario where SRAA-on hurts.
- The IEPPA+SRAA path is the published top performer on kk1204 (3.7s / max_err 2.4e-14). Reusing the same code path is the lowest-risk choice.
- SRAA off (plain IEPPA) is the *fallback* in WI-2's status-1 recovery path: if SRAA produces NaN or status≠0 inside the warm-start, the shim re-runs once with `accelerate = false` before bailing to "skip warm-start, run Newton from λ=0".

### Q3: Warm-start budget — **`K_warm = 8`, fixed**

See the contraction-rate derivation above. Not adaptive on gap because (a) the gap-to-Newton-basin-radius mapping isn't a clean threshold (basin depends on local `H` curvature, which Newton itself measures), (b) every dynamic-budget scheme adds branching that needs its own tests, and (c) 8 iters of IEPPA+SRAA on K=9 is `<0.5s` (per spec table extrapolation; ieppa+sraa needs 10 iters for full convergence at 3.7s on K=20 with n=1M; per-iter cost there is ~0.4s; on K=9 n=200K it is ~0.05s). The budget is dwarfed by Newton's own per-iter cost on these sizes.

### Q4: What if IEPPA already converges (max_err < 1e-4) inside K_warm?

**Decision: always run Newton afterwards (option (a)).** Justification:
1. Newton's polishing cost from a converged warm start is ~1-3 iters (basin entry → quadratic floor in 1-2 steps). Per-iter cost on stepstone is ~0.05s; total Newton tail is ~0.1s. Negligible.
2. Newton from a converged IEPPA `λ` is a *correctness check by construction*: if Newton's first-iter gradient is `<tol_abs`, it exits at iter 0 with `status = RK_OK` and the same primal weights IEPPA produced. The user gets a quality-checked Newton-KL result with no measurable cost. If somehow Newton's gradient is *not* tiny (e.g., target rounding mismatch between ieppa's per-cat residual and newton's per-(k,j-1) gradient), Newton polishes; the user still wins.
3. Option (b) "adaptive handoff: skip Newton if IEPPA converged" branches the API surface for negligible savings and creates two distinct user-facing convergence paths to test. Avoid.

### Q5: Failure modes

| Failure | Detection | Recovery |
|---|---|---|
| IEPPA SRAA produces NaN / `status = RK_ERR_*` inside warm-start | `IEPPAResult.base.status != RK_OK` from shim | Retry once with `accelerate = false` (plain IEPPA, same K_warm). If that also fails, **skip warm-start entirely**: log a one-line warning, set `n_warmstart_iters_used = 0`, run Newton from `λ = 0`. Newton's existing best-iterate-fallback handles the rest. |
| IEPPA bound infeasibility (max_weight=3 + skew, no feasible primal exists) | `IEPPAResult.base.status` returns the bounds error code | Same as above: skip warm-start, run Newton from `λ=0`. Newton itself doesn't enforce bounds (it's the dual optimiser only); the bounds-infeasibility surface is the IEPPA cell-mode invariant, not Newton's concern. The user sees a Newton result with the unmodified `max_err` characterisation. |
| Conversion produces non-finite `λ` (lf overflow during warm-start before SRAA quenches) | `std::isfinite(lam[a])` check on conversion output | Skip warm-start, log warning, Newton from `λ=0`. |
| Newton from warm-started `λ` has `status != RK_OK` after `max_iter` | Caller-side: existing post-Newton path handles this exactly as today. | No change — Newton's own NOCONV/best-iterate logic handles it. |

In all four failure modes, the behaviour reduces to "Newton from `λ = 0`", which is exactly today's master behaviour. The warm-start is *strictly additive*: it can only improve quality, never regress it (modulo the warm-start's own wall cost).

### Q6: Public API impact — none

`harvest(method="newton_kl", ...)` signature is unchanged. AUTO routing in `src/c_api.cpp:182` is unchanged. The only user-visible changes are:
1. Two new fields on `NewtonCalibResult` (`n_warmstart_iters_used`, `warmstart_max_err_at_handoff`) — additive.
2. A NEWS.md bullet (mandatory in same commit as WI-2 wiring per discipline §4).

### Q7: Field naming — rename `n_homotopy_levels_used` → `n_warmstart_iters_used`

**Decision: rename in place (option (a)).** Justification:
- The field was added by Epic-B's TH-1a (`5355f69`); Epic-B is BLOCKED and its semantics did not ship to users. The R-side surfacing (LM-1c) was *deliberately deferred* per the homotopy investigation's "Outstanding scope" section: *"`lm_mu_final` and `n_homotopy_levels_used` R-side surfacing deferred until a method actually needs it"*. So no R-side consumer exists; no `r_bridge` SEXP-pack code currently reads the field.
- Grep for `n_homotopy_levels_used` in the worktree (verified during planning) finds it only in `src/newton_calib.cpp` (zero writes — the homotopy code that wrote it was discarded), `src/newton_calib.hpp` (the field declaration), the spec doc, the homotopy plan, and the homotopy investigation. No live reader exists outside the source files this epic touches anyway.
- Option (b) "deprecate-in-place + add new field" leaves a permanently-zero stale field on a publicly-allocated struct. Worse hygiene; same risk profile.

The rename is a single-commit edit on `src/newton_calib.hpp` (rename declaration), `src/newton_calib.cpp` (rename in `NewtonCalibResult res;` initialisation paths). Done in WI-2 alongside the warm-start wiring (algorithmic-with-code colocation per discipline §4).

Add `double warmstart_max_err_at_handoff = 0.0;` as a second new field. Rationale: it gives reviewers and users a one-glance signal of "did the warm-start actually pre-converge or not?" — useful for debugging without re-running with `verbose`.

---

## Forbidden

Carried forward from rev 2 homotopy plan where applicable, plus IEPPA-specific entries. Italicised items are *not* applicable to this epic and are listed only to pre-empt review re-raising.

- **No target homotopy of any form.** Epic-B BLOCKED on this; do not reintroduce.
- **No touching the IEPPA coordinate-descent loop body** (`src/ieppa.cpp:637-690` and surrounding sweep machinery). The shim only changes IEPPA's *external* affordances (one new lf-capture parameter on `ieppa_solve`, defaulted to `nullptr` so existing callers are bit-identical).
- **No SRAA disabled by default** on the warm-start. SRAA off is the recovery path only.
- **No loosening of any test gate.** T2 stays at `<1e-4`. T8 stays at `<1e-3`. T7 (new) lands at `<1e-4`. NO relaxation to mask convergence shortfall.
- **No exposing `K_warm` via `CalibState` or `harvest()`** in this epic. Hard-coded constant. Tunability is a follow-up.
- **No "skip Newton if IEPPA converged" branch.** Always run Newton after warm-start.
- **No homotopy-bounds reuse.** `ieppa.homotopy_levels` is a different concept (bounds homotopy on `max_weight`); leave it alone. The warm-start runs at the user's `max_weight` directly.
- **No public API changes.** `harvest()` signature, AUTO routing in `src/c_api.cpp:182` and `src/r_bridge.cpp:425` — all unchanged.
- **No NEWS.md update deferred to a later ticket.** Algorithmic change colocates with code in WI-2 per discipline §4.
- **No skipping pre-commit hooks.**
- **No T2 amendment.** WI-3 explicitly forbids touching T2; only ADDs T7 + T8 (T8 is the K=20 severe-skew unit test mirrored from Epic-B's design).
- **No bundling tickets.** WI-0, WI-1, WI-2, WI-3, WI-4, WI-5a, WI-5b are each one bead. WI-2 is "warm-start wiring + NEWS.md + field rename" — atomic only because the rename is meaningless without the wiring.
- **No touching files outside:** `src/newton_calib.cpp`, `src/newton_calib.hpp`, `src/ieppa.cpp` (one new lf-capture parameter on `ieppa_solve`), `src/ieppa.hpp` (one new function signature), `tests/testthat/test-newton-kl.R`, `benchmarks/newton_kl_bench.R`, `benchmarks/results/`, `NEWS.md`, `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md` (WI-0 only). NO touches to `src/r_bridge.cpp`, `src/c_api.cpp`, `R/*.R`.
- *Not applicable: `T_eps` schedules, `lm_mu` per-ε reset semantics, `T_eps` feasibility guard, ε-schedule tail-jump risk, `n_homotopy_levels_used == N_EPS` brittle assertion.* All Epic-B concerns; this epic does no homotopy.
- *Not applicable: "rank deficiency in `H` is fixed by warm start".* It isn't. Warm start lands `λ` inside the basin; rank deficiency can still strike if the LM solver oversteps from the warm-start. The LM-Newton's own best-iterate-restoration handles that, exactly as today.

---

## Plan Steps

Each ticket is independently revertible. The sequence is strict: WI-0 → WI-1 → WI-2 → WI-3 → WI-4 → WI-5a → WI-5b. No parallel paths in this epic (the parallel `LM-1c` r_bridge surfacing from rev 2 is deferred — no R-side consumer needs it, and it is properly Epic-D scope if ever raised).

---

### WI-0 — Spec amendment

- **File:** `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md`
- **Objective:** Add an "IEPPA warm-start" subsection after the "Target homotopy" subsection (which Epic-B's TH-0 already landed). Document: rationale (rank deficiency robust to coordinate descent, basin entry by construction), `K_warm = 8` derivation, lf → λ conversion math, SRAA default, failure-fallback behaviour, the two new `NewtonCalibResult` fields, AUTO routing impact (none).
- **Constraints:** Spec amendment only. No code touched. No NEWS.md. Conventional commit `docs(spec): Newton-KL — add IEPPA warm-start section`.
- **Body sketch (skeletal — controller expands at ticket creation):**
  - Subsection header.
  - Math block: lf → λ conversion (copied verbatim from §Mechanism).
  - Justification block: linear-rate-contraction argument for K_warm = 8.
  - Field documentation: `n_warmstart_iters_used`, `warmstart_max_err_at_handoff`, deprecation of `n_homotopy_levels_used`.
  - Cross-reference to the homotopy investigation as the "why not target homotopy" explanation.

---

### WI-1 — IEPPA lf-capture orchestration shim (no Newton wiring yet)

- **Files:** `src/ieppa.hpp`, `src/ieppa.cpp`, `src/newton_calib.cpp` (new helper function — does NOT call into Newton inner), `src/newton_calib.hpp` (forward declare warmstart helper if exposed at TU level; otherwise file-static).
- **Objective:** Land the orchestration shim (`newton_kl_ieppa_warmstart`) and the IEPPA lf-capture affordance. **No behavioural change to the Newton solver in this ticket** — `run_newton_inner(T, max_iter)` is still called with `lam` initialised to zero. The shim is dead code at end of WI-1 (called by no one). Build clean. All existing tests pass bit-identically.
- **Constraints:** `ieppa_solve` signature change must default-nullptr the new lf-capture parameter so all existing callers compile and behave bit-identically. Confirm this with grep before commit: *all* call sites of `ieppa_solve` (the `harvest` path through `r_bridge.cpp`, AUTO routing through `c_api.cpp`, any internal callers) take the no-arg default. Build clean (`R CMD INSTALL --preclean .`); zero new warnings.
- **Body sketch:**
  1. Edit `src/ieppa.hpp`: change `IEPPAResult ieppa_solve(CalibState& state);` to `IEPPAResult ieppa_solve(CalibState& state, std::vector<double>* lf_capture = nullptr);`.
  2. Edit `src/ieppa.cpp`: at the function's exit paths (convergence return + recovery returns + early-exit returns), if `lf_capture != nullptr`, write `*lf_capture = lf;` (copy current lf vector).
  3. Edit `src/newton_calib.cpp`: add `static WarmstartResult newton_kl_ieppa_warmstart(const CalibState& st, int K_warm);` definition with body per §Mechanism Q1. Include the lf → λ conversion. Include the SRAA-fallback recovery logic per Q5.
  4. Add a TU-local unit-test affordance: a small inline assertion (gated on `#ifndef NDEBUG`) that verifies `var(u_newton(lam) - u_ieppa(lf)) < 1e-12` on the shim's output, run on a tiny synthetic fixture (n=8, K=2, nj=3) at TU init time. Or — cleaner — drop this in favour of a proper testthat unit test as part of WI-3 (T6). Implementer's choice; the testthat path is recommended.
  5. Conventional commit: `feat(newton_kl,ieppa): add IEPPA lf-capture shim for Newton-KL warm-start (no wiring yet)`.

---

### WI-2 — Warm-start wiring + field rename + NEWS.md

- **Files:** `src/newton_calib.cpp`, `src/newton_calib.hpp`, `NEWS.md`. Single atomic commit.
- **Objective:** Rename `n_homotopy_levels_used` → `n_warmstart_iters_used`; add `warmstart_max_err_at_handoff`; wire the warm-start shim before the Newton inner call at `src/newton_calib.cpp:371`. Update NEWS.md.
- **Constraints:** Bit-identical Newton inner behaviour for `K_warm = 0` (i.e., if someone hard-codes `kNewtonKLWarmstartIters = 0`, the path collapses to today's master). Build clean.
- **Body sketch:**
  1. `src/newton_calib.hpp`: rename `int n_homotopy_levels_used = 0;` to `int n_warmstart_iters_used = 0;`. Add `double warmstart_max_err_at_handoff = 0.0;`.
  2. `src/newton_calib.cpp`: above the lambda definition, add `constexpr int kNewtonKLWarmstartIters = 8;`. Comment cites contraction-rate derivation and pins `8` as default.
  3. Replace the single `run_newton_inner(T, max_iter);` at line 371 with:
     ```cpp
     // Warm-start via IEPPA coordinate descent (8 sweeps, SRAA on).
     // Lands λ inside Newton's quadratic basin on overlapping-margin fixtures
     // where bare LM-Newton drifts (e.g., stepstone K=9). Strictly additive:
     // failure modes fall back to λ=0 with no behaviour change. See spec.
     const auto warm = newton_kl_ieppa_warmstart(st, kNewtonKLWarmstartIters);
     if (warm.status == RK_OK && warm.lam.size() == static_cast<size_t>(n_lam)) {
         std::copy(warm.lam.begin(), warm.lam.end(), lam.begin());
         res.n_warmstart_iters_used = warm.n_iters;
         res.warmstart_max_err_at_handoff = warm.max_err;
         // Refresh the Newton inner's local g_curr / Z_curr / u_max_curr
         // so the first iter's convergence check sees the warm λ.
         g_curr = eval_dual(lam, Z_curr, u_max_curr);
     }
     // else: lam is still zero, n_warmstart_iters_used = 0; Newton runs cold.
     run_newton_inner(T, max_iter);
     ```
  4. NEWS.md (same commit): add bullet under `## Newton-KL calibration`:
     > * `method="newton_kl"` now warm-starts via 8 iterations of IEPPA coordinate descent (with SRAA-m acceleration) on the original targets before entering the LM-damped Newton inner. Lands `λ` inside Newton's quadratic basin on overlapping-margin fixtures (e.g., stepstone K=9) where bare LM-Newton drifts past the optimum. Strictly additive: warm-start failure falls back to cold start with no behaviour change. New result fields: `n_warmstart_iters_used`, `warmstart_max_err_at_handoff`. The Epic-B `n_homotopy_levels_used` field is renamed to `n_warmstart_iters_used` (no R-side consumer existed).
  5. Build clean.
  6. Conventional commit: `feat(newton_kl): IEPPA warm-start before LM-Newton inner; field rename`.

---

### WI-3 — T6 + T7 + T8 regression tests

- **Files:** `tests/testthat/test-newton-kl.R`.
- **Objective:** Three new tests. T6 = lf → λ conversion correctness (the structural gate). T7 = stepstone K=9 + warm-start audit. T8 = K=20 severe-skew unit test.
- **Constraints:** **DO NOT modify T2.** T2's `<1e-4` gate stays. The warm-start makes T2 pass at `<1e-4` from `2.8e-4` (the entire point of this epic). Only ADD T6, T7, T8.
- **Body sketch:**
  - **T6 (lf → λ conversion).** Synthetic n=64 K=3 nj=3 fixture, simple targets. Run the warm-start shim directly (skip the public R API; call via a thin testthat helper that exposes `newton_kl_ieppa_warmstart` for testing — alternative: assert via the IEPPA solver and the Newton solver independently, computing each side's `u_i` and asserting `var(u_newton - u_ieppa) < 1e-10`). Asserts: warm-start produces `λ` such that recovered weights from Newton at this `λ` match IEPPA's recovered weights to within 1e-10 relative.
  - **T7 (stepstone K=9 + audit).** Existing stepstone fixture from T2. Run `harvest(method="newton_kl")`. Assert `max_err < 1e-4`, `status == 0`, `n_warmstart_iters_used >= 1`, `warmstart_max_err_at_handoff < 1e-2` (sanity — handoff was inside basin). The assertion `n_warmstart_iters_used >= 1` is loose-and-correct (loose per Epic-B reviewer F7's earlier guidance against brittle `==N` asserts).
  - **T8 (K=20 severe-skew unit).** n=10000, K=20, nj=5, target {0.6, 0.2, 0.1, 0.07, 0.03}, `max_weight = 3`. Assert `status == 0` (or `1` if bounds violation expected at this n; implementer determines empirically) AND `max_err < 1e-3`. Catches regressions where warm-start fails to break the K=20 drift.
  - Conventional commit: `test(newton_kl): T6 conversion + T7 stepstone audit + T8 K=20 severe-skew`.

---

### WI-4 — Verify test suite

- **Files:** none (verification only).
- **Objective:** Run the full testthat suite and confirm zero failures. Specific gate: T2 (stepstone) passes at `<1e-4`.
- **Body:**
  - `Rscript -e "testthat::test_file('tests/testthat/test-newton-kl.R')"` — T1, T2, T3, T4, T5, T6, T7, T8 all PASS.
  - `Rscript -e "testthat::test_local('.')"` — FAIL=0.
  - **Halt criterion:** if T2 still fails (warm-start did not break the basin floor), STOP — escalate to BLOCKED. The warm-start either works or it doesn't; there is no rev-2-style mid-course adjustment available. (See §Risks for what BLOCKED looks like and what the next epic would be.)
  - No commit on this ticket — verification only.

---

### WI-5a — kk1204 benchmark

- **Files:** `benchmarks/results/newton_kl_kk1204.csv`.
- **Objective:** Run kk1204 severe-skew benchmark with `OMP_NUM_THREADS=1`. Capture wall, iters (Newton + warm-start separately), `max_err`, `n_warmstart_iters_used`, `warmstart_max_err_at_handoff`, `lm_mu_final`, status.
- **Body:**
  - `OMP_NUM_THREADS=1 Rscript benchmarks/newton_kl_bench.R` (existing script; if it doesn't capture warm-start fields, edit it in this ticket to do so — same-commit).
  - Save CSV.
  - Conventional commit: `bench(newton_kl): kk1204 with IEPPA warm-start`.

---

### WI-5b — Verdict + investigation doc + close epic

- **Files:** `docs/investigations/2026-05-01-newton-kl-ieppa-warmstart-result.md` (new); possibly `NEWS.md` (verdict bullet only).
- **Objective:** Write the result document. Decide GATE_MET / PARTIAL / BLOCKED. Close epic appropriately.
- **Decision rules:**
  - **GATE_MET** (kk1204 severe-skew wall `<3s` ∧ `max_err <1e-4` ∧ stepstone `<1e-4`): commit investigation doc + verdict NEWS.md bullet, close Epic-C as resolved. *Honest assessment: GATE_MET is plausible on quality (warm-start solves the structural problem) but unlikely on wall (see §Cost). Expected outcome is PARTIAL.*
  - **PARTIAL** (quality OK on T2/T7/T8 but wall ≥3s on kk1204, OR severe-skew `max_err ∈ [1e-4, 1e-3]`): commit `fix(newton_kl): IEPPA warm-start partial — wall budget unmet` analogue. **Close Epic-C** with the warm-start as the *correct* algorithm; file Epic-D ("Newton-KL wall-time tuning" or "AUTO routing on quality verdict") as a follow-up. PARTIAL is the planned successful closure mode.
  - **BLOCKED** (T2 stays `>1e-4` even with warm-start): the basin-warm-start hypothesis itself is wrong. The investigation doc explains why; epic closes BLOCKED. The follow-up at this point is Epic-E: ESCAPE-only AUTO routing (route severe-skew K≥5 to `ieppa+sraa` directly, leaving Newton-KL for moderate-skew K≥5). This is the rev 2 homotopy plan's contingency, copy-pasted forward.

---

## Cost Estimate

Honest, source-derived. Numbers cited from the spec table at `docs/superpowers/specs/2026-05-01-newton-kl-calibration-design.md:9-21` and the LM rev 2 measurements.

| Component | kk1204 severe-skew (K=20, n=1M) | stepstone (K=9, n=200K) |
|---|---|---|
| IEPPA+SRAA warm-start (8 sweeps) | ≈ 8/10 × 3.7s = **3.0s** | ≈ 8/10 × ~0.3s = **0.24s** |
| LM-Newton tail (1-3 iters from inside basin) | ≈ 1-3 × 0.6s = **0.6 — 1.8s** | ≈ 1-3 × 0.05s = **0.05 — 0.15s** |
| **Total** | **3.6 — 4.8s** | **≈ 0.3 — 0.4s** |

**Honest gate:** kk1204 < 3s. **Realistic outcome:** PARTIAL on wall (3.6-4.8s vs 3s gate). The warm-start is paying for IEPPA's already-3.7s cost and adding Newton on top.

The original Newton-KL spec gate (<2s on kk1204) is dead under any approach that uses IEPPA as a warm-start. The path to <2s would require either (a) IEPPA with fewer sweeps that still land Newton in the basin (need to test K_warm = 4 in a follow-up), or (b) a different first-order pre-conditioner (out of scope; the only one we have that works on these fixtures is IEPPA).

**Quality gate:** stepstone `<1e-4`. **Realistic outcome:** GATE_MET. IEPPA+SRAA already hits 1.13e-4 on stepstone alone; warm-start to Newton from there is mathematically a polish step. Expected stepstone result: `max_err < 1e-7` (Newton converges to its tol_abs from inside the basin).

**Quality gate (T8):** K=20 severe-skew unit, `<1e-3`. **Realistic outcome:** GATE_MET. The K=20 severe-skew failure under bare Newton is a basin-escape problem; warm-start prevents it.

---

## Risks & Open Questions

- **Risk: warm-start lands `λ` inside the basin but Newton's first step still drifts.** Mitigation: existing best-iterate-restoration in `run_newton_inner` (pre-existing TH-1a code at `src/newton_calib.cpp:357-365`) catches this. Worst case: Newton returns its warm-started `λ` unchanged, with quality identical to IEPPA-alone (1.13e-4 on stepstone) — still better than master's 2.8e-4.
- **Risk: IEPPA+SRAA fails on a problem where Newton-KL would have worked.** Mitigation: failure-fallback recovery (Q5) reverts to `λ = 0` and runs Newton as today. Strictly additive.
- **Risk: lf-capture parameter on `ieppa_solve` lands a subtle copy-cost regression.** Mitigation: `lf_capture` defaults to nullptr; the cost is one if-check at function exit. Confirm via WI-1's existing-tests-pass-bit-identically gate — if any IEPPA bench regresses by even 1%, dig into it.
- **Risk: `kNewtonKLWarmstartIters = 8` is wrong.** Mitigation: the constant is derived from a contraction-rate calculation, not pulled from a hat. WI-5a's benchmark CSV captures `n_warmstart_iters_used`; if early-exit fires at 3-4 sweeps consistently, future Epic-D can lower the budget. If post-warm-start Newton consistently runs >5 iters (bad warm-up), Epic-D raises it. The follow-up tuning is small-scope.
- **Open: should T6 use the public R API or a TU-local C++ test harness?** Recommended path: testthat-only via the public API (run IEPPA via `harvest(method="ieppa")` + Newton via `harvest(method="newton_kl")` and compare recovered weights). This is the cleanest scope but doesn't directly probe the lf → λ conversion. Implementer chooses; if testthat-only is too coarse, a tiny standalone C++ unit-test file `tests/testcpp/test_newton_kl_warmstart.cpp` would be the alternative.
- **Open: AUTO routing impact.** This epic does NOT change AUTO routing. If WI-5b lands GATE_MET, AUTO's `K ≥ 5 ∧ M_cell/n ≥ 0.9 → newton_kl` rule becomes more aggressive in practice (newton_kl is now reliable); no code change needed. If WI-5b lands PARTIAL on wall, AUTO routing logic stays as-is (newton_kl is still correct, just not faster than ieppa+sraa on K=20). If BLOCKED, the rev 2 homotopy plan's ESCAPE patch (severe-skew gate routing K=20 to ieppa+sraa) becomes the next epic.

---

## Reviewer Concerns from Epic-B That DO NOT Apply

Pre-empting the plan-review-gate (rev 2 homotopy reviewers caught these; they should not re-raise on Epic-C):

| Concern | Why not applicable here |
|---|---|
| C1 (mechanism handwaving — claims without proof) | §Mechanism includes the exact lf → λ conversion, the algebraic argument that recovered weights are bit-identical (constant column shift cancels in LSE), and the contraction-rate justification for `K_warm = 8`. No "this should work because basin" hand-waving — IEPPA is *empirically known* to converge on overlapping margins (1.13e-4 on stepstone today). |
| C3 (T2 gate loosening) | §Forbidden explicitly forbids it; WI-3 explicitly forbids touching T2. |
| C4 (iter-accumulation contract) | Not relevant — we don't accumulate iter counts across an outer loop. The Newton inner runs once with its full budget; the warm-start uses IEPPA's own counter (one number). |
| C7 (NEWS.md update timing) | WI-2 colocates NEWS.md update with the wiring commit per discipline §4. |
| C10 (per-method internal homotopy) | We don't use per-method homotopy at all. SRAA is acceleration (different concept). `ieppa.homotopy_levels` is bounds-homotopy (different concept again). Both are inert from this epic's perspective. |
| C12 (stagnation criterion plan/ticket disagreement) | No stagnation criterion exists in this plan. Newton uses its existing `tol_abs + max_iter` + best-iterate-restoration. |
| F2 (ε schedule tail jump) | No ε schedule. |
| F3 (`lm_mu` per-ε reset) | No per-ε anything. `lm_mu` is initialised to `1e-3` exactly as today (no change). |
| F4 (T_eps feasibility guard) | No T_eps. |
| F7 (`n_homotopy_levels_used == N` brittle assertion) | T7 uses `>= 1`, not `== K_warm`, exactly because of this concern. |
| F8 (target-skew gate in c_api.cpp) | This epic doesn't touch `c_api.cpp`. AUTO-routing gate is deferred to Epic-D/E if WI-5b lands BLOCKED. |
| S6 (`c_api.cpp:182` ternary, not fall-through) | Out of scope; `c_api.cpp` is in §Forbidden's untouched-files list. |

---

## Out of Scope

- AUTO routing changes (deferred to Epic-D/E based on WI-5b verdict).
- R-side surfacing of `n_warmstart_iters_used`, `warmstart_max_err_at_handoff`, `lm_mu_final` via `r_bridge.cpp` (the prior LM-1c parallel ticket). Filed as Epic-D scope only if WI-5b GATE_MET / PARTIAL and a user-facing tuning loop ever needs it.
- `K_warm` parameter exposed via `CalibState` or `harvest()`. Tunability is a follow-up.
- SIMD / multithreading optimisation of either IEPPA's sweep or Newton's accumulation.
- Different first-order pre-conditioners (raking-warmstart, sinkhorn-warmstart, etc.). IEPPA is chosen because (a) it works empirically on overlapping margins, (b) it has SRAA acceleration already, (c) it has the same `CalibState` ABI as Newton.
- Different K_warm values per fixture / dynamic adjustment based on K. Hard-coded constant; revisit only on bench-driven evidence in a follow-up.
- Touching the IEPPA coordinate-descent loop body for any reason (cosmetic, performance, or otherwise). Off-limits in this epic per §Forbidden.
