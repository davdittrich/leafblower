## Conflict Detection Report

Ingest set: 44 classified docs (32 SPEC, 11 DOC, 1 PRD; 0 ADR, 0 UNKNOWN; all
`confidence: high`, all `manifest_override: true`, none `locked`). Mode: new.
Precedence applied: SPEC > PRD > DOC, then document date (later supersedes earlier on the
same subject), per the ingest brief. Cross-ref cycle detection: run, no cycles found.

Third run. This replaces the previous report. The PRD that run 2 recorded as missing has
been located at `tasks/prd-leafblower-core.md` — it sits outside the
`docs/superpowers/specs/` and `docs/methods/` convention, which is why directory-based
discovery missed it. `requirements.md` is no longer empty by construction: nine
requirement entries are extracted from it.

Two further run-2 warnings were verified against live source by the orchestrator between
runs and are closed: the `omega_mode_id` R/Python wiring and the `kNCatsTotalMax` value.
Both are recorded below as INFO with the verified values, and must not be re-raised.

### BLOCKERS (0)

None. There is no ADR-class document and no `locked` document in the set, so
LOCKED-vs-LOCKED contradiction is structurally impossible here. The corpus is safe to
route once the single warning is acknowledged.

### WARNINGS (1)

[WARNING] docs/raking.md §8.2 attributes the unimplemented iEPPA outer proximal-point
  algorithm to ORIS
  Found: docs/raking.md §8.2 "The ORIS Solver (renamed from iEPPA)" describes ORIS as
    wrapping "a proximal point framework around an inner solver", executing an "Entropic
    Proximal Step", a "Dual Block Coordinate Descent Subsolver" and "Inexact Stopping
    Criteria", and claims it "outperforms both stabilized Dykstra's algorithms and premium
    commercial linear programming solvers like Gurobi". §12 repeats the attribution.
  Found: docs/methods/oris.md states the exact opposite as the reason the rename happened:
    the paper's "headline contribution — an outer inexact entropic-proximal-point loop …
    is NOT implemented here", it is "mathematically inert" at `C = 0`, and citing
    arXiv:2011.14312 "over-claims a contribution the code does not contain".
  Found: docs/superpowers/specs/2026-05-30-oris-rename-design.md §1 gives the same
    rationale and §8.10 lists `docs/raking.md` explicitly as a LIVE doc in rename scope.
  Found: tasks/prd-leafblower-core.md § US-005 is where the "paper-faithful iEPPA
    (Chu-Liang-Toh-Yang 2022, arXiv:2011.14312) at C=0" framing originates (FR-11, FR-12
    and the FR-14 Bregman/`tolRb` inner stop encode the outer loop). The PRD is the oldest
    document in the corpus and is superseded on this point; docs/raking.md is not.
  Impact: The rename appears to have been applied to this file as a name substitution that
    preserved the very claim the rename exists to repudiate. Precedence resolves the truth
    cleanly (SPEC > DOC; and DOC-vs-DOC broken by the rename spec), but the false
    capability claim is still sitting in a live document that a reader or a downstream
    generator will pick up.
  → Rewrite docs/raking.md §8.2 and §12 to describe the paper's algorithm under its own
    name, or delete the passage; do not leave an ORIS-labelled description of an
    unimplemented outer loop. Already ticketed as `leafblower-05ha`.

### INFO (25)

[INFO] Resolved (was WARNING): the referenced PRD has been located
  Note: run 2 recorded that `2026-04-23-ieppa-faithful-design.md` §3.1/§4.3/§10/§11 cite a
    PRD (§ US-001 AC, § US-005, new § US-005b, "FR entries updated to match", "Non-Goals:
    as in current PRD plus …") that was absent from the ingest set. The document is
    `tasks/prd-leafblower-core.md` — outside the `docs/` convention that discovery walked.
    It is now classified (`type: PRD`, `confidence: high`) and extracted into
    `intel/requirements.md` as nine entries. The requirement layer the SPEC corpus assumes
    as its base is no longer missing.

[INFO] The PRD is heavily superseded on the solver roster — record, do not resolve
  Found: tasks/prd-leafblower-core.md is "Draft v3", dated 2026-04-18 — the OLDEST
    document in a corpus running to 2026-08-14 — and describes a TWO-algorithm package:
    § 6 declares `typedef enum { RK_ALG_AUTO = 0, RK_ALG_IEPPA = 1, RK_ALG_LBFGSB = 2 }
    rk_algorithm_t`; US-006 and FR-20 … FR-28 specify L-BFGS-B in full (link selection,
    dual objective, 2-loop recursion, Wolfe line search); US-007 auto-routes between
    `RK_ALG_IEPPA` and `RK_ALG_LBFGSB` on `complexity > 500000 || max_weight < 3.0 ||
    min_weight > 0`; US-001 maps `method="rake"` and `method="nr"` to `"lbfgsb"` with a
    named warning; FR-34 maps anesrake's `choosemethod` the same way; § 11 KPIs include an
    L-BFGS-B row (100K rows, 5 margins < 2 s).
  Note: All of the above is WITHDRAWN or SUPERSEDED, and none of it needs resolution —
    precedence puts SPEC above PRD and document date reinforces it.
    `docs/superpowers/specs/2026-08-14-removed-solver-slots-supersession.md` records
    `RK_ALG_LBFGSB = 2` as removed and the slot as permanently reserved
    (`src/leafblower.h:44`); the package now ships eight solvers (AUTO 0, ORIS 1, RAKING 3,
    SINKHORN 4, CHEBYSHEV 5, GREG 6, ORIS_SOFT 8, GREENKHORN 9, LOGIT 10, NEWTON_KL 11);
    routing is the three-way rule in
    `2026-05-01-newton-kl-calibration-design.md` (K >= 5, `M_cell/n >= 0.9`, `target_skew`
    vs 5.0). This PRD is the ORIGIN of the slot-2 hole. The autumn `"rake"`/`"nr"` synonyms
    still need a destination, but the PRD no longer supplies one. Recorded in
    `intel/requirements.md` § Superseded / withdrawn items 1–4. The Deville-Sarndal logit
    DISTANCE survives in the separate `logit` solver (`RK_ALG_LOGIT = 10`,
    `src/logit_calib.cpp`) — a different solver, not a survival of US-006.

[INFO] PRD § 5 Non-Goals: four disclaimed items have since shipped or been worked
  Found: tasks/prd-leafblower-core.md § 5 disclaims for v1 — "SQUAREM acceleration — not
    ported from autumn"; "Bounded IPF (water-filling) as a named `method` — not exposed";
    "`auto_collapse` — sparse level merging not in v1 (raises informative error)";
    "`add_na_proportion` — NA-as-category not in v1 (raises informative error)".
  Note: SQUAREM SqS3 shipped via `2026-04-28-raking-squarem-design.md` (and has since been
    superseded in turn by SRAA-m); water-filling is the live ORIS `bounds_mode="unit"`
    per-cell path (`docs/methods/00-overview.md`, `2026-04-24-ieppa-speed-convergence-
    bounds-design.md`); later work exists on `auto_collapse` and `add_na_proportion`.
    These are superseded SCOPE, not violated requirements — a non-goal that was later
    delivered is a scope expansion, and the US-001 acceptance criterion that
    `auto_collapse`/`add_na_proportion` "raise an informative error" is stale with it.
    Two § 5 non-goals are NOT superseded and remain in force: CPU-only (no GPU) and
    single-threaded (the FR-9/§ 7 reentrancy contract is a caller-side guarantee, not
    internal parallelism, and single-thread BLAS is required for deterministic R<->Python
    parity at rtol=1e-6).

[INFO] PRD § 7 Build System `-O3` flag is superseded
  Found: tasks/prd-leafblower-core.md § 7 specifies
    `src/Makevars.in: PKG_CXXFLAGS = @CXXFLAGS_STD@ -I. -O3 -DSTRICT_R_HEADERS` and a
    `PKG_SOURCES` list of all `.cpp` files.
  Note: The R build now deliberately sets NO `-O` level of its own — CRAN's
    `tools:::.check_make_vars` rejects `-O` flags in `PKG_CXXFLAGS`, and R supplies the
    user/site level through `$(CXXFLAGS)` in `$(ALL_CXXFLAGS)`; a user wanting `-O3` sets
    it in `~/.R/Makevars`. Only `python/CMakeLists.txt` sets `-O3`. The `PKG_SOURCES` list
    is decorative — R auto-globs `src/*.cpp` — whereas the Python build does NOT glob and
    needs each new `src/*.cpp` added to its explicit `CORE_SOURCES` list, which is the
    inverse of the mechanism FR-38 describes (though the outcome, no `r_bridge.cpp` in the
    Python build, is the same).

[INFO] PRD US-005 is partially superseded: the paper-faithfulness framing, not the solver
  Found: tasks/prd-leafblower-core.md § US-005 requires "the paper-faithful iEPPA
    algorithm (Chu-Liang-Toh-Yang 2022, arXiv:2011.14312) at C=0", with FR-11/FR-12
    specifying the outer entropic-proximal loop (`epsilon` fixed at 0.05, outer proximal
    center `w^0`, `outer_max_iter = 50`) and FR-14 specifying the inner Bregman stopping
    rule `D(w~, w^k)/(1 + normU) < tolRb` with `normU = max_weight`, `tolRp` adaptive and
    `tolRb = 1/outer_iter^1.1`.
  Note: That outer loop is not implemented and is mathematically inert at `C = 0` — the
    stated reason for the `ieppa` -> `oris` rename. FR-11/FR-12/FR-14 are historical. The
    durable half of US-005 — cell compression with `M_cell <= min(n, prod cat_counts)`,
    log-space Sinkhorn factors with the 700 clip, the capacity clamp
    `X[c] = clamp(X_tilde[c]*W[c], L_c, U_c)`, `errRp < tol_abs` convergence, the
    `S_j < 1e-15*W` empty-category test (never exact `== 0`), NA/OOV skip, and the
    every-inner-step bounds INVARIANT — is corroborated by `docs/methods/oris.md` and
    `docs/methods/00-overview.md` and is extracted as
    `REQ-us005-oris-capacity-constrained-solver`.

[INFO] PRD US-005b projection geometry is superseded; the routing requirement is not
  Found: tasks/prd-leafblower-core.md § US-005b specifies `method="raking"` as classical
    cyclic IPF (Deming-Stephan 1940; Csiszar 1975) with ADDITIVE Dykstra box projection
    (Boyle-Dykstra 1986) plus hyperplane projection.
  Note: `docs/superpowers/specs/2026-04-27-raking-bregman-dykstra-design.md` replaces the
    Euclidean corrections with multiplicative KL-Bregman projections at cell level, with
    `Xc <= 0` and `scale <= 0` guards. `docs/raking_bounds.md` (DOC) independently
    diagnoses the metric mixing that motivated the change. The geometry moves; the
    requirement that `RK_ALG_RAKING` exist, that `method="raking"` route to
    `raking_solve`, and that the pre-rev2 iEPPA tests act as a regression guard, stands
    unchanged — and it is exactly the "§ US-005b" the SPEC corpus cites.

[INFO] PRD states the medium-scale performance target twice with different numbers
  Found: tasks/prd-leafblower-core.md § 1 Success Criteria: "Performance — medium |
    100K rows, 5 margins (3–5 cats each) < 1 s".
  Found: the same document's § 11 KPI table: "L-BFGS-B convergence | 100K rows, 5 margins
    < 2 s, `max_error < 1e-6` | `test-lbfgsb.R` Phase 1 gate".
  Note: Not treated as a competing acceptance variant, because both statements were
    written against the removed `lbfgsb` solver and neither has a live measuring artefact
    (`test-lbfgsb.R` is void). Recorded in `REQ-us003-large-scale-performance` as a target
    that must be re-benchmarked against a live solver before being used as a gate. The
    large-scale target (1M rows, 20 margins, `max_weight = 3`, < 30 s) is stated
    consistently in § 1, § 3, § 10 Phase 2 and § 11 and carries forward intact.

[INFO] Resolved (was WARNING): `omega_mode_id` R/Python wiring is fixed
  Note: run 2 warned that `2026-06-02-oris-iterate-change-omega-design.md` §3.7 recorded a
    silent R<->Python divergence — C default 2, R `parse_sor` defaulting to 1L, and the
    Python `_parse_sor` returning a 6-tuple that DROPPED `omega_mode_id` so `_bindings.cpp`
    never forwarded it (tickets e18t.9 / e18t.10). Verified against live source by the
    orchestrator between runs: all three default sites now agree at 2 —
    `src/types.hpp:74`, `src/c_api.cpp:213`, `R/harvest.R:1077` — and Python forwards it
    via `python/leafblower/_harvest.py:593` and `_bindings.cpp:140-142`. Closed; do not
    re-raise.

[INFO] Resolved (was WARNING): `kNCatsTotalMax` is 2048
  Note: run 2 warned that `2026-04-25-calibration-solvers-design.md` gives the constant as
    2048 in its § 2 code block and 8192 in its § 12 Resolved Design Questions row — a 4x
    disagreement on a user-visible input limit. Verified against live source:
    `src/calib_validate.hpp:10` defines it as **2048**. The § 12 row is the losing
    statement and is a spec-internal transcription error, not a live design question.
    Closed; do not re-raise.

[INFO] Resolved (was BLOCKER): `grake` / slot 7 is removed, not live
  Note: Four specs through 2026-04-29 treat `grake` as a shipped method —
    `2026-04-25-calibration-solvers-design.md` (§3 `#define RK_ALG_GRAKE 7`, Goal 4
    `method='grake'`, §7 algorithm, §11 acceptance criterion A4, files
    `src/grake.cpp/.hpp`), `2026-04-27-sinkhorn-a1-fix-ieppa-admm-method.md`
    (`case RK_ALG_GRAKE:` in `select_solver_objective`),
    `2026-04-29-greenkhorn-solver.md` (`alg_names` pins `"grake"` at index 7) and
    `2026-04-29-chebyshev-greg-fix.md` (Out of Scope: "grake separate fix").
    `docs/superpowers/specs/2026-08-14-removed-solver-slots-supersession.md` supersedes all
    four: `RK_ALG_GRAKE = 7` was removed pre-release by commit `9a67891` (2026-04-30,
    "no release, no ABI constraint"); `src/grake.cpp`/`src/grake.hpp` do not exist; slot 7
    is permanently reserved; NO grake acceptance criterion (A4 included) carries forward.
    The live `grake_norm` METRIC (`src/greg.cpp:153`, `src/logit_calib.cpp:558`, packed at
    `src/c_api.cpp:74/:110/:426`) is a different object and is unaffected — a
    `grep grake` in the tree returns the metric, not the solver.

[INFO] Resolved (was WARNING): the `lbfgsb` removal is on record
  Note: Four specs target L-BFGS-B — `2026-04-18-bounded-convergence-fix-design.md`
    ("Fix 2", link selection `exponential = !std::isfinite(U)`, the `use_logit` guard),
    `2026-04-20-algo-selection-design.md` (its entire subject is the L-BFGS-B vs iEPPA
    routing contour), `2026-04-24-convergence-metrics-sor-design.md` §6/§8 (edits to
    `src/lbfgsb_solver.{cpp,hpp}`), and `2026-04-29-ieppa-alm-soft-capacity.md`
    (`alm_lambda`/`alm_mu` protected as "NOT dead (used by lbfgsb)") — and, now,
    `tasks/prd-leafblower-core.md` US-006 / FR-20 … FR-28, which is where the solver was
    specified in the first place. The 2026-08-14 record confirms removal: slot 2 annotated
    at `src/leafblower.h:44` (`/* 2 = removed (was RK_ALG_LBFGSB) */`), doc drift purged in
    commit `7fa211c`, residual header note at `src/leafblower.h:65`. One live consequence
    is carried into `constraints.md` rather than left implicit: with lbfgsb gone, the
    sum-to-n ALM that justified reserving `alm_lambda`/`alm_mu` no longer exists, so those
    two fields must be re-checked against current source before reuse or deletion.

[INFO] Resolved (was WARNING): Epic-K `method="cp"` is a withdrawn proposal, not pending
  work
  Note: `2026-05-02-epic-k-cp-productionization-design.md` (rev 6) specifies
    `harvest(method="cp")` with `RK_ALG_CP = 12`, `src/cp_calib.{hpp,cpp}`, a 47-element
    result list, CP-only `accelerate=TRUE`, and an 8-test suite, decomposed into work units
    K-1…K-7. The 2026-08-14 record traces the port to commit `00a3f10` (2026-05-03) and its
    revert to `3fac1d6` the same day; `grep -c 'RK_ALG_CP' src/leafblower.h` returns 0 and
    slot 12 was never occupied. The spec is retained in `constraints.md` with its three
    entries re-titled "(WITHDRAWN PROPOSAL — never landed)"; its T1…T8 gates do not carry
    forward. This is a scope decision, not a defect: if Epic-K is deliberately revived,
    slot 12 becomes available and the supersession record must be amended.

[INFO] Resolved (was WARNING): `RK_ALG_RAKING` is 3, not 2
  Note: `2026-04-27-sinkhorn-a1-fix-ieppa-admm-method.md` Part 1 gives
    "(RK_ALG_IEPPA=1, RK_ALG_RAKING=2, RK_ALG_SINKHORN=4, etc.)" — later in date than the
    2026-04-23 definition, so pure date resolution would have adopted the wrong value.
    `2026-08-14-removed-solver-slots-supersession.md` states explicitly: "one spec in the
    corpus gives `RK_ALG_RAKING` as 2. That is wrong — raking is 3, and 2 is the reserved
    `lbfgsb` hole." This agrees with `2026-04-23-ieppa-faithful-design.md` §4.2,
    `docs/methods/raking.md` and `docs/methods/00-overview.md`. The 04-27 parenthetical is a
    transcription error; the authoritative enum is recorded in `constraints.md`.

[INFO] Exactly one PRD, zero ADRs — there is no locked-decision layer
  Note: Of 44 docs, 32 are SPEC, 11 are DOC and 1 is PRD; none is `locked` and none
    carries a per-doc `precedence` override. LOCKED-vs-LOCKED contradiction is therefore
    structurally impossible in this set, and no competing PRD acceptance variants exist
    because there is only one PRD. `decisions.md` is empty of decision entries by
    construction and says so rather than being back-filled from lower-precedence material.
    Two documents act as decision records without being ADRs: the 2026-08-14 supersession
    record (Status: Accepted, commit-traced, self-declared authoritative) and the PRD's
    § 12 Open Questions table (OQ-1 … OQ-5, three of which still stand — fixed
    `epsilon = 0.05`, `log_fn` callback for verbose output, Windows CRAN binary as a
    non-goal).

[INFO] Cross-ref cycle detection: clean
  Note: The directed graph over the classifications' `cross_refs` fields is acyclic with
    maximum depth 2. The only doc->doc edges are
    `docs/methods/sinkhorn.md -> docs/methods/oris.md`, the rename spec's fan-out
    (`2026-05-30-oris-rename-design.md` -> `2026-04-23-ieppa-faithful-design.md`,
    `docs/methods/oris.md`, `docs/methods/00-overview.md`) and the supersession record's
    fan-out (`2026-08-14-removed-solver-slots-supersession.md` -> the four grake specs, the
    Epic-K spec and `docs/methods/00-overview.md`). The newly added PRD contributes two
    outbound edges (`README.md`, `DESCRIPTION`), both to non-documents, and is referenced
    by `2026-04-23-ieppa-faithful-design.md` in prose only — its `cross_refs` array does
    not name it, so it adds no edge and no cycle. Every other `cross_refs` entry points at
    source files or at documents outside the ingest set (`docs/investigations/**`,
    `docs/superpowers/plans/**`, `docs/superpowers/derivations/oris_structure_map.md`,
    `research/**`). The method docs contain mutual markdown links in their bodies
    (oris.md <-> sinkhorn.md <-> 00-overview.md), but these are navigational, not
    derivation dependencies, and create no synthesis loop.

[INFO] Auto-resolved: AUTO routing supersession chain, later spec wins
  Note: tasks/prd-leafblower-core.md § US-007 (`complexity = n * sum(cat_counts)` vs
    500000, plus `max_weight < 3.0` or `min_weight > 0`, choosing between iEPPA and
    L-BFGS-B) -> 2026-04-20-algo-selection-design.md (route to L-BFGS-B inside a
    benchmarked 1.2x iso-contour on complexity x tolerance) ->
    2026-04-23-ieppa-faithful-design.md §4.2 ("AUTO returns RK_ALG_IEPPA always, no
    heuristic exceptions") -> 2026-05-01-newton-kl-calibration-design.md (three-way rule on
    K >= 5, `M_cell/n >= 0.9` and `target_skew` vs 5.0, implemented identically in
    `c_api.cpp` and `r_bridge.cpp`). The 2026-05-01 rule is recorded in `constraints.md` as
    current; the three earlier rules are marked superseded. The codebase map's
    ARCHITECTURE.md describes the 2026-05-01 rule, corroborating the resolution. The first
    two links are doubly dead: their routing target (L-BFGS-B) was itself removed.

[INFO] Auto-resolved: convergence-configuration supersession chain
  Note: 2026-04-24-convergence-metrics-sor-design.md (default `pct = 0.001`, single
    five-value `criterion` enum, `pct_change` result field) ->
    2026-04-25-convergence-redesign.md (orthogonal `metric` + `rule`, default
    `max_err + improvement + 0.001`, `pct` REDEFINED to `l1_weight + plateau`,
    `pct_change` renamed `l1_weight_change`) -> 2026-04-25-calibration-solvers-design.md §9
    (ieppa default metric -> `kl`) -> 2026-04-27-sinkhorn-a1-fix (sinkhorn default metric
    -> `kl`; `convergence_objective` -> `convergence_solver_objective`, now reporting the
    solver's mathematical objective) -> 2026-04-28-convergence-status-design.md
    (`RK_ERR_BUDGET = 4`, `RK_ERR_STALL = 5`, weight-KL stall for the flat loop). Each link
    is dated and non-contradictory once ordered; only the final state is carried forward.
    The PRD's four base codes (`RK_OK=0`, `RK_ERR_NOCONV=1`, `RK_ERR_INFEAS=2`,
    `RK_ERR_BADARG=3`) are EXTENDED by this chain, not replaced.

[INFO] Auto-resolved: acceleration supersession chain (six superseded mechanisms)
  Note: The longest chain in the corpus, all on the same subject — accelerating the
    ORIS/raking/greenkhorn fixed point.
    (1) APVA / joint Anderson (2026-04-24-ieppa-speed-convergence-bounds-design.md §5.2)
    — empirically failed: 487 of 500 steps rejected on kk1204 by capacity-clamp
    non-smoothness.
    (2) Tang 2024 primal-dual with Newton on constraint duals (§5.3) — SHELVED, marked
    "DO NOT IMPLEMENT as-is" after a unanimous review failure (K-count mismatch, Delta not
    runtime-knowable, scope overrun).
    (3) Halpern mixing (§5.4) — promoted to primary, then absent from all later specs.
    (4) SQUAREM SqS3 (2026-04-28-raking-squarem-design.md) plus its obs-level-alpha
    geometry fix (2026-04-28-squarem-geometry-fix.md).
    (5) SRAA-m replacing SQUAREM (2026-04-29-i0am-sraa-acceleration.md), then the global
    safeguard + revert-to-best (…-global-safeguard.md), then the adaptive-sort +
    outer-revert fix (…-correct-all-scales.md).
    (6) SRAA-m in log-factor space for ORIS/oris_soft
    (2026-04-30-ieppa-sraa-acceleration.md).
    Current state = (5)+(6). Superseded entries are retained in `constraints.md` and
    marked, because several carry measured negative results that a future plan must not
    re-discover. Note that SQUAREM at step (4) is itself the item the PRD § 5 disclaimed as
    a non-goal.

[INFO] Auto-resolved: two same-date SRAA specs give different root causes but compose
  Note: 2026-04-29-i0am-sraa-global-safeguard.md attributes the K=9 basin escape
    (2.12e-3 vs plain 1.57e-3) to a too-loose LOCAL safeguard and fixes it with a global
    quality floor; 2026-04-29-i0am-sraa-correct-all-scales.md attributes the SAME
    regression to fixed-vs-adaptive sort creating two distinct fixed points and states
    that "no per-step safeguard can distinguish the two basins in early iterations". The
    ordering is inferable only from content: the latter adds its constants "after
    `kSRAARestartGamma`" and declares plateau gating unnecessary "with adaptive sort",
    i.e. it builds on the former. Both mechanisms coexist in `sraa.hpp`; no revocation is
    implied. Recorded so the contradictory diagnoses are not mistaken for a live dispute.

[INFO] Auto-resolved: ORIS optimal-omega supersession chain (two recorded negatives)
  Note: mj1p.2 global spectral theta2 (NO-GO — global residual ratio -> 1 on pinned cells,
    omega -> 1.9, slower than fixed omega) ->
    2026-05-31-oris-free-subspace-omega-design.md (free-subspace MARGINAL-RESIDUAL theta2;
    NO-SHIP — 500-iter budget vs 140 fixed on bounded stepstone, because `R2_free` is
    algebraically identical to the global residual once clamped mass folds into the free
    target) -> 2026-06-02-oris-iterate-change-omega-design.md (free-coordinate
    ITERATE-CHANGE theta2 with block-root cadence recovery; feasibility-agnostic).
    `docs/methods/oris.md` documents the last as shipped, and the wiring is confirmed live
    (see the `omega_mode_id` INFO entry). The two negatives are preserved deliberately —
    both failure mechanisms are re-derivable traps.

[INFO] Auto-resolved: raking Dykstra geometry, SPEC over DOC
  Note: docs/raking_bounds.md (DOC) diagnoses the Euclidean/multiplicative metric mixing
    and prescribes Bregman-Dykstra with obs-level code (`std::vector<double> q(st.n)`,
    `w[i]`); docs/superpowers/specs/2026-04-27-raking-bregman-dykstra-design.md (SPEC,
    higher precedence) implements the identical mathematics at CELL level (`p[c]` over
    `ct.M_cell`, `L_cell[c]`/`U_cell[c]`) with the additional `Xc <= 0` and `scale <= 0`
    guards. SPEC wins on the level; the DOC is retained as the rationale record.

[INFO] Not a conflict: SOR direction differs by solver, by design
  Note: `2026-04-24-convergence-metrics-sor-design.md` §3 specifies SOR as adaptive
    UNDER-relaxation (omega starts at 1.0, adapts downward, floor `omega_min = 0.3`) to
    damp oscillation, and `2026-04-27-raking-bregman-dykstra-design.md` wires raking with
    `eff_omega <= 1`; docs/methods/oris.md specifies OVER-relaxation (omega > 1, ceiling
    1.8). These coexist: under-relaxation is a damping response to bound-induced
    oscillation, over-relaxation is a rate accelerator on the free subspace. Both are
    confined to the globally convergent window omega in (0, 2), and ORIS's `alpha` damping
    keeps the NET exponent `alpha*omega` inside it.

[INFO] Not a conflict: the `ieppa` -> `oris` rename
  Note: Per the ingest brief, `ieppa`/`ieppa_soft` and `oris`/`oris_soft` denote the same
    solvers throughout, including every occurrence of "iEPPA" in
    `tasks/prd-leafblower-core.md` (US-005, US-007, FR-11 … FR-19, § 6, § 8, § 9).
    `2026-05-30-oris-rename-design.md` freezes the enum VALUES (1 and 8), so identifier
    drift across the corpus carries no semantic difference. All constraint and requirement
    entries use the identifier that appears in their source document, with the mapping
    stated.

[INFO] ABI tripwire values drift across the corpus with no final record
  Note: `EXPECTED_RK_PARAMS_BYTES` appears as 152 (2026-04-24-convergence-metrics-sor),
    220 (2026-04-25-convergence-redesign, "verified on Linux x86_64"), 224 and 232
    (2026-04-29-ieppa-alm-soft-capacity), while `EXPECTED_RK_RESULT_BYTES` appears as 448
    and 480. The codebase map records the current values as 264 and 536. This is expected
    behaviour, not a conflict: every spec instructs the implementer to MEASURE `sizeof()`
    at implementation time and hard-code the observed value beside the `static_assert`. No
    spec claims to state the final size, so nothing needs resolving — but no document in
    the set records the current values either. The PRD's § 6 `rk_params_t` / `rk_result_t`
    layouts are the v1 seed of both structs, not a competing ABI claim.

[INFO] Superseded in full: the 2026-04-18 bounded-convergence fix
  Note: `2026-04-18-bounded-convergence-fix-design.md` rewrites the then-`ieppa` solver
    around `dykstra_solve`, deleting `bcd_sweep`, `bregman_dist`, the EPPA outer loop and
    the post-convergence projection loop. Five days later,
    `2026-04-23-ieppa-faithful-design.md` renames that entire hybrid to `raking` and
    replaces `ieppa` with a new cell-compressed algBCD solver, so the earlier spec's
    subject no longer exists under that name. Its second half (the L-BFGS-B logit-link
    fix) is additionally void, the solver having been removed. Retained in
    `constraints.md` only where its content survives in the renamed solver. It shares its
    2026-04-18 date with the PRD, and the two agree — both predate every change that
    matters.
