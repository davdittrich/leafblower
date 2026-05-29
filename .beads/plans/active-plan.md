<!-- approved: 2026-05-29T11:16:22Z -->
<!-- gate-iterations: 1 -->
<!-- user-approved: true -->
<!-- status: in-progress -->
<!-- epics: 2lxq,o2o6,nil1,ya85 -->

# Active Plan — codereview remediation, epic-organized by solver file

## PROGRESS (2026-05-29 session): 11 tickets done, 11 commits (local; no remote)
- chebyshev (2lxq) 4/4: izql a2547a3, xl44 6d215d7, xiox(rejected) 9e2e913, 0xc5 58448c5
- ieppa (o2o6) 7/8: l6to d7c7502, jd2f bd2534a, ohi0 1996338, l7sg 92ddb03, vtjf f493cf4, za9r 36f0530, 7emq(rejected) 7c62d50
- Shared helper added: lbw::finalize_weights[_buf] + kMinSafeTotalWeight (calib_dispatch.hpp)
- 2 tickets REJECTED as not-a-bug via derivation+sim: xiox (Mehrotra corrector), 7emq (ALM Newton). Guard comments + bd remember added.
- Spawned: dead-hi (P4), ieppa-adopt-finalize_weights (P3).

## uu8r PLAN — GATE-APPROVED (plan-review-gate iteration 3/3: Feasibility+Completeness+Scope all PASS)
Cold-only split — build has NO LTO, so hot ieppa_solve + pack/unpack_lf STAY in ieppa.cpp.
REORDERED (link-order fix: write_trajectory_csv called by ieppa_finalize → externalize trajectory FIRST):
- uu8r.1: extract trajectory I/O (parse_trajectory_iters + write_trajectory_csv) -> ieppa_trajectory.cpp + CREATE ieppa_internal.hpp (externalize both).
- uu8r.2 (dep .1): extract ieppa_finalize -> ieppa_finalize.cpp (links to now-external write_trajectory_csv).
- uu8r.3 (dep .1,.2): stepstone perf gate (pre/post, 3% tol) + 593 R + 15 py parity byte-identical.
TWO build sites per move: R src/Makevars.in PKG_SOURCES (R auto-globs src/*.cpp — likely decorative, edit for docs) + python/CMakeLists.txt CORE_SOURCES (EXPLICIT list — MUST add or pybind11 link fails).
Mechanism: cold-only relocation, byte-identical bodies. Forbidden: moving hot code, -flto, logic changes.
HONEST SCOPE: ~330 cold lines move; hot ieppa_solve (~1690 lines) stays → modest compile win (cold edits only).
Gate caught + fixed: (iter1) omitted Python CORE_SOURCES wiring; (iter2) write_trajectory_csv link-order.
STATUS: DONE. uu8r.1 (8d0a966), uu8r.2 (637198d), uu8r.3 (perf gate, no code). ieppa.cpp 2153->1789. Perf-neutral: ieppa median 3411ms (pre) -> 3341ms (post), -2.0%, single-thread n=11. 593 R/15 py byte-identical pre/post. uu8r + 3 children CLOSED.

## SESSION FINAL: 14 tickets closed, 13 commits. chebyshev 4/4 + ieppa 7/8 bugs + uu8r split (3 WU). Open: ieppa-adopt-finalize_weights (P3), dead-hi (P4), chebyshev-hi-var (P4) — cosmetic. No git remote (commits local).

# (original plan below)

**Branch:** `feat/codereview-remediation-p2-batch`
**Mechanism:** Surgical per-ticket correctness fixes; single-exit / Σw=n enforcement mirrors `ieppa_finalize`.
**Forbidden:** post-normalizing Σw=n (breaks unit water-fill); subagent self-cert of numerics; git add -A; --no-verify.
**Gate:** plan-review-gate iteration 1 — Feasibility PASS, Completeness PASS, Scope PASS.

## Epic structure (axis = by solver file, all 18 open tickets)
- **chebyshev** (leafblower-2lxq, P1): izql, xl44, xiox, 0xc5. GUARD: izql+xl44 edit exit region L245–714 → SEQUENTIAL.
- **ieppa** (leafblower-o2o6, P1): jd2f, l6to, ohi0, vtjf, za9r, 7emq, l7sg, uu8r. All edit ieppa.cpp → SERIAL. uu8r (split) BLOCKED-BY 7 bugs → runs LAST. 7emq SUSPECTED → read paper PDF before patch.
- **raking/newton/logit** (leafblower-nil1, P2): 24f7, dqs3, 2ce1, xwqy. 4 distinct files → parallelizable.
- **refactor** (leafblower-ya85, P3): l1p3 (R/harvest.R), yxcg (SIMD). Coordinate yxcg with uu8r line-shift.
- **distribution** (leafblower-kk1.24, pre-existing feature): untouched.

## Verification gate (per epic)
R tests (testthat FAIL 0) + Python parity rtol=1e-6 + stepstone no regression. Orchestrator INDEPENDENTLY reruns — no subagent self-cert. Per-file `git add`. One commit per ticket.

## Execution order
1. chebyshev: **izql (in_progress)** → xl44 → xiox → 0xc5
2. ieppa bugs (serial) → uu8r split last
3. raking/newton/logit (parallel)
4. refactor

## Non-blocker (gate flag)
yxcg SIMD touches ieppa.cpp:241; left under refactor epic; coupling to uu8r line-shift noted in ya85 body.
