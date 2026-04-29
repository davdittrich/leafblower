# SRAA-m Global Safeguard Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix SRAA-m trajectory drift on K≥3 problems by replacing the per-step local safeguard with a global quality-floor safeguard tracking best_err_seen, plus revert-to-best on stall.

**Architecture:** Two changes to src/sraa.hpp only: (1) replace `err_AA ≤ err_plain` with `err_AA ≤ best_err_seen * 1.001` (global quality floor); (2) add revert-to-best after kSRAAStallWindow=15 steps without improvement. No R changes, no solver changes.

**Tech Stack:** C++17 (sraa.hpp header-only template), R testthat3.

---

## Plan Header

- **Mechanism:** Global best_err_seen quality floor + revert-to-best on stall (sraa.hpp)
- **Forbidden:** R-side changes; solver changes (greenkhorn.hpp/sinkhorn.hpp); silent acceptance bypassing best_err_seen tracking; std::swap on X_best (must use copy to keep X_best valid)
- **Audit:** Strict TDD RED→GREEN gate via T_sraa_global testthat case on K=4 overlapping-margin problem; per-file compile gate `R CMD INSTALL --preclean .`; full `devtools::test()` to confirm no regression on K=2 / K=9 cases.

---

## Spec Reference

- Spec: [`docs/superpowers/specs/2026-04-29-i0am-sraa-global-safeguard.md`](../specs/2026-04-29-i0am-sraa-global-safeguard.md)
- Bug evidence: K=9 stepstone harvest with `accelerate=TRUE` produces `max_err = 2.12e-3` while `accelerate=FALSE` produces `max_err = 1.57e-3` (35% worse). Per-step local safeguard `err_AA ≤ err_plain` admits monotone trajectory drift on K≥3 problems.
- Source of truth: `/home/dd/Gemini/leafblower/src/sraa.hpp` (header-only template).

---

## Task 1: Add T_sraa_global test (RED)

**Mechanism:** New testthat3 case in `tests/testthat/` exercising `harvest()` with `method="greenkhorn"` and `accelerate=TRUE/FALSE` on a K=4 overlapping-margin synthetic frame; assert `max_error(AA) <= max_error(plain) * 1.001 + 1e-10`.

**Forbidden:** Modifying sraa.hpp in this task; weakening the assertion to a noop; using a K=2 DGP that masks the trajectory-drift mechanism.

**Audit:** Run the new test against unmodified sraa.hpp and confirm it FAILS (RED gate). Capture exact `me_aa`, `me_plain` from the failure label.

### Steps

- [ ] **1.1 Create the new test file**

  Create `/home/dd/Gemini/leafblower/tests/testthat/test-sraa-global.R` unconditionally.
  T_sraa_global is an SRAA-specific global-safeguard regression test — it belongs in its
  own file, not appended to generic harvest or calibration tests.

  ```bash
  ls /home/dd/Gemini/leafblower/tests/testthat/test-sraa-global.R 2>/dev/null && echo "EXISTS" || echo "CREATING"
  ```

- [ ] **1.2 Add the T_sraa_global test case**

  Create the file with the following content:

  ```r
  test_that("T_sraa_global: greenkhorn+SRAA max_err <= plain on K=4 overlapping-margin problem", {
    set.seed(5); n <- 3000L
    df <- data.frame(
      a = factor(sample(letters[1:4], n, TRUE)),
      b = factor(sample(LETTERS[1:3], n, TRUE)),
      c = factor(sample(c("x","y"),   n, TRUE)),
      d = factor(sample(c("M","F"),   n, TRUE))
    )
    tgt <- list(
      a = setNames(c(0.4, 0.3, 0.2, 0.1), letters[1:4]),
      b = setNames(c(0.5, 0.3, 0.2),      LETTERS[1:3]),
      c = c(x=0.6, y=0.4),
      d = c(M=0.45, F=0.55)
    )
    r_aa    <- suppressWarnings(harvest(df, tgt, method="greenkhorn", accelerate=TRUE,
                                         max_iterations=500L, attach_weights=FALSE))
    r_plain <- suppressWarnings(harvest(df, tgt, method="greenkhorn", accelerate=FALSE,
                                         max_iterations=500L, attach_weights=FALSE))
    me_aa    <- attr(r_aa,    "result")$max_error
    me_plain <- attr(r_plain, "result")$max_error
    expect_lte(me_aa, me_plain * 1.001 + 1e-10,
      label=sprintf("SRAA K=4 (%.2e) must not exceed plain (%.2e)", me_aa, me_plain))
  })
  ```

  If creating a new file, the file header is:

  ```r
  # tests/testthat/test-sraa-global.R
  # T_sraa_global: SRAA-m global safeguard regression test (K=4 overlapping margins)
  ```

  No `library(testthat)` line is required — devtools::test() and R CMD check load testthat automatically.

- [ ] **1.3 Compile gate (no source changes expected, but ensure baseline still installs)**

  ```bash
  cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
  ```

  Expected last line: `* DONE (leafblower)` (or equivalent install-success marker). If install fails, halt and diagnose before continuing.

- [ ] **1.4 RED gate: run the new test against unmodified sraa.hpp**

  ```bash
  cd /home/dd/Gemini/leafblower && Rscript -e 'devtools::test(filter="sraa-global")' 2>&1 | tail -30
  ```

  Expected output: `T_sraa_global` FAILS with a message of the form:

  ```
  Failure (test-sraa-global.R:NN:5): T_sraa_global: greenkhorn+SRAA max_err <= plain on K=4 overlapping-margin problem
  me_aa not less than or equal to me_plain * 1.001 + 1e-10.
  SRAA K=4 (X.XXe-XX) must not exceed plain (Y.YYe-YY)
  ```

  Where `X.XXe-XX > Y.YYe-YY * 1.001`. Record the two numbers in the commit message body. If the test PASSES on the unmodified sraa.hpp, the DGP does not exhibit the bug — halt and revisit DGP parameters before proceeding (do NOT relax the assertion).

- [ ] **1.5 Commit the failing test**

  ```bash
  cd /home/dd/Gemini/leafblower && git add tests/testthat/test-sraa-global.R && git status
  ```

  Expected: shows the new test file staged.

  ```bash
  cd /home/dd/Gemini/leafblower && git commit -m "$(cat <<'EOF'
  test(sraa): add T_sraa_global RED gate for K=4 trajectory drift

  Demonstrates the SRAA-m per-step local safeguard admits monotone
  drift to suboptimal fixed points on K>=3 problems. Failing assertion
  records me_aa vs me_plain so the GREEN gate (Task 2) can confirm
  the global-best safeguard restores parity with plain greenkhorn.
  EOF
  )"
  ```

  Verify:

  ```bash
  cd /home/dd/Gemini/leafblower && git log -1 --stat
  ```

  Expected: single commit, single new file under `tests/testthat/`.

---

## Task 2: Implement global safeguard in sraa.hpp (GREEN)

**Mechanism:** Add `kSRAAGlobalEps`, `kSRAAStallWindow` constants; add `best_err_seen`, `X_best`, `stall_count` fields to `SRAAState`; allocate `X_best` and reset `best_err_seen`/`stall_count` in `init()`; in `clear()` reset only `stall_count` (preserve `best_err_seen` and `X_best`); replace the local safeguard in step 10 with the global quality floor; introduce a `track_best` lambda called at every return in `sraa_step`, including revert-to-best after `kSRAAStallWindow` consecutive non-improving accepts.

**Forbidden:** Resetting `best_err_seen` or `X_best` inside `clear()`; using `std::swap` on `X_best` (must `X = state.X_best;` copy so X_best stays valid for future reverts); skipping `track_best` on any return path; placing the constants outside the existing `kSRAA*` cluster (must keep cluster cohesion); extending `kSRAAMaxM` or any unrelated state.

**Audit:** GREEN gate is `devtools::test(filter="sraa-global")` PASS, plus full `devtools::test()` confirming no regression on existing K=2 and K=9 SRAA tests.

### Steps

- [ ] **2.1 Re-read the full current `sraa.hpp`**

  ```bash
  cd /home/dd/Gemini/leafblower && wc -l src/sraa.hpp && grep -n "^" src/sraa.hpp | head -1
  ```

  Then read the entire file once via the Read tool to confirm exact line numbers for: the `kSRAA*` constant cluster, the `SRAAState` struct, `init()`, `clear()`, and all return points in `sraa_step`. Do NOT proceed to edits until each insertion site is identified by exact surrounding text.

- [ ] **2.2 Add the two new constants**

  Insert immediately after the existing `kSRAARestartGamma` line (line 19 in the current file):

  ```cpp
  static constexpr double kSRAAGlobalEps   = 1e-3;  // 0.1% slack on global safeguard
  static constexpr int    kSRAAStallWindow = 15;    // steps without improvement -> revert
  ```

  Style: match the column alignment of the existing `kSRAA*` block (constants align on `=`, comment aligned past the value).

- [ ] **2.3 Add three new `SRAAState` fields**

  Inside `struct SRAAState`, after the existing `aa_accepted_count` field (or grouped with the other persistent state), add:

  ```cpp
  double best_err_seen = std::numeric_limits<double>::infinity();  // NOT reset by clear()
  std::vector<double> X_best;   // iterate at best_err_seen; allocated in init()
  int stall_count = 0;           // reset by clear()
  ```

  Add `#include <limits>` and `#include <vector>` near the top of the file if not already present (verify by grepping):

  ```bash
  cd /home/dd/Gemini/leafblower && grep -nE "^#include <(limits|vector)>" src/sraa.hpp
  ```

  Add only those headers that are missing.

- [ ] **2.4 Update `init()` — allocate X_best and reset best_err_seen / stall_count**

  Inside the existing `try { ... }` block in `init()`, AFTER `F_cur.assign(...)` and `scratch.assign(...)`, AND BEFORE the `clear()` call, insert:

  ```cpp
  X_best.assign(M, 0.0);
  best_err_seen = std::numeric_limits<double>::infinity();
  stall_count   = 0;
  ```

  Rationale: `clear()` (called immediately after) resets `stall_count` again, but explicit assignment here documents the post-init contract. `best_err_seen` and `X_best` are NOT touched by `clear()`, so init MUST own their first-use initialization.

- [ ] **2.5 Update `clear()` — reset stall_count only**

  Replace the current single-line `clear()` with:

  ```cpp
  void clear() {
      head = 0; count = 0; has_prev = false; prev_resid_norm = 0.0;
      stall_count = 0;
      // best_err_seen and X_best NOT reset - quality floor persists across restarts
  }
  ```

  Verify by grep that no other site mutates `best_err_seen`:

  ```bash
  cd /home/dd/Gemini/leafblower && grep -n "best_err_seen\|X_best\|stall_count" src/sraa.hpp
  ```

  Expected (after this step): three fields declared, three references in `init()`, one reference (`stall_count`) in `clear()`.

- [ ] **2.6 Add `track_best` lambda at the top of `sraa_step`**

  Inside `sraa_step`, BEFORE step 1 (i.e. immediately after the function signature opens and `err_plain` is first computed via `f_eval(state.F_cur)`, but BEFORE any `return` statement), add:

  ```cpp
  // track_best: call before every return, after X is set to accepted iterate.
  // - Records the lowest accepted_err ever seen (monotone quality floor).
  // - Reverts X to the best iterate after kSRAAStallWindow consecutive non-improving accepts.
  // - Uses copy (X = state.X_best) so X_best remains valid for future reverts.
  auto track_best = [&](double accepted_err) {
      if (accepted_err < state.best_err_seen) {
          state.best_err_seen = accepted_err;
          state.X_best        = X;   // O(M) copy - X_best stays valid for future reverts
          state.stall_count   = 0;
      } else {
          state.stall_count++;
      }
      // Revert-to-best on stall (copy not swap - keeps X_best valid)
      if (state.stall_count >= kSRAAStallWindow &&
          state.best_err_seen < std::numeric_limits<double>::infinity()) {
          X = state.X_best;  // O(M) copy
          state.clear();     // resets history + stall_count; preserves best_err_seen + X_best
      }
  };
  ```

  Placement note: the lambda captures `X`, `state` by reference. It must be declared after `X` and `state` are in scope (they are function parameters, so anywhere inside the function body works) and before any `return`. Place it directly after the line `double err_plain = f_eval(state.F_cur);` to keep all returns within scope.

- [ ] **2.7 Replace step 10 safeguard predicate**

  Current (line 183):

  ```cpp
  if (err_AA <= err_plain) {
  ```

  Replace with:

  ```cpp
  if (err_AA <= state.best_err_seen * (1.0 + kSRAAGlobalEps)) {
  ```

  Verify uniqueness before edit:

  ```bash
  cd /home/dd/Gemini/leafblower && grep -nc "if (err_AA <= err_plain)" src/sraa.hpp
  ```

  Expected: `1` (single match). If >1, disambiguate with surrounding context in the Edit call.

- [ ] **2.8 Insert `track_best(...)` calls before EVERY return in `sraa_step`**

  There are six return points. For each, insert a single line `track_best(<err>);` immediately BEFORE the existing `return {...};` (and AFTER the `std::swap(X, ...)` that sets X to the accepted iterate). Exact list:

  | # | Path                          | Site (current line) | After swap                    | Insert                       | Return                            |
  |---|-------------------------------|---------------------|-------------------------------|------------------------------|-----------------------------------|
  | 1 | Restart (norm_k > gamma^2)    | ~95                 | `std::swap(X, state.F_cur);`  | `track_best(err_plain);`     | `return {false, 1, err_plain};`   |
  | 2 | count < kSRAAMinCount         | ~121                | `std::swap(X, state.F_cur);`  | `track_best(err_plain);`     | `return {false, 1, err_plain};`   |
  | 3 | LDLT failure                  | ~156                | `std::swap(X, state.F_cur);`  | `track_best(err_plain);`     | `return {false, 1, err_plain};`   |
  | 4 | NaN guard on err_AA           | ~178                | `std::swap(X, state.F_cur);`  | `track_best(err_plain);`     | `return {false, 1, err_plain};`   |
  | 5 | AA accept (step 10 then)      | ~184                | `std::swap(X, state.scratch);`| `track_best(err_AA);`        | `return {true, 2, err_AA};`       |
  | 6 | AA reject (step 10 else)      | ~191                | `std::swap(X, state.F_cur);`  | `track_best(err_plain);`     | `return {false, 2, err_plain};`   |

  Double-revert safety: paths 1, 3, 4 call `state.clear()` immediately before the swap, which sets `stall_count = 0`. After `track_best` increments, `stall_count` is at most 1, well below `kSRAAStallWindow = 15`, so the revert branch cannot fire on the same call that just cleared. This is the intended contract — verify with grep:

  ```bash
  cd /home/dd/Gemini/leafblower && grep -nE "state\.clear\(\);|std::swap\(X|track_best\(|return \{" src/sraa.hpp
  ```

  Expected pattern at each return site: `state.clear();` (paths 1/3/4 only) → `std::swap(X, ...);` → `track_best(...);` → `return {...};`.

- [ ] **2.9 Per-file compile gate**

  ```bash
  cd /home/dd/Gemini/leafblower && R CMD INSTALL --preclean . 2>&1 | tail -3
  ```

  Expected last line: `* DONE (leafblower)`.

  If compilation fails, halt. Common causes: missing `<limits>`/`<vector>` include; lambda capture of out-of-scope `X` (move lambda decl down); accidental use of `std::swap(X, state.X_best)` (must be `X = state.X_best;`). Diagnose root cause before retrying — do NOT band-aid by removing the failing line.

- [ ] **2.10 GREEN gate: T_sraa_global must PASS**

  ```bash
  cd /home/dd/Gemini/leafblower && Rscript -e 'devtools::test(filter="sraa-global")' 2>&1 | tail -10
  ```

  Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 1 ]` for the sraa-global filter, with `T_sraa_global` listed as PASS. If still RED, halt and diagnose — do NOT relax the assertion or add slack to `kSRAAGlobalEps`.

- [ ] **2.11 Full regression gate: no other test may regress**

  ```bash
  cd /home/dd/Gemini/leafblower && Rscript -e 'devtools::test()' 2>&1 | tail -20
  ```

  Expected: `[ FAIL 0 | WARN ? | SKIP ? | PASS N ]` where N is the full pre-existing PASS count plus 1 (for T_sraa_global). Pay specific attention to any K=2 or K=9 SRAA test (e.g. test names mentioning `sraa`, `accelerate`, `stepstone`) — these characterize the prior behavior and MUST still pass. If any test regresses:
  - Capture the failing test name and its assertion.
  - Trace the failure to one of: (a) `X_best` aliasing/lifetime bug, (b) `track_best` not called on a return, (c) `clear()` accidentally resetting `best_err_seen`.
  - Halt and fix the root cause; do NOT skip or `expect_warning` away the regression.

- [ ] **2.12 Optional empirical confirmation (K=9 stepstone)**

  If a stepstone benchmark script exists (look in `bench/`, `benchmarks/`, `scripts/`):

  ```bash
  cd /home/dd/Gemini/leafblower && ls bench/ benchmarks/ scripts/ 2>/dev/null | grep -iE "stepstone|k9|sraa" || echo "no stepstone bench script"
  ```

  If found, run it and confirm K=9 `max_err` for `accelerate=TRUE` is now within `1.001 *` the plain figure (was `2.12e-3` vs `1.57e-3` before fix). Record the new pair in the commit body. If no script exists, skip this step — the testthat gate is the authoritative pass/fail.

- [ ] **2.13 Commit the implementation**

  ```bash
  cd /home/dd/Gemini/leafblower && git add src/sraa.hpp && git status
  ```

  Expected: only `src/sraa.hpp` staged (compiled artifacts `src/leafblower.so`, `src/ieppa.o` MUST NOT be staged — they are pre-existing modifications unrelated to this task and tracked separately).

  ```bash
  cd /home/dd/Gemini/leafblower && git commit -m "$(cat <<'EOF'
  fix(sraa): replace per-step local safeguard with global quality floor

  Per-step safeguard (err_AA <= err_plain) admitted monotone trajectory
  drift to suboptimal fixed points on K>=3 problems (K=9 stepstone:
  max_err 2.12e-3 vs plain 1.57e-3, 35% worse).

  Changes (sraa.hpp only):
  - kSRAAGlobalEps=1e-3, kSRAAStallWindow=15 (new constants)
  - SRAAState: best_err_seen, X_best, stall_count (new fields)
  - init() allocates X_best, resets best_err_seen and stall_count
  - clear() resets stall_count only; preserves best_err_seen + X_best
  - sraa_step step 10: err_AA <= best_err_seen * 1.001 (global floor)
  - track_best lambda called before every return; reverts X to X_best
    after 15 consecutive non-improving accepts (copy, not swap)

  Closes T_sraa_global RED gate from prior commit. Full devtools::test()
  confirms K=2 and K=9 prior tests still pass.
  EOF
  )"
  ```

  Verify:

  ```bash
  cd /home/dd/Gemini/leafblower && git log -2 --oneline
  ```

  Expected: two new commits — first the test (Task 1), second the fix (Task 2), in that order.

---

## Done Criteria (must all hold)

1. `tests/testthat/test-sraa-global.R` exists and contains `T_sraa_global`.
2. `R CMD INSTALL --preclean .` succeeds with `* DONE (leafblower)`.
3. `devtools::test(filter="sraa-global")` reports `PASS 1, FAIL 0`.
4. `devtools::test()` reports `FAIL 0` overall, with PASS count = prior PASS count + 1.
5. `git log` shows two commits in order: `test(sraa): ...` then `fix(sraa): ...`.
6. Working tree is clean for `tests/testthat/` and `src/sraa.hpp` (only pre-existing modifications to `.wolf/`, `graphify-out/`, `src/ieppa.o`, `src/leafblower.so` remain).

## Halt Conditions

- **SPEC_FAILURE** if T_sraa_global passes against unmodified sraa.hpp at step 1.4 — DGP does not reproduce the bug; do not soften the assertion to make it RED.
- **SPEC_FAILURE** if any pre-existing test regresses at step 2.11 — the global safeguard must be a strict improvement; do not whitelist or skip a regression.
- **SPEC_FAILURE** if `R CMD INSTALL --preclean .` fails after the sraa.hpp edits — fix the root cause, do not retry blindly.
