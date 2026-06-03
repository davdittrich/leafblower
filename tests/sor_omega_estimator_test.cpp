// Standalone unit tests for src/sor_omega_estimator.hpp
// Compile: g++ -std=c++17 -I<repo>/src -o /tmp/sor_test tests/sor_omega_estimator_test.cpp
// Run:     /tmp/sor_test
// Exit 0 = all PASS; non-zero = at least one FAIL.
#include <cstdio>
#include <cmath>
#include <vector>
#include <cassert>
#include "sor_omega_estimator.hpp"

using namespace lbw;

static int failures = 0;

static void check(bool ok, const char* name) {
    if (ok) {
        std::printf("PASS  %s\n", name);
    } else {
        std::printf("FAIL  %s\n", name);
        failures++;
    }
}

// ── Test 1: omega_from_theta2 basic properties ────────────────────────────────
static void test_omega_from_theta2() {
    // theta2=0 -> omega = 2/(1+sqrt(1))  = 2/2 = 1.0
    double o0 = omega_from_theta2(0.0, kSorProdCeiling);
    check(std::abs(o0 - 1.0) < 1e-12, "omega_from_theta2: theta2=0 -> 1.0");

    // theta2 < 0 -> returns ceiling unchanged
    double o_neg = omega_from_theta2(-0.5, kSorProdCeiling);
    check(std::abs(o_neg - kSorProdCeiling) < 1e-12,
          "omega_from_theta2: theta2<0 -> ceiling");

    // theta2 very close to 1 -> large omega, should be clamped to ceiling
    double o_hi = omega_from_theta2(1.0 - 1e-10, kSorProdCeiling);
    check(std::abs(o_hi - kSorProdCeiling) < 1e-9,
          "omega_from_theta2: theta2~1 clamped to ceiling");

    // omega_from_theta2 with a ceiling of 0.1 should trigger kSorOmegaMin clamp
    // theta2=0 gives omega=1.0 but ceiling=0.1 < kSorOmegaMin=0.3, so min(1.0,0.1)=0.1
    // then max(0.3, 0.1) = 0.3
    double o_low = omega_from_theta2(0.0, 0.1);
    check(std::abs(o_low - kSorOmegaMin) < 1e-12,
          "omega_from_theta2: clamped to kSorOmegaMin when ceiling < floor");
}

// Helper: make all-unpinned bool vector
static std::vector<bool> all_free(int M) {
    return std::vector<bool>(M, false);
}

// Helper: make X with all elements equal to val
static std::vector<double> make_X(int M, double val) {
    return std::vector<double>(M, val);
}

// ── Test 2: Gate 3 warm-up — first call returns 1.0 and sets snapshot ─────────
static void test_gate3_warmup() {
    SorOmegaEstimator est;
    int M = 4;
    std::vector<double> X = {1.0, 2.0, 3.0, 4.0};
    auto pins = all_free(M);

    double omega = est.update(X, pins, M);
    check(std::abs(omega - 1.0) < 1e-12,
          "gate3_warmup: first call returns 1.0");
    check(est.snap_taken,
          "gate3_warmup: snap_taken set after first call");
    check(est.X_snapshot == X,
          "gate3_warmup: snapshot equals X after first call");
    check(est.n_warmup_fb == 1,
          "gate3_warmup: n_warmup_fb incremented");
}

// ── Test 3: Gate 4 tiny denominator — S_dX_prev=0 returns 1.0 ────────────────
static void test_gate4_tiny_denom() {
    SorOmegaEstimator est;
    int M = 4;
    auto pins = all_free(M);

    // First call (gate 3): initialises snapshot = {1,1,1,1}, S_dX_prev=inf
    std::vector<double> X1 = make_X(M, 1.0);
    est.update(X1, pins, M);

    // Manually force S_dX_prev below kResidFloor to trigger gate 4
    est.S_dX_prev = kResidFloor / 2.0;   // < kResidFloor

    std::vector<double> X2 = make_X(M, 1.5);
    double omega = est.update(X2, pins, M);
    check(std::abs(omega - 1.0) < 1e-12,
          "gate4_tiny_denom: returns 1.0 when S_dX_prev < kResidFloor");
    check(est.n_conv_fb >= 1,
          "gate4_tiny_denom: n_conv_fb incremented");
}

// ── Test 4: Gates 7+8 — two consecutive decreasing S_dX → omega > 1.0 ────────
static void test_gates78_recovery() {
    SorOmegaEstimator est;
    int M = 4;
    auto pins = all_free(M);

    // Gate 3: seed snapshot (X_snapshot = {0,0,0,0})
    std::vector<double> X0 = make_X(M, 0.0);
    est.update(X0, pins, M);

    // Set up a known S_dX_prev so we control the ratio.
    // We want ratio = S_dX_new / S_dX_prev < 1 (decreasing).
    // S_dX_prev must be > kResidFloor.
    // X1 differs from X0 by 1 per cell: S_dX = M * 1^2 = 4.
    // Manually override S_dX_prev to something larger so ratio < 1.
    est.S_dX_prev = 1000.0;  // large denominator -> ratio will be tiny
    // Also ensure prev_decreasing=false so gate 9 won't fire on first informative call.
    est.prev_decreasing = false;

    std::vector<double> X1 = make_X(M, 1.0);
    double o1 = est.update(X1, pins, M);
    // S_dX = 4, ratio = 4/1000 = 0.004 -> theta2 = 0.004^(1/10) ~ 0.631
    // omega = 2/(1+sqrt(1-0.631)) ~ 2/(1+0.607) ~ 1.24 — well above 1.0
    check(o1 > 1.0,
          "gates78_recovery: first informative call returns omega > 1.0");

    // Second decreasing call: snap is now X1, so we advance further.
    // S_dX_prev is now ~4 (from prior call). Make X2 close to X1 so S_dX < 4.
    std::vector<double> X2 = {1.1, 1.1, 1.1, 1.1};
    // Set S_dX_prev explicitly larger to ensure ratio < 1
    est.S_dX_prev = 100.0;
    est.prev_decreasing = true;  // was decreasing last time
    double o2 = est.update(X2, pins, M);
    // S_dX = 4*(0.1^2) = 0.04, ratio = 0.04/100 = 4e-4, ratio<1 -> gate 7+8
    // EMA updated -> omega should still be > 1
    check(o2 > 1.0,
          "gates78_recovery: second consecutive decreasing call returns omega > 1.0");
}

// ── Test 5: Gate 9 ordering — oscillation fires BEFORE gate 6 ────────────────
// Critical: decrease then increase (sign flip) must return damped omega,
// NOT reset to 1.0 as gate 6 would do.
static void test_gate9_before_gate6() {
    SorOmegaEstimator est;
    int M = 4;
    auto pins = all_free(M);

    // Gate 3 warm-up
    std::vector<double> X0 = make_X(M, 0.0);
    est.update(X0, pins, M);

    // Seed a non-trivial omega (> 1) by simulating a decreasing step through gate 7+8
    // Set prev_decreasing=false, S_dX_prev large -> first informative call -> omega>1
    est.S_dX_prev      = 1000.0;
    est.prev_decreasing = false;
    std::vector<double> X1 = make_X(M, 1.0);  // S_dX=4, ratio=4/1000 tiny
    double o1 = est.update(X1, pins, M);
    // o1 > 1.0 now confirmed by test 4; omega_current = o1
    check(o1 > 1.0, "gate9_ordering: pre-condition omega > 1.0 after decreasing step");

    // Now simulate oscillation: prev_decreasing was true (S_dX went down),
    // but next step S_dX goes UP (sign flip). Gate 9 must fire first.
    // Construct X2 where S_dX > S_dX_prev (S_dX_prev was set by last update = 4).
    // X2 = X1 + 10 per cell -> dx from X_snapshot(X1) = 10 -> S_dX = 400 > 4.
    std::vector<double> X2 = make_X(M, 11.0);
    // Confirm prev_decreasing is still true from last gate-7+8 call
    check(est.prev_decreasing, "gate9_ordering: prev_decreasing is true before sign-flip call");

    double o2 = est.update(X2, pins, M);

    // Gate 9 fires (sign flip): damped = max(0.3, o1 * 0.7) which is < o1 but > 0.3
    // Gate 6 would set omega=1.0 instead
    double expected_damped = std::max(kSorOmegaMin, o1 * kSorOscillationDamp);
    bool damped_correctly  = std::abs(o2 - expected_damped) < 1e-10;
    bool not_reset_to_one  = std::abs(o2 - 1.0) > 1e-10;  // gate 6 would give 1.0
    check(damped_correctly && not_reset_to_one,
          "gate9_ordering: sign-flip returns damped omega (gate9 before gate6)");
}

// ── Test 6: Gate 6 — non-sign-flip increase resets omega to 1.0 ──────────────
static void test_gate6_ratio_ge1() {
    SorOmegaEstimator est;
    int M = 4;
    auto pins = all_free(M);

    // Warm-up
    std::vector<double> X0 = make_X(M, 0.0);
    est.update(X0, pins, M);

    // Prime: prev_decreasing = false so no sign flip possible
    est.S_dX_prev      = 1000.0;
    est.prev_decreasing = false;

    // Call with S_dX > S_dX_prev (ratio >= 1): no sign flip since prev_decreasing=false
    std::vector<double> X1 = make_X(M, 20.0);  // S_dX = 4*400 = 1600 > 1000
    int n_before = est.n_resid_grew;
    double omega = est.update(X1, pins, M);

    check(std::abs(omega - 1.0) < 1e-12,
          "gate6_ratio_ge1: ratio>=1 returns omega=1.0");
    check(est.n_resid_grew == n_before + 1,
          "gate6_ratio_ge1: n_resid_grew incremented");
    check(est.consec_up == 1,
          "gate6_ratio_ge1: consec_up incremented to 1");
}

// ── Test 7: Gate 10 latch — repeated ratio>=1 eventually latches ──────────────
static void test_gate10_latch() {
    // To latch: need kSorLatchTrips=3 cooldown trips.
    // Each trip requires kSorLatchTrips=3 consecutive gate-6 fires.
    // After each trip, there are kSorCooldown=5 cooldown checks (gate 2b).
    // We drive the estimator through the full sequence.
    SorOmegaEstimator est;
    int M = 2;
    auto pins = all_free(M);

    // Warm-up
    est.update(make_X(M, 0.0), pins, M);

    // Helper: inject N consecutive ratio>=1 steps (no sign flip)
    // Uses large X increments to ensure S_dX > S_dX_prev consistently.
    auto inject_ratio_ge1 = [&](int N) {
        for (int i = 0; i < N; i++) {
            double base = (double)(i + 1) * 100.0;
            // Force no sign-flip by ensuring prev_decreasing = false
            est.prev_decreasing = false;
            est.S_dX_prev       = 1.0;  // tiny denom so ratio = S_dX/1 >> 1
            est.update(make_X(M, base), pins, M);
        }
    };

    // Helper: drain cooldown checks
    auto drain_cooldown = [&]() {
        for (int i = 0; i < kSorCooldown + 2; i++) {
            est.update(make_X(M, 0.0), pins, M);
        }
    };

    // Trip 1: 3 gate-6 fires -> gate 10 -> cooldown
    inject_ratio_ge1(kSorLatchTrips);
    check(est.cooldown_trips == 1, "gate10_latch: cooldown_trips=1 after first trip");
    check(!est.latched, "gate10_latch: not latched after first trip");
    drain_cooldown();

    // Trip 2
    inject_ratio_ge1(kSorLatchTrips);
    check(est.cooldown_trips == 2, "gate10_latch: cooldown_trips=2 after second trip");
    check(!est.latched, "gate10_latch: not latched after second trip");
    drain_cooldown();

    // Trip 3 — should latch
    inject_ratio_ge1(kSorLatchTrips);
    check(est.cooldown_trips == 3, "gate10_latch: cooldown_trips=3 after third trip");
    check(est.latched, "gate10_latch: latched after kSorLatchTrips trips");
}

// ── Test 8: Gate 2 latched — any call returns 1.0 ────────────────────────────
static void test_gate2_latched() {
    SorOmegaEstimator est;
    int M = 4;
    auto pins = all_free(M);

    // Manually force latch
    est.latched = true;

    // Warm-up call should still return 1.0
    double o1 = est.update(make_X(M, 1.0), pins, M);
    check(std::abs(o1 - 1.0) < 1e-12,
          "gate2_latched: returns 1.0 when latched (first call)");

    // Any subsequent call also returns 1.0
    double o2 = est.update(make_X(M, 2.0), pins, M);
    check(std::abs(o2 - 1.0) < 1e-12,
          "gate2_latched: returns 1.0 when latched (second call)");
}

int main() {
    std::printf("=== SorOmegaEstimator unit tests ===\n");
    test_omega_from_theta2();
    test_gate3_warmup();
    test_gate4_tiny_denom();
    test_gates78_recovery();
    test_gate9_before_gate6();
    test_gate6_ratio_ge1();
    test_gate10_latch();
    test_gate2_latched();
    std::printf("===================================\n");
    if (failures == 0) {
        std::printf("ALL PASS (%d tests)\n", 0);  // counted by check() calls
    } else {
        std::printf("%d FAILURE(S)\n", failures);
    }
    return failures > 0 ? 1 : 0;
}
