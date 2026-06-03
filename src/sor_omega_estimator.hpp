#pragma once
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>

namespace lbw {

// Constants — all local to this header (mirrors oris.cpp verbatim).
// Sources: oris.cpp:97, 104, 107-111, 442; calib_dispatch.hpp:304
static constexpr int    kErrCheckInterval    = 10;
static constexpr double kSorEmaAlpha         = 0.2;
static constexpr double kSorProdCeiling      = 1.8;     // production omega ceiling (NOT 1.99)
static constexpr double kSorOscillationDamp  = 0.7;     // sign-flip damp factor
static constexpr double kResidFloor          = 1e-12;   // denominator guard for gate 4
static constexpr int    kSorCooldown         = 5;       // gate-10 cooldown window
static constexpr int    kSorLatchTrips       = 3;       // gate-10 latch threshold
static constexpr double kPinTol              = 1e-9;    // relative pinning tolerance
static constexpr double kSorOmegaMin         = 0.3;     // omega floor (sor_cfg.omega_min default)
static constexpr double kMinSafeTotalWeight  = 1e-100;  // from calib_dispatch.hpp:304
static constexpr int    kSkOmegaBurnin       = 20;      // sinkhorn-local burnin; NOT from sor_cfg.burnin

/// Optimal SOR omega from spectral radius estimate theta2 in [0,1).
/// theta2 < 0 means uninformative — returns ceiling unchanged.
/// Result clamped to [kSorOmegaMin, ceiling].
/// Mirrors oris.cpp:431-435.
inline double omega_from_theta2(double theta2, double ceiling) {
    if (theta2 < 0.0) return ceiling;
    double omega = 2.0 / (1.0 + std::sqrt(1.0 - theta2));
    return std::max(kSorOmegaMin, std::min(omega, ceiling));
}

/// Free-subspace iterate-change theta2 estimator for SOR over-relaxation.
/// Mirrors mode-2 v2 global block in oris.cpp:1684-1813 (e18t.8).
///
/// Observable: S_dX = sum_{c:!pinned} (X[c] - X_snapshot[c])^2
///   -> rho(M_II)^(2*I) per check-interval I.
/// I-th root recovers per-sweep spectral radius; single global omega applied.
///
/// Gate order (INVARIANT): 2 -> 2b -> 3 -> 4 -> 9 -> 6 -> 7+8 -> 10.
/// Gate 9 (oscillation damp) MUST precede gate 6 (ratio>=1) — oris.cpp:1750.
///
/// Caller contract: only call update() when
///   iter >= burnin AND iter % kErrCheckInterval == 0.
struct SorOmegaEstimator {
    // --- Persistent state (mirrors v2_* variables in oris.cpp) ---
    std::vector<double> X_snapshot;
    double S_dX_prev       = std::numeric_limits<double>::infinity();
    double theta2_ema      = -1.0;   // EMA of per-check theta2; -1 = uninformative
    double omega_current   = 1.0;    // last issued omega (for gate-9 damp of prior value)
    int    consec_up       = 0;      // consecutive S_dX non-decreases (gate 10)
    int    cooldown_left   = 0;      // remaining cooldown checks (gate 10)
    int    cooldown_trips  = 0;      // total cooldown trips this solve (gate 10 latch)
    bool   snap_taken      = false;  // guards warm-up: first post-burnin check inits snapshot
    bool   latched         = false;  // permanent omega=1 latch (gate 2)
    bool   prev_decreasing = false;  // S_dX trend for oscillation gate 9

    // --- Observability counters ---
    int    n_pinned_fb   = 0;    // gate-2 fallback count
    int    n_warmup_fb   = 0;    // gate-3 fallback count
    int    n_conv_fb     = 0;    // gate-4 fallback count
    int    n_resid_grew  = 0;    // gate-6 fallback count
    int    n_monotone_cd = 0;    // gate-10 cooldown count
    double omega_sum     = 0.0;  // accumulates gate-7+8 omega values only
    int    omega_n       = 0;    // count of gate-7+8 emissions

    /// Reset all state to initial values.
    void reset(int M = 0) {
        X_snapshot.assign(M, 0.0);
        S_dX_prev      = std::numeric_limits<double>::infinity();
        theta2_ema     = -1.0;
        omega_current  = 1.0;
        consec_up      = 0;
        cooldown_left  = 0;
        cooldown_trips = 0;
        snap_taken     = false;
        latched        = false;
        prev_decreasing = false;
        n_pinned_fb    = 0;
        n_warmup_fb    = 0;
        n_conv_fb      = 0;
        n_resid_grew   = 0;
        n_monotone_cd  = 0;
        omega_sum      = 0.0;
        omega_n        = 0;
    }

    /// update(): call at iter % kErrCheckInterval == 0 AND iter >= burnin.
    ///
    /// @param X         current cell weights (M_cell elements)
    /// @param is_pinned per-cell pinning mask (true = at bound)
    /// @param M_cell    number of cells
    /// @param ceiling   omega cap (use kSorProdCeiling for production)
    /// @returns recommended omega (1.0 = identity / conservative)
    double update(const std::vector<double>& X,
                  const std::vector<bool>& is_pinned,
                  int M_cell,
                  double ceiling = kSorProdCeiling) {
        const int M = M_cell;

        // Ensure snapshot buffer is sized.
        if ((int)X_snapshot.size() != M) X_snapshot.assign(M, 0.0);

        // --- Compute S_dX and n_free in one pass (oris.cpp:1686-1694) ---
        double S_dX = 0.0;
        int    n_free = 0;
        for (int c = 0; c < M; c++) {
            if (!is_pinned[c]) {
                double dx = X[c] - X_snapshot[c];
                S_dX += dx * dx;
                n_free++;
            }
        }

        // Gate 2: free-set empty OR permanently latched -> omega=1, reset EMA.
        // (oris.cpp:1697-1704)
        if (n_free == 0 || latched) {
            theta2_ema = -1.0;
            X_snapshot = X;
            S_dX_prev  = S_dX;
            n_pinned_fb++;
            return _emit(1.0);
        }

        // Gate 2b: cooldown active -> omega=1, decrement, reset EMA.
        // (oris.cpp:1707-1715)
        if (cooldown_left > 0) {
            theta2_ema = -1.0;
            cooldown_left--;
            n_monotone_cd++;
            X_snapshot = X;
            S_dX_prev  = S_dX;
            return _emit(1.0);
        }

        // Gate 3: warm-up — first post-burnin check: initialize snapshot, hold omega=1.
        // (oris.cpp:1718-1726)
        if (!snap_taken) {
            X_snapshot    = X;
            snap_taken    = true;
            theta2_ema    = -1.0;
            S_dX_prev     = std::numeric_limits<double>::infinity();
            n_warmup_fb++;
            return _emit(1.0);
        }

        // Gate 4: finiteness + tiny denominator or tiny free mass.
        // (oris.cpp:1729-1743)
        {
            double free_mass = 0.0;
            for (int c = 0; c < M; c++)
                if (!is_pinned[c]) free_mass += X[c];
            if (!std::isfinite(S_dX) ||
                S_dX_prev < kResidFloor ||
                free_mass < kMinSafeTotalWeight) {
                theta2_ema = -1.0;
                n_conv_fb++;
                X_snapshot = X;
                S_dX_prev  = S_dX;
                return _emit(1.0);
            }
        }

        // --- Gates 9, 6, 7+8, 10; lag update (oris.cpp:1745-1811) ---
        bool dX_decreasing = (S_dX < S_dX_prev);
        bool dX_sign_flip  = !dX_decreasing && prev_decreasing;

        // Gate 9: oscillation — damp prior omega * factor, reset EMA, advance lag.
        // MUST come BEFORE gate 6 ratio check. (oris.cpp:1749-1764)
        if (dX_sign_flip) {
            double damped  = std::max(kSorOmegaMin, omega_current * kSorOscillationDamp);
            theta2_ema     = -1.0;
            prev_decreasing = dX_decreasing;
            X_snapshot     = X;
            S_dX_prev      = S_dX;
            return _emit(damped);
        }

        double ratio    = S_dX / S_dX_prev;
        double omega_out;

        // Gate 6: iterate-change grew -> reset EMA, omega=1, consec_up++.
        // (oris.cpp:1769-1773)
        if (ratio >= 1.0) {
            theta2_ema = -1.0;
            consec_up++;
            n_resid_grew++;
            omega_out = 1.0;
        } else {
            // Gates 7+8: I-th root cadence recovery, EMA, formula.
            // (oris.cpp:1775-1791)
            consec_up = 0;
            double theta2 = std::pow(ratio, 1.0 / static_cast<double>(kErrCheckInterval));
            theta2 = std::min(std::max(theta2, 0.0), 1.0 - 1e-9);
            if (theta2_ema < 0.0) {
                theta2_ema = theta2;   // seed first informative sample
            } else {
                theta2_ema = kSorEmaAlpha * theta2 +
                             (1.0 - kSorEmaAlpha) * theta2_ema;
            }
            omega_out = omega_from_theta2(theta2_ema, ceiling);
            omega_sum += omega_out;
            omega_n++;
        }

        // Gate 10: monotone latch — cooldown and permanent latch.
        // NOTE: oris.cpp uses hardcoded 3 here matching kSorLatchTrips. (oris.cpp:1794-1805)
        if (consec_up >= kSorLatchTrips) {
            omega_out      = 1.0;
            cooldown_left  = kSorCooldown;
            cooldown_trips++;
            consec_up      = 0;
            theta2_ema     = -1.0;
            n_monotone_cd++;
            if (cooldown_trips >= kSorLatchTrips) {
                latched = true;
            }
        }

        // Update lag. (oris.cpp:1807-1811)
        prev_decreasing = dX_decreasing;
        X_snapshot      = X;
        S_dX_prev       = S_dX;
        return _emit(omega_out);
    }

private:
    double _emit(double omega) noexcept {
        omega_current = omega;
        return omega;
    }
};

}  // namespace lbw
