#pragma once
// calib_dispatch.hpp — shared CalibMetric / CalibRule helpers (WU-B).
// Eliminates triplicated switch blocks across ieppa, raking, and lbfgsb_solver.
//
// Semantics:
//   select_metric — map CalibMetric enum → the pre-computed scalar value.
//   apply_rule    — apply a CalibRule stopping condition, updating prev in-place.
//
// apply_rule behaviour (all rules require tol > 0 and finite curr):
//   THRESHOLD  : curr < tol
//   IMPROVEMENT: |curr - prev| / prev < tol
//                Special: curr <= 1e-15 → trivially converged (metric at machine zero).
//                Skip (false) when prev is non-finite or prev <= 1e-15 (first check).
//   PLATEAU    : !(curr < prev * (1 - tol))  — i.e. metric did not improve by ≥tol fraction.
//                Skip (false) when prev is non-finite (first check).
//
// The caller is responsible for passing the correct local variable names for each
// metric; lbfgsb passes pct_change for L1_WEIGHT because it is a batch solver.

#include "types.hpp"
#include <cmath>
#include <limits>

namespace lbw {

// Select the active metric value from the pre-computed metric set.
// The caller provides all six metric values in canonical order; only the
// one matching `metric` is used.
inline double select_metric(
    CalibMetric metric,
    double max_err,
    double mean_err,
    double kl,
    double chi2,
    double grake_norm,
    double l1_weight) noexcept
{
    switch (metric) {
        case CalibMetric::MAX_ERR:    return max_err;
        case CalibMetric::MEAN_ERR:   return mean_err;
        case CalibMetric::KL:         return kl;
        case CalibMetric::CHI2:       return chi2;
        case CalibMetric::GRAKE_NORM: return grake_norm;
        case CalibMetric::L1_WEIGHT:  return l1_weight;
    }
    return max_err;  // unreachable; exhaustive enum silences -Wreturn-type
}

// Apply the stopping rule.
// `prev` is read and updated by this function; initialise to +infinity before
// the first iteration.  Returns false on the first call when prev==+inf.
//
// Note: lbfgsb_solver does NOT use this helper because it is a batch solver
// with no per-iteration baseline; it retains its own THRESHOLD-only logic.
inline bool apply_rule(
    CalibRule rule,
    double    curr,
    double&   prev,   // in-out: updated to curr on return when finite
    double    tol) noexcept
{
    if (!std::isfinite(curr) || tol <= 0.0) return false;

    bool converged = false;
    switch (rule) {
        case CalibRule::THRESHOLD:
            converged = (curr < tol);
            break;

        case CalibRule::IMPROVEMENT: {
            // When curr is at machine-zero the metric cannot improve further;
            // treat relative change as 0 (trivially converged).
            if (curr <= 1e-15) {
                converged = true;
            } else if (std::isfinite(prev) && prev > 1e-15) {
                const double rel = std::fabs(curr - prev) / prev;
                converged = (rel < tol);
            }
            // else: first check (prev==inf) or prev near-zero → skip
            break;
        }

        case CalibRule::PLATEAU:
            // Converge when curr did NOT drop by at least tol fraction vs prev.
            // Skip on first check (prev==inf); once prev is finite always evaluate.
            if (std::isfinite(prev)) {
                converged = !(curr < prev * (1.0 - tol));
            }
            break;
    }

    prev = curr;   // always update so next call has a valid baseline
    return converged;
}

} // namespace lbw
