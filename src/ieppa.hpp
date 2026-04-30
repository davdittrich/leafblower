#pragma once
#include "types.hpp"

namespace lbw {

struct IEPPAResult {
    CalibResult base;
    // ── iEPPA-specific extras ──
    int    M_cell                       = 0;   // compression info
    int    n_cap_active                 = 0;   // cells with W[c] != 1 at convergence
    int    n_xcur_writes_per_iter_linear = 0;  // 0 outside linear path; counter for P1.1 RED test
    double min_alpha_seen               = 1.0; // min alpha over all sweeps
    double final_alpha                  = 1.0; // alpha at solver exit
    int    n_bounds_violated            = 0;   // cell-mode diagnostic: count of w_i outside bounds
    int    n_bounds_clamped             = 0;   // unit-mode action: count of w_i clamped
    // ── Overlay diagnostics ──
    int    homotopy_levels_used         = 0;
    double homotopy_final_factor        = 1.0;
    int    greedy_sweeps_taken          = 0;
    double eta_final                    = 0.0;
    // ── SOR diagnostics ──
    double sor_min_omega                = 1.0;
    int    sor_n_damped                 = 0;
    // ── iEPPA internal metrics ──
    double best_objective_seen          = 0.0;
    double marginal_kl_at_iter          = 0.0;
    // ── ALM diagnostics (ieppa_soft only; zero elsewhere) ──
    double alm_capacity_mu_final        = 0.0;
    int    alm_n_growth_events          = 0;
    double alm_max_dual_norm            = 0.0;
    double alm_sum_drift                = 0.0;
    // ── End iEPPA-specific extras ──
};

// Faithful iEPPA (paper-faithful algBCD at C=0). See
// docs/superpowers/specs/2026-04-23-ieppa-faithful-design.md.
IEPPAResult ieppa_solve(CalibState& state);

} // namespace lbw
