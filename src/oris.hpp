#pragma once
#include <vector>
#include "types.hpp"

namespace lbw {

struct ORISResult {
    CalibResult base;
    // ── ORIS-specific extras ──
    int    M_cell                       = 0;   // compression info
    int    n_cap_active                 = 0;   // cells with W[c] != 1 at convergence
    int    n_xcur_writes_per_iter_last = 0;  // 0 outside linear path; counter for P1.1 RED test
    // Infeasibility-damping diagnostics. FLAT-PATH ONLY: compute_alpha() runs only in
    // the non-accelerated BCD loop, so under accelerate=TRUE (SRAA) both stay 1.0 by
    // design — 1.0 means "damping not applicable on this path", NOT "no damping needed"
    // (CR-B6/y2ks.6: damping shown inert under SRAA; the accept/reject safeguard is the
    // step control). Do not interpret 1.0 as a convergence signal on the SRAA path.
    double min_alpha_seen               = 1.0; // min alpha over all sweeps (flat path)
    double final_alpha                  = 1.0; // alpha at solver exit (flat path)
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
    double sor_omega_mean               = 1.0;
    int    sor_any_latched              = 0;
    int    sor_n_pinned_fb              = 0;
    int    sor_n_warmup_fb              = 0;
    int    sor_n_conv_fb                = 0;
    int    sor_n_resid_grew             = 0;
    int    sor_n_monotone_cd            = 0;
    // ── ORIS internal metrics ──
    double best_objective_seen          = 0.0;
    double marginal_kl_at_iter          = 0.0;
    // ── SRAA-m diagnostics (Anderson Acceleration; 0 when accelerate=FALSE) ──
    int    aa_accepted_count            = 0;   // cumulative AA-accepted super-steps this solve
    // ── ALM diagnostics (oris_soft only; zero elsewhere) ──
    double alm_capacity_mu_final        = 0.0;
    int    alm_n_growth_events          = 0;
    double alm_max_dual_norm            = 0.0;
    double alm_sum_drift                = 0.0;
    // ── SRAA scheduler-demotion flag ──
    // Set TRUE iff SRAA-m acceleration was requested (st.accelerate) AND
    // greedy scheduler was demoted to round-robin (greedy is incompatible
    // with SRAA's fixed-point geometry). FALSE in all other cases.
    bool   sraa_demoted                 = false;
    // ── End ORIS-specific extras ──
};

// ORIS (Over-Relaxed Iterative Scaling; paper-faithful algBCD at C=0). See
// docs/superpowers/specs/2026-04-23-ieppa-faithful-design.md (historical; solver renamed to ORIS).
//
// Optional `lf_capture` (default nullptr): if non-null, on solve exit the
// best-iterate log-Sinkhorn factors `lf` are written to *lf_capture. The
// snapshot mirrors the W_best best-iterate tracking (NOT the trajectory's
// final lf). Used by Newton warm-start (WI-2). Default nullptr keeps every
// existing caller bit-identical.
ORISResult oris_solve(CalibState& state,
                      std::vector<double>* lf_capture = nullptr);

} // namespace lbw
