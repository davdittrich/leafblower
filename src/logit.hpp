#pragma once
#include <cmath>
#include <algorithm>

namespace lbw {

struct LinkFn {
    bool   exponential;     // true = exp link; false = logit link
    double L;               // min_weight
    double U;               // max_weight
    // Deville-Sarndal (1992) logit scale factor: A = (U-L)/((U-1)*(1-L))
    // Governs how fast F(u) transitions from L to U; only valid when !exponential.
    double logit_scale;

    explicit LinkFn(double min_weight, double max_weight)
        : L(min_weight), U(max_weight)
    {
        // Use exponential link only when max_weight is infinite
        exponential = !std::isfinite(U);
        logit_scale = exponential ? 0.0 : (U - L) / ((U - 1.0) * (1.0 - L));
    }

    // Clamp exp(x) to prevent IEEE 754 overflow.
    // 700.0: exp(700) ≈ 1.01e304 < DBL_MAX ≈ 1.8e308.
    // Clamp-safety: both F(u) and H(u) use safe_exp(logit_scale*u) for the
    // same clamped value e, so the identity H'(u) = F(u) is preserved algebraically
    // even when the clamp fires. L-BFGS-B gradient correctness is maintained.
    static double safe_exp(double x) {
        return std::exp(std::min(x, 700.0));
    }

    // F(u): link function mapping dual variable to weight
    double F(double u) const {
        if (exponential) {
            return safe_exp(u);
        }
        double e = safe_exp(logit_scale * u);
        return (L * (U - 1.0) + U * (1.0 - L) * e) / ((U - 1.0) + (1.0 - L) * e);
    }

    // F'(u): derivative of link function
    double dF(double u) const {
        if (exponential) {
            return safe_exp(u);
        }
        double f = F(u);
        double ls = (U - L) / ((U - 1.0) * (1.0 - L));
        return ls * (f - L) * (U - f) / (U - L);
    }

    // H(u): antiderivative of F(u).
    // Logit branch: H(0) = 0 by construction (constant of integration chosen).
    // Exp branch: H(u) = exp(u); H(0) = 1 (additive constant irrelevant for optimization).
    // Logit (Deville-Sarndal 1992):
    //   H(u) = L*u + (U-L)/logit_scale * ln(((U-1)+(1-L)*exp(logit_scale*u)) / (U-L))
    double H(double u) const {
        if (exponential) {
            return safe_exp(u);
        }
        double e = safe_exp(logit_scale * u);
        double num = (U - 1.0) + (1.0 - L) * e;
        return L * u + (U - L) / logit_scale * std::log(num / (U - L));
    }

    // FH(u): F and H from a single safe_exp call.
    // Halves transcendental evaluations in the Wolfe inner loop.
    struct FHResult { double F; double H; };
    FHResult FH(double u) const {
        if (exponential) {
            double e = safe_exp(u);
            return {e, e};
        }
        double e = safe_exp(logit_scale * u);
        double denom = (U - 1.0) + (1.0 - L) * e;
        double f = (L * (U - 1.0) + U * (1.0 - L) * e) / denom;
        double h = L * u + (U - L) / logit_scale * std::log(denom / (U - L));
        return {f, h};
    }

    // F and H from pre-computed e = exp(logit_scale * u).
    // Precondition: !exponential and e == exp(logit_scale * u).
    // Callers MUST check fn.exponential and fall back to FH() when true;
    // calling these with exponential==true produces incorrect results.
    double F_from_e(double e) const {
        return (L * (U - 1.0) + U * (1.0 - L) * e) /
               ((U - 1.0) + (1.0 - L) * e);
    }
    double H_from_e(double e, double u) const {
        double num = (U - 1.0) + (1.0 - L) * e;
        return L * u + (U - L) / logit_scale * std::log(num / (U - L));
    }
};

} // namespace lbw
