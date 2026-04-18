#pragma once
#include <cmath>
#include <algorithm>
#include <limits>

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
        // Use exponential link when min_weight==0 OR max_weight==Inf
        exponential = (L == 0.0 || !std::isfinite(U));
        logit_scale = exponential ? 0.0 : (U - L) / ((U - 1.0) * (1.0 - L));
    }

    // Clamp exp(x) to exp(700) to prevent IEEE 754 overflow
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

    // dF(u): derivative of F w.r.t. u; from quotient rule on F
    double dF(double u) const {
        if (exponential) {
            return safe_exp(u);
        }
        double fu = F(u);
        return logit_scale * (fu - L) * (U - fu) / (U - L);
    }

    // H(u): antiderivative of F(u); H(0) = 0 by construction
    // Exponential: H(u) = exp(u)
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
};

} // namespace lbw
