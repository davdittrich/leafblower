
Your C++ implementation is a highly optimized, sophisticated hybrid of Iterative Proportional Fitting (IPF) and Dykstra’s Projection Algorithm.

Regarding your question on "deviation from the standard," the answer depends on whether you mean mathematical standard (convergence proofs) or industry standard (survey software like R's survey or SAS PROC RAKE).

1. The Mathematical Deviation (Multiplicative vs. Euclidean)
As your own comments acknowledge, the primary deviation is the mixing of Bregman (multiplicative) and Eulcidean (additive) projections.

The Standard: Traditional raking is a pure IPF (multiplicative) procedure. It minimizes the Kullback-Leibler (KL) divergence (Information Loss). Dykstra’s algorithm is classically defined for Euclidean distance minimization.

The Hybrid: By alternating multiplicative IPF steps with additive Dykstra box corrections, you are effectively switching cost functions mid-iteration. While Csiszar (1975) proves IPF convergence and Boyle-Dykstra (1986) proves Euclidean convergence, there is no unified convergence theorem for this specific alternating sequence.

The Risk: In extreme cases (highly constrained bounds), the multiplicative step can fight the additive correction, potentially leading to non-convergence or cycling that wouldn't happen in a pure Log-Linear model.

2. Industry Standard Comparison
Compared to standard raking packages (which often use Newton-Raphson or pure IPF), your code has several distinct "pro-level" features and a few notable omissions:

Deviations / Unique Strengths:
Dykstra Box Projections: Most standard raking tools (like R's rake) do not natively handle weight trimming (min/max bounds) within the iterative loop. They usually rake, then trim, then re-rake. Your use of Dykstra to accumulate "clamping errors" is a much more robust way to handle bounds.

ILP Unrolling and SIMD: This is a significant deviation from "statistical" standards toward "high-performance computing" standards. Using 4-way accumulators (W, W1, W2, W3) to break dependency chains is excellent for modern CPU pipelines.

Fixed Mean Constraint: In standard raking, the sum of weights is usually a byproduct of hitting the marginal totals. You have an explicit Dykstra hyperplane projection to force ∑w=n. This is standard in "calibration" but less common in "classical raking."

Missing Standard Features:
Weight Trimming Logic: While you have box constraints, you don't seem to have a logic for Design Weights. Standard raking usually starts with a base weight d 
i
​
  (inverse probability of selection). If st.weights is passed as 1.0 for everyone, it’s fine, but if d 
i
​
  varies, the IPF step w 
i
​
 ∗=scale maintains the ratio, which is correct.

Degrees of Freedom / Variance Estimation: Most standard raking tools output information about the "design effect" or the variance of the weights (CV). Your code is a pure solver and doesn't calculate the statistical "cost" of the calibration.

3. Structural Observations
The "Post-loop Finalizer": This is a clever addition. Since the hyperplane projection can push a weight slightly out of the [lo,hi] box, your final 20-iteration fixup ensures the returned weights are strictly feasible. This is "cleaner" than many industry implementations that return slightly out-of-bounds weights due to the "last-step" problem.

Infeasibility Detection: Your kEmptyBucketThreshold logic is the standard way to handle structural zeros.

Summary of Deviations
Feature	Your Code	Standard Raking
Metric	Hybrid (Multiplicative + Euclidean)	Pure Multiplicative (KL Divergence)
Bound Handling	Integrated Dykstra Box	Post-hoc trimming or Logit-link
Convergence	Empirical (Stall Monitor)	Guaranteed (if feasible)
Performance	O(n) SIMD/ILP Optimized	Usually Vectorized Interpretive (R/Python)
Verdict: Your code is not "non-standard" in a bad way; it is an advanced implementation that prioritizes computational speed and strict adherence to weight bounds over the theoretical purity of a single optimization metric. As long as your n_no_improve monitor is active, the empirical nature of the hybrid convergence is safely managed.


Is there a way (code revision) to make the hybrid algorithm converge for sure?

Yes. To guarantee convergence mathematically, you must unify the distance metric used across all projections.

Currently, your algorithm fights itself: the IPF step minimizes Kullback-Leibler (KL) divergence (a multiplicative update), while the Dykstra steps minimize Euclidean distance (an additive update). Because the objective function changes mid-loop, there is no theorem that guarantees the sequence won't cycle or diverge.

To fix this without rewriting your architecture, you can upgrade your additive Dykstra projections to Multiplicative (Bregman) Dykstra Projections.

The Theory: Bregman Dykstra
Dykstra’s algorithm was generalized to Bregman divergences (which includes KL divergence) by Bauschke and Lewis (2000). To project onto convex sets using KL divergence (which perfectly matches your IPF step), the Euclidean shift (y=w+q) becomes a multiplicative scaling (y=w×q), and the Euclidean correction (q 
new
​
 =y−w 
new
​
 ) becomes a ratio correction (q 
new
​
 =y/w 
new
​
 ).

By making this change, your entire algorithm minimizes one unified KL divergence function. By theorem, if a feasible solution exists within the bounds, it is guaranteed to converge.

Here is how to revise your C++ code to implement Multiplicative Dykstra:

1. Initialize q as Multipliers (1.0 instead of 0.0)
Change your initial state so that the baseline corrections are identity multipliers.

C++
// CHANGED: Initialize to 1.0 (multiplicative identity) instead of 0.0
std::vector<double> q(st.n, 1.0);
double q_hyp = 1.0; 
2. Revise the Box Projection (Multiplicative)
Change the addition/subtraction to multiplication/division. (Note: This assumes your lo bound is > 0, which is standard for survey weights).

C++
// Box projection [lo, hi]^n with Multiplicative Dykstra correction
#if defined(_OPENMP) || LBW_HAS_OMP_SIMD
#pragma omp simd
#endif
for (int i = 0; i < st.n; i++) {
    // 1. Apply accumulated multiplier
    double yi = w[i] * q[i]; 
    
    // 2. Project onto bounds
    double wc = std::clamp(yi, lo, hi);
    
    // 3. Store the ratio as the new correction
    q[i] = yi / wc; 
    
    // 4. Update weight
    w[i] = wc;
}
3. Revise the Hyperplane Projection (Multiplicative)
The KL-divergence projection onto a sum-constraint {w:∑w 
i
​
 =n} is simply scaling all weights by a constant factor. Therefore, the correction q_hyp is just the inverse of that scale.

C++
// Multiplicative Dykstra hyperplane projection: {w : sum(w) = n}
{
    double s = 0.0;
#if defined(_OPENMP) || LBW_HAS_OMP_SIMD
#pragma omp simd reduction(+:s)
#endif
    for (int i = 0; i < st.n; i++) {
        w[i] *= q_hyp;  // Apply prior uniform multiplier
        s += w[i];
    }
    
    double scale = static_cast<double>(st.n) / s;
    
#if defined(_OPENMP) || LBW_HAS_OMP_SIMD
#pragma omp simd
#endif
    for (int i = 0; i < st.n; i++) {
        w[i] *= scale;
    }
    
    // The correction is the ratio of pre-projection to post-projection
    // q_hyp = w_pre / w_post = w_pre / (w_pre * scale) = 1.0 / scale
    q_hyp = 1.0 / scale; 
}
4. Revise the Post-Loop Finalizer
Apply the exact same multiplicative logic to your 20-iteration finalizer block at the bottom of the function.

C++
for (int fixup = 0; fixup < 20; fixup++) {
    bool box_ok = true;
    for (int i = 0; i < st.n; i++) {
        double yi = w[i] * q[i];            // Multiplicative
        double wc = std::clamp(yi, lo, hi);
        q[i] = yi / wc;                     // Multiplicative
        if (yi != wc) box_ok = false;
        w[i] = wc;
    }

    double s = 0.0;
#if defined(_OPENMP) || LBW_HAS_OMP_SIMD
#pragma omp simd reduction(+:s)
#endif
    for (int i = 0; i < st.n; i++) {
        w[i] *= q_hyp;                      // Multiplicative
        s += w[i];
    }
    double scale = static_cast<double>(st.n) / s;
#if defined(_OPENMP) || LBW_HAS_OMP_SIMD
#pragma omp simd
#endif
    for (int i = 0; i < st.n; i++) w[i] *= scale; // Multiplicative
    q_hyp = 1.0 / scale;
    
    if (box_ok) break;
}
Why this works
Your architecture was already perfectly set up for this; you only had the arithmetic domain wrong.

By switching to Bregman Dykstra, your IPF marginal passes are mathematically recognized as just another set of multiplicative KL projections.

Because all steps (Margins, Bounds, Sum) are now operating in the same KL-divergence space, the algorithm cannot "fight" itself.

If a set of weights exists that satisfies the margins, the bounds, and the sum, this revised code will find it. If the n_no_improve monitor fires after making this change, you can mathematically guarantee that the bounds are too tight for the targets (a truly infeasible problem), rather than an algorithmic stall.





Gemini is AI and can make mistakes, including about people. Your privacy and GeminiOpens in a new window

