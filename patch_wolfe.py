import re
import sys
import os
import tempfile
import shutil

# F3 derivation trace (verified 2026-05-03 against phi_from_u canonical reference):
#
#   phi(lam) = T·lam − Σ_i d_i·H(u_i) + lam_reg·r(lam) + (mu/2)·r(lam)^2
#   where r(lam) = Σ_i d_i·F(u_i) − n     (ALM constraint residual)
#         F = H'                          (logit.hpp: H is antiderivative of F)
#         u_i(lam)  = Σ_k lam[off[k]+g_k(i)]
#         du_i/dα  = du[i]                (precomputed via compute_du)
#
# Directional derivative along search dir:
#   dr/dα      = Σ_i d_i · F'(u_i) · du[i] = Σ_i d_i · fn.dF(u_i) · du[i]   ≡ sum_dw
#   d/dα[lam_reg·r + (mu/2)·r^2] = (lam_reg + mu·r) · dr/dα = alm_scale · sum_dw
#
# Sign: phi is MAXIMIZED (slope_0 > 0). The non-ALM loop computes
#   slope -= d[i] * F_i * du[i]   (negative because phi has − Σ d·H term)
# The ALM penalty is ADDED to phi (line 87 of phi_from_u: obj += λ·r + (μ/2)·r²),
# so its directional derivative ADDS to slope:
#   slope += alm_scale * sum_dw
# This matches phi_from_u's gradient assembly at line 84:
#   grad[off[k]+g] += alm_scale * d[i] * fn.dF(u[i])
# (directional derivative = grad · dir collapses to alm_scale · Σ d_i · dF · du_i).
# Verdict: arithmetic CORRECT in sign and units.

target_file = 'src/lbfgsb_solver.cpp'

with open(target_file, 'r') as f:
    content = f.read()

# Replace in wolfe_zoom (loop variable: i)
wolfe_zoom_pattern = r"""        for \(int i = 0; i < st\.n; i\+\+\) \{
            double Fi, Hi;
            if \(fn\.exponential\) \{ auto fh = fn\.FH\(u_work\[i\]\); Fi = fh\.F; Hi = fh\.H; \}
            else                \{ Fi = fn\.F_from_e\(e_vec\[i\]\); Hi = fn\.H_from_e\(e_vec\[i\], u_work\[i\]\); \}
            phi_trial -= d\[i\] \* Hi;
            slope    -= d\[i\] \* Fi \* du\[i\];
        \}"""
wolfe_zoom_replacement = """        double sum_w = 0.0;
        double sum_dw = 0.0;
        for (int i = 0; i < st.n; i++) {
            double Fi, Hi;
            if (fn.exponential) { auto fh = fn.FH(u_work[i]); Fi = fh.F; Hi = fh.H; }
            else                { Fi = fn.F_from_e(e_vec[i]); Hi = fn.H_from_e(e_vec[i], u_work[i]); }
            phi_trial -= d[i] * Hi;
            slope    -= d[i] * Fi * du[i];

            if (st.alm_mu > 0.0 || st.alm_lambda != 0.0) {
                sum_w  += d[i] * Fi;
                sum_dw += d[i] * fn.dF(u_work[i]) * du[i];
            }
        }
        if (st.alm_mu > 0.0 || st.alm_lambda != 0.0) {
            double residual = sum_w - static_cast<double>(st.n);
            double alm_scale = st.alm_lambda + st.alm_mu * residual;
            phi_trial += st.alm_lambda * residual + (st.alm_mu / 2.0) * residual * residual;
            slope += alm_scale * sum_dw;
        }"""

content, n1 = re.subn(wolfe_zoom_pattern, wolfe_zoom_replacement, content)

if n1 == 0:
    print(f"ERROR: wolfe_zoom pattern not found in {target_file} — already patched or source changed",
          file=sys.stderr)
    sys.exit(1)
if n1 > 1:
    print(f"ERROR: wolfe_zoom pattern matched {n1} times — ambiguous patch, refusing to apply",
          file=sys.stderr)
    sys.exit(1)

# Replace in wolfe_line_search (loop variable: j)
wolfe_ls_pattern = r"""        for \(int j = 0; j < st\.n; j\+\+\) \{
            double Fj, Hj;
            if \(fn\.exponential\) \{ auto fh = fn\.FH\(u_work\[j\]\); Fj = fh\.F; Hj = fh\.H; \}
            else                \{ Fj = fn\.F_from_e\(e_vec\[j\]\); Hj = fn\.H_from_e\(e_vec\[j\], u_work\[j\]\); \}
            phi_trial -= d\[j\] \* Hj;
            slope    -= d\[j\] \* Fj \* du\[j\];
        \}"""
wolfe_ls_replacement = """        double sum_w = 0.0;
        double sum_dw = 0.0;
        for (int j = 0; j < st.n; j++) {
            double Fj, Hj;
            if (fn.exponential) { auto fh = fn.FH(u_work[j]); Fj = fh.F; Hj = fh.H; }
            else                { Fj = fn.F_from_e(e_vec[j]); Hj = fn.H_from_e(e_vec[j], u_work[j]); }
            phi_trial -= d[j] * Hj;
            slope    -= d[j] * Fj * du[j];

            if (st.alm_mu > 0.0 || st.alm_lambda != 0.0) {
                sum_w  += d[j] * Fj;
                sum_dw += d[j] * fn.dF(u_work[j]) * du[j];
            }
        }
        if (st.alm_mu > 0.0 || st.alm_lambda != 0.0) {
            double residual = sum_w - static_cast<double>(st.n);
            double alm_scale = st.alm_lambda + st.alm_mu * residual;
            phi_trial += st.alm_lambda * residual + (st.alm_mu / 2.0) * residual * residual;
            slope += alm_scale * sum_dw;
        }"""

content, n2 = re.subn(wolfe_ls_pattern, wolfe_ls_replacement, content)

if n2 == 0:
    print(f"ERROR: wolfe_line_search pattern not found in {target_file} — already patched or source changed",
          file=sys.stderr)
    sys.exit(1)
if n2 > 1:
    print(f"ERROR: wolfe_line_search pattern matched {n2} times — ambiguous patch, refusing to apply",
          file=sys.stderr)
    sys.exit(1)

# Atomic write via temp file + rename
dirpath = os.path.dirname(os.path.abspath(target_file))
with tempfile.NamedTemporaryFile(mode='w', dir=dirpath, delete=False, suffix='.tmp') as tmp:
    tmp.write(content)
    tmp_path = tmp.name

shutil.move(tmp_path, target_file)
print(f"Patched {target_file} (2 substitutions: wolfe_zoom + wolfe_line_search)")
