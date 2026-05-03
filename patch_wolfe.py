import re
import sys
import os
import tempfile
import shutil

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
