import re

with open('src/lbfgsb_solver.cpp', 'r') as f:
    content = f.read()

# Replace in wolfe_zoom
def patch_wolfe(code):
    pattern = r"""        for \(int i = 0; i < st\.n; i\+\+\) \{
            double Fi, Hi;
            if \(fn\.exponential\) \{ auto fh = fn\.FH\(u_work\[i\]\); Fi = fh\.F; Hi = fh\.H; \}
            else                \{ Fi = fn\.F_from_e\(e_vec\[i\]\); Hi = fn\.H_from_e\(e_vec\[i\], u_work\[i\]\); \}
            phi_trial -= d\[i\] \* Hi;
            slope    -= d\[i\] \* Fi \* du\[i\];
        \}"""
    replacement = """        double sum_w = 0.0;
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
    return re.sub(pattern, replacement, code)

content = patch_wolfe(content)

# Replace in wolfe_line_search (variables j instead of i)
def patch_wolfe_line_search(code):
    pattern = r"""        for \(int j = 0; j < st\.n; j\+\+\) \{
            double Fj, Hj;
            if \(fn\.exponential\) \{ auto fh = fn\.FH\(u_work\[j\]\); Fj = fh\.F; Hj = fh\.H; \}
            else                \{ Fj = fn\.F_from_e\(e_vec\[j\]\); Hj = fn\.H_from_e\(e_vec\[j\], u_work\[j\]\); \}
            phi_trial -= d\[j\] \* Hj;
            slope    -= d\[j\] \* Fj \* du\[j\];
        \}"""
    replacement = """        double sum_w = 0.0;
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
    return re.sub(pattern, replacement, code)

content = patch_wolfe_line_search(content)

with open('src/lbfgsb_solver.cpp', 'w') as f:
    f.write(content)
