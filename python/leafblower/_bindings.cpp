#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "leafblower.h"
#include <vector>
#include <string>
#include <cstring>

namespace py = pybind11;

// GIL is held throughout rk_calibrate() call.
// py_log_trampoline casts ctx to PyObject* (a callable) and invokes it.
static void py_log_trampoline(const char* msg, void* ctx) {
    py::gil_scoped_acquire gil;
    PyObject* callable = reinterpret_cast<PyObject*>(ctx);
    if (callable && PyCallable_Check(callable)) {
        PyObject* str = PyUnicode_FromString(msg);
        if (str) {
            PyObject* res = PyObject_CallOneArg(callable, str);
            if (!res) PyErr_Clear();  // log failure must not abort calibration
            Py_XDECREF(res);
            Py_DECREF(str);
        }
    }
}

PYBIND11_MODULE(_leafblower, m) {
    m.doc() = "leafblower C API bindings";

    m.def("calibrate",
        [](int n, int K,
           py::array_t<double, py::array::c_style | py::array::forcecast> weights_np,
           std::vector<py::array_t<int32_t, py::array::c_style | py::array::forcecast>> group_ids_list,
           std::vector<int> cat_counts_list,
           std::vector<py::array_t<double, py::array::c_style | py::array::forcecast>> targets_list,
           py::dict params_dict,
           py::object log_callable)
        -> py::tuple {
            // Validate dtype and layout
            if (weights_np.dtype().kind() != 'f' || weights_np.itemsize() != 8)
                throw py::value_error("weights must be float64");
            if (weights_np.ndim() != 1 || weights_np.size() != n)
                throw py::value_error("weights must be 1D array of length n");

            // Copy weights (so we can return a new array without aliasing)
            std::vector<double> weights_copy(weights_np.data(), weights_np.data() + n);

            // Build group_ids pointers
            std::vector<const int32_t*> gid_ptrs(K);
            // group_ids_list elements outlive rk_calibrate() — raw pointer capture is safe.
            // forcecast may allocate a temporary; that temporary lives in group_ids_list[k].
            for (int k = 0; k < K; k++) {
                auto& arr = group_ids_list[k];
                if (arr.size() != (size_t)n)
                    throw py::value_error("group_ids[k] length must equal n");
                gid_ptrs[k] = arr.data();
            }

            // Build targets pointers
            std::vector<const double*> tgt_ptrs(K);
            for (int k = 0; k < K; k++) {
                tgt_ptrs[k] = targets_list[k].data();
            }

            // Build params
            rk_params_t p;
            rk_params_init(&p);
            if (params_dict.contains("min_weight"))
                p.min_weight = params_dict["min_weight"].cast<double>();
            if (params_dict.contains("max_weight"))
                p.max_weight = params_dict["max_weight"].cast<double>();
            if (params_dict.contains("inner_max_iter"))
                p.inner_max_iter = params_dict["inner_max_iter"].cast<int>();
            if (params_dict.contains("outer_max_iter"))
                p.outer_max_iter = params_dict["outer_max_iter"].cast<int>();
            if (params_dict.contains("tol_abs"))
                p.tol_abs = params_dict["tol_abs"].cast<double>();
            if (params_dict.contains("verbose"))
                p.verbose = params_dict["verbose"].cast<int>();
            if (params_dict.contains("algorithm"))
                p.algorithm = (rk_algorithm_t)params_dict["algorithm"].cast<int>();
            if (params_dict.contains("epsilon"))
                p.epsilon = params_dict["epsilon"].cast<double>();
            if (params_dict.contains("lbfgs_m"))
                p.lbfgs_m = params_dict["lbfgs_m"].cast<int>();
            if (params_dict.contains("bounds_mode"))
                p.bounds_mode = (rk_bounds_mode_t)params_dict["bounds_mode"].cast<int>();
            // Convergence config (WU-G)
            if (params_dict.contains("pct_tol"))
                p.pct_tol = params_dict["pct_tol"].cast<double>();
            if (params_dict.contains("absolute_tol"))
                p.absolute_tol = params_dict["absolute_tol"].cast<double>();
            if (params_dict.contains("metric"))
                p.metric = params_dict["metric"].cast<int>();
            if (params_dict.contains("rule"))
                p.rule = params_dict["rule"].cast<int>();
            if (params_dict.contains("stop_when"))
                p.stop_when = params_dict["stop_when"].cast<int>();
            // SOR config (WU-F)
            if (params_dict.contains("sor_enabled"))
                p.sor_enabled = params_dict["sor_enabled"].cast<int>();
            if (params_dict.contains("sor_auto"))
                p.sor_auto = params_dict["sor_auto"].cast<int>();
            if (params_dict.contains("sor_omega_init"))
                p.sor_omega_init = params_dict["sor_omega_init"].cast<double>();
            if (params_dict.contains("sor_omega_min"))
                p.sor_omega_min = params_dict["sor_omega_min"].cast<double>();
            if (params_dict.contains("sor_omega_fixed"))
                p.sor_omega_fixed = params_dict["sor_omega_fixed"].cast<double>();
            if (params_dict.contains("sor_burnin"))
                p.sor_burnin = params_dict["sor_burnin"].cast<int>();
            // Homotopy (PY-1)
            if (params_dict.contains("homotopy_levels"))
                p.homotopy.n_levels = params_dict["homotopy_levels"].cast<int>();
            if (params_dict.contains("homotopy_start_factor"))
                p.homotopy.start_factor = params_dict["homotopy_start_factor"].cast<double>();
            if (params_dict.contains("homotopy_end_factor"))
                p.homotopy.end_factor = params_dict["homotopy_end_factor"].cast<double>();
            if (params_dict.contains("homotopy_budget_p"))
                p.homotopy.budget_split_p = params_dict["homotopy_budget_p"].cast<double>();
            // Scheduler / eta (PY-1)
            if (params_dict.contains("scheduler"))
                p.scheduler = (rk_scheduler_t)params_dict["scheduler"].cast<int>();
            if (params_dict.contains("eta_mode"))
                p.eta_mode = (rk_eta_mode_t)params_dict["eta_mode"].cast<int>();
            if (params_dict.contains("eta_start"))
                p.eta_start = params_dict["eta_start"].cast<double>();
            if (params_dict.contains("eta_end"))
                p.eta_end = params_dict["eta_end"].cast<double>();
            if (params_dict.contains("eta_schedule_power"))
                p.eta_schedule_power = params_dict["eta_schedule_power"].cast<double>();
            // Method-specific (PY-1)
            if (params_dict.contains("capacity_penalty"))
                p.capacity_penalty = params_dict["capacity_penalty"].cast<double>();
            if (params_dict.contains("newton_tsvd_ratio"))
                p.newton_tsvd_ratio = params_dict["newton_tsvd_ratio"].cast<double>();
            // SRAA / ALM (PY-2)
            if (params_dict.contains("accelerate"))
                p.accelerate = params_dict["accelerate"].cast<int>();
            if (params_dict.contains("alm_penalty"))
                p.alm_penalty = params_dict["alm_penalty"].cast<double>();

            // Wire log callback if verbose and callable provided
            PyObject* callable_ptr = nullptr;
            if (!log_callable.is_none() && p.verbose > 0) {
                callable_ptr = log_callable.ptr();
                p.log_fn  = py_log_trampoline;
                p.log_ctx = callable_ptr;
            }

            rk_result_t result;
            int rc = rk_calibrate(n, K, weights_copy.data(),
                                  gid_ptrs.data(),
                                  cat_counts_list.data(),
                                  tgt_ptrs.data(),
                                  &p, &result);

            // Return (status, weights_out_copy, result_dict)
            // weights_out is a NEW ndarray — never a view into input
            py::array_t<double> weights_out(n);
            std::memcpy(weights_out.mutable_data(), weights_copy.data(), (size_t)n * sizeof(double));

            py::dict result_dict;
            result_dict["status"]           = result.status;
            result_dict["iterations"]       = result.iterations;
            result_dict["max_error"]        = result.max_error;
            result_dict["algorithm_used"]   = (int)result.algorithm_used;
            result_dict["message"]          = std::string(result.message);
            result_dict["n_bounds_violated"] = result.n_bounds_violated;
            result_dict["n_bounds_clamped"]  = result.n_bounds_clamped;
            // Extended result fields (WU-G)
            result_dict["mean_error"]          = result.mean_error;
            result_dict["kl"]                  = result.kl;
            result_dict["chi2"]                = result.chi2;
            result_dict["l1_weight_change"]    = result.l1_weight_change;
            result_dict["grake_norm"]          = result.grake_norm;
            result_dict["convergence_metric"]  = result.convergence_metric;
            result_dict["convergence_rule"]    = result.convergence_rule;
            result_dict["convergence_tol"]     = result.convergence_tol;
            result_dict["convergence_iter"]               = result.convergence_iter;
            result_dict["convergence_solver_objective"]   = result.convergence_solver_objective;
            result_dict["convergence_minimized_metric"]   = result.convergence_minimized_metric;
            result_dict["best_error"]          = result.best_error;
            result_dict["best_iter"]           = result.best_iter;
            result_dict["sor_min_omega"]       = result.sor_min_omega;
            result_dict["sor_n_damped"]        = result.sor_n_damped;

            return py::make_tuple(rc, weights_out, result_dict);
        },
        py::arg("n"), py::arg("K"), py::arg("weights"),
        py::arg("group_ids"), py::arg("cat_counts"), py::arg("targets"),
        py::arg("params") = py::dict(),
        py::arg("log_callable") = py::none(),
        "Calibrate survey weights in-place. Returns (status, weights_copy, result_dict)."
    );
}
