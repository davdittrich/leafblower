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
    PyObject* callable = reinterpret_cast<PyObject*>(ctx);
    if (callable && PyCallable_Check(callable)) {
        PyObject* str = PyUnicode_FromString(msg);
        if (str) {
            PyObject* res = PyObject_CallOneArg(callable, str);
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
                                  (const double**)tgt_ptrs.data(),
                                  &p, &result);

            // Return (status, weights_out_copy, result_dict)
            // weights_out is a NEW ndarray — never a view into input
            py::array_t<double> weights_out(n);
            std::memcpy(weights_out.mutable_data(), weights_copy.data(), n * sizeof(double));

            py::dict result_dict;
            result_dict["status"]         = result.status;
            result_dict["iterations"]     = result.iterations;
            result_dict["max_error"]      = result.max_error;
            result_dict["algorithm_used"] = (int)result.algorithm_used;
            result_dict["message"]        = std::string(result.message);

            return py::make_tuple(rc, weights_out, result_dict);
        },
        py::arg("n"), py::arg("K"), py::arg("weights"),
        py::arg("group_ids"), py::arg("cat_counts"), py::arg("targets"),
        py::arg("params") = py::dict(),
        py::arg("log_callable") = py::none(),
        "Calibrate survey weights in-place. Returns (status, weights_copy, result_dict)."
    );
}
