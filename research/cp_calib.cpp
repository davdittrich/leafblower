// Epic-J WU-1 stub
#include "cp_calib.hpp"

CPResult cp_calibrate(
    const int n_row,
    const int n_col,
    const std::vector<int>& p,
    const std::vector<int>& j,
    const std::vector<double>& x,
    const std::vector<double>& b,
    const std::vector<double>& d,
    const std::vector<double>& lo,
    const std::vector<double>& hi,
    const int max_iterations,
    const bool capture_trace
) {
  CPResult res;
  res.weights.assign(n_row, 1.0);
  res.status_code = 99;
  res.status_msg = "WU-1 stub";
  res.iterations = 0;
  res.wall_time_ms = 0.0;
  return res;
}
