// Epic-J WU-1 stub
#ifndef LEAFBLOWER_RESEARCH_CP_CALIB_HPP_
#define LEAFBLOWER_RESEARCH_CP_CALIB_HPP_

#include <vector>
#include <string>

struct CPResult {
  std::vector<double> weights;
  int status_code;
  std::string status_msg;
  int iterations;
  double wall_time_ms;
};

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
);

#endif
