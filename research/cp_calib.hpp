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
  std::vector<double> trace_data;
};

CPResult cp_calibrate(
    const int n_row,
    const int n_col,
    const int* p,
    const int* j,
    const double* x,
    const double* b,
    const double* d,
    const double* lo,
    const double* hi,
    const int max_iterations,
    const bool capture_trace
);

#endif
