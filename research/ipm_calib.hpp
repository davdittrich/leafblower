// Epic-J WU-5: Interior-Point Newton (IPM) calibration
#ifndef LEAFBLOWER_RESEARCH_IPM_CALIB_HPP_
#define LEAFBLOWER_RESEARCH_IPM_CALIB_HPP_

#include <vector>
#include <string>

struct IPMResult {
  std::vector<double> weights;
  int status_code;
  std::string status_msg;
  int iterations;       // total inner Newton steps performed
  double wall_time_ms;
  std::vector<double> trace_data;  // 8 columns per row, see Sec 3.2 diagnostic table
};

IPMResult ipm_calibrate(
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
