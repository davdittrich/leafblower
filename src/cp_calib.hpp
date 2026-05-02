#ifndef LEAFBLOWER_CP_CALIB_HPP_
#define LEAFBLOWER_CP_CALIB_HPP_

#include <string>
#include "types.hpp"

namespace lbw {

struct CpCalibResult {
    CalibResult base;
    int         n_cells              = 0;
    std::string algorithm_requested;
    std::string algorithm_used;
    double      A_norm_estimate      = 0.0;
    int         n_power_iter         = 0;
    double      final_theta          = 0.0;
    double      final_tau            = 0.0;
    double      final_sigma          = 0.0;
    bool        fell_back_to_pdhg    = false;
    double      wall_time_ms         = 0.0;
};

CpCalibResult cp_calibrate(CalibState& st);

} // namespace lbw

#endif
