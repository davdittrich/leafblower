#include "newton_calib.hpp"
#include "calib_dispatch.hpp"

namespace lbw {

NewtonCalibResult newton_calibrate(CalibState& st) {
    NewtonCalibResult res;
    res.base.status = RK_ERR_BADARG;  // Stub: newton_calib.cpp not yet implemented (N2)
    return res;
}

} // namespace lbw
