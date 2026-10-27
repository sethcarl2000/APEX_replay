#ifndef APEX_S2_h
#define APEX_S2_h

#include <APEX/units.h>
// stdlib 
#include <limits> 

namespace APEX
{
namespace S2 
{   
    /// some misc. information for the S2 scintillator

    // number of paddles
    constexpr int N_paddles = 16; 

    // TDC resolution
    constexpr double TDC_resolution = 0.5 * units::ns; 

    // paddle width (X-direction)
    constexpr double paddle_width_X = 13.975 * units::cm;

    /// @brief represents an invalid time measurement 
    const double nulltime = std::numeric_limits<double>::quiet_NaN(); 

    /// @brief Get X-position of paddle, in DET coordinate system. 
    /// @param is_RHRS 'true' for RHRS, 'false' for LHRS
    /// @param paddle paddle idx [0-15]
    /// @return x-position for paddle
    inline double get_paddle_X(bool is_RHRS, int paddle) {
        return ((double)paddle)*S2::paddle_width_X  + (is_RHRS ? -121.413*units::cm : -116.913*units::cm);
    }

    /// @return Y-position of S2. 
    inline double get_Y() { return 0.*units::m; }

    /// @brief Z-position of S2. 
    /// @param is_RHRS 'true' for RHRS, 'false' for LHRS
    /// @return z-position of S2 face
    inline double get_Z(bool is_RHRS) { return is_RHRS ? 3.3098*units::m : 3.1790*units::m; };

    /// @brief deterrmine whether a reported coinc. hit is actually a coinc hit
    /// @param is_RHRS 'true' if RHRS, 'false' if LHRS
    /// @param paddle paddle of S2 hit [0-15]
    /// @param rawtime_pmt_R rawtime reported the Right (Left) PMT TDC
    /// @return calibrated realtime  
    double compute_realtime_pmtR(bool is_RHRS, int paddle, double rawtime_pmt_R);

    /// @brief deterrmine whether a reported coinc. hit is actually a coinc hit
    /// @param is_RHRS 'true' if RHRS, 'false' if LHRS
    /// @param paddle paddle of S2 hit [0-15]
    /// @param rawitme_pmt_L rawtime reported the Right (Left) PMT TDC
    /// @return calibrated realtime  
    double compute_realtime_pmtL(bool is_RHRS, int paddle, double rawtime_pmt_L);

}
}


#endif