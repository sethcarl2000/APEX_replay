#ifndef APEX_Cerenkov_h
#define APEX_Cerenkov_h

#include <APEX/units.h>
// stdlib 
#include <limits> 

namespace APEX
{
namespace Cerenkov
{   
    /// some misc. information for the S2 scintillator

    // number of paddles
    constexpr int N_pmts = 10; 

    // TDC resolution
    constexpr double TDC_resolution = 0.5 * units::ns; 

    /// @brief represents an invalid time measurement 
    const double nulltime = std::numeric_limits<double>::quiet_NaN(); 
 
    /// @return Z-position of Cerenkov PMT mirrors, TRANSPORT coordinates 
    inline double Z(bool is_RHRS) {
        return is_RHRS ? 180.*units::cm : 181.*units::cm; 
    }

}
}


#endif