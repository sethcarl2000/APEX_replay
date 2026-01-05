#include "TapexVDCNamespace.h"
#include <stdexcept> 
#include <limits> 
#include <sstream>
#include <cmath>  

namespace {
    //NaN to return in case function fails
    constexpr double ret_nan = std::numeric_limits<double>::quiet_NaN(); 
}

//___________________________________________________________________________________________________________________________
double VDC::WirePos(const bool is_RHRS, const int plane, const int wire_num) 
{
    //check that the plane and wire number are valid
    if ((plane < 0 || plane >= n_planes) || (wire_num < 0 || wire_num >= n_wires)) {
        std::ostringstream oss; 
        oss << "in <VDC::" << __func__ << ">: Plane (" << plane << ") or wire number (" << wire_num << ") is invalid."; 
        throw std::invalid_argument(oss.str()); 
        return ret_nan; 
    }

    //now, return the wire position. note that wire 0 is at the highest u/v position, and the last wire is at the lowest 
    return (is_RHRS ? R_plane_wire0[plane] : L_plane_wire0[plane]) - wire_spacing*((double)wire_num);
}
//___________________________________________________________________________________________________________________________
int VDC::WireNum(const bool is_RHRS, const int plane, const double position)
{
    //check that the plane and wire number are valid
    if (plane < 0 || plane >= n_planes) {
        std::ostringstream oss; 
        oss << "in <VDC::" << __func__ << ">: Plane (" << plane << ") is invalid."; 
        throw std::invalid_argument(oss.str()); 
        return -1; 
    }

    double offset = is_RHRS ? R_plane_wire0[plane] : L_plane_wire0[plane];
    int wire_num = (int)std::round((offset - position)/wire_spacing);
    
    if (wire_num > n_wires) { return n_wires; }
    else { if (wire_num < 0) return 0; }

    return wire_num; 
}
//___________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________
