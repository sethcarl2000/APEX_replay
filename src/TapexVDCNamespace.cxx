#include "TapexVDCNamespace.h"
#include <stdexcept> 
#include <limits> 
#include <sstream>
#include <cmath>  

namespace {
    //NaN to return in case function fails
    constexpr double ret_nan = std::numeric_limits<double>::quiet_NaN(); 

    constexpr double uvw_to_xyz[9] = {  
        +0.500000,  +0.500000, -0.707107, 
        -0.707107,  +0.707107,  0.,
        +0.500000,  +0.500000, +0.707107 
    }; 

    constexpr double xyz_to_uvw[9] = {  
        +0.500000, -0.707107,  +0.500000, 
        +0.500000, +0.707107,  +0.500000,
        -0.707107,  0.,        +0.707107 
    }; 
}

//___________________________________________________________________________________________________________________________
double VDC::WirePos(const bool is_RHRS, const int plane, const int wire_num) 
{
#ifdef RANGE_CHECKS
    //check that the plane and wire number are valid
    if ((plane < 0 || plane >= n_planes) || (wire_num < 0 || wire_num >= n_wires)) {
        std::ostringstream oss; 
        oss << "in <VDC::" << __func__ << ">: Plane (" << plane << ") or wire number (" << wire_num << ") is invalid."; 
        throw std::invalid_argument(oss.str()); 
        return ret_nan; 
    }
#endif 
    //now, return the wire position. note that wire 0 is at the highest u/v position, and the last wire is at the lowest 
    return (is_RHRS ? R_plane_wire0[plane] : L_plane_wire0[plane]) - wire_spacing*((double)wire_num);
}
//___________________________________________________________________________________________________________________________
int VDC::WireNum(const bool is_RHRS, const int plane, const double position)
{
#ifdef RANGE_CHECKS
    //check that the plane and wire number are valid
    if (plane < 0 || plane >= n_planes) {
        std::ostringstream oss; 
        oss << "in <VDC::" << __func__ << ">: Plane (" << plane << ") is invalid."; 
        throw std::invalid_argument(oss.str()); 
        return -1; 
    }
#endif 
    double offset = is_RHRS ? R_plane_wire0[plane] : L_plane_wire0[plane];
    int wire_num = (int)std::round((offset - position)/wire_spacing);
    
    if (wire_num > n_wires) { return n_wires; }
    else { if (wire_num < 0) return 0; }

    return wire_num; 
}
//___________________________________________________________________________________________________________________________
TVector3 VDC::UVW_to_XYZ(const TVector3& X)
{
    return TVector3(
        X[0]*uvw_to_xyz[0] + X[1]*uvw_to_xyz[1] + X[2]*uvw_to_xyz[2], 
        X[0]*uvw_to_xyz[3] + X[1]*uvw_to_xyz[4] + X[2]*uvw_to_xyz[5], 
        X[0]*uvw_to_xyz[6] + X[1]*uvw_to_xyz[7] + X[2]*uvw_to_xyz[8] 
    ); 
}
//___________________________________________________________________________________________________________________________
TVector3 VDC::XYZ_to_UVW(const TVector3& X)
{
    return TVector3(
        X[0]*xyz_to_uvw[0] + X[1]*xyz_to_uvw[1] + X[2]*xyz_to_uvw[2], 
        X[0]*xyz_to_uvw[3] + X[1]*xyz_to_uvw[4] + X[2]*xyz_to_uvw[5], 
        X[0]*xyz_to_uvw[6] + X[1]*xyz_to_uvw[7] + X[2]*xyz_to_uvw[8] 
    ); 
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
