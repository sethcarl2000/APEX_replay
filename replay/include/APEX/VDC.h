
#ifndef APEX_VDC_H
#define APEX_VDC_H

#include <APEX/units.h>

namespace APEX
{
namespace VDC 
{   
    
/// number of VDC planes 
constexpr int kN_planes = 4;

/// number of wires per plane 
constexpr int kN_wires_per_plane = 368; 

/// resolution of TDC (s)
constexpr double kTDC_resolution = 0.5*units::ns; 

/// spacing between VDC wires (perp. distance) (m)
constexpr double kWireSpacing = 4.2426*units::mm;

/// @return real-time of VDC hit (s) 
double RealTime( const bool is_RHRS, 
            const int plane, 
            const int wireNum, 
            const double rawTime );

/// @return position of VDC wire in uvw coordinates (m)
double WirePos( const bool is_RHRS, 
            const int  plane, 
            const int  wireNum );

/// @return number of closest wire 
int WireNum( const bool   is_RHRS, 
            const int    plane, 
            const double x );

/// @return w-position (vertical) of the plane (m)
double w( const bool is_RHRS, const int plane );

/// w-values (coordinate normal-vertical to plane) of each RHRS vdc plane (m)
constexpr double kW_RHRS[kN_planes] = { 0.*units::m, 0.026*units::m, 0.333245*units::m, 0.359245*units::m };
/// w-values (coordinate normal-vertical to plane) of each RHRS vdc plane (m)
constexpr double kW_LHRS[kN_planes] = { 0.*units::m, 0.026*units::m, 0.335344*units::m, 0.361444*units::m };

/// @brief name of the namespace as a const char array
constexpr char kNamespaceName[] = "APEX::VDC";

};
}

#endif