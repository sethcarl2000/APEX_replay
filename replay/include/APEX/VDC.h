
#ifndef APEX_VDC_H
#define APEX_VDC_H

#include <APEX/replay.h>
#include <APEX/units.h>
// stdlib
#include <array>
#include <string> 

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

// VDC chamber (low or high)
enum class Chamber { Low, High }; 

// VDC plane
enum class Plane { U1, V1, U2, V2 }; 


/// @tparam VDC Chamber C  
/// @return 'true' if low-chamber (Chamber::Low) 'false' otherwise (Chamber::High)
template<Chamber C> inline constexpr bool is_loChamber() { return C == Chamber::Low; } 

/// @tparam P VDC plane 
/// @return 'true' if low-chamber (U1 || V1) else otherwise (U2 || V2)
template<Plane P> inline constexpr bool is_loChamber() {
    if constexpr ((P == Plane::U1) || (P == Plane::V1)) { return true; } else { return false; } 
}

template<Chamber C> inline constexpr Plane get_U() {
    if constexpr (is_loChamber<C>()) return Plane::U1; 
    return Plane::U2; 
}
template<Chamber C> inline constexpr Plane get_V() {
    if constexpr (is_loChamber<C>()) return Plane::V1; 
    return Plane::V2; 
}

// template expression that can be called with any VDC::Plane as a template argument
template<typename F> concept PlaneCallable = 
requires(F& f) {
    f.template operator()<Plane::U1>();
    f.template operator()<Plane::V1>();
    f.template operator()<Plane::U2>();
    f.template operator()<Plane::V2>();
} || requires(const F& f) {
    f.template operator()<Plane::U1>();
    f.template operator()<Plane::V1>();
    f.template operator()<Plane::U2>();
    f.template operator()<Plane::V2>();
}; 

// Call the expression 'f' once for each VDC plane. 
template<PlaneCallable F> void for_each_plane(F& f) {
    f.template operator()<Plane::U1>();
    f.template operator()<Plane::V1>();
    f.template operator()<Plane::U2>();
    f.template operator()<Plane::V2>();
};
// Call the expression 'f' once for each VDC plane. 
template<PlaneCallable F> void for_each_plane(const F& f) {
    f.template operator()<Plane::U1>();
    f.template operator()<Plane::V1>();
    f.template operator()<Plane::U2>();
    f.template operator()<Plane::V2>();
};

// template expression that can be called with any VDC::Plane as a template argument
template<typename F> concept ChamberCallable = requires(F& f) {
    f.template operator()<Chamber::Low>();
    f.template operator()<Chamber::High>();
} || requires(const F& f) {
    f.template operator()<Chamber::Low>();
    f.template operator()<Chamber::High>();
}; 

// Call the expression 'f' once for each VDC chamber. 
template<ChamberCallable F> void for_each_chamber(F& f) {
    f.template operator()<Chamber::Low>();
    f.template operator()<Chamber::High>();
};
// Call the expression 'f' once for each VDC chamber. 
template<ChamberCallable F> void for_each_chamber(const F& f) {
    f.template operator()<Chamber::Low>();
    f.template operator()<Chamber::High>();
};

// number of VDC plane 
template<Plane P> inline constexpr int number() {
    if constexpr (P == Plane::U1) return 0; 
    if constexpr (P == Plane::V1) return 1; 
    if constexpr (P == Plane::U2) return 2; 
    /* only option left is V2 */  return 3; 
}

/// @brief Get Names of each plane 
/// @tparam P VDC Plane
/// @return std::string_view of plane name
template<Plane P> std::string plane_name() {
    static const std::array<std::string_view, VDC::kN_planes> names{
        "U1", "V1", "U2", "V2"
    }; 
    return std::string{names[VDC::number<P>()]}; 
}

/// @brief Get Names of each Chamber 
/// @tparam C VDC Chamber
/// @return std::string_view of chamber name
template<Chamber C> std::string chamber_name() {
    if constexpr (C == Chamber::Low) { return "Low"; } else { return "High"; }
}

template<Plane P> inline constexpr Chamber get_Chamber() {
    if constexpr (VDC::is_loChamber<P>()) { return Chamber::Low; } else { return Chamber::High; }
}

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

template<Plane P> inline double w(bool is_RHRS) {
    return w(is_RHRS, VDC::number<P>() );
}

/// w-values (coordinate normal-vertical to plane) of each RHRS vdc plane (m)
constexpr double kW_RHRS[kN_planes] = { 0.*units::m, 0.026*units::m, 0.333245*units::m, 0.359245*units::m };
/// w-values (coordinate normal-vertical to plane) of each RHRS vdc plane (m)
constexpr double kW_LHRS[kN_planes] = { 0.*units::m, 0.026*units::m, 0.335344*units::m, 0.361444*units::m };

/// @brief name of the namespace as a const char array
constexpr char kNamespaceName[] = "APEX::VDC";


}
}

#endif