
#ifndef ApexVDC_H
#define ApexVDC_H

namespace ApexVDC 
{   
    /// number of VDC planes 
    constexpr int kN_planes = 4;

    /// number of wires per plane 
    constexpr int kN_wires_per_plane = 368; 

    /// resolution of TDC (s)
    constexpr double kTDC_resolution = 0.5e-9; 

    /// spacing between VDC wires (perp. distance) (m)
    constexpr double kWireSpacing = 4.2426e-3;

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
    constexpr double kW_RHRS[ApexVDC::kN_planes] = { 0, 0.026, 0.333245, 0.359245 };
    /// w-values (coordinate normal-vertical to plane) of each RHRS vdc plane (m)
    constexpr double kW_LHRS[ApexVDC::kN_planes] = { 0, 0.026, 0.335344, 0.361444 };

    /// @brief name of the namespace as a const char array
    constexpr char kNamespaceName[] = "ApexHRS";

};


#endif