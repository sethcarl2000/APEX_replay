
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

};


#endif