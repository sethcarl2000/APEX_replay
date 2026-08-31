#ifndef APEX_units_H
#define APEX_units_H

namespace APEX
{
namespace units 
{

    //length
    constexpr double m   = 1.; 
    constexpr double cm  = m/100.; 
    constexpr double mm  = m/1000.;

    //energy
    constexpr double GeV = 1.; 
    constexpr double MeV = GeV/1000.; 

    //time 
    constexpr double s   = 1.;
    constexpr double ns  = s/1e9; 

    //angles 
    constexpr double rad = 1.;
    constexpr double mrad = rad/1000.; 
    constexpr double deg  = rad * (2 * 3.14159265359 / 360.); 

}
}

#endif 