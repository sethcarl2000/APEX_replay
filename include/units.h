#ifndef units_H
#define units_H

namespace units {

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
};

#endif 