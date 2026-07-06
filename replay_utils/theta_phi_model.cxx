#ifndef theta_phi_model_H
#define theta_phi_model_H


#include "replay_utils.h"
#include <ApexVDCTrack.h> 
#include <math.h> 
#include <ApexUtils.h> 


namespace replay_utils
{

double Phi_model(const ApexVDC::Track& track) { 
    
    return 
    -0.231*APEX::square(track.S2_y()) 
    +0.2532*track.S2_y() 
    +0.539e-3; 
};
double  Theta_model(const ApexVDC::Track& track) { 
    
    return 0.109648*track.S2_x(); 
};

};

#endif