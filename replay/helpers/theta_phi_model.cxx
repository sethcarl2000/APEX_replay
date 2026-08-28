
#include <APEX/replay/helpers.h>
#include <APEX/VDC/Track.h> 
#include <math.h> 
#include <APEX/utils.h> 


namespace APEX
{
namespace replay
{
namespace helpers
{

double Phi_model(const VDC::Track& track) { 
    
    return 
    -0.231*utils::intpow<2>(track.S2_y()) 
    +0.2532*track.S2_y() 
    +0.539e-3; 
};
double  Theta_model(const VDC::Track& track) { 
    
    return 0.109648*track.S2_x(); 
};

}
}
}