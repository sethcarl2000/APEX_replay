
// APEX
#include <APEX/replay/helpers.h>
#include <APEX/EventHandler.h>
#include <APEX/S2Hit.h> 
// ROOT
#include <ROOT/RVec.hxx>
// stdlib
#include <cmath> 

namespace APEX
{
namespace replay
{
namespace helpers
{

ROOT::RVec<EventHandler> gen_coinc_window_events(
    double dt_min, 
    double dt_max, 
    double beam_current, 
    unsigned int run_number,
    const ROOT::RVec<S2Hit>& R_s2_hits,
    const ROOT::RVec<S2Hit>& L_s2_hits
)
{   
    ROOT::RVec<EventHandler> coinc_events; 

    for (const auto& R_hit : R_s2_hits) {
        for (const auto& L_hit : L_s2_hits) {

            double dt = R_hit.Time() - L_hit.Time(); 

            if ( dt_min < dt && dt < dt_max ) 
                coinc_events.push_back(EventHandler(false, beam_current, run_number, &R_hit, &L_hit)); 
        }
    }

    return coinc_events; 
}   

}
}
} 
