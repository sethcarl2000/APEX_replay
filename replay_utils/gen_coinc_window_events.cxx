#ifndef gen_coinc_events_H
#define gen_coinc_events_H

#include "replay_utils.h"
#include <TapexEventHandler.h>
#include <TapexS2Hit.h> 
#include <ROOT/RVec.hxx>
#include <math.h> 

namespace replay_utils
{

ROOT::RVec<TapexEventHandler> gen_coinc_window_events(
    double dt_min, 
    double dt_max, 
    double beam_current, 
    unsigned int run_number,
    const ROOT::RVec<TapexS2Hit>& R_s2_hits,
    const ROOT::RVec<TapexS2Hit>& L_s2_hits
)
{   
    ROOT::RVec<TapexEventHandler> coinc_events; 

    for (const auto& R_hit : R_s2_hits) {
        for (const auto& L_hit : L_s2_hits) {

            double dt = R_hit.Time() - L_hit.Time(); 

            if ( dt_min < dt && dt < dt_max ) 
                coinc_events.push_back(TapexEventHandler(false, beam_current, run_number, &R_hit, &L_hit)); 
        }
    }

    return coinc_events; 
}   

}; 

#endif