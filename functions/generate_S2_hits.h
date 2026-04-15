#ifndef generate_S2_hits_H
#define generate_S2_hits_H

// APEX headers
#include <TapexS2Hit.h> 
#include <TapexEventHandler.h>
#include "../run_parameters.h"
// ROOT headers
#include <ROOT/RVec.hxx>
// stdlib headers
#include <cmath> 

ROOT::VecOps::RVec<TapexS2Hit> generate_S2_hits(
    const bool is_RHRS,
    const ROOT::RVec<double>& PMT_R, 
    const ROOT::RVec<double>& PMT_L )
{
    //generate a vector of all coinc s2-paddle hits
    ROOT::RVec<TapexS2Hit> coinc_hits{}; 
    
    int last_paddle(-999); 
    
    for (int p=0; p<TapexS2Hit::N_paddles(); p++) { 
      
        TapexS2Hit hit( is_RHRS, p, PMT_R[p], PMT_L[p] );
        
        //check if it's a coinc-hit
        if ( !hit.IsCoinc() ) continue; 
        
        //check if this hit is a 'twin hit' with a previous hit
        if (coinc_hits.size()>0) { 

            auto& last_hit = coinc_hits.back(); 
        
            if ( hit.Paddle() - last_hit.Paddle() ==1  &&
                std::fabs(hit.Time() - last_hit.Time()) < run_parameters::kCUT_S2_twinHit_timeErr ) { 
        
                //these two hits are 'twins', i.e., likely caused by the same particle 
                last_hit.Make_twinHit( &hit ); 
            } 
        }
        
        coinc_hits.push_back( hit );    
    }
    
    return coinc_hits; 
}


#endif