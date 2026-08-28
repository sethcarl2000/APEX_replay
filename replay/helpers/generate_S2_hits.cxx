
// APEX headers
#include <APEX/replay/helpers.h>
#include <APEX/S2Hit.h> 
#include <APEX/EventHandler.h>
#include <APEX/run_parameters.h>
// ROOT headers
#include <ROOT/RVec.hxx>
// stdlib headers
#include <cmath> 

namespace APEX
{
namespace replay
{
namespace helpers
{

ROOT::VecOps::RVec<S2Hit> generate_S2_hits(
    const bool is_RHRS,
    const ROOT::RVec<double>& PMT_R, 
    const ROOT::RVec<double>& PMT_L )
{
    //generate a vector of all coinc s2-paddle hits
    ROOT::RVec<S2Hit> coinc_hits{}; 
    
    int last_paddle(-999); 
    
    for (int p=0; p<S2Hit::N_paddles(); p++) { 
      
        S2Hit hit( is_RHRS, p, PMT_R[p], PMT_L[p] );
        
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

}
}
}