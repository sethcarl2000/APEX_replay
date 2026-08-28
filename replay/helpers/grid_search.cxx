
// APEX 
#include <APEX/replay/helpers.h>
#include <APEX/EventCounter.h>
#include <APEX/VDC/HitGroup.h> 
#include <APEX/EventHandler.h> 
#include <APEX/Utils.h> 
#include <APEX/run_parameters.h>
// stdlib
#include <cmath> 
#include <cstdio> 

namespace APEX
{
namespace replay
{
namespace helpers
{

//#define DEBUG_GRID

/////////////////////////////////////////////////////////////////////////////
double grid_search( 
    const EventHandler& evt, 
    VDC::HitGroup& group, 
    double m1, //slope 
    double m2, 
    double x_Lo,
    double x_Hi, 
    double TAU_sigma, 
    double TAU_buffer) 
{
    using APEX::square;

#ifdef DEBUG_GRID
    printf("<%s> in body\n", __func__); 
#endif 

    double eta(0.); 
    
    double m_min = std::min( m1, m2 ); 
    double m_max = std::max( m1, m2 ); 
    
    double m_avg = 0.5*(m1 + m2); 
    
    double x_span = x_Hi - x_Lo; 
    
    double x_avg = (x_Lo + x_Hi)/2.; 

    #ifdef DEBUG_GRID
    printf("<%s> ready to loop over hits\n", __func__); 
#endif 

    //loop over all hits
    for (int h=0; h<group.Nhits(); h++) { 

#ifdef DEBUG_GRID
        printf("<%s> hit %i / %i\n", __func__, h, group.Nhits()-1); 
#endif 
        double x = group.WirePos(h); 
        
        double tau_Lo; 
        
        if ( std::fabs(x-x_avg) < x_span/2. ) { //this wire-pos is between x_Hi & _Lo
        
            tau_Lo 
                = evt.Drift_T( 0., m_avg ) 
                - TAU_buffer; 
            
        } else {
        
            tau_Lo 
                = evt.Drift_T( m_min*( x-(x<x_avg ? x_Lo : x_Hi) ), m_avg ) 
                - TAU_buffer; 
            
            //cut hits that are too far away
            if ( tau_Lo > run_parameters::kVDC_max_realTime ) continue; 
        }
        
        double tau_Hi 
            = evt.Drift_T( m_max*( x-(x<x_avg ? x_Hi : x_Lo) ), m_avg ) 
            + TAU_buffer; 
            
        /*cout << TString::Format(" dv %2.1f mm, DT = %3.1f ns", 
        (x-x_avg)*1e3, (tau_Hi-tau_Lo)*1e9 ) << endl; */ 
        
        double tau = group.Time(h); 
        
        //tau is below lowest guess
        if (tau < tau_Lo) { 
            eta += std::exp( -0.5*square((tau-tau_Lo)/TAU_sigma) ); 
            continue; 
        }

        //'goldilocks'-zone, between hi- & lo-predictions
        if (tau < tau_Hi) { 
            eta += 1.; 
            continue; 
        }
        
        //tau is above highest guess
        eta += std::exp( -0.5*square((tau-tau_Hi)/TAU_sigma) ); 

    }//for (int h=0; h<group->Nhits(); 

    return eta; 
};
//////////////////////////////////////////////////////////////////////////////

}
}
}