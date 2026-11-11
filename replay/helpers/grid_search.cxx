
// APEX 
#include <APEX/replay/helpers.h>
#include <APEX/VDC/HitGroup.h> 
#include <APEX/EventHandler.h> 
#include <APEX/S2Event.h>
#include <APEX/utils.h> 
#include <APEX/units.h>
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
template<VDC::Plane P> double grid_search( 
    S2Event& evt, 
    const VDC::HitGroup<P>& group, 
    double m1, //slope 
    double m2, 
    double x_Lo,
    double x_Hi, 
    double TAU_sigma, 
    double TAU_buffer, 
    double TAU_offset
) 
{
    double eta{0.}; 
    
    double m_min = std::min( m1, m2 ); 
    double m_max = std::max( m1, m2 ); 
    
    double m_avg = 0.5*(m1 + m2); 
    
    double x_span = x_Hi - x_Lo; 
    
    double x_avg = (x_Lo + x_Hi)/2.; 

#ifdef DEBUG_GRID
    printf("<%s> ready to loop over hits\n", __func__); 
#endif 

    //loop over all hits
    for (const auto& hit : group) { 

        double x = hit.WirePos();
        
        //find min / max drift-distance
        double w_min = m_min * std::min( std::fabs(x-x_Lo), std::fabs(x-x_Hi) ); 
        double w_max = m_max * std::max( std::fabs(x-x_Lo), std::fabs(x-x_Hi) ); 
 
#ifdef DEBUG_GRID
        std::printf(" w_min %+3.1f mm, w_max: %+3.1f mm",
            w_min/units::mm, w_max/units::mm    
        ); 
#endif
        // skip any hits above a max. reasonable drift-distance 
        if (w_min > run_parameters::max_drift_distance) {
#ifdef DEBUG_GRID
        std::printf(" out-of-range\n"); 
#endif        
            continue; 
        }
        double tau_min = evt.Drift_T(w_min, m_avg, 0) + TAU_offset - TAU_buffer; 
        double tau_max = evt.Drift_T(w_max, m_avg, 0) + TAU_offset + TAU_buffer; 
        
        double tau = hit.Time(); 
        
        double d_eta; 
        //tau is below lowest guess
        if (tau_min < tau && tau < tau_max) { 
            d_eta = 1.; 
        } else {
            d_eta = std::exp( -0.5*utils::intpow<2>( 
                std::min( 
                    std::fabs(tau-tau_min), 
                    std::fabs(tau-tau_max) 
                )/TAU_sigma 
            )); 
        }

        eta += d_eta; 

#ifdef DEBUG_GRID
        std::printf(" tau: %+6.1f, tau_min: %+6.1f,  tau_max: %+6.1f, d_eta %.3f\n",
            tau/units::ns, tau_min/units::ns, tau_max/units::ns, d_eta
        ); 
#endif
    }

    return eta; 
};
//////////////////////////////////////////////////////////////////////////////
// explicit template instantiations 
template double grid_search( 
    S2Event&, 
    const VDC::HitGroup<VDC::Plane::U1>&, 
    double, double, double, double, double, double, double
); 
template double grid_search( 
    S2Event&, 
    const VDC::HitGroup<VDC::Plane::V1>&, 
    double, double, double, double, double, double, double
); 
template double grid_search( 
    S2Event&, 
    const VDC::HitGroup<VDC::Plane::U2>&, 
    double, double, double, double, double, double, double
); 
template double grid_search( 
    S2Event&, 
    const VDC::HitGroup<VDC::Plane::V2>&, 
    double, double, double, double, double, double, double
); 

}
}
}