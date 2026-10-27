#include <APEX/S2.h>
// stdlib
#include <cmath> 

namespace APEX 
{
namespace S2
{

constexpr double RHRS_PMT_R[16] 
= {2553.56, 2547.79, 2550.14, 2549.62, 2552.32, 2549.35, 2554.73, 2552.29, 2550.83, 2547.63, 2551.22, 2551.11, 2552.26, 2549.53, 2548.21, 2551.16};

constexpr double RHRS_PMT_L[16] 
= {2551.36, 2550.84, 2552.93, 2552.42, 2553.03, 2553.02, 2552.42, 2553.06, 2555.18, 2553.52, 2554.33, 2552.95, 2549.92, 2550.57, 2551.83, 2548.26};

constexpr double LHRS_PMT_R[16] 
= {2817.96, 2760.24, 2763.8, 2760.99, 2764.22, 2764.18, 2765.38, 2767.82, 2764.86, 2763.6, 2764.22, 2758.28, 2759.5, 2759.49, 2763.37, 2821.76};

constexpr double LHRS_PMT_L[16] 
= {2759.06, 2820.8, 2815.81, 2819.58, 2821.29, 2818.89, 2821.84, 2821.3, 2821.87, 2822.54, 2822.36, 2821.5, 2820.98, 2821.42, 2820.59, 2759.86};

//__________________________________________________________________________________________________________________
double compute_realtime_pmtR(bool is_RHRS, int paddle, double pmt_rawtime)
{
    // this is beyond the upper limit for the maximum reasonable PMT time. 
    if (std::fabs(pmt_rawtime) > 1e5) return S2::nulltime; 

    double shift = is_RHRS ? RHRS_PMT_R[paddle] : LHRS_PMT_R[paddle]; 

    // we subtract the rawtime to reflect the fact that the TDC counts 'up' from the trigger, going backwards in time. 
    // i.e., a singal arriving close time to the trigger is a small rawtime, and a signal arriving long before the trigger 
    // is a large rawtime. 
    return ( shift - pmt_rawtime )*TDC_resolution; 
}
//__________________________________________________________________________________________________________________
double compute_realtime_pmtL(bool is_RHRS, int paddle, double pmt_rawtime)
{    
    // this is beyond the upper limit for the maximum reasonable PMT time. 
    if (std::fabs(pmt_rawtime) > 1e5) return S2::nulltime; 
    
    double shift = is_RHRS ? RHRS_PMT_L[paddle] : LHRS_PMT_L[paddle]; 

    // we subtract the rawtime, to reflect the fact that the TDC counts 'up' from the trigger, going backwards in time. 
    // i.e., a singal arriving close time to the trigger is a small rawtime, and a signal arriving long before the trigger 
    // is a large rawtime. 
    return ( shift - pmt_rawtime )*TDC_resolution; 
}
//__________________________________________________________________________________________________________________

}
}