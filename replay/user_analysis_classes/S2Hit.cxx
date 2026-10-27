
#include <APEX/S2Hit.h>

#include <cmath> 
#include <limits> 

namespace {
    
    
  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN(); 

  inline bool is_nan(double x) { return x != x; }
}

namespace APEX
{

////////////////////////////////////////////////////////////////////////////////////
S2Hit::S2Hit( bool is_RHRS, int paddle, double T_pmtR, double T_pmtL ) 
  : fis_RHRS{is_RHRS}, fPaddle{paddle}, fRawTime_pmtR{T_pmtR}, fRawTime_pmtL{T_pmtL}
{
  fRealTime_pmtR = S2::compute_realtime_pmtR(fis_RHRS, fPaddle, fRawTime_pmtR); 
  fRealTime_pmtL = S2::compute_realtime_pmtL(fis_RHRS, fPaddle, fRawTime_pmtL); 
}

}