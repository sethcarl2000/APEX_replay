
#include <APEX/S2Hit.h>

#include <cmath> 
#include <array>
#include <limits> 


namespace APEX
{

namespace {
    
  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN(); 

  constexpr std::array<double,S2::N_paddles> pmt_dt_offsets_RHRS{
    -5.42,
    -6.15,
    -5.86,
    -6.06,
    -1.97,
    -7.94,
    3.99,
    -3.40,
    -9.27,
    -12.52,
    -7.21,
    -5.01,
    4.21,
    -4.30,
    -9.40,
    -8.33
  };

  constexpr std::array<double, S2::N_paddles> pmt_dt_offsets_LHRS{
      -115.00,
      -121.29,
      -104.60,
      -117.71,
      -114.71,
      -110.07,
      -113.13,
      -104.64,
      -113.92,
      -117.92,
      -116.11,
      -126.19,
      -122.52,
      -123.05,
      -113.68,
      -115.00
  };

  // this is the effective 'velocity' of scintillation photons in the S2 paddles, in units of cm/TDC-step 
  constexpr double v_effective = (1.2)*units::cm; 

  inline bool is_nan(double x) { return x != x; }
}

////////////////////////////////////////////////////////////////////////////////////
S2Hit::S2Hit( bool is_RHRS, int paddle, double T_pmtR, double T_pmtL ) 
  : fis_RHRS{is_RHRS}, fPaddle{paddle}, fRawTime_pmtR{T_pmtR}, fRawTime_pmtL{T_pmtL}
{
  fRealTime_pmtR = S2::compute_realtime_pmtR(fis_RHRS, fPaddle, fRawTime_pmtR); 
  fRealTime_pmtL = S2::compute_realtime_pmtL(fis_RHRS, fPaddle, fRawTime_pmtL); 
}
//__________________________________________________________________________________
double S2Hit::Y() const 
{
  double pmt_diff = fRawTime_pmtR - fRawTime_pmtL;
 
  double offset = Is_RightArm() ? pmt_dt_offsets_RHRS[Paddle()] : pmt_dt_offsets_LHRS[Paddle()]; 
 
  return (pmt_diff - offset)/v_effective; 
}

}