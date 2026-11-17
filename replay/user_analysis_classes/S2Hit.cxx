
#include <APEX/S2Hit.h>
#include <APEX/S2.h>

#include <cmath> 
#include <array>
#include <limits> 


namespace APEX
{

namespace {
    
  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN(); 

  constexpr std::array<double,S2::N_paddles> pmt_dt_offsets_RHRS{
    -2.52,
    -2.90,
    -2.88,
    -2.95,
    -1.00,
    -3.93,
    1.97,
    -1.73,
    -4.68,
    -6.27,
    -3.64,
    -2.54,
    1.99,
    -2.26,
    -4.77,
    -4.23
  };

  constexpr std::array<double, S2::N_paddles> pmt_dt_offsets_LHRS{
    0.00,
    -60.64,
    -52.44,
    -58.91,
    -57.40,
    -55.11,
    -56.60,
    -52.37,
    -57.00,
    -58.97,
    -58.08,
    -63.06,
    -61.23,
    -61.44,
    -56.79,
    0.00
  };

  // this is the effective 'velocity' of scintillation photons in the S2 paddles, in units of cm/TDC-step 
  constexpr double v_effective = (4.664)*units::cm/units::ns; 

  inline bool is_nan(double x) { return x != x; }
}

////////////////////////////////////////////////////////////////////////////////////
S2Hit::S2Hit( bool is_RHRS, int paddle, double T_pmtR, double T_pmtL, double ADC_R, double ADC_L ) 
  : fis_RHRS{is_RHRS}, fPaddle{paddle}, fRawTime_pmtR{T_pmtR}, fRawTime_pmtL{T_pmtL}, fADC_pmtR{ADC_R}, fADC_pmtL{ADC_L}
{
  fRealTime_pmtR = S2::compute_realtime_pmtR(fis_RHRS, fPaddle, fRawTime_pmtR); 
  fRealTime_pmtL = S2::compute_realtime_pmtL(fis_RHRS, fPaddle, fRawTime_pmtL); 
}
//__________________________________________________________________________________
double S2Hit::Y() const 
{
  double pmt_diff = fRawTime_pmtR - fRawTime_pmtL;
 
  double offset = fis_RHRS ? pmt_dt_offsets_RHRS[Paddle()] : pmt_dt_offsets_LHRS[Paddle()]; 
 
  return (pmt_diff - offset)*S2::TDC_resolution * v_effective; 
}

}