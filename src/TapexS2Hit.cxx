
#include <TapexS2Hit.h>

#include <cmath> 
#include <limits> 

namespace {
    
    const Double_t Right_R[16] 
    = {2551.36, 2550.84, 2552.93, 2552.42, 2553.03, 2553.02, 2552.42, 2553.06, 2555.18, 2553.52, 2554.33, 2552.95, 2549.92, 2550.57, 2551.83, 2548.26};

    const Double_t Right_L[16] 
    = {2553.56, 2547.79, 2550.14, 2549.62, 2552.32, 2549.35, 2554.73, 2552.29, 2550.83, 2547.63, 2551.22, 2551.11, 2552.26, 2549.53, 2548.21, 2551.16};

    const Double_t Left_R[16] 
    = {2759.06, 2820.8, 2815.81, 2819.58, 2821.29, 2818.89, 2821.84, 2821.3, 2821.87, 2822.54, 2822.36, 2821.5, 2820.98, 2821.42, 2820.59, 2759.86};

    const Double_t Left_L[16] 
    = {2817.96, 2760.24, 2763.8, 2760.99, 2764.22, 2764.18, 2765.38, 2767.82, 2764.86, 2763.6, 2764.22, 2758.28, 2759.5, 2759.49, 2763.37, 2821.76};

  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN(); 

  inline bool is_nan(double x) { return x != x; }
}

////////////////////////////////////////////////////////////////////////////////////
TapexS2Hit::TapexS2Hit( bool arm, int paddle, double T_pmtL, double T_pmtR ) {
  
  f_isRightArm =arm; 
  fPaddle      =paddle; 
  
  fRawTime_pmtL =T_pmtL; 
  fRawTime_pmtR =T_pmtR; 
  
  fTime = Compute_RealTime(); 
  
  fIsCoinc = !is_nan(fTime); 

  //get the location of the paddle
  fZ = f_isRightArm ? 3.3098 : 3.1790 ; 
  
  fX = ((double)fPaddle)*fPaddleWidth_X  + (f_isRightArm ? -1.21413 : -1.16913); 
  
  fY = 0; 
}
double TapexS2Hit::DeltaT_raw() const
{ 

  return (fRawTime_pmtR - fRawTime_pmtL)/fTDC_resolution; 
}
double TapexS2Hit::Compute_RealTime() {
  
  //check to make sure that the PMTs registered non-null times
  if (std::fabs(fRawTime_pmtL) > 1e7 || 
      std::fabs(fRawTime_pmtR) > 1e7) return kNaN;  
  
  if (f_isRightArm) { //RHRS 
    
    fRealTime_pmtR = ( Right_R[fPaddle] - fRawTime_pmtR )*fTDC_resolution;
    fRealTime_pmtL = ( Right_L[fPaddle] - fRawTime_pmtL )*fTDC_resolution; 
  
  } else      { 
    
    fRealTime_pmtR = ( Left_R[fPaddle]  - fRawTime_pmtR )*fTDC_resolution; 
    fRealTime_pmtL = ( Left_L[fPaddle]  - fRawTime_pmtL )*fTDC_resolution; 
  }
    
  //check to make sure the raw times agree 
  if (std::fabs(fRealTime_pmtR-fRealTime_pmtL) > 7.5e-9) return kNaN; 
  
  return 0.5*(fRealTime_pmtR + fRealTime_pmtL); 
}
void TapexS2Hit::Make_twinHit( TapexS2Hit *neighbor ) { 
  
  fRawTime_pmtL = (this->fRawTime_pmtL + neighbor->fRawTime_pmtL)/2.;
  fRawTime_pmtR = (this->fRawTime_pmtR + neighbor->fRawTime_pmtR)/2.; 

  fTime = 0.5*(fTime + neighbor->fTime); 
  
  fX += 0.5*fPaddleWidth_X; 
  
  f_isTwinHit = true; 
} 