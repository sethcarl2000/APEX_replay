#include "TapexS2Hit.h"
#include "TapexS2PMTOffsets.h"
#include <cmath> 

//____________________________________________________________________________________________________________________________
TapexS2Hit::TapexS2Hit( bool arm, int paddle, double T_pmtL, double T_pmtR ) {
  
  f_isRightArm =arm; 
  fPaddle      =paddle; 
  
  fRawTime_pmtL =T_pmtL; 
  fRawTime_pmtR =T_pmtR; 
  
  fTime = Compute_RealTime(); 
  
  fIsCoinc = (fTime > -1e29); 

  //get the location of the paddle
  fZ = f_isRightArm ? 3.3098 : 3.1790 ; 
  
  fX = ((double)fPaddle)*fPaddleWidth_X  + (f_isRightArm ? -1.21413 : -1.16913); 
  
  fY = 0; 
}
//____________________________________________________________________________________________________________________________
double TapexS2Hit::Compute_RealTime() {
  
  //check to make sure that the PMTs registered non-null times
  if (std::abs(fRawTime_pmtL) > 1e7 || 
      std::abs(fRawTime_pmtR) > 1e7) return -1e30; 
  
  if (f_isRightArm) { //RHRS 
    
    fRealTime_pmtR = ( S2_PMT::Right_R[fPaddle] - fRawTime_pmtR )*fTDC_resolution;
    fRealTime_pmtL = ( S2_PMT::Right_L[fPaddle] - fRawTime_pmtL )*fTDC_resolution; 
  
  } else      { 
    
    fRealTime_pmtR = ( S2_PMT::Left_R[fPaddle]  - fRawTime_pmtR )*fTDC_resolution; 
    fRealTime_pmtL = ( S2_PMT::Left_L[fPaddle]  - fRawTime_pmtL )*fTDC_resolution; 
  }
    
  //check to make sure the raw times agree 
  
  if (std::abs(fRealTime_pmtR-fRealTime_pmtL) > 7.5e-9) return -1e30; 
  
  return 0.5*(fRealTime_pmtR + fRealTime_pmtL); 
}
//____________________________________________________________________________________________________________________________
void TapexS2Hit::Make_twinHit( TapexS2Hit *neighbor ) { 
  
  fTime = 0.5*(fTime + neighbor->Time()); 
  
  fX += 0.5*fPaddleWidth_X; 
  
  f_isTwinHit = true; 
} 
//____________________________________________________________________________________________________________________________

ClassImp(TapexS2Hit); 