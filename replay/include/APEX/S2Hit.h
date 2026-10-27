#ifndef APEX_S2Hit_H 
#define APEX_S2Hit_H 

#include <APEX/utils.h>
#include <APEX/units.h> 
#include <APEX/S2.h> 
// stdlib
#include <cmath>
#include <optional>
#include <memory> 

namespace APEX 
{



////////////////////////////////////////////////////////////////////////////////////
class S2Hit { 
  
  bool   fis_RHRS; 
  int    fPaddle; 
  double fRawTime_pmtL; 
  double fRawTime_pmtR; 
  double fRealTime_pmtR;
  double fRealTime_pmtL;

  /// @brief If this is a twin hit, then this will be the twin. 
  std::shared_ptr<S2Hit> fMyTwin{nullptr}; 

 public: 

  //maximum difference between L/R pmt for this hit to be considered a coinc. hit 
  static constexpr double coinc_timing_cut = 7.5 * units::ns; 

  S2Hit( bool is_RHRS, int paddle, double T_pmtR, double T_pmtL ); 

  /// @return 'true' if this S2 hit is coinsidered a coinc. hit, 'false' if not. 
  bool   IsCoinc()      const { 
    // if either PMT time is NAN, return 'false'
    if ((fRealTime_pmtL != fRealTime_pmtL) || (fRealTime_pmtR != fRealTime_pmtR)) return false; 
  
    // if both PMT times are valid, reutrn 'true' if they're in each other's coinc-window 
    return std::fabs(fRealTime_pmtR - fRealTime_pmtL) < coinc_timing_cut; 
  }

  bool   Is_RightArm()  const { return fis_RHRS; } 
  double Time()         const { return 0.5*(fRealTime_pmtR + fRealTime_pmtL); } 
  int    Paddle()       const { return fPaddle; }

  double Time_PMT_R() const { return fRealTime_pmtR; }
  double Time_PMT_L() const { return fRealTime_pmtL; }

  double DeltaT_raw()   const { return (fRawTime_pmtR - fRawTime_pmtL)/S2::TDC_resolution; }
  
  double X()            const { return S2::get_paddle_X(fis_RHRS,fPaddle); }
  double Y()            const { return S2::get_Y(); }
  double Z()            const { return S2::get_Z(fis_RHRS); }
  
  double PaddleWidth() const { return S2::paddle_width_X; }
  
  //this hit is now a 'twin'-hit, i.e., both this paddle (and its neighbor) were 
  // likely triggerd by the same particle. merge it with its neighbor. 
  void AddTwin(S2Hit twin) {  
    fMyTwin = std::make_shared<S2Hit>(std::move(twin));    
  }

  S2Hit* GetTwin() {
    return fMyTwin.get(); 
  }

  bool Is_twinHit() const { return fMyTwin != nullptr; }
  
  //number of S2 segments
  static constexpr int N_paddles() { return S2::N_paddles; }
};
/////////////////////////////////////////////////////////////////////////////////

}

#endif 