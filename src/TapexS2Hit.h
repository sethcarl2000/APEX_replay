#ifndef TapexS2Hit_H 
#define TapexS2Hit_H 

#include <TObject.h> 

////////////////////////////////////////////////////////////////////////////////////
class TapexS2Hit { 
  
 public: 
  /***
   *      Tracks information for S2m hits of either arm, including the ability 
   *      to convert between raw/realtime. 
   ***/
  TapexS2Hit( bool arm=true, int paddle=-1, double T_pmtL=-1e30, double T_pmtR=-1e30 ); 
  
  virtual ~TapexS2Hit() {/*noop*/}; 
  
  bool   IsCoinc()      const { return fIsCoinc; }
  bool   Is_RightArm()  const { return f_isRightArm; } 
  double Time()         const { return fTime; }
  int    Paddle()       const { return fPaddle; }
  
  double DeltaT_raw()   const 
  { return (fRawTime_pmtR - fRawTime_pmtL)/fTDC_resolution; }
  
  double X()            const { return fX; }
  double Y()            const { return fY; }
  double Z()            const { return fZ; }
  
  double PaddleWidth() const { return fPaddleWidth_X; }
  
  
  //this hit is now a 'twin'-hit, i.e., both this paddle (and its neighbor) were 
  // likely triggerd by the same particle. merge it with its neighbor. 
  void Make_twinHit(TapexS2Hit *neighbor); 

  bool Is_twinHit() const { return f_isTwinHit; }
  
  //number of S2 segments
  static constexpr int N_paddles() { return 16; }
  
 private: 
  bool   f_isRightArm; 
  int    fPaddle; 
  bool   fIsCoinc; 
  double fTime;
  double fRawTime_pmtL; 
  double fRawTime_pmtR; 
  double fRealTime_pmtR;
  double fRealTime_pmtL;
  
  bool f_isTwinHit=false; 
  
  double fX, fY, fZ; 

  
  static constexpr double fPaddleWidth_X = 0.13975; //in m
  
  static constexpr double fTDC_resolution = 0.5e-9; 
  
  double Compute_RealTime(); //const; 
  
  ClassDef(TapexS2Hit,1); 
};
/////////////////////////////////////////////////////////////////////////////////

#endif 