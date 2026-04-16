#ifndef TapexEventHandler_H 
#define TapexEventHandler_H 

#include <TapexS2Hit.h>
#include <ROOT/RVec.hxx>

class TapexEventHandler { 
  
 public: 
  /*** 
   *    Keeps event variables like beam current, drift-parameters, etc. 
   * 
   ***/ 
  
  TapexEventHandler( bool arm=true, 
		 double beamCurrent=0.,
		 int runNumber=-1, 
		 const TapexS2Hit *fHit_R=0, 
		 const TapexS2Hit *fHit_L=0  ); 
  
  virtual ~TapexEventHandler() {}; 
  
  void SetActiveArm(bool arm) { f_activeArm=arm; }
  
  bool ActiveArm() const { return f_activeArm; }
  
  const TapexS2Hit* GetS2Hit(bool is_RHRS) const { return is_RHRS ? fS2Hit_Right : fS2Hit_Left; }
    
  //un-blurred drift function
  double Drift_X( double tau, double slope, int derivative=0 ) const; 
  double Drift_T( double x,   double slope, int derivative=0 ) const;  

  int GenUniqueTrackID(bool is_RHRS) { return is_RHRS ? fUniqueTrackCounter_RHRS++ : fUniqueTrackCounter_LHRS++; }; 
    
  double GetBeamCurrent() const { return fBeamCurrent; }
  
  double Get_tauSigma() const; 
  
  bool Is_nullBeamCurrent() const { return f_isNullBeamCurrent; }
  
 private: 
  bool f_activeArm; 
  const TapexS2Hit *fS2Hit_Right; 
  const TapexS2Hit *fS2Hit_Left; 

  int fUniqueTrackCounter_RHRS{0};
  int fUniqueTrackCounter_LHRS{0};

  //un-blurred drift function
  double Drift_T_raw( const double x, 
		      const double *par, 
		      const int derivative=0 ) const;  
  
  
  double f_tauSigma; 
  
  int    fRunNumber; 
  double fBeamCurrent; 
  double fTimeStamp; 
  
  bool   f_isNullBeamCurrent;  //sometimes, the beam-current reading is null 

  //each TapexEventHandler instance must have its own copy of this; as these parameters
  // are dependent on beam current, and whether a patricular run is pre-or-post-
  // VDC-fix. For the Right arm, however, the parameters do not change. 
  double fParams_L[5][5]; 
  
  //linearly interpolate between points Y (each with x-values X)
  // this assumes the size of the array to be 
  double Interpolate( const double x, 
		      const double *X, 
		      const double *Y, 
		      const int nPts ) const; 
  
  ClassDef(TapexEventHandler,1);
}; 

#endif 