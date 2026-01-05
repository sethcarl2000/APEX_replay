#ifndef TapexEventHandler_H 
#define TapexEventHandler_H

#include <TObject.h> 
#include "TapexS2Hit.h"

////////////////////////////////////////////////////////////////////////////////////
/***
 * @class TapexEventHandler
 * 
 * @brief Class which manages drift model and other global event info    
 * 
 * This class is responsible for providing access to the VDC drift model for both arms, 
 * as well as other global event info, like expected VDC wire-TDC rms values, as a function of beam current
 * 
 ***/
class TapexEventHandler : public TObject { 
  
 public: 
  
  TapexEventHandler( bool arm=true, 
		 double beamCurrent=0.,
		 int runNumber=-1, 
		 TapexS2Hit *fHit_R=0, 
		 TapexS2Hit *fHit_L=0  ); 
  
  ~TapexEventHandler() {}; 
  
  void SetActiveArm(bool arm) { f_activeArm=arm; }
  
  bool ActiveArm() const { return f_activeArm; }
  
  TapexS2Hit* GetS2Hit() { return f_activeArm ? fS2Hit_Right : fS2Hit_Left; }
    
  //un-blurred drift function
  double Drift_X( double tau, double slope, int derivative=0 ) const; 
  double Drift_T( double x,   double slope, int derivative=0 ) const;  
    
  double GetBeamCurrent() const { return fBeamCurrent; }
  
  double Get_tauSigma()   const; 
  
  bool Is_nullBeamCurrent() const { return f_isNullBeamCurrent; }
  
 private: 
  bool f_activeArm; 
  TapexS2Hit *fS2Hit_Right; 
  TapexS2Hit *fS2Hit_Left; 
  
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
		      const int nPts   ) const; 
  
  
  //we apply a gaussian blur to the drift model, to get rid of sharp edges (which would cause problems when
  // we take first- and second-derivates. so, we need to perform a numerical gaussian blur: 
  static constexpr int gausBlur_nSamples = 25; 
  static constexpr double gausBlur_sigma =0.20e-3; //this will be in units of mm

  static constexpr const double dP_spread = 1./((double)gausBlur_nSamples); 

  double z_spread[TapexEventHandler::gausBlur_nSamples]; 


  ClassDef(TapexEventHandler,0);
}; 
/////////////////////////////////////////////////////////////////////////////////

#endif 