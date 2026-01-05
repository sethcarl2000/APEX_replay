/***
 * @class TapexVDCHit.h
 * 
 * @brief Class which contains one VDC wire TDC hit   
 * 
 * This class contains and processes a single VDC wire hit
 * 
 ***/

#ifndef TapexVDCHit_H 
#define TapexVDCHit_H

#include <TObject.h> 
#include "TapexEventHandler.h"

class TapexVDCHit : public TObject { 
 
 public: 
  //TapexVDCHit(); 
  TapexVDCHit( int plane=-999, int wire=-1, double rawTime=-1e30, TapexEventHandler *event=nullptr );
  
  ~TapexVDCHit() {/*noop*/}; 
  
  //void FillHit(int plane, double wire, double rawTime); 
  double Time() const { return fRealTime; }
  double wPos() const { return fWirePos; } 
  int    wNum() const { return fWireNum; } 

  void SetArm(bool arm) { f_isRightArm=arm; }
  bool IsRightArm() const { return f_isRightArm; } 
  
  double GetRawTime() const { return fRawTime; }
  void   SetRawTime(const double rawTime); 
  
  double W()     const { return fW; } 
  
  int    Plane() const { return fPlane; } 
  
  static double WireSpacing() { return 4.2426e-3; } 
  
  static double RawTime( const bool   is_RightArm, 
			 const int    plane, 
			 const int    wireNum, 
			 const double realTime ); 
  
  static int    WireNum( const bool   is_RightArm, 
			 const int    plane, 
			 const double wirePos ); 
  
  static double WirePos( const bool   is_RightArm, 
			 const int    plane, 
			 const int    wireNum ); 

 private: 
  TapexEventHandler *fEvent; 
  int    fPlane;
  double fRawTime; 
  int    fWireNum; 
  double fWirePos; 
  double fRealTime;  
  bool   f_isRightArm=true; 
  
  double fW; 
  static constexpr double fWireSpacing = 4.2426e-3; 
  
  double GetRealTime( double rawTime ) const; 
  double GetWirePos( int wire )        const; 
  double GetWireNum( double pos  )     const; 
  
  ClassDef(TapexVDCHit,0);
}; 

#endif 