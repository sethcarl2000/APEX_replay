#ifndef ApexVDCHit_H
#define ApexVDCHit_H

#include <TapexEventHandler.h> 
#include <ApexVDC.h> 

namespace ApexVDC 
{

/////////////////////////////////////////////////////////////////////////////////
class Hit { 

 public: 
  //Hit(); 
  Hit( int plane=-999, 
	   double wire=-1e30, 
	   double rawTime=-1e30, 
	   TapexEventHandler *event=0 );
  
  ~Hit() {/*noop*/}; 
  
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
  

 private: 
  TapexEventHandler *fEvent; 
  int    fPlane;
  double fRawTime; 
  int    fWireNum; 
  double fWirePos; 
  double fRealTime;  
  bool   f_isRightArm=true; 
  
  double fW; 
  
  ClassDef(Hit,0);
}; 

};

#endif