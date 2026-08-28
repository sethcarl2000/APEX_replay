#ifndef APEX_VDC_Hit_h
#define APEX_VDC_Hit_h

#include <APEX/EventHandler.h> 
#include <APEX/VDC.h> 

namespace APEX 
{
namespace VDC
{

/////////////////////////////////////////////////////////////////////////////////
class Hit { 

 public: 
  //Hit(); 
  Hit( int plane=-999, 
	   double wire=-1e30, 
	   double rawTime=-1e30, 
	   const EventHandler *event=0 );
    
  virtual ~Hit() {}; 

  Hit& operator=(const Hit& rhs) = default; 

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
  const EventHandler *fEvent; 
  int    fPlane;
  double fRawTime; 
  int    fWireNum; 
  double fWirePos; 
  double fRealTime;  
  bool   f_isRightArm=true; 
  
  double fW; 
  
}; 

}
}

#endif