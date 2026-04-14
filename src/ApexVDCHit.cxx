#include <ApexVDCHit.h> 
#include <stdexcept> 
#include <TString.h> 

namespace ApexVDC {

Hit::Hit( int plane, 
		  double wire, 
		  double rawTime, 
		  TapexEventHandler *event ) { 
  
  fEvent = event; 
  
  f_isRightArm = fEvent->ActiveArm(); 
  
  fPlane    = plane; 
  fRawTime  = rawTime; 
  fWireNum  = (int)wire; 
  fWirePos  = ApexVDC::WirePos( f_isRightArm, fPlane, fWireNum );
  fRealTime = ApexVDC::RealTime( f_isRightArm, fPlane, fWireNum, fRawTime ); 
    
  fW = ApexVDC::w(f_isRightArm, fPlane);
}
void   Hit::SetRawTime(const double rawTime)  { 
  
  fRawTime  = rawTime; 
  fRealTime = ApexVDC::RealTime( f_isRightArm, fPlane, fWireNum, fRawTime ) - fEvent->GetS2Hit()->Time(); 
} 

}; 
////////////////////////////////////////////////////////////////////////////////  