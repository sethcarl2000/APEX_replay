#include <APEX/VDC/Hit.h> 

namespace APEX
{
namespace VDC 
{

////////////////////////////////////////////////////////////////////////////////  
Hit::Hit( int plane, 
		  double wire, 
		  double rawTime, 
		  const EventHandler *event ) { 
  
  fEvent = event; 
  
  f_isRightArm = fEvent->ActiveArm(); 
  
  fPlane    = plane; 
  fRawTime  = rawTime; 
  fWireNum  = (int)wire; 
  fWirePos  = VDC::WirePos( f_isRightArm, fPlane, fWireNum );
  fRealTime = VDC::RealTime( f_isRightArm, fPlane, fWireNum, fRawTime ); 
    
  fW = VDC::w(f_isRightArm, fPlane);
}
////////////////////////////////////////////////////////////////////////////////  
void   Hit::SetRawTime(const double rawTime)  { 
  
  fRawTime  = rawTime; 
  fRealTime = VDC::RealTime( f_isRightArm, fPlane, fWireNum, fRawTime ) - fEvent->GetS2Hit(f_isRightArm)->Time(); 
} 

} 
}
////////////////////////////////////////////////////////////////////////////////  