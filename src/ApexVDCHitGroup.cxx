
#include <ApexVDCHitGroup.h> 
#include <stdexcept> 
#include <TString.h> 
#include <ApexUtils.h> 

using APEX::kNaN_double;
using APEX::kNull_int; 

namespace ApexVDC
{

void HitGroup::AddHit( double wire, double rawTime ) 
{   
  fHits.push_back(Hit(fPlane,wire,rawTime)); 
}
double HitGroup::WirePos( unsigned int h ) const 
{   
  if (h>=Nhits()) {
    throw std::logic_error(Form("<%s::%s> hit index %u out-of-range; max is %u", 
        kNamespaceName, __func__, 
        Nhits()-1, h
    ));
    return kNaN_double; 
  }
  return fHits[h].wPos();
}
int    HitGroup::WireNum( unsigned int h ) const 
{ 
  if (h>=Nhits()) {
    throw std::logic_error(Form("<%s::%s> hit index %u out-of-range; max is %u", 
        kNamespaceName, __func__, 
        Nhits()-1, h
    ));
    return APEX::kNull_int; 
  }
  
  return fHits[h].wNum(); 
}
double HitGroup::Time( unsigned int h )    const 
{ 
  if (h>=Nhits()) {
    throw std::logic_error(Form("<%s::%s> hit index %u out-of-range; max is %u", 
        kNamespaceName, __func__, 
        Nhits()-1, h
    ));
    return kNaN_double; 
  }
  
  return fHits[h].Time(); 
}
int    HitGroup::FirstWire() const { return WireNum(0); }
double HitGroup::LoEdge()    const { return WirePos(Nhits()-1); }
double HitGroup::HiEdge()    const { return WirePos(0); }

};