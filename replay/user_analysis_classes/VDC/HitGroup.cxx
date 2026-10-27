
#include <APEX/VDC/HitGroup.h> 
#include <stdexcept> 
#include <TString.h> 
#include <APEX/utils.h> 

namespace APEX
{
namespace VDC
{

double HitGroup::WirePos( unsigned int h ) const 
{   
  if (h>=Nhits()) {
    throw std::logic_error(Form("<%s::%s> hit index %u out-of-range; max is %u", 
        kNamespaceName, __func__, 
        Nhits()-1, h
    ));
    return utils::kNaN; 
  }
  return fHits[h].WirePos();
}
int    HitGroup::WireNum( unsigned int h ) const 
{ 
  if (h>=Nhits()) {
    throw std::logic_error(Form("<%s::%s> hit index %u out-of-range; max is %u", 
        kNamespaceName, __func__, 
        Nhits()-1, h
    ));
    return -1; 
  }
  
  return fHits[h].WireNum(); 
}
double HitGroup::Time( unsigned int h )    const 
{ 
  if (h>=Nhits()) {
    throw std::logic_error(Form("<%s::%s> hit index %u out-of-range; max is %u", 
        kNamespaceName, __func__, 
        Nhits()-1, h
    ));
    return utils::kNaN; 
  }
  
  return fHits[h].Time(); 
}
int    HitGroup::FirstWire() const { return fHits.front().WireNum(); }
double HitGroup::LoEdge()    const { return fHits.back().WirePos(); }
double HitGroup::HiEdge()    const { return fHits.front().WirePos(); }

};
}