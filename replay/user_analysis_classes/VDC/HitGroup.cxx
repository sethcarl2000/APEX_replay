
#include <APEX/VDC/HitGroup.h> 
#include <stdexcept> 
#include <TString.h> 
#include <APEX/utils.h> 

namespace APEX
{
namespace VDC
{

//______________________________________________________________________________________________________________________________________
HitGroup::HitGroup(std::span<Hit> hits)
  : fis_RHRS{hits.begin()->Is_RHRS()}, fPlane{hits.begin()->Plane()}, fHits(hits.begin(), hits.end())
{}
//______________________________________________________________________________________________________________________________________
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
//______________________________________________________________________________________________________________________________________
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
//______________________________________________________________________________________________________________________________________
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
//______________________________________________________________________________________________________________________________________
int    HitGroup::FirstWire() const { return fHits.front().WireNum(); }
//______________________________________________________________________________________________________________________________________
double HitGroup::LoEdge()    const { return fHits.back().WirePos(); }
//______________________________________________________________________________________________________________________________________
double HitGroup::HiEdge()    const { return fHits.front().WirePos(); }
//______________________________________________________________________________________________________________________________________

};
}