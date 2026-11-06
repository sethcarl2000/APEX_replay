
#include <APEX/VDC/HitGroup.h> 
#include <stdexcept> 
#include <TString.h> 
#include <APEX/utils.h> 

namespace APEX
{
namespace VDC
{

//______________________________________________________________________________________________________________________________________
template<Plane P> HitGroup<P>::HitGroup(std::span<const Hit<P>> hits)
  : fis_RHRS{hits.begin()->Is_RHRS()}, fHits(hits.begin(), hits.end())
{}
//______________________________________________________________________________________________________________________________________
template<Plane P> double HitGroup<P>::WirePos( unsigned int h ) const 
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
template<Plane P> int    HitGroup<P>::WireNum( unsigned int h ) const 
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
template<Plane P> double HitGroup<P>::Time( unsigned int h )    const 
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
//______________________________________________________________________________________________________________________________________

};
}