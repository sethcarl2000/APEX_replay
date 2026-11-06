#ifndef APEX_VDC_HitCluster_H 
#define APEX_VDC_HitCluster_H 

#include <APEX/ObjCollection.h> 
#include <APEX/VDC.h> 
#include <APEX/VDC/HitGroup.h>

namespace APEX
{
class S2Event; 

namespace VDC
{

template<VDC::Plane P> class HitCluster {
private: 

  friend class ObjCollection<HitCluster<P>>; 
  object_id id{0}; 

  //ref to parent group 
  HitGroup<P> fGroup; 

  double fIntercept; 
  double fEta_score;   

public: 

  HitCluster(
    const HitGroup<P>& group, 
    double intercept, 
    double eta
  ) : fGroup{group}, fIntercept{intercept}, fEta_score{eta} {}; 
  
  HitGroup<P>& GetGroup() { return fGroup; } 
  const HitGroup<P>& GetGroup() const { return fGroup; } 
  
  void   SetIntercept(double x) { fIntercept=x; }
  double GetIntercept() const { return fIntercept; }

  void    SetEta(double eta) { fEta_score=eta; }
  double  GetEta() const { return fEta_score; }

  bool Is_RHRS() const { return fGroup.Is_RHRS(); }

  // Get unique ID
  object_id ID() const { return id; }
  bool operator==(const HitCluster<P>& rhs) const { return rhs.id == id; }
};
// explicit instantiations 
template class HitCluster<Plane::U1>; 
template class HitCluster<Plane::V1>; 
template class HitCluster<Plane::U2>; 
template class HitCluster<Plane::V2>; 

}
}

#endif 