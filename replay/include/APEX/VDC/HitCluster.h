#ifndef APEX_VDC_HitCluster_H 
#define APEX_VDC_HitCluster_H 

#include <APEX/ObjRef.h>
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
  unique_id id{0}; 

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

  // Update the ptr to the S2Event, in case it's moved in memory
  void SetS2Event(S2Event* evt) { fEvent=evt; } 

  bool Is_RHRS() const { return fGroup.Is_RHRS(); }

  // Get unique ID
  unique_id ID() const { return id; }
  bool operator==(const HitCluster<P>& rhs) const { return rhs.id == id; }
};

}
}

#endif 