#ifndef APEX_VDC_ChamberPair_H 
#define APEX_VDC_ChamberPair_H 

#include <APEX/ObjCollection.h>
#include <APEX/VDC.h> 
#include <APEX/VDC/HitCluster.h> 

namespace APEX
{
class S2Event; 

namespace VDC 
{

template<VDC::Chamber C> class ChamberPair { 
private: 

  friend class ObjCollection<ChamberPair<C>>; 
  object_id id{0}; 

  using ClustU = HitCluster<VDC::get_U<C>()>; 
  using ClustV = HitCluster<VDC::get_V<C>()>; 

  // groups this cluster is identified with 
  ClustU fClust_U; 
  ClustV fClust_V; 

public: 
 
  ChamberPair( 
        const ClustU& clust_U, 
        const ClustV& clust_V
  ) : fClust_U{clust_U}, fClust_V{clust_V} {}; 
  
  double u() const { return fClust_V.GetIntercept(); }
  double v() const { return fClust_U.GetIntercept(); }
  
  void Set_u(double u) { fClust_V.SetIntercept(u); }
  void Set_v(double v) { fClust_U.SetIntercept(v); }
  
  void Get_uv(double &u, double &v) {
    u=fClust_V.GetIntercept(); 
    v=fClust_U.GetIntercept(); 
  };
  
  void Set_uv(double u, double v) {
    fClust_V.SetIntercept(u); 
    fClust_U.SetIntercept(v); 
  };

  ClustU& GetCluster_U() { return fClust_U; }
  ClustV& GetCluster_V() { return fClust_V; }

  const ClustU& GetCluster_U() const { return fClust_U; }
  const ClustV& GetCluster_V() const { return fClust_V; }
  
  bool Is_loChamber() const { return VDC::is_loChamber<C>(); }

  bool Is_RHRS() const { return fClust_U.Is_RHRS(); }
  
  // Get unique ID
  object_id ID() const { return id; }
  bool operator==(const ChamberPair<C>& rhs) const { return rhs.id == id; }
};
// explicit instantiations 
template class ChamberPair<Chamber::Low>; 
template class ChamberPair<Chamber::High>; 

}; 
} 


#endif 