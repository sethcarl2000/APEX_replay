#ifndef APEX_VDC_HitCluster_H 
#define APEX_VDC_HitCluster_H 

#include <APEX/VDC.h> 
#include <APEX/VDC/HitGroup.h> 

namespace APEX
{
namespace VDC
{

class HitCluster {
public: 
  HitCluster(HitGroup *group=nullptr, 
	      const double intercept=0, 
	      const double eta=0); 

  virtual ~HitCluster() {}; 
  
  HitGroup* GetGroup(); 
    
  double Intercept() const { return fIntercept; }
  double Eta()       const { return fEta_score; }
  
 private: 
  //ptr to parent group 
  APEX::VDC::HitGroup *fGroup; 
  double fIntercept; 
  double fEta_score; 
  
};

}
}

#endif 