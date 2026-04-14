#ifndef ApexVDCHitCluster_H 
#define ApexVDCHitCluster_H 

#include <ApexVDC.h> 
#include <ApexVDCHitGroup.h> 

namespace ApexVDC
{

class HitCluster {
public: 
  HitCluster(ApexVDC::HitGroup *group=nullptr, 
	      const double intercept=0, 
	      const double eta=0); 
  
  HitGroup* GetGroup() { return fGroup; }
  
  double Intercept() const { return fIntercept; }
  double Eta()       const { return fEta_score; }
  
 private: 
  //ptr to parent group 
  ApexVDC::HitGroup *fGroup; 
  double fIntercept; 
  double fEta_score; 
  
  ClassDef(HitCluster,1); 
};

}; 

#endif 