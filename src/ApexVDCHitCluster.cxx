
#include <ApexVDCHitCluster.h> 

namespace ApexVDC
{

////////////////////////////////////////////////////////////////////////////////  
HitCluster::HitCluster( ApexVDC::HitGroup *group, 
			  const double intercept, 
			  const double eta ) { 
  fGroup     =group; 
  fIntercept =intercept; 
  fEta_score =eta; 
}

}; 

ClassImp(ApexVDC::HitCluster);
