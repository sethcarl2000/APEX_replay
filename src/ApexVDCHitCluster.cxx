
#include <ApexVDCHitCluster.h> 
#include <stdexcept> 

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
//________________________________________________________________________________________________
HitGroup* HitCluster::GetGroup()
{
  if (fGroup) return fGroup; 

  throw std::logic_error("<ApexVDC::HitGroup::GetGroup>: requested ptr for group which was null."); 
  return nullptr; 
}
//________________________________________________________________________________________________
//________________________________________________________________________________________________
//________________________________________________________________________________________________
}; 

