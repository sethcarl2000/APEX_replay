
#include <APEX/VDC/HitCluster.h> 
#include <stdexcept> 

namespace APEX
{
namespace VDC
{

////////////////////////////////////////////////////////////////////////////////  
HitCluster::HitCluster( HitGroup *group, 
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

  throw std::logic_error("<APEX::VDC::HitGroup::GetGroup>: requested ptr for group which was null."); 
  return nullptr; 
}
//________________________________________________________________________________________________
//________________________________________________________________________________________________
//________________________________________________________________________________________________

}; 
}
