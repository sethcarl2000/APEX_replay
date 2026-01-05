#include "TapexVDCChamberPair.h"
#include "TapexVDCNamespace.h"
#include <stdexcept> 
#include <sstream> 
#include <cmath> 

//____________________________________________________________________________________________________________________________
TapexVDCChamberPair::TapexVDCChamberPair(bool is_loChamber, TapexVDCHitCluster* clust_u, TapexVDCHitCluster* clust_v, int unique_id)
    : fUnique_ID{unique_id}
{
    //check if the pointer to either cluster is null
    if ((clust_u==nullptr) || (clust_v==nullptr)) {
        std::ostringstream oss; 
        oss << "in <TapexChamberPair::" << __func__ << ">: pointer to u or v cluster is null."; 
        throw std::invalid_argument(oss.str()); 
        return; 
    }

    //get the intercepts
    fu = clust_u->Intercept(); 
    fv = clust_v->Intercept(); 
}

//____________________________________________________________________________________________________________________________
double TapexVDCChamberPair::ClosestWirePos_Lo( double x ) const 
{   
  int plane = f_isLoChamber ? 0 : 2; 
  
  double plane_offset = fGroup_U->IsRightArm() ? VDC::R_plane_wire0[plane] : VDC::L_plane_wire0[plane]; 

  double wire_num = (double)std::ceil( (plane_offset - x)/VDC::wire_spacing );  
    
  return VDC::WirePos(fGroup_U->IsRightArm(), plane, wire_num);  
}
//____________________________________________________________________________________________________________________________
double TapexVDCChamberPair::ClosestWirePos( const double x ) const { 
  
  int plane = f_isLoChamber ? 0 : 2; 
  
  const bool is_RHRS = fGroup_U->IsRightArm(); 
  
  return VDC::WirePos(is_RHRS, plane, VDC::WireNum(is_RHRS, plane, x));  
}
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
//____________________________________________________________________________________________________________________________
