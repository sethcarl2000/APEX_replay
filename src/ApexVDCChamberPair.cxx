
#include <ApexVDCChamberPair.h>
#include <algorithm> 

namespace ApexVDC
{

///////////////////////////////////////////////////////////////////////////////
ChamberPair::ChamberPair( bool is_loChamber,
			    double u,
			    double v,
			    ApexVDC::HitGroup *Group_U,
			    ApexVDC::HitGroup *Group_V, 
			    int unique_id ) {
  
  f_isLoChamber = is_loChamber; 
  fu = u; 
  fv = v; 
  fGroup_U = Group_U;
  fGroup_V = Group_V; 
  
  fUnique_ID = unique_id; 
}
ChamberPair::ChamberPair( bool is_loChamber, 
			    HitCluster *clust_u,
			    HitCluster *clust_v, 
			    int unique_id ) { 

  f_isLoChamber = is_loChamber; 
  
  fu = clust_u->Intercept(); 
  fv = clust_v->Intercept(); 
  
  fGroup_U = clust_v->GetGroup(); 
  fGroup_V = clust_u->GetGroup(); 
  
  fUnique_ID = unique_id; 
}
void ChamberPair::Remove_track(Track* track)
{
  std::remove_if( fTracks.begin(), fTracks.end(), [track](const Track* rhs){ return track==rhs; }); 
}


}; 