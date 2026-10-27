#ifndef APEX_S2Event_H 
#define APEX_S2Event_H 

#include <APEX/utils.h>
#include <APEX/units.h> 
#include <APEX/S2.h> 
#include <APEX/VDC.h>
#include <APEX/VDC/Hit.h>
#include <APEX/VDC/HitGroup.h>
#include <APEX/VDC/HitCluster.h>
#include <APEX/VDC/ChamberPair.h>
#include <APEX/VDC/Track.h>
#include <ROOT/RVec.hxx>
// stdlib
#include <vector>

namespace APEX 
{

/// @brief This class is a S2-centered event; all analysis processes (from raw data to track reconstruction) proceed around a single S2 hit (or twin hit.)
class S2Event { 
private: 

  bool fis_RHRS; 

  // the S2 hit this event is based around 
  S2Hit fS2Hit; 

  // list of all hit-groups for each plane
  std::vector<VDC::HitGroup> fHitGroups[VDC::kN_planes]; 

  // list of all hit clusters for each VDC plane
  std::vector<VDC::HitCluster> fHitClusters[VDC::kN_planes]; 

  // list of all chamber pairs for each VDC chamber 
  std::vector<VDC::ChamberPair> fChamberPairs_lo, fChamberPairs_hi; 

  // list of all VDC tracks
  std::vector<VDC::Track> fTracks; 

public: 

  S2Event(bool is_RHRS, S2Hit hit); 

  void make_vdc_hit_groups(
    int plane, 
    const ROOT::VecOps::RVec<double>& wire, 
    const ROOT::VecOps::RVec<double>& rawtime 
  );

  std::vector<VDC::HitGroup>& get_groups(int plane) { return fHitGroups[plane]; }
  
  std::vector<VDC::HitCluster>& get_hit_clusters(int plane) { return fHitClusters[plane]; }
  
  std::vector<VDC::ChamberPair>& get_chamber_pairs(bool is_lo_chamber) { return is_lo_chamber ? fChamberPairs_lo : fChamberPairs_hi; }
  
  std::vector<VDC::Track>& get_tracks() { return fTracks; } 

  //get the timing of this S2 event. 
  double S2Time() const; 

  bool Is_RHRS() const { return fis_RHRS; }

};
/////////////////////////////////////////////////////////////////////////////////

}

#endif 