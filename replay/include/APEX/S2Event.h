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
#include <APEX/EventHandler.h>
#include <ROOT/RVec.hxx>
// stdlib
#include <vector>
#include <memory> 

namespace APEX 
{

/// @brief This class is a S2-centered event; all analysis processes (from raw data to track reconstruction) proceed around a single S2 hit (or twin hit.)
class S2Event { 
private: 

  bool fis_RHRS; 

  // the S2 hit this event is based around 
  S2Hit fS2Hit; 

  // id of this event, unique per-arm 
  std::size_t id{0};

  // we wrap this event handler to pass it to some legacy functions that use it 
  std::shared_ptr<EventHandler> fEventHandler{nullptr}; 

  // list of all hit-groups for each plane
  ROOT::VecOps::RVec<VDC::HitGroup> fHitGroups[VDC::kN_planes]; 

  // list of all hit clusters for each VDC plane
  ROOT::VecOps::RVec<VDC::HitCluster> fHitClusters[VDC::kN_planes]; 

  // list of all chamber pairs for each VDC chamber 
  ROOT::VecOps::RVec<VDC::ChamberPair> fChamberPairs_lo, fChamberPairs_hi; 

  // list of all VDC tracks
  ROOT::VecOps::RVec<VDC::Track> fTracks; 

public: 

  S2Event(bool is_RHRS, S2Hit hit, const EventHandler* event_handler, std::size_t _id=0); 

  void make_vdc_hit_groups(
    int plane, 
    const ROOT::VecOps::RVec<double>& wire, 
    const ROOT::VecOps::RVec<double>& rawtime 
  );

  ROOT::VecOps::RVec<VDC::HitGroup>& get_groups(int plane) { return fHitGroups[plane]; }
  
  ROOT::VecOps::RVec<VDC::HitCluster>& get_hit_clusters(int plane) { return fHitClusters[plane]; }
  
  ROOT::VecOps::RVec<VDC::ChamberPair>& get_chamber_pairs(bool is_lo_chamber) { return is_lo_chamber ? fChamberPairs_lo : fChamberPairs_hi; }
  
  ROOT::VecOps::RVec<VDC::Track>& get_tracks() { return fTracks; } 

  //get the timing of this S2 event. 
  double S2Time() const; 

  /// @return unique ID of this event. unique among all events _for this arm_ but _not_ between arms (i.e, both arms will have an event with 'id = 0')
  std::size_t GetID() const { return id; }

  EventHandler& GetEventHandler() { return *fEventHandler.get(); } 

  bool Is_RHRS() const { return fis_RHRS; }

};
/////////////////////////////////////////////////////////////////////////////////

}

#endif 