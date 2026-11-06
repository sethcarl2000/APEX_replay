#ifndef APEX_S2Event_H 
#define APEX_S2Event_H 

#include <APEX/utils.h>
#include <APEX/ObjCollection.h>
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
#include <cstddef> 
#include <limits> 
#include <map> 
#include <stdexcept> 
#include <tuple> 
#include <concepts> 


namespace APEX 
{

// concept where type must be one of the analysis classes
template<typename T> concept AnaClass = 
  std::same_as<T, VDC::HitGroup<VDC::Plane::U1>> || // VDC Hit Groups 
  std::same_as<T, VDC::HitGroup<VDC::Plane::V1>> || 
  std::same_as<T, VDC::HitGroup<VDC::Plane::U2>> || 
  std::same_as<T, VDC::HitGroup<VDC::Plane::V2>> || 
  std::same_as<T, VDC::HitCluster<VDC::Plane::U1>> || // VDC Hit Clusters
  std::same_as<T, VDC::HitCluster<VDC::Plane::V1>> || 
  std::same_as<T, VDC::HitCluster<VDC::Plane::U2>> || 
  std::same_as<T, VDC::HitCluster<VDC::Plane::V2>> || 
  std::same_as<T, VDC::ChamberPair<VDC::Chamber::Low>>  || // VDC Chamber Pairs
  std::same_as<T, VDC::ChamberPair<VDC::Chamber::High>> || 
  std::same_as<T, VDC::Track>; // VDC Track
  

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

  std::tuple<

    ObjCollection<VDC::HitGroup<VDC::Plane::U1>>,  // U1 Hit groups
    ObjCollection<VDC::HitGroup<VDC::Plane::V1>>,  // V1 
    ObjCollection<VDC::HitGroup<VDC::Plane::U2>>,  // U2 
    ObjCollection<VDC::HitGroup<VDC::Plane::V2>>,  // V2 
  
    ObjCollection<VDC::HitCluster<VDC::Plane::U1>>, // U1 Hit clusters
    ObjCollection<VDC::HitCluster<VDC::Plane::V1>>, // V1 
    ObjCollection<VDC::HitCluster<VDC::Plane::U2>>, // U2 
    ObjCollection<VDC::HitCluster<VDC::Plane::V2>>, // V2 

    ObjCollection<VDC::ChamberPair<VDC::Chamber::Low>>,  // Low-chamber Chamber pairs
    ObjCollection<VDC::ChamberPair<VDC::Chamber::High>>, // Hi-chamber 
    
    ObjCollection<VDC::Track> //tracks 

  > fObjects; 

public: 

  S2Event(bool is_RHRS, S2Hit hit, const EventHandler* event_handler, std::size_t _id=0); 

  /// @return Ref to S2 hit this event is based around
  S2Hit& GetS2Hit() { return fS2Hit; }
  
  /// @return Ref to S2 hit this event is based around
  const S2Hit& GetS2Hit() const { return fS2Hit; }

  /// @brief Return size of current event's collection of a given type 
  /// @tparam T Type of object 
  /// @return size in current collection 
  template<AnaClass T> std::size_t size() const noexcept {
    return std::get<ObjCollection<T>>( fObjects ).size(); 
  }

  /// @tparam T Type of object 
  /// @return 'true' if current collection of that object type is empty, 'false' otherwise 
  template<AnaClass T> bool empty() const noexcept {
    return this->size<T>() < 1; 
  }

  /// @brief Add object of type 'T' to collection of objects for that type 
  /// @tparam T Type of object
  /// @param Args&&...args constructor arguments
  /// @return size of collection, after new object added
  template<AnaClass T, typename... Args> std::size_t emplace(Args&&... args) {
    auto& collection = std::get<ObjCollection<T>>( fObjects ); 
    return collection.emplace(std::forward<Args>(args)...); 
    return this->size<T>(); 
  }

  /// @brief Get collection of all objects of given type
  /// @tparam T type of object to fetch 
  /// @return ObjCollection<T> of all objects currently in event 
  template<AnaClass T> const ObjCollection<T>& GetCollection() const {
    return std::get<ObjCollection<T>>(fObjects); 
  }

  /// @brief Get collection of all objects of given type
  /// @tparam T type of object to fetch 
  /// @return ObjCollection<T> of all objects currently in event 
  template<AnaClass T> ObjCollection<T>& GetCollection() {
    return std::get<ObjCollection<T>>(fObjects); 
  }

  /// @brief Try to remove object with same ID as 'item'
  /// @tparam T Type of object
  /// @param item intsance 
  /// @return 'true' if an item is removed, 'false' otherwise 
  template<AnaClass T> bool remove_object(T& item) {
    auto& collection = std::get<ObjCollection<T>>( fObjects ); 
    return collection.remove_item(item); 
  }

  /// @brief Item collection with certain ID 
  /// @return 
  template<AnaClass T> T& get(std::size_t id) {
    auto& collection = std::get<ObjCollection<T>>(fObjects); 
    auto find_it = collection.get(id); 
    if (find_it == collection.end()) {
      throw std::out_of_range("in <get>: tried to access object using ref that doesn't exist"); 
    }
    return *find_it; 
  }

  /// @brief Compute the drift time for a given drift-distance, slope, (and derivative order).
  /// @param x absolute drift-distance from midplane of VDC drift chamber
  /// @param slope slope of current VDC track (affects some drift-parameters)
  /// @param derivative which derivative to compute, 0 (no derivative 1-st or 2-nd order. 
  /// @return drift time from track intercept with wire's 'vertical plane' 
  double Drift_T( double x,   double slope, int derivative=0 ) const;  

  //get the timing of this S2 event. 
  double S2Time() const; 

  void SetID(std::size_t _id) { id=_id; }
  /// @return unique ID of this event. unique among all events _for this arm_ but _not_ between arms (i.e, both arms will have an event with 'id = 0')
  std::size_t GetID() const { return id; }

  EventHandler& GetEventHandler() { return *fEventHandler.get(); } 

  double Get_tauSigma() const { return fEventHandler->Get_tauSigma(); } 

  bool Is_RHRS() const { return fis_RHRS; }

  /// Do this to make sure that each track keeps track of this object's address in memory (in case RDataFrame moves it)
  void RefereshEventPointers() {
    auto& tracks = std::get<ObjCollection<VDC::Track>>( fObjects ); 
    for (auto& track : tracks) track.SetEvent(this); 
  }

};
/////////////////////////////////////////////////////////////////////////////////

}

#endif 


/*
  std::size_t n_groups(int plane)   const { return fHitGroups[plane].size(); }
  std::size_t n_clusters(int plane) const { return fHitClusters[plane].size(); }
  std::size_t n_chamber_pairs(bool is_loChamber) const { return is_loChamber ? fChamberPairs_lo.size() : fChamberPairs_hi.size(); }
  std::size_t n_tracks() const { return fTracks.size(); }

  template<typename... Args> void add_group(int plane, Args&&... args) {
    //check plane range 
    if (plane<0 || plane>=VDC::kN_planes)
      throw std::out_of_range("in <"+std::string{__func__}+"> VDC plane index '"+std::to_string(plane)+"' out-of-range"); 
  
    fHitGroups[plane].emplace(std::forward<Args>(args)...); 
  }

  template<typename... Args> void add_cluster(int plane, Args&&... args) {
    //check plane range 
    if (plane<0 || plane>=VDC::kN_planes)
      throw std::out_of_range("in <"+std::string{__func__}+"> VDC plane index '"+std::to_string(plane)+"' out-of-range"); 
  
    fHitClusters[plane].emplace(std::forward<Args>(args)...); 
  }

  template<typename... Args> void add_chamber_pair(bool is_loChamber, Args&&... args) {
    auto& my_map = is_loChamber ? fChamberPairs_lo : fChamberPairs_hi; 
    my_map.emplace(std::forward<Args>(args)...); 
  }

  template<typename... Args> void add_track(Args&&... args) {
    fTracks.emplace(std::forward<Args>(args)...); 
  }  */ 
