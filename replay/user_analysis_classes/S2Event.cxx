
#include <APEX/S2Event.h>
#include <APEX/units.h>
#include <APEX/run_parameters.h>
// stdlib
#include <stdexcept> 
#include <string> 

namespace APEX
{

//_________________________________________________________________________________________________________
S2Event::S2Event(bool is_RHRS, S2Hit hit, const EventHandler* evt_handler, std::size_t _id) 
    : fis_RHRS{is_RHRS}, fS2Hit{hit}, id{_id}
{
    if (!fS2Hit.IsCoinc()) {
        throw std::invalid_argument(
            "in <"+std::string{__func__}+">: S2 hit is not valid!"
        ); 
    }

    if (!evt_handler) {
        throw std::invalid_argument(
            "in <"+std::string{__func__}+">: event handler passed is null"
        ); 
    }

    fEventHandler = std::make_shared<EventHandler>(
        Is_RHRS() ? replay::kRHRS : replay::kLHRS, 
        evt_handler->GetBeamCurrent(),
        evt_handler->GetRunNumber(), 
        Is_RHRS() ? &fS2Hit : evt_handler->GetS2Hit(false), 
        Is_RHRS() ? evt_handler->GetS2Hit(true)  : &fS2Hit   
    ); 
}
//_________________________________________________________________________________________________________
double S2Event::S2Time() const
{
    if (fS2Hit.Is_twinHit()) {
        // if it's a twin hit, then take the average of the two hits 
        return 0.5*( fS2Hit.Time() + fS2Hit.GetTwin()->Time() );  
    } else {
        return fS2Hit.Time(); 
    }
}
//_________________________________________________________________________________________________________
double S2Event::Drift_T(double x, double slope, int derivative) const 
{
    return fEventHandler->Drift_T(x,slope,derivative); 
}
/*
//_________________________________________________________________________________________________________
VDC::HitGroup& S2Event::get(int plane, ObjRef<VDC::HitGroup> ref)
{
    if (plane<0 || plane>=VDC::kN_planes)
      throw std::out_of_range("in <"+std::string{__func__}+"> VDC plane index '"+std::to_string(plane)+"' out-of-range");

    auto& my_map = fHitGroups[plane]; 
    auto itor = my_map.find(ref); 

    if (itor == fHitGroups[plane].end()) {
        throw std::out_of_range(
            "in <S2Event::"+std::string{__func__}+">: Ref to hit group is invalid."
        ); 
    }
    return itor->second; 
}
//_________________________________________________________________________________________________________
VDC::HitCluster& S2Event::get(int plane, ObjRef<VDC::HitCluster> ref)
{
    if (plane<0 || plane>=VDC::kN_planes)
      throw std::out_of_range("in <"+std::string{__func__}+"> VDC plane index '"+std::to_string(plane)+"' out-of-range"); 
    
    auto& my_map = fHitClusters[plane]; 
    auto itor = my_map.find(ref); 

    if (itor == my_map.end()) {
        throw std::out_of_range(
            "in <S2Event::"+std::string{__func__}+">: Ref to hit cluster is invalid."
        ); 
    }
    return itor->second; 
}
//_________________________________________________________________________________________________________
VDC::ChamberPair& S2Event::get(bool is_loChamber, ObjRef<VDC::ChamberPair> ref)
{
    auto& my_map = is_loChamber ? fChamberPairs_lo : fChamberPairs_hi;  
    auto itor = my_map.find(ref); 

    if (itor == my_map.end()) {
        throw std::out_of_range(
            "in <S2Event::"+std::string{__func__}+">: Ref to hit cluster is invalid."
        ); 
    }
    return itor->second; 
}
//_________________________________________________________________________________________________________
VDC::Track& S2Event::get(ObjRef<VDC::Track> ref)
{
    auto& my_map = fTracks; 
    auto itor = my_map.find(ref); 

    if (itor == my_map.end()) {
        throw std::out_of_range(
            "in <S2Event::"+std::string{__func__}+">: Ref to hit cluster is invalid."
        ); 
    }
    return itor->second; 
}
//_________________________________________________________________________________________________________
//_________________________________________________________________________________________________________
*/ 

}
