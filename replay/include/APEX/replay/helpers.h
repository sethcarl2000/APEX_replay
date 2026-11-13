#ifndef APEX_replay_helpers_h
#define APEX_replay_helpers_h

// APEX headers
#include <APEX/units.h>
#include <APEX/replay.h> 
#include <APEX/ObjCollection.h>
#include <APEX/replay/EventCounter.h> 
#include <APEX/EventHandler.h>
#include <APEX/ReactVertex.h>
#include <APEX/S2Event.h> 
#include <APEX/S2Hit.h>
#include <APEX/PidManager.h>
// ---  VDC headers
#include <APEX/VDC.h> 
#include <APEX/VDC/Track.h>
#include <APEX/VDC/ChamberPair.h>
#include <APEX/VDC/HitGroup.h>
// --- replay header
#include <APEX/replay.h> 
#include <RDFNodeAccumulator.h>
// ROOT headers
#include <ROOT/RVec.hxx>
// stdlib
#include <vector> 
#include <cstddef> 
#include <cstddef> 


namespace APEX
{
namespace replay
{
namespace helpers
{


/// @brief Form 'VDC::HitGroup's from groups of VDC hits
/// @param evt Apex event handler
/// @param p the VDC plane in question, [0,..,3].
/// @param h_rawtime a collection of all VDC wire rawtimes 
/// @param h_wire a collection of all VDC wire numbers
/// @return a vector of all valid VDC hit groups 
template<VDC::Plane P> void group_vdc_hits( 
    S2Event& event,
    const ROOT::RVec<double>& h_wire,
    const ROOT::RVec<double>& h_rawtime 
);

template<VDC::Plane P> double grid_search( 
    S2Event& evt,  
    const VDC::HitGroup<P>& group, 
    double m1, //slope 
    double m2, 
    double x_Lo,
    double x_Hi, 
    double TAU_sigma  =9.*units::ns, 
    double TAU_buffer =20.*units::ns, 
    double TAU_offset =0.*units::ns
);

/// @brief The VDC chamber to iterate through. 
/// @tparam C chamber to process (VDC::Chamber::High or VDC::Chamber::Low)
/// @param evt S2 Event to base this around
template<VDC::Chamber C> void gen_pairs( 
    S2Event& evt
);


/// @brief generate Raw VDC tracks from chamber-pairs
/// @param evt S2 event to process 
void gen_rawtracks( 
    S2Event& evt
); 


void Compute_trackError(VDC::Track& trk);

void compute_trackdata(VDC::Track &trk); 

ROOT::RVec<EventHandler> gen_coinc_events(
    double dt_center,  //the 'central' TR - TL for this event
    double dt_cut,     //the cut for events too far away from this value
    double beam_current, 
    unsigned int run_number,
    const ROOT::RVec<S2Hit>& R_s2_hits,
    const ROOT::RVec<S2Hit>& L_s2_hits
);

ROOT::RVec<EventHandler> gen_coinc_window_events(
    double dt_min, 
    double dt_max, 
    double beam_current, 
    unsigned int run_number,
    const ROOT::RVec<S2Hit>& R_s2_hits,
    const ROOT::RVec<S2Hit>& L_s2_hits
);




void gen_pid_data(EArmMode is_RHRS, RDFNodeAccumulator& rna, const PidManager* pid_manager);


/// @brief Generates react vertex
/// @param rna RDFNodeAccumulator to add branches to 
/// @param vtx const ptr to ReactVertex obj (arm specific!)
void gen_react_vertex(RDFNodeAccumulator& rna, ReactVertex* vtx);

ROOT::VecOps::RVec<S2Hit> generate_S2_hits(
    EArmMode arm_mode,
    const ROOT::RVec<double>& PMT_R, 
    const ROOT::RVec<double>& PMT_L 
);

template<typename T> concept Comparable = std::is_default_constructible_v<T> || 
    requires(T& t1, T& t2) {  
        {t1 <= t2} -> std::same_as<bool>; 
    }; 

template<Comparable T> struct Bound { 
    T min{}, max{}; 

    T clamp(const T& val) const { return (val <= max) ? ((min <= val) ? val : min) : max; }

    bool within(const T& val) const { return (min <= val) && (val <= max); }
}; 

struct pair_candidate { 
    // intercept
    double u,v;

    // eta of each intercept
    double eta_u,eta_v; 

    double min_eta() const { return std::min(eta_u,eta_v); }

    bool operator<(const pair_candidate& rhs) const { return this->min_eta() < rhs.min_eta(); }

    // 'boundary' of this cluster. 
    Bound<double> boundary_u, boundary_v; 
    std::size_t id_U,id_V; 
};

/// @brief Search for cluster-pairs in each VDC chamber
/// @tparam C VDC Chamber
/// @param event S2 Event (containing 'VDC::HitGroup' collections)
/// @param min_cluster_gap Minimum gap between valid clusters, for them to not be combined 
/// @param min_eta Minimum eta (per-plane) for a choice of intercepts to be considered a cluster 
template<VDC::Chamber C> std::vector<pair_candidate> chamber_grid_search( 
    S2Event& event,
    double min_cluster_gap,
    double min_eta
); 

/// @brief generate 'raw' tracks, and check which ones are compatible with the given S2 hit. 
/// @param event S2 Event to generate tracks from
void generate_all_tracks(S2Event& event); 


template<VDC::Plane P> struct PlaneData { ROOT::VecOps::RVec<double> wire, rawtime; };
    
using RawVDCData = std::tuple<
    PlaneData<VDC::Plane::U1>,
    PlaneData<VDC::Plane::V1>,
    PlaneData<VDC::Plane::U2>,
    PlaneData<VDC::Plane::V2>
>;

/// @brief Generates VDC tracks from an 'S2Event' object
ROOT::RVec<S2Event> generate_vdc_tracks_from_S2Event( 
    const bool is_RHRS, 
    const ROOT::RVec<EventHandler>& event_handlers,
    const ROOT::RVec<S2Hit>& S2_hits, 
    const RawVDCData& VDC_data, 
    EventCounter_atomic &nPass_1group, 
    EventCounter_atomic &nPass_1pair, 
    EventCounter_atomic &nPass_1rawTrack,
    EventCounter_atomic &nPass_1refinedTrack, 
    EventCounter_atomic &total_S2_hits,
    EventCounter_atomic &total_refined_tracks,
    bool keep_empty_events=false
);

//
// Refine tracks from the hi-chamber using newton's method
//   
void refine_track( 
    VDC::Track& trk, 
    const int nCycles=10, 
    double sigma=25e-9 
); 

double Phi_model(const VDC::Track& track);
double Theta_model(const VDC::Track& track);

}
}
}


#endif
