#ifndef APEX_replay_helpers_h
#define APEX_replay_helpers_h

// APEX headers
#include <APEX/replay.h> 
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
#include <atomic> 


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
    double TAU_sigma  =9e-9, 
    double TAU_buffer =20e-9 
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

/// @brief Generates VDC tracks from an 'S2Event' object
/// @param is_RHRS
/// @param node_in input RDF node
/// @param n_pass_1group EventCounter representing the number of events which reconstructed at least 1 group 
void generate_vdc_tracks_from_S2Event( 
    EArmMode arm_mode, 
    RDFNodeAccumulator& rna, 
    EventCounter_atomic &nPass_1group, 
    EventCounter_atomic &nPass_1pair, 
    EventCounter_atomic &nPass_1rawTrack,
    EventCounter_atomic &nPass_1refinedTrack
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
