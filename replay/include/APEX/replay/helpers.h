#ifndef APEX_replay_helpers_h
#define APEX_replay_helpers_h

// APEX headers
#include <APEX/EventHandler.h>
#include <APEX/ReactVertex.h>
#include <APEX/S2Hit.h>
#include <APEX/PidManager.h>
// ---  VDC headers
#include <APEX/VDC/Track.h>
#include <APEX/VDC/ChamberPair.h>
#include <APEX/VDC/HitGroup.h>
// --- replay header
#include <APEX/replay.h> 
#include <RDFNodeAccumulator.h>
// ROOT headers
#include <ROOT/RVec.hxx>


namespace APEX
{
namespace replay
{
namespace helpers
{


double grid_search( 
    const EventHandler& evt, 
    VDC::HitGroup& group, 
    double m1, //slope 
    double m2, 
    double x_Lo,
    double x_Hi, 
    double TAU_sigma  =9e-9, 
    double TAU_buffer =20e-9 );

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

ROOT::RVec<VDC::ChamberPair> gen_pairs( 
    const EventHandler& evt, 
    ROOT::RVec<VDC::HitGroup>& gVec_U, 
    ROOT::RVec<VDC::HitGroup>& gVec_V, 
    bool is_LoChamber 
);


void gen_pid_data(const bool is_RHRS, RDFNodeAccumulator& rna, const PidManager* pid_manager);

ROOT::RVec<VDC::Track> gen_rawtracks( 
    EventHandler& evt, 
    ROOT::RVec<VDC::ChamberPair>& pairs_Lo, 
    ROOT::RVec<VDC::ChamberPair>& pairs_Hi 
); 

/// @brief Generates react vertex
/// @param rna RDFNodeAccumulator to add branches to 
/// @param vtx const ptr to ReactVertex obj (arm specific!)
void gen_react_vertex(RDFNodeAccumulator& rna, ReactVertex* vtx);

ROOT::VecOps::RVec<S2Hit> generate_S2_hits(
    const bool is_RHRS,
    const ROOT::RVec<double>& PMT_R, 
    const ROOT::RVec<double>& PMT_L 
);

/// @brief Generates VDC tracks for RHRS / LHRS
/// @param is_RHRS
/// @param node_in input RDF node
/// @param n_pass_1group EventCounter representing the number of events which reconstructed at least 1 group 
void generate_vdc_tracks( 
    const bool is_RHRS, 
    RDFNodeAccumulator& rna, 
    EventCounter_RPtr &nPass_1group, 
    EventCounter_RPtr &nPass_1pair, 
    EventCounter_RPtr &nPass_1rawTrack,
    EventCounter_RPtr &nPass_1refinedTrack
);

/// @brief Form 'VDC::HitGroup's from groups of VDC hits
/// @param evt Apex event handler
/// @param p the VDC plane in question, [0,..,3].
/// @param h_rawtime a collection of all VDC wire rawtimes 
/// @param h_wire a collection of all VDC wire numbers
/// @return a vector of all valid VDC hit groups 
ROOT::RVec<VDC::HitGroup> group_vdc_hits ( 
    const EventHandler& evt, 
    int   p, 
    const ROOT::RVec<double>& h_wire,
    const ROOT::RVec<double>& h_rawtime 
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