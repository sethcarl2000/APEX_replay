#ifndef APEX_replay_helpers_h
#define APEX_replay_helpers_h

#include <APEX/replay.h> 
#include <RDFNodeAccumulator.h>
// ROOT headers
#include <ROOT/RVec.hxx>

class TapexEventHandler; 
class TapexReactVertex; 
class TapexS2Hit; 
class ApexPidManager; 
namespace ApexVDC { 
    class Track;
    class ChamberPair;
    class HitGroup; 
}; 


namespace APEX
{
namespace replay
{
namespace helpers
{


double grid_search( 
    const TapexEventHandler& evt, 
    ApexVDC::HitGroup& group, 
    double m1, //slope 
    double m2, 
    double x_Lo,
    double x_Hi, 
    double TAU_sigma  =9e-9, 
    double TAU_buffer =20e-9 );

void Compute_trackError(ApexVDC::Track& trk);

void compute_trackdata(ApexVDC::Track &trk); 

ROOT::RVec<TapexEventHandler> gen_coinc_events(
    double dt_center,  //the 'central' TR - TL for this event
    double dt_cut,     //the cut for events too far away from this value
    double beam_current, 
    unsigned int run_number,
    const ROOT::RVec<TapexS2Hit>& R_s2_hits,
    const ROOT::RVec<TapexS2Hit>& L_s2_hits
);

ROOT::RVec<TapexEventHandler> gen_coinc_window_events(
    double dt_min, 
    double dt_max, 
    double beam_current, 
    unsigned int run_number,
    const ROOT::RVec<TapexS2Hit>& R_s2_hits,
    const ROOT::RVec<TapexS2Hit>& L_s2_hits
);

ROOT::RVec<ApexVDC::ChamberPair> gen_pairs( 
    const TapexEventHandler& evt, 
    ROOT::RVec<ApexVDC::HitGroup>& gVec_U, 
    ROOT::RVec<ApexVDC::HitGroup>& gVec_V, 
    bool is_LoChamber 
);


void gen_pid_data(const bool is_RHRS, RDFNodeAccumulator& rna, const ApexPidManager* pid_manager);

ROOT::RVec<ApexVDC::Track> gen_rawtracks( 
    TapexEventHandler& evt, 
    ROOT::RVec<ApexVDC::ChamberPair>& pairs_Lo, 
    ROOT::RVec<ApexVDC::ChamberPair>& pairs_Hi 
); 

/// @brief Generates react vertex
/// @param rna RDFNodeAccumulator to add branches to 
/// @param vtx const ptr to TapexReactVertex obj (arm specific!)
void gen_react_vertex(RDFNodeAccumulator& rna, TapexReactVertex* vtx);

ROOT::VecOps::RVec<TapexS2Hit> generate_S2_hits(
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

/// @brief Form 'ApexVDC::HitGroup's from groups of VDC hits
/// @param evt Apex event handler
/// @param p the VDC plane in question, [0,..,3].
/// @param h_rawtime a collection of all VDC wire rawtimes 
/// @param h_wire a collection of all VDC wire numbers
/// @return a vector of all valid VDC hit groups 
ROOT::RVec<ApexVDC::HitGroup> group_vdc_hits ( 
    const TapexEventHandler& evt, 
    int   p, 
    const ROOT::RVec<double>& h_wire,
    const ROOT::RVec<double>& h_rawtime 
);

//
// Refine tracks from the hi-chamber using newton's method
//   
void refine_track( 
    ApexVDC::Track& trk, 
    const int nCycles=10, 
    double sigma=25e-9 
); 

double Phi_model(const ApexVDC::Track& track);
double Theta_model(const ApexVDC::Track& track);

}
}
}


#endif