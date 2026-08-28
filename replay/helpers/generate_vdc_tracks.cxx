
// APEX headers
#include <APEX/run_parameters.h>
#include <APEX/EventCounter.h> 
#include <APEX/EventHandler.h> 
#include <APEX/VDC/HitGroup.h> 
#include <APEX/VDC/Track.h>
//
//functions 
#ifdef DEBUG_TRACK
//if the 'DEBUG_TRACK' flag is set, then propagate the debug macro down to sub-routines
#define DEBUG_GROUP
#define DEBUG_GRID
#define DEBUG_PAIR
#define DEBUG_RAW_TRACK
#endif
//
#include <APEX/replay/helpers.h>
//
// ROOT headers
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <TString.h> 
//
// stdlib headers
#include <vector>
#include <string> 
#include <algorithm> 

namespace APEX
{
namespace replay
{
namespace helpers
{

namespace { 
    /// @return 'true' if RVec iss non-empty  
    template<typename T> bool RVec_not_empty(const ROOT::RVec<T>& v) { return !v.empty(); }

    /// @brief delete an element from a collection  
    template<typename T> void remove_element(T& elem, ROOT::RVec<T>& collection) {
        collection.erase(
            std::remove_if( collection.begin(), collection.end(), [&elem](const T& rhs) { return (&elem)==(&rhs); } ),
            collection.end()
        ); 
    }
}

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
) 
{
    using namespace std;            

    using RVecD = ROOT::RVec<double>; 
    using namespace ROOT::VecOps; 

    //feed this script your node, with generated S2-hits, ready for tracking data. 

    //it'll spit out refined tracks for the plane you tell it to. 

    string arm = is_RHRS ? "R" : "L" ; 

    vector<string> plane_name = { "u1", "v1", "u2", "v2" }; 

    vector<string> branch_rawtime; 
    vector<string> branch_wire; 

    for (int p=0; p<4; p++) { 
        
        string rawtime = arm+".vdc."+plane_name[p]+".rawtime"; 
        
        string wire    = arm+".vdc."+plane_name[p]+".wire"; 
        
        branch_rawtime.push_back( rawtime ); 
        branch_wire   .push_back( wire );   
    }
        
    const char branch_event_vec[] = "coinc_events"; 

    //name of the branch which points to the (single) event manager we want to use for this event
    string branch_event = arm + "_event_handler"; 

    //pick the first event manager

    //first, make sure there is at least 1 event for this event (this cut should have already been made, but its good to be safe).
    rna.Filter(RVec_not_empty<EventHandler>, {branch_event_vec}); 

    rna.DefineIfMissing(branch_event, [is_RHRS](RVec<EventHandler>& events)
    {
        auto& event = events.front();  
        event.SetActiveArm(is_RHRS); 
        return event; 
    }, {branch_event_vec}); 
    
    
    vector<string> branch_group; 

    for (int p=0; p<4; p++) { 

        //create a name for this group
        branch_group.push_back( (string)arm+"_groups_"+plane_name[p] ); 
        
        //form the groups
        rna.Define(branch_group[p], [p](const EventHandler& evt, const RVecD& wire, const RVecD& time)
            { 
                return group_vdc_hits(evt,p,wire,time); 
            },
            {branch_event, branch_wire[p], branch_rawtime[p]});
        
        //require that at least 1 group was successfully formed 
        rna.Filter( RVec_not_empty<VDC::HitGroup>, {branch_group[p]} );
    } 

    //record the number of events that form at least 1 group 
    nPass_1group = rna.Count(); 
     
    //create lo-chamber pairs
    string br_pairs_lo = arm+"_pairs_LoChamber";
    rna.Define(br_pairs_lo, [](
        const EventHandler& evt, 
        RVec<VDC::HitGroup>& vec_gU, 
        RVec<VDC::HitGroup>& vec_gV) { 
            return gen_pairs(evt, vec_gU,vec_gV, true); 
        }, { branch_event, branch_group[0], branch_group[1] });
        
    rna.Filter( RVec_not_empty<VDC::ChamberPair>, {br_pairs_lo}); 
        
    //create hi-chamber pairs
    string br_pairs_hi = arm+"_pairs_HiChamber";
    rna.Define(arm+"_pairs_HiChamber", [](
        const EventHandler& evt, 
        RVec<VDC::HitGroup>& vec_gU, 
        RVec<VDC::HitGroup>& vec_gV) { 
            return gen_pairs(evt, vec_gU,vec_gV, false);
        }, {branch_event, branch_group[2], branch_group[3] });
        
    rna.Filter( RVec_not_empty<VDC::ChamberPair>, {br_pairs_hi}); 
    
    nPass_1pair = rna.Count(); 

    
    //generate raw tracks
    string br_tracks_raw = arm+"_tracks_raw";
    rna.Define(br_tracks_raw, 
        gen_rawtracks, 
        {branch_event, arm+"_pairs_LoChamber", arm+"_pairs_HiChamber"}); 
        
    rna.Filter(RVec_not_empty<VDC::Track>, {br_tracks_raw}); 
    
    nPass_1rawTrack = rna.Count(); 

    
    const double CUT_goodPoints_error = 40e-9; 
    
    int arm_int_index = is_RHRS ? 0 : 1; 
    const double CUT_ph_min = run_parameters::CUT_ph_min[arm_int_index];
    const double CUT_ph_max = run_parameters::CUT_ph_max[arm_int_index];

    const double CUT_th_min = run_parameters::CUT_th_min[arm_int_index];
    const double CUT_th_max = run_parameters::CUT_th_max[arm_int_index];
        
    const double x_param_offset = is_RHRS ? run_parameters::TRK_R_xParam_offset : run_parameters::TRK_L_xParam_offset; 

    const double dt_offset = is_RHRS ? run_parameters::TRK_R_Dt_offset : run_parameters::TRK_L_Dt_offset; 
    const double dt_cut    = is_RHRS ? run_parameters::TRK_R_CUT_Dt : run_parameters::TRK_L_CUT_Dt; 

    //refine track candidates
    string br_refined = arm+"_tracks_refined"; 
    rna.Define(br_refined, [CUT_th_min, CUT_th_max, CUT_ph_min, CUT_ph_max, x_param_offset, dt_offset,dt_cut] ( ROOT::RVec<VDC::Track>& tracks ) 
        { 
            RVec<VDC::Track> refined_tracks; refined_tracks.reserve(tracks.size());
                
            const double tau_sigma = tracks.at(0).GetEvent()->Get_tauSigma(); 
            
            for (auto& trk : tracks) {
                
                refine_track(trk, 20, 25e-9); 
                
                double err_Theta = trk.Theta() - Theta_model(trk); 
                double err_Phi   = trk.Phi()   - Phi_model(trk); 
                        
                if ( err_Theta < CUT_th_min || 
                    err_Theta > CUT_th_max || 
                    err_Phi   < CUT_ph_min || 
                    err_Phi   > CUT_ph_max ) { //cut this track 
                
                    continue; // remove_element( trk, tracks );
                } 
                
                //find out how many 'good' points each plane has for this track
                compute_trackdata( trk ); 
                
                if (trk.Get_Eta() < run_parameters::TRK_CUT_Eta || 
                    trk.Get_nGoodPoints(0) < run_parameters::TRK_CUT_nGoodPts_min_perPlane || 
                    trk.Get_nGoodPoints(1) < run_parameters::TRK_CUT_nGoodPts_min_perPlane || 
                    trk.Get_nGoodPoints(2) < run_parameters::TRK_CUT_nGoodPts_min_perPlane || 
                    trk.Get_nGoodPoints(3) < run_parameters::TRK_CUT_nGoodPts_min_perPlane || 
                    trk.Get_nGoodPoints()  < run_parameters::TRK_CUT_nGoodPts_min) 
                {
                    continue; 
                }
                
                //do some basic checks
                double xParam = trk.xParam();
                double Dt     = trk.T0(); 

                if ( std::fabs(xParam-x_param_offset) > run_parameters::TRK_CUT_xParam || std::fabs(Dt-dt_offset) > dt_cut ) 
                { 
                    //discard this track
                    continue;  
                } 

                //compute track error
                Compute_trackError( trk ); 
        
                //add this track to the list of refined tracks
                refined_tracks.push_back( trk );
            }
            return refined_tracks; 
        }, {arm+"_tracks_raw"});
        
    rna.Filter(RVec_not_empty<VDC::Track>, {br_refined}); 

    rna.DefineOutput(arm+"_tracks_S2time", [is_RHRS]( ROOT::RVec<VDC::Track>& tracks )
    {
      ROOT::RVec<double> time; time.reserve( tracks.size() ); 
      for (auto& trk : tracks) {
	time.push_back( trk.GetEvent()->GetS2Hit(is_RHRS)->Time() );
      }
      return time; 
    }, {arm+"_tracks_refined"}); 
    
    nPass_1refinedTrack = rna.Count(); 
    return; 
}

}
}
}
