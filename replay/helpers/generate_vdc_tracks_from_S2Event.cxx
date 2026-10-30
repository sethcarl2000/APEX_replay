
// APEX headers
#include <APEX/replay.h> 
#include <APEX/run_parameters.h>
#include <APEX/EventHandler.h> 
#include <APEX/S2Event.h>
#include <APEX/S2Hit.h>
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
#include <numeric> 

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
void generate_vdc_tracks_from_S2hit( 
    EArmMode arm_mode, 
    RDFNodeAccumulator& rna, 
    EventCounter_RPtr &nPass_1group, 
    EventCounter_RPtr &nPass_1pair, 
    EventCounter_RPtr &nPass_1rawTrack,
    EventCounter_RPtr &nPass_1refinedTrack
) 
{
    const bool is_RHRS = (arm_mode == kRHRS); 

    using RVecD = ROOT::RVec<double>; 
    using namespace ROOT::VecOps; 

    //feed this script your node, with generated S2-hits, ready for tracking data. 

    //it'll spit out refined tracks for the plane you tell it to. 

    std::string arm = is_RHRS ? "R" : "L" ; 

    std::vector<std::string> plane_name = { "u1", "v1", "u2", "v2" }; 

    std::vector<std::string> branch_rawtime; 
    std::vector<std::string> branch_wire; 

    for (int p=0; p<4; p++) { 
        
        std::string rawtime = arm+".vdc."+plane_name[p]+".rawtime"; 
        
        std::string wire    = arm+".vdc."+plane_name[p]+".wire"; 
        
        branch_rawtime.push_back( rawtime ); 
        branch_wire   .push_back( wire );   
    }
        
    const char branch_event_vec[] = "coinc_events"; 

    //name of the branch which points to the (single) event manager we want to use for this event
    std::string branch_event = arm + "_event_handler"; 

    //pick the first event manager

    //first, make sure there is at least 1 event for this event (this cut should have already been made, but its good to be safe).
    rna.Filter(RVec_not_empty<EventHandler>, {branch_event_vec}); 


    

    std::string br_S2Event_array = arm+"_S2Events"; 
    
    //check to make sure this branch doesn't already exist
    if (rna.IsBranchDefined(br_S2Event_array)) {
        throw std::runtime_error(
            "in <"+std::string{__func__}+">:\n"
            "   it looks like this function has been called more than once for the " + (is_RHRS? "RHRS":"LHRS") + ", which is not allowed.\n"
            "   The branch '" + br_S2Event_array + "' has already beeen defined, which should not be possible if it was only called once-per-arm."
        );
        return; 
    }
    
    rna.Define(br_S2Event_array, [is_RHRS](const RVec<EventHandler>& events, const RVec<S2Hit>& S2_hits)
    {   
        // use the first event in the list as a 'template' for all following events
        // it's really arbitrary which event we use; what matters is that we have an event we can pass down
        // which keeps track of the run number and beam current information. 
        const auto& evt_template = events.front(); 

        // define an array of 'S2Events', which will be the basis of all our data. 
        RVec<S2Event> S2_events; S2_events.reserve(S2_hits.size()); 

        for (const auto& S2_hit : S2_hits) {

            EventHandler temp(
                is_RHRS ? kRHRS : kLHRS, 
                evt_template.GetBeamCurrent(), 
                evt_template.GetRunNumber(),
                is_RHRS ? &S2_hit : evt_template.GetS2Hit(false), 
                is_RHRS ? evt_template.GetS2Hit(true)  : &S2_hit
            ); 

            S2_events.emplace_back(
                is_RHRS,
                S2_hit, 
                &temp
            ); 
        }

        return S2_events; 

    }, {branch_event_vec, arm+"_S2_hits"}); 

    ///now, group VDC hits

    std::vector<std::string> branch_group; 

    for (int p=0; p<VDC::kN_planes; p++) { 

        //form the groups
        rna.Overwrite(br_S2Event_array, [p](RVec<S2Event>& S2_events, const RVecD& wire, const RVecD& rawtime)
            {   
                // go through the list of S2 events, make VDC hit groups for each one 
                for (auto& S2_event : S2_events) {

                    auto& group_vec = S2_event.get_groups(p); 

                    group_vec = group_vdc_hits(
                        S2_event.GetEventHandler(),
                        p, 
                        wire, 
                        rawtime
                    );
                }

                return S2_events; 
            },
            {br_S2Event_array, branch_wire[p], branch_rawtime[p]});
    } 

    //sort out events that don't have at least one group
    rna.Overwrite(br_S2Event_array, [](RVec<S2Event>& S2_events){

        // sort events by whether or not they have at least 1 group. 
        // the 'sort' aglorithm is analogous to the '<' comparator, and it uses 
        // the set-theory convention that two elements are equivalent in a poset when a < b && b < a. 
        // since we want to sort by the yes-no criteria that an event has a minimum number of groups, 
        // our criteria becomes a < b if 'a' has min. N. groups, and b < a if 'b' has min N. groups. 
        // this way, only those events _without_ the min N. groups survive.
        S2_events.erase(std::remove_if(S2_events.begin(), S2_events.end(), [](S2Event& e){
            for (int p=0; p<VDC::kN_planes; p++) { if (e.get_groups(p).empty()) return true; }
            return false;  
        }), S2_events.end()); 

        return S2_events; 
    }, {br_S2Event_array}); 

    // Filter out events that don't have any S2 events left 
    rna.Filter(RVec_not_empty<S2Event>, {br_S2Event_array}); 

    //record the number of events that form at least 1 group 
    nPass_1group = rna.Count(); 
     
    rna.Overwrite(br_S2Event_array, [](RVec<S2Event>& S2_events){

        static constexpr bool kLowChamber{true}, kHiChamber{false}; 

        // go through the list of S2 events, make VDC hit groups for each one 
        for (auto& S2_event : S2_events) {

            auto& pairs_lo = S2_event.get_chamber_pairs(kLowChamber); 

            pairs_lo = gen_pairs(
                S2_event.GetEventHandler(), 
                S2_event.get_groups(0), 
                S2_event.get_groups(1),
                kLowChamber
            ); 

            auto& pairs_hi = S2_event.get_chamber_pairs(kHiChamber); 

            pairs_hi = gen_pairs(
                S2_event.GetEventHandler(), 
                S2_event.get_groups(2), 
                S2_event.get_groups(3),
                kHiChamber
            ); 
        }

        //now, sort out events that don't have at least 1 pair. 
        S2_events.erase(std::remove_if(S2_events.begin(), S2_events.end(), [](S2Event& e){
            if (e.get_chamber_pairs(true) .empty()) { return true; }
            if (e.get_chamber_pairs(false).empty()) { return true; }
            return false; 
        }), S2_events.end()); 

        return S2_events; 

    }, {br_S2Event_array}); 

    // Filter out events that don't have any S2 events left 
    rna.Filter(RVec_not_empty<S2Event>, {br_S2Event_array}); 
    
    // count the number of CODA events which have at least 1 chamber-pair for both the hi and lo chambers
    nPass_1pair = rna.Count(); 

    // Now, generate 'Raw' tracks
    rna.Overwrite(br_S2Event_array, [](RVec<S2Event>& S2_events){

        static constexpr bool kLowChamber{true}, kHiChamber{false}; 

        // go through the list of S2 events, make VDC hit groups for each one 
        for (auto& S2_event : S2_events) {
            
            auto& tracks = S2_event.get_tracks(); 

            tracks = gen_rawtracks(
                S2_event.GetEventHandler(),
                S2_event.get_chamber_pairs(kLowChamber), 
                S2_event.get_chamber_pairs(kHiChamber)
            ); 
        }

        //now, sort out events that don't have at least 1 raw track. 
        S2_events.erase(std::remove_if(S2_events.begin(), S2_events.end(), [](S2Event& e){
            if (e.get_tracks().empty()) { return true; }
            return false; 
        }), S2_events.end()); 

        return S2_events; 

    }, {br_S2Event_array}); 
    
    
    // Filter out events that don't have any S2 events left 
    rna.Filter(RVec_not_empty<S2Event>, {br_S2Event_array}); 
    
    // count all CODA events which have at least 1 'raw' track left 
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
    rna.Overwrite(br_S2Event_array, [=](RVec<S2Event>& S2_events) { 
        
        for (auto& S2_event : S2_events) {

            const double tau_sigma = S2_event.GetEventHandler().Get_tauSigma(); 
            
            RVec<VDC::Track> refined_tracks; refined_tracks.reserve(S2_event.get_tracks().size()); 

            for (auto& trk : S2_event.get_tracks()) {    

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
                refined_tracks.emplace_back( trk );
            }

            //overwrite the list of 'raw' tracks with the new, 'refined' tracks
            S2_event.get_tracks() = refined_tracks; 
        }


        //now, sort out events that don't have at least 1 refined track. 
        S2_events.erase(std::remove_if(S2_events.begin(), S2_events.end(), [](S2Event& e){
            if (e.get_tracks().empty()) { return true; }
            return false; 
        }), S2_events.end()); 

        return S2_events; 

    }, {br_S2Event_array});

    std::string br_refined = arm+"_tracks_refined"; 
    rna.Define(br_refined, [](RVec<S2Event>& S2_events){
        
        std::size_t n_tracks{0}; 
        for (auto& evt : S2_events) n_tracks += evt.get_tracks().size();  

        RVec<VDC::Track> tracks; tracks.reserve(n_tracks);  

        std::size_t S2hit_id{0}; 
        for (auto& evt : S2_events) {

            auto& tracks_evt = evt.get_tracks(); 
            for (auto& track : tracks_evt) track.Set_S2_id(S2hit_id); 

            tracks.insert(tracks.end(), tracks_evt.begin(), tracks_evt.end());
            
            ++S2hit_id; 
        } 
        
        return tracks; 
        
    }, {br_S2Event_array}); 

    // Filter out events that don't have any S2 events left 
    rna.Filter(RVec_not_empty<VDC::Track>, {br_refined}); 

    // keep track of the ID of the S2 hit associated with this track 
    rna.Define(arm+"_S2hit_id", [](const RVec<VDC::Track>& trks){
        RVec<std::size_t> ids; ids.reserve(trks.size());
        for (const auto& t : trks) ids.emplace_back(t.Get_S2_id()); 
        return ids; 
    }, {br_refined}); 

    // count all CODA events which have at least 1 'refined' track left 
    nPass_1refinedTrack = rna.Count(); 
    
    return; 
}

}
}
}
