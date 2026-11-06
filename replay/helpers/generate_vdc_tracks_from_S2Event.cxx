
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
#include <tuple> 

namespace APEX
{
namespace replay
{
namespace helpers
{

//#define DEBUG_S2Event

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


    /// @brief Remove all S2 event which have an empty collection of type 'T' 
    /// @tparam T type of collection to measure 
    /// @param events list of S2 Events 
    template<typename T> void remove_and_sort(ROOT::VecOps::RVec<S2Event>& S2_events)
    {
#ifdef DEBUG_S2Event
        Info("remove_and_sort", "N. events start: %zu", S2_events.size());
#endif
        // this shifts all events to the left if they aren't empty, and to the right if they are.
        // returns an iterator pointing to the boundary between the haves and have-nots 
        auto new_end = std::remove_if(
            S2_events.begin(), 
            S2_events.end(), 
            [](const S2Event& evt){ return evt.empty<T>(); }
        ); 

        // now, delete all events which are empty 
        S2_events.erase(new_end, S2_events.end()); 

        // refresh the pointers of all the S2 hits 
        for (auto& event : S2_events) event.RefereshEventPointers(); 
#ifdef DEBUG_S2Event
        Info("remove_and_sort", "N. events end: %zu", S2_events.size());
#endif
    }

    template<VDC::Plane P> struct PlaneData { ROOT::VecOps::RVec<double> wire, rawtime; };
    
    using RawVDCData = std::tuple<
        PlaneData<VDC::Plane::U1>,
        PlaneData<VDC::Plane::V1>,
        PlaneData<VDC::Plane::U2>,
        PlaneData<VDC::Plane::V2>
        >; 

}

/// @brief Generates VDC tracks for RHRS / LHRS
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
) 
{   
    /// increment an event counter by 1 
    auto increment = [](EventCounter_atomic& counter) { return counter.fetch_add(1, std::memory_order_relaxed); };

    if (arm_mode == kBothArms) {
        throw std::invalid_argument(
            "in <"+std::string{__func__}+">: called with arm_mode=kBothArms, which is not valid (must be kRHRS or kLHRS)" 
        ); 
        return; 
    }
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


    ///now, group VDC hits

    const double CUT_goodPoints_error = 40e-9; 
    
    int arm_int_index = is_RHRS ? 0 : 1; 
    const double CUT_ph_min = run_parameters::CUT_ph_min[arm_int_index];
    const double CUT_ph_max = run_parameters::CUT_ph_max[arm_int_index];

    const double CUT_th_min = run_parameters::CUT_th_min[arm_int_index];
    const double CUT_th_max = run_parameters::CUT_th_max[arm_int_index];
        
    const double x_param_offset = is_RHRS ? run_parameters::TRK_R_xParam_offset : run_parameters::TRK_L_xParam_offset; 

    const double dt_offset = is_RHRS ? run_parameters::TRK_R_Dt_offset : run_parameters::TRK_L_Dt_offset; 
    const double dt_cut    = is_RHRS ? run_parameters::TRK_R_CUT_Dt : run_parameters::TRK_L_CUT_Dt; 


    //form the groups 
    rna.Overwrite(br_S2Event_array, [&](
        const RVec<EventHandler>& events, 
        const RVec<S2Hit>& S2_hits,
        const RVecD& wire_u1, const RVecD& rawtime_u1,
        const RVecD& wire_v1, const RVecD& rawtime_v1,
        const RVecD& wire_u2, const RVecD& rawtime_u2,
        const RVecD& wire_v2, const RVecD& rawtime_v2
    )
    {
        
        // use the first event in the list as a 'template' for all following events
        // it's really arbitrary which event we use; what matters is that we have an event we can pass down
        // which keeps track of the run number and beam current information. 
        const auto& evt_template = events.front(); 

        // define an array of 'S2Events', which will be the basis of all our data. 
        RVec<S2Event> S2_events; S2_events.reserve(S2_hits.size()); 

        std::size_t id{0}; 


        for (const auto& S2_hit : S2_hits) {

            EventHandler temp(
                is_RHRS ? kRHRS : kLHRS, 
                evt_template.GetBeamCurrent(), 
                evt_template.GetRunNumber(),
                is_RHRS ? &S2_hit : evt_template.GetS2Hit(false), 
                is_RHRS ? evt_template.GetS2Hit(true)  : &S2_hit
            ); 

            S2Event event(is_RHRS,S2_hit,&temp); 
                
            event.SetID(id); ++id;  

            S2_events.emplace_back( std::move(event) ); 
        }

#ifdef DEBUG_S2Event
        Info("Event Loop", "N. S2 events: %zu", S2_events.size());
#endif
        ///////////////////////////////////////////////////////////////////////////////////
        ///
        ///     Create VDC Hits &  VDC HitGroups
        ///
        RawVDCData data{
            {wire_u1, rawtime_u1},
            {wire_v1, rawtime_v1},
            {wire_u2, rawtime_u2},
            {wire_v2, rawtime_v2}
        };

        VDC::for_each_plane([&data,&S2_events]<VDC::Plane P>(){

#ifdef DEBUG_S2Event
        Info("Event Loop", "Creating Groups for VDC Plane %i", VDC::number<P>());
#endif
            for (auto& S2_event : S2_events) {
                // do this to make sure each track has the correct event pointer. 
                S2_event.RefereshEventPointers(); 

                const auto& d = std::get<PlaneData<P>>( data );
                // form VDC::HitGroup objects for this plane
                group_vdc_hits<P>(S2_event, d.wire, d.rawtime); 
            }
            // remove any S2 hits which don't have at least 1 group 
            remove_and_sort<VDC::HitGroup<P>>(S2_events); 
        });

        // if there are no events left, then carry on 
        if (S2_events.empty()) return S2_events; 

#ifdef DEBUG_S2Event
        Info("Event Loop", "%zu events have at least 1 group per plane", S2_events.size());
#endif

        // count the number of events we still have 
        increment( nPass_1group ); 

#ifdef DEBUG_S2Event
        Info("Event Loop", "%llu events so far surviving '1 group' cut", static_cast<EventCounter>(nPass_1group));
#endif
        ///////////////////////////////////////////////////////////////////////////////////
        ///
        ///     Create VDC HitClusters &  VDC ChamberPairs
        ///
        VDC::for_each_chamber([&S2_events]<VDC::Chamber C>(){

            //if (S2_events.empty()) return; 
#ifdef DEBUG_S2Event
        Info("Event Loop", "making pairs for chamber 'VDC %s'", (C == VDC::Chamber::Low ? "Low" : "High"));
#endif

            // gen chamber pairs for each event
            for (auto& S2_event : S2_events) {

#ifdef DEBUG_S2Event
        Info("Event Loop", "<event %zu> Refreshing event pointers", S2_event.GetID());
#endif
                // do this to make sure each track has the correct event pointer. 
                S2_event.RefereshEventPointers(); 

#ifdef DEBUG_S2Event
        Info("Event Loop", "<event %zu> Making pairs for event", S2_event.GetID());
#endif
                gen_pairs<C>(S2_event); 
            }
            
#ifdef DEBUG_S2Event
        Info("Event Loop", "Refresh & Sort (%zu total events)", S2_events.size());
#endif
            // remove the events that didn't have any pairs generated
            remove_and_sort<VDC::ChamberPair<C>>(S2_events); 
        });

        // if there are no events left, then carry on 
        if (S2_events.empty()) return S2_events; 
        // count the number of events we still have 
        increment( nPass_1pair ); 

#ifdef DEBUG_S2Event
        Info("Event Loop", "%llu events so far surviving '1 pair' cut", static_cast<EventCounter>(nPass_1pair));
#endif
        ///////////////////////////////////////////////////////////////////////////////////
        ///
        ///     generate 'raw' tracks for each S2 event
        ///
        for (auto& S2_event : S2_events) { 
            // do this to make sure each track has the correct event pointer. 
            S2_event.RefereshEventPointers(); 

            gen_rawtracks(S2_event); 
        }
        remove_and_sort<VDC::Track>(S2_events); 

        // if there are no events left, then carry on 
        if (S2_events.empty()) return S2_events; 
        // count the number of events we still have 
        increment( nPass_1rawTrack ); 

#ifdef DEBUG_S2Event
        Info("Event Loop", "%llu events so far surviving '1 raw track' cut", static_cast<EventCounter>(nPass_1rawTrack));
#endif

        ///////////////////////////////////////////////////////////////////////////////////
        ///
        ///     generate 'refined' tracks for each S2 event
        ///
        for (auto& S2_event : S2_events) {
            // do this to make sure each track has the correct event pointer. 
            S2_event.RefereshEventPointers(); 

            const double tau_sigma = S2_event.GetEventHandler().Get_tauSigma(); 
            
            RVec<VDC::Track> refined_tracks; refined_tracks.reserve(S2_event.size<VDC::Track>()); 

            auto& tracks = S2_event.GetCollection<VDC::Track>(); 

            //first, 'refine' the tracks. 
            for (auto& track : tracks) {
                refine_track(track, 20, 25.*units::ns); 
            }

            //now, check to make it has good angular agreement. 
            tracks.remove_if([=](const VDC::Track& trk){
                
                double err_Theta = trk.Theta() - Theta_model(trk); 
                double err_Phi   = trk.Phi()   - Phi_model(trk); 
                        
                if ( err_Theta < CUT_th_min || 
                     err_Theta > CUT_th_max || 
                     err_Phi   < CUT_ph_min || 
                     err_Phi   > CUT_ph_max ) { //cut this track 
                    
                    return true; //remove this track 
                } 
                return false; //keep this track 
            });

            // now, compute the 'track data' 
            for (auto& trk : tracks) {
                compute_trackdata(trk); 
            }

            // see if there are enough 'good points' for each plane 
            tracks.remove_if([=](const VDC::Track& trk){
                
                // see if the track's 'Eta' is good enough 
                if (trk.GetEta() < run_parameters::TRK_CUT_Eta) return true; 

                // see if each VDC Plane has the min. number of hits 
                if (trk.Get_nGoodPoints() < run_parameters::TRK_CUT_nGoodPts_min) true; 

                // see if there are enough hits overall 
                bool discard=false; 
                VDC::for_each_plane([&]<VDC::Plane P>(){
                    if (trk.Get_nGoodPoints<P>() < run_parameters::TRK_CUT_nGoodPts_min_perPlane) { discard=true; }
                });
                if (discard) return true; 

                return false; 
            });

            // see if this track's projection onto the S2 matches up well enough 
            tracks.remove_if([=](const VDC::Track& trk){    
                
                //do some basic checks
                double xParam = trk.xParam();
                double Dt     = trk.T0(); 

                if ( std::fabs(xParam-x_param_offset) > run_parameters::TRK_CUT_xParam || std::fabs(Dt-dt_offset) > dt_cut ) 
                { 
                    //discard this track
                    return true;  
                } 
                return false; 
            }); 

            //compute track errors
            for (auto& trk : tracks) {
                Compute_trackError( trk ); 
            }
        }
        //now, sort out events that don't have at least 1 refined track. 
        remove_and_sort<VDC::Track>(S2_events); 

        // if there are no events left, then carry on 
        if (S2_events.empty()) return S2_events; 
        // count the number of events we still have 
        increment( nPass_1refinedTrack );

#ifdef DEBUG_S2Event
        Info("Event Loop", "%llu events so far surviving '1 refined track' cut", static_cast<EventCounter>(nPass_1refinedTrack));
#endif
        return S2_events; 
    },
    {   branch_event_vec, 
        arm+"_S2_hits",
        arm+".vdc.u1.wire", arm+".vdc.u1.rawtime",
        arm+".vdc.v1.wire", arm+".vdc.v1.rawtime",
        arm+".vdc.u2.wire", arm+".vdc.u2.rawtime",
        arm+".vdc.v2.wire", arm+".vdc.v2.rawtime"
    });

    // Filter out events that don't have any S2 events left 
    rna.Filter(RVec_not_empty<S2Event>, {br_S2Event_array}); 
    
    return; 
}

}
}
}
