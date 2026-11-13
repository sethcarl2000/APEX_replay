
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
        std::printf("<remove and sort>: N. items for each event:<remove and sort>: ");
        for (const auto& S2_event : S2_events) std::printf("%2zu  ", S2_event.size<T>()); 
        std::printf("\n"); 
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
        std::printf("<remove and sort>: N. events end: %zu\n", S2_events.size());
#endif
    }

     

}

/// @brief Generates VDC tracks for RHRS / LHRS
/// @param is_RHRS
/// @param node_in input RDF node
/// @param n_pass_1group EventCounter representing the number of events which reconstructed at least 1 group 
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
    EventCounter_atomic &total_tracks,
    bool keep_empty_events
) 
{   
    /// increment an event counter by 1 
    auto increment = [](EventCounter_atomic& counter, EventCounter step=1) { 
        return counter.fetch_add(step, std::memory_order_relaxed); 
    };

    using RVecD = ROOT::RVec<double>; 
    using namespace ROOT::VecOps; 

    ///now, group VDC hits
    // define an array of 'S2Events', which will be the basis of all our data. 
    RVec<S2Event> S2_events; S2_events.reserve(S2_hits.size()); 
    
    if (event_handlers.empty()) return S2_events; 

    // use the first event in the list as a 'template' for all following events
    // it's really arbitrary which event we use; what matters is that we have an event we can pass down
    // which keeps track of the run number and beam current information. 
    const auto& evt_template = event_handlers.front(); 

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

    // count the total number of tracks
    increment( total_S2_hits, S2_hits.size() ); 

#ifdef DEBUG_S2Event
    std::printf("<Event Loop>: N. S2 events: %zu\n", S2_events.size());
#endif
    
#ifdef DEBUG_S2Event
        std::printf("<Event Loop>: New Event ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
#endif
    VDC::for_each_plane([&VDC_data,&S2_events,keep_empty_events]<VDC::Plane P>(){ 
        
        if (S2_events.empty()) return; 

#ifdef DEBUG_S2Event
        std::printf("<Event Loop>: Creating Groups for VDC Plane %i\n", VDC::number<P>());
#endif
        for (auto& S2_event : S2_events) {
            // do this to make sure each track has the correct event pointer. 
            S2_event.RefereshEventPointers(); 

            const auto& d = std::get<PlaneData<P>>( VDC_data );
            // form VDC::HitGroup objects for this plane
            group_vdc_hits<P>(S2_event, d.wire, d.rawtime); 
        }
#ifdef DEBUG_S2Event
        std::printf("<Event Loop>: remove & sort (%s)\n", VDC::plane_name<P>().c_str());
#endif
        // remove any S2 hits which don't have at least 1 group 
        if (!keep_empty_events) remove_and_sort<VDC::HitGroup<P>>(S2_events); 
    });

    // if there are no events left, then carry on 
    if (S2_events.empty()) return S2_events; 

#ifdef DEBUG_S2Event
    std::printf("<Event Loop>: %zu events have at least 1 group per plane\n", S2_events.size());
#endif

    // count the number of events we still have 
    increment( nPass_1group ); 

#ifdef DEBUG_S2Event
    std::printf("<Event Loop>: %llu events so far surviving '1 group' cut\n", static_cast<EventCounter>(nPass_1group));
#endif
    ///////////////////////////////////////////////////////////////////////////////////
    ///
    ///     Create VDC HitClusters &  VDC ChamberPairs
    ///
    VDC::for_each_chamber([&S2_events,keep_empty_events]<VDC::Chamber C>(){

        if (S2_events.empty()) return; 

#ifdef DEBUG_S2Event
    std::printf("<Event Loop>: making pairs for chamber 'VDC %s'\n", (C == VDC::Chamber::Low ? "Low" : "High"));
#endif

        // gen chamber pairs for each event
        for (auto& S2_event : S2_events) {

#ifdef DEBUG_S2Event
    std::printf("<Event Loop>: <event %zu> Refreshing event pointers\n", S2_event.GetID());
#endif
            // do this to make sure each track has the correct event pointer. 
            S2_event.RefereshEventPointers(); 

#ifdef DEBUG_S2Event
    std::printf("<Event Loop>: <event %zu> Making pairs for event\n", S2_event.GetID());
#endif
            chamber_grid_search<C>(S2_event, 
                run_parameters::min_cluster_gap, run_parameters::kCUT_minEta
            ); 
        }
        
#ifdef DEBUG_S2Event
        std::printf("<Event Loop>: Refresh & Sort (Chamber %s)\n", VDC::chamber_name<C>().c_str());
#endif
        // remove the events that didn't have any pairs generated
        if (!keep_empty_events) remove_and_sort<VDC::ChamberPair<C>>(S2_events); 
    });

    // if there are no events left, then carry on 
    if (S2_events.empty()) return S2_events; 
    // count the number of events we still have 
    increment( nPass_1pair ); 

#ifdef DEBUG_S2Event
    std::printf("<Event Loop>: %llu events so far surviving '1 pair' cut\n", static_cast<EventCounter>(nPass_1pair));
#endif

    ///////////////////////////////////////////////////////////////////////////////////
    ///
    ///     generate tracks for each S2 event
    ///
    for (auto& S2_event : S2_events) { 
        // do this to make sure each track has the correct event pointer. 
        S2_event.RefereshEventPointers(); 
#ifdef DEBUG_S2Event
    std::printf("<Event Loop>: Gen rawtracks\n");
#endif
        generate_all_tracks(S2_event); 
    }

    // at this point, we want to keep events that have 1 (and only 1) track compatible with the S2 hit. 
    auto new_end = std::remove_if(
        S2_events.begin(), 
        S2_events.end(), 
        [](const S2Event& event){ return event.size<VDC::Track>() != 1; }
    );

    S2_events.erase(
        new_end, 
        S2_events.end()
    );

#ifdef DEBUG_S2Event
    std::printf("<Event Loop>: Remove & sort rawtracks\n");
#endif

    // if there are no events left, then carry on 
    if (S2_events.empty()) return S2_events; 

    // count the number of events we still have 
    increment( nPass_1rawTrack ); 
    increment( nPass_1refinedTrack ); 

    // count the total number of tracks
    increment( total_tracks, S2_events.size() ); 

#ifdef DEBUG_S2Event
        std::printf("<Event Loop>: %llu events so far surviving '1 refined track' cut\n", static_cast<EventCounter>(nPass_1refinedTrack));
#endif
    return S2_events; 
}

}
}
}
