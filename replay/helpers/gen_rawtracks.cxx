
//APEX headers
#include <APEX/replay/helpers.h>
#include <APEX/run_parameters.h>
#include <APEX/EventHandler.h>
#include <APEX/VDC/Track.h> 
#include <APEX/VDC/ChamberPair.h> 
#include <APEX/utils.h> 
//ROOT headers
#include <ROOT/RVec.hxx>
//stdlib headers
#include <math.h> 
#include <algorithm> 
#include <vector> 

namespace APEX
{
namespace replay
{
namespace helpers
{

void gen_rawtracks( 
    S2Event& evt
) 
{             
#ifdef DEBUG_RAW_TRACK
    std::printf("<%s>: in body\n", __func__); 
#endif

    // VDC chambers 
    constexpr VDC::Chamber Hi = VDC::Chamber::High; 
    constexpr VDC::Chamber Lo = VDC::Chamber::Low; 

    const bool is_RHRS = evt.Is_RHRS(); 
    
    const int arm_int_index = is_RHRS ? 0 : 1;  
    const double CUT_ph_min = run_parameters::CUT_ph_min[arm_int_index];
    const double CUT_ph_max = run_parameters::CUT_ph_max[arm_int_index];

    const double CUT_th_min = run_parameters::CUT_th_min[arm_int_index];
    const double CUT_th_max = run_parameters::CUT_th_max[arm_int_index];

    using Pair_lo = VDC::ChamberPair<Lo>; 
    using Pair_hi = VDC::ChamberPair<Hi>; 

    const auto& pairs_lo = evt.GetCollection<Pair_lo>(); 

    evt.GetCollection<VDC::Track>().reserve( 
        evt.size<Pair_lo>() * evt.size<Pair_hi>()
    );

    for (const auto& pair_lo : evt.GetCollection<Pair_lo>()) { 
        for (const auto& pair_hi : evt.GetCollection<Pair_hi>()) { 
            
            //add a new track to the list of tracks
            //tracks.emplace_back( &evt, &pLo, &pHi, evt.GenUniqueTrackID(is_RHRS) ); 
            
            VDC::Track track(&evt, pair_lo, pair_hi);

            double err_Theta = track.Theta() - Theta_model( track ); 
            double err_Phi   = track.Phi()   - Phi_model( track ); 
            
            //if this track has a good angular match, keep it. 
            if ( err_Theta  < CUT_th_min ||
                 err_Theta  > CUT_th_max ||
                 err_Phi    < CUT_ph_min ||
                 err_Phi    > CUT_ph_max ) { continue; }
            
            //check to make sure each plane passes the min-eta cut
            bool pass_minEtaCut=true; 
            
            VDC::for_each_plane([&evt, &track]<VDC::Plane P>(){ 
            
                double eta = grid_search<P>( 
                    evt, 
                    track.GetGroup<P>(), 
                    track.Slope<P>(), 
                    track.Slope<P>(), 
                    track.GetIntercept<P>() - run_parameters::kGridSpacing*2.,
                    track.GetIntercept<P>() + run_parameters::kGridSpacing*2.,  
                    9.*units::ns, 
                    25.*units::ns 
                ); 
                
                track.SetEta<P>(eta);
            });

            if (!pass_minEtaCut) { continue; }

            //add this track to the list 
            evt.emplace<VDC::Track>( track ); 

        } // for (auto& pair_hi : evt.GetCollection<Pair_hi>())
    } // for (auto& pair_lo : evt.GetCollection<Pair_lo>())
    
    //if we didn't find any tracks, quit
    if (evt.empty<VDC::Track>()) return; 
        
    /// @brief Scan all tracks that share a given cluster, only keep the track with the highest 'eta' 
    VDC::for_each_chamber([&evt]<VDC::Chamber C>(){ 

        for (auto& pair : evt.GetCollection<VDC::ChamberPair<C>>()) {

            //check to seem how many tracks match this cluster
            auto& all_tracks = evt.GetCollection<VDC::Track>(); 

            std::vector<VDC::Track*> pair_tracks; pair_tracks.reserve(all_tracks.size()); 

            std::size_t n_tracks_that_match{0}; 
            for (auto& track : evt.GetCollection<VDC::Track>()) {

                const auto& other_pair = track.GetPair<C>(); 
            
                if (other_pair == pair) pair_tracks.emplace_back(&track); 
            }

            if (pair_tracks.empty()) { 
                // no tracks to compare here 
                continue; 
            } 
            
            VDC::Track* best_track = &(*all_tracks.begin()); 
            double best_eta = best_track->GetEta(); 

            //find the track with the highest eta
            for (auto track : pair_tracks) { 

                if (track->GetEta() > best_eta) {
                    best_track = track; 
                    best_eta = track->GetEta(); 
                }
            }
            
            //now, delete all tracks besides the best one
            for (auto track : pair_tracks) {
                if (track != best_track) {
                    evt.remove_object(*track); 
                }
            }
        } // for (auto& pair : evt.GetCollection<VDC::ChamberPair<C>>())
    }); 

    //cout << "Size after pruning = " << tracks.size() << endl; 
    return; 
}; 

}
}
}