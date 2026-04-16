#ifndef gen_rawtracks_H
#define gen_rawtracks_H

//APEX headers
#include "../run_parameters.h"
#include "grid_search.h"
#include "theta_phi_model.h"
#include <TapexEventHandler.h>
#include <ApexVDCTrack.h> 
#include <ApexVDCChamberPair.h> 
#include <ApexUtils.h> 
//ROOT headers
#include <ROOT/RVec.hxx>
//stdlib headers
#include <math.h> 
#include <algorithm> 
#include <vector> 

ROOT::RVec<ApexVDC::Track> gen_rawtracks( 
    TapexEventHandler& evt, 
    ROOT::RVec<ApexVDC::ChamberPair>& pairs_Lo, 
    ROOT::RVec<ApexVDC::ChamberPair>& pairs_Hi ) 
{             
#ifdef DEBUG_RAW_TRACK
    std::printf("<%s>: in body\n", __func__); 
#endif

    using APEX::square; 

    const bool is_RHRS = evt.ActiveArm(); 
    
    const int arm_int_index = is_RHRS ? 0 : 1;  
    const double CUT_ph_min = run_parameters::CUT_ph_min[arm_int_index];
    const double CUT_ph_max = run_parameters::CUT_ph_max[arm_int_index];

    const double CUT_th_min = run_parameters::CUT_th_min[arm_int_index];
    const double CUT_th_max = run_parameters::CUT_th_max[arm_int_index];
    
    ROOT::RVec<ApexVDC::Track> tracks; 

    //returns ptr to track with given id 
    //__________________________________________________________________________________________________________________
    auto get_track_ptr = [&tracks](int id) {
#ifdef DEBUG_RAW_TRACK
        printf("<gen_rawtracks::get_track_ptr> searching for track %i...", id); 
#endif
        auto it = std::find_if(tracks.begin(), tracks.end(), [id](const ApexVDC::Track& rhs){ return rhs.GetID()==id; });
        if (it == tracks.end()) {
#ifdef DEBUG_RAW_TRACK
            printf("not found, returning nullptr\n"); 
#endif    
            return (ApexVDC::Track*)nullptr;
        }
#ifdef DEBUG_RAW_TRACK
        printf("found\n"); 
#endif    
        return (ApexVDC::Track*)it; 
    };
    //__________________________________________________________________________________________________________________
    
    //deletes track with given id 
    //__________________________________________________________________________________________________________________
    auto delete_track = [&tracks, &get_track_ptr](int id) {
#ifdef DEBUG_RAW_TRACK
        printf("<gen_rawtracks::delete_track> list of track ids: ");
        for (const auto& trk : tracks) printf("%i ", trk.GetID()); 
        printf("\n"); 
        printf("<gen_rawtracks::delete_track> deleting track %i\n", id); 
#endif
        auto trk = get_track_ptr(id); 
        if (!trk) return; 
        trk->GetPair_Lo()->Remove_track(id);
        trk->GetPair_Hi()->Remove_track(id);
        tracks.erase(
            std::remove_if(tracks.begin(), tracks.end(), [id](const ApexVDC::Track& rhs){ return rhs.GetID()==id; }),
            tracks.end()
        );
#ifdef DEBUG_RAW_TRACK
        printf("<gen_rawtracks::delete_track> list of track ids: ");
        for (const auto& trk : tracks) printf("%i ", trk.GetID()); 
        printf("\n"); 
#endif
        return; 
    };
    //__________________________________________________________________________________________________________________
    
    tracks.reserve(pairs_Hi.size()*pairs_Lo.size());
    for (int pH=0; pH<pairs_Hi.size(); pH++) { 
        for (int pL=0; pL<pairs_Lo.size(); pL++) { 
            
            auto& pLo = pairs_Lo.at(pL); 
            auto& pHi = pairs_Hi.at(pH); 
            
            //add a new track to the list of tracks
            tracks.emplace_back( &evt, &pLo, &pHi, evt.GenUniqueTrackID(is_RHRS) ); 
            
            ApexVDC::Track& track = tracks.back(); 
#ifdef DEBUG_RAW_TRACK
            printf("<gen_rawtracks> created test track %i (pairs lo-%i / hi-%i)\n", track.GetID(), pL,pH); 
#endif 

            double err_Theta = track.Theta() - Theta_model( track ); 
            double err_Phi   = track.Phi()   - Phi_model( track ); 
            
            //if this track has a good angular match, keep it. 
            if ( err_Theta < CUT_th_min ||
                err_Theta  > CUT_th_max ||
                err_Phi    < CUT_ph_min ||
                err_Phi    > CUT_ph_max ) { delete_track(track.GetID()); continue; }
            
            //check to make sure each plane passes the min-eta cut
            bool pass_minEtaCut=true; 
            
            double net_eta(0.); 
            
            for (int p=0; p<4; p++) { 
            
                double eta = grid_search( 
                    evt, 
                    *track.GetGroup(p), 
                    track.Slope(p), 
                    track.Slope(p), 
                    track.Intercept(p) -run_parameters::kGridSpacing*2.,
                    track.Intercept(p) +run_parameters::kGridSpacing*2.,  9e-9, 25e-9 
                ); 
                
                track.Set_Eta( p, eta ); 
                
                if (eta < run_parameters::kCUT_minEta) { pass_minEtaCut=false; break; }

            }//for (int p=0; p<4; p++) 
            
            if (!pass_minEtaCut) { delete_track(track.GetID()); continue; }

#ifdef DEBUG_RAW_TRACK
            printf("<gen_rawtracks> kept track %i; passed basic cuts\n", track.GetID()); 
#endif      
        }//for (int pL=0; pL<pairs_Lo.size(); pL++) 
    }//for (int pH=0; pH<pairs_Hi.size(); pH++) 
    
    //if we didn't find any tracks, quit
    if (tracks.size()<1) return {}; 
        
    //check to make sure that tracks aren't 'sharing' clusters
    
    //cout << "Size before pruning = " << tracks.size() << endl; 
    
    auto delete_shared_tracks = [&tracks, &get_track_ptr, &delete_track](ApexVDC::ChamberPair *pair) { 
#ifdef DEBUG_RAW_TRACK 
        printf("<gen_rawtracks::delete_shared_tracks> "
            "in body. event tracks ids (%zi total): { ", 
            tracks.size()); 
        for (auto& track : tracks) { printf("%i ", track.GetID()); }
        printf("}\n");

        printf("<gen_rawtracks::delete_shared_tracks> "
            "track ids in pair (%i total): { ", 
            pair->N_tracks()); 
        for (int t=0; t<pair->N_tracks(); t++) { printf("%i ", pair->GetTrackID(t)); }
        printf("}\n");
#endif 
        if (pair->N_tracks() <= 1) return; 
        
        //ApexVDC::Track* best_track = get_track_ptr(pair->GetTrackID(0));  
        int best_track_id = pair->GetTrackID(0); 
        double best_eta = get_track_ptr(best_track_id)->Get_Eta(); 

        //find the track with the highest eta
        for (int t=0; t<pair->N_tracks(); t++) { 
/*#ifdef DEBUG_RAW_TRACK 
            printf("<gen_rawtracks::delete_shared_tracks> "
                "tracks & id's (%zi total): ", 
                tracks.size()); 
            for (auto& track : tracks) { printf("%i ", track.GetID()); }
            printf("\n"); 
#endif*/  
            auto track = get_track_ptr(pair->GetTrackID(t));
            if (!track) {
                throw std::logic_error(Form(
                    "<gen_rawtracks::delete_shared_tracks>"
                    " track id %i (pair index %i) was not found in list of tracks",
                    pair->GetTrackID(t), t
                )); 
                return; 
            }
            if (track->Get_Eta() > best_eta) {
                best_track_id = track->GetID(); 
                best_eta = track->Get_Eta(); 
            }
        }
        
        //now, delete all tracks besides the best one
        int t=0; 
        while (pair->N_tracks() > 1) {
            
            int id = pair->GetTrackID(t); 
            if (id == best_track_id) { ++t; } else { delete_track(id); }; 
        }
    
    }; 
    
    for (int p=0; p<pairs_Lo.size(); p++) { 
        
        delete_shared_tracks( &pairs_Lo.at(p) ); 
    }
    
    for (int p=0; p<pairs_Hi.size(); p++) {  
    
        delete_shared_tracks( &pairs_Hi.at(p) ); 
    }
    //cout << "Size after pruning = " << tracks.size() << endl; 
    
    return tracks; 
}; 



#endif