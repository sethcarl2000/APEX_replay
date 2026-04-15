#ifndef gen_rawtracks_H
#define gen_rawtracks_H

//APEX headers
#include "run_parameters.h"
#include "functions/grid_search.h"
#include <TapexEventHandler.h>
#include <ApexVDCTrack.h> 
#include <ApexVDCChamberPair.h> 
#include <ApexUtils.h> 
//ROOT headers
#include <ROOT/RVec.hxx>
//stdlib headers
#include <math.h> 
#include <algorithm> 

ROOT::RVec<ApexVDC::Track> gen_rawtracks( 
    const TapexEventHandler& evt, 
    ROOT::RVec<ApexVDC::ChamberPair>& pairs_Lo, 
    ROOT::RVec<ApexVDC::ChamberPair>& pairs_Hi ) 
{             
    using APEX::square; 

    auto Phi_model   = [](const ApexVDC::Track& track) { 
        
        return 
        -0.231*square(track.S2_y()) 
        +0.2532*track.S2_y() 
        +0.539e-3; 
    };
    auto Theta_model = [](const ApexVDC::Track& track) { 
        
        return 0.109648*track.S2_x(); 
    };

    const bool is_RHRS = evt.ActiveArm(); 
    
    const double CUT_ph_min = is_RHRS ? -0.012 : -0.018;
    const double CUT_ph_max = is_RHRS ?  0.012 :  0.008;   

    const double CUT_th_min = is_RHRS ? -0.018 :-0.020;
    const double CUT_th_max = is_RHRS ?  0.018 : 0.018; 
        
    const double CUT_xParam = 2.00; 
    
    ROOT::RVec<ApexVDC::Track> tracks; 
    
    for (int pH=0; pH<pairs_Hi.size(); pH++) { 
        for (int pL=0; pL<pairs_Lo.size(); pL++) { 
            
            auto pLo = pairs_Lo.at(pL); 
            auto pHi = pairs_Hi.at(pH); 
            
            //track for this chamber pair
            ApexVDC::Track track( &evt, &pLo, &pHi ); 
            
            double err_Theta = track.Theta() - Theta_model( track ); 
            double err_Phi   = track.Phi()   - Phi_model( track ); 
            
            //if this track has a good angular match, keep it. 
            if ( err_Theta < CUT_th_min ||
                err_Theta  > CUT_th_max ||
                err_Phi    < CUT_ph_min ||
                err_Phi    > CUT_ph_max ) { continue; }
            
            //check to make sure each plane passes the min-eta cut
            bool pass_minEtaCut=true; 
            
            double net_eta(0.); 
            
            for (int p=0; (p<4 && pass_minEtaCut); p++) { 
            
                double eta = grid_search( 
                    evt, 
                    *track.GetGroup(p), 
                    track.Slope(p), 
                    track.Slope(p), 
                    track.Intercept(p) -run_parameters::kGridSpacing*2.,
                    track.Intercept(p) +run_parameters::kGridSpacing*2.,  9e-9, 25e-9 
                ); 
                
                track.Set_Eta( p, eta ); 
                
                if (eta < run_parameters::kCUT_minEta) pass_minEtaCut=false;
            }//for (int p=0; p<4; p++) 
            
            if (!pass_minEtaCut) { continue; }
            
            tracks.push_back( track ); 	

        }//for (int pL=0; pL<pairs_Lo.size(); pL++) 
    }//for (int pH=0; pH<pairs_Hi.size(); pH++) 
    
    //if we didn't find any tracks, quit
    if (tracks.size()<1) return {}; 
        
    //check to make sure that tracks aren't 'sharing' clusters
    
    //cout << "Size before pruning = " << tracks.size() << endl; 
    
    auto delete_shared_tracks = [&tracks](ApexVDC::ChamberPair *pair) { 
        
        if (pair->N_tracks() <= 1) return; 
                
        ApexVDC::Track* best_track = pair->GetTrack(0);  
        
        //find the track with the highest eta
        for (int t=0; t<pair->N_tracks(); t++) { 
            
            ApexVDC::Track* new_track = pair->GetTrack(t); 
            
            if ( new_track->Get_Eta() > best_track->Get_Eta() ) 
            best_track = new_track; 
        }
        
        //now, delete all other tracks
        for (int t=0; t<pair->N_tracks();) { 
            
            ApexVDC::Track* test_track = pair->GetTrack(t); 
            
            if ( best_track == test_track ) { t++; continue; } 
            
            //remove this track from the overall track vector
            std::remove_if( 
                tracks.begin(), tracks.end(), 
                [test_track](const ApexVDC::Track& rhs){ return &rhs == test_track; } 
            ); 
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