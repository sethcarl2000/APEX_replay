#ifndef generate_vdc_tracks_H
#define generate_vdc_tracks_H

// APEX headers
#include "group_vdc_hits.h"
#include "grid_search.h"
#include "gen_rawtracks.h"
#include "refine_track.h"
#include "gen_pairs.h"
#include "compute_trackdata.h"
#include "compute_track_error.h"
#include "theta_phi_model.h"
#include "../run_parameters.h"
#include <EventCounter.h> 
#include <TapexEventHandler.h> 
#include <ApexVDCHitGroup.h> 
// ROOT headers
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
// stdlib headers
#include <vector>
#include <string> 
#include <algorithm> 

namespace { 
    /// @return 'true' if RVec iss non-empty  
    template<typename T> bool RVec_not_empty(const ROOT::RVec<T>& v) { return !v.empty(); }

    /// @brief delete an element from a collection  
    template<typename T> void remove_element(T& elem, ROOT::RVec<T>& collection) {
        std::remove_if( collection.begin(), collection.end(), [&elem](const T& rhs) { return (&elem)==(&rhs); } );
    }
}

/// @brief Generates VDC tracks for RHRS / LHRS
/// @param is_RHRS
/// @param node_in input RDF node
/// @param n_pass_1group EventCounter representing the number of events which reconstructed at least 1 group 
ROOT::RDF::RNode generate_vdc_tracks( 
    const bool is_RHRS, 
    ROOT::RDF::RNode inNode, 
    EventCounter &nPass_1group, 
    EventCounter &nPass_1pair, 
    EventCounter &nPass_1rawTrack,
    EventCounter &nPass_1refinedTrack
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
        
    string branch_event = (string)TString("event_"+arm); 
    
    vector<string> branch_group; 
    
    //Now that we have eliminated non-coinc events, we can proceed with the 
    // analysis. (Starting with the right arm)
    inNode = inNode
        .Define(branch_event, [is_RHRS](TapexEventHandler *event) 
            { event->SetActiveArm(is_RHRS); return event; }, {"event"}); 
    
    for (int p=0; p<4; p++) { 
        
        branch_group.push_back( (string)"groups_"+arm+"_"+plane_name[p] ); 
        
        inNode = inNode
            .Define(branch_group[p].data(), [p](const TapexEventHandler& evt, const RVecD& wire, const RVecD& time)
            { 
                return group_vdc_hits(evt,p,wire,time); 
            },
            {branch_event, branch_wire[p].data(), branch_rawtime[p].data()});
    } 

    //genrate groups
    auto nEvents_1group = inNode  
        .Filter( RVec_not_empty<ApexVDC::HitGroup>, {branch_group[0].data()} )
        .Filter( RVec_not_empty<ApexVDC::HitGroup>, {branch_group[1].data()} )
        .Filter( RVec_not_empty<ApexVDC::HitGroup>, {branch_group[2].data()} )
        .Filter( RVec_not_empty<ApexVDC::HitGroup>, {branch_group[3].data()} ); 
    
    
    //generate pairs
    auto nEvents_1pair = nEvents_1group
        
        //create lo-chamber pairs
        .Define("pairs_"+arm+"_LoChamber", [](const TapexEventHandler& evt, 
            RVec<ApexVDC::HitGroup>& vec_gU, 
            RVec<ApexVDC::HitGroup>& vec_gV) 
            { return gen_pairs(evt, vec_gU,vec_gV, true); }, 
            { branch_event, branch_group[0].data(), branch_group[1].data() })
        
        .Filter( RVec_not_empty<ApexVDC::ChamberPair>, {"pairs_"+arm+"_LoChamber"}) 
        
        //create hi-chamber pairs
        .Define("pairs_"+arm+"_HiChamber", [](const TapexEventHandler& evt, 
            RVec<ApexVDC::HitGroup>& vec_gU, 
            RVec<ApexVDC::HitGroup>& vec_gV) 
            { return gen_pairs(evt, vec_gU,vec_gV, false); }, 
            {branch_event, 
            branch_group[2].data(), 
            branch_group[3].data() })
        
        .Filter( RVec_not_empty<ApexVDC::ChamberPair>, {"pairs_"+arm+"_HiChamber"}); 
        

    auto nEvents_1rawTrack = nEvents_1pair 
        
        //generate raw tracks
        .Define("tracks_"+arm+"_raw", gen_rawtracks, 
            {branch_event, "pairs_"+arm+"_LoChamber", "pairs_"+arm+"_HiChamber"}) 
        
        .Filter(RVec_not_empty<ApexVDC::Track>, {"tracks_"+arm+"_raw"}); 
    
    
    const double CUT_goodPoints_error = 40e-9; 
    
    int arm_int_index = is_RHRS ? 0 : 1; 
    const double CUT_ph_min = run_parameters::CUT_ph_min[arm_int_index];
    const double CUT_ph_max = run_parameters::CUT_ph_max[arm_int_index];

    const double CUT_th_min = run_parameters::CUT_th_min[arm_int_index];
    const double CUT_th_max = run_parameters::CUT_th_max[arm_int_index];
        
    auto nEvents_1refinedTrack = nEvents_1rawTrack
        
        //refine track candidates
        .Define("tracks_"+arm+"_refined", [
                        CUT_th_min, CUT_th_max, 
                        CUT_ph_min, CUT_ph_max]
            ( ROOT::RVec<ApexVDC::Track>& tracks ) { 
            
            ROOT::RVec<ApexVDC::Track> refined_tracks; 
                
            const double tau_sigma = tracks.at(0).GetEvent()->Get_tauSigma(); 
            
            for (int t=0; t<tracks.size(); t++) {
                
                auto& trk = tracks.at(t); 
                
                refine_track(trk, 20, 25e-9); 
            
        #if DEBUG
                bool is_nan=false; 
                
                for (int p=0; p<4; p++) 
                if ( trk->Intercept(0)!=trk->Intercept(0) ) is_nan=true; 
                
                if ( trk->T0()!=trk->T0() ) is_nan=true; 
                
                if (is_nan) { cout << "NAN intercept!!!!!!!!!!!" << endl; }
                else        { cout << "intercept exists." << endl; }
                
                cout << " intercepts = { " ; 
                for (int p=0; p<4; p++) { 
                cout << TString::Format("%0.3f ", trk->Intercept(p)); 
                } cout << "}" << endl; 
                cout << TString::Format( "Theta=%0.3f, Phi=%0.3f", 
                            trk->Theta(), trk->Phi() ) << endl; 
                
        #endif 
                
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
                
                //do some basic checks
                double xParam = trk.xParam();
                double Dt     = trk.T0(); 
                        
                if (trk.Get_Eta() < run_parameters::TRK_CUT_Eta || 
                    trk.Get_nGoodPoints(0) < run_parameters::TRK_CUT_nGoodPts_min_perPlane || 
                    trk.Get_nGoodPoints(1) < run_parameters::TRK_CUT_nGoodPts_min_perPlane || 
                    trk.Get_nGoodPoints(2) < run_parameters::TRK_CUT_nGoodPts_min_perPlane || 
                    trk.Get_nGoodPoints(3) < run_parameters::TRK_CUT_nGoodPts_min_perPlane || 
                    trk.Get_nGoodPoints()  < run_parameters::TRK_CUT_nGoodPts_min ||
                    std::fabs(xParam) > run_parameters::TRK_CUT_xParam ||
                    std::fabs(Dt)     > run_parameters::TRK_CUT_Dt        ) { 
                
                    //delete this track 		  
                    continue; 
                } 
                
                //compute track error
                Compute_trackError( trk ); 

                //add this track to the list of refined tracks
                refined_tracks.push_back( trk );
            }
            return tracks; 
        }, {"tracks_"+arm+"_raw"})
        
        .Filter(RVec_not_empty<ApexVDC::Track>, {"tracks_"+arm+"_refined"}); 
    
    
    nPass_1group        = nEvents_1group   .Count(); 
    nPass_1pair         = nEvents_1pair    .Count(); 
    nPass_1rawTrack     = nEvents_1rawTrack.Count(); 
    nPass_1refinedTrack = nEvents_1refinedTrack.Count(); 
    
    return nEvents_1refinedTrack; 
}



#endif