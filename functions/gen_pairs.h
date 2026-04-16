#ifndef gen_pairs_H
#define gen_pairs_H

#include <ApexVDCHitGroup.h>
#include <ApexVDCHitCluster.h>
#include <ApexVDCChamberPair.h> 
#include <ApexVDCTrack.h> 
#include "grid_search.h"
#include "../run_parameters.h"
#include <math.h> 
#include <TVector3.h> 
#include <stdio.h> 

//#define DEBUG_PAIR

ROOT::RVec<ApexVDC::ChamberPair> gen_pairs( 
    const TapexEventHandler& evt, 
    ROOT::RVec<ApexVDC::HitGroup>& gVec_U, 
    ROOT::RVec<ApexVDC::HitGroup>& gVec_V, 
    bool is_LoChamber ) 
{ 
#ifdef DEBUG_PAIR
#endif
    
#ifdef DEBUG_PAIR
    std::printf("<%s>: in body\n", __func__); 
#endif
    //use by tracks to identify pairs later
    int pair_unique_id(0); 
    
    //constants
    bool is_RHRS = evt.ActiveArm(); 
    
    //this involves angular prediction based on typical S2-angles. 
    const double uFix = 0.026/std::sqrt(2);
    
    const double wVDC = gVec_U.at(0).W(); 

    const double M_Theta = is_RHRS ? 0.1096 : 0.1096 ; 
    const double M_Phi   = is_RHRS ? 0.242  : 0.242  ; 

#ifdef DEBUG_PAIR
    std::printf("<%s>: fetching s2 hit info\n", __func__); 
#endif
    
    const double Z_fPoint_Theta = evt.GetS2Hit(is_RHRS)->Z() - 1./M_Theta;     
    const double Z_fPoint_Phi   = evt.GetS2Hit(is_RHRS)->Z() - 1./M_Phi; 

#ifdef DEBUG_PAIR
    std::printf("<%s>: in body\n", __func__); 
#endif
    
    //____________________________________________________________________________________________________________________    
    auto Guess_slopes = [uFix, wVDC, 
            Z_fPoint_Theta, 
            Z_fPoint_Phi]( const double u, 
                    const double v, 
                    double &slope_u, 
                    double &slope_v ) { 
        TVector3 r_xyz 
        = ApexVDC::Track::Rotate_uvw_to_xyz( TVector3( u-uFix, v, wVDC ) ); 
        
        TVector3 S_uvw 
        = ApexVDC::Track::Rotate_xyz_to_uvw( TVector3( r_xyz[0]/(r_xyz[2]-Z_fPoint_Theta), 
                            r_xyz[1]/(r_xyz[2]-Z_fPoint_Phi),
                            1.) );       
        slope_u = S_uvw.Z() / S_uvw.X(); 
        slope_v = S_uvw.Z() / S_uvw.Y(); 
    }; 
    //____________________________________________________________________________________________________________________    
    

    //____________________________________________________________________________________________________________________    
    auto Guess_s2x = [evt, uFix, wVDC, Z_fPoint_Theta, is_RHRS]( const double u, const double v ) 
    { 
        TVector3 r_xyz = ApexVDC::Track::Rotate_uvw_to_xyz( TVector3( u-uFix, v, wVDC ) ); 
        return r_xyz.x()*(evt.GetS2Hit(is_RHRS)->Z() - Z_fPoint_Theta)/(r_xyz.Z() - Z_fPoint_Theta); 		 
    }; 
    //____________________________________________________________________________________________________________________
    
    const double CUT_vu_diff_min = is_LoChamber 
        ? (is_RHRS ? -83e-3 : -85e-3) 
        : (is_RHRS ? -95e-3 : -96e-3); 
    
    const double CUT_vu_diff_max = is_LoChamber 
        ? (is_RHRS ?  58e-3 :  58e-3)
        : (is_RHRS ?  70e-3 :  70e-3); 
    
    //____________________________________________________________________________________________________________________
    auto Find_clusters = [&evt, &Guess_slopes, CUT_vu_diff_min, CUT_vu_diff_max] ( 
        bool is_Uplane, 
        ROOT::RVec<ApexVDC::HitGroup>& groups 
    ) { 
    
        //how many wires over are we gonna look for tracks?
        const int   wire_buffer = 30; 
        
        //if clusters are closer than this, then join them together
        const int min_cluster_gap = 30; 
        
        ROOT::RVec<ApexVDC::HitCluster> clust_all; 
        
        
        for (int g=groups.size()-1; g>=0; g+=-1) { 
            //iterate backwards thru groups
            auto& group = groups.at(g); 
            
            //find span of our groupings
            int span
                = (int)( group.Span()/run_parameters::kGridSpacing ) + 2*wire_buffer + 1; 
                
            double x0 = group.LoEdge() - ((double)wire_buffer)*run_parameters::kGridSpacing; 
            
            for (int ix=0; ix<span; ix++) { 
                
                double x = x0 + ((double)ix)*run_parameters::kGridSpacing; 
                
                double y_min = x + (is_Uplane ? CUT_vu_diff_min : -CUT_vu_diff_min); 
                double y_max = x + (is_Uplane ? CUT_vu_diff_max : -CUT_vu_diff_max); 
                
                double m1, m2, m_dummy;  
                
                if (is_Uplane) { 
                    Guess_slopes( y_min,x, m_dummy,m1 ); 
                    Guess_slopes( y_max,x, m_dummy,m2 ); 
                } else         { 
                    Guess_slopes( x,y_min, m1,m_dummy ); 
                    Guess_slopes( x,y_max, m2,m_dummy ); 
                } 
                
                double Eta = grid_search( evt, group, 
                                m1,m2, 
                                x - run_parameters::kGridSpacing, 
                                x + run_parameters::kGridSpacing, 
                                9e-9, 
                                20e-9 ); 
                
                if (Eta > run_parameters::kCUT_minEta) 
                    clust_all.push_back( ApexVDC::HitCluster( &group, x, Eta ) ); 
                
            }//for (int ix=0; ix<span; ix++) 
        }//for (int g=0; g<groups.size(); g++) 
        
        
        //no clusters found! 
        if (clust_all.size()<1) return clust_all; 
        
        //now, make sure that that clusters aren't bunched up together
        ROOT::RVec<ApexVDC::HitCluster> clust_keep; 
        
        //now, prune the clusters to decide which hits to get rid of
        auto& bestClust = clust_all.front(); 
                
        double prev_x = bestClust.Intercept(); 
        
        for (int c=0; c<clust_all.size(); c++) { 
                
            auto clust = clust_all.at(c); 
            
            //how far is this cluster from the last one? 
            int gap = (int)TMath::Nint( (clust.Intercept()-prev_x)/run_parameters::kGridSpacing ); 
            
#if 0 //DEBUG	
            cout << TString::Format("c:%3i, int:%2.4f, eta:%1.3f, gap:%1.4f", 
                        c, 
                        clust->Intercept(), 
                        clust->Eta(), 
                        (clust->Intercept()-prev_x)/run_parameters::kGridSpacing ) << endl;
#endif

            if (gap > min_cluster_gap) { 
                
                clust_keep.push_back( bestClust ); 
                
                bestClust = clust; //start a new potential cluster grouping
            
            } else  { //this cluster is close to the last one
                
                if (clust.Eta() > bestClust.Eta()) //the new best cluster
                    bestClust = clust; 
                }
            prev_x = clust.Intercept(); 
        }
        clust_keep.push_back( bestClust ); 

#if 0 //DEBUG 
        cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~ kept ("<<clust_keep.size()<<"):"<<endl; 
        
        for (int c=0; c<clust_keep.size(); c++) 
        cout << TString::Format("c:%3i, int:%2.4f, eta:%1.3f", 
                    c, 
                    clust_keep.at(c)->Intercept(), 
                    clust_keep.at(c)->Eta() ) << endl; 
        
        cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl; 
#endif 
        
        return clust_keep; 
    }; 
    //____________________________________________________________________________________________________________________
    
    //for each possible pairing of groups, check to see if they
    // might be allowed to form a ApexVDC::ChamberPair, and, if so, go ahead and generate it. 
    ROOT::RVec<ApexVDC::ChamberPair> pairs; 
    
    ROOT::RVec<ApexVDC::HitCluster> clusters_u = Find_clusters( false, gVec_V ); 
    
    ROOT::RVec<ApexVDC::HitCluster> clusters_v = Find_clusters( true, gVec_U ); 
    
#ifdef DEBUG_PAIR
    std::printf("<%s>: beginning to pair clusters\n", __func__); 
#endif           

    
    //now, see if any of these clusters might be valid pairs
    for (int cv=0; cv<clusters_v.size(); cv++) {    
        auto& clust_v = clusters_v.at(cv); 
        
        for (int cu=0; cu<clusters_u.size(); cu++) {  
            auto& clust_u = clusters_u.at(cu); 

#ifdef DEBUG_PAIR
    std::printf("<%s>: cluster pair %i / %i\n", __func__, cv, cu); 
#endif           
            double vu_diff = clust_v.Intercept() - clust_u.Intercept(); 
            
            if ( vu_diff + 2.*run_parameters::kGridSpacing < CUT_vu_diff_min ||
                vu_diff - 2.*run_parameters::kGridSpacing > CUT_vu_diff_max  ) continue; 
            
            //now check to see if these clusters ACTUALLY match
            
            double m_u,m_v; 
            Guess_slopes( clust_u.Intercept(), clust_v.Intercept(), m_u, m_v );

#ifdef DEBUG_PAIR
    std::printf("<%s>: starting grid search (u) group nhits: %p\n ", __func__, clust_u.GetGroup()); 
#endif         
            double Eta_u = grid_search( 
                evt, *clust_u.GetGroup(), 
                m_u,m_u, 
                clust_u.Intercept() -run_parameters::kGridSpacing*2., 
                clust_u.Intercept() +run_parameters::kGridSpacing*2., 
                9e-9,    //TAU_sigma
                20e-9 
            ); //TAU_buffer 
            
            if (Eta_u < run_parameters::kCUT_minEta) continue; 

#ifdef DEBUG_PAIR
    std::printf("<%s>: starting grid search (v)\n", __func__); 
#endif         
            double Eta_v = grid_search( 
                evt, *clust_v.GetGroup(), 
                m_v,m_v, 
                clust_v.Intercept() -run_parameters::kGridSpacing*2., 
                clust_v.Intercept() +run_parameters::kGridSpacing*2., 
                9e-9,    //TAU_sigma
                20e-9 
            ); //TAU_buffer 
            
            if (Eta_v < run_parameters::kCUT_minEta) continue; 
            
            pair_unique_id++; //this will be used by tracks later
#ifdef DEBUG_PAIR
    std::printf("<%s>: keeping pair!\n", __func__); 
#endif         
            pairs.push_back( ApexVDC::ChamberPair(is_LoChamber, &clust_u, &clust_v, pair_unique_id++) ); 

        } //for (int cu=0; cu<clusters_u.size(); cu++) 
    } //for (int cv=0; cv<clusters_v.size(); cv++)
    
    return pairs; 
}; 

#endif