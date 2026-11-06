
// APEX
#include <APEX/replay/helpers.h>
#include <APEX/run_parameters.h>
#include <APEX/VDC/HitGroup.h>
#include <APEX/VDC/HitCluster.h>
#include <APEX/VDC/ChamberPair.h> 
#include <APEX/VDC/Track.h> 
// ROOT
#include <TVector3.h> 
// stdlib 
#include <cmath> 
#include <cstdio> 

namespace APEX
{
namespace replay
{
namespace helpers
{

namespace {
    //this involves angular prediction based on typical S2-angles. 
    const double uFix = 26.*units::mm/std::sqrt(2);
}; 

//#define DEBUG_PAIR

template<VDC::Chamber C> void gen_pairs( 
    S2Event& evt) 
{ 
#ifdef DEBUG_PAIR
#endif
    
#ifdef DEBUG_PAIR
    std::printf("<%s>: in body\n", __func__); 
#endif
    
    static constexpr bool is_loChamber = VDC::is_loChamber<C>(); 

    // these are the enums corresponding to the U and V planes 
    constexpr VDC::Plane U = VDC::get_U<C>(); 
    constexpr VDC::Plane V = VDC::get_V<C>(); 

#ifdef DEBUG_PAIR
    std::printf("<gen_pairs>: VDC Planes: %s / %s\n", VDC::plane_name<U>().c_str(), VDC::plane_name<V>().c_str()); 
#endif

    //constants
    bool is_RHRS = evt.Is_RHRS(); 
    
    // we're using the U-plane (lower plane) as the base w-value (the +w-axis is normal-up from each VDC plane)
    const double wVDC = VDC::w<U>(is_RHRS); 

    const double M_Theta = is_RHRS ? 0.1096 : 0.1096 ; 
    const double M_Phi   = is_RHRS ? 0.242  : 0.242  ; 

#ifdef DEBUG_PAIR
    std::printf("<%s>: fetching s2 hit info\n", __func__); 
#endif
    
    const double z_S2 = evt.GetS2Hit().Z(); 
    const double Z_fPoint_Theta = z_S2 - 1./M_Theta;     
    const double Z_fPoint_Phi   = z_S2 - 1./M_Phi; 

    //____________________________________________________________________________________________________________________    
    auto Guess_slopes = [wVDC, 
            Z_fPoint_Theta, 
            Z_fPoint_Phi]( 
                const double u, 
                const double v, 
                double &slope_u, 
                double &slope_v 
    ) { 
        TVector3 r_xyz = VDC::Track::Rotate_uvw_to_xyz( TVector3( u-uFix, v, wVDC ) ); 
        
        TVector3 S_uvw = VDC::Track::Rotate_xyz_to_uvw( TVector3( 
            r_xyz.X()/(r_xyz.Z()-Z_fPoint_Theta), 
            r_xyz.Y()/(r_xyz.Z()-Z_fPoint_Phi),
            1.
        ));       

        slope_u = S_uvw.Z() / S_uvw.X(); 
        slope_v = S_uvw.Z() / S_uvw.Y(); 
    }; 
    //____________________________________________________________________________________________________________________    
    

    const double CUT_vu_diff_min = is_loChamber 
        ? (is_RHRS ? -83.*units::mm : -85.*units::mm) 
        : (is_RHRS ? -95.*units::mm : -96.*units::mm); 
    
    const double CUT_vu_diff_max = is_loChamber 
        ? (is_RHRS ?  58.*units::mm :  58.*units::mm)
        : (is_RHRS ?  70.*units::mm :  70.*units::mm); 
    
    //____________________________________________________________________________________________________________________
    auto find_clusters = [&evt, &Guess_slopes, CUT_vu_diff_min, CUT_vu_diff_max]<VDC::Plane P>( 
        const ObjCollection<VDC::HitGroup<P>>& groups 
    ) { 
        static constexpr bool is_Uplane = ((P==VDC::Plane::U1) || (P==VDC::Plane::U2));

        //how many wires over are we gonna look for tracks?
        const int   wire_buffer = 30; 
        
        //if clusters are closer than this, then join them together
        const int min_cluster_gap = 30; 
        
        ROOT::RVec<VDC::HitCluster<P>> clust_all; 
        clust_all.reserve(groups.size());
        
        for (const auto& group : groups) { 
            
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
                                9.*units::ns, 
                                20.*units::ns 
                            ); 
                    
                if (Eta > run_parameters::kCUT_minEta) 
                    clust_all.emplace_back( group, x, Eta ); 
                
            }//for (int ix=0; ix<span; ix++) 
        }//for (int g=0; g<groups.size(); g++) 
        
        //no clusters found! 
        if (clust_all.empty()) { 
#ifdef DEBUG_PAIR
            std::printf("<gen_pairs: find_clusters %s>: No clusters found\n", VDC::plane_name<P>().c_str()); 

#endif  
            return; 
        }
        
        //now, make sure that that clusters aren't bunched up together
        ROOT::RVec<VDC::HitCluster<P>> clust_keep; 
        clust_keep.reserve(clust_all.size()); 
        
        //now, prune the clusters to decide which hits to get rid of
        auto best_clust = &clust_all.front(); 
                
        double prev_x = best_clust->GetIntercept(); 
        
        for (auto& clust : clust_all) { 
               
            //how far is this cluster from the last one? 
            int gap = (int)std::round( (clust.GetIntercept()-prev_x)/run_parameters::kGridSpacing ); 
            
            if (gap > min_cluster_gap) { 
                
                evt.emplace<VDC::HitCluster<P>>( *best_clust ); 
                
                best_clust = &clust; //start a new potential cluster grouping
            
            } else  { //this cluster is close to the last one
                
                //the new best cluster
                if (clust.GetEta() > best_clust->GetEta()) best_clust = &clust; 
            }
            
            prev_x = clust.GetIntercept(); 
            
        }

        // add this cluster to the running list for the event 
        evt.emplace<VDC::HitCluster<P>>( *best_clust ); 
#ifdef DEBUG_PAIR
        std::printf("<gen_pairs: find_clusters %s>: N. 'raw' clusters: %zu\n", 
            VDC::plane_name<P>().c_str(), 
            evt.size<VDC::HitCluster<P>>()
        );
#endif
    }; 
    //____________________________________________________________________________________________________________________
    
    //for each possible pairing of groups, check to see if they
    // might be allowed to form a VDC::ChamberPair, and, if so, go ahead and generate it. 
    ROOT::RVec<VDC::ChamberPair<C>> pairs;  
#ifdef DEBUG_PAIR
    std::printf("<%s>: Finding Clusters <%s>\n", __func__, VDC::plane_name<U>().c_str()); 
#endif
    find_clusters( evt.GetCollection<VDC::HitGroup<U>>() ); 
#ifdef DEBUG_PAIR
    std::printf("<%s>: Finding Clusters <%s>\n", __func__, VDC::plane_name<V>().c_str()); 
#endif
    find_clusters( evt.GetCollection<VDC::HitGroup<V>>() ); 


    //now, see if any of these clusters might be valid pairs
    for (auto& clust_U : evt.GetCollection<VDC::HitCluster<U>>()) {    
        for (auto& clust_V : evt.GetCollection<VDC::HitCluster<V>>()) {

            double u = clust_V.GetIntercept(); 
            double v = clust_U.GetIntercept(); 
#ifdef DEBUG_PAIR
            std::printf("<gen_pairs>: <U/V: %zu/%zu>: cluster interecepts: u = %+6.1f mm   v = %+6.1f mm    v-u = %+.1f mm\n", 
                clust_U.ID(), clust_V.ID(), u/units::mm, v/units::mm, (v-u)/units::mm
            );
#endif
            double vu_diff = v - u; 
            
            if ( 
                vu_diff + 2.*run_parameters::kGridSpacing < CUT_vu_diff_min ||
                vu_diff - 2.*run_parameters::kGridSpacing > CUT_vu_diff_max  ) continue; 
            
            //now check to see if these clusters ACTUALLY match
            
            double m_u,m_v; 
            Guess_slopes( u, v, m_u, m_v );

#ifdef DEBUG_PAIR
            std::printf("<%s>: starting grid search (u) group nhits: %i\n ", __func__, clust_V.GetGroup().Nhits()); 
#endif         
            double Eta_u = grid_search( 
                evt, clust_V.GetGroup(), 
                m_u,m_u, 
                u - run_parameters::kGridSpacing*2., 
                u + run_parameters::kGridSpacing*2., 
                9.*units::ns,    //TAU_sigma
                20.*units::ns 
            ); //TAU_buffer 
#ifdef DEBUG_PAIR
            std::printf("<%s>: Eta:  %5.3f\n ", __func__, Eta_u); 
#endif             
            if (Eta_u < run_parameters::kCUT_minEta) continue; 

#ifdef DEBUG_PAIR
    std::printf("<%s>: starting grid search (v) group nhits: %i\n ", __func__, clust_U.GetGroup().Nhits()); 
#endif        
            double Eta_v = grid_search( 
                evt, clust_U.GetGroup(), 
                m_v,m_v, 
                v - run_parameters::kGridSpacing*2., 
                v + run_parameters::kGridSpacing*2., 
                9.*units::ns,    //TAU_sigma
                20.*units::ns 
            ); //TAU_buffer 
#ifdef DEBUG_PAIR
            std::printf("<%s>: Eta:  %5.3f\n ", __func__, Eta_v); 
#endif             
            if (Eta_v < run_parameters::kCUT_minEta) continue; 
            
#ifdef DEBUG_PAIR
            std::printf("<%s>: keeping pair!\n", __func__); 
#endif         
            evt.emplace<VDC::ChamberPair<C>>( clust_U, clust_V ); 

        } //for (int cu=0; cu<clusters_u.size(); cu++) 
    } //for (int cv=0; cv<clusters_v.size(); cv++)

}; 

template void gen_pairs<VDC::Chamber::Low> (S2Event&); 
template void gen_pairs<VDC::Chamber::High>(S2Event&); 

}
}
}