
// APEX
#include <APEX/replay/helpers.h>
#include <APEX/run_parameters.h>
#include <APEX/VDC/HitGroup.h>
#include <APEX/VDC/HitCluster.h>
#include <APEX/VDC/ChamberPair.h> 
#include <APEX/VDC/Track.h> 
#include <APEX/utils.h>
// ROOT
#include <TVector3.h> 
// stdlib 
#include <cmath> 
#include <cstdio> 
#include <vector> 

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

template<VDC::Chamber C> std::vector<pair_candidate> chamber_grid_search( 
    S2Event& evt, 
    const double min_cluster_gap, 
    const double min_eta
) { 
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
    const double VDC_w = 0.5*( VDC::w<U>(is_RHRS) + VDC::w<V>(is_RHRS) ); 

    const double CUT_vu_diff_min = is_loChamber 
        ? (is_RHRS ? -83.*units::mm : -85.*units::mm) 
        : (is_RHRS ? -95.*units::mm : -96.*units::mm); 
    
    const double CUT_vu_diff_max = is_loChamber 
        ? (is_RHRS ?  58.*units::mm :  58.*units::mm)
        : (is_RHRS ?  70.*units::mm :  70.*units::mm); 


    // how much could the 'true' intercept of a track 'hang over' the edge of a group? 
    // this is a reasonable upper-limit.
    const double overlap_buffer = ((double)3.)*VDC::kWireSpacing; 

    
    // grid-spacing in u and v 
    static constexpr double dx = run_parameters::kGridSpacing; 

    static constexpr double eta_threshold = run_parameters::kCUT_minEta; 

    // get position of S2 hit
    const auto& S2_hit = evt.GetS2Hit(); 
    TVector3 S2_int_xyz( S2_hit.X(), S2_hit.Y(), S2_hit.Z() ); 

    const TVector3 S2_int_uvw = VDC::Track::Rotate_xyz_to_uvw( S2_int_xyz );

    const double S2_u{S2_int_uvw.X()}, S2_v{S2_int_uvw.Y()}, S2_w{S2_int_uvw.Z()}; 
    
    std::vector<pair_candidate> all_pairs; 

    // do brute-force grid search with all groups 
    for (const auto& group_U : evt.GetCollection<VDC::HitGroup<U>>()) {
        for (const auto& group_V : evt.GetCollection<VDC::HitGroup<V>>()) {
            
            //first, check if these groups overlap at all 
            double vu_diff_min = group_U.LoEdge() - group_V.HiEdge() - overlap_buffer; 
            double vu_diff_max = group_U.HiEdge() - group_V.LoEdge() + overlap_buffer; 

            // if these groups don't overlap at all, then continue 
            if (!(vu_diff_min < CUT_vu_diff_max || vu_diff_max > CUT_vu_diff_min)) continue;   

            //do a brute-force grid-search
            const int n_steps_u = 1 + static_cast<int>((2.*overlap_buffer + group_V.Span())/dx); 
            const int n_steps_v = 1 + static_cast<int>((2.*overlap_buffer + group_U.Span())/dx); 
            
            std::vector<pair_candidate> pairs; pairs.reserve(n_steps_u * n_steps_v); 

            /// keep pair candidate (create new cluster, or combine into exsisting cluster)
            //_______________________________________________________________________________________________________________
            auto keep_candidate = [&pairs, min_cluster_gap](double u, double v, double eta_u, double eta_v) {

                pair_candidate cand; 
                cand.u = u; 
                cand.v = v; 
                cand.eta_u = eta_u; 
                cand.eta_v = eta_v; 

                for (auto& pair : pairs) {

                    //is this candidate within this cluster's boundary? 
                    if (pair.boundary_u.within(cand.u) && pair.boundary_v.within(cand.v)) {
                        
                        // this pair is in the 'boundary' of an existing cluster. 
                        // check if this new pair is a better candidate then the existing one. 
                        if (cand.min_eta() > pair.min_eta()) {

                            pair.u = cand.u; 
                            pair.v = cand.v; 

                            pair.eta_u = cand.eta_u; 
                            pair.eta_v = cand.eta_v; 
                        }

                        // set the new boundary of this cluster
                        pair.boundary_u.min = std::min( pair.boundary_u.min, cand.u - min_cluster_gap ); 
                        pair.boundary_u.max = std::max( pair.boundary_u.max, cand.u + min_cluster_gap ); 

                        pair.boundary_v.min = std::min( pair.boundary_v.min, cand.v - min_cluster_gap ); 
                        pair.boundary_v.max = std::max( pair.boundary_v.max, cand.v + min_cluster_gap ); 

                        // exit (we've found a matching cluster)
                        return; 
                    }
                }      

                // if this doesn't match any existing clusters, then create a new one. 
                cand.boundary_u = Bound<double>{
                    u - run_parameters::min_cluster_gap,
                    u + run_parameters::min_cluster_gap
                }; 
                cand.boundary_v = Bound<double>{
                    v - run_parameters::min_cluster_gap,
                    v + run_parameters::min_cluster_gap
                }; 
                pairs.emplace_back( std::move(cand) );
                return; 
            }; 
            //_______________________________________________________________________________________________________________
            
            for (int iu=0; iu<n_steps_u; iu++) {
                double u = group_V.LoEdge() - overlap_buffer + ((double)iu)*dx; 
            
                for (int iv=0; iv<n_steps_v; iv++) {
                    double v = group_U.LoEdge() - overlap_buffer + ((double)iv)*dx; 
                    
                    // is this u,v pair feasible? 
                    if (v-u - overlap_buffer > CUT_vu_diff_max || v-u + overlap_buffer < CUT_vu_diff_min) continue; 

                    // distance to S2 
                    double dist_to_S2 = std::sqrt( 
                        utils::intpow<2>(S2_u - u) + 
                        utils::intpow<2>(S2_v - v) + 
                        utils::intpow<2>(S2_w - VDC_w) 
                    ); 

                    // time of flight for this cluster (form the midpoint of the two VDC planes to the S2 hit)
                    double ToF = dist_to_S2/units::speed_of_light_vacuum; 

                    // get slopes
                    double dw_du = (S2_w - VDC::w<V>(is_RHRS)) / (S2_u - u); 
                    double dw_dv = (S2_w - VDC::w<U>(is_RHRS)) / (S2_v - v); 
                    
                    // grid-search in 'u' and 'v'
                    double eta_u = grid_search(
                        evt, group_V, dw_du,dw_du, u - dx/2., u + dx/2., 20.*units::ns, 0., -ToF
                    );
                    double eta_v = grid_search(
                        evt, group_U, dw_dv,dw_dv, v - dx/2., v + dx/2., 20.*units::ns, 0., -ToF
                    ); 
#ifdef DEBUG_PAIR
                    std::printf("u,v %+5.1f mm, %+5.1f mm,  dw/du: %.3f,  dw/dv: %.3f,  eta_u/eta_v:  %.2f/%.2f\n",
                        u/units::mm, v/units::mm, dw_du, dw_dv, eta_u,eta_v
                    );
#endif
                    if (eta_u > min_eta && eta_v > min_eta) {

                        keep_candidate(u,v, eta_u,eta_v); 
                    }
                }
            }

            // now, for all the cluster-candidates which remain, fill their info into the event. 
            for (auto & pair : pairs) {
                evt.emplace<VDC::ChamberPair<C>>(
                    VDC::HitCluster<U>(group_U, pair.v, pair.eta_v), 
                    VDC::HitCluster<V>(group_V, pair.u, pair.eta_u)
                ); 
            }

            // remove this later, it's only for debug purposes
        } // for (const auto& group_V ... 
    } // for (const auto& group_U ... 

    return {};
}; 

template std::vector<pair_candidate> chamber_grid_search<VDC::Chamber::Low> (S2Event&, double,double); 
template std::vector<pair_candidate> chamber_grid_search<VDC::Chamber::High>(S2Event&, double,double); 

}
}
}