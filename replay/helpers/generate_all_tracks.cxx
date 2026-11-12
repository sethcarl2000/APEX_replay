#include <APEX/replay/helpers.h>
#include <APEX/VDC/Track.h>
#include <APEX/S2Hit.h>
#include <APEX/S2Event.h>
#include <APEX/run_parameters.h>

namespace APEX
{
namespace replay
{
namespace helpers
{

void generate_all_tracks(S2Event& S2_event)
{
    const auto& S2_hit = S2_event.GetS2Hit(); 

    const bool is_RHRS = S2_event.Is_RHRS(); 

    const double CUT_ph_min = is_RHRS ? run_parameters::CUT_ph_min_R : run_parameters::CUT_ph_min_L; 
    const double CUT_ph_max = is_RHRS ? run_parameters::CUT_ph_max_R : run_parameters::CUT_ph_max_L;

    const double CUT_th_min = is_RHRS ? run_parameters::CUT_th_min_R : run_parameters::CUT_th_min_L;
    const double CUT_th_max = is_RHRS ? run_parameters::CUT_th_max_R : run_parameters::CUT_th_max_L;
   
    const double x_param_offset = is_RHRS ? run_parameters::TRK_R_xParam_offset : run_parameters::TRK_L_xParam_offset; 

    const double dt_offset = is_RHRS ? run_parameters::TRK_R_Dt_offset : run_parameters::TRK_L_Dt_offset; 
    const double dt_cut    = is_RHRS ? run_parameters::TRK_R_CUT_Dt : run_parameters::TRK_L_CUT_Dt; 

    constexpr auto Lo = VDC::Chamber::Low; 
    constexpr auto Hi = VDC::Chamber::High; 

    // pre-size the array of all tracks
    S2_event.reserve<VDC::Track>( 
        S2_event.size<VDC::ChamberPair<Lo>>() * S2_event.size<VDC::ChamberPair<Hi>>() 
    ); 

    static constexpr double rawtrack_mult = run_parameters::TRK_rawtrack_cut_multiplier; 

    // create all raw tracks
    for (auto& pair_lo : S2_event.GetCollection<VDC::ChamberPair<Lo>>()) {
        for (auto& pair_hi : S2_event.GetCollection<VDC::ChamberPair<Hi>>()) {

            VDC::Track trk( &S2_event, pair_lo, pair_hi );

            double err_Theta = trk.Theta() - Theta_model(trk); 
            double err_Phi   = trk.Phi()   - Phi_model(trk); 
                
            if ( err_Theta < CUT_th_min * rawtrack_mult || 
                    err_Theta > CUT_th_max * rawtrack_mult || 
                    err_Phi   < CUT_ph_min * rawtrack_mult || 
                    err_Phi   > CUT_ph_max * rawtrack_mult 
            ) continue; 

            //do some basic checks
            double xParam = trk.xParam();
            double error_y = trk.S2_y() - trk.GetEvent()->GetS2Hit().Y(); 

            if (std::fabs(xParam-x_param_offset) > run_parameters::TRK_CUT_xParam * rawtrack_mult ||  // track x-alignment with paddle
                std::fabs(error_y) > run_parameters::TRK_S2_y_CUT * rawtrack_mult                     // track y-alignment with paddle
            ) continue; 

            S2_event.emplace<VDC::Track>( std::move(trk) ); 
        }
    }

    // now, we can sort through our tracks, and see which we should reject. 
    // do this to make sure each track has the correct event pointer. 
    S2_event.RefereshEventPointers(); 

    const double tau_sigma = S2_event.GetEventHandler().Get_tauSigma(); 
    
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

    // see if this track's projection onto the S2 matches up well enough 
    tracks.remove_if([=](const VDC::Track& trk){    
        
        //do some basic checks
        double xParam = trk.xParam();
        double Dt     = trk.T0(); 
        double error_y = trk.S2_y() - trk.GetEvent()->GetS2Hit().Y(); 

        if ( 
            std::fabs(xParam-x_param_offset) > run_parameters::TRK_CUT_xParam ||  // track x-alignment with paddle
            std::fabs(error_y) > run_parameters::TRK_S2_y_CUT ||                  // track y-alignment with paddle
            std::fabs(Dt-dt_offset) > dt_cut                                      // track t-agreement with paddle
        ) 
        { 
            //discard this track
            return true;  
        } 
        return false; 
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

    //compute track errors (for tracks surviving all cuts)
    for (auto& trk : tracks) {
        Compute_trackError( trk ); 
    }
}

}
}
}
