
// APEX
#include <APEX/replay/helpers.h>
#include <APEX/VDC/Track.h>
#include <APEX/VDC.h> 
#include <APEX/utils.h> 
// ROOT
#include <ROOT/RVec.hxx>
// stdlib
#include <cmath> 
#include <span> 
#include <array> 
#include <algorithm> 

namespace APEX
{
namespace replay
{
namespace helpers
{

//
// Refine tracks from the hi-chamber using newton's method
//   
void refine_track( 
    VDC::Track& trk, 
    const int nCycles, 
    double sigma) 
{ 
    /// @return 'true' if x is nan
    auto is_nan = [](double x) { return x != x; }; 

    static constexpr int nudge_size = VDC::kN_planes + 1; 

    //_________________________________________________________________________________________
    //check to see if a nudge is 'reasonable' i.e., not 'inf' or 'nan'
    auto Check_nudge = [&is_nan](std::span<const double, nudge_size> nudge) { 
        
        const double maxNudge_X = 5e-2; //no nudge should ever be this large
        
        const double maxNudge_T = 3e-7; 
        
        for (int p=0; p<4; p++) { 
            if (is_nan(nudge[p]) || std::fabs(nudge[p]) > maxNudge_X) return false; 
        }
        
        if (is_nan(nudge[4]) || std::fabs(nudge[4]) > maxNudge_T) return false; 
        return true; 
    }; 
    //_________________________________________________________________________________________
            
    //random-walk track minimizaiton
    const double GRAD_momentum =0.50; 
    const double GRAD_step0    =0.05; 
    const double GRAD_exponent =0.33; 

    //cout << "refining track..." << endl; 
    const double scale_T =3e-9;   
    const double scale_X =2.2e-3; 

    const double nudge_multiplier =1.; 

    const double Chi_cutoff  =6.; 

    static constexpr  double sigma_decay =0.950;
    static constexpr  double measure_sigma   =5.*units::ns; //sigma*2.; 

    double final_eta[4] = {0.}; 

    const double s2 = sigma*sigma; 

    ROOT::RVec<double> tVec; 

    std::array<double,nudge_size> nudge, F, Deriv; 
    std::array<double,4> J_ii, J_4i, objective_eta; 
    double J_44; 

    for (uint c=0; c<nCycles; c++) { 

        //this will be how we 'nudge' the parameters
        
        J_44 =0.; 
        // reset all arrays to 0 
        for (auto arr : {&nudge, &F, &Deriv}) arr->fill(0.);  
        for (auto arr : {&J_ii, &J_4i, &objective_eta}) arr->fill(0.);  
            
        VDC::for_each_plane([&]<VDC::Plane P>(){

            const auto& group = trk.GetGroup<P>(); 

            constexpr int p = VDC::number<P>(); 

            for (uint h=0; h<group.Nhits(); h++) { 
                    
                /*/draw event picture ///////////////////////////////////////////////
                h2d->Fill( group->WirePos(h)-int0[p], group->Time(h) ); 	  
                //*/////////////////////////////////////////////////////////////////
                
                double Chi = trk.Get_T_model<P>(group.WirePos(h)) + trk.T0() -  group.Time(h); 
                
                if ( std::fabs(Chi) > Chi_cutoff*sigma ) continue; 
                
                double Eta = std::exp( -0.5*utils::intpow<2>(Chi/sigma) ); 
                
                double m = trk.Slope<P>(); 
                
                F[p] += -m * Eta * Chi * trk.Get_T_model<P>(group.WirePos(h),1); 
                
                F[4] +=  Eta * Chi; 

                J_ii[p] 
                += (Chi * trk.Get_T_model<P>(group.WirePos(h),2)  
                + utils::intpow<2>(trk.Get_T_model<P>(group.WirePos(h),1)))*m*m*Eta; 
                
                J_4i[p] 
                +=  -Eta * m * trk.Get_T_model<P>(group.WirePos(h),1); 
                
                J_44 += Eta; 
                
                objective_eta[p] 
                += std::exp( -0.5*utils::intpow<2>(Chi/measure_sigma) ); 
                
            }
            
        }); //for (int p=0; p<4; p++) 
        
        /*///////////////////////////////////////////////////////////////////////////
        canv->Print("~/ftp-dump/plots/test_track_hi_Dt.gif+33");
        canv->~TCanvas(); 
        //*//////////////////////////////////////////////////////////////////////////
        
        sigma *= sigma_decay; 
        
        for (int p=0; p<4; p++) { 

            J_44  +=  J_4i[p]*(-J_4i[p]/J_ii[p]); 
            F[4]  +=  F[p]   *(-J_4i[p]/J_ii[p]); 
        }
        
        nudge[4] = -nudge_multiplier*F[4]/J_44; 
        
        for (int p=0; p<4; p++) { 

            nudge[p] = -nudge_multiplier*(F[p] - J_4i[p]*nudge[4])/J_ii[p]; 
        }
        
        //check to see if this nudge is reasonalbe
        if (!Check_nudge(nudge)) break; 
        
        trk.Nudge_params( nudge ); 
        
        trk.SetEta( objective_eta ); 
    }
};   

}
}
}