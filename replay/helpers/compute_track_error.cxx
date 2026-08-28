
#include <APEX/replay/helpers.h>
#include <APEX/run_parameters.h>
#include <APEX/VDC/Track.h> 
#include <APEX/utils.h> 
#include <math.h> 
#include <Math/QuantFuncMathCore.h> 

namespace APEX
{
namespace replay
{
namespace helpers
{

void Compute_trackError(VDC::Track& trk) 
{
    using namespace APEX; 

    const double ERR_tau_sigma = trk.GetEvent()->Get_tauSigma(); 

    double param_error[5]; unsigned int total_pts(0); 

    for (int p=0; p<4; p++) { 
        
        double S0(0.), S1(0.); 

        VDC::HitGroup* group = trk.GetGroup(p); 
        
        for (int h=0; h<group->Nhits(); h++) { 
                
            double err = trk.Get_T_model( p, group->WirePos(h) ) + trk.T0() - group->Time(h); 

            if ( std::fabs(err)/ERR_tau_sigma > 4. ) continue; 

            total_pts++; 
            
            S0 += utils::intpow<2>( trk.Get_T_model(p,group->WirePos(h),1) );

            S1 += trk.Get_T_model(p, group->WirePos(h), 1) * trk.Get_T_model(p, group->WirePos(h), 2); 
        }
        S0 = std::sqrt(S0); 
        
        double v_square(0.); 
        //sample all possible values of i
        const double dp = 1./((double)run_parameters::kGausIntPoints); 
        double cum_dp=0.; 
        for (int i=0; i<run_parameters::kGausIntPoints; i++) { 

            //cout << z[i] << " " << endl; 
            double zi = ROOT::Math::normal_quantile( (cum_dp += dp) - 0.5*dp, 1. );

            double v2 = ERR_tau_sigma*zi*S0 / ( ERR_tau_sigma*zi*S1/S0 + S0*S0 ); 

            v_square += utils::intpow<2>( v2/trk.Slope(p) ) * dp ; 
        }
        
        param_error[p] = std::sqrt(v_square); 
        //cout << TString::Format("  -  %3.3f", param_error[p]*1e3) << flush; 
        
    }// for (int p=0; p<4; p++) 

    //error of T0 
    param_error[4] = ERR_tau_sigma/std::sqrt( (double)total_pts ); 

    trk.Set_Errors( param_error ); 
}

}
}
}
