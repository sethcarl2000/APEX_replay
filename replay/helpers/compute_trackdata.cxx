

//544444444444444/"
//"""""
// - muon's comment 

#include "APEX_replay_helpers.h"
#include <run_parameters.h>
#include <ApexVDCTrack.h> 
#include <math.h> 

namespace APEX
{
namespace replay
{
namespace helpers
{

void compute_trackdata(ApexVDC::Track &trk) 
{ 

#if DEBUG 
    cout << "Compute_trackData():"; 
    //check if track exists 
    cout << "track = " << (trk ? "exists." : "does not exist!!!!!!!!" ) << endl; 
#endif 

    //good-point group
    const double CUT_goodPoint = 15e-9; 
    const double measure_sigma = run_parameters::TRK_measureSigma; 

    for (int p=0; p<4; p++) { 
        
        auto group = trk.GetGroup(p); 
        
        double eta(0.); 
        double RMS(0.); 
        int    nGoodPoints(0); 
        
#if DEBUG
        cout << TString::Format("plane %2i ", p);   
#endif 
        
        for (int h=0; h<group->Nhits(); h++) { 

            double err = std::fabs( trk.Get_T_model(p,group->WirePos(h)) + trk.T0() - group->Time(h) ); 

            if (err > CUT_goodPoint) continue; 

            //now, lets use this point (it's good, according to our cut)
            nGoodPoints++; 

            eta += TMath::Exp( -0.5*TMath::Power(err/measure_sigma, 2) ); 
#if DEBUG 
            cout 
                << TString::Format("h:%3i  wp=%3.3f t-T=%3.3f t-T0=%3.3f g-Time=%3.3f err=%3.3f",
                            h, 
                            group->WirePos(h), 
                            1e9*trk->Get_T_model(p,group->WirePos(h)),
                            1e9*trk->T0(),
                            1e9*group->Time(h), 
                            1e9*err )            << endl; 
#endif 

            RMS += err*err; 
        }
        
#if DEBUG
        if ( eta==eta ) { cout << TString::Format("Eta= %0.3f", eta) << endl; }
        else            { cout << "Eta= NaN !!!!!!!!!!!!!!!!!!!!!!!" << endl; } 
#endif 
        
        trk.Set_Eta( p, eta ); 
        
        trk.Set_nGoodPoints( p, nGoodPoints ); 
        
        if (nGoodPoints > 0) { 

            trk.Set_RMS( p, TMath::Sqrt( RMS/((double)nGoodPoints) ) ); 

        } else { trk.Set_RMS( p, -1e30 ); }
        
    }//for (int p=0; p<4; p++) 

    return;  
};

} 
}
}