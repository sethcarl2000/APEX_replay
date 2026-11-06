

//544444444444444/"
//"""""
// - muon's comment 

#include <APEX/units.h>
#include <APEX/replay/helpers.h>
#include <APEX/run_parameters.h>
#include <APEX/VDC/Track.h> 
#include <math.h> 

namespace APEX
{
namespace replay
{
namespace helpers
{

void compute_trackdata(VDC::Track &trk) 
{ 

#if DEBUG 
    cout << "Compute_trackData():"; 
    //check if track exists 
    cout << "track = " << (trk ? "exists." : "does not exist!!!!!!!!" ) << endl; 
#endif 

    //good-point group
    const double CUT_goodPoint = 15.*units::ns; 
    const double measure_sigma = run_parameters::TRK_measureSigma; 

    VDC::for_each_plane([&]<VDC::Plane P>(){ 
        
        auto& group = trk.GetGroup<P>(); 
        
        double eta(0.); 
        double RMS(0.); 
        int    nGoodPoints(0); 
        
#if DEBUG
        cout << TString::Format("plane %2i ", p);   
#endif 
        
        for (int h=0; h<group.Nhits(); h++) { 

            double err = std::fabs( trk.Get_T_model<P>(group.WirePos(h)) + trk.T0() - group.Time(h) ); 

            if (err > CUT_goodPoint) continue; 

            //now, lets use this point (it's good, according to our cut)
            nGoodPoints++; 

            eta += std::exp( -0.5*utils::intpow<2>(err/measure_sigma) ); 
#if DEBUG 
            printf("h:%3i  wp=%3.3f t-T=%3.3f t-T0=%3.3f g-Time=%3.3f err=%3.3f",
                    h, 
                    group->WirePos(h), 
                    1e9*trk->Get_T_model(p,group->WirePos(h)),
                    1e9*trk->T0(),
                    1e9*group->Time(h), 
                    1e9*err 
                ); 
#endif 

            RMS += err*err; 
        }
        
#if DEBUG
        if ( eta==eta ) { cout << TString::Format("Eta= %0.3f", eta) << endl; }
        else            { cout << "Eta= NaN !!!!!!!!!!!!!!!!!!!!!!!" << endl; } 
#endif 
        
        trk.SetEta<P>( eta ); 
        
        trk.Set_nGoodPoints<P>( nGoodPoints ); 
        
    });//for (int p=0; p<4; p++) 

    return;  
};

} 
}
}