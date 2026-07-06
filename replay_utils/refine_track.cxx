#ifndef refine_track_H
#define refine_track_H

#include "replay_utils.h"
#include <ApexVDCTrack.h>
#include <math.h> 
#include <ROOT/RVec.hxx>
#include <ApexUtils.h> 

namespace replay_utils
{

//
// Refine tracks from the hi-chamber using newton's method
//   
void refine_track( 
    ApexVDC::Track& trk, 
    const int nCycles, 
    double sigma) 
{ 
    using APEX::square; 

    /// @return 'true' if x is nan
    auto is_nan = [](double x) { return x != x; }; 

    //_________________________________________________________________________________________
    //check to see if a nudge is 'reasonable' i.e., not 'inf' or 'nan'
    auto Check_nudge = [&is_nan](const double nudge[5]) { 
        
        const double maxNudge_X = 5e-2; //no nudge should ever be this large
        
        const double maxNudge_T = 3e-7; 
        
        for (int p=0; p<4; p++) { 
            if (is_nan(nudge[p])) return false; 
            if (std::fabs(nudge[p]) > maxNudge_X) return false; 
        }
        
        if (is_nan(nudge[4])) return false; 
        if (std::fabs(nudge[4]) > maxNudge_T) return false; 
        
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

    const double sigma_decay =0.950;
    const double measure_sigma   =5e-9; //sigma*2.; 

    ///////////////////////////////////////////////////////////////////////////////
    double int0[5]; 
    for (int p=0; p<4; p++) int0[p] = trk.Intercept(p); 
    //TString plane_name[4] = {"U1","V1","U2","V2"}; 
    int0[4] = trk.T0(); 
    ///////////////////////////////////////////////////////////////////////////////

    double final_eta[4] = {0.}; 

    const double s2 = sigma*sigma; 

    ROOT::RVec<double> tVec; 

    for (uint c=0; c<nCycles; c++) { 
        
        double objective_eta[4] = {0.};
        //this will be how we 'nudge' the parameters
        double nudge[5] = {0.}; 
        
        double Deriv[5] = {0.}; 
        
        //first-derivatives
        double F[5] = {0.}; 
        
        double J_ii[4] = {0.}; 
        double J_4i[4] = {0.}; 
        double J_44(0.); 
    
        for (int p=0; p<4; p++) { 

            /*/draw event picture ///////////////////////////////////////////////////
            canv->cd( iCanvas );  iCanvas++; 
            auto h2d = new TH2D(TString::Format("vdcPlane_%i",p), "", 
                                200, -25e-3, 25e-3, 
                                200, VDC_min_realTime, VDC_max_realTime); 
            //*//////////////////////////////////////////////////////////////////////

            ApexVDC::HitGroup* group = trk.GetGroup(p); 

            for (uint h=0; h<group->Nhits(); h++) { 
                    
                /*/draw event picture ///////////////////////////////////////////////
                h2d->Fill( group->WirePos(h)-int0[p], group->Time(h) ); 	  
                //*/////////////////////////////////////////////////////////////////
                
                double Chi = trk.Get_T_model(p, group->WirePos(h)) +  trk.T0() -  group->Time(h); 
                
                if ( std::fabs(Chi) > Chi_cutoff*sigma ) continue; 
                
                double Eta = TMath::Exp( -0.5*square(Chi/sigma) ); 
                
                double m = trk.Slope(p); 
                
                F[p] += -m * Eta * Chi * trk.Get_T_model(p,group->WirePos(h),1); 
                
                F[4] +=  Eta * Chi; 

                J_ii[p] 
                += (Chi * trk.Get_T_model(p,group->WirePos(h),2)  
                + square(trk.Get_T_model(p,group->WirePos(h),1)))*m*m*Eta; 
                
                J_4i[p] 
                +=  -Eta * m * trk.Get_T_model(p,group->WirePos(h),1); 
                
                J_44 += Eta; 
                
                objective_eta[p] 
                += TMath::Exp( -0.5*square(Chi/measure_sigma) ); 
                
            }//for (uint h=0; h<trk->pGet_N(p); 

            /*/////////////////////////////////////////////////////////////////////////
            h2d->SetMarkerStyle(kOpenCircle); 
            h2d->SetMarkerSize(0.50); 
            h2d->DrawCopy(); 

            auto ltx = new TLatex; 

            ltx->DrawLatex( -23e-3, VDC_max_realTime, 
                    plane_name[p] ); 
                    
            const int nPts = 100; 

            double max_X = 14e-3; 
            double X[nPts]; 
            double T[nPts]; 

            double dx = 2.*max_X/((double)nPts-1); 

            double xx=(trk->Intercept(p)-int0[p])-max_X; 
            for (int i=0; i<nPts; i++) { 
                X[i] = xx;  xx += dx; 
                T[i] = trk->Get_T_model( p, xx+int0[p] ) + trk->T0(); 
            }	
            auto g = new TGraph(nPts, X, T); 
            g->SetLineColor(kRed); 
            g->Draw("SAME"); 

            h2d->~TH2D(); 
            ltx->~TLatex(); 
            //*////////////////////////////////////////////////////////////////////////
            
        }//for (int p=0; p<4; p++) 
        
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
        
        trk.Set_Eta( objective_eta ); 
    }
};   

};

#endif