
#include <TCanvas.h> 
#include <TH2D.h>
#include "AnalyzeAllData.cxx"
#include <cmath> 

int draw_coinc_time()
{

  //
  AnalyzeAllData ana(20); 
  

  double ns = 1e-9; 

  //__________________________________________________________________________________________
  RDataframeUpdateFcn scale_time = [ns](ROOT::RDF::RNode df) {
    auto df_out = df
      .Redefine("S2R_S2L_dt_center", [ns](double t) { return t/ns; }, {"S2R_S2L_dt_center"})
      .Redefine("S2R_S2L_dt_sigma",  [ns](double t) { return t/ns; }, {"S2R_S2L_dt_sigma"});

    return df_out;
  }; 
  //__________________________________________________________________________________________

  //__________________________________________________________________________________________
  RDataframeUpdateFcn log_sigma = [ns](ROOT::RDF::RNode df) {
    auto df_out = df
      .Define("log10_S2R_S2L_dt_sigma",  [ns](double t) { return std::log10(std::fabs(t)); }, {"S2R_S2L_dt_sigma"});

    return df_out;
  }; 
  //__________________________________________________________________________________________

  
  auto hist_coinc = new TH2D("h_coinc", "T_{R} - T_{L} coinc time;run number;coinc time (ns)",
			     200, 3800,4800,
			     200, -100,+100);
  
  hist_coinc->SetDirectory(0); 
  
  ana.Fill_TH2D(hist_coinc, "run_number", "S2R_S2L_dt_center", &scale_time, "meta_data");


  auto hist_sigma = new TH2D("h_sigma", "T_{R} - T_{L} log_{10}[ #sigma ];run number;log_{10}( #sigma[R-L](s) )",
			     200, 3800,4800,
			     200, -12,-5);
  ana.Fill_TH2D(hist_sigma, "run_number", "log10_S2R_S2L_dt_sigma", &log_sigma, "meta_data"); 

  
  new TCanvas; 
  hist_coinc->SetMarkerStyle(kOpenCircle);
  hist_coinc->SetMarkerSize(0.70); 
  hist_coinc->Draw();

  new TCanvas;
  hist_sigma->SetMarkerStyle(kOpenCircle);
  hist_sigma->SetMarkerSize(0.70); 
  hist_sigma->Draw();

  
  return 0; 
}
