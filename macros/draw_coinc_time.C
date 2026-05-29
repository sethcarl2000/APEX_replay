
#include <TCanvas.h> 
#include <TH2D.h>
#include "AnalyzeAllData.cxx"

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

  auto hist_coinc = new TH2D("h_coinc", "T_{R} - T_{L} coinc time;run number;coinc time (ns)",
			     200, 3800,4800,
			     200, -100,+100);

  hist_coinc->SetDirectory(0); 

  ana.Fill_TH2D(hist_coinc, "run_number", "S2R_S2L_dt_center", &scale_time, "meta_data");

  hist_coinc->SetMarkerStyle(kOpenCircle);
  hist_coinc->SetMarkerSize(0.70); 
  hist_coinc->Draw();

  return 0; 
}
