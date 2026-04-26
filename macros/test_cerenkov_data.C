
#include <ROOT/RDataFrame.hxx>
#include <TH1D.h>
#include <math.h>
#include <TCanvas.h> 
#include <TStyle.h> 
#include <TLegend.h> 

using namespace std;
using namespace ROOT::VecOps;

namespace {

  constexpr int n_paddles=10;
}

int test_cerenkov_data(const char* path_infile)
{

  struct CerPMT { int paddle;  double adc, tdc; }; 
  
  ROOT::RDF::RNode df = ROOT::RDataFrame("T", path_infile);
  
  
  
  df = df
    .Define("paddles", [](const RVec<double>& adc, const RVec<double>& tdc)
  {
    RVec<CerPMT> paddles; paddles.reserve(n_paddles);

    for (int i=0; i<n_paddles; i++) {

      paddles.push_back({ i, adc.at(i), tdc.at(i) }); 
    }
    return paddles; 
  }, {"L.cer.a_c", "L.cer.t_c"})
    
    .Define("good_hits", [](const RVec<CerPMT>& paddles) {
      RVec<double> adcs;
      for (const auto& paddle : paddles) {
	if (log(fabs(paddle.tdc)) < 0) {
	  adcs.push_back(paddle.adc);
	}
      } return adcs;
    }, {"paddles"})
    
    .Define("bad_hits", [](const RVec<CerPMT>& paddles) {
      RVec<double> adcs;
      for (const auto& paddle : paddles) {
	if (log(fabs(paddle.tdc)) > 0) {
	  adcs.push_back(paddle.adc);
	}
      } return adcs;
    }, {"paddles"});
  
  
  const double min_adc = -1000;
  const double max_adc = +15000; 
  
  auto hist_good = (TH1D*)df.Histo1D<RVec<double>>({
      "h_good", "Cerenkov ADC values;ADC (arb. units)",
      200, min_adc, max_adc},
    "good_hits" )->Clone(); 
      
  auto hist_bad = (TH1D*)df.Histo1D<RVec<double>>({
      "h_bad", "Cerenkov ADC values;ADC (arb. units)",
      200, min_adc, max_adc},
    "bad_hits" )->Clone(); 

  new TCanvas;
  gStyle->SetOptStat(0);
  gPad->SetLogy(1);
  
  auto legend = new TLegend;
  
  
  hist_bad->SetLineColor(kRed); 
  hist_bad->SetFillColor(kRed);
  hist_bad->SetFillStyle(3004);
  hist_bad->Draw();
  legend->AddEntry(hist_bad,  "null");  

  hist_good->SetLineWidth(2); 
  hist_good->SetLineColor(kBlack);
  hist_good->Draw("SAME");
  legend->AddEntry(hist_good, "valid");  
  
  legend->SetHeader("Paddle TDC time");
  legend->Draw(); 
  return 0;
}
