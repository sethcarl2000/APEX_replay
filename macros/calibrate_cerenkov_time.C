
// APEX
#include <APEX/S2.h>
#include <APEX/Cerenkov.h>
#include <APEX/units.h> 
#include <APEX/utils.h>
// ROOT
#include <TF1.h> 
#include <TH1D.h> 
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFitResult.h> 
#include <TFitResultPtr.h>
#include <TCanvas.h> 
#include <TString.h> 
#include <TGraph.h> 
#include <TAxis.h> 
#include <TLine.h> 
// stdlib
#include <string> 
#include <vector> 
#include <map> 
#include <memory> 
#include <fstream> 
#include <algorithm> 

template<typename T> using Rv = ROOT::VecOps::RVec<T>; 

void calibrate_cerenkov_time(
    std::string infile_path="data-v2/replay.4016.root", 
    std::string outgraphic_path="plots/Cerenkov_pmt_offsets",
    std::string outfile_path="data/csv/pmt_dt_offsets.txt"
)
{
    ROOT::EnableImplicitMT(); 

    ROOT::RDataFrame df("track_data", infile_path);

    using namespace APEX; 

    auto canv = new TCanvas; 

    auto f_gaus = std::make_unique<TF1>("f_gaus", "gaus(0)", -1., +1.); 
    f_gaus->SetLineStyle(kDashed); 
        
    struct S2_pmt_data { double pmt_dt, S2_y; }; 

    TCanvas* canv_all_pmts = new TCanvas; 
    canv_all_pmts->Divide(2,1); 


    constexpr double min_cerenkov_adc = 1000.; 

    TCanvas* canv_all_tdc_adc = new TCanvas; 
    canv_all_tdc_adc->Divide(2,1); 

    TCanvas* canv_error = new TCanvas; 
    canv_error->Divide(2,1); 

    // open the output file 
    std::fstream outfile(outfile_path, std::ios::out | std::ios::trunc); 
        
    // do right-arm first
    for (const bool is_RHRS : {true, false}) {

        std::string arm = is_RHRS ? "R" : "L"; 

        const double pmt_err_min = is_RHRS ? -1280*units::ns : -1400*units::ns; 
        const double pmt_err_max = is_RHRS ? -1240*units::ns : -1370*units::ns; 

        std::vector<double> pmt_offsets(Cerenkov::N_pmts, 0.); 
        
        for (int i_pmt=0; i_pmt<Cerenkov::N_pmts; i_pmt++) {

            auto df_out = df    

                .Define("ToF_cer_S2", [](Rv<double> t_S2, Rv<double> t_cer){ 
                    Rv<double> ToF; ToF.reserve(t_S2.size()); 
                    for (int i=0; i<t_S2.size(); i++) ToF.emplace_back(t_S2[i] - t_cer[i]);
                    return ToF; 
                }, {arm+"_S2_dt", arm+"_Cerenkov_t"})

            // fix typo in replay file
            .Define("t_cer_adc", [i_pmt, min_cerenkov_adc](
                Rv<int> cer_pmt,
                Rv<double> cer_ADC, 
                Rv<double> ToF
            ){
                    
                Rv<double> ADC; 
                for (std::size_t i=0; i<cer_pmt.size(); i++) { 
                    int pmt = cer_pmt[i];
                    if (pmt != i_pmt) continue; 
                    if (cer_ADC[i] < min_cerenkov_adc) continue; 

                    for (std::size_t j=0; j<ToF.size(); j++) {
                        ADC.emplace_back( cer_ADC[i] ); 
                    }
                }
                return ADC; 

            }, {arm+"_cer_paddle", arm+"_cer_adc", "ToF_cer_S2"})

                // fix typo in replay file
                .Define("t_cer_minus_t_S2", [i_pmt, min_cerenkov_adc](
                    Rv<int> cer_pmt,
                    Rv<double> cer_TDC, 
                    Rv<double> cer_ADC,
                    Rv<double> S2_TDC, 
                    Rv<double> ToF
                ){        
                    Rv<double> t_S2_minus_t_cer;
                    for (std::size_t i=0; i<cer_pmt.size(); i++) { 
                        int pmt = cer_pmt[i];
                        if (pmt != i_pmt) continue; 
                        if (cer_ADC[i] < min_cerenkov_adc) continue; 

                        for (std::size_t j=0; j<ToF.size(); j++) {
                            t_S2_minus_t_cer.emplace_back(
                                ((cer_TDC[i]+ToF[j]) - S2_TDC[j])/units::ns
                            ); 
                        }
                    }
                    return t_S2_minus_t_cer; 

                }, {arm+"_cer_paddle", arm+"_cer_time",arm+"_cer_adc", arm+"_S2_hit_time", "ToF_cer_S2"}); 

            auto hist_tdc_adc = df_out
                .Histo2D<Rv<double>>({
                    "h_pmt_dt", 
                    Form("%sHRS Cerenkov, PMT %i ADC vs TDC;Paddle (T_{cer} + ToF_{cer->S2}) - T_{S2} (ns);Paddle ADC", (is_RHRS ? "R" : "L"), i_pmt), 
                    75, pmt_err_min/units::ns, pmt_err_max/units::ns,
                    75, 0.,   1e4
                }, "t_cer_minus_t_S2", "t_cer_adc"); 

            auto hist_p = df_out
                .Histo1D<Rv<double>>({
                    "h_pmt_dt", 
                    Form("%s Cerenkov, PMT %i ADC vs TDC;Paddle (T_{cer} + ToF_{cer->S2}) - T_{S2} (ns);", (is_RHRS ? "RHRS" : "LHRS"), i_pmt), 
                    75, pmt_err_min/units::ns, pmt_err_max/units::ns}, 
                    "t_cer_minus_t_S2"
                ); 
            
            //save each 'file' to a pdf
            auto save_canv_as_pdf = [canv](std::string path){
                static std::map<std::string, std::size_t> pages{};  
                
                int page; 
                auto find_pdf = pages.find(path);
                if (find_pdf == pages.end()) { //no pages added for this pdf (yet), so make a new one 
                    pages[path] = 1 ;
                    page = 1; 
                } else { // this pdf already has pages, so we just add to the page count
                    auto& this_page = find_pdf->second; 
                    page = (++this_page); 
                }

                switch (page) {

                    //first page 
                    case 1 : 
                        canv->Print(std::string{path + "["}.c_str(), "pdf"); break; 
                    
                    //last page
                    case (Cerenkov::N_pmts*2) : 
                        canv->Print(std::string{path + "]"}.c_str(), "pdf"); break; 

                    //any other page
                    default : 
                        canv->Print(path.c_str(), "pdf"); 
                }
            }; 

            canv->cd(); 
            canv->Clear(); 

            hist_p->DrawCopy(); 
            
            auto xax = hist_p->GetXaxis(); 
            
            const double radius = 4.; 
            int max_bin = hist_p->GetMaximumBin(); 
            double center = xax->GetBinCenter(max_bin); 
            double max    = hist_p->GetMaximum(); 

            f_gaus->SetRange(center-radius, center+radius); 

            f_gaus->SetParameter(0, max);
            f_gaus->SetParameter(1, center); 
            f_gaus->SetParameter(2, radius/2.); 

            auto fit_result = hist_p->Fit("f_gaus", "R N Q S"); 

            if (static_cast<int>(fit_result)==0 && fit_result->IsValid()) {
                
                double offset = fit_result->Parameter(1);
                auto line = new TLine(offset,0., offset,hist_p->GetMaximum()); 
                line->SetLineStyle(kDashed); 
                line->SetLineColor(kRed); 
                line->Draw();  

                pmt_offsets[i_pmt] = offset; 

                double sigma  = fit_result->Parameter(2); 

                std::printf("%sHRS paddle: %i offset: %.1f, sigma: %.1f\n", arm.c_str(), i_pmt, offset, sigma); 
                f_gaus->DrawClone("SAME"); 
            }

            canv->Modified(); 
            canv->Update(); 

            if (!outgraphic_path.empty())
                save_canv_as_pdf(outgraphic_path + "_"+arm+"HRS"+ ".pdf"); 

            
            canv->cd(); 
            canv->Clear(); 

            hist_tdc_adc->DrawCopy("col"); 

            canv->Modified(); 
            canv->Update(); 
            
            if (!outgraphic_path.empty())
                save_canv_as_pdf(outgraphic_path + "_"+arm+"HRS"+ ".pdf"); 

        }
        
        auto df_all = df

            .Define("ToF_cer_S2", [](Rv<double> t_S2, Rv<double> t_cer){ 
                Rv<double> ToF; ToF.reserve(t_S2.size()); 
                for (int i=0; i<t_S2.size(); i++) ToF.emplace_back(t_S2[i] - t_cer[i]);
                return ToF; 
            }, {arm+"_S2_dt", arm+"_Cerenkov_t"})

            // fix typo in replay file
            .Define("t_cer_adc", [&pmt_offsets](
                Rv<int> cer_pmt,
                Rv<double> cer_ADC, 
                Rv<double> ToF
            ){
                    
                Rv<double> ADC; 
                for (std::size_t i=0; i<cer_pmt.size(); i++) { 
                    int pmt = cer_pmt[i];
                    double offset = pmt_offsets[pmt];

                    for (std::size_t j=0; j<ToF.size(); j++) {
                        ADC.emplace_back( cer_ADC[i] ); 
                    }
                }
                return ADC; 

            }, {arm+"_cer_paddle", arm+"_cer_adc", "ToF_cer_S2"})

            // fix typo in replay file
            .Define("t_cer_minus_t_S2", [&pmt_offsets](
                Rv<int> cer_pmt,
                Rv<double> cer_TDC, 
                Rv<double> S2_TDC, 
                Rv<double> ToF){
                    
                Rv<double> t_S2_minus_t_cer;
                for (std::size_t i=0; i<cer_pmt.size(); i++) { 
                    int pmt = cer_pmt[i];
                    double offset = pmt_offsets[pmt];

                    for (std::size_t j=0; j<ToF.size(); j++) {
                        t_S2_minus_t_cer.emplace_back(
                            ((cer_TDC[i]+ToF[j]) - S2_TDC[j])/units::ns - offset
                        ); 
                    }
                }
                return t_S2_minus_t_cer; 

            }, {arm+"_cer_paddle", arm+"_cer_time", arm+"_S2_hit_time", "ToF_cer_S2"}); 

        auto hist_p_all = df_all
            .Histo1D<Rv<double>>({
                "h_pmt_dt", 
                Form("%sHRS Cerenkov, T_{cer} - T_{S2} (all PMTs, ToF corrected);(T_{cer} + ToF_{cer->S2}) - T_{S2} (ns);", (is_RHRS ? "RHRS" : "LHRS")), 
                100, -20., +20.}, 
                "t_cer_minus_t_S2"
            ); 

        auto hist_tdc_adc = df_all
            .Histo2D<Rv<double>>({
                "h_pmt_dt", 
                Form("%sHRS Cerenkov, T_{cer} - T_{S2} (all PMTs, ToF corrected) vs paddle ADC;Paddle (T_{cer} + ToF_{cer->S2}) - T_{S2} (ns);Paddle ADC", (is_RHRS ? "R" : "L")), 
                50, -20., +20.,
                50, 0., 1e4
            }, "t_cer_minus_t_S2", "t_cer_adc"); 

        canv_all_pmts->cd(is_RHRS ? 1 : 2); 

        auto h = (TH1D*)hist_p_all->Clone(std::string{arm+"my_cpy"}.c_str());
        h->Draw();

        double center, sigma; 
        utils::fit_gaus_to_hist(h, 5., center, sigma, false); 
        std::printf("%sHRS offset: %.3f, sigma: %.3f\n", arm.c_str(), center, sigma); 
    
        canv_all_tdc_adc->cd(is_RHRS ? 1 : 2); 
        hist_tdc_adc->DrawCopy("lego"); 
    }
}