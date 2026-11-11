
// APEX
#include <APEX/S2.h>
#include <APEX/units.h> 
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

void fit_S2_pmt_dt(
    std::string infile_path="data-v2/replay.4016.root", 
    std::string outgraphic_path="plots/S2_pmt_dt",
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

    // effective scintillator photon speed, in TDC tick / cm
    constexpr double v_effective = (5./3.)*units::cm/units::ns; 

    TCanvas* canv_all_paddles = new TCanvas; 
    canv_all_paddles->Divide(2,1); 


    TCanvas* canv_error = new TCanvas; 
    canv_error->Divide(2,1); 

    // open the output file 
    std::fstream outfile(outfile_path, std::ios::out | std::ios::trunc); 
        
    // do right-arm first
    for (const bool is_RHRS : {true, false}) {


        std::string arm = is_RHRS ? "R" : "L"; 

        double pmt_dt_range[2] = {
            is_RHRS ? -30. : -135.,
            is_RHRS ? +15. : -95.
        }; 
        
        std::vector<double> pmt_offsets(S2::N_paddles, 0.); 
        
        for (int paddle=0; paddle<S2::N_paddles; paddle++) {

            auto df_out = df    
                .Define("pmt_data", [paddle, v_effective](
                    Rv<int> t_paddle, 
                    Rv<double> t_right_pmt_tdc, 
                    Rv<double> t_left_pmt_tdc, 
                    Rv<double> t_y){
                    Rv<S2_pmt_data> data; data.reserve(t_paddle.size());
                    
                    for (std::size_t i=0; i<t_paddle.size(); i++) {
                        if (t_paddle[i] == paddle) {
                            double pmt_dt = t_right_pmt_tdc[i] - t_left_pmt_tdc[i];
                            data.emplace_back(
                                ((pmt_dt)*APEX::S2::TDC_resolution*v_effective - t_y[i])/units::cm, t_y[i]/units::cm
                            ); 
                        }
                    }
                    return data; 

                }, {arm+"_S2_paddle", arm+"_S2_right_pmt_tdc", arm+"_S2_left_pmt_tdc", arm+"_S2_y"})

                .Filter([](Rv<S2_pmt_data> data){ return !data.empty(); }, {"pmt_data"})

                .Define("S2_pmt_dt_raw", [](const Rv<S2_pmt_data>& v){ 
                    Rv<double> out; 
                    out.reserve(v.size()); 
                    for(const auto& d : v) out.emplace_back(d.pmt_dt); 
                    return out; 
                }, {"pmt_data"}); 

            auto hist_p = df_out
                .Histo1D<Rv<double>>({
                    "h_pmt_dt", 
                    Form("%s S2, Paddle %i;TDC PMT_{R} - TDC PMT_{L} - 2/3*y_{S2};", (is_RHRS ? "RHRS" : "LHRS"), paddle), 
                    200, pmt_dt_range[0], pmt_dt_range[1]}, 
                    "S2_pmt_dt_raw"
                ); 
            
            canv->cd(); 
            canv->Clear(); 

            hist_p->DrawCopy(); 
            
            auto xax = hist_p->GetXaxis(); 
            
            const double radius = 7.5; 
            int max_bin = hist_p->GetMaximumBin(); 
            double center = xax->GetBinCenter(max_bin); 
            double max    = hist_p->GetMaximum(); 

            f_gaus->SetRange(center-radius, center+radius); 

            f_gaus->SetParameter(0, max);
            f_gaus->SetParameter(1, center); 
            f_gaus->SetParameter(2, radius/2.); 

            auto fit_result = hist_p->Fit("f_gaus", "R N S"); 

            if (static_cast<int>(fit_result)==0 && fit_result->IsValid()) {
                
                double offset = fit_result->Parameter(1);
                auto line = new TLine(offset,0., offset,hist_p->GetMaximum()); 
                line->SetLineStyle(kDashed); 
                line->SetLineColor(kRed); 
                line->Draw();  

                pmt_offsets[paddle] = offset; 

                double sigma  = fit_result->Parameter(2); 

                std::printf("offset: %.1f, sigma: %.1f\n", offset, sigma); 
                f_gaus->DrawClone("SAME"); 
            }

            canv->Modified(); 
            canv->Update(); 

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
                    case (S2::N_paddles) : 
                        canv->Print(std::string{path + "]"}.c_str(), "pdf"); break; 

                    //any other page
                    default : 
                        canv->Print(path.c_str(), "pdf"); 
                }
            }; 
            save_canv_as_pdf(outgraphic_path + "_"+arm+"HRS"+ ".pdf"); 
        }

        //make the output
        outfile << 
            "\n\n"
            "std::array<double> pmt_dt_offsets_"+arm+"{\n";
        for (int i=0; i<S2::N_paddles; i++) {
            char buff[100]; std::sprintf(buff, "    %.2f%s\n", pmt_offsets[i], (i<15?",":""));  
            outfile << buff;
        }      
        outfile << "}\n"; 


        auto df_all = df
            
            .Define("pmt_data_corrected", [&pmt_offsets](Rv<int> t_paddle, Rv<double> t_pmt_dt, Rv<double> t_y){
                Rv<S2_pmt_data> data; data.reserve(t_y.size()); 
                for (std::size_t i=0; i<t_y.size(); i++) {
                    data.emplace_back(
                        t_pmt_dt[i]*units::ns - pmt_offsets[t_paddle[i]], 
                        t_y[i]
                    );
                }
                return data; 
            }, {arm+"_S2_paddle", arm+"_S2_pmt_dt", arm+"_S2_y"})

            .Define("S2_y", [](const Rv<S2_pmt_data>& v){ 
                Rv<double> out; 
                out.reserve(v.size()); 
                for(const auto& d : v) out.emplace_back(d.S2_y/units::cm); 
                return out; 
            }, {"pmt_data_corrected"})

            .Define("S2_pmt_dt_raw_corrected", [](const Rv<S2_pmt_data>& v){ 
                Rv<double> out; 
                out.reserve(v.size()); 
                for(const auto& d : v) out.emplace_back(d.pmt_dt); 
                return out; 
            }, {"pmt_data_corrected"})

           .Define("S2_y_error", [v_effective](const Rv<S2_pmt_data>& v){ 
                Rv<double> out; 
                out.reserve(v.size()); 
                for(const auto& d : v) out.emplace_back(d.S2_y/units::cm - d.pmt_dt/v_effective); 
                return out; 
            }, {"pmt_data_corrected"});
            
            
        auto histo2d = df_all
            
            .Histo2D<Rv<double>>({
                "h_y_dt_corr", 
                Form("%s: Corrected PMT-DT times; PMT_{R} TDC - PMT_{L} TDC;y_{S2} (cm)", (is_RHRS?"RHRS":"LHRS")), 
                100, -30., +30.,
                150, -30., +30.
            }, "S2_pmt_dt_raw_corrected", "S2_y"); 
        

        auto histo1d = df_all

            .Histo1D<Rv<double>>({
                "h_y_dt_corr", 
                Form("%s: Corrected PMT-DT times;y_{S2} - y(#Delta T) (cm)", (is_RHRS?"RHRS":"LHRS")), 
                100, -16.,+16.
            }, "S2_y_error"); 

        
        canv_all_paddles->cd((is_RHRS ? 1 : 2)); 
        histo2d->DrawCopy("col"); 
        
        auto line = new TLine(-30.*v_effective,-30, +30*v_effective,30); 
        line->SetLineColor(kRed);
        line->SetLineStyle(kDashed); 
        line->Draw(); 

        canv_error->cd((is_RHRS ? 1 : 2));
        histo1d->DrawCopy("col"); 

        f_gaus->SetRange(-16, +16); 

        f_gaus->SetParameter(0, histo1d->GetMaximum());
        f_gaus->SetParameter(1, 0.); 
        f_gaus->SetParameter(2, 2.);
        
        auto fit_result = histo1d->Fit("f_gaus", "R N S");
        
        if (static_cast<int>(fit_result)==0 && fit_result->IsValid()) {
        
            std::printf("sigma %s: %.2f cm\n", (is_RHRS?"RHRS":"LHRS"), std::fabs(f_gaus->GetParameter(2))); 
            
            f_gaus->DrawClone("SAME"); 
        }
        
            
    }

    outfile.close(); 
}