
// APEX
#include "gaus_integrate.h"
#include <APEX/S2.h>
#include <APEX/Cerenkov.h>
#include <APEX/units.h> 
#include <APEX/utils.h>
// ROOT
#include <TF1.h> 
#include <TH1D.h> 
#include <Math/RootFinder.h>
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
#include <Math/ProbFunc.h>
// stdlib
#include <string> 
#include <limits> 
#include <vector> 
#include <map> 
#include <memory> 
#include <fstream> 
#include <algorithm> 


template<typename T> using Rv = ROOT::VecOps::RVec<T>; 

std::pair<double,double> get_paddle_pos(bool is_RHRS, int i)
{
    double x_offset = is_RHRS ? (-106.0*APEX::units::cm) : (-103.7*APEX::units::cm); 
    double x_pace   = is_RHRS ? (+42.5*APEX::units::cm)  : (+43.9*APEX::units::cm); 

    double x = x_offset + x_pace*((double)(i/2));   
    double y = ((i%2)==0) ? +0.150 : -0.150; 
    return {x,y}; 
}; 

/// @brief Skew-normal distribution
/// @param x absolute position 
/// @param A height 
/// @param x0 offset 
/// @param l length scale (sigma)
/// @param alpha skewness parameter. when alpha=0, this is a regular normal dist. 
double skewnormal_dist(double x, double A, double x0, double l, double alpha){
    x = (x - x0)/l; 
    return A * std::exp( -0.5*x*x ) * ROOT::Math::normal_cdf(alpha * x); 
};

double find_skewnormal_mode(double x0, double l, double alpha)
{
    // for simplicity, we do all these root finding calculations with x>0, as the skew-normal is 
    // invariant under reflections w/r/t alpha. i.e. 
    // f_skewnormal(x, alpha) = f_skewnormal(-x, -alpha); 
    double sign_alpha = alpha/std::fabs(alpha); 
    double abs_alpha = std::fabs(alpha); 
    // first, estimate the mode 
    // i got these formulas from: https://en.wikipedia.org/wiki/Skew_normal_distribution

    auto f_skewnorm_derivative = [abs_alpha](double x){ 
        double arg = x*abs_alpha; 
        return 
            abs_alpha * std::exp( -0.5*arg*arg )/(std::sqrt(2.*3.14159265)*ROOT::Math::normal_cdf(arg)) - x; 
    }; 
    ROOT::Math::Functor1D functor(f_skewnorm_derivative); 
    ROOT::Math::RootFinder rf(ROOT::Math::RootFinder::kBRENT); 
    rf.SetFunction(functor, -0.5, 1.5); 
    double mode; 
    if (rf.Solve(100, 1e-12, 1e-8)) {
        mode = rf.Root(); 
    } else {
        Error(__func__, "Something went wrong trying to find root of skew-normal dist."); 
        return std::numeric_limits<double>::quiet_NaN(); 
    }
    // scale & shift the mode we found 
    mode *= sign_alpha; 
    return (mode*l) + x0;   
} 

double mirror_width_X(bool is_RHRS) {
    return is_RHRS ? 68.5*APEX::units::cm : 63.5*APEX::units::cm; 
}
double mirror_width_Y(bool is_RHRS) {
    return is_RHRS ? 47.0*APEX::units::cm : 47.0*APEX::units::cm; 
}

/// @brief Determine if track is likely associated with a particular PMT. 
/// @param is_RHRS RHRS=true, LHRS=false
/// @param i_pmt pmt number [0,9]
/// @param x track x-intercept with Cerenkov mirror plane
/// @param y track y-intercept with Cerenkov mirror plane
/// @return 'true' if track's intercept is consistent with a given Cerenkov PMT mirror, false otherwise. 
bool is_track_associated(bool is_RHRS, int i_pmt, double x, double y)
{
    auto [x_mirror, y_mirror] = get_paddle_pos(is_RHRS, i_pmt); 

    //see if it agrees in the x-direction
    if ( std::fabs(x - x_mirror) > 0.5*mirror_width_X(is_RHRS) ) return false; 
    if ( std::fabs(y - y_mirror) > 0.5*mirror_width_Y(is_RHRS) ) return false; 

    return true; 
}

void calibrate_cerenkov_time(
    std::string infile_path="data-v2/replay-NEW.4016.root", 
    std::string outgraphic_path="plots/Cerenkov_pmt_offsets_TEST",
    std::string outfile_path="data/csv/cerenkov_pmt_offsets-TEST.txt"
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


    constexpr double min_cerenkov_adc = 0.; 

    TCanvas* canv_all_tdc_adc = new TCanvas; 
    canv_all_tdc_adc->Divide(2,1); 

    TCanvas* canv_error = new TCanvas; 
    canv_error->Divide(2,1); 

    struct Track {
        double x_cer, y_cer, ToF, S2_time;
    };

    
    const int n_fit_bins = 40; 
    const double fit_radius = 20.; 
    const double dx = 2.*fit_radius / ((double)n_fit_bins);
    auto skewnormal_integ = [dx](double *X, double *par){
        return gaus_integrate([par](double x){
            return skewnormal_dist(x, par[0], par[1], par[2], par[3]); 
        }, X[0] - dx/2., X[0] + dx/2.); 
    };
    auto f_skewnormal = new TF1("f_skewnormal", skewnormal_integ, 
        -fit_radius, 
        +fit_radius,
        4
    ); 

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

                .Define("associated_tracks", [is_RHRS,i_pmt](
                    Rv<double> x_trk,
                    Rv<double> y_trk,
                    Rv<double> S2_tdc, 
                    Rv<double> S2_dt,
                    Rv<double> t_cer
                ){
                    ROOT::RVec<Track> tracks; tracks.reserve(x_trk.size()); 
                    for (std::size_t t=0; t<x_trk.size(); t++) {
                        if (is_track_associated(is_RHRS,i_pmt, x_trk[t], y_trk[t])) {
                            tracks.emplace_back(
                                x_trk[t], y_trk[t], 
                                (S2_dt[t] - t_cer[t]), 
                                S2_tdc[t]
                            ); 
                        }
                    }
                    return tracks; 

                }, {arm+"_Cerenkov_x", arm+"_Cerenkov_y", arm+"_S2_hit_time", arm+"_S2_dt", arm+"_Cerenkov_t"})


                .Define("ToF_cer_S2", [](Rv<double> t_S2, Rv<double> t_cer){ 
                    Rv<double> ToF; ToF.reserve(t_S2.size()); 
                    for (int i=0; i<t_S2.size(); i++) ToF.emplace_back(t_S2[i] - t_cer[i]);
                    return ToF; 
                }, {arm+"_S2_dt", arm+"_Cerenkov_t"})

            // fix typo in replay file
            .Define("t_cer_adc", [i_pmt, min_cerenkov_adc](
                Rv<int> cer_pmt,
                Rv<double> cer_ADC, 
                Rv<Track> tracks
            ){
                    
                Rv<double> ADC; 
                for (std::size_t i=0; i<cer_pmt.size(); i++) { 
                    int pmt = cer_pmt[i];
                    if (pmt != i_pmt) continue; 
                    if (cer_ADC[i] < min_cerenkov_adc) continue; 

                    for (std::size_t j=0; j<tracks.size(); j++) {
                        ADC.emplace_back( cer_ADC[i] ); 
                    }
                }
                return ADC; 

            }, {arm+"_cer_paddle", arm+"_cer_adc", "associated_tracks"})

                // fix typo in replay file
                .Define("t_cer_minus_t_S2", [i_pmt, min_cerenkov_adc](
                    Rv<int> cer_pmt,
                    Rv<double> cer_TDC, 
                    Rv<double> cer_ADC,
                    Rv<Track> tracks
                ){        
                    Rv<double> t_S2_minus_t_cer;
                    for (std::size_t i=0; i<cer_pmt.size(); i++) { 
                        int pmt = cer_pmt[i];
                        if (pmt != i_pmt) continue; 
                        if (cer_ADC[i] < min_cerenkov_adc) continue; 

                        for (std::size_t j=0; j<tracks.size(); j++) {
                            t_S2_minus_t_cer.emplace_back(
                                ((cer_TDC[i]+tracks[j].ToF) - tracks[j].S2_time)/units::ns
                            ); 
                        }
                    }
                    return t_S2_minus_t_cer; 

                }, {arm+"_cer_paddle", arm+"_cer_time",arm+"_cer_adc", "associated_tracks"}); 

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
                n_fit_bins, -fit_radius, +fit_radius}, 
                "t_cer_minus_t_S2"
            ); 

        auto hist_tdc_adc = df_all
            .Histo2D<Rv<double>>({
                "h_pmt_dt", 
                Form("#splitline{%sHRS Cerenkov, T_{cer} - T_{S2}}{(all PMTs, ToF corrected) vs paddle ADC};Paddle (T_{cer} + ToF_{cer->S2}) - T_{S2} (ns);Paddle ADC", (is_RHRS ? "R" : "L")), 
                n_fit_bins, -fit_radius, +fit_radius,
                40, 0., 1e4
            }, "t_cer_minus_t_S2", "t_cer_adc"); 

        canv_all_pmts->cd(is_RHRS ? 1 : 2); 

        auto h = (TH1D*)hist_p_all->Clone(std::string{arm+"my_cpy"}.c_str());
        h->Draw();

        // fill the output file 
        outfile << "\n\n\n"


            "\\*________________________________________________________________________\n"
            "\n"
            "   Cerenkov PMT corrections - "+arm+"HRS\n"
            "\n"
            "   Basically to get a 'reconstructed' time (a maximum-likelihood estimate (MLE) of the\n"
            "   'true' time of intercept of the track with the PMT mirror),\n"
            "   for a PMT TDC time: $t_{TDC}$ and pmt ADC $adc$, we have: \n"
            "   read from pmt 'i_pmt', you want to add:\n"
            "\n"
            "       t_{MLE} = t_{TDC} - offset[i_pmt] - ADC_corr(adc)\n"
            "\n"
            "   Where:\n"
            "       offset[i_pmt] is read from the following table, per-PMT\n"
            "       ADC_corr(adc) is the value 't_mode', interpolated from the table below (see below for description of t_mode.)\n"

            "  ________________________________________________________________________\n"\
            "*\\\n\n"
            "\\\\ "+arm+"HRS Cerenkov PMT offsets (indexed [0-9]):\n"
            "\n"
            "{\n"; 
        for (int pmt=0; pmt<Cerenkov::N_pmts; pmt++) {
            outfile << Form("   %+7.3f*units::ns%s\n", 
                pmt_offsets[pmt], 
                (pmt<Cerenkov::N_pmts-1?",":"")
            ); 
        }
        outfile << "}\n"; 
    

        double center, sigma; 
        utils::fit_gaus_to_hist(h, 5., center, sigma, false); 
        std::printf("%sHRS offset: %.3f, sigma: %.3f\n", arm.c_str(), center, sigma); 
    
        auto c = canv_all_tdc_adc->cd(is_RHRS ? 1 : 2); 
        c->SetLeftMargin(0.15); 
        c->SetTopMargin(0.15); 
        hist_tdc_adc->DrawCopy("col"); 

        auto xax = hist_tdc_adc->GetXaxis(); 
        const double dx = xax->GetBinWidth(1);
        // integrate the skew-normal distribution for each bin  

        double x0, l, alpha; 

        outfile << 
            "\n\n\n\n"
            "\\*________________________________________________________________________\n"
            "\n"
            "   Skew-normal shift:\n"
            "   \n"
            "   The way this PDF is defined is in the following way:\n"
            "   \n"
            "       f(t | t0,l,a) =  N((t-t0)/l) * CDF(a * (t-t0)/l)\n"
            "   \n"
            "   Where N(.) is a normalized gaus. PDF, and CDF(t) is: int_{-\\infty}^{t} N(t') dt'\n"
            "   The following list is given as: {ADC, t_mode, RMS, t0, l, a}, where:\n"
            "       ADC     - Cerenkov PMT ADC amplitude for this fit\n"
            "       t_mode  - the t-point at which the PDF f(t) is maximal. coincides with t0 iff a==0\n"
            "       RMS     - the variance: <t^2> - <t>^2. Be careful, as this is _not_ centered on t_mode\n"
            "       t0      - the 'shift' of the distribution. when a=0, this is the gaus. mean\n"
            "       l       - the 'width' of the distribution. when a=0, this is the gaus. rms\n"
            "       a       - skewness parameter, this parametrizes how 'skew' the distribution is.\n"
            "                     when a=0, this is a regular gaus. dist\n"
            "  ________________________________________________________________________\n"
            "*\\\n\n"
            "{\n"; 


        auto yax = hist_tdc_adc->GetYaxis(); 
        auto slice_canv = new TCanvas; 
        for (int iy=1; iy<=yax->GetNbins(); iy++) {

            double adc = yax->GetBinCenter(iy); 

            if (adc < min_cerenkov_adc) continue; 

            auto slice = hist_tdc_adc->ProjectionX("slice", iy,iy ); 
            slice->SetDirectory(nullptr); 
            slice->SetBit(kMustCleanup); 
            slice->ResetBit(kCanDelete); 

            slice->SetTitle(Form("%sHRS, ADC=[%.0f,%.0f];Paddle (T_{cer} + ToF_{cer->S2}) - T_{S2} (ns)", 
                arm.c_str(), yax->GetBinLowEdge(iy), yax->GetBinUpEdge(iy)
            ));
            slice_canv->Clear(); 

            //slice->Draw(); 
            
            // set guesses for all parameters 

            // if it's the first pass, esimate valid parameters. 
            if (iy==1) {
                x0 = xax->GetBinCenter(slice->GetMaximumBin()); 
                l  = slice->GetRMS(); 
                alpha = 0.; 
            }

            f_skewnormal->SetParameter(0, slice->GetMaximum()*dx); //height 
            f_skewnormal->SetParameter(1, x0); //center 
            f_skewnormal->SetParameter(2, l); // width 
            f_skewnormal->SetParameter(3, alpha); //skewness 

            auto fit_result = slice->Fit("f_skewnormal", "R S");
            double skewnormal_mode; 
            if (static_cast<int>(fit_result)==0 && fit_result->IsValid()) {
            
                std::printf("sigma %s: %.2f cm\n", (is_RHRS?"RHRS":"LHRS"), std::fabs(f_skewnormal->GetParameter(2))); 
                f_skewnormal->DrawClone("SAME"); 

                x0 = fit_result->Parameter(1); 
                l  = fit_result->Parameter(2); 
                alpha = fit_result->Parameter(3); 

                skewnormal_mode = find_skewnormal_mode(x0, l, alpha); 
                if (skewnormal_mode != skewnormal_mode) {
                    Error(__func__, "Something went wrong trying to find the root"); 
                } else {
                    auto line = new TLine(skewnormal_mode,0., skewnormal_mode,slice->GetMaximum()); 
                    line->SetLineStyle(kDashed); 
                    line->SetLineColor(kRed); 
                    line->Draw(); 
                
                    std::printf("%sHRS, slice [%.0f,%.0f] skewnormal x0: %.3f, l: %.3f, alpha: %.3f mode: %.3f\n", 
                        arm.c_str(), yax->GetBinLowEdge(iy),yax->GetBinUpEdge(iy), x0, l, alpha, 
                        skewnormal_mode
                    ); 
                }

                // compute the variance 
                double d2 = alpha*alpha/(1. + alpha*alpha); 
                double RMS = std::sqrt( l*l*( 1. - 2.*d2/3.14159265 ) ); 

                outfile << Form("    {%5.0f, %+8.4f*units::ns, %+8.4f*units::ns, %+8.4f*units::ns, %+8.4f*units::ns, %+8.8e}%s\n",
                    adc, skewnormal_mode, RMS, x0, l, alpha, (iy==yax->GetNbins()?"":",")
                );
            }

            slice_canv->Modified(); 
            slice_canv->Update(); 


            slice_canv->Print(Form("plots/%sHRS_cerenkov_adc_slice.pdf%s", arm.c_str(), (iy==1?"[":""))); 
        }
        outfile << "}\n"; 
        slice_canv->Print(Form("plots/%sHRS_cerenkov_adc_slice.pdf]", arm.c_str()));     
    }
    outfile.close(); 
}