#include <APEX/units.h>
// ROOT 
#include <TH2D.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <APEX/Cerenkov.h>
#include <TString.h> 
#include <TCanvas.h> 
#include <TBox.h> 
// stdlib
#include <string> 
#include <cmath> 
#include <utility> 


class PdfWriter {
private: 
    std::map<std::string, TCanvas*> pdfs{};  

public: 
    PdfWriter() = default; 

    void write_pdf_page(std::string path, TCanvas* canv) {
     
        canv->cd(); 

        auto it = pdfs.find(path); 
        // open new PDF
        if (it==pdfs.end()) {

            pdfs[path] = canv;
            canv->Print(std::string{path + "["}.c_str(), "pdf"); 

        } else {
            
            canv->Print(path.c_str(), "pdf"); 
        }
    }

    void ClosePdfs() { 
        for (const auto& [path, canv]: pdfs) {
            canv->Print(std::string{path + "]"}.c_str());
        }
    }
}; 

std::pair<double,double> get_paddle_pos(bool is_RHRS, int i)
{
    double x_offset = is_RHRS ? (-106.0*APEX::units::cm) : (-103.7*APEX::units::cm); 
    double x_pace   = is_RHRS ? (+42.5*APEX::units::cm)  : (+43.9*APEX::units::cm); 

    double x = x_offset + x_pace*((double)(i/2));   
    double y = ((i%2)==0) ? +0.150 : -0.150; 
    return {x,y}; 
}; 

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
    if ( 0.5*std::fabs(x - x_mirror) > 0.5*mirror_width_X(is_RHRS) ) return false; 
    if ( 0.5*std::fabs(y - y_mirror) > 0.5*mirror_width_Y(is_RHRS) ) return false; 

    return true; 
}

void draw_all_cerenkov_paddles(std::string path_infile, std::string path_outgraphic="plots/cerenkov_paddles_associated")
{
    ROOT::EnableImplicitMT(); 

    ROOT::RDataFrame df("track_data", path_infile); 

    auto canv = new TCanvas; 
    PdfWriter pdf_writer; 

    for (const bool is_RHRS : {true, false}) {

        std::string arm = is_RHRS ? "R" : "L"; 

        for (int i=0; i<APEX::Cerenkov::N_pmts; i++) {

            auto hist = df  

                .Filter([i](ROOT::RVec<int> pnum){ 
                    return ROOT::VecOps::Filter(pnum, [i](int p){ return p==i; }).empty() == false; 
                }, {arm+"_cer_paddle"})

                .Define(arm+"_Cerenkov_xy", [i,is_RHRS](ROOT::RVec<double> x_trk, ROOT::RVec<double> y_trk){
                    ROOT::RVec<std::pair<double,double>> pairs; pairs.reserve(x_trk.size()); 
                    for (int t=0; t<x_trk.size(); t++) {
                        if (is_track_associated(is_RHRS, i, x_trk[t],y_trk[t])) 
                            pairs.emplace_back( x_trk[t], y_trk[t] );
                    } 
                    return pairs; 
                }, {arm+"_Cerenkov_x", arm+"_Cerenkov_y"})

                .Redefine(arm+"_Cerenkov_x", [](ROOT::RVec<std::pair<double,double>> v){
                    ROOT::RVec<double> t; t.reserve(v.size());
                    for (const auto& [x,y] : v) t.emplace_back(x); 
                    return t; 
                }, {arm+"_Cerenkov_xy"})

                .Redefine(arm+"_Cerenkov_y", [](ROOT::RVec<std::pair<double,double>> v){
                    ROOT::RVec<double> t; t.reserve(v.size());
                    for (const auto& [x,y] : v) t.emplace_back(y); 
                    return t; 
                }, {arm+"_Cerenkov_xy"})

                .Histo2D({
                    "h", 
                    Form("%sHRS, pmt %i - All tracks projected onto cerenkov plane;x_{Cerenkov} (m);y_{Cerenkov} (m)",arm.c_str(),i),
                    200, -1.1, 0.9, 
                    200, -0.16, 0.20},
                    arm+"_Cerenkov_x", arm+"_Cerenkov_y"
                );

            canv->Clear(); 

            auto [pad_x, pad_y] = get_paddle_pos(is_RHRS, i); 

    
            hist->DrawCopy("col"); 

            const double paddle_x_width = mirror_width_X(is_RHRS); 
            const double paddle_y_width = mirror_width_Y(is_RHRS); 

            auto box = new TBox(
                pad_x - paddle_x_width/2., pad_y - paddle_y_width/2., 
                pad_x + paddle_x_width/2., pad_y + paddle_y_width/2.
            ); 
            box->SetFillStyle(0);
            box->SetLineColor(kRed); 
            box->Draw();

            canv->Modified();
            canv->Update(); 

            pdf_writer.write_pdf_page(path_outgraphic+"_"+arm+"HRS.pdf", canv); 
        }
    }

    pdf_writer.ClosePdfs(); 
}