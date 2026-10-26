#include <APEX/VDC.h>
#include <APEX/VDC/Hit.h>
#include <APEX/replay/helpers.h>
#include <RDFNodeAccumulator.h>
// ROOT
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TH2D.h> 
#include <ROOT/RVec.hxx>
#include <TCanvas.h> 
#include <THStack.h> 
#include <TLegend.h> 
#include <TStyle.h> 
// stdlib
#include <string> 
#include <algorithm> 

namespace {
    //Set units 
    constexpr double ns = 1e-9; 
}

template<typename T> using Rv = ROOT::VecOps::RVec<T>; 

void construct_all_VDC_hits(std::string path_infile)
{
    ROOT::EnableImplicitMT(); 
    ROOT::RDataFrame df("T", path_infile); 

    ROOT::RDF::Experimental::AddProgressBar(df); 

    RDFNodeAccumulator rna(df); 

    const bool kRHRS_ = true; 
    const bool kLHRS_ = false; 

    struct TimeRL { double tR, tL; }; 

    std::string plane_names[APEX::VDC::kN_planes] = {"u1", "v1", "u2", "v2"};  

    for (bool is_RHRS : {kRHRS_, kLHRS_}) {
        
        std::string arm = is_RHRS ? "R" : "L"; 

        rna.Define(arm+"_S2_hits", [is_RHRS](const Rv<double>& pmt_r, const Rv<double>& pmt_l)
        {
            Rv<APEX::S2Hit> hits; hits.reserve(APEX::S2Hit::N_paddles());  

            for (int i=0; i<APEX::S2Hit::N_paddles(); i++) {
                auto hit = APEX::S2Hit(is_RHRS, i, pmt_l[i], pmt_r[i]); 
                if (hit.IsCoinc()) hits.emplace_back(hit); 
            }
            return hits; 
        }, {arm+".s2.rt", arm+".s2.lt"}); 

        for (int p=0; p<APEX::VDC::kN_planes; p++) {

            std::string plane = plane_names[p]; 

            auto arm_mode = is_RHRS ? APEX::replay::kRHRS : APEX::replay::kLHRS;  

            rna.Define("VDC_hits_"+arm+"_"+plane_names[p], 
                [arm_mode,p](const Rv<double>& wire, const Rv<double>& rawtime){

                    return APEX::replay::helpers::generate_VDC_hits(arm_mode,p,wire,rawtime); 

                }, {arm+".vdc."+plane_names[p]+".wire", arm+".vdc."+plane_names[p]+".rawtime"}
            ); 
        }
        rna.Define(arm+"_VDC_hits",
            [arm](
                const Rv<APEX::VDC::Hit>& u1, 
                const Rv<APEX::VDC::Hit>& v1, 
                const Rv<APEX::VDC::Hit>& u2, 
                const Rv<APEX::VDC::Hit>& v2){
                    
                //get total number of hits
                std::size_t all{0}; 
                for (const auto& pts : {u1,v1,u2,v2}) all += pts.size(); 

                Rv<APEX::VDC::Hit> all_hits; all_hits.reserve(all); 
                for (const auto& pts : {u1,v1,u2,v2})
                    all_hits.insert(all_hits.end(), pts.cbegin(), pts.cend()); 

                return all_hits;

            }, {"VDC_hits_"+arm+"_u1", "VDC_hits_"+arm+"_v1", "VDC_hits_"+arm+"_u2", "VDC_hits_"+arm+"_v2"}
        );

        rna.Define(arm+"_VDC_hits_realtime",
            [arm](const Rv<APEX::VDC::Hit>& hits){
                Rv<double> realtime; realtime.reserve(hits.size()); 
                for (const auto& hit : hits) realtime.emplace_back(hit.Time()); 
                return realtime; 

            }, {arm+"_VDC_hits"}
        );

        rna.Define(arm+"_S2rel_VDC_hits_realtime",
            [arm](const Rv<APEX::VDC::Hit>& vdc_hits, const Rv<APEX::S2Hit>& s2_hits){
                Rv<double> realtime; realtime.reserve(vdc_hits.size() * s2_hits.size()); 
                for (const auto& vdc_hit : vdc_hits) 
                    for (const auto& s2_hit : s2_hits)
                        realtime.emplace_back(vdc_hit.Time() - s2_hit.Time()); 

                return realtime; 

            }, {arm+"_VDC_hits", arm+"_S2_hits"}
        );
    
        rna.Define(arm+"_vdc_nhits", [](const Rv<APEX::VDC::Hit>& hits){
            return static_cast<double>(hits.size()); 
        }, {arm+"_VDC_hits"}); 
        
        rna.Define(arm+"_S2_nhits", [](const Rv<APEX::S2Hit>& hits){
            return static_cast<double>(hits.size()); 
        }, {arm+"_S2_hits"}); 
        
        rna.Define(arm+"_S2_paddle", [](const Rv<APEX::S2Hit>& hits){
            Rv<double> paddles; paddles.reserve(hits.size()); 
            for (const auto& hit : hits) paddles.push_back(static_cast<double>(hit.Paddle())); 
            return paddles; 
        }, {arm+"_S2_hits"}); 
        
        rna.Define(arm+"_S2_PMTR_minus_PMTL", [](const Rv<APEX::S2Hit>& hits){
            Rv<double> ret; ret.reserve(hits.size()); 
            for (const auto& hit : hits) ret.push_back(hit.Time_PMT_R() - hit.Time_PMT_L()); 
            return ret; 
        }, {arm+"_S2_hits"}); 
    }

    auto h_all_R = rna.Get()
        .Histo1D({"h_all_R", "Timing of all VDC hits - RHRS;;", 
            70, -100e-9, +500e-9},
            "R_VDC_hits_realtime"
        ); 
    auto h_all_L = rna.Get()
        .Histo1D({"h_all_L", "Timing of all VDC hits - LHRS;;", 
            70, -100e-9, +500e-9},
            "L_VDC_hits_realtime"
        ); 

    auto h_all_S2rel_R = rna.Get()
        .Histo1D({"h_all_S2rel_R", "Timing of all VDC hits, rel to S2 hits - RHRS;;", 
            70, -100e-9, +500e-9},
            "R_S2rel_VDC_hits_realtime"
        ); 
    auto h_all_S2rel_L = rna.Get()
        .Histo1D({"h_all_S2rel_L", "Timing of all VDC hits, rel to S2 hits - LHRS;;", 
            70, -100e-9, +500e-9},
            "L_S2rel_VDC_hits_realtime"
        ); 
    

    auto h_R_S2_reltime = rna.Get()
        .Histo2D({"h_R_S2_reltime", "rel. PMT times for S2 paddles - RHRS;paddle ind.;PMT-R - PMT-L (s)",
            16, -0.5, 16-0.5, 
            41, -10e-9, +10e-9
        }, "R_S2_paddle", "R_S2_PMTR_minus_PMTL"); 
    auto h_L_S2_reltime = rna.Get()
        .Histo2D({"h_L_S2_reltime", "rel. PMT times for S2 paddles - LHRS;paddle ind.;PMT-R - PMT-L (s)",
            16, -0.5, 16-0.5, 
            41, -10e-9, +10e-9
        }, "L_S2_paddle", "L_S2_PMTR_minus_PMTL"); 



    auto h_R_S2_nhits = rna.Get()
        .Histo1D({"h_R_S2_nhits", "Number of Coinc S2 hits - RHRS;n. R+L PMT coinc hits / event;", 
            17, -0.5, 16+0.5},
            "R_S2_nhits"
        ); 
    auto h_L_S2_nhits = rna.Get()
        .Histo1D({"h_L_S2_nhits", "Number of Coinc S2 hits - LHRS;n. R+L PMT coinc hits / event;", 
            17, -0.5, 16+0.5},
            "L_S2_nhits"
        ); 


    TCanvas* canv; 
    
    canv = new TCanvas; 
    h_all_R->DrawCopy(); 
    canv = new TCanvas; 
    h_all_L->DrawCopy(); 


    canv = new TCanvas; 
    h_all_S2rel_R->DrawCopy(); 
    canv = new TCanvas; 
    h_all_S2rel_L->DrawCopy(); 


    canv = new TCanvas; 
    h_R_S2_nhits->DrawCopy(); 
    canv = new TCanvas; 
    h_L_S2_nhits->DrawCopy(); 

    canv = new TCanvas; 
    h_R_S2_reltime->DrawCopy("colz"); 
    canv = new TCanvas; 
    h_L_S2_reltime->DrawCopy("colz"); 
    
}