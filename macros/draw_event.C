// EventDisplayKernel
#include <EDK.hpp>
// APEX stuff
#include <APEX/units.h>
#include <APEX/replay.h>
#include <APEX/replay/helpers.h>
#include <APEX/run_parameters.h>
// Replay classes
#include <APEX/EventHandler.h>
#include <APEX/VDC.h>
#include <APEX/VDC/Hit.h>
#include <APEX/VDC/HitGroup.h>
// ROOT 
#include <TH1D.h> 
#include <TGraph.h>
#include <TVirtualPad.h> 
#include <TStyle.h> 
#include <TAxis.h> 
#include <TBox.h> 
#include <ROOT/RVec.hxx>
// stdlib
#include <string> 
#include <span> 
#include <vector> 



namespace {
    constexpr bool kRHRS{true}, kLHRS{false}; 
    // left arm 
    constexpr bool is_RHRS = kLHRS;  
}

/// draw VDC window 
void drawVDC(); 

void draw_event(std::string path_infile="data-v2/decode.4016.root")
{
    using namespace APEX; 

    auto& kernel = EDK::Instance(); 

    ROOT::RDataFrame df("T", path_infile); 

    auto df_out = df    
        .Define("event_number", [](ULong64_t e_num){ return static_cast<UInt_t>(e_num); }, {"fEvtHdr.fEvtNum"}); 

    kernel.SetDataFrame(df_out); 
    kernel.SetEventNumberBranch("event_number"); 
    //kernel.SetTreeAndFile("T", path_infile); 

    kernel.AddUserSlider("min_VDC_time", -400., 0., 401, "%3.0f ns"); 
    kernel.AddUserSlider("max_VDC_time", 0., +800., 801, "%3.0f ns"); 

    kernel.AddUserWindow("VDC", 1200, 400); 
    kernel.AddDrawFunction("VDC", "VDC", drawVDC, EDK::Frequency::Type::kEachEvent); 

    kernel.LaunchApp(); 
}

template<typename T> using Rv = ROOT::VecOps::RVec<T>; 

//__________________________________________________________________________________________________
void drawVDC()
{
    auto& kernel = EDK::Instance(); 

    using namespace APEX;
    

    std::vector<Rv<double>> wire, rawtime; 

    const double min_VDC = kernel.GetParameter<double>("min_VDC_time"); 
    const double max_VDC = kernel.GetParameter<double>("max_VDC_time"); 

    static const std::vector<std::string> plane_names{
        "u1", "v1", "u2", "v2"
    }; 
    static const std::string arm = is_RHRS ? "R" : "L"; 
    static const replay::EArmMode mode = is_RHRS ? replay::kRHRS : replay::kLHRS; 

    auto S2_hits_RHRS = replay::helpers::generate_S2_hits(
        replay::kRHRS, 
        kernel.GetData<Rv<double>>("R.s2.rt"),
        kernel.GetData<Rv<double>>("R.s2.lt")
    );
    auto S2_hits_LHRS = replay::helpers::generate_S2_hits(
        replay::kLHRS, 
        kernel.GetData<Rv<double>>("L.s2.rt"),
        kernel.GetData<Rv<double>>("L.s2.lt")
    );


    if (kernel.IsActive() && (S2_hits_RHRS.empty() || S2_hits_LHRS.empty())) return; 

    std::shared_ptr<EventHandler> dummy_event;  

    double beam_current = kernel.GetData<double>("hac_bcm_average"); 
    
    if (kernel.IsActive()) 
        dummy_event = std::make_shared<EventHandler>(mode, 
            beam_current, 4016,
            &S2_hits_RHRS.front(),
            &S2_hits_LHRS.front()
        ); 

    TCanvas* canv{nullptr}; 
    static TVirtualPad* pad[] = {nullptr, nullptr, nullptr, nullptr}; 
    if (kernel.IsActive()) {
        canv = kernel.GetCanvas(); 
        if (!pad[0]){
            canv->Divide(1,4., 0.01,0.0); 
            for (int i=0; i<APEX::VDC::kN_planes; i++) pad[i] = canv->cd(4-i); 
        }
    }

    for (int plane=0; plane<APEX::VDC::kN_planes; plane++) {

        const auto& pname = plane_names[plane]; 
        //wire   .emplace_back(kernel.GetData<Rv<double>>(arm+".vdc."+pname+".wire")); 
        //rawtime.emplace_back(kernel.GetData<Rv<double>>(arm+".vdc."+pname+".rawtime")); 

        const auto& wire    = kernel.GetData<Rv<double>>(arm+".vdc."+pname+".wire"); 
        const auto& rawtime = kernel.GetData<Rv<double>>(arm+".vdc."+pname+".rawtime"); 

        if (!kernel.IsActive()) continue; 
        
        auto mpad = pad[plane];
        mpad->cd();  
        mpad->SetTopMargin(0.);
        mpad->SetBottomMargin(0.); 
        auto frame = mpad->DrawFrame(
            0.,min_VDC, ((double)APEX::VDC::kN_wires_per_plane)-1.,max_VDC 
        ); 
        frame->SetTitle(";wire num;VDC time (ns)"); 

        Rv<double> pts_wire, pts_time; 

        //make all hits
        TGraph *g; 
        for (std::size_t i=0; i<wire.size(); i++) {
            
            APEX::VDC::Hit hit(is_RHRS, plane, (int)wire[i], rawtime[i], 0.); 
        
            pts_wire.emplace_back(hit.WireNum()); 
            pts_time.emplace_back(hit.Time()/units::ns);
        }
        g = new TGraph(wire.size(), pts_wire.data(), pts_time.data()); 
        g->SetMarkerStyle(kOpenCircle);  
        kernel.Draw(g, "SAME P"); 

        for (auto& S2_hit : (is_RHRS ? S2_hits_RHRS : S2_hits_LHRS)) {

            EventHandler event(
                mode, beam_current, 4016,
                is_RHRS ? &S2_hit : dummy_event->GetS2Hit(false), 
                is_RHRS ? dummy_event->GetS2Hit(true)  : &S2_hit
            );

            auto hit_groups = replay::helpers::group_vdc_hits(*dummy_event, plane, wire, rawtime); 

            for (const auto& group : hit_groups) {
                
                double t_S2 = S2_hit.Time()/units::ns;

                auto box = new TBox(
                    ((double)group.FirstWire()) - ((double)run_parameters::kGroup_max_gap)/2., 
                    t_S2 + run_parameters::kVDC_min_realTime/units::ns,
                    ((double)group.LastWire())  + ((double)run_parameters::kGroup_max_gap)/2.,
                    t_S2 + run_parameters::kVDC_max_realTime/units::ns
                ); 
                box->SetFillStyle(0);
                kernel.Draw(box, "SAME"); 
            }
        }

       
    }

}