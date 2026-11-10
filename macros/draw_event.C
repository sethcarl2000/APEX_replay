// EventDisplayKernel
#include <EDK.hpp>
// APEX stuff
#include <APEX/units.h>
#include <APEX/replay.h>
#include <APEX/replay/helpers.h>
#include <APEX/run_parameters.h>
// Replay classes
#include <APEX/EventHandler.h>
#include <APEX/S2Event.h>
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
#include <TLine.h> 
#include <ROOT/RVec.hxx>
// stdlib
#include <string> 
#include <span> 
#include <vector> 



namespace {
    constexpr bool kRHRS{true}, kLHRS{false}; 
    // left arm 
    constexpr bool is_RHRS = kRHRS;  

    constexpr int run_number = 4016; 

    constexpr double NaN = std::numeric_limits<double>::quiet_NaN(); 
}

/// @brief Draw VDC window
/// @param min_wire 
/// @param max_wire 
void drawVDC(const bool draw_whole); 

void drawVDC_whole() { drawVDC(true); }
void drawVDC_part()  { drawVDC(false); }


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

    kernel.AddUserSlider("VDC_window_size", 0., (double)APEX::VDC::kN_wires_per_plane, APEX::VDC::kN_wires_per_plane, "%.0f"); 
    kernel.AddUserSlider("VDC_window_center", 0., (double)APEX::VDC::kN_wires_per_plane, APEX::VDC::kN_wires_per_plane, "%.0f"); 

    kernel.AddUserSlider("min_VDC_time", -400., 0., 401, "%3.0f ns"); 
    kernel.AddUserSlider("max_VDC_time", 0., +800., 801, "%3.0f ns"); 

    kernel.AddUserWindow("VDC_whole", 1200, 900); 
    kernel.AddDrawFunction("VDC_whole", "VDC_whole", drawVDC_whole, EDK::Frequency::Type::kEachEvent); 

    kernel.AddUserWindow("VDC_part", 1200, 900); 
    kernel.AddDrawFunction("VDC_part", "VDC_part", drawVDC_part, EDK::Frequency::Type::kEachEvent); 

    kernel.LaunchApp(); 
}

template<typename T> using Rv = ROOT::VecOps::RVec<T>; 

//__________________________________________________________________________________________________
void drawVDC(const bool draw_whole)
{
    auto& kernel = EDK::Instance(); 

    using namespace APEX;
    
    double window_center = std::round( kernel.GetParameter<double>("VDC_window_center") ); 
    double window_size   = std::round( kernel.GetParameter<double>("VDC_window_size") ); 

    // clamp values
    const int max_wire = APEX::VDC::kN_wires_per_plane-1; 
    const int wire_min = std::max<int>( std::min<int>(window_center - window_size/2., max_wire), 0 ); 
    const int wire_max = std::max<int>( std::min<int>(window_center + window_size/2., max_wire), 0 ); 

    std::vector<Rv<double>> wire, rawtime; 



    const double min_VDC = kernel.GetParameter<double>("min_VDC_time"); 
    const double max_VDC = kernel.GetParameter<double>("max_VDC_time"); 

    static const std::vector<std::string> plane_names{
        "u1", "v1", "u2", "v2"
    }; 
    static const std::string arm = is_RHRS ? "R" : "L"; 
    static const replay::EArmMode mode = is_RHRS ? replay::kRHRS : replay::kLHRS; 

    using RVecD = Rv<double>; 

    APEX::replay::helpers::RawVDCData VDC_data{
        { kernel.GetData<Rv<double>>(arm+".vdc.u1.wire"), kernel.GetData<Rv<double>>(arm+".vdc.u1.rawtime") },
        { kernel.GetData<Rv<double>>(arm+".vdc.v1.wire"), kernel.GetData<Rv<double>>(arm+".vdc.v1.rawtime") },
        { kernel.GetData<Rv<double>>(arm+".vdc.u2.wire"), kernel.GetData<Rv<double>>(arm+".vdc.u2.rawtime") },
        { kernel.GetData<Rv<double>>(arm+".vdc.v2.wire"), kernel.GetData<Rv<double>>(arm+".vdc.v2.rawtime") }
    };


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
    
    if (!kernel.IsActive()) return; 

    dummy_event = std::make_shared<EventHandler>(mode, 
        beam_current, run_number,
        &S2_hits_RHRS.front(),
        &S2_hits_LHRS.front()
    ); 

    APEX::replay::EventCounter_atomic d1,d2,d3,d4; 

    const auto & S2_hits = is_RHRS ? S2_hits_RHRS : S2_hits_LHRS; 
    Rv<S2Event> S2_events = APEX::replay::helpers::generate_vdc_tracks_from_S2Event(
        is_RHRS, 
        ROOT::RVec<EventHandler>{ *dummy_event }, 
        S2_hits, 
        VDC_data, 
        d1,d2,d3,d4, true
    ); 


    TH1F *plane_frames[APEX::VDC::kN_planes]; 

    TCanvas* canv{nullptr}; 
    TVirtualPad* pad[] = {nullptr, nullptr, nullptr, nullptr}; 
    if (kernel.IsActive()) {
        canv = kernel.GetCanvas(); 
        canv->Clear(); 
        canv->Divide(1,4, 0.01,0.00); 
        for (int i=0; i<APEX::VDC::kN_planes; i++) pad[i] = canv->cd(4-i); 
    }

    APEX::VDC::for_each_plane([&]<APEX::VDC::Plane P>(){

        const int pad_num = APEX::VDC::number<P>(); 
        const auto pname = APEX::VDC::plane_name<P>(); 
        
        auto mpad = pad[pad_num];
        mpad->cd();  
        mpad->SetTopMargin(0.);
        mpad->SetBottomMargin(0.);
        gStyle->SetOptStat(0); 

        std::printf("drawing frame for plane %s\n", pname.c_str()); 
        TH2D* frame; 
        double xmin, xmax; 
        if (draw_whole) {
            xmin = APEX::VDC::WirePos(is_RHRS, APEX::VDC::number<P>(), max_wire);
            xmax = APEX::VDC::WirePos(is_RHRS, APEX::VDC::number<P>(), 0); 
        } else {
            xmin = APEX::VDC::WirePos(is_RHRS, APEX::VDC::number<P>(), wire_max);
            xmax = APEX::VDC::WirePos(is_RHRS, APEX::VDC::number<P>(), wire_min);
        }
        frame = new TH2D(Form("frame_%s",pname.c_str()), "", 
            100, xmin, xmax,
            100, min_VDC, max_VDC 
        ); 
        
        kernel.Draw(frame); 

        if (draw_whole) {
            double center = kernel.GetParameter<double>("VDC_window_center"); 
            double window_size = kernel.GetParameter<double>("VDC_window_size"); 
            //draw a box on each VDC plane corresponding to the 'window' being drawn
            auto box = new TBox(
                APEX::VDC::WirePos(is_RHRS, APEX::VDC::number<P>(), wire_max), min_VDC,
                APEX::VDC::WirePos(is_RHRS, APEX::VDC::number<P>(), wire_min), max_VDC
            ); 

            box->SetLineColor(kRed); 
            box->SetLineStyle(kDashed); 
            box->SetFillStyle(0); 
            kernel.Draw(box); 
        }

        Rv<double> pts_wire, pts_time; 
        const auto& plane_data = std::get<APEX::replay::helpers::PlaneData<P>>( VDC_data ); 

        for (std::size_t i=0; i<plane_data.wire.size(); i++) {

            APEX::VDC::Hit<P> hit(is_RHRS, static_cast<int>(plane_data.wire[i]), plane_data.rawtime[i], 0.); 

            pts_wire.emplace_back( hit.WirePos() ); 
            pts_time.emplace_back( hit.Time()/units::ns );
        }
        //make all hits
        TGraph *g; 
        g = new TGraph(pts_wire.size(), pts_wire.data(), pts_time.data()); 
        g->SetMarkerStyle(kOpenCircle);  
        kernel.Draw(g, "SAME P"); 

        for (auto& S2_event : S2_events) {

            double t_S2 = S2_event.S2Time()/units::ns; 
            std::printf(" t_S2: %.1f ns\n", t_S2); 

            for(const auto& group : S2_event.GetCollection<APEX::VDC::HitGroup<P>>()) {

                auto box = new TBox(
                    ((double)group.LoEdge())  - APEX::VDC::kWireSpacing*((double)run_parameters::kGroup_max_gap)/2., 
                    t_S2 + run_parameters::kVDC_min_realTime/units::ns,
                    ((double)group.HiEdge())  + APEX::VDC::kWireSpacing*((double)run_parameters::kGroup_max_gap)/2.,
                    t_S2 + run_parameters::kVDC_max_realTime/units::ns
                ); 
                box->SetFillStyle(0);
                kernel.Draw(box, "SAME"); 
            }
        }
    }); 

    APEX::VDC::for_each_chamber([&]<APEX::VDC::Chamber C>(){

        for (auto& S2_event : S2_events) {

            double t_S2 = S2_event.S2Time();
            auto& pairs = S2_event.GetCollection<APEX::VDC::ChamberPair<C>>(); 

            int plane_num_U = APEX::VDC::number<APEX::VDC::get_U<C>()>();
            int plane_num_V = APEX::VDC::number<APEX::VDC::get_V<C>()>();

            double t_min = (t_S2 + APEX::run_parameters::kVDC_min_realTime)/units::ns; 
            double t_max = (t_S2 + APEX::run_parameters::kVDC_max_realTime)/units::ns;

            std::printf("n. pairs: %zu  t_min/max: %.1f ns / %.1f ns\n", pairs.size(), t_min, t_max); 

            for (const auto& pair : pairs) {

                pad[plane_num_U]->cd(); 

                TLine *line; 
                line = new TLine(pair.v(),t_min, pair.v(),t_max); 
                line->SetLineColor(kRed); 
                line->SetLineStyle(kDashed); 
                kernel.Draw(line); 

                pad[plane_num_V]->cd(); 

                line = new TLine(pair.u(),t_min, pair.u(),t_max); 
                line->SetLineColor(kRed); 
                line->SetLineStyle(kDashed); 
                kernel.Draw(line); 
            }
        }
    }); 


}