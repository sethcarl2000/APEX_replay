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

    template<APEX::VDC::Plane P> constexpr bool is_U_plane() {
        return (P == APEX::VDC::Plane::U1 || P == APEX::VDC::Plane::U2); 
    }
}

/// @brief Draw VDC window
/// @param min_wire 
/// @param max_wire 
void drawVDC(const bool draw_whole); 

void drawVDC_whole() { drawVDC(true); }
void drawVDC_part()  { drawVDC(false); }

void drawS2(); 

void drawGridSearch(); 

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

    // '1' if low chamber, '0' otherwise 
    kernel.AddUserSlider("is_loChamber", 0., 1., 2, "%.0f"); 
    kernel.AddUserSlider("group_U_id", 0., 5., 6, "%.0f"); 
    kernel.AddUserSlider("group_V_id", 0., 5., 6, "%.0f"); 


    //kernel.AddUserWindow("VDC_whole", 1200, 900); 
    //kernel.AddDrawFunction("VDC_whole", "VDC_whole", drawVDC_whole, EDK::Frequency::Type::kEachEvent); 

    kernel.AddUserWindow("VDC_part", 1200, 900); 
    kernel.AddDrawFunction("VDC_part", "VDC_part", drawVDC_part, EDK::Frequency::Type::kEachEvent); 

    kernel.AddUserWindow("VDC_grid_search", 900, 900); 
    kernel.AddDrawFunction("VDC_grid_search", "VDC_grid_search", drawGridSearch, EDK::Frequency::Type::kEachEvent); 

    kernel.AddUserWindow("S2", 1200, 500); 
    kernel.AddDrawFunction("S2", "S2", drawS2, EDK::Frequency::Type::kEachEvent); 

    kernel.LaunchApp(); 
}

template<typename T> using Rv = ROOT::VecOps::RVec<T>; 

TGraph* draw_track(APEX::S2Event& evt, double slope, double intercept, double t_offset); 

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

    const APEX::VDC::Chamber selected_chamber = 
        (static_cast<int>(kernel.GetParameter<double>("is_loChamber")) == 0) 
        ? APEX::VDC::Chamber::Low : APEX::VDC::Chamber::High; 


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

    const std::size_t id_U = static_cast<std::size_t>(kernel.GetParameter<double>("group_U_id")); 
    const std::size_t id_V = static_cast<std::size_t>(kernel.GetParameter<double>("group_V_id")); 

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

    auto grid_search_results_lo = replay::helpers::chamber_grid_search<APEX::VDC::Chamber::Low>( S2_events.front() );
    auto grid_search_results_hi = replay::helpers::chamber_grid_search<APEX::VDC::Chamber::High>( S2_events.front() );

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

            const auto& S2_hit = S2_event.GetS2Hit(); 

            TVector3 S2_int_xyz( S2_hit.X(), S2_hit.Y(), S2_hit.Z() ); 

            const TVector3 S2_int_uvw = APEX::VDC::Track::Rotate_xyz_to_uvw( S2_int_xyz );

            const double S2_u{S2_int_uvw.X()}, S2_v{S2_int_uvw.Y()}, S2_w{S2_int_uvw.Z()}; 

            double t_S2 = S2_event.S2Time(); 
            std::printf(" t_S2: %.1f ns\n", t_S2/units::ns); 

            for(const auto& group : S2_event.GetCollection<APEX::VDC::HitGroup<P>>()) {

                auto box = new TBox(
                    ((double)group.LoEdge())  - APEX::VDC::kWireSpacing*((double)run_parameters::kGroup_max_gap)/2., 
                    (t_S2 + run_parameters::kVDC_min_realTime)/units::ns,
                    ((double)group.HiEdge())  + APEX::VDC::kWireSpacing*((double)run_parameters::kGroup_max_gap)/2.,
                    (t_S2 + run_parameters::kVDC_max_realTime)/units::ns
                ); 
                
                box->SetFillStyle(0);
                box->SetLineColor(kBlack); 
                box->SetLineStyle(kSolid); 
                
                std::size_t id; 
                if (selected_chamber == APEX::VDC::get_Chamber<P>()) {
                    if (selected_chamber == APEX::VDC::Chamber::Low) {
                        id = (P == APEX::VDC::Plane::U1) ? id_U : id_V; 
                    } else {
                        id = (P == APEX::VDC::Plane::U2) ? id_U : id_V; 
                    }

                    if (group.ID() == id) {

                        box->SetFillColor(kRed); 
                        box->SetFillStyle(3004);
                            
                        auto& gs_data = (APEX::VDC::is_loChamber<P>()) ? 
                            grid_search_results_lo : 
                            grid_search_results_hi; 
                            
                        double best_u, best_v; 
                        double best_min_eta=0.; 
                        for (const auto& dat : gs_data) {

                            double min_eta = std::min(dat.eta_u,dat.eta_v); 

                            if (min_eta > best_min_eta) {
                                best_min_eta = min_eta; 
                                best_u = dat.u; 
                                best_v = dat.v; 
                            }
                        }

                        // distance to S2 
                        double ToF = std::sqrt( 
                            utils::intpow<2>(S2_u - best_u) + 
                            utils::intpow<2>(S2_v - best_v) + 
                            utils::intpow<2>(S2_w - APEX::VDC::w<P>(is_RHRS)) 
                        )/units::speed_of_light_vacuum; 

                        double best_x = is_U_plane<P>() ? best_v : best_u; 

                        double slope = is_U_plane<P>() 
                            ? (S2_w - APEX::VDC::w<P>(is_RHRS))/(S2_v - best_x)
                            : (S2_w - APEX::VDC::w<P>(is_RHRS))/(S2_u - best_x); 

                        auto graph = draw_track(S2_event, slope, best_x, t_S2 - ToF); 
                        graph->SetLineColor(kRed); 
                        graph->SetLineWidth(2); 
                        kernel.Draw(graph, "SAME"); 
                    } 
                }
                kernel.Draw(box, "SAME"); 

            } // for (const auto& group ... 

            // draw tracks 
            std::printf("n. tracks: %zu\n", S2_event.size<APEX::VDC::Track>()); 
            static const std::size_t n_draw_pts = 100; 
            S2_event.RefereshEventPointers(); 
            for (const auto& track : S2_event.GetCollection<APEX::VDC::Track>()) {
                
                double slope = track.Slope<P>(); 
                double intercept = track.GetIntercept<P>(); 

                double x_min = intercept - APEX::run_parameters::max_drift_distance/slope; 
                double x_max = intercept + APEX::run_parameters::max_drift_distance/slope; 

                const double dx = (x_max - x_min)/((double)n_draw_pts-1); 
                double x = x_min; 

                std::vector<double> pts_x(n_draw_pts), pts_t(n_draw_pts); 
                for (int i=0; i<n_draw_pts; i++) {

                    pts_x[i] = x; 
                    pts_t[i] = (track.Get_T_model<P>(x) + t_S2)/units::ns; 
                    
                    x += dx; 
                }
                auto graph = new TGraph(n_draw_pts, pts_x.data(), pts_t.data()); 
                kernel.Draw(graph, "SAME");  
            
            } // for (const auto& track ... 

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

//__________________________________________________________________________________________________
void drawGridSearch()
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

    const APEX::VDC::Chamber selected_chamber = 
        (static_cast<int>(kernel.GetParameter<double>("is_loChamber")) == 0) 
        ? APEX::VDC::Chamber::Low : APEX::VDC::Chamber::High; 


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

    const std::size_t id_U = static_cast<std::size_t>(kernel.GetParameter<double>("group_U_id")); 
    const std::size_t id_V = static_cast<std::size_t>(kernel.GetParameter<double>("group_V_id")); 

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

    // take just the first S2 event
    auto& event = S2_events.front(); 

    auto grid_search_results_lo = replay::helpers::chamber_grid_search<APEX::VDC::Chamber::Low>( event );
    auto grid_search_results_hi = replay::helpers::chamber_grid_search<APEX::VDC::Chamber::High>( event );

    const double overlap_buffer = ((double)3.)*APEX::VDC::kWireSpacing; 
    // grid-spacing in u and v 
    static constexpr double dx = run_parameters::kGridSpacing; 


    APEX::VDC::for_each_chamber([&]<APEX::VDC::Chamber C>(){

        if (C != selected_chamber) return; 

        constexpr auto U = APEX::VDC::get_U<C>(); 
        constexpr auto V = APEX::VDC::get_V<C>(); 

        for (const auto& group_U : event.GetCollection<APEX::VDC::HitGroup<U>>()) {
            for (const auto& group_V : event.GetCollection<APEX::VDC::HitGroup<V>>()) {
                
                if (group_U.ID() != id_U || group_V.ID() != id_V) continue; 

                //do a brute-force grid-search
                const int n_steps_u = 1 + static_cast<int>((2.*overlap_buffer + group_V.Span())/dx); 
                const int n_steps_v = 1 + static_cast<int>((2.*overlap_buffer + group_U.Span())/dx); 

                const double u_min = group_V.LoEdge() - overlap_buffer; 
                const double v_min = group_U.LoEdge() - overlap_buffer; 
                
                const double u_max = group_V.HiEdge() + overlap_buffer; 
                const double v_max = group_U.HiEdge() + overlap_buffer; 

                if (group_U.ID() == id_U || group_V.ID() == id_V) {

                    auto hist = new TH2D(
                        "h_grid_search", 
                        "Chamber Grid-search: min( #eta_{u}, #eta_{v} );u;v;min( #eta_{u}, #eta_{v} )", 
                        n_steps_u, u_min-dx/2.,u_max+dx/2.,
                        n_steps_v, v_min-dx/2.,v_max+dx/2.
                    );
                    auto& gs_data = (APEX::VDC::is_loChamber<C>()) ? 
                        grid_search_results_lo : 
                        grid_search_results_hi; 

                    double best_u, best_v; 
                    double best_min_eta=0.; 
                    for (const auto& dat : gs_data) {

                        double min_eta = std::min(dat.eta_u,dat.eta_v); 

                        if (min_eta > best_min_eta) {
                            best_min_eta = min_eta; 
                            best_u = dat.u; 
                            best_v = dat.v; 
                        }

                        if (dat.id_U == id_U && dat.id_V == id_V) hist->Fill(dat.u, dat.v, min_eta); 
                    }

                    auto canv = kernel.GetCanvas(); 
                    canv->SetLeftMargin(0.125); canv->SetRightMargin(0.125); 
                    kernel.Draw(hist, "colz"); 
                }
            }   
        }
    });


}
//__________________________________________________________________________________________________
void drawS2()
{
    auto& kernel = EDK::Instance(); 

    using namespace APEX;

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

    const std::size_t id_U = static_cast<std::size_t>(kernel.GetParameter<double>("group_U_id")); 
    const std::size_t id_V = static_cast<std::size_t>(kernel.GetParameter<double>("group_V_id")); 

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

        
    //draw the S2 plane & paddle
    const double S2_x_range[2] = { 
        S2::get_paddle_X(is_RHRS,0)  - S2::paddle_width_X/2., 
        S2::get_paddle_X(is_RHRS,15) + S2::paddle_width_X/2.
    };

    static constexpr double paddle_width_y = 43.2*units::cm; 
    const double S2_y_range[2] = { -paddle_width_y/2., +paddle_width_y/2. };

    // buffer around total size of S2 
    double buffer = 1.1; 

    auto hist = new TH2D("h_S2", "S2 plane;x_{S2};y_{S2}", 
        100,  
        S2_x_range[0] - (S2_x_range[1]-S2_x_range[0])*(buffer - 1.)/2., 
        S2_x_range[1] + (S2_x_range[1]-S2_x_range[0])*(buffer - 1.)/2.,

        100,  
        S2_y_range[0] - (S2_y_range[1]-S2_y_range[0])*(buffer - 1.)/2., 
        S2_y_range[1] + (S2_y_range[1]-S2_y_range[0])*(buffer - 1.)/2.
    );

    kernel.Draw(hist); 

    for (int i=0; i<APEX::S2::N_paddles; i++) {

        auto box = new TBox(
            S2::get_paddle_X(is_RHRS,i) - S2::paddle_width_X/2., -paddle_width_y/2., 
            S2::get_paddle_X(is_RHRS,i) + S2::paddle_width_X/2., +paddle_width_y/2.
        ); 
        box->SetFillStyle(0); 
        kernel.Draw(box); 
    }

    double S2_y_rms = 2.5*units::cm; 
    for (const auto& S2_event : S2_events) {

        const auto& hit = S2_event.GetS2Hit(); 

        double S2_x = S2::get_paddle_X(is_RHRS,hit.Paddle()); 

        std::printf("S2 hit: (%.1f cm, %.1f cm)\n", hit.X()/units::cm, hit.Y()/units::cm); 

        auto box = new TBox(
            S2_x - S2::paddle_width_X/2., hit.Y() - S2_y_rms, 
            S2_x + S2::paddle_width_X/2., hit.Y() + S2_y_rms
        ); 
        box->SetLineColor(kBlack);
        box->SetFillStyle(0); 
        box->SetLineStyle(kDashed); 
        kernel.Draw(box); 

        // draw track 'crosshair' 
        TLine* line; 
        line = new TLine(
            S2_x - S2::paddle_width_X/2., 
            hit.Y(), 
            S2_x + S2::paddle_width_X/2., 
            hit.Y()
        ); 
        line->SetLineWidth(1); 
        line->SetLineColor(kBlack); 
        kernel.Draw(line); 

        line = new TLine(
            S2_x, 
            hit.Y() - S2_y_rms, 
            S2_x, 
            hit.Y() + S2_y_rms            
        ); 
        line->SetLineWidth(1); 
        line->SetLineColor(kBlack); 
        kernel.Draw(line); 
    }

    const double track_line_width = 5.*units::cm; 

    //now, draw tracks
    for (const auto& S2_event : S2_events) {
        for (const auto& track : S2_event.GetCollection<APEX::VDC::Track>()) {
            //plot the S2 x & y 
            double S2_x = track.S2_x(); 
            double S2_y = track.S2_y(); 
    
            // draw track 'crosshair' 
            TLine* line; 
            line = new TLine(
                S2_x, 
                S2_y-track_line_width/2., 
                S2_x, 
                S2_y+track_line_width/2.
            ); 
            line->SetLineWidth(2); 
            line->SetLineColor(kBlack); 
            kernel.Draw(line); 

            line = new TLine(
                S2_x-track_line_width/2., 
                S2_y, 
                S2_x+track_line_width/2., 
                S2_y
            ); 
            line->SetLineWidth(2); 
            line->SetLineColor(kBlack); 
            kernel.Draw(line); 
        }
    }
    
}

//__________________________________________________________________________________________________
TGraph* draw_track(APEX::S2Event& evt, double slope, double intercept, double t_offset)
{
    using namespace APEX;

    double x_min = intercept - APEX::run_parameters::max_drift_distance/slope;
    double x_max = intercept + APEX::run_parameters::max_drift_distance/slope;

    static constexpr int n_draw_pts = 100; 
    const double dx = (x_max - x_min)/((double)n_draw_pts-1); 
    std::vector<double> pts_x(n_draw_pts), pts_t(n_draw_pts); 

    for (int i=0; i<n_draw_pts; i++) {
        
        double x = x_min + ((double)i)*dx; 
        double w = slope * std::fabs(x - intercept); 
        
        double tau = evt.Drift_T(w, slope) + t_offset; 

        pts_x[i] = x; 
        pts_t[i] = tau/units::ns; 

        std::printf(" x,tau: %.1f, %.1f\n", x, tau); 
    }

    return new TGraph(n_draw_pts, pts_x.data(), pts_t.data()); 
}