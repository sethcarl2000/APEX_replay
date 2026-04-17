// APEX headers
//#define DEBUG_TRACK
#include "include/RDFNodeAccumulator.h"
#include "run_parameters.h"
#include "functions/generate_vdc_tracks.h"
#include "functions/generate_S2_hits.h"
#include "functions/fit_gaus_to_hist.h"
#include "functions/gen_coinc_events.h"
#include "functions/gen_react_vertex.h"
#include <TapexReactVertex.h> 
#include <TapexS2Hit.h> 
#include <EventCounter.h> 
#include <ArmMode.h> 
#include "include/units.h" 
// ROOT headers
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h> 
#include <TH1D.h> 
#include <TString.h>
#include <TError.h>  
#include <TStopwatch.h> 
// stdlib headers
#include <vector> 
#include <string> 
#include <stdio.h> 
#include <iostream> 
#include <cmath> 
#include <map> 
#include <stdexcept> 

namespace {
  constexpr bool kRHRS{true}, kLHRS{false};

  template<typename T> bool vec_is_nonempty(const ROOT::RVec<T>& v){ return !v.empty(); }

  const std::map<std::string, ArmMode::Bit> kArmModeOpt {
    {"RHRS", ArmMode::kRHRS},
    {"LHRS", ArmMode::kLHRS},
    {"both", ArmMode::kBoth} 
  }; 
}
//____________________________________________________________________________________________________________________________
int vdc_track_replay(
  std::string path_infile, 
  std::string path_outfile, 
  int run_number, 
  std::string arm_mode_str="both", //valid options are 'both', 'RHRS', "LHRS"
  ULong64_t max_entries=0
)
{

  using namespace std; 
  using namespace ROOT::VecOps; 
  using RVecD = ROOT::RVec<double>; 
  using namespace units; 

  std::printf("<%s> " 
    "In funciton call body.\n"
    "--input file      %s\n"
    "--output file     %s\n"
    "--run nubmer      %i\n"
    "--event cap       %s\n",
    __func__,
    path_infile.data(), 
    path_outfile.data(), 
    run_number, 
    (max_entries > 0 ? Form("%lli",max_entries) : "no event cap") 
  ); 
  

  //get which arm we're dealing with 
  auto arm_it = kArmModeOpt.find(arm_mode_str); 
  if (arm_it == kArmModeOpt.end()) {
    Error(__func__, "Arm mode option passed is invalid: %s", arm_mode_str.c_str());
    return -1; 
  }
  //
  const ArmMode::Bit arm_mode = arm_it->second; 

  //these are just for readability 
  auto RHRS_active      = [arm_mode](){ return (bool)(arm_mode & ArmMode::kRHRS); };
  auto LHRS_active      = [arm_mode](){ return (bool)(arm_mode & ArmMode::kLHRS); };
  auto both_arms_active = [arm_mode](){ return (bool)(arm_mode == ArmMode::kBoth); };
  

  //find the central-momentum of both arms using EPICS vars
  ROOT::EnableImplicitMT(); 
  ROOT::RDataFrame d_E("E", path_infile.data());

  //check if epics tree 'E' is empty
  if ((ULong64_t)*d_E.Count() < 1) {
    Error(__func__,"Epics tree 'E' of path \"%s\" is empty!",path_infile.data());
    return -1; 
  }

  const double momentum_RHRS
    = d_E.Histo1D({"","",200,-1,-1},"HacR_D1_P0rb")->GetMean() * (GeV / MeV);
  
  const double momentum_LHRS
    = d_E.Histo1D({"","",200,-1,-1},"HacL_D1_P0rb")->GetMean() * (GeV / MeV);

  Info(__func__, 
    "Momentum reported by RHRS/LHRS: %.1f / %.1f (MeV)", momentum_RHRS, momentum_LHRS
  ); 

  //initialize the react vertex
  //create the 'TReactVertex' object, which computes the reaction vertex,
  // using raster information, and known v-wire positions. 
  TapexReactVertex *R_react_vertex, *L_react_vertex;
  

  const bool single_threadding = max_entries > 0 || (!run_parameters::kEnableMT);

  //check status of multithreadding d
  if (!single_threadding) {
  
    if (!ROOT::IsImplicitMTEnabled()) ROOT::EnableImplicitMT();
    Info(__func__, "Multithreadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize());
    
  } else {

    if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT();
    Info(__func__, "Multithreadding is disabled.");
  }
  
  RDFNodeAccumulator rna(ROOT::RDataFrame("T", path_infile.data())); 

  //try to construct the RHRS react vertex handler
  if (RHRS_active()) {
    printf("<%s> generating RHRS react vertex handler...\n",__func__); std::cout << std::flush; 
    
    try {
      R_react_vertex = new TapexReactVertex( kRHRS, path_infile, rna.Get() );
    } catch (const std::exception& e) {
      Error(__func__, "Something went wrong trying to construct the RHRS react-vertex handler.\n what(): %s", e.what()); 
      return -1; 
    }

    printf("<%s> done.\n",__func__); std::cout << std::flush; 
  } 
  
  //try to reconstruct the LHRS react vertex handler 
  if (LHRS_active()) {
    printf("<%s> generating LHRS react vertex handler...\n",__func__); std::cout << std::flush; 
    
    try {
      L_react_vertex = new TapexReactVertex( kLHRS, path_infile, rna.Get() );
    } catch (const std::exception& e) {
      Error(__func__, "Something went wrong trying to construct the LHRS react-vertex handler.\n what(): %s", e.what()); 
      return -1; 
    }
    
    printf("<%s> done.\n",__func__); std::cout << std::flush; 
  } 

  //get the total number of events
  const EventCounter total_events_in_file = *rna.Get().Count(); 

  //max number of events to run 
  const EventCounter total_events = (max_entries > 0) ? min<Long64_t>( max_entries, total_events_in_file ) : total_events_in_file; 

  std::printf("Processing %lli events.\n", total_events); 

  if (single_threadding) {
    rna = rna.Get().Range(0, total_events);
  }

  //create S2 hits
  if (RHRS_active()) {
    rna.Define("R_S2_hits", [](const RVec<double>& R_pmt, const RVec<double>& L_pmt)
      {
          return generate_S2_hits(kRHRS, R_pmt, L_pmt); 
      }, {"R.s2.rt", "R.s2.lt"});
    rna.Filter(vec_is_nonempty<TapexS2Hit>, {"R_S2_hits"});
  }
  
  if (LHRS_active()) {
    rna.Define("L_S2_hits", [](const RVec<double>& R_pmt, const RVec<double>& L_pmt)
      {
          return generate_S2_hits(kLHRS, R_pmt, L_pmt); 
      }, {"L.s2.rt", "L.s2.lt"}); 
    rna.Filter(vec_is_nonempty<TapexS2Hit>, {"L_S2_hits"}); 
  }
  //add logic here to process both S2 hits
  
  EventCounter_RPtr rptr_nPass_coinc = rna.Count(); 

  //now, we will generate the coinc-time.
  //this histogram will be the 'subtraction' of all coinc-times  
  double dt_sigma, dt_center; 

  EventCounter_RPtr rptr_nPass_1event; 

  if (both_arms_active()) {

    auto hist_dt = rna.Get()
      
      .Define("dt", [](const RVec<TapexS2Hit>& R_hits, const RVec<TapexS2Hit>& L_hits)
      {
        RVecD dt; dt.reserve(R_hits.size() * L_hits.size()); 

        for (const auto& R_hit : R_hits) {
          for (const auto& L_hit : L_hits) {
            dt.push_back( R_hit.Time() - L_hit.Time() );
          }
        }

        return dt; 
      
      }, {"R_S2_hits", "L_S2_hits"}) 
      
      .Histo1D<RVecD>({"h_dt", "coinc_times;s;", 200, -75*ns, +75*ns }, "dt");

    fit_gaus_to_hist( (TH1D*)hist_dt->Clone("test_h"), 20*ns, dt_center, dt_sigma);

    std::printf(
      " central coinc-time (S2_r - S2_l): %.2f ns\n"
      " central coinc-time sigma: %.3f ns (cut is %.3f stddev-s from center)\n",
      dt_center/ns, 
      dt_sigma/ns, run_parameters::kS2_coinc_sigma_cut
    );  

  }
  
  //create TeventHandler objets 
  switch (arm_mode) {

    //if both arms are active, then generate events based on coincidences between the left & right arms. 
    case ArmMode::kBoth : { //_________________________________________________________________________ 

      rna = rna.Get()
        .Define("coinc_events", [dt_center, dt_sigma](
          double beam_current, 
          unsigned int run_number, 
          const RVec<TapexS2Hit>& R_hits, 
          const RVec<TapexS2Hit>& L_hits )
        {
          return gen_coinc_events(
            dt_center, dt_sigma * run_parameters::kS2_coinc_sigma_cut, 
            beam_current, 
            run_number, 
            R_hits, 
            L_hits
          );
        }, {"hac_bcm_average","fEvtHdr.fRun","R_S2_hits","L_S2_hits"}) 
        
        //make a cut on events with at least 1 coincidence event
        .Filter(vec_is_nonempty<TapexEventHandler>, {"coinc_events"}); 
      
      break; 

    } //_________________________________________________________________________

    // create events for just the right arm
    case ArmMode::kRHRS : { //___________________________________________________

      rna = rna.Get() 
        .Define("coinc_events", []( 
          double beam_current, 
          unsigned int run_number, 
          const RVec<TapexS2Hit>& hits)
        {
          RVec<TapexEventHandler> events; events.reserve(hits.size());  
          for (const auto& hit : hits) {
            events.push_back(TapexEventHandler(kRHRS, beam_current, run_number, &hit, nullptr));  
          }
          return events; 
        }, {"hac_bcm_average","fEvtHdr.fRun","R_S2_hits"}); 

      break; 
    } //_________________________________________________________________________

    
    // create events for just the left arm
    case ArmMode::kLHRS : { //___________________________________________________

      rna.Define("coinc_events", []( 
          double beam_current, 
          unsigned int run_number, 
          const RVec<TapexS2Hit>& hits) 
        {
          RVec<TapexEventHandler> events; events.reserve(hits.size());  
          for (const auto& hit : hits) {
            events.push_back(TapexEventHandler(kLHRS, beam_current, run_number, nullptr, &hit));  
          }
          return events; 
        }, {"hac_bcm_average","fEvtHdr.fRun","L_S2_hits"});

      /*/if this parameter is true, we will only keep events which have one (and only one) coincidence between both arms 
      if (run_parameters::kKillMultiCoincEvents) {
        rna = rna.Get().Filter([](const RVec<TapexEventHandler>& vec){ return vec.size() != 1; }, {"coinc_events"}); 
      }*/ 
      break; 
    } //_________________________________________________________________________

    default : { 
      Error(__func__, "unsupported arm-mode %s, (it should not be possible to get here..)", arm_mode_str.c_str()); 
      return -1; 
      break;
    }
  
  }
  rptr_nPass_1event = rna.Count(); 

  EventCounter_RPtr nPass_1group_R, nPass_1pair_R, nPass_1raw_R, nPass_1ref_R; 
  EventCounter_RPtr nPass_1group_L, nPass_1pair_L, nPass_1raw_L, nPass_1ref_L; 

  //first, generate tracks in the right arm 
  if (RHRS_active()) {
    generate_vdc_tracks(kRHRS, rna, 
      nPass_1group_R, 
      nPass_1pair_R, 
      nPass_1raw_R, 
      nPass_1ref_R
    ); 
  }

  //now, do the left-arm tracks
  if (LHRS_active()) {
    generate_vdc_tracks(kLHRS, rna, 
      nPass_1group_L, 
      nPass_1pair_L, 
      nPass_1raw_L, 
      nPass_1ref_L
    );
  }

  //let's do some react-vertex calculations
  if (RHRS_active()) { gen_react_vertex(rna, R_react_vertex); }
  if (LHRS_active()) { gen_react_vertex(rna, L_react_vertex); }

  //let's define some output branches 
  //_________________________________________________________________________________________________________________
  auto define_output_from_track = [&rna](string trk_branch, string output, double (ApexVDC::Track::*method)() const) 
  {
    rna.DefineOutput(output, [method](const ROOT::RVec<ApexVDC::Track>& tracks)
        { 
          RVec<double> vals; vals.reserve(tracks.size()); 
          for (const auto& trk : tracks) vals.push_back( (trk.*method)() ); 
          return vals; 
        }, {trk_branch}); 
  };
  //_________________________________________________________________________________________________________________

  string track_branch, arm; 
  if (RHRS_active()) {
    track_branch = "R_tracks_refined"; 
    arm = "R";
    define_output_from_track(track_branch, arm+"_x_fp", &ApexVDC::Track::FP_x);
    define_output_from_track(track_branch, arm+"_y_fp", &ApexVDC::Track::FP_y);
    define_output_from_track(track_branch, arm+"_dxdz_fp", &ApexVDC::Track::dx_dz);
    define_output_from_track(track_branch, arm+"_dydz_fp", &ApexVDC::Track::dy_dz);
    rna.AddBranchToOutput("R_position_vtx");
    rna.DefineOutput("R_y_BPM", [](double bpma_y, double bpmb_y){ return (bpma_y + bpmb_y)/2; }, {"Rrb.BPMA.y","Rrb.BPMB.y"});  
  }

  if (LHRS_active()) {
    track_branch = "L_tracks_refined"; 
    arm = "L";
    define_output_from_track(track_branch, arm+"_x_fp", &ApexVDC::Track::FP_x);
    define_output_from_track(track_branch, arm+"_y_fp", &ApexVDC::Track::FP_y);
    define_output_from_track(track_branch, arm+"_dxdz_fp", &ApexVDC::Track::dx_dz);
    define_output_from_track(track_branch, arm+"_dydz_fp", &ApexVDC::Track::dy_dz);
    rna.AddBranchToOutput("L_position_vtx");
    rna.DefineOutput("L_y_BPM", [](double bpma_y, double bpmb_y){ return (bpma_y + bpmb_y)/2; }, {"Lrb.BPMA.y","Lrb.BPMB.y"});  
  } 

  if (both_arms_active()) {
    rna.DefineOutput("y_BPM", [](double R_y, double L_y){ return (R_y + L_y)/2.; }, {"R_y_BPM","L_y_BPM"}); 
  }

  
  //rna.AddBranchToOutput("y_BPM"); 
  
  TStopwatch timer; 

  //turn this OFF; have it throw an exception instead. 
  rna.SetAbortOnError(false);

  try {  
    rna.Snapshot("track_data", path_outfile); 
  } catch (const std::exception& e) {
    Error(__func__, "Something went wrong trying to make a snapshot.\n what(): %s", e.what());
    return -1; 
  }
  double elapsed  = timer.RealTime(); 
  double cpu_time = timer.CpuTime(); 
  
  EventCounter n_pass_1s2hit   = *rptr_nPass_coinc; 
  EventCounter n_pass_coinc    = *rptr_nPass_1event; 

  EventCounter n_pass_1group_R, n_pass_1pair_R, n_pass_1raw_R, n_pass_1ref_R; 
  EventCounter n_pass_1group_L, n_pass_1pair_L, n_pass_1raw_L, n_pass_1ref_L; 
  
  if (RHRS_active()) {
    n_pass_1group_R = *nPass_1group_R; 
    n_pass_1pair_R  = *nPass_1pair_R; 
    n_pass_1raw_R   = *nPass_1raw_R; 
    n_pass_1ref_R   = *nPass_1ref_R; 
  } 
  if (LHRS_active()) {
    n_pass_1group_L = *nPass_1group_L; 
    n_pass_1pair_L  = *nPass_1pair_L; 
    n_pass_1raw_L   = *nPass_1raw_L;
    n_pass_1ref_L   = *nPass_1ref_L; 
  }

  std::printf(
    "N. events with at least: --------------------------------------------\n"
    "   1 S2 hit in each arm:       %10llu  (%5.1f %%)\n"
    "   1 coinc with LHRS & RHRS:   %10llu  (%5.1f %%)\n",
    n_pass_1s2hit,    100.*((double)n_pass_1s2hit)/((double)total_events), 
    n_pass_coinc,     100.*((double)n_pass_coinc)/((double)total_events)
  ); 
  if (RHRS_active()) std::printf(
    " ~~ Right arm ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    "   1 hit group                 %10llu  (%5.1f %%)\n"
    "   1 hit group pair            %10llu  (%5.1f %%)\n"
    "   1 raw track                 %10llu  (%5.1f %%)\n"
    "   1 refined track             %10llu  (%5.1f %%)\n",
    n_pass_1group_R,  100.*((double)n_pass_1group_R)/((double)total_events),
    n_pass_1pair_R,   100.*((double)n_pass_1pair_R)/((double)total_events),
    n_pass_1raw_R,    100.*((double)n_pass_1raw_R)/((double)total_events),
    n_pass_1ref_R,    100.*((double)n_pass_1ref_R)/((double)total_events)
  ); 
  if (LHRS_active()) std::printf(
    " ~~ Left arm ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    "   1 hit group                 %10llu  (%5.1f %%)\n"
    "   1 hit group pair            %10llu  (%5.1f %%)\n"
    "   1 raw track                 %10llu  (%5.1f %%)\n"
    "   1 refined track             %10llu  (%5.1f %%)\n",
    n_pass_1group_L,  100.*((double)n_pass_1group_L)/((double)total_events),
    n_pass_1pair_L,   100.*((double)n_pass_1pair_L)/((double)total_events),
    n_pass_1raw_L,    100.*((double)n_pass_1raw_L)/((double)total_events),
    n_pass_1ref_L,    100.*((double)n_pass_1ref_L)/((double)total_events)
  ); 
  std::printf(
    "---------------------------------------------------------------------\n"
  );

  std::printf(
    "Real time: %.3f s  ( %.6f ms/event )\n"
    "Cpu time:  %.3f s  ( %.6f ms/event )\n",
    elapsed, 1000.*elapsed/((double)total_events),
    cpu_time, 1000.*cpu_time/((double)total_events)
  ); 

  std::cout << "exiting..." << std::endl;  

  return 0; 
}
