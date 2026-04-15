// APEX headers
#include "include/RDFNodeAccumulator.h"
#include "run_parameters.h"
#include "functions/generate_vdc_tracks.h"
#include "functions/generate_S2_hits.h"
#include "functions/fit_gaus_to_hist.h"
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
  Long64_t max_entries=-1
)
{
  printf("<%s> " 
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
  
  using namespace std; 
  using namespace ROOT::VecOps; 
  using RVecD = ROOT::RVec<double>; 
  using namespace units; 

  //get which arm we're dealing with 
  auto arm_it = kArmModeOpt.find(arm_mode_str); 
  if (arm_it == kArmModeOpt.end()) {
    Error(__func__, "Arm mode option passed is invalid: %s", arm_mode_str.c_str());
    return -1; 
  }
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
  TapexReactVertex react_vertex( kLHRS, path_infile );


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

  //get the total number of events
  const Long64_t total_events_in_file = *rna.Get().Count(); 

  //max number of events to run 
  const Long64_t total_events = (max_entries > 0) ? min<Long64_t>( max_entries, total_events_in_file ) : total_events_in_file; 

  printf("Processing %lli events.\n", total_events); 

  if (single_threadding) {
    rna = rna.Get().Range(0, total_events);
  }

  //create S2 hits
  if (RHRS_active()) {
    rna.Define("R_S2_hits", [](const RVec<double>& R_pmt, const RVec<double>& L_pmt)
      {
          return generate_S2_hits(kRHRS, R_pmt, L_pmt); 
      }, {"R.s2.rt", "R.s2.lt"});
  }
  
  if (LHRS_active()) {
    rna.Define("L_S2_hits", [](const RVec<double>& R_pmt, const RVec<double>& L_pmt)
      {
          return generate_S2_hits(kLHRS, R_pmt, L_pmt); 
      }, {"L.s2.rt", "L.s2.lt"}); 
  }
  //add logic here to process both S2 hits

  rna = rna.Get()
    .Filter(vec_is_nonempty<TapexS2Hit>, {"R_S2_hits"}) 
    .Filter(vec_is_nonempty<TapexS2Hit>, {"L_S2_hits"}); 

  //now, we will generate the coinc-time.
  //this histogram will be the 'subtraction' of all coinc-times  
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
      
      .Histo1D<RVecD>({"h_dt", "coinc_times;s;", 200, -60*ns, +60*ns }, "dt");

      //hist_dt->DrawCopy(); 
      double dt_center, dt_sigma;
      fit_gaus_to_hist( (TH1D*)hist_dt->Clone("test_h"), 20*ns, dt_center, dt_sigma, true);
      
  }
  
  EventCounter rptr_nPass_coinc = rna.Count(); 

  /*/first, generate tracks in the right arm 
  EventCounter nPass_1group_R, nPass_1pair_R, nPass_1raw_R, nPass_1real_R; 

  rna = generate_vdc_tracks(kRHRS, rna.Get(), 
    nPass_1group_R, 
    nPass_1pair_R, 
    nPass_1raw_R, 
    nPass_1real_R
  ); 

  //now, do the left-arm tracks
  EventCounter nPass_1group_L, nPass_1pair_L, nPass_1raw_L, nPass_1real_L; 

  rna = generate_vdc_tracks(kLHRS, rna.Get(), 
    nPass_1group_L, 
    nPass_1pair_L, 
    nPass_1raw_L, 
    nPass_1real_L
  );*/ 

  TStopwatch timer; 

  Long64_t n_pass_coinc = *rptr_nPass_coinc; 

  printf(
    "N. events with at least: --------------------------------------------\n"
    "   1 S2 hit in each arm:    %10llu  (%5.1f %%)\n"
    "---------------------------------------------------------------------\n",
    n_pass_coinc, 100.*((double)n_pass_coinc)/((double)total_events) 
  ); 

  double elapsed  = timer.RealTime(); 
  double cpu_time = timer.CpuTime(); 

  printf(
    "Real time: %.3f s  ( %.6f ms/event )\n"
    "Cpu time:  %.3f s  ( %.6f ms/event )\n",
    elapsed, 1000.*elapsed/((double)total_events),
    cpu_time, 1000.*cpu_time/((double)total_events)
  ); 

  return 0; 


}


