// APEX headers
#include "include/RDFNodeAccumulator.h"
#include "run_parameters.h"
#include "functions/generate_vdc_tracks.h"
#include "functions/generate_S2_hits.h"
#include <TapexReactVertex.h> 
#include <TapexS2Hit.h> 
#include <EventCounter.h> 
// ROOT headers
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h> 
#include <TString.h>
#include <TError.h>  
// stdlib headers
#include <vector> 
#include <string> 
#include <stdio.h> 
#include <iostream> 
#include <cmath> 

#define DISABLE_EPICS

using namespace std; 
using namespace ROOT::VecOps; 

namespace {
  constexpr bool kRHRS{true}, kLHRS{false};

  template<typename T> bool vec_is_nonempty(const ROOT::RVec<T>& v){ return !v.empty(); }
}

int vdc_track_replay(
  string path_infile, 
  string path_outfile, 
  int run_number, 
  Long64_t max_entries=-1
)
{
  Info(__func__, 
    "In funciton call body.\n"
    "--input file      %s\n"
    "--output file     %s\n"
    "--run nubmer      %i\n"
    "--event cap       %s\n",
    path_infile.data(), 
    path_outfile.data(), 
    run_number, 
    (max_entries > 0 ? Form("%lli",max_entries) : "no event cap") 
  ); 

#ifndef DISABLE_EPICS
  //find the central-momentum of both arms using EPICS vars
  ROOT::EnableImplicitMT(); 
  ROOT::RDataFrame d_E("E", path_infile.data());

  //check if epics tree 'E' is empty
  if ((ULong64_t)*d_E.Count() < 1) {
    Error(__func__,"Epics tree 'E' of path \"%s\" is empty!",path_infile.data());
    return -1; 
  }
    
  const double momentum_RHRS
    = d_E.Histo1D({"","",200,-1,-1},"HacR_D1_P0rb")->GetMean();
  
  const double momentum_LHRS
    = d_E.Histo1D({"","",200,-1,-1},"HacL_D1_P0rb")->GetMean();

  Info(__func__, 
    "Momentum reported by RHRS/LHRS: %.f / %.f (MeV)", momentum_RHRS, momentum_LHRS
  ); 
#endif 

  //initialize the react vertex
  //create the 'TReactVertex' object, which computes the reaction vertex,
  // using raster information, and known v-wire positions. 
  TapexReactVertex react_vertex( kLHRS, path_infile );


  //check status of multithreadding 
  if (run_parameters::kEnableMT && max_entries < 0 ) {
  
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


  //create S2 hits
  rna.Define("R_S2_hits", [](const RVec<double>& R_pmt, const RVec<double>& L_pmt)
    {
        return generate_S2_hits(kRHRS, R_pmt, L_pmt); 
    }, {"R.s2.rt", "R.s2.lt"});
  
  rna.Define("L_S2_hits", [](const RVec<double>& R_pmt, const RVec<double>& L_pmt)
    {
        return generate_S2_hits(kLHRS, R_pmt, L_pmt); 
    }, {"L.s2.rt", "L.s2.lt"}); 

  //add logic here to process both S2 hits

  rna = rna.Get()
    .Filter(vec_is_nonempty<TapexS2Hit>, {"R_S2_hits"}) 
    .Filter(vec_is_nonempty<TapexS2Hit>, {"L_S2_hits"}); 
  
  EventCounter nPass_coinc = rna.Count(); 

  //first, generate tracks in the right arm 
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
  ); 

  printf("n. events with at least 1 coinc in each arm: %llu\n", *nPass_coinc); 

  return 0; 


}


