
#include <APEX/replay.h>
#include <APEX/utils.h>
#include <APEX/units.h>
#include <APEX/replay/helpers.h>
#include <APEX/replay/Manager.h>
#include <APEX/run_parameters.h>
#include <APEX/ReactVertex.h>
#include <APEX/S2Hit.h>
#include <APEX/EventHandler.h>
#include <APEX/PidManager.h>
#include <APEX/replay/EventCounter.h> 
#include <APEX/replay/UniqueIDGenerator.h> 
// ROOT headers 
#include <ROOT/RVec.hxx>
#include <TStopwatch.h>
#include <TError.h>
// stdlib headers
#include <stdexcept> 
#include <cstdio> 
#include <cmath> 
#include <sstream> 
#include <string> 
#include <iostream>
#include <functional> 
#include <concepts> 

namespace {

    template<typename T> std::function<bool(const ROOT::VecOps::RVec<T>&)> vector_isnt_empty() {
        return [](const ROOT::RVec<T>& vec) { return !vec.empty(); };
    };

    /// @brief Add a constant value as a branch. 
    /// @tparam T type of the value
    /// @param rna the RNA to add it to 
    /// @param name the name of the new branch 
    /// @param val the value of the constnat  
    template<typename T> void define_constant_valued_output(RDFNodeAccumulator& rna, const std::string& name, const T& val) {
      
      //check to make sure this value can be copied
      static_assert(std::is_copy_constructible_v<T>, "type of constant 'T' is not copy constructible!"); 
      
      //now, add this value as a branch 
      rna.DefineOutput(name, [val](){ return val; }, {});       
    }

};

namespace APEX
{
namespace replay
{

/// concept requiruing that a callable function takes a VDC::Track as argument, and returns type 'T' 
template<typename T, typename F> concept TrackCallable = requires(F&& expr, const VDC::Track& t){
  { std::forward<F>(expr)(t) } -> std::same_as<T>; 
};

/// @brief Add an output as a member field of a track branch
/// @tparam T Return type of expression
/// @tparam F type of expression 
/// @param rna RDF node accumulator
/// @param event_branch name of event branch
/// @param name name of output branch 
/// @param expr expression to return result 
template<typename T, typename F> 
  requires TrackCallable<T,F>  
void define_track_output(RDFNodeAccumulator& rna, std::string event_branch, std::string name, F&& expr)
{
  //check if this output already exists
  if (rna.IsBranchDefined(name)) {
    throw std::invalid_argument("<define_track_output>: branch '"+name+"' already exists"); 
  }

  rna.DefineOutput(name, [&](ROOT::RVec<S2Event>& events){

    std::size_t n_tracks{0}; for (const auto& evt : events) n_tracks += evt.size<VDC::Track>(); 

    ROOT::RVec<T> out; out.reserve(n_tracks); 

    for (auto& evt : events) {
      evt.RefereshEventPointers(); 
      for (const auto& trk : evt.GetCollection<VDC::Track>()) {
        out.emplace_back( std::forward<F>(expr)(trk) ); 
      } 
    }
  
    return out; 

  }, {event_branch}); 

}


//______________________________________________________________________________________________________________
Manager::Manager()
{

    //this will choose the maximum number of threads available on the machine (unless the user overwrites this later)
    SetMaxNThreads(0);
    SetMaxEventsToProcess(0);  
}
//______________________________________________________________________________________________________________
void Manager::SetMaxNThreads(unsigned int n_threads)
{
    //get the number of cpus 
    unsigned int max_cpus = utils::get_available_cpus(); 

    if (n_threads < 1) {
        fN_threads = max_cpus; 
    } else {
        fN_threads = std::min( n_threads, max_cpus ); 
    }
}
//______________________________________________________________________________________________________________
void Manager::SetMaxEventsToProcess(Long64_t max_events)
{
    if (max_events < 1) {
        fMaxEventsToProcess = 1e9; 
    } else {
        fMaxEventsToProcess = max_events; 
    }
}
//______________________________________________________________________________________________________________
void Manager::SetArmMode(const std::string& mode)
{
    if (mode=="both") { fArmMode = kBothArms; return; }
    if (mode=="RHRS") { fArmMode = kRHRS;     return; }
    if (mode=="LHRS") { fArmMode = kLHRS;     return; }

    //if we got here, the user entered an invalid mode 
    std::ostringstream oss; 

    oss << "in <Manager::"<<__func__<<">: invalid arm mode given: '"<< mode <<"', must be 'RHRS', 'LHRS', or 'both'."; 
    throw std::runtime_error(oss.str());
}
//______________________________________________________________________________________________________________
int Manager::Process(const std::string& path_input, const std::string& stem_output, int rawfile_number, int segment_number)
{
  using namespace units; 
  using ROOT::VecOps::RVec; 
  
    //assemble the output file path
    std::ostringstream oss; 
    oss << stem_output << "." << fRunNumber << ".rawfile-" << rawfile_number << ".segment-" << segment_number << ".root"; 

    const std::string path_output = oss.str(); 

    Info(__func__, 
    "In funciton call body.\n"
    "--input file      %s\n"
    "--output file     %s\n"
    "--run nubmer      %i\n"
    "--event cap       %s\n",
    path_input.data(), 
    path_output.data(), 
    fRunNumber, 
    (fMaxEventsToProcess > 0 ? Form("%lli",fMaxEventsToProcess) : "no event cap") 
  ); 
  
  //these are just for readability 
  const bool is_RHRS_active = (fArmMode & kRHRS);
  const bool is_LHRS_active = (fArmMode & kLHRS);
  const bool are_both_arms_active = (fArmMode == kBothArms);
  
  //find the central-momentum of both arms using EPICS vars
  //ROOT::EnableImplicitMT(fN_threads);   

 
  //initialize the react vertex
  //create the 'TReactVertex' object, which computes the reaction vertex,
  // using raster information, and known v-wire positions. 
  ReactVertex *R_react_vertex, *L_react_vertex;

  const bool single_threadding = (run_parameters::kEnableMT == false) || fN_threads==1;

  //check status of multithreadding d
  if (!single_threadding) {
  
    if (!ROOT::IsImplicitMTEnabled()) ROOT::EnableImplicitMT(fN_threads);
    Info(__func__, "Multithreadding is enabled. Thread pool size: %i", ROOT::GetThreadPoolSize());
    
  } else {

    if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT();
    Info(__func__, "Multithreadding is disabled.");
  }
  
  RDFNodeAccumulator rna(ROOT::RDataFrame("T", path_input.data())); 

  rna.SetAbortOnError(false); 
  
  //try to construct the RHRS react vertex handler
  if (is_RHRS_active) {
    printf("<%s> generating RHRS react vertex handler...\n",__func__); std::cout << std::flush; 
    
    try {
      R_react_vertex = new ReactVertex( kRHRS, path_input, rna.Get() );
    } catch (const std::exception& e) {

      Error(__func__, "Something went wrong trying to construct the RHRS react-vertex handler.\n what(): %s", e.what()); 
      return -1; 
    }

    printf("<%s> done.\n",__func__); std::cout << std::flush; 
  } 
  
  //try to reconstruct the LHRS react vertex handler 
  if (is_LHRS_active) {
    printf("<%s> generating LHRS react vertex handler...\n",__func__); std::cout << std::flush; 
    
    try {
      L_react_vertex = new ReactVertex( kLHRS, path_input, rna.Get() );
    } catch (const std::exception& e) {

      Error(__func__, "Something went wrong trying to construct the LHRS react-vertex handler.\n what(): %s", e.what()); 
      return -1; 
    }
    
    printf("<%s> done.\n",__func__); std::cout << std::flush; 
  } 

  //get the total number of events
  const EventCounter total_events_in_file = *rna.Get().Count(); 

  if (total_events_in_file < 1) {
    Warning(__func__, "No events found in decode file.");
    return 0;
  } 
   
  //max number of events to run 
  const EventCounter total_events = (fMaxEventsToProcess > 0) ? std::min( fMaxEventsToProcess, (Long64_t)total_events_in_file ) : total_events_in_file; 

  Info(__func__, "Processing %lli events.\n", total_events); 

  if (single_threadding) {
    rna = rna.Get().Range(0, total_events);
  }

  //create S2 hits
  if (is_RHRS_active) {
    rna.Define("R_S2_hits", [](const RVec<double>& R_pmt, const RVec<double>& L_pmt)
      {
        return helpers::generate_S2_hits(kRHRS, R_pmt, L_pmt); 
      }, {"R.s2.rt", "R.s2.lt"});
    rna.Filter(vector_isnt_empty<S2Hit>(), {"R_S2_hits"});
  }
  
  if (is_LHRS_active) {
    rna.Define("L_S2_hits", [](const RVec<double>& R_pmt, const RVec<double>& L_pmt)
      {
        return helpers::generate_S2_hits(kLHRS, R_pmt, L_pmt); 
      }, {"L.s2.rt", "L.s2.lt"}); 
    rna.Filter(vector_isnt_empty<S2Hit>(), {"L_S2_hits"}); 
  }
  //add logic here to process both S2 hits
  
  EventCounter_RPtr rptr_nPass_coinc = rna.Count(); 
  
  EventCounter_RPtr rptr_nPass_1event; 

  
  //create TeventHandler objets 
  switch (fArmMode) {

    //if both arms are active, then generate events based on coincidences between the left & right arms. 
    case kBothArms : { //_________________________________________________________________________ 

      rna = rna.Get()
        .Define("coinc_events", [/*dt_center, dt_sigma*/](
          double beam_current, 
          unsigned int fRunNumber, 
          const RVec<S2Hit>& R_hits, 
          const RVec<S2Hit>& L_hits )
        {
          return helpers::gen_coinc_window_events(
            run_parameters::min_tr_tl, 
            run_parameters::max_tr_tl,
            beam_current, 
            fRunNumber, 
            R_hits, 
            L_hits 
          );
        }, {"hac_bcm_average","fEvtHdr.fRun","R_S2_hits","L_S2_hits"}) 
        
        //make a cut on events with at least 1 coincidence event
        .Filter(vector_isnt_empty<EventHandler>(), {"coinc_events"}); 
      
      break; 

    } //_________________________________________________________________________

    // create events for just the right arm
    case kRHRS : { //___________________________________________________

      rna = rna.Get() 
        .Define("coinc_events", []( 
          double beam_current, 
          unsigned int fRunNumber, 
          const RVec<S2Hit>& hits)
        {
          RVec<EventHandler> events; events.reserve(hits.size());  
          for (const auto& hit : hits) {
            events.push_back(EventHandler(kRHRS, beam_current, fRunNumber, &hit, nullptr));  
          }
          return events; 
        }, {"hac_bcm_average","fEvtHdr.fRun","R_S2_hits"}); 

      break; 
    } //_________________________________________________________________________

    
    // create events for just the left arm
    case kLHRS : { //___________________________________________________

      rna.Define("coinc_events", []( 
          double beam_current, 
          unsigned int fRunNumber, 
          const RVec<S2Hit>& hits) 
        {
          RVec<EventHandler> events; events.reserve(hits.size());  
          for (const auto& hit : hits) {
            events.push_back(EventHandler(kLHRS, beam_current, fRunNumber, nullptr, &hit));  
          }
          return events; 
        }, {"hac_bcm_average","fEvtHdr.fRun","L_S2_hits"});

      /*/if this parameter is true, we will only keep events which have one (and only one) coincidence between both arms 
      if (run_parameters::kKillMultiCoincEvents) {
        rna = rna.Get().Filter([](const RVec<EventHandler>& vec){ return vec.size() != 1; }, {"coinc_events"}); 
      }*/ 
      break; 
    } //_________________________________________________________________________

    default : { 
      Error(__func__, "unsupported arm-mode (it should not be possible to get here..)"); 
      return -1; 
      break;
    }
  
  }
  rptr_nPass_1event = rna.Count(); 

  //add a unique event id to each event that passes our S2-hit cuts
  UniqueIDGenerator unique_id_generator;

  rna.DefineOutput("event_id",   [&unique_id_generator](){
    return unique_id_generator.GetID();
  }, {}); 
  
  PidManager pid_manager; 
  
  EventCounter_atomic nPass_1group_R{0}, nPass_1pair_R{0}, nPass_1raw_R{0}, nPass_1ref_R{0}; 
  EventCounter_atomic nPass_1group_L{0}, nPass_1pair_L{0}, nPass_1raw_L{0}, nPass_1ref_L{0}; 

  //first, generate tracks in the right arm 
  if (is_RHRS_active) {
    //generate tracks
    helpers::generate_vdc_tracks_from_S2Event(kRHRS, rna, 
      nPass_1group_R, 
      nPass_1pair_R, 
      nPass_1raw_R, 
      nPass_1ref_R
    ); 
    //generate pid data
    helpers::gen_pid_data(kRHRS, rna, &pid_manager); 
  }

  //now, do the left-arm tracks
  if (is_LHRS_active) {
    //geenrate tracks
    helpers::generate_vdc_tracks_from_S2Event(kLHRS, rna, 
      nPass_1group_L, 
      nPass_1pair_L, 
      nPass_1raw_L, 
      nPass_1ref_L
    );
    //generate pid data
    helpers::gen_pid_data(kLHRS, rna, &pid_manager); 
  }

  //let's do some react-vertex calculations
  if (is_RHRS_active) { helpers::gen_react_vertex(rna, R_react_vertex); }
  if (is_LHRS_active) { helpers::gen_react_vertex(rna, L_react_vertex); }

  for (const bool is_RHRS : {true, false}) {

    const bool is_LHRS = !is_RHRS; 

    if (is_RHRS && (!is_RHRS_active)) continue; 
    if (is_LHRS && (!is_LHRS_active)) continue; 
 
    std::string arm = is_RHRS ? "R" : "L";
    std::string event_branch = arm+"_S2Events";
    
    define_track_output<double>(rna, event_branch, arm+"_x_fp",    [](const VDC::Track& t){ return t.FP_x(); });
    define_track_output<double>(rna, event_branch, arm+"_y_fp",    [](const VDC::Track& t){ return t.FP_y(); });
    define_track_output<double>(rna, event_branch, arm+"_dxdz_fp", [](const VDC::Track& t){ return t.FP_dx_dz(); });
    define_track_output<double>(rna, event_branch, arm+"_dydz_fp", [](const VDC::Track& t){ return t.FP_dy_dz(); });

    define_track_output<double>(rna, event_branch, arm+"_S2_x",    [](const VDC::Track& t){ return t.S2_x(); });
    define_track_output<double>(rna, event_branch, arm+"_S2_y",    [](const VDC::Track& t){ return t.S2_y(); });
    define_track_output<double>(rna, event_branch, arm+"_S2_dt",   [](const VDC::Track& t){ return t.S2_t(); });
    
    // S2 hit information 
    define_track_output<double>(rna, event_branch, arm+"_S2_xParam",    [](const VDC::Track& t){ return t.xParam(); });
    define_track_output<double>(rna, event_branch, arm+"_S2_hit_time",  [](const VDC::Track& t){ return t.GetEvent()->S2Time(); });
    define_track_output<double>(rna, event_branch, arm+"_S2_hit_PMT_dt_raw",[](const VDC::Track& t){ return t.GetEvent()->GetS2Hit().DeltaT_raw(); });
    define_track_output<int>   (rna, event_branch, arm+"_S2_hit_paddle",[](const VDC::Track& t){ return t.GetEvent()->GetS2Hit().Paddle(); });
    define_track_output<std::size_t>(rna, event_branch, arm+"_S2_hit_id", [](const VDC::Track& t){ return t.GetEvent()->GetID(); });

    // Cerenkov information 
    define_track_output<double>(rna, event_branch, arm+"_Cerenkov_x",    [](const VDC::Track& t){ return t.Cerenkov_x(); });
    define_track_output<double>(rna, event_branch, arm+"_Cerenkov_y",    [](const VDC::Track& t){ return t.Cerenkov_y(); });
    define_track_output<double>(rna, event_branch, arm+"_Cerenkov_t",    [](const VDC::Track& t){ return t.Cerenkov_t(); });

    // Raw Cerenkov data 
    rna.DefineOutput(arm+"_cer_a_c", [](const RVec<double>& v){ return v; }, {arm+".cer.a_c"});
    rna.DefineOutput(arm+"_cer_t_c", [](const RVec<double>& v){ return v; }, {arm+".cer.t_c"});

    //define_output_from_track(track_branch, arm+"_Eta", &VDC::Track::Get_Eta);
    rna.AddBranchToOutput(arm+"_position_vtx");
    rna.DefineOutput(arm+"_y_BPM", [](double bpma_y, double bpmb_y){ return (bpma_y + bpmb_y)/2; }, {arm+"rb.BPMA.y", arm+"rb.BPMB.y"});  
  }

  if (are_both_arms_active) {
    rna.DefineOutput("y_BPM", [](double R_y, double L_y){ return (R_y + L_y)/2.; }, {"R_y_BPM","L_y_BPM"}); 
  }

  //define some misc. information 
  define_constant_valued_output<int>(rna, "run_number", fRunNumber); 
  define_constant_valued_output<int>(rna, "rawfile_number", rawfile_number);
  define_constant_valued_output<int>(rna, "segment_number", segment_number);
  
  rna.DefineOutput("beam_current", [](double beam_current){
    return beam_current;
  }, {"hac_bcm_average"}); 
  
  TStopwatch timer; 

  //turn this OFF; have it throw an exception instead. 
  //rna.SetAbortOnError(false);

  try {  
    rna.Snapshot("track_data", path_output); 
  } catch (const std::exception& e) {
    Error(__func__, "Something went wrong trying to make a snapshot.\n what(): %s", e.what());
    return -1; 
  }
  double elapsed  = timer.RealTime(); 
  double cpu_time = timer.CpuTime(); 

  EventCounter n_pass_1s2hit   = *rptr_nPass_coinc; 
  EventCounter n_pass_coinc    = *rptr_nPass_1event; 

  EventCounter 
    n_pass_1group_R{nPass_1group_R},
    n_pass_1pair_R{nPass_1pair_R},
    n_pass_1raw_R{nPass_1raw_R},
    n_pass_1ref_R{nPass_1ref_R}, 

    n_pass_1group_L{nPass_1group_L},
    n_pass_1pair_L{nPass_1pair_L},
    n_pass_1raw_L{nPass_1raw_L},
    n_pass_1ref_L{nPass_1ref_L}; 

  std::printf(
    "N. events with at least: --------------------------------------------\n"
    "   1 S2 hit in each arm:       %10llu  (%5.1f %%)\n"
    "   1 coinc with LHRS & RHRS:   %10llu  (%5.1f %%)\n",
    n_pass_1s2hit,    100.*((double)n_pass_1s2hit)/((double)total_events), 
    n_pass_coinc,     100.*((double)n_pass_coinc)/((double)total_events)
  ); 
  if (is_RHRS_active) std::printf(
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
  if (is_LHRS_active) std::printf(
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

  //check to make sure we've kept at least 1 event
  if (nPass_1ref_L > 0) {

    Info(__func__, "Acceptable replay completed for file: %s", path_output.c_str());
    
  } else {

    Info(__func__, "No events kept in replay file: %s", path_output.c_str());

    utils::PathCheckStatus status;

    utils::remove_file_from_disk(path_output, &status);

    if (!status) {

      Error(__func__, "Something went wrong trying to delete replay file. message: %s", status.message.c_str());
      return -1; 
      
    }

  }

  std::printf(
    "Elapsed time:  %8.3f s  ( %.6f ms/raw event )\n"
    "Cpu time:      %8.3f s  ( %.6f ms/raw event )\n"
    "exiting...",
    elapsed, 1000.*elapsed/((double)total_events),
    cpu_time, 1000.*cpu_time/((double)total_events)
  ); 

  return 0; 
}
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________


}
}
