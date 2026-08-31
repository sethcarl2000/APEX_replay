
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

namespace {

    std::string here(const char* const func) {
        static constexpr char classname[] = "Manager"; 
        std::ostringstream oss;
        oss << "<"<<classname<<"::"<<func<<">:";
        return oss.str();  
    }

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
  using RVecD = ROOT::VecOps::RVec<double>; 

  const bool kRHRS_bool = true;
  const bool kLHRS_bool = false; 
  
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
  ROOT::EnableImplicitMT(fN_threads);   

 
  //initialize the react vertex
  //create the 'TReactVertex' object, which computes the reaction vertex,
  // using raster information, and known v-wire positions. 
  ReactVertex *R_react_vertex, *L_react_vertex;

  const bool single_threadding = !run_parameters::kEnableMT;

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

  /*
  //now, we will generate the coinc-time.
  //this histogram will be the 'subtraction' of all coinc-times  
  double dt_sigma, dt_center; 

  if (are_both_arms_active) {

    auto hist_dt = rna.Get()
      
      .Define("dt", [](const RVec<S2Hit>& R_hits, const RVec<S2Hit>& L_hits)
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

    try {
      utils::fit_gaus_to_hist( (TH1D*)hist_dt->Clone("test_h"), 20*ns, dt_center, dt_sigma);
    } catch (const std::exception& e) {
      
      //if the fit fails, we will just make it so that there is no coinc-cut. 
      Warning(__func__, "Fitting gaussian to coincidence peak failed; A TR-TL time-coincidence cut will not be applied for this run.\nn. events in hist: %.2e\n what(): %s", hist_dt->Integral(), e.what()); 
      dt_center = utils::kNaN;
      dt_sigma  = utils::kNaN;  
    }

    std::printf(
      " central coinc-time (S2_r - S2_l): %.2f ns\n"
      " central coinc-time sigma: %.3f ns (cut is %.3f stddev-s from center)\n",
      dt_center/ns, 
      dt_sigma/ns, run_parameters::kS2_coinc_sigma_cut
    );  

    }*/
  
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
  
  EventCounter_RPtr nPass_1group_R, nPass_1pair_R, nPass_1raw_R, nPass_1ref_R; 
  EventCounter_RPtr nPass_1group_L, nPass_1pair_L, nPass_1raw_L, nPass_1ref_L; 

  //first, generate tracks in the right arm 
  if (is_RHRS_active) {
    //generate tracks
    helpers::generate_vdc_tracks(kRHRS, rna, 
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
    helpers::generate_vdc_tracks(kLHRS, rna, 
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

  //let's define some output branches 
  //_________________________________________________________________________________________________________________
  auto define_output_from_track = [&rna](std::string trk_branch, std::string output, double (VDC::Track::*method)() const) 
  {
    rna.DefineOutput(output, [method](const ROOT::RVec<VDC::Track>& tracks)
        { 
          RVec<double> vals; vals.reserve(tracks.size()); 
          for (const auto& trk : tracks) vals.push_back( (trk.*method)() ); 
          return vals; 
        }, {trk_branch}); 
  };
  //_________________________________________________________________________________________________________________

  std::string track_branch, arm; 

  if (is_RHRS_active) {

    track_branch = "R_tracks_refined"; 
    arm = "R";
    define_output_from_track(track_branch, arm+"_x_fp",     &VDC::Track::FP_x);
    define_output_from_track(track_branch, arm+"_y_fp",     &VDC::Track::FP_y);
    define_output_from_track(track_branch, arm+"_dxdz_fp",  &VDC::Track::FP_dx_dz);
    define_output_from_track(track_branch, arm+"_dydz_fp",  &VDC::Track::FP_dy_dz);
    
    define_output_from_track(track_branch, arm+"_s2_x",     &VDC::Track::S2_x);
    define_output_from_track(track_branch, arm+"_s2_y",     &VDC::Track::S2_y);

    //define_output_from_track(track_branch, arm+"_S2_x_param", &VDC::Track::xParam);
    define_output_from_track(track_branch, arm+"_S2_dt",    &VDC::Track::T0);

    //define_output_from_track(track_branch, arm+"_Eta", &VDC::Track::Get_Eta);

    rna.DefineOutput(arm+"_n_points", [](const RVec<VDC::Track>& tracks){
      RVec<int> ret; ret.reserve(tracks.size()); 
      for (const auto& trk : tracks) {
	  ret.push_back( trk.Get_nGoodPoints() ); 
      }
      return ret; 
    }, {track_branch});

    
    //separation between PMT times for our S2 hit
    rna.DefineOutput(arm+"_S2_pmt_dt", [](const RVec<VDC::Track>& tracks){
      RVecD ret; ret.reserve(tracks.size()); 
      for (const auto& trk : tracks) {
        
        const auto s2_hit = trk.GetEvent()->GetS2Hit(trk.IsRightArm());
        ret.push_back(s2_hit->DeltaT_raw());
      }
      return ret; 
    }, {track_branch});

    rna.DefineOutput(arm+"_S2_paddle", [](const RVec<VDC::Track>& tracks){
      RVec<int> ret; ret.reserve(tracks.size()); 
      for (const auto& trk : tracks) {
        
        const auto s2_hit = trk.GetEvent()->GetS2Hit(trk.IsRightArm());
        ret.push_back(s2_hit->Paddle());
      }
      return ret; 
    }, {track_branch});

    //separation between PMT times for our S2 hit
    rna.DefineOutput(arm+"_S2_is_twin_hit", [](const RVec<VDC::Track>& tracks){
      RVec<bool> ret; ret.reserve(tracks.size()); 
      for (const auto& trk : tracks) {
        
        const auto s2_hit = trk.GetEvent()->GetS2Hit(trk.IsRightArm());
        ret.push_back(s2_hit->Is_twinHit());
      }
      return ret; 
    }, {track_branch});

    rna.AddBranchToOutput("R_position_vtx");
    rna.DefineOutput("R_y_BPM", [](double bpma_y, double bpmb_y){ return (bpma_y + bpmb_y)/2; }, {"Rrb.BPMA.y","Rrb.BPMB.y"});  
  }

  if (is_LHRS_active) {
    track_branch = "L_tracks_refined"; 
    arm = "L";
    define_output_from_track(track_branch, arm+"_x_fp", &VDC::Track::FP_x);
    define_output_from_track(track_branch, arm+"_y_fp", &VDC::Track::FP_y);
    define_output_from_track(track_branch, arm+"_dxdz_fp", &VDC::Track::FP_dx_dz);
    define_output_from_track(track_branch, arm+"_dydz_fp", &VDC::Track::FP_dy_dz);
    
    define_output_from_track(track_branch, arm+"_S2_x", &VDC::Track::S2_x);
    define_output_from_track(track_branch, arm+"_S2_y", &VDC::Track::S2_y);

    //define_output_from_track(track_branch, arm+"_S2_x_param", &VDC::Track::xParam);
    define_output_from_track(track_branch, arm+"_S2_dt", &VDC::Track::T0);

    //define_output_from_track(track_branch, arm+"_Eta", &VDC::Track::Get_Eta);

    rna.DefineOutput(arm+"_n_points", [](const RVec<VDC::Track>& tracks){
      RVec<int> ret; ret.reserve(tracks.size()); 
      for (const auto& trk : tracks) {
	ret.push_back( trk.Get_nGoodPoints() ); 
      }
      return ret; 
    }, {track_branch});
    
    //separation between PMT times for our S2 hit
    rna.DefineOutput(arm+"_S2_pmt_dt", [](const RVec<VDC::Track>& tracks){
      RVecD ret; ret.reserve(tracks.size()); 
      for (const auto& trk : tracks) {
        
        const auto s2_hit = trk.GetEvent()->GetS2Hit(trk.IsRightArm());
        ret.push_back(s2_hit->DeltaT_raw());
      }
      return ret; 
    }, {track_branch});

    rna.DefineOutput(arm+"_S2_paddle", [](const RVec<VDC::Track>& tracks){
      RVec<int> ret; ret.reserve(tracks.size()); 
      for (const auto& trk : tracks) {
        
        const auto s2_hit = trk.GetEvent()->GetS2Hit(trk.IsRightArm());
        ret.push_back(s2_hit->Paddle());
      }
      return ret; 
    }, {track_branch});
    
    //separation between PMT times for our S2 hit
    rna.DefineOutput(arm+"_S2_is_twin_hit", [](const RVec<VDC::Track>& tracks){
      RVec<bool> ret; ret.reserve(tracks.size()); 
      for (const auto& trk : tracks) {
        
        const auto s2_hit = trk.GetEvent()->GetS2Hit(trk.IsRightArm());
        ret.push_back(s2_hit->Is_twinHit());
      }
      return ret; 
    }, {track_branch});

    rna.AddBranchToOutput("L_position_vtx");
    rna.DefineOutput("L_y_BPM", [](double bpma_y, double bpmb_y){ return (bpma_y + bpmb_y)/2; }, {"Lrb.BPMA.y","Lrb.BPMB.y"});  
  } 

  if (are_both_arms_active) {
    rna.DefineOutput("y_BPM", [](double R_y, double L_y){ return (R_y + L_y)/2.; }, {"R_y_BPM","L_y_BPM"}); 
  }

  //define some misc. information 
  define_constant_valued_output(rna, "run_number", fRunNumber); 
  define_constant_valued_output(rna, "rawfile_number", rawfile_number);
  define_constant_valued_output(rna, "segment_number", segment_number);
  
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

  EventCounter n_pass_1group_R, n_pass_1pair_R, n_pass_1raw_R, n_pass_1ref_R; 
  EventCounter n_pass_1group_L, n_pass_1pair_L, n_pass_1raw_L, n_pass_1ref_L; 
  
  if (is_RHRS_active) {
    n_pass_1group_R = *nPass_1group_R; 
    n_pass_1pair_R  = *nPass_1pair_R; 
    n_pass_1raw_R   = *nPass_1raw_R; 
    n_pass_1ref_R   = *nPass_1ref_R; 
  } 
  if (is_LHRS_active) {
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
  if (n_pass_1ref_L > 0) {

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
    "Cpu time:  %.3f s  ( %.6f ms/raw event )\n"
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
