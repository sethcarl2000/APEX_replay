#ifndef vdc_track_replay_C
#define vdc_track_replay_C

//argparse
#include <include/argparse.hpp>
// APEX headers
// replay-utils lib
#include <replay_utils.h>
#include "replay_utils/DataLog.h"
//#define DEBUG_TRACK
#include "include/RDFNodeAccumulator.h"
#include "run_parameters.h"
#include <TapexReactVertex.h> 
#include <TapexS2Hit.h> 
#include <EventCounter.h> 
#include <ArmMode.h>
#include <ApexPidManager.h>
#include <UniqueIDGenerator.h>
#include "include/units.h" 
// ROOT headers
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/RSnapshotOptions.hxx>
#include <TFile.h> 
#include <TH1D.h> 
#include <TString.h>
#include <TError.h>  
#include <TStopwatch.h> 
#include <TDatime.h> 
// stdlib headers
#include <vector> 
#include <sstream> 
#include <string> 
#include <stdio.h> 
#include <iostream> 
#include <cmath> 
#include <map> 
#include <stdexcept> 
#include <fstream>
#include <stdlib.h> 
#include <sstream> 
#include <thread> 
#include <limits> 

namespace {
  constexpr bool kRHRS{true}, kLHRS{false};

  //wrap a parameter in a callable object, so you can feed it to an rdataframe
  template<typename T> std::function<T(void)> wrap_value(const T& value) {
    return (std::function<T(void)>)[&value]() { return value; };  
  };
  
  template<typename T> bool vec_is_nonempty(const ROOT::RVec<T>& v){ return !v.empty(); }

  const std::map<std::string, ArmMode::Bit> kArmModeOpt {
    {"RHRS", ArmMode::kRHRS},
    {"LHRS", ArmMode::kLHRS},
    {"both", ArmMode::kBoth} 
  }; 

  constexpr double kNaN = std::numeric_limits<double>::quiet_NaN(); 
}

/// @brief manages the conversion of 'raw' (decoded) data, and outputs VDC tracks
/// @param path_infile path to decoded data file
/// @param path_outfile path to output .root file
/// @param run_number the run number
/// @param rawfile_number the number of the rawfile to run 
/// @param segment_number the number of the segment to run 
/// @param arm_mode 'both' = replay both arms in coinc. mode. 'RHRS'/'LHRS' replay the R or L arm only
/// @param max_entries 0 = process all events, use multithreadding. n (>0) = process 'n' events, run in single-threadding mode 
/// @param n_threads number of threads to use in multi-threadding mode
int vdc_track_replay(
  const std::string path_infile, 
  const std::string path_outfile, 
  const int run_number,
  const int rawfile_number, 
  const int segment_number,
  const std::string arm_mode_str="both", //valid options are 'both', 'RHRS', "LHRS"
  const ULong64_t max_entries=0,
  const size_t n_threads=1
);


int main(int argc, char* argv[])
{
  //parse arguments
  argparse::ArgumentParser program("vdc_track_replay"); 

  program.add_argument("--run-number")
    .required()
    .help("index of run")
    .scan<'i', int>()
    .metavar("N")
    .nargs(1);

  program.add_argument("--rawfile-number")
    .help("rawfile number")
    .scan<'i', int>()
    .default_value(0)
    .metavar("R")
    .nargs(1);

  program.add_argument("--segment-number")
    .help("segment number")
    .scan<'i', int>()
    .default_value(0)
    .metavar("S")
    .nargs(1);

  
  program.add_argument("--arm-mode")
    .help("active-arm mode to operate in")
    .metavar("MODE")
    .choices("RHRS","LHRS","both")
    .default_value("both")
    .nargs(1);

  program.add_argument("--max-entries")
    .help("maximum raw-events to process. 0 => all events")
    .metavar("N")
    .scan<'i',unsigned int>()
    .default_value(0U)
    .nargs(1);

  program.add_argument("--n-threads")
    .help("number of threads to use in multi-threadding mode")
    .metavar("N")
    .scan<'i',int>()
    .default_value((int)std::thread::hardware_concurrency())
    .nargs(1);

  
  program.add_argument("stem_output")
    .required()
    .help("stem to output file destination. format: 'stem_output.[run].rawfile-[rawfile-num].seg-[seg-index].root")
    .nargs(1);

  program.add_argument("path_input")
    .required()
    .help("absolute path to input file")
    .nargs(1);
  
  try {
    program.parse_args(argc, argv);
  } catch (const std::exception& e) {
    std::cerr << e.what() << "\n"; 
    std::cerr << program; 
    return -1; 
  }
  
  std::string path_input          = program.get<std::string>("path_input"); 
  std::string stem_output         = program.get("stem_output");
  const int run_number            = program.get<int>("--run-number");
  const int rawfile_number        = program.get<int>("--rawfile-number");
  const int segment_number        = program.get<int>("--segment-number");
  std::string arm_mode_str        = program.get("--arm-mode");
  const unsigned int max_entries  = program.get<unsigned int>("--max-entries");
  const int n_threads             = program.get<int>("--n-threads");

  try {
      
    std::printf(
		"<vdc_track_replay>: There are %zi files to process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
		path_input.size()
		);
    
    std::string path_output = Form("%s.%i.raw-num-%i.seg-%i.root",
				   stem_output.data(),
				   run_number,
				   rawfile_number,
				   segment_number
				   );

    printf(
	   "<vdc_track_replay>: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
	   "                     Processing rawfile %i, segment %i/%zi\n"
	   "                     Input file '%s' ...\n"
	   "                     Output file '%s' ...\n",
	   rawfile_number, segment_number, path_input.size(), 
	   path_input.data(),
	   path_output.data() 
	   );
    
    int ret = vdc_track_replay(
			       path_input, 
			       path_output, 
			       run_number, 
			       rawfile_number, 
			       segment_number, 
			       arm_mode_str, 
			       max_entries, 
			       n_threads
			       );
    
    if (ret < 0) {
      Warning("vdc_tracK_replay", "Bad return code on replay %i/%i/%i (<0)", run_number, rawfile_number, segment_number);
    }
    
    
    printf(
	   "<replay_run_series>: done. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
	   ); 
    
  } catch (const std::exception& e) {
    std::cerr <<
      "<replay_run_series>: exception caught executing replay.\nwhat():" << e.what() << "\n"; 
    return -1; 
  }
  
  return 0; 
}


//________________________________________________________________________________________________________
int vdc_track_replay(
  const std::string path_infile, 
  const std::string path_outfile, 
  const int run_number,
  const int rawfile_number, 
  const int segment_number,
  const std::string arm_mode_str, //valid options are 'both', 'RHRS', "LHRS"
  const ULong64_t max_entries,
  const size_t n_threads
)
{
  using namespace replay_utils; 

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
  
  //make it so that, whenever we exit, this function is called (it will log all our data)
  auto& data_log = replay_utils::DataLog::instance(); 

  //register some data with the log file 
  data_log.Append(path_infile, '|');
  data_log.Append(path_outfile, '|');
  data_log.Append(run_number, '|');
  data_log.Append(rawfile_number, '|');
  data_log.Append(segment_number, '|');


  //get which arm we're dealing with 
  auto arm_it = kArmModeOpt.find(arm_mode_str); 
  if (arm_it == kArmModeOpt.end()) {
    ostringstream oss;
    for (const auto& [key, val] : kArmModeOpt) oss << "'" << key << "' "; 
    Error(__func__, "Arm mode option passed is invalid: '%s'. options are: %s",
	  arm_mode_str.c_str(),
	  oss.str().c_str()
	  );
    data_log.Append("invalid-arm-mode"); 
    data_log.BadExit(); 
    return -1; 
  }
  //
  const ArmMode::Bit arm_mode = arm_it->second; 

  //these are just for readability 
  auto RHRS_active      = [arm_mode](){ return (bool)(arm_mode & ArmMode::kRHRS); };
  auto LHRS_active      = [arm_mode](){ return (bool)(arm_mode & ArmMode::kLHRS); };
  auto both_arms_active = [arm_mode](){ return (bool)(arm_mode == ArmMode::kBoth); };
  
  //find the central-momentum of both arms using EPICS vars
  ROOT::EnableImplicitMT(n_threads); 
  ROOT::RDataFrame d_E("E", path_infile.data());

  //check if epics tree 'E' is empty
  if ((ULong64_t)*d_E.Count() < 1) {
    Error(__func__,"Epics tree 'E' of path \"%s\" is empty!",path_infile.data());
    data_log.Append("no-epics-data"); 
    data_log.BadExit(); 
    return -1; 
  }

  const double momentum_RHRS = *d_E.Mean("HacR_D1_P0rb") * (GeV / MeV);  
  const double momentum_LHRS = *d_E.Mean("HacL_D1_P0rb") * (GeV / MeV);

  Info(__func__, 
    "Momentum reported by RHRS/LHRS: %.1f / %.1f (MeV)", momentum_RHRS, momentum_LHRS
  ); 

  bool min_momentum_met=true; 
  if (RHRS_active() && (momentum_RHRS < run_parameters::min_momentum)) {
    data_log.Append("RHRS-momentum-too-low"); 
    min_momentum_met=false; 
  } 
  if (LHRS_active() && (momentum_LHRS < run_parameters::min_momentum)) {
    data_log.Append("LHRS-momentum-too-low"); 
    min_momentum_met=false; 
  } 
  if (!min_momentum_met) {
    Warning(__func__, "One or both arms failed to meet minimum spectrometer momentum threshold (%.1f MeV/c).", run_parameters::min_momentum); 
    //data_log.BadExit(); 
    //return -1; 
  }

  //initialize the react vertex
  //create the 'TReactVertex' object, which computes the reaction vertex,
  // using raster information, and known v-wire positions. 
  TapexReactVertex *R_react_vertex, *L_react_vertex;
  

  const bool single_threadding = max_entries > 0 || (!run_parameters::kEnableMT);

  //check status of multithreadding d
  if (!single_threadding) {
  
    if (!ROOT::IsImplicitMTEnabled()) ROOT::EnableImplicitMT(n_threads);
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
      data_log.Append("bad-react-vtx-R"); 
      data_log.BadExit();   
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
      data_log.Append("bad-react-vtx-L");   
      data_log.BadExit();   
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

    try {
      fit_gaus_to_hist( (TH1D*)hist_dt->Clone("test_h"), 20*ns, dt_center, dt_sigma);
    } catch (const std::exception& e) {
      
      //if the fit fails, we will just make it so that there is no coinc-cut. 
      Warning(__func__, "Fitting gaussian to coincidence peak failed; A TR-TL time-coincidence cut will not be applied for this run.\nn. events in hist: %.2e\n what(): %s", hist_dt->Integral(), e.what()); 
      data_log.Append("coinc-fit-failed");   
      dt_center = kNaN;
      dt_sigma  = kNaN;  
    }

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
          return gen_coinc_window_events(
            run_parameters::min_tr_tl, 
            run_parameters::max_tr_tl,
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
      data_log.Append("invalid-arm-mode");
      data_log.BadExit(); 
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
  
  ApexPidManager pid_manager; 
  
  EventCounter_RPtr nPass_1group_R, nPass_1pair_R, nPass_1raw_R, nPass_1ref_R; 
  EventCounter_RPtr nPass_1group_L, nPass_1pair_L, nPass_1raw_L, nPass_1ref_L; 

  //first, generate tracks in the right arm 
  if (RHRS_active()) {
    //generate tracks
    generate_vdc_tracks(kRHRS, rna, 
      nPass_1group_R, 
      nPass_1pair_R, 
      nPass_1raw_R, 
      nPass_1ref_R
    ); 
    //generate pid data
    gen_pid_data(kRHRS, rna, &pid_manager); 
  }

  //now, do the left-arm tracks
  if (LHRS_active()) {
    //geenrate tracks
    generate_vdc_tracks(kLHRS, rna, 
      nPass_1group_L, 
      nPass_1pair_L, 
      nPass_1raw_L, 
      nPass_1ref_L
    );
    //generate pid data
    gen_pid_data(kLHRS, rna, &pid_manager); 
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
    define_output_from_track(track_branch, arm+"_dxdz_fp", &ApexVDC::Track::FP_dx_dz);
    define_output_from_track(track_branch, arm+"_dydz_fp", &ApexVDC::Track::FP_dy_dz);
    
    define_output_from_track(track_branch, arm+"_s2_x", &ApexVDC::Track::S2_x);
    define_output_from_track(track_branch, arm+"_s2_y", &ApexVDC::Track::S2_y);

    define_output_from_track(track_branch, arm+"_S2_x_param", &ApexVDC::Track::xParam);
    define_output_from_track(track_branch, arm+"_S2_dt", &ApexVDC::Track::T0);

    //separation between PMT times for our S2 hit
    rna.DefineOutput(arm+"_S2_pmt_dt", [](const RVec<ApexVDC::Track>& tracks){
      RVecD ret; ret.reserve(tracks.size()); 
      for (const auto& trk : tracks) {
        
        const auto s2_hit = trk.GetEvent()->GetS2Hit(trk.IsRightArm());
        ret.push_back(s2_hit->DeltaT_raw());
      }
      return ret; 
    }, {track_branch});

    rna.DefineOutput(arm+"_S2_paddle", [](const RVec<ApexVDC::Track>& tracks){
      RVec<int> ret; ret.reserve(tracks.size()); 
      for (const auto& trk : tracks) {
        
        const auto s2_hit = trk.GetEvent()->GetS2Hit(trk.IsRightArm());
        ret.push_back(s2_hit->Paddle());
      }
      return ret; 
    }, {track_branch});

    //separation between PMT times for our S2 hit
    rna.DefineOutput(arm+"_S2_is_twin_hit", [](const RVec<ApexVDC::Track>& tracks){
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

  if (LHRS_active()) {
    track_branch = "L_tracks_refined"; 
    arm = "L";
    define_output_from_track(track_branch, arm+"_x_fp", &ApexVDC::Track::FP_x);
    define_output_from_track(track_branch, arm+"_y_fp", &ApexVDC::Track::FP_y);
    define_output_from_track(track_branch, arm+"_dxdz_fp", &ApexVDC::Track::FP_dx_dz);
    define_output_from_track(track_branch, arm+"_dydz_fp", &ApexVDC::Track::FP_dy_dz);
    
    define_output_from_track(track_branch, arm+"_S2_x", &ApexVDC::Track::S2_x);
    define_output_from_track(track_branch, arm+"_S2_y", &ApexVDC::Track::S2_y);

    define_output_from_track(track_branch, arm+"_S2_x_param", &ApexVDC::Track::xParam);
    define_output_from_track(track_branch, arm+"_S2_dt", &ApexVDC::Track::T0);

    //separation between PMT times for our S2 hit
    rna.DefineOutput(arm+"_S2_pmt_dt", [](const RVec<ApexVDC::Track>& tracks){
      RVecD ret; ret.reserve(tracks.size()); 
      for (const auto& trk : tracks) {
        
        const auto s2_hit = trk.GetEvent()->GetS2Hit(trk.IsRightArm());
        ret.push_back(s2_hit->DeltaT_raw());
      }
      return ret; 
    }, {track_branch});

    rna.DefineOutput(arm+"_S2_paddle", [](const RVec<ApexVDC::Track>& tracks){
      RVec<int> ret; ret.reserve(tracks.size()); 
      for (const auto& trk : tracks) {
        
        const auto s2_hit = trk.GetEvent()->GetS2Hit(trk.IsRightArm());
        ret.push_back(s2_hit->Paddle());
      }
      return ret; 
    }, {track_branch});
    
    //separation between PMT times for our S2 hit
    rna.DefineOutput(arm+"_S2_is_twin_hit", [](const RVec<ApexVDC::Track>& tracks){
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

  if (both_arms_active()) {
    rna.DefineOutput("y_BPM", [](double R_y, double L_y){ return (R_y + L_y)/2.; }, {"R_y_BPM","L_y_BPM"}); 

    //defined as (S2_r - S2_l)/sigma **for the hit associated with this specific track**
    
    /*rna.DefineOutput("R_track_dt", [dt_sigma, dt_center](RVec<ApexVDC::Track>& tracks)
    {
      RVec<double> dt; dt.reserve(tracks.size());
      for (auto& trk : tracks) {
	dt.push_back( (trk.GetEvent()->Get_Dt()-dt_center)/dt_sigma ); 
      }
      return dt; 
    }, {"R_tracks_refined"});

    
    rna.DefineOutput("L_track_dt", [dt_sigma, dt_center](RVec<ApexVDC::Track>& tracks)
    {
      RVec<double> dt; dt.reserve(tracks.size());
      for (auto& trk : tracks) {
	dt.push_back( (trk.GetEvent()->Get_Dt()-dt_center)/dt_sigma ); 
      }
      return dt; 
    }, {"L_tracks_refined"});*/ 
      
  }

  //define some misc. information 
  rna.DefineOutput("run_number", wrap_value(run_number), {}); 
  rna.DefineOutput("rawfile_number", wrap_value(rawfile_number), {}); 
  rna.DefineOutput("segment_number", wrap_value(segment_number), {}); 
  rna.DefineOutput("S2R_S2L_dt_center", wrap_value(dt_center), {});
  rna.DefineOutput("S2R_S2L_dt_sigma",  wrap_value(dt_sigma), {});

  
  rna.DefineOutput("beam_current", [](double beam_current){
    return beam_current;
  }, {"hac_bcm_average"}); 
  
  
  //rna.AddBranchToOutput("y_BPM"); 
  
  TStopwatch timer; 

  //turn this OFF; have it throw an exception instead. 
  rna.SetAbortOnError(false);

  try {  
    rna.Snapshot("track_data", path_outfile); 
  } catch (const std::exception& e) {
    Error(__func__, "Something went wrong trying to make a snapshot.\n what(): %s", e.what());
    data_log.Append("snapshot-threw-exception");
    data_log.BadExit(); 
    return -1; 
  }
  double elapsed  = timer.RealTime(); 
  double cpu_time = timer.CpuTime(); 
  
  data_log.Append("main-tree-success"); 

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
    "Real time: %.3f s  ( %.6f ms/raw event )\n"
    "Cpu time:  %.3f s  ( %.6f ms/raw event )\n",
    elapsed, 1000.*elapsed/((double)total_events),
    cpu_time, 1000.*cpu_time/((double)total_events)
  ); 

  //now, let's add some meta-data
  ROOT::RDataFrame meta_df(1); // <- create one 'event' (one entry in the tree for each file).  
  RDFNodeAccumulator rna_meta(meta_df);

  rna_meta.DefineOutput("run_number",        wrap_value(run_number), {});
  rna_meta.DefineOutput("rawfile_number",    wrap_value(rawfile_number), {});
  rna_meta.DefineOutput("segment_number",    wrap_value(segment_number), {});
  rna_meta.DefineOutput("S2R_S2L_dt_center", wrap_value(dt_center), {});
  rna_meta.DefineOutput("S2R_S2L_dt_sigma",  wrap_value(dt_sigma), {});

  rna_meta.DefineOutput("n_events_total", wrap_value(total_events), {}); 
  rna_meta.DefineOutput("n_events_s2hit", wrap_value(n_pass_1s2hit), {});   
  rna_meta.DefineOutput("n_events_coinc", wrap_value(n_pass_coinc), {}); 

  rna_meta.DefineOutput("RHRS_momentum",   wrap_value(momentum_RHRS), {});
  rna_meta.DefineOutput("LHRS_momentum",   wrap_value(momentum_LHRS), {});
  
  if (RHRS_active()) {
    rna_meta.DefineOutput("n_events_1group_R", wrap_value(n_pass_1group_R), {}); 
    rna_meta.DefineOutput("n_events_1pair_R",  wrap_value(n_pass_1pair_R), {}); 
    rna_meta.DefineOutput("n_events_1raw_R",   wrap_value(n_pass_1raw_R), {}); 
    rna_meta.DefineOutput("n_events_1ref_R",   wrap_value(n_pass_1ref_R), {});
  }
  if (LHRS_active()) {
    rna_meta.DefineOutput("n_events_1group_L", wrap_value(n_pass_1group_L), {}); 
    rna_meta.DefineOutput("n_events_1pair_L",  wrap_value(n_pass_1pair_L), {}); 
    rna_meta.DefineOutput("n_events_1raw_L",   wrap_value(n_pass_1raw_L), {}); 
    rna_meta.DefineOutput("n_events_1ref_L",   wrap_value(n_pass_1ref_L), {}); 
  }

  //record the date, and the amount of cpu time used 
  TDatime dtime;
  std::string date{ dtime.AsSQLString() }; 
  rna_meta.DefineOutput("run_process_date", wrap_value(date), {});
  rna_meta.DefineOutput("total_cpu_seconds", wrap_value(cpu_time), {});

  std::cout << "saving meta-data tree to file... \n" << std::flush;

  try {

    //use this option to make sure the file is not wiped & overwritten, but only added to. 
    ROOT::RDF::RSnapshotOptions opts;
    opts.fMode = "UPDATE";
    
    rna_meta.Snapshot("meta_data", path_outfile, &opts); 

  } catch (const std::exception& e) {
    Error(__func__, "Something went wrong trying to make the meta-data snapshot.\n what(): %s", e.what());
    data_log.Append("meta-data-tree-fail"); 
    data_log.BadExit(); 
    return -1; 
  }
  data_log.Append("meta-data-tree-success");   
  data_log.GoodExit(); 
    
  std::cout << "exiting..." << std::endl;
  return 0; 
}


#endif
