#include <TROOT.h>
#include <TStopwatch.h>
#include <TFile.h> 
#include <TTree.h> 
#include <TError.h>
#include <TString.h>

#include "vdc_track_replay.C" 

#include <vector>
#include <utility> 
#include <string>
#include <stdio.h> 

/// this macro is meant to compile the 'vdc_track_replay' macro, and process the
/// seires of runs which are input. 

/// @brief replays a series of runs 
/// @param run_number run number
/// @param rawfile_number index of the rawfile for this run. first rawfile is '0'
/// @param stem_output stem of the output file to make
/// @param paths_input list of input paths for the run
/// @param mode mode to run the arms in 
int replay_run_series(
  const int run_number,
  const int rawfile_number, 
  std::string stem_output, 
  const std::vector<std::string>& paths_input, 
  std::string mode )
{
  using namespace std; 
  
  printf(
    "<replay_run_series>: There are %zi files to process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
    paths_input.size()
  );

  int i_section=0; 
  for (const auto& path_input : paths_input) {
    
    string path_output = Form("%s.%i.raw-num-%i.seg-%i.root",
			      stem_output.data(),
			      run_number,
			      rawfile_number,
			      i_section 
			      );

    printf(
      "<replay_run_series>: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
      "                     Processing rawfile %i, segment %i/%zi\n"
      "                     Input file '%s' ...\n"
      "                     Output file '%s' ...\n",
      rawfile_number, i_section, paths_input.size(), 
      path_input.data(),
      path_output.data() 
    );
    
    vdc_track_replay(path_input, path_output, run_number, rawfile_number, i_section, mode, 0);

    ++i_section; 
  }

 
  
  printf(
    "<replay_run_series>: done. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
  ); 
  
  return 0; 
}
