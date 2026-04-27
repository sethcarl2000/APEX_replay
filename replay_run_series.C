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


int replay_run_series(
  const int run_number,
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
    
    string path_output = Form("%s.%i.part-%i.root",
			      stem_output.data(),
			      run_number,
			      i_section 
			      );

    printf(
      "<replay_run_series>: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
      "                     Processing sub-part %i/%zi\n"
      "                     Input file '%s' ...\n"
      "                     Output file '%s' ...\n",
      i_section, paths_input.size(), 
      path_input.data(),
      path_output.data() 
    );
    
    vdc_track_replay(path_input, path_output, run_number, i_section, mode, 0);

    ++i_section; 
  }

 
  
  printf(
    "<replay_run_series>: done. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
  ); 
  
  return 0; 
}
