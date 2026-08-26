#ifndef replay_APEX_coda_file_C
#define replay_APEX_coda_file_C

#include "decode/count_physics_events.h"
#include "decode/init_APEX_decoders.h" 
#include "decode/decode_APEX_coda_file.h"
#include "decode/check_existence_and_redability.h" 
// stdlib headers
#include <string>
#include <stdexcept> 

constexpr char* 


void replay_APEX_coda_file(const int run_number,
			   const int rawfile_number,
			   const std::string& path_output_stem,
			   const std::string& path_scratch_directory)
{

  //initialize decoding apps
  init_APEX_decoders();  
  
  //check for existence of the output file 
  //
  // .... implement later

  
  
  //construct the path to the coda file (input) 
  std::ostringstream oss;
  oss << "/cache/halla/apex/raw/apex_" << run_number << ".dat." << rawfile_number; 
  
  const std::string path_coda_file = oss.str(); 


  const std::string path_decode_file = path_scratch_directory + "decode.root";

  //first, we must count the number of events 
  Long64_t n_physics_events = count_physics_events(path_coda_file); 
  
}


#endif
