
#include <APEX/decode.h>
#include <APEX/utils.h> 
// Podd headers
#include "THaRun.h"
#include "THaEvent.h" 
#include "THaAnalyzer.h"
// ROOT headers
#include <TError.h> 
// stdlib headers
#include <memory>
#include <iostream>
#include <fstream> 
#include <string>
#include <sstream> 
#include <stdexcept> 

namespace APEX
{
namespace decode
{

Long64_t count_physics_events(const std::string& path_input, std::string path_output_dump)
{ 
  using Status = utils::PathCheckStatus; 
  
  std::cout << "in <" << __func__ << ">: in body\n";

  bool status_ok=true;  

  //create the 'analyzer' 
  auto analyzer = std::make_unique<THaAnalyzer>(); 
  
  //create the 'event'
  auto event = std::make_unique<THaEvent>(); 
  
  //create the 'run' 
  auto run = std::make_unique<THaRun>(path_input.c_str(), "CODA input file");
  
  //only count physics events
  analyzer->SetCountMode( THaAnalyzer::kCountPhysics ); 

  //set the minimum verbosity for output to stdout 
  analyzer->SetVerbosity(0);


  
  if (path_output_dump.empty()) {

    Status status; 
    path_output_dump = utils::get_slurm_scratch_directory(&status);

    if (!status) {
      path_output_dump = "data/test_slurm_dir";  
      Warning(__func__, "We are NOT in a slurm job, so we will use path: %s", path_output_dump.c_str());   
    }
  }
  path_output_dump += "/dump.root"; 
  
  //if the scratch file already exists, that's bad!
  if (utils::is_path_regular_file(path_output_dump)) {

    Error(__func__, "File already exists under output dump path: %s", path_output_dump.c_str());
    return -1; 
  }
  
  
  //give the analyzer adummy output file. we'll put it on the scratch disk, unless we aren't in a slurm-job.
  analyzer->SetOutFile(path_output_dump.c_str()); 
  
  //now, actually perform the analysis. 
  Long64_t count = analyzer->Process( run.get() );
  
  if (count < 0) {

    std::ostringstream oss;
    oss << "in <"<<__func__<<">: something went wrong trying to count the number of physics events."
      " File: '" << path_input << "'";
    
    
    throw std::runtime_error(oss.str());
    return -1; 
  }

  std::cout << "in <" << __func__ << ">: exiting, "<<count<<" physics events counted.\n"; 

  //now, let's delete the 'dump' file:
  Status status; 
  utils::remove_file_from_disk(path_output_dump, &status); 

  if (!status) {
    Error(__func__, "Something went wrong trying to remove the 'dump' output file. Message: %s", status.message.c_str());
    return -1;
  }
  
  return count; 
}

}
}
