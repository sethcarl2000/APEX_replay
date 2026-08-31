
#include "count_physics_events.C"

// ROOT headers
#include <TError.h> 
// stdlib headers
#include <string>
#include <fstream> 
#include <stdexcept> 

void create_eventcount_csv(const std::string& path_outfile, int min_run=3076, int max_run=5007)
{
  //try to open the output file 
  std::fstream outfile(path_outfile, std::ios::out | std::ios::trunc);

  if (!outfile.is_open()) {
    Error(__func__, "Unable to open csv output file: '%s'", path_outfile.c_str());
    return;
  }


  



}
