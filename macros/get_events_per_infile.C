
#include "replay_paths.h" 

#include <ROOT/RDataFrame.hxx>
#include <TString.h> 

#include <fstream> 
#include <string> 
#include <stdexcept> 

int get_events_per_infile(const char* path_output)
{
  using namespace std; 
  
  
  fstream outfile(path_output, ios::out | ios::trunc);   

  
  int run=3800; 
    
  for (const auto& path : replay_paths::list) {

    bool is_good=true;
    ULong64_t count=0;

    try { 

      std::printf("trying to open run %i - %i...", run, run+24);
      std::cout << std::flush; 
      
      ROOT::RDataFrame df("track_data", path); 
      
      count = *df.Count(); 

      std::printf("counted %llu events.\n", count);  
      
    } catch (const std::exception& e) {

      is_good=false;
      count=0; 
      Error("get_events_per_infile", "Something went wrong with this file");
      
    }
    outfile << Form("  {%4i, %4i, %8llu, %5s, \"%s\"},\n",
		    run,
		    run+24,
		    count,
		    (is_good?"true":"false"), path.c_str()
		    ); 

    run = run +25; 
  }

  outfile.close();

  return 0; 
}
