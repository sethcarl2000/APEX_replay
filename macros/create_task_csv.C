
// ROOT
#include <TError.h>
#include <TString.h> 
// stdlib
#include <fstream> 
#include <cstdio> 
#include <iostream> 
#include <string> 
#include <sstream>
#include <vector> 

struct Rawfile { unsigned int rawfile_number, n_events; };

struct Run { std::vector<Rawfile> rawfiles; unsigned int run_number, n_events=0; }; 

int create_task_csv(
    const std::string& path_input, 
    const std::string& path_output, 
    unsigned int first_run, 
    unsigned int last_run, 
    unsigned int max_events_per_task=7e6
) 
{   
    //try to parse csv. this is a quick & dirty way to do it, but it's what we have. 
    std::ifstream infile(path_input);

    if (!infile.is_open()) {
        Error(__func__, "Input file could not be opened: %s", path_input.c_str());
        return -1; 
    }

    std::fstream outfile(path_output, std::ios::out | std::ios::trunc); 

    if (!outfile.is_open()) {
        Error(__func__, "Output file could not be opened: %s", path_output.c_str());
        return -1; 
    }

    std::string line;
    
    std::vector<Run> runs{}; 

    Run* current_run=nullptr;

    //strip header
    std::getline(infile, line);
    
    while (std::getline(infile, line)) {

        //try to parse line 
        std::istringstream iss(line); 

        unsigned int run_number, rawfile_number, n_events; 
        
        // get data from the csv
        iss >> run_number >> rawfile_number >> n_events; 
        
        // add this rawifle to the run (or make a new run)
        if (!current_run || (current_run->run_number != run_number)) {

            Run new_run; 
            new_run.rawfiles.emplace_back(Rawfile{rawfile_number, n_events}); 
            new_run.run_number = run_number; 
            new_run.n_events   = n_events; 
	    
            runs.emplace_back(new_run); 

            current_run = &runs.back(); 

        } else {

 	    current_run->rawfiles.emplace_back(Rawfile{rawfile_number, n_events}); 
	  
	    current_run->n_events += n_events; 
        }
    } 

    //now, we're done with the infile. 
    infile.close(); 

    unsigned int n_events_in_task=0; 
    unsigned int task_id=0; 

    //headers for our csv 
    outfile << "task-id|run-number|rawfile-number\n"; 

    std::vector<Run> task_runs{}; 

    auto new_task = [&task_runs, &task_id, &n_events_in_task, &outfile]() {
        
        //make a new task 
        for (const auto& task_run : task_runs) {
            for (const auto& rawfile : task_run.rawfiles) {

                outfile << task_id << " " << task_run.run_number << " " << rawfile.rawfile_number << "\n"; 
            }
        } 

        task_runs.clear(); 
        ++task_id; 
        n_events_in_task =0; 	
    };

    for (const auto& run : runs) {
      
      if (run.run_number < first_run) continue; 
      
      if (run.run_number > last_run) break;  
      
      if (n_events_in_task + run.n_events > max_events_per_task) new_task();
	      
      task_runs.emplace_back(run); 
      n_events_in_task += run.n_events; 

    }
        
    
    if (!task_runs.empty()) new_task(); 
    
    outfile.close(); 
    return 0; 
}
