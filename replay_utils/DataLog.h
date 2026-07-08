#ifndef DataLog_h
#define DataLog_h

#include "replay_utils.h"
#include <TError.h> 

#include <fstream>
#include <sstream> 
#include <stdlib.h> 

namespace replay_utils
{

//meyers singleton class, which logs data to an output file. 
class DataLog {
private: 

    //keeps track of whether execution is aborted unexpectedly. 
    bool unexpected_exit{true};
    std::ostringstream replay_log; 

public: 

    //delete copy constructor & assignment operator (so no copies can be made)
    DataLog(const DataLog&) = delete; 
    DataLog& operator=(const DataLog&) = delete; 


    static DataLog& instance() { 
        static DataLog inst; 
        return inst; 
    } 

    
    template<typename T> void Append(const T& data, char delim=',') {
        replay_log << data << delim; 
    }; 

    void GoodExit() {
        unexpected_exit=false;
        replay_log << "good-exit\n"; 
    }
    void BadExit() {
        unexpected_exit=false;
        replay_log << "bad-exit\n";
    }

private:

    //upon program termination, log information to the logfile. 
    static void WriteInfo() {
        
        //get instance of DataLog
        auto& inst = instance(); 
        
        const char here[] = "logdata::WriteInfo"; 
        
        const char* path_apex_volatile = std::getenv("PATH_APEX_VOLATILE"); 
        
        std::string path_logfile; 

        if (!path_apex_volatile) {
        Error(here, "unable to parse environment variable PATH_APEX_VOLATILE, log file will be written in 'misc' directory"); 
            path_logfile = "misc/"; 
        } else {
            path_logfile = std::string{ path_apex_volatile } + "/slurm/replay_log/"; 
        }
        path_logfile += "replay_"; 

        //get the slurm array / job id
        const char* slurm_job_id = std::getenv("SLURM_ARRAY_JOB_ID"); 
        if (!slurm_job_id) {
            Warning(here, "variable SLURM_ARRAY_JOB_ID is null");
        } else {
            path_logfile += std::string{ slurm_job_id }; 
        }
        path_logfile += "_"; 

        const char* slurm_array_id = std::getenv("SLURM_ARRAY_TASK_ID"); 
        if (!slurm_array_id) {
            Warning(here, "variable SLURM_ARRAY_TASK_ID is null");
        } else {
            path_logfile += std::string{ slurm_array_id }; 
        }
        path_logfile += ".log";

        //now, try and and open the file 
        std::fstream logfile(path_logfile, std::ios::out | std::ios::app);  

        if (logfile.is_open()) {
            
            logfile << inst.replay_log.str();  
            
            //was this exit unexpected? 
            if (inst.unexpected_exit) {
                Warning(here, "unexpected exit / interrupted execution might have occured (neither GoodExit() nor BadExit() were called..)"); 
                logfile << "unexpected-exit\n"; 
            }
            logfile.close();                 

            Info(here, "Wrote log data to file '%s'", path_logfile.c_str()); 
            //inst.replay_log.clear(); 
            return; 

        } else {

            Error(here, "Unable to open logfile under path '%s'", path_logfile.c_str()); 
            return; 
        }
    }    

    //l;;;;;;;;;;;we
    // - Muon's comment 5 Jul 26 
    //constructor
    DataLog() {
        //tell the OS to execute this fcn when the program exits 
        std::atexit( DataLog::WriteInfo );
    }
};

};

#endif
