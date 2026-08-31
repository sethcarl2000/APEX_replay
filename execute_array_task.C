// APEX specific headers
#include <APEX/utils.h>
#include <APEX/utils/ErrorHandler.h>
#include <APEX/decode.h>
#include <APEX/replay/Manager.h>
// ROOT
#include <TStopwatch.h> 
#include <TError.h> 
// stdlib
#include <stdexcept> 
#include <iostream> 
#include <cstdio> 
#include <string>
#include <sstream>

namespace 
{
    // preferred size of decode 'chunks'. large CODA files's decoded output will be split into .root files with this many events.  
    constexpr Long64_t event_chunk_size = 200e3;

    // the maximum number of events to permit in one CODA file. this is to prevent a CODA file with 201k events from being split into two files with 200k and 1k events, respectively. 
    constexpr Long64_t max_chunk_size = 250e3;
}

void execute_array_task(const std::string& path_task_list, const std::string& output_directory)
{
    using namespace APEX; 

    //so, we're going to try and replay the given file. 

    //connect our custom error handler to ROOT's error handler interface 
    utils::ErrorHandler::Connect(); 

    //so, first, let's try to figure out what 'task' is ours: 
    unsigned int task_id; 
    auto task_id_result = utils::get_env_variable_int("SLURM_ARRAY_TASK_ID"); 

    if (!task_id_result) {

        task_id =0; 
        Error(__func__, "could not fetch 'SLURM_ARRAY_TASK_ID'! we must not be in a slurm task. choosign task ID=%u", task_id); 
    } else {

        task_id = task_id_result.value(); 
    }

    Info(__func__, "in body. SLURM_ARRAY_TASK_ID: %u. Fetching tasks...", task_id);

    //get list of rawfiles we've been tasked to process
    const auto rawfile_list = utils::get_task_rawfile_list(path_task_list, task_id); 

    if (rawfile_list.empty()) {

        Error(__func__, "No rawfiles find assigned to task %u. task list file: %s", task_id, path_task_list.c_str()); 
        return;
    }

    Info(__func__, "Task id %u has %zi rawfiles to process:", task_id, rawfile_list.size()); 

    std::printf("   run-number/rawfile-number\n");
    
    for (const auto& rawfile : rawfile_list) {
        std::printf("   %4u/%u\n", rawfile.run_number, rawfile.rawfile_number); 
    }

    //initialize the decoders
    decode::init_decoders(); 

    //so, let's loop over all rawfiles. 

    for (const auto rawfile : rawfile_list) {

        const unsigned int run_number     = rawfile.run_number; 
        const unsigned int rawfile_number = rawfile.rawfile_number; 

        Info(__func__, "Processing rawfile: run/rawfile: %u/%u ", run_number, rawfile_number); 

        const std::string path_cache = utils::get_env_variable_string("PATH_APEX_CACHE").value_or("null"); 
        
        std::ostringstream oss;

        oss << path_cache << "/apex_" << run_number << ".dat." << rawfile_number;
        const std::string path_coda_file = oss.str(); 
        
        
        //first, let's check the path of the file. 
        auto result = utils::is_path_regular_file(path_coda_file); 

        Info(__func__, "checking if path '%s' is regular file...\n", path_coda_file.c_str()); 

        if (!result) {
            Error(__func__, "coda file's path is not regular file: %s\n"
                "   message: %s", 
                path_coda_file.c_str(),
                result.message.c_str()
            ); 
            return; 
        }

        //now, let's count the number of physics CODA events to process. 
        Long64_t n_CODA_events; 
            
        Info(__func__, "Counting physics events..."); 
        try {

        n_CODA_events = decode::count_physics_events(path_coda_file); 
        
        } catch (const std::exception& e) {

            Error(__func__, "Error caught trying to count physics events.\n"
                "   CODA path:  %s\n"
                "   what:       %s",
                path_coda_file.c_str(),
                e.what()
            );
        return; 

        //tyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
        // - muon's comment 29 aug 26
        }
        Info(__func__, "Counted %lli physics events.", n_CODA_events); 


        //initialize the repaly manager
        replay::Manager replay_mgr; 
        
        replay_mgr.SetRunNumber(run_number); 
        
        replay_mgr.SetMaxNThreads(0); 

        replay_mgr.SetArmMode("both"); 

        
        //now, let's process these raw events. 
        Long64_t n_CODA_events_processed =0;
        int segment_number =0; 
        
        //get the decode file name
        std::string path_decode;
        
        using Status = utils::PathCheckStatus; 
        Status status; 

        path_decode = utils::get_slurm_scratch_directory(&status); 

        if (!status) {
        
        Error(__func__, "We must not be in a slurm job! exception caught: %s", status.message.c_str());
        path_decode = "data/test_slurm_dir"; 

        }
        
        path_decode = path_decode + "/decode.root"; 
        
        //now, let's loop over all 'event chunks' 
        while (n_CODA_events_processed < n_CODA_events) {

            Long64_t first_event = n_CODA_events_processed; 
            Long64_t last_event  = n_CODA_events_processed + event_chunk_size; 

            //if there's only a few events left, let's just process all remaining events. 
            if (n_CODA_events - first_event <= max_chunk_size) {
                last_event = first_event + max_chunk_size; 
            }

            Info(__func__, "Processing decode.\n"
                "   input path:     %s\n"
                "   output path:    %s\n"
                "   event range:    [%lli, %lli]",
                path_coda_file.c_str(),
                path_decode.c_str(),
                first_event, last_event
            );

            //now, process the decode. 
            try {

                decode::decode_coda_file(path_coda_file, path_decode, first_event, last_event); 
            
            } catch (const std::exception& e) {

                Error(__func__, "Error processing decode of coda file.\n"
                    "   what: %s", e.what()
                );
                return; 
            }
            Info(__func__, "Done with decode phase."); 

        //let's try to check how big the 'decode' file is (and that is exists)...
        utils::PathCheckStatus status; 

        double decode_size_MB = utils::get_file_size_KB(path_decode, &status) / (1024.);

        if (!status || !decode_size_MB) {
        Error(__func__, "Something went wrong checking size of decode file: '%s', size: %.2f MB. message: %s",
            path_decode.c_str(),
            decode_size_MB,
            status.message.c_str()
            );
        return;
        }

        Info(__func__, "decode file is %.2f MB", decode_size_MB); 
        

            auto replay_stem = output_directory + "/replay"; 

            //so, if all is well, then we just need to replay the file. 
            try {

                replay_mgr.Process(path_decode, replay_stem, rawfile_number, segment_number); 
            
            } catch (const std::exception& e) {

                Error(__func__, "Error replaying decoded file.\n"
                    "   what: %s", e.what()
                );
                return; 
            }

        Info(__func__, "Done with replay phase. trying to delete the 'decode.root' file..."); 

        status.message = "";
            //try to delete decode file 
            utils::remove_file_from_disk(path_decode, &status);

        if (!status) {
        Error(__func__, "Something went wrong trying to delete decode file.\n"
            "   message: %s", status.message.c_str()
            );
        return; 
        }
            //iterate segment number / event range 

            n_CODA_events_processed += event_chunk_size; 
            ++segment_number; 
        }

        Info(__func__, "Done processing rawfile: %u/%u ", run_number, rawfile_number); 
    }


}
