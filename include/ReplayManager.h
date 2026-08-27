
#include <APEX_replay.h>
// ROOT headers
#include <RtypesCore.h> 
// stdlib headers
#include <string> 

namespace APEX
{
namespace replay
{

class ReplayManager {
private: 

    // if this is <= 0, then we process all 
    Long64_t fMaxEventsToProcess{-1}; 

    // if this is '0', then use all available cpus concurrently 
    unsigned int fN_threads{0}; 

    int fArmMode{kBothArms}; 
    int fRunNumber{0}; 

public: 
    //default ctor. 
    ReplayManager(); 

    //we don't want anyone to make copies!!! (this deletes the copy-ctor)
    ReplayManager(const ReplayManager&) = delete; 
    ReplayManager& operator=(const ReplayManager&) = delete; 

    /// @brief Process a replay. throws std::runtime_error and quits on failure. 
    /// @param path_input absolute path to the input 'decoded' .root file. 
    /// @param stem_output 'stem' of the output repaly file. The final file output will be under the path [stem_output].[run-num].rawfile-[num].segment-[num].root 
    /// @param rawfile_number index of the rawfile (starts at 0)
    /// @param segment_number index of the rawfile 'segment' for rawfiles that had to be split into smaller pieces (starts at 0)
    void Process(const std::string& path_input, const std::string& stem_input, int rawfile_number, int segment_number);


    //setters / getters

    /// @brief Set the max. number of threads. true maximum is capped at number of threads available in slurm job / on machine. 
    /// @param n_threads max n. threads to run concurrently. 0 = use all available concurrent threads 
    void SetMaxNThreads(unsigned int n_threads=0);

    /// @brief Set max number of events to process. This can only be done in single-thread mode, and if its called, it sets the max number of threads to 1. 
    /// @param max_events 
    void SetMaxEventsToProcess(Long64_t max_events); 

    /// @brief Set the 'arm mode' to use. 
    /// @param mode valid options are 'RHRS', 'LHRS', and 'both' 
    void SetArmMode(const std::string& mode); 
}; 


}
}
