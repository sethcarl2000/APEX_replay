#ifndef slurm_info_h
#define slurm_info_h

#include <string> 
#include <cstdlib>

namespace slurm_info {
    
/// @return If we're in a job array, SLURM_ARRAY_JOB_ID, otherwise SLURM_JOB_ID. If we're not in a slurm job, return -1
inline static int job_id() { 
    //check if we're in an array. 
    const char* array_job_id_env = std::getenv("SLURM_ARRAY_JOB_ID");
    if (array_job_id_env) {
        return std::stoi( std::string{array_job_id_env} );
    } else {
        const char* job_id_env = std::getenv("SLURM_JOB_ID");
        if (job_id_env) {
            return std::stoi( std::string{job_id_env} );  
        } else {
            return -1; 
        }
    }        
}

/// @return If we're in a job array, SLURM_ARRAY_TASK_ID. Otherwise, -1
inline static int array_task_id() { 
    //check if we're in an array. 
    const char* array_task_id_env = std::getenv("SLURM_ARRAY_TASK_ID");
    if (array_task_id_env) {
        return std::stoi( std::string{array_task_id_env} );
    } else {
        return -1; 
    }        
}

};

#endif