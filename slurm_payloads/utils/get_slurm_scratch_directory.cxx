// APEX headers
#include <APEX/utils.h>
// stdlib headers
#include <filesystem>
#include <stdexcept> 
#include <sstream> 
#include <cstdio> 


namespace APEX
{
namespace utils
{

std::string get_slurm_scratch_directory(PathCheckStatus* status)
{

  //so, we're going to use what we know about ifarm's slurm disk formatting. 
    
  // when in a slurm job, and you've asked for some space on a 'scratch disk' with '--gres:disk...', this is where slurm makes a new directory for you:
    
  // /scracth/slurm/[SLURM_JOB_ID]

  auto slurm_job_id = get_env_variable_int("SLURM_JOB_ID");

  if (!slurm_job_id) {
    std::ostringstream oss; 
    
    oss << "in <APEX::utils::" << __func__ <<">: failed to fetch slurm job id! we must not be in a slurm job."; 

    if (status == nullptr) { 
      throw std::runtime_error(oss.str()); 
    } else {
      status->message = oss.str();
    } 
    return ""; 
  }
  
  char buff[100]; 
  std::sprintf(buff, "/scratch/slurm/%i", slurm_job_id.value()); 

  if (status) status->message = "success"; 
  
  return std::string{buff}; 
    
}
    
}
}
