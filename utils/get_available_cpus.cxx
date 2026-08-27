#include <APEX_utils.h>
// stdlib headers 
#include <thread> 
#include <string> 
#include <stdexcept> 
#include <sstream> 


namespace APEX
{
namespace utils
{

unsigned int get_available_cpus()
{
    //try to fetch the env variable.
    auto n_cpus_slurm = get_env_variable_int("SLURM_CPUS_PER_TASK"); 

    if (!n_cpus_slurm) {

        // we **are not** in a slurm job. just use all the cpus we got. 
        return std::thread::hardware_concurrency(); 
    
    } else {

        return n_cpus_slurm.value(); 
    }
    
}

}
}