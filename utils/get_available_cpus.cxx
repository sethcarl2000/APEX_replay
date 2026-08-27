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

int get_available_cpus()
{
    //try to fetch the env variable.
    auto n_cpus_slurm = get_env_variable("SLURM_CPUS_PER_TASK"); 

    if (!n_cpus_slurm) {

        // we **are not** in a slurm job. just use all the cpus we got. 
        return std::thread::hardware_concurrency(); 
    
    } else {

        // we **are** in a slurm job
        int n_cpus;

        try { 
            n_cpus = std::stoi(n_cpus_slurm.value()); 
        }
        catch (const std::exception& e) {
            
            //report this failure (this should not have happened!)
            std::ostringstream oss; 
            oss << "in <APEX::utils::"<<__func__<<">: exception caught trying to convert env value 'SLURM_CPUS_PER_TASK' to an integer.\n"
                " SLURM_CPUS_PER_TASK: '"<<n_cpus_slurm.value()<<"'\n"
                " what(): "<< e.what(); 
            
            throw std::runtime_error(oss.str()); 
            return -1; 
        }
        
        return n_cpus; 
    }
    
}

}
}