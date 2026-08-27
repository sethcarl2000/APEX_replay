#ifndef APEX_utils_h
#define APEX_utils_h

#include <string> 
#include <optional> 

namespace APEX
{
namespace utils
{

/// @brief Fetch an environment variable, return value as string  
/// @param param name of the parameter
/// @return returned value for parameter if not-null, or nullopt_t if param is null. 
std::optional<std::string> get_env_variable(const std::string& param); 

/// @return the number of available cpus (wether in a slurm job or not)
int get_available_cpus(); 

}
}



#endif