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
std::optional<std::string> get_env_variable_string(const std::string& param); 

/// @brief Fetch an environment variable, return value as int  
/// @param param name of the parameter
/// @return returned value for parameter if not-null. nullopt_t if param is null, or failed to be converted to an integer. 
std::optional<int> get_env_variable_integer(const std::string& param); 


/// @return the number of available cpus (wether in a slurm job or not)
unsigned int get_available_cpus(); 

}
}



#endif