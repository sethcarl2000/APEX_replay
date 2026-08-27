#include <APEX_utils.h>
// stdlib headers
#include <cstdlib> 
#include <mutex> 


namespace APEX
{
namespace utils
{

std::optional<std::string> get_env_variable_string(const std::string& param)
{
    //try to fetch the env variable.
    
    //use a mutex here, in case two threads try to call this at once (and some funny business might result...)
    static std::mutex mut; 

    mut.lock();
    const char* const value_ptr = std::getenv(param.c_str()); 
    mut.unlock(); 

    //return null if we didn't get a value. 
    if (!value_ptr) { return std::nullopt; }

    return std::string{value_ptr}; 
}

}
}