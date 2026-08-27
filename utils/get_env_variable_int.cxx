#include <APEX_utils.h>
// stdlib headers
#include <cstdlib> 
#include <mutex> 
#include <stdexcept> 
#include <iostream> 

namespace APEX
{
namespace utils
{

std::optional<int> get_env_variable_int(const std::string& param)
{
    //try to fetch the env variable.
    
    //use a mutex here, in case two threads try to call this at once (and some funny business might result...)
    static std::mutex mut; 

    mut.lock();
    const char* const value_ptr = std::getenv(param.c_str()); 
    mut.unlock(); 

    //return null if we didn't get a value. 
    if (!value_ptr) { return std::nullopt; }

    //now, try to convert it to an int.
    int value; 
    try {
        
        value = std::stoi(std::string{value_ptr}); 
    
    } catch (const std::exception& e) {

        std::cerr << "Warning in <APEX::utils::"<<__func__<<">: exception caught trying to convert env var '"<<param<<"'='"<<value_ptr<<"' to integer.\n"
        "   what(): " << e.what() << "\n"; 
        return std::nullopt; 
    }
    
    return value; 
}

}
}