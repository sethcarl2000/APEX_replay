#include <APEX/utils.h>
// stdlib headers
#include <cstdlib> 
#include <mutex> 
#include <stdexcept> 
#include <iostream> 
#include <filesystem> 
#include <system_error> 

namespace APEX
{
namespace utils
{

PathCheckStatus is_path_directory(const std::string& path_str)
{
    //claude opus 5 showed me how to use the <filesystem> library to do this. this is more-or-less copied from a code snipet it gave me 
    // -seth, 28 aug 2026

    //use this for the filesystem namespace
    namespace fs = std::filesystem; 

    //try to fetch the env variable.
    PathCheckStatus result{"null"}; 

    auto path = fs::path{path_str}; 

    std::error_code error_code; 

    //check the status
    auto status = fs::status(path, error_code); 

    //quit if there was an error code
    if (error_code) {
        result.message = "status threw error: " + error_code.message(); 
        return result; 
    }

    //check if the path exists
    if (!fs::exists(status)) {
        result.message = "path does not exist"; 
        return result; 
    }

    //check if its a regular file 
    if (!fs::is_directory(status)) {
        result.message = "path is not a directory"; 
        return result; 
    }

    result.message = "success"; 
    return result; 
}

}
}