// APEX headers
#include <APEX/utils.h>
// ROOT headers
#include <TError.h> 
// stdlib headers
#include <filesystem>
#include <stdexcept> 
#include <sstream> 
#include <cstdio> 
#include <fstream> 


namespace APEX
{
namespace utils
{

std::vector<Rawfile> get_task_rawfile_list(const std::string& path, unsigned int task_id)
{
    //try to open the file 
    std::fstream infile(path, std::ios::in); 

    if (!infile.is_open()) {

        Error(__func__, "Unable to open file: %s", path.c_str());   
        return {}; 
    }

    std::vector<Rawfile> rawfiles{}; 

    std::string line; 

    while (std::getline(infile, line)) {

        std::istringstream iss(line); 

        unsigned int task, run_number, rawfile_number; 

        iss >> task >> run_number >> rawfile_number; 

        // add all rawfiles given to this task 
        if (task == task_id) { rawfiles.emplace_back(run_number, rawfile_number); } 
    }

    return rawfiles; 
}
    
}
}
