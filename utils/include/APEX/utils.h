#ifndef APEX_utils_h
#define APEX_utils_h

// ROOT headers
#include <RtypesCore.h>
// stdlib headers
#include <string> 
#include <optional> 
#include <limits> 

class TH1; 

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
std::optional<int> get_env_variable_int(const std::string& param); 

/// @return the number of available cpus (wether in a slurm job or not)
unsigned int get_available_cpus(); 
  
/// @return a 'FileCheck' 
struct PathCheckStatus { 

  std::string message{""};

  operator bool() const { return message=="success"; }; 
};

/// @brief Check if the given path 1.) exists and 2.) is a regular file
PathCheckStatus is_path_regular_file(const std::string& path); 

/// @brief Check if the given path 1.) exists and 2.) is a directory 
PathCheckStatus is_path_directory(const std::string& path); 

/// @brief Remove a file from disk. throws exception if it doesn't succeed (unless non-null status is passed, then error message is given to 'status' and no expcetion is thrown). 
/// @param path path to file to be removed
  void remove_file_from_disk(const std::string& path, PathCheckStatus* status=nullptr);

/// @brief Based on ifarm's slurm configuration
/// @return path to the root directory of this slurm job's /scratch disk allocation. if it detects we are **not** in a slurm job / array, or the directory is missing, then throws a std::runtime_exception. 
std::string get_slurm_scratch_directory(PathCheckStatus* status=nullptr);

/// @brief get the size of the file, in KB.
/// @param path path to file
/// @param status status object. if the status object is not null, and the size check fails, then the status will be populated with the error message. If the check fails, and the status is null, then an exception is thrown 
double get_file_size_KB(const std::string& path, PathCheckStatus* status=nullptr);
  
/// @brief Attempts to fit a gaussian to a histogram, with a constant background
/// @param hist histogram to fit
/// @param radius radius of fit
/// @param center center of gaus, to be optimized by fit (value of arg. overwritten)
/// @param sigma sigma of gaus, to be optimized by fit (value of arg. overwritten) 
/// @param do_draw if true, then the fit will be drawn 
void fit_gaus_to_hist( 
    TH1 *hist, 
    double radius, 
    double &center, 
    double &sigma, 
    bool do_draw=false
); 


/// @brief Raise value 'x' to integer power 'N'
/// @tparam N power. only valid for N = 2,3,4.  
/// @param x argument  
/// @return x^N 
template <int N> double intpow(double x); 


constexpr double kNaN = std::numeric_limits<double>::quiet_NaN(); 

//null integer value
constexpr int kNaN_int = -99999999;


}
}



#endif
