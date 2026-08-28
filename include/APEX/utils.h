#ifndef APEX_utils_h
#define APEX_utils_h


#include <APEX/EventCounter.h>
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

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN(); 


class UniqueIDGenerator { 
  
 public: 
  /*** 
   *    Generates Unique ID each time it's called in thread-safe manner
   * 
   ***/

  UniqueIDGenerator() : fCounter{0} {};
  virtual ~UniqueIDGenerator() {}; 
  
  /// @return ID number which is unique across all threads 
  Long64_t GetID() { return fCounter.fetch_add(1, std::memory_order_relaxed); }

  //delete copy / move constructors
  UniqueIDGenerator(const UniqueIDGenerator&) = delete;
  UniqueIDGenerator& operator=(const UniqueIDGenerator&) = delete;
  
  
private:

  std::atomic<EventCounter> fCounter; 
  
}; 


}
}



#endif