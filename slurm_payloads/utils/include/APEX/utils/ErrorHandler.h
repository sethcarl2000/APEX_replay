#ifndef APEX_utils_ErrorHandler_h
#define APEX_utils_ErrorHandler_h

// ROOT
#include <TError.h>
// stdlib
#include <iostream>
#include <cstdlib>
#include <string>
#include <map> 
#include <chrono> 
#include <mutex> 


namespace APEX
{
namespace utils
{

//meyers singleton class for handling errors
class ErrorHandler {
public: 
    struct level_and_name_t { int level; std::string name; }; 
 
private:

  //this is initialized upon creation of this error handler
  std::chrono::time_point<std::chrono::system_clock> fStartPoint; 
  
  std::vector<level_and_name_t> fLevels;

  /// @brief true == print error messages to stdout as well
  bool fPrintToStdout{true};

  int fMinPrintLevel{0}; 

  //make constructor private 
  ErrorHandler(); 
  ~ErrorHandler() = default; 
 
public: 

  /// @brief Get reference to the error handler instance
  static ErrorHandler& Instance(); 

  /// @brief Initialize the error handler, and connect it ROOT's TError interface. 
  static void Connect() { Instance(); }

  /// @return whether to also print to stdout 
  bool PrintToStdout() const { return fPrintToStdout; }

  /// @brief if 'true' then also print to stdout
  void SetPrintToStdout(bool val) { fPrintToStdout=val; }

  /// @brief Given a duration, format the time in the following way: [hh:mm:ss.sss]
  /// @param duration a duration 
  /// @return 
  std::string FormattedTimestamp(const std::chrono::system_clock::time_point& time) const;
  
  //delete copy construction & assignment operator
  ErrorHandler(const ErrorHandler&) = delete;
  ErrorHandler* operator=(const ErrorHandler&) = delete;


  const std::vector<level_and_name_t>& GetLevelPrefixes() const { return fLevels; }

  int GetMinPrintLevel() const { return fMinPrintLevel; }

  void SetMinPrintLevel(int level) { fMinPrintLevel=level; }

};

static void ProcessErrorReport(int level, bool abort, const char* location, const char* msg);
  
}
}
    
#endif

