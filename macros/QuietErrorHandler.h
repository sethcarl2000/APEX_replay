#ifndef QuietErrorHandler_h
#define QuietErrorHandler_h

#include <TError.h>
#include <iostream>
#include <cstdlib>
#include <string>
#include <map> 

//meyers singleton class for handling errors
class QuietErrorHandler {
 private:

  struct level_and_name_t { int level; std::string name; }; 
  
  static std::vector<level_and_name_t> fLevels;
  
  static int fMinPrintLevel; 

  //make constructor / destructor private 
  QuietErrorHandler() {
    SetErrorHandler( &QuietErrorHandler::ErrorFcn );
  }
  ~QuietErrorHandler() {}; 
  
 public: 

  //delete copy construction & assignment operator
  QuietErrorHandler(const QuietErrorHandler&) = delete;
  QuietErrorHandler& operator=(const QuietErrorHandler&) = delete;

  static QuietErrorHandler& Instance() {
    static QuietErrorHandler inst;
    return inst;
  }
  
  static void ErrorFcn(int level, bool abort, const char* location, const char* msg) {
    //decide whether or not to report this error
    std::string name = ""; 
    
    for (const auto& level_name : fLevels) {
      if (level <= level_name.level) { name = level_name.name; break; } 
    }
    
    if (level >= fMinPrintLevel)
      std::cerr << name <<
	" in <"<< (location != nullptr ? location : "unspecified location") <<
	">: " << msg << "\n";
    
    if (abort) std::abort();
  }

  void SetMinPrintLevel(int level) { fMinPrintLevel=level; }

};

int QuietErrorHandler::fMinPrintLevel = 0; 

std::vector<QuietErrorHandler::level_and_name_t> QuietErrorHandler::fLevels = std::vector<QuietErrorHandler::level_and_name_t>{
  { kPrint, "" },
  { kInfo, "Info" },
  { kWarning, "Warning" },
  { kError, "Error" },
  { kBreak, "\n *** Break ***\n" },
  { kFatal, "\n *** Break ***\nFatal Error" }
};
  
    
#endif

