
#include <APEX/utils/ErrorHandler.h>
// stdlib
#include <cstdio>
#include <iostream> 
#include <sstream> 
#include <cmath> 


//nmmmmmmmmmmmmmioooooooooooooooooooooooooooooooooooooooooo 
// - muon's comment 28 aug 26


namespace APEX
{
namespace utils
{

//__________________________________________________________________________________________________________________________
ErrorHandler& ErrorHandler::Instance() 
{
    static ErrorHandler instance; 
    return instance;
}
//__________________________________________________________________________________________________________________________
ErrorHandler::ErrorHandler()
{
    std::cout << "in ctor" << std::endl; 

    fLevels = std::vector<level_and_name_t>{
        { kPrint,   "" },
        { kInfo,    "Info" },
        { kWarning, "Warning" },
        { kError,   "Error" },
        { kBreak,   "Error (Break)" },
        { kFatal,   "Fatal Error" }
    };

    fStartPoint = std::chrono::system_clock::now(); 

    //signal that the clock has started
    Info(__func__, "Error handler has been initialized; message timestamps counting from this initialization point."); 

    //tell ROOT to call this function when throwing an error. 
    SetErrorHandler( &ProcessErrorReport );
}
//__________________________________________________________________________________________________________________________
void ProcessErrorReport(int level, bool abort, const char* location, const char* msg) 
{
    
    static std::mutex mut_cout, mut_cerr; 

    auto& handler = ErrorHandler::Instance(); 
    //decide whether or not to report this error
    const std::string* name{nullptr}; 

    //quit if this error isn't important enough
    if (level < handler.GetMinPrintLevel()) return; 

    for (const auto& level_name : handler.GetLevelPrefixes()) {
        if (level <= level_name.level) { name = &level_name.name; break; } 
    }

    std::ostringstream oss; 

    //print the timestamp 
    const auto now = std::chrono::system_clock::now(); 
    oss << "[" << handler.FormattedTimestamp(now) << "] "; 

    //print the error level, location and message 
    oss << *name << " from <" << location << ">: " << msg; 

    const bool print_to_stdout = handler.PrintToStdout(); 

    if (print_to_stdout) {
        mut_cout.lock(); 
        std::cout << oss.str() << "\n"; 
        mut_cout.unlock(); 
    }

    if (level >= kWarning) {
        mut_cerr.lock(); 
        std::cout << oss.str() << "\n"; 
        mut_cerr.unlock(); 
    }

    // in these cases, abort execution of the program. 
    if (level >= kFatal) { std::abort(); } else { if (level >= kBreak) std::exit(1); }
    return; 

}
//__________________________________________________________________________________________________________________________
std::string ErrorHandler::FormattedTimestamp(const std::chrono::system_clock::time_point& time) const
{
    std::chrono::duration<double, std::milli> duration{ time - fStartPoint }; 

    //duration from start point in ms
    double seconds = duration.count() / 1e3; 

    int hours   = (int)std::floor(seconds / 3600); 
    int minutes = (int)std::floor(seconds / 60) - hours*60; 
    seconds     = seconds - (hours*3600.) - (minutes*60.); 

    char buff[15]; 
    
    std::sprintf(buff, "%02i:%02i:%07.4f", hours, minutes, seconds); 

    return std::string{buff}; 
} 

//__________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________


}
}