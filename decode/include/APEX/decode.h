#ifndef APEX_decode_h
#define APEX_decode_h

// ROOT headers
#include <RtypesCore.h> 
// stdlib headers
#include <string> 

namespace APEX
{
namespace decode
{

//initialze detector objects needed to decode APEX's CODA files, and hand them over to 'gHaApps'. 
void init_decoders();

/// @brief Attempt to decode CODA file, and output result in ROOT file
/// @param path_input absolute path to CODA file. 
/// @param path_output absolute path to output .root file. 
/// @param first_event first physics event to process (0 is the first physics event in **this** CODA file, **not** the first physics event in the run!!!!) 
/// @param last_event last physics event to process, indexed the same way as 'first_event' 
/// @param path_odef path to output variable definitions
void process_coda_file(const std::string& path_input, const std::string& path_output, Long64_t first_event=0, Long64_t last_event=1e7, const std::string& path_odef="outDefs/full_replay.odef");

/// @brief Count the number of physics events in a CODA file 
/// @param path_coda absolute path to CODA file 
/// @return number of physics triggers counted in CODA file. If there was a processing error, '0' is returned (and a std::runtime_exception is thrown). 
Long64_t count_physics_events(const std::string& path_coda); 

}
}

#endif