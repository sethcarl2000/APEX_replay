#ifndef APEX_replay_h
#define APEX_replay_h

// ROOT headers
#include <RtypesCore.h> 
#include <ROOT/RDataFrame.hxx>
// stdlib headers
#include <string> 

namespace APEX
{
namespace replay
{

// determines both arms. 
enum EArmMode : int { kRHRS=1<<0, kLHRS=1<<1 }; 
constexpr int kBothArms = kRHRS | kLHRS; 

/// @brief simple rename for a RResultPtr which can be used to count the number of events passing a given analysis step 
using EventCounter = ULong64_t; 
using EventCounter_RPtr = ROOT::RDF::RResultPtr<EventCounter>; 

}
}

#endif