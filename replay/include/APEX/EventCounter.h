#ifndef APEX_EventCounter_H
#define APEX_EventCounter_H

#include <ROOT/RResultPtr.hxx> 

namespace APEX
{
    
/// @brief simple rename for a RResultPtr which can be used to count the number of events passing a given analysis step 
using EventCounter = ULong64_t; 
using EventCounter_RPtr = ROOT::RDF::RResultPtr<ULong64_t>; 

}

#endif