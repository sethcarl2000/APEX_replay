#ifndef EventCounter_H
#define EventCounter_H

#include <ROOT/RResultPtr.hxx> 

/// @brief simple rename for a RResultPtr which can be used to count the number of events passing a given analysis step 
using EventCounter = ULong64_t; 
using EventCounter_RPtr = ROOT::RDF::RResultPtr<ULong64_t>; 

#endif