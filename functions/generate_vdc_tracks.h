#ifndef generate_vdc_tracks_H
#define generate_vdc_tracks_H

// APEX headers
#include "group_vdc_hits.h"
#include "run_parameters.h"
#include <EventCounter.h> 
// ROOT headers
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
// stdlib headers
#include <vector>
#include <string> 

/// @brief Generates VDC tracks for RHRS / LHRS
/// @param is_RHRS
/// @param node_in input RDF node
/// @param n_pass_1group EventCounter representing the number of events which reconstructed at least 1 group 
ROOT::RDF::RNode generate_vdc_tracks( 
    const bool is_RHRS, 
    ROOT::RDF::RNode inNode, 
    EventCounter &nPass_1group, 
    EventCounter &nPass_1pair, 
    EventCounter &nPass_1rawTrack ) 
{
    using namespace std;            

    //feed this script your node, with generated S2-hits, ready for tracking data. 

    //it'll spit out refined tracks for the plane you tell it to. 

    string arm = is_RHRS ? "R" : "L" ; 

    vector<string> plane_name = { "u1", "v1", "u2", "v2" }; 

    vector<string> branch_rawtime; 
    vector<string> branch_wire; 

    for (int p=0; p<4; p++) { 
        
        string rawtime = arm+".vdc."+plane_name[p]+".rawtime"; 
        
        string wire    = arm+".vdc."+plane_name[p]+".wire"; 
        
        branch_rawtime.push_back( rawtime ); 
        branch_wire   .push_back( wire );   
    }


        
    
}



#endif