#ifndef generate_vdc_tracks_H
#define generate_vdc_tracks_H

// APEX headers
#include "group_vdc_hits.h"
#include "run_parameters.h"
#include <EventCounter.h> 
// ROOT headers
#include <ROOT/RVec.hxx>
// stdlib headers
#include <vector>
#include <string> 

/// @brief Generates VDC tracks for RHRS / LHRS
/// @param evt Apex event handler
/// @param node_in input RDF node
/// @param n_pass_1group EventCounter representing the number of events which reconstructed at least 1 group 
ROOT::RDF::RNode generate_vdc_tracks(
    TapexEventHandler* evt, 
    ROOT::RDF::RNode node_in, 
    EventCounter& n_pass_1group  
)
{
    using namespace std; 
    using rvecd = ROOT::VecOps::RVec<double>; 

    std::vector<ROOT::RDF::RNode> nodes{node_in}; 
    
    const bool is_RHRS = evt->ActiveArm(); 

    string arm = is_RHRS ? "R" : "L"; 
    const vector<string> plane_name{"U1","V1","U2","V2"}; 
    
    vector<string> branch_rawtime, branch_wire; 

    //group hits
    for (int p=0; p<4; p++) { 
        
        string br_rawtime = arm+".vdc."+plane_name[p]+".rawtime"; 
        
        string br_wire    = arm+".vdc."+plane_name[p]+".wire"; 

        auto new_node = nodes.back()

            .Define(Form("%s_%s_groups",arm.c_str(),plane_name[p].c_str()), [p, evt](const rvecd& rawtime, const rvecd& wire)
            {
                return group_vdc_hits(evt, p, rawtime, wire); 
            }, {br_rawtime.c_str(), br_wire.c_str()}); 

        nodes.push_back(new_node); 
    }
    
    auto output_node = nodes.back(); 

    n_pass_1group = output_node.Count(); 

    return output_node; 
}



#endif