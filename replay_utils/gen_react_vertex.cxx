#ifndef gen_react_vertex_H
#define gen_react_vertex_H 

#include <replay_utils/interface.h>
#include <TapexReactVertex.h> 
#include "../include/RDFNodeAccumulator.h"
#include <string> 


namespace replay_utils
{

/// @brief Generates react vertex
/// @param rna RDFNodeAccumulator to add branches to 
/// @param vtx const ptr to TapexReactVertex obj (arm specific!)
void gen_react_vertex(RDFNodeAccumulator& rna, TapexReactVertex* vtx)
{
    const bool is_RHRS = vtx->IsRHRS(); 

    //names for BPM branches
    
    const std::string raster_name = is_RHRS ? "Rrb.Raster2" : "Lrb.Raster2"; 
    const std::string rb_name = is_RHRS ? "Rrb" : "Lrb"; 

    rna.Define((is_RHRS ?"R_position_vtx":"L_position_vtx"), [vtx](double current_x, double current_y, double bpma_y, double bpmb_y) 
        {
            return vtx->Compute_reactVertex( current_x, current_y, (bpma_y+bpmb_y)/2 );

        }, {raster_name+".rawcur.x", raster_name+".rawcur.y", rb_name+".BPMA.y", rb_name+".BPMB.y"});

}

}; 

#endif