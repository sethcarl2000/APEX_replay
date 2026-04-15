#ifndef group_vdc_hits_H
#define group_vdc_hits_H

// APEX headers
#include <ApexVDCHit.h> 
#include <ApexVDCHitGroup.h> 
#include <TapexEventHandler.h> 
#include "../run_parameters.h"
// ROOT headers
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx> 
// stdlib headers
#include <vector> 
#include <cmath> 

/// @brief Form 'ApexVDC::HitGroup's from groups of VDC hits
/// @param evt Apex event handler
/// @param p the VDC plane in question, [0,..,3].
/// @param h_rawtime a collection of all VDC wire rawtimes 
/// @param h_wire a collection of all VDC wire numbers
/// @return a vector of all valid VDC hit groups  
ROOT::RVec<ApexVDC::HitGroup> group_vdc_hits ( 
                                            const TapexEventHandler& evt, 
                                            int   p, 
                                            const ROOT::RVec<double>& h_rawtime, 
                                            const ROOT::RVec<double>& h_wire
                                            ) 
{
#if DEBUG
  cout << "DEBUG: Vdc valid real-time range: ["
  << VDC_min_realTime << ", " << VDC_max_realTime << "]" << endl;
  int n_validHits(0);
  cout << "Is right arm? " << (evt->ActiveArm() ? "true" : "false") << "\n";
  cout << " plane " << p << endl;
  cout << " S2 hit time: " << evt->GetS2Hit()->Time()*1e9 << endl; 
#endif
  
  //return no groups if there aren't enough hits to make any 
  if (h_rawtime.size() < run_parameters::kGroup_min_hits) return {};

  int n_hits = h_rawtime.size(); 
  
  ROOT::RVec<ApexVDC::HitGroup> group_vec{}; 

  //cBeg->push_back((int)wire[0]);   
  ApexVDC::HitGroup group(p); 

  int w_prev = -1; 
  
  //the clust_points vector stores each hit as a 2-vector;
  // the [0] element is the wirePos,
  // the [1] element is the T_TDC (corrected tdc time of that hit)
  for (int h=0; h<n_hits; h++) { 
    
    ApexVDC::Hit hit( p, (int)std::round(h_wire[h]), h_rawtime[h], &evt );  
    
    //vdc timing cut, these times won't ever be useful for a coinc track 
    if (hit.Time() > run_parameters::kVDC_max_realTime || 
        hit.Time() < run_parameters::kVDC_min_realTime ) {
      #if DEBUG
        cout << "killed! " << endl; 
      #endif 
      continue; 
    }
    
#if DEBUG
    cout << "kept!" << endl; 
    n_validHits++; 
#endif       
    //if wPrev =-1, then no wires have yet been added
    int gap = w_prev<0 ? 0 : hit.wNum() - w_prev  -1; 

    if (gap > run_parameters::kGroup_max_gap) { // for adjacent wires, 

      //does this cluster already have enough hits? 
      if (group.Nhits() >= run_parameters::kGroup_min_hits) { 

        group_vec.push_back( group );  
      } 

      //otherwise, make a new group anyway.
      group.Clear();  

    }//if (span/gap too big) 
      
    //if span is NOT too big, then add this hit to the current cluster

    group.AddHit( hit );   
    w_prev = hit.wNum(); 
    
  }//for (int i=0; h<nHtis; h++) 

  //the last hit is reached, and it's not over-span
  if (group.Nhits() >= run_parameters::kGroup_min_hits)  { 
    
    group_vec.push_back( group ); 

  } 
  
#if DEBUG 
  cout << TString::Format( "group_hits() -  total hits:%3i, valid-hits:%3i, groups made:%2i",
          (int)h_wire.size(), 
          n_validHits, 
          (int)groupVec.size() ) << endl; 
#endif
  
  return group_vec;       
}; 

#endif 