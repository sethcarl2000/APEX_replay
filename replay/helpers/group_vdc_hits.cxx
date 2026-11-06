
// APEX headers
#include <APEX/replay/helpers.h>
#include <APEX/VDC/Hit.h> 
#include <APEX/VDC/HitGroup.h> 
#include <APEX/EventHandler.h> 
#include <APEX/run_parameters.h>
#include <APEX/units.h>
// ROOT headers
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx> 
// stdlib headers
#include <vector> 
#include <cmath> 

namespace APEX
{
namespace replay
{
namespace helpers
{

//#define DEBUG_GROUP

/// @brief Form 'VDC::HitGroup's from groups of VDC hits
/// @param event Apex S2-event handler
/// @param plane the VDC plane in question, [0,..,3].
/// @param h_wire a collection of all VDC wire numbers
/// @param h_rawtime a collection of all VDC wire rawtimes 
/// @return a vector of all valid VDC hit groups  
template<VDC::Plane P> void group_vdc_hits( 
  S2Event& evt, 
  const ROOT::RVec<double>& wire,
  const ROOT::RVec<double>& rawtime 
) 
{
  constexpr int plane = VDC::number<P>(); 

#ifdef DEBUG_GROUP
    const std::string plane_names[] = {"u1", "v1", "u2", "v2"}; 
    printf(
      "<%s>: grouping VDC hits for %s, VDC plane %s\n",
      __func__, (evt.ActiveArm() ? "RHRS" : "LHRS"), plane_names[plane].c_str() 
    );
#endif        

  //return no groups if there aren't enough hits to make any 
  if (rawtime.size() < run_parameters::kGroup_min_hits) return;

  const bool is_RHRS = evt.Is_RHRS();  
  const double t_S2 = evt.S2Time(); 
  
  ROOT::RVec<VDC::HitGroup<P>> group_vec{}; 

  //make a buffer of all valid hits
  const auto nhits_total = wire.size(); 
  std::vector<VDC::Hit<P>> buffer; buffer.reserve(nhits_total); 

  // make buffer of all valid VDC hits
  for (std::size_t i=0; i<nhits_total; i++) {

    VDC::Hit<P> hit(is_RHRS, static_cast<int>(wire[i]), rawtime[i], t_S2); 

    if (auto t = hit.Time(); 
      t > run_parameters::kVDC_min_realTime && 
      t < run_parameters::kVDC_max_realTime) 
    {
      buffer.emplace_back( std::move(hit) ); 
    }
  }

  //quit if there's not enough hits to form at least 1 hit. 
  if (buffer.size() < run_parameters::kGroup_min_hits) return; 

  //now, we can make groups. 
  int prev_wire = buffer.front().WireNum(); 

  std::size_t ind_begin=0; 
  std::size_t ind_end=ind_begin; 

  //add a new group by copying events from the current buffer 
  auto add_group = [&](){ 
    auto hits = std::span<const VDC::Hit<P>>{buffer.begin()+ind_begin, buffer.begin()+ind_end}; 
    evt.emplace<VDC::HitGroup<P>>(hits);
  };

  //keep iterating until we've gone through all hits
  for (size_t i=0; i<buffer.size(); i++) {
     
    const auto& hit = buffer[i]; 

#ifdef DEBUG_GROUP 
    printf(
      "plane %i, hit %zu/%zu\n"
      "   wire:       %i\n"
      "   time:       %.1f ns\n"
      "   raw time:   %.f\n",
      plane, ind_end, buffer.size(),
      hit.WireNum(), 
      hit.Time()/units::ns, 
      hit.Rawtime()
    );
#endif        

    int gap = i==ind_begin ? 0 : hit.WireNum() - prev_wire -1; 
    if (gap > run_parameters::kGroup_max_gap) { // if the gap is too big to be part of one group 

      //does this cluster already have enough hits? 
      size_t n_hits_in_group = ind_end - ind_begin + 1; 
      if (n_hits_in_group >= run_parameters::kGroup_min_hits) {

        add_group(); 
      }     
      
      //otherwise, make a new group anyway. bring the 'begin' ptr to up to meet this hit (it's potentially the start of a new group)
      ind_begin = i;  

    }//if (span/gap too big) 

    //get ready to re-start this loop with the next hisdfst
    ind_end = i;
    prev_wire = hit.WireNum(); 
  }

  //We've reached the end, but we still might have enough hits to make a new group.  
  size_t n_hits_in_group = ind_end - ind_begin + 1; 
  if (n_hits_in_group >= run_parameters::kGroup_min_hits) {

    add_group(); 
  }     

}; 
//////////////////////////////////////////////////////////////////////////////
// explicit template instantiations 
template void group_vdc_hits<VDC::Plane::U1>( 
  S2Event&, 
  const ROOT::RVec<double>&,
  const ROOT::RVec<double>& 
);
template void group_vdc_hits<VDC::Plane::V1>( 
  S2Event&, 
  const ROOT::RVec<double>&,
  const ROOT::RVec<double>& 
);
template void group_vdc_hits<VDC::Plane::U2>( 
  S2Event&, 
  const ROOT::RVec<double>&,
  const ROOT::RVec<double>& 
);
template void group_vdc_hits<VDC::Plane::V2>( 
  S2Event&, 
  const ROOT::RVec<double>&,
  const ROOT::RVec<double>& 
);


}
}
}