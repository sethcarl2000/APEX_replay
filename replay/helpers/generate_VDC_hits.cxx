
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

template<typename T> using Rv = ROOT::VecOps::RVec<T>; 

//#define DEBUG_GROUP

/// @brief generate a list of all VDC hits for a given arm & plane
/// @param arm_mode which arm is currently in use
/// @param p VDC plane index [0-3]
/// @param h_wire list of VDC wire indices [0-367]
/// @param h_rawtime raw TDC time for each wire
/// @return an RVec<VDC::Hit> of all hits with reasonable times. Sorted by ascending wire-number. 
ROOT::VecOps::RVec<VDC::Hit> generate_VDC_hits(
    EArmMode arm_mode,
    int   plane, 
    const Rv<double>& h_wire,
    const Rv<double>& h_rawtime 
)
{
#ifdef DEBUG_GROUP
  using namespace units; 
  std::cout << "DEBUG: Vdc valid real-time range: ["
  << run_parameters::kVDC_min_realTime/ns << ", " << run_parameters::kVDC_max_realTime/ns << "]" << std::endl;
  int n_validHits(0);
  std::cout << "Is right arm? " << (evt.ActiveArm() ? "true" : "false") << "\n";
  std::cout << " plane " << p << std::endl;
  std::cout << " S2 hit time: " << evt.GetS2Hit(evt.ActiveArm())->Time()*1e9 << std::endl; 
#endif
  
  //hits are stored by descending wire-number by default. let's flip the order. 
  const std::size_t n_raw = h_wire.size(); 
  
  Rv<VDC::Hit> hits; hits.reserve(n_raw); 

  const bool is_RHRS = (arm_mode == kRHRS); 

  double default_offset = is_RHRS ? run_parameters::RHRS_S2_reftime : run_parameters::LHRS_S2_reftime; 

  const double t_min = is_RHRS ? run_parameters::RHRS_VDC_min_realtime : run_parameters::LHRS_VDC_min_realtime; 
  const double t_max = is_RHRS ? run_parameters::RHRS_VDC_max_realtime : run_parameters::LHRS_VDC_max_realtime; 
  
  for (std::size_t i=0; i<n_raw; i++) {

    VDC::Hit hit(
      is_RHRS, 
      plane, 
      static_cast<int>(h_wire[i]), 
      h_rawtime[i], 
      default_offset
    ); 

    //if this hit's time is reasonable, then keep it. 
    if (double t = hit.Time(); (t > t_min && t < t_max))
      hits.emplace_back(hit); 
  }

  return hits; 
}; 

}
}
}