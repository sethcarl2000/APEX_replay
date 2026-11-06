#ifndef APEX_VDC_HitGroup_H
#define APEX_VDC_HitGroup_H

#include <APEX/ObjCollection.h>
#include <APEX/VDC.h> 
#include <APEX/VDC/Hit.h> 
// stdlib 
#include <span> 

namespace APEX 
{

namespace VDC
{

template<VDC::Plane P> class HitGroup {
private: 

  friend class ObjCollection<HitGroup<P>>; 
  object_id id{0}; 

  bool fis_RHRS; 
  std::vector<Hit<P>> fHits{}; 


public:   

  using iterator = std::vector<Hit<P>>::iterator; 
  using const_iterator = std::vector<Hit<P>>::const_iterator; 

  HitGroup(bool is_RHRS) : fis_RHRS{is_RHRS} {}; 
  
  // fill out a hit group from a span hits (make copies of all underlying hits). 
  HitGroup(std::span<const Hit<P>> hits); 

  template<typename... Args> void AddHit(Args&&... args) { 
    fHits.emplace_back(std::forward<Args>(args)...); 
  } 
  
  void AddHit( int wire, double rawtime, double T_offset=0.) {
    fHits.emplace_back(fis_RHRS, wire, rawtime, T_offset);  
  }; 
  
  int Nhits() const { return fHits.size(); }
    
  iterator begin() { return fHits.begin(); }
  iterator end()   { return fHits.end(); }

  const_iterator begin() const { return fHits.begin(); }
  const_iterator end()   const { return fHits.end(); }

  const_iterator cbegin() const { return fHits.cbegin(); }
  const_iterator cend()   const { return fHits.cend(); }


  double WirePos( unsigned int h ) const; 
  int    WireNum( unsigned int h ) const; 
  double    Time( unsigned int h ) const; 
  
  Hit<P>& GetHit( unsigned int h ) { return fHits.at(h); } 

  // First wire in the hit-group 
  int    FirstWire() const { return fHits.front().WireNum(); }
  // Last wire in the hit-group  
  int    LastWire()  const { return fHits.back().WireNum(); }
  
  // Low edge
  double LoEdge()    const { return fHits.back().WirePos(); }
   
  // hi edge
  double HiEdge()    const { return fHits.front().WirePos(); }
  
  double Span()     const { return HiEdge()-LoEdge(); }
  
  bool Is_RHRS() const { return fis_RHRS; }
  
  double W()        const { return fHits.front().W(); }

  //delete all hits
  void Clear() { fHits.clear(); }

  // Get unique ID
  object_id ID() const { return id; }
  bool operator==(const HitGroup<P>& rhs) const { return rhs.id == id; }

private: 
  
}; 
// explicit instantiations 
template class HitGroup<Plane::U1>; 
template class HitGroup<Plane::V1>; 
template class HitGroup<Plane::U2>; 
template class HitGroup<Plane::V2>; 

}
}
 
#endif 