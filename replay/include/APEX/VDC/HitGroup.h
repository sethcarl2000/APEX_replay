#ifndef APEX_VDC_HitGroup_H
#define APEX_VDC_HitGroup_H

#include <APEX/VDC.h> 
#include <APEX/VDC/Hit.h> 

namespace APEX 
{
namespace VDC
{

class HitGroup {
private: 
  bool fis_RHRS; 
  std::vector<Hit> fHits; 
  int fPlane; 


public:   
  
  HitGroup(bool is_RHRS, int plane) : fis_RHRS{is_RHRS}, fPlane{plane} {}; 
  
  void AddHit( const Hit& hit ) { fHits.push_back(hit); } 
  
  void AddHit( int wire, double rawtime, double T_offset=0.) {
    fHits.emplace_back(fis_RHRS, fPlane, wire, rawtime, T_offset);  
  }; 
  
  int Nhits() const { return fHits.size(); }
    
  double WirePos( unsigned int h ) const; 
  int    WireNum( unsigned int h ) const; 
  double    Time( unsigned int h ) const; 
  
  Hit& GetHit( unsigned int h ) { return fHits.at(h); } 
  
  int    FirstWire()  const; 
  
  double LoEdge()   const; 
  double HiEdge()   const; 
  
  double Span()     const { return HiEdge()-LoEdge(); }
  
  bool IsRightArm() const { return fHits.front().Is_RHRS(); }
  
  double W()        const { return fHits.front().W(); }

  //delete all hits
  void Clear() { fHits.clear(); }
  
private: 
  
}; 

}
}
 
#endif 