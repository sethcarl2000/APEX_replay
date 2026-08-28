#ifndef ApexVDCHitGroup_H
#define ApexVDCHitGroup_H

#include <ApexVDC.h> 
#include <ApexVDCHit.h> 

namespace ApexVDC 
{

class HitGroup {

 public:   
  
  HitGroup(int plane=-1) : fPlane{plane} {}; 
  
  virtual ~HitGroup() {}; 

  void AddHit( const Hit& hit ) { fHits.push_back(hit); } 
  
  void AddHit( double wire, double rawTime ); 
  
  int Nhits() const { return fHits.size(); }
    
  double WirePos( unsigned int h ) const; 
  int    WireNum( unsigned int h ) const; 
  double    Time( unsigned int h ) const; 
  
  Hit& GetHit( unsigned int h ) { return fHits.at(h); } 
  
  int    FirstWire()  const; 
  
  double LoEdge()   const; 
  double HiEdge()   const; 
  
  double Span()     const { return HiEdge()-LoEdge(); }
  
  bool IsRightArm() const { return fHits.front().IsRightArm(); }
  
  double W()        const { return fHits.front().W(); }

  //delete all hits
  void Clear() { fHits.clear(); }
  
private: 
  std::vector<Hit> fHits; 
  int fPlane; 

}; 

}; 
#endif 