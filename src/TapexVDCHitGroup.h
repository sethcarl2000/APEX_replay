/***
 * @class TapexVDCHitGroup.h
 * 
 * @brief Class which contains one VDC wire TDC hit group   
 * 
 * This class contains and processes a 'hit group', which is a loose collection of wire-hits
 * which are in a certain proximity to one another. A hit group may or may not eventually be converted
 * to a 'legitimage' cluster. But, starting out, it is not assumed that one (or any) clusters exist 
 * within the group 
 * 
 ***/

#ifndef TapexVDCHitGroup_H 
#define TapexVDCHitGroup_H

#include <TObject.h> 
#include <vector> 
#include "TapexVDCHit.h"

class TapexVDCHitGroup : public TObject { 
  
 public:   
  TapexVDCHitGroup(int plane=-1) : fPlane{plane} {};
  
  ~TapexVDCHitGroup() {}; 
  
  void AddHit(const TapexVDCHit& hit) { fHits.push_back(hit); } 
  
  void AddHit( int wire, double rawtime ) { fHits.emplace_back(fPlane,wire,rawtime); }; 
  
  unsigned int Nhits() const { return fHits.size(); }

  double WirePos( unsigned int h ) const; 
  int    WireNum( unsigned int h ) const; 
  double    Time( unsigned int h ) const; 
  
  TapexVDCHit& GetHit( unsigned int h );

  int    FirstWire() const { return fHits.front().wNum(); } 
  
  double LoEdge()    const { return fHits.back().wPos(); } 
  double HiEdge()    const { return fHits.front().wPos(); } 
  
  double Span()     const { return HiEdge()-LoEdge(); }
  
  bool IsRightArm() const { return fHits.at(0).IsRightArm(); }
  
  double W()        const { return fHits.at(0).W(); }
  
  
private: 
  std::vector<TapexVDCHit> fHits; 
  int fPlane; 
  
  const double kNull_double=-1e30; 
  const int    kNull_int   =-999; 

  ClassDef(TapexVDCHitGroup,0); 
};

#endif 