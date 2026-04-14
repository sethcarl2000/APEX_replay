
#include "TvdcHit.h"
#include "THitGroup.h"
#include "../def_apex.h"
#include "../function_lib.C"

ClassImp(THitGroup); 

THitGroup::THitGroup() {/*THitGroup*/}; 
THitGroup::THitGroup(int plane) { 

  fPlane = plane; 
}
int    THitGroup::Nhits() { return fHits.size(); }
int    THitGroup::AddHit( TvdcHit *hit ) {
  
  fHits.push_back( hit ); 
  return Nhits(); 
}
int    THitGroup::AddHit( double wire, double rawTime ) { 
  
  auto hit = new TvdcHit(fPlane,wire,rawTime); 
  AddHit( hit ); 
  return Nhits();  
}
double THitGroup::WirePos( int h ) {
   
  if (h>=Nhits()) return -1e30;   
  return fHits.at(h)->wPos(); 
}
int    THitGroup::WireNum( int h ) { 
  
  if (h>=Nhits()) return -99999;   
  return fHits.at(h)->wNum(); 
} 
double THitGroup::Time( int h ) { 
  
  if (h>=Nhits()) return -1e30;   
  return fHits.at(h)->Time(); 
} 
double THitGroup::LoEdge() { return WirePos(Nhits()-1); }
double THitGroup::HiEdge() { return WirePos(0); }
int    THitGroup::FirstWire() { return WireNum(0); }
double THitGroup::Span()   { return HiEdge()-LoEdge(); }

THitGroup::~THitGroup() {/*default destructor*/
  
  //delete constituent hits
  for (int h=0; h<Nhits(); h++) fHits.at(h)->~TvdcHit();   
}

