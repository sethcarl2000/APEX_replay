#include "TObject.h"

class TvdcHit; 

class THitGroup : public TObject { 
 private: 
  std::vector<TvdcHit*> fHits; 
  int fPlane;
  
 public: 
  THitGroup(); 
  ~THitGroup(); 
  
  THitGroup(int plane); 
  
  int AddHit(TvdcHit *hit); 
  int AddHit(double wire,double rawTime); 
  
  int Nhits();
  int FirstWire(); 
  double LoEdge(); 
  double HiEdge(); 
  double Span(); 
  double WirePos( int h ); 
  double Time( int h ); 
  int WireNum( int h ); 
  
  ClassDef(THitGroup,1);
}; 
  
  
