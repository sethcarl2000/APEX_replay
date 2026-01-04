#ifndef Podd_THRSTrack_h_
#define Podd_THRSTrack_h_

//////////////////////////////////////////////////////////////////////////
//
// TOpticsModel
//
//////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TObject.h"

using namespace std; 

using RVecD = ROOT::RVec<double>;

class THRSTrack : public TObject {
  
 public:
  
  enum EVDCCoordType { kUVW=0, kDetector, kTransport, kFocalPlane }; 
  
  TVDCTrack() {};  
  ~TVDCTrack();
  
  void Set_Intercept(const int p,const double intercept);  
  void Set_Intercepts(const double intercept[4]); 
  
  double Get_Intercept(const int p) const;
  void   Get_Intercepts(double (&intercepts)[4]); 
  
  double Get_x   (const EVDCCoordType coord) const; 
  double Get_y   (const EVDCCoordType coord) const; 
  double Get_dxdz(const EVDCCoordType coord) const; 
  double Get_dydz(const EVDCCoordType coord) const; 
  
 private:

  //I want to make these as light-weight as possible, so they take up little data
  // in the full-replay.
  double fIntercept[4];
  //timing of track (at S2 plane, specifically)
  double fTime; 
    
  ClassDef(THRSTrack,0);   
};

#endif
