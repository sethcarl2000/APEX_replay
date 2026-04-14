#ifndef ApexVDCTrack_H
#define ApexVDCTrack_H 

//APEX headers 
#include <ApexVDC.h> 
#include <ApexVDCChamberPair.h> 
#include <TapexEventHandler.h> 
#include <math.h> 
//ROOT headers
#include <TVector3.h>

namespace ApexVDC 
{
    
class Track { 
  
 public: 
  // This is the track object, which handles self-refinement, 
  // as well as acting as a container for its constituent points. 
  Track( TapexEventHandler *event=0, 
	     ChamberPair *pLo=0, 
	     ChamberPair *pHi=0 ); 
  
  ~Track() {
    
    //cout << "Track:: deleting track..."<<endl; 
    
    if (fPair_Lo) fPair_Lo->Remove_track( this ); 
    if (fPair_Hi) fPair_Hi->Remove_track( this ); 
  }; //noop
  
  void SetPair_Hi( ChamberPair *pHi );  
  void SetPair_Lo( ChamberPair *pLo );  
  
  ChamberPair* GetPair_Hi() { return fPair_Hi; }
  ChamberPair* GetPair_Lo() { return fPair_Lo; }
  
  void SetEvent( TapexEventHandler *evt ) { fEvent=evt; }
  TapexEventHandler *GetEvent() { return fEvent; }
  
  bool IsRightArm() const { return f_isRightArm; }
  
  void Set_uv_Hi( const double u, const double v ); 
  
  void Set_uv_Lo( const double u, const double v );
  
  double T0() const { return fT0; }
  void   Set_T0(const double T0) { fT0=T0; }
    
  //returns the agreement of this hit with the S2-paddle
  double xParam() const; 
  
  
  //adjust the params by this amount
  void Nudge_params( double nudge[5] ); 
  
  void Set_params( const double vLo, 
		   const double uLo,
		   const double vHi, 
		   const double uHi, 
		   const double T0 ); 
  
  void Set_params( const double params[5] ); 
  
  ApexVDC::HitGroup *GetGroup(int plane) { return fGroup[plane]; }
    
  double Slope_u(); 
  double Slope_v(); 
  
  //per-plane operations
  void   Set_intercept(int plane, double x); 
  double Intercept(int plane)       const;

  //this is for monte-carlo processing, and it avoids needing actual hit data
  void Set_S2int_angles( const double s2x, 
			 const double s2y, 
			 const double theta, 
			 const double phi ); 
  
  int    Nhits(int plane)           const;
  
  double Tau(int plane, int h)      const; 
  double WirePos(int plane, int h)  const; 
  double Slope(int plane)           const; 
  
  double Get_T_model(int plane, double x, int derivative=0) const; 
  
  double ToF(int plane, double x)   const; 
  
  void UpdateTrackInfo();
    
  //overall track parameters
  double FP_x()  const; 
  double FP_y()  const; 
  

  //optics data ~~~~~~~~~~~~~~~~~~~~~~~~~~ (Computation handled by THRS class)
  //track target coordinates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //placeholder functions 
  double Transport_x() const { return FP_x(); } 
  double Transport_y() const { return FP_y(); } 
  double dx_dz()       const { return std::tan(Theta()); } 
  double dy_dz()       const { return std::tan(Phi()); } 
    
  void Set_targetCoords( const double tg_y, 
			 const double tg_theta, 
			 const double tg_phi, 
			 const double tg_dp ) 
  { f_tg_y     =tg_y;
    f_tg_theta =tg_theta; 
    f_tg_phi   =tg_phi; 
    f_tg_dp    =tg_dp; } 
  
  //all angles in radians 
  double Get_tg_y()     const { return f_tg_y; } 
  double Get_tg_theta() const { return f_tg_theta; } 
  double Get_tg_phi()   const { return f_tg_phi; } 
  
  
  //track fp coordinates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void Set_focalPlaneCoords( const double fp_y, 
			     const double fp_theta, 
			     const double fp_phi ) 
  { f_fp_y     =fp_y;
    f_fp_theta =fp_theta; 
    f_fp_phi   =fp_phi; } 
  
  //all angles in radians 
  double Get_fp_y()     const { return f_fp_y; } 
  double Get_fp_theta() const { return f_fp_theta; } 
  double Get_fp_phi()   const { return f_fp_phi; } 
    
  
  double GetTimeAtZ(const double z) const; 
  
  
  double S2_x()  const; 
  double S2_y()  const; 
  
  double Theta() const; 
  double Phi()   const; 
  
  //set/get track errors. 
  void Set_Errors( double errors[5] ); 
  
  double Error_intercept(int plane) const { return fIntercept_ERR[plane]; }
  double Error_T0()                 const { return fT0_ERR; }
  
  //a general figure which is a stand-in for the error of all track intercepts
  double Error_allIntercept()       const; 
    
  
  double Error_Theta()              const { return fTheta_ERR; }
  double Error_Phi()                const { return fPhi_ERR; }
  
  
  //
  void Set_goodPointGroup( int plane, ApexVDC::HitGroup *goodPointGroup ) {
    f_goodPointGroup[plane] = goodPointGroup; 
  }
  
  ApexVDC::HitGroup *Get_goodPointGroup( int plane ) { 
    return f_goodPointGroup[plane]; 
  }
  
  void   Set_Eta(int plane, double eta) { fPlane_Eta[plane]=eta; }
  void   Set_Eta(double eta[4]); 
  double Get_Eta(int plane) const       { return fPlane_Eta[plane]; }
  double Get_Eta()          const; //for all planes
  
  void   Set_RMS(int plane, double RMS) { fPlane_RMS[plane]=RMS; }
  double Get_RMS(int plane) const       { return fPlane_RMS[plane]; }
  double Get_RMS()          const; //for all planes
  
  void Set_nGoodPoints(int plane, int pts) { fPlane_goodPoints[plane]=pts; }
  int  Get_nGoodPoints(int plane) const    { return fPlane_goodPoints[plane]; }
  int  Get_nGoodPoints()          const; //for all planes
    
  
  //this flag is used for monte-carlo testing
  void Set_isGoodTrack(bool isGood) { f_isGoodTrack =isGood; } 
  bool IsGoodTrack() const { return f_isGoodTrack; }
    
  static TVector3 Rotate_uvw_to_xyz( const TVector3& vec ); 
  static TVector3 Rotate_xyz_to_uvw( const TVector3& vec ); 
  
  static void Compute_Theta_Phi( const bool arm, 
				 const double intercepts[4], 
				 double &Theta, 
				 double &Phi ); 
    
  enum ECoordType { kDetector=0, kTransport, kFocalPlane };  

  TVector3 ComputeIntercept_w(const double w) const; 
  
  TVector3 ComputeIntercept_z(const double z) const;   
  
 private: 
  TapexEventHandler *fEvent;
  bool f_isRightArm; 
  
  //does this track have associated hits/hitgroups? if not make sure we don't try 
  // to access them. 
  bool f_hasVDCdata; 
  
  ChamberPair *fPair_Hi; 
  ChamberPair *fPair_Lo; 

  /// VDC plane w-coordinate (normal-up from vdc plane) (m)
  const double* fW; 
  double fIntercept[4]; 
  
  double fPlane_RMS[4]; //per-plane RMS
  double fPlane_Eta[4]; //per-plane Eta
  int    fPlane_goodPoints[4]; 
    
  double fTheta;
  double fPhi; 
  
  //errors
  double fIntercept_ERR[4] = {-1e30}; 
  double fT0_ERR    =-1e30; 
  
  
  double fTheta_ERR =-1e30;
  double fPhi_ERR   =-1e30; 
  
  TVector2 f_FPInt_xy_ERR = TVector2(-1e30,-1e30); 
  
  TVector3 f_S2Int_xyz;   //intercept with S2-plane
  TVector3 f_FPInt_xyz;   //intercept with focal plane (defined by z=0)
  
  double fEta; 
  
  const double fC = 2.99e8; 
  
  
  //target & fp coordinates
  double f_fp_y;     //focal-plane coords
  double f_fp_theta; 
  double f_fp_phi; 
  
  double f_tg_y;     //target coords
  double f_tg_theta;
  double f_tg_phi; 
  double f_tg_dp; 
    
  
  //intercepts with the S2, in u-v coords
  double fS2_u, fS2_v; 

  //track speed (direction) in 
  double fC_u, fC_v; 
  
  double fT0=0.; 
  
  ApexVDC::HitGroup *fGroup[4]; 
  
  //this tracks the 'good points' which agree well enough with the final track
  ApexVDC::HitGroup *f_goodPointGroup[4]; 
  
  double fSlope_u; 
  double fSlope_v; 
  
  bool f_isGoodTrack;  
  
  ClassDef(Track,0); 
}; 



};

#endif