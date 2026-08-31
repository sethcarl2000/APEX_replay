#ifndef APEX_VDC_Track_H
#define APEX_VDC_Track_H 

//APEX headers 
#include <APEX/VDC.h> 
#include <APEX/VDC/ChamberPair.h> 
#include <APEX/EventHandler.h> 
#include <math.h> 
//ROOT headers
#include <TVector3.h>

namespace APEX 
{
namespace VDC
{
    
class Track { 
  
 public: 
  // This is the track object, which handles self-refinement, 
  // as well as acting as a container for its constituent points. 
  Track( 
    const EventHandler *event=nullptr, 
    ChamberPair *pLo=0, 
    ChamberPair *pHi=0,
    int unique_id=-1
  ); 
  
  inline ~Track() {
    //cout << "Track:: deleting track..."<<endl; 
    
    //if (fPair_Lo) fPair_Lo->Remove_track(GetID()); 
    //if (fPair_Hi) fPair_Hi->Remove_track(GetID()); 
  }; //noop
  
  void SetPair_Hi( ChamberPair *pHi );  
  void SetPair_Lo( ChamberPair *pLo );  
  
  inline ChamberPair* GetPair_Hi() { return fPair_Hi; }
  inline ChamberPair* GetPair_Lo() { return fPair_Lo; }
  
  inline void SetEvent( EventHandler *evt ) { fEvent=evt; }
  inline const EventHandler *GetEvent() const { return fEvent; }
  
  inline bool IsRightArm() const { return f_isRightArm; }
  
  void Set_uv_Hi( const double u, const double v ); 
  
  void Set_uv_Lo( const double u, const double v );
  
  inline double T0() const { return fT0; }
  inline void   Set_T0(const double T0) { fT0=T0; }
    
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
  
  inline APEX::VDC::HitGroup *GetGroup(int plane) { return fGroup[plane]; }
    
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
    
  //overall track parameters (in focal-plane coordinates)
  double FP_x()     const; 
  double FP_y()     const; 
  double FP_dx_dz() const;
  double FP_dy_dz() const; 

  //optics data ~~~~~~~~~~~~~~~~~~~~~~~~~~ (Computation handled by THRS class)
  //track target coordinates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //placeholder functions 
  inline double Transport_x() const { return FP_x(); } 
  inline double Transport_y() const { return FP_y(); } 
  inline double dx_dz()       const { return std::tan(Theta()); } 
  inline double dy_dz()       const { return std::tan(Phi()); } 
    
  inline void Set_targetCoords( const double tg_y, 
			 const double tg_theta, 
			 const double tg_phi, 
			 const double tg_dp ) 
  { f_tg_y     =tg_y;
    f_tg_theta =tg_theta; 
    f_tg_phi   =tg_phi; 
    f_tg_dp    =tg_dp; } 
  
  //all angles in radians 
  inline double Get_tg_y()     const { return f_tg_y; } 
  inline double Get_tg_theta() const { return f_tg_theta; } 
  inline double Get_tg_phi()   const { return f_tg_phi; } 
  
  
  //track fp coordinates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  inline void Set_focalPlaneCoords( const double fp_y, 
			     const double fp_theta, 
			     const double fp_phi ) 
  { f_fp_y     =fp_y;
    f_fp_theta =fp_theta; 
    f_fp_phi   =fp_phi; } 
  
  //all angles in radians 
  inline double Get_fp_y()     const { return f_fp_y; } 
  inline double Get_fp_theta() const { return f_fp_theta; } 
  inline double Get_fp_phi()   const { return f_fp_phi; } 
    
  
  double GetTimeAtZ(const double z) const; 
  
  
  double S2_x()  const; 
  double S2_y()  const; 
  
  double Theta() const; 
  double Phi()   const; 
  
  //set/get track errors. 
  void Set_Errors( double errors[5] ); 
  
  inline double Error_intercept(int plane) const { return fIntercept_ERR[plane]; }
  inline double Error_T0()                 const { return fT0_ERR; }
  
  //a general figure which is a stand-in for the error of all track intercepts
  double Error_allIntercept()       const; 
    
  
  inline double Error_Theta()              const { return fTheta_ERR; }
  inline double Error_Phi()                const { return fPhi_ERR; }
  
  
  //
  inline void Set_goodPointGroup( int plane, APEX::VDC::HitGroup *goodPointGroup ) {
    f_goodPointGroup[plane] = goodPointGroup; 
  }
  
  inline APEX::VDC::HitGroup *Get_goodPointGroup( int plane ) { 
    return f_goodPointGroup[plane]; 
  }
  
  inline void   Set_Eta(int plane, double eta) { fPlane_Eta[plane]=eta; }
  void   Set_Eta(double eta[4]); 
  inline double Get_Eta(int plane) const       { return fPlane_Eta[plane]; }
  double Get_Eta()          const; //for all planes
  
  inline void   Set_RMS(int plane, double RMS) { fPlane_RMS[plane]=RMS; }
  inline double Get_RMS(int plane) const       { return fPlane_RMS[plane]; }
  double Get_RMS()          const; //for all planes
  
  inline void Set_nGoodPoints(int plane, int pts) { fPlane_goodPoints[plane]=pts; }
  inline int  Get_nGoodPoints(int plane) const    { return fPlane_goodPoints[plane]; }
  int  Get_nGoodPoints()          const; //for all planes
    
  
  //this flag is used for monte-carlo testing
  inline void Set_isGoodTrack(bool isGood) { f_isGoodTrack =isGood; } 
  inline bool IsGoodTrack() const { return f_isGoodTrack; }
    
  static TVector3 Rotate_uvw_to_xyz( const TVector3& vec ); 
  static TVector3 Rotate_xyz_to_uvw( const TVector3& vec ); 
  
  static void Compute_Theta_Phi( const bool arm, 
				 const double intercepts[4], 
				 double &Theta, 
				 double &Phi ); 
    
  enum ECoordType { kDetector=0, kTransport, kFocalPlane };  

  TVector3 ComputeIntercept_w(const double w) const; 
  
  TVector3 ComputeIntercept_z(const double z) const;   

  //use the uniqueness of each track's memory address to define a unique instance of a track
  inline bool operator==(const APEX::VDC::Track& rhs) { return (this)==(&rhs); }
  
  /// @return unique ID of this track 
  int GetID() const { return fID; }

  inline bool operator==(const Track& rhs) const { return (rhs.IsRightArm()==IsRightArm()) && (rhs.GetID()==GetID()); }

 private: 
  const EventHandler *fEvent;
  bool f_isRightArm; 
    
  //unique id assigned to this track. 
  int fID; 

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
  
  APEX::VDC::HitGroup *fGroup[4]; 
  
  //this tracks the 'good points' which agree well enough with the final track
  APEX::VDC::HitGroup *f_goodPointGroup[4]; 
  
  double fSlope_u; 
  double fSlope_v; 
  
  bool f_isGoodTrack;  
  
}; 

}
}

#endif
