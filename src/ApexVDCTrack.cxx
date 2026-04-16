
//APEX headers
#include <ApexVDCTrack.h>
#include <ApexVDC.h> 
#include <ApexUtils.h> 
//ROOT headers
#include <TString.h> 
#include <TMatrixD.h> 
#include <TVectorD.h> 
//stdlib headers
#include <math.h> 
#include <limits> 
#include <stdexcept> 

using APEX::kNaN_double;
using APEX::kNull_int; 

namespace {
    
    //rotates uvw to xyz
    constexpr double 
        R_xu( 0.500000),  R_xv(0.500000), R_xw(-0.707107), 
        R_yu(-0.707107),  R_yv(0.707107), R_yw( 0.000000),
        R_zu( 0.500000),  R_zv(0.500000), R_zw( 0.707107); 

    //rotatess xyz to uvw 
    constexpr double 
        R_ux( 0.500000),  R_uy(-0.707107), R_uz( 0.500000), 
        R_vx( 0.500000),  R_vy( 0.707107), R_vz( 0.500000),
        R_wx(-0.707107),  R_wy( 0.000000), R_wz( 0.707107); 


    constexpr double ARRAY_rotate_uvw_to_xyz[9] = {  
         0.500000,  0.500000, -0.707107, 
        -0.707107,  0.707107,  0.000000,
         0.500000,  0.500000,  0.707107 
    }; 

    constexpr double ARRAY_rotate_xyz_to_uvw[9] = {  
         0.500000, -0.707107,  0.500000, 
         0.500000,  0.707107,  0.500000,
        -0.707107,  0.000000,  0.707107 
    }; 

    constexpr char report_prefix[] = "ApexVDC::Track"; 

    //speed of light
    constexpr double kC = 2.99e8; 
}

namespace ApexVDC
{

//______________________________________________________________________________________________________________
void Track::UpdateTrackInfo() { 
  
  //check to see if this track has VDC data. 
  if ( !fPair_Lo || !fPair_Hi ) { f_hasVDCdata=false; } else { f_hasVDCdata=true; }
    
  //tell the point-pairs about their actual slope-values
  double wSep = fW[2]-fW[0]; //the U2-U1 & V2-V1 separations are the same
  
  if (f_hasVDCdata) { 
    //if we DO have vdc data, we want to keep track of the 'pair' intercepts as 
    // best we can.
    fIntercept[0] = fPair_Lo->v(); 
    fIntercept[1] = fPair_Lo->u(); 
    fIntercept[2] = fPair_Hi->v(); 
    fIntercept[3] = fPair_Hi->u(); 
  } 
  
  double 
    uLo(fIntercept[1]),
    vLo(fIntercept[0]),
    uHi(fIntercept[3]),
    vHi(fIntercept[2]); 
  
  fSlope_u = wSep/(uHi-uLo); 
  fSlope_v = wSep/(vHi-vLo);  
  
  if (f_hasVDCdata) { 
    fPair_Hi->SetSlope_uv( fSlope_u, fSlope_v ); 
    fPair_Lo->SetSlope_uv( fSlope_u, fSlope_v ); 
  }
  
  f_S2Int_xyz = ComputeIntercept_z( fEvent->GetS2Hit(f_isRightArm)->Z() ); 
  
  //we need to find the 'speed' of the track in the u & v directions. 
  //the 'true' speed of the track thru the VDC is effectivley c, but 
  //since the track isn't travelling parallel to the 'u' or 'v' axes, we 
  //need to adjust accordingly. this will be used to compute Time-of-Flight (ToF)
  //for the track at any given wire crossing. 
  TVector3 S2Int_uvw =  Track::Rotate_xyz_to_uvw( f_S2Int_xyz ); 
  
  fS2_u = S2Int_uvw[0]; 
  fS2_v = S2Int_uvw[1]; 
    
  TVector3 S_uvw  = TVector3( 1./fSlope_u, 1./fSlope_v, 1. ).Unit(); 
  
  //these are the 'speeds' of the particle in the 'U' or 'V' direction, assuming 
  // the overall speed of the particle is 'c'. 
  fC_u = kC * S_uvw[0]; 
  fC_v = kC * S_uvw[1]; 
  
  //compute focal plane intercept
  f_FPInt_xyz = ComputeIntercept_z(0.); 
    
  TVector3 S_xyz = Track::Rotate_uvw_to_xyz( S_uvw ); 
  
  fTheta = std::atan( S_xyz.X()/S_xyz.Z() ); 
  fPhi   = std::atan( S_xyz.Y()/S_xyz.Z() ); 
}
//______________________________________________________________________________________________________________
Track::Track( const TapexEventHandler *event, 
		      ChamberPair *pLo, 
		      ChamberPair *pHi,
          int unique_id ) : fID{unique_id} 
{ 
  
    if (!event) {
      throw std::logic_error("<ApexVDC::Track::Track> event handler ptr is null"); 
      return; 
    }

    fEvent = event; 
            
    f_isRightArm = fEvent->ActiveArm(); 
    
    fW = f_isRightArm ? ApexVDC::kW_RHRS : ApexVDC::kW_LHRS; 
    
    //if we're making a new monte-carlo track, skip this step. 
    if ( !pLo || !pHi ) { f_hasVDCdata=false; } else { f_hasVDCdata=true; }
    
    if (f_hasVDCdata) { 
        
        fPair_Lo = pLo;   fPair_Lo->Add_track(GetID());
        fPair_Hi = pHi;   fPair_Hi->Add_track(GetID()); 
        
        fGroup[0] = (ApexVDC::HitGroup*)fPair_Lo->GetGroup_U(); 
        fGroup[1] = (ApexVDC::HitGroup*)fPair_Lo->GetGroup_V();
        
        fGroup[2] = (ApexVDC::HitGroup*)fPair_Hi->GetGroup_U(); 
        fGroup[3] = (ApexVDC::HitGroup*)fPair_Hi->GetGroup_V(); 
    
        UpdateTrackInfo(); 
    }
  
}


//______________________________________________________________________________________________________________
void Track::Set_Errors( double errors[5] ) { 
  
  for (int p=0; p<4; p++) fIntercept_ERR[p] =errors[p]; 

  fT0_ERR = errors[4];     
}
//______________________________________________________________________________________________________________
void   Track::Set_Eta(double eta[4]) { 

  for (int p=0; p<4; p++) Set_Eta(p, eta[p]); 
} 
//______________________________________________________________________________________________________________
double Track::Get_Eta() const { 
  
  double eta(0.); 
  for (int p=0; p<4; p++) eta += Get_Eta(p); 
  return eta; 
}
//______________________________________________________________________________________________________________
int    Track::Get_nGoodPoints() const { 
  
  int nPts(0); 
  for (int p=0; p<4; p++) nPts += Get_nGoodPoints(p); 
  return nPts; 
} 
//______________________________________________________________________________________________________________
double Track::Get_RMS() const { 
  
  double RMS(0.); 
  for (int p=0; p<4; p++) RMS += 0.25*TMath::Power( Get_RMS(p), 2. ); 
  return TMath::Sqrt( RMS ); 
} 
//______________________________________________________________________________________________________________
double Track::Error_allIntercept() const { 
  
  double RMS(0.); 
  for (int p=0; p<4; p++) RMS += 0.25*TMath::Power( Error_intercept(p), 2. ); 
  return TMath::Sqrt( RMS ); 
}
//______________________________________________________________________________________________________________
void Track::Set_S2int_angles( double s2x, 
				  double s2y, 
				  double theta, 
				  double phi )   { 
  
  //this is an alternate way to initialize the track without VDC data 
  // (used for monte-carlo processing, in which a track's position is created 
  //  before the corresponding VDC hits are created). 
  
  if ( !fPair_Lo || !fPair_Hi ) { f_hasVDCdata=false; } else { f_hasVDCdata=true; }
    
  f_S2Int_xyz = TVector3( s2x, s2y, fEvent->GetS2Hit(f_isRightArm)->Z() ); 
  
  TVector3 S_xyz( TMath::Tan(theta), TMath::Tan(phi), 1. ); 
  
  TVector3 S_uvw = Track::Rotate_xyz_to_uvw( S_xyz ).Unit(); 
  
  fSlope_u = S_uvw.Z()/S_uvw.X(); 
  fSlope_v = S_uvw.Z()/S_uvw.Y(); 
    
  f_FPInt_xyz = TVector3( f_S2Int_xyz - S_xyz*f_S2Int_xyz.Z() ); 
    
  auto S2u = Track::Rotate_xyz_to_uvw( f_S2Int_xyz ); 
  auto FPu = Track::Rotate_xyz_to_uvw( f_FPInt_xyz ); 
    
  fIntercept[0] = ( fW[0]-FPu.Z() )/fSlope_v  +  FPu.Y(); 
  
  fIntercept[1] = ( fW[1]-FPu.Z() )/fSlope_u  +  FPu.X();
  
  fIntercept[2] = ( fW[2]-FPu.Z() )/fSlope_v  +  FPu.Y(); 
  
  fIntercept[3] = ( fW[3]-FPu.Z() )/fSlope_u  +  FPu.X();

  Track::Compute_Theta_Phi( f_isRightArm, 
				fIntercept, 
				fTheta, fPhi ); 
  
  //fTheta = std::atan( S_xyz.X()/S_xyz.Z() ); 
  //fPhi   = std::atan( S_xyz.Y()/S_xyz.Z() ); 
  
  fS2_u = S2u[0]; 
  fS2_v = S2u[1]; 
    
  fC_u = kC * S_uvw[0]; 
  fC_v = kC * S_uvw[1]; 
  
  UpdateTrackInfo(); 
}
//______________________________________________________________________________________________________________
void Track::SetPair_Hi( ChamberPair *pHi ) { 
  
  //if the old pair exists, tell it that we're breaking up with it
  // (it's not you, it's me..) 
  if (fPair_Hi) fPair_Hi->Remove_track(GetID());
  
  fPair_Hi = pHi; fPair_Hi->Add_track(GetID()); 
  
  fGroup[2] = fPair_Hi->GetGroup_U(); 
  fGroup[3] = fPair_Hi->GetGroup_V(); 
  UpdateTrackInfo(); 
}
//______________________________________________________________________________________________________________
void Track::SetPair_Lo( ChamberPair *pLo ) { 
  
  //if the old pair exists, tell it that we're breaking up with it
  // (it's not you, it's me..) 
  if (fPair_Lo) fPair_Lo->Remove_track(GetID());
  
  fPair_Lo = pLo; fPair_Lo->Add_track(GetID()); 
  
  fGroup[0] = fPair_Lo->GetGroup_U(); 
  fGroup[1] = fPair_Lo->GetGroup_V();
  UpdateTrackInfo(); 
} 
//______________________________________________________________________________________________________________
void Track::Set_uv_Hi( const double u, 
			   const double v ) { 
  if(fPair_Hi) { fPair_Hi->Set_uv( u,v ); } 
  else         { fIntercept[2]=v; fIntercept[3]=u; }
  UpdateTrackInfo(); 
}
//______________________________________________________________________________________________________________
void Track::Set_uv_Lo( const double u, 
			   const double v ) { 
  if(fPair_Lo) { fPair_Lo->Set_uv( u,v ); } 
  else         { fIntercept[0]=v; fIntercept[1]=u; }
  UpdateTrackInfo(); 
}
//______________________________________________________________________________________________________________
double Track::Slope_u() { return fSlope_u; }
//______________________________________________________________________________________________________________
double Track::Slope_v() { return fSlope_v; }

//______________________________________________________________________________________________________________
void   Track::Set_intercept(int plane, double x) { 
  
  if (plane==0) { if (fPair_Lo) { fPair_Lo->Set_v(x); } else { fIntercept[0]=x; } }
  if (plane==1) { if (fPair_Lo) { fPair_Lo->Set_u(x); } else { fIntercept[1]=x; } }
  
  if (plane==2) { if (fPair_Hi) { fPair_Hi->Set_v(x); } else { fIntercept[2]=x; } }
  if (plane==3) { if (fPair_Hi) { fPair_Hi->Set_u(x); } else { fIntercept[3]=x; } }
  
  UpdateTrackInfo(); 
}
//______________________________________________________________________________________________________________
double Track::Intercept(int plane) const { 
  
  return fIntercept[plane]; 
}
//______________________________________________________________________________________________________________
int    Track::Nhits(int plane) const { 
  
  if (!f_hasVDCdata) { 
    throw std::logic_error(Form("<%s::%s> "
        "this track does not have vdc data, illegal call", 
        report_prefix,__func__
    ));
    return kNull_int; 
  }
  if (!f_hasVDCdata || !fGroup[plane]) { 
    throw std::logic_error(Form("<%s::%s> "
        "pointer for ApexVDC::HitGroup for plane %i invalid", 
        report_prefix,__func__, plane
    ));
    return kNull_int; 
  }
  
  return fGroup[plane]->Nhits(); 
}
//______________________________________________________________________________________________________________
double Track::Tau(int plane, int h) const 
{    
    if (!f_hasVDCdata) { 
        throw std::logic_error(Form("<%s::%s> "
            "this track does not have vdc data, illegal call", 
            report_prefix,__func__
        ));
        return kNaN_double; 
    }
    if (!fGroup[plane]) { 
        throw std::logic_error(Form("<%s::%s> "
            "pointer for ApexVDC::HitGroup for plane %i invalid", 
            report_prefix,__func__, plane
        ));
        return kNaN_double; 
    }
    if (h >= Nhits(plane) || h<0 ) { 
        throw std::logic_error(Form("<%s::%s> "
            "hit index %i is invalid. range is [0, %i]", 
            report_prefix,__func__, h, Nhits(plane)-1
        ));
        return kNaN_double; 
    }
    
    return fGroup[plane]->Time(h); 
}
//______________________________________________________________________________________________________________
double Track::WirePos(int plane, int h) const 
{   
    if (!f_hasVDCdata) { 
        throw std::logic_error(Form("<%s::%s> "
            "this track does not have vdc data, illegal call", 
            report_prefix,__func__
        ));
        return kNaN_double; 
    }
    if (!fGroup[plane]) { 
        throw std::logic_error(Form("<%s::%s> "
            "pointer for ApexVDC::HitGroup for plane %i invalid", 
            report_prefix,__func__, plane
        ));
        return kNaN_double; 
    }
    if (h >= Nhits(plane) || h<0 ) { 
        throw std::logic_error(Form("<%s::%s> "
            "hit index %i is invalid. range is [0, %i]", 
            report_prefix,__func__, h, Nhits(plane)-1
        ));
        return kNaN_double; 
    }
    
    return fGroup[plane]->WirePos(h); 
}
//______________________________________________________________________________________________________________
double Track::Slope(int plane) const { 
  
    if (plane==0 || plane==2) { return fSlope_v; }
    if (plane==1 || plane==3) { return fSlope_u; }
    
    throw std::logic_error(Form("<%s::%s> "
        "invalid plane index: %i", 
        report_prefix,__func__, plane
    ));
    return kNaN_double; 
}
//______________________________________________________________________________________________________________
TVector3  Track::Rotate_uvw_to_xyz( const TVector3& v ) 
{   
    TVector3 r_ret; 
    
    int i_mat=0;
    for (int i=0; i<3; i++) 
        for (int j=0; j<3; j++) 
            r_ret[i] += ARRAY_rotate_uvw_to_xyz[i_mat++]*v[j];

    return r_ret; 
}
//______________________________________________________________________________________________________________
TVector3  Track::Rotate_xyz_to_uvw( const TVector3& v ) 
{     
    TVector3 r_ret; 
    
    int i_mat=0;
    for (int i=0; i<3; i++) 
        for (int j=0; j<3; j++) 
            r_ret[i] += ARRAY_rotate_xyz_to_uvw[i_mat++]*v[j];

    return r_ret; 
}
//______________________________________________________________________________________________________________
void Track::Compute_Theta_Phi( 
    const bool is_RHRS,  
    const double intercepts[4], 
    double &Theta, 
    double &Phi ) 
{ //static method

  double wSep = ApexVDC::w( is_RHRS, 2 ) - ApexVDC::w( is_RHRS, 0 );
  
  //compute slopes
  double m_v = wSep/(intercepts[2]-intercepts[0]); 
  double m_u = wSep/(intercepts[3]-intercepts[1]); 

  TVector3 S_xyz = Track::Rotate_uvw_to_xyz( TVector3( 1./m_u, 1./m_v, 1. ) ); 
  
  Theta = std::atan( S_xyz.X()/S_xyz.Z() ); 
  Phi   = std::atan( S_xyz.Y()/S_xyz.Z() ); 
  
}  
//______________________________________________________________________________________________________________
TVector3 Track::ComputeIntercept_w(const double w) const {
    
  return TVector3( fIntercept[1] + ( w - fW[1] )/fSlope_u, 
		   fIntercept[0] + ( w - fW[0] )/fSlope_v, w ); 
}
//______________________________________________________________________________________________________________
TVector3 Track::ComputeIntercept_z(const double z) const {
  
  //get two arbitrary ponints on the track so we can interpolate to find 'z'
  TVector3 S  = Track::Rotate_uvw_to_xyz( TVector3( 1./fSlope_u, 1./fSlope_v, 1. ) );  
  
  TVector3 r0 = Track::Rotate_uvw_to_xyz( ComputeIntercept_w(0.) ); 
  
  double t    = ( z - r0.Z() )/S.Z();   
  
  return TVector3( r0  +  S*t ); 
}
//______________________________________________________________________________________________________________
double Track::ToF(int plane, double x) const { 
  
  //U-plane (v-coord)
  if (plane==0 || plane==2) { return (fS2_v - x)/fC_v; }
  //V-plane (u-coord)
  else                      { return (fS2_u - x)/fC_u; }
}
//______________________________________________________________________________________________________________
double Track::GetTimeAtZ(const double z) const {
  
  //computes the time which this track intercets some z-plane
  // (z in transport coords)
  // note that this time is given relative to the time of the S2-hit used to make
  // this track.
  
  //intercept with this plane, in uvw (vdc) coordinates
  TVector3 intercept_uvw
    = Track::Rotate_xyz_to_uvw( ComputeIntercept_z(z) ); 
  
  double u = intercept_uvw[0];
  
  return (u - fS2_u)/fC_u;
}
//______________________________________________________________________________________________________________
double Track::Get_T_model(int plane, double v, int derivative) const { 
  
  double slope = Slope(plane); 
  
  double w = slope*( v - Intercept(plane) ); 
  
  double Tau_model = fEvent->Drift_T( w,slope,derivative ); 
  
  //time-of-flight does not contribute to the derivatives of T
  if (derivative==0) Tau_model += -ToF( plane, v );   
  
  return Tau_model; 
}
//______________________________________________________________________________________________________________
double Track::xParam() const { 
  
  return  
    ( f_S2Int_xyz.X() - fEvent->GetS2Hit(f_isRightArm)->X() ) 
    / fEvent->GetS2Hit(f_isRightArm)->PaddleWidth();     
}
//______________________________________________________________________________________________________________
void Track::Nudge_params( double nudge[5] ) { 
  
  double u,v; 
  
  if ( !fPair_Lo || !fPair_Hi || !f_hasVDCdata ) { //is this a vdc-data track?
    
    for (int p=0; p<4; p++) fIntercept[p] += nudge[p]; 
    
  } else                                         { 
    
    fPair_Lo->Get_uv( u,v ); 
    fPair_Lo->Set_uv( u + nudge[1], 
		      v + nudge[0] );
    
    fPair_Hi->Get_uv( u,v ); 
    fPair_Hi->Set_uv( u + nudge[3], 
		      v + nudge[2] ); 
  } 
   
  fT0 += nudge[4]; 

  UpdateTrackInfo(); 
}
//______________________________________________________________________________________________________________
void   Track::Set_params( const double vLo, 
			      const double uLo,
			      const double vHi, 
			      const double uHi, 
			      const double T0 ) { 
  
  if (fPair_Lo && fPair_Hi && f_hasVDCdata) {  
    
    fPair_Lo->Set_uv( uLo,vLo );
    fPair_Hi->Set_uv( uHi,vHi ); 
  
  } else  { 
    
    fIntercept[0] =vLo;
    fIntercept[1] =uLo;
    fIntercept[2] =vHi;
    fIntercept[3] =uHi;
  }
  
  fT0 = T0;  

  UpdateTrackInfo(); 
}
//______________________________________________________________________________________________________________
void   Track::Set_params( const double params[5] ) { 
  
    Set_params( 
        params[0], 
        params[1], 
        params[2], 
        params[3], 
        params[4] 
    ); 
}
//______________________________________________________________________________________________________________
double Track::FP_x()  const { return f_FPInt_xyz.X(); }
//______________________________________________________________________________________________________________
double Track::FP_y()  const { return f_FPInt_xyz.Y(); }

//______________________________________________________________________________________________________________
double Track::S2_x()  const { return f_S2Int_xyz.X(); }
//______________________________________________________________________________________________________________
double Track::S2_y()  const { return f_S2Int_xyz.Y(); }

//______________________________________________________________________________________________________________
double Track::Theta() const { return fTheta; }
//______________________________________________________________________________________________________________
double Track::Phi()   const { return fPhi; }
////////////////////////////////////////////////////////////////////////////////////


};