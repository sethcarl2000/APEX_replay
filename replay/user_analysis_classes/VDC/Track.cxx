
//APEX/ headers
#include <APEX/S2Event.h>
#include <APEX/VDC/Track.h>
#include <APEX/VDC.h> 
#include <APEX/utils.h> 
#include <APEX/units.h>
//ROOT headers
#include <TString.h> 
#include <TMatrixD.h> 
#include <TVectorD.h> 
//stdlib headers
#include <math.h> 
#include <limits> 
#include <stdexcept> 



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

    constexpr char report_prefix[] = "APEX::VDC::Track"; 

    //speed of light
    constexpr double kC = APEX::units::speed_of_light_vacuum; 
}

namespace APEX
{
namespace VDC
{

//______________________________________________________________________________________________________________
Track::Track( 
  const S2Event* event, 
  const VDC::ChamberPair<Chamber::Low>& pair_lo, 
  const VDC::ChamberPair<Chamber::High>& pair_hi)
  : fPair_lo{pair_lo}, fPair_hi{pair_hi}, fS2Event{event}
{ 
  if (!fS2Event) {
    throw std::logic_error("<APEX::VDC::Track (ctor)> event handler ptr is null"); 
    return; 
  }

  fT0 = fS2Event->S2Time();   

  /// fill the errors / n.goodponts with null data 
  InvalidateData(); 
}
//______________________________________________________________________________________________________________
void Track::InvalidateData()
{
  // already null 
  if (fErrors.front() != fErrors.front()) return; 
  
  // fill errors with null values
  std::ranges::fill(fErrors.begin(), fErrors.end(), std::numeric_limits<double>::quiet_NaN()); 

  // fill n. goodpoints with null values
  std::ranges::fill(fN_goodpoints_per_plane.begin(), fN_goodpoints_per_plane.end(), -1); 
}
//______________________________________________________________________________________________________________
bool Track::Is_RHRS() const { 
  if (!fS2Event) {
    throw std::logic_error("<APEX::VDC::Track::"+std::string{__func__}+">: evt ptr is null"); 
  }
  return fS2Event->Is_RHRS(); 
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
TVector3 Track::ComputeIntercept_w(const double w) const {
    
  return TVector3( 
    GetIntercept<Plane::V1>() + ( w - VDC::w<Plane::V1>(Is_RHRS()) )/Slope_u(), 
    GetIntercept<Plane::U1>() + ( w - VDC::w<Plane::U1>(Is_RHRS()) )/Slope_v(), 
    w 
  ); 
}
//______________________________________________________________________________________________________________
TVector3 Track::ComputeIntercept_z(const double z) const {
  
  //get two arbitrary ponints on the track so we can interpolate to find 'z'
  TVector3 S  = Track::Rotate_uvw_to_xyz( TVector3(1./Slope_u(), 1./Slope_v(), 1.) );  
  
  TVector3 r0 = Track::Rotate_uvw_to_xyz( ComputeIntercept_w(0.) ); 
  
  double t    = ( z - r0.Z() )/S.Z();   
  
  return TVector3( r0  +  S*t ); 
}
//______________________________________________________________________________________________________________
TVector3 Track::Get_S2_intercept() const
{
  return ComputeIntercept_z( fS2Event->GetS2Hit().Z() ); 
}
//______________________________________________________________________________________________________________
template<Plane P> double Track::ToF(double pos) const { 
  
  auto S2_uvw = Track::Rotate_xyz_to_uvw( Get_S2_intercept() ); 

  auto velocity = S2_uvw.Unit() * units::speed_of_light_vacuum; 
  
  if constexpr ((P==Plane::V1) || (P==Plane::V2)) { //U-plane (v-coord)
    return (S2_uvw.x() - pos)/velocity.x(); 
  } else {  //V-plane (u-coord)
    return (S2_uvw.y() - pos)/velocity.y(); 
  }
}
//______________________________________________________________________________________________________________
template<Plane P> double Track::Get_T_model(double x, int derivative) const 
{   
  double slope = Slope<P>(); 
  
  //drift-distance 
  double w = slope*( x - GetIntercept<P>() ); 
  
  double Tau_model = fS2Event->Drift_T( w,slope,derivative ); 
  
  //time-of-flight does not contribute to the derivatives of T
  if (derivative==0) Tau_model += -ToF<P>(x);   
  
  return Tau_model; 
}
//______________________________________________________________________________________________________________
double Track::xParam() const 
{ 
  return ( Get_S2_intercept().X() - fS2Event->GetS2Hit().X() ) / S2::paddle_width_X;     
}
//______________________________________________________________________________________________________________
void Track::Nudge_params(std::span<const double, VDC::kN_planes + 1> nudge) 
{ 
  /// Nudge each VDC intercept
  VDC::for_each_plane([nudge,this]<VDC::Plane P>(){

    this->SetIntercept<P>( this->GetIntercept<P>() + nudge[VDC::number<P>()] ); 
  }); 

  /// nudge the timing offset 
  fT0 += nudge.back(); 
}
// Explicit template instantiations 
//______________________________________________________________________________________________________________
template ChamberPair<Chamber::Low>&  Track::GetPair();
template ChamberPair<Chamber::High>& Track::GetPair();
//______________________________________________________________________________________________________________
template const ChamberPair<Chamber::Low>&  Track::GetPair() const;
template const ChamberPair<Chamber::High>& Track::GetPair() const;
//______________________________________________________________________________________________________________
template void Track::Set_uv<Chamber::Low>(double,double); 
template void Track::Set_uv<Chamber::High>(double,double); 
//______________________________________________________________________________________________________________
template VDC::HitCluster<Plane::U1>& Track::GetCluster(); 
template VDC::HitCluster<Plane::V1>& Track::GetCluster(); 
template VDC::HitCluster<Plane::U2>& Track::GetCluster(); 
template VDC::HitCluster<Plane::V2>& Track::GetCluster(); 
//______________________________________________________________________________________________________________
template const VDC::HitCluster<Plane::U1>& Track::GetCluster() const; 
template const VDC::HitCluster<Plane::V1>& Track::GetCluster() const; 
template const VDC::HitCluster<Plane::U2>& Track::GetCluster() const; 
template const VDC::HitCluster<Plane::V2>& Track::GetCluster() const; 
//______________________________________________________________________________________________________________
template VDC::HitGroup<Plane::U1>& Track::GetGroup(); 
template VDC::HitGroup<Plane::V1>& Track::GetGroup(); 
template VDC::HitGroup<Plane::U2>& Track::GetGroup(); 
template VDC::HitGroup<Plane::V2>& Track::GetGroup(); 
//______________________________________________________________________________________________________________
template const VDC::HitGroup<Plane::U1>& Track::GetGroup() const; 
template const VDC::HitGroup<Plane::V1>& Track::GetGroup() const; 
template const VDC::HitGroup<Plane::U2>& Track::GetGroup() const; 
template const VDC::HitGroup<Plane::V2>& Track::GetGroup() const; 
//______________________________________________________________________________________________________________
template double Track::Slope<Plane::U1>() const; 
template double Track::Slope<Plane::V1>() const; 
template double Track::Slope<Plane::U2>() const; 
template double Track::Slope<Plane::V2>() const; 
//______________________________________________________________________________________________________________
template double Track::GetIntercept<Plane::U1>() const; 
template double Track::GetIntercept<Plane::V1>() const; 
template double Track::GetIntercept<Plane::U2>() const; 
template double Track::GetIntercept<Plane::V2>() const; 
//______________________________________________________________________________________________________________
template void Track::SetIntercept<Plane::U1>(double); 
template void Track::SetIntercept<Plane::V1>(double); 
template void Track::SetIntercept<Plane::U2>(double); 
template void Track::SetIntercept<Plane::V2>(double); 
//______________________________________________________________________________________________________________
template double Track::Tau<Plane::U1>(int) const; 
template double Track::Tau<Plane::V1>(int) const; 
template double Track::Tau<Plane::U2>(int) const; 
template double Track::Tau<Plane::V2>(int) const; 
//______________________________________________________________________________________________________________
template double Track::WirePos<Plane::U1>(int) const; 
template double Track::WirePos<Plane::V1>(int) const; 
template double Track::WirePos<Plane::U2>(int) const; 
template double Track::WirePos<Plane::V2>(int) const; 
//______________________________________________________________________________________________________________
template double Track::ToF<Plane::U1>(double) const; 
template double Track::ToF<Plane::V1>(double) const; 
template double Track::ToF<Plane::U2>(double) const; 
template double Track::ToF<Plane::V2>(double) const; 
//______________________________________________________________________________________________________________
template double Track::Get_T_model<Plane::U1>(double,int) const; 
template double Track::Get_T_model<Plane::V1>(double,int) const; 
template double Track::Get_T_model<Plane::U2>(double,int) const; 
template double Track::Get_T_model<Plane::V2>(double,int) const; 
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________
//______________________________________________________________________________________________________________


}
}