#include "TapexVDCTrack.h"
#include "TapexVDCNamespace.h"
#include <cmath> 

namespace {
    //in meters/second 
    constexpr double speed_of_light = 2.99e8; 
}

//________________________________________________________________________________________________________
TapexVDCTrack::TapexVDCTrack(TapexEventHandler* event, TapexVDCChamberPair* pair_lo, TapexVDCChamberPair* pair_hi)
    : fEvent{event}, fPair_lo{pair_lo}, fPair_hi{pair_hi}, f_hasVDCdata{false}
{
    //check if all pointers are valid. if not, this track does not have valid data, and we don't need to process it 
    // right now until valid data is provided. 
    if (!(fEvent && fPair_lo && fPair_hi)) return; 

    f_isRightArm = fEvent->ActiveArm(); 

    //add this track to the list of tracks in these pairs
    fPair_lo->Add_track(this); 
    fPair_hi->Add_track(this); 

    //get pointers to the TapexVDCHitGroup's 
    fGroup[1] = fPair_lo->GetGroup_U(); 
    fGroup[0] = fPair_lo->GetGroup_V(); 
    fGroup[3] = fPair_hi->GetGroup_U(); 
    fGroup[2] = fPair_hi->GetGroup_V(); 

    UpdateTrackInfo(); 
    return; 
}
//________________________________________________________________________________________________________
TapexVDCTrack::~TapexVDCTrack() 
{    
    //remove this track from the list of tracks in these TapexVDCChamberPair objects 
    if (fPair_lo) fPair_lo->Remove_track(this); 
    if (fPair_hi) fPair_hi->Remove_track(this); 
};  
//________________________________________________________________________________________________________
void TapexVDCTrack::SetPair_Lo( TapexVDCChamberPair *new_pair ) { 
  
  //if the old pair exists, tell it that we're breaking up with it
  // (it's not you, it's me..) 
  if (fPair_lo) fPair_lo->Remove_track( this );
  
  new_pair->Add_track( this ); 
  
  fGroup[0] = new_pair->GetGroup_U(); 
  fGroup[1] = new_pair->GetGroup_V();

  fPair_lo = new_pair; 
  UpdateTrackInfo(); 
} 
//________________________________________________________________________________________________________
void TapexVDCTrack::SetPair_Hi( TapexVDCChamberPair *new_pair ) { 
  
  //if the old pair exists, tell it that we're breaking up with it
  // (it's not you, it's me..) 
  if (fPair_hi) fPair_hi->Remove_track( this );
  
  new_pair->Add_track( this ); 
  
  fGroup[2] = new_pair->GetGroup_U(); 
  fGroup[3] = new_pair->GetGroup_V();
  
  fPair_hi = new_pair; 
  UpdateTrackInfo();
} 
//________________________________________________________________________________________________________
void TapexVDCTrack::UpdateTrackInfo() 
{
    if (!f_hasVDCdata) return; 
    
    if (!(fPair_lo && fPair_hi)) {
        f_hasVDCdata=false; 
        return; 
    }

    double u1 = fPair_lo->u(); 
    double v1 = fPair_lo->v(); 
    double u2 = fPair_hi->u(); 
    double v2 = fPair_hi->v(); 
    
    fIntercept[0] = v1;
    fIntercept[1] = u1; 
    fIntercept[2] = v2; 
    fIntercept[3] = u2; 

    double plane_w_sep = f_isRightArm 
        ? VDC::R_plane_w[2] - VDC::R_plane_w[0] 
        : VDC::L_plane_w[2] - VDC::L_plane_w[0]; 

    fSlope_u = plane_w_sep/(u2 - u1); 
    fSlope_v = plane_w_sep/(v2 - v1); 

    //compute the electron's speed (in m/s) in the U and V directions
    TVector3 direction = TVector3( 1./fSlope_u, 1./fSlope_v, 1. ).Unit(); 
    fC_u = speed_of_light * direction.x(); 
    fC_v = speed_of_light * direction.y(); 
}
//________________________________________________________________________________________________________
TVector3 TapexVDCTrack::ComputeIntercept_w(const double w) const
{
    double w_u1 = f_isRightArm ? VDC::R_plane_w[0] : VDC::L_plane_w[0]; 
    double w_v1 = f_isRightArm ? VDC::R_plane_w[1] : VDC::L_plane_w[1]; 
    
    return TVector3(
        fIntercept[1] + (w - w_u1)/fSlope_u, 
        fIntercept[0] + (w - w_v1)/fSlope_v, 
        w
    ); 
}
//________________________________________________________________________________________________________
TVector3 TapexVDCTrack::ComputeIntercept_z(const double z) const
{
    TVector3 S  = VDC::UVW_to_XYZ( TVector3(1./fSlope_u, 1./fSlope_v, 1.) );
    
    TVector3 R0 = VDC::UVW_to_XYZ( ComputeIntercept_w(0.) ); 
    
    double t = ( z - R0.z() )/S.z(); 

    return R0 + (S*t); 
}
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________
//________________________________________________________________________________________________________