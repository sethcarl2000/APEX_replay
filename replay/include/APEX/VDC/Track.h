#ifndef APEX_VDC_Track_H
#define APEX_VDC_Track_H 

// APEX headers 
#include <APEX/ObjCollection.h> 
#include <APEX/Cerenkov.h>
#include <APEX/VDC.h> 
#include <APEX/VDC/ChamberPair.h> 
#include <APEX/EventHandler.h> 
// ROOT headers
#include <TVector3.h>
// stdlib
#include <cmath> 
#include <stdexcept> 
#include <span> 
#include <algorithm>
#include <array> 
#include <numeric> 

namespace APEX 
{
class S2Event; 

namespace VDC
{
    
class Track { 
private: 

  // we need to friend this so the 'id' field is accessible 
  friend class ObjCollection<Track>; 

  ChamberPair<Chamber::Low>  fPair_lo; //lower chamber pair
  ChamberPair<Chamber::High> fPair_hi; //upper chamber pair 

  const S2Event* fS2Event{nullptr}; // S2Event 

  double fT0; // time offset rel. to S2 time 

  // errors for each intercept / time. the inidces are: 
  /// [0] <= error of U1 intercept
  /// [1] <= error of V1 intercept
  /// [2] <= error of U2 intercept
  /// [3] <= error of V2 intercept
  /// [4] <= error of T0 (time shift rel. to S2 hit)
  std::array<double, VDC::kN_planes + 1> fErrors; 

  std::array<int, VDC::kN_planes> fN_goodpoints_per_plane; 

  /// @brief Any time an intercept / timing offset is shifted, we have to nullify the meta-data (like n. goodpoints and errors)
  void InvalidateData(); 

public: 
  // This is the track object, which handles self-refinement, 
  // as well as acting as a container for its constituent points. 
  Track( 
    const S2Event* event,
    const ChamberPair<Chamber::Low>&  pair_low, 
    const ChamberPair<Chamber::High>& pair_hi
  ); 
  
  /// @brief Get chamber pair 
  /// @tparam C VDC Chamber 
  /// @return the ChamberPair<C> this track owns
  template<Chamber C> ChamberPair<C>& GetPair() {
    if constexpr (VDC::is_loChamber<C>()) { return fPair_lo; } else { return fPair_hi; }
  }

  /// @brief Get chamber pair (const)
  /// @tparam C VDC Chamber 
  /// @return the ChamberPair<C> this track owns
  template<Chamber C> const ChamberPair<C>& GetPair() const {
    if constexpr (VDC::is_loChamber<C>()) { return fPair_lo; } else { return fPair_hi; }
  }
  
  ChamberPair<Chamber::Low>&  GetPair_Lo() { return GetPair<Chamber::Low>(); }
  ChamberPair<Chamber::High>& GetPair_Hi() { return GetPair<Chamber::High>(); }

  const ChamberPair<Chamber::Low>&  GetPair_Lo() const { return GetPair<Chamber::Low>(); }
  const ChamberPair<Chamber::High>& GetPair_Hi() const { return GetPair<Chamber::High>(); }
  
  void SetEvent(const S2Event* ptr) { fS2Event=ptr; }
  const S2Event* GetEvent() const { return fS2Event; }

  bool Is_RHRS() const; 
  
  /// @brief Set u & v intercept for Low/Hi-chamber. Recall that the 'U' plane determines the 'v' intercept; and the 'V' planes determines the 'u' intercept. 
  /// @tparam C VDC Chamber
  /// @param u u-intercept
  /// @param v v-intercept
  template<Chamber C> void Set_uv(double u, double v) {
    GetPair<C>().Set_uv(u,v); 
    InvalidateData(); 
  }
  
  /// @brief Set u & v intercept for Low-chamber. Recall that the 'U' plane determines the 'v' intercept; and the 'V' planes determines the 'u' intercept. 
  /// @param u u-intercept
  /// @param v v-intercept
  void Set_uv_Lo(double u, double v) { Set_uv<Chamber::Low>(u,v); }

  /// @brief Set u & v intercept for High-chamber. Recall that the 'U' plane determines the 'v' intercept; and the 'V' planes determines the 'u' intercept. 
  /// @param u u-intercept
  /// @param v v-intercept
  void Set_uv_Hi(double u, double v) { Set_uv<Chamber::High>(u,v); } 
  
  double T0() const { return fT0; }
  void   Set_T0(const double T0) { fT0=T0; InvalidateData(); }
    
  //returns the agreement of this hit with the S2-paddle
  double xParam() const; 
  
  //adjust the params by this amount
  void Nudge_params(std::span<const double, VDC::kN_planes + 1> nudge); 

  /// @brief This track's hit-cluster corresponding to VDC Plane P
  /// @tparam P VDC Plane
  /// @return Ref. to Hit Cluster
  template<Plane P> VDC::HitCluster<P>& GetCluster() {
    static constexpr Chamber Ch = VDC::get_Chamber<P>(); 
    auto& pair = GetPair<Ch>(); 
    if constexpr (P == VDC::get_U<Ch>()) { return pair.GetCluster_U(); } else { return pair.GetCluster_V(); }
  }

  /// @brief This track's hit-cluster corresponding to VDC Plane P
  /// @tparam P VDC Plane
  /// @return Ref. to Hit Cluster
  template<Plane P> const VDC::HitCluster<P>& GetCluster() const {
    static constexpr Chamber Ch = VDC::get_Chamber<P>(); 
    auto& pair = GetPair<Ch>(); 
    if constexpr (P == VDC::get_U<Ch>()) { return pair.GetCluster_U(); } else { return pair.GetCluster_V(); }
  }

  /// @brief This track's hit-group corresponding to VDC Plane P
  /// @tparam P VDC Plane
  /// @return Ref. to Hit Group
  template<Plane P> VDC::HitGroup<P>& GetGroup() { 
    return GetCluster<P>().GetGroup(); 
  }

  /// @brief This track's hit-group corresponding to VDC Plane P
  /// @tparam P VDC Plane
  /// @return Ref. to Hit Group
  template<Plane P> const VDC::HitGroup<P>& GetGroup() const { 
    return GetCluster<P>().GetGroup(); 
  }

  /// @return this track's (dw/du) slope in the VDC's UVW coordinates 
  double Slope_u() const {
    const auto& clust_V1 = GetCluster<Plane::V1>(); 
    const auto& clust_V2 = GetCluster<Plane::V2>(); 
    double w_sep = VDC::w<Plane::V2>(Is_RHRS()) - VDC::w<Plane::V1>(Is_RHRS()); 
    return w_sep/(clust_V2.GetIntercept() - clust_V1.GetIntercept()); 
  }; 

  /// @return this track's (dw/dv) slope in the VDC's UVW coordinates 
  double Slope_v() const {
    const auto& clust_U1 = GetCluster<Plane::U1>(); 
    const auto& clust_U2 = GetCluster<Plane::U2>(); 
    double w_sep = VDC::w<Plane::U2>(Is_RHRS()) - VDC::w<Plane::U1>(Is_RHRS()); 
    return w_sep/(clust_U2.GetIntercept() - clust_U1.GetIntercept()); 
  }; 

  /// @brief Get VDC slope (dw/du or dw/dv) for Plane P 
  /// @tparam P VDC Plane
  /// @return slope (dw/du or dw/dv)
  template<Plane P> double Slope() const {
    if constexpr ((P == Plane::U1) || (P == Plane::U2)) { return Slope_v(); } else { return Slope_u(); }
  }
   
  /// @tparam P VDC Plane 
  /// @return Track intercept, in VDC UVW coords  
  template<Plane P> double GetIntercept() const { return GetCluster<P>().GetIntercept(); } 
   
  /// @tparam P VDC Plane 
  /// @return Track intercept, in VDC UVW coords  
  template<Plane P> void SetIntercept(double x) { return GetCluster<P>().SetIntercept(x); } 

  /// @brief Set 'Eta'
  /// @tparam P VDC Plane  
  template<Plane P> void SetEta(double eta) { GetCluster<P>().SetEta(eta); }

  /// @brief Set 'Eta' for each plane
  void SetEta(std::span<const double,VDC::kN_planes> _data) { 
    VDC::for_each_plane([_data,this]<VDC::Plane P>(){
      SetEta<P>( _data[VDC::number<P>()] );
    });
  }

  /// @tparam P VDC Plane 
  /// @return 'Eta' 
  template<Plane P> double GetEta() const { return GetCluster<P>().GetEta(); }

  /// @return Sum of 'Eta' for all planes 
  double GetEta() const { 
    double sum=0.; 
    VDC::for_each_plane([&sum,this]<VDC::Plane P>(){
      sum += GetCluster<P>().GetEta(); 
    }); 
    return sum; 
  }

  /// @brief Set errors for each VDC plane intercept + timing offset
  void SetErrors(std::span<const double, VDC::kN_planes+1> _err) { std::ranges::copy(_err, fErrors.begin()); } 

  /// @brief Get intercept error est. for given VDC plane
  /// @tparam P VDC Plane
  /// @return Est. of intercept RMS
  template<Plane P> double GetInterceptError() const {
    return fErrors[VDC::number<P>()]; 
  }

  /// @return Est. of T0 RMS
  double GetT0Error() const {
    return fErrors.back(); 
  }

  /// @brief Set the number of goodpoints for a given plane
  /// @tparam P VDC Plane 
  /// @param n_goodpoints n. of 'goodpoints' 
  template<Plane P> void Set_nGoodPoints(int n_goodpoints) {
    fN_goodpoints_per_plane[VDC::number<P>()] = n_goodpoints; 
  }

  /// @brief Get the number of goodpoints for a given plane
  /// @tparam P VDC Plane 
  /// @param n_goodpoints Number of 'goodpoints'
  template<Plane P> int Get_nGoodPoints() const {
    return fN_goodpoints_per_plane[VDC::number<P>()];
  }

  /// @return Number of 'goodpoints' across all 4 planes
  int Get_nGoodPoints() const {
    return std::accumulate(fN_goodpoints_per_plane.begin(), fN_goodpoints_per_plane.end(), 0); 
  }
  
  template<Plane P> double Tau(int h)       const { return GetGroup<P>().Time(h); } 
  template<Plane P> double WirePos(int h)   const { return GetGroup<P>().WirePos(h); } 
  
  static TVector3 Rotate_uvw_to_xyz( const TVector3& vec ); 
  static TVector3 Rotate_xyz_to_uvw( const TVector3& vec ); 
  
  TVector3 ComputeIntercept_w(const double w) const; 
  
  TVector3 ComputeIntercept_z(const double z) const;   

  /// @return The 3-vector position of this track's intercept with the S2 plane. 
  TVector3 Get_S2_intercept() const; 

  /// @brief Compute the time-of-fligth between a certain VDC coordinate (v-pos or u-pos), and this track's intercept with the S2 plane. 
  /// @tparam P VDC plane 
  /// @param pos position in VDC-coordinates
  /// @return time-of-flight between this track's intercept with the S2-plane, and a plane defined by a given VDC coordinate
  template<Plane P> double ToF(double pos) const; 

  template<Plane P> double Get_T_model(double x, int derivative=0) const; 
  
  //overall track parameters (in focal-plane coordinates)
  double Transport_x()     const { return ComputeIntercept_z(0.).X(); }
  double Transport_y()     const { return ComputeIntercept_z(0.).Y(); } 
  double Transport_dx_dz() const {
    auto S = Track::Rotate_uvw_to_xyz(TVector3(1./Slope_u(), 1./Slope_v(), 1.)); 
    return S.x() / S.z(); 
  };
  double Transport_dy_dz() const {
    auto S = Track::Rotate_uvw_to_xyz(TVector3(1./Slope_u(), 1./Slope_v(), 1.)); 
    return S.y() / S.z(); 
  };

  double Theta()  const { return Transport_dx_dz(); }
  double Phi()    const { return Transport_dy_dz(); } 

  //overall track parameters (in focal-plane coordinates)
  double FP_x()     const { return Transport_x(); }
  double FP_y()     const { return Transport_y(); } 
  double FP_dx_dz() const { return Transport_dx_dz() - Transport_x()/6.; }
  double FP_dy_dz() const { return Transport_dy_dz(); }

private:
  std::size_t fS2_id; 
public: 
  void Set_S2_id(std::size_t id) { fS2_id=id; }
  std::size_t Get_S2_id() const { return fS2_id; }
  
    
  /// @return X-intercept of this track with the S2 plane, TRANSPORT coordinates
  double S2_x()  const { return Get_S2_intercept().X(); }
  /// @return Y-intercept of this track with the S2 plane, TRANSPORT coordinates
  double S2_y()  const { return Get_S2_intercept().Y(); } 
  /// @return time of this track's arrival at the S2 plane. 
  double S2_t()  const { return fT0; }


  /// @return X-intercept of this track with the Cerenkov plane, TRANSPORT coordinates
  double Cerenkov_x() const {
    return ComputeIntercept_z(Cerenkov::Z(Is_RHRS())).X(); 
  }; 
  /// @return Y-intercept of this track with the Cerenkov plane, TRANSPORT coordinates
  double Cerenkov_y() const {
    return ComputeIntercept_z(Cerenkov::Z(Is_RHRS())).Y(); 
  };  
  /// @return time-of-intercept of this track with the Cerenkov plane (rel. to S2 hit)
  double Cerenkov_t() const {

    double dist_to_S2 = 
      (Get_S2_intercept() - ComputeIntercept_z(Cerenkov::Z(Is_RHRS()))).Mag(); 

    return fT0 - (dist_to_S2/units::speed_of_light_vacuum); 
  }; 


private: 
  object_id id;
public: 
  object_id ID() const { return id; }
  bool operator==(const Track& rhs) const { return this->id == rhs.id; }
}; 

}
}

#endif
