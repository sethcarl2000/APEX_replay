/***
 * @class TapexVDCChamberPair.h
 * 
 * @brief Class which represents a pair of clusters in two partnered VDC planes  
 * 
 * 
 * 
 ***/

#ifndef TapexVDCChamberPair_H 
#define TapexVDCChamberPair_H

#include "TapexVDCHitGroup.h"
#include "TapexVDCHitCluster.h"
#include <TObject.h> 
#include <vector> 

//we can't include the header file for this, to prevent circular header include statements,
//but we need it in a few function signatures
class TapexVDCTrack;

class TapexVDCChamberPair : public TObject { 
  
 public: 
  TapexVDCChamberPair(bool is_loChamber=true,
                      double u=0,
                      double v=0,
                      TapexVDCHitGroup *Group_U=nullptr, 
                      TapexVDCHitGroup *Group_V=nullptr,
                      int unique_id=-1)
    : fu{u}, fv{v}, fGroup_U{Group_U}, fGroup_V{Group_V}, fUnique_ID{unique_id} {}; 
  
  TapexVDCChamberPair( bool is_loChamber, 
		TapexVDCHitCluster *clust_u, 
		TapexVDCHitCluster *clust_v,
		int unique_id=-1); 
  
  ~TapexVDCChamberPair() {/*noop*/}; 
  
  double u() const { return fu; }
  double v() const { return fv; }
  
  void Set_u(double u) { fu=u; }
  void Set_v(double v) { fv=v; }
  
  void Get_uv(double &u, double &v) { u=fu; v=fv; }
  void Set_uv(double u, double v)   { fu=u; fv=v; }
  
  TapexVDCHitGroup* GetGroup_U() { return fGroup_U; }
  TapexVDCHitGroup* GetGroup_V() { return fGroup_V; }
  
  bool Is_loChamber() const { return f_isLoChamber; }
  
  void SetSlope_uv( double mu, double mv )   { fSlope_u=mu; fSlope_v=mv; }
  void GetSlope_uv( double &mu, double &mv ) { mu=fSlope_u; mv=fSlope_v; }
  
  //gets the wire wire position (truncated downward)
  double ClosestWirePos_Lo( double x ) const; 
  double ClosestWirePos( const double x ) const; 
  
  int Get_ID() const { return fUnique_ID; }
  
  int N_tracks() const { return fTracks.size(); }
  
  void Add_track   ( TapexVDCTrack *track ) { fTracks.push_back( track ); }
  
  void Remove_track( TapexVDCTrack *track ); 
  
  TapexVDCTrack* GetTrack( unsigned int h ) { return fTracks.at(h); } 
  
 private: 
  int fUnique_ID; 
  std::vector<TapexVDCTrack*> fTracks; 
  //this will be used so that tracks can tell if they're using the same clusters
  
  bool f_isLoChamber; 
  TapexVDCHitGroup *fGroup_U; 
  TapexVDCHitGroup *fGroup_V; 
  double fu; 
  double fv; 
  
  double fSlope_u; 
  double fSlope_v; 
  
  ClassDef(TapexVDCChamberPair,0); 
};

#endif 