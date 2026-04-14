#ifndef ApexVDCChamberPair_H 
#define ApexVDCChamberPair_H 

#include <ApexVDC.h> 
#include <ApexVDCHitGroup.h> 
#include <ApexVDCHitCluster.h> 

namespace ApexVDC
{
    
class Track; 

class ChamberPair { 
  
 public: 
  ChamberPair( bool is_loChamber=true,
		double u=0,
		double v=0,
	        HitGroup *Group_U=0, 
	        HitGroup *Group_V=0,
		int unique_id=-1); 
  
  ChamberPair( bool is_loChamber, 
		HitCluster *clust_u, 
		HitCluster *clust_v,
		int unique_id=-1); 

  
  double u() const { return fu; }
  double v() const { return fv; }
  
  void Set_u(double u) { fu=u; }
  void Set_v(double v) { fv=v; }
  
  void Get_uv(double &u, double &v) { u=fu; v=fv; }
  void Set_uv(double u, double v)   { fu=u; fv=v; }
  
  HitGroup* GetGroup_U() { return fGroup_U; }
  HitGroup* GetGroup_V() { return fGroup_V; }
  
  bool Is_loChamber() const { return f_isLoChamber; }
  
  void SetSlope_uv( double mu, double mv )   { fSlope_u=mu; fSlope_v=mv; }
  void GetSlope_uv( double &mu, double &mv ) { mu=fSlope_u; mv=fSlope_v; }
  
  double ClosestWirePos_Lo( double x ) const; 
  double ClosestWirePos( const double x ) const; 
  
  int Get_ID() const { return fUnique_ID; }
  
  int N_tracks() const { return fTracks.size(); }
  
  void Add_track   ( Track *track ) { fTracks.push_back( track ); }
  
  void Remove_track( Track *track ); 
  
  Track* GetTrack( unsigned int h ) { return fTracks.at(h); } 
  
 private: 
  int fUnique_ID; 
  std::vector<Track*> fTracks; 
  //this will be used so that tracks can tell if they're using the same clusters
  
  bool f_isLoChamber; 
  HitGroup *fGroup_U; 
  HitGroup *fGroup_V; 
  double fu; 
  double fv; 
  
  double fSlope_u; 
  double fSlope_v; 
  
  ClassDef(ChamberPair,0); 
};

}; 


#endif 