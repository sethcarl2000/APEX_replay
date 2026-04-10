#ifndef Podd_THaTrack_h_
#define Podd_THaTrack_h_

//////////////////////////////////////////////////////////////////////////
//
// THaTrack
//
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TVector3.h"
#include "THaPIDinfo.h"
#include <cstring>   // for memset

#include "THaVDCPointPair.h"
#include "THaVDCPoint.h"
#include "THaVDCCluster.h"

class TClonesArray;
class THaTrackingDetector;
class THaCluster;
class THaTrackID;

class THaTrack : public TObject {

public:

  // Bits for fType
  enum {
    kHasDet        = BIT(0),  // Detector coordinates set
    kHasFP         = BIT(1),  // Focal plane coordinates set
    kHasRot        = BIT(2),  // Rotating TRANSPORT coordinates set
    kHasTarget     = BIT(3),  // Target coordinates reconstructed
    kHasVertex     = BIT(4)   // Vertex reconstructed
  };

  // Default constructor
  THaTrack()
    : TObject(),
    fX(kBig), fY(kBig), fTheta(kBig), fPhi(kBig), fP(kBig),
    fDX(kBig), fDY(kBig), fDTheta(kBig), fDPhi(kBig),
    fRX(kBig), fRY(kBig), fRTheta(kBig), fRPhi(kBig),
    fTX(kBig), fTY(kBig), fTTheta(kBig), fTPhi(kBig), fDp(kBig),
    fPvect(kBig,kBig,kBig), fVertex(kBig,kBig,kBig),
    fVertexError(kBig,kBig,kBig),
    fPathl(kBig), fToff(kBig), fTime(kBig), fdTime(kBig), fBeta(kBig), fdBeta(kBig),
    fNclusters(0), fPIDinfo(0), fCreator(0), fIndex(-1), fTrkNum(0),
    fID(0), fFlag(0), fType(0), fChi2(kBig), fNDoF(0),
    fDedx(kBig), fEnergy(kBig),
    fNPMT(0), fBetaChi2(kBig), fFPTime(kBig),
    fGoodPlane3(0), fGoodPlane4(0),
    fUV12X(kBig),fUV12Y(kBig),fUV12PX(kBig),fUV12PY(kBig),
    fUV21X(kBig),fUV21Y(kBig),fUV21PX(kBig),fUV21PY(kBig)
  { memset(fClusters,0,kMAXCL*sizeof(THaCluster*)); }

  // Constructor with fp coordinates
  // FIXME: this really should be setting detector coordinates
  THaTrack( Double_t x, Double_t y, Double_t theta, Double_t phi,
	    THaTrackingDetector* creator=0, THaTrackID* id=0,
	    THaPIDinfo* pid=0 )
    : TObject(),
      fX(x), fY(y), fTheta(theta), fPhi(phi), fP(kBig),
      fDX(kBig), fDY(kBig), fDTheta(kBig), fDPhi(kBig),
      fRX(kBig), fRY(kBig), fRTheta(kBig), fRPhi(kBig),
      fTX(kBig), fTY(kBig), fTTheta(kBig), fTPhi(kBig), fDp(kBig),
      fPvect(kBig,kBig,kBig), fVertex(kBig,kBig,kBig),
      fVertexError(kBig,kBig,kBig),
      fPathl(kBig), fToff(kBig), fTime(kBig), fdTime(kBig), fBeta(kBig), fdBeta(kBig),
      fNclusters(0), fPIDinfo(pid), fCreator(creator), fIndex(-1), fTrkNum(0),
      fID(id), fFlag(0), fType(kHasFP), fChi2(kBig), fNDoF(0),
      fDedx(kBig), fEnergy(kBig),
      fNPMT(0), fBetaChi2(kBig), fFPTime(kBig),
      fGoodPlane3(0), fGoodPlane4(0)
  {
    memset(fClusters,0,kMAXCL*sizeof(THaCluster*));
    if(pid) pid->Clear();
  }

  virtual ~THaTrack();

  Int_t             AddCluster( THaCluster* c );
  void              Clear( Option_t* opt="" );
  THaTrackingDetector* GetCreator()    const { return fCreator; }
  Int_t             GetNclusters()     const { return fNclusters; }
  Int_t             GetIndex()         const { return fIndex; }
  THaCluster*       GetCluster( Int_t i )    { return fClusters[i]; }
  UInt_t            GetFlag()          const { return fFlag; }
  UInt_t            GetType()          const { return fType; }
  THaTrackID*       GetID()            const { return fID; }
  Int_t             GetTrkNum()        const { return fTrkNum; }

  Double_t          GetError()         const { return fError; }
  Double_t          GetP()             const { return fP; }
  Double_t          GetPhi()           const { return fPhi; }
  THaPIDinfo*       GetPIDinfo()       const { return fPIDinfo; }
  Double_t          GetTheta()         const { return fTheta; }
  Double_t          GetX()             const { return fX; }
  Double_t          GetY()             const { return fY; }
  Double_t          GetX( Double_t z ) const { return fX + z*fTheta; }
  Double_t          GetY( Double_t z ) const { return fY + z*fPhi; }

  Double_t          GetChi2()          const { return fChi2; }
  Int_t             GetNDoF()          const { return fNDoF; }

  Double_t          GetDX()            const { return fDX; }
  Double_t          GetDY()            const { return fDY; }
  Double_t          GetDTheta()        const { return fDTheta; }
  Double_t          GetDPhi()          const { return fDPhi; }
  Double_t          GetRX()            const { return fRX; }
  Double_t          GetRY()            const { return fRY; }
  Double_t          GetRTheta()        const { return fRTheta; }
  Double_t          GetRPhi()          const { return fRPhi; }
  Double_t          GetToff()          const { return fToff; }
  Double_t          GetTX()            const { return fTX; }
  Double_t          GetTY()            const { return fTY; }
  Double_t          GetTTheta()        const { return fTTheta; }
  Double_t          GetTPhi()          const { return fTPhi; }
  Double_t          GetDp()            const { return fDp; }
  Double_t          GetLabPx()         const { return fPvect.X(); }
  Double_t          GetLabPy()         const { return fPvect.Y(); }
  Double_t          GetLabPz()         const { return fPvect.Z(); }
  Double_t          GetVertexX()       const { return fVertex.X(); }
  Double_t          GetVertexY()       const { return fVertex.Y(); }
  Double_t          GetVertexZ()       const { return fVertex.Z(); }
  Double_t          GetAltVertexX()    const { return fAltVertex.X(); }
  Double_t          GetAltVertexY()    const { return fAltVertex.Y(); }
  Double_t          GetAltVertexZ()    const { return fAltVertex.Z(); }
  Double_t          GetPathLen()       const { return fPathl; }

  TVector3&         GetPvect()               { return fPvect; }
  TVector3&         GetVertex()              { return fVertex; }
  TVector3&         GetVertexError()         { return fVertexError; }
  TVector3&         GetAltVertex()           { return fAltVertex; }
  TVector3&         GetAltVertexError()      { return fAltVertexError; }

  Double_t          GetTime()          const { return fTime; } // at refplane (s)
  Double_t          GetdTime()         const { return fdTime; } // (s)
  Double_t          GetBeta()          const { return fBeta; } // from scint.
  Double_t          GetdBeta()         const { return fdBeta; }
  Double_t          GetDedx()          const { return fDedx; }
  Double_t          GetEnergy()        const { return fEnergy; }
  Int_t             GetNPMT()          const { return fNPMT; }
  Double_t          GetBetaChi2()      const { return fBetaChi2; }
  Double_t          GetFPTime()        const { return fFPTime; }
  Int_t             GetGoodPlane3()    const { return fGoodPlane3; }
  Int_t             GetGoodPlane4()    const { return fGoodPlane4; }
  
  bool              HasDet()           const { return (fType&kHasDet); }
  bool              HasFP()            const { return (fType&kHasFP); }
  bool              HasRot()           const { return (fType&kHasRot); }
  bool              HasTarget()        const { return (fType&kHasTarget); }
  bool              HasVertex()        const { return (fType&kHasVertex); }

  void              Print( Option_t* opt="" ) const;

  void              Set( Double_t x, Double_t y, Double_t theta, Double_t phi )
  { fX = x; fY = y; fTheta = theta; fPhi = phi; fType |= kHasFP; }
  void              SetR( Double_t x, Double_t y,
			  Double_t theta, Double_t phi );
  void              SetD( Double_t x, Double_t y,
			  Double_t theta, Double_t phi );
  void              SetTarget( Double_t x, Double_t y,
			       Double_t theta, Double_t phi );
  
  void              SetPathLen( Double_t pathl ) { fPathl = pathl; /* meters */ }
  void              SetTime( Double_t time )     { fTime = time; /* seconds */ }
  void              SetToff( Double_t toff )     { fToff = toff; /* seconds */ }
  void              SetdTime( Double_t dt )      { fdTime = dt; /* seconds */ }
  void              SetBeta( Double_t beta )     { fBeta = beta; }
  void              SetdBeta( Double_t db )      { fdBeta = db; }
  void              SetDedx(Double_t dedx)       { fDedx = dedx; }
  void              SetEnergy(Double_t energy)   { fEnergy = energy; }
  void              SetNPMT(Int_t npmt)          { fNPMT = npmt; }
  void              SetBetaChi2(Double_t betachi2) { fBetaChi2 = betachi2; }
  void              SetFPTime(Double_t fptime)     { fFPTime = fptime; }
  void              SetGoodPlane3(Int_t gdplane3)  { fGoodPlane3 = gdplane3; }
  void              SetGoodPlane4(Int_t gdplane4)  { fGoodPlane4 = gdplane4; }

  void              SetClustNums(THaVDCPointPair* thePair);
  void              SetError( Double_t Error) { fError = Error; }
  void              SetChi2( Double_t chi2, Int_t ndof ) { fChi2=chi2; fNDoF=ndof; }

  void              SetID( THaTrackID* id )   { fID   = id; }
  void              SetFlag( UInt_t flag )    { fFlag = flag; }
  void              SetType( UInt_t flag )    { fType = flag; }
  void              SetMomentum( Double_t p ) { fP    = p; }
  void              SetDp( Double_t dp )      { fDp   = dp; }
  void              SetTrkNum( Int_t n )      { fTrkNum = n; }
  void              SetCreator( THaTrackingDetector* d ) { fCreator = d; }
  void              SetIndex( Int_t idx )                { fIndex   = idx; }
  void              SetPIDinfo( THaPIDinfo* pid )        { fPIDinfo = pid; }
  void              SetPvect( const TVector3& pvect )    { fPvect   = pvect; }
  void              SetVertex( const TVector3& vert )
  { fVertex = vert; fType |= kHasVertex; }
  void              SetVertex( Double_t x, Double_t y, Double_t z )
  { fVertex.SetXYZ( x, y, z ); fType |= kHasVertex; }
  void              SetVertexError( const TVector3& err )
  { fVertexError = err; }
  void              SetVertexError( Double_t x, Double_t y, Double_t z )
  { fVertexError.SetXYZ( x, y, z ); }
  void              SetAltVertex( const TVector3& vert )
  { fAltVertex = vert; fType |= kHasVertex; }
  void              SetAltVertex( Double_t x, Double_t y, Double_t z )
  { fAltVertex.SetXYZ( x, y, z ); fType |= kHasVertex; }
  void              SetAltVertexError( const TVector3& err )
  { fAltVertexError = err; }
  void              SetAltVertexError( Double_t x, Double_t y, Double_t z )
  { fAltVertexError.SetXYZ( x, y, z ); }
  
  void              SetUVXY(Double_t UV12X, Double_t UV12Y, Double_t UV12PX, Double_t UV12PY, Double_t UV21X, Double_t UV21Y, Double_t UV21PX, Double_t UV21PY)
  { fUV12X = UV12X;
    fUV12Y = UV12Y;
    fUV12PX = UV12PX;
    fUV12PY = UV12PY;
    fUV21X = UV21X;
    fUV21Y = UV21Y;
    fUV21PX = UV21PX;
    fUV21PY = UV21PY;
  }

  void              SetUVTP(Double_t UV12Theta, Double_t UV12Phi, Double_t UV21Theta, Double_t UV21Phi)
  { fUV12Theta = UV12Theta;
    fUV12Phi = UV12Phi;
    fUV21Theta = UV21Theta;
    fUV21Phi = UV21Phi;
  }

  virtual Bool_t    IsSortable() const { return kTRUE; }
  virtual Int_t	    Compare(const TObject* obj) const;

protected:

  enum { kMAXCL = 4 };

  // Focal plane coordinates (TRANSPORT system projected to z=0)
  Double_t          fX;              // x position in TRANSPORT plane (m)
  Double_t          fY;              // y position in TRANSPORT plane (m)
  Double_t          fTheta;          // Tangent of TRANSPORT Theta (x')
  Double_t          fPhi;            // Tangent of TRANSPORT Phi (y')
  Double_t          fP;              // Track momentum (GeV)

  // coordinates in the detector system
  Double_t fDX;     // x position in DCS
  Double_t fDY;     // y position in DCS
  Double_t fDTheta; // Tangent of DCS Theta
  Double_t fDPhi;   // Tangent of DCS Phi

  // coordinates in the rotated TRANSPORT system
  Double_t fRX;     // x position in focal plane (m)
  Double_t fRY;     // y position in focal plane (m)
  Double_t fRTheta; // Tangent of TRANSPORT Theta (x')
  Double_t fRPhi;   // Tangent of TRANSPORT Phi (y')

  // reconstructed coordinates at the target in the target TRANSPORT system (z=0)
  Double_t fTX;     // x position at target (m)
  Double_t fTY;     // y position at target (m)
  Double_t fTTheta; // Tangent of TRANSPORT Theta (out-of-plane angle) at target
  Double_t fTPhi;   // Tangent of TRANSPORT Phi (in-plane angle) at target
  Double_t fDp;     // dp/p_center -- fractional change in momentum

  TVector3          fPvect;  // Momentum vector at target in lab system (GeV)
  TVector3          fVertex; // Vertex location in lab (m) valid if fHasVertex
  TVector3          fVertexError; // Uncertainties in fVertex coordinates.
  TVector3          fAltVertex; // Alternative Vertex location in lab (m); calculated using reconstructed tg_y directly 
  TVector3          fAltVertexError; // Uncertainties in fAltVertex coordinates.

  Double_t          fPathl;  // pathlength from target (z=0) (meters)

  Double_t          fToff;   // average of timing offset in 4 planes (s)
  Double_t          fTime;   // time of track at focal plane (s)
  Double_t          fdTime;  // uncertainty in fTime
  Double_t          fBeta;   // beta of track
  Double_t          fdBeta;  // uncertainty in fBeta

  // Status variables and objects related to this track
  Int_t             fNclusters;      //! Number of clusters
  THaCluster*       fClusters[kMAXCL]; //! Clusters of this track
  THaPIDinfo*       fPIDinfo;        //! Particle ID information for this track
  THaTrackingDetector* fCreator;     //! Detector creating this track
  Int_t             fIndex;          // Track index (-1 = none, 0 = first, etc.)
  Int_t             fTrkNum;         // Track number (0 = unassigned)
  Int_t             fU1Num;           // Num of clust in U1 plane (-1 = unused)
  Int_t             fV1Num;           // Num of clust in V1 plane (-1 = unused)
  Int_t             fU2Num;           // Num of clust in U2 plane (-1 = unused)
  Int_t             fV2Num;           // Num of clust in V2 plane (-1 = unused)
  Double_t          fError;

  THaTrackID*       fID;     //! Track identifier
  UInt_t            fFlag;   // General status flag (for use by tracking det.)
  UInt_t            fType;   // Flag indicating which vectors reconstructed

  Double_t          fChi2;   // goodness of track fit
  Int_t             fNDoF;   // number of hits on the track contributing to chi2

  Double_t          fDedx;   // dEdX from hodoscopes
  Double_t          fEnergy; // Energy from calorimeter
  // Needed for "prune" select best track method
  Int_t             fNPMT;       // Number of PMTs hit in track
  Double_t          fBetaChi2;   // (reduced) Chisq of fit on Beta
  Double_t          fFPTime;     // Focal Plane time (same as fTime?)
  Int_t             fGoodPlane3; // Track hit a plane 3 paddle
  Int_t             fGoodPlane4; // Track hit a plane 4 paddle
  
  Double_t          fUV12X;     // X first chamber
  Double_t          fUV12Y;     // Y first chamber
  Double_t          fUV12PX;    // X in first chamber projected from second chamber
  Double_t          fUV12PY;    // Y in first chamber projected from second chamber
  Double_t          fUV12Theta; // Theta in first chamber
  Double_t          fUV12Phi;   // Phi in first chamber

  Double_t          fUV21X;     // X second chamber
  Double_t          fUV21Y;     // Y second chamber
  Double_t          fUV21PX;    // X in second chamber projected from first chamber
  Double_t          fUV21PY;    // Y in second chamber projected from first chamber
  Double_t          fUV21Theta; // Theta in secopnd chamber
  Double_t          fUV21Phi;   // Phi in second chamber


  
  static const Double_t kBig;

  ClassDef(THaTrack,5)       // A generic particle track
};

//__________________ inlines __________________________________________________
inline
void THaTrack::SetD( Double_t x, Double_t y, Double_t theta, Double_t phi )
{
  // set the coordinates in the rotated focal-plane frame
  fDX = x;
  fDY = y;
  fDTheta = theta;
  fDPhi = phi;
  fType |= kHasDet;
}

//_____________________________________________________________________________
inline
void THaTrack::SetR( Double_t x, Double_t y, Double_t theta, Double_t phi )
{
  // set the coordinates in the rotated focal-plane frame
  fRX = x;
  fRY = y;
  fRTheta = theta;
  fRPhi = phi;
  fType |= kHasRot;
}

//_____________________________________________________________________________
inline
void THaTrack::SetTarget( Double_t x, Double_t y, Double_t theta, Double_t phi )
{
  // set the coordinates in the target frame
  fTX = x;
  fTY = y;
  fTTheta = theta;
  fTPhi = phi;
  fType |= kHasTarget;
}

#endif
