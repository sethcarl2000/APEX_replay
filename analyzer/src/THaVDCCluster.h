#ifndef Podd_THaVDCCluster_h_
#define Podd_THaVDCCluster_h_

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// THaVDCCluster                                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include <utility>
#include <vector>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"


class THaVDCHit;
class THaVDCPlane;
class THaVDCPointPair;
class THaTrack;

namespace VDC {
  struct FitCoord_t {
    FitCoord_t( Double_t _x, Double_t _y, Double_t _w = 1.0, Int_t _s = 1 )
      : x(_x), y(_y), w(_w), s(_s) {}
    Double_t x, y, w;
    Int_t s;
  };
  struct FitCoordT_t {
    FitCoordT_t( THaVDCHit* _hit,Double_t _x, Double_t _t, Double_t _w = 1.0, Int_t _s = 1 )
      : hit(_hit), x(_x), t(_t), w(_w), s(_s) {}
    THaVDCHit* hit;
    Double_t x, t, w;
    Int_t s;    
    };

  //  FitCoordT_t( Double_t _x, Double_t _t, Double_t _w = 1.0, Int_t _s = 1 )
  //  VDC::Vhit_t    fHits;              // Hits associated w/this cluster
  
  typedef std::pair<Double_t,Int_t>  chi2_t;
  typedef THaVDCPointPair VDCpp_t;
  typedef std::vector<THaVDCHit*> Vhit_t;
  typedef std::vector<FitCoord_t> Vcoord_t;
  typedef std::vector<FitCoordT_t> VcoordT_t;

  extern const Double_t kBig;

  inline chi2_t operator+( chi2_t a, const chi2_t& b ) {
    a.first  += b.first;
    a.second += b.second;
    return a;
  };
}

class THaVDCCluster : public TObject {

public:

  THaVDCCluster( THaVDCPlane* owner = 0 );
  //  virtual ~THaVDCCluster() {}

  ~THaVDCCluster() {
    delete fmin;
    delete func;
  }

  //  enum EMode { kSimple, kWeighted, kT0 };
  enum EMode { kSimple, kWeighted};
  enum EFitMode {kLinear, kTwoParam, kThreeParam, kT0}; // type of fit to be used for cluster fitting: kLinear is two parameter fit done via standard formulae for linear regression, kTwoParam is a two parameter fit performed via minuit, kThreeParam is a three parameter fit performed via minuit, kT0 is a three parameter if performed via standard formulae
  enum EConvMode {kNormal, kToffapp}; //  controls intial conversion from time to distance: kNormal uses unadjusted time, kToffapp takes into account approximation of t0

  virtual void   AddHit( THaVDCHit* hit );
  virtual void   EstTrackParameters();
  virtual void   ConvertTimeToDist();
  virtual void   ConvertTimeToDist(EConvMode convFit);
  virtual void   FitTrack( EMode mode = kWeighted,  EFitMode modeFit = kLinear);
  virtual void   ClearFit();
  virtual void   CalcChisquare(Double_t& chi2, Int_t& nhits) const;
  VDC::chi2_t    CalcDist();    // calculate global track to wire distances

  // TObject functions redefined
  virtual void   Clear( Option_t* opt="" );
  virtual Int_t  Compare( const TObject* obj ) const;
  virtual Bool_t IsSortable() const        { return kTRUE; }
  virtual void   Print( Option_t* opt="" ) const;

  //Get and Set Functions
  THaVDCHit*     GetHit(Int_t i)     const { return fHits[i]; }
  THaVDCPlane*   GetPlane()          const { return fPlane; }
  Int_t          GetSize ()          const { return fHits.size(); }
  Double_t       GetSlope()          const { return fSlope; }
  Double_t       GetLocalSlope()     const { return fLocalSlope; }
  Double_t       GetSigmaSlope()     const { return fSigmaSlope; }
  Double_t       GetIntercept()      const { return fInt; }
  Double_t       GetSigmaIntercept() const { return fSigmaInt; }
  THaVDCHit*     GetPivot()          const { return fPivot; }
  Int_t          GetPivotWireNum()   const;
  Double_t       GetTimeCorrection() const { return fTimeCorrection; }
  Double_t       GetT0()             const { return fT0; }
  VDC::VDCpp_t*  GetPointPair()      const { return fPointPair; }
  THaTrack*      GetTrack()          const { return fTrack; }
  Int_t          GetTrackIndex()     const;
  Int_t          GetTrkNum()         const { return fTrkNum; }
  Int_t          GetClsNum()         const { return fClsNum; }
  Double_t       GetSigmaT0()        const { return fSigmaT0; }
  Double_t       GetChi2()           const { return fChi2; }
  Double_t       GetNDoF()           const { return fNDoF; }
  Bool_t         IsFitOK()           const { return fFitOK; }
  Bool_t         IsUsed()            const { return (fTrack != 0); }

  void           SetPlane( THaVDCPlane* plane )     { fPlane = plane; }
  void           SetIntercept( Double_t intercept ) { fInt = intercept; }
  void           SetSlope( Double_t slope)          { fSlope = slope;}
  void           SetPivot( THaVDCHit* piv)          { fPivot = piv; }
  void           SetTimeCorrection( Double_t dt )   { fTimeCorrection = dt; }
  void           SetPointPair( VDC::VDCpp_t* pp )   { fPointPair = pp; }
  void           SetTrack( THaTrack* track );
  void           SetClsNum( Int_t clsnum)           { fClsNum = clsnum; }
  void           SetT0_app_adj( Double_t t0_app_adj) {fT0_app_adj = t0_app_adj; }
  

protected:

  // minimiser
  ROOT::Minuit2::Minuit2Minimizer* fmin;

  ROOT::Math::Functor* func;

  const Int_t NPara = 3;
  
  VDC::Vhit_t    fHits;              // Hits associated w/this cluster
  THaVDCPlane*   fPlane;             // Plane the cluster belongs to
  VDC::VDCpp_t*  fPointPair;         // Lower/upper combo we're assigned to
  THaTrack*      fTrack;             // Track the cluster belongs to
  Int_t          fTrkNum;            // Number of the track using this cluster
  Int_t          fClsNum;            // Number of cluster (-1 = unused)


  //Track Parameters
  Double_t       fSlope;             // Current best estimate of actual slope
  Double_t       fLocalSlope;        // Fitted slope, from FitTrack()
  Double_t       fSigmaSlope;        // Error estimate of fLocalSlope from fit
  Double_t       fInt, fSigmaInt;    // Intercept and error estimate
  Double_t       fLocalInt, fSigmaLocalInt;    // Local Intercept and error estimate    
  Double_t       fT0 = 0, fSigmaT0;      // Fitted common timing offset and error
  Double_t       fT0_app;      // estimate common timing offset
  Double_t       fT0_app_adj = 6.7e-9;  // adjustment between approximation and 'true' timing offset (imperical)
  Double_t       fT0_fake = 0.0;           // fake offset added at start of fitting and removed to get result (200 ns)
  THaVDCHit*     fPivot;             // Pivot - hit with smallest drift time
  //FIXME: in the code, this is used as a distance correction!!
  Double_t       fTimeCorrection;    // correction to be applied when fitting
				     // drift times
  Bool_t         fFitOK;             // Flag indicating that fit results valid
  Double_t       fChi2;              // chi2 for the cluster (using fSlope)
  Double_t       fNDoF;              // NDoF in local chi2 calculation
  Int_t          fClsBeg;	     // Starting wire number
  Int_t          fClsEnd;            // Ending wire number

  // Workspace for fitting routines
  VDC::Vcoord_t  fCoord;             // coordinates to be fit
  VDC::VcoordT_t  fCoordT;            // alternative coordinates to be fit (time for each hit inlclded)

  void   CalcLocalDist();     // calculate the local track to wire distances

  void   FitSimpleTrack( Bool_t weighted = false );
  void   FitNLTrack();        // Non-linear 3-parameter fit

  VDC::chi2_t CalcChisquare( Double_t slope, Double_t icpt, Double_t d0 ) const;
  //  VDC::chi2_t CalcChisquareTwoParam( Double_t slope, Double_t icpt, Double_t d0 ) const;
  void   DoCalcChisquare( Double_t& chi2, Int_t& nhits,
			  Double_t slope, bool do_print = false ) const;
  
  void   FitTwoParamTrack( Bool_t weighted = false );
  //  VDC::chi2_t CalcChisquareTwoParam( Double_t slope, Double_t icpt, Double_t d0 ) const;
  VDC::chi2_t TwoParamFit( Double_t& slope, Double_t& icpt);
  Double_t fcn_2P(const Double_t* par);
  
  void   FitThreeParamTrack( Bool_t weighted = false );
  //  VDC::chi2_t CalcChisquareTwoParam( Double_t slope, Double_t icpt, Double_t d0 ) const;
  VDC::chi2_t ThreeParamFit( Double_t& slope, Double_t& icpt, Double_t& t, Double_t& slopeErr, Double_t& icptErr, Double_t& tError);
  Double_t fcn_3P(const Double_t* par);
  
  void   Linear3DFit( Double_t& slope, Double_t& icpt, Double_t& d0 ) const;
  Int_t  LinearClusterFitWithT0();

  ClassDef(THaVDCCluster,0)          // A group of VDC hits
};

//////////////////////////////////////////////////////////////////////////////

#endif
