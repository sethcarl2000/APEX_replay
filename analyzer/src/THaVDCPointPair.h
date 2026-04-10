#ifndef Podd_THaVDCPointPair_h_
#define Podd_THaVDCPointPair_h_

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// THaVDCPointPair                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "THaVDCCluster.h"   // for chi2_t

class THaVDCPoint;
class THaTrack;

// X and Y Differences between point pair in both chambers and timing offset difference
struct PointPair_Diff {
  Double_t X1Diff; // X difference in lower chamber
  Double_t Y1Diff; // Y difference in lower chamber
  Double_t X2Diff; // X difference in upper chamber
  Double_t Y2Diff; // Y difference in upper chamber
  Double_t t0Diff;  // t0 difference

};


class THaVDCPointPair : public TObject {

public:
  THaVDCPointPair( THaVDCPoint* lp, THaVDCPoint* up, Double_t spacing, Double_t XRes, Double_t YRes)
    : fLowerPoint(lp), fUpperPoint(up), fSpacing(spacing), fXRes(XRes), fYRes(YRes),fError(1e38),
      fStatus(0) {}
  virtual ~THaVDCPointPair() {}

  void            Analyze();
  void            Associate( THaTrack* track );
  VDC::chi2_t     CalcChi2() const;
  virtual Int_t   Compare( const TObject* ) const;
  Double_t        GetError()   const { return fError; }
  THaVDCPoint*    GetLower()   const { return fLowerPoint; }
  THaVDCPoint*    GetUpper()   const { return fUpperPoint; }
  Double_t        GetSpacing() const { return fSpacing; }
  Int_t           GetStatus()  const { return fStatus; }
  THaTrack*       GetTrack()   const;
  Bool_t          HasUsedCluster() const;
  virtual Bool_t  IsSortable() const { return kTRUE; }
  virtual void    Print( Option_t* opt="" ) const;
  void            Release();
  void            SetStatus( Int_t i ) { fStatus = i; }
  void            Use();


  Double_t CalcError(Bool_t &Pass);

  
  /*static PointPair_Diff CalcErrorEst( THaVDCPoint* lowerPoint,
			     THaVDCPoint* upperPoint,
				Double_t spacing);
  */

  static void CalcXY( THaVDCPoint* lowerPoint,
		      THaVDCPoint* upperPoint,
		      Double_t spacing, Double_t &UV12X, Double_t &UV12Y, Double_t &UV12PX, Double_t &UV12PY, Double_t &UV21X, Double_t &UV21Y, Double_t &UV21PX, Double_t &UV21PY);


  static void GetTP(THaVDCPoint* here, THaVDCPoint* there, Double_t &UV12Theta, Double_t &UV12Phi,  Double_t &UV21Theta, Double_t &UV21Phi);

  static Double_t GetProjectedDistance( THaVDCPoint* here,
					THaVDCPoint* there,
					Double_t spacing, Bool_t &Pass, Double_t XRes = 0.009, Double_t YRes = 0.006); // default values for x and y resolution here established from data

  static Double_t GetProjectedDistance( THaVDCPoint* here,
					THaVDCPoint* there,
					Double_t spacing, Double_t XRes = 0.009, Double_t YRes = 0.006); // default values for x and y resolution here established from data


  static void GetProjectedXY( THaVDCPoint* here,
			      THaVDCPoint* there,
			      Double_t spacing,
			      Double_t &X,
			      Double_t &Y,
			      Double_t &PX,
			      Double_t &PY);

  // get sum of t0 differences between all planes 
  static Double_t GetTimeOffsetDifference(THaVDCPoint* here,
					  THaVDCPoint* there);

  // third parameter to check if all to diff combinations fall within 3 sigma of each other
  static Double_t GetTimeOffsetDifference(THaVDCPoint* here,
					  THaVDCPoint* there, Bool_t &Pass);


protected:

  THaVDCPoint*    fLowerPoint;  // Lower UV point
  THaVDCPoint*    fUpperPoint;  // Upper UV point
  Double_t        fSpacing;     // Spacing between lower and upper chambers [m]
  Double_t        fXRes;        // Resolution of X-PX [m]
  Double_t        fYRes;        // Resolution of Y-PY [m]
  Double_t        fError;       // Goodness of match between the points
  Bool_t        fPass;        // check if Bool passes all cuts
  Int_t           fStatus;      // Status flag

private:
  THaVDCPointPair();

  ClassDef(THaVDCPointPair,0)     // Pair of lower/upper VDC points
};

//////////////////////////////////////////////////////////////////////////////

#endif
