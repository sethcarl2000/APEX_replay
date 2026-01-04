#ifndef Podd_ApexUtils_h_
#define Podd_ApexUtils_h_

//////////////////////////////////////////////////////////////////////////
//
// ApexUtils
//
// Some misc. utilities which will be used by various scripts
//
//////////////////////////////////////////////////////////////////////////

#include <limits>
#include "TVectorD.h"
#include "TMatrixD.h"
#include "ROOT/RVec.hxx"
#include "TGraph.h"
#include "TFitResult.h"
#include "TVector3.h"
#include "TH1.h"
#include "TPad.h"
#include "OpticsVectors.h"

namespace APEX {
  
  enum EFPCoordinate { kY=0, kDxdz, kDydz }; 
  enum ETargetCoord { kX_tg=0, kY_tg, kDxdz_tg, kDydz_tg };
  
  //constexpr double min_determinant = 1e-200; 
  TVectorD SolveLinSystem(TMatrixD A, TVectorD B);

  //i tried to templatize this but it was too much of a pain
  TVectorD TVectorD_cast(const ROOT::RVec<double>  &inVec);
  TVectorD TVectorD_cast(const std::vector<double> &inVec);

  ROOT::RVec<double> RVecD_cast(const TVectorD &inVec); 

  template <typename T> T* FindTypeInPad(TPad *pad); 

  TFitResult* Draw_gausFit( TH1 *hist, 
			    double radius, 
			    double center=-1e30, 
			    int color=2 ); 
  
  namespace Coord { 
    enum ETransformType { kHall_to_Target=0, kTarget_to_Hall }; 
    
    void RotateTV3   (TVector3 &r, ETransformType type, bool isRHRS);
    void TransformTV3(TVector3 &r, ETransformType type, bool isRHRS);
    
    void Transform_Xtg(TXtg &Xtg,
		       ETransformType type,
		       bool isRHRS); 
  };

  void SetTGraphBounds(TGraph *g, const double xmin, const double ymin, const double xmax, const double ymax); 

  //computes polynomial
  //template for vector-objects
  double ComputePol(const ROOT::RVec<double>  *a,     const double x);
  double ComputePol(const std::vector<double> *a,     const double x);
  double ComputePol(const double *a,                  const double x,  const int order); 
  double ComputePol(const ROOT::RVec<double>  &a,     const double x); 

  bool hasNAN_RVecD(const ROOT::RVec<double> &v);

  
  
  //some RVecD Handlers
  double Length( const ROOT::RVec<double> &r1 );
  double Dot   ( const ROOT::RVec<double> &r1, const ROOT::RVec<double> &r2 ); 
  ROOT::RVec<double> Unit( const ROOT::RVec<double> &v ); 
  
  
  //some NaN declarations
  static constexpr double kNAN_double = std::numeric_limits<double>::quiet_NaN();  
  
  namespace Graphic {
    void VLine(TH1 *hist, double x, unsigned int color=1, unsigned int style=1); 
    void HLine(TH1 *hist, double y, unsigned int color=1, unsigned int style=1); 
  };
  
}; 
#endif
