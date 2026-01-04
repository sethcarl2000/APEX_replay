#ifndef Podd_TOpticsModel_h_
#define Podd_TOpticsModel_h_

//////////////////////////////////////////////////////////////////////////
//
// TOpticsModel
//
//////////////////////////////////////////////////////////////////////////

//#include "THaPhysicsModule.h"
#include <vector>
#include <map>
#include <memory>
#include <array>
#include "ROOT/RVec.hxx"
#include "TString.h" 
#include "TObject.h"
#include "ApexUtils.h"
#include "TNPoly.h"
#include "TXMap.h"
#include "TVector3.h"
#include <fstream> 
#include "OpticsVectors.h"
#include "OpticsPolynomials.h"
#include "TOpticsRays.h"

//_________________________________________________________________________________
class TOpticsModel : public TObject {
  
public:
  
  //this class will be ported to the apex_lib library, and will be responsible
  // for handling all optics calculations.
  TOpticsModel() {};
  TOpticsModel(const int  nDoF,
	       const bool isRHRS=true,
	       const int  order=0); 
  
  ~TOpticsModel();

  //functions to handle adding poly elements (done when filling in data from files)
  void Add_polyElement   (const ROOT::RVec<int>    &ePows,  
			  const ROOT::RVec<double> &ePoly,      
			  APEX::EFPCoordinate  coord);  
  
  void Add_fwdPolyElement(const ROOT::RVec<int>    &ePows,
			  //const ROOT::RVec<double> &ePoly,
			  double coeff,
			  APEX::ETargetCoord coord);
  
  
  std::unique_ptr<TXMap> Get_Xmap(const double x_fp) const; 
  
  //get x-tg, but only for a fixed value in the total span of possible values. 
  //ROOT::RVec<double> Get_Xtg(const ROOT::RVec<double> &X_fp) const; 
  TXtg Get_Xtg(const TXfp &Xfp) const; 
  
  
  TNPoly* Get_poly(const APEX::EFPCoordinate coord) { return &fPolys.at(coord); }

  
  int Parse_file         (TString inFile_path);
  int Parse_forwardTensor(TString inFile_path);
  
  int Parse_reverseTensor(TString inFile_path);
  

  std::unique_ptr<TOpticsRays> Get_opticsRays() const; 
  
  //TXfp Compute_Xfp(TXtg Xtg) const;
  
  //this is the 'r0' reaction point, used to find a first solution when converging to a
  // track.
  void Set_r0(double &x0, double &y0);
  
private:
  
  
  static const unsigned int fDoF_fp =4;
  static const unsigned int fDoF_tg =5; 
  
  //this should probably be read in from a db later on, or be able to change on a per-element basis
  const unsigned int f_elem_polyDegree=4; 

  unsigned int fNelems,fnDoF,fnDoF_fp,fPolySize;
  bool f_isRHRS; 
  unsigned int fOrder;

  std::array<TNPoly,TOpticsModel::fDoF_fp-1> fPolys;
  
  //these polynomials map from fp-coords to target-coords, but only for a limited
  // range! this is why the 'solver' is necessary.
  std::array<TFPoly, 2> fPolys_forward;
  TVector3 f_R0; //this is used as the (approximate) react-point for the forward poly
  
  std::array<TRPoly, 4> fPolys_reverse;
  bool f_init_reversePolys; 
  
  //z-position of the sieve-plane, in target-coordinates
  double fZ_sv; 
  
  
  const TString className="TOpticsModel";
  
  
  const std::map<const APEX::EFPCoordinate,const TString> fCoordName = {
    {APEX::kY,    "Y_fp"},
    {APEX::kDxdz, "dx/dz_fp"},
    {APEX::kDydz, "dy/dz_fp"}
  };
  
  ClassDef(TOpticsModel,0); 
};

#endif
