#ifndef Podd_TapexReactVertex_h_
#define Podd_TapexReactVertex_h_

//////////////////////////////////////////////////////////////////////////
//
// TapexReactVertex
// Handles computation of the react vertex. 
//
//////////////////////////////////////////////////////////////////////////

//#include "THaPhysicsModule.h"
#include "TObject.h"
#include "TVector3.h"
#include "TVector2.h"
#include "RMatrixD.h"
#include "TString.h" 
#include <ROOT/RDataFrame.hxx>
#include <vector> 
#include <map> 
#include <functional> 

//_________________________________________________________________________________
class TapexReactVertex { 

public:
  
  TapexReactVertex() {};
  
  enum EOpticsTarget { kNone=0, kV1,kV2,kV3, kH1,kH2,kH3,kH4 }; 
  
  TapexReactVertex(bool isRHRS,
	       TString path_decode_epics,
         ROOT::RDF::RNode df,
	       TString target="");

  virtual ~TapexReactVertex() {}; 

  /// @brief computes react vertex x & y 
  /// @param current_x raw current from Raster2 (x): [R/L]rb.Raster2.rawcur.x 
  /// @param current_y raw current from Raster2 (y): [R/L]rb.Raster2.rawcur.y 
  /// @param y_BPM average of the 2 BPM readouts for this event
  /// @return a 3-vector whose x & y components are the best-estimate for x_hcs and y_hcs (m)
  TVector3 Compute_reactVertex(double current_x, double current_y, double y_BPM) const; 
  
  TVector2 Get_beamCenter() const; 

  TVector2 Get_rasterAmplitude() const { return fRaster_amplitude; } 
  
  struct OpticsWire_t {
    TString name;
    bool isVertical; 
    double x,y,z;
  };  

  //not sure if i'll ever need these methods, but they might be useful 
  bool IsWireMode() const { return f_isWireMode; }
  OpticsWire_t Get_wire() const; 

  bool IsRHRS() const { return fis_RHRS; }
  
private:

  TString fTargetName;
  bool f_isWireMode;  //this is true when a valid optic-wire name is passed as 'target' 
  bool f_hasData;    //this is true when no valid file-name is passed. 
  
  OpticsWire_t fWire;

  const std::map<TString,TapexReactVertex::OpticsWire_t> fWireMap = {
    {"V1", {.name="V1", .isVertical=true,  .x=-3.225e-3, .y= 0e-3,     .z=-196.214e-3}},
    {"V2", {.name="V2", .isVertical=true,  .x=-0.725e-3, .y= 0e-3,     .z=   3.786e-3}},
    {"V3", {.name="V3", .isVertical=true,  .x= 1.725e-3, .y= 0e-3,     .z= 203.786e-3}},
    {"H1", {.name="H1", .isVertical=false, .x= 0e-3,     .y=-7.660e-3, .z=-241.249e-3}},
    {"H2", {.name="H2", .isVertical=false, .x= 0e-3,     .y=-2.676e-3, .z= -91.204e-3}},
    {"H3", {.name="H3", .isVertical=false, .x= 0e-3,     .y= 2.324e-3, .z= 108.750e-3}},
    {"H4", {.name="H4", .isVertical=false, .x= 0e-3,     .y= 7.328e-3, .z= 258.746e-3}}
  }; 
    
  TVector2 fBeamCenter;  //center of beam at target (z_HCS = 0) 
  TVector2 fBeam_dXdz;

  //slope of the beam in the x-direction (HCS)
  double fBeam_dxdz;
  
  //slope of the beam in the x-direction (HCS)
  double fBeam_dydz;

  //the beam's x_hcs position at the z_hcs = 0 plane (m)
  double fBeam_x0; 
  
  //the beam's x_hcs position at the z_hcs = 0 plane (m)
  double fBeam_y0; 

  //this function performs the raster phase correction.
  // signature: (y_rast, y_BPM) -> y_hcs_phase_correction (m)
  std::function<double(double,double)> fRasterPhaseCorrection; 

  //the computed phase correction to be applied to y_hcs (m)
  double fYhcs_phase_correction;

  //mean values of raw raster
  double fMeanCurrent_x, fMeanCurrent_y; 

  bool fis_RHRS; 

  RMatrixD fMatrix_rast; //channel (dc) to meter (dX) conversion
  
  RMatrixD fMatrix_BPMA; //correction matrix for BPM readouts
  ROOT::RVec<double> fR0_BPMA, fR_BPMA;
  
  RMatrixD fMatrix_BPMB; //correction matrix for BPM readouts
  ROOT::RVec<double> fR0_BPMB, fR_BPMB; 

  
  TVector2 fRast_avg; 

  //the x-y span of the raster, throughout the run (in HCS)
  TVector2 fRaster_amplitude; 
  
};
//_________________________________________________________________________________

#endif
