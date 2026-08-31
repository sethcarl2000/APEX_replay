//*-- Author :    Seth Hall   8-Oct-24

//////////////////////////////////////////////////////////////////////////
//     
// ReactVertex
//
// Computes the reaction-vertex 
//     
//////////////////////////////////////////////////////////////////////////

// APEX 
#include <APEX/ReactVertex.h>
#include <RMatrixD.h> 
// ROOT
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TVector2.h>
// stdlib 
#include <memory> 
#include <string> 
#include <math.h> 
#include <stdexcept> 

using namespace std;
//using namespace APEX/; 

using RVecD = ROOT::RVec<double>; 

namespace {
  /// z-positions of BPMs (relative to apex scattering chamber)
  
  /// z-position of BPMA (m)
  const double kBPMA_z = -7.354 + 1.0537;

  /// z-position of BPMB (m)
  const double kBPMB_z = -2.215 + 1.0537;


  //when deciding where to make vertex cuts, reject BPM reports outside of this range
  constexpr double kMax_y_BPM = +6.00e-3;
  constexpr double kMin_y_BPM =  0.00e-3;

  constexpr double kBPM_buffer = 1e-3; 
}

namespace APEX 
{

//_____________________________________________________________________________
ReactVertex::ReactVertex(
			 replay::EArmMode arm_mode,
  TString path_decode_epics,
  ROOT::RDF::RNode dT, 
  TString target)
  : fis_RHRS{(arm_mode == replay::kRHRS)},
    fTargetName(target),
    f_isWireMode(false),
    f_hasData(false), 
    fRaster_amplitude(TVector2(0,0))
{ 
  //these constants are hard-coded for now... it's probably better to have them loaded in from a file...
  //the 'df' expects a node of type '
  
  //initialize the wire as null, unless its specified otherwise
  fWire = OpticsWire_t({.name="null",.isVertical=false,.x=0,.y=0,.z=0}); 
  
  if (target != "") {
    //if the creator of this object told it that this is an optic-wire run (target != ""),
    // then we will search for whichever wire it wanted to use
    if (auto findWire = fWireMap.find(target); findWire == fWireMap.end()) {
      //wire-name given is invalid
      Warning("ReactVertex",
	      "Target-name passed: '%s', but this name does not match any wire-name. Proceeding without any optic-wire selected..",
	      target.Data());
    } else {
      //assign the wire to the matching name
      f_isWireMode=true;
      
      fWire = findWire->second; 
      Info("ReactVertex",
	   "Engaging optics-wire mode.. (Target name = %s)", fWire.name.Data() );
    }
  }

  //raster matrices 
  if (fis_RHRS) { //RHRS 
    fMatrix_rast = RMatrixD(2,2,{
       0.,         0.,
      -2.418e-8, -3.234e-7
    });
  } else      { //LHRS
    fMatrix_rast = RMatrixD(2,2,{
       0.,         0.,
      -3.500e-8, -4.533e-7
    });
  }
  
  //bpm matrix A
  fMatrix_BPMA = RMatrixD(2,2,{
    -1.044e0,  1.707e-1,
     2.449e-2, 1.132e0
  });
  
  //bpm offset A 
  fR0_BPMA     = {1.185e-3, 4.223e-4};

  //bpm matrix B 
  fMatrix_BPMB = RMatrixD(2,2,{
    -2.833e-1, -3.360e-2,
     2.781e-3,  8.627e-1
  });
  
  //bpm offset B 
  fR0_BPMB     = {1.314e-3, 1.713e-3}; 

  fR_BPMA = {0,0};
  fR_BPMB = {0,0};
  
  fBeam_dXdz  = TVector2(0,0);
  fBeamCenter = TVector2(0,0);
  
  fRast_avg   = TVector2(0,0);

  
  //do some basic checks before we start

  //you can launch ReactVertex in 'no-decode' mode, in which you're just using it
  // basically to store optics wire-data.
  if (path_decode_epics=="") return;
    
  auto file = unique_ptr<TFile>(new TFile(path_decode_epics.Data()));
  
  if (!file || file->IsZombie()) {
    Error("ReactVertex()", "Pointer to file \"%s\" invalid. Check path?",
	  path_decode_epics.Data());
    return;
  }
  
  
  if (!file->IsOpen() || file->IsZombie()) {
    Error("ReactVertex()", "File \"%s\" is zombie / is not open.", path_decode_epics.Data());
    return;
  }
    
  //now, check to make sure that the right trees are present
  if (!file->GetListOfKeys()->Contains("E")) {
    Error("ReactVertex",
	  "File \"%s\" does not contain Epics tree (E). React-vertex & raster set to 0,0",
	  path_decode_epics.Data()); 
    return;    
  } 
  
  ROOT::RDataFrame df_epics("E", path_decode_epics.Data()); 

  //only choose events with non-zero beam-current readings
  auto dE = df_epics.Filter([](double bcm){return bcm>1;}, {"hac_bcm_average"}); 

  const string bpma_x_epics{"IPM1H04A.XPOS"};
  const string bpma_y_epics{"IPM1H04A.YPOS"};

  const string bpmb_x_epics{"IPM1H04B.XPOS"};
  const string bpmb_y_epics{"IPM1H04B.YPOS"};

  //get the average epics BPM-reading for all events
  RVecD fR_BPMA_raw = { *dE.Mean(bpma_x_epics)/1e3, *dE.Mean(bpma_y_epics)/1e3 };
  RVecD fR_BPMB_raw = { *dE.Mean(bpmb_x_epics)/1e3, *dE.Mean(bpmb_y_epics)/1e3 };

  fR_BPMA = fMatrix_BPMA*fR_BPMA_raw + fR0_BPMA;
  fR_BPMB = fMatrix_BPMB*fR_BPMB_raw + fR0_BPMB;

  //compute the slope of the beam
  fBeam_dxdz = (fR_BPMB[0]-fR_BPMA[0])/(kBPMB_z-kBPMA_z); 
  fBeam_dydz = (fR_BPMB[1]-fR_BPMA[1])/(kBPMB_z-kBPMA_z); 

  //compute the beam position at z_hcs = 0; 
  fBeam_x0 = fR_BPMB[0] - fBeam_dxdz*kBPMB_z;
  fBeam_y0 = fR_BPMB[1] - fBeam_dydz*kBPMB_z;
  
  
  //now, check to make sure that the right trees are present
  if (!file->GetListOfKeys()->Contains("T")) {
    Error("ReactVertex",
	  "File \"%s\" does not contain CODA tree (T). Avg. Raster set to 0,0",
	  path_decode_epics.Data()); 
    return;    
  } 
  
  const string arm = fis_RHRS ? "R" : "L";
  const string raster_name = arm+"rb.Raster2";
  const string name_bpma = arm+"rb.BPMA"; 
  const string name_bpmb = arm+"rb.BPMB"; 
  

  const double max_current_x = *dT.Max(raster_name+".rawcur.x");
  const double min_current_x = *dT.Min(raster_name+".rawcur.x");
  
  const double max_current_y = *dT.Max(raster_name+".rawcur.y");
  const double min_current_y = *dT.Min(raster_name+".rawcur.y");

  fMeanCurrent_x = ( max_current_x + min_current_x )/2.;
  fMeanCurrent_y = ( max_current_y + min_current_y )/2.;

  const double y_hcs_min = fMatrix_rast.get(1,1)*(max_current_y - fMeanCurrent_y) + fBeam_y0; 
  const double y_hcs_max = fMatrix_rast.get(1,1)*(min_current_y - fMeanCurrent_y) + fBeam_y0; 

  //distance between min and max y_hcs for this run
  const double y_hcs_amplitude = std::fabs( y_hcs_max - y_hcs_min );

  auto dt_raster_param = dT

    .Define("raster_parameter", [max_current_y, min_current_y](double y_rawcur)
    {
      return (max_current_y-y_rawcur)/(max_current_y-min_current_y);
    }, {raster_name+".rawcur.y"}) 
    
    .Define("y_BPM", [](double a_y, double b_y)
    {
      return (a_y + b_y)/2; 
    }, {name_bpma+".y", name_bpmb+".y"})
    
    .Filter([](double y_BPM){ return (y_BPM > kMin_y_BPM) && (y_BPM < kMax_y_BPM);}, {"y_BPM"}); 


  //bpm-y value of events with highest raster (lowest y-hcs)
  const double min_rast_y_bpm = *dt_raster_param 
      //the events with the smallest raster; what is their average y-bpm value? 
      .Filter([](double rast_param){ return rast_param < 0.01; }, {"raster_parameter"})
      .Mean("y_BPM"); 

  const double max_rast_y_bpm = *dt_raster_param
      //the events with the smallest raster; what is their average y-bpm value? 
      .Filter([](double rast_param){ return rast_param > 0.99; }, {"raster_parameter"})
      .Mean("y_BPM"); 
  
  printf("<%s> min/max y_bpm (for phase offset)  %5.2f mm / %5.2f mm\n", __func__, min_rast_y_bpm*1e3, max_rast_y_bpm*1e3); 

  //we need to reverse-engineer some informaiton to 'fix' the y-raster. 
  //  this constant here (0.070302882883/2) is the measured ratio: (raster-timing-delay / 1-rast.-period).
  //  computed by looking at the graphs of h-wire runs, we can look at the 'peaks' when we plot y-rast vs. y-BPM  
  //                  /\<=|
  //      ___________/  \__________________________ <= bpm
  //                      |=>/\
  //      __________________/  \___________________ => bpm 
  //
  //                   |------|  <= 2*raster-timing-delay
  //
  //      |----------------------------------------| <= total raster amplitude
  //  
  const double y_hcs_phase_correction = y_hcs_amplitude * (fis_RHRS ? 0.0684620982311 : 0.070302882883) / 2.; 
  
  //____________________________________________________________________________________________
  fRasterPhaseCorrection = 
  [ y_hcs_phase_correction, 
    min_current_y,
    max_current_y, 
    min_rast_y_bpm, 
    max_rast_y_bpm](double current_y, double y_BPM)
  {
    double raster_param = (max_current_y - current_y) / (max_current_y - min_current_y); 

    double phase_line_value = min_rast_y_bpm + (max_rast_y_bpm - min_rast_y_bpm) * raster_param; 

    if (y_BPM > phase_line_value) {
      return -y_hcs_phase_correction;
    } else {
      return +y_hcs_phase_correction;
    }
  }; 
  //____________________________________________________________________________________________
  
  /*
  //now, we can compute the average raster position, which we will need later to
  // compute the beam-position on a per-event basis. 
  bool is_mt_enabled = ROOT::IsImplicitMTEnabled(); 
//  if (is_mt_enabled) ROOT::DisableImplicitMT(); 

  auto findAvgRast = [&dT](TString rastName)
  {
    double rastMin(1e30), rastMax(-1e30);

    auto avg = dT
      .Range(0,20e3)
      .Define("rast", [&rastMin,&rastMax](double rast)
      {
	if (rast < rastMin) rastMin=rast;
	if (rast > rastMax) rastMax=rast;
	return rast;
      }, {rastName.Data()})
      .Histo1D({"", "", 200, -1,-1}, "rast")->GetMean(); 

    std::array<double,2> rastSpan = {rastMin,rastMax}; 
    return rastSpan; 
  };
  //Now, find the average-raster values
  std::array<double,2> xRast = findAvgRast( TString(is_RHRS?"R":"L")+"rb.Raster2.rawcur.x" );
  std::array<double,2> yRast = findAvgRast( TString(is_RHRS?"R":"L")+"rb.Raster2.rawcur.y" );
  
  fRast_avg = TVector2( 0.5*(xRast[0]+xRast[1]),
			0.5*(yRast[0]+yRast[1]) );
  
  fRaster_amplitude = TVector2( xRast[1]-xRast[0],
				yRast[1]-yRast[0] ); 
			
    */ 
  Info("ReactVertex", "Done with initial react-point calculations");

  //reset this to its previous value
  //if (is_mt_enabled) ROOT::EnableImplicitMT(); 
  
  f_hasData = true; 
}
//_____________________________________________________________________________
TVector3 ReactVertex::Compute_reactVertex(double current_x, double current_y, double y_BPM) const
{
  //NOTE: this is given in hall-coordinates (HCS), in meters, from the APEX/
  // scattering-chamber center. 
  if (!f_hasData) {
    throw std::logic_error("<ReactVertex::Compute_reactVertex> "
	    "not created with valid run-data, cannot compute react-vertex."
    );
    return TVector3(0,0,0);
  }
  double x_hcs, y_hcs; 
  
  //this accounts for the fact that the y-raster was not calibrated correctly for these runs. (m) 
  const double y_correction = fis_RHRS ? -2.0559325e-3 : +1.72835e-3;

  //compute y_hcs using bpms, raster
  RVecD currents{ current_x - fMeanCurrent_x, current_y - fMeanCurrent_y };

  RVecD raster_offset = fMatrix_rast*currents; 

  x_hcs = raster_offset[0] + fBeam_x0; 
  y_hcs = raster_offset[1] + fBeam_y0; 
  
  //now, add corrections for y_hcs
  y_hcs += fRasterPhaseCorrection(current_y, y_BPM) + y_correction; 

  return TVector3( x_hcs, y_hcs, 0. );
  /* 
  //compute the beam's offset from the run-average position
  RVecD rast = {
    rastX - fRast_avg.X(),
    rastY - fRast_avg.Y()
  };
  
  RVecD beam_offset = fMatrix_rast * rast; 

  TVector3 react_point( fBeamCenter.X() + beam_offset[0],
			fBeamCenter.Y() + beam_offset[1],
			0. ); 

  //check if we're running in optics (wire) mode
  if (f_isWireMode) {

    //vertical or horizontal wire? 
    if (fWire.isVertical) {
      //V-wire
      react_point.SetX( fWire.x );
      react_point.SetZ( fWire.z );

    } else {
      //H-wire
      react_point.SetY( fWire.y );
      react_point.SetZ( fWire.z );
    }
  }    
   */ 
}
//_____________________________________________________________________________
TVector2 ReactVertex::Get_beamCenter() const
{
  if (!f_hasData) {
    Error("Get_beamCenter",
	  "ReactVertex not created with valid run-data, cannot compute beam-center.");
    return TVector2(0,0);
  }

  //NOTE: this is given in hall-coordinates (HCS), in meters, from the APEX/
  // scattering-chamber center. 
  
  return fBeamCenter; 
}
//_____________________________________________________________________________
ReactVertex::OpticsWire_t ReactVertex::Get_wire() const {

  if (!f_isWireMode) {
    Warning("Get_wire", "Wire requested, but not in wire-mode, so wire returned is null.");
  }
  return fWire;
}
//_____________________________________________________________________________


}
