// APEX
#include <APEX/decode.h> 
//#include <APEX/decode/Shower.h>
// Podd Headers
#include <THaRun.h>
#include <THaAnalyzer.h>
#include <THaEvent.h>
#include <THaHRS.h>
#include <THaVDC.h>
#include <THaScintillator.h>
#include <THaCherenkov.h>
#include <THaShower.h>
#include <THaIdealBeam.h>
#include <THaRasteredBeam.h>
#include <THaRaster.h>
#include <THaBPM.h>
#include <THaGlobals.h>
// ROOT headers
#include <TList.h> 
// stdlib headers 
#include <string> 
#include <cstdio>

namespace APEX
{
namespace decode
{

//initialize apps needed to decode raw APEX data
void init_decoders()
{
  //add the left HRS & associated detectors
  THaHRS *RHRS = new THaHRS("R", "Right HRS");
  RHRS->AutoStandardDetectors( kFALSE ); //add only our custom detectors

  auto R_vdc = new THaVDC("vdc", "VDC - Right"); 
  R_vdc->SetBit(THaVDC::kDecodeOnly);
  RHRS->AddDetector( R_vdc ); 
  
  //add other HRS detectors
  
  // The 'TapexShower' class (src: /test_lib) is a user-made
  // class which inherits from THaShower, and exists only to overload the 'Decode()'
  // method, for debug purposes.

  //decode showers only
  RHRS->AddDetector( new THaShower("ps",  "pre-shower EMC") ); 
  RHRS->AddDetector( new THaShower("sh",  "shower EMC")     ); 
  RHRS->AddDetector( new THaScintillator("s2",  "S2 Scintillator"));  
  //*/
  //decode cerenkov only
  RHRS->AddDetector( new THaCherenkov ("cer", "Gas Cherenkov counter"));
  
  
  gHaApps->Add( RHRS );   

  
  //add the left HRS & associated detectors
  THaHRS *LHRS = new THaHRS("L", "Right HRS");
  LHRS->AutoStandardDetectors( kFALSE ); //add only our custom detectors
  
  auto L_vdc = new THaVDC("vdc", "VDC - Left"); 
  L_vdc->SetBit(THaVDC::kDecodeOnly);
  LHRS->AddDetector( L_vdc ); 

  //decode shower only
  LHRS->AddDetector( new THaShower ("prl1",  "pre-shower EMC")  ); 
  LHRS->AddDetector( new THaShower ("prl2",  "shower EMC")      ); 
  LHRS->AddDetector( new THaScintillator("s2",  "S2 Scintillator"));  
  //*/
  //decode cerenkov only
  LHRS->AddDetector( new THaCherenkov ("cer", "Gas Cherenkov counter"));
    
  gHaApps->Add( LHRS ); 
  

  //add beam
  gHaApps->Add( new THaIdealBeam("ib", "Ideal beam") ); 
  
  //Left Rastered beam
  auto Lrb = new THaRasteredBeam("Lrb", "Rastered beam - Left"); 

  Lrb->AddDetector( new THaRaster("Raster2", "Downstream raster") );
  Lrb->AddDetector( new THaBPM   ("BPMA",    "Beam position monitor A") );
  Lrb->AddDetector( new THaBPM   ("BPMB",    "Beam position monitor B") );
  
  gHaApps->Add( Lrb ); 

  
  //Right Rastered beam
  auto Rrb = new THaRasteredBeam("Rrb", "Rastered beam - Right"); 
  
  Rrb->AddDetector( new THaRaster("Raster2", "Downstream raster") );
  Rrb->AddDetector( new THaBPM   ("BPMA",    "Beam position monitor A") );
  Rrb->AddDetector( new THaBPM   ("BPMB",    "Beam position monitor B") );
  
  gHaApps->Add( Rrb );
}

}
}
