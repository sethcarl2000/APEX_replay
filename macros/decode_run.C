
#include "../include/def_apex.h"

#include <fstream>

using namespace std;

#define DECODE_ONLY true
#define VERBOSITY 1

void decode_run( TString inFile_path
		 ="/cache/mss/halla/apex/raw/apex_4016.dat.0",
		 TString outFile_path
		 ="decode.4016.root",
		 TString path_outDef
		 ="outDefs/optics.odef",
		 ULong64_t event_first =0, 
		 ULong64_t event_last  =1e7 )
{
  const char* const here = "decode_run.C"; 
  
  Info(here, "inFile  path = \"%s\"",inFile_path.Data());
  Info(here, "outFile path = \"%s\"",outFile_path.Data());
  Info(here, "var defs     = \"%s\"",path_outDef.Data());   
  
  TString path_cuts   = "coinc.cuts";
  
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
  RHRS->AddDetector( new THaApexShower  ("ps",  "pre-shower EMC") ); 
  RHRS->AddDetector( new THaApexShower  ("sh",  "shower EMC")     ); 
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
  LHRS->AddDetector( new THaApexShower ("prl1",  "pre-shower EMC")  ); 
  LHRS->AddDetector( new THaApexShower ("prl2",  "shower EMC")      ); 
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

  
  
  //now, setup 
  auto analyzer = new THaAnalyzer;    
    
  auto event = new THaEvent; 
  
  auto run = new THaRun( inFile_path ); 
  
  //debug - limit number of events
  run->SetEventRange( event_first, event_last ); 
  
  //feed some important parameters to the 'analyzer' object
  analyzer->SetEvent( event ); 
  analyzer->SetOutFile( outFile_path );
  analyzer->SetOdefFile( path_outDef );
  analyzer->SetCutFile( path_cuts ); 
  
  analyzer->SetSummaryFile( "summary_example.log" ); 
  
  analyzer->SetCompressionLevel( 0 ); 
  
  Info(here,"Processing run"); 
  
  analyzer->SetOutFile( outFile_path );

  analyzer->SetVerbosity( VERBOSITY ); 
  
  analyzer->Process( run ); 
  
  Info(here,"Done processing run"); 
}
