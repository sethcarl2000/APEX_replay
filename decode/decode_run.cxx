
#include "../include/argparse.hpp"
#include "../include/def_apex.h"

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
// src headers
#include <THaApexShower.h>

#include <fstream>
#include <stdexcept>
#include <string> 


#define DECODE_ONLY true
#define VERBOSITY 1

void decode_run( 
  std::string inFile_path,
  std::string outFile_path,
  std::string path_outDe,
  unsigned long event_first, 
  unsigned long event_last 
); 

int main(int argc, char* argv[])
{
  argparse::ArgumentParser program("decode_run", "1.0"); 

  program.add_description("This program decodes a CODA file and writes the output to a ROOT file. "
      "The user can specify the input and output file paths, as well as the range of events to process. "
      "The program uses the THaRun class to handle the CODA file and the THaAnalyzer class to perform the decoding."
  );

  program.add_argument("--input-path")
      .help("Input CODA file path")
      .required();

  program.add_argument("--output-path")
      .help("Output file path")
      .required();

  program.add_argument("--outDef-path")
      .help("rel. path to outDef file")
      .default_value("outDefs/optics.odef");
  
  program.add_argument("--event-first")
      .help("first event to process")
      .scan<'i', unsigned long>()
      .default_value(0UL);

  program.add_argument("--event-last")
      .help("last event to process")
      .scan<'i', unsigned long>()
      .default_value(1000000UL);

  try {
      program.parse_args(argc, argv);
  }
  catch (const std::runtime_error& err) {
      std::cerr << err.what() << std::endl;       
      std::cerr << program << std::endl;
      return -1; 
  }

  decode_run(
      program.get<std::string>("--input-path"),
      program.get<std::string>("--output-path"),
      program.get<std::string>("--outDef-path"),
      program.get<unsigned long>("--event-first"),
      program.get<unsigned long>("--event-last")
  );

  return 0;
}

void decode_run( std::string inFile_path,
		 std::string outFile_path,
		 std::string path_outDef,
		 unsigned long event_first, 
		 unsigned long event_last )
{
  const char* const here = "decode_run"; 
  
  Info(here, "inFile  path = \"%s\"",inFile_path.c_str());
  Info(here, "outFile path = \"%s\"",outFile_path.c_str());
  Info(here, "var defs     = \"%s\"",path_outDef.c_str());   
  
  std::string path_cuts   = "coinc.cuts";
  
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
  
  auto run = new THaRun( inFile_path.c_str() ); 
  
  //debug - limit number of events
  run->SetEventRange( event_first, event_last ); 
  
  //feed some important parameters to the 'analyzer' object
  analyzer->SetEvent( event ); 
  analyzer->SetOutFile( outFile_path.c_str() );
  analyzer->SetOdefFile( path_outDef.c_str() );
  analyzer->SetCutFile( path_cuts.c_str() ); 
  
  analyzer->SetSummaryFile( "summary_example.log" ); 
  
  analyzer->SetCompressionLevel( 0 ); 
  
  Info(here,"Processing run"); 
  
  analyzer->SetOutFile( outFile_path.c_str() );

  analyzer->SetVerbosity( VERBOSITY ); 
  
  analyzer->Process( run ); 
  
  Info(here,"Done processing run"); 

  return; 
}
