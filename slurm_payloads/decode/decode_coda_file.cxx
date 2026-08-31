#include <APEX/decode.h> 

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
// stdlib headers
#include <fstream>
#include <stdexcept>
#include <string> 
#include <cstdio>
#include <memory>

namespace APEX
{
namespace decode
{

void decode_coda_file(const std::string& path_infile,
		      const std::string& path_outfile,
		      Long64_t first_event, 
		      Long64_t last_event, 
		      const std::string& path_odef)
{  
  Info(__func__, "starting decode:\n"
    "  input file:   %s\n"
    "  output file:  %s\n"
    "  odef:         %s\n"
    "  event range:  [%lli, %lli]\n",
    path_infile.c_str(),
    path_outfile.c_str(),
    path_odef.c_str(),
    first_event, last_event
  ); 

  std::string path_cuts   = "coinc.cuts";
  
  auto analyzer = std::make_unique<THaAnalyzer>(); 
  
  auto event    = std::make_unique<THaEvent>(); 
  
  auto run      = std::make_unique<THaRun>( path_infile.c_str() ); 
  
  //debug - limit number of events
  run->SetEventRange( first_event, last_event ); 
  
  //feed some important parameters to the 'analyzer' object
  analyzer->SetEvent( event.get() ); 
  analyzer->SetOutFile( path_outfile.c_str() );
  analyzer->SetOdefFile( path_odef.c_str() );
  analyzer->SetCutFile( path_cuts.c_str() ); 
  
  analyzer->SetSummaryFile( "summary_example.log" ); 
  
  analyzer->SetCompressionLevel( 0 ); 
  
  analyzer->SetOutFile( path_outfile.c_str() );

  analyzer->SetVerbosity( 0 ); 
  
  analyzer->Process( run.get() );

  std::printf(__func__, "done\n"); 
}

}
}
