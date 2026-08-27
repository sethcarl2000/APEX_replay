#include <APEX_decode.h> 

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

#define DECODE_ONLY true
#define VERBOSITY 1

namespace APEX
{
namespace decode
{

void process_coda_file(const std::string& path_infile,
			   const std::string& path_outfile,
			   ULong64_t first_event =0, 
			   ULong64_t last_event  =1e7, 
			   const std::string& path_odef ="outDefs/optics.odef" )
{  
  std::printf("in <%s>: starting decode:\n"
    "  input file:   %s\n"
    "  output file:  %s\n"
    "  odef:         %s\n"
    "  first event:  %ull\n"
    "  last event:   %ull\n",
    __func__,
    path_infile.c_str(),
    path_outfile.c_str(),
    path_odef.c_str(),
    first_event,
    last_event
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

  analyzer->SetVerbosity( VERBOSITY ); 
  
  analyzer->Process( run.get() );

  std::printf("in <%s>: done\n", __func__); 
}

}
}
