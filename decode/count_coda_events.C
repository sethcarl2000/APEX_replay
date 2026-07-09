#include "THaRun.h"

#include <memory>
#include <iostream>
#include <fstream> 
#include <string>

using namespace std;

void count_coda_events(std::string path_input, std::string path_output)
{
  bool status_ok=true;
  
  auto run = make_unique<THaRun>(path_input.c_str(), "CODA input file");
  run->SetDataRequired(0);
  
  auto st = run->Init();
  if( st != THaRunBase::READ_OK ) {
    cerr << "<count_coda_events>:  Error initializing" << endl;
    status_ok=false; 
  }

  st = run->Open();
  if( st != THaRunBase::READ_OK ) {
    cerr << "<count_coda_events>:  Error opening" << endl;
    status_ok=false; 
  }
  
  ULong64_t iev = 0;
  while( (st = run->ReadEvent()) == THaRunBase::READ_OK ) { ++iev; }
  
  cout << "<count_coda_events>:  Number of events read = " << iev << endl;
  
  cout << "<count_coda_events>:  Finished file scan, final status = ";
  switch( st ) {
  case THaRunBase::READ_EOF:
    cout << "EOF"; break;
  case THaRunBase::READ_ERROR:
    status_ok=false; 
    cout << "ERROR"; break;
  case THaRunBase::READ_FATAL:
    status_ok=false; 
    cout << "FATAL"; break;
  default:
    status_ok=false; 
    cout << "UNKNOWN? = " << st; break;
  }
  cout << endl;
  if( st != THaRunBase::READ_EOF ) {
    status_ok=false; 
    cerr << "<count_coda_events>:  Error reading ev " << iev << endl;
  }

  st = run->Close();
  if( st != THaRunBase::READ_OK ) {
    status_ok=false; 
    cerr << "<count_coda_events>:  Error closing?" << endl;
  }
  cout << "<count_coda_events>:  All successful" << endl;

  std::fstream outfile(path_output.c_str(), std::ios::out | std::ios::trunc);

  if (status_ok) {
    outfile << iev << "\n";
  } else {
    outfile << "bad-read\n";
  }
  outfile.close();

  return; 
}
