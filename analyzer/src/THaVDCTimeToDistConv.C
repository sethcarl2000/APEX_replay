///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// THaVDCTimeToDistConv                                                      //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaVDCTimeToDistConv.h"

using namespace std;

ClassImp(VDC::TimeToDistConv)

namespace VDC {

//_____________________________________________________________________________
TimeToDistConv::TimeToDistConv( UInt_t npar )
  : fNparam(npar), fDriftVel(kBig), fIsSet(false)
{
  // Constructor
}

//_____________________________________________________________________________
void TimeToDistConv::SetDriftVel( Double_t v )
{
  fDriftVel = v;
  if( fNparam == 0 )
    fIsSet = true;
}

//_____________________________________________________________________________
Int_t TimeToDistConv::SetParameters( const vector<double>& )
{
  if( fNparam == 0 )
    fIsSet = true;
  return 0;
}


  //_____________________________________________________________________________
Int_t TimeToDistConv::SetAngleParameters( Double_t R, Double_t Theta )
{
  if( fNparam == 0 )
    fIsSet = true;
  return 0;
}


// Int_t TimeToDistConv::PrintParameters(){
  
  
//   }
  

Int_t TimeToDistConv::SetLookupParams( std::vector<Double_t> Table, Int_t NBins, Double_t Low, Double_t R, Double_t Theta, Double_t M1, Double_t M2){

  fIsSet = true;
  return 0;
  
}
  


  
} // namespace VDC

////////////////////////////////////////////////////////////////////////////////
