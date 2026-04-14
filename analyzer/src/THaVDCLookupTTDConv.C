///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// THaVDCLookupTTDConv                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaVDCLookupTTDConv.h"
#include "TError.h"
#include <vector>

ClassImp(VDC::LookupTTDConv)

using namespace std;

namespace VDC {

//_____________________________________________________________________________
  LookupTTDConv::LookupTTDConv() : TimeToDistConv(2)

{
  // constuctor

  LNumBins = 0;

  LowTime = kBig;
  RCorr = kBig;
  Theta0 = kBig;
  
  
  
}

//_____________________________________________________________________________
Double_t LookupTTDConv::ConvertTimeToDist( Double_t time, Double_t tanTheta,
					     Double_t* ddist) const
{
  // Drift Velocity in m/s
  // time in s
  // Return m

  if( !fIsSet ) {
    Error( "VDC::LookupTTDConv::ConvertTimeToDist", "Parameters not set. "
	   "Fix database." );
    return kBig;
  }


  
  double dist = GetLookupVal(time); // get lookup table value for central angle
  
   
  const double R = RCorr; // radius of electric field of wire
  //  const double* fslope0 = &(par[1]); // central angle
  const double fslope0 = Theta0; // central 
  
  const double fTheta0 = TMath::ATan(fslope0);
  //  const double fTheta0 = 1/(*fslope0);
  /* const double fTheta0 = (*fslope0); */


  // slope in analyzer is 1/m where m is gradient of slope
  tanTheta = TMath::ATan(1.0/tanTheta);

  Double_t SecTheta = (1/TMath::Cos(tanTheta));
  Double_t CosTheta = TMath::Cos(tanTheta);

  Double_t SecTheta0 = (1/TMath::Cos(fTheta0));
  

  if (dist >= R*SecTheta0 ) {
    dist += R*(SecTheta -  SecTheta0);

  } else if (dist < R*SecTheta0 ) { 
    
    dist *= (SecTheta /  SecTheta0 );
  }
		 
  return dist;
}

 


Double_t LookupTTDConv::GetLookupVal( Double_t time) const{

  time -= LowTime;
  
  Int_t bin_no = time/(bin_res);

  Double_t dist = 0.0;

  if(time < 0.0){
    // extension for times smaller than zero
    dist = LTable[0] + M2*time;
    }
  else if(bin_no > LNumBins-2){
    // extension for times greater than max
    dist = LTable[LNumBins-1] + M1*(time - ((LNumBins-1)*bin_res)); 
  }
  else{
    dist = LTable[bin_no] + std::fmod(time,bin_res)*(LTable[bin_no+1]-LTable[bin_no]);
  }

  return dist;  
}


  
  //_____________________________________________________________________________
  Int_t LookupTTDConv::SetLookupParams( std::vector<Double_t> Table, Int_t NBins, Double_t Low, Double_t R, Double_t Theta, Double_t m1, Double_t m2)
{

  LTable = Table;     // Lookup table
  LNumBins = NBins;   // number of entries in table
  LowTime = Low;      // Smallest time for plane
  RCorr = R;          // Distance at which angular correction changes (nature of E field)
  Theta0 = Theta;     // central angle used in correction
  M1 = m1;
  M2 = m2;

  
  fIsSet = true;
  return 0;
}


  
  Int_t LookupTTDConv::PrintParameters () const
{

  std::cout << "LNumBins = " << LNumBins << std::endl;
  std::cout << "LowTime = " << LowTime << std::endl;
  std::cout << "RCorr = " << RCorr << std::endl;
  std::cout << "Theta0 = " << Theta0 << std::endl;

  std::cout << "Printing Lookup table: " << std::endl;

  Int_t i = 1;
  for( auto val : LTable){
    std::cout << val << " ";
    if(i%10==0){
      std::cout << std::endl;
    }    
    i++;
  }
  
  return 0;


}


  
//_____________________________________________________________________________
// void LookupTTDConv::SetDefaultParam()
// {
//   // Set some reasonable defaults for the polynomial coefficients and
//   // drift time uncertainty. Applicable to Hall A VDCs.

//   fA1tdcCor[0] = 2.12e-3;
//   fA1tdcCor[1] = 0.0;
//   fA1tdcCor[2] = 0.0;
//   fA1tdcCor[3] = 0.0;
//   fA2tdcCor[0] = -4.20e-4;
//   fA2tdcCor[1] =  1.3e-3;
//   fA2tdcCor[2] = 1.06e-4;
//   fA2tdcCor[3] = 0.0;

//   fdtime    = 4.e-9; // 4ns -> 200 microns
// }

} //namespace VDC

///////////////////////////////////////////////////////////////////////////////
