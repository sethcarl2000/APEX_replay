///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// THaVDCAnalyticTTDConv                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaVDCAnalyticTTDConv.h"
#include "TError.h"

ClassImp(VDC::AnalyticTTDConv)

using namespace std;

namespace VDC {

//_____________________________________________________________________________
AnalyticTTDConv::AnalyticTTDConv() : TimeToDistConv(9), fdtime(0)

{
  // Constructor

  for( Int_t i=0; i<4; ++i ) {
    fA1tdcCor[i] = fA2tdcCor[i] = kBig;
  }
}

//_____________________________________________________________________________
Double_t AnalyticTTDConv::CalcDist( Double_t time, Double_t tanTheta,
					     Double_t* ddist) const
{
  // Drift Velocity in m/s
  // time in s
  // Return m

  if( !fIsSet ) {
    Error( "VDC::AnalyticTTDConv::ConvertTimeToDist", "Parameters not set. "
	   "Fix database." );
    return kBig;
  }

//    printf("Converting Drift Time to Drift Distance!\n");

  // Find the values of a1 and a2 by evaluating the proper polynomials
  // a = A_3 * x^3 + A_2 * x^2 + A_1 * x + A_0
  Double_t a1 = 0.0, a2 = 0.0;

  tanTheta = 1.0 / tanTheta;

  for (Int_t i = 3; i >= 1; i--) {
    a1 = tanTheta * (a1 + fA1tdcCor[i]);
    a2 = tanTheta * (a2 + fA2tdcCor[i]);
  }
  a1 += fA1tdcCor[0];
  a2 += fA2tdcCor[0];

  Double_t dist = fDriftVel * time;
  Double_t unc  = fDriftVel * fdtime;  // watch uncertainty in the timing
  if (dist < 0) {
    // something screwy is going on
  } else if (dist < a1 ) {
    //    dist = fDriftVel * time * (1 + 1 / (a1/a2 + 1));
    dist *= ( 1 + a2 / a1);
    unc *=  ( 1 + a2 / a1);
  }  else {
    dist +=  a2;
  }

  if (ddist) *ddist = unc;

  return dist;
}

//_____________________________________________________________________________
Double_t AnalyticTTDConv::ConvertTimeToDist( Double_t time, Double_t tanTheta,
					     Double_t* ddist) const
{
  // Drift Velocity in m/s
  // time in s
  // Return m

  if( !fIsSet ) {
    Error( "VDC::AnalyticTTDConv::ConvertTimeToDist", "Parameters not set. "
	   "Fix database." );
    return kBig;
  }

  //    printf("Converting Drift Time to Drift Distance!\n");

  double dist = CalcDist(time, tanTheta, ddist); // get lookup table value for central angle
  
   
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

  
//_____________________________________________________________________________
Double_t AnalyticTTDConv::GetParameter( UInt_t i ) const
{
  // Get i-th parameter

  switch(i) {
  case 0:
  case 1:
  case 2:
  case 3:
    return fA1tdcCor[i];
  case 4:
  case 5:
  case 6:
  case 7:
    return fA2tdcCor[i-4];
  case 8:
    return fdtime;
  }
  return kBig;
}

//_____________________________________________________________________________
Int_t AnalyticTTDConv::SetAngleParameters( Double_t R, Double_t Theta)
{
  // Set coeffecients for angular correction

  
  RCorr = R;          // Distance at which angular correction changes (nature of E field)
  Theta0 = Theta;     // central angle used in correction
    
  return 0;
}


  
//_____________________________________________________________________________
Int_t AnalyticTTDConv::SetParameters( const vector<double>& parameters )
{
  // Set coefficients of a1 and a2 4-th order polynomial and uncertainty
  // of drift time measurement
  // 0-3: A1
  // 4-7: A2
  // 8: sigma_time (s)

  if( (UInt_t)parameters.size() < fNparam )
    return -1;

  for( size_t i = 0; i<4; ++i ) {
    fA1tdcCor[i] = parameters[i];
    fA2tdcCor[i] = parameters[i+4];
  }
  fdtime = parameters[8];

  fIsSet = true;
  return 0;
}

//_____________________________________________________________________________
// void AnalyticTTDConv::SetDefaultParam()
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
