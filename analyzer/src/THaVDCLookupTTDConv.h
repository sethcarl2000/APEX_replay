#ifndef PODD_VDC_LookupTTDConv
#define PODD_VDC_LookupTTDConv

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// THaVDCLookupTTDConv                                                     //
//                                                                           //
// Uses a drift velocity (um/ns) to convert time (ns) into distance (cm)     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaVDCTimeToDistConv.h"
#include "TMath.h"
#include <iostream>

namespace VDC {

  class LookupTTDConv : public TimeToDistConv {

  public:
    //    LookupTTDConv();
    LookupTTDConv();


    virtual ~LookupTTDConv() {}

    virtual Double_t ConvertTimeToDist( Double_t time, Double_t tanTheta,
				        Double_t* ddist=0 ) const;

    //    virtual Double_t ConvertDistToTime( Double_t dist, Double_t tanTheta) const = 0;
    
    Double_t GetLookupVal( Double_t time) const;

    Double_t GetLookupTime( Double_t dist) const;
    
    //    virtual Int_t    SetParameters( std::vector<Double_t> Table, Int_t NBins, Double_t R, Double_t Theta);

    virtual Int_t SetLookupParams( std::vector<Double_t> Table, Int_t NBins, Double_t Low, Double_t R, Double_t Theta, Double_t m1, Double_t m2);

    virtual Int_t PrintParameters() const;

protected:

    // Lookup table
    std::vector<Double_t> LTable;
    Int_t LNumBins; // number of bins in lookup table
    Double_t LowTime; // Smallest time for plane
    
    // angular parameter

    Double_t RCorr; //angular correction parameter: distance where correction shifts form from prop to time to constant
    Double_t Theta0; // angular correction central angle 


    // extension corrections gives values of distances for extreme times (less than zero or larger time than in table)

    Double_t M1 = 0.0; // gradient for times less than zero
    Double_t M2 = 0.0; // gradient for time greater than max
    
    const Double_t bin_res = 0.5e-9; // TDC resolution of 0.5 ns
    
    //    Double_t fdtime;      // uncertainty in the measured time

    ClassDef(LookupTTDConv,0)   // VDC Lookup TTD Conv class
  };
}

////////////////////////////////////////////////////////////////////////////////

#endif
