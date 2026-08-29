#ifndef APEX_decode_Shower_h_
#define APEX_decode_Shower_h_

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// APEX::Shower                                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <THaPidDetector.h>
#include <vector>

class TClonesArray;

namespace APEX
{
namespace decode
{

class Shower : public THaPidDetector {
  
public:
  Shower( 
    const char* name, 
    const char* description = "",
    THaApparatus* a = nullptr 
  ) : THaPidDetector(name,description,a) {};

  // default ctor
  Shower() : THaPidDetector() {};

  ~Shower() override;

  void       Clear( Option_t* ="" ) override;
  Int_t      Decode( const THaEvData& ) override;
  Int_t      CoarseProcess( TClonesArray& tracks ) override;
  Int_t      FineProcess( TClonesArray& tracks ) override;
  
          Int_t      GetNclust() const { return fNclust; }
          Int_t      GetNhits() const  { return fNhits; }
          Float_t    GetE() const      { return fE; }
          Float_t    GetX() const      { return fX; }
          Float_t    GetY() const      { return fY; }

protected:

  // Mapping (see also fDetMap)
  std::vector< std::vector<UShort_t> > fChanMap{}; // Logical channel numbers
                                                 // for each detector map module

  // Configuration
  Int_t      fNclublk{0};   // Max. number of blocks composing a cluster
  Int_t      fNrows{0};     // Number of rows

  // Geometry
  Float_t*   fBlockX{nullptr};    // [fNelem] x positions (cm) of block centers
  Float_t*   fBlockY{nullptr};    // [fNelem] y positions (cm) of block centers

  // Calibration
  Float_t*   fPed{nullptr};       // [fNelem] Pedestals for each block
  Float_t*   fGain{nullptr};      // [fNelem] Gains for each block

  // Other parameters
  Float_t    fEmin{0.};      // Minimum energy for a cluster center

  // Per-event data
  Int_t      fNhits{0};     // Number of hits
  Float_t*   fA{nullptr};         // [fNelem] Array of ADC amplitudes of blocks
  Float_t*   fA_p{nullptr};       // [fNelem] Array of ADC minus pedestal values of blocks
  Float_t*   fA_c{nullptr};       // [fNelem] Array of corrected ADC amplitudes of blocks
  Float_t    fAsum_p{0.};    // Sum of blocks ADC minus pedestal values
  Float_t    fAsum_c{0.};    // Sum of blocks corrected ADC amplitudes
  Int_t      fNclust{0};    // Number of clusters
  Float_t    fE{0.};         // Energy (MeV) of main cluster
  Float_t    fX{0.};         // x position (cm) of main cluster
  Float_t    fY{0.};         // y position (cm) of main cluster
  Int_t      fMult{0};      // Number of blocks in main cluster
  Int_t*     fNblk{nullptr};      // [fNclublk] Numbers of blocks composing main cluster
  Float_t*   fEblk{nullptr};      // [fNclublk] Energies of blocks composing main cluster

  void           DeleteArrays();
  Int_t  ReadDatabase( const TDatime& date ) override;
  Int_t  DefineVariables( EMode mode = kDefine ) override;

  //ClassDefOverride(Shower,1); 
  
};

}
}

////////////////////////////////////////////////////////////////////////////////

#endif