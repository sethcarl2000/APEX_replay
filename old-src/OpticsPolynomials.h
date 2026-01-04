#ifndef Podd_OpticsPolynomials_h_
#define Podd_OpticsPolynomials_h_

//////////////////////////////////////////////////////////////////////////
//
// OpticsPolynomials
// Basically, these are polynomials that mostly act as vectors, but can also
// compute things when handed TXtg or TXfp classes. 
//
//////////////////////////////////////////////////////////////////////////

//#include "THaPhysicsModule.h"
#include <vector>
#include <array>
#include <algorithm>
#include "TObject.h"
#include "OpticsVectors.h"

//____________________________________________________________________________________
class TRPoly : public TObject {

public: 
  
  struct RPolyElem_t {
    double coeff;
    std::array<int,5> powers;
  }; 
  
  TRPoly() : fElems{} {};
  
  //vector-like operations
  void Clear() { fElems.clear(); }
  
  void Push_back(RPolyElem_t elem) { fElems.push_back(elem); } 
  
  unsigned int Size() const { return (unsigned int)fElems.size(); }

  
  double Eval(const TXtg &Xtg) const;
  
  double df_dxj    (const TXtg &Xtg, unsigned int j) const; 
  double df_dxj_dxk(const TXtg &Xtg, unsigned int j, unsigned int k) const; 

  
private:

  std::vector<RPolyElem_t> fElems; 
  
  ClassDef(TRPoly,0); 
}; 
//____________________________________________________________________________________
class TFPoly : public TObject {

public: 
  
  struct FPolyElem_t {
    double coeff;
    std::array<int,4> powers;
  }; 
  
  TFPoly() : fElems{} {};
  
  //vector-like operations
  void Clear() { fElems.clear(); }
  
  void Push_back(FPolyElem_t elem) { fElems.push_back(elem); } 
  
  unsigned int Size() const { return (unsigned int)fElems.size(); }
  
  double Eval(const TXfp &Xfp) const;
  
private:

  std::vector<FPolyElem_t> fElems; 
  
  ClassDef(TFPoly,0); 
}; 
//____________________________________________________________________________________

#endif
