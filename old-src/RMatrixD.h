#ifndef Podd_RMatrixD_h_
#define Podd_RMatrixD_h_

//////////////////////////////////////////////////////////////////////////
//
// RMatrixD
//
// A jankly little lin-algebra class which is able to work (implicitly)
// with the ROOT::RVec<double> class, 'cause they're a lot more convenient
// to work with than TVectorD's are. 
//
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "TObject.h"
#include "ROOT/RVec.hxx"

using RVecD = ROOT::RVec<double>; 

//Basically just a wrapper for TMatrixD
//________________________________________________________________________________
class RMatrixD : public TObject {
public:
  
  RMatrixD(const unsigned int nr=1, const unsigned int nc=1,
	   const double init=0.); 

  RMatrixD(const unsigned int nr, const unsigned int nc,
	   const ROOT::RVec<double> &array);
  
  RMatrixD(const unsigned int nr, const unsigned int nc,
	   const ROOT::RVec<double> *array);
  
  //copy constructor
  //RMatrixD(const RMatrixD &mat); 
  
  ~RMatrixD() { fElems.clear(); } 

  ROOT::RVec<double> Solve(const ROOT::RVec<double> &B) const;

  //element-wise access
  double& at(unsigned int i,unsigned int j);
  double  at(unsigned int i,unsigned int j) const; 
    
  //multiplication by ROOT::RVec<double>
  ROOT::RVec<double> operator*(const ROOT::RVec<double> &rhs) const; 

  //adding two matrices
  RMatrixD           operator+(const RMatrixD &rhs) const; 

  //same as at(i,j), but without bounds-checking
  double get(unsigned int i, unsigned int j) const
  { return fElems[GetNCols()*i + j]; }
  
  unsigned int GetNCols() const { return fnCols; }
  unsigned int GetNRows() const { return fnRows; }

  void Print(); 
  
  bool &ReportSingular() {return f_reportSingular;}

  ROOT::RVec<double> Data() const { return fElems; }
  
private:
  
  unsigned int fnCols, fnRows; 
  bool f_isSquare; 
  
  // (default==true) 
  // controls wheter or not an error message is printed when a singular matrix is encountered
  bool f_reportSingular; 
  
  ROOT::RVec<double> fElems; 

  ClassDef(RMatrixD,0)
}; 
//________________________________________________________________________________
  

#endif
