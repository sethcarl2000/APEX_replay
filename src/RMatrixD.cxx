//*-- Author :    Seth Hall;  22 Jul 2024

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// RMatrixD
//
// This is a very bare-bones Lin. algebra class, meant to deal primarily
// with square-matricies. Specifically, it's designed to work with
// the ROOT::RVecD framework. 
//
//////////////////////////////////////////////////////////////////////////

#include "TString.h"
#include "ROOT/RVec.hxx"
#include "RMatrixD.h"

using namespace std;
using namespace ROOT::VecOps; 

using RVecD = RVec<double>; 

//_______________________________________________________________________________
//_______________________________________________________________________________
RMatrixD::RMatrixD(const unsigned int nr,const unsigned int nc, const double init)
  : fnCols(nc),
    fnRows(nr),
    f_isSquare(nr==nc), 
    f_reportSingular(true)
{
  fElems = RVecD(nr*nc,init); 
}
//_______________________________________________________________________________
RMatrixD::RMatrixD(const unsigned int nr,const unsigned int nc, const RVecD &array)
  : fnCols(nc),
    fnRows(nr),
    f_isSquare(nr==nc)
{
  if (nr*nc != array.size()) {
    Error("RMatrixD(const RVecD *array)",
	  "Array size does not match given matrix dims! Initialzed as all zeros.");
    fElems = RVecD(nr*nc,0.); 
  } else { 
    fElems = RVecD(array);
  }
}
//_______________________________________________________________________________
RMatrixD::RMatrixD(const unsigned int nr,const unsigned int nc, const RVecD *array)
  : fnCols(nc),
    fnRows(nr),
    f_isSquare(nr==nc)
{
  if (nr*nc != array->size()) {
    Error("RMatrixD(const RVecD *array)",
	  "Array size does not match given matrix dims! Initialzed as all zeros.");
    fElems = RVecD(nr*nc,0.); 
  } else { 
    fElems = RVecD(*array);
  }
}
//_______________________________________________________________________________
/*RMatrixD::RMatrixD(const RMatrixD &mat)
  : TObject(),
    fnCols(mat.GetNCols()),
    fnRows(mat.GetNRows()),
    f_isSquare(fnCols==fnRows),
    f_reportSingular(true),
    fElems(mat.Data())
{
  //copy-constructor
  } */
//_______________________________________________________________________________
double& RMatrixD::at(unsigned int i, unsigned int j)
{
  if (i>=GetNRows() || j>=GetNCols()) {
    Error("at","Invalid element access attempt: (%i,%i), last elem is (%i,%i)",
	  i,j, GetNRows()-1,GetNCols()-1);
    return *(new double(std::nan("1")));
  }
  return fElems[fnCols*i + j]; 
}
//_______________________________________________________________________________
double RMatrixD::at(unsigned int i, unsigned int j) const
{
  if (i>=GetNRows() || j>=GetNCols()) {
    Error("at","Invalid element access attempt: (%i,%i), last elem is (%i,%i)",
	  i,j, GetNRows()-1,GetNCols()-1);
    return *(new double(std::nan("1")));
  }
  return fElems[GetNCols()*i + j]; 
}
//_______________________________________________________________________________
RVecD RMatrixD::operator*(const RVecD &rhs) const
{
  //check vector size
  if (rhs.size() != GetNCols()) {
    Error("operator*","Input vector does not match matrix col. size.");
    return {};
  }
  
  RVecD out(fnRows,0.);
  
  for (unsigned int i=0; i<GetNRows(); i++) 
    for (unsigned int j=0; j<GetNCols(); j++) out[i] += get(i,j) * rhs[j];
  
  return out; 
} 
//_______________________________________________________________________________
RMatrixD RMatrixD::operator+(const RMatrixD &rhs) const
{
  //check matrix sizes
  if (rhs.GetNCols() != GetNCols() ||
      rhs.GetNRows() != GetNRows()) {
    Error("operator+ (Matrix addition)","Matrices are not of matching dimension (%ix%i vs %ix%i).",
	  (int)GetNRows(),    (int)GetNCols(),
	  (int)rhs.GetNRows(),(int)rhs.GetNCols()
	  );
    return RMatrixD();
  }

  //make a copy of this matrix
  RMatrixD sum = *this;
  
  for (unsigned int i=0; i<GetNRows(); i++) 
    for (unsigned int j=0; j<GetNCols(); j++) sum.at(i,j) += rhs.at(i,j);
  
  return sum; 
} 
//_______________________________________________________________________________
RVecD RMatrixD::Solve(const RVecD &B) const
{
  //check if we're trying to call this when this matrix is not square
  if (!f_isSquare) {
    Error("Solve", "Trying to invert non-square matrix!");
    return {};
  }

  const unsigned int N = GetNRows(); 
  
  //takes the NxN matrix A, and the N-vector B, and solves it. 
  //performs numerical LU factorization 
  
  //ASSUMES MATRIX IS NONSINGULAR!!!!
  
  //check to make sure all our inputs are of a consistent size
  if (N != (unsigned int)B.size()) {
      Error("Solve",
	    "Matrix R^(%ux%u) and vector R^(%u) do not match!",
	    N,N,(unsigned int)B.size()); 
    return {}; 
  }
  
  //the U-matrix starts as a copy of the 'A' input-matrix
  RMatrixD U = RMatrixD(*this);  
  //we initialize the L-matrix with all zeros
  RMatrixD L(N,N, 0.);  
  
  //Create the U (upper triangular) and L (lower triangular) matrices
  for (unsigned int ii=0; ii<N; ii++) { 
    //loop over a (successivley smaller) sub-matrix
    double a_00 = U.at(ii,ii); 
    for (unsigned int i=ii+1; i<N; i++) { 
      double a_i0 = U.at(i,ii); 
      L.at(i,ii) = a_i0/a_00; 
      for (unsigned int j=ii; j<N; j++) { 
        U.at(i,j) += ( -a_i0/a_00 ) * U.at(ii,j); 
      }
    }
  }
  double det = 1.; 
  for (unsigned int ii=0; ii<N; ii++) { L.at(ii,ii)=1.; det *= U.at(ii,ii); }

  //det is NaN 
  if ( det != det ) {
    if (f_reportSingular) Error("Solve", "Determinant is NaN."); 
    return {}; 
  }
  
  RVecD y(N,0.);

  //solve the system Ly = b
  for (unsigned int i=0; i<N; i++) { 
    y[i] = B[i]; 
    for (unsigned int j=0; j<i; j++) y[i] += -L.at(i,j)*y[j];
  }

  RVecD X(N,0.); 
  
  //now solve Ux = y
  for (int i=N-1; i>=0; i--) { 
    X[i] = y[i]; 
    for (int j=N-1; j>i; j--) { 
      X[i]  +=  - U.at(i,j) * X[j]; 
    }
    X[i] *= 1./U.at(i,i); 
  }
  return X; 
}
//_______________________________________________________________________________
//_______________________________________________________________________________
//_______________________________________________________________________________
void RMatrixD::Print()
{
  //print the elems of this matrix (for debugging purposes)
  cout << TString::Format("Matrix: %ix%i",GetNRows(),GetNCols()) << endl; 

  for (unsigned int i=0; i<fnRows; i++) {
    for (unsigned int j=0; j<fnCols; j++) 
      cout << TString::Format(" %+0.3e", this->at(i,j) );
    cout << endl;
  } 
}
//_______________________________________________________________________________

ClassImp(RMatrixD)
