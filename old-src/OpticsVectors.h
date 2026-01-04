#ifndef Podd_OpticsVectors_h_
#define Podd_OpticsVectors_h_

//////////////////////////////////////////////////////////////////////////
//
// OpticsVectors
// The trajectory-variables explicitly used in optics calculations.
// not much to see here. 
//
//////////////////////////////////////////////////////////////////////////

//#include "THaPhysicsModule.h"
#include <vector>
#include <array>
#include <algorithm>
#include "TObject.h"
#include "ROOT/RVec.hxx"

//____________________________________________________________________________________
class TXtg : public TObject {
  std::array<double,5> arr;
  //std::vector<double> arr; 
  
public:
  //default constructor
  TXtg() : arr{} {}; 
  
  TXtg(double _x, double _y, double _dxdz, double _dydz, double _dp) : arr{}
  { arr[0]=_x; arr[1]=_y; arr[2]=_dxdz; arr[3]=_dydz; arr[4]=_dp; }
        
  //conversion constructor wtih vector
  TXtg(const std::vector<double> &v);
        
  //return an RVecD 
  //std::array<double,5> Data()      const { return arr; }  
  ROOT::RVec<double>   Get_RVecD() const
  { return ROOT::RVec<double>( arr.data(), arr.data()+5 ); }  
  
  //const & not-const access to member elements
  double & x()       { return arr[0]; }//________x
  double   x() const { return arr[0]; }
  double & y()       { return arr[1]; }//________y
  double   y() const { return arr[1]; }
  double & dxdz()       { return arr[2]; }//__dxdz
  double   dxdz() const { return arr[2]; }
  double & dydz()       { return arr[3]; }//__dydz
  double   dydz() const { return arr[3]; }
  double & dp()       { return arr[4]; }//______dp
  double   dp() const { return arr[4]; }
    
  //these are here so we can access this vec analogously to an RVecD 
  double & at(unsigned int i);
  double  get(unsigned int i) const; 
  
  void operator+=(const ROOT::RVec<double> &rhs); 

  //arithmetic operators
  TXtg operator+(const TXtg &rhs) const;
  TXtg operator-(const TXtg &rhs) const;
  TXtg operator*(double scale)    const;
  
  
  ClassDef(TXtg,0)
};
//____________________________________________________________________________________  
class TXfp : public TObject {
  std::array<double,4> arr;
  //std::vector<double> arr;
  
public:
  //default constructor
  TXfp() : arr{} {}; 
  
  TXfp(double _x, double _y, double _dxdz, double _dydz) : arr{}
  { arr[0]=_x; arr[1]=_y; arr[2]=_dxdz; arr[3]=_dydz; }
  
  //conversion constructor wtih vector
  TXfp(const std::vector<double> &v);
        
  //return an RVecD 
  //std::array<double,4> Data()      const { return arr; }  
  ROOT::RVec<double>   Get_RVecD() const
  { return ROOT::RVec<double>( arr.data(), arr.data()+4 ); }  
  
  //const & not-const access to member elements
  double & x()       { return arr[0]; }//________x
  double   x() const { return arr[0]; }
  double & y()       { return arr[1]; }//________y
  double   y() const { return arr[1]; }
  double & dxdz()       { return arr[2]; }//__dxdz
  double   dxdz() const { return arr[2]; }
  double & dydz()       { return arr[3]; }//__dydz
  double   dydz() const { return arr[3]; }
    
  //these are here so we can access this vec analogously to an RVecD 
  double & at(unsigned int i);
  double  get(unsigned int i) const; 
  
  double Distance(const TXfp &rhs) const; 

  ClassDef(TXfp,0)
};
//____________________________________________________________________________________

#endif
