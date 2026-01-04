//*-- Author :    Seth Hall   09-Jul-24

//////////////////////////////////////////////////////////////////////////
//     
// OpticsVectors
// both of these classes are basically wrappers for std::arrays, with few
// extra little featuers.
//    My thinking in making these is that it will be helpful to have these
// classes defined explcitily, rahter than passing aroudn a bunch of
// std::array<double,5> and std::array<double,4> classes, if only for the
// sake of readability. 
//     
//////////////////////////////////////////////////////////////////////////


#include "OpticsVectors.h"
#include "ROOT/RVec.hxx"
#include <vector>
#include <cmath> 
#include <algorithm>
#include <cmath> 
#include <array>

using namespace std; 

//_____________________________________________________________________________
TXtg::TXtg(const std::vector<double> &v) : arr{}
{
  if (v.size() != 5) {
    Error("TXtg (vector constructor)", "input vector wrong size (%i), should be 5.",
	  (int)v.size());
    return;
  }
  std::copy( v.begin(), v.begin()+5, arr.begin() ); 
}
//_____________________________________________________________________________
void TXtg::operator+=(const ROOT::RVec<double> &rhs)
{
  //this is explicitly for use with the TOpticsRays::Find_Xtg method.
  if (rhs.size() != 5) {
    Error("operator+=(RVecD)", "input vector wrong size (%i), should be 4.",
	  (int)rhs.size());
    return;
  }
  for (unsigned int i=0; i<5; i++) arr[i] += rhs.at(i); 
}
//_____________________________________________________________________________
TXtg TXtg::operator+(const TXtg &rhs) const
{
  return TXtg( arr[0] + rhs.get(0),
	       arr[1] + rhs.get(1),
	       arr[2] + rhs.get(2),
	       arr[3] + rhs.get(3),
	       arr[4] + rhs.get(4) ); 
}
//_____________________________________________________________________________
TXtg TXtg::operator-(const TXtg &rhs) const
{
  return TXtg( arr[0] - rhs.get(0),
	       arr[1] - rhs.get(1),
	       arr[2] - rhs.get(2),
	       arr[3] - rhs.get(3),
	       arr[4] - rhs.get(4) ); 
}
//_____________________________________________________________________________
TXtg TXtg::operator*(double scale) const
{
  return TXtg( arr[0]*scale,
	       arr[1]*scale,
	       arr[2]*scale,
	       arr[3]*scale,
	       arr[4]*scale ); 
}
//_____________________________________________________________________________
TXfp::TXfp(const std::vector<double> &v) : arr{}
{
  if (v.size() != 4) {
    Error("TXfp (vector constructor)", "input vector wrong size (%i), should be 4.",
	  (int)v.size());
    return;
  }
  std::copy( v.begin(), v.begin()+4, arr.begin() ); 
}
//_____________________________________________________________________________
double& TXtg::at(unsigned int i)
{
  if (i>=5) {
    Error("at", "Index given (i=%i) out-of-range; should be i=[0,4].", (int)i);
    return *(new double(std::nan("1")));
  }
  return arr[i]; 
}
//_____________________________________________________________________________
double& TXfp::at(unsigned int i)
{
  if (i>=4) {
    Error("at", "Index given (i=%i) out-of-range; should be i=[0,3].", (int)i);
    return *(new double(std::nan("1")));
  }
  return arr[i]; 
}
//_____________________________________________________________________________
double TXtg::get(unsigned int i) const 
{
  if (i>=5) {
    Error("at", "Index given (i=%i) out-of-range; should be i=[0,3].", (int)i);
    return *(new double(std::nan("1")));
  }
  return arr[i]; 
}
//_____________________________________________________________________________
double TXfp::get(unsigned int i) const 
{
  if (i>=4) {
    Error("at", "Index given (i=%i) out-of-range; should be i=[0,3].", (int)i);
    return *(new double(std::nan("1")));
  }
  return arr[i]; 
}
//_____________________________________________________________________________
double TXfp::Distance(const TXfp &rhs) const
{
  double ret=0;//arr[0];
  for (uint i=0; i<4; i++) ret += pow( get(i) - rhs.get(i), 2); // - rhs.at(i), 2 );
  return sqrt(ret); 
}
ClassImp(TXfp); 
ClassImp(TXtg);

