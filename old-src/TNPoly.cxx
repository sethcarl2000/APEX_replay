//*-- Author :    Seth Hall   09-Jul-24

//////////////////////////////////////////////////////////////////////////
//     
// TNPoly
// This is the backbone of the 'reverse optics' model.
// effectivley, all it is is a polynomial, with a appendable list of elements,
// defined on some R^n space.
// In the case of the reverse-optics, it will be defined on an R^4 space,
//  + using 3 polynomials in conjunction lets you make a f:R^4->R^3 map.
//  
//     
//////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iostream>
#include "TNPoly.h"
#include "ROOT/RVec.hxx"
#include "TString.h"

using namespace ROOT::VecOps;
using RVecD = RVec<double>; 

using namespace std;

//_____________________________________________________________________________
TNPoly::TNPoly(const int nDoF, const int order) :
  fnDoF(nDoF)
{  
  if (order>0) {
    fOrder = order;
    AutoConstructPoly(order,nDoF);
  }
}
//_____________________________________________________________________________
TNPoly::~TNPoly() {/*destructor*/};
//_____________________________________________________________________________
RVecD TNPoly::Eval_noCoeff(const RVecD &X) const
{
  //check to make sure that the input vector is the right size
  if ((int)X.size() != fnDoF) {
    Error("TNPoly::Eval_noCoeff",
	  "Size of input vector (%i) does not match poly nDoF (%i)",
	  (int)X.size(), fnDoF);
    return {};
  }

  RVecD vals; vals.reserve(Get_nElems());
  
  for (const TNPolyElem &elem : fElems) {

    double val=1.;
    for (int i=0; i<fnDoF; i++) {
      if (elem.powers.at(i)>0) {
	val *= pow( X.at(i), elem.powers.at(i) );
      }
    }
    //check to see if this is the one element for which all pows. are zero
    vals.push_back( val );
  }

  return vals;
}
//_____________________________________________________________________________
double TNPoly::Eval(const RVecD &coeff, const RVecD &X) const
{
  if (coeff.size() != Get_nElems()) { 
    Error("Eval(RVecD,RVecD)",
	  "Wrong number of coeffs. given: expected %i, recieved %i.",
	  Get_nElems(),
	  (int)coeff.size() );
    return -1e30;
  }
        
  if ((int)X.size() != fnDoF) { 
    Error("Eval(RVecD,RVecD)", "Size of input vector (%i) does not match poly nDoF (%i)",
	  (int)X.size(), (int)fnDoF);
    return -1e30;
  }
      
  //now, actually evaluate
  double val=0.; 

  //for each element, raise each val in X to the right power, the multiply
  // by the coressponding coeff.
  for (unsigned int i=0; i<Get_nElems(); i++) { const auto elem = fElems[i]; 
    
    double elem_val=1.;
    
    for (int d=0; d<fnDoF; d++) {
      if (elem.powers.at(d)>0) { 
	elem_val *= pow( X.at(d), elem.powers.at(d) );
      }
    }
    //check to see if this is the one element for which all pows. are zero
    val += elem_val * coeff.at(i); 
	
  }//for (i=0; i<Get_nElems(); i++) 
      
  return val;
}
//_____________________________________________________________________________
double TNPoly::Eval(const RVecD &X) const 
{
  //same as above, but use the coefficients 'hard-coded' to each element
  
  if ((int)X.size() != fnDoF) { 
    Error("Eval(RVecD)", "Size of input vector (%i) does not match poly nDoF (%i)",
	  (int)X.size(), (int)fnDoF);
    return -1e30;
  }
      
  //now, actually evaluate
  double val=0.; 

  //for each element, raise each val in X to the right power, the multiply
  // by the coressponding coeff.
  for (const auto &elem : fElems) {
    
    double elem_val=1.;
    
    for (int d=0; d<fnDoF; d++) {
      if (elem.powers.at(d)>0) { 
	elem_val *= pow( X.at(d), elem.powers.at(d) );
      }
    }
    //check to see if this is the one element for which all pows. are zero
    val += elem_val * elem.coeff; 
	
  }//for (i=0; i<Get_nElems(); i++) 
      
  return val;
}
//_____________________________________________________________________________
RVecD TNPoly::Gradient(const RVecD &coeff, const RVecD &X) const
{
  if (coeff.size() != Get_nElems()) { 
    Error("TNPoly::Eval()",
	  "Wrong number of coeffs. given; expected %i, recieved %i.",
	  Get_nElems(), (int)coeff.size() );
    return {};
  }
  
  if ((int)X.size() != fnDoF) { 
    Error("TNPoly::Eval_noCoeff",
	  "Size of input vector (%i) does not match poly nDoF (%i)",
	  (int)X.size(), fnDoF);
    return {};
  }

  RVecD grad(fnDoF,0.); 
  
  //for each element, raise each val in X to the right power, the multiply
  // by the coressponding coeff. 
  for (unsigned int i=0; i<Get_nElems(); i++) { auto elem = fElems[i];
    
    RVecD grad_elem(fnDoF,1.); 

    for (int g=0; g<fnDoF; g++) { 

      int pow_deriv = elem.powers[g];
      if (pow_deriv<1) {
        grad_elem[g] = 0.; continue; 
      } else { 
        grad_elem[g] *= ((double)pow_deriv)*pow(X[g], pow_deriv-1);
      } 

      //now, deal with the elems we are NOT taking the deriv of 
      for (int d=0; d<fnDoF; d++) { if (d==g) continue; 
        if (elem.powers.at(d)>0) { 
          grad_elem[g] *= pow( X.at(d), elem.powers.at(d) );
        }
      }
    }

    /* std::cout << "pows {"; 
    for (int p    : elem->powers) std::cout << p << " "; 
    std::cout << "} grad += {";
    for (double g : grad_elem) std::cout << TString::Format("% -*.4f  ",4,g); 
    std::cout << TString::Format("}  grad += %-*.4f",5,coeff->at(i)) << endl;  */
     

    grad += coeff.at(i) * grad_elem; 
  }

  /* std::cout << "grad = {";
  for (double g : grad) std::cout << TString::Format("% -*.4f  ",4,g); 
  std::cout << "}" << endl;  */
  
  //for (i=0; i<Get_nElems(); i++) 
  return grad;
}
//_____________________________________________________________________________
RVecD TNPoly::Get_elemPoly(unsigned int i) const
{
  if (i >= Get_nElems()) {
    Warning("TNPoly::Get_elemPoly",
	    "Invalid element-index requested; %i, max is %i",i,Get_nElems()-1);
    return {};
  }
  return fElems.at(i).poly;
}
//_____________________________________________________________________________
RVec<int> TNPoly::Get_elemPowers(unsigned int i) const
{
  if (i >= Get_nElems()) {
    Warning("TNPoly::Get_elemPowers",
	    "Invalid element-index requested; %i, max is %i",i,Get_nElems()-1);
    return {};
  }
  return fElems.at(i).powers;
}
//_____________________________________________________________________________
double TNPoly::Get_elemCoeff(unsigned int i) const
{
  if (i >= Get_nElems()) {
    Warning("Get_elemCoeff",
	    "Invalid element-index requested; %i, max is %i",i,Get_nElems()-1);
    return {};
  }
  return fElems.at(i).coeff;
}
//_____________________________________________________________________________
void TNPoly::Add_element(const RVec<int> &pows, const RVec<double> &poly,
			 double coefficient)
{
  if (pows.size()!=fnDoF) {
    Error("Add_element", "Num. of powers (%i) incorrect size for this Poly's DoF (%i)",
	  (int)pows.size(), (int)fnDoF);
    return;
  }
  fElems.push_back( TNPolyElem(pows,poly,coefficient) ); 
}
//_____________________________________________________________________________
void TNPoly::AutoConstructPoly(const int max_power, const int nDoF)
{
  RVec<int> powers(max_power,0); 
  
  while (1) {
    //create a new poly element with our current 'powers'
    //auto elem = TNPolyElem;
    RVec<int> pows; 
    
    for (int x=0; x<nDoF; x++) {
      int myPow=0;
      for (int xPow : powers) if (xPow==x) myPow++;
      
      pows.push_back(myPow);
    }
          
    /*cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; 
      for (int p=0; p<max_power; p++) { 
      for (int d=0; d<nDoF+1; d++) cout << (powers.at(p)==d?1:0) << " ";
      cout << endl; 
      }*/ 

    /*cout << "elem : "; 
      for (int d=0; d<nDoF; d++) cout << elem->powers.at(d) << ","; 
      cout << endl; */ 

    fElems.push_back( TNPolyElem(pows,{}) );
    
    //now, increase to the next set of powers
    for (int it1=powers.size()-1; it1>=0; it1--) { 
      //go back up and reset higher powers
            
      if (powers.at(it1)<nDoF) {
	powers.at(it1) += 1; 
	for (unsigned int it2=it1+1; it2<powers.size(); it2++) {
	  powers.at(it2) = powers.at(it1);
	}
	break; 
      }
      
      if (it1==0) {
	//cout << "polynomial size = " << elems.size() << endl; 
	return; 
      }
            
    }//for (int it1=powers.size()-1 ... [power-escalating group]        
  }//while (1)
}
//_____________________________________________________________________________

ClassImp(TNPoly)

