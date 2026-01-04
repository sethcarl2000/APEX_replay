//*-- Author :    Seth Hall   15-Dec-24

//////////////////////////////////////////////////////////////////////////
//     
// TOpticsRays
// responsible for computing all possible target trajectories associated with
// a single X-fp vector (i.e. 'rays'). 
//
//////////////////////////////////////////////////////////////////////////

// \/\/ comment from senior dev Muon
//
//  aqQ1`B               V;PPPPPPPPPPPPPPPPPPPPL;................./
//

#include "TOpticsRays.h"
#include "OpticsPolynomials.h"
#include "OpticsVectors.h"
#include "RMatrixD.h" 
#include "ApexUtils.h"
#include <vector>
#include <cmath> 
#include <array>

using namespace std; 
using uint = unsigned int; 

//_________________________________________________________________________________________
TXfp TOpticsRays::Compute_Xfp(const TXtg &Xtg) const
{
  //assuming that poly's have been initialized
  return TXfp( fPolys_reverse[0]->Eval(Xtg),
	       fPolys_reverse[1]->Eval(Xtg),
	       fPolys_reverse[2]->Eval(Xtg),
	       fPolys_reverse[3]->Eval(Xtg) ); 
}
//_________________________________________________________________________________________
double TOpticsRays::Find_Xtg(const TXfp &Xfp,
			     TXtg &Xtg,
			     uint max_iterations,
			     double min_error) const
{
  //converge to a value of Xtg which most closely matches the given Xfp.
  double rms = std::nan("1"); 
  
  for (uint i=0; i<max_iterations; i++) {
    
    
    //the gradients of each coordinate (Xfp), with respect to each coordinate (of Xtg)
    double dfi_dxj[4][5];
    
    for (uint i=0; i<4; i++)
      for (uint j=0; j<5; j++)
	dfi_dxj[i][j] = fPolys_reverse[i]->df_dxj(Xtg, j);
    
    
    double fi_minus_zi[4];
    //compute difference between model (fi) and actual (Xfp)
    TXfp fi = Compute_Xfp(Xtg);
    
    for (uint i=0; i<4; i++)
      fi_minus_zi[i] = fi.at(i) - Xfp.get(i); 
    
    
    RVecD Fx(5);
    
    for (uint j=0; j<5; j++)
      for (uint i=0; i<4; i++)
	Fx[j] += fi_minus_zi[i]*dfi_dxj[i][j];
    
    
    //now, compute the Jacobian (starts out as a 5x5 matrix with all 0-s)
    RMatrixD Jx(5,5, 0.); Jx.ReportSingular()=false; 
    
    for (uint j=0; j<5; j++)
      for (uint k=0; k<5; k++)
	for (uint i=0; i<4; i++)
	  Jx.at(j,k)
	    += dfi_dxj[i][j]*dfi_dxj[i][k]
	    +  fi_minus_zi[i]*fPolys_reverse[i]->df_dxj_dxk(Xtg, j, k); 
    
    //now, actually finding our new iteration is simple:
    RVecD Jx_Fx = Jx.Solve(Fx);
    
    if (Jx_Fx.size()<5) break; 

    Xtg += -Jx_Fx; 
    
    //check how far off our new iteration is
    rms = Compute_Xfp(Xtg).Distance(Xfp); 
    
    if (rms < min_error) break; 
  }
  
  return rms; 
} 
//_________________________________________________________________________________________
void TOpticsRays::Make_rays(const TXfp &Xfp, const TXtg &Xtg)
{
  fRays.clear(); 

  const double s = 0.100e-3; //the distance (in Xtg coords) between steps
  
  const double min_error = s/25;
  const double max_error = s/10;
  
  unsigned int pivotElem(0);
  
  fRays.push_back(Xtg); 
  
  //__________________________________________________________________________
  //
  // finds the next point.
  //
  // Returns the 'error' of the solver (correction step)
  // returns NaN if next point could not be found.
  //
  auto findNextPoint = [this,&pivotElem,&Xfp,&min_error](TXtg &Xtg, double s)
  {
    RMatrixD J_sub(4,4, 0.); //this is the jacobian, with one column 'moved' to the rhs
    J_sub.ReportSingular()=false; 
    RVecD J_pivot(4, 0.); 
    
    for (uint i=0; i<4; i++) { uint jj=0; 
      for (uint j=0; j<5; j++) { //loop over all 5 Xtg-coordinates, EXCEPT the pivot-element
	
	if (j==pivotElem) { //if j==pivotElem, fill out the pivot-vector
	  
	  J_pivot[i]     = this->fPolys_reverse[i]->df_dxj(Xtg, j);
	  
	} else            {
	  
	  J_sub.at(i,jj) = this->fPolys_reverse[i]->df_dxj(Xtg, j);
	  jj++; 
	}
      }
    }//for (uint i=0; i<4; i++)
    
    //solve this linear system
    RVecD dXtg = J_sub.Solve(J_pivot); 
    
    if (dXtg.size()==0) return std::nan("1"); 
    
    //now, put our 'guessed' element back in
    dXtg.insert( dXtg.begin()+pivotElem, -1. );
    
    dXtg = APEX::Unit(dXtg); //<- helper-function which normalizes an RVecD
    
    Xtg += s*dXtg; 
    
    return this->Find_Xtg(Xfp, Xtg, 50, min_error); 
  };
  
  //Info("Make_rays", "Starting.."); 
  
  TXtg Xtg_i = Xtg; 
  //find rays searching 'backward' from our starting place (Xtg)
  for (uint t=0; t<f_maxRays; t++) {
    if ( findNextPoint(Xtg_i, -s) > max_error ) break;
    fRays.push_back(Xtg_i);
  }
  
  //Info("Make_rays", "%i rays made (reverse)", (int)fRays.size() ); 
  
  //now, let's do something funky to reverse the order of all these rays
  uint i0(0), i1(fRays.size()-1);
  while (i0<i1) {
    TXtg temp = fRays[i1];
    fRays[i1] = fRays[i0];
    fRays[i0] = temp; 
    i0++;
    i1--;
  }
  
  Xtg_i = Xtg;
  
  //now, create rays in the 'forward' sense
  for (uint t=0; t<f_maxRays; t++) {
    if ( findNextPoint(Xtg_i,  s) > max_error ) break;
    fRays.push_back(Xtg_i);
  }

  //Info("Make_rays", "%i rays made (reverse+forward)", (int)fRays.size() ); 
  
  f_dInd = 1./((double)fRays.size()-1);

  //Info("Make_rays", "%i rays created.", (int)fRays.size()); 
}
//_________________________________________________________________________________________
TXtg TOpticsRays::Get_Xtg(double index_double) const
{
  if (fRays.size()<1) {
    Error("Get_Xtg", "No rays made. Must call TOpticsRays::Make_rays() first.");
    return TXtg();
  }
  
  int index = (int)floor(index_double/f_dInd); 
  
  if (index<0 || index>=(int)fRays.size()-1) {
    Error("Get_Xtg", "Index out-of-range (ind=%f). Must be 0 <= ind < 1.", index_double);
    return TXtg();
  }
  
  double modulo = index_double/f_dInd - index;

  const TXtg &ray1 = fRays.at(index);
  const TXtg &ray2 = fRays.at(index+1);
  
  return ray1 + (ray2 - ray1)*modulo; 
}
//_________________________________________________________________________________________
void TOpticsRays::Cast_Xtg_onto_target(const TXtg &Xtg,
				       double &x, double &y, double z) const
{
  TVector3 dir( Xtg.dxdz(), Xtg.dydz(), 1.        );
  TVector3 r0 ( Xtg.x(),    Xtg.y(),    f_Z_sieve ); 
  
  APEX::Coord::RotateTV3   ( dir, APEX::Coord::kTarget_to_Hall, f_isRHRS );
  APEX::Coord::TransformTV3( r0,  APEX::Coord::kTarget_to_Hall, f_isRHRS );

  x = r0.X() + (dir.X()/dir.Z())*( z - r0.Z() );
  y = r0.Y() + (dir.Y()/dir.Z())*( z - r0.Z() );
}
//_________________________________________________________________________________________

// Additional comments from senior dev muon
// 
//  w33333333333333333333333333333333333333333333333333333333333333333333333333333333333
//  3333333333

ClassImp(TRPoly);

