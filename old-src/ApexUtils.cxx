//*-- Author :    Ole Hansen   04-Dec-03

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// ApexUtils
//
// helper functions used by the APEX software.
// for now, all members will be defined within the APEX namespace. 
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <memory>
#include "TVectorD.h"
#include "TMatrixD.h"
#include "ApexUtils.h"
#include "RMatrixD.h" 
#include "ROOT/RVec.hxx"
#include "TFitResult.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TPad.h"
#include "TLatex.h"
#include "TClass.h"
#include "TString.h"
#include "TVector3.h"
#include "TLine.h"
#include "TClass.h"
#include "OpticsVectors.h"

using namespace std;
using namespace ROOT::VecOps; 

using RVecD = RVec<double>; 

//_____________________________________________________________________________
TVectorD APEX::SolveLinSystem(TMatrixD A, TVectorD B)
{ 
  //takes the NxN matrix A, and the N-vector B, and solves it. 
  //performs numerical LU factorization 
  
  //ASSUMES MATRIX IS NONSINGULAR!!!!
  const int N = B.GetNoElements(); 

  TVectorD X(N); 

  //check to make sure all our inputs are of a consistent size
  if (A.GetNcols() != N || 
      A.GetNrows() != N ) { 
    Error("APEX::SolveLinSystem",
	  "Mismatch in matrix/vector dims: A => %ix%i, B => %i. A must be square "
	  "and match the dims. of B.",
          A.GetNrows(), A.GetNcols(), N); 
    return X; 
  }
  
  //the U-matrix starts as a copy of the 'A' input-matrix
  TMatrixD U(A);  
  //we initialize the L-matrix with all zeros
  TMatrixD L(N,N);  for (int i=0; i<N; i++) for (int j=0; j<N; j++)   L(i,j)=0.; 

  //Create the U (upper triangular) and L (lower triangular) matrices
  for (int ii=0; ii<N; ii++) { 
    //loop over a (successivley smaller) sub-matrix
    double a_00 = U(ii,ii); 
    for (int i=ii+1; i<N; i++) { 
      double a_i0 = U(i,ii); 
      L(i,ii) = a_i0/a_00; 
      for (int j=ii; j<N; j++) { 
        U(i,j) += ( -a_i0/a_00 ) * U(ii,j); 
      }
    }
  }
  double det = 1.; 
  for (int ii=0; ii<N; ii++) { L(ii,ii)=1.; det *= U(ii,ii); }

  //det is NaN 
  if ( det != det ) {
    Error("fit_points::solveLinSystem()", "Determinant is NaN."); 
    return X; 
  }
  
  TVectorD y(N);

  //solve the system Ly = b
  for (int i=0; i<N; i++) { 
    y(i) = B(i); 
    for (int j=0; j<i; j++) y(i) += -L(i,j)*y(j);
  }
  
  //now solve Ux = y
  for (int i=N-1; i>=0; i--) { 
    X(i) = y(i); 
    for (int j=N-1; j>i; j--) { 
      X(i)  +=  -U(i,j) * X(j); 
    }
    X(i) *= 1./U(i,i); 
  }
  return X; 
}
//_____________________________________________________________________________
TVectorD APEX::TVectorD_cast(const RVec<double> &inVec)
{
  return TVectorD(inVec.size(),inVec.data()); 
}
//_____________________________________________________________________________
TVectorD APEX::TVectorD_cast(const vector<double> &inVec)
{
  return TVectorD(inVec.size(),inVec.data()); 
}
//_____________________________________________________________________________
RVec<double> APEX::RVecD_cast(const TVectorD &inVec)
{
  RVec<double> vec(*inVec.GetMatrixArray()); 
  return vec; 
}
//_____________________________________________________________________________
template <typename T> T* APEX::FindTypeInPad(TPad *pad)
{
  if (pad) {
    Error("ApexUtils::FindTypeInPad<T>",
	  "ptr to TPad is null! is there no canvas?");
    return 0; 
  }
  T *out=0;
  
  auto list = gPad->GetListOfPrimitives();
    
  for (TObject *obj : *list) 
    if (obj->IsA() == TClass::GetClass<T>()) out = (T*)obj;
    
  return out;
}
//_____________________________________________________________________________
TFitResult* APEX::Draw_gausFit( TH1 *hist,
                                double radius,
                                double center,
                                int color )
{
  //center the fit aruond the max. value of the histogram
  if (center < -1e29) { 
    center = hist->GetBinCenter( hist->GetMaximumBin() ); 
  }
  
  TF1 *gausFit = new TF1( "gausFit", "gaus(0) + [3]", 
			  -radius+center, 
			  radius+center ); 
  
  TAxis *xAxis = hist->GetXaxis(); 
  
  Double_t base = 0.;
  
  base += 0.5 * hist->GetBinContent( 1 ); 
  base += 0.5 * hist->GetBinContent( xAxis->GetNbins() ); 
    
  Double_t A = hist->GetMaximum() - base; 
  
  if (color >=0)
    cout << "A,base = " << A << ", " << base << endl; 
  
  gausFit->SetParameter( 0, A ); 
  gausFit->SetParameter( 1, center ); 
  gausFit->SetParameter( 2, radius/2 ); 
  gausFit->SetParameter( 3, base ); 
  
  TFitResultPtr fitPtr = hist->Fit("gausFit", "L N R S Q"); 
  
  double x0     = fitPtr->Parameter(1); 
  double x0_err = fitPtr->ParError(1); 
  
  double sig     = fitPtr->Parameter(2); 
  double sig_err = fitPtr->ParError(2); 
  
  A    = fitPtr->Parameter(0); 
  base = fitPtr->Parameter(3); 
  
  double yMax = hist->GetMaximum(); 
  double xMin = xAxis->GetBinCenter(1); 
  
  if (color>=0) {
    
    gausFit->SetLineColor(color); 
    gausFit->Draw("SAME");
    
    cout << "Parameters = p0( A ), p1( mean ), p2( sigma ), p3( const. )" << endl; 
  
    TF1 *gaus_copy = new TF1( (*gausFit) ); //copy the fit, to draw another one
    
  
    gaus_copy->SetRange( xAxis->GetXmin(), xAxis->GetXmax() ); 
    
    gaus_copy->SetLineStyle(kDotted);
    gaus_copy->Draw("SAME"); 
  
    cout << "sigma = " << sig << " +/- " << sig_err << endl; 
    cout << "mean  = " << x0  << " +/- " << x0_err  << endl; 
    cout << "A,base = " << A << ", " << base << endl; 
  
    TLatex *ltx = new TLatex; 
  
    ltx->DrawLatex(xMin*0.90, yMax, 
		   TString::Format("#bar{x} = %0.4e",x0)+
		   TString::Format(" #pm %0.4e",x0_err) ); 
  
    ltx->DrawLatex(xMin*0.90, yMax*0.92, 
		   TString::Format("#sigma = %0.4e",sig)+
		   TString::Format(" #pm %0.4e",sig_err) ); 
  }  

  return fitPtr.Get(); 
}
//_______________________________________________________________________________
void APEX::Coord::RotateTV3(TVector3 &r, Coord::ETransformType type, bool isRHRS) 
{
  //angles of TCS systems, offset from HCS centerline, in rad
  const double theta = isRHRS ? -0.093654368 : 0.093759087; 

  if (type==Coord::kHall_to_Target) { 
    r.RotateY(-theta);
    r.RotateZ(TMath::Pi()/2.); 
  } else {
    r.RotateZ(-TMath::Pi()/2.); 
    r.RotateY(theta); 
  }
}
//_______________________________________________________________________________
void APEX::Coord::TransformTV3(TVector3 &r, Coord::ETransformType type, bool isRHRS)
{
  TVector3 D0 = isRHRS 
    ? TVector3(-1.101e-3, -3.885e-3, 0.)
    : TVector3(-1.301e-3,  6.672e-3, 0.);

  if (type==Coord::kHall_to_Target) {
    Coord::RotateTV3(r,type,isRHRS);
    r += -D0;  
  } else {
    r +=  D0;
    Coord::RotateTV3(r,type,isRHRS);
  }
}
//_______________________________________________________________________________
void APEX::Coord::Transform_Xtg(TXtg &Xtg,
				APEX::Coord::ETransformType type,
				bool isRHRS)
{
  RVecD D0 = isRHRS 
    ? RVecD{ -1.101e-3, -3.885e-3, 0.7946 }
    : RVecD{ -1.301e-3,  6.672e-3, 0.7958 };

  RVecD dir = { Xtg.dxdz(), Xtg.dydz(), 1. }; 
  RVecD r0  = { Xtg.x(),    Xtg.y(),    0. }; 
  
  if (type==Coord::kTarget_to_Hall) { 
    RMatrixD A = isRHRS
      ? RMatrixD(3,3, { 0.000000,  0.995608,  0.093622,
		       -1.000000,  0.000000,  0.000000,
		       -0.000000, -0.093622,  0.995608 })
      : RMatrixD(3,3, { 0.000000,  0.995618, -0.093518,
		       -1.000000,  0.000000,  0.000000,
		        0.000000,  0.093518,  0.995618 }); 
    dir = A*dir; 

    r0  = A*(r0 + D0); 
    
  } //else {}
  
  
  
}
//_______________________________________________________________________________
void APEX::SetTGraphBounds(TGraph *g, const double xmin, const double ymin, const double xmax, const double ymax)
{
  //this is the best way I found to do this, which is frankly very irritating. ergo, I have this method now. 
  //you can have a TGraph drawn with these specific bounds using: 

  // using namespace APEX;
  // { 
  //    auto gr = new TGraph(...[your data]...);
  //    SetTGraphBounds(gr, xmin,ymin,...);
  //    gr->Draw("A L");  // A must be used; here 'L' means it will be drawn with lines. See TGraphPainter for more options
  // }


  if (!g) { Error("APEX::SetTGraphBounds", "TGraph is null!"); return; }  

  if (xmin>xmax) { Error("SetTGraphBounds", "Bad x-range {min,max} = {% .3e,% .3e}", xmin,xmax); return; }
  if (ymin>ymax) { Error("SetTGraphBounds", "Bad y-range {min,max} = {% .3e,% .3e}", ymin,ymax); return; }

  g->GetXaxis()->SetLimits(xmin, xmax); 
  g->GetHistogram()->SetMinimum(ymin);
  g->GetHistogram()->SetMaximum(ymax); 
}
//_______________________________________________________________________________
double APEX::ComputePol(const std::vector<double> *A, const double x)
{
  //T must be a std::vector<double>-like container of the coefficients. 
  // this function uses the sense: 
  // y = coeffs[0] + coeffs[1]*x + coeffs[2]*x*x + ...
  double val=0.; 
  for (auto a=A->end()-1; a>=A->begin(); a--) val = val*x + *a;
  return val; 
}
//_______________________________________________________________________________
double APEX::ComputePol(const ROOT::RVec<double> *A, const double x)
{
  double val=0.; 
  for (auto a=A->end()-1; a>=A->begin(); a--) val = val*x + *a;
  return val; 
}
//_______________________________________________________________________________
double APEX::ComputePol(const double *A, const double x, const int order)
{
  double val=0.; 
  for (int i=order; i>=0; i--) val = val*x + A[i];
  return val; 
}
//_______________________________________________________________________________
double APEX::ComputePol(const ROOT::RVec<double> &A, const double x)
{
  double val=0.; 
  for (int i=A.size()-1; i>=0; i--) val = val*x + A[i]; 
  return val; 
}
//_______________________________________________________________________________
bool APEX::hasNAN_RVecD(const RVecD &v)
{
  for (const double &x : v) if (x!=x) return true;
  return false; 
}
//_______________________________________________________________________________
double APEX::Length(const RVecD &v)
{
  double mod(0.);
  for (const double &x : v) mod += pow(x,2);
  return sqrt(mod); 
}
//_______________________________________________________________________________
double APEX::Dot(const RVecD &v1, const RVecD &v2)
{
  if (v1.size()!=v2.size()) {
    Error("Dot", "RVecD size mismatch v1.size()=%i, v2.size()=%i.",
	  (int)v1.size(), (int)v2.size());
    return -1e30;
  }
  
  double val(0.);
  for (unsigned int i=0; i<v1.size(); i++) val += v1[i]*v2[i];
  return val;
}
//_______________________________________________________________________________
RVecD APEX::Unit(const RVecD &v)
{
  double mod(0.);
  for (const double &x : v) mod += pow(x,2); 
  return v/sqrt(mod); 
}
//_______________________________________________________________________________
void APEX::Graphic::VLine(TH1 *hist, double x, unsigned int color, unsigned int style) 
{
  if (!hist) { Error("VLine", "hist passed is null"); return; }

  auto line = unique_ptr<TLine>(new TLine); 

  double y0,y1; 
  //is 2D-histogram
  if (hist->InheritsFrom(TH2::Class())) {
    y0 = hist->GetYaxis()->GetXmin(); 
    y1 = hist->GetYaxis()->GetXmax(); 
  
  } else { //is a 1D histogram

    y0 = min<double>( hist->GetMinimum(), 0. );
    y1 = hist->GetXaxis()->GetXmax();  
  }
  line->SetLineColor(color); line->SetLineStyle(style); 

  line->DrawLine(x,y0, x,y1); 
}
//_______________________________________________________________________________
void APEX::Graphic::HLine(TH1 *hist, double y, unsigned int color, unsigned int style) 
{
  if (!hist) { Error("HLine", "hist passed is null"); return; }

  auto line = unique_ptr<TLine>(new TLine); 

  double x0,x1; 
  x0 = hist->GetXaxis()->GetXmin(); 
  x1 = hist->GetXaxis()->GetXmax(); 
  
  line->SetLineColor(color); line->SetLineStyle(style); 

  line->DrawLine(x0,y, x1,y); 
}
