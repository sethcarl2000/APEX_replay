#include "def_apex.h"
#include "function_lib.h"
#include <fstream>
#include <typeinfo>
#include <limits> 

class function_lib;

using namespace std; 

//////////////////////////////////////////////////////////////////////////////////
//
//  These are a bunch of little functions which mostly have to do with the 
//  geometry of the HRSs for the APEX setup. 
//  
//  a more future-proof way (if the detectors are in a different position)
//  is to actually ask the Analyzer's detector objects and classes where 
//  the detectors are. However, this requires booting up & going through an
//  analyzer run, if only for a few events. I'm too lasy to do that each time
//  i want to run some analysis on replay data, ergo this code. 
//  
//  - Seth 16 sept. 23
//
/////////////////////////////////////////////////////////////////////////////////
template <typename T>
T* FindTypeInPad()
{
  if (!gPad) {
    Error("function_lib::FindTypeInPad<T>()",
	  "ptr 'gPad' is null! is there no canvas?");
    return 0; 
  }
  T *out=0; 
  
  auto list = gPad->GetListOfPrimitives();
  
  for (TObject *obj : *list) 
    if (obj->IsA() == TClass::GetClass<T>()) out = (T*)obj;
  
  return out;
}     
void setArm(bool armChoice) { 
  //choices are either kArm_Right (true) or kArm_Left (false)
  isRHRS = armChoice;
  return; 
} 

Bool_t checkArm() { return isRHRS; }   
void xyz_to_uvw(Double_t x, 
		Double_t y, 
		Double_t z, 
		Double_t &u, 
		Double_t &v, 
		Double_t &w)
{
  TArrayD uData(9);
  
  Double_t a = TMath::Power(2, -1);
  Double_t b = TMath::Power(2, -0.5);
    
  uData[0] =  a; uData[1] = -b; uData[2] =  a; //xyz -> uvw
  uData[3] =  a; uData[4] =  b; uData[5] =  a;
  uData[6] = -b; uData[7] =  0; uData[8] =  b;
  
  TMatrixD U(3, 3); 
  U.SetMatrixArray( uData.GetArray() ); 
  
  TVectorD r(3); r[0] =x; r[1] =y; r[2] =z; 
  
  TVectorD rPrime = U*r; 
  
  u = rPrime[0];
  v = rPrime[1];
  w = rPrime[2];

  return; 
}
void uvw_to_xyz(Double_t u, 
		Double_t v, 
		Double_t w, 
		Double_t &x, 
		Double_t &y, 
		Double_t &z)
{
  TArrayD uData(9);
  
  Double_t a = TMath::Power(2, -1);
  Double_t b = TMath::Power(2, -0.5);
  
  uData[0] =  a; uData[1] =  a; uData[2] = -b; //uvw -> xyz
  uData[3] = -b; uData[4] =  b; uData[5] =  0;
  uData[6] =  a; uData[7] =  a; uData[8] =  b;
  
  TMatrixD U(3, 3); 
  U.SetMatrixArray( uData.GetArray() ); 
  
  TVectorD r(3); r[0] =u; r[1] =v; r[2] =w; 
  
  TVectorD rPrime = U*r; 
  
  x = rPrime[0];
  y = rPrime[1];
  z = rPrime[2];

  return; 
}
TVector3 vec_xyz_uvw( TVector3 xyz )
{
  Double_t 
    x( xyz.X() ),
    y( xyz.Y() ),
    z( xyz.Z() );  
  
  Double_t u,v,w; 
    
  xyz_to_uvw(x,y,z, u,v,w); 
  
  TVector3 uvw( u,v,w ); 
  
  return uvw; 
}
TVector3 vec_uvw_xyz( TVector3 uvw )
{
  Double_t 
    u( uvw.X() ),
    v( uvw.Y() ),
    w( uvw.Z() );  
  
  Double_t x,y,z; 
    
  uvw_to_xyz(u,v,w, x,y,z); 
  
  TVector3 xyz( x,y,z ); 
  
  return xyz; 
}
Double_t wirePos(int p, Int_t wire) 
{
  if (isRHRS) {
    
    Double_t v0[4] = {0.77852, 0.77852, 1.02793, 1.02793}; 
    return v0[p] - dWire*((double)wire);
    
  } else { 
    
    Double_t v0[4] = {0.77852, 0.77852, 1.02718, 1.02718};
    return v0[p] - dWire*((double)wire);
    
  }
}
Int_t wireNum(int p, Double_t pos) 
{ 
  if (isRHRS) {
    Double_t v0[4] = {0.77852, 0.77852, 1.02793, 1.02793}; 
    int num = (int)TMath::Nint( (v0[p]-pos)/dWire ); 
    
    return num;
  
  } else { 
    Double_t v0[4] = {0.77852, 0.77852, 1.02718, 1.02718};
    int num = (int)TMath::Nint( (v0[p]-pos)/dWire ); 
    
    return num;
  
  }
}
Double_t vdc_Toff(int plane)
{
  //this works via = 
  Double_t Rt_off[4] = {48.1e-9, 48.3e-9, 38.8e-9, 41.5e-9}; 
  Double_t Lt_off[4] = {36.6e-9, 36.3e-9, 29.5e-9, 29.6e-9};  
  
  if(isRHRS) { return Rt_off[plane]; }
  else       { return Lt_off[plane]; }
}
Double_t vdc_w(int plane) {
  
  if (isRHRS) { 
    
    Double_t my_w[4] = { 0, 0.026, 0.333245, 0.359245 };
    return my_w[plane]; 
    
  } else {
    
    Double_t my_w[4] = { 0, 0.026, 0.335344, 0.361444 };
    return my_w[plane]; 
  
  }

}
Double_t wL() { return 0.5*(vdc_w(0)+vdc_w(1)); }
Double_t wU() { return 0.5*(vdc_w(2)+vdc_w(3)); }
Double_t wSep() { return vdc_w(2)-vdc_w(0); }

Bool_t s2_isCoinc(Double_t rt[16], Double_t lt[16], int pad) 
{
  
  if (TMath::Abs(rt[pad]) > 1e10 || TMath::Abs(lt[pad]) > 1e10) 
    { return false; } 
  
  Float_t sep = rt[pad] - lt[pad] + (!checkArm())*55;
  
  if (TMath::Abs(sep) > 25) 
    { return false; }
  
  return true; 
}
TVector2 leastSquares_points( vector<TVector2>  p_in,
			      vector<TVector3> &p_out,
			      double y_error ) 
{
  //take a set of input points (data) ordered by 
  //
  // p_in = {xi, yi, dyi}         ~~where:
  // 
  //  xi = data-point x
  //  yi = data-point y
  // dyi = y-uncertainty of i-th point
  
  // Then, given x-points in 'predict', it will return those predictions 
  //       in the same format, given that the x-values have been filled out. 
  
  
  const int N0 = p_in .size(); 
  const int N1 = p_out.size(); 
  double n = (double)N0; 
  
  Double_t Sxy(0), Sx(0), Sy(0), Sxx(0); 
  
  //do least-square linear fit, return R^2 error
  for (int j=0; j<N0; j++) { 
    
    Double_t x(p_in[j].X()), y(p_in[j].Y()); 
    
    Sx  += x;
    Sy  += y;
    Sxx += x*x; 
    Sxy += x*y;
  }
  
  double denom = n*Sxx - Sx*Sx;
    
  double m = ( n*Sxy - Sx*Sy )/( n*Sxx - Sx*Sx );
  
  double b = (Sy - m*Sx)/n;
  
  
  for (int j=0; j<N1; j++) { 
    
    double X  =  p_out[j].X(); 
    double Y  =  m*X + b; 
    
    double dY = 0; 
    
    for (int jj=0; jj<N0; jj++) { 
      
      double x = p_in[jj].X(); 
      
      double dY_dyi = ((n*x - Sx)/denom)*( X - Sx/n )  +  1.0/n; 
      
      dY += TMath::Power( dY_dyi, 2.0 ); 
          
    }
        
    dY = y_error * TMath::Sqrt( dY ); 
    
    //cout << X << "  " << Y << "  " << dY << endl; 
    
    p_out[j] = TVector3( X, Y, dY ); 
  }  
  
  return TVector2( m, b ); 
}
TVector3 leastSquares_quad( vector<TVector2> pts ) 
{ 
  double Sx[4] = {0,0,0,0}; double Sy[3] = {0,0,0};
  
  for (int h=0; h<pts.size(); h++ ) {
    
    double x = pts[h].X();
    double y = pts[h].Y(); 
    
    Sx[0] += x;
    Sx[1] += x*x; 
    Sx[2] += x*x*x; 
    Sx[3] += x*x*x*x; 
    
    Sy[0] += y;
    Sy[1] += y*x;
    Sy[2] += y*x*x; 
  }
  
  TArrayD Ad(9); 
  
  Ad[0]=Sx[3];  Ad[1]=Sx[2];  Ad[2]=Sx[1];
  Ad[3]=Sx[2];  Ad[4]=Sx[1];  Ad[5]=Sx[0];
  Ad[6]=Sx[1];  Ad[7]=Sx[0];  Ad[8]=((double)pts.size());
  
  TMatrixD A( 3,3 ); A.SetMatrixArray( Ad.GetArray() ); 
  
  TVectorD bVec(3); 
  bVec[0]=Sy[2]; bVec[1]=Sy[1]; bVec[2]=Sy[0];
  
  bVec = ( A.Invert() )*bVec; 
  
  return TVector3( bVec[0],bVec[1],bVec[2] ); 
}
TVector3 leastSquares_quad_error( vector<TVector3> pts, TVector3 &error ) 
{
  //the input vector<TVector3> pts:
  //   pts.at(i)[0] - x
  //   pts.at(i)[1] - y
  //   pts.at(i)[2] - dy (y-error of this point)
  
  double Sx[4] = {0,0,0,0}; double Sy[3] = {0,0,0};
  
  for (int h=0; h<pts.size(); h++ ) {
    
    double x = pts[h].X();
    double y = pts[h].Y(); 
    
    Sx[0] += x;
    Sx[1] += x*x; 
    Sx[2] += x*x*x; 
    Sx[3] += x*x*x*x; 
    
    Sy[0] += y;
    Sy[1] += y*x;
    Sy[2] += y*x*x; 
  }
  
  TArrayD Ad(9); 
  
  Ad[0]=Sx[3];  Ad[1]=Sx[2];  Ad[2]=Sx[1];
  Ad[3]=Sx[2];  Ad[4]=Sx[1];  Ad[5]=Sx[0];
  Ad[6]=Sx[1];  Ad[7]=Sx[0];  Ad[8]=((double)pts.size());
  
  TMatrixD A( 3,3 ); A.SetMatrixArray( Ad.GetArray() ); 
  
  TVectorD bVec(3); 
  bVec[0]=Sy[2]; bVec[1]=Sy[1]; bVec[2]=Sy[0];
  
  A.Invert(); 
  
  bVec = A*bVec; 
  
  double a(bVec[0]), b(bVec[1]), c(bVec[2]); 
  
  TVectorD xi_vec(3); xi_vec[0]=0; xi_vec[1]=0; xi_vec[2]=0; 
  
  double *inv_a = A.GetMatrixArray(); 
  
  double dA_dyi(0), dX_dyi(0), dY_dyi(0);
  
  for (int h=0; h<pts.size(); h++) { 
    
    TVectorD bVec_i(3); 
    
    double x = pts[h].X();
    
    double dy = pts[h].Z(); 
    
    bVec_i[0] = x*x;
    bVec_i[1] = x;
    bVec_i[2] = 1.0; 
    
    bVec_i = A*bVec_i; 
    
    double dA = bVec_i[0]; 
    
    double dX = (bVec_i[0]*b/a  -  bVec_i[1])/(2.0*a);
    
    double dY = bVec_i[2]  -  b*bVec_i[1]/(2.0*a)  +  b*b*bVec_i[0]/(4.0*a*a); 
    
    dA_dyi += TMath::Power( dA * dy, 2 ); 
    dX_dyi += TMath::Power( dX * dy, 2 );
    dY_dyi += TMath::Power( dY * dy, 2 ); 
  }
  
  //compute error vec
  
  error = TVector3( TMath::Sqrt(dA_dyi),
		    TMath::Sqrt(dX_dyi),
		    TMath::Sqrt(dY_dyi) ); 
  
  //returning { A, x0, y0 }
  return TVector3( a, -b/(2.0*a), c - (b*b)/(4.0*a) ); 
  
}
Double_t drift_X(Double_t tau) 
{
  Double_t V_drift, B_drift, x0_drift, t0_drift;  
  
  if (checkArm()==kArm_Right) {
    V_drift  = V_drift_R; 
    B_drift  = B_drift_R;  
    x0_drift = x0_drift_R; 
    t0_drift = t0_drift_R;
  } else                      { 
    V_drift  = V_drift_L; 
    B_drift  = B_drift_L;  
    x0_drift = x0_drift_L; 
    t0_drift = t0_drift_L;
  }
  
  if ( tau > t0_drift ) { return V_drift*tau + B_drift;   }
  else                  { return (x0_drift/t0_drift)*tau; }
}
double drift_X_new( double tau, double slope ) 
{
  return drift_X(tau); //temporary substitution
  
  if (tau < 50e-9) {  //low-t region
    
    double dw=0; 
    
    for (int i=0; i<=lowT_order; i++) { 
      for (int j=0; j<=lowT_order; j++) { 
	dw +=  param_lowT[i][j] * TMath::Power(tau/tauScale,i)*TMath::Power(slope-1.5,j); 
      }
    }
    return dw; 
  }
  
  //'linear'-t region
  return  k0 + kT*tau + kM*slope; 
}
Double_t drift_T(Double_t x) 
{  
  Double_t V_drift, B_drift, x0_drift, t0_drift;  
  
  if (checkArm()==kArm_Right) {
    V_drift  = V_drift_R; 
    B_drift  = B_drift_R;  
    x0_drift = x0_drift_R; 
    t0_drift = t0_drift_R;
  } else                      { 
    V_drift  = V_drift_L; 
    B_drift  = B_drift_L;  
    x0_drift = x0_drift_L; 
    t0_drift = t0_drift_L;
  }
    
  x = TMath::Abs(x);
  if ( x > x0_drift ) { return (x-B_drift)/V_drift;   } 
  else                { return (t0_drift/x0_drift)*x; }
}
Double_t s2_t_vdc(Double_t t_s2) 
{
  Double_t T_s = (checkArm()) ? 580.23e-9 : 596.73e-9; 
  
  return t_s2 = T_s - t_s2;
}
double drift_T_new( double dw, double slope ) 
{
  return drift_T(dw); //temporary substition 

  
  dw = TMath::Abs(dw);
  
  if (dw < 3.5e-3) {//low-x region
    
    double tau=0; 
    
    for (int i=0; i<=lowT_order; i++) { 
      for (int j=0; j<=lowT_order; j++) { 
	tau +=  param_lowX[i][j] * TMath::Power(dw/dwScale,i)*TMath::Power(slope-1.5,j); 
      }
    }
    return tau; 
  }
  
  //'linear' region
  return  ( dw -k0 -kM*slope )/kT; //note this is the inverse of the above formula
}
Double_t s2_x(int pad) {
  //return the x-pos of the i-th paddle [0-15]
  
  Double_t pad_width = 0.13975; //in m
    
  if (isRHRS) {
    Double_t pad_x0 = -1.21413;
    
    return pad_x0 + ((double)pad)*pad_width;
  
  } else { 
    Double_t pad_x0 = -1.16913;
    
    return pad_x0 + ((double)pad)*pad_width;
  } 
  
}
Double_t s2_v(int pad) {
  //all these values are calculated by applyin the xyz_to_uvw to the values of the aapove functions.
  Double_t pad_vWidth = 6.9875e-2; 
    
  if (isRHRS) {
    Double_t pad_v0 = -1.21413;
    return pad_v0 + ((double)pad)*pad_vWidth;
  } else      {
    Double_t pad_v0 = -1.16913;
    return pad_v0 + ((double)pad)*pad_vWidth;
  }
}
Double_t s2_padWidth() { return 0.13975; } 
Double_t s2_x_offset() 
{ 
  if (isRHRS) { return -0.166; }
  else        { return -0.121; }
}   
Double_t s2_z() { 
  if (isRHRS) { return 3.3098; }
  else        { return 3.1790; }
}
Double_t s2_realTime(bool isPMT_R, int paddle, Double_t rawTDC[16]) {
  
  //choices for isPMT_R = 
  // kPMT_R (true)
  // kPMT_L (true)
  

 Double_t Right_R[16] 
   = {2551.36, 2550.84, 2552.93, 2552.42, 2553.03, 2553.02, 2552.42, 2553.06, 2555.18, 2553.52, 2554.33, 2552.95, 2549.92, 2550.57, 2551.83, 2548.26};

Double_t Right_L[16] 
  = {2553.56, 2547.79, 2550.14, 2549.62, 2552.32, 2549.35, 2554.73, 2552.29, 2550.83, 2547.63, 2551.22, 2551.11, 2552.26, 2549.53, 2548.21, 2551.16};

Double_t Left_R[16] 
  = {2759.06, 2820.8, 2815.81, 2819.58, 2821.29, 2818.89, 2821.84, 2821.3, 2821.87, 2822.54, 2822.36, 2821.5, 2820.98, 2821.42, 2820.59, 2759.86};

Double_t Left_L[16] 
  = {2817.96, 2760.24, 2763.8, 2760.99, 2764.22, 2764.18, 2765.38, 2767.82, 2764.86, 2763.6, 2764.22, 2758.28, 2759.5, 2759.49, 2763.37, 2821.76};

  
  
  Double_t resolution = 0.5e-9; //resolution of tdc, half a ns.
  
  //re-do the s2m TDC calculations (without the time-walk corrections)
  Double_t RHRS_R[16] = {1390.9, 1530.88, 1532.87, 1531.96, 1532.17, 1531.76, 1530.96, 1530.9, 1532.62, 1530.46, 1531.27, 1529.89, 1526.86, 1527.51, 1528.87, 1387.6}; 
  
  Double_t RHRS_L[16] = {1393.1, 1527.83,  1530.08,  1529.16,  1531.46,  1528.09,  1533.27,  1530.13,  1528.27,  1524.57,  1528.16,  1528.05,  1529.2,  1526.47,  1525.25,  1390.5};
  
  Double_t LHRS_R[16] = {1565.6, 1576.54,  1572.45,  1575.92,  1577.13,  1574.53,  1577.68,  1576.04,  1577.11,  1578.08,  1578.3,  1577.44,  1577.32,  1577.86,  1576.93,  1570.6};
  
  Double_t LHRS_L[16] = {1624.5,  1515.98,  1520.44,  1517.33,  1520.06,  1519.82,  1521.22,  1522.56,  1520.1,  1519.14,  1520.16,  1514.22,  1515.84,  1515.93,  1519.71, 1632.5};
    
  Double_t correct_RHRS[16] = {0, 140.5, 140.4, 140, 139.6, 139.2, 139, 138.3, 137.9, 137.4, 137.4, 137.4, 137.4, 137.4, 137.5, -0.2}; 

  Double_t correct_LHRS[16] = {0, -50.8, -49.9, -50.2, -50.7, -50.9, -50.7, -51.8, -51.3, -51, -50.6, -50.6, -50.2, -50.1, -50.2, 4.2};
  
  if (checkArm()==kArm_Right) { //RHRS 
    if (isPMT_R) {  
      return (Right_R[paddle] - rawTDC[paddle])*resolution; 
    } else    { 
      return (Right_L[paddle] - rawTDC[paddle])*resolution; 
    }
  } else      { 
    if (isPMT_R) { 
      return (Left_R[paddle]  - rawTDC[paddle])*resolution;
    } else    { 
      return (Left_L[paddle]  - rawTDC[paddle])*resolution; 
    }
  }
  
}
Double_t s2_coincTime(Double_t tR[16], Double_t tL[16], int pad) 
{    
  if (!s2_isCoinc(tR, tL, pad)) return -1e38;
  
  Double_t timeR = s2_realTime(kPMT_R, pad, tR);
  Double_t timeL = s2_realTime(kPMT_L, pad, tL);
  
  return 0.5*( timeR + timeL ); 
}
double s2_paddleTime( double rt, double lt, int pad, bool arm) 
{
  if (TMath::Abs(rt) > 1e10 || TMath::Abs(lt) > 1e10) return -1e30; 
  
  Float_t sep = rt - lt + (!arm)*55;
  
  if (TMath::Abs(sep) > 25) return -1e30; 
  
  Double_t Right_R[16] 
   = {2551.36, 2550.84, 2552.93, 2552.42, 2553.03, 2553.02, 2552.42, 2553.06, 2555.18, 2553.52, 2554.33, 2552.95, 2549.92, 2550.57, 2551.83, 2548.26};

Double_t Right_L[16] 
  = {2553.56, 2547.79, 2550.14, 2549.62, 2552.32, 2549.35, 2554.73, 2552.29, 2550.83, 2547.63, 2551.22, 2551.11, 2552.26, 2549.53, 2548.21, 2551.16};

Double_t Left_R[16] 
  = {2759.06, 2820.8, 2815.81, 2819.58, 2821.29, 2818.89, 2821.84, 2821.3, 2821.87, 2822.54, 2822.36, 2821.5, 2820.98, 2821.42, 2820.59, 2759.86};

Double_t Left_L[16] 
  = {2817.96, 2760.24, 2763.8, 2760.99, 2764.22, 2764.18, 2765.38, 2767.82, 2764.86, 2763.6, 2764.22, 2758.28, 2759.5, 2759.49, 2763.37, 2821.76};

  Double_t resolution = 0.5e-9; //resolution of tdc, half a ns.
  
  //re-do the s2m TDC calculations (without the time-walk corrections)
  Double_t RHRS_R[16] = {1390.9, 1530.88, 1532.87, 1531.96, 1532.17, 1531.76, 1530.96, 1530.9, 1532.62, 1530.46, 1531.27, 1529.89, 1526.86, 1527.51, 1528.87, 1387.6}; 
  
  Double_t RHRS_L[16] = {1393.1, 1527.83,  1530.08,  1529.16,  1531.46,  1528.09,  1533.27,  1530.13,  1528.27,  1524.57,  1528.16,  1528.05,  1529.2,  1526.47,  1525.25,  1390.5};
  
  Double_t LHRS_R[16] = {1565.6, 1576.54,  1572.45,  1575.92,  1577.13,  1574.53,  1577.68,  1576.04,  1577.11,  1578.08,  1578.3,  1577.44,  1577.32,  1577.86,  1576.93,  1570.6};
  
  Double_t LHRS_L[16] = {1624.5,  1515.98,  1520.44,  1517.33,  1520.06,  1519.82,  1521.22,  1522.56,  1520.1,  1519.14,  1520.16,  1514.22,  1515.84,  1515.93,  1519.71, 1632.5};
    
  Double_t correct_RHRS[16] = {0, 140.5, 140.4, 140, 139.6, 139.2, 139, 138.3, 137.9, 137.4, 137.4, 137.4, 137.4, 137.4, 137.5, -0.2}; 
  
  Double_t correct_LHRS[16] = {0, -50.8, -49.9, -50.2, -50.7, -50.9, -50.7, -51.8, -51.3, -51, -50.6, -50.6, -50.2, -50.1, -50.2, 4.2};
 
  double Rpmt_REAL, Lpmt_REAL; 
  
  if (arm==kArm_Right) { //RHRS 
    
    Rpmt_REAL = (Right_R[pad] - rt)*resolution; 
    Lpmt_REAL = (Right_L[pad] - lt)*resolution; 
  
  } else      { 
    
    Rpmt_REAL = (Left_R[pad]  - rt)*resolution; 
    Lpmt_REAL = (Left_L[pad]  - lt)*resolution; 
  }
  
  return 0.5*(Rpmt_REAL + Lpmt_REAL); 
}
Bool_t s2_coincCut(Double_t Rtr[16], Double_t Rtl[16], 
		   Double_t Ltr[16], Double_t Ltl[16],
		   int &pR, 
		   int &pL,
		   Double_t &time_R, 
		   Double_t &time_L,
		   const Double_t stdv_cutoff=2.5) //how many stdvs to cut off?
{
  //return the fist s2m coinc time found (if any) 
  
  bool start_arm = checkArm(); 
  
  for (int r=0; r<16; r++) { 
    for (int l=0; l<16; l++) { 
      
      setArm(kArm_Right); 
      if (!s2_isCoinc( Rtr,Rtl, r )) continue; 
      Double_t time_R_test = s2_coincTime( Rtr,Rtl, r ); 
      
      setArm(kArm_Left); 
      if (!s2_isCoinc( Ltr,Ltl, l)) continue; 
      Double_t time_L_test = s2_coincTime( Ltr,Ltl, l ); 
      
      Double_t T_coinc = time_R_test - time_L_test; 
      
      Double_t stDev = TMath::Abs( RUN_Tcoinc - T_coinc )/RUN_Tcoinc_sigma; 
      
      if ( stDev < stdv_cutoff ) {
	
	time_R = time_R_test;
	time_L = time_L_test; 
	pR = r;
	pL = l;
	
	setArm(start_arm); 
	return true; 
      }
      
    }// for (int l=0; 
  }// for (int r=0; 	
  
  setArm(start_arm); 
  return false; 
}
void mlp_cluster_coarse(Double_t t_lower[11], 
			Double_t t_upper[11],
			Int_t wLower,
			Int_t wUpper,
			Double_t &int1,
			Double_t &int2)
{ 
  
  if (!mlp_coarse) {
    TFile *f = new TFile("mlp/mlp_cluster_coarse.root", "READ");
    mlp_coarse = (TMultiLayerPerceptron*)f->Get("mlp_cluster"); 
  }
  
  bool right=true;
  Double_t w_sep = 0.333245; //w-separation between like vdc planes 
  Double_t time_const = 400e-9; 
  Double_t inputs[23];
  
  inputs[0] = ( wirePos(2, wUpper) - wirePos(0, wLower) )/w_sep;
  
  for (int j=0; j<11; j++) {
    inputs[1+j]  = (TMath::Abs(t_lower[10-j])>1) ? 0 : 1. - (t_lower[10-j]/time_const);
    inputs[12+j] = (TMath::Abs(t_upper[10-j])>1) ? 0 : 1. - (t_upper[10-j]/time_const);
  }
  
  //for (int j=0; j<23; j++) {cout << i << " " << inputs[j] << endl;}
  
  int1 = mlp_coarse->Evaluate(0, inputs);
  int2 = mlp_coarse->Evaluate(1, inputs); 
  
  int1 = (5*dWire*int1) + wirePos(0, wLower);
  int2 = (5*dWire*int2) + wirePos(2, wUpper); 
  
  return; 
}
TGraph* th2d_to_graph(TH2D *h, 
		      Double_t xMin, 
		      Double_t xMax,
		      Double_t fit_radius=-1,//restrict the fit to a certain radius
		      Double_t sigma_guess=-1, 
		      bool is_errorGraph=false)//try this sigma as a guess  
{
  //takes a 2d histogram, with a specified x-range, and fits a gaus dist to each
  // y-bin. then, adds a point on a tgraph corresponding to each point
  
  const int nBins = h->GetNbinsX(); 
  
  TAxis *xAxis = (TAxis*)h->GetXaxis(); 
  TAxis *yAxis = (TAxis*)h->GetYaxis(); 
  
  int bin0 = xAxis->FindBin( xMin ); if (bin0<1)     return 0; 
  int bin1 = xAxis->FindBin( xMax ); if (bin1>nBins) return 0;
    
  const int nPoints = bin1 - bin0 + 1;
  Double_t 
    px[nPoints], 
    py[nPoints],
    p_sigma[nPoints]; 
    
  Double_t y0, y1; 
  y0 = yAxis->GetXmin();
  y1 = yAxis->GetXmax();
  
  TF1 *gausFit = new TF1("gausFit", "gaus", y0, y1); 
  
  Double_t mean_sigma(0); 
  
  //start taking sections
  for (int bin=bin0; bin<=bin1; bin++) { 
    
    TH1D *proj = h->ProjectionY("proj", bin, bin); 
    
    Double_t sigma = proj->GetRMS(); 
    if (sigma_guess>0) sigma=sigma_guess; //if we want a custom guess for sigma 
    
    Double_t mean  = yAxis->GetBinCenter( proj->GetMaximumBin() ); 
    Double_t A     = proj->GetMaximum(); 
    
    
    Double_t fitRange_low = (mean-fit_radius > y0) ? mean-fit_radius : y0 ; 
    Double_t fitRange_hi  = (mean+fit_radius < y1) ? mean+fit_radius : y1 ; 
    
    
    if (fit_radius>0) {//restric the fit radius about the mean
      gausFit = new TF1("gausFit", "gaus", 
			fitRange_low,
			fitRange_hi ); 
    }      
    
    gausFit->SetParameters(A, mean, sigma); 
    
    proj->Fit("gausFit", "Q N R L"); 
    //Q means no printed output, N means don't plot, R means only fit 'in-range'
    
    px[bin-bin0]      = xAxis->GetBinCenter( bin ); 
    py[bin-bin0]      = gausFit->GetParameter(1); //the 'mean'
    
    p_sigma[bin-bin0] = gausFit->GetParameter(2); 
    
    mean_sigma += gausFit->GetParameter(2); 
    
    proj->~TH1D(); //if we don't delete this here, root gets irritated 
    
    if (fit_radius>0) { gausFit->~TF1(); } //same case here. 
  }
  
  cout << "mean sigma = " << mean_sigma/((float)nPoints) << endl; 

  
  TGraph *graph_y     = new TGraph(nPoints, px, py); 
  
  TGraph *sigma_graph = new TGraph(nPoints, px, p_sigma); 
  
  if (is_errorGraph) { return sigma_graph; }
  else               { return graph_y;     } 
}
double get_Z_value( TH1D *h, double centile, bool fromLeft=true ) 
{ 
  //get the z-value, along the x-axis of h. 
  //z-value, defined for some centile, is when 'centile'=fraction of a plot is 
  // to the left of that z-value. 

  auto axis = h->GetXaxis(); 
  
  const int nBins = axis->GetNbins(); 
  
  double dx = axis->GetBinWidth(1); 
  
  double sum(0.), Target(h->Integral()*centile); 

  int bin=0; double dBin; 
  
  while ( bin<=nBins && sum<Target ) {
    
    bin++; dBin = h->GetBinContent(bin); 
    sum += dBin; 
  }
  
  return axis->GetBinCenter(bin)  +  dx*(Target-(sum-dBin))/dBin;
}
TGraph* th2d_to_graph_percentile(TH2D *h, 
				 Double_t xMin, 
				 Double_t xMax,
				 Double_t percentile=0.500)
{
  //takes a 2d histogram, with a specified x-range, and fits a gaus dist to each
  // y-bin. then, adds a point on a tgraph corresponding to each point
  
  TAxis *xAxis = (TAxis*)h->GetXaxis(); 
  TAxis *yAxis = (TAxis*)h->GetYaxis(); int nBinsY = yAxis->GetNbins(); 
  
  int bin0 = xAxis->FindBin( xMin ); 
  int bin1 = xAxis->FindBin( xMax );  
  
  const int nPoints = bin1 - bin0 + 1;
  Double_t 
    px[nPoints], 
    py[nPoints]; 
  
  Double_t y0, y1; 
  y0 = yAxis->GetXmin();
  y1 = yAxis->GetXmax();
  
  TF1 *gausFit = new TF1("gausFit", "gaus", y0, y1); 
  
  Double_t mean_sigma(0); 
  
  //start taking sections
  for (int bx=bin0; bx<=bin1; bx++) { 
    
    TH1D *proj = h->ProjectionY("proj", bx, bx); 
    
    double net = proj->Integral(); 
    double sum = 0; 
    
    double x = xAxis->GetBinCenter(bx); 
    
    int by=1; 
    while  ( by <= nBinsY  &&  sum/net < percentile ) {
      sum += proj->GetBinContent(by);
      by++; 
    }
    
    py[bx-bin0]  =  yAxis->GetBinCenter(by); 
    px[bx-bin0]  =  xAxis->GetBinCenter(bx); 
    
    proj->~TH1D(); //if we don't delete this here, root gets irritated 
  }
    
  TGraph *graph = new TGraph(nPoints, px, py); 
  
  return graph; 
}
Double_t vdc_realTime(int p, int w, Double_t rawTime) 
{
  /*/vdc wire offsets
  
  Double_t offset = (isRHRS) 
    ? R_abs[p][w]-1.04e-9 : L_abs[p][w]+3.20e-9; 
    
    return (-0.5e-9)*rawTime + offset; */ return 1.; 
}
Double_t vdc_rawTime(int p, int w, Double_t realTime)
{
  /*/convert 'real' time back into 'raw' time (for monte-carlo generation) 
  Double_t offset = (isRHRS) 
    ? R_abs[p][w]-1.04e-9 : L_abs[p][w]+3.20e-9; 
  
  return ( realTime - offset )/(-0.5e-9); */ return 1.; 
}
TVector3 compute_focalPlane_intercept(Double_t intercept[4], int axis=kUVW) 
{
  Double_t Sv = (intercept[2]-intercept[0])/wSep(); 
  Double_t Su = (intercept[3]-intercept[1])/wSep(); 
  
  Double_t v0 = intercept[0] - vdc_w(0)*Sv; 
  Double_t u0 = intercept[1] - vdc_w(1)*Su; 
  
  TVector3 r0_uvw( u0,v0,0.0 );  // = new TVector3( u0,v0,0. ); 
  
  if (axis==kUVW) { 
    return r0_uvw; //return the point in UVW coords (default) 
  } else          { 
    return vec_uvw_xyz(r0_uvw); //return the point in XYZ coords
  }
}
TVector3 compute_s2_intercept(Double_t intercept[4], int axis=kXYZ) 
{ 
  
  TVector3 r0_uvw = compute_focalPlane_intercept(intercept); 
  
  Double_t Sv = (intercept[2]-intercept[0])/wSep(); 
  Double_t Su = (intercept[3]-intercept[1])/wSep(); 
  
  TVector3  S_uvw( Su,Sv,1. );  //  = new TVector3( Su,Sv,1. ); 
  
  TVector3 r0_xyz = vec_uvw_xyz( r0_uvw ); 
  TVector3 S_xyz  = vec_uvw_xyz( S_uvw  ); 
    
  Double_t tz = (s2_z() - r0_xyz.Z())/S_xyz.Z();  
  
  TVector3 rs_xyz( r0_xyz + S_xyz*tz );   
  
  if (axis==kXYZ) { //return the point in XYZ coords (default) 
    return rs_xyz;
  } else          { //return the point in UVW coords 
    return vec_xyz_uvw(rs_xyz); 
  }
}
Double_t compute_s2_xParam(Double_t intercept[4], int pad) 
{
  //a track's intercept with the S2, normalized to the width of a paddle. 
  
  TVector3 rs_xyz = compute_s2_intercept(intercept); 
  
  return (rs_xyz.X()-s2_x(pad))/s2_padWidth();
}
Double_t compute_ToF(int plane,   //0&2=U-planes, 1&3=V-planes
		     Double_t vi, //wire pos of ToF to evaluate
		     Double_t intercept[4])
{
  //compute flight time between S2 & and wirePos 'vi'
  
  TVector3 r0_uvw = compute_focalPlane_intercept(intercept); 
  TVector3 rs_uvw = vec_xyz_uvw( compute_s2_intercept(intercept) ); 
  
  Double_t Tf_total = (rs_uvw-r0_uvw).Mag()/3e8; 
  
  Double_t ToF; 
  if (plane==0 || plane==2) { //if this is the U1 or U2 plane
    ToF = ( rs_uvw.Y()-vi )/( rs_uvw.Y()-r0_uvw.Y() ); 
  } else         {
    ToF = ( rs_uvw.X()-vi )/( rs_uvw.X()-r0_uvw.X() ); 
  }
  
  return ToF*Tf_total; 
}
void draw_track_times(TCanvas *c, 
		      vector<TVector2> points[4], 
		      Double_t intercept[4], 
		      Double_t T_s2)
{
  TH2D *hist2d[4]; 
  
  TString planeName[4] = {"U1", "V1", "U2", "V2"}; 
    
  c = new TCanvas; 
  
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0);
  
  c->Divide( 2,2, 0,0 ); 
  
  int iCanvas=0;
  
  
  Double_t mv = wSep()/(intercept[2]-intercept[0]); 
  Double_t mu = wSep()/(intercept[3]-intercept[1]);   
  
  
  
  for (int p=0; p<4; p++) {
    
    Double_t m, dtdv; bool is_Uplane; 
    if (p==0 || p==2) { 
      m = mv; is_Uplane=true;  
    } 
    else              { 
      m = mu; is_Uplane=false; 
    }
    
    hist2d[p] = new TH2D("hist2d_"+TString::Itoa(p,10), "", 
			 450, -22, 22,
			 450, -50, 500); 
    
    for (int h=0; h<points[p].size(); h++) { 
      
      TVector2 pt = points[p].at(h); 
      
      Double_t vi    = pt.X();
      Double_t T_tdc = pt.Y(); 
      
      hist2d[p]->Fill( (vi-intercept[p])*1e3, T_tdc*1e9 ); 
    }
    
    
    iCanvas++; 
    c->cd(iCanvas); 
      
    hist2d[p]->SetMarkerStyle(kOpenCircle); 
    //hist2d[p]->GetYaxis()->SetNdivisions(1); 
    hist2d[p]->SetMarkerSize(0.45); 
    hist2d[p]->DrawCopy();
    
    //draw drift-time prediction lines
    
    Double_t T_flt[5] 
      = { intercept[p] - 18e-3/m,
	  intercept[p] - 3e-3/m,
	  intercept[p],
	  intercept[p] + 3e-3/m,
	  intercept[p] + 18e-3/m }; 
    
    for (int j=0; j<5; j++) { 
      T_flt[j] = compute_ToF(is_Uplane, T_flt[j], intercept); 
    }
    
    
    
    //drift-lines
    TLine *mLine[4] 
      = { new TLine(-18/m, T_s2*1e9 + (20.32*18-25.56) -1e9*T_flt[0], 
		    -3/m,  T_s2*1e9 + (20.32*3 -25.56) -1e9*T_flt[1]),
	  
	  new TLine(-3/m,  T_s2*1e9 + (11.8*3)         -1e9*T_flt[1], 
		    0,     T_s2*1e9 + (0)              -1e9*T_flt[2]),
	  
	  new TLine(0,     T_s2*1e9 + (0)              -1e9*T_flt[2], 
		    3/m,   T_s2*1e9 + (11.8*3)         -1e9*T_flt[3]),
	  
	  new TLine(3/m,   T_s2*1e9 + (11.8*3)         -1e9*T_flt[3], 
		    18/m,  T_s2*1e9 + (20.32*18-25.56) -1e9*T_flt[4]) }; 
    
    for (int j=0; j<4; j++) mLine[j]->Draw("SAME"); 
    
    //Draw text which labels each cluster
    TText *tex = new TText;
    
    tex->DrawTextNDC(0.5, 0.925, planeName[p]); 
    

  }
  
  return; 
}
void draw_vLine( TH2 *hist, Double_t x, Int_t color=kBlack, Int_t style=1) 
{
  //WARNING!!! this only works for a th2. use draw_vLine_th1() for 1-d hist. 
  Double_t yMin = hist->GetYaxis()->GetXmin();
  Double_t yMax = hist->GetYaxis()->GetXmax();
  
  TLine *l = new TLine(x, yMin, x, yMax); 
  
  l->SetLineColor(color); 
  l->SetLineStyle(style); // 1 = solid line
  l->Draw("SAME"); 
  return; 
}	
void draw_vLine_th1(TH1 *hist, Double_t x, Int_t color=kBlack, Int_t style=1) 
{
  //does the same as above, but works for 1d-histos
  Double_t yMin = 0;
  Double_t yMax = hist->GetMaximum(); 
  
  TLine *l = new TLine(x, yMin, x, yMax); 
  
  l->SetLineColor(color); 
  l->SetLineStyle(style); // 1 = solid
  l->Draw("SAME"); 
  return; 
}	
void draw_hLine(TH1 *hist, Double_t y, Int_t color=kBlack, Int_t style=1) 
{
  Double_t xMin = hist->GetXaxis()->GetXmin();
  Double_t xMax = hist->GetXaxis()->GetXmax();
  
  TLine *l = new TLine(xMin, y, xMax, y); 
  
  l->SetLineColor(color); 
  l->SetLineStyle(style); 
  l->Draw("SAME"); 
  return; 
}	
void draw_box( double x0, 
	       double y0, 
	       double x1, 
	       double y1, int color=kBlack, int style=1 ) 
{
  //draws a box. 
  
  auto line = new TLine; 
  
  line->SetLineColor(color); 
  line->SetLineStyle(style); 

  line->DrawLine(x0,y0, x1,y0); 
  line->DrawLine(x1,y0, x1,y1); 
  line->DrawLine(x1,y1, x0,y1); 
  line->DrawLine(x0,y1, x0,y0); 
  
  return; 
}
void draw_eventPicture( TCanvas *canv, 
			int nHits[4],
			Double_t hit_rawTime[4][maxNum_hits],
			Double_t hit_wire[4][maxNum_hits],
			Double_t tR[16], //s2 rawTime (R pmt) 
			Double_t tL[16], //s2 rawTime (L pmt)
			Int_t nTrack,
			Double_t track_int[4][maxNum_newTrack],
			Double_t track_time[maxNum_newTrack],
			Double_t T_s2=-1e30,
			const Double_t vdc_minTime=-75e-9,
			const Double_t vdc_maxTime= 700e-9)
{
  
  bool draw_coincTime=false;  
  if (T_s2 > -1e29) { //then draw coinctime
    draw_coincTime=true;
  }
  
  
  cout << "Drawnig event picture." << endl; 
  cout << "TDC hits in each plane = {"
       << nHits[0]<<", "<<nHits[1]<<", "<<nHits[2]<<", "<<nHits[3]<<"}"<<endl; 
  cout << "Number of tracks to draw = " << nTrack << endl; 
  cout << "VDC plane TDC window drawn = { " 
       << vdc_minTime*1e9 << " ns,  " << vdc_maxTime*1e9 << " ns }" << endl; 
  cout << "Drawing dotted line over S2 coinc time? " 
       << ((draw_coincTime) ? "yes." : "no.") << endl; 
  
  
  TH2D *hist2d[4]; 
  
  TString planeName[4] = {"U1", "V1", "U2", "V2"}; 
  
  canv = new TCanvas("eventPicture", "Event Drawnig", 0,0, 1350, 500); 
    
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0);
  canv->Divide( 1,5, 0,0 ); 
    
  TLine *mLine; 
  
  int iCanvas=0;
  
  for (int p=0; p<4; p++) {
    
    
    hist2d[p] = new TH2D("vdcPlane_"+TString::Itoa(p,10), "", 
			 368, wirePos(p,367)-(dWire/2), wirePos(p,0)+(dWire/2),
			 450, -75,  vdc_maxTime*1e9); 
    
    for (int h=0; h<nHits[p]; h++) { 
      
      int wire = (int)hit_wire[p][h];
      
      Double_t T_tdc = vdc_realTime(p, wire, hit_rawTime[p][h]); 
      
      hist2d[p]->Fill( wirePos(p, wire), T_tdc*1e9 ); 
    }
    
    canv->cd(5-iCanvas); 
    iCanvas++; 
    
    hist2d[p]->SetMarkerStyle(kOpenCircle); 
    hist2d[p]->GetYaxis()->SetNdivisions(1); 
    hist2d[p]->SetMarkerSize(0.45); 
    hist2d[p]->DrawCopy();
    
    //Draw text which labels each cluster
    TText *tex = new TText;
    
    tex->DrawText(0.5, 0.925, planeName[p]); 
  }//for (int p=0; p<4; p++) { 
  
  for (int p=0; p<4; p++) { 
    canv->cd(5-p);
    hist2d[p]->DrawCopy(); 
  }
  
  //draw S2
  canv->cd(1); 
  TH2D *h_s2 = new TH2D("h_s2", "", 
			16, 
			s2_x(0)-(s2_padWidth()/2), 
			s2_x(15)+(s2_padWidth()/2), 
			35, -35, 35); 
  
  for (int p=0; p<16; p++) { 
    
    if (!s2_isCoinc( tR,tL, p )) continue; 
    
    h_s2->Fill( s2_x(p), s2_coincTime( tR,tL, p )*1e9 ); 
  }
  
  h_s2->GetXaxis()->SetNdivisions(15); 
  h_s2->GetYaxis()->SetNdivisions(1); 
  h_s2->Draw("box"); 
  
  //draw coinc time
  if (draw_coincTime) { 
    draw_hLine( h_s2, T_s2*1e9, kBlue, kDotted ); 
  }
  
  int oneTrack = -1; 
  if (nTrack < 0) { //user has asked to draw just one track
    
    oneTrack = -nTrack-1; 
    nTrack = maxNum_newTrack; 
  }
  
  //now do track 
  for (int t=0; t<nTrack; t++) { 
    
    if ( oneTrack >= 0 && t != oneTrack ) continue; //draw only 1 track?  
    
    Double_t T_track = track_time[t]; 
    
    Double_t tInt[4]; for (int p=0; p<4; p++) tInt[p] = track_int[p][t]; 
    
    iCanvas=0; 
    
    for (int p=0; p<4; p++) { 
      
      Double_t m; bool is_Uplane; 
      if (p==0 || p==2) { 
	m = wSep()/(tInt[2]-tInt[0]); 
	is_Uplane=true;  
      } 
      else              { 
	m = wSep()/(tInt[3]-tInt[1]);   
	is_Uplane=false; 
      }
      
      //draw tracks
      Double_t p_dw[5] 
	= { -16e-3, -3e-3, 0.0, 3e-3, 16e-3 }; 
      
      const int nPts=5; 
      
      Double_t p_dv[nPts]; 
      
      Double_t p_time[nPts];
      
      for (int j=0; j<nPts; j++) { 
	//cout << p_dv[j] << " " << p_time[j] << endl; 
	
	p_dv[j] = p_dw[j]/m + track_int[p][t]; 
	
	p_time[j] = 
	  T_track + drift_T(p_dw[j]) 
	  - compute_ToF(is_Uplane, p_dv[j], tInt); 
	//cout << "vi,T = " << p_dv[j] << ", " << (p_time[j]*1e9) << endl; 
      }
      
      canv->cd(5-iCanvas); 
      iCanvas++;
      
      for (int j=0; j<nPts-1; j++) { 
		
	mLine = new TLine(p_dv[j],   p_time[j]*1e9, 
			  p_dv[j+1], p_time[j+1]*1e9); 
	mLine->SetLineColor(kRed); 
	mLine->Draw("SAME"); 
	//mLine->~TLine(); 
      }
    }//for (int p=0; p<4; p++) { 
    
    //draw S2
    canv->cd(1); 
    
    TVector3 s2_int = compute_s2_intercept( tInt ); 
    
    mLine = new TLine(s2_int.X(), T_track*1e9-20, 
		      s2_int.X(), T_track*1e9+20);  
    mLine->SetLineColor(kRed); mLine->Draw("SAME"); 
    
    mLine = new TLine(s2_int.X()-0.15, T_track*1e9, 
		      s2_int.X()+0.15, T_track*1e9);  
    mLine->SetLineColor(kRed); mLine->Draw("SAME"); 
    
  }//for int p=0; p<4; p++) { 
  
  return; 
}
void draw_ePic_oneTrack( TCanvas *canv, 
			 int nHits[4],
			 Double_t hit_rawTime[4][maxNum_hits],
			 Double_t hit_wire[4][maxNum_hits],
			 Double_t tR[16], //s2 rawTime (R pmt) 
			 Double_t tL[16], //s2 rawTime (L pmt)
			 Double_t track_index,
			 Double_t track_int[4][maxNum_newTrack],
			 Double_t track_time[maxNum_newTrack],
			 Double_t T_s2=-1e30,
			 const Double_t vdc_minTime=-75e-9,
			 const Double_t vdc_maxTime= 700e-9 ) 
{
  //execute eventPicture just as above, but only run 1 track. 
  draw_eventPicture( canv, nHits, hit_rawTime, hit_wire, tR, tL, 
		     -track_index-1,
		     track_int, track_time, T_s2, vdc_minTime, vdc_maxTime );
    
  return;
}
	
void draw_track_TDC( TCanvas *canv, 
		     int nHits[4],
		     Double_t hit_rawTime[4][maxNum_hits],
		     Double_t hit_wire[4][maxNum_hits],
		     Double_t t_int[4],
		     Double_t t_time,
		     double Dt,
		     double xParam,
		     const Double_t vdc_minTime= -15e-9,
		     const Double_t vdc_maxTime=  300e-9, int t_index=-1)
{
  
  bool draw_coincTime=false;  
    
  cout << "Drawnig event picture." << endl; 
  cout << "TDC hits in each plane = {"
       << nHits[0]<<", "<<nHits[1]<<", "<<nHits[2]<<", "<<nHits[3]<<"}"<<endl; 
  cout << "VDC plane TDC window drawn = { " 
       << vdc_minTime*1e9 << " ns,  " << vdc_maxTime*1e9 << " ns }" << endl; 
  
  
  TH2D *hist2d[4]; 
  
  TString planeName[4] = {"U1", "V1", "U2", "V2"}; 
  
  canv = new TCanvas; 
    
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0);
  canv->Divide( 2,2, 0,0 ); 
  
  TString 
    title = "Run: "+TString::Itoa(RUN_runNumber,10);
  
  title += ((checkArm()) ? "-R" : "-L"); 
    
  if (t_index > -1) {
    title += "     track index: " + TString::Itoa(t_index,10); 
  }
  
  canv->SetTitle(title); 
      
  TLine *mLine; 
  
  int iCanvas=0;
  
  const float max_dv = 20e-3;
  
  cout << "track time = " << t_time << endl; 
  
  
  double slope;
    
  const double max_dw = 18e-3; 
  
  const int nFitPts =25; 
  
  TGraph *gFit[4];   
   
  cout << "intercept = {" 
       << t_int[0] << ", " 
       << t_int[1] << ", " 
       << t_int[2] << ", " 
       << t_int[3] << "}" << endl; 
  
  for (int p=0; p<4; p++) {
    
    double plane_time = t_time - compute_ToF(p, t_int[p], t_int); 
    
    cout << "plane " << p << " time = " << plane_time << endl; 
    
    hist2d[p] = new TH2D("vdcPlane_"+TString::Itoa(p,10), "", 
			 200, t_int[p]-max_dv, t_int[p]+max_dv, 
			 450, 
			 (plane_time+vdc_minTime)*1e9,  
			 (plane_time+vdc_maxTime)*1e9 ); 
        
    
    for (int h=0; h<nHits[p]; h++) { 
      
      int wire = (int)hit_wire[p][h];
      
      Double_t T_tdc = vdc_realTime(p, wire, hit_rawTime[p][h]); 
      
      hist2d[p]->Fill( wirePos(p, wire), T_tdc*1e9 ); 
    }
    
    canv->cd(p+1); 
    
    hist2d[p]->SetMarkerStyle(kOpenCircle); 
    hist2d[p]->GetYaxis()->SetNdivisions(1); 
    hist2d[p]->SetMarkerSize(0.45); 
    hist2d[p]->DrawCopy();
    
    //Draw text which labels each cluster
    TLatex *tex = new TLatex;
    
    tex->DrawLatex(t_int[p]-max_dv*0.90, (plane_time + vdc_maxTime*0.9)*1e9, 
		   planeName[p]); 
    
    //now, draw the fitted track
    if (p==0 || p==2) { slope = wSep()/(t_int[2]-t_int[0]); }
    else              { slope = wSep()/(t_int[3]-t_int[1]); }
    
    cout << "slope " << slope << endl; 
    
    iCanvas++;

    double p_v[nFitPts], p_Tau[nFitPts]; 
    for (int h=0; h<nFitPts; h++) { 
      
      double dw = max_dw*( 1.0 - 2.0*(((double)h)/((double)nFitPts-1)) ); 
      
      p_v[h]   = dw/slope + t_int[p]; 
      
      p_Tau[h] = drift_T(dw) + (t_time-compute_ToF(p, p_v[h], t_int)); 
      
      p_Tau[h] *= 1e9; //convert to ns
      
      //cout << p << " | " << p_v[h] << "  " << p_Tau[h] << endl; 
      
    }//for (int h=0; h<nFitPts; h++) { 
    
    gFit[p] = new TGraph(nFitPts, p_v, p_Tau); 
    
    gFit[p]->SetLineColor(kRed); 
    gFit[p]->Draw("SAME"); 
    
  }//for (int p=0; p<4; p++) { 
  
  
   
  return; 
}    
TVector2 draw_gausFit( TH1 *hist, 
		       Double_t radius, 
		       Double_t center=-1e30, 
		       Int_t color=kRed) 
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

  return TVector2( x0, sig ); 
}
TVector2 draw_gausFit_fixedBase( TH1 *hist, 
				 double radius, 
				 double base=0,
				 double center=-1e30, 
				 Int_t color=kRed) 
{
  //center the fit aruond the max. value of the histogram
  if (center < -1e29) { 
    center = hist->GetBinCenter( hist->GetMaximumBin() ); 
  }
  
  TF1 *gausFit = new TF1( "gausFit", "gaus(0) + [3]", 
			  -radius+center, 
			   radius+center ); 
  
  TAxis *xAxis = hist->GetXaxis(); 
  
  Double_t A = hist->GetMaximum() - base;
  
  if (color >= 0)
    cout << "A,base = " << A << ", " << base << endl; 
  
  gausFit->SetParameter( 0, A ); 
  gausFit->SetParameter( 1, center ); 
  gausFit->SetParameter( 2, radius/2 ); 
  gausFit->FixParameter( 3, base ); 
  
  TFitResultPtr fitPtr = hist->Fit("gausFit", "L N R S Q B"); 
  
  double x0     = fitPtr->Parameter(1); 
  double x0_err = fitPtr->ParError(1); 
  
  double sig     = fitPtr->Parameter(2); 
  double sig_err = fitPtr->ParError(2); 
  
  A = fitPtr->Parameter(0); 
  
  double yMax = hist->GetMaximum(); 
  double xMin = xAxis->GetBinCenter(1); 
  
  if (color >= 0) {
    gausFit->SetLineColor(color); 
    gausFit->Draw("SAME");
    
    TF1 *gaus_copy = new TF1( (*gausFit) ); //copy the fit, to draw another one
    
    gaus_copy->SetRange( xAxis->GetXmin(), xAxis->GetXmax() ); 
    
    gaus_copy->SetLineStyle(kDotted);
    gaus_copy->Draw("SAME"); 
  
    cout << "sigma = " << sig << " +/- " << sig_err << endl; 
    cout << "mean  = " << x0  << " +/- " << x0_err  << endl; 
    
    cout << "A,base = " << A << ", " << base << endl; 
    
    TLatex *ltx = new TLatex; 
  
    ltx->DrawLatex(xMin*0.90, yMax, 
		   TString::Format("#bar{x} = %0.3e",x0)+
		   TString::Format(" #pm %0.1e",x0_err) ); 
    
    ltx->DrawLatex(xMin*0.90, yMax*0.92, 
		   TString::Format("#sigma = %0.3e",sig)+
		   TString::Format(" #pm %0.1e",sig_err) ); 
  }
  return TVector2(x0,sig); 
}
double getMedian( TH1D *h )
{
  //find the median of the histogram (copied this code from the root forum)
  const int n = h->GetXaxis()->GetNbins(); 
  double x[n]; 
  double y[n]; 
  h->GetXaxis()->GetCenter(x); 
  
  for (int i=1; i<=n; i++) y[i] = h->GetBinContent(i); 
  
  return TMath::Median(n, x, y); 
}
double StdDev_fixedMean( TH1D *h, double mean=0. ) 
{
  const int n = h->GetXaxis()->GetNbins(); 
  double x[n]; 
  
  double variance(0.); 
  h->GetXaxis()->GetCenter( x ); 
  
  for (int i=1; i<=n; i++) 
    variance += h->GetBinContent(i)*TMath::Power( x[i]-mean, 2 ); 
  
  return TMath::Sqrt( variance/h->Integral() ); 
}
void SetBranchAddress_VDC( TTree *tree, 
			   TString varName, 
			   Double_t (*b_address)[4][maxNum_hits] )
{ 
  
  TString Arm;
  if (checkArm()==kArm_Right) { Arm = "R"; }
  else                        { Arm = "L"; } 
  
  tree->SetBranchAddress( Arm+".vdc.u1."+varName, &(*(b_address))[0] );
  tree->SetBranchAddress( Arm+".vdc.v1."+varName, &(*(b_address))[1] );
  tree->SetBranchAddress( Arm+".vdc.u2."+varName, &(*(b_address))[2] );
  tree->SetBranchAddress( Arm+".vdc.v2."+varName, &(*(b_address))[3] );
  
  return; 
}
void SetBranchAddress_VDC_Ndata( TTree *tree, 
				 TString varName, 
				 Int_t (*b_address)[4] )
{ 
  TString Arm;
  if (checkArm()==kArm_Right) { Arm = "R"; }
  else                        { Arm = "L"; } 
  
  
  tree->SetBranchStatus( Arm+".vdc.*."+varName, 1 ); 
  
  tree->SetBranchAddress( "Ndata."+Arm+".vdc.u1."+varName, &(*(b_address))[0] );
  tree->SetBranchAddress( "Ndata."+Arm+".vdc.v1."+varName, &(*(b_address))[1] );
  tree->SetBranchAddress( "Ndata."+Arm+".vdc.u2."+varName, &(*(b_address))[2] );
  tree->SetBranchAddress( "Ndata."+Arm+".vdc.v2."+varName, &(*(b_address))[3] );
  
  return; 
}
bool Solve_symmetric_system( double *A_arr, double *b, const int N=5  ) { 
  
  //takes the NxN matrix A, and the N-vector B, and solves it. 
  
  //performs numerical LU factorization 
  // A must be an array of size N*(N+1)/2 
  //   (i.e., the number of unique elements in a symmetric, NxN matrix)
  
  //ASSUMES MATRIX IS NONSINGULAR!!!!
  
  const double det_cutoff = 1e-80; 
  
  const int n_unique_elem = N*(N+1)/2; 
  
  double U[N][N]; int elem(-1); 
  
  for (int i=0; i<N; i++) {  
    for (int j=i; j<N; j++) { elem++; 
      
      U[i][j] = A_arr[elem]; 
      U[j][i] = A_arr[elem]; 
    }
  }
  
  double det = 1.; 
  
  //create the L-matrix
  double L[N][N];  
  for (int i=0; i<N; i++) 
for (int j=0; j<N; j++) { L[i][j] =0.; }
  	
  for (int Lj=0; Lj<N; Lj++) { 
    
    double a_ii = U[Lj][Lj]; 
          
    for (int Li=Lj; Li<N; Li++) { 
            
      if (Li==Lj) { L[Li][Lj] =1.; continue; }
      
      double m_ij = U[Li][Lj] / a_ii; 
      
      for (int Uj=Lj+1; Uj<N; Uj++) { 
	
	      U[Li][Uj] += -m_ij * U[Lj][Uj];
      } 
      
      L[Li][Lj] = m_ij; 
      
      
    }
    det *= U[Lj][Lj]; 
  }
  //get diagonal elements
  
  //check if matrix is invertable
  cout << "Det(A) = " << det << endl; 
  if (TMath::Abs(det) < det_cutoff || det != det ) return false; 
  
  double y[N]; 
  
  //solve the system Ly = b
  for (int i=0; i<N; i++) { 
    
    y[i] = b[i]; 
    
    for (int j=0; j<i; j++) { 
      
      y[i]  +=  -L[i][j] * y[j]; 
    }
  }
  
  for (int i=N-1; i>=0; i+=-1) { 
    
    b[i] = y[i]; 
    
    for (int j=N-1; j>i; j+=-1) { 
      
      b[i]  +=  -U[i][j] * b[j]; 
    }
    b[i] *= 1./U[i][i]; 
  }
  
  return true; 
}
double compute_minVariance( TH1 *h_input ) 
{
  TH1 *h = (TH1*)h_input->Clone(); 
  
  TAxis *axis = h->GetXaxis(); 
  
  h->Scale( ((double)axis->GetNbins()) / h->Integral() ); 
  
  const int N_Bins = axis->GetNbins();   double N_H = (double)N_Bins; 
  
  //normalize the histogram ( scale it to represent dP/dz )
  double dx = (axis->GetXmax() - axis->GetXmin()) / (N_H-1.); 
  
  vector<double> phi(N_Bins,0.); 
  
  double Phi_i(0.); 
  
  for (int n=0; n<N_Bins; n++) {
    
    phi[n] = Phi_i;
    
    double dPhi_dx = (1.0/h->GetBinContent(n+1)) - 1;  
    
    Phi_i += dPhi_dx*dx; 
  }

  double X[N_Bins]; axis->GetCenter(X); 
  
  double slope_offset
    = (phi[N_Bins-1]-phi[0]) /  (axis->GetXmax()-axis->GetXmin()); 
  
  double square_error(0.); 
  
  for (int n=0; n<N_Bins; n++) {
    
    square_error += pow( phi[n] - slope_offset*(X[n]-X[0]), 2 ); 
    //phi[n] += -phi_avg;
    phi[n] += -slope_offset*(X[n]-X[0]);
    phi[n] += X[n]; 
  }

  /*
  new TCanvas; 
  auto g = new TGraph(N_Bins, X, phi.data()); 
  //g->GetYaxis()->SetRangeUser( -1, 1 ); 
  g->Draw(); 
  
  auto func = new TF1("func", "x", axis->GetXmin(), axis->GetXmax());
  func->SetLineStyle(kDotted); 
  func->Draw("SAME"); //*/ 
  
  return square_error / N_H;
} 
vector<double> TVectorD_to_vector( const TVectorD vec ) {
  //create a vector from this array
  auto arr = vec.GetMatrixArray();
  
  const int size = vec.GetNoElements();
  
  vector<double> myVec( arr, arr+size );
  
  return myVec;
}
int function_lib(int runNum)
{
  
  
  RUN_runNumber = runNum; 
  
  //is this run a VDC-voltage-affected run?
  RUN_is_VDC_affected = (RUN_runNumber < first_goodRun); 
  
  
  //debug - for drift-model generation
  if (RUN_runNumber == 4208) { 
    min_slope = -0.42; 
    max_slope =  0.42; 
    
    min_tau = -15e-9; 
    max_tau = 380e-9;
    
    cutLine = 0; //-33.803e-9; 
  }
  if (RUN_runNumber == 3930) { 
    min_slope = -0.40; 
    max_slope =  0.30; 
    
    min_tau = -10e-9; 
    max_tau = 320e-9;
    
    cutLine = 0; //-17e-9/0.55; 
  }
  if (RUN_runNumber == 4329) { 
    min_slope = -0.40; 
    max_slope =  0.30; 
    
    min_tau = -10e-9; 
    max_tau = 310e-9; 
    
    cutLine = 0; //-33.962e-9;
  }


  
  //find the beam current associated with this run
  fstream datFile;
    
  datFile.open("affected_VDC_runs.dat", ios::in); 
    
  TString iStr(" ");
  
  bool found_data=false; 
  
  while (iStr != "") { 
      
    datFile >> iStr;
    Int_t    iRun     = atoi( iStr.Data() ); 
      
    datFile >> iStr; 
    Double_t iCurrent = atof( iStr.Data() );
            
    datFile >> iStr; 
    Double_t iT_coinc = atof( iStr.Data() ); 
      
    datFile >> iStr; 
    Double_t iT_coincSigma = atof( iStr.Data() );  
      
    if (iRun == runNum) {
      found_data=true; 
      RUN_beamCurrent  = iCurrent;
      RUN_Tcoinc       = iT_coinc/1e9;
      RUN_Tcoinc_sigma = iT_coincSigma/1e9;
	
      break; 
    }
  }
  datFile.close(); 
  
  if (!found_data) { 
    cout << "FATAL ERROR !!! run number '" << runNum << "' not found!!!!" << endl;  
    return -1; //return fatal error 
  }
    
  //print some data about this run
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; 
  
  cout << "processing run " << RUN_runNumber << endl;
  cout << "RUN_beamCurrent=("<< RUN_beamCurrent <<" uA)" << endl
       << "RUN_Tcoinc(R-L)=("<< RUN_Tcoinc*1e9 <<" ns)" << endl 
       << "RUN_Tcoinc_sigma=("<< RUN_Tcoinc_sigma*1e9 <<" ns)" << endl
       << "is vdc-affected? = "<< ((RUN_is_VDC_affected)?"yes":"no") << endl; 
  
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; 
      
  
  if (RUN_is_VDC_affected) {
    
    //set new drift parameters
    if (checkArm()==kArm_Right) { 
      
      std::copy( &p_lowT_right[0][0], &p_lowT_right[0][0]+5*5, &param_lowT[0][0] ); 
      std::copy( &p_lowX_right[0][0], &p_lowX_right[0][0]+5*5, &param_lowX[0][0] ); 
      
      
    } else { 
      
      if (RUN_runNumber==4329) {
	
	std::copy( &p_lowT_4329[0][0], &p_lowT_4329[0][0]+5*5, &param_lowT[0][0] ); 
	std::copy( &p_lowX_4329[0][0], &p_lowX_4329[0][0]+5*5, &param_lowX[0][0] ); 
 	
	k0 =1.030e-3; 
	kT =46931.9; 
	kM =1.060e-4; 
	
	cout << "run number 4329 params used" << endl; 
      }
      
      if (RUN_runNumber==3926) {
	
	std::copy( &p_lowT_3926[0][0], &p_lowT_3926[0][0]+5*5, &param_lowT[0][0] ); 
	std::copy( &p_lowX_3926[0][0], &p_lowX_3926[0][0]+5*5, &param_lowX[0][0] ); 
	
	k0 =0.885e-3; 
	kT =47433.8;
	kM =3.344e-4; 
      }
      
    }
    
    //interpolate based on beam current
    double I_lo(33.2),     I_hi(52.8); //beam current in uA 
      
    double V_lo(4.738e4),  V_hi(4.693e4);  
    double B_lo(1.00e-3),  B_hi(0.91e-3); 
      
    double x_lo(3.56e-3),  x_hi(3.67e-3);
    double t_lo(54.0e-9),  t_hi(58.7e-9); 
    
    V_drift_L  = ( (V_hi-V_lo) / (I_hi-I_lo) )*( RUN_beamCurrent - I_lo ) + V_lo; 
    
    B_drift_L  = ( (B_hi-B_lo) / (I_hi-I_lo) )*( RUN_beamCurrent - I_lo ) + B_lo; 
      
    x0_drift_L = ( (x_hi-x_lo) / (I_hi-I_lo) )*( RUN_beamCurrent - I_lo ) + x_lo; 
      
    t0_drift_L = ( (t_hi-t_lo) / (I_hi-I_lo) )*( RUN_beamCurrent - I_lo ) + t_lo; 
    
  } 
  
  return 1; //return 'good'
}
