
#include <TapexEventHandler.h>
#include <cmath> 
#include <TMath.h> 

namespace {
    ////////////////////////////////////////////////////////////////////////////////  
    double fDrift_param_L_BC_Lo[5][5] 
    = { {3.562e-3, 3.517e-3, 3.587e-3, 3.403e-3, 3.342e-3}, 
        {4.287e-8, 4.200e-8, 4.087e-8, 3.834e-8, 3.554e-8}, 
        {2.114e-5, 2.095e-5, 2.077e-5, 2.036e-5, 2.033e-5}, 
        {4.735e-3, 4.652e-3, 4.547e-3, 4.271e-3, 3.654e-3}, 
        {8.510e-9, 8.510e-9, 8.510e-9, 8.510e-9, 8.510e-9} }; 

    /*= { {2.718e-3, 2.658e-3, 2.530e-3, 2.299e-3, 2.349e-3}, 
        {3.106e-8, 2.968e-8, 2.844e-8, 2.879e-8, 2.950e-8},
        {2.201e-5, 2.152e-5, 2.098e-5, 2.038e-5, 2.046e-5},
        {5.042e-3, 4.709e-3, 3.965e-3, 4.430e-3, 4.501e-3}, 
        {6.603e-9, 4.276e-9, 1.888e-9, 1.872e-9, 3.195e-9} }; */ 

    double fDrift_param_L_BC_Hi[5][5] 
    = { {3.173e-3, 3.280e-3, 3.347e-3, 3.500e-3, 3.621e-3}, 
        {4.218e-8, 4.205e-8, 4.156e-8, 4.164e-8, 4.128e-8},
        {2.105e-5, 2.101e-5, 2.087e-5, 2.091e-5, 2.083e-5}, 
        {4.514e-3, 4.546e-3, 4.543e-3, 4.584e-3, 4.605e-3}, 
        {1.325e-8, 1.550e-8, 1.800e-8, 1.950e-8, 2.150e-8} }; 

    //    {5.750e-9, 8.000e-9, 1.050e-8, 1.200e-8, 1.400e-8} }; */ 


    double fParams_R[5][5] 
    = { {3.774e-3, 3.774e-3, 3.774e-3, 3.774e-3, 3.774e-3}, 
        {3.618e-8, 3.618e-8, 3.618e-8, 3.618e-8, 3.618e-8},
        {2.024e-5, 2.024e-5, 2.024e-5, 2.024e-5, 2.024e-5},
        {3.864e-3, 3.864e-3, 3.864e-3, 3.864e-3, 3.864e-3},
        {1.431e-8, 1.431e-8, 1.431e-8, 1.431e-8, 1.431e-8} }; 

    //first run which is 'fixed' 
    const int FIRST_FIXED_RUN = 4619; 

    const int N_PARAMS =5;
    const int Npts_SLOPE =5; 
    
    double slopeNodes_R[5] = { 1.144, 1.363, 1.473, 1.582, 1.868 };
    double slopeNodes_L[5] = { 1.155, 1.336, 1.419, 1.509, 1.824 }; 
    
    double f_BEAM_CURRENT[2] = {17., 53.};

    double fDrift_param_R[5][5];

    const int blur_nSamples = 25; 
    const double blur_sigma =0.20e-3; 

    const double beamCurrent_DEFAULT_VALUE = 30.; 

    const double zQuantSpread[25] = { -2.053749, 
                    -1.554774, 
                    -1.281552, 
                    -1.080319, 
                    -0.915365, 
                    -0.772193, 
                    -0.643345, 
                    -0.524401, 
                    -0.412463, 
                    -0.305481, 
                    -0.201893, 
                    -0.100434, 
                    0.000000, 
                    0.100434, 
                    0.201893, 
                    0.305481, 
                    0.412463, 
                    0.524401, 
                    0.643345, 
                    0.772193, 
                    0.915365, 
                    1.080319, 
                    1.281552, 
                    1.554774, 
                    2.053749  }; 


    double z_spread[25]; 
    double dP_spread; 
}

TapexEventHandler::TapexEventHandler( bool arm, 
			      double beamCurrent,
			      int runNumber, 
			      const TapexS2Hit *fHit_R, 
			      const TapexS2Hit *fHit_L ) { 
  f_activeArm  =arm; 
  
  fRunNumber = runNumber; 
  
  fBeamCurrent = beamCurrent; 
  
  
  if (fBeamCurrent < -1e3) { 
    f_isNullBeamCurrent=true; 
    fBeamCurrent = beamCurrent_DEFAULT_VALUE; 
  }
  
  fS2Hit_Right = fHit_R; 
  fS2Hit_Left  = fHit_L;  

  //develop drift-parameters for left arm
  for (int s=0; s<Npts_SLOPE; s++) { 
    for (int p=0; p<N_PARAMS; p++) { 
      
      if (fRunNumber < FIRST_FIXED_RUN) {
	double f_PARAM_INTERPOLATE[2] = { fDrift_param_L_BC_Lo[p][s], 
					  fDrift_param_L_BC_Hi[p][s] }; 
	
	fParams_L[p][s] = Interpolate( beamCurrent, 
				       f_BEAM_CURRENT, 
				       f_PARAM_INTERPOLATE, 2 ); 
	
      } else { fParams_L[p][s]=fParams_R[p][s]; }
    }
  }
  //this is to blur the model, to get rid of 'sharp points
  dP_spread = 1./((double)blur_nSamples);  double cum_P(dP_spread*0.50); 
  
  for (int i=0; i<blur_nSamples; i++) { 
    
    z_spread[i] = blur_sigma*zQuantSpread[i]; 
    cum_P += dP_spread;
  }
  
} 
double f_SIGMA_INTERPOLATE[2] = { 5e-9, 5e-9 }; 
double f_SIGMA_FIXED = 5e-9; 

double TapexEventHandler::Get_Dt() const
{
  return GetS2Hit(true)->Time() - GetS2Hit(false)->Time(); 
}

double TapexEventHandler::Get_tauSigma() const { 
  
  if (fRunNumber < FIRST_FIXED_RUN && ActiveArm()==false ) { 
    
    return Interpolate( fBeamCurrent, f_BEAM_CURRENT, f_SIGMA_INTERPOLATE, 2 ); 
  } 
  return f_SIGMA_FIXED; 
}
double TapexEventHandler::Drift_X( double tau, double slope, int derivative ) const { 
  
  double par[N_PARAMS];
  
  if (ActiveArm()==true) {
    
    for (int p=0; p<N_PARAMS; p++) 
      par[p] = Interpolate( slope, slopeNodes_R, fParams_R[p], Npts_SLOPE ); 
    
  } else                 { 
    
    for (int p=0; p<N_PARAMS; p++) 
      par[p] = Interpolate( slope, slopeNodes_L, fParams_L[p], Npts_SLOPE ); 
  }
  
  double x1(par[0]), y1(par[1]), m(par[2]), b(par[3]), y0(par[4]); 
  
  
  if (tau<y0) return 0; 
  
  if (tau<y0+y1) 
    return 0.5*( -b + sqrt( b*b + 4.*(tau-y0)*(x1*x1 + b*x1)/y1 )  ); 
  
  return ( tau - (y0+y1) )/m + x1; 
}
double TapexEventHandler::Drift_T_raw( const double x, 
				   const double *par, 
				   const int derivative ) const {
  
  double x1(par[0]), y1(par[1]), m(par[2]), b(par[3]), y0(par[4]); 
  
  double absX = std::fabs(x); 
  
  if (derivative==0) 
    return absX < x1  
      ? y1*(x*x + b*absX)/(x1*x1 + b*x1)  +  y0 
      : m*(absX - x1)  +  y1  +  y0; 
  
  if (derivative==1) 
    return absX < x1 
      ? y1*(2.*x + TMath::Sign(b,x) )/(x1*x1 + b*x1) 
      : TMath::Sign(m,x); 
  
  //derivative==2
  return absX < x1  ?  y1*( 2. )/(x1*x1 + b*x1)  :  0.; 
}
double TapexEventHandler::Drift_T( double x, double slope, int derivative ) const { 
  
  double par[N_PARAMS];
  
  if (ActiveArm()==true) {
    
    for (int p=0; p<N_PARAMS; p++) 
      par[p] = Interpolate( slope, slopeNodes_R, fParams_R[p], Npts_SLOPE ); 
    
  } else                 {    
    
    for (int p=0; p<N_PARAMS; p++) 
      par[p] = Interpolate( slope, slopeNodes_L, fParams_L[p], Npts_SLOPE ); 
  }

  //blur higher derivatives
  const double h=blur_sigma*0.6;

  if (std::fabs(x) < par[0]+blur_sigma*5.) {
        
    double Tau(0.); 
        
    if (derivative >0) {

      for (int i=0; i<blur_nSamples; i++) { 
      	Tau 
	  += Drift_T_raw( x+h+z_spread[i], par, derivative-1 )
	  -  Drift_T_raw( x-h+z_spread[i], par, derivative-1 ); 
      }
      return Tau*dP_spread/(2.*h); 
      
    } else             {
      
      for (int i=0; i<blur_nSamples; i++) 
	Tau 
	  += Drift_T_raw( x+z_spread[i], par ); 
      
      return Tau*dP_spread; 
    }
  } else  { 
    
    return Drift_T_raw( x, par, derivative );
  }    
  
}
double TapexEventHandler::Interpolate( const double x, 
				   const double *X, 
				   const double *Y, 
				   const int nPts   ) const {

  //interpolate between the points in X & Y. 'nPts' should be the n. of points
  //                                          in BOTH the X & Y arrays.
  // if 'x' is out of the bounds of X, then it just returns the value of the 
  // closest extreme value of Y. 
  int zone(0); 
  
  //find which 'zone' we're in
  while ( zone < nPts ) { 
    
    if (x < X[zone]) break; zone++; 
  }
  
  //we're in zone=0, which means we're out-of-bounds, to the left. 
  if (zone==0)    return Y[0];
  
  //we're in zone=nPts, which means we're out-of-bounds, to the right. 
  if (zone==nPts) return Y[nPts-1]; 
  
  //if were in a valid 'zone', then interpolate
  return Y[zone-1] + (Y[zone]-Y[zone-1])*( (x - X[zone-1])/(X[zone]-X[zone-1]) ); 
 
}
