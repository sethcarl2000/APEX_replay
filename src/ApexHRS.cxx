
#include <TROOT.h>
#include "TMath.h"
#include "ApexHRS.h"
#include <fstream>



#define DEBUG false

//some hard-coded APEX values, that do not change for any run-number
const int VDC_nPlanes = 4; 
const int VDC_nWires  = 368; 

//w-value of each plane (equivalent to det-coords z-value)
const double VDC_w_R[VDC_nPlanes] = { 0, 0.026, 0.333245, 0.359245 };
const double VDC_w_L[VDC_nPlanes] = { 0, 0.026, 0.335344, 0.361444 };

//positions of first-wires
const double VDC_wire0_R[VDC_nPlanes] = { 0.77852, 0.77852, 1.02793, 1.02793 }; 
const double VDC_wire0_L[VDC_nPlanes] = { 0.77852, 0.77852, 1.02718, 1.02718 }; 

double THRS::VDC::offset_RHRS[THRS::VDC::fN_PLANES][THRS::VDC::fN_WIRES];
double THRS::VDC::offset_LHRS[THRS::VDC::fN_PLANES][THRS::VDC::fN_WIRES];

THRS::VDC::InitStatus THRS::VDC::status_RHRS[4] = { THRS::VDC::NotInit }; 
THRS::VDC::InitStatus THRS::VDC::status_LHRS[4] = { THRS::VDC::NotInit }; 


map<THRS::MatrixType,THRS::TMatrix*> THRS::fMatrixMap_R{ 
  {THRS::Mat_fp_y,     new THRS::TMatrix(THRS::Mat_fp_y)},
  {THRS::Mat_tan_rho,  new THRS::TMatrix(THRS::Mat_tan_rho)},
  {THRS::Mat_fp_phi,   new THRS::TMatrix(THRS::Mat_fp_phi)},
  {THRS::Mat_tg_y,     new THRS::TMatrix(THRS::Mat_tg_y)},
  {THRS::Mat_tg_theta, new THRS::TMatrix(THRS::Mat_tg_theta)},
  {THRS::Mat_tg_phi,   new THRS::TMatrix(THRS::Mat_tg_phi)},
  {THRS::Mat_tg_dp,    new THRS::TMatrix(THRS::Mat_tg_dp)} }; 
    
map<THRS::MatrixType,THRS::TMatrix*> THRS::fMatrixMap_L{ 
  {THRS::Mat_fp_y,     new THRS::TMatrix(THRS::Mat_fp_y)},
  {THRS::Mat_tan_rho,  new THRS::TMatrix(THRS::Mat_tan_rho)},
  {THRS::Mat_fp_phi,   new THRS::TMatrix(THRS::Mat_fp_phi)},
  {THRS::Mat_tg_y,     new THRS::TMatrix(THRS::Mat_tg_y)},
  {THRS::Mat_tg_theta, new THRS::TMatrix(THRS::Mat_tg_theta)},
  {THRS::Mat_tg_phi,   new THRS::TMatrix(THRS::Mat_tg_phi)},
  {THRS::Mat_tg_dp,    new THRS::TMatrix(THRS::Mat_tg_dp)} }; 

THRS::HRS_initStatus THRS::OpticsStatus_RHRS = THRS::Init_notDone;
THRS::HRS_initStatus THRS::OpticsStatus_LHRS = THRS::Init_notDone;


////////////////////////////////////////////////////////////////////////////////////
void THRS::Parse_opticsData( const bool arm, const TString path_DB ) { 
  
  map<TString,MatrixType> typeNames{ 
    {"t",Mat_tan_rho}, 
    {"y",Mat_fp_y}, 
    {"p",Mat_fp_phi}, 
    {"D",Mat_tg_dp}, 
    {"T",Mat_tg_theta}, 
    {"P",Mat_tg_phi}, 
    {"Y",Mat_tg_y}     }; 
     

  const int N_exponent =3; 
  const int N_poly     =7; 
  
  const double minVal = 1e-30; 
  
  const TString here = "THRS::Parse_opticsData(TString)"; 
    
  cout << TString::Format("%s => parsing data from \"%s\"...",
			  here.Data(),
			  path_DB.Data()) 
       << flush; 
  
  fstream db; db.open( path_DB.Data(), ios::in );   
  
  //failed to open file! 
  if (!db) { 
    THRS::ReportError(here,"File could not be opened!"); 
    return; 
  }
    
  TString str; 
    
  while ( !db.eof() ) { 
    
    //let's find which element we're working with..
    db >> str; 
    
    //we've reached the end of the file
    if (str=="eof") { 
      if (arm==THRS::Arm_Right) { THRS::OpticsStatus_RHRS = Init_ok; }
      else                      { THRS::OpticsStatus_LHRS = Init_ok; }
      break; 
    }
    
    const TString name = str; 
    
    //check to see if the one-letter 'name' is actually a valid matrix name    
    auto findName = typeNames.find(name); 
    
    if (findName == typeNames.end()) {  
      
      if (arm==THRS::Arm_Right) { THRS::OpticsStatus_RHRS = Init_error; }
      else                      { THRS::OpticsStatus_LHRS = Init_error; }
      
      THRS::ReportError(here,"Invalid type read from optics file"); 
      return; 
    }
    
    MatrixType type = findName->second; 
    
    
    auto elem = new TMatrixElement(type); 
    
    for (int p=0; p<N_exponent; p++) { db >> str; elem->fExp[p] = atoi(str); }
    
    //now, get the polynomials
    for (int p=0; p<N_poly; p++) { 
      
      db >> str; double val = stod((string)str); 
      
      elem->fPoly.push_back(val); 
    }
    
    //check for 0-values on the end of the polynomial
    for (int p=elem->fPoly.size()-1; p>=0; p--) { 
      
      double val = elem->fPoly.at(p);  
      
      if (TMath::Abs(val) < minVal) { elem->fPoly.erase( elem->fPoly.begin()+p ); }
      else { break; }
    }
    
    
    //this extra 'int' on the end of each line, idk what it's for
    //UPDATE: I think this last power is powers in (fExp); 
    db >> str; elem->fExp[3] = atoi(str); 
        
    //check size of the polynomial
    if (elem->fPoly.size()==0) { 
      cout << here << " => zero poly. elems. found." << endl; 
      elem->~TMatrixElement(); 
      continue; 
    }			    
    
    auto matrix = (arm==THRS::Arm_Right) 
      ? THRS::fMatrixMap_R.find(type)->second 
      : THRS::fMatrixMap_L.find(type)->second; 
      
    matrix->fElems.push_back( elem ); 
    //now, add this element to its respective matrix
    
    //check which arm we're looking at
  }  
  
  //fill out what 'powers' we need for each computation of target variables
  
}
bool THRS::PrintStatus_optics() { 
  
  bool good=true; 
  
  map<MatrixType,TString> matrix_names{ 
    {Mat_tan_rho,  "tan_rho"},
    {Mat_fp_y,     "fp_y"},
    {Mat_fp_phi,   "fp_phi"}, 
    {Mat_tg_dp,    "tg_dp"}, 
    {Mat_tg_theta, "tg_theta"}, 
    {Mat_tg_phi,   "tg_phi"}, 
    {Mat_tg_y,     "tg_y"} }; 
  
  const TString here = "THRS::Print_matrixStatus()"; 
  
  //check all matricies that should exist
  vector<MatrixType> matrices{ 
    Mat_fp_y, 
    Mat_tan_rho, 
    Mat_fp_phi, 
    Mat_tg_y,
    Mat_tg_theta,
    Mat_tg_phi,
    Mat_tg_dp };  
  
  cout << here << " => Checking status of optics matricies.." << endl; 
  
  cout << "RHRS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; 
  for (UInt_t m=0; m<matrices.size(); m++) { 
    
    auto m_it = THRS::fMatrixMap_R.find( matrices.at(m) ); 
    
    auto matrix = m_it->second; 
    
    TString name = matrix_names.find(matrices.at(m))->second; 
    
    cout << TString::Format("Matrix '%s', nElems=%2i", 
			    name.Data(), 
			    (int)matrix->fElems.size() ) << endl; 
  }
  
  cout << "overall status = "; 
  switch (THRS::OpticsStatus_RHRS) 
    { 
    case THRS::Init_ok      : cout << "Ok" << endl; break; 
    case THRS::Init_notDone : cout << "Not initialized" << endl; good=false; break; 
    case THRS::Init_error   : cout << "Error!" << endl;          good=false; break; 
    }
  
  
  cout << "LHRS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; 
  for (UInt_t m=0; m<matrices.size(); m++) { 
    
    auto m_it = THRS::fMatrixMap_L.find( matrices.at(m) ); 
    
    auto matrix = m_it->second; 
    
    TString name = matrix_names.find(matrices.at(m))->second; 
    
    cout << TString::Format("Matrix '%s', nElems=%2i", 
			    name.Data(), 
			    (int)matrix->fElems.size() ) << endl; 
  }
  
  cout << "overall status = "; 
  switch (THRS::OpticsStatus_LHRS) 
    { 
    case THRS::Init_ok      : cout << "Ok" << endl; break; 
    case THRS::Init_notDone : cout << "Not initialized" << endl; good=false; break; 
    case THRS::Init_error   : cout << "Error!" << endl;          good=false; break; 
    }
  
  return good; 
}
////////////////////////////////////////////////////////////////////////////////////
double THRS::TMatrixElement::Evaluate_poly( const double x_fp ) const { 
  
  const TString here = "THRS::TMatrixElement::Evaluate_poly()"; 
    
  double val=0.0; 
  
  for (int p=fPoly.size()-1; p>=0; p--) { val = val*x_fp + fPoly.at(p); }
  
  return val; 
}
////////////////////////////////////////////////////////////////////////////////////
double THRS::TMatrix::Evaluate( const double fp_x, 
				const double fp_y[fMatrix_MAX_EXP], 
				const double fp_theta[fMatrix_MAX_EXP], 
				const double fp_phi[fMatrix_MAX_EXP]    ) const { 
  
#if DEBUG
  cout << "THRS::TMatrix::Evaluate() - matrix n.=" << fType << endl; 
  
  cout << TString::Format("inputs = [x=%0.3e, y=%0.3e, th=%0.3e, ph=%0.3e]", 
			  fp_x, fp_y[0], fp_theta[0], fp_phi[0]) << endl; 
#endif   

  bool is_fpMatrix = fElems.size()==1; 
  
  double val = 0.;
  
  for (UInt_t i=0; i<fElems.size(); i++) { 
    
    TMatrixElement *elem = fElems.at(i); 
    
    //compute the x-fp polynomial of this element
    double Ei = elem->Evaluate_poly( fp_x ); 
    
    if (is_fpMatrix) return Ei; 
    
#if DEBUG 
    cout << TString::Format("elem-(%2i), exp=[%i,%i,%i,%i], poly=[%0.3e], ", 
			    (int)i, 
			    elem->fExp[0],elem->fExp[1],elem->fExp[2],elem->fExp[3],
			    Ei); 
#endif 
    
    //if this ISN'T a fp-matrx, then we will add the powers of the other vars
    if (elem->fExp[0]>0) Ei *= fp_theta[elem->fExp[0]-1]; 
    if (elem->fExp[1]>0) Ei *= fp_y    [elem->fExp[1]-1]; 
    if (elem->fExp[2]>0) Ei *= fp_phi  [elem->fExp[2]-1]; 
    
    //if (elem->fExp[3]>0) Ei *= TMath::Abs(fp_theta[elem->fExp[3]-1]); 
    
#if DEBUG 
    cout << TString::Format("M=[%0.3e]", Ei) << endl; 
#endif 
    
    val += Ei;
  }

#if DEBUG 
    cout << TString::Format("final val = %0.3e, ", val) << endl; 
#endif 
    
  
  return val; 
}
double THRS::TMatrix::Evaluate( const double fp_x ) const { 
  
  //same as above, but for single-element fp-matrices
  const TString here = "THRS::Evaluate()";   
  
  if (fElems.size()!=1) { 
    THRS::ReportError(here,"This is not a FP matrix!"); 
    return -1e30; 
  }
    
  return fElems.at(0)->Evaluate_poly( fp_x );
}
void THRS::Compute_trackOptics( const bool arm, 
				const double x_tr, 
				const double y_tr, 
				const double dxdz,
				const double dydz, 
				double &tg_y,
				double &tg_dxdz, 
				double &tg_dydz, 
				double &tg_dp    ) { 
  
  //compute the 'focal plane' coordinates of the track, as well as the target 
  // coordinates. 
  // These computations are informed mainly by John Williamson's thesis
  // 'APEX (A'Epxeriment): The Search for a Dark Photon at Jefferson Lab' (2022), 
  // specifically, see section 3.5.5 as well as section 4.5.  
  
  const TString here = "THRS::Compute_trackOptics()"; 
  
  //check optics status 
  HRS_initStatus optics_status = arm 
    ? THRS::OpticsStatus_RHRS : THRS::OpticsStatus_LHRS; 
  
  switch (optics_status) 
    { 
    case Init_ok     : break; //init status ok 

    case Init_notDone: 
      THRS::ReportError(here,"Optics not initalized."); 
      return; 
    case Init_error  : 
      THRS::ReportError(here,"Problem with optics init."); 
      return; 
    default          :
      THRS::ReportError(here,"Unknown error with optics init."); 
      return; 
    }
  
  //get optics matrix-map 
  map<MatrixType,TMatrix*> matMap = arm ? THRS::fMatrixMap_R : THRS::fMatrixMap_L; 
  
  double fp_x = x_tr; 
  
  double fp_y[TMatrix::fMatrix_MAX_EXP+1]; 
  double fp_theta[TMatrix::fMatrix_MAX_EXP+1]; 
  double fp_phi[TMatrix::fMatrix_MAX_EXP+1]; 
    
  //focal-plane y
  //note: these 'fp-matrices' are only polynomials in powers of 'fp_x'
  fp_y[0] = y_tr - matMap.find( Mat_fp_y )->second->Evaluate( fp_x ); 
  
  //get tan(rho), which is the 'local central ray' of particles
  double tan_rho = matMap.find( Mat_tan_rho )->second->Evaluate( fp_x ); 
  
  double cos_rho = 1./TMath::Sqrt( 1. + tan_rho*tan_rho ); 
  
  //now, compute fp_theta & fp_phi (first order)
  fp_theta[0] = (dxdz + tan_rho)/(1. - dxdz*tan_rho); 
    
  fp_phi[0] 
    = ( dydz - matMap.find( Mat_fp_phi )->second->Evaluate( fp_x ) )
    / (1. - dxdz*tan_rho) / cos_rho;   
    
  //we've now evaluated all of our focal-plane coords. let's find their higher-pow's
  for (int p=2; p<=TMatrix::fMatrix_MAX_EXP; p++) { 
    
    fp_y[p-1]     = pow( fp_y[0],     p ); 
    fp_theta[p-1] = pow( fp_theta[0], p ); 
    fp_phi[p-1]   = pow( fp_phi[0],   p ); 
  }
  
  //now, we're ready to evaluate target-coords. 
  //target_y
  tg_y    = matMap.find( Mat_tg_y     )->second->Evaluate( fp_x, 
							   fp_y, 
							   fp_theta, 
							   fp_phi ); 
  //target_theta
  tg_dxdz = matMap.find( Mat_tg_theta )->second->Evaluate( fp_x, 
							   fp_y, 
							   fp_theta, 
							   fp_phi ); 
  //target_phi
  tg_dydz = matMap.find( Mat_tg_phi   )->second->Evaluate( fp_x, 
							   fp_y, 
							   fp_theta, 
							   fp_phi ); 
  //target_dp
  tg_dp   = matMap.find( Mat_tg_dp    )->second->Evaluate( fp_x, 
							   fp_y, 
							   fp_theta, 
							   fp_phi ); 
}
////////////////////////////////////////////////////////////////////////////////////
void THRS::ReportError( const TString location, 
			const TString message ) { 
  
  cerr << endl 
       << TString::Format("ERROR in %s: %s", location.Data(), message.Data()) 
       << endl; 
  return; 
}
////////////////////////////////////////////////////////////////////////////////////
void THRS::VDC::Parse_offsets( const TString path_DB ) { 
  
  const TString here = "THRS::VDC::Parse_offsets()"; 
    
  cout << TString::Format("%s => parsing data from \"%s\"...",
			  here.Data(),
			  path_DB.Data()) 
       << flush; 
  
  fstream db; db.open( path_DB.Data(), ios::in );   
  
  //failed to open file! 
  if (!db) { 
    THRS::ReportError(here,"File could not be opened!"); 
    return; 
  }
    
  TString str; 
  
  bool arm; 
  
  while ( !db.eof() ) { 
    
    //check which arm we're looking at
    db >> str; 
    
    if (str!="R" && str!="L") { 
      THRS::ReportError(here,TString::Format("Bad arm-name given! \"%s\"", 
					     str.Data()) ); 
      return; 
    }
    
    if (str=="R") { arm=THRS::Arm_Right; } else { arm=THRS::Arm_Left; }
    
    //check which plane we're looking at
    db >> str; 
    
    int plane = atoi(str); 
    
    if (plane>=VDC::fN_PLANES || plane<0) { 
      THRS::ReportError(here,"Bad plane number given!"); 
      return; 
    }
    
    //read offsets
    int wire=0; db >> str; 
    while (str != "eol") { 
      
      if (arm==THRS::Arm_Right) { 
	VDC::offset_RHRS[plane][wire] = stod((string)str); 
      } else                    { 
	VDC::offset_LHRS[plane][wire] = stod((string)str); 
      }
      
      wire++; 
      db >> str; 
    }
    
    if (wire > VDC::fN_WIRES || wire < VDC::fN_WIRES) { 
      THRS::ReportError(here,
			TString::Format("bad num of wires! p-%i, n.wires=%i",
					plane, wire)); 
      
      if (arm==THRS::Arm_Right) { VDC::status_RHRS[plane] = VDC::Error; }
      else                      { VDC::status_LHRS[plane] = VDC::Error; }
      
    } else { 
      
      //if we've gotten here, then we've parsed the offsets successfully
      if (arm==THRS::Arm_Right) { VDC::status_RHRS[plane] = VDC::Ok; }
      else                      { VDC::status_LHRS[plane] = VDC::Ok; } 
    }
    
    if (plane==VDC::fN_PLANES-1) break; 
  }
  
  /*for (int p=0; p<VDC::fN_PLANES; p++) { 
    
    copy( begin(R_abs[p]), end(R_abs[p]), begin(VDC::offset_RHRS[p]) ); 
    copy( begin(L_abs[p]), end(L_abs[p]), begin(VDC::offset_LHRS[p]) ); 
    }*/ 
  
  
  
  cout << "done" << endl; 
  return; 
}
////////////////////////////////////////////////////////////////////////////////////
bool THRS::VDC::Print_status() const { 
  
  //print status for offset-parsing for all planes. 
  const TString here = "THRS::VDC::Print_status()"; 
  
  cout << here << " => Parsing status for all planes:" << endl; 
  
  bool allGood=true; 
  
  cout << "HRS - Right~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; 
  for (int p=0; p<VDC::fN_PLANES; p++) { 
    
    cout << TString::Format("plane %i, status = ",p); 
    
    switch (VDC::status_RHRS[p]) { 
    case VDC::Ok      : cout << "Ok" << endl;              break; 
    case VDC::NotInit : cout << "Not initialized" << endl; allGood=false; break;  
    case VDC::Error   : cout << "Error!" << endl;          allGood=false; break; 
    default           : cout << "Unknown status!" << endl; allGood=false; break; 
    }
  } 
  cout << "HRS - Left~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; 
  for (int p=0; p<VDC::fN_PLANES; p++) { 
    
    cout << TString::Format("plane %i, status = ",p); 
    
    switch (VDC::status_LHRS[p]) { 
    case VDC::Ok      : cout << "Ok" << endl;              break; 
    case VDC::NotInit : cout << "Not initialized" << endl; allGood=false; break;  
    case VDC::Error   : cout << "Error!" << endl;          allGood=false; break; 
    default           : cout << "Unknown status!" << endl; allGood=false; break; 
    }
  } 
    
  cout << "Overall initialization status: " << (allGood?"good.":"bad!") << endl; 
  
  return allGood; 
}
////////////////////////////////////////////////////////////////////////////////////
double THRS::VDC::WirePos( const bool arm, 
			   const int  plane, 
			   const int  wireNum ) { 
  
  const TString here = "THRS::VDC::WirePos()"; 
  
  //check if plane is valid
  if (plane >= THRS::VDC::fN_PLANES || plane < 0) {
    
    string msg 
      = (string)TString::Format("second arg. out-of-range: planes=[0,%i], given %i",
				THRS::VDC::fN_PLANES, 
				plane); 
    
    THRS::ReportError( here, msg ); return -1e30; 
  }
  
  //check if wire is valid
  if (wireNum >= THRS::VDC::fN_WIRES || wireNum < 0) {
    
    string msg 
      = (string)TString::Format("third arg. out-of-range: wires=[0,%i], given %i",
				THRS::VDC::fN_WIRES, 
				wireNum); 
  
    THRS::ReportError( here, msg ); return -1e30; 
  }
  
  //the position of the closest wire (with a lower u- or v-value) 
  double x_offset = (arm==THRS::Arm_Right) 
    ? VDC_wire0_R[plane]
    : VDC_wire0_L[plane]; 
     
  return x_offset  -  THRS::VDC::wireSpacing*((double)wireNum);    
}
////////////////////////////////////////////////////////////////////////////////  
int    THRS::VDC::WireNum( const bool   arm, 
			   const int    plane, 
			   const double x ) { 
  
  const TString here = "THRS::VDC::WireNum()"; 
  
  //check if plane is valid
  if (plane >= THRS::VDC::fN_PLANES || plane < 0) {
    
    string msg 
      = (string)TString::Format("second arg. out-of-range: planes=[0,%i], given %i",
				THRS::VDC::fN_PLANES, 
				plane); 
    
    THRS::ReportError( here, msg ); return -999; 
  }
  
  double x_offset = (arm==THRS::Arm_Right) 
    ? VDC_wire0_R[plane]
    : VDC_wire0_L[plane]; 
  
  return (int)TMath::Nint( (x_offset - x)/THRS::VDC::wireSpacing );  
}
double THRS::VDC::RealTime( const bool arm, 
			    const int plane, 
			    const int wireNum, 
			    const double rawTime ) {
  
  const TString here = "THRS::VDC::RealTime()"; 
       
  //check if plane is valid
  if (plane >= THRS::VDC::fN_PLANES || plane < 0) {
    
    string msg 
      = (string)TString::Format("second arg. out-of-range: planes=[0,%i], given %i",
				THRS::VDC::fN_PLANES, 
				plane); 
    
    THRS::ReportError( here, msg ); return -1e30; 
  }
  
  //check if it's been initalized
  if ( (arm==THRS::Arm_Right && VDC::status_RHRS[plane] != VDC::Ok) ||
       (arm==THRS::Arm_Left  && VDC::status_LHRS[plane] != VDC::Ok) ) { 
    
    THRS::ReportError( here, "VDC offsets not correctly initialized!" ); 
    return -1e30; 
  }
  
  if (wireNum >= THRS::VDC::fN_WIRES || wireNum < 0) {
    
    string msg 
      = (string)TString::Format("third arg. out-of-range: wires=[0,%i], given %i",
				THRS::VDC::fN_WIRES, 
				wireNum); 
    
    THRS::ReportError( here, msg ); return -1e30; 
  }
  //vdc wire offsets
  
  Double_t offset = (arm==THRS::Arm_Right) 
    ? VDC::offset_RHRS[plane][wireNum]
    : VDC::offset_LHRS[plane][wireNum]; 
    
  return (THRS::VDC::TDC_resolution * rawTime)  +  offset; 
}
double THRS::VDC::RawTime( const bool arm, 
			   const int plane, 
			   const int wireNum, 
			   const double realTime ) {
  
  const TString here = "THRS::VDC::RealTime()"; 
  
  //check if plane is valid
  if (plane >= THRS::VDC::fN_PLANES || plane < 0) {
    
    string msg 
      = (string)TString::Format("second arg. out-of-range: planes=[0,%i], given %i",
				THRS::VDC::fN_PLANES, 
				plane); 
    
    THRS::ReportError( here, msg ); return -1e30; 
  }
  //check if it's been initalized
  if ( (arm==THRS::Arm_Right && VDC::status_RHRS[plane] != VDC::Ok) ||
       (arm==THRS::Arm_Left  && VDC::status_LHRS[plane] != VDC::Ok) ) { 
    
    THRS::ReportError( here, "VDC offsets not correctly initialized!" ); 
    return -1e30; 
  }
  
  
  if (wireNum >= THRS::VDC::fN_WIRES || wireNum < 0) {
    
    string msg 
      = (string)TString::Format("third arg. out-of-range: wires=[0,%i], given %i",
				THRS::VDC::fN_WIRES, 
				wireNum); 
    
    THRS::ReportError( here, msg ); return -1e30; 
  }
  //vdc wire offsets
  
  Double_t offset = (arm==THRS::Arm_Right) 
    ? VDC::offset_RHRS[plane][wireNum]
    : VDC::offset_LHRS[plane][wireNum]; 
    
  return (realTime  - offset)/THRS::VDC::TDC_resolution; 
}
////////////////////////////////////////////////////////////////////////////////////
double THRS::VDC::w( const bool arm, const int plane ) { 
  
  const string here = "THRS::VDC::w()"; 
  
  //check if plane is valid
  if (plane >= THRS::VDC::fN_PLANES || plane < 0) {
    
    string msg 
      = (string)TString::Format("second arg. out-of-range: planes=[0,%i], given %i",
				THRS::VDC::fN_PLANES, 
				plane); 
    
    THRS::ReportError( here, msg ); return -1e30; 
  }
  
  return (arm==THRS::Arm_Right) 
    ? VDC_w_R[plane] 
    : VDC_w_L[plane]; 
}
////////////////////////////////////////////////////////////////////////////////////
TS2Hit::TS2Hit( bool arm, int paddle, double T_pmtL, double T_pmtR ) {
  
  f_isRightArm =arm; 
  fPaddle      =paddle; 
  
  fRawTime_pmtL =T_pmtL; 
  fRawTime_pmtR =T_pmtR; 
  
  fTime = Compute_RealTime(); 
  
  fIsCoinc = (fTime > -1e29); 

  //get the location of the paddle
  fZ = f_isRightArm ? 3.3098 : 3.1790 ; 
  
  fX = ((double)fPaddle)*fPaddleWidth_X  + (f_isRightArm ? -1.21413 : -1.16913); 
  
  fY = 0; 
}
double TS2Hit::Compute_RealTime() {
  
  //check to make sure that the PMTs registered non-null times
  if (TMath::Abs(fRawTime_pmtL) > 1e7 || 
      TMath::Abs(fRawTime_pmtR) > 1e7) return -1e30; 
  
  if (f_isRightArm) { //RHRS 
    
    fRealTime_pmtR = ( Right_R[fPaddle] - fRawTime_pmtR )*fTDC_resolution;
    fRealTime_pmtL = ( Right_L[fPaddle] - fRawTime_pmtL )*fTDC_resolution; 
  
  } else      { 
    
    fRealTime_pmtR = ( Left_R[fPaddle]  - fRawTime_pmtR )*fTDC_resolution; 
    fRealTime_pmtL = ( Left_L[fPaddle]  - fRawTime_pmtL )*fTDC_resolution; 
  }
    
  //check to make sure the raw times agree 
  
  if (TMath::Abs(fRealTime_pmtR-fRealTime_pmtL) > 7.5e-9) return -1e30; 
  
  return 0.5*(fRealTime_pmtR + fRealTime_pmtL); 
}
void TS2Hit::Make_twinHit( TS2Hit *neighbor ) { 
  
  fTime = 0.5*(fTime + neighbor->Time()); 
  
  fX += 0.5*fPaddleWidth_X; 
  
  f_isTwinHit = true; 
} 
/////////////////////////////////////////////////////////////////////////////////
TvdcHit::TvdcHit( int plane, 
		  double wire, 
		  double rawTime, 
		  TEventHandler *event ) { 
  
  fEvent = event; 
  
  f_isRightArm = fEvent->ActiveArm(); 
  
  fPlane    = plane; 
  fRawTime  = rawTime; 
  fWireNum  = (int)wire; 
  fWirePos  = THRS::VDC::WirePos( f_isRightArm, 
				  fPlane,
				  fWireNum );
  
  fRealTime 
    = THRS::VDC::RealTime( f_isRightArm, fPlane, fWireNum, fRawTime )  -  fEvent->GetS2Hit()->Time();
    
  fW = f_isRightArm ? VDC_w_R[fPlane] : VDC_w_L[fPlane]; 
}
void   TvdcHit::SetRawTime(const double rawTime)  { 
  
  fRawTime  = rawTime; 
  fRealTime = GetRealTime( fRawTime ) - fEvent->GetS2Hit()->Time(); 
} 
double TvdcHit::GetWirePos( int wire )      const {
  
  if (f_isRightArm) {
    
    Double_t v0[4] = {0.77852, 0.77852, 1.02793, 1.02793}; 
    return v0[fPlane] - fWireSpacing*((double)wire);
    
  } else { 
    
    Double_t v0[4] = {0.77852, 0.77852, 1.02718, 1.02718};
    return v0[fPlane] - fWireSpacing*((double)wire);
  }
}
double TvdcHit::GetRealTime(double rawTime) const {
  
  //vdc wire offsets
  
  Double_t offset = (f_isRightArm) 
    ? R_abs[fPlane][fWireNum]-1.04e-9 : L_abs[fPlane][fWireNum]+3.20e-9; 
    
  return (-0.5e-9)*rawTime + offset; 
}
double TvdcHit::RawTime( const bool is_RightArm, 
			 const int plane, 
			 const int wireNum, 
			 const double realTime ) { 
  
  double offset = is_RightArm ? R_abs[plane][wireNum] : L_abs[plane][wireNum]; 
  
  return (realTime - offset)/(-0.5e-9); 
}
int    TvdcHit::WireNum( const bool   is_RightArm, 
			 const int    plane, 
			 const double wirePos ) { 
  
  //the position of the closest wire (with a lower u- or v-value) 
  double x_offset; 
  if (is_RightArm) { x_offset = plane<2 ? 0.77852 : 1.02793; }
  else             { x_offset = plane<2 ? 0.77852 : 1.02718; }
  
  return (int)TMath::Nint( (x_offset - wirePos)/TvdcHit::WireSpacing() );  
}
double TvdcHit::WirePos( const bool is_RightArm, 
			 const int  plane, 
			 const int  wireNum ) { 
  
  //the position of the closest wire (with a lower u- or v-value) 
  double x_offset; 
  if (is_RightArm) { x_offset = plane<2 ? 0.77852 : 1.02793; }
  else             { x_offset = plane<2 ? 0.77852 : 1.02718; }
      
  return x_offset  -  TvdcHit::WireSpacing()*((double)wireNum);    
}
////////////////////////////////////////////////////////////////////////////////  
THitGroup::THitGroup(int plane) { 
    
  fPlane=plane; 
}  
THitGroup::~THitGroup() {  
  //delete all constituent hits
  for (unsigned int h=0; h<fHits.size(); h++) fHits.at(h)->~TvdcHit(); 
}
void THitGroup::AddHit( double wire, double rawTime ) { 
  
  fHits.push_back( new TvdcHit(fPlane,wire,rawTime) ); 
}
double THitGroup::WirePos( unsigned int h ) const { 
  
  if (h>=Nhits()) return kNull_double; 
  
  return fHits.at(h)->wPos();
}
int    THitGroup::WireNum( unsigned int h ) const { 
  
  if (h>=Nhits()) return kNull_int; 
  
  return fHits.at(h)->wNum(); 
}
double THitGroup::Time( unsigned int h )    const { 
  
  if (h>=Nhits()) return kNull_double; 
  
  return fHits.at(h)->Time(); 
}
int    THitGroup::FirstWire() const { return WireNum(0); }
double THitGroup::LoEdge()    const { return WirePos(Nhits()-1); }
double THitGroup::HiEdge()    const { return WirePos(0); }
////////////////////////////////////////////////////////////////////////////////  
THitCluster::THitCluster( THitGroup *group, 
			  const double intercept, 
			  const double eta ) { 
  fGroup     =group; 
  fIntercept =intercept; 
  fEta_score =eta; 
}
///////////////////////////////////////////////////////////////////////////////
TChamberPair::TChamberPair( bool is_loChamber,
			    double u,
			    double v,
			    THitGroup *Group_U,
			    THitGroup *Group_V, 
			    int unique_id ) {
  
  f_isLoChamber = is_loChamber; 
  fu = u; 
  fv = v; 
  fGroup_U = Group_U;
  fGroup_V = Group_V; 
  
  fUnique_ID = unique_id; 
}
TChamberPair::TChamberPair( bool is_loChamber, 
			    THitCluster *clust_u,
			    THitCluster *clust_v, 
			    int unique_id ) { 

  f_isLoChamber = is_loChamber; 
  
  fu = clust_u->Intercept(); 
  fv = clust_v->Intercept(); 
  
  fGroup_U = clust_v->GetGroup(); 
  fGroup_V = clust_u->GetGroup(); 
  
  fUnique_ID = unique_id; 
}
void TChamberPair::Remove_track( TObject* track ) { 
  
  //cout << "TChamberPair:: removing track.."<< endl; 
  
  auto trk_it = find( begin(fTracks), end(fTracks), track ); 
  
  if (trk_it==fTracks.end()) {
    cout << "WARNING: track not found in TChamberPair" << endl; 
    return; 
  }
  
  fTracks.erase( trk_it ); 
}

  
double TChamberPair::ClosestWirePos_Lo( double x ) const { 
  
  //the position of the closest wire (with a lower u- or v-value) 
  
  double vOffR[4] = {0.77852, 0.77852, 1.02793, 1.02793}; 
  double vOffL[4] = {0.77852, 0.77852, 1.02718, 1.02718};
  
  int fPlane = f_isLoChamber ? 0 : 2; 
  
  double planeOffset = fGroup_U->IsRightArm() ? vOffR[fPlane] : vOffL[fPlane]; 
    
  double wireNum = (double)ceil( (planeOffset - x)/TvdcHit::WireSpacing() );  
    
  return planeOffset  -  TvdcHit::WireSpacing()*wireNum;    
}
double TChamberPair::ClosestWirePos( const double x ) const { 
  
  double vOffR[4] = {0.77852, 0.77852, 1.02793, 1.02793}; 
  double vOffL[4] = {0.77852, 0.77852, 1.02718, 1.02718};
  
  int fPlane = f_isLoChamber ? 0 : 2; 
  
  double planeOffset = fGroup_U->IsRightArm() ? vOffR[fPlane] : vOffL[fPlane]; 
    
  double wireNum = TMath::Nint( (planeOffset - x)/TvdcHit::WireSpacing() );  
    
  return planeOffset  -  TvdcHit::WireSpacing()*wireNum;    
}
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


using namespace std;

TEventHandler::TEventHandler( bool arm, 
			      double beamCurrent,
			      int runNumber, 
			      TS2Hit *fHit_R, 
			      TS2Hit *fHit_L ) { 
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

double TEventHandler::Get_tauSigma() const { 
  
  if (fRunNumber < FIRST_FIXED_RUN && ActiveArm()==false ) { 
    
    return Interpolate( fBeamCurrent, f_BEAM_CURRENT, f_SIGMA_INTERPOLATE, 2 ); 
  } 
  return f_SIGMA_FIXED; 
}
double TEventHandler::Drift_X( double tau, double slope, int derivative ) const { 
  
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
double TEventHandler::Drift_T_raw( const double x, 
				   const double *par, 
				   const int derivative ) const {
  
  double x1(par[0]), y1(par[1]), m(par[2]), b(par[3]), y0(par[4]); 
  
  double absX = TMath::Abs(x); 
  
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
double TEventHandler::Drift_T( double x, double slope, int derivative ) const { 
  
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

  if (TMath::Abs(x) < par[0]+blur_sigma*5.) {
        
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
double TEventHandler::Interpolate( const double x, 
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
//////////////////////////////////////////////////////////////////////////////////
void TvdcTrack::UpdateTrackInfo() { 
  
  //check to see if this track has VDC data. 
  if ( !fPair_Lo || !fPair_Hi ) { f_hasVDCdata=false; } else { f_hasVDCdata=true; }
    
  //tell the point-pairs about their actual slope-values
  double wSep = fW[2]-fW[0]; //the U2-U1 & V2-V1 separations are the same
  
  if (f_hasVDCdata) { 
    //if we DO have vdc data, we want to keep track of the 'pair' intercepts as 
    // best we can.
    fIntercept[0] = fPair_Lo->v(); 
    fIntercept[1] = fPair_Lo->u(); 
    fIntercept[2] = fPair_Hi->v(); 
    fIntercept[3] = fPair_Hi->u(); 
  } 
  
  double 
    uLo(fIntercept[1]),
    vLo(fIntercept[0]),
    uHi(fIntercept[3]),
    vHi(fIntercept[2]); 
  
  fSlope_u = wSep/(uHi-uLo); 
  fSlope_v = wSep/(vHi-vLo);  
  
  if (f_hasVDCdata) { 
    fPair_Hi->SetSlope_uv( fSlope_u, fSlope_v ); 
    fPair_Lo->SetSlope_uv( fSlope_u, fSlope_v ); 
  }
  
  f_S2Int_xyz = ComputeIntercept_z( fEvent->GetS2Hit()->Z() ); 
  
  //we need to find the 'speed' of the track in the u & v directions. 
  //the 'true' speed of the track thru the VDC is effectivley c, but 
  //since the track isn't travelling parallel to the 'u' or 'v' axes, we 
  //need to adjust accordingly. this will be used to compute Time-of-Flight (ToF)
  //for the track at any given wire crossing. 
  TVector3 S2Int_uvw =  TvdcTrack::Rotate_xyz_to_uvw( f_S2Int_xyz ); 
  
  fS2_u = S2Int_uvw[0]; 
  fS2_v = S2Int_uvw[1]; 
    
  TVector3 S_uvw  = TVector3( 1./fSlope_u, 1./fSlope_v, 1. ).Unit(); 
  
  //these are the 'speeds' of the particle in the 'U' or 'V' direction, assuming 
  // the overall speed of the particle is 'c'. 
  fC_u = fC * S_uvw[0]; 
  fC_v = fC * S_uvw[1]; 
  
  //compute focal plane intercept
  f_FPInt_xyz = ComputeIntercept_z(0.); 
    
  TVector3 S_xyz = TvdcTrack::Rotate_uvw_to_xyz( S_uvw ); 
  
  fTheta = TMath::ATan( S_xyz.X()/S_xyz.Z() ); 
  fPhi   = TMath::ATan( S_xyz.Y()/S_xyz.Z() ); 
}
TvdcTrack::TvdcTrack( TEventHandler *event, 
		      TChamberPair *pLo, 
		      TChamberPair *pHi ) { 
  
  fEvent = event; 
  
  f_isRightArm = fEvent->ActiveArm(); 
  
  if (f_isRightArm) { //
    
    copy( begin(VDC_w_R), end(VDC_w_R), begin(fW) );  
  
  } else                   {
    
    copy( begin(VDC_w_L), end(VDC_w_L), begin(fW) );   
  }
  
  //if we're making a new monte-carlo track, skip this step. 
  if ( !pLo || !pHi ) { f_hasVDCdata=false; } else { f_hasVDCdata=true; }
  
  if (f_hasVDCdata) { 
    
    fPair_Lo = pLo;   fPair_Lo->Add_track( this );
    fPair_Hi = pHi;   fPair_Hi->Add_track( this ); 
    
    fGroup[0] = (THitGroup*)fPair_Lo->GetGroup_U(); 
    fGroup[1] = (THitGroup*)fPair_Lo->GetGroup_V();
    
    fGroup[2] = (THitGroup*)fPair_Hi->GetGroup_U(); 
    fGroup[3] = (THitGroup*)fPair_Hi->GetGroup_V(); 
  
    UpdateTrackInfo(); 
  }
  
}

const double 
R_xu( 0.500000),  R_xv(0.500000), R_xw(-0.707107), 
  R_yu(-0.707107),  R_yv(0.707107), R_yw( 0.0),
  R_zu( 0.500000),  R_zv(0.500000), R_zw(0.707107); 

void TvdcTrack::Set_Errors( double errors[5] ) { 
  
  for (int p=0; p<4; p++) fIntercept_ERR[p] =errors[p]; 

  fT0_ERR = errors[4]; 
  
  /*/compute Theta-error
    TVector3 S_xyz 
    = TvdcTrack::Rotate_uvw_to_xyz( TVector3(1./fSlope_u, 1./fSlope_v, 1.) ); 
  
  double sx(S_xyz[0]), sy(S_xyz[1]), sz(S_xyz[2]); 
  
  double d_Theta = 0.; 
  
  double d_arctan = 1./( 1. + (sx*sx)/(sz*sz) ); 
  
  double Dw = f_isRightArm 
    ? VDC_w_R[2]-VDC_w_R[0] 
    : VDC_w_L[2]-VDC_w_L[0]; //separation between vdc planes
  
  //U1 - plane
  d_Theta
    +=  TMath::Power( fIntercept_ERR[0]
		      * d_arctan * (sz*R_xv/Dw - sx*R_zv/Dw)/(sz*sz), 2 ); 
  
  //U2 - plane
  d_Theta 
    +=  TMath::Power( fIntercept_ERR[2]
		      * d_arctan * (sz*R_xv/Dw - sx*R_zv/Dw)/(sz*sz), 2 ); 
  
  //V1 - plane
  d_Theta 
    +=  TMath::Power( fIntercept_ERR[1]
		      * d_arctan * (sz*R_xu/Dw - sx*R_zu/Dw)/(sz*sz), 2 ); 
  
  //V2 - plane
  d_Theta 
    +=  TMath::Power( fIntercept_ERR[3]
		      * d_arctan * (sz*R_xu/Dw - sx*R_zu/Dw)/(sz*sz), 2 ); 
  
  
		      fTheta_ERR = TMath::Sqrt( d_Theta ); //*/ 
    
}
void   TvdcTrack::Set_Eta(double eta[4]) { 

  for (int p=0; p<4; p++) Set_Eta(p, eta[p]); 
} 
double TvdcTrack::Get_Eta() const { 
  
  double eta(0.); 
  for (int p=0; p<4; p++) eta += Get_Eta(p); 
  return eta; 
}
int    TvdcTrack::Get_nGoodPoints() const { 
  
  int nPts(0); 
  for (int p=0; p<4; p++) nPts += Get_nGoodPoints(p); 
  return nPts; 
} 
double TvdcTrack::Get_RMS() const { 
  
  double RMS(0.); 
  for (int p=0; p<4; p++) RMS += 0.25*TMath::Power( Get_RMS(p), 2. ); 
  return TMath::Sqrt( RMS ); 
} 
double TvdcTrack::Error_allIntercept() const { 
  
  double RMS(0.); 
  for (int p=0; p<4; p++) RMS += 0.25*TMath::Power( Error_intercept(p), 2. ); 
  return TMath::Sqrt( RMS ); 
}
void TvdcTrack::Set_S2int_angles( double s2x, 
				  double s2y, 
				  double theta, 
				  double phi )   { 
  
  //this is an alternate way to initialize the track without VDC data 
  // (used for monte-carlo processing, in which a track's position is created 
  //  before the corresponding VDC hits are created). 
  
  if ( !fPair_Lo || !fPair_Hi ) { f_hasVDCdata=false; } else { f_hasVDCdata=true; }
    
  f_S2Int_xyz = TVector3( s2x, s2y, fEvent->GetS2Hit()->Z() ); 
  
  TVector3 S_xyz( TMath::Tan(theta), TMath::Tan(phi), 1. ); 
  
  TVector3 S_uvw = TvdcTrack::Rotate_xyz_to_uvw( S_xyz ).Unit(); 
  
  fSlope_u = S_uvw.Z()/S_uvw.X(); 
  fSlope_v = S_uvw.Z()/S_uvw.Y(); 
    
  f_FPInt_xyz = TVector3( f_S2Int_xyz - S_xyz*f_S2Int_xyz.Z() ); 
    
  auto S2u = TvdcTrack::Rotate_xyz_to_uvw( f_S2Int_xyz ); 
  auto FPu = TvdcTrack::Rotate_xyz_to_uvw( f_FPInt_xyz ); 
    
  fIntercept[0] = ( fW[0]-FPu.Z() )/fSlope_v  +  FPu.Y(); 
  
  fIntercept[1] = ( fW[1]-FPu.Z() )/fSlope_u  +  FPu.X();
  
  fIntercept[2] = ( fW[2]-FPu.Z() )/fSlope_v  +  FPu.Y(); 
  
  fIntercept[3] = ( fW[3]-FPu.Z() )/fSlope_u  +  FPu.X();

  TvdcTrack::Compute_Theta_Phi( f_isRightArm, 
				fIntercept, 
				fTheta, fPhi ); 
  
  //fTheta = TMath::ATan( S_xyz.X()/S_xyz.Z() ); 
  //fPhi   = TMath::ATan( S_xyz.Y()/S_xyz.Z() ); 
  
  fS2_u = S2u[0]; 
  fS2_v = S2u[1]; 
    
  fC_u = fC * S_uvw[0]; 
  fC_v = fC * S_uvw[1]; 
  
  UpdateTrackInfo(); 
}
void TvdcTrack::SetPair_Hi( TChamberPair *pHi ) { 
  
  //if the old pair exists, tell it that we're breaking up with it
  // (it's not you, it's me..) 
  if (fPair_Hi) fPair_Hi->Remove_track( this );
  
  fPair_Hi = pHi; fPair_Hi->Add_track( this ); 
  
  fGroup[2] = fPair_Hi->GetGroup_U(); 
  fGroup[3] = fPair_Hi->GetGroup_V(); 
  UpdateTrackInfo(); 
}
void TvdcTrack::SetPair_Lo( TChamberPair *pLo ) { 
  
  //if the old pair exists, tell it that we're breaking up with it
  // (it's not you, it's me..) 
  if (fPair_Lo) fPair_Lo->Remove_track( this );
  
  fPair_Lo = pLo; fPair_Lo->Add_track( this ); 
  
  fGroup[0] = fPair_Lo->GetGroup_U(); 
  fGroup[1] = fPair_Lo->GetGroup_V();
  UpdateTrackInfo(); 
} 
void TvdcTrack::Set_uv_Hi( const double u, 
			   const double v ) { 
  if(fPair_Hi) { fPair_Hi->Set_uv( u,v ); } 
  else         { fIntercept[2]=v; fIntercept[3]=u; }
  UpdateTrackInfo(); 
}
void TvdcTrack::Set_uv_Lo( const double u, 
			   const double v ) { 
  if(fPair_Lo) { fPair_Lo->Set_uv( u,v ); } 
  else         { fIntercept[0]=v; fIntercept[1]=u; }
  UpdateTrackInfo(); 
}
double TvdcTrack::Slope_u() { return fSlope_u; }
double TvdcTrack::Slope_v() { return fSlope_v; }
  
void   TvdcTrack::Set_intercept(int plane, double x) { 
  
  if (plane==0) { if (fPair_Lo) { fPair_Lo->Set_v(x); } else { fIntercept[0]=x; } }
  if (plane==1) { if (fPair_Lo) { fPair_Lo->Set_u(x); } else { fIntercept[1]=x; } }
  
  if (plane==2) { if (fPair_Hi) { fPair_Hi->Set_v(x); } else { fIntercept[2]=x; } }
  if (plane==3) { if (fPair_Hi) { fPair_Hi->Set_u(x); } else { fIntercept[3]=x; } }
  
  UpdateTrackInfo(); 
}
double TvdcTrack::Intercept(int plane) const { 
  
  return fIntercept[plane]; 
}
int    TvdcTrack::Nhits(int plane) const { 
  
  if (!f_hasVDCdata) { 
    THRS::ReportError("TvdcTrack::Nhits(int)",
		      "Non-vdc track does not have TDC data.");
    return -999; 
  }
  
  if (!f_hasVDCdata || !fGroup[plane]) { 
    THRS::ReportError("TvdcTrack::Nhits(int)","Pointer for group invalid.");
    return -999; 
  }
  return fGroup[plane]->Nhits(); 
}
double TvdcTrack::Tau(int plane, int h) const {  
  
  if (!f_hasVDCdata) { 
    THRS::ReportError("TvdcTrack::Tau(int,int)",
		      "Non-vdc track does not have TDC data.");
    return -1e30; 
  }
  
  if (!fGroup[plane]) { 
    THRS::ReportError("TvdcTrack::Tau(int,int)","Pointer for group invalid.");
    return -1e30; 
  }
  
  if (h >= Nhits(plane) || h<0 ) { 
    THRS::ReportError("TvdcTrack::Tau(int,int)","Out-of-range hit requested."); 
    return -1e30; 
  }
  
  return fGroup[plane]->Time(h); 
}
double TvdcTrack::WirePos(int plane, int h) const {  
  
  if (!f_hasVDCdata) { 
    THRS::ReportError("TvdcTrack::WirePos(int,int)",
		      "Non-vdc track does not have TDC data.");
    return -1e30; 
  }
  
  if (!f_hasVDCdata || !fGroup[plane]) { 
    THRS::ReportError("TvdcTrack::WirePos(int,int)","Pointer for group invalid.");
    return -1e30; 
  }
  
  if (h >= Nhits(plane) || h<0 ) { 
    THRS::ReportError("TvdcTrack::WirePos(int,int)","Out-of-range hit requested."); 
    return -1e30; 
  }
  
  return fGroup[plane]->WirePos(h); 
}
double TvdcTrack::Slope(int plane)          const { 
  
  if (plane==0 || plane==2) { return fSlope_v; }
  if (plane==1 || plane==3) { return fSlope_u; }
  
  THRS::ReportError("TvdcTrack::Slope(int)","Invalid plane-num given."); 
  return -1e30; 
}

double ARRAY_rotate_uvw_to_xyz[9] = {  0.500000,  0.500000, -0.707107, 
				      -0.707107,  0.707107,  0.,
				       0.500000,  0.500000,  0.707107 }; 

double ARRAY_rotate_xyz_to_uvw[9] = {  0.500000, -0.707107,  0.500000, 
				       0.500000,  0.707107,  0.500000,
				      -0.707107,  0.,        0.707107 }; 


TVector3  TvdcTrack::Rotate_uvw_to_xyz( const TVector3 vec ) { 
  
  double ARRAY_vec[3] = { vec[0], vec[1], vec[2] }; 
  
  TVectorD r(3, ARRAY_vec); 
  
  TMatrixD U(3,3, ARRAY_rotate_uvw_to_xyz); 
  
  r = U*r; 

  return TVector3( r[0], r[1], r[2] ); 
}
TVector3  TvdcTrack::Rotate_xyz_to_uvw( const TVector3 vec ) { 
  
  double ARRAY_vec[3] = { vec[0], vec[1], vec[2] }; 
  
  TVectorD r(3, ARRAY_vec); 
  
  TMatrixD U(3,3, ARRAY_rotate_xyz_to_uvw); 
  
  r = U*r; 

  return TVector3( r[0], r[1], r[2] ); 
}
void TvdcTrack::Compute_Theta_Phi( const bool arm,  
				   const double intercepts[4], 
				   double &Theta, 
				   double &Phi ) { //static method
  
  double wSep = arm 
    ? VDC_w_R[2]-VDC_w_R[0] 
    : VDC_w_L[2]-VDC_w_L[0] ; 
  
  //compute slopes
  double m_v = wSep/(intercepts[2]-intercepts[0]); 
  double m_u = wSep/(intercepts[3]-intercepts[1]); 

  TVector3 S_xyz = TvdcTrack::Rotate_uvw_to_xyz( TVector3( 1./m_u, 1./m_v, 1. ) ); 
  
  Theta = TMath::ATan( S_xyz.X()/S_xyz.Z() ); 
  Phi   = TMath::ATan( S_xyz.Y()/S_xyz.Z() ); 
  
}  
TVector3 TvdcTrack::ComputeIntercept_w(const double w) const {
    
  return TVector3( fIntercept[1] + ( w - fW[1] )/fSlope_u, 
		   fIntercept[0] + ( w - fW[0] )/fSlope_v, w ); 
}
TVector3 TvdcTrack::ComputeIntercept_z(const double z) const {
  
  //get two arbitrary ponints on the track so we can interpolate to find 'z'
  TVector3 S  = TvdcTrack::Rotate_uvw_to_xyz( TVector3( 1./fSlope_u, 1./fSlope_v, 1. ) );  
  
  TVector3 r0 = TvdcTrack::Rotate_uvw_to_xyz( ComputeIntercept_w(0.) ); 
  
  double t    = ( z - r0.Z() )/S.Z();   
  
  return TVector3( r0  +  S*t ); 
}
double TvdcTrack::ToF(int plane, double x) const { 
  
  //U-plane (v-coord)
  if (plane==0 || plane==2) { return (fS2_v - x)/fC_v; }
  //V-plane (u-coord)
  else                      { return (fS2_u - x)/fC_u; }
}
double TvdcTrack::GetTimeAtZ(const double z) const {
  
  //computes the time which this track intercets some z-plane
  // (z in transport coords)
  // note that this time is given relative to the time of the S2-hit used to make
  // this track.
  
  //intercept with this plane, in uvw (vdc) coordinates
  TVector3 intercept_uvw
    = TvdcTrack::Rotate_xyz_to_uvw( ComputeIntercept_z(z) ); 
  
  double u = intercept_uvw[0];
  
  return (u - fS2_u)/fC_u;
}
double TvdcTrack::Get_T_model(int plane, double v, int derivative) const { 
  
  double slope = Slope(plane); 
  
  double w = slope*( v - Intercept(plane) ); 
  
  double Tau_model = fEvent->Drift_T( w,slope,derivative ); 
  
  //time-of-flight does not contribute to the derivatives of T
  if (derivative==0) Tau_model += -ToF( plane, v );   
  
  return Tau_model; 
}
double TvdcTrack::xParam() const { 
  
  return  
    ( f_S2Int_xyz.X() - fEvent->GetS2Hit()->X() ) 
    / fEvent->GetS2Hit()->PaddleWidth();     
}
void TvdcTrack::Nudge_params( double nudge[5] ) { 
  
  double u,v; 
  
  if ( !fPair_Lo || !fPair_Hi || !f_hasVDCdata ) { //is this a vdc-data track?
    
    for (int p=0; p<4; p++) fIntercept[p] += nudge[p]; 
    
  } else                                         { 
    
    fPair_Lo->Get_uv( u,v ); 
    fPair_Lo->Set_uv( u + nudge[1], 
		      v + nudge[0] );
    
    fPair_Hi->Get_uv( u,v ); 
    fPair_Hi->Set_uv( u + nudge[3], 
		      v + nudge[2] ); 
  } 
   
  fT0 += nudge[4]; 

  UpdateTrackInfo(); 
}
void   TvdcTrack::Set_params( const double vLo, 
			      const double uLo,
			      const double vHi, 
			      const double uHi, 
			      const double T0 ) { 
  
  if (fPair_Lo && fPair_Hi && f_hasVDCdata) {  
    
    fPair_Lo->Set_uv( uLo,vLo );
    fPair_Hi->Set_uv( uHi,vHi ); 
  
  } else  { 
    
    fIntercept[0] =vLo;
    fIntercept[1] =uLo;
    fIntercept[2] =vHi;
    fIntercept[3] =uHi;
  }
  
  fT0 = T0;  

  UpdateTrackInfo(); 
}
void   TvdcTrack::Set_params( const double params[5] ) { 
  
  Set_params( params[0], 
	      params[1], 
	      params[2], 
	      params[3], 
	      params[4] ); 

}
double TvdcTrack::FP_x()  const { return f_FPInt_xyz.X(); }
double TvdcTrack::FP_y()  const { return f_FPInt_xyz.Y(); }

double TvdcTrack::S2_x()  const { return f_S2Int_xyz.X(); }
double TvdcTrack::S2_y()  const { return f_S2Int_xyz.Y(); }

double TvdcTrack::Theta() const { return fTheta; }
double TvdcTrack::Phi()   const { return fPhi; }
////////////////////////////////////////////////////////////////////////////////////


ClassImp(THRS); 

ClassImp(TS2Hit);

ClassImp(TvdcHit);

ClassImp(THitGroup);

ClassImp(THitCluster); 

ClassImp(TChamberPair); 

ClassImp(TEventHandler); 

ClassImp(TvdcTrack); 
