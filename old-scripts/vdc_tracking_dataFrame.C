#include "def_apex.h"
#include "function_lib.C"
#include <fstream>

//#include "class/TvdcHit.h"
//#include "class/TvdcHit.C"
#include "ApexHRS.h"
#include <functional>
#include <string>

class vdc_tracking_dataFrame;

using namespace std; 

using RVecD  = ROOT::VecOps::RVec<double>;
using RNode  = ROOT::RDF::RNode; 


//these are used for track error-estimation
double dP; 
double cum_P; 

const int err_nSamples = 25; 

double z[err_nSamples]; 

bool choose_goodOrBad=true;  //true==look at 'good' tracks; false==bad tracks 

double f_model_( double *x, double *par ) { 
  
  auto dummy_s2R = new TS2Hit(true, 0,200,200); 
  auto dummy_s2L = new TS2Hit(false,0,200,200); 
    
  auto dummy_event = new TEventHandler(false, 53., dummy_s2R, dummy_s2L); 
  
  return dummy_event->Drift_T(x[0], 1.4, 2 );
}  
//#endif 

#define IS_MONTECARLO true

void vdc_tracking_dataFrame(int maxEntries=5e4, 
			    bool isRightArm=true, //true=Right, false=Left 
			    int runNum=4712, 
			    TString path_outFile="",  //path to out file (if used)
			    bool is_monteCarlo=false, 
			    TString path_inFile="replay/out")
{
  const TString here = "vdc_tracking_singleArm"; 
  

  //check EPICS tree
  ROOT::EnableImplicitMT();
  
  ////////////////////////////////////////
  ////////////////////////////////////////
  //
  vector<string> branch_S2; 
  
  branch_S2.push_back( is_monteCarlo ? "R_s2_rt" : "R.s2.rt" );
  branch_S2.push_back( is_monteCarlo ? "R_s2_lt" : "R.s2.lt" ); 
  branch_S2.push_back( is_monteCarlo ? "L_s2_rt" : "L.s2.rt" ); 
  branch_S2.push_back( is_monteCarlo ? "L_s2_lt" : "L.s2.lt" );  
  
  vector<string> plane_name = { "u1", "v1", "u2", "v2" }; 
  
  vector<string> branch_rawtime; 
  vector<string> branch_wire; 
  
  for (int p=0; p<4; p++) { 
    
    string rawtime = is_monteCarlo 
      ? "L_vdc_"+plane_name[p]+"_rawtime" 
      : "L.vdc."+plane_name[p]+".rawtime"; 
    
    string wire    = is_monteCarlo 
      ? "L_vdc_"+plane_name[p]+"_wire" 
      : "L.vdc."+plane_name[p]+".wire"; 
    
    branch_rawtime.push_back( rawtime ); 
    branch_wire   .push_back( wire );   
  }
  
  cout << "searching for both-arm coincidence peak..." << flush; 
  
  path_inFile += "."+TString::Itoa(runNum, 10)+".root";
  
  
  /***
   *    Takes Raw S2-m data from the decoded-data tree, and converts them to 
   *     the objects called 'TS2Hit*', stored in an RVec.  
  ***/ 
  auto Generate_S2_hits = []( bool arm, RVecD PMT_R, RVecD PMT_L ) { 
    
    //generate a vector of all coinc s2-paddle hits
    ROOT::RVec<TS2Hit*> coincHits; 
  
    for (int p=0; p<16; p++) { 
    
      auto hit = new TS2Hit( arm, p, PMT_R[p], PMT_L[p] );
      
      if ( hit->IsCoinc() ) { coincHits.push_back( hit ); }
      else                  { hit->~TS2Hit();             } 
    }
    
    return coincHits; 
  }; 
  
  
  /*** 
   *    Performs a similar function to the above, but compares the 
   *     T_R - T_L times for all coinc hits, so that the 
   *     variables tR_tL_diff/sigma (below) can be filled; 
   *     these variables are the basis for making both-arm S2 coinc-cuts 
   *     for all events during the full analysis.
   ***/ 
  auto get_dt = [&Generate_S2_hits]( RVecD RHRS_rt, RVecD RHRS_lt,
				     RVecD LHRS_rt, RVecD LHRS_lt ) { 
    
    //compare all s2m-coinc hits from both arms
    auto hits_RHRS = Generate_S2_hits( true,  RHRS_rt, RHRS_lt ); 
    auto hits_LHRS = Generate_S2_hits( false, LHRS_rt, LHRS_lt ); 
  
    RVecD dt; 
    
    for (int iR=0; iR<hits_RHRS.size(); iR++) { 
      for (int iL=0; iL<hits_LHRS.size(); iL++) { 
      	
	dt.push_back(  hits_RHRS.at(iR)->Time() 
		      -hits_LHRS.at(iL)->Time()  ); 
      }
    }
    
    return dt; 
  };
  
  
  //find the S2-coincidence peak
  ROOT::RDataFrame d_T("T", path_inFile.Data() ); 
  
  auto h_dt_temp = d_T
    .Range(0,5e4)
    .Define("dt", get_dt, branch_S2)
    .Histo1D({"", "dt", 200, -60e-9, 60e-9}, {"dt"}); 
  
  //h_dt_temp->DrawCopy(); 
  
  TVector2 dt_fitResult = draw_gausFit( (TH1D*)h_dt_temp->Clone(), 
					6e-9, -1e30, -1 ); 
  
  const double tR_tL_diff  = dt_fitResult.X(); 
  const double tR_tL_sigma = dt_fitResult.Y(); 
  
  const double CUT_S2Sep_stdDev = 6; 
  
  
  cout << "done." << endl; 
  
  cout << TString::Format("T_R - T_L offset,sigma = %2.2f ns, %2.3f ns", 
			  tR_tL_diff*1e9, tR_tL_sigma*1e9 ) << endl; 
  //auto get_RL_coincHits = [&Generate_S2_hits,
  
  
  /*** 
   *     Create a TEventHandler object for each both-arm S2-coinc, 
   *      using the variables above (tR_tL_diff/sigma) to make cuts
   ***/
  auto Gen_S2_coincHits 
    = [tR_tL_diff, tR_tL_sigma, CUT_S2Sep_stdDev,
       momentum_R, momentum_L]
    ( double beamCurrent, double beamEnergy, 
      ROOT::RVec<TS2Hit*> hits_R,
      ROOT::RVec<TS2Hit*> hits_L ) {
    
    ROOT::RVec<TEventHandler*> coincEvents; 
    
    for (int iR=0; iR<hits_R.size(); iR++) { 
      for (int iL=0; iL<hits_L.size(); iL++) { 
	
	double TR_TL = hits_R.at(iR)->Time() - hits_L.at(iL)->Time(); 
	
	if ( TMath::Abs(TR_TL-tR_tL_diff)/tR_tL_sigma < CUT_S2Sep_stdDev ) 
	  
	  //this is a bona-fide both-arm S2-coinc; let's keep it
	  coincEvents.push_back( new TEventHandler(false,
						   beamEnergy,
						   momentum_L, momentum_R,
						   beamCurrent, 
						   hits_R.at(iR), 
						   hits_L.at(iL) ) );  
		
      }//for (int iL=0; iL<hits_L.size(); iL++) { 
    }
    
    return coincEvents;
  }; 
  
  
  //
  //     Create 'hit groups' for a particular plane/arm.      
  //
  
  auto  group_hits = [] ( TEventHandler *evt, //true = RHRS, false = LHRS
			  int   p, //the plane we're working with
			  RVecD h_wire, 
			  RVecD h_time ) {
    
    const unsigned int clust_minHits=1; //minimum number of hits in a group
    const int maxGap = 3;      //maximum empty wire-gap in a group
    
    ROOT::RVec<THitGroup*> groupVec; 
    
    int nHits = h_time.size(); 
    if (h_time.size() < clust_minHits) return groupVec;
    
    //cBeg->push_back((int)wire[0]);   
    auto group = new THitGroup(p); 
    
    int wPrev = -1; 
    
    //the clust_points vector stores each hit as a 2-vector;
    // the [0] element is the wirePos,
    // the [1] element is the T_TDC (corrected tdc time of that hit)
    
    for (int h=0; h<nHits; h++) { 
      
      //vdc timing cut, these times won't ever be useful for a coinc track 
      auto hit = new TvdcHit( p, 
			      h_wire[h], 
			      h_time[h],  evt ); 
      
      if (hit->Time() > VDC_max_realTime || 
	  hit->Time() < VDC_min_realTime ) {
	
	hit->~TvdcHit(); //throw out this hit 
	continue; 
      }
      
      //if wPrev =-1, then no wires have yet been added
      int gap = wPrev<0 ? 0 : hit->wNum() - wPrev  -1; 
            
      if (gap > maxGap) { //for adjacent wires, 
	
	//does this cluster already have enough hits? 
	if (group->Nhits() >= clust_minHits) { 
	
	  groupVec.push_back( group );  
	  
	} else { group->~THitGroup(); }
	   
	//otherwise, make a new group anyway.
	group = new THitGroup(p);  
	
      }//if (span/gap too big) 
      //if span is NOT too big, then add this hit to the current cluster
      
      group->AddHit( hit );   wPrev = hit->wNum(); 
      
    }//for (int i=0; h<nHtis; h++) 
    
    //the last hit is reached, and it's not over-span
    if (group->Nhits() >= clust_minHits)  { 
      
      groupVec.push_back( group ); 
      
    } else  { 
      
      group->~THitGroup(); 
    }
    
    return groupVec; 
  }; 
  
  
  
  ///
  //        
  //
  ///
  auto Check_interWire = []( TEventHandler *evt, 
			     TChamberPair *pair, 
			     double &Eta_u,
			     double &Eta_v, 
			     double sigma =20e-9 ) {
    
    //get the wire above & below this guess
    double u,v;    pair->Get_uv( u,v ); 
    
    double mu,mv;  pair->GetSlope_uv( mu,mv ); 
    
    double guessBuffer = 0.125*dWire; 
    
    for (int p=0; p<2; p++) { 
      
      bool is_UPlane = (p==0); 
      
      double guess_Lo, m;
      THitGroup *group = new THitGroup; 
      
      if (is_UPlane) { //U-plane
	
	m = mv; 
	guess_Lo = pair->ClosestWirePos_Lo( v ) - guessBuffer; 
	group    = (THitGroup*)pair->GetGroup_U();  
	
      } else         { //V-plane
	
	m = mu; 
	guess_Lo = pair->ClosestWirePos_Lo( u ) - guessBuffer;  
	group    = (THitGroup*)pair->GetGroup_V(); 
      }
      
      double guess_Hi = guess_Lo + dWire + 2.*guessBuffer; 
      
      double eta(0.); 
      
      //loop over all hits
      for (int h=0; h<group->Nhits(); h++) { 
	
	double x   = group->WirePos(h); 
	
	double tau_a = evt->Drift_T( m*(x-guess_Lo), m ); 
	double tau_b = evt->Drift_T( m*(x-guess_Hi), m ); 
	
	double tau_Lo = TMath::Min( tau_a, tau_b ); 
	double tau_Hi = TMath::Max( tau_a, tau_b ); 
		
	//cut hits that are too far away
	if ( tau_Lo > VDC_max_realTime ) continue; 
	
        double tau = group->Time(h); 
	
	//tau is below lowest guess
	if (tau < tau_Lo) { 
	  
	  eta += TMath::Exp( -0.5*TMath::Power((tau-tau_Lo)/sigma, 2) ); 
	  continue; 
	}
	
	//'goldilocks'-zone, between hi- & lo-predictions
	if (tau < tau_Hi) { 
	  
	  eta += 1.; 
	  continue; 
	}
	
	//tau is above highest guess
	eta += TMath::Exp( -0.5*TMath::Power((tau-tau_Hi)/sigma, 2) ); 
	
      }//for (int h=0; h<group->Nhits(); 
      
      if (is_UPlane) { Eta_v =eta; } else { Eta_u =eta; }
      
    }//for (int p=0; p<2; p++) 
    
  };
  

  
  /*** 
   *    Take a pair of U1 & U2 groups, and generate 'chamber pairs' 
   *     (TChamberPair*) which are effectivley the stubs of new track candidates. 
   ***/   
  const double overlapBuffer = dWire*2.0; 
  const double uFix = 0.026/TMath::Sqrt(2);
  const double CUT_minEta      = 1.0;              //min eta of one plane
  const double CUT_minEta_both = 3.9;    //min eta of both planes
    
  auto Gen_pairs_loChamber = [&Check_interWire, overlapBuffer, uFix, 
			      CUT_minEta, CUT_minEta_both] 
    ( TEventHandler *evt, 
      ROOT::RVec<THitGroup*> gVec_U, 
      ROOT::RVec<THitGroup*> gVec_V ) { 
    
    //constants
    bool is_RightArm = evt->ActiveArm(); 
        
    const double group_xParamCut  = 10.; 
    const double gridSearch_CoarseSpacing = 0.5*dWire; 
    
    const double vu_Lo_diff_min = -0.082;
    const double vu_Lo_diff_max =  0.057; 
    
    const double CUT_vu_diff_min = is_RightArm ? -83e-3 : -85e-3; 
    const double CUT_vu_diff_max = is_RightArm ?  58e-3 :  58e-3; 
    
    const double mPhi_X = 0.1072; //a 'slope' which models the disperson of tracks
    // from the +z axis
    const double mPhi_Y = 0.0031; // same, but for the y-direction 
    
    const double Phi_X_Error = 0.020; 
    const double Phi_Y_Error = 0.011; 
    
    
    auto S2Hit = (TS2Hit*)evt->GetS2Hit(); 
    
    //for each possible pairing of groups, check to see if they
    // might be allowed to form a TChamberPair, and, if so, go ahead and generate it. 
    ROOT::RVec<TChamberPair*> pairs; 
    
    RVecD vTest; 
        
    //the 's-point' is an imaginary point far away from the plane of the S2 
    // at which the nominal track trajectories converge. 
    double S_point = S2Hit->Z() - 1./mPhi_X; 
    
    //loop over all groups in both planes
    for (int iU=0; iU<gVec_U.size(); iU++) { 
      for (int iV=0; iV<gVec_V.size(); iV++) { 
	
	//get the U & V groups
	auto gU = gVec_U.at(iU);
	auto gV = gVec_V.at(iV); 
	
	//check to make sure these clusters have an acceptable overlap
	// (this is a very tolerant cut, but is still useful!) 
	// u-v 
	double vu_diff_max = gU->HiEdge()-gV->LoEdge() + 2.*overlapBuffer; 
	// u-v 
	double vu_diff_min = gU->LoEdge()-gV->HiEdge() - 2.*overlapBuffer; 
	
	if ( vu_diff_max < CUT_vu_diff_min ) continue; 
	
	if ( vu_diff_min > CUT_vu_diff_max ) continue; 
	
	//check to make sure that an OK s2 intercept is possible
	TVector3 r_xMax( gV->HiEdge() + overlapBuffer - uFix,
			 gU->HiEdge() + overlapBuffer,  gU->W() ); 
	
	//rotate the vector, so it's now in xyz-coords
	r_xMax = vec_uvw_xyz( r_xMax ); 
	
	//now approximate the s2 intercept of the cluster's most extreme point
	Double_t X_s2_max = 
	  r_xMax.X()*(S2Hit->Z() - S_point)/(r_xMax.Z() - S_point) 
	  + Phi_X_Error*(S2Hit->Z() - r_xMax.Z()); 
	
	if ( ( X_s2_max - S2Hit->X() )/S2Hit->PaddleWidth() < -group_xParamCut ) 
	  continue; 
	
	//now, compute the maximum y-intercept
	TVector3 r_xMin( gV->LoEdge() - overlapBuffer - uFix,
			 gU->LoEdge() - overlapBuffer,  gU->W() ); 
	
	//rotate the vec.
	r_xMin = vec_uvw_xyz( r_xMin ); 
	
	//now approximate the S2 intercept, like before
	Double_t X_s2_min = 
	  r_xMin.X()*(S2Hit->Z() - S_point)/(r_xMin.Z() - S_point)
	  - Phi_X_Error*(S2Hit->Z() - r_xMax.Z()); 
	
	if ( ( X_s2_min - S2Hit->X() )/S2Hit->PaddleWidth() >  group_xParamCut ) 
	  continue; 
	
	//now, check for potential tracks
		
	//how many wires over are we gonna look for tracks?
	const int wire_buffer = 2; 
	
	//find span of our groupings
	int span_u 
	  = (int)( gV->Span()/dWire ) + 2*wire_buffer + 1; 
	
	int span_v 
	  = (int)( gU->Span()/dWire ) + 2*wire_buffer + 1; 
	
	
	double u0 
	  = gV->LoEdge() - ((double)wire_buffer)*dWire + 0.5*dWire; 
	
	double v0 
	  = gU->LoEdge() - ((double)wire_buffer)*dWire + 0.5*dWire; 
		
	
	auto testPair =  new TChamberPair( true, u0,v0, gU,gV ); 
	
	double 
	  best_Eta_u(CUT_minEta), 
	  best_Eta_v(CUT_minEta); 
	
	double best_u, best_v; 
	
	double best_Eta(-999.); 
	
	//now, loop over all inter-wire positions within these groupings
	for (int iu=0; iu<span_u; iu++) { 
	  double u = u0 + ((double)iu)*dWire; 
	  
	  for (int iv=0; iv<span_v; iv++) { 
	    double v = v0 + ((double)iv)*dWire;
	    
	    if ( v-u < CUT_vu_diff_min ||
		 v-u > CUT_vu_diff_max ) continue; 
	    
	    
	    testPair->Set_uv( u, v ); 
	    
	    //make a guess for the slopes at this point
	    TVector3 r0( u-0.026/TMath::Sqrt(2), v, gU->W() ); 
	    r0 = vec_uvw_xyz( r0 ); 
	    
	    double Phi_x_guess =
	      mPhi_X * r0.X() * (S2Hit->Z() - S_point)/(r0.Z() - S_point); 
 
	    	    
	    float Phi_y_guess = mPhi_Y*r0.X(); 
	    
	    TVector3 S( Phi_x_guess, Phi_y_guess, 1.0 ); 
	    
	    S = vec_xyz_uvw( S ); 
	    
	    Double_t m_u = S.Z()/S.X(); 
	    Double_t m_v = S.Z()/S.Y(); 
	    
	    
	    testPair->SetSlope_uv( m_u, m_v ); 
	    
	    
	    double Eta_u, Eta_v; 
	    	    
	    Check_interWire( evt, testPair, Eta_u, Eta_v ); 
	    
	    //if ( Eta_u < 1e-3 || Eta_v < 1e-3 ) continue; 
	    //cout << Eta_u << " " << Eta_v << endl; 
	    	    
	    if ( Eta_u         < CUT_minEta ||
		 Eta_v         < CUT_minEta || 
		 Eta_u + Eta_v < CUT_minEta_both ) continue; 
	    
	    //is this inter-wire guess the best found so far? 
	    if ( Eta_u+Eta_v > best_Eta ) { 
	      
	      best_u = u; 
	      best_v = v;  
	      
	      best_Eta = Eta_u+Eta_v; 
	      
	      best_Eta_u = Eta_u; 
	      best_Eta_v = Eta_v; 
	    }
	    	    
	  }//for (int iv=0; 
	}//for (int iu=0; 
	
	
	if (best_Eta_u > CUT_minEta && 
	    best_Eta_v > CUT_minEta ) {
	  
	  testPair->Set_uv( best_u,best_v ); 
	  
	  pairs.push_back( testPair ); 
	  
	} else { testPair->~TChamberPair(); }  
		
      }//for (int iV=0; 
    }//for (int iU=0; 

    //return vTest; 
    return pairs; 
  };  
  
  
  auto Phi_model   = [](TvdcTrack *track) { 
    
    return 
    -0.231*TMath::Power(track->S2_y(),2) 
    +  0.2532*track->S2_y() 
    +  0.539e-3; 
  };
  auto Theta_model = [](TvdcTrack *track) { 
    
    return 0.109648*track->S2_x(); 
  };
    
  
  const double CUT_ph_min = -0.018; 
  const double CUT_ph_max =  0.008; 
  
  const double CUT_th_min = -0.0200; 
  const double CUT_th_max =  0.0175; 
  
  
  auto Gen_pairs_hiChamber = [&Check_interWire, 
			      overlapBuffer, 
			      uFix, CUT_minEta, CUT_minEta_both,
			      &Phi_model,   CUT_ph_min, CUT_ph_max,
			      &Theta_model, CUT_th_min, CUT_th_max] 
    ( TEventHandler *evt, 
      ROOT::RVec<TChamberPair*> pairs_Lo, 
      ROOT::RVec<THitGroup*> gVec_U, 
      ROOT::RVec<THitGroup*> gVec_V ) { 
    
    
    bool is_RightArm = evt->ActiveArm(); 
    
    const double CUT_xParam = 2.00; 
    
    const double CUT_vu_diff_min = is_RightArm ? -95e-3 : -96e-3;  
    const double CUT_vu_diff_max = is_RightArm ?  70e-3 :  70e-3; 
    
    ROOT::RVec<TvdcTrack*> tracks; 
            
    auto S2Hit = (TS2Hit*)evt->GetS2Hit(); 
    
    
    //first, check upper-chamber groups, to see which match (if any) 
    ROOT::RVec<TChamberPair*> pairs_Hi; 
    
    for (int iU=0; iU<gVec_U.size(); iU++) { 
      for (int iV=0; iV<gVec_V.size(); iV++) { 
	
	//get the U & V groups
	auto gU = gVec_U.at(iU);
	auto gV = gVec_V.at(iV); 
	
	//check to make sure these clusters have an acceptable overlap
	// (this is a very tolerant cut, but is still useful!) 
	// u-v 
	double vu_diff_max = gU->HiEdge()-gV->LoEdge() + 2.*overlapBuffer; 
	// u-v 
	double vu_diff_min = gU->LoEdge()-gV->HiEdge() - 2.*overlapBuffer; 
	
	if ( vu_diff_max < CUT_vu_diff_min ) continue;
	
	if ( vu_diff_min > CUT_vu_diff_max ) continue;  
	
	//acceptable agreement... proceed using this pair
	
	pairs_Hi.push_back( new TChamberPair( false, //false == hiChamber 
					      0.,0., 
					      gU,gV ) ); 
      }
    }
  
    //no possible pairs found in the upper chamber!! return the empty vec. 
    if (pairs_Hi.size()<1) return tracks; 
    
    
    for (int pH=0; pH<pairs_Hi.size(); pH++) { 
      for (int pL=0; pL<pairs_Lo.size(); pL++) { 

	auto pHi = pairs_Hi.at(pH); 
	auto pLo = pairs_Lo.at(pL); 
	
	auto track = new TvdcTrack( evt, pLo, pHi ); 
	
	//check highest possible s2 X-intercept
	
	track->Set_uv_Hi( pHi->GetGroup_V()->HiEdge() + overlapBuffer, 
			  pHi->GetGroup_U()->HiEdge() + overlapBuffer ); 
	
	
	if (track->xParam() < -CUT_xParam) { 
	  
	  track->~TvdcTrack(); continue; 
	}
	
	
	//check lowest possible s2 X-intercept
	
	track->Set_uv_Hi( pHi->GetGroup_V()->LoEdge() - overlapBuffer, 
			  pHi->GetGroup_U()->LoEdge() - overlapBuffer ); 
	
	//cout << "min = " << track->xParam() << endl; 
	
	if (track->xParam() > CUT_xParam) { 
	  
	  track->~TvdcTrack(); continue; 
	}
		
	//now, we can finally scan for new inter-wire positions in the upper plane
	
	const int wire_buffer = 2; 
	
	//find span of our groupings
	int span_u 
	  = (int)( pHi->GetGroup_V()->Span()/dWire ) + 2*wire_buffer + 1; 
	
	int span_v 
	  = (int)( pHi->GetGroup_U()->Span()/dWire ) + 2*wire_buffer + 1; 
		
	double u0 
	  = pHi->GetGroup_V()->LoEdge() - ((double)wire_buffer)*dWire + 0.5*dWire; 
	
	double v0 
	  = pHi->GetGroup_U()->LoEdge() - ((double)wire_buffer)*dWire + 0.5*dWire; 
	
	double 
	  best_Eta_u(-999), 
	  best_Eta_v(-999); 
	
	double best_u, best_v; 
	
	double best_Eta(0.); 
		
	//now, loop over all inter-wire positions within these groupings
	for (int iu=0; iu<span_u; iu++) { 
	  double u = u0 + ((double)iu)*dWire; 
	  
	  for (int iv=0; iv<span_v; iv++) { 
	    double v = v0 + ((double)iv)*dWire;
	    
	    //check the intra-plane agreement of this guess
	    if ( v-u < CUT_vu_diff_min ||
		 v-u > CUT_vu_diff_max ) continue; 
	    	    
	    const double error_theta = 8e-3; 
	    const double error_phi   = 12e-3; 
	    
	    track->Set_uv_Hi( u,v ); 
	    	    
	    //check to make sure this pair has a plausible angle
	    double modelErr; 
	    
	    modelErr = track->Theta() - Theta_model( track ); 
	    
	    if (modelErr-error_theta > CUT_th_max || 
		modelErr+error_theta < CUT_th_min  ) continue; 
	    
	    modelErr = track->Phi()   - Phi_model( track ); 
	    
	    if (modelErr-error_phi   > CUT_ph_max || 
		modelErr+error_phi   < CUT_ph_min  ) continue; 
	    
	    
	    double Eta_u, Eta_v; 
	    
	    Check_interWire( evt, track->GetPair_Hi(), Eta_u, Eta_v ); 
	    	    
	    if ( Eta_u         < CUT_minEta ||
		 Eta_v         < CUT_minEta || 
		 Eta_u + Eta_v < CUT_minEta_both ) continue; 
	    
	    //is this inter-wire guess the best found so far? 
	    if ( Eta_u+Eta_v > best_Eta ) { 
	      
	      best_u = u; 
	      best_v = v;  
	      
	      best_Eta = Eta_u+Eta_v; 
	      
	      best_Eta_u = Eta_u; 
	      best_Eta_v = Eta_v; 
	    }
	    
	  }//for (int iv=0; 
	}//for (int iu=0; 
	
	
	if (best_Eta_u > CUT_minEta && 
	    best_Eta_v > CUT_minEta ) {
	  
	  track->Set_uv_Hi( best_u,best_v ); 
	  
	  tracks.push_back( track ); 
	  
	} else { 
	  
	  track->~TvdcTrack();
	}  
	
	
      }//for (int pH=0; pH<pairs_Hi.size(); 
    }//for (int pL=0; pL<pairs_Lo.size(); 
       
    return tracks; 
  }; 


  const double d_Theta_d_Xfp = 0.16745; 
  
  //
  // refine tracks from the hi-chamber using newton's method
  //
  //   
  auto Refine_track = [CUT_th_max, CUT_th_min,
		       CUT_ph_max, CUT_ph_min]( TvdcTrack *trk, 
						const int nCycles=10, 
						double sigma=40e-9 ) { 
    
    //random-walk track minimizaiton
    const double GRAD_momentum = 0.50; 
    const double GRAD_step0    = 0.05; 
    const double GRAD_exponent = 0.33; 
    
    //cout << "refining track..." << endl; 
    
    const double scale_T = 3e-9; 
    const double scale_X = 2.2e-3; 
    
    const double nudge_multiplier = 1.; 
    
    const double Chi_cutoff = 6.; 
    
    const double sigma_decay = 0.9;
    const double min_sigma   = sigma*2.; 
    
    ///////////////////////////////////////////////////////////////////////////////
    double int0[5]; 
    for (int p=0; p<4; p++) int0[p] = trk->Intercept(p); 
    TString plane_name[4] = {"U1","V1","U2","V2"}; 
    int0[4] = trk->T0(); 
    //*/////////////////////////////////////////////////////////////////////////////
    
    const double s2 = sigma*sigma; 
    
    ROOT::RVec<double> tVec; 
    
    double objective_eta0(0.); 
    
    for (UInt_t c=0; c<nCycles; c++) { 
      
      //this measures if 'eta' is actually increasing
      double objective_eta1(0.); 
      
      //this will be how we 'nudge' the parameters
      double nudge[5] = {0.}; 
      
      double Deriv[5] = {0.}; 
      
      /*///////////////////////////////////////////////////// DRAW MINIMIZE
      auto canv = new TCanvas("ePic", "Event draw", 0,0, 900, 500); 
      canv->Divide( 2,2 ); 
      
      gStyle->SetOptStat(0); 
      int iCanvas=1; 
      //*////////////////////////////////////////////////////////////////////////
      
      //first-derivatives
      double F[5] = {0.}; 
      
      double J_ii[4] = {0.}; 
      double J_4i[4] = {0.}; 
      double J_44(0.); 
      
      for (int p=0; p<4; p++) { 
	
	/*/draw event picture ///////////////////////////////////////////////////
	canv->cd( iCanvas );  iCanvas++; 
	auto h2d = new TH2D(TString::Format("vdcPlane_%i",p), "", 
			    200, -25e-3, 25e-3, 
			    200, VDC_min_realTime, VDC_max_realTime); 
	//*//////////////////////////////////////////////////////////////////////
	
	auto group = trk->GetGroup(p); 
	
	for (UInt_t h=0; h<group->Nhits(); h++) { 
	  	  
	  /*/draw event picture ///////////////////////////////////////////////
	  h2d->Fill( group->WirePos(h)-int0[p], group->Time(h) ); 	  
	  //*/////////////////////////////////////////////////////////////////
	  
	  double Chi 
	    = trk->Get_T_model(p, group->WirePos(h))  
	    +  trk->T0() 
	    -  group->Time(h); 
	  
	  if ( TMath::Abs(Chi) > Chi_cutoff*sigma ) continue; 
	  
	  double Eta = TMath::Exp( -0.5*TMath::Power(Chi/sigma,2) ); 
	  
	  double m = trk->Slope(p); 
	  
	  F[p] += -m * Eta * Chi * trk->Get_T_model(p,group->WirePos(h),1); 
	  
	  F[4] +=  Eta * Chi; 
	
	  J_ii[p] 
	    += (Chi * trk->Get_T_model(p,group->WirePos(h),2)  
		+ TMath::Power(trk->Get_T_model(p,group->WirePos(h),1),2))*m*m*Eta; 
	  
	  J_4i[p] 
	    +=  -Eta * m * trk->Get_T_model(p,group->WirePos(h),1); 
	  
	  J_44 += Eta; 
	  
	  objective_eta1 += TMath::Exp( -0.5*TMath::Power(Chi/min_sigma,2) ); 
	  
	}//for (UInt_t h=0; h<trk->pGet_N(p); 
	
	/*/////////////////////////////////////////////////////////////////////////
	h2d->SetMarkerStyle(kOpenCircle); 
	h2d->SetMarkerSize(0.50); 
	h2d->DrawCopy(); 
	
	auto ltx = new TLatex; 
	
	ltx->DrawLatex( -23e-3, VDC_max_realTime, 
			plane_name[p] ); 
    	
	const int nPts = 100; 
	
	double max_X = 14e-3; 
	double X[nPts]; 
	double T[nPts]; 
	
	double dx = 2.*max_X/((double)nPts-1); 
	
	double xx=(trk->Intercept(p)-int0[p])-max_X; 
	for (int i=0; i<nPts; i++) { 
	  X[i] = xx;  xx += dx; 
	  T[i] = trk->Get_T_model( p, xx+int0[p] ) + trk->T0(); 
	}	
	auto g = new TGraph(nPts, X, T); 
	g->SetLineColor(kRed); 
	g->Draw("SAME"); 
	
	h2d->~TH2D(); 
	ltx->~TLatex(); 
	//*////////////////////////////////////////////////////////////////////////
		
      }//for (int p=0; p<4; p++) 
      
      /*///////////////////////////////////////////////////////////////////////////
      canv->Print("~/ftp-dump/plots/test_track_hi_Dt.gif+33");
      canv->~TCanvas(); 
      //*//////////////////////////////////////////////////////////////////////////
      
      /*cout << TString::Format("eta-obj = %7.7f (delta = %7.7e)", 
	objective_eta1,objective_eta1-objective_eta0); */ 
      
      //if (objective_eta1 < objective_eta0) break; 
      objective_eta0 = objective_eta1; 
          
      sigma *= sigma_decay; 
      
      
      for (int p=0; p<4; p++) { 
	
	J_44  +=  J_4i[p]*(-J_4i[p]/J_ii[p]); 
	F[4]  +=  F[p]   *(-J_4i[p]/J_ii[p]); 
      }
      
      nudge[4] = -nudge_multiplier*F[4]/J_44; 
      
      for (int p=0; p<4; p++) { 
	
	nudge[p] = -nudge_multiplier*(F[p] - J_4i[p]*nudge[4])/J_ii[p]; 
      }
      
      trk->Nudge_params( nudge ); 
      
    }
    
    trk->SetEta( objective_eta0 ); 
  };   
  
  
  
  
  dP = 1./((double)err_nSamples); 
  
  cum_P = dP; 
  
  for (int i=0; i<err_nSamples; i++) { 
    
    z[i] = ROOT::Math::normal_quantile( cum_P-0.5*dP, 1. );  
    cum_P += dP;
  }
  
    
  auto Compute_trackError = [](TvdcTrack *trk, 
			       double ERR_tau_sigma =9e-9) {
    
    double param_error[5]; unsigned int total_pts(0); 
    
    for (int p=0; p<4; p++) { 
      
      double S0(0.), S1(0.); 
    
      auto group = trk->GetGroup(p); 
      
      for (int h=0; h<group->Nhits(); h++) { 
      	
	double err 
	  = trk->Get_T_model( p, group->WirePos(h) )
	  + trk->T0()  
	  - group->Time(h); 
	
	if ( TMath::Abs(err)/ERR_tau_sigma > 4. ) continue; 
	
	total_pts++; 
	
	S0 += TMath::Power( trk->Get_T_model(p,group->WirePos(h),1), 2 ); 
	
	S1 
	  += trk->Get_T_model(p,group->WirePos(h),1) 
	  *  trk->Get_T_model(p,group->WirePos(h),2); 
      }
      S0 = TMath::Sqrt(S0); 
      
      double v_square(0.); 
      //sample all possible values of i
      for (int i=0; i<err_nSamples; i++) { 
	
	//cout << z[i] << " " << endl; 
	
	double v2 = ERR_tau_sigma*z[i]*S0 / ( ERR_tau_sigma*z[i]*S1/S0 + S0*S0 ); 
	
	v_square += TMath::Power( v2/trk->Slope(p), 2) * dP ; 
      }
      
      param_error[p] = TMath::Sqrt(v_square); 
    
      //cout << TString::Format("  -  %3.3f", param_error[p]*1e3) << flush; 
      
    }// for (int p=0; p<4; p++) 
    
    //error of T0 
    param_error[4] = ERR_tau_sigma/TMath::Sqrt( (double)total_pts ); 
    
    trk->Set_Errors( param_error ); 
  }; 
  
  
  const int CUT_plane_goodPts = 2; 
  
  //now actually perform the analysis chain
  auto nEvents_coinc = d_T
    .Range(0,maxEntries)

#if IS_MONTECARLO 
    //MONTE-CARLO  -- check how many 'good pts' per plane
    .Filter([CUT_plane_goodPts](ROOT::RVec<int> nPts_good, 
				ROOT::RVec<int> nPts_all) { 
	
	for (int p=0; p<4; p++) { 
	  if (nPts_good.at(p) < CUT_plane_goodPts) 
	    return false; 
	  
	  if (nPts_all.at(p)-nPts_good.at(p) < CUT_plane_goodPts) 
	    return false; 
	}
	return true; }, {"MONTE_nGoodHits", "MONTE_nAllHits"}) 
#endif 
    
    //find all valid S2 hits for the Rigth arm
    .Define("S2_hits_RHRS", [&Generate_S2_hits](RVecD rt, RVecD lt)
	    { return Generate_S2_hits(true,  rt,lt); }, 
	    {branch_S2[0].data(), branch_S2[1].data()})
    
    //find all valid S2 hits for the Left arm
    .Define("S2_hits_LHRS", [&Generate_S2_hits](RVecD rt, RVecD lt)
	    { return Generate_S2_hits(false, rt,lt); }, 
	    {branch_S2[2].data(), branch_S2[3].data()})
    
    //Find all possible both-arm S2 coincidences
    .Define("coincEventVec", Gen_S2_coincHits, 
	    {"hac_bcm_average","S2_hits_RHRS","S2_hits_LHRS"})
    
    //Make a cut on Coincs, we only allow exactly 1-per-event
    .Filter([](ROOT::RVec<TEventHandler*> v)
	    { return v.size()==1; },
	    {"coincEventVec"})
    
    //get the event handler (so we don't have to keep passing it as a vec)
    .Define("event", [](ROOT::RVec<TEventHandler*> v) 
	    { return v.at(0); }, 
	    {"coincEventVec"}); 
    
 
  vector<string> branch_group; 
  
  //Now that we have eliminated non-coinc events, we can proceed with the 
  // analysis. (Starting with the right arm)
  
  for (int p=0; p<4; p++) { 
    
    branch_group.push_back( (string)"groups_L_"+plane_name[p] ); 
    
    nEvents_coinc = nEvents_coinc
      
      .Define(branch_group[p].data(), [&group_hits,p](TEventHandler *evt, 
						      RVecD wire, RVecD time)
	      { return group_hits(evt,p,wire,time); },
	      {"event", branch_wire[p].data(), branch_rawtime[p].data()});
    
  }
  
  
  auto nEvents_1group = nEvents_coinc 
    .Filter([](ROOT::RVec<THitGroup*> v)
	    {return v.size()>0;}, {branch_group[0].data()})
    .Filter([](ROOT::RVec<THitGroup*> v)
	    {return v.size()>0;}, {branch_group[1].data()})
    .Filter([](ROOT::RVec<THitGroup*> v)
	    {return v.size()>0;}, {branch_group[2].data()})
    .Filter([](ROOT::RVec<THitGroup*> v)
	    {return v.size()>0;}, {branch_group[3].data()}); 
  
  
  
  /*auto h1d_test = nEvents_1group
    
    .Define("hit_sep", [](ROOT::RVec<THitGroup*> gVec) { 
	
	RVecD v; 
	
	for (int g=0; g<gVec.size(); g++) { 
	  
	  auto group = gVec.at(g); 
	  
	  int lastWire = -1; double lastTime; 
	  for (int h=0; h<group->Nhits(); h++) { 
	    
	    if (group->WireNum(h) == lastWire) { 
	      
	      v.push_back( (group->Time(h) - lastTime)*1e9 ); 
	    }   
	    lastTime = group->Time(h); 
	    lastWire = group->WireNum(h); 
	  }	  
	}
	return v; 
      }, {"groups_L_u1"})
    
    .Histo1D({"", "T_{i} - T_{i+1} (ns)", 200, -100, 100}, "hit_sep"); 
  

  h1d_test->DrawCopy(); 
  return; //*/ 
  
 
  
  auto nEvents_1pair = nEvents_1group
    
    //create lo-chamber pairs
    .Define("pairs_L_loChamber", Gen_pairs_loChamber, 
	    { "event", 
		branch_group[0].data(), 
		branch_group[1].data() })
    
    .Filter([](ROOT::RVec<TChamberPair*> v)
	    { return v.size()>0; }, {"pairs_L_loChamber"}); 
  
  
  //'refine' the tracks; 
  const double CUT_track_minEta = 6.5; //chosen arbitrarily, will change later
    
  auto nEvents_1track = nEvents_1pair 
    
    //generate raw tracks
    .Define("tracks_L_raw", Gen_pairs_hiChamber, 
	    { "event",
	      "pairs_L_loChamber", 
	      branch_group[2].data(), 
	      branch_group[3].data() })
    
    .Filter([](ROOT::RVec<TvdcTrack*> v) 
	    { return v.size()>0; }, {"tracks_L_raw"})
    
    //refine tracks
    .DefineSlotEntry("tracks_L", [CUT_track_minEta, 
				  &Refine_track, &Compute_trackError]
		     (unsigned int RDF_slot,
		      ULong64_t    RDF_entry, 
		      ROOT::RVec<TvdcTrack*> tracks) { 
		       
	      double best_eta(-1.); int best_track=0; 
	      
	      //if (RDF_entry != 111) return tracks; 
	      
	      for (int t=0; t<tracks.size();) { 
		
		Refine_track( tracks.at(t) ); 
		
		double eta = tracks.at(t)->GetEta(); 
		
		//delete this track
		if ( eta < CUT_track_minEta ) { 
		  
		  tracks.at(t)->~TvdcTrack(); 
		  tracks.erase( tracks.begin()+t ); 
		  continue; 
		}
				
		Compute_trackError( tracks.at(t) ); 
				
		//keep only best track
		if (eta > best_eta) { 
		  best_track = t; 
		  best_eta   = eta; 
		}
		 		
		t++; 
	      }
	      
	      /*/this keeps only the best track
	      ROOT::RVec<TvdcTrack*> tVec; 
	      
	      tVec.push_back( tracks.at(best_track) ); 
	      	      
	      return tVec; //*/ 
	      return tracks; }, {"tracks_L_raw"}); 
    
  
  
#if IS_MONTECARLO 
    
  ///MONTE-CARLO
  const double CUT_maxTrackErr    = 0.8e-3; 
  const double CUT_trackOffsetErr = 20e-9; 
  
  auto nEvents_trueTrack = nEvents_1track 
    
    .Define("nGoodTracks", [CUT_maxTrackErr, CUT_trackOffsetErr]
	    (ROOT::RVec<TvdcTrack*> tracks, 
	     RVecD good_track, 
	     RVecD accident_track) { 
	      
	      int nGoodTracks(0); 
	      
	      for (int t=0; t<tracks.size(); t++) { 
				
		auto trk = tracks.at(t); 
		
		bool match_goodTrack=true; 
		bool match_accidentTrack=true; 
				
		for (int p=0; p<4; p++) { 
		  
		  double Err;
		  
		  Err = TMath::Abs(trk->Intercept(p) - good_track[p]); 
		  if (Err > CUT_maxTrackErr) match_goodTrack=false; 
		  
		  Err = TMath::Abs(trk->Intercept(p) - accident_track[p]); 
		  if (Err > CUT_maxTrackErr) match_accidentTrack=false; 
		}
		
		/*double T_err; 
		
		T_err = TMath::Abs(trk->T0()-good_track[4]); 
		if (T_err > CUT_trackOffsetErr) match_goodTrack=false; 
		
		T_err = TMath::Abs(trk->T0()-accident_track[4]); 
		if (T_err > CUT_trackOffsetErr) match_accidentTrack=false; */ 
				
		if (match_goodTrack || match_accidentTrack) {
		  trk->Set_isGoodTrack(true); 
		} else { 
		  trk->Set_isGoodTrack(false); 
		}
		
		if (trk->IsGoodTrack()) nGoodTracks++; 
	      }
	      return nGoodTracks; 
	      
	    }, {"tracks_L", "MONTE_track_params", 
		"MONTE_track_params_accident"})
    
    //forcal plane intercepts
    .Define("tracks_s2_x", [](ROOT::RVec<TvdcTrack*> tracks, 
			      int nGoodTracks) { 
	      RVecD v; 
	      
	      //I need the column 'nGoodTracks' as an (irritating) workaround. 
	      //Basically, I need the tracks to 'remember' whether I set
	      //  them as 'good' or 'bad.' the 'volatile' keyword basically tells 
	      //  the compiler "hey, i know it looks like i dont actually USE 
	      //  'dummy_int', but it's important, so don't ignore it." 
	      
	      for (int t=0; t<tracks.size(); t++) 
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad)  
		  v.push_back(tracks.at(t)->S2_x()); 
	      
	      return v; }, {"tracks_L", "nGoodTracks"})
    
    .Define("tracks_s2_y", [](ROOT::RVec<TvdcTrack*> tracks, 
			      int nGoodTracks) { 
	      RVecD v; 
	      
	      for (int t=0; t<tracks.size(); t++)	  
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad) 
		  v.push_back(tracks.at(t)->S2_y()); 
	      
	      return v; 
	    }, {"tracks_L", "nGoodTracks"})
    
    .Define("tracks_theta", [](ROOT::RVec<TvdcTrack*> tracks, 
					       int nGoodTracks) { 
	      RVecD v; 
	      
	      for (int t=0; t<tracks.size(); t++)	  
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad) 
		  v.push_back(tracks.at(t)->Theta()); 
	      
	      return v; 
	    }, {"tracks_L", "nGoodTracks"})
    
    .Define("tracks_phi", [](ROOT::RVec<TvdcTrack*> tracks,  
			     int nGoodTracks) { 
	      RVecD v; 
	      
	      for (int t=0; t<tracks.size(); t++)	  
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad) 
		  v.push_back(tracks.at(t)->Phi()); 
	      
	      return v; 
	    }, {"tracks_L", "nGoodTracks"}); 
  
  
  
  
#endif
  
  
  //define some parameters to print
  auto h2d_Theta_s2x = nEvents_trueTrack 
    .Histo2D({"theta_s2x", "#Theta vs X_{S2};X_{S2} (m);#Theta (rad)", 
	  200, -1.2, 1.2, 200, -0.15, 0.15}, "tracks_s2_x", "tracks_theta"); 
    
  /*auto h2d_Theta_s2y = nEvents_1track 
    .Histo2D({"theta_s2y", "#Theta vs Y_{FP};Y_{FP} (m);#Theta (rad)", 
	  200, -0.2, 0.3, 200, -0.2, 0.2}, 
	  "tracks_s2_y", "tracks_theta"); //*/ 
  
  
  
  /*auto h2d_Phi_s2x = nEvents_1track 
    .Histo2D({"Phi_fpx", "#Phi vs X_{FP};X_{FP} (m);#Phi (rad)", 
	  200, -1.2, 1.0, 100, -0.05, 0.06}, 
	  "tracks_s2_x", "tracks_phi"); */ 
  
  auto h2d_Phi_s2y = nEvents_trueTrack 
    .Histo2D({"phi_s2y", "#Phi vs Y_{S2};Y_{S2} (m);#Phi (rad)", 
	  200, -0.325, 0.325, 200, -0.15, 0.15}, "tracks_s2_y", "tracks_phi"); //*/ 
  
  
  

  auto h2d_xParam_Dt = nEvents_trueTrack 
    
    .Define("tracks_xParam", [](ROOT::RVec<TvdcTrack*> tracks)  
	    { 
	      ROOT::RVec<double> xParam; 
	      
	      for (int t=0; t<tracks.size(); t++) { 
		
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad)
		  xParam.push_back( tracks.at(t)->xParam() ); 
	      }
	      return xParam; }, {"tracks_L"})
    
    .Define("tracks_Dt", [](ROOT::RVec<TvdcTrack*> tracks)  
	    { 
	      ROOT::RVec<double> xParam; 
	      
	      for (int t=0; t<tracks.size(); t++) { 
		
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad)
		  xParam.push_back( tracks.at(t)->T0() ); 
	      }
	      return xParam; }, {"tracks_L"})
    
    .Histo2D({"","track agreement with coinc-hit;xParam;T_{track}-T_{S2} (s)", 
	  200, -4.5, 4.5, 
	  200, -100e-9, 100e-9}, "tracks_xParam", "tracks_Dt"); 
  
  
  
  auto h2d_trackErr = nEvents_1track 
    
    .Define("err_actual", [](RVecD monte_track, 
			     ROOT::RVec<TvdcTrack*> tracks) 
	    { 
	      RVecD err; 
	      
	      for (int t=0; t<tracks.size(); t++) { auto trk = tracks.at(t); 
		
		for (int p=0; p<4; p++) { 
		  
		  err.push_back( monte_track[p] - trk->Intercept(p) ); 
		}
	      }
	      
	      return err; 
	    }, {"MONTE_track_params", "tracks_L"}) 
    
    .Define("err_guess", [](ROOT::RVec<TvdcTrack*> tracks) 
	    { 
	      RVecD err; 
	      
	      for (int t=0; t<tracks.size(); t++) { auto trk = tracks.at(t); 
		for (int p=0; p<4; p++) { 
		  
		  err.push_back( trk->GetError_intercept(p) ); 
		}
	      }
	      
	      return err; 
	    }, {"tracks_L"})
    
    .Histo2D({"", ";error - model (m);error - actual (m)", 
	  7, 0.13e-3, 0.20e-3, 45, -0.6e-3, 0.6e-3}, "err_guess","err_actual"); //*/ 
  
  //.Histo1D({"", ";error - model (m);", 200, 0, 0.6e-3}, "err_guess"); 
  
  
  
  auto h2d_w_tau = nEvents_1track 
    .Define("w", [](ROOT::RVec<TvdcTrack*> tracks) { 
	
	RVecD w; 
	
	for (int t=0; t<tracks.size(); t++) { 
	  
	  auto trk = tracks.at(t); 
	  
	  for (int p=0; p<4; p++) { 
	    
	    for (int h=0; h<trk->Nhits(p); h++) { 
	      
	      w.push_back( TMath::Abs( ( trk->WirePos(p,h)
					-trk->Intercept(p) )*trk->Slope(p) ) ); 
	    }
	  }
	}
	
	return w; 
      }, {"tracks_L"})
    
    .Define("tau", [](ROOT::RVec<TvdcTrack*> tracks) { 
	
	RVecD tau; 
	
	for (int t=0; t<tracks.size(); t++) { 
	  
	  auto trk = tracks.at(t); 
	  
	  for (int p=0; p<4; p++) { 
	    
	    for (int h=0; h<trk->Nhits(p); h++) { 
	      
	      tau.push_back( trk->Tau(p,h) ); 
	      /*- (trk->T0()+trk->Get_T_model(p,trk->WirePos(p,h)))*/
	    }
	  }
	}
	
	return tau; 
      }, {"tracks_L"})
    
    .Histo2D({"w_tau", ";drift-dist (m);drift-time (s)", 
	  200, 0, 20e-3, 200, -20e-9, 380e-9}, "w", "tau"); 
  
  
  auto h1d_eta = nEvents_trueTrack
    //define some parameters
    .Define("tracks_eta", [](ROOT::RVec<TvdcTrack*> tracks)  
	    { 
	      ROOT::RVec<double> eta;
	      
	      for (int t=0; t<tracks.size(); t++) { 
		
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad)
		  eta.push_back( tracks.at(t)->GetEta() ); 
	      }
	      return eta; 
	    }, {"tracks_L"})
    
    .Histo1D({"eta", "Eta distribution", 200, 0, 25}, "tracks_eta"); 
    
  
  auto h1d_dv = nEvents_1track 
    
    .Define("tracks_dv", [](ROOT::RVec<TvdcTrack*> tracks) 
	    {
	      ROOT::RVec<double> dv_vec; 
	      
	      for (int t=0; t<tracks.size(); t++) { 
		
		for (int p=0; p<4; p++) { 
		  
		  auto pair 
		    = p<2 ? tracks.at(t)->GetPair_Lo() 
		    : tracks.at(t)->GetPair_Hi(); 
		  		  
		  double dv 
		    = tracks.at(t)->Intercept(p) 
		    - pair->ClosestWirePos( tracks.at(t)->Intercept(p) ); 
		  
		  dv_vec.push_back( dv*1e3 ); 
		}
	      }
	      
	      return dv_vec; }, {"tracks_L"})
    
    .Histo1D({"h_dv", "dv-distribution (mm)", 125, -0.5e3*dWire, 0.5e3*dWire}, 
	     "tracks_dv"); //*/ 
  
  
  
  /*auto h1d_Theta = nEvents_1track
    .Define("theta_test", [&Phi_model](ROOT::RVec<TvdcTrack*> tracks) { 
	RVecD v; 
	
	double nudge[5] = { dWire, -dWire, -dWire, dWire, 0. }; 
	
	for (int i=0; i<tracks.size(); i++) {
	  
	  auto trk = tracks.at(i); 
	  
	  double old_theta = trk->Phi(); 
	  
	  double uHi,vHi; trk->GetPair_Hi()->Get_uv( uHi,vHi ); 
	  
	  double uLo,vLo; trk->GetPair_Lo()->Get_uv( uLo,vLo ); 
	  
	  uHi = trk->GetPair_Hi()->ClosestWirePos_Lo( uHi ); 
	  vHi = trk->GetPair_Hi()->ClosestWirePos_Lo( vHi ); 
	  
	  uLo = trk->GetPair_Lo()->ClosestWirePos_Lo( uLo ); 
	  vLo = trk->GetPair_Lo()->ClosestWirePos_Lo( vLo ); 
	  
	  trk->Set_params( vLo+dWire/2., 
			   uLo+dWire/2., 
			   vHi+dWire/2.,
			   uHi+dWire/2., 0. ); 
	  
	  double new_theta = trk->Phi(); 
	  
	  v.push_back( old_theta - new_theta  ); 
	}
	
	return v; 
      }, {"tracks_L"})
    
    .Histo1D({"h_tmodel", "#Phi_{grid-search} -  #Phi_{actual} (rad)", 
    200, -0.03, 0.03}, "theta_test"); */ 
  
  
  /*auto h2d_Test = nEvents_1track 
    
    .Define("int_u1", [](ROOT::RVec<TvdcTrack*> tracks) 
	    { return tracks.at(0)->Intercept(0); }, {"tracks_L"})
    .Define("int_v1", [](ROOT::RVec<TvdcTrack*> tracks) 
	    { return tracks.at(0)->Intercept(1); }, {"tracks_L"})
    
    .Histo2D({"test", "Real tracks;U1-intercept (m);V1-intercept(m);", 
	  200, -1.6, 1.6, 200, -1.6, 1.6}, "int_u1", "int_v1"); 

	  h2d_Test->DrawCopy("col"); 
	  
	  
	  
	  auto h2d_S2_xy = nEvents_1track 
	  .Histo2D({"S2_xy", "S2-intercept - X&Y;X_{S2} (m);Y_{S2}", 
    200, -1.2, 1.0, 50, -0.15, 0.30}, "tracks_s2_x", "tracks_s2_y");  */ 
  
    
  TString canv_title = TString::Format("run %4i, %5e evts,",
				       runNum,(double)maxEntries); 
  
#if IS_MONTECARLO 
  canv_title += choose_goodOrBad ? " GOOD tracks" : " BAD tracks"; 
#endif 
  
  canv_title += isRightArm ? "Right arm" : "Left arm"; 
  
  auto stopwatch = new TStopwatch; 
  
  int nEventsPass_trueTrack 
    = *( nEvents_trueTrack 
	 .Filter([](int NGT) {return NGT>0;}, {"nGoodTracks"}) ).Count();  
  
  int nEventsPass_track = *nEvents_1track.Count(); 
  int nEventsPass_pair  = *nEvents_1pair .Count();
  int nEventsPass_group = *nEvents_1group.Count(); 
  int nEventsPass_coinc = *nEvents_coinc .Count();
    
  double net_time = stopwatch->RealTime(); 
    
  
  auto Phi_model_draw   = [](double *S2_y, double *offset) { 
    
    return 
    -0.231*TMath::Power(S2_y[0],2) 
    +  0.2532*S2_y[0]
    +  0.539e-3  + offset[0]; 
  };
  auto Theta_model_draw = [](double *S2_x, double *offset) { 
    
    return 0.109648*S2_x[0] + offset[0]; 
  };
    
  auto canv_1 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);   //*/ 
  

  auto allTrack_err = nEvents_trueTrack
    .Define("trk_err", [](ROOT::RVec<TvdcTrack*> tracks, 
			  RVecD mTrk_good, 
			  RVecD mTrk_accident) { 
	      RVecD err; 
	      
	      for (int t=0; t<tracks.size(); t++) { 
		//for (int p=0; p<4; p++) { 
		err.push_back( tracks.at(t)->T0()-mTrk_good[4] );
		err.push_back( tracks.at(t)->T0()-mTrk_accident[4] );
		  //}
	      }
	      return err; }, {"tracks_L", 
		"MONTE_track_params",
		"MONTE_track_params_accident"})
    
    .Histo1D({"trk_err", "Track T-offset error (s)", 200, -100e-9, 100e-9}, "trk_err"); 

  allTrack_err->DrawCopy(); 
  
  draw_gausFit( (TH1D*)allTrack_err->Clone(), 40e-9 ); 
	  
  
  /*gStyle->SetOptTitle(0); 
    
  auto nTracks_good = nEvents_trueTrack
    .Histo1D({"nGood", "good tracks", 11, -0.5, 10.5}, "nGoodTracks"); 
  
  nTracks_good->Scale( 1./nTracks_good->Integral() ); 
  nTracks_good->GetYaxis()->SetRangeUser( 0., 1. ); 
  nTracks_good->DrawCopy("HIST"); 
  
  auto nTracks_all = nEvents_trueTrack
    .Define("nAllTracks", [](ROOT::RVec<TvdcTrack*> tracks) 
	    { return tracks.size(); }, {"tracks_L"})
    .Histo1D({"nAll", "all tracks", 11, -0.5, 10.5}, "nAllTracks"); 
  
  nTracks_all->Scale( 1./nTracks_all->Integral() ); 
  nTracks_all->SetLineColor(kRed); 
  nTracks_all->DrawCopy("HIST SAME"); 
  
  gPad->BuildLegend(); //*/ 
  
  
  //h2d_xParam_Dt->DrawCopy("col"); 
  
  
  /*auto canv_3 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);   
  
  h2d_Theta_s2x->DrawCopy("col"); 
    
  auto fThetaCut = new TF1("thetaCut", Theta_model_draw, 
			   -1.2, 1.2, 1);
  
  fThetaCut->SetLineStyle(kDashed); 
  
  fThetaCut->SetParameter(0,CUT_th_min); 
  fThetaCut->DrawCopy("SAME"); 
  
  fThetaCut->SetParameter(0,CUT_th_max); 
  fThetaCut->DrawCopy("SAME"); //*/ 

  
  /*auto canv_4 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);
  
  h2d_Phi_s2y->DrawCopy("col"); 
  
  auto fPhiCut = new TF1("phiCut", Phi_model_draw,
			 -0.32, 0.32, 1);
  
  fPhiCut->SetLineStyle(kDashed); 
  
  fPhiCut->SetParameter(0,CUT_ph_min); 
  fPhiCut->DrawCopy("SAME"); 
  
  fPhiCut->SetParameter(0,CUT_ph_max); 
  fPhiCut->DrawCopy("SAME");            //*/ 
  
  
  /*h2d_trackErr->DrawCopy("surf"); 
  
  auto xAxis = h2d_trackErr->GetXaxis(); 
  
  
  double err_X[xAxis->GetNbins()]; xAxis->GetCenter( err_X );  
  double err_Y[xAxis->GetNbins()]; 

  for (int b=1; b<=xAxis->GetNbins(); b++) { 
    
    auto proj = h2d_trackErr->ProjectionY("proj", b,b); 
    
    auto fit = draw_gausFit_fixedBase( proj, 0.6e-3, 0, 0, -1 );   
    
    err_Y[b-1] = TMath::Abs(fit.Y()); 
  }
  
  
  auto err_graph = new TGraph(xAxis->GetNbins(), err_X, err_Y); 
  
  auto lineFit = new TF1("fitLine", "x", 0, 3e-4); 
  lineFit->SetLineStyle(kDashed); 
  lineFit->Draw(); 
  
  err_graph->GetXaxis()->SetRangeUser( 0, 3e-4 );
  err_graph->GetYaxis()->SetRangeUser( 0, 3e-4 );
  
  err_graph->SetTitle("Error predictor;actual (m);real (m)"); 
 
  
  err_graph->SetMarkerStyle(kOpenCircle); 
  err_graph->Draw("SAME P"); //*/ 
    
    
  
  //h2d_Phi_fpx->DrawCopy("col"); 
  
  //
  //h2d_Theta_fpx->DrawCopy("col");
  
  /*auto g = th2d_to_graph( (TH2D*)h2d_Theta_fpx->Clone(), 
			  -0.50, 0.50, 0.03 ); 
  g->SetMarkerStyle(kOpenCircle); 
  g->SetMarkerSize(0.5); 
  g->Draw("SAME P");
  
  auto fit = new TF1("fit", "[0]*x + [1]", -0.45, 0.40); 
  
  g->Fit("fit"); //*/ 
  
  
  /*auto canv_2 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);   
   
  h1d_eta->DrawCopy(); */ 

  /*auto h1d_nGoodTracks = nEvents_trueTrack 
    .Histo1D({"NGT", "N. 'good' tracks (Err < 0.75 mm)", 5, -0.5, 4.5}, 
	     "nGoodTracks"); 
  
  h1d_nGoodTracks->Scale( 1./h1d_nGoodTracks->Integral() ); 
  
  h1d_nGoodTracks->GetYaxis()->SetRangeUser(0, 1.); 
  h1d_nGoodTracks->DrawCopy("E"); 
    
  //h2d_xParam_Dt->DrawCopy("col"); //*/ 
  
  
  /*auto canv_3 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);   
  
  h1d_dv->GetYaxis()->SetRangeUser(0., h1d_dv->GetMaximum()*1.2); 
  h1d_dv->DrawCopy("E"); //*/ 
  
    
  /*auto canv_4 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);   
  
  h1d_eta->DrawCopy(); //*/ 
  
  
  cout << TString::Format(" time taken = %3.1f (%8.8f ms/entry)", 
			  net_time, 
			  1e3*net_time/((double)maxEntries) ) << endl; 
  
  
  cout << "coinc events           " << nEventsPass_coinc << endl; 
  
  
  cout << "events with 1 group    " << nEventsPass_group 
       << TString::Format(" (%0.4f)", 
			  ((float)nEventsPass_group)/((float)nEventsPass_coinc) ) 
       << endl; 
 
  cout << "events with 1 pair     " << nEventsPass_pair 
       << TString::Format(" (%0.4f)", 
			  ((float)nEventsPass_pair)/((float)nEventsPass_coinc) ) 
       << endl; 
  
  cout << "events with 1 track    " << nEventsPass_track
       << TString::Format(" (%0.4f)", 
			  ((float)nEventsPass_track)/((float)nEventsPass_coinc) ) 
       << endl; 
  
  cout << "events with 'true' track found  " << nEventsPass_trueTrack
       << TString::Format(" (%0.4f)", 
			  ((float)nEventsPass_trueTrack)/((float)nEventsPass_coinc) ) 
       << endl; 
  
  
 
  return; 
}

#if 0 

//using newton's method
auto Refine_track = [CUT_th_max, CUT_th_min,
		       CUT_ph_max, CUT_ph_min]( TvdcTrack *trk, 
						const int nCycles=10, 
						double sigma=35e-9 ) { 
    
    //random-walk track minimizaiton
    const double GRAD_momentum = 0.50; 
    const double GRAD_step0    = 0.05; 
    const double GRAD_exponent = 0.33; 
    
    //cout << "refining track..." << endl; 
    
    const double scale_T = 3e-9; 
    const double scale_X = 2.2e-3; 
    
    const double nudge_multiplier = 1.; 
    
    const double Chi_cutoff = 6.; 
    
    const double sigma_decay = 1.;
    const double min_sigma   = sigma*2.; 
    
    ///////////////////////////////////////////////////////////////////////////////
    double int0[5]; 
    for (int p=0; p<4; p++) int0[p] = trk->Intercept(p); 
    TString plane_name[4] = {"U1","V1","U2","V2"}; 
    int0[4] = trk->T0(); 
    //*/////////////////////////////////////////////////////////////////////////////
    
    const double s2 = sigma*sigma; 
    
    ROOT::RVec<double> tVec; 
    
    double objective_eta0(0.); 
    
    for (UInt_t c=0; c<nCycles; c++) { 
      
      //this measures if 'eta' is actually increasing
      double objective_eta1(0.); 
      
      //this will be how we 'nudge' the parameters
      double nudge[5] = {0.}; 
      
      double Deriv[5] = {0.}; 
      
      /*///////////////////////////////////////////////////// DRAW MINIMIZE
      auto canv = new TCanvas("ePic", "Event draw", 0,0, 900, 500); 
      canv->Divide( 2,2 ); 
      
      gStyle->SetOptStat(0); 
      int iCanvas=1; 
      //*////////////////////////////////////////////////////////////////////////
      
      //first-derivatives
      double F[5] = {0.}; 
      
      double J_ii[4] = {0.}; 
      double J_4i[4] = {0.}; 
      double J_44(0.); 
      
      for (int p=0; p<4; p++) { 
	
	/*/draw event picture ///////////////////////////////////////////////////
	canv->cd( iCanvas );  iCanvas++; 
	auto h2d = new TH2D(TString::Format("vdcPlane_%i",p), "", 
			    200, -25e-3, 25e-3, 
			    200, VDC_min_realTime, VDC_max_realTime); 
	//*//////////////////////////////////////////////////////////////////////
	
	auto group = trk->GetGroup(p); 
	
	for (UInt_t h=0; h<group->Nhits(); h++) { 
	  	  
	  /*/draw event picture ///////////////////////////////////////////////
	  h2d->Fill( group->WirePos(h)-int0[p], group->Time(h) ); 	  
	  //*/////////////////////////////////////////////////////////////////
	  
	  double Chi 
	    = trk->Get_T_model(p, group->WirePos(h))  
	    +  trk->T0() 
	    -  group->Time(h); 
	  
	  if ( TMath::Abs(Chi) > Chi_cutoff*sigma ) continue; 
	  
	  double Eta = TMath::Exp( -0.5*TMath::Power(Chi/sigma,2) ); 
	  
	  double m = trk->Slope(p); 
	  
	  F[p] += -m * Eta * Chi * trk->Get_T_model(p,group->WirePos(h),1); 
	  
	  F[4] +=  Eta * Chi; 
	
	  J_ii[p] 
	    += (Chi * trk->Get_T_model(p,group->WirePos(h),2)  
		+ TMath::Power(trk->Get_T_model(p,group->WirePos(h),1),2))*m*m*Eta; 
	  
	  J_4i[p] 
	    +=  -Eta * m * trk->Get_T_model(p,group->WirePos(h),1); 
	  
	  J_44 += Eta; 
	  
	  objective_eta1 += TMath::Exp( -0.5*TMath::Power(Chi/min_sigma,2) ); 
	  
	}//for (UInt_t h=0; h<trk->pGet_N(p); 
	
	/*/////////////////////////////////////////////////////////////////////////
	h2d->SetMarkerStyle(kOpenCircle); 
	h2d->SetMarkerSize(0.50); 
	h2d->DrawCopy(); 
	
	auto ltx = new TLatex; 
	
	ltx->DrawLatex( -23e-3, VDC_max_realTime, 
			plane_name[p] ); 
    	
	const int nPts = 100; 
	
	double max_X = 14e-3; 
	double X[nPts]; 
	double T[nPts]; 
	
	double dx = 2.*max_X/((double)nPts-1); 
	
	double xx=(trk->Intercept(p)-int0[p])-max_X; 
	for (int i=0; i<nPts; i++) { 
	  X[i] = xx;  xx += dx; 
	  T[i] = trk->Get_T_model( p, xx+int0[p] ) + trk->T0(); 
	}	
	auto g = new TGraph(nPts, X, T); 
	g->SetLineColor(kRed); 
	g->Draw("SAME"); 
	
	h2d->~TH2D(); 
	ltx->~TLatex(); 
	//*////////////////////////////////////////////////////////////////////////
		
      }//for (int p=0; p<4; p++) 
      
      /*///////////////////////////////////////////////////////////////////////////
      canv->Print("~/ftp-dump/plots/test_track_hi_Dt.gif+33");
      canv->~TCanvas(); 
      //*//////////////////////////////////////////////////////////////////////////
      
      /*cout << TString::Format("eta-obj = %7.7f (delta = %7.7e)", 
	objective_eta1,objective_eta1-objective_eta0); */ 
      
      //if (objective_eta1 < objective_eta0) break; 
      objective_eta0 = objective_eta1; 
          
      sigma *= sigma_decay; 
      
      
      for (int p=0; p<4; p++) { 
	
	J_44  +=  J_4i[p]*(-J_4i[p]/J_ii[p]); 
	F[4]  +=  F[p]   *(-J_4i[p]/J_ii[p]); 
      }
      
      nudge[4] = -nudge_multiplier*F[4]/J_44; 
      
      for (int p=0; p<4; p++) { 
	
	nudge[p] = -nudge_multiplier*(F[p] - J_4i[p]*nudge[4])/J_ii[p]; 
      }
      
      trk->Nudge_params( nudge ); 
      
      /*cout << TString::Format("   change (T)= %5.7f ns", (trk->T0()-int0[4])*1e9)
	<< endl; */ 
    }
    
    trk->SetEta( objective_eta0 ); 
  };   

 auto Refine_track = []( TvdcTrack *trk, 
			  const int nCycles=8, 
			  double sigma=100e-9 ) { 
    
    //random-walk track minimizaiton
    const double GRAD_momentum = 0.50; 
    const double GRAD_step0    = 0.05; 
    const double GRAD_exponent = 0.33; 
    
    cout << "refining track..." << endl; 
    
    const double scale_T = 3e-9; 
    const double scale_X = 2.2e-3; 
    
    const double Chi_cutoff = sigma*6.; 
    
    const double sigma_decay = 0.85;
    const double min_sigma   = sigma*2.; 
    
    ///////////////////////////////////////////////////////////////////////////////
    double int0[5]; 
    for (int p=0; p<4; p++) int0[p] = trk->Intercept(p); 
    TString plane_name[4] = {"U1","V1","U2","V2"}; 
    int0[4] = trk->T0(); 
    //*/////////////////////////////////////////////////////////////////////////////
    
    double eta_prev(0.); 
    
    const double s2 = sigma*sigma; 
    
    double B[5] = {0.}; 

    ROOT::RVec<double> tVec; 
    
    double objective_eta0(0.); 
    
    for (UInt_t c=0; c<nCycles; c++) { 
      
      //this measures if 'eta' is actually increasing
      double objective_eta1(0.); 
      
      //this will be how we 'nudge' the parameters
      double nudge[5]; 
      
      double Deriv[5] = {0.}; 
      
      ////////////////////////////////////////////////////// DRAW MINIMIZE
      auto canv = new TCanvas("ePic", "Event draw", 0,0, 900, 500); 
      
      canv->Divide( 2,2 ); 
      
      gStyle->SetOptStat(0); 
      int iCanvas=1; 
      //*////////////////////////////////////////////////////////////////////////
      
      //first-derivatives
      double F[5] = {0.}; 
      
      double J_ii[4] = {0.}; 
      double J_4i[4] = {0.}; 
      double J_44(0.); 
      
      for (int p=0; p<4; p++) { 
	
	//draw event picture ///////////////////////////////////////////////////
	canv->cd( iCanvas );  iCanvas++; 
	
	auto h2d = new TH2D(TString::Format("vdcPlane_%i",p), "", 
			    200, -25e-3, 25e-3, 
			    200, VDC_min_realTime, VDC_max_realTime); 
	//*//////////////////////////////////////////////////////////////////////
	
	auto group = trk->GetGroup(p); 
	
	for (UInt_t h=0; h<group->Nhits(); h++) { 
	  	  
	  //draw event picture ///////////////////////////////////////////////
	  h2d->Fill( group->WirePos(h)-int0[p], group->Time(h) ); 	  
	  //*/////////////////////////////////////////////////////////////////
	  
	  //our current model for Tau
	  double T_mod[3] 
	    = { trk->Get_T_model(p,group->WirePos(h),0),   //0-th derivative
		trk->Get_T_model(p,group->WirePos(h),1),   //1-st derivative 
		trk->Get_T_model(p,group->WirePos(h),2) }; //2-nd derivative 
	  
	  double Chi 
	    = T_mod[0]
	    +  trk->T0() 
	    -  group->Time(h); 
	  
	  if ( TMath::Abs(Chi) > Chi_cutoff ) continue; 
		  
	  double Chi2 = TMath::Power( Chi/sigma, 2 ); 
		  
	  double Eta = TMath::Exp( -0.5*Chi2 ); 
	  
	  double m = trk->Slope(p); 
	  
	  F[p] +=  m * Eta * Chi * T_mod[1]; 
	  
	  F[4] +=  -Eta * Chi; 
	  
	  J_ii[p] 
	    += Eta*( Chi2*T_mod[1]  
		     -  T_mod[1] 
		     -  T_mod[2]/T_mod[1]*Chi )*T_mod[1]*m*m; 
	  
	  J_4i[p] 
	    +=  Eta*( 1.  -  Chi2 )*T_mod[1]*m; 
	  
	  J_44 += Eta*( Chi2  +  1. ); 
	  
	  objective_eta1 += TMath::Exp( -0.5*TMath::Power(Chi/min_sigma,2) ); 
	  
	}//for (UInt_t h=0; h<trk->pGet_N(p); 
	
	//////////////////////////////////////////////////////////////////////////
	h2d->SetMarkerStyle(kOpenCircle); 
	h2d->SetMarkerSize(0.50); 
	h2d->DrawCopy(); 
	
	auto ltx = new TLatex; 
	
	ltx->DrawLatex( -23e-3, VDC_max_realTime, 
			plane_name[p] ); 
    	
	const int nPts = 100; 
	
	double max_X = 14e-3; 
	double X[nPts]; 
	double T[nPts]; 
	
	double dx = 2.*max_X/((double)nPts-1); 
	
	double xx=(trk->Intercept(p)-int0[p])-max_X; 
	for (int i=0; i<nPts; i++) { 
	  X[i] = xx;  xx += dx; 
	  T[i] = trk->Get_T_model( p, xx+int0[p] ) + trk->T0(); 
	}	
	auto g = new TGraph(nPts, X, T); 
	g->SetLineColor(kRed); 
	g->Draw("SAME"); 
	
	h2d->~TH2D(); 
	ltx->~TLatex(); 
	//*////////////////////////////////////////////////////////////////////////
		
      }//for (int p=0; p<4; p++) 
      
      ////////////////////////////////////////////////////////////////////////////
      canv->Print("plots/test_track.gif+100");
      canv->~TCanvas(); 
      //*//////////////////////////////////////////////////////////////////////////
      
      cout << TString::Format("eta-obj = %7.7f (delta = %7.7e)", 
			      objective_eta1,objective_eta1-objective_eta0); 
      
      //if (objective_eta1 < objective_eta0) break; 
      objective_eta0 = objective_eta1; 
          
      //sigma *= sigma_decay; 
      
      
      //nudge the track
      for (int p=0; p<4; p++) { 
	F[p]    *= -1/s2; 
	J_ii[p] *= 1/s2; 
	J_4i[p] *= 1/s2; 
      }
      F[4] *= -1/s2; 
      J_44 *= 1/s2; 
      

      for (int p=0; p<4; p++) { 
	
	J_44 += J_4i[p]*(-J_4i[p]/J_ii[p]); 
	F[4] += F[p]   *(-J_4i[p]/J_ii[p]); 
      }
            
      nudge[4] = F[4]/J_44; 
      
      for (int p=0; p<4; p++) { nudge[p] = (F[p] - J_4i[p]*nudge[4])/J_ii[p]; }
      
      trk->Nudge_params( nudge ); 
      
      cout << TString::Format("   change (T)= %5.7f ns", (trk->T0()-int0[4])*1e9)
	   << endl; 
    }
    
    //noop
    return 1.; 
 };   



auto Refine_track = []( TvdcTrack *trk, 
			int nCycles=60, 
			const double sigma=50e-9 ) { 
  
    //random-walk track minimizaiton
    const double GRAD_momentum = 0.50; 
    const double GRAD_step0    = 0.05; 
    const double GRAD_exponent = 0.33; 
    
    cout << "refining track..." << endl; 
    
    const double scale_T = 3e-9; 
    const double scale_X = 2.2e-3; 
    
    const double Chi_cutoff = sigma*6.; 
    
    ///////////////////////////////////////////////////////////////////////////////
    double int0[5]; 
    for (int p=0; p<4; p++) int0[p] = trk->Intercept(p); 
    TString plane_name[4] = {"U1","V1","U2","V2"}; 
    int0[4] = trk->T0(); 
    ///////////////////////////////////////////////////////////////////////////////
    
    const double s2 = sigma*sigma; 
    
    double B[5] = {0.}; 

    ROOT::RVec<double> tVec; 
    
    for (UInt_t c=0; c<nCycles; c++) { 
      
      //get random direction vector
      double rVec[5]; // = { gRandom->Gaus() }; 
      
      //normalize this random vector
      double modulo(0.); 
      for (int j=0; j<5; j++) { 
	rVec[j] = gRandom->Gaus(); 
	modulo += rVec[j]*rVec[j]; 
      }
      modulo = TMath::Sqrt( modulo ); 
      
      for (int j=0; j<5; j++) { rVec[j] *= 1./modulo; }
            
      
      //this will be how we 'nudge' the parameters
      double nudge[5]; 
      
      double Deriv[5] = {0.}; 
      
      ////////////////////////////////////////////////////// DRAW MINIMIZE
      auto canv = new TCanvas("ePic", "Event draw", 0,0, 900, 500); 
      
      canv->Divide( 2,2 ); 
      
      gStyle->SetOptStat(0); 
      int iCanvas=1; 
      //////////////////////////////////////////////////////////////////////////
      
      for (int p=0; p<4; p++) { 
	
	//draw event picture ///////////////////////////////////////////////////
	canv->cd( iCanvas );  iCanvas++; 
	
	auto h2d = new TH2D(TString::Format("vdcPlane_%i",p), "", 
			    200, -25e-3, 25e-3, 
			    200, VDC_min_realTime, VDC_max_realTime); 
	////////////////////////////////////////////////////////////////////////
	
	auto group = trk->GetGroup(p); 
	
	for (UInt_t h=0; h<group->Nhits(); h++) { 
	  	  
	  //draw event picture ///////////////////////////////////////////////
	  h2d->Fill( group->WirePos(h)-int0[p], group->Time(h) ); 	  
	  /////////////////////////////////////////////////////////////////////	  
	  
	  double Chi 
	    = trk->Get_T_model(p, group->WirePos(h))  
	    +  trk->T0() 
	    -  group->Time(h); 
	  
	  if ( TMath::Abs(Chi) > Chi_cutoff ) continue; 
	  
	  double Eta = TMath::Exp( -0.5*TMath::Power(Chi/sigma, 2) ); 
	  
	  double m = trk->Slope(p); 
	  
	  Deriv[p] += Eta * Chi * m * trk->Get_T_model(p,group->WirePos(h),1); 
	  	  
	  Deriv[4] += Eta * Chi; 
	}//for (UInt_t h=0; h<trk->pGet_N(p); 
		
	//////////////////////////////////////////////////////////////////////////
	h2d->SetMarkerStyle(kOpenCircle); 
	h2d->SetMarkerSize(0.50); 
	h2d->DrawCopy(); 
	
	auto ltx = new TLatex; 
	
	ltx->DrawLatex( -23e-3, VDC_max_realTime, 
			plane_name[p] ); 
    
	
	const int nPts = 100; 
	
	double max_X = 14e-3; 
	double X[nPts]; 
	double T[nPts]; 
	
	double dx = 2.*max_X/((double)nPts-1); 
	
	double xx=(trk->Intercept(p)-int0[p])-max_X; 
	for (int i=0; i<nPts; i++) { 
	  X[i] = xx;  xx += dx; 
	  T[i] = trk->Get_T_model( p, xx+int0[p] ); 
	}	
	auto g = new TGraph(nPts, X, T); 
	g->SetLineColor(kRed); 
	g->Draw("SAME"); 
	
	h2d->~TH2D(); 
	ltx->~TLatex(); 
	//////////////////////////////////////////////////////////////////////////
		
      }//for (int p=0; p<4; p++) 
            
      ////////////////////////////////////////////////////////////////////////////
      canv->Print("~/ftp-dump/plots/test_track.gif+10");
      canv->~TCanvas(); 
      ////////////////////////////////////////////////////////////////////////////
      
      for (int p=0; p<4; p++) Deriv[p] *= scale_X/s2; 
      
      Deriv[4] *= -scale_T/s2; 
      
      
      double K = GRAD_step0*TMath::Power(c+1,-GRAD_exponent); 
      
      cout << TString::Format(" cycle =%3i eta =%i "
      //now, nudge the parameters
      for (int i=0; i<5; i++) { 
	
	double scale = i<4 ? scale_X : scale_T; 
	
	B[i] = GRAD_momentum*B[i]  +  K * Deriv[i] * scale; 
	
	nudge[i] = B[i]; 
	
	if (i<4) {
	  cout << TString::Format("%7.5f, ", 1e3*(trk->Intercept(i)-int0[i]) );
	} else   { 
	  cout << TString::Format("%7.5f, ", 1e9*(trk->T0()-int0[i]) );
	}
      }      
      cout << endl; 
      
      trk->Nudge_params( nudge ); 
    }
    
    //noop
    return 1.; 
  };   
    
    auto Refine_track = []( TvdcTrack *trk, 
			    const int nCycles=8, 
			    double sigma=50e-9 ) { 
      
      //random-walk track minimizaiton
    const double GRAD_momentum = 0.50; 
    const double GRAD_step0    = 0.05; 
    const double GRAD_exponent = 0.33; 
    
    cout << "refining track..." << endl; 
    
    const double scale_T = 3e-9; 
    const double scale_X = 2.2e-3; 
    
    const double Chi_cutoff = sigma*6.; 
    
    const double sigma_decay = 0.85;
    const double min_sigma   = sigma*2.; 
    
    ///////////////////////////////////////////////////////////////////////////////
    double int0[5]; 
    for (int p=0; p<4; p++) int0[p] = trk->Intercept(p); 
    TString plane_name[4] = {"U1","V1","U2","V2"}; 
    int0[4] = trk->T0(); 
    //*/////////////////////////////////////////////////////////////////////////////
    
    double eta_prev(0.); 
    
    const double s2 = sigma*sigma; 
    
    double B[5] = {0.}; 

    RVecD tVec; 
    
    double objective_eta0(0.); 
    
    for (UInt_t c=0; c<nCycles; c++) { 
      
      //this measures if 'eta' is actually increasing
      double objective_eta1(0.); 
      
      //this will be how we 'nudge' the parameters
      double nudge[5]; 
      
      double Deriv[5] = {0.}; 
      
      ////////////////////////////////////////////////////// DRAW MINIMIZE
      auto canv = new TCanvas("ePic", "Event draw", 0,0, 900, 500); 
      
      canv->Divide( 2,2 ); 
      
      gStyle->SetOptStat(0); 
      int iCanvas=1; 
      //*////////////////////////////////////////////////////////////////////////
      
      //first-derivatives
      double F[5] = {0.}; 
      
      double J_ii[4] = {0.}; 
      double J_4i[4] = {0.}; 
      double J_44(0.); 
      
      for (int p=0; p<4; p++) { 
	
	//draw event picture ///////////////////////////////////////////////////
	canv->cd( iCanvas );  iCanvas++; 
	
	auto h2d = new TH2D(TString::Format("vdcPlane_%i",p), "", 
			    200, -25e-3, 25e-3, 
			    200, VDC_min_realTime, VDC_max_realTime); 
	//*//////////////////////////////////////////////////////////////////////
	
	auto group = trk->GetGroup(p); 
	
	for (UInt_t h=0; h<group->Nhits(); h++) { 
	  	  
	  //draw event picture ///////////////////////////////////////////////
	  h2d->Fill( group->WirePos(h)-int0[p], group->Time(h) ); 	  
	  //*/////////////////////////////////////////////////////////////////
	  
	  double Chi 
	    = trk->Get_T_model(p, group->WirePos(h))  
	    +  trk->T0() 
	    -  group->Time(h); 
	  
	  if ( TMath::Abs(Chi) > Chi_cutoff ) continue; 
	  
	  double Eta = TMath::Exp( -0.5*TMath::Power(Chi/sigma,2) ); 
	  
	  double m = trk->Slope(p); 
	  
	  F[p] += -m * Eta * Chi * trk->Get_T_model(p,group->WirePos(h),1); 
	  
	  F[4] +=  Eta * Chi; 
	
	  J_ii[p] 
	    += (Chi * trk->Get_T_model(p,group->WirePos(h),2)  
		+ TMath::Power(trk->Get_T_model(p,group->WirePos(h),1),2))*m*m*1.; 
	  
	  J_4i[p] 
	    +=  1. * m * trk->Get_T_model(p,group->WirePos(h),1); 
	  
	  J_44 += 1.; 
	  
	  objective_eta1 += TMath::Exp( -0.5*TMath::Power(Chi/min_sigma,2) ); 
	  
	}//for (UInt_t h=0; h<trk->pGet_N(p); 
	
	//////////////////////////////////////////////////////////////////////////
	h2d->SetMarkerStyle(kOpenCircle); 
	h2d->SetMarkerSize(0.50); 
	h2d->DrawCopy(); 
	
	auto ltx = new TLatex; 
	
	ltx->DrawLatex( -23e-3, VDC_max_realTime, 
			plane_name[p] ); 
    	
	const int nPts = 100; 
	
	double max_X = 14e-3; 
	double X[nPts]; 
	double T[nPts]; 
	
	double dx = 2.*max_X/((double)nPts-1); 
	
	double xx=(trk->Intercept(p)-int0[p])-max_X; 
	for (int i=0; i<nPts; i++) { 
	  X[i] = xx;  xx += dx; 
	  T[i] = trk->Get_T_model( p, xx+int0[p] ) + trk->T0(); 
	}	
	auto g = new TGraph(nPts, X, T); 
	g->SetLineColor(kRed); 
	g->Draw("SAME"); 
	
	h2d->~TH2D(); 
	ltx->~TLatex(); 
	//*////////////////////////////////////////////////////////////////////////
		
      }//for (int p=0; p<4; p++) 
      
      ////////////////////////////////////////////////////////////////////////////
      canv->Print("plots/test_track.gif+100");
      canv->~TCanvas(); 
      //*//////////////////////////////////////////////////////////////////////////
      
      cout << TString::Format("eta-obj = %7.7f (delta = %7.7e)", 
			      objective_eta1,objective_eta1-objective_eta0); 
      
      //if (objective_eta1 < objective_eta0) break; 
      objective_eta0 = objective_eta1; 
          
      //sigma *= sigma_decay; 
      
      
      for (int p=0; p<4; p++) { J_44 += -J_4i[p]*J_4i[p]/J_ii[p]; }
      
      nudge[4] = F[4]/J_44; 
      
      for (int p=0; p<4; p++) { nudge[p] = -(F[p] - J_4i[p]*nudge[4])/J_ii[p]; }
      
      trk->Nudge_params( nudge ); 
      
      cout << TString::Format("   change (T)= %5.7f ns", (trk->T0()-int0[4])*1e9)
	   << endl; 
    }
    
    //noop
    return 1.; 
  };   

    
    
    
 
#endif 
