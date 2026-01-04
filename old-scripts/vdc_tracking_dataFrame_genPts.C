#include "def_apex.h"
#include "function_lib.C"
#include <fstream>

//#include "class/TvdcHit.h"
//#include "class/TvdcHit.C"
#include "lib/NewTrack.h"
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
  
  auto dummy_event = new TEventHandler(false, 4329, 53., dummy_s2R, dummy_s2L); 
  
  return dummy_event->Drift_T(x[0], 1.4, 2 );
}  
//#endif 

#define OUTPUT_U_PLANE true  //
#define IS_MONTECARLO false   //
#define RUN_LIMITED false     //
#define MAKE_OUTPUT true     //
#define DEBUG false           //


//track cuts
double TRK_CUT_Dt     = 27.5e-9; 
double TRK_CUT_xParam = 1.0; 
double TRK_CUT_Eta    = 4.250; 
double TRK_measureSigma = 9e-9; 
int    TRK_CUT_nGoodPts_min = 2; //good points per plane

double CUT_minEta     = 2.95;  //min eta of one plane (during grid-searching). 
                               // this is analogous to the num. of points found
  


RNode generate_vdcTracks( bool is_RightArm, RNode inNode, 			  
			  ROOT::RDF::RResultPtr<unsigned long long> &nPass_1group, 
			  ROOT::RDF::RResultPtr<unsigned long long> &nPass_1pair, 
			  ROOT::RDF::RResultPtr<unsigned long long> &nPass_1rawTrack
			  ) { 
  
  //feed this script your node, with generated S2-hits, ready for tracking data. 
  
  //it'll spit out refined tracks for the plane you tell it to. 
  
  string arm = is_RightArm ? "R" : "L" ; 
  
  vector<string> plane_name = { "u1", "v1", "u2", "v2" }; 
  
  vector<string> branch_rawtime; 
  vector<string> branch_wire; 
  
  for (int p=0; p<4; p++) { 
    
    string rawtime = IS_MONTECARLO 
      ? arm+"_vdc_"+plane_name[p]+"_rawtime" 
      : arm+".vdc."+plane_name[p]+".rawtime"; 
    
    string wire    = IS_MONTECARLO
      ? arm+"_vdc_"+plane_name[p]+"_wire" 
      : arm+".vdc."+plane_name[p]+".wire"; 
    
    branch_rawtime.push_back( rawtime ); 
    branch_wire   .push_back( wire );   
  }
    
  //
  //     Create 'hit groups' for a particular plane/arm.      
  //
  auto  group_hits = [] ( TEventHandler *evt, //true = RHRS, false = LHRS
			  int   p, //the plane we're working with
			  RVecD h_wire, 
			  RVecD h_time ) {
#if DEBUG 
    int n_validHits(0); 
    cout << (evt->ActiveArm() ? "arm=Right" : "arm=Left") << endl;
#endif 

    const unsigned int clust_minHits=2; //minimum number of hits in a group
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
      
#if DEBUG 
      cout << TString::Format("wire=%3i, rawtime=%5.1f, time=%5.1f   ",
			      (int)h_wire[h], 
			      h_time[h], 
			      hit->Time()*1e9 ); 
#endif
      
      
      if (hit->Time() > VDC_max_realTime || 
	  hit->Time() < VDC_min_realTime ) {

#if DEBUG
	cout << "killed! " << endl; 
#endif 
	
	hit->~TvdcHit(); //throw out this hit 
	continue; 
      }
      
#if DEBUG
      cout << "kept!" << endl; 
      n_validHits++; 
#endif 
      
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
    
#if DEBUG 
    cout << TString::Format( "group_hits() -  total hits:%3i, valid-hits:%3i, groups made:%2i", 
			     (int)h_wire.size(), 
			     n_validHits, 
			     (int)groupVec.size() ) << endl; 
#endif 
    
    return groupVec; 
        
  }; 

  
  
  /////////////////////////////////////////////////////////////////////////////
  auto Grid_search = [] ( TEventHandler *evt, 
			  THitGroup *group, 
			  double m1, //slope 
			  double m2, 
			  double x_Lo,
			  double x_Hi, 
			  double TAU_sigma  =9e-9, 
			  double TAU_buffer =20e-9) {
	
    double eta(0.); 
	
    double m_min = TMath::Min( m1, m2 ); 
    double m_max = TMath::Max( m1, m2 ); 
    
    double m_avg = 0.5*(m1 + m2); 
    
    double x_span = x_Hi - x_Lo; 
	
    double x_avg = (x_Lo + x_Hi)/2.; 

    //loop over all hits
    for (int h=0; h<group->Nhits(); h++) { 
      
      double x = group->WirePos(h); 
      
      double tau_Lo; 
      
      if ( TMath::Abs(x-x_avg) < x_span/2. ) { //this wire-pos is between x_Hi & _Lo
	
	tau_Lo 
	  = evt->Drift_T( 0., m_avg ) 
	  - TAU_buffer; 
	
      } else {
	
	tau_Lo 
	  = evt->Drift_T( m_min*( x-(x<x_avg ? x_Lo : x_Hi) ), m_avg ) 
	  - TAU_buffer; 
	
	//cut hits that are too far away
	if ( tau_Lo > VDC_max_realTime ) continue; 
      }
      
      double tau_Hi 
	= evt->Drift_T( m_max*( x-(x<x_avg ? x_Hi : x_Lo) ), m_avg ) 
	+ TAU_buffer; 
      
      /*cout << TString::Format(" dv %2.1f mm, DT = %3.1f ns", 
	(x-x_avg)*1e3, (tau_Hi-tau_Lo)*1e9 ) << endl; */ 
      
      double tau = group->Time(h); 
      
      //tau is below lowest guess
      if (tau < tau_Lo) { 

	eta += TMath::Exp( -0.5*TMath::Power((tau-tau_Lo)/TAU_sigma, 2) ); 
	continue; 
      }

      //'goldilocks'-zone, between hi- & lo-predictions
      if (tau < tau_Hi) { 

	eta += 1.; 
	continue; 
      }
      
      //tau is above highest guess
      eta += TMath::Exp( -0.5*TMath::Power((tau-tau_Hi)/TAU_sigma, 2) ); 

    }//for (int h=0; h<group->Nhits(); 

    return eta; 
  };
  //////////////////////////////////////////////////////////////////////////////
  
  
  // 
  //    Take a pair of U1 & U2 groups, and generate 'chamber pairs' 
  //     (TChamberPair*) which are effectivley the stubs of new track candidates. 
  //   
  const double gridSpacing = dWire/3.; 
    
  auto Gen_pairs = [&Grid_search, gridSpacing] 
    ( TEventHandler *evt, 
      ROOT::RVec<THitGroup*> gVec_U, 
      ROOT::RVec<THitGroup*> gVec_V, 
      bool is_LoChamber ) { 
    
    //use by tracks to identify pairs later
    int pair_unique_id(0); 
    
    //constants
    bool is_RightArm = evt->ActiveArm(); 
    
    //this involves angular prediction based on typical S2-angles. 
    const double uFix = 0.026/TMath::Sqrt(2);
    
    const double wVDC = gVec_U.at(0)->W(); 

    const double M_Theta = is_RightArm ? 0.1096 : 0.1096 ; 
    const double M_Phi   = is_RightArm ? 0.242  : 0.242  ; 
    
    const double Z_fPoint_Theta = evt->GetS2Hit()->Z() - 1./M_Theta;     
    const double Z_fPoint_Phi   = evt->GetS2Hit()->Z() - 1./M_Phi; 
    
    //////////////////////////////////////////////////////////////////////////////
    auto Guess_slopes = [uFix, wVDC, 
			 Z_fPoint_Theta, 
			 Z_fPoint_Phi]( const double u, 
					const double v, 
					double &slope_u, 
					double &slope_v ) { 
      TVector3 r_xyz 
      = TvdcTrack::Rotate_uvw_to_xyz( TVector3( u-uFix, v, wVDC ) ); 
	
      TVector3 S_uvw 
      = TvdcTrack::Rotate_xyz_to_uvw( TVector3( r_xyz[0]/(r_xyz[2]-Z_fPoint_Theta), 
						r_xyz[1]/(r_xyz[2]-Z_fPoint_Phi),
						1.) );       
      slope_u = S_uvw.Z() / S_uvw.X(); 
      slope_v = S_uvw.Z() / S_uvw.Y(); 
    }; 
        
    //////////////////////////////////////////////////////////////////////////////
    auto Guess_s2x = [evt, uFix, wVDC, 
		      Z_fPoint_Theta]( const double u, 
				       const double v ) { 
      TVector3 r_xyz 
      = TvdcTrack::Rotate_uvw_to_xyz( TVector3( u-uFix, v, wVDC ) ); 
      
      return 
      r_xyz.X()*(evt->GetS2Hit()->Z() - Z_fPoint_Theta)
      /(r_xyz.Z()            - Z_fPoint_Theta); 		 
    }; 
    
    const double CUT_vu_diff_min = is_LoChamber 
    ? (is_RightArm ? -83e-3 : -85e-3) 
    : (is_RightArm ? -95e-3 : -96e-3); 
    
    const double CUT_vu_diff_max = is_LoChamber 
    ? (is_RightArm ?  58e-3 :  58e-3)
    : (is_RightArm ?  70e-3 :  70e-3); 
    
    //////////////////////////////////////////////////////////////////////////////
    auto Find_clusters = [&evt,
			  gridSpacing,
			  &Grid_search, 
			  &Guess_slopes,
			  CUT_vu_diff_min, 
			  CUT_vu_diff_max] ( bool is_Uplane, 
					     ROOT::RVec<THitGroup*> groups ) { 
      
      //how many wires over are we gonna look for tracks?
      const int   wire_buffer = 25; 
      
      //if clusters are closer than this, then join them together
      const int min_cluster_gap = 30; 
      
      ROOT::RVec<THitCluster*> clust_all; 
      
#if 0 //DEBUG      
      cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; 
#endif 
      
      for (int g=groups.size()-1; g>=0; g+=-1) { 
	
	//iterate backwards thru groups
	auto group = groups.at(g); 
	
	//find span of our groupings
	int span
	  = (int)( group->Span()/gridSpacing ) + 2*wire_buffer + 1; 
        
	double x0 = group->LoEdge() - ((double)wire_buffer)*gridSpacing; 
	
	for (int ix=0; ix<span; ix++) { 
	  
	  double x = x0 + ((double)ix)*gridSpacing; 
	  
	  double y_min = x + (is_Uplane ? CUT_vu_diff_min : -CUT_vu_diff_min); 
	  double y_max = x + (is_Uplane ? CUT_vu_diff_max : -CUT_vu_diff_max); 
	  
	  double m1, m2, m_dummy;  
	  
	  if (is_Uplane) { 
	    Guess_slopes( y_min,x, m_dummy,m1 ); 
	    Guess_slopes( y_max,x, m_dummy,m2 ); 
	  } else         { 
	    Guess_slopes( x,y_min, m1,m_dummy ); 
	    Guess_slopes( x,y_max, m2,m_dummy ); 
	  } 
	  
	  double Eta = Grid_search( evt, group, 
				    m1,m2, 
				    x -gridSpacing, 
				    x +gridSpacing, 
				    9e-9, 
				    20e-9 ); 
	  
	  if (Eta > CUT_minEta) 
	    clust_all.push_back( new THitCluster( group, x, Eta ) ); 
	  
	}//for (int ix=0; ix<span; ix++) 
      }//for (int g=0; g<groups.size(); g++) 
      
      
      //no clusters found! 
      if (clust_all.size()<1) return clust_all; 
      
      //now, make sure that that clusters aren't bunched up together
      ROOT::RVec<THitCluster*> clust_keep; 
      
      //now, prune the clusters to decide which hits to get rid of
      auto bestClust = clust_all.at(0); 
            
      double prev_x   = bestClust->Intercept(); 
      
      for (int c=0; c<clust_all.size(); c++) { 
		
	auto clust = clust_all.at(c); 
	
	//how far is this cluster from the last one? 
	int gap = (int)TMath::Nint( (clust->Intercept()-prev_x)/gridSpacing ); 
	
#if 0 //DEBUG	
	cout << TString::Format("c:%3i, int:%2.4f, eta:%1.3f, gap:%1.4f", 
				c, 
				clust->Intercept(), 
				clust->Eta(), 
				(clust->Intercept()-prev_x)/gridSpacing ) << endl;
#endif

	if (gap > min_cluster_gap) { 
	  
	  clust_keep.push_back( bestClust ); 
	  
	  bestClust = clust; //start a new potential cluster grouping
	  
	} else  { //this cluster is close to the last one
	  
	  if (clust->Eta() > bestClust->Eta()) //the new best cluster
	    bestClust = clust; 
	}
	prev_x = clust->Intercept(); 
      }
      clust_keep.push_back( bestClust ); 

#if 0 //DEBUG 
      cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~ kept ("<<clust_keep.size()<<"):"<<endl; 
      
      for (int c=0; c<clust_keep.size(); c++) 
	cout << TString::Format("c:%3i, int:%2.4f, eta:%1.3f", 
				c, 
				clust_keep.at(c)->Intercept(), 
				clust_keep.at(c)->Eta() ) << endl; 
      
      cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl; 
#endif 
      
      return clust_keep; 
    }; 
    
    //for each possible pairing of groups, check to see if they
    // might be allowed to form a TChamberPair, and, if so, go ahead and generate it. 
    ROOT::RVec<TChamberPair*> pairs; 
    
    ROOT::RVec<THitCluster*> clusters_u 
    = Find_clusters( false, 
		     gVec_V ); 
    
    ROOT::RVec<THitCluster*> clusters_v
    = Find_clusters( true, 
		     gVec_U ); 
    

    //now, see if any of these clusters might be valid pairs
    for (int cv=0; cv<clusters_v.size(); cv++) {    
      auto clust_v = clusters_v.at(cv); 
      
      for (int cu=0; cu<clusters_u.size(); cu++) {  
	auto clust_u = clusters_u.at(cu); 
	
	double vu_diff = clust_v->Intercept() - clust_u->Intercept(); 
	
	if ( vu_diff + 2.*gridSpacing < CUT_vu_diff_min ||
	     vu_diff - 2.*gridSpacing > CUT_vu_diff_max  ) continue; 
	
	//now check to see if these clusters ACTUALLY match
	
	double m_u,m_v; 
	Guess_slopes( clust_u->Intercept(), 
		      clust_v->Intercept(), m_u, m_v );
	
	double Eta_u = Grid_search( evt, clust_u->GetGroup(), 
				    m_u,m_u, 
				    clust_u->Intercept() -gridSpacing*2., 
				    clust_u->Intercept() +gridSpacing*2., 
				    9e-9,    //TAU_sigma
				    20e-9 ); //TAU_buffer 
	
	if (Eta_u < CUT_minEta) continue; 
	
	double Eta_v = Grid_search( evt, clust_v->GetGroup(), 
				    m_v,m_v, 
				    clust_v->Intercept() -gridSpacing*2., 
				    clust_v->Intercept() +gridSpacing*2., 
				    9e-9,    //TAU_sigma
				    20e-9 ); //TAU_buffer 
	
	if (Eta_v < CUT_minEta) continue; 
	
	pair_unique_id++; //this will be used by tracks later
	
	pairs.push_back( new TChamberPair( is_LoChamber, clust_u,clust_v, 
					   pair_unique_id++ ) ); 
      }
    }
    
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
    

#define ANGLE_CUT_OFFSET 25e-3 //25e-3 

  
  const double CUT_ph_min = -0.018   -ANGLE_CUT_OFFSET; 
  const double CUT_ph_max =  0.008   +ANGLE_CUT_OFFSET; 
  
  const double CUT_th_min = -0.0200  -ANGLE_CUT_OFFSET; 
  const double CUT_th_max =  0.0175  +ANGLE_CUT_OFFSET; 
  
  
  auto Gen_rawTracks = [&Grid_search, gridSpacing,
			&Phi_model,   CUT_ph_min,  CUT_ph_max,
			&Theta_model, CUT_th_min,  CUT_th_max] 
    ( TEventHandler *evt, 
      ROOT::RVec<TChamberPair*> pairs_Lo, 
      ROOT::RVec<TChamberPair*> pairs_Hi ) { 
    
    bool is_RightArm = evt->ActiveArm(); 
    
    const double CUT_xParam = 2.00; 
    
    ROOT::RVec<TvdcTrack*> tracks; 
    
    //no possible pairs found in the upper chamber!! return the empty vec. 
    for (int pH=0; pH<pairs_Hi.size(); pH++) { 
      for (int pL=0; pL<pairs_Lo.size(); pL++) { 
	
	auto pLo = pairs_Lo.at(pL); 
	auto pHi = pairs_Hi.at(pH); 
	
	auto track = new TvdcTrack( evt, pLo, pHi ); 
	
	double err_Theta = track->Theta() - Theta_model( track ); 
	double err_Phi   = track->Phi()   - Phi_model( track ); 
	
	//if this track has a good angular match, keep it. 
	if ( err_Theta < CUT_th_min ||
	     err_Theta > CUT_th_max ||
	     err_Phi   < CUT_ph_min ||
	     err_Phi   > CUT_ph_max ) { track->~TvdcTrack(); continue; }
	
	//check to make sure each plane passes the min-eta cut
	bool pass_minEtaCut=true; 
	
	double net_eta(0.); 
	
	for (int p=0; (p<4 && pass_minEtaCut); p++) { 
	  
	  double eta = Grid_search( evt, track->GetGroup(p), 
				    track->Slope(p), track->Slope(p), 
				    track->Intercept(p) -gridSpacing*2.,
				    track->Intercept(p) +gridSpacing*2.,
				    9e-9, 25e-9 ); 
	  
	  track->Set_Eta( p, eta ); 
	  
	  if (eta < CUT_minEta) pass_minEtaCut=false;
	}//for (int p=0; p<4; p++) 
	
	if (!pass_minEtaCut) { track->~TvdcTrack(); continue; }
	
	tracks.push_back( track ); 	
	
      }//for (int pL=0; pL<pairs_Lo.size(); pL++) 
    }//for (int pH=0; pH<pairs_Hi.size(); pH++) 
    
    //if we didn't find any tracks, quit
    if (tracks.size()<1) return tracks; 
        
    //check to make sure that tracks aren't 'sharing' clusters
    
    //cout << "Size before pruning = " << tracks.size() << endl; 
    
    auto Delete_shared_tracks = [&tracks](TChamberPair *pair) { 
    
      if (pair->N_tracks() <= 1) return; 
            
      auto best_track = (TvdcTrack*)pair->GetTrack(0);  
      
      //find the track with the highest eta
      for (int t=0; t<pair->N_tracks(); t++) { 
	
	auto new_track = (TvdcTrack*)pair->GetTrack(t); 
	
	if ( new_track->Get_Eta() > best_track->Get_Eta() ) 
	  best_track = new_track; 
      }
      
      //cout << "best eta = " << best_track->Get_Eta() << endl; 
            
      //now, delete all other tracks
      for (int t=0; t<pair->N_tracks();) { 
	
	/*cout << TString::Format("trk:%2i size:%2i  ",
	  t,pair->N_tracks()) << flush; */ 
	
	auto test_track = (TvdcTrack*)pair->GetTrack(t); 
	
	if ( best_track == test_track ) { 
	  //cout << "Best-found!" << endl; 
	  t++; continue; } 
	
	//cout << " removing track..." << flush; 
	
	//note that the call of ~TvdcTrack(); lets our pairs know to remove this 
	// track from their internal vectors. 
	//((TvdcTrack*)pair->GetTrack(t))->~TvdcTrack(); 
	
	test_track->~TvdcTrack(); 
	
	//remove this track from the overall track vector
	tracks.erase( find(begin(tracks), end(tracks), test_track) ); 
      }
    }; 
    
    for (int p=0; p<pairs_Lo.size(); p++) { 
      
      /*cout << TString::Format("p:%2i, nTrks:%2i",
	p, pairs_Lo.at(p)->N_tracks()) << endl; */ 
      
      Delete_shared_tracks( pairs_Lo.at(p) ); 
    }
    
    for (int p=0; p<pairs_Hi.size(); p++) {  
      
      /*cout << TString::Format("p:%2i, nTrks:%2i",
	p, pairs_Hi.at(p)->N_tracks()) << endl; */ 
      
      Delete_shared_tracks( pairs_Hi.at(p) ); 
    }
    
    //cout << "Size after pruning = " << tracks.size() << endl; 
    
    return tracks; 
  }; 
  
  
    //
  // Refine tracks from the hi-chamber using newton's method
  //   
  auto Refine_track = []( TvdcTrack *trk, 
			  const int nCycles=15, 
			  double sigma=25e-9 ) { 
    
    //random-walk track minimizaiton
    const double GRAD_momentum =0.50; 
    const double GRAD_step0    =0.05; 
    const double GRAD_exponent =0.33; 
    
    //cout << "refining track..." << endl; 
    const double scale_T =3e-9;   
    const double scale_X =2.2e-3; 
    
    const double nudge_multiplier =1.; 
    
    const double Chi_cutoff  =6.; 
    
    const double sigma_decay =1.00;
    const double measure_sigma   =9e-9; //sigma*2.; 
    
    ///////////////////////////////////////////////////////////////////////////////
    double int0[5]; 
    for (int p=0; p<4; p++) int0[p] = trk->Intercept(p); 
    TString plane_name[4] = {"U1","V1","U2","V2"}; 
    int0[4] = trk->T0(); 
    //*/////////////////////////////////////////////////////////////////////////////
    
    double objective_eta[4] = {0.}; 
    
    const double s2 = sigma*sigma; 
    
    ROOT::RVec<double> tVec; 
    
    for (UInt_t c=0; c<nCycles; c++) { 
      
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
	  
	  if (c==nCycles-1) 
	    objective_eta[p] 
	      += TMath::Exp( -0.5*TMath::Power(Chi/measure_sigma,2) ); 
	  
	  
	  
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
    
    trk->Set_Eta( objective_eta ); 
    
  };   
  
  
  auto Compute_trackError = [](TvdcTrack *trk) {
    
    const double ERR_tau_sigma = trk->GetEvent()->Get_tauSigma(); 
    
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
  
    
  
  auto Compute_trackData = [](TvdcTrack *trk) { 
    
    //good-point group
    const double CUT_goodPoint = 40e-9; 
    const double measure_sigma = TRK_measureSigma; 
    
    for (int p=0; p<4; p++) { 
      
      auto group = trk->GetGroup(p); 
      
      double eta(0.); 
      double RMS(0.); 
      int    nGoodPoints(0); 
      
      for (int h=0; h<group->Nhits(); h++) { 
	
	double err 
	  = TMath::Abs( trk->Get_T_model(p,group->WirePos(h)) + trk->T0()
			- group->Time(h) ); 
	
	if (err > CUT_goodPoint) continue; 
	
	//now, lets use this point (it's good, according to our cut)
	nGoodPoints++; 
	
	eta += TMath::Exp( -0.5*TMath::Power(err/measure_sigma, 2) ); 
	
	RMS += err*err; 
      }
      
      trk->Set_Eta( p, eta ); 
      
      trk->Set_nGoodPoints( p, nGoodPoints ); 
      
      if (nGoodPoints > 0) { 
	
	trk->Set_RMS( p, TMath::Sqrt( RMS/((double)nGoodPoints) ) ); 
	
      } else { trk->Set_RMS( p, -1e30 ); }
      
    }//for (int p=0; p<4; p++) 
    
    return;  
  };
    
  
  
  string branch_event = (string)TString("event_"+arm); 
  
  vector<string> branch_group; 
  
  //Now that we have eliminated non-coinc events, we can proceed with the 
  // analysis. (Starting with the right arm)
  inNode = inNode
    .Define(branch_event, [is_RightArm](TEventHandler *event) 
	    { event->SetActiveArm(is_RightArm); return event; }, {"event"}); 
  
  for (int p=0; p<4; p++) { 
    
    branch_group.push_back( (string)"groups_"+arm+"_"+plane_name[p] ); 
    
    inNode = inNode
      
      .Define(branch_group[p].data(), [&group_hits,p](TEventHandler *evt, 
						      RVecD wire, RVecD time)
	      { return group_hits(evt,p,wire,time); },
	      {branch_event, branch_wire[p].data(), branch_rawtime[p].data()});
  }
  
  
  //genrate groups
  auto nEvents_1group = inNode  
    .Filter([](ROOT::RVec<THitGroup*> v)
	    {return v.size()>0;}, {branch_group[0].data()})
    .Filter([](ROOT::RVec<THitGroup*> v)
	    {return v.size()>0;}, {branch_group[1].data()})
    .Filter([](ROOT::RVec<THitGroup*> v)
	    {return v.size()>0;}, {branch_group[2].data()})
    .Filter([](ROOT::RVec<THitGroup*> v)
	    {return v.size()>0;}, {branch_group[3].data()}); 
  
  
  //generate pairs
  auto nEvents_1pair = nEvents_1group
    
    //create lo-chamber pairs
    .Define("pairs_"+arm+"_LoChamber", [&Gen_pairs]
	    (TEventHandler *evt, 
	     ROOT::RVec<THitGroup*> vec_gU, 
	     ROOT::RVec<THitGroup*> vec_gV) 
	    
	    { return Gen_pairs(evt, vec_gU,vec_gV, true); }, 
	    {branch_event, 
		branch_group[0].data(), 
		branch_group[1].data() })
    
    .Filter([](ROOT::RVec<TChamberPair*> v) 
	    { return v.size()>0; }, {"pairs_"+arm+"_LoChamber"}) 
    
    .Define("pairs_"+arm+"_HiChamber", [&Gen_pairs]
	    (TEventHandler *evt, 
	     ROOT::RVec<THitGroup*> vec_gU, 
	     ROOT::RVec<THitGroup*> vec_gV) 
	    
	    { return Gen_pairs(evt, vec_gU,vec_gV, false); }, 
	    {branch_event, 
		branch_group[2].data(), 
		branch_group[3].data() })
    
    .Filter([](ROOT::RVec<TChamberPair*> v) 
	    { return v.size()>0; }, {"pairs_"+arm+"_HiChamber"}); 
    

  auto nEvents_1rawTrack = nEvents_1pair 
    
    //generate raw tracks
    .Define("tracks_"+arm+"_raw", Gen_rawTracks, 
	    {branch_event,"pairs_"+arm+"_LoChamber","pairs_"+arm+"_HiChamber"}) 
    
    .Filter([](ROOT::RVec<TvdcTrack*> v) { return v.size(); }, 
	    {"tracks_"+arm+"_raw"}); 
  
  nPass_1group    = nEvents_1group   .Count(); 
  nPass_1pair     = nEvents_1pair    .Count(); 
  nPass_1rawTrack = nEvents_1rawTrack.Count(); 
  
  return nEvents_1rawTrack; 
}
bool OUTPUT_ARM; 

void vdc_tracking_dataFrame_genPts( TString path_inFile ="replay/out.4329.root", 
				    bool in_OUTPUT_ARM=false,
				    TString path_outFile
				    ="driftModel/data/test_RDF.root", 
				    long int maxEntries=15e3 )
{
  ////////////////////////////////////////
  ////////////////////////////////////////
  //
  vector<string> branch_S2; 
  
  OUTPUT_ARM = in_OUTPUT_ARM; 
    
  bool is_monteCarlo = IS_MONTECARLO; 
  
  
  branch_S2.push_back( is_monteCarlo ? "R_s2_rt" : "R.s2.rt" );
  branch_S2.push_back( is_monteCarlo ? "R_s2_lt" : "R.s2.lt" ); 
  branch_S2.push_back( is_monteCarlo ? "L_s2_rt" : "L.s2.rt" ); 
  branch_S2.push_back( is_monteCarlo ? "L_s2_lt" : "L.s2.lt" );  
  
  
  cout << "searching for both-arm coincidence peak..." << flush; 
  
  /***
   *    Takes Raw S2-m data from the decoded-data tree, and converts them to 
   *     the objects called 'TS2Hit*', stored in an RVec.  
   ***/ 
  const double CUT_twinHit_timeErr =5e-9; 
  
  auto Generate_S2_hits = [CUT_twinHit_timeErr]
    ( bool arm, RVecD PMT_R, RVecD PMT_L ) { 
    
    //generate a vector of all coinc s2-paddle hits
    ROOT::RVec<TS2Hit*> coincHits; 
    
    int last_paddle(-999); 
    
    for (int p=0; p<16; p++) { 
      
      auto hit = new TS2Hit( arm, p, PMT_R[p], PMT_L[p] );
      
      //check if it's a coinc-hit
      if ( !hit->IsCoinc() ) { hit->~TS2Hit(); continue; }
      
      if (coincHits.size()>0) { 
	auto lastHit = coincHits.at(coincHits.size()-1); 
		
	if ( hit->Paddle() - lastHit->Paddle() ==1  &&
	     TMath::Abs(hit->Time() - lastHit->Time()) 
	     < CUT_twinHit_timeErr ) { 
	  
	  //these two hits are 'twins', i.e., likely caused by the same particle 
	  lastHit->Make_twinHit( hit ); 
	  
	  hit->~TS2Hit(); continue; 
	  //his hit & the last one are connected
	} 
      }
      
      coincHits.push_back( hit ); 
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
    = [tR_tL_diff, tR_tL_sigma, CUT_S2Sep_stdDev]
    ( double beamCurrent,  
      int    runNumber, 
      ROOT::RVec<TS2Hit*> hits_R,
      ROOT::RVec<TS2Hit*> hits_L ) {

    ROOT::RVec<TEventHandler*> coincEvents; 
    
    for (int iR=0; iR<hits_R.size(); iR++) { 
      for (int iL=0; iL<hits_L.size(); iL++) { 
	
	double TR_TL = hits_R.at(iR)->Time() - hits_L.at(iL)->Time(); 
	
	if ( TMath::Abs(TR_TL-tR_tL_diff)/tR_tL_sigma < CUT_S2Sep_stdDev ) 
	  
	  //this is a bona-fide both-arm S2-coinc; let's keep it
	  coincEvents.push_back( new TEventHandler(false, 
						   runNumber, 
						   beamCurrent, 
						   hits_R.at(iR), 
						   hits_L.at(iL) ) );  
	
      }//for (int iL=0; iL<hits_L.size(); iL++) { 
    }//for (int iR=0; ...
    
    
    
    return coincEvents;
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
    
  
  const double d_Theta_d_Xfp = 0.16745; 
  
  dP    = 1./((double)err_nSamples); 
  cum_P = dP; 
  
  for (int i=0; i<err_nSamples; i++) { 
    
    z[i] = ROOT::Math::normal_quantile( cum_P-0.5*dP, 1. );  
    cum_P += dP;
  }
  
  
#if !RUN_LIMITED 
  ROOT::EnableImplicitMT(); cout << "running all events..." << endl; 
#else 
  cout << TString::Format("running %i events (single-threading)...", 
			  (int)maxEntries) << endl; 
#endif 
  
  const int CUT_plane_goodPts = 2; 
  
  //now actually perform the analysis chain
  auto nEvents_coinc = d_T
    
#if RUN_LIMITED
    .Range(0,10e3)
#endif 
    
#if IS_MONTECARLO 
    
    //MONTE-CARLO  -- check how many 'good pts' per plane
    .Filter([CUT_plane_goodPts]( ROOT::RVec<int> nPts_good, 
				 ROOT::RVec<int> nPts_all ) { 
	      
	      for (int p=0; p<4; p++) { 
		if (nPts_good.at(p) < CUT_plane_goodPts) 
		  return false; 
		
		/* if (nPts_all.at(p)-nPts_good.at(p) < CUT_plane_goodPts) 
		   return false; */ 
	      }
	      return true; }, {"MONTE_nGoodHits", "MONTE_nAllHits"}) 
#endif 
        
    //find all valid S2 hits for the Rigth arm
    .Define("S2_hits_RHRS", [&Generate_S2_hits](RVecD rt, RVecD lt)
	    { return Generate_S2_hits(true,  rt,lt); }, 
	    { branch_S2[0].data(), branch_S2[1].data() })
    
    //find all valid S2 hits for the Left arm
    .Define("S2_hits_LHRS", [&Generate_S2_hits](RVecD rt, RVecD lt)
	    { return Generate_S2_hits(false, rt,lt); }, 
	    { branch_S2[2].data(), branch_S2[3].data() })
    
    
    //Find all possible both-arm S2 coincidences
    .Define("coincEventVec", Gen_S2_coincHits, 
	    {"hac_bcm_average","fEvtHdr.fRun","S2_hits_RHRS","S2_hits_LHRS"})
    
    //Make a cut on Coincs, we only allow exactly 1-per-ev    
    .Filter([](ROOT::RVec<TEventHandler*> v)
	    { return v.size()==1;}, {"coincEventVec"})
    
    
    //get the event handler (so we don't have to keep passing it as a vec)
    .Define("event", [](ROOT::RVec<TEventHandler*> v) 
	    { 
	      v.at(0)->SetActiveArm(OUTPUT_ARM); //start with the right arm
	      return v.at(0); }, {"coincEventVec"}); 
  
  //these 'RResultPtr's will tell us about 
  ROOT::RDF::RResultPtr<unsigned long long> nPass_1group; 
  ROOT::RDF::RResultPtr<unsigned long long> nPass_1pair; 
  ROOT::RDF::RResultPtr<unsigned long long> nPass_1rawTrack;  
  
  auto track_node = generate_vdcTracks( OUTPUT_ARM, 
					nEvents_coinc, 
					nPass_1group, 
					nPass_1pair, 
					nPass_1rawTrack ); 
  
  auto stopwatch = new TStopwatch; 
    
  TString arm = OUTPUT_ARM ? "R" : "L";   
  
  string branch_rawTrack = (string)TString("tracks_"+arm+"_raw");
  
  vector<string> output_columns = { "w0","vi","tau","int_lo","int_hi","slope" }; 
  
  double CUT_tau_error = 80e-9; 
  
  auto output_node = track_node 
    
    .Define("best_track", [](ROOT::RVec<TvdcTrack*> tracks) { 
	
	auto best_track = tracks.at(0); 

	for (int t=0; t<tracks.size(); t++) {
	  if (tracks.at(t)->Get_Eta() > best_track->Get_Eta()) 
	    best_track = tracks.at(t); 
	}
	return best_track; }, {branch_rawTrack})
    
    .Define("vdc_hits", [CUT_tau_error]
	    (TvdcTrack* trk) {
	      
	      ROOT::RVec<TvdcHit*> hits; 
	      
	      for (int p=0; p<4; p++) { 
		
#if OUTPUT_U_PLANE 
		if ( !(p==0 || p==2)) continue; 
#else 
		if ( !(p==1 || p==3)) continue; 
#endif 
		
		auto group = trk->GetGroup(p);
		
		for (int h=0; h<group->Nhits(); h++) { 
		  
		  double tau_err 
		    = trk->Get_T_model(p,group->WirePos(h))+trk->T0() 
		    - group->Time(h); 
		  
		  if ( TMath::Abs(tau_err) > CUT_tau_error ) continue; 
		  
		  hits.push_back( group->GetHit(h) ); 
		}
	      }
	      
	      return hits; }, {"best_track"}) 
    
    //check to make sure each plane has 3 'good' hits
    .Filter([](ROOT::RVec<TvdcHit*> hits) { 
	
	int pHi(0), pLo(0); 
	
	for (int h=0; h<hits.size(); h++) { 

	  int p = hits.at(h)->Plane(); 
	  
	  if (p>=2) { pHi++; } else { pLo++; }
	  
	}
	if (pHi>=3 && pLo>=3) return true; 
	
	return false; }, {"vdc_hits"})
    
    .Define("w0", [](ROOT::RVec<TvdcHit*> hits) { 
	
	RVecD v; 
	for (int h=0; h<hits.size(); h++) 
	  v.push_back( (hits.at(h)->Plane()<2 ? -wSep()/2. : wSep()/2.) );  
	
	return v; }, {"vdc_hits"}) 
    
    .Define("vi", [](ROOT::RVec<TvdcHit*> hits) { 
	
	RVecD v; 
	for (int h=0; h<hits.size(); h++) 
	  v.push_back( hits.at(h)->wPos() ); 
	
	return v; }, {"vdc_hits"}) 
    
    .Define("tau", [](ROOT::RVec<TvdcHit*> hits) { 
	
	RVecD v; 
	for (int h=0; h<hits.size(); h++) 
	  v.push_back( hits.at(h)->Time() ); 
	
	return v; }, {"vdc_hits"}) 

    .Define("int_lo", [](TvdcTrack* trk) { 
	
	return OUTPUT_U_PLANE 
	? trk->Intercept(0) 
	: trk->Intercept(1); }, {"best_track"}) 
    
    .Define("int_hi", [](TvdcTrack* trk) { 
	
	return OUTPUT_U_PLANE 
	? trk->Intercept(2) 
	: trk->Intercept(3); }, {"best_track"}) 
    
    .Define("slope", [](TvdcTrack* trk) { 
	
	return OUTPUT_U_PLANE 
	? trk->Slope(0) 
	: trk->Slope(1); },     {"best_track"})
    
#if MAKE_OUTPUT 
    
    .Snapshot( "drift_points",
	       path_outFile.Data(), 
	       output_columns ); 
#else 

  ; 
#endif   

  int count_1refTrack = (int)*track_node.Count(); 
  
  int count_coinc     = (int)*nEvents_coinc.Count(); 
  int count_1group    = (int)*nPass_1group; 
  int count_1pair     = (int)*nPass_1pair; 
  int count_1rawTrack = (int)*nPass_1rawTrack; 
  
  double net_time = stopwatch->RealTime(); 
  
  
  cout << TString::Format("nCoinc: %i, passed...", count_coinc) << endl; 
  
  cout << 
    TString::Format(" 1-group:     %6i (%0.4f)", 
		    count_1group, 
		    ((double)count_1group)/((double)count_coinc))     << endl; 
  
  cout << 
    TString::Format(" 1-pair:      %6i (%0.4f)", 
		    count_1pair, 
		    ((double)count_1pair)/((double)count_coinc))      << endl; 
  
  cout << 
    TString::Format(" 1-raw-track: %6i (%0.4f)", 
		    count_1rawTrack, 
		    ((double)count_1rawTrack)/((double)count_coinc))  << endl; 
  
  cout << 
    TString::Format("Time taken: %3.2f (%2.4f ms/rawEvent)",
		    net_time, 
		    1e3*net_time/((double)maxEntries))                << endl; 
  
  
  return; 
  
}
