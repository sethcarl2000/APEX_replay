#include "def_apex.h"
#include "function_lib.C"
#include <fstream>

//#include "class/TvdcHit.h"
//#include "class/TvdcHit.C"
#include "lib/NewTrack.h"
#include <functional>
#include <string>
#include <cmath>

class vdc_tracking_dataFrame;

using namespace std; 

using RVecD  = ROOT::VecOps::RVec<double>;
using RNode  = ROOT::RDF::RNode; 

//path which this code is intended to be executed from
const string PATH_EXECUTE = "/work/halla/apex/disk1/sethhall/analyzer/apex-install/analyzer-APEX/scripts/"; 



//                           //
#define IS_MONTECARLO false  //
#define RUN_LIMITED false    //
#define DEBUG false          //
#define MAKE_OUTPUT true     //
//                           //
//these are used for track error-estimation
double dP; 
double cum_P; 

const int err_nSamples = 25; 

double z[err_nSamples]; 

bool choose_goodOrBad=true;  //true==look at 'good' tracks; false==bad tracks 

double f_model_( double *x, double *par ) { 
  
  auto dummy_s2R = new TS2Hit(true, 0,200,200); 
  auto dummy_s2L = new TS2Hit(false,0,200,200); 
  
  auto dummy_event = new TEventHandler(false, 53., 4329, dummy_s2R, dummy_s2L); 
  
  return dummy_event->Drift_T(x[0], 1.4, 2 );
}  
//#endif 

//track cuts
const double TRK_CUT_Dt     = 40e-9; 
const double TRK_CUT_xParam = 1.55; 
const double TRK_CUT_Eta    = 3.750; 
const double TRK_measureSigma = 9e-9; 
const int    TRK_CUT_nGoodPts_min_perPlane = 1; //good points per plane
const int    TRK_CUT_nGoodPts_min          = 8;   

const double CUT_minEta     = 1.95;//min eta of one plane (during grid-searching). 
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
#if 0 //DEBUG 
    int n_validHits(0); 
    cout << (evt->ActiveArm() ? "arm=Right" : "arm=Left") << endl;
#endif 

    const unsigned int clust_minHits=2; //minimum number of hits in a group
    const int maxGap = 4;      //maximum empty wire-gap in a group
    
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
      
#if 0 //DEBUG 
      cout << TString::Format("wire=%3i, rawtime=%5.1f, time=%5.1f   ",
			      (int)h_wire[h], 
			      h_time[h], 
			      hit->Time()*1e9 ); 
#endif
      
      
      if (hit->Time() > VDC_max_realTime || 
	  hit->Time() < VDC_min_realTime ) {

#if 0 //DEBUG
	cout << "killed! " << endl; 
#endif 
	
	hit->~TvdcHit(); //throw out this hit 
	continue; 
      }
      
#if 0 //DEBUG
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
    
#if 0 //DEBUG 
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
  const double gridSpacing = dWire/10.; 
    
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
      const int   wire_buffer = 30; 
      
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
      
  const double CUT_ph_min = is_RightArm ? -0.012 : -0.018;
  const double CUT_ph_max = is_RightArm ?  0.012 :  0.008;   
  
  const double CUT_th_min = is_RightArm ? -0.018 :-0.020;
  const double CUT_th_max = is_RightArm ?  0.018 : 0.018; 
    
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
			  const int nCycles=10, 
			  double sigma=25e-9 ) { 
    
    
    //check to see if a nudge is 'reasonable' i.e., not 'inf' or 'nan'
    auto Check_nudge = [](const double nudge[5]) { 
      
      const double maxNudge_X = 5e-2; //no nudge should ever be this large
      
      const double maxNudge_T = 3e-7; 
      
      for (int p=0; p<4; p++) { 
	if ( nudge[p] != nudge[p] )            return false; 
	if (TMath::Abs(nudge[p]) > maxNudge_X) return false; 
      }
      
      if (nudge[4] != nudge[4])              return false; 
      if (TMath::Abs(nudge[4]) > maxNudge_T) return false; 
      
      return true; 
    }; 
          
    //random-walk track minimizaiton
    const double GRAD_momentum =0.50; 
    const double GRAD_step0    =0.05; 
    const double GRAD_exponent =0.33; 
    
    //cout << "refining track..." << endl; 
    const double scale_T =3e-9;   
    const double scale_X =2.2e-3; 
    
    const double nudge_multiplier =1.; 
    
    const double Chi_cutoff  =6.; 
    
    const double sigma_decay =0.950;
    const double measure_sigma   =9e-9; //sigma*2.; 
    
    ///////////////////////////////////////////////////////////////////////////////
    double int0[5]; 
    for (int p=0; p<4; p++) int0[p] = trk->Intercept(p); 
    TString plane_name[4] = {"U1","V1","U2","V2"}; 
    int0[4] = trk->T0(); 
    //*/////////////////////////////////////////////////////////////////////////////
    
    double final_eta[4] = {0.}; 
    
    const double s2 = sigma*sigma; 
    
    ROOT::RVec<double> tVec; 
    
    for (UInt_t c=0; c<nCycles; c++) { 
      
      double objective_eta[4] = {0.};
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
      	
      sigma *= sigma_decay; 
      
      for (int p=0; p<4; p++) { 
	
	J_44  +=  J_4i[p]*(-J_4i[p]/J_ii[p]); 
	F[4]  +=  F[p]   *(-J_4i[p]/J_ii[p]); 
      }
      
      nudge[4] = -nudge_multiplier*F[4]/J_44; 
      
      for (int p=0; p<4; p++) { 
	
	nudge[p] = -nudge_multiplier*(F[p] - J_4i[p]*nudge[4])/J_ii[p]; 
      }
      
      //check to see if this nudge is reasonalbe
      if (!Check_nudge(nudge)) break; 
      
      trk->Nudge_params( nudge ); 
      
      trk->Set_Eta( objective_eta ); 
    }
    
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
    
#if DEBUG 
    cout << "Compute_trackData():"; 
    //check if track exists 
    cout << "track = " << (trk ? "exists." : "does not exist!!!!!!!!" ) << endl; 
#endif 

    //good-point group
    const double CUT_goodPoint = 40e-9; 
    const double measure_sigma = TRK_measureSigma; 
    
    for (int p=0; p<4; p++) { 
      
      auto group = trk->GetGroup(p); 
      
      double eta(0.); 
      double RMS(0.); 
      int    nGoodPoints(0); 
      
#if DEBUG
      cout << TString::Format("plane %2i ", p);   
#endif 
      
      for (int h=0; h<group->Nhits(); h++) { 
	
	double err 
	= TMath::Abs( trk->Get_T_model(p,group->WirePos(h)) + trk->T0()
			- group->Time(h) ); 
	
	if (err > CUT_goodPoint) continue; 
	
	//now, lets use this point (it's good, according to our cut)
	nGoodPoints++; 
	
	eta += TMath::Exp( -0.5*TMath::Power(err/measure_sigma, 2) ); 
#if DEBUG 
	cout 
	  << TString::Format("h:%3i  wp=%3.3f t-T=%3.3f t-T0=%3.3f g-Time=%3.3f err=%3.3f",
			     h, 
			     group->WirePos(h), 
			     1e9*trk->Get_T_model(p,group->WirePos(h)),
			     1e9*trk->T0(),
			     1e9*group->Time(h), 
			     1e9*err )            << endl; 
#endif 
	
	RMS += err*err; 
      }
      
#if DEBUG
      if ( eta==eta ) { cout << TString::Format("Eta= %0.3f", eta) << endl; }
      else            { cout << "Eta= NaN !!!!!!!!!!!!!!!!!!!!!!!" << endl; } 
#endif 
      
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
  
  
  const double CUT_goodPoints_error = 40e-9; 
  
  auto nEvents_1refinedTrack = nEvents_1rawTrack
    
    //refine track candidates
    .Define("tracks_"+arm+"_refined", [ &Compute_trackError, 
					&Refine_track, 
					&Theta_model,  CUT_th_min, CUT_th_max, 
					&Phi_model,    CUT_ph_min, CUT_ph_max, 
					&Compute_trackData ]
	    ( ROOT::RVec<TvdcTrack*> tracks ) { 
	      
	      ROOT::RVec<TvdcTrack*> refined_tracks; 

	      const double tau_sigma = tracks.at(0)->GetEvent()->Get_tauSigma(); 
	      
	      for (int t=0; t<tracks.size(); t++) {
		
		auto trk = tracks.at(t); 
		
		Refine_track(trk, 20, 25e-9); 
		
#if DEBUG
		bool is_nan=false; 
		
		for (int p=0; p<4; p++) 
		  if ( trk->Intercept(0)!=trk->Intercept(0) ) is_nan=true; 
		
		if ( trk->T0()!=trk->T0() ) is_nan=true; 
		
		if (is_nan) { cout << "NAN intercept!!!!!!!!!!!" << endl; }
		else        { cout << "intercept exists." << endl; }
		
		cout << " intercepts = { " ; 
		for (int p=0; p<4; p++) { 
		  cout << TString::Format("%0.3f ", trk->Intercept(p)); 
		} cout << "}" << endl; 
		cout << TString::Format( "Theta=%0.3f, Phi=%0.3f", 
					 trk->Theta(), trk->Phi() ) << endl; 
		
#endif 
		
		double err_Theta = trk->Theta() - Theta_model(trk); 
		double err_Phi   = trk->Phi()   - Phi_model(trk); 
				
		if ( err_Theta < CUT_th_min || 
		     err_Theta > CUT_th_max || 
		     err_Phi   < CUT_ph_min || 
		     err_Phi   > CUT_ph_max ) { //cut this track 
		  
		  //delete this track 		  
		  trk->~TvdcTrack(); continue; 
		} 
		
		//find out how many 'good' points each plane has for this track
		Compute_trackData( trk ); 
		
		//do some basic checks
		double xParam = trk->xParam();
		double Dt     = trk->T0(); 
				
		if (trk->Get_Eta() < TRK_CUT_Eta || 
		    trk->Get_nGoodPoints(0) < TRK_CUT_nGoodPts_min_perPlane || 
		    trk->Get_nGoodPoints(1) < TRK_CUT_nGoodPts_min_perPlane || 
		    trk->Get_nGoodPoints(2) < TRK_CUT_nGoodPts_min_perPlane || 
		    trk->Get_nGoodPoints(3) < TRK_CUT_nGoodPts_min_perPlane || 
		    trk->Get_nGoodPoints()  < TRK_CUT_nGoodPts_min ||
		    TMath::Abs(xParam) > TRK_CUT_xParam ||
		    TMath::Abs(Dt)     > TRK_CUT_Dt        ) { 
		  
		  //delete this track 		  
		  trk->~TvdcTrack(); continue; 
		} 
		
		//compute track error
		Compute_trackError( trk ); 
		
		//maybe some cuts here for the track error? 
		
		refined_tracks.push_back( trk ); 
	      }
	      return refined_tracks; }, {"tracks_"+arm+"_raw"})
    
    .Filter([](ROOT::RVec<TvdcTrack*> v) 
	    { return v.size()>0; }, {"tracks_"+arm+"_refined"}); 
  
  
  nPass_1group    = nEvents_1group   .Count(); 
  nPass_1pair     = nEvents_1pair    .Count(); 
  nPass_1rawTrack = nEvents_1rawTrack.Count(); 
  
  return nEvents_1refinedTrack; 
}
bool is_monteCarlo; 

RNode generate_track_output( bool is_RightArm, RNode track_node, 
			     vector<string> &out_branches ) { 

  string branch_tr = is_RightArm ? "tracks_R_refined" : "tracks_L_refined"; 
  
  TString arm = is_RightArm ? "R" : "L"; 
  
  out_branches.push_back( (string)(arm+"_tr_fp_x") );
  
  track_node = track_node  
    .Define( out_branches.at(out_branches.size()-1),  
	     [](ROOT::RVec<TvdcTrack*> tr) 
	     { RVecD v; 
	       for (int t=0; t<tr.size(); t++) 
		 v.push_back( tr.at(t)->FP_x() );
	       return v; }, {branch_tr}); 
    
  
  out_branches.push_back( (string)(arm+"_tr_fp_y") );
  
  track_node = track_node
    .Define( out_branches.at(out_branches.size()-1),  
	     [](ROOT::RVec<TvdcTrack*> tr) 
	     { RVecD v; 
	       for (int t=0; t<tr.size(); t++) 
		 v.push_back( tr.at(t)->FP_y() ); 
	       return v; }, {branch_tr}); 


  out_branches.push_back( (string)(arm+"_tr_theta") );
  
  track_node = track_node
    .Define( out_branches.at(out_branches.size()-1),  
	     [](ROOT::RVec<TvdcTrack*> tr) 
	     { RVecD v; 
	       for (int t=0; t<tr.size(); t++) 
		 v.push_back( tr.at(t)->Theta() ); 
	       return v; }, {branch_tr}); 
  
  
  out_branches.push_back( (string)(arm+"_tr_phi") );
  
  track_node = track_node 
    .Define( out_branches.at(out_branches.size()-1),  
	     [](ROOT::RVec<TvdcTrack*> tr) 
	     { RVecD v; 
	       for (int t=0; t<tr.size(); t++) 
		 v.push_back( tr.at(t)->Phi() ); 
	       return v; }, {branch_tr}); 
  
  
  string plane_name[4] = { "U1","V1","U2","V2" }; 
  
  const double sigma_eta = 9e-9; 
  const double sigma_cutoff = 40e-9; 
  
  for (int p=0; p<4; p++) { 
    
    out_branches.push_back( (string)(arm+"_tr_"+plane_name[p]+"_intercept" ) );
    out_branches.push_back( (string)(arm+"_tr_"+plane_name[p]+"_intercept_error" ) );
    out_branches.push_back( (string)(arm+"_tr_"+plane_name[p]+"_nGoodPoints" ) );
    out_branches.push_back( (string)(arm+"_tr_"+plane_name[p]+"_Eta" ) );
    
    track_node = track_node 
      
      .Define( (string)arm+"_tr_"+plane_name[p]+"_intercept", 
	      [p](ROOT::RVec<TvdcTrack*> tr) 
	      { RVecD v; 
		for (int t=0; t<tr.size(); t++) 
		  v.push_back( tr.at(t)->Intercept(p) ); 
		return v; }, {branch_tr})
      
      .Define( (string)arm+"_tr_"+plane_name[p]+"_intercept_error", 
	      [p](ROOT::RVec<TvdcTrack*> tr) 
	      { RVecD v; 
		for (int t=0; t<tr.size(); t++) 
		  v.push_back( tr.at(t)->Error_intercept(p) ); 
		return v; }, {branch_tr})
      
      .Define( (string)arm+"_tr_"+plane_name[p]+"_nGoodPoints", 
	      [p](ROOT::RVec<TvdcTrack*> tr) 
	      { ROOT::RVec<int> v; 
		for (int t=0; t<tr.size(); t++) 
		  v.push_back( tr.at(t)->Get_nGoodPoints(p) ); 
		return v; }, {branch_tr})
      
      .Define( (string)arm+"_tr_"+plane_name[p]+"_Eta", 
	      [p, sigma_eta](ROOT::RVec<TvdcTrack*> tr) 
	      { RVecD v; 
		for (int t=0; t<tr.size(); t++) 
		  v.push_back( tr.at(t)->Get_Eta(p) ); 
		return v; }, {branch_tr}); 
    
  }
  
  out_branches.push_back( (string)(arm+"_tr_Eta") );
  
  track_node = track_node 
    .Define( out_branches.at(out_branches.size()-1),  
	     [](ROOT::RVec<TvdcTrack*> tr) 
	     { RVecD v; 
	       for (int t=0; t<tr.size(); t++) 
		 v.push_back( tr.at(t)->Get_Eta() ); 
	       return v; }, {branch_tr}); 
  
  
  out_branches.push_back( (string)(arm+"_tr_RMS") );
  
  track_node = track_node 
    .Define( out_branches.at(out_branches.size()-1),  
	     [](ROOT::RVec<TvdcTrack*> tr) 
	     { RVecD v; 
	       for (int t=0; t<tr.size(); t++) 
		 v.push_back( tr.at(t)->Get_RMS() ); 
	       return v; }, {branch_tr}); 
  
  
  out_branches.push_back( (string)(arm+"_tr_nGoodPoints") );
  
  track_node = track_node 
    .Define( out_branches.at(out_branches.size()-1),  
	     [](ROOT::RVec<TvdcTrack*> tr) 
	     { RVecD v; 
	       for (int t=0; t<tr.size(); t++) 
		 v.push_back( tr.at(t)->Get_nGoodPoints() ); 
	       return v; }, {branch_tr}); 

  
  return track_node; 
}

bool isRightArm; 

const char* PATH_VOLATILE_DIR = {"/volatile/halla/apex/full_replay"}; 

void vdc_tracking_dataFrame_fork( int maxEntries=5e4, 
				  bool IN_isRightArm=true, //true=RHRS, false=LHRS 
				  int runNum=4712,
				  TString path_inFile
				  ="$VOLATILE_DIR/production/decode.4329.root",
				  TString path_outFile="" )
{
  ////////////////////////////////////////
  ////////////////////////////////////////
  //
  vector<string> branch_S2; 
  
  
  
  isRightArm = IN_isRightArm; 
  
  TString STRING_replay
         ="   ";
  
  cout << "vdc_tracking_dataFrame_fork.C -- executing replay script.~~~~~~~~~~~~~~~~~~~~~~~" 
       << endl; 
  
  cout << "   input file:  \"" <<path_inFile<< "\"" << endl; 
  cout << "   output file: \"" <<path_outFile<< "\"" << endl; 
  cout << "   test-parameters (ignore if full replay):\n" 
       << TString::Format( "   arm=%s, run=%i", 
			   (isRightArm ? "R-hrs":"L-hrs"),
			   runNum ) << endl; 
  
  cout << "   monte-carlo run? " 
       << (IS_MONTECARLO ? "true" : "false") << endl; 
  
  cout << "   debug active?    " 
       << (DEBUG ? "true" : "false") << endl; 
  
  cout << "   run all events?  "; 
  if (RUN_LIMITED==true) { 
    
    cout << TString::Format( "false ( max raw-events to run: %0.2e )",
			     (double)maxEntries ) << endl; 
  } else            { 
    cout << "true" << endl; 
  }
  
  cout << "   create output?   " 
       << (MAKE_OUTPUT ? "true" : "false") << endl; 
  
  
  //print some track-cut parameters
  cout << "track-cut parameters ******************************************"<<endl; 
  
  cout << TString::Format("  cut => |T_trk-T_coinc| ......... <  %0.2f ns", 
			  TRK_CUT_Dt*1e9 ) << endl; 
  
  cout << TString::Format("  cut => |xParam| ................ <  %0.2f", 
			  TRK_CUT_xParam ) << endl; 
    
  cout << TString::Format("  cut => track Eta ............... >  %0.3f (measure sigma=%0.3f ns)", 
			  TRK_CUT_Eta, TRK_measureSigma*1e9 ) << endl; 
  
  cout << TString::Format("  cut => GoodPts/track ........... >= %i", 
			  TRK_CUT_nGoodPts_min ) << endl; 
  
  cout << TString::Format("  cut => GoodPts/plane ........... >= %i", 
			  TRK_CUT_nGoodPts_min_perPlane ) << endl; 
  
  cout << TString::Format("  cut => Eta/plane (grid-search) . >  %0.3f", 
			  CUT_minEta ) << endl; 
    
  cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" 
       << endl; 
  
  is_monteCarlo = IS_MONTECARLO; 
  
  branch_S2.push_back( is_monteCarlo ? "R_s2_rt" : "R.s2.rt" );
  branch_S2.push_back( is_monteCarlo ? "R_s2_lt" : "R.s2.lt" ); 
  branch_S2.push_back( is_monteCarlo ? "L_s2_rt" : "L.s2.rt" ); 
  branch_S2.push_back( is_monteCarlo ? "L_s2_lt" : "L.s2.lt" );  
  
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
  
  
  ROOT::RDataFrame d_T_S2coinc("T", path_inFile.Data() ); 
  
  auto nTotalEvents = d_T_S2coinc.Count(); 
  
  const int min_events    = 1500; 
  const int min_coincHits = 3000; 
  
  if ( (int)*nTotalEvents < min_events ) { 
    
    cout << "ERROR! too few events to process this decode-file.\n"
	 << TString::Format(" events %i (min events: %i)\n", 
			    (int)*nTotalEvents, 
			    min_events ) 
	 << " quitting replay script..." << endl; 
    return; 
    
  } else { cout << "Total number of events: " << (int)*nTotalEvents << endl; }
  
  
  cout << "searching for both-arm coincidence peak..." << flush; 
  
  int nEvents_coincTest = (int)TMath::Min( (double)*nTotalEvents, 30e3 ); 
  
  auto h_dt_temp = d_T_S2coinc
    .Range(0,nEvents_coincTest)
    .Define("dt", get_dt, branch_S2)
    .Histo1D({"", "dt", 100, -60e-9, 60e-9}, {"dt"}); 
  

  if (h_dt_temp->Integral() < min_coincHits) { 
    
    cout << "ERROR! too few S2-hits in both arms in decode-file.\n"
	 << " perhaps one arm is inactive?\n"
	 << TString::Format(" S2-hits %i (min events: %i)\n", 
			    (int)h_dt_temp->Integral(), 
			    min_events ) 
	 << " quitting replay script..." << endl; 
    return; 
  }
    
    
    
  
  //h_dt_temp->DrawCopy(); 
  
  TVector2 dt_fitResult = draw_gausFit( (TH1D*)h_dt_temp->Clone(), 
					6e-9, -1e30, -1 ); 
  
  const double tR_tL_diff  = dt_fitResult.X(); 
  const double tR_tL_sigma = dt_fitResult.Y(); 

  const double CUT_S2Sep_stdDev = 6; 
  
  cout << TString::Format("TR-TL avg.= %0.2f  sigma= %0.2f", 
			  tR_tL_diff*1e9, 
			  tR_tL_sigma*1e9) << endl; 
			  
  
  
  /*** 
   *     Create a TEventHandler object for each both-arm S2-coinc, 
   *      using the variables above (tR_tL_diff/sigma) to make cuts
   ***/
  auto Gen_S2_coincHits 
    = [tR_tL_diff, tR_tL_sigma, CUT_S2Sep_stdDev]
    ( double beamCurrent,  
      int runNumber, 
      ROOT::RVec<TS2Hit*> hits_R,
      ROOT::RVec<TS2Hit*> hits_L ) {

    ROOT::RVec<TEventHandler*> coincEvents; 
    
    for (int iR=0; iR<hits_R.size(); iR++) { 
      for (int iL=0; iL<hits_L.size(); iL++) { 
	
	double TR_TL = hits_R.at(iR)->Time() - hits_L.at(iL)->Time(); 
	
	if ( TMath::Abs(TR_TL-tR_tL_diff)/tR_tL_sigma < CUT_S2Sep_stdDev ) 
	  
	  //this is a bona-fide both-arm S2-coinc; let's keep it
	  coincEvents.push_back( new TEventHandler(false, 
						   beamCurrent, 
						   runNumber, 
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
  
    
  
#if RUN_LIMITED   
  cout << TString::Format("running %i events (single-threading)...", 
			  maxEntries) << endl; 
#else 
  ROOT::EnableImplicitMT(); 
  cout << "running all events (mulit-threadding enabled)..." << endl; 
#endif
  
  
  ROOT::RDataFrame d_T("T", path_inFile.Data() ); 
    
  
  const int CUT_plane_goodPts = 2; 
    
  //now actually perform the analysis chain
  auto nEvents_coinc = d_T
    
#if RUN_LIMITED
    .Range(0,maxEntries)
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
	    { return v.at(0); }, {"coincEventVec"}); 
  
  
  //these 'RResultPtr's will tell us about 
  ROOT::RDF::RResultPtr<unsigned long long> R_nPass_1group; 
  ROOT::RDF::RResultPtr<unsigned long long> R_nPass_1pair; 
  ROOT::RDF::RResultPtr<unsigned long long> R_nPass_1rawTrack;  
    
  ROOT::RDF::RResultPtr<unsigned long long> L_nPass_1group; 
  ROOT::RDF::RResultPtr<unsigned long long> L_nPass_1pair; 
  ROOT::RDF::RResultPtr<unsigned long long> L_nPass_1rawTrack;  
  

  bool node_arm = isRightArm ? kArm_Right : kArm_Left; 
  
  auto trackNode_right = generate_vdcTracks( kArm_Right, 
					     nEvents_coinc, 
					     R_nPass_1group, 
					     R_nPass_1pair, 
					     R_nPass_1rawTrack ); 
  
  //make sure there's at least 1 refined-track in the right arm 
  trackNode_right = trackNode_right
    
    .Filter([](ROOT::RVec<TvdcTrack*> tracks) 
	    { return tracks.size()>0; }, {"tracks_R_refined"}); 
  
  
  auto trackNode_left = generate_vdcTracks( kArm_Left, 
					    trackNode_right, 
					    L_nPass_1group, 
					    L_nPass_1pair, 
					    L_nPass_1rawTrack ); 
  
  //make sure the left-arm found 1 track
  trackNode_left = trackNode_left 
    
    .Filter([](ROOT::RVec<TvdcTrack*> tracks) 
	    { return tracks.size()>0; }, {"tracks_L_refined"}); 
  
  
  vector<string> out_branches; 
  out_branches.push_back( "L_cerenkov_a_c" ); //.......Left HRS info
  out_branches.push_back( "L_cerenkov_t_c" ); 
  out_branches.push_back( "L_preShower_a_c" );  
  out_branches.push_back( "L_shower_a_c" ); 
  out_branches.push_back( "L_S2_paddle_realTime" );
  out_branches.push_back( "L_S2_paddle_num" );
  out_branches.push_back( "L_S2_isTwinHit" ); 
  out_branches.push_back( "L_S2_pmtDiff_raw" ); 
  out_branches.push_back( "R_cerenkov_a_c" ); //.......Right HRS info
  out_branches.push_back( "R_cerenkov_t_c" ); 
  out_branches.push_back( "R_preShower_a_c" );  
  out_branches.push_back( "R_shower_a_c" ); 
  out_branches.push_back( "R_S2_paddle_realTime" );
  out_branches.push_back( "R_S2_paddle_num" );
  out_branches.push_back( "R_S2_isTwinHit" ); 
  out_branches.push_back( "R_S2_pmtDiff_raw" ); 
  out_branches.push_back( "beamPos_BPMA" );//..........BPMA 
  out_branches.push_back( "beamPos_BPMB" );//..........BPMB
  out_branches.push_back( "beam_current" );//..........beam current 
  out_branches.push_back( "run_number" );//............run number 
  out_branches.push_back( "TR_TL_diff" ); 
  out_branches.push_back( "fEvtHdr.fEvtTime" );//......event 'time' 
  
  //define some output branches for the left-arm
  auto output_node = trackNode_left 
    
    //left arm 
    .Define("L_S2_paddle_realTime", [](TEventHandler *event) 
	    { return event->GetS2Hit()->Time(); }, {"event_L"}) 
    .Define("L_S2_paddle_num", [](TEventHandler *event) 
	    { return event->GetS2Hit()->Paddle(); }, {"event_L"}) 
    .Define("L_S2_pmtDiff_raw", [](TEventHandler *event) 
	    { return event->GetS2Hit()->DeltaT_raw(); }, {"event_L"}) 
    .Define("L_S2_isTwinHit", [](TEventHandler *event) 
	    { return event->GetS2Hit()->Is_twinHit(); }, {"event_L"}) 
    
    .Define("L_cerenkov_t_c",  [](RVecD v){ return v; }, {"L.cer.t_c"}) 
    .Define("L_cerenkov_a_c",  [](RVecD v){ return v; }, {"L.cer.a_c"}) 
    .Define("L_preShower_a_c", [](RVecD v){ return v; }, {"L.prl1.a_c"}) 
    .Define("L_shower_a_c",    [](RVecD v){ return v; }, {"L.prl2.a_c"}) 
     
    //right arm 
    .Define("R_S2_paddle_realTime", [](TEventHandler *event) 
	    { return event->GetS2Hit()->Time(); }, {"event_R"}) 
    .Define("R_S2_paddle_num", [](TEventHandler *event) 
	    { return event->GetS2Hit()->Paddle(); }, {"event_R"}) 
    .Define("R_S2_pmtDiff_raw", [](TEventHandler *event) 
	    { return event->GetS2Hit()->DeltaT_raw(); }, {"event_R"}) 
    .Define("R_S2_isTwinHit", [](TEventHandler *event) 
	    { return event->GetS2Hit()->Is_twinHit(); }, {"event_R"})
    
    .Define("R_cerenkov_t_c",  [](RVecD v){ return v; }, {"R.cer.t_c"}) 
    .Define("R_cerenkov_a_c",  [](RVecD v){ return v; }, {"R.cer.a_c"}) 
    .Define("R_preShower_a_c", [](RVecD v){ return v; }, {"R.ps.a_c"}) 
    .Define("R_shower_a_c",    [](RVecD v){ return v; }, {"R.sh.a_c"})
    
    .Define("beam_current",    [](double x){ return x>-1e5 ? x : -1e30; }, 
	    {"hac_bcm_average"})
    
    .Define("run_number",      [](int    x){ return x; }, {"fEvtHdr.fRun"}) 
    .Define("event_number",    [](UInt_t x){ return x; }, {"fEvtHdr.fEvtNum"})
    .Define("TR_TL_diff",      [tR_tL_diff](){ return tR_tL_diff; })
    
    .Define("beamPos_BPMA",   [](double x,double y,double z) 
	    { return TVector3(x,y,z); },{"Rrb.BPMA.x","Rrb.BPMA.y","Rrb.BPMA.z"}) 
    .Define("beamPos_BPMB",   [](double x,double y,double z) 
	    { return TVector3(x,y,z); },{"Rrb.BPMB.x","Rrb.BPMB.y","Rrb.BPMB.z"});  
  
  
  //get track output branches
  output_node = generate_track_output( kArm_Right, output_node, out_branches ); 
    
  output_node = generate_track_output( kArm_Left,  output_node, out_branches ); 
  
  
  auto stopwatch = new TStopwatch; 
  
  //number of null beam current events
  auto nonNullBC = d_T
    .Filter([](double bcm) { return bcm>-1e29; }, {"hac_bcm_average"}); 
  
  auto h1d_beamCurrent = nonNullBC
    .Histo1D({"beamCurrent", "", 200, 0, 60}, "hac_bcm_average"); 
  
  
  //print output data to a file
#if MAKE_OUTPUT
  output_node.Snapshot( "track_data", path_outFile.Data(), out_branches ); 
#endif 

  int count_nonNullBC  = (int)*nonNullBC.Count(); 
    
  int Lcount_1refTrack = (int)*trackNode_left.Count(); 
  
  int Lcount_1group    = (int)*L_nPass_1group; 
  int Lcount_1pair     = (int)*L_nPass_1pair; 
  int Lcount_1rawTrack = (int)*L_nPass_1rawTrack; 
  
  
  int Rcount_1refTrack = (int)*trackNode_right.Count(); 
  
  int Rcount_1group    = (int)*R_nPass_1group; 
  int Rcount_1pair     = (int)*R_nPass_1pair; 
  int Rcount_1rawTrack = (int)*R_nPass_1rawTrack; 
  
  
  int count_coinc      = (int)*nEvents_coinc.Count(); 
  
  double net_time = stopwatch->RealTime(); 
    
  
  cout << TString::Format("nCoinc: %i, passed...", count_coinc) << endl; 
  
  cout << "Right arm (fractions relative to total coinc count) ~~~~~~~~~~~~~~~~~~~" 
       << endl; 
  
  cout << 
    TString::Format(" 1-group:     %7i (%0.4f)", 
		    Rcount_1group, 
		    ((double)Rcount_1group)/((double)count_coinc))    << endl; 
  
  cout << 
    TString::Format(" 1-pair:      %7i (%0.4f)", 
		    Rcount_1pair, 
		    ((double)Rcount_1pair)/((double)count_coinc))     << endl; 
  
  cout << 
    TString::Format(" 1-raw-track: %7i (%0.4f)", 
		    Rcount_1rawTrack, 
		    ((double)Rcount_1rawTrack)/((double)count_coinc)) << endl; 
  
  cout << 
    TString::Format(" 1-ref-track: %7i (%0.4f)", 
		    Rcount_1refTrack, 
		    ((double)Rcount_1refTrack)/((double)count_coinc)) << endl
       << endl; 
  
  cout << "Left arm (fractions relative to R-events with track) ~~~~~~~~~~~~~~~~~" 
       << endl; 
  
  
  cout << 
    TString::Format(" 1-group:     %7i (%0.4f)", 
		    Lcount_1group, 
		    ((double)Lcount_1group)/((double)Rcount_1refTrack))    << endl; 
  
  cout << 
    TString::Format(" 1-pair:      %7i (%0.4f)", 
		    Lcount_1pair, 
		    ((double)Lcount_1pair)/((double)Rcount_1refTrack))     << endl; 
  
  cout << 
    TString::Format(" 1-raw-track: %7i (%0.4f)", 
		    Lcount_1rawTrack, 
		    ((double)Lcount_1rawTrack)/((double)Rcount_1refTrack)) << endl; 
  
  cout << 
    TString::Format(" 1-ref-track: %7i (%0.4f)", 
		    Lcount_1refTrack, 
		    ((double)Lcount_1refTrack)/((double)Rcount_1refTrack)) << endl
       << endl; 
  

  
  //print beam-current information
  cout << 
    TString::Format("--Average beam current: %2.2f \n--(fraction of events with non-null beam-current = %0.4f)",
		    h1d_beamCurrent->GetMean(), 
		    ((double)count_nonNullBC)/((double)*nTotalEvents) ) << endl; 
  
  
  
  cout << 
    TString::Format("Time taken: %3.2f (%2.4f ms/rawEvent)",
		    net_time, 
		    1e3*net_time/((double)maxEntries)) << endl; 
  
  
  
#if MAKE_OUTPUT 
  //create csv data
  cout << TString::Format("writing data to csv file \"%s\"...", 
			  PATH_VOLATILE_DIR ) << flush; 
  
  fstream out_CSV;
  
  out_CSV.open( TString::Format("%s/csv/info.%i.csv",
				PATH_VOLATILE_DIR, 
				runNum).Data(),     ios::out | ios::app ); 
  

  
  
  out_CSV << TString::Format( "%i,%0.2f,%0.4f,%0.3f,%0.2f,%i,%i,%i,%i,%i,%i",
			      runNum, 
			      h1d_beamCurrent->GetMean(),
			      ((double)count_nonNullBC)/((double)*nTotalEvents),
			      1e9*tR_tL_diff,
			      1e9*tR_tL_sigma*CUT_S2Sep_stdDev, 
			      (int)*nTotalEvents, 
			      count_coinc,   
			      Rcount_1group,
			      Rcount_1refTrack,
			      Lcount_1group,
			      Lcount_1refTrack  ) << endl; 
  out_CSV.close(); 

  cout << "done" << endl; 
  
  return; 
#endif 
  

#if !IS_MONTECARLO 
  
  //non-monte-carlo data-processing
  
  
  
  
  string branch_track = isRightArm ? "tracks_R_refined" : "tracks_L_refined"; 
  
  auto trackNode_test = isRightArm ? trackNode_right : trackNode_left; 
  
  
  auto h2d_eta = trackNode_test 
    
    .Define("eta",         [](ROOT::RVec<TvdcTrack*> tracks) 
	    { RVecD out; 
	      for (int t=0; t<tracks.size(); t++) 
		out.push_back( tracks.at(t)->Get_Eta() ); 
	      
	      return out; }, {branch_track})
    
    .Define("nGoodPoints", [](ROOT::RVec<TvdcTrack*> tracks) 
	    { RVecD out; 
	      for (int t=0; t<tracks.size(); t++) 
		out.push_back( tracks.at(t)->Get_nGoodPoints() ); 
	      
	      return out; }, {branch_track})
        
    .Histo2D({"eta","Eta of all tracks", 29, -0.5, 28.5, 200, 0, 28}, 
	     "nGoodPoints", "eta"); 
  
  
  
  auto h2d_xParam_Dt = trackNode_test
    
    .Define("xParam", [](ROOT::RVec<TvdcTrack*> tracks) 
	    { RVecD out; 
	      for (int t=0; t<tracks.size(); t++) 
		out.push_back( tracks.at(t)->xParam() ); 
	      
	      return out; }, {branch_track}) 
    
    .Define("Dt",      [](ROOT::RVec<TvdcTrack*> tracks) 
	    { RVecD out; 
	      for (int t=0; t<tracks.size(); t++) 
		out.push_back( tracks.at(t)->T0() ); 
	      
	      return out; }, {branch_track})
    
    .Histo2D({"xp_dt", ";x-Param;T_{S2}-T_{Track} (s)", 
	  200, -4.5, 4.5, 200, -100e-9, 100e-9}, "xParam", "Dt"); 
  
  
  
  auto h2d_phi_s2y = trackNode_test 
        
    /*  .Define("theta",      [](ROOT::RVec<TvdcTrack*> tracks) 
	    { RVecD out; 
	      for (int t=0; t<tracks.size(); t++) 
		out.push_back( tracks.at(t)->Theta() ); 
	      
	      return out; }, {branch_track})
    
    .Define("s2_x",      [](ROOT::RVec<TvdcTrack*> tracks) 
	    { RVecD out; 
	      for (int t=0; t<tracks.size(); t++) 
		out.push_back( tracks.at(t)->S2_x() ); 
	      
	      return out; }, {branch_track})
    
    .Histo2D({"theta_s2x", "#Theta vs X_{S2};X_{S2} (m);#Theta (rad)", 
    200, -1.2, 1.2, 200, -0.15, 0.15}, "s2_x", "theta"); */
    
    .Define("phi_err",  [&Phi_model](ROOT::RVec<TvdcTrack*> tracks) 
	    { RVecD out; 
	      for (int t=0; t<tracks.size(); t++) 
		out.push_back( tracks.at(t)->Phi() 
			       - Phi_model( tracks.at(t) ) ); 
	      
	      return out; }, {branch_track})
    
    .Histo1D({"theta_err", "#Phi-model error (rad)", 200, -50e-3, 50e-3}, 
	     "phi_err"); 
  

  
  
  auto h2d_theta_s2x = trackNode_test 
        
    /*.Define("phi",      [](ROOT::RVec<TvdcTrack*> tracks) 
	    { RVecD out; 
	      for (int t=0; t<tracks.size(); t++) 
		out.push_back( tracks.at(t)->Phi() ); 
	      
	      return out; }, {branch_track})
	      
    .Define("s2_y",      [](ROOT::RVec<TvdcTrack*> tracks) 
	    { RVecD out; 
	      for (int t=0; t<tracks.size(); t++) 
		out.push_back( tracks.at(t)->S2_y() ); 
	      
	      return out; }, {branch_track})
    
    .Histo2D({"phi_s2y", "#Phi vs Y_{S2};Y_{S2} (m);#Phi (rad)", 
	200, -0.30, 0.30, 200, -0.15, 0.15}, "s2_y", "phi"); */ 
    
    .Define("theta_err",  [&Theta_model](ROOT::RVec<TvdcTrack*> tracks) 
	    { RVecD out; 
	      for (int t=0; t<tracks.size(); t++) 
		out.push_back( tracks.at(t)->Theta() 
			       - Theta_model( tracks.at(t) ) ); 
	      
	      return out; }, {branch_track})
    
    .Histo1D({"theta_err", "#Theta-model error (rad)", 200, -50e-3, 50e-3}, 
	     "theta_err");
  
  
  
  auto Compute_dv = [](TvdcTrack *trk, int p) { 
    
    bool arm = trk->IsRightArm(); 
    
    int wire       = TvdcHit::WireNum( arm, p, trk->Intercept(p) ); 
    
    double wirePos = TvdcHit::WirePos( arm, p, wire ); 
    
    return trk->Intercept(p) - wirePos; 
  }; 
    
  
  auto h2d_slope_Dt = trackNode_test 
    
    .Define("slope",   [](ROOT::RVec<TvdcTrack*> tracks) 
	    { 
	      RVecD out; 
	      for (int t=0; t<tracks.size(); t++) {
		
		out.push_back( tracks.at(t)->Slope(0) ); 
		out.push_back( tracks.at(t)->Slope(1) ); 

		} return out; }, {branch_track})  
		
    .Define("Dt",  [](ROOT::RVec<TvdcTrack*> tracks) 
	    { 
	      RVecD out; 
	      for (int t=0; t<tracks.size(); t++) {
		
		out.push_back( tracks.at(t)->T0() ); 
		out.push_back( tracks.at(t)->T0() ); 

	      } return out; }, {branch_track})
    
    //  .Histo1D({"Dt", "Dt (s)", 150, -35e-9, 35e-9}, "Dt"); 
    
    .Histo2D({"slope_dt", ";dw/dx (slope);Dt (s)", 
	  60, 1.00, 1.85, 85, -30e-9, 30e-9}, "slope", "Dt");  
  
  
  
  auto h1d_dv = trackNode_test 
    
    .Define("dv",      [&Compute_dv](ROOT::RVec<TvdcTrack*> tracks) 
	    { 
	      RVecD out; 
	      for (int t=0; t<tracks.size(); t++) 
		for (int p=0; p<4; p++) 
		  out.push_back( Compute_dv(tracks.at(t),p) ); 
	      
	      return out; }, {branch_track})
    
    /*.Define("slope", [](ROOT::RVec<TvdcTrack*> tracks) 
	    { 
	      RVecD out; 
	      for (int t=0; t<tracks.size(); t++) 
		for (int p=0; p<4; p++) 
		  out.push_back( tracks.at(t)->Slope(p) ); 
		
	      return out; }, {branch_track})
        
    .Histo2D({"dv", ";dv (m);slope", 
      100, -dWire/2., dWire/2, 100, 1.000, 1.800 }, "dv","slope");   */ 
    
    .Histo1D({"dv", ";dv (m);", 100, -dWire/2, dWire/2}, "dv"); 
	     
  
  auto h2d_w_tau = trackNode_test 
    
    .Define("w", [](ROOT::RVec<TvdcTrack*> tracks) 
	    { 
	      RVecD out; 
	      for (int t=0; t<tracks.size(); t++) { auto trk = tracks.at(t); 
		for (int p=0; p<4; p++) { 
		  auto group = trk->GetGroup(p); 
		  for (int h=0; h<group->Nhits(); h++) { 
		    
		    double w = group->WirePos(h) - trk->Intercept(p); 
		    
		    w = TMath::Abs(w*trk->Slope(p)); 
		    
		    out.push_back( w ); 
		  }
		}
	      }
	      return out; }, {branch_track})
    
    /*.Define("tau", [](ROOT::RVec<TvdcTrack*> tracks) 
	    { 
	      RVecD out; 
	      for (int t=0; t<tracks.size(); t++) { auto trk = tracks.at(t); 
		for (int p=0; p<4; p++) { 
		  auto group = trk->GetGroup(p); 
		  for (int h=0; h<group->Nhits(); h++) { 
		    out.push_back( group->Time(h) ); 
		  }
		}
	      }
	      return out; }, {branch_track})
    
    .Histo2D({"w_tau", ";drift-dist (m);drift-time (s);", 
    200, 0, 20e-3, 200, 0, -15e-9, 300e-9}, "w", "tau"); */ 
    
    .Histo1D({"w", ";drift-dist (m);", 200, 0, 20e-3}, "w"); 
    
    
  auto h1d_nGoodPts_plane = trackNode_test
    
    .Define("nGoodPts_plane", [](ROOT::RVec<TvdcTrack*> tracks) 
	    { 
	      RVecD out; 
	      for (int t=0; t<tracks.size(); t++) { auto trk = tracks.at(t); 
		for (int p=0; p<4; p++) 
		  out.push_back( trk->Get_nGoodPoints(p) ); 
		  
	      }
	      return out; }, {branch_track})
    
    .Histo1D({"ngp", "Number of good points per plane", 16, -0.5, 15.5}, 
	     "nGoodPts_plane"); 
    
	    

  auto h1d_tau_error = trackNode_test 
    
    .Define("tau_err", [](ROOT::RVec<TvdcTrack*> tracks) 
	    { 
	      RVecD out; 
	      for (int t=0; t<tracks.size(); t++) { auto trk = tracks.at(t); 
		for (int p=0; p<4; p++) { 
		  auto group = trk->GetGroup(p); 
		  for (int h=0; h<group->Nhits(); h++) { 
		    double tau_err 
		      = trk->Get_T_model(p,group->WirePos(h)) + trk->T0() 
		      - group->Time(h); 
		    
		    out.push_back( tau_err ); 
		  }
		}
	      }
	      return out; }, {branch_track})
    
    .Histo1D({"tau_err", "Tau-error (s)", 200, -100e-9, 100e-9}, "tau_err"); 
    
  TString canv_title = TString::Format("run %4i, %2.2e events,",
				       runNum, (double)maxEntries); 
  
  canv_title += isRightArm ? "Right arm" : "Left arm"; 
  
  
  auto canv_1 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);  
  
  /*h2d_slope_Dt->DrawCopy("col"); 

  auto graph = th2d_to_graph( (TH2D*)h2d_slope_Dt->Clone(), 1.16, 1.70 ); 
  
  graph->Draw("SAME"); */   
  const double dv_error = compute_minVariance( (TH1D*)h1d_dv->Clone() ); 
  
  h1d_dv->SetTitle(TString::Format("DV-error = %3.5f mm", 
				   TMath::Sqrt(dv_error)*1e3)); 
  
  h1d_dv->GetYaxis()->SetRangeUser(0, h1d_dv->GetMaximum()*1.20); 
  h1d_dv->DrawCopy("E");
    
  return; 
  
  //h2d_xParam_Dt->DrawCopy("col"); 
  
  //h2d_theta_s2x->DrawCopy(); 
  
  auto canv_2 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);  
  
  
  //draw x-param plot
  h2d_xParam_Dt->DrawCopy("col"); 
  
  draw_vLine( (TH2D*)h2d_xParam_Dt->Clone(), -TRK_CUT_xParam ); 
  draw_vLine( (TH2D*)h2d_xParam_Dt->Clone(),  TRK_CUT_xParam ); 
    
  draw_hLine( (TH2D*)h2d_xParam_Dt->Clone(), -TRK_CUT_Dt ); 
  draw_hLine( (TH2D*)h2d_xParam_Dt->Clone(),  TRK_CUT_Dt ); //*/ 
  
  //h2d_phi_s2y->DrawCopy(); 
  //
  //h2d_eta->DrawCopy("col"); //*/
  /*h2d_slope_Dt->DrawCopy("E"); 
  
    draw_gausFit_fixedBase( (TH1D*)h2d_slope_Dt->Clone(), 10e-9 ); */ 
  
  auto canv_3 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);  
  
  h1d_tau_error->DrawCopy("E"); 
  
  draw_vLine_th1( (TH1D*)h1d_tau_error->Clone(), -40e-9 ); 
  draw_vLine_th1( (TH1D*)h1d_tau_error->Clone(),  40e-9 ); 
  
  draw_gausFit_fixedBase( (TH1D*)h1d_tau_error->Clone(), 10e-9 );  

  
  
  auto canv_4 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);  
  
  h2d_w_tau->DrawCopy("E"); 
    
  
  
  auto canv_5 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);  
  
  h1d_nGoodPts_plane->DrawCopy(); 

#else 

  
  const double CUT_monteError    = 0.75e-3; 
  const double CUT_monte_T_error = 15e-9; 
  
  const double goodPts_tauCut = 40e-9; 
  
  
  auto monte_check = track_node 
    
    .Define("nGoodTracks", [CUT_monteError, 
			    CUT_monte_T_error]( ROOT::RVec<TvdcTrack*> tracks, 
						RVecD params_good,
						RVecD params_accident ) { 
	      int nGoodTracks(0); 
	      
	      for (int t=0; t<tracks.size(); t++) { 
		
		bool match_good     =true; 
		bool match_accident =true; 
		
		auto trk = tracks.at(t); trk->Set_isGoodTrack(true); 
		
		for (int p=0; p<4; p++) { 
		  
		  double Err; 
		  Err = TMath::Abs(params_good[p]-trk->Intercept(p)); 
		  if (Err > CUT_monteError) match_good    =false;  
		    
		  Err = TMath::Abs(params_accident[p]-trk->Intercept(p)); 
		  if (Err > CUT_monteError) match_accident=false;
		}
		
		double Err_T; 
		Err_T = TMath::Abs(trk->T0()-params_good[4]); 
		if (Err_T > CUT_monte_T_error) match_good     =false; 
		
		Err_T = TMath::Abs(trk->T0()-params_accident[4]); 
		if (Err_T > CUT_monte_T_error) match_accident =false; 
		
		if (match_good || match_accident) { trk->Set_isGoodTrack(true);  }
		else                              { trk->Set_isGoodTrack(false); }
		
		if (trk->IsGoodTrack()) nGoodTracks++; 
	      }
	      return nGoodTracks; 
	      
	    }, {"tracks_L_refined", "MONTE_track_params",
		"MONTE_track_params_accident"})
    
    .Define("nBadTracks", [](int nGoodTracks, 
			     ROOT::RVec<TvdcTrack*> v)
	    { return v.size() - nGoodTracks; }, 
	    {"nGoodTracks", "tracks_L_refined"})
    
    .Define("tracks_L_eta", [](ROOT::RVec<TvdcTrack*> tracks, 
			       int nGoodTracks) 
	    { RVecD eta; 
	      
	      for (int t=0; t<tracks.size(); t++) 
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad)
		  eta.push_back( tracks.at(t)->Get_Eta() ); 
	      
	      return eta; }, {"tracks_L_refined", "nGoodTracks"})
    
    
    //forcal plane intercepts
    .Define("tracks_s2_x", [](ROOT::RVec<TvdcTrack*> tracks, 
			      int nGoodTracks) { 
	      RVecD v; 
	      //I need the column 'nGoodTracks' as an (irritating) workaround. 
	      //Basically, I need the tracks to 'remember' whether I set
	      //  them as 'good' or 'bad.' 
	      
	      for (int t=0; t<tracks.size(); t++) 
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad)  
		  v.push_back(tracks.at(t)->S2_x()); 
	      
	      return v; }, {"tracks_L_refined", "nGoodTracks"})
    
    .Define("tracks_s2_y", [](ROOT::RVec<TvdcTrack*> tracks, 
			      int nGoodTracks) { 
	      RVecD v; 
	      
	      for (int t=0; t<tracks.size(); t++)	  
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad) 
		  v.push_back(tracks.at(t)->S2_y()); 
	      
	      return v; }, {"tracks_L_refined", "nGoodTracks"})
    
    .Define("tracks_theta", [](ROOT::RVec<TvdcTrack*> tracks, 
			       int nGoodTracks) { 
	      RVecD v; 
	      
	      for (int t=0; t<tracks.size(); t++)	  
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad) 
		  v.push_back(tracks.at(t)->Theta()); 
	      
	      return v; }, {"tracks_L_refined", "nGoodTracks"})
    
    .Define("tracks_phi", [](ROOT::RVec<TvdcTrack*> tracks,  
			     int nGoodTracks) { 
	      RVecD v; 
	      
	      for (int t=0; t<tracks.size(); t++)	  
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad) 
		  v.push_back(tracks.at(t)->Phi()); 
	      
	      return v; }, {"tracks_L_refined", "nGoodTracks"}) 
    
    .Define("tracks_xParam", [](ROOT::RVec<TvdcTrack*> tracks, 
				int dummy) { 
	      RVecD v; 
	      
	      for (int t=0; t<tracks.size(); t++)	  
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad) 
		  v.push_back(tracks.at(t)->xParam()); 
	      
	      return v; }, {"tracks_L_refined", "nGoodTracks"}) 
    
    .Define("tracks_Dt", [](RVecD track_param, 
			    ROOT::RVec<TvdcTrack*> tracks, 
			    int dummy) { 
	      RVecD v; 
	      
	      for (int t=0; t<tracks.size(); t++)	  
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad) 
		  v.push_back(tracks.at(t)->T0() - track_param[4]); 
	      
	      return v; }, 
	    {"MONTE_track_params", "tracks_L_refined", "nGoodTracks"}) 
    
    .Define("tracks_ERR_int_guess", [](ROOT::RVec<TvdcTrack*> tracks, 
					 int dummy) { 
	      RVecD v; 
	      
	      //cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; 
	      
	      for (int t=0; t<tracks.size(); t++) 
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad) 
		  for (int p=0; p<4; p++) { 
		    
		    /*cout << TString::Format("err-guess(%1i) = %2.4f", 
					    p,tracks.at(t)->Error_intercept(p)*1e3)
					    << endl; */ 
		    v.push_back(tracks.at(t)->Error_intercept(p)); 
		  }
	      
	      return v; }, {"tracks_L_refined", "nGoodTracks"}) 
    
    .Define("tracks_ERR_int_actual", [](ROOT::RVec<TvdcTrack*> tracks, 
					RVecD track_params, 
					int dummy) { 
	      RVecD v; 
	      
	      for (int t=0; t<tracks.size(); t++) 
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad) 
		  for (int p=0; p<4; p++)
		    v.push_back( tracks.at(t)->Intercept(p) - track_params[p] ); 
	      
	      return v; }, {"tracks_L_refined", "MONTE_track_params", "nGoodTracks"})
    
    .Define("tracks_ERR_Theta_guess",[](ROOT::RVec<TvdcTrack*> tracks, 
					int dummy) { 
	      
	      RVecD v; 
	      
	      for (int t=0; t<tracks.size(); t++) 
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad) 
		  v.push_back( tracks.at(t)->Error_Theta() ); 
	      
	      return v; }, {"tracks_L_refined", "nGoodTracks"}) 
    
    .Define("tracks_ERR_Theta_actual",[](ROOT::RVec<TvdcTrack*> tracks, 
					 RVecD track_params, 
					 int dummy) { 
	      
	      RVecD v; 
	      
	      for (int t=0; t<tracks.size(); t++) {
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad) { 
		  
		  double params[4]; 
		  for (int p=0; p<4; p++) params[p]=track_params[p]; 
		  
		  double Theta=0., Phi=0.; 
		  
		  TvdcTrack::Compute_Theta_Phi( kArm_Left, 
						params, 
						Theta, 
						Phi );  
		  
		  v.push_back( tracks.at(t)->Theta() - Theta ); 
		}
	      }
	      
	      return v; }, {"tracks_L_refined", "MONTE_track_params", "nGoodTracks"}); 		  
  
  
  TString canv_title = TString::Format("run %4i, %1.3e evts,",
				       runNum,(double)maxEntries); 
  
  canv_title += isRightArm ? "Right arm" : "Left arm"; 
  

  canv_title += " Monte-Carlo"; 
  canv_title += choose_goodOrBad ? " GOOD tracks" : " BAD tracks"; 
  
  auto Phi_model_draw   = [](double *S2_y, double *offset) { 
    
    return -0.231*TMath::Power(S2_y[0],2) 
    +  0.2532*S2_y[0]
    +  0.539e-3  + offset[0]; 
  };
  
  auto Theta_model_draw = [](double *S2_x, double *offset) { 
    
    return 0.109648*S2_x[0] + offset[0]; 
  };
    
  
  auto canv_1 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);  
  
  /*auto h1d_Theta_err = monte_check 
    
    .Histo1D({"Theta_err", "#Theta-error prediction;rad;", 150, 0.2e-3, 0.6e-3}, 
	     "tracks_ERR_Theta_guess"); 
  
  h1d_Theta_err->DrawCopy(); 
  
  return; //*/ 
  
  
  
  
  auto h2d_err = monte_check 
    
    .Histo2D({"err_guess", "Track #Theta-error;guess (rad);actual (rad)", 
	  9, 0.275e-3, 0.500e-3, 
	  600, -0.8e-3, 0.8e-3}, 
      "tracks_ERR_Theta_guess", "tracks_ERR_Theta_actual"); 
  
  const int nBinsX = h2d_err->GetXaxis()->GetNbins(); 
  
  double X[nBinsX];  h2d_err->GetXaxis()->GetCenter( X ); 
  double Y[nBinsX]; 
  double Y_err[nBinsX]; 
  
  double X_err[nBinsX]; 
  
  for (int bin=1; bin<=nBinsX; bin++) { 
    
    auto proj = h2d_err->ProjectionY("proj",bin,bin); 
    
    double sigma = StdDev_fixedMean( proj, 0. ); 
    
    Y    [bin-1] = sigma; 
    Y_err[bin-1] = sigma/TMath::Sqrt(proj->Integral()); 
    
    X_err[bin-1] = 0.; 
  }
  
  auto errGraph = new TGraphErrors( nBinsX, 
				    X, Y, 
				    X_err, Y_err ); 
  
  auto fit_ = new TF1("line_fit", "x", 0.2e-3, 1.0e-3);
  
  
  errGraph->SetMarkerStyle(kOpenCircle); 
  
  errGraph->GetXaxis()->SetRangeUser(0.15e-3,0.60e-3); 
  errGraph->GetYaxis()->SetRangeUser(0.15e-3,0.60e-3); 
  
  errGraph->Draw("A L"); 
  
  
  fit_->SetLineWidth(1); 
  fit_->SetLineStyle(kDashed); 
    
  fit_->Draw("SAME"); 
  
    
  
  //h1d_err_guess->DrawCopy(); 
    
#if 0 
  
  auto h2d_eta_ngp = nEvents_1refineTrack 
    
    .Define("tracks_L_nGoodPts", [](ROOT::RVec<TvdcTrack*> tracks, 
				    RVecD dummy) { 
	      ROOT::RVec<int> ngp; 
	      
	      for (int t=0; t<tracks.size(); t++) 
		if (tracks.at(t)->IsGoodTrack()==choose_goodOrBad) 
		  ngp.push_back( tracks.at(t)->Get_nGoodPts() ); 
	      
	      return ngp; }, {"tracks_L_refined","tau_error"}) //*/ 
    
    .Histo2D({"eta_ngp", "N.good pts vs Eta;nGoodPts;Eta", 
	  26, -0.5, 25.5, 
	  200, 0, 25}, "tracks_L_nGoodPts", "tracks_L_eta"); //*/ 
  
  /*.Histo1D({"", "track tau-RMS (cutoff = 40ns);(s);",  100, 0, 15e-9}, 
    "tau_error"); //*/  
      
  h2d_eta_ngp->DrawCopy("col"); 
  
  return; 
  
    
  gStyle->SetOptTitle(0); 
  
  auto h1d_nGood = monte_check 
    .Histo1D({"nGood", "N. good", 6, -0.5, 5.5}, "nGoodTracks"); 
  
  auto h1d_nBad = monte_check 
    .Histo1D({"nBad", "N. bad",   6, -0.5, 5.5}, "nBadTracks"); 
  
  h1d_nGood->Scale( 1./h1d_nGood->Integral() ); 
  h1d_nGood->DrawCopy("HIST"); 
  
  h1d_nBad->Scale( 1./h1d_nBad->Integral() ); 
  h1d_nBad->SetLineColor(kRed); h1d_nBad->SetLineStyle(kDashed); 
  h1d_nBad->DrawCopy("HIST SAME"); 
  
  gPad->BuildLegend(); //*/ 
    
  
  return; //*/ 
  
  /*auto h1d_etaGood = monte_check 
    
    .Define("eta_good", [](ROOT::RVec<TvdcTrack*> tracks, 
			   int dummy)
	    { RVecD eta; 
	      
	      for (int t=0; t<tracks.size(); t++) 
		if (tracks.at(t)->IsGoodTrack()==true) 
		  eta.push_back( tracks.at(t)->Get_Eta() ); 
	      
	      return eta; }, {"tracks_L_refined", "nGoodTracks"})     
    
    .Histo1D({"etagood", "good tracks", 200, 0, 25}, "eta_good"); 
  
  h1d_etaGood->DrawCopy(); 
  
  
  
  auto h1d_etaBad = monte_check 
    
    .Define("eta_bad", [](ROOT::RVec<TvdcTrack*> tracks, 
			  int dummy)
	    { RVecD eta; 
	      
	      for (int t=0; t<tracks.size(); t++) 
		if (tracks.at(t)->IsGoodTrack()==false)
		  eta.push_back( tracks.at(t)->Get_Eta() ); 
	      
	      return eta; }, {"tracks_L_refined", "nGoodTracks"})     
    
    .Histo1D({"etabad", "bad tracks", 200, 0, 25}, "eta_bad"); 
  
  h1d_etaBad->SetLineColor(kRed); 
  h1d_etaBad->DrawCopy("SAME"); 
    
  gPad->BuildLegend(); 
  
  
  
  /*auto h2d_xParam_Dt = monte_check 
    
    .Histo2D({"xP_Dt", ";x-param;T_{track}-T_{S2}", 
	  200, -4.5, 4.5, 200, -100e-9, 100e-9}, "tracks_xParam", "tracks_Dt"); 
	  
	  
	  h2d_xParam_Dt->DrawCopy("col"); //*/ 
  
  
 
  
  /*auto h1d_tauErr = monte_check 
    .Histo1D({"tauErr", "Tau error (s)", 200, -100e-9, 100e-9}, 
	     "tau_error"); 
  
	     h1d_tauErr->DrawCopy(); return; //*/ 
    
  
  
  /*auto h2d_Theta_s2x = monte_check
    .Histo2D({"theta_s2x", "#Theta vs X_{S2};X_{S2} (m);#Theta (rad)", 
	  200, -1.2, 1.2, 200, -0.15, 0.15}, "tracks_s2_x", "tracks_theta"); 
  
  h2d_Theta_s2x->DrawCopy("col"); 
    
  auto fThetaCut = new TF1("thetaCut", Theta_model_draw, 
			   -1.2, 1.2, 1);
  
  fThetaCut->SetLineStyle(kDashed); 
  
  fThetaCut->SetParameter(0,CUT_th_min); 
  fThetaCut->DrawCopy("SAME"); 
  
  fThetaCut->SetParameter(0,CUT_th_max); 
  fThetaCut->DrawCopy("SAME"); //*/ 

  
  
  /*auto canv_2 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);
  
  auto h2d_Phi_s2y = monte_check
    .Histo2D({"phi_s2y", "#Phi vs Y_{S2};Y_{S2} (m);#Phi (rad)", 
	  200, -0.325, 0.325, 200, -0.15, 0.15}, "tracks_s2_y", "tracks_phi"); 
  
  h2d_Phi_s2y->DrawCopy("col"); 
  
  auto fPhiCut = new TF1("phiCut", Phi_model_draw,
			 -0.32, 0.32, 1);
  
  fPhiCut->SetLineStyle(kDashed); 
  
  fPhiCut->SetParameter(0,CUT_ph_min); 
  fPhiCut->DrawCopy("SAME"); 
  
  fPhiCut->SetParameter(0,CUT_ph_max); 
  fPhiCut->DrawCopy("SAME");     //*/ 
  
  
  auto canv_2 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);   
  
  auto h2d_goodBad = monte_check 
    .Histo2D({"nGood_nBad", ";N. good;N. bad", 6, -0.5, 5.5, 6, -0.5, 5.5}, 
	     "nGoodTracks","nBadTracks"); 
  
  h2d_goodBad->Scale( 1./h2d_goodBad->Integral() ); 
  h2d_goodBad->DrawCopy("colz"); //*/
    
#endif 


#endif 
}
