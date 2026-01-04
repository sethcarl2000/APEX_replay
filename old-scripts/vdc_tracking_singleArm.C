#include <TROOT.h>
#include "def_apex.h"
#include "function_lib.C"
#include <fstream>

//#include "class/TvdcHit.h"
//#include "class/TvdcHit.C"
#include "apex_lib/ApexHRS.h"
#include <functional>
#include <string>
#include <cmath>

using namespace std; 
using namespace APEX; 

using RVecD = ROOT::VecOps::RVec<double>;
using RNode = ROOT::RDF::RNode; 
using uint  = unsigned int; 

//this function adds branches to the output branch list
  
//path which this code is intended to be executed from
const string PATH_EXECUTE = "/work/halla/apex/disk1/sethhall/analyzer/apex-install/analyzer-APEX/scripts/"; 

//                           //
#define IS_MONTECARLO false  //
#define RUN_LIMITED true     //
#define DEBUG false          //
#define MAKE_OUTPUT true     //
//                           //
//these are used for track error-estimation
double dP; 
double cum_P; 

const int err_nSamples = 25; 

double z[err_nSamples]; 

bool choose_goodOrBad=true;  //true==look at 'good' tracks; false==bad tracks 

//#endif 

//track cuts
const double TRK_CUT_Dt     = 40e-9; 
const double TRK_CUT_xParam = 1.55; 
const double TRK_CUT_Eta    = 3.750; 
const double TRK_measureSigma = 5e-9; 
const int    TRK_CUT_nGoodPts_min_perPlane = 2; //good points per plane
const int    TRK_CUT_nGoodPts_min          = 12;   

const double CUT_minEta     = 1.950; //min eta of one plane (during grid-searching). 
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
    cout << "DEBUG: Vdc valid real-time range: ["
	 << VDC_min_realTime << ", " << VDC_max_realTime << "]" << endl;
    int n_validHits(0);
    cout << "Is right arm? " << (evt->ActiveArm() ? "true" : "false") << "\n";
    cout << " plane " << p << endl;
    cout << " S2 hit time: " << evt->GetS2Hit()->Time()*1e9 << endl; 
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
      
#if DEBUG 
      cout << TString::Format("wire=%3i, rawtime=%5.f, time=%5.1f  time(check)=%5.1f",
			      (int)h_wire[h], 
			      h_time[h], 
			      hit->Time()*1e9,
			      THRS::VDC::RealTime(evt->ActiveArm(),
						  p,
						  h_wire[h],
						  h_time[h])*1e9 ); 
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
    const double measure_sigma   =5e-9; //sigma*2.; 
    
    ///////////////////////////////////////////////////////////////////////////////
    double int0[5]; 
    for (int p=0; p<4; p++) int0[p] = trk->Intercept(p); 
    TString plane_name[4] = {"U1","V1","U2","V2"}; 
    int0[4] = trk->T0(); 
    //*/////////////////////////////////////////////////////////////////////////////
    
    double final_eta[4] = {0.}; 
    
    const double s2 = sigma*sigma; 
    
    ROOT::RVec<double> tVec; 
    
    for (uint c=0; c<nCycles; c++) { 
      
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
	
	for (uint h=0; h<group->Nhits(); h++) { 
	  	  
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
	  
	}//for (uint h=0; h<trk->pGet_N(p); 
	
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


//used to transfer target-coordinate information to process optics
struct TTargetCoords {
  
  TTargetCoords() {}; 
  TTargetCoords( const double y,
		 const double dxdz, 
		 const double dydz, 
		 const double dp ) 
    : tg_y(y), tg_dxdz(dxdz), tg_dydz(dydz), tg_dp(dp) {}; 
  
  ~TTargetCoords(); 
  
  double tg_y; 
  double tg_dxdz; 
  double tg_dydz; 
  double tg_dp; 
}; 


RNode generate_optics_output( bool is_RightArm, RNode track_node, 
			      vector<string> &out_branches, 
			      THRS *hrs ) { 
  //compute & output optics data
  string branch_tr 
    = is_RightArm ? "track_R_refined" : "tracks_L_refined"; 
  
  
  TString arm = is_RightArm ? "R" : "L"; 
  
  
  return track_node; 
}

struct TTrackCoords { 
  
  //this used to pass track coords back and forth between nodes in the function 
  // below
  TTrackCoords(const double ix, 
	       const double iy, 
	       const double idxdz, 
	       const double idydz) 
    : x(ix), y(iy), dxdz(idxdz), dydz(idydz) {}; 
  
  ~TTrackCoords() {}; 
  
  double x, y, dxdz, dydz; 
}; 

RNode generate_track_output( bool is_RightArm, RNode track_node, 
			     vector<string> &out_branches, 
			     THRS *hrs ) { 
  
  string branch_tr = is_RightArm ? "tracks_R_refined" : "tracks_L_refined"; 
  
  TString arm = is_RightArm ? "R" : "L"; 
  
  //this defines 4 positional variables, which can be defined for transport & 
  // detector coordinate systems.
  auto Make_coordBranches = [arm, &out_branches](const TString coordName, 
						 RNode node) { 
    
    string br_coord = (string)TString(arm+"_tr_coords_"+coordName); 
    
    out_branches.push_back( (string)TString(arm+"_tr_"+coordName+"_x") );
    
    node = node  
    .Define( out_branches.at(out_branches.size()-1),  
	     [](ROOT::RVec<TTrackCoords*> coords) { 
	       RVecD v; 
	       for (TTrackCoords *c: coords) v.push_back( c->x );
	       return v; }, {br_coord}); 
    
    out_branches.push_back( (string)TString(arm+"_tr_"+coordName+"_y") );
    
    node = node
    .Define( out_branches.at(out_branches.size()-1),  
	     [](ROOT::RVec<TTrackCoords*> coords) { 
	       RVecD v; 
	       for (TTrackCoords *c: coords) v.push_back( c->y );
	       return v; }, {br_coord}); 
    
    out_branches.push_back( (string)TString(arm+"_tr_"+coordName+"_th") );
    
    node = node
    .Define( out_branches.at(out_branches.size()-1),  
	     [](ROOT::RVec<TTrackCoords*> coords) { 
	       RVecD v; 
	       for (TTrackCoords *c: coords) v.push_back( c->dxdz );
	       return v; }, {br_coord}); 
    
    out_branches.push_back( (string)TString(arm+"_tr_"+coordName+"_ph") );
    
    node = node
    .Define( out_branches.at(out_branches.size()-1),  
	     [](ROOT::RVec<TTrackCoords*> coords) { 
	       RVecD v; 
	       for (TTrackCoords *c: coords) v.push_back( c->dydz );
	       return v; }, {br_coord}); 
    
    return node; 
  }; 
        
  string br_transport = (string)arm+"_tr_coords_tra"; 
  
  //define the 'transport' coordinates (z-axis normal to S2-plane)
  track_node = track_node
    .Define( (string)arm+"_tr_coords_tra", 
	     [](ROOT::RVec<TvdcTrack*> tr) 
	     { ROOT::RVec<TTrackCoords*> coords; 
	       for (int t=0; t<tr.size(); t++) { 
		 auto trk = tr.at(t); 
		 
		 coords.push_back( new TTrackCoords( trk->FP_x(), 
						     -trk->FP_y(), 
						     TMath::Tan(trk->Theta()), 
						     -TMath::Tan(trk->Phi()) ) ); 
	       } return coords; }, {branch_tr}); 
  
  //define the positional parameters for this position
  track_node = Make_coordBranches( "tra", track_node ); 
  
  //define the 'transport' coordinates (z-axis normal to S2-plane)
  track_node = track_node
    .Define( (string)arm+"_tr_coords_det", 
	     [](ROOT::RVec<TvdcTrack*> tr) 
	     { ROOT::RVec<TTrackCoords*> coords; 
	       for (int t=0; t<tr.size(); t++) { 
		 
		 auto trk = tr.at(t); 
		 
		 TVector3 S( TMath::Tan(trk->Theta()), 
			     TMath::Tan(trk->Phi()), 
			     1. ); 
		 
		 S.RotateY( TMath::Pi()/4. ); 
		 
		 //now, compute x & y 
		 TVector3 r0( trk->Intercept(1) - 0.026/trk->Slope(1), 
			      trk->Intercept(0), 
			      0. ); 
		 
		 r0.RotateZ( -TMath::Pi()/4. ); 
		 
		 coords.push_back( new TTrackCoords( r0.X(), 
						     -r0.Y(), 
						     S.X()/S.Z(), 
						     -S.Y()/S.Z() ) ); 
	       } return coords; }, {branch_tr}); 
  
    
  track_node = Make_coordBranches( "det", track_node ); 
    

  
  //now, we can compute optics information
  string branch_tgCoords 
    = is_RightArm ? "track_R_tgCoords" : "tracks_L_tgCoords"; 
    
  //first, we define a TTargetCoords object which will store data for us. 
  //this obj will store the results of the optics computations. 
  track_node = track_node

    
    //now, compute target coor
    .Define(branch_tgCoords, 
	    [hrs, is_RightArm](ROOT::RVec<TTrackCoords*> coords_det, 
			       ROOT::RVec<TTrackCoords*> coords_tra) { 
	      
	      ROOT::RVec<TTargetCoords*> tgCoords; 
	      
	      for (uint t=0; t<coords_det.size(); t++) { 
	      	
		double TG_y, TG_dxdz, TG_dydz, TG_dp; 

		auto c_det = coords_det.at(t); 
		auto c_tra = coords_tra.at(t); 
				
		hrs->Compute_trackOptics( is_RightArm, 
					  c_tra->x, 
					  c_tra->y,
					  c_det->dxdz, 
					  c_det->dydz, 
					  TG_y,
					  TG_dxdz, 
					  TG_dydz, 
					  TG_dp ); 
		
		tgCoords.push_back( new TTargetCoords( TG_y, 
						       TG_dxdz,
						       TG_dydz,
						       TG_dp ) ); 
		
		/*cout << TString::Format("y_tg new=(%3.3e) old=(%3.3e) dif=(%3.3e)", 
					TG_y, 
					tg_y.at(t), 
					TG_y-tg_y.at(t)) << endl; */ 
		
	      } 
	      return tgCoords; }, 
	    { TString(arm+"_tr_coords_det").Data(), 
	      TString(arm+"_tr_coords_tra").Data() }); 
  
  //write the target coordinates as output branches
  out_branches.push_back( (string)TString(arm+"_tr_tg_y") );
  
  track_node = 
    track_node.Define(out_branches.at(out_branches.size()-1), 
		      [](ROOT::RVec<TTargetCoords*> tgCoords) {
			RVecD v; 
			for (TTargetCoords *t: tgCoords) v.push_back( t->tg_y ); 
			return v; }, {branch_tgCoords}); 
  
  out_branches.push_back( (string)TString(arm+"_tr_tg_th") );
  
  track_node = 
    track_node.Define(out_branches.at(out_branches.size()-1), 
		      [](ROOT::RVec<TTargetCoords*> tgCoords) {
			RVecD v; 
			for (TTargetCoords *t: tgCoords) v.push_back( t->tg_dxdz ); 
			return v; }, {branch_tgCoords}); 
  
  
  out_branches.push_back( (string)TString(arm+"_tr_tg_ph") );
  
  track_node = 
    track_node.Define(out_branches.at(out_branches.size()-1), 
		      [](ROOT::RVec<TTargetCoords*> tgCoords) {
			RVecD v; 
			for (TTargetCoords *t: tgCoords) v.push_back( t->tg_dydz ); 
			return v; }, {branch_tgCoords}); 
  
  
  out_branches.push_back( (string)TString(arm+"_tr_dp") );
  
  track_node = 
    track_node.Define(out_branches.at(out_branches.size()-1), 
		      [](ROOT::RVec<TTargetCoords*> tgCoords) {
			RVecD v; 
			for (TTargetCoords *t: tgCoords) v.push_back( t->tg_dp ); 
			return v; }, {branch_tgCoords});   
  
  //now write plane-wise parameters
  TString plane_name[4] = { "U1","V1","U2","V2" }; 
  
  const double sigma_eta = 9e-9; 
  const double sigma_cutoff = 40e-9; 
  
  for (int p=0; p<4; p++) { 
    
    out_branches.push_back((string)TString(arm+"_tr_"+plane_name[p]+"_intercept"));
    out_branches.push_back((string)TString(arm+"_tr_"+plane_name[p]+"_intercept_error"));
    out_branches.push_back((string)TString(arm+"_tr_"+plane_name[p]+"_nGoodPoints"));
    out_branches.push_back((string)TString(arm+"_tr_"+plane_name[p]+"_Eta"));
    
    track_node = track_node 
      
      .Define( (string)TString(arm+"_tr_"+plane_name[p]+"_intercept"), 
	      [p](ROOT::RVec<TvdcTrack*> tr) 
	      { RVecD v; 
		for (TvdcTrack *trk: tr) 
		  v.push_back( trk->Intercept(p) ); 
		return v; }, {branch_tr})
      
      .Define( (string)TString(arm+"_tr_"+plane_name[p]+"_intercept_error"), 
	      [p](ROOT::RVec<TvdcTrack*> tr) 
	      { RVecD v; 
		for (TvdcTrack *trk: tr) 
		  v.push_back( trk->Error_intercept(p) ); 
		return v; }, {branch_tr})
      
      .Define( (string)TString(arm+"_tr_"+plane_name[p]+"_nGoodPoints"), 
	      [p](ROOT::RVec<TvdcTrack*> tr) 
	      { ROOT::RVec<int> v; 
		for (TvdcTrack *trk: tr) 
		  v.push_back( trk->Get_nGoodPoints(p) ); 
		return v; }, {branch_tr})
      
      .Define( (string)TString(arm+"_tr_"+plane_name[p]+"_Eta"), 
	      [p, sigma_eta](ROOT::RVec<TvdcTrack*> tr) 
	      { RVecD v; 
		for (TvdcTrack *trk: tr) 
		  v.push_back( trk->Get_Eta(p) ); 
		return v; }, {branch_tr}); 
    
  }
  
  out_branches.push_back( (string)(arm+"_tr_Eta") );
  
  track_node = track_node 
    .Define( out_branches.at(out_branches.size()-1),  
	     [](ROOT::RVec<TvdcTrack*> tr) 
	     { RVecD v; 
	       for (TvdcTrack *trk: tr) 
		 v.push_back( trk->Get_Eta() ); 
	       return v; }, {branch_tr}); 
  
  
  out_branches.push_back( (string)(arm+"_tr_RMS") );
  
  track_node = track_node 
    .Define( out_branches.at(out_branches.size()-1),  
	     [](ROOT::RVec<TvdcTrack*> tr) 
	     { RVecD v; 
	       for (TvdcTrack *trk: tr) 
		 v.push_back( trk->Get_RMS() ); 
	       return v; }, {branch_tr}); 
  
  
  out_branches.push_back( (string)(arm+"_tr_nGoodPoints") );
  
  track_node = track_node 
    .Define( out_branches.at(out_branches.size()-1),  
	     [](ROOT::RVec<TvdcTrack*> tr) 
	     { RVecD v; 
	       for (TvdcTrack *trk: tr) 
		 v.push_back( trk->Get_nGoodPoints() ); 
	       return v; }, {branch_tr}); 

  out_branches.push_back( (string)(arm+"_tr_dt") );
  
  track_node = track_node 
    .Define( out_branches.at(out_branches.size()-1),  
	     [](ROOT::RVec<TvdcTrack*> tr) 
	     { RVecD v; 
	       for (TvdcTrack *trk: tr) 
		 v.push_back( trk->T0() ); 
	       return v; }, {branch_tr}); 

  
  //find which cerenkov-cell a track is connected to
  auto find_cerenkovCell = [is_RightArm](TTrackCoords *trk) {

    const double z_cer = is_RightArm ?  1.800  :  1.810;
    const double x0    = is_RightArm ? -0.9495 : -0.9188;

    const double cell_width = 0.3829;

    TVector3 r_cer( trk->x + trk->dxdz*z_cer,
		    trk->y + trk->dydz*z_cer,
		    z_cer ); 
    
    //find x-number
    int cell = (int)TMath::Nint( (r_cer.X() - x0)*(2./cell_width) ); 

    if (r_cer.Y() > -1.612e-2) cell += 1;
    
    if (cell<0) cell=0;
    if (cell>9) cell=8;

    return cell;
  }; 
  
  //do a bit of PiD calculations
  auto get_cerenkov_adc = [&find_cerenkovCell,
			   is_RightArm](TvdcTrack *trk,
					RVecD cer_adc) {
    
    int cell = find_cerenkovCell( new TTrackCoords( trk->FP_x(), 
						    trk->FP_y(), 
						    TMath::Tan(trk->Theta()), 
						    TMath::Tan(trk->Phi()) ) ); 
    double adc = cer_adc.at(cell) + (is_RightArm ? 7.838e3:0); 
    
    //add the other-pmt
    if (cell%2==0) { adc += cer_adc.at(cell+1); }
    else           { adc += cer_adc.at(cell-1); }
    
    return adc; 
  }; 

  
  out_branches.push_back( (string)(arm+"_tr_n") );

  //This is a really dumb way to measure the size of an RVec, but
  // when i try to do it the more sensible way
  // (as in the '#if 0/#endif' block below), I get a segfault. go figure. 
  track_node = track_node
    .Define( out_branches.at(out_branches.size()-1),
	     [&get_cerenkov_adc](ROOT::RVec<TvdcTrack*> tr) 
	     { int size=0; 
	       for (TvdcTrack *trk: tr) size++;  
	       return size; }, {branch_tr});   
#if 0
  track_node = track_node
    .Define( out_branches.at(out_branches.size()-1),
	     [](ROOT::RVec<TvdcTrack*> tr) 
	     { return (int)tr.size(); }, {branch_tr});
#endif


  
  out_branches.push_back( (string)(arm+"_tr_cer_a_c") );
  
  track_node = track_node
    .Define( out_branches.at(out_branches.size()-1),
	     [&get_cerenkov_adc](ROOT::RVec<TvdcTrack*> tr,
				 RVecD cer_adc) 
	     { RVecD v; 
	       for (TvdcTrack *trk: tr) 
		 v.push_back( get_cerenkov_adc(trk,cer_adc) ); 
	       return v; }, {branch_tr, TString(arm+".cer.a_c").Data()}); 
  
  return track_node; 
}

bool isRightArm; 

const TString PATH_Execute = "/work/halla/apex/disk1/sethhall/analyzer/scripts";

const char* PATH_VOLATILE_DIR = {"/volatile/halla/apex/full_replay"}; 

template <typename F>
void add_output_branch( RNode &output_node,
			vector<string> &output_branches,
			const string output,
			const F* func,
			const vector<string> inputs ) {
  
  output_branches.push_back( output );
  
  output_node.Define( output, (*func), inputs );  
  return; 
} 


void vdc_tracking_singleArm( int maxEntries=3, 
			     bool IN_isRightArm=true, //true=RHRS, false=LHRS 
			     int runNum=3782,
			     TString path_inFile
			     ="$VOLATILE_DIR/optics/decode.3782.root",
			     TString path_outFile="replay/test/replay.3782.root",
			     TString optics_target="")
{

  
  
  const TString here = "vdc_tracking_singleArm()"; 
  
  //find the central-momentum of both arms using EPICS vars
  ROOT::EnableImplicitMT(); 
  ROOT::RDataFrame d_E("E", path_inFile.Data() );

  //check if epics tree 'E' is empty
  if ((ULong64_t)*d_E.Count() < 1) {
    Error(here,"Epics tree 'E' of path \"%s\" is empty!",path_inFile.Data());
    return; 
  }
    
  const double momentum_RHRS
    = d_E.Histo1D({"","",200,-1,-1},"HacR_D1_P0rb")->GetMean();
  
  const double momentum_LHRS
    = d_E.Histo1D({"","",200,-1,-1},"HacL_D1_P0rb")->GetMean();


  //check to see if our selected arm has zero central-momentum
  if ( (IN_isRightArm  && momentum_RHRS<0.05) ||
       (!IN_isRightArm && momentum_LHRS<0.05) ) {

    Warning(here,
	    "%s chosen, but central momentum for this run is %0.3f! (<0.05).",
	    IN_isRightArm ? "RHRS" : "LHRS",
	    IN_isRightArm ? momentum_RHRS : momentum_LHRS);
  }
  
  ROOT::DisableImplicitMT(); 

  auto hrs = new THRS; 
  
  //parse optics offsets
  const TString path_optics_R = PATH_Execute + "/DB/R.opticsMatrix.dat"; 
  const TString path_optics_L = PATH_Execute + "/DB/L.opticsMatrix.dat"; 
  
  hrs->Parse_opticsData(THRS::Arm_Right, path_optics_R); 
  hrs->Parse_opticsData(THRS::Arm_Left,  path_optics_L); 
  
  //check the initialization of the optics
  if (!hrs->PrintStatus_optics()) { 
    
    Error(here, "Optics not successfully initialized!");
    return; 
  }
  
  auto vdc = new THRS::VDC; 
  
  //parse VDC offsets
  const TString path_vdc_offsets_R = PATH_Execute + "/DB/R.vdc.offsets.dat"; 
  const TString path_vdc_offsets_L = PATH_Execute + "/DB/L.vdc.offsets.dat"; 
  
  vdc->Parse_offsets(path_vdc_offsets_R); 
  vdc->Parse_offsets(path_vdc_offsets_L); 
  
  
  //check initialization status of VDC-offsets
  if (!vdc->Print_status()) { 
    
    Error(here, "VDC not successfully initialized!"); 
    return; 
  }
  
  /*
  cout << "realtime = " << THRS::VDC::RealTime( THRS::Arm_Right, 0, 203, 2469.0 )*1e9 << endl;
  return; 
  
  cout << "Offsets: L:\n";
  for (uint p=0; p<4; p++) {
    for (uint w=0; w<10; w++) {
      cout << THRS::VDC::RealTime( THRS::Arm_Left, p, w, 2147 ) << endl;
    }
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"; 
  }
  cout << endl;

  cout << "Offsets: R:\n";
  for (uint p=0; p<4; p++) {
    for (uint w=0; w<10; w++) {
      cout << THRS::VDC::RealTime( THRS::Arm_Right, p, w, 2147 ) << endl;
    }
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"; 
  }
  cout << endl;
  return;
  */
  
  //create the 'TReactVertex' object, which computes the reaction vertex,
  // using raster information, and known v-wire positions. 
  TReactVertex reactVertex( IN_isRightArm,
			    path_inFile );
  
  isRightArm = IN_isRightArm; 
  
  vector<string> branch_S2; 
  
  TString STRING_replay
    ="   ";
  
  cout << "vdc_tracking_singleArm.C -- executing replay script.~~~~~~~~~~~~~~~~~~~~~~~" 
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
  
  if (IN_isRightArm) { 
    branch_S2.push_back( is_monteCarlo ? "R_s2_rt" : "R.s2.rt" );
    branch_S2.push_back( is_monteCarlo ? "R_s2_lt" : "R.s2.lt" ); 
  } else          {
    branch_S2.push_back( is_monteCarlo ? "L_s2_rt" : "L.s2.rt" );
    branch_S2.push_back( is_monteCarlo ? "L_s2_lt" : "L.s2.lt" ); 
  }
    
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
  cout << "single-threading mode, "; 
#else 
  ROOT::EnableImplicitMT(); 
  cout << "multi-threading mode, ";
#endif
  
  
  ROOT::RDataFrame d_T( "T", path_inFile.Data() ); 

  auto nTotalEvents = d_T.Count(); 

  
  cout << TString::Format("running %i events...",
			  RUN_LIMITED
			  ? (int)TMath::Min( (int)*nTotalEvents, maxEntries )
			  : (int)*nTotalEvents ) << endl; 
  
  
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
    .Define("S2_hits", [&Generate_S2_hits, IN_isRightArm](RVecD rt, RVecD lt)
    {
#if DEBUG
      cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

      cout << "S2 hits: (arm=" << (IN_isRightArm ? "R":"L") << ") \n";

      cout << "paddle raw times:\n";
      for (int p=0; p<16; p++)
	cout << TString::Format("%2i,   %5.f,%5.f\n", p, rt[p],lt[p]); 
            
      auto S2hits = Generate_S2_hits(IN_isRightArm, rt,lt);
      for (const TS2Hit* hit : S2hits)
	cout << TString::Format(" pad,time = %2i, % 3.1f",
				hit->Paddle(),
				hit->Time()*1e9) << endl; 
#endif      
      return Generate_S2_hits(IN_isRightArm,  rt,lt);
    }, 
      { branch_S2[0].data(), branch_S2[1].data() }) 
    
    //Find all possible both-arm S2 coincidences
    .Define("eventVec", [IN_isRightArm](ROOT::RVec<TS2Hit*> S2_hits, 
					double              beamCurrent, 
					uint                runNumber ) 
	    { ROOT::RVec<TEventHandler*> eVec; 
	      
	      for (uint i=0; i<S2_hits.size(); i++) 
		eVec.push_back( new TEventHandler( IN_isRightArm, 
						   beamCurrent, 
						   runNumber, 
						   S2_hits.at(i), 
						   S2_hits.at(i) ) ); 
	      return eVec; }, 
	    {"S2_hits", "hac_bcm_average", "fEvtHdr.fRun"})
    
    //Make a cut on Coincs, we only allow exactly 1-per-ev    
    .Filter([](ROOT::RVec<TEventHandler*> v)
	    { return v.size()>0;}, {"eventVec"})
    
    //get the event handler (so we don't have to keep passing it as a vector)
    .Define("event", [](ROOT::RVec<TEventHandler*> v) 
	    { return v.at(0); }, {"eventVec"}); 
  
  
  //these 'RResultPtr's will tell us about 
  ROOT::RDF::RResultPtr<unsigned long long> nPass_1group; 
  ROOT::RDF::RResultPtr<unsigned long long> nPass_1pair; 
  ROOT::RDF::RResultPtr<unsigned long long> nPass_1rawTrack;  
  
  bool node_arm = isRightArm ? kArm_Right : kArm_Left; 
  
  auto trackNode = generate_vdcTracks( isRightArm, 
				       nEvents_coinc, 
				       nPass_1group, 
				       nPass_1pair, 
				       nPass_1rawTrack ); 
  

  //this is to translate from target-coordinates to hall-coordinates
  const TVector3 D0_R(-1.101e-3, -3.885e-3, 0./*0.7946*/);
  const TVector3 D0_L(-1.301e-3,  6.672e-3, 0./*0.7958*/); 

  enum ECoordTranslate { kTarget_to_Hall=0, kHall_to_Target }; 

  auto Rotate_hcs_tcs = [](TVector3 &R,
			   const bool arm,
			   const ECoordTranslate type ) {
    
    const double theta0 = (arm==true? -5.372 : 5.366)*TMath::Pi()/180.;  
    
    if (type==kHall_to_Target) { //Hall-Coords to Target-coords
      R.RotateY(-theta0);
      R.RotateZ(TMath::Pi()/2.);
    } else                     { //Target-coords to hall-coords
      R.RotateZ(-TMath::Pi()/2.); 
      R.RotateY(theta0);
    }
  };

  
  auto Translate_coordinates = [D0_R,D0_L,
				&Rotate_hcs_tcs](const TVector3 RR,
						 const bool arm,
						 const ECoordTranslate type ) { 
    TVector3 R( RR ); 
    
    if (type==kHall_to_Target) {
      
      Rotate_hcs_tcs( R, arm, kHall_to_Target );    
      R +=  -(arm ? D0_R : D0_L); 
      
    } else                     {
    
      R +=   (arm ? D0_R : D0_L); 
      Rotate_hcs_tcs( R, arm, kTarget_to_Hall );
    }

    return R; 
  };

  
  TString arm = isRightArm ? "R" : "L"; 
  
  vector<string> out_branches; 
  
  out_branches.push_back( (string)TString(arm+"_cer_a_c") ); //.......Left HRS info
  out_branches.push_back( (string)TString(arm+"_cer_t_c") ); 
  out_branches.push_back( isRightArm?"R_ps_a_c":"L_prl1_a_c" );  
  out_branches.push_back( isRightArm?"R_sh_a_c":"L_prl2_a_c" ); 
  out_branches.push_back( "S2_paddle_realTime" );
  out_branches.push_back( "S2_paddle_num" );
  out_branches.push_back( "S2_isTwinHit" ); 
  out_branches.push_back( "S2_pmtDiff_raw" ); 
  out_branches.push_back( "Raster2_current_x" );
  out_branches.push_back( "Raster2_current_y" );
  out_branches.push_back( "r_BPMA" );
  out_branches.push_back( "r_BPMB" ); 
  out_branches.push_back( "beam_current" );//..........beam current
  out_branches.push_back( "RHRS_p" ); 
  out_branches.push_back( "LHRS_p" ); 
  out_branches.push_back( "beam_MeV" );   //...........beam energy
  out_branches.push_back( "event_number" );//..........event 'number' 
  out_branches.push_back( "react_vertex" ); 
  out_branches.push_back( "react_vertex_TCS" );
  //out_branches.push_back( "avg_react_vertex" ); 
  
  TString rb_name  = isRightArm ? "Rrb" : "Lrb"; 
  
  //define some output branches for the arm we've chosen
  auto output_node = trackNode
    
    .Define("S2_paddle_realTime", [](TEventHandler *event) 
    { return event->GetS2Hit()->Time(); },       {"event"}) 
    .Define("S2_paddle_num", [](TEventHandler *event) 
    { return event->GetS2Hit()->Paddle(); },     {"event"}) 
    .Define("S2_pmtDiff_raw", [](TEventHandler *event) 
    { return event->GetS2Hit()->DeltaT_raw(); }, {"event"}) 
    .Define("S2_isTwinHit", [](TEventHandler *event) 
    { return event->GetS2Hit()->Is_twinHit(); }, {"event"}) 

    .Define("LHRS_p", [momentum_LHRS](){ return momentum_LHRS; }, {})
    .Define("RHRS_p", [momentum_RHRS](){ return momentum_RHRS; }, {})
    
    .Define((string)TString(arm+"_cer_t_c"),
	    [](RVecD v){ return v; }, {(string)arm+".cer.t_c"})
    
    .Define((string)TString(arm+"_cer_a_c"),
	    [IN_isRightArm](RVecD v){
	      RVecD outVec; 
	      for (double x : v) outVec.push_back( x );
	      return outVec; },
	    {(string)arm+".cer.a_c"})
    
    .Define(isRightArm?"R_ps_a_c":"L_prl1_a_c",
	    [](RVecD v){ return v; }, {(isRightArm ? "R.ps.a_c":"L.prl1.a_c")})
    
    .Define(isRightArm?"R_sh_a_c":"L_prl2_a_c",
	    [](RVecD v){ return v; }, {(isRightArm ? "R.sh.a_c":"L.prl2.a_c")}) 

    .Define("r_BPMA",
	    [](double x, double y){ return TVector2(x,y); },
	    {(string)arm+"rb.BPMA.x", (string)arm+"rb.BPMA.y"}) 
    .Define("r_BPMB",
	    [](double x, double y){ return TVector2(x,y); },
	    {(string)arm+"rb.BPMB.x", (string)arm+"rb.BPMB.y"}) 
    
    .Define("beam_current",    [](double x){ return x>-1e5 ? x : -1e30; }, 
	    {"hac_bcm_average"})
    
    .Define("beam_MeV", [](double x){ return x; }, {"HALLA_p"}) 
    
    .Define("Raster2_current_y", [](double x){ return x; },
	    {TString(arm+"rb.Raster2.rawcur.y").Data()})     
    .Define("Raster2_current_x", [](double x){ return x; },
	    {TString(arm+"rb.Raster2.rawcur.x").Data()}) 
    
    //if we're in vertical wire mode, then compute the reaction vertex. 
    .Define("react_vertex", [&reactVertex](double rast_Xcurr,
					   double rast_Ycurr)
    { return reactVertex.Compute_reactVertex( rast_Xcurr, rast_Ycurr ); },
	    {TString(arm+"rb.Raster2.rawcur.x").Data(),
	     TString(arm+"rb.Raster2.rawcur.y").Data()}) //*/ 
    
    .Define("react_vertex_TCS", [IN_isRightArm](TVector3 r)
    { APEX::Coord::TransformTV3(r,APEX::Coord::kHall_to_Target,IN_isRightArm);
      return r; 
    }, {"react_vertex"})
    
    //.Define("beam_center", [&reactVertex](){ return reactVertex.Get_beamCenter(); }, {}) 
    .Define("event_number",    [](uint x){ return x; }, {"fEvtHdr.fEvtNum"});

  
  //get track output branches
  output_node 
    = generate_track_output( isRightArm, output_node, out_branches, hrs ); 
    
  
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
    
  int count_1refTrack  = (int)*trackNode.Count(); 
    
  int count_1group     = (int)*nPass_1group; 
  int count_1pair      = (int)*nPass_1pair; 
  int count_1rawTrack  = (int)*nPass_1rawTrack; 
    
  int count_coinc      = (int)*nEvents_coinc.Count(); 
  
  double net_time = stopwatch->RealTime(); 
  
  
  
  
  
  cout << TString::Format("nCoinc: %i, passed...", count_coinc) << endl; 
  
  cout << (isRightArm?"Right":"Left") 
       << " arm (fractions relative to total coinc count) ~~~~~~~~~~~~~~~~~~~" 
       << endl; 
  
  cout << 
    TString::Format(" 1-group:     %7i (%0.4f)", 
		    count_1group, 
		    ((double)count_1group)/((double)count_coinc))    << endl; 
  
  cout << 
    TString::Format(" 1-pair:      %7i (%0.4f)", 
		    count_1pair, 
		    ((double)count_1pair)/((double)count_coinc))     << endl; 
  
  cout << 
    TString::Format(" 1-raw-track: %7i (%0.4f)", 
		    count_1rawTrack, 
		    ((double)count_1rawTrack)/((double)count_coinc)) << endl; 
  
  cout << 
    TString::Format(" 1-ref-track: %7i (%0.4f)", 
		    count_1refTrack, 
		    ((double)count_1refTrack)/((double)count_coinc)) << endl
       << endl; 
  
  //print beam-current information
  cout << 
    TString::Format("--Average beam current: %2.2f \n--(fraction of events with non-null beam-current = %0.4f)",
		    h1d_beamCurrent->GetMean(), 
		    ((double)count_nonNullBC)/((double)*nTotalEvents) ) << endl; 
    
  cout << 
    TString::Format( "Time taken: %3.2f (%2.4f ms/rawEvent)",
		     net_time, 
		     1e3*net_time/((double)maxEntries) ) << endl; 
  

#if MAKE_OUTPUT
  return;
#endif 
  
#if !IS_MONTECARLO 
  
  //non-monte-carlo data-processing
  TString canv_title = TString::Format("run %4i, %2.2e events,",
				       runNum, (double)maxEntries); 
  
  canv_title += isRightArm ? "Right arm" : "Left arm"; 

  //cute litte function to make a new TCanvas with the properties we want
  auto new_TCanvas = [canv_title](bool showStats=false)
  {
    new TCanvas;
    gPad->SetTitle(canv_title);
    gStyle->SetOptStat(showStats);
    return;
  };
  
  
  string branch_track = isRightArm ? "tracks_R_refined" : "tracks_L_refined"; 
  
  const double tauCut_center =  0.;
  const double tauCut_radius = 20e-9; 

  const double w_bound = 22.5e-3; 
  
  
  auto h2d_xParam_Dt = trackNode
    
    .Define("xParam", [](ROOT::RVec<TvdcTrack*> tracks) 
    {
      RVecD out; out.reserve(tracks.size()); 
      for (TvdcTrack *trk : tracks) out.push_back( trk->xParam()*0.13975 ); 
      return out;
    }, {branch_track}) 
    
    .Define("Dt",      [](ROOT::RVec<TvdcTrack*> tracks) 
    {
      RVecD out; out.reserve(tracks.size()); 
      for (TvdcTrack *trk : tracks) out.push_back( trk->T0()*1e9 + 14.92 ); 
      return out;
    }, {branch_track})
    
    .Histo2D({"xp_dt", ";x_tr - x_s2 (m);T_{S2}-T_{Track} (ns)", 
	200, -2.5*0.14, 2.5*0.14, 200, -40, 40}, "xParam", "Dt"); 

  new_TCanvas();
  auto h2d_xparam_clone = (TH1D*)h2d_xParam_Dt->Clone();
  h2d_xparam_clone->Draw("col"); 
  return; 
  
  
  auto h1d_tau_error = output_node

    .Define("tau_error", [](ROOT::RVec<TvdcTrack*> tracks)
    {
      ROOT::RVec<double> out; 
      for (TvdcTrack *trk : tracks) { 
	for (int p=0; p<4; p++) { 
	  auto group = trk->GetGroup(p); 
	  for (int h=0; h<group->Nhits(); h++) {
	    
	    double w = fabs((group->WirePos(h) - trk->Intercept(p))/trk->Slope(p)); 
	    
	    double err
	      = group->Time(h)
	      - (trk->T0() + trk->Get_T_model(p,group->WirePos(h)));
	    
	    
	    out.push_back( err ); 
	  }
	}
      }
      return out;      
    }, {branch_track})

    .Histo1D({"", "tau-error", 200, -45e-9, 45e-9}, "tau_error"); 

  
  new_TCanvas(); 
  h1d_tau_error->DrawCopy();

  
  auto test_wireEff = output_node
    
    .Define("w_hits", [tauCut_center,tauCut_radius](ROOT::RVec<TvdcTrack*> tracks)
    {
      const double dWire = THRS::VDC::wireSpacing; 
      const int max_wires = 7; 
      
      ROOT::RVec<pair<double,bool>> out;  
      for (TvdcTrack *trk : tracks) { 
	for (int p=0; p<4; p++) { 

	  auto pair = p>=2 ? trk->GetPair_Hi() : trk->GetPair_Lo(); 

	  //get closest wire pos
	  double v0 = pair->ClosestWirePos(trk->Intercept(p)); 
	  
	  int wire0 = THRS::VDC::WireNum(trk->IsRightArm(), p, v0); 
	    	  
	  double vHi = v0 + ((double)max_wires)*dWire;

	  map<int,bool> hitMap; 
	  
	  //build our map (which we will use to keep track of hits
	  for (int wire=max_wires+wire0; wire>=-max_wires+wire0; wire--)   
	    if (wire>0 && wire<=368)
	      hitMap[wire] = false;
	  	  
	  
	  auto group = trk->GetGroup(p); 
	  for (int h=0; h<group->Nhits(); h++) {
	    
	    int wire = group->WireNum(h);

	    auto findWire = hitMap.find(wire);

	    //wire is not in map (skip); 
	    if (findWire==hitMap.end()) continue;
	      
	    //compute w	    
	    double err
	      = group->Time(h)
	      - (trk->T0() + trk->Get_T_model(p,group->WirePos(h)));
	    
	    //mark this wire as having a valid hit
	    if (fabs(err)<tauCut_radius) findWire->second = true; 
	  }
	  
	  //now check to see which wires have hits found
	  for (auto it = hitMap.begin(); it != hitMap.end(); ++it) {
	    
	    //wire pos
	    double v = THRS::VDC::WirePos(trk->IsRightArm(), p, it->first); 
	    
	    double w = fabs( (v - trk->Intercept(p))*trk->Slope(p) );
	    
	    out.push_back( {w, it->second} ); 
	  }
	  
	}//for (int p=0; p<4; p++) 
      }//for (TvdcTrack *trk : tracks)
      
      return out;
    }, {branch_track})

    .Define("hits_found", [](ROOT::RVec<pair<double,bool>> hits)
    {
      RVecD out; 
      for (pair <double,bool> hit : hits) if (hit.second) out.push_back( hit.first ); 
      return out;
    }, {"w_hits"}) 
    
    .Define("hits_all", [](ROOT::RVec<pair<double,bool>> hits)
    {
      RVecD out; 
      for (pair <double,bool> hit : hits) out.push_back( hit.first ); 
      return out;
    }, {"w_hits"}); 
  
  
  new_TCanvas(); 

  const int nBins = 125;
  const double w_max = 22.5e-3; 
  
  auto h1d_wHits_all = test_wireEff
    .Histo1D({"", "all", nBins, 0, w_max}, "hits_all");
  
  h1d_wHits_all->SetLineColor(kRed); 
  h1d_wHits_all->GetYaxis()->SetRangeUser( 0, h1d_wHits_all->GetMaximum()*1.20 ); 
  h1d_wHits_all->DrawCopy(); 
  
  auto h1d_wHits_found = test_wireEff
    .Histo1D({"", "found", nBins, 0, w_max}, "hits_found");
  
  h1d_wHits_found->DrawCopy("SAME"); 


  new_TCanvas();
  auto h1d_eff = new TH1D("wireEff", "wire-#eta vs X_{drift};m;", nBins, 0, w_max); 

  for (int bin=1; bin<=h1d_eff->GetXaxis()->GetNbins(); bin++) {

    if (h1d_wHits_all->GetBinContent(bin) < 1.) continue;
    
    h1d_eff->SetBinContent(bin,
			   h1d_wHits_found->GetBinContent(bin)/
			   h1d_wHits_all  ->GetBinContent(bin) );
    
    h1d_eff->SetBinError(bin,
			 sqrt( h1d_wHits_found->GetBinContent(bin) +
			       h1d_wHits_all  ->GetBinContent(bin) )
			 /h1d_wHits_all->GetBinContent(bin) ); 
  }

  new TFile(TString::Format("hist_vdcEfficiency_%i_%s.root",
			    runNum,isRightArm?"R":"L"), "RECREATE");
  
  auto h1d_eff_save = (TH1D*)h1d_eff->Clone(); 
  
  h1d_eff_save->GetYaxis()->SetRangeUser( 0, h1d_eff->GetMaximum()*1.20 );
  h1d_eff_save->SetLineColor(kBlue); 
  h1d_eff_save->Draw("HIST L"); 

  h1d_eff_save->Write(); 
  
  
  
  auto Compute_dv = [](TvdcTrack *trk, int p) { 
    
    bool arm = trk->IsRightArm(); 
    
    int wire       = TvdcHit::WireNum( arm, p, trk->Intercept(p) ); 
    
    double wirePos = TvdcHit::WirePos( arm, p, wire ); 
    
    return trk->Intercept(p) - wirePos; 
  }; 
  
  auto h1d_dv = output_node
    
    .Define("dv", [&Compute_dv](ROOT::RVec<TvdcTrack*> tracks) 
	    { 
	      RVecD out; 
	      for (TvdcTrack *trk : tracks)
		for (int p=0; p<4; p++) 
		  out.push_back( Compute_dv(trk,p) ); 
	      
	      return out;
	    }, {branch_track})
    
    .Histo1D({"dv", ";dv (m);", 100, -dWire/2, dWire/2}, "dv"); 
  
  new_TCanvas(); 

  h1d_dv->GetYaxis()->SetRangeUser( 0, h1d_dv->GetMaximum()*1.20 ); 
  double dv_error = compute_minVariance( (TH1D*)h1d_dv->Clone() );

  h1d_dv->SetTitle(TString::Format("dv = %0.5f mm;x_{track}-x_{p.wire} (m);",
				   sqrt(dv_error)*1e3)); 
  h1d_dv->DrawCopy("E");
    
  return; 
  

  
  
  //cerenkov analysis

  //basic containter for cerenkov hit information
  const int CER_nCells = 10; 
  
  struct TCerenkovHit {

    TCerenkovHit() {};
    TCerenkovHit( const double itr_x,
		  const double itr_y,
		  const double itr_time,
		  const int    icellNum,
		  const double iadc,
		  const double itime )
      : tr_x(itr_x),
	tr_y(itr_y),
	tr_time(itr_time),
	cell(icellNum),
	adc(iadc),
	time(itime) {};

    ~TCerenkovHit(); 
    
    double tr_x,tr_y,tr_time,adc,time;
    int cell; 
  };

      
  
  
  //find which cerenkov-cell a track is connected to
  auto find_cerenkovCell = [IN_isRightArm](TTrackCoords *trk) {

    const double z_cer = IN_isRightArm ?  1.800  :  1.810;
    const double x0    = IN_isRightArm ? -0.9495 : -0.9188;

    const double cell_width = 0.3829;

    TVector3 r_cer( trk->x + trk->dxdz*z_cer,
		    trk->y + trk->dydz*z_cer,
		    z_cer ); 
    
    //find x-number
    int cell = (int)TMath::Nint( (r_cer.X() - x0)*(2./cell_width) ); 

    if (r_cer.Y() > -1.612e-2) cell += 1;
    
    if (cell<0) cell=0;
    if (cell>9) cell=8;

    return cell;
  }; 


  /*
  struct TShowerInfo {

    enum EDetType { kPreShower=0, kShower };
    
    TShowerInfo() {};
    TShowerInfo(const bool isRHRS, 
		RVecD adc_sh,
		RVecD adc_ps ) : f_isRHRS(isRHRS) {
      
      fZ_preShower = isRHRS ? 3.500 : 4.580;
      fZ_shower    = isRHRS ? 3.640 : 4.770;
      
      fAdc_sh = adc_sh;
      fAdc_ps = adc_ps; 
    }
    
    ~TShowerInfo();

    
    //test pos-model
    int Get_closestBlockId( const TVector3 pos,
			    const EDetType type ) const { 

      const int nRows = f_isRHRS
	? (type==kPreShower ? 24 : 15)
	: (type==kPreShower ? 17 : 17); 
            
      const int nCols = f_isRHRS
	? (type==kPreShower ? 2 : 5)
	: (type==kPreShower ? 2 : 2); 
      
      const double cell_width = f_isRHRS
	? (type==kPreShower ?  0.1011 :  0.1513)
	: (type==kPreShower ? -0.1476 : -0.1480); 
      
      const double cell_x0    = f_isRHRS
	? (type==kPreShower ? -1.2267 : -1.0800)
	: (type==kPreShower ?  0.9829 :  1.0169); 

      int cell; 
      
      //now, compute the cell using y-pos
      if (f_isRHRS) { }
      else          { 
	
	int cell = (int)TMath::Nint( (pos.X()-cell_x0)/cell_width ); 
	
	const double cell_yEdge = type==kPreShower ? -3.044e-2 : 2.328e-2;

	if (pos.Y() < cell_yEdge) {            }; 
	
	//start here. 
      }
            
      //for the L-shower (prl2), blocks 2 & 3 are swapped. 
      if (!isRHRS && !is_preShower) { 
	if (b_num==2) { b_num=3; }
	else          { if (b_num==3) b_num=2; }
      }

      //once the block gets to the halfway point, you get to the next row
      b_num = b_num % n_rows;
    
      return cell_x0 - cell_width*((double)b_num);    
      return -1000; 
    }; 
  
    
    TVector3 Get_track_at_det(TvdcTrack *trk, EDetType type) const {

      auto coord = new TTrackCoords( trk->FP_x(),
				     trk->FP_y(),
				     TMath::Tan(trk->Theta()),
				     TMath::Tan(trk->Phi()) );
      
      //project the track-coords onto the det-plane
      double z = type==kPreShower ? fZ_preShower : fZ_shower;
      
      return TVector3( coord->x + coord->dxdz*z,
		       coord->y + coord->dydz*z, z );
    }

    //get the adc of the track at the detector
    void Get_track_adc( TvdcTrack *trk,
			double &adc_ps,
			double &adc_sh ) const {
      
      TVector3 pos = Get_trackInt_at_det( trk, kPreShower );
    }

    
    
  private:
    
    double fZ_shower, fZ_preShower;
    RVecD fAdc_sh, fAdc_ps;
    bool f_isRHRS; 
    }; */ 
		
  
  auto get_cerenkov_adc = [&find_cerenkovCell](TvdcTrack *trk,
					       RVecD cer_adc) {
    
    int cell = find_cerenkovCell( new TTrackCoords( trk->FP_x(), 
						    trk->FP_y(), 
						    TMath::Tan(trk->Theta()), 
						    TMath::Tan(trk->Phi()) ) ); 
    double adc = cer_adc.at(cell); 
    
    //add the other-pmt
    if (cell%2==0) { adc += cer_adc.at(cell+1); }
    else           { adc += cer_adc.at(cell-1); }
    
    return adc; 
  }; 
    

  const double det_z = isRightArm ? 1.800 : 1.810;
  
  auto node_cerenkov = output_node

    .Define("cer_hits",
	    [det_z, &find_cerenkovCell](ROOT::RVec<TvdcTrack*> tracks,
					RVecD cer_tdc,
					RVecD cer_adc, 
					double t_s2){
	      
	      ROOT::RVec<TCerenkovHit*> hits;
	      for (TvdcTrack *trk : tracks) {
		
		TVector3 S( TMath::Tan( trk->Theta() ),
			    TMath::Tan( trk->Phi() ), 1.);
		
		TVector3 r0( trk->FP_x(),
			     trk->FP_y(), 0. );
		
		int cell = find_cerenkovCell( new TTrackCoords( r0.X(),
								r0.Y(),
								S.X(),
								S.Y() ) );

		//compute this track's intercept witht the cerenkov-plane
		TVector3 r_cer( r0.X() + S.X()*det_z,
				r0.Y() + S.Y()*det_z,
				det_z );

		//now, compute the time-to-flight to the cerenkov-plane
		S = S.Unit() * 2.99e8; //we can now use this vector to compute speed

		double ToF = fabs( (r_cer.X() - r0.X())/S.X() );
		
		hits.push_back( new TCerenkovHit( r_cer.X(),
						  r_cer.Y(),
						  t_s2 + ToF,
						  cell,
						  cer_adc.at(cell),
						  cer_tdc.at(cell) ) );

	      }//for (TvdcTrack *trk : tracks) {
	      return hits; },
	    { branch_track,
	      TString(arm+".cer.t_c").Data(),
	      TString(arm+".cer.a_c").Data(),
	      "S2_paddle_realTime" });

  
  auto h2d_cerenkov_timing = node_cerenkov

    .Define("t_offset", [](ROOT::RVec<TCerenkovHit*> hits) {
      RVecD outVec;
      for (TCerenkovHit *hit : hits) outVec.push_back( hit->time - hit->tr_time );
      return outVec; }, {"cer_hits"})

    .Define("hit_cell", [](ROOT::RVec<TCerenkovHit*> hits) {
      ROOT::RVec<int> outVec;
      for (TCerenkovHit *hit : hits) outVec.push_back( hit->cell );
      return outVec; }, {"cer_hits"})
    
    .Histo2D({"", "timing vs cell;cell no.;timing offset",
	CER_nCells, -0.5, ((double)CER_nCells)-0.5,
	60, -1.45e-6, -1.34e-6}, "hit_cell", "t_offset");
  
  
  new TCanvas;
  gStyle->SetOptStat(0);

  h2d_cerenkov_timing->DrawCopy("col");
  
  for (int c=1; c<=CER_nCells; c++) {
    auto proj = h2d_cerenkov_timing->ProjectionY("proj", c,c);

    if (proj->Integral()<5.) {
      cout << TString::Format("cell %i insufficient stats (%i)",
			      c,(int)proj->Integral() ) << endl; 
      continue;
    }
    
    auto fitResult = proj->Fit("gaus", "S Q N");
    
    cout << TString::Format("cell,avg,gaus-average = {%i, %+0.4e, %+0.4e}",
			    c,
			    proj->GetMean(),
			    fitResult->Parameter(1) ) << endl;
  }
    

  return; 
  
  
  auto h1d_cerenkov = output_node
    
    .Define((string)TString(arm+"_tr_cer_adc"),
	    [&get_cerenkov_adc] (ROOT::RVec<TvdcTrack*> tracks,
				 RVecD cer_adc) {
	      RVecD out;
	      for (TvdcTrack *trk : tracks) {

		double adc = get_cerenkov_adc(trk, cer_adc);
		
		out.push_back( adc + (isRHRS ? 7.838e3 : 0) );
	      }
	      return out; }, {branch_track, "cerenkov_a_c"})
    
    .Histo1D({"", "cerenkov-adc associated with track", 200, -500,10e3},
	     TString(arm+"_tr_cer_adc").Data());
  
  new TCanvas;
  gPad->SetTitle(canv_title);  
  gStyle->SetOptStat(0);

  h1d_cerenkov->DrawCopy();
  return; 
	  
      
      

    
  
  


  
  auto h2d_eta = trackNode
    
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
  
  
  
  
  
  auto h2d_phi_s2y = trackNode 
        
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
  

  
  
  auto h2d_theta_s2x = trackNode
        
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
  
  
  
  
  auto h2d_slope_Dt = trackNode
    
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
  
  
  	     
  
  auto h2d_w_tau = trackNode
    
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
    
    
  auto h1d_nGoodPts_plane = trackNode
    
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
    
	       

  auto h1d_int_error = trackNode 
    
    .Define("intercept_err", [](ROOT::RVec<TvdcTrack*> tracks) 
	    { 
	      RVecD out; 
	      for (int t=0; t<tracks.size(); t++) { auto trk = tracks.at(t); 
		for (int p=0; p<4; p++) 
		  out.push_back( trk->Error_intercept(p) ); 
	      }
	      return out; }, {branch_track})
    
    .Histo1D({"int_err", "Intercept-error (m)", 200, 0, 0.3e-3}, "intercept_err"); 
  
  
  
  
  
  
  
  auto diff_RVecD = [](RVecD v1, RVecD v2) { 
    
    RVecD vCompare; 
    
    for (double x1: v1) 
      for (double x2: v2) 
	vCompare.push_back( x1-x2 ); 
    
    return vCompare; 
  }; 


  const double Cerenkov_z = 1.81; 
  
  auto h1d_cerenkov_timeCheck = output_node
    
    .Define("cerenkov_Dt", [Cerenkov_z](ROOT::RVec<TvdcTrack*> track_vec,
					double t_s2,
					RVecD cer_vec) {
      RVecD v;
      
      for (TvdcTrack *trk: track_vec) { 
	for (double t_cer: cer_vec)   {
	  
	  if (abs(t_cer)>1e20) continue; 
	  
	  v.push_back( t_cer - (trk->GetTimeAtZ(Cerenkov_z) + t_s2) );
	}
      }
      return v; },
	    { TString("tracks_"+arm+"_refined").Data(),
	      "S2_paddle_realTime",
	      TString(arm+".cer.t_c").Data() })

    .Define("cerenkov_pmtNum", [](ROOT::RVec<TvdcTrack*> track_vec,
				  double t_s2,
				  RVecD cer_vec) {
      RVecD v;
      
      for (TvdcTrack *trk: track_vec) { 
	for (uint n=0; n<cer_vec.size(); n++)   {
	  
	  double t_cer = cer_vec.at(n); 
	  
	  if (abs(t_cer)>1e20) continue; 
	  
	  v.push_back( n );
	}
      }
      return v; },
	    { TString("tracks_"+arm+"_refined").Data(),
	      "S2_paddle_realTime",
	      TString(arm+".cer.t_c").Data() })

    
    .Histo2D({"cer-dt", "pmt-number;t_{cer} - t_{track} (s)",
	10, -0.5, 9.5, 200, -1.40e-6, -1.32e-6},
      "cerenkov_pmtNum", "cerenkov_Dt");

  
  auto canv_1 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);  
    
  h1d_cerenkov_timeCheck->DrawCopy("col");
  
  return; 
	    
  

  /*/find the difference between all possible pairs of values in two RVecD's 
    
  auto h2d_coordCheck = output_node 
    
    .Define("err_x", diff_RVecD, 
	    {TString(arm+"_tr_tg_dxdz").Data(),TString(arm+".tr.tg_th").Data()})
    .Define("err_y", diff_RVecD, 
	    {TString(arm+"_tr_tg_dp").Data(),TString(arm+".tr.tg_dp").Data()})
    
    .Histo2D({"coord_err", ";x-err (new-old) (m);y-err (new-old) (m)", 
	  200, -0.1,0.1, 200, -0.1,0.1}, "err_x", "err_y"); 
  
  
  h2d_coordCheck->DrawCopy("col");   
	
  return;*/ 

  

  /*auto canv_split = new TCanvas("c", "target coord error", 1200, 800);  
  gStyle->SetOptStat(0); 

  canv_split->Divide(2,2); 
  
  auto h2d_coordCheck__y  = output_node
    .Histo2D({"check-y", "y", 200, -0.05,0.05, 200, -0.05,0.05},
	     TString(arm+".tr.tg_y").Data(),
	     TString(arm+"_tr_tg_y").Data() ); 
  
  canv_split->cd(1); 
  h2d_coordCheck__y->DrawCopy("col");   
  
  auto h2d_coordCheck__th = output_node
    .Histo2D({"check-th", "dx/dz", 200, -0.05,0.05, 200, -0.05,0.05},
	     TString(arm+".tr.tg_th").Data(),
	     TString(arm+"_tr_tg_dxdz").Data() ); 

  canv_split->cd(2); 
  h2d_coordCheck__th->DrawCopy("col"); 
  
  auto h2d_coordCheck__ph = output_node
    .Histo2D({"check-ph", "dy/dz", 200, -0.05,0.05, 200, -0.05,0.05},
	     TString(arm+".tr.tg_ph").Data(),
	     TString(arm+"_tr_tg_dydz").Data() ); 
  
  canv_split->cd(3); 
  h2d_coordCheck__ph->DrawCopy("col"); 

  auto h2d_coordCheck__dp = output_node
    .Histo2D({"check-dp", "#delta p", 200, -0.05,0.05, 200, -0.05,0.05},
	     TString(arm+".tr.tg_dp").Data(),
	     TString(arm+"_tr_tg_dp").Data() ); 

  canv_split->cd(4); 
  h2d_coordCheck__dp->DrawCopy("col"); 
  
  return; */
  
  
  	   
  
  auto canv_2 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);  
    
  //draw x-param plot
  h2d_xParam_Dt->DrawCopy("col"); 
  
  draw_vLine( (TH2D*)h2d_xParam_Dt->Clone(), -TRK_CUT_xParam ); 
  draw_vLine( (TH2D*)h2d_xParam_Dt->Clone(),  TRK_CUT_xParam ); 
    
  draw_hLine( (TH2D*)h2d_xParam_Dt->Clone(), -TRK_CUT_Dt ); 
  draw_hLine( (TH2D*)h2d_xParam_Dt->Clone(),  TRK_CUT_Dt ); //*/ 
  
  
  
  /*h2d_slope_Dt->DrawCopy("col"); 

  auto graph = th2d_to_graph( (TH2D*)h2d_slope_Dt->Clone(), 1.16, 1.70 ); 
  
  graph->Draw("SAME"); */   
  
  h1d_dv->SetTitle(TString::Format("DV-error = %3.5f mm", 
				   TMath::Sqrt(dv_error)*1e3)); 
  
  h1d_dv->GetYaxis()->SetRangeUser(0, h1d_dv->GetMaximum()*1.20); 
  h1d_dv->DrawCopy("E");
    
  
  //h2d_xParam_Dt->DrawCopy("col"); 
  
  //h2d_theta_s2x->DrawCopy(); 
  
  
  //h2d_phi_s2y->DrawCopy(); 
  //
  //h2d_eta->DrawCopy("col"); //*/
  /*h2d_slope_Dt->DrawCopy("E"); 
    
    draw_gausFit_fixedBase( (TH1D*)h2d_slope_Dt->Clone(), 10e-9 ); */ 
  
  auto canv_3 = new TCanvas; 
  gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);  
  
  
  auto canv_4 = new TCanvas; 
  //gStyle->SetOptStat(0); 
  gPad->SetTitle(canv_title);  
  
  h1d_int_error->DrawCopy();   
  
  return; 
  
  
  
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
