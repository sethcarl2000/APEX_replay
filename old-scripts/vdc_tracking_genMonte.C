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


double f_model_( double *x, double *par ) { 
  
  auto dummy_s2R = new TS2Hit(true, 0,200,200); 
  auto dummy_s2L = new TS2Hit(false,0,200,200); 
    
  auto dummy_event = new TEventHandler(false, 53., dummy_s2R, dummy_s2L); 
  
  return dummy_event->Drift_T(x[0], 1.4, 2 );
}  
void vdc_tracking_genMonte(int maxEntries=5e4,
			   int runNum=4329, 
			   TString path_outFile //path to out file (if used)
			   ="replay/monte-carlo/test", 
			   //="replay/monte-carlo/multiTrack_sig-9ns_wireEff-80",  
			   TString path_inFile="replay/out")
{
  ////////////////////////////////////////
  ////////////////////////////////////////
  //    
  
  cout << "processing RDataFrame" << endl; 
  
  path_inFile  += "."+TString::Itoa(runNum, 10)+".root";
  path_outFile += "."+TString::Itoa(runNum, 10)+".root"; 
  
  
  /***
   *    Takes Raw S2-m data from the decoded-data tree, and converts them to 
   *     the objects called 'TS2Hit*', stored in an RVec.  
  ***/ 
  auto Generate_S2_hits = []( bool arm, RVecD PMT_R, RVecD PMT_L ) { 
    
    //generate a vector of all coinc s2-paddle hits
    ROOT::VecOps::RVec<TS2Hit*> coincHits; 
  
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
  
  
  //ROOT::EnableImplicitMT(); 
  
  
  ROOT::RDataFrame d_T("T", path_inFile.Data() ); 
  
  auto h_dt_temp = d_T
    .Range(0,5e4)
    .Define("dt", get_dt, {"R.s2.rt","R.s2.lt","L.s2.rt","L.s2.lt"})
    .Histo1D({"", "dt", 200, -60e-9, 60e-9}, {"dt"}); 
  
  //h_dt_temp->DrawCopy(); 
  
  TVector2 dt_fitResult = draw_gausFit( (TH1D*)h_dt_temp->Clone(), 
					6e-9, -1e30, -1 ); 
  
  const double tR_tL_diff  = dt_fitResult.X(); 
  const double tR_tL_sigma = dt_fitResult.Y(); 
  
  const double CUT_S2Sep_stdDev = 6; 
    
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
      ROOT::VecOps::RVec<TS2Hit*> hits_R,
      ROOT::VecOps::RVec<TS2Hit*> hits_L ) {
    
    ROOT::VecOps::RVec<TEventHandler*> coincEvents; 
    
    for (int iR=0; iR<hits_R.size(); iR++) { 
      for (int iL=0; iL<hits_L.size(); iL++) { 
	
	double TR_TL = hits_R.at(iR)->Time() - hits_L.at(iL)->Time(); 
	
	if ( TMath::Abs(TR_TL-tR_tL_diff)/tR_tL_sigma < CUT_S2Sep_stdDev ) 
	  
	  //this is a bona-fide both-arm S2-coinc; let's keep it
	  coincEvents.push_back( new TEventHandler(false, 
						   beamCurrent, 
						   hits_R.at(iR), 
						   hits_L.at(iL) ) );  
	
      }//for (int iL=0; iL<hits_L.size(); iL++) { 
    }
    
    return coincEvents;
  }; 
  
  
  //now actually perform the analysis chain
  auto Phi_model   = [](double S2_y) { 
    
    return 
    -0.231*TMath::Power(S2_y,2) 
    +  0.2532*S2_y 
    +  0.539e-3; 
  };
  auto Theta_model = [](double S2_x) { 
    
    return 0.109648*S2_x; 
  };
 
  auto Gen_monteTrack = [&Theta_model, &Phi_model](TEventHandler *evt) {
    
    const double M_S2_x_max 
    = evt->GetS2Hit()->X() + evt->GetS2Hit()->PaddleWidth()*0.5; 
    
    const double M_S2_x_min
    = evt->GetS2Hit()->X() - evt->GetS2Hit()->PaddleWidth()*0.5; 
	
    const double M_S2_y_max =  0.17; 
    const double M_S2_y_min = -0.10; 
    
    const double M_ph_min = -0.018;//  + 0.5e-3; 
    const double M_ph_max =  0.008;//  - 0.5e-3; 
	
    const double M_th_min = -0.0200;// + 0.5e-3; 
    const double M_th_max =  0.0175;// - 0.5e-3; 
        
    const double M_vu_Lo_diff_min = -0.082;
    const double M_vu_Lo_diff_max =  0.057; 
        
    
    auto track = new TvdcTrack( evt ); 
    
    while (1) { 
      
      double S2_x 
	= M_S2_x_min + (M_S2_x_max - M_S2_x_min)*gRandom->Rndm(); 
	
      double S2_y =
	+ M_S2_y_min + (M_S2_y_max - M_S2_y_min)*gRandom->Rndm(); 
	
      double Theta = Theta_model( S2_x ) + 
	M_th_min + (M_th_max - M_th_min)*gRandom->Rndm(); 
	
      double Phi   = Phi_model( S2_y ) + 
	M_ph_min + (M_ph_max - M_ph_min)*gRandom->Rndm(); 
      
      
      track->Set_S2int_angles( S2_x, S2_y, Theta, Phi ); 
      //the '0' on the end is the offset
      
      double vu_Lo_diff = track->Intercept(0)-track->Intercept(1); 
      
      if (M_vu_Lo_diff_min < vu_Lo_diff &&
	  M_vu_Lo_diff_max > vu_Lo_diff ) break; 
        
    }
          
    return track; 
  }; 
    
  
    
    

  auto Gen_monteHits = []( TvdcTrack *trk,
			   int plane ) { 
    
    auto group = new THitGroup(plane); 
    //generate monte-carlo TDC hits for each plane
    bool arm = trk->GetEvent()->ActiveArm(); 
    
    double intercept = trk->Intercept(plane); 
        
    double T_s2 = trk->GetEvent()->GetS2Hit()->Time(); 
    
    const double M_tau_sigma = 9e-9; //sigma of TDC hits
    const double M_w_max = 20e-3; //max drift-dist
    
    const double M_efficiency = 0.75; //prob of an in-range wire having a hit
    
    double slope = trk->Slope(plane); 
    
    int wireNum    = TvdcHit::WireNum( arm, plane, intercept +M_w_max/slope ) - 2; 
        
    double wirePos = TvdcHit::WirePos( arm, plane, wireNum ); 
    
    while( (wirePos-intercept)*slope > -M_w_max ) { 
      
      wirePos += -dWire; 
      wireNum += 1; 
          
      double w = TMath::Abs( slope*(wirePos-intercept) ); 
      
      double prob_hit = M_efficiency/( TMath::Exp(1400.*(w - 14.2e-3)) + 1. ); 
      
      //random check to see if this hit should be generated
      if ( gRandom->Rndm() > prob_hit ) continue; 
      
      
      double realTime
	= trk->Get_T_model( plane, wirePos ) 
	+ T_s2 
	+ gRandom->Gaus()*M_tau_sigma; 
      
      double rawTime = TvdcHit::RawTime( arm, 
					 plane, 
					 wireNum, 
					 realTime ); 
      
      auto hit = new TvdcHit( plane, 
			      (int)wireNum, 
			      rawTime, 
			      trk->GetEvent() );
            
      group->AddHit( hit ); 
      
    }
    return group; 
  }; 
  
  auto Gen_accidentTrack 
    = [&Theta_model, &Phi_model]( ROOT::RVec<THitGroup*> groups,
				  TvdcTrack *trk, 
				  TEventHandler *evt ) { 
    
    //because of the way the NewTrack library is written, it's
    //much easier to do track-operations when you have groups/pairs to deal with
    
    auto pair_Lo = new TChamberPair( true, 
				     trk->Intercept(1), 
				     trk->Intercept(0), 
				     groups.at(0), groups.at(1) ); 
    
    auto pair_Hi = new TChamberPair( false, 
				     trk->Intercept(3), 
				     trk->Intercept(2), 
				     groups.at(2), groups.at(3) ); 
    
    auto new_track = new TvdcTrack( evt, pair_Lo, pair_Hi ); 
    
    
    //now, that we have a 'fully-formed' track, we can screw with it! 
    
    double old_params[5] = { trk->Intercept(0), 
			     trk->Intercept(1), 
			     trk->Intercept(2), 
			     trk->Intercept(3), 
			     trk->T0() }; 
        
    const double ERR_T = 125e-9;
    const double ERR_x = 15.*dWire;     
        
    const double M_ph_min = -0.0180; 
    const double M_ph_max =  0.0080;
	
    const double M_th_min = -0.0200;
    const double M_th_max =  0.0175;
    
    while (1) { 
      
      double new_params[5]; 
      
      
      double nudgeT = ERR_T*( 1. - 2.*gRandom->Rndm() ); 
      
      new_params[4] = old_params[4] + nudgeT; 
            
      for (int p=0; p<4; p++) { 
	
	double nudgeX = ERR_x*( 1. - 2.*gRandom->Rndm() ); 
	
	//make sure this guess isn't too close to the ACTUAL 'correct' track
	while (TMath::Abs(nudgeT) < 60e-9 && 
	       TMath::Abs(nudgeX) < 1.1*dWire) {
	  
	  nudgeX = ERR_x*( 1. - 2.*gRandom->Rndm() ); 
	}
	new_params[p] = old_params[p] + nudgeX; 
      }
              
      //now, see if this accidental track is valid. 
      new_track->Set_params( new_params ); 
	
      //don't let this track have valid angles!!
      double theta_err = new_track->Theta() - Theta_model( new_track->S2_x() ); 
      double phi_err   = new_track->Phi()   - Phi_model( new_track->S2_y() ); 
      
      if ( theta_err > M_th_min &&
	   theta_err < M_th_max &&
	   phi_err   > M_ph_min &&
	   phi_err   < M_ph_max ) break; //keep this track
    }
      
    return new_track; 
  };   
  
  auto Gen_allHits = []( TvdcTrack *trk_good, 
			 TvdcTrack *trk_bad   ) { 
    
    bool arm = trk_good->GetEvent()->ActiveArm(); 
        
    ROOT::RVec<THitGroup*> groups_sorted; 
    
    for (int p=0; p<4; p++) { 
      
      ROOT::RVec<TvdcHit*> all_hits; 
      
      auto good_group = trk_bad->GetGroup(p); 
      //paradoxically, the 'bad' track has access to the 'good' groups 
      
      for (int h=0; h<good_group->Nhits(); h++) 
	all_hits.push_back( good_group->GetHit(h) ); 
        
      
      //generate monte-carlo TDC hits for each plane
      
      double intercept = trk_bad->Intercept(p); 
        
      double T_s2 = trk_bad->T0() + trk_bad->GetEvent()->GetS2Hit()->Time(); 
      
      
      const double M_tau_sigma = 9e-9; //sigma of TDC hits
      const double M_w_max = 20e-3; //max drift-dist
    
      const double M_efficiency = 0.75; //prob of an in-range wire having a hit
    
      double slope = trk_bad->Slope(p); 
    
      int wireNum    = TvdcHit::WireNum( arm, p, intercept +M_w_max/slope ) - 2; 
      
      double wirePos = TvdcHit::WirePos( arm, p, wireNum ); 
    
      while( (wirePos-intercept)*slope > -M_w_max ) { 
      
	wirePos += -dWire; 
	wireNum += 1; 
	
	double w = TMath::Abs( slope*(wirePos-intercept) ); 
	
	double prob_hit = M_efficiency/( TMath::Exp(1400.*(w - 14.2e-3)) + 1. ); 
      
	//random check to see if this hit should be generated
	if ( gRandom->Rndm() > prob_hit ) continue; 
	
	double realTime
	  = trk_bad->Get_T_model( p, wirePos ) 
	  + T_s2 
	  + gRandom->Gaus()*M_tau_sigma; 
	
	double rawTime = TvdcHit::RawTime( arm, 
					   p, 
					   wireNum, 
					   realTime ); 
      
	auto hit = new TvdcHit( p, 
				(int)wireNum, 
				rawTime, 
				trk_bad->GetEvent() );
            
	all_hits.push_back( hit ); 
      }
      
      
      //sort hits
      auto group_sorted = new THitGroup(p); 
      
      int elem = 0; 
      
      int    lo_wire    = 999; 
      double hi_rawtime = -1e30; 
      int    elem_lowest = 0; 
      
      int it_cutoff = 15; 
      int it(0); 
      
      while (all_hits.size() > 0 && it_cutoff > it) { 
	
	it++; 
	 
	int lo_wire(999); double hi_rawtime(-1e30);  
	int lo_wire_index(-1); 
	
	for (int h=0; h<all_hits.size(); h++) { 
	  
	  if (  (all_hits.at(h)->wNum()  < lo_wire) || 
		(all_hits.at(h)->wNum() == lo_wire &&
		 all_hits.at(h)->GetRawTime() > hi_rawtime)  ) { 
	    
	    //either this is the lowest wire found,
	    // OR, it's TIED for the lowest wire, with a larger rawtime. 
	    
	    //add this hit to the sorted group, and remove it from the 
	    // 'unsorted' vector
	    
	    lo_wire    = all_hits.at(h)->wNum(); 
	    hi_rawtime = all_hits.at(h)->GetRawTime(); 
	    
	    lo_wire_index = h; 
	  } 
	  
	  	  
	}//for (int h=0; h<all_hits.size(); h++) 
		
	group_sorted->AddHit( all_hits.at(lo_wire_index) ); 
	
	all_hits.erase( begin(all_hits)+lo_wire_index ); 
		
      }//while (all_hits.size() > 0) 
      
      
      //make sure this hit isn't 'too close' to the last one.
      //  The TDC has this funny property that there's a 'dead-zone'
      //  of about 70ns (in the LHRS) between hits...
      //  therefore, if there are 2 hits on a single wire closer than 70ns, 
      //  we 'push' the higher-time hit up so it's at that 70ns-level. 
      int lastWire = -999;
      double lastRawTime = 1e30; 
      
      double min_buffer = 140; //corresponds to ~70ns in real-time
      
      for (int h=0; h<group_sorted->Nhits(); h++) {
	
	if ( group_sorted->WireNum(h) == lastWire  &&
	     lastRawTime - group_sorted->GetHit(h)->GetRawTime() < min_buffer ) { 
	  
	  group_sorted->GetHit(h)->SetRawTime( lastRawTime - min_buffer ); 
	   
	} 
	
	lastWire    = group_sorted->WireNum(h); 
	lastRawTime = group_sorted->GetHit(h)->GetRawTime(); 
      }    
      
      groups_sorted.push_back( group_sorted ); 
      //now, 
      
    }//for (int p=0; p<4; p++) 
    
    return groups_sorted; 
  }; 
  
  
  TString BRANCH_R_pmtR = "R_s2_rt"; 
  TString BRANCH_R_pmtL = "R_s2_lt"; 
  TString BRANCH_L_pmtR = "L_s2_rt"; 
  TString BRANCH_L_pmtL = "L_s2_lt"; 
  
  TString BRANCH_vdc_rawtime[4] = { "L_vdc_u1_rawtime", 
				    "L_vdc_v1_rawtime", 
				    "L_vdc_u2_rawtime", 
				    "L_vdc_v2_rawtime" };  
  
  TString BRANCH_vdc_wire[4]    = { "L_vdc_u1_wire", 
				    "L_vdc_v1_wire", 
				    "L_vdc_u2_wire", 
				    "L_vdc_v2_wire" };  
  
    
  vector<string> out_branches; 
  
  out_branches.push_back( BRANCH_R_pmtR.Data() ); 
  out_branches.push_back( BRANCH_R_pmtL.Data() ); 
  out_branches.push_back( BRANCH_L_pmtR.Data() ); 
  out_branches.push_back( BRANCH_L_pmtL.Data() ); 
  
  for (int p=0; p<4; p++) { 
    
    out_branches.push_back( BRANCH_vdc_rawtime[p].Data() ); 
    out_branches.push_back( BRANCH_vdc_wire[p].Data() ); 
  }
  
  
  ROOT::EnableImplicitMT(); 
  
  
  auto nEvents_coinc = d_T
    
    //.Range(0,maxEntries)
    
    //find all valid S2 hits for the Rigth arm
    .Define("S2_hits_RHRS", [&Generate_S2_hits](RVecD rt, RVecD lt)
	    { return Generate_S2_hits(true,  rt,lt); }, 
	    {"R.s2.rt","R.s2.lt"})
    
    //find all valid S2 hits for the Left arm
    .Define("S2_hits_LHRS", [&Generate_S2_hits](RVecD rt, RVecD lt)
	    { return Generate_S2_hits(false, rt,lt); }, 
	    {"L.s2.rt","L.s2.lt"})
    
    //Find all possible both-arm S2 coincidences
    .Define("coincEventVec", Gen_S2_coincHits, 
	    {"hac_bcm_average","S2_hits_RHRS","S2_hits_LHRS"})
    
    .Filter([](ROOT::VecOps::RVec<TEventHandler*> v) 
	    { return v.size()==1; }, {"coincEventVec"})
    
    .Define("event", [](ROOT::VecOps::RVec<TEventHandler*> v) 
	    { return v.at(0); },    {"coincEventVec"})
    
    //generate the monte-carlo track
    .Define("track", Gen_monteTrack, {"event"})
    
    //generate monte hit_groups
    .Define("groups_good", [&Gen_monteHits](TvdcTrack *trk) 
	    { 
	      ROOT::VecOps::RVec<THitGroup*> groups; 
	      
	      for (int p=0; p<4; p++) { 
		
		groups.push_back( Gen_monteHits(trk,p) ); 
	      }
	      return groups;
	      
	    }, {"track"})
          
    //check to make sure each group has at least 1 hit
    .Filter([](ROOT::RVec<THitGroup*> groups) {
	
	for (int p=0; p<4; p++) {
	  
	  if ( groups.at(p)->Nhits() < 1 ) return false; 
	}
	return true; }, {"groups_good"})
    
    //generate accidental track
    .Define("track_accidental", Gen_accidentTrack, 
	    {"groups_good","track","event"})
    
    //generate 'full' groups
    .Define("groups", Gen_allHits, 
	    {"track","track_accidental"}) 
    
    
    //now, we're ready to define the output data. 
    
    .Define(BRANCH_R_pmtR.Data(), [](RVecD inBranch) 
	    { return inBranch; }, {"R.s2.rt"})
    .Define(BRANCH_R_pmtL.Data(), [](RVecD inBranch) 
	    { return inBranch; }, {"R.s2.lt"})
    .Define(BRANCH_L_pmtR.Data(), [](RVecD inBranch) 
	    { return inBranch; }, {"L.s2.rt"})
    .Define(BRANCH_L_pmtL.Data(), [](RVecD inBranch) 
	    { return inBranch; }, {"L.s2.lt"}) 
        
    //monte-carlo metadata
    //count the number of good hits
    .Define("MONTE_nGoodHits", [](ROOT::RVec<THitGroup*> groups) { 
	
	ROOT::RVec<int> nHits; 
	
	for (int p=0; p<4; p++) { nHits.push_back( groups.at(p)->Nhits() ); } 
	
	return nHits; }, {"groups_good"})

    
    .Define("MONTE_nAllHits", [](ROOT::RVec<THitGroup*> groups) { 
	
	ROOT::RVec<int> nHits; 
	
	for (int p=0; p<4; p++) { nHits.push_back( groups.at(p)->Nhits() ); } 
	
	return nHits; }, {"groups"})

    
    .Define("MONTE_track_params", [](TvdcTrack *trk) { 
	
	RVecD par; 
	
	for (int p=0; p<4; p++) { par.push_back( trk->Intercept(p) ); }
	
	par.push_back( 0. ); 
	
	return par; }, {"track"}) 
    
    .Define("MONTE_track_params_accident", [](TvdcTrack *trk) { 
	
	RVecD par; 
	
	for (int p=0; p<4; p++) { par.push_back( trk->Intercept(p) ); }
	
	par.push_back( trk->T0() ); 
	
	return par; }, {"track_accidental"}); 


  
  out_branches.push_back( "MONTE_nGoodHits" ); 
  out_branches.push_back( "MONTE_nAllHits" ); 
  out_branches.push_back( "MONTE_track_params" ); 
  out_branches.push_back( "MONTE_track_params_accident" ); 
  
  out_branches.push_back( "hac_bcm_average" ); 
  
  
  //fill the vdc data
  for (int p=0; p<4; p++) { 
    
    nEvents_coinc = nEvents_coinc 
      .Define(BRANCH_vdc_rawtime[p].Data(), 
	      [p](ROOT::VecOps::RVec<THitGroup*> groups) {
		
		RVecD rawtime; 
		
		auto g = groups.at(p); 
		
		for (int h=0; h<g->Nhits(); h++) 
		  rawtime.push_back( g->GetHit(h)->GetRawTime() ); 
		
		return rawtime; }, {"groups"})
      
      .Define(BRANCH_vdc_wire[p].Data(), 
	      [p](ROOT::VecOps::RVec<THitGroup*> groups) {
	      
		RVecD rawtime; 
		
		auto g = groups.at(p); 
		
		for (int h=0; h<g->Nhits(); h++) 
		  rawtime.push_back( g->GetHit(h)->wNum() ); 
		
		return rawtime; }, {"groups"});
    
  }
  
  
  
  /*auto h1d_test = nEvents_1group
    
    .Define("hit_sep", [](RVecD wire, 
			  RVecD rawtime) { 
	      
	      RVecD v; 
	      
	      cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl; 
	      
	      int lastWire = -1; double lastTime; 
	      for (int h=0; h<group->Nhits(); h++) { 
		
		cout << TString::Format(" h %1i   w %3i", 
					h, group->WireNum(h)) << endl;  
		
		if (group->WireNum(h) == lastWire) { 
		  
		  v.push_back( (group->Time(h) - lastTime)*1e9 ); 
	    
		  cout << (group->Time(h)-lastTime)*1e9 << endl; 
		}   
		lastTime = group->Time(h); 
		lastWire = group->WireNum(h); 
	      }	  
	      
	      return v;}, {"L_vdc_u1_wire","L_vdc_u1_rawtime"})
    
    .Histo1D({"", "T_{i} - T_{i+1} (ns)", 200, -100, 100}, "hit_sep"); 
  
  
  h1d_test->DrawCopy(); 
  return; //*/ 
  
 
  

  
  /*auto h2d_test = nEvents_coinc 
    
    .Define("x_err", []( TvdcTrack *trk_good, 
			 TvdcTrack *trk_bad ) { 
	      	
	      return trk_bad->S2_y(); 
	    }, {"track","track_accidental"})
    
    .Define("T_err", []( TvdcTrack *trk_good, 
			 TvdcTrack *trk_bad ) { 
	      	
	      return trk_bad->Phi(); 
	    }, {"track","track_accidental"})
    
    .Histo2D({"", ";Y_{S2} (m);#Phi (rad)", 
	  200, -0.3, 0.3, 200, -0.15, 0.15}, "x_err","T_err"); 
  
	  h2d_test->DrawCopy("col"); */ 
  
  
	      
	      
  
  
  auto stopwatch = new TStopwatch; 
  
  auto nTracks  = nEvents_coinc.Count(); 
  auto snapshot = nEvents_coinc.Snapshot( "T", path_outFile.Data(), out_branches ); 
  
  int nTracks_count = *nTracks; 
  
  double net_time = stopwatch->RealTime(); 
  
  cout << TString::Format(" total time = %5.5f s (%5.5f ms/trk)", 
			  net_time, 1e3*net_time/((double)nTracks_count) )
       << endl; 
  
  cout << TString::Format("Total tracks reconstructed: %i", nTracks_count) << endl; 
  
  
  
  /*auto h2d_Test = nEvents_coinc 
    
    .Define("w", [](ROOT::VecOps::RVec<THitGroup*> groups,
		    TvdcTrack *trk) { 
	      
	      RVecD w; 
	      
	      for (int p=0; p<4; p++) { 
		
		auto group = groups.at(p); 
		
		for (int h=0; h<group->Nhits(); h++) { 
		  
		  w.push_back( trk->Slope(p)*( group->WirePos(h) - 
					       trk->Intercept(p) ) ); 
		}
	      }
	      
	      return w; 
	    }, {"groups","track"})
        
    .Histo1D({"test",";Monte-carlo drift-dist (m);", 
	  200, -25e-3, 25e-3}, "w"); 
  
  
	  h2d_Test->DrawCopy(); */     
}

