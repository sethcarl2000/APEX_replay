#ifndef gen_pid_data_H
#define gen_pid_data_H

#include <TapexEventHandler.h>
#include "../include/RDFNodeAccumulator.h"
#include <TapexS2Hit.h>
#include <ROOT/RVec.hxx>
#include <math.h>
#include <string> 

struct CerenkovPMT_t { int paddle; double adc, tdc; };

namespace {

  constexpr int n_cerenkov_paddles = 10; 
  
}; 

void gen_pid_data(const bool is_RHRS, RDFNodeAccumulator& rna)
{
  const std::string arm = is_RHRS ? "R" : "L"; 
  
  using rvecd = ROOT::VecOps::RVec<double>; 

  
  rna.Define(arm+"_cerenkov_paddles", [](const rvecd& adc, const rvecd& tdc) {
      
      //add only the hits that have non-null TDC times
      ROOT::RVec<CerenkovPMT_t> paddles; //paddles.reserve(adc.size()); 

      for (int i=0; i<n_cerenkov_paddles; i++) {

	double t = tdc.at(i); 
	
	if ( std::fabs(t) > 1. ) continue;

	paddles.push_back({ i, adc.at(i), tdc.at(i) });
      }
      
      return paddles;
      
    }, {arm+".cer.a_c", arm+".cer.t_c"});

  rna.DefineOutput(arm+"_cer_paddle", [](const ROOT::RVec<CerenkovPMT_t>& hits) {
    
    ROOT::RVec<int> v; v.reserve(hits.size());
    for (const auto& hit : hits) { v.push_back(hit.paddle); }
    return v; 
  }, {arm+"_cerenkov_paddles"}); 
  
  
  rna.DefineOutput(arm+"_cer_time", [](const ROOT::RVec<CerenkovPMT_t>& hits) {
    
    rvecd v; v.reserve(hits.size());
    for (const auto& hit : hits) { v.push_back(hit.tdc); }
    return v; 
  }, {arm+"_cerenkov_paddles"}); 
  
  
  rna.DefineOutput(arm+"_cer_adc", [](const ROOT::RVec<CerenkovPMT_t>& hits) {
      
    rvecd v; v.reserve(hits.size());
    for (const auto& hit : hits) { v.push_back(hit.adc); }
    return v; 
  }, {arm+"_cerenkov_paddles"}); 


  //now, lets define the pre-shower & shower
  rna.DefineOutput(arm+"_ps_adc", [](const rvecd& v){ return v; },
		   {is_RHRS ? "R.ps.a_c" : "L.prl1.a_c"}); 
  
  rna.DefineOutput(arm+"_sh_adc", [](const rvecd& v){ return v; },
		   {is_RHRS ? "R.sh.a_c" : "L.prl2.a_c"}); 
  
}   



#endif
