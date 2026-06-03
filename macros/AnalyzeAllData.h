#ifndef AnalyzeAllData_h
#define AnalyzeAllData_h

//APEX headers
#include "MetadataFetcher.h"
//root headers
#include <TH2D.h>
#include <ROOT/RDataFrame.hxx> 
#include <ROOT/RDF/HistoModels.hxx>
//stdlib headers
#include <vector>
#include <map> 
#include <string> 
#include <functional> 

using RDataframeUpdateFcn = std::function<ROOT::RDF::RNode(ROOT::RDF::RNode)>; 

class AnalyzeAllData {
 private: 

  //verbosity of script 
  int fVerbose{0};
  
  int fNThreads{1};

  static const std::vector<std::string> fPathList; 
  
  //add contents of 'source' into 'target' 
  void StackHistograms(TH2D* target, TH2D* source); 
  
  static constexpr char kClassName[] = "AnalyzeAllData"; 

 public: 

  /// @brief construct AnalyzeAllData class.
  /// @param n_threads number of threads to use. '0' means use all available threads. 
  AnalyzeAllData(int n_threads=1) : fNThreads{n_threads} {}; 
  
  MetadataFetcher MakeMetadataFetcher(std::string branch) const; 
  
  void SetVerbosity(int v) { fVerbose=v; }
  
  /// @brief given a TH2 ptr, and a function that takes an RDF as input and reutnrs a lazy TH2 result ptr, the histogram will be filled. 
  /// @param h_params list of parameters to initialize the TH2 with. 
  /// @param branch_x branch to fill on x-axis of histogram
  /// @param branch_y branch to fill on y-axis of histogram
  /// @param fcn function that defines x and y branches. if none is provided, then no new branches are added to the df (they must already be present!)
  /// @param target_tree TTree to examine. can be either 'track_data' or 'meta_data'. 
  TH2D* Make_TH2D(const ROOT::RDF::TH2DModel& h_params, std::string branch_x, std::string branch_y, const RDataframeUpdateFcn *fcn=nullptr, std::string target_tree="track_data"); 
  
}; 


#endif
