#ifndef AnalyzeAllData_h
#define AnalyzeAllData_h

//APEX headers
#include "MetadataFetcher.h"
#include "QuietErrorHandler.h"
#include "replay_paths.h"
//root headers
#include <TH2D.h>
#include <TH1D.h> 
#include <ROOT/RDataFrame.hxx> 
#include <ROOT/RDF/HistoModels.hxx>
//stdlib headers
#include <vector>
#include <memory> 
#include <map> 
#include <string> 
#include <functional> 

using RDataframeUpdateFcn = std::function<ROOT::RDF::RNode(ROOT::RDF::RNode)>; 

class AnalyzeAllData {
 private: 

  //verbosity of script 
  int fVerbose{0};
  
  int fNThreads{1};

  int fMinRun, fMaxRun; 
  
  std::vector<replay_paths::replay_segment> fSegmentList; 
  
  //add contents of 'source' into 'target' 
  void StackHistograms(TH2D* target, TH2D* source); 
  void StackHistograms(TH1D* target, TH1D* source); 

  bool fAddMetadata{false}; 
  std::unique_ptr<MetadataFetcher> fMetadata{nullptr}; 
  
  static constexpr char kClassName[] = "AnalyzeAllData"; 

  ROOT::RDF::RNode Add_metadata_to_node(ROOT::RDF::RNode node) const; 
  
  std::string fDataset_identifier; 

 public: 

  /// @brief construct AnalyzeAllData class.
  /// @param n_threads number of threads to use. '0' means use all available threads, '1' for single-threadded execution. 
  /// @param verbose verbosity level: 0 -> no output; 1 -> [default] print msg at start and end of analysis loop. 2 -> print status for each file 
  /// @param dataset_identifier identifier for dataset to use. see 'replay_paths.h' for valid options.  
  AnalyzeAllData(int n_threads, int verbose=1, std::string dataset_identifier="ifarm-all");
  
  std::unique_ptr<MetadataFetcher> MakeMetadataFetcher(std::string branch) const; 
  
  void SetVerbosity(int v) { fVerbose=v; }

  //call this before using meta-data branches
  void AddMetadata(); 
  
  /// @brief given a TH1D constructor, and a function that adds/redefines RDF branches, the histogram will be filled.
  /// @param h_params list of parameters to initialize the TH2D with. 
  /// @param branch_x branch to fill on x-axis of histogram
  /// @param branch_y branch to fill on y-axis of histogram
  /// @param fcn function that defines x and y branches. if none is provided, then no new branches are added to the df (they must already be present!)
  /// @param target_tree TTree to examine. can be either 'track_data' or 'meta_data'. 
  TH2D* Make_TH2D(const ROOT::RDF::TH2DModel& h_params, std::string branch_x, std::string branch_y, const RDataframeUpdateFcn *fcn=nullptr, std::string target_tree="track_data"); 

  
  /// @brief given a TH1D constructor, and a function that adds/redefines RDF branches, a filled TH1D will be returned 
  /// @param h_params list of parameters to initialize the TH1D with. 
  /// @param branch_x branch to fill on x-axis of histogram
  /// @param branch_y branch to fill on y-axis of histogram
  /// @param fcn function that defines x and y branches. if none is provided, then no new branches are added to the df (they must already be present!)
  /// @param target_tree TTree to examine. can be either 'track_data' or 'meta_data'. 
  TH1D* Make_TH1D(const ROOT::RDF::TH1DModel& h_params, std::string branch_x, const RDataframeUpdateFcn *fcn=nullptr, std::string target_tree="track_data"); 
  
}; 


#endif
