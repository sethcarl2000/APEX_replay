#ifndef AnalyzeAllData_cxx
#define AnalyzeAllData_cxx

#include "AnalyzeAllData.h"
#include "replay_paths.h" 
#include "SlurmInterface.h" 

#include <ROOT/RResultPtr.hxx> 
#include <TString.h>
#include <TAxis.h>
#include <TError.h> 
#include <DllImport.h>
#include <TFile.h>
#include <TTree.h> 

#include <stdexcept>
#include <cstdio>
#include <iostream> 
#include <algorithm> 
#include <thread> 
#include <mutex> 

#include <RtypesCore.h> 

namespace {
  const std::string kTestPath = "/volatile/halla/apex/full_replay/production/replay-4175-4199.root";
};
//__________________________________________________________________________________________
AnalyzeAllData::AnalyzeAllData(int n_threads, int verbose, bool test_only)
  : fNThreads{n_threads},
    fVerbose{verbose},
    fMinRun{replay_paths::min_run}, fMaxRun{replay_paths::max_run}
{
  //
  if (test_only) {
    fMinRun = 4175;
    fMaxRun = 4199; 
    fPathList.push_back( kTestPath );
  } else {
    fPathList = replay_paths::list;
  }
  
  //silence errors, unless they are fatal. 
  auto& err_handler = QuietErrorHandler::Instance();
  err_handler.SetMinPrintLevel(kBreak); 
}
//__________________________________________________________________________________________
ROOT::RDF::RNode AnalyzeAllData::Add_metadata_to_node(ROOT::RDF::RNode node) const
{
  if (!fAddMetadata) return node;  

  if (!fMetadata) {
    throw std::logic_error("in <AnalyzeAllData::Add_metadata_to_node>:"
			   " fMetadata ptr is null!"); 
    return node; 
  }
  
  return node
    
    .Define("metadata", [this](int run, int segment, int rawfile)
    {
      return fMetadata->Get(run,segment,rawfile); 
    }, {"run_number","segment_number","rawfile_number"})

    //if this value is nan, then there is no metadata for this event
    .Filter([](const metadata_branch_t& md)
    {
      return md.dt_sigma == md.dt_sigma;
    }, {"metadata"})
    
    .Define("S2R_S2L_dt_sigma", [](const metadata_branch_t& md){ return md.dt_sigma; },
	    {"metadata"})
	    
    .Define("S2R_S2L_dt_center", [](const metadata_branch_t& md){ return md.dt_center; },
	    {"metadata"});    
}
//__________________________________________________________________________________________
void AnalyzeAllData::AddMetadata()
{
  fAddMetadata=true;

  fMetadata = std::move( MakeMetadataFetcher("null") );  
}
//__________________________________________________________________________________________
TH2D* AnalyzeAllData::Make_TH2D(const ROOT::RDF::TH2DModel& hmod, std::string branch_x, std::string branch_y, const RDataframeUpdateFcn *fcn, std::string target_tree)
{  
  //make the histogram
  auto hist = new TH2D(hmod.fName,   hmod.fTitle,
		       hmod.fNbinsX, hmod.fXLow, hmod.fXUp,
		       hmod.fNbinsY, hmod.fYLow, hmod.fYUp
		       );
  
  hist->SetDirectory(0); 

  auto xax = hist->GetXaxis(); 
  auto yax = hist->GetYaxis(); 

  size_t n_threads = fNThreads;

  //auto-determine the number of threads. note that this still works outside a slurm-job
  if (n_threads==1) {
    if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT();
    if (fVerbose>=1) 
      std::printf("<%s::%s>: running in single-thread mode.\n",kClassName,__func__); 
  } else {
    //auto-detect the number of threads we should use 
    if (n_threads==0) {
      auto& slurm_interface = SlurmInterface::Instance();
      n_threads = slurm_interface.get_n_threads(); 
    }
    ROOT::EnableImplicitMT(n_threads);   
    if (fVerbose>=1) 
      std::printf("<%s::%s>: running with %zi threads.\n",kClassName,__func__, n_threads); 
  }
  
  //now, we will decide how much work to give each thread.
  //to do this, we will first find how many events each thread has.
  ULong64_t total_events=0;
  for (const auto& segment : replay_paths::segments) {
    if (segment.is_good) total_events += segment.n_events; 
  }

  ULong64_t n_events_per_thread = total_events / ((ULong64_t)n_threads);

  std::mutex stack_mutex; 
  
  std::vector<std::thread> threads;
  threads.reserve(n_threads); 
  
  for (int t=0; t<n_threads; t++) {

    std::string thread_name{ Form("<thread [%2i/%zi]>",t+1,n_threads) }; 
    
    size_t start=t; 
    
    ULong64_t n_events_thread = 0; 
    
    std::vector<replay_paths::replay_segment> thread_segments{}; 

    const auto& all_segments = replay_paths::segments; 

    size_t end=start; 
    while (end < replay_paths::segments.size()) {

      const auto& seg = all_segments[end]; 
      end += n_threads; 
      
      if (!seg.is_good) continue;

      n_events_thread += seg.n_events;

      thread_segments.push_back( seg ); 
    }
    
    if (fVerbose>=2) {
      std::printf("%s: to-do list (%llu events total):\n",
		  thread_name.c_str(), n_events_thread);
      
      for (const auto& seg : thread_segments) {
	std::printf("  runs %i - %i (%llu events)\n",
		    seg.min_run, seg.max_run, seg.n_events); 
      }
    }
    
    //now, we can launch the threads.
    threads.emplace_back([&stack_mutex, t, thread_name, thread_segments, hist,
			 hmod, branch_x, branch_y, fcn, target_tree, this]
    {
      
      TH2D *hist_thread = new TH2D(Form("hthread_t%i",t), hmod.fTitle,
				   hmod.fNbinsX, hmod.fXLow, hmod.fXUp,
				   hmod.fNbinsY, hmod.fYLow, hmod.fYUp
				   );
      hist_thread->SetDirectory(0); 
      
      if (fVerbose>=2) std::cout << thread_name << ": starting loop\n"; 

      size_t i_task=1;
      for (const auto& seg : thread_segments) {
	
	if (fVerbose>=2) std::printf("%s: processing task [%2zi/%zi]; runs %i - %i\n",
				     thread_name.c_str(), i_task++,thread_segments.size(),
				     seg.min_run, seg.max_run); 
	
	TH2D* sub_hist=nullptr; 
	
	try {      

	  ROOT::RDF::RNode out_node = ROOT::RDataFrame(target_tree, seg.path);      
	  
	  //there are no new branches to add 
	  if (fcn == nullptr) {
      
	    sub_hist = (TH2D*)out_node.Histo2D({Form("h_%i",t), "",
		hmod.fNbinsX, hmod.fXLow, hmod.fXUp,
		hmod.fNbinsY, hmod.fYLow, hmod.fYUp
	      }, branch_x, branch_y)->Clone(Form("h_clone_%i",t));
	    
	  } else {
	
	    //add new branches
	    out_node = (*fcn)(out_node); 
	    sub_hist = (TH2D*)out_node.Histo2D({Form("h_%i",t), "",
		hmod.fNbinsX, hmod.fXLow, hmod.fXUp,
		hmod.fNbinsY, hmod.fYLow, hmod.fYUp
	      }, branch_x, branch_y)->Clone(Form("h_clone_%i",t));
	  }

	  sub_hist->SetDirectory(0); 

	} catch (const std::exception& e) {

	  if (fVerbose>=2)
	    std::printf("%s: error encountered; file skipped.\n"
			" ~~ what(): %s\n", thread_name.c_str(), e.what());  
	  continue; 
	}
    
	    
	if (sub_hist==nullptr) {
      
	  if (fVerbose>=2) std::printf("sub-hist is null; file skipped.\n");
	  continue; 
	}
	//tell ROOT we want to manage the memory for this ourselves
	//sub_hist->SetDirectory(0);

	//stack histograms
	if (fVerbose>=2) {
	  std::printf("sub-hist has %.0f entries. adding to main result...\n", sub_hist->Integral());
	}
	
	StackHistograms(hist_thread, sub_hist); 
	
	delete sub_hist;
	
	
      }//for (const auto& segment : thread_segments)

      stack_mutex.lock();
      if (fVerbose>=2) { std::cout << thread_name << ": done. stacking final result.\n"; }
      StackHistograms(hist, hist_thread);
      stack_mutex.unlock();

      delete hist_thread; 
    }); 
    
    
  }//for (int t=0; t<n_threads; t++) 

  for (auto& thread : threads) thread.join();  
  
  return hist; 
  /*  
  if (fNThreads != 1) {
    std::cout << "In <"<<kClassName<<"::"<<__func__<<">: using " << n_threads << " threads to process files.\n";  
  } else {
    if (fVerbose>=1)
      std::cout << "In <"<<kClassName<<"::"<<__func__<<">: executing in single-thread mode\n";
  }
  

  if (fVerbose>=1) {
    std::printf("in<%s::%s>: starting loop over all %zi files...\n", kClassName,__func__, fPathList.size());
  }
  for (size_t i=0; i<fPathList.size(); i++) {    
    
    const auto& path = fPathList[i];
    
    if (fVerbose>=2){
      std::printf(" ~~ processing file: %2zi/%zi '%s'...\n"
		  " ~~ ", i+1,fPathList.size(), path.c_str());  
      std::cout << std::flush; 
    }

    TH2D* sub_hist = nullptr; 

    try {      

      ROOT::RDF::RNode out_node = ROOT::RDataFrame(target_tree, path);

      if (fAddMetadata) 
	out_node = Add_metadata_to_node(out_node); 
      
      
      //there are no new branches to add 
      if (fcn == nullptr) {
      
	sub_hist = (TH2D*)out_node.Histo2D({Form("h_%zi",i), "",
	    xax->GetNbins(), xax->GetXmin(), xax->GetXmax(), 
	    yax->GetNbins(), yax->GetXmin(), yax->GetXmax()
	  }, branch_x, branch_y)->Clone(Form("h_clone_%zi",i));

      } else {
	
	//add new branches
	out_node = (*fcn)(out_node); 
	sub_hist = (TH2D*)out_node.Histo2D({Form("h_%zi",i), "",
	    xax->GetNbins(), xax->GetXmin(), xax->GetXmax(), 
	    yax->GetNbins(), yax->GetXmin(), yax->GetXmax()
	  }, branch_x, branch_y)->Clone(Form("h_clone_%zi",i));
      }
      sub_hist->SetDirectory(0); 

    } catch (const std::exception& e) {

      if (fVerbose>=2)
	std::printf("error encountered; file skipped.\n"
		    " ~~ what(): %s\n", e.what());  
      continue; 
    }
    
    
    if (sub_hist==nullptr) {
      
      if (fVerbose>=2) std::printf("sub-hist is null; file skipped.\n");
      continue; 
    }
    //tell ROOT we want to manage the memory for this ourselves
    //sub_hist->SetDirectory(0);

    //stack histograms
    if (fVerbose>=2) {
      std::printf("sub-hist has %.0f entries. adding to main result...", sub_hist->Integral());
      std::cout << std::flush; 
    }
    
    StackHistograms(hist, sub_hist); 

    delete sub_hist;
    
    if (fVerbose>=2) std::cout << "done." << std::endl; 
  }
  */ 
  return hist; 
}
//__________________________________________________________________________________________
TH1D* AnalyzeAllData::Make_TH1D(const ROOT::RDF::TH1DModel& hmod, std::string branch_x, const RDataframeUpdateFcn *fcn, std::string target_tree)
{

  //make metadata branches  
  
  //make the histogram
  auto hist = new TH1D(hmod.fName,   hmod.fTitle,
		       hmod.fNbinsX, hmod.fXLow, hmod.fXUp
		       );
  
  hist->SetDirectory(0); 
  
  //don't define any new branches
  auto xax = hist->GetXaxis(); 
  
  if (fNThreads != 1) {
    if (fVerbose>=1) std::cout << "In <"<<kClassName<<"::"<<__func__<<">: using " << fNThreads << " threads to process files.\n";  
    ROOT::EnableImplicitMT(fNThreads);
  } else {
    if (fVerbose>=1) std::cout << "In <"<<kClassName<<"::"<<__func__<<">: executing in single-thread mode\n";
    if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT();
  }
  
  
  if (fVerbose>=1) std::printf("in<%s::%s>: starting loop over all %zi files...\n", kClassName,__func__, fPathList.size()); 
  for (size_t i=0; i<fPathList.size(); i++) {    
    
    const auto& path = fPathList[i];

    if (fVerbose>=2) {
      std::printf(" ~~ processing file: %2zi/%zi '%s'...\n"
		  " ~~ ", i+1,fPathList.size(), path.c_str());  
      std::cout << std::flush; 
    }
    
    TH1D* sub_hist = nullptr; 

    try {

      ROOT::RDF::RNode out_node = ROOT::RDataFrame(target_tree, path);

      if (fAddMetadata) 
	out_node = Add_metadata_to_node(out_node); 

      //there are no new branches to add 
      if (fcn == nullptr) {
	
	sub_hist = (TH1D*)out_node.Histo1D({Form("h_%zi",i), "",
	    xax->GetNbins(), xax->GetXmin(), xax->GetXmax()
	  }, branch_x)->Clone(Form("h_clone_%zi",i));
	
      } else {
	
	//add new branches
	out_node = (*fcn)(out_node); 
	sub_hist = (TH1D*)out_node.Histo1D({Form("h_%zi",i), "",
	    xax->GetNbins(), xax->GetXmin(), xax->GetXmax()
	  }, branch_x)->Clone(Form("h_clone_%zi",i));

      }
      sub_hist->SetDirectory(0); 

    } catch (const std::exception& e) {

      if (fVerbose>=2) 
	std::printf("error encountered; file skipped.\n ~~ what(): %s\n", e.what());
      
      continue; 
    }
    
    
    if (sub_hist==nullptr) {
      if (fVerbose>=2) std::printf("sub-hist is null; file skipped.\n");
      continue; 
    }
    //tell ROOT we want to manage the memory for this ourselves
    //sub_hist->SetDirectory(0);

    //stack histograms
    if (fVerbose>=2) { 
      std::printf("sub-hist has %.0f entries. adding to main result...", sub_hist->Integral());
      std::cout << std::flush; 
    }
    StackHistograms(hist, sub_hist); 

    delete sub_hist;
    
    if (fVerbose>=2) std::cout << "done." << std::endl; 
  }

  return hist; 
}
//__________________________________________________________________________________________
void AnalyzeAllData::StackHistograms(TH2D* target, TH2D* source)
{
  if (target==nullptr) {
    Error(__func__, "target (output) histogram is null!");
    return;
  }

  if (source==nullptr) {
    Error(__func__, "source (input) histogram is null!");
    return;
  }

  auto xax = source->GetXaxis();
  auto yax = source->GetYaxis();

  for (int bx=1; bx<=xax->GetNbins(); bx++) {
    for (int by=1; by<=yax->GetNbins(); by++) {

      target->Fill(xax->GetBinCenter(bx),
		   yax->GetBinCenter(by),
		   source->GetBinContent(bx,by)
		   );
    }
  }
  
}

//__________________________________________________________________________________________
void AnalyzeAllData::StackHistograms(TH1D* target, TH1D* source)
{
  if (target==nullptr) {
    Error(__func__, "target (output) histogram is null!");
    return;
  }

  if (source==nullptr) {
    Error(__func__, "source (input) histogram is null!");
    return;
  }
  

  auto xax = source->GetXaxis();
  
  for (int bx=1; bx<=xax->GetNbins(); bx++) {
    
    target->Fill(xax->GetBinCenter(bx),
		 source->GetBinContent(bx)
		 );
  }
  
}
//__________________________________________________________________________________________
std::unique_ptr<MetadataFetcher> AnalyzeAllData::MakeMetadataFetcher(std::string branch) const 
{

  auto fetcher = std::make_unique<MetadataFetcher>(branch, fMinRun, fMaxRun); 

  if (fVerbose>=1) 
    std::printf("in<%s::%s>: starting loop over all %zi files...\n",
		kClassName,__func__, fPathList.size()); 
  
  //we need to do this in single-threading mode
  if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT(); 
  
  int ind=0; 
  for (size_t i=0; i<fPathList.size(); i++) {    

    const auto& path = fPathList[i]; 

    if (fVerbose>=2) { 
      std::printf(" ~~ processing file: %2zi/%zi '%s'...\n"
		  " ~~ ", i+1,fPathList.size(), path.c_str());  
      std::cout << std::flush; 
    }
    
    try {
            
      ROOT::RDataFrame df_meta("meta_data", path);
      
      auto md_subresult = *df_meta
	
	.Define("md", [&fetcher,&ind](int run, int segment, int rawfile,
				      double dt, double sigma)
	{
	  if (fetcher->fRunIndex[run - fetcher->fMinRun] < 0) {
	    fetcher->fRunIndex[run - fetcher->fMinRun] = ind;
	  }
	  
	  ++ind;
	  
	  return metadata_branch_t{ run, segment, rawfile, dt, sigma };

	}, {"run_number", "segment_number", "rawfile_number",
		 "S2R_S2L_dt_center", "S2R_S2L_dt_sigma"})
	
	.Take<metadata_branch_t>("md"); 


      if (md_subresult.empty()) {
	
	std::printf("no meta-data to add; file skipped.\n");
	continue; 
      }

      if (fVerbose >= 2) {
	std::printf("%zi results to add...", md_subresult.size());
	std::cout << std::flush;
      }
      
      //copy this result into the overall vector
      fetcher->fData.insert(std::end(fetcher->fData),
			    std::begin(md_subresult),
			    std::end(md_subresult)
			    );

      
      if (fVerbose >= 2) std::printf("done. overall size: %zi\n", fetcher->fData.size());
      
    
    } catch (const std::exception& e) {

      if (fVerbose >= 2) 
	std::printf("error encountered; file skipped.\n"
		    " ~~ what(): %s\n", e.what());  
      continue; 
    }
  }
  
  return fetcher; 
} 
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________
//__________________________________________________________________________________________


#endif
