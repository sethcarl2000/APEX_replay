#ifndef AnalyzeAllData_cxx
#define AnalyzeAllData_cxx

#include "AnalyzeAllData.h"
#include "replay_paths.h" 

#include <ROOT/RResultPtr.hxx> 
#include <TString.h>
#include <TAxis.h>
#include <TError.h> 
#include <DllImport.h>

#include <stdexcept>
#include <cstdio>
#include <iostream> 
#include <algorithm> 


const std::vector<std::string> AnalyzeAllData::fPathList = replay_paths::list; 

//__________________________________________________________________________________________
//__________________________________________________________________________________________
void AnalyzeAllData::AddMetadata()
{
  fAddMetadata=true;

  fMetadata = std::move( MakeMetadataFetcher("null") );  
}
//__________________________________________________________________________________________
TH2D* AnalyzeAllData::Make_TH2D(const ROOT::RDF::TH2DModel& hmod, std::string branch_x, std::string branch_y, const RDataframeUpdateFcn *fcn, std::string target_tree)
{

  //make metadata branches
  
  
  //make the histogram
  auto hist = new TH2D(
		       hmod.fName,   hmod.fTitle,
		       hmod.fNbinsX, hmod.fXLow, hmod.fXUp,
		       hmod.fNbinsY, hmod.fYLow, hmod.fYUp
		       );
  
  hist->SetDirectory(0); 
  
  if (fNThreads != 1) {
    if (fVerbose>=1) 
      std::cout << "In <"<<kClassName<<"::"<<__func__<<">: using " << fNThreads << " threads to process files.\n";  
    ROOT::EnableImplicitMT(fNThreads);
  } else {
    if (fVerbose>=1)
      std::cout << "In <"<<kClassName<<"::"<<__func__<<">: executing in single-thread mode\n";
    if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT();
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

    ROOT::RDF::RResultPtr<TH2D> hist_rptr; 
    
    //don't define any new branches
    auto xax = hist->GetXaxis(); 
    auto yax = hist->GetYaxis(); 

    TH2D* sub_hist = nullptr; 

    try {      

      ROOT::RDataFrame df(target_tree, path);

      //there are no new branches to add 
      if (fcn == nullptr) {
      
	sub_hist = (TH2D*)df.Histo2D({Form("h_%zi",i), "",
	    xax->GetNbins(), xax->GetXmin(), xax->GetXmax(), 
	    yax->GetNbins(), yax->GetXmin(), yax->GetXmax()
	  }, branch_x, branch_y)->Clone(Form("h_clone_%zi",i));

	sub_hist->SetDirectory(0); 
	
      } else {
      
	//add new branches
	auto df_out = (*fcn)(df); 
	sub_hist = (TH2D*)df_out.Histo2D({Form("h_%zi",i), "",
	    xax->GetNbins(), xax->GetXmin(), xax->GetXmax(), 
	    yax->GetNbins(), yax->GetXmin(), yax->GetXmax()
	  }, branch_x, branch_y)->Clone(Form("h_clone_%zi",i));

	sub_hist->SetDirectory(0); 
      }

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

  return hist; 
}
//__________________________________________________________________________________________
TH1D* AnalyzeAllData::Make_TH1D(const ROOT::RDF::TH1DModel& hmod, std::string branch_x, const RDataframeUpdateFcn *fcn, std::string target_tree)
{

  //make metadata branches
  
  
  //make the histogram
  auto hist = new TH1D(
		       hmod.fName,   hmod.fTitle,
		       hmod.fNbinsX, hmod.fXLow, hmod.fXUp
		       );
  
  hist->SetDirectory(0); 
  
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
    
    //don't define any new branches
    auto xax = hist->GetXaxis(); 
    auto yax = hist->GetYaxis(); 

    TH1D* sub_hist = nullptr; 

    try {
      ROOT::RDataFrame df(target_tree, path);
      
      //there are no new branches to add 
      if (fcn == nullptr) {
	
	sub_hist = (TH1D*)df.Histo1D({Form("h_%zi",i), "",
	    xax->GetNbins(), xax->GetXmin(), xax->GetXmax()
	  }, branch_x)->Clone(Form("h_clone_%zi",i));

	sub_hist->SetDirectory(0); 
	
      } else {
      
	//add new branches
	auto df_out = (*fcn)(df); 
	sub_hist = (TH1D*)df_out.Histo1D({Form("h_%zi",i), "",
	    xax->GetNbins(), xax->GetXmin(), xax->GetXmax()
	  }, branch_x)->Clone(Form("h_clone_%zi",i));

	sub_hist->SetDirectory(0); 
      }

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

  auto fetcher = std::make_unique<MetadataFetcher>(branch, replay_paths::min_run, replay_paths::max_run); 

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
