#ifndef AnalyzeAllData_cxx
#define AnalyzeAllData_cxx

#include "AnalyzeAllData.h"

#include <ROOT/RResultPtr.hxx> 
#include <TString.h>
#include <TAxis.h>
#include <TError.h> 

#include <stdexcept>
#include <cstdio>
#include <iostream> 

  
//list of all APEX full-replay paths
const std::vector<std::string> AnalyzeAllData::fPathList = std::vector<std::string>{
  "/volatile/halla/apex/full_replay/production/replay-3800-3824.root",
  "/volatile/halla/apex/full_replay/production/replay-3825-3849.root",
  "/volatile/halla/apex/full_replay/production/replay-3850-3874.root",
  "/volatile/halla/apex/full_replay/production/replay-3875-3899.root",
  "/volatile/halla/apex/full_replay/production/replay-3900-3924.root",
  "/volatile/halla/apex/full_replay/production/replay-3925-3949.root",
  "/volatile/halla/apex/full_replay/production/replay-3950-3974.root",
  "/volatile/halla/apex/full_replay/production/replay-3975-3999.root",
  "/volatile/halla/apex/full_replay/production/replay-4000-4024.root",
  "/volatile/halla/apex/full_replay/production/replay-4025-4049.root",
  "/volatile/halla/apex/full_replay/production/replay-4050-4074.root",
  "/volatile/halla/apex/full_replay/production/replay-4075-4099.root",
  "/volatile/halla/apex/full_replay/production/replay-4100-4124.root",
  "/volatile/halla/apex/full_replay/production/replay-4125-4149.root",
  "/volatile/halla/apex/full_replay/production/replay-4150-4174.root",
  "/volatile/halla/apex/full_replay/production/replay-4175-4199.root",
  "/volatile/halla/apex/full_replay/production/replay-4200-4224.root",
  "/volatile/halla/apex/full_replay/production/replay-4225-4249.root",
  "/volatile/halla/apex/full_replay/production/replay-4250-4274.root",
  "/volatile/halla/apex/full_replay/production/replay-4275-4299.root",
  "/volatile/halla/apex/full_replay/production/replay-4300-4324.root",
  "/volatile/halla/apex/full_replay/production/replay-4325-4349.root",
  "/volatile/halla/apex/full_replay/production/replay-4350-4374.root",
  "/volatile/halla/apex/full_replay/production/replay-4375-4399.root",
  "/volatile/halla/apex/full_replay/production/replay-4400-4424.root",
  "/volatile/halla/apex/full_replay/production/replay-4425-4449.root",
  "/volatile/halla/apex/full_replay/production/replay-4450-4474.root",
  "/volatile/halla/apex/full_replay/production/replay-4475-4499.root",
  "/volatile/halla/apex/full_replay/production/replay-4500-4524.root",
  "/volatile/halla/apex/full_replay/production/replay-4525-4549.root",
  "/volatile/halla/apex/full_replay/production/replay-4550-4574.root",
  "/volatile/halla/apex/full_replay/production/replay-4575-4599.root",
  "/volatile/halla/apex/full_replay/production/replay-4600-4624.root",
  "/volatile/halla/apex/full_replay/production/replay-4625-4649.root",
  "/volatile/halla/apex/full_replay/production/replay-4650-4674.root",
  "/volatile/halla/apex/full_replay/production/replay-4675-4699.root",
  "/volatile/halla/apex/full_replay/production/replay-4700-4724.root",
  "/volatile/halla/apex/full_replay/production/replay-4725-4749.root",
  "/volatile/halla/apex/full_replay/production/replay-4750-4774.root",
  "/volatile/halla/apex/full_replay/production/replay-4775-4799.root",
  "/volatile/halla/apex/full_replay/production/replay-4800-4824.root",
  "/volatile/halla/apex/full_replay/production/replay-4825-4849.root",
  "/volatile/halla/apex/full_replay/production/replay-4850-4874.root",
  "/volatile/halla/apex/full_replay/production/replay-4875-4899.root",
  "/volatile/halla/apex/full_replay/production/replay-4900-4924.root",
  "/volatile/halla/apex/full_replay/production/replay-4925-4949.root",
  "/volatile/halla/apex/full_replay/production/replay-4950-4974.root",
  "/volatile/halla/apex/full_replay/production/replay-4975-4999.root"
};


//______________________________________________________________________________________________________
void AnalyzeAllData::Fill_TH2D(TH2D* hist, std::string branch_x, std::string branch_y, const std::function<ROOT::RDF::RNode(ROOT::RDF::RNode)> *fcn, std::string target_tree)
{
  
  
  if (!hist) {
    Error(__func__, "hist ptr passed is null!");
    return;
  }
  
  
  if (fNThreads != 1) {
    std::cout << "In <"<<kClassName<<"::"<<__func__<<">: using " << fNThreads << " threads to process files.\n";  
    ROOT::EnableImplicitMT(fNThreads);
  } else {
    std::cout << "In <"<<kClassName<<"::"<<__func__<<">: executing in single-thread mode\n";
    if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT();
  }
  
  
  std::printf("in<%s::%s>: starting loop over all %zi files...\n", kClassName,__func__, fPathList.size()); 
  for (size_t i=0; i<fPathList.size(); i++) {
     
    
    const auto& path = fPathList[i];

    std::printf(" ~~ processing file: %2zi/%zi '%s'...\n"
		" ~~ ", i+1,fPathList.size(), path.c_str());  
    std::cout << std::flush; 
    
    ROOT::RDataFrame df(target_tree, path);

    ROOT::RDF::RResultPtr<TH2D> hist_rptr; 
    
    //don't define any new branches
    auto xax = hist->GetXaxis(); 
    auto yax = hist->GetYaxis(); 

    TH2D* sub_hist = nullptr; 

    try {      
      //there are no new branches to add 
      if (fcn == nullptr) {
      
	sub_hist = (TH2D*)df.Histo2D<double>({Form("h_%zi",i), "",
	    xax->GetNbins(), xax->GetXmin(), xax->GetXmax(), 
	    yax->GetNbins(), yax->GetXmin(), yax->GetXmax()
	  }, branch_x, branch_y)->Clone(Form("h_clone_%zi",i));

      } else {
      
	//add new branches
	auto df_out = (*fcn)(df); 
	sub_hist = (TH2D*)df_out.Histo2D<double>({Form("h_%zi",i), "",
	    xax->GetNbins(), xax->GetXmin(), xax->GetXmax(), 
	    yax->GetNbins(), yax->GetXmin(), yax->GetXmax()
	  }, branch_x, branch_y)->Clone(Form("h_clone_%zi",i));
  
      }

    } catch (const std::exception& e) {

      std::printf("error encountered; file skipped.\n"
		  " ~~ what(): %s\n", e.what());  
      continue; 
    }
    
    
    if (sub_hist==nullptr) {
      std::printf("sub-hist is null; file skipped.\n");
      continue; 
    }
    //tell ROOT we want to manage the memory for this ourselves
    //sub_hist->SetDirectory(0);

    //stack histograms

    std::printf("sub-hist has %.0f entries. adding to main result...", sub_hist->Integral());
    std::cout << std::flush; 

    StackHistograms(hist, sub_hist); 

    //delete sub_hist; 
    
    std::cout << "done." << std::endl; 
  }
  
}
//_______________________________________________________________________________________________________
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

      target->Fill(
		   xax->GetBinCenter(bx),
		   yax->GetBinCenter(by),
		   source->GetBinContent(bx,by)
		   );
    }
  }
  
}
//_______________________________________________________________________________________________________
//_______________________________________________________________________________________________________
//_______________________________________________________________________________________________________
//_______________________________________________________________________________________________________
//_______________________________________________________________________________________________________
//_______________________________________________________________________________________________________
//_______________________________________________________________________________________________________



#endif
