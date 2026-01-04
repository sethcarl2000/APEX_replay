#ifndef Podd_TrdfHelper_h_
#define Podd_TrdfHelper_h_

//////////////////////////////////////////////////////////////////////////
//
// TrdfHelper
//
//////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <string>
#include <functional>

using RNode = ROOT::RDF::RNode; 

class TrdfHelper {
  
public:

  TrdfHelper(RNode *node, vector<string> branches={})
    : fNode(node), fBranches(branches) {}; 
  
  template <typename F> void Define(const TString out,
				    const *F,
				    const vector<TString> cols);

  template <typename T> void Redefine(const TString bold,
				      const vecotr<TString> bnew); 

  void Snapshot(const TString treeName, const TString path);
  
private:

  RNode *fNode;
  vector<string> fBranches; 
  
  ClassDef(TrdfHelper,0)   //Example physics module
};

#endif
