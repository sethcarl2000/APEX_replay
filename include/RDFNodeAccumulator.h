#ifndef RDFNodeAccumulator_h_
#define RDFNodeAccumulator_h_

#include <vector>
#include <string> 
#include <stdexcept> 
#include <sstream>
#include <cstdlib> 
#include <iostream> 
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>

//A small helper class which is meant to avoid some of the awkward syntax assocaited with RDataFrame creation. 
class RDFNodeAccumulator {
public: 
    enum EStatus {
        kGood    =0,
        kWarning =1, 
        kError =-1
    }; 

private: 
    ROOT::RDF::RNode fNode;     
    std::ostringstream fErrorMsg; 
    EStatus fStatus{kGood};
    bool fAbortOnError{true}; 


    // list of branches to include in an output file. 
    std::vector<std::string> fOutputBranches; 

public: 
    RDFNodeAccumulator(ROOT::RDF::RNode start); 

    //constructor that seeks to construct RDataFrame object (without one needing to be manually constructed and passed)
    RDFNodeAccumulator(const char* tree_name, const char* path_infile); 

    ~RDFNodeAccumulator(); 

    //define a new branch, with lambda func. 'expression', and inputs 'inputs' 
    template<typename F> void Define(const char* new_branch, F expression, const std::vector<std::string>& inputs); 

    //overwrite a previously defined branch, or define it for the first time if it doesn't exist 
    template<typename F> void Overwrite(const char* new_branch, F expression, const std::vector<std::string>& inputs); 

    //same as above, but only define if this column is not yet defined in the dataframe
    template<typename F> void DefineIfMissing(const char* new_branch, F expression, const std::vector<std::string>& inputs);

    //apply a filter to this analysis chain; 
    template<typename F> void Filter(F expression, const std::vector<std::string>& inputs); 

    //check if a branch is defined in the dataframe as it currently exists      
    bool IsBranchDefined(std::string branch); 

    //return a count at the current last node
    ROOT::RDF::RResultPtr<ULong64_t> Count() { return fNode.Count(); }

    //bool IsBranchDefined(const char* branch) const { string str(branch); return IsBranchDefined(str); }

    //these do exactly the same as the functions above with the same names, but allow the input branch name to be a
    // std::string
    template<typename F> void Define(std::string new_branch, F expression, const std::vector<std::string>& inputs);
    
    template<typename F> void DefineIfMissing(std::string new_branch, F expression, const std::vector<std::string>& inputs);

    template<typename F> void Overwrite(std::string new_branch, F expression, const std::vector<std::string>& inputs);

    bool IsBranchDefined(const char* branch); 

    ROOT::RDF::RNode& Get() { return fNode; }

    EStatus GetStatus() const { return fStatus; }
    std::string GetErrorMsg() const { return fErrorMsg.str(); }

    //called when we want to print the stored error message and abort. 
    void PrintMsgAndAbort() {
        std::cerr << "<RDFNodeAccumulator>: Aborted execution. Error message: " << GetErrorMsg() << std::endl; 
        std::abort(); 
    }

    //add a pre-existing branch to the output
    void AddBranchToOutput(std::string branch) { fOutputBranches.push_back(branch); }

    //define a branch specifically for the output. 
    template<typename F> void DefineOutput(std::string new_branch, F expression, const std::vector<std::string>& inputs);
    
    void Snapshot(std::string tree_name, std::string path_outfile); 

    //set whether an abort should be called when an error is encountered. 
    void SetAbortOnError(bool _val) { fAbortOnError=_val; }

    //assignment operator 
    ROOT::RDF::RNode& operator=(const ROOT::RDF::RNode& node) {
        fNode = node; 
        return fNode; 
    }
    ClassDef(RDFNodeAccumulator,1); 
};

//now, we get actual class definitions. 

//__________________________________________________________________________________________________________________________________
RDFNodeAccumulator::RDFNodeAccumulator(ROOT::RDF::RNode start) : fNode{ start } {/* noop */}

//__________________________________________________________________________________________________________________________________
RDFNodeAccumulator::RDFNodeAccumulator(const char* tree_name, const char* path_infile)
    : fNode{ROOT::RDataFrame(tree_name, path_infile)}
{
    /*/try to construct RDataFrame
    try {

        ROOT::RDataFrame df(tree_name, path_infile); 
        fNode = df;

    } catch (const std::exception& e) {

        fErrorMsg << "in <RDFNodeAccumulator(constructor)>: Something went wrong trying to construct RDataFrame.\n"
            " -- File:   '" << path_infile << "'\n"
            " -- Tree:   '" << tree_name << "'\n"
            " -- what(): " << e.what();
         
        PrintMsgAndAbort(); 
        return; 
    }*/ 
}

//__________________________________________________________________________________________________________________________________
RDFNodeAccumulator::~RDFNodeAccumulator() {/* noop */}

//__________________________________________________________________________________________________________________________________
template<typename F> void RDFNodeAccumulator::Define(const char* new_branch, F expression, const std::vector<std::string>& inputs) 
{
    //check if the status is ok. if not, then abort. 
    if (GetStatus() != kGood) {
        //fErrorMsg << "\n - <Define>: branch '" << new_branch << "' cannot be added, status is not 'kGood'";
        if (fAbortOnError) { PrintMsgAndAbort(); } else { return; }
    }

    //check if this branch is already defined. if so, error! 
    if (IsBranchDefined(new_branch)) {
        //Error("RDFNodeAccumulator::Define", "exception caught while trying to define new branch '%s'.\n -- what(): %s", new_branch, e.what());
        fErrorMsg << "\n - <Define>: branch '" << new_branch << "' is already defined."
                     "\n             use 'DefineIfMissing' to skip definition if already defined, or 'Overwrite' to overwrite preexisting definition.";
        fStatus = kError; 
        if (fAbortOnError) { PrintMsgAndAbort(); } else { return; }
    }

    try {        
        fNode = fNode.Define(new_branch, expression, inputs); 
    
    } catch (const std::exception& e) {
        //Error("RDFNodeAccumulator::Define", "exception caught while trying to define new branch '%s'.\n -- what(): %s", new_branch, e.what());
        fErrorMsg << "\n - <Define>: exception caught defining branch '" << new_branch << "'. what(): " << e.what(); 
        fStatus = kError; 
        if (fAbortOnError) { PrintMsgAndAbort(); } else { return; }
    }
}

//__________________________________________________________________________________________________________________________________
template<typename F> void RDFNodeAccumulator::Overwrite(const char* new_branch, F expression, const std::vector<std::string>& inputs) {
        
    if (GetStatus() != kGood) {
        fErrorMsg << "\n - <Overwrite>: branch '" << new_branch << "' cannot be added, status is not ok.";
        if (fAbortOnError) { PrintMsgAndAbort(); } else { return; }
    }

    try {
        if (IsBranchDefined(new_branch)) { 
            fNode = fNode.Redefine(new_branch, expression, inputs); 
        } else {
            fNode = fNode.Define(new_branch, expression, inputs); 
        }
        
    } catch (const std::exception& e) {
        //Error("RDFNodeAccumulator::Define", "exception caught while trying to define new branch '%s'.\n -- what(): %s", new_branch, e.what());
        fErrorMsg << "\n - <Overwrite>: exception caught defining branch '" << new_branch << "'. "
                     "\n                what(): " << e.what(); 
        fStatus = kError; 
        if (fAbortOnError) { PrintMsgAndAbort(); } else { return; }
    }
}
//__________________________________________________________________________________________________________________________________
template<typename F> void RDFNodeAccumulator::DefineIfMissing(const char* new_branch, F expression, const std::vector<std::string>& inputs) 
{
    if (!IsBranchDefined(new_branch)) Define(new_branch, expression, inputs); 
}
//__________________________________________________________________________________________________________________________________
template<typename F> void RDFNodeAccumulator::Filter(F expression, const std::vector<std::string>& inputs)
{
    fNode = fNode.Filter(expression, inputs); 
} 
//__________________________________________________________________________________________________________________________________
bool RDFNodeAccumulator::IsBranchDefined(std::string branch) 
{
    for (const std::string& column : fNode.GetColumnNames()) { if (branch == column) return true; }
    return false; 
}
//__________________________________________________________________________________________________________________________________

//these functions are the same as above, but with a different signature

//__________________________________________________________________________________________________________________________________
template<typename F> void RDFNodeAccumulator::Define(std::string new_branch, F expression, const std::vector<std::string>& inputs)
{
    Define(new_branch.c_str(), expression, inputs); 
} 
//__________________________________________________________________________________________________________________________________
template<typename F> void RDFNodeAccumulator::Overwrite(std::string new_branch, F expression, const std::vector<std::string>& inputs) 
{
    Overwrite(new_branch.c_str(), expression, inputs); 
}
//__________________________________________________________________________________________________________________________________
template<typename F> void RDFNodeAccumulator::DefineIfMissing(std::string new_branch, F expression, const std::vector<std::string>& inputs) 
{
    DefineIfMissing(new_branch.c_str(), expression, inputs); 
}
//__________________________________________________________________________________________________________________________________
bool RDFNodeAccumulator::IsBranchDefined(const char* branch)
{
    std::string branch_str(branch); 
    return IsBranchDefined(branch_str); 
}
//__________________________________________________________________________________________________________________________________
template<typename F> void RDFNodeAccumulator::DefineOutput(std::string new_branch, F expression, const std::vector<std::string>& inputs)
{
    fOutputBranches.push_back(new_branch); 
    Define(new_branch, expression, inputs);  
}
//__________________________________________________________________________________________________________________________________
void RDFNodeAccumulator::Snapshot(std::string tree_name, std::string path_outfile)
{
    if (fOutputBranches.empty()) {
        fErrorMsg << "\n - <Snapshot>: no output branches! you may define them using 'DefineOutput' or 'AddBranchToOutput'.\n";
        fStatus = kError; 
        if (fAbortOnError) { PrintMsgAndAbort(); }
        return; 
    }

    //first, let's scan to be sure that we have all the branches we need. 
    std::vector<std::string> missing_branches{}; 
    for (const auto& out_br : fOutputBranches) {
        if (!IsBranchDefined(out_br.c_str())) missing_branches.push_back(out_br); 
    }

    if (!missing_branches.empty()) {
        fErrorMsg << 
        "\n"
        "- <Snapshot>: tried to write to output file, but the following branches are missing:\n"; 
        for (auto& br : missing_branches) fErrorMsg << "               " << br << "\n";
        
        fStatus = kError; 
        if (fAbortOnError) { PrintMsgAndAbort(); }
        return; 
    }

    try {
        fNode.Snapshot(tree_name, path_outfile, fOutputBranches); 
    } catch (const std::exception& e) {
        fErrorMsg << 
        "\n"
        "- <Snapshot>: something unexpected went wrong trying to make a snapshot.\n" 
        "              what(): "<< e.what() <<"\n";
        
        fStatus = kError; 
        if (fAbortOnError) { PrintMsgAndAbort(); }
        return; 
    }

}
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________________

ClassImp(RDFNodeAccumulator); 

#endif 