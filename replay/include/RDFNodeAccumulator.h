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
#include <ROOT/RSnapshotOptions.hxx> 

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
    inline RDFNodeAccumulator(ROOT::RDF::RNode start) : fNode{ start } {/* noop */};

    //__________________________________________________________________________________________________________________________________
    //constructor that seeks to construct RDataFrame object (without one needing to be manually constructed and passed)
    inline RDFNodeAccumulator(const char* tree_name, const char* path_infile)
        : fNode{ROOT::RDataFrame(tree_name, path_infile)}
    {};

    inline ~RDFNodeAccumulator() {/* noop */}; 
    
    //__________________________________________________________________________________________________________________________________
    template<typename F> inline void Define(const char* new_branch, F expression, const std::vector<std::string>& inputs) 
    {
        //check if the status is ok. if not, then abort. 
        if (GetStatus() != kGood) {
            //fErrorMsg << "\n - <Define>: branch '" << new_branch << "' cannot be added, status is not 'kGood'";
            PrintMsgAndAbort();
            return;
        }

        //check if this branch is already defined. if so, error! 
        if (IsBranchDefined(new_branch)) {
            //Error("RDFNodeAccumulator::Define", "exception caught while trying to define new branch '%s'.\n -- what(): %s", new_branch, e.what());
            fErrorMsg << "\n - <Define>: branch '" << new_branch << "' is already defined."
                        "\n             use 'DefineIfMissing' to skip definition if already defined, or 'Overwrite' to overwrite preexisting definition.";
            fStatus = kError; 
            PrintMsgAndAbort();
            return;
        }

        try {        
            fNode = fNode.Define(new_branch, expression, inputs); 
        
        } catch (const std::exception& e) {
            //Error("RDFNodeAccumulator::Define", "exception caught while trying to define new branch '%s'.\n -- what(): %s", new_branch, e.what());
            fErrorMsg << "\n - <Define>: exception caught defining branch '" << new_branch << "'. what(): " << e.what(); 
            fStatus = kError; 
            PrintMsgAndAbort();
            return;
        }
    }

    //__________________________________________________________________________________________________________________________________
    //overwrite a previously defined branch, or define it for the first time if it doesn't exist 
    template<typename F> inline void Overwrite(const char* new_branch, F expression, const std::vector<std::string>& inputs) 
    {
            
        if (GetStatus() != kGood) {
            fErrorMsg << "\n - <Overwrite>: branch '" << new_branch << "' cannot be added, status is not ok.";
            PrintMsgAndAbort();
            return;
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
            PrintMsgAndAbort();
            return;
        }
    }

    //__________________________________________________________________________________________________________________________________
    //same as above, but only define if this column is not yet defined in the dataframe
    template<typename F> void DefineIfMissing(const char* new_branch, F expression, const std::vector<std::string>& inputs) 
    {
        if (!IsBranchDefined(new_branch)) Define(new_branch, expression, inputs); 
    }

    //__________________________________________________________________________________________________________________________________
    //apply a filter to this analysis chain; 
    template<typename F> inline void Filter(F expression, const std::vector<std::string>& inputs)
    {
        fNode = fNode.Filter(expression, inputs); 
    } 

    //__________________________________________________________________________________________________________________________________
    //check if a branch is defined in the dataframe as it currently exists      
    inline bool IsBranchDefined(std::string branch) 
    {
        for (const std::string& column : fNode.GetColumnNames()) { if (branch == column) return true; }
        return false; 
    }

    //__________________________________________________________________________________________________________________________________
    //return a count at the current last node
    inline ROOT::RDF::RResultPtr<ULong64_t> Count() { return fNode.Count(); }

    //__________________________________________________________________________________________________________________________________
    //these do exactly the same as the functions above with the same names, but allow the input branch name to be a
    // std::string
    template<typename F> inline void Define(std::string new_branch, F expression, const std::vector<std::string>& inputs)
    {
        Define(new_branch.c_str(), expression, inputs); 
    } 

    //__________________________________________________________________________________________________________________________________
    template<typename F> inline void DefineIfMissing(std::string new_branch, F expression, const std::vector<std::string>& inputs)
    {
        Overwrite(new_branch.c_str(), expression, inputs); 
    }

    //__________________________________________________________________________________________________________________________________
    template<typename F> inline void Overwrite(std::string new_branch, F expression, const std::vector<std::string>& inputs)
    {
        DefineIfMissing(new_branch.c_str(), expression, inputs); 
    }
    

    //__________________________________________________________________________________________________________________________________
    inline ROOT::RDF::RNode& Get() { return fNode; }

    //__________________________________________________________________________________________________________________________________
    inline EStatus GetStatus() const { return fStatus; }
    
    //__________________________________________________________________________________________________________________________________
    inline std::string GetErrorMsg() const { return fErrorMsg.str(); }

    //__________________________________________________________________________________________________________________________________
    //called when we want to print the stored error message and abort. 
    inline void PrintMsgAndAbort() {
        if (fAbortOnError) {
            std::cerr << "<RDFNodeAccumulator>: Aborted execution. Error message: " << GetErrorMsg() << std::endl; 
            std::abort(); 
        } else {
            throw std::logic_error(GetErrorMsg()); 
            return; 
        }
    }

    //__________________________________________________________________________________________________________________________________
    //add a pre-existing branch to the output
    inline void AddBranchToOutput(std::string branch) { fOutputBranches.push_back(branch); }

    //__________________________________________________________________________________________________________________________________
    //define a branch specifically for the output. 
    template<typename F> inline void DefineOutput(std::string new_branch, F expression, const std::vector<std::string>& inputs)
    {
        fOutputBranches.push_back(new_branch); 
        Define(new_branch, expression, inputs);  
    }
  
    //__________________________________________________________________________________________________________________________________
    inline void Snapshot(std::string tree_name, std::string path_outfile, ROOT::RDF::RSnapshotOptions* opts=nullptr)
    {
        if (fOutputBranches.empty()) {
            fErrorMsg << "\n - <Snapshot>: no output branches! you may define them using 'DefineOutput' or 'AddBranchToOutput'.\n";
            fStatus = kError; 
            PrintMsgAndAbort();
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
            PrintMsgAndAbort();
            return;
        }

        try {
        
        if (opts == nullptr) {
            fNode.Snapshot(tree_name, path_outfile, fOutputBranches);
        } else {
            fNode.Snapshot(tree_name, path_outfile, fOutputBranches, *opts);
        }
        } catch (const std::exception& e) {
            fErrorMsg << 
            "\n"
            "- <Snapshot>: something unexpected went wrong trying to make a snapshot.\n" 
            "              what(): "<< e.what() <<"\n";
            
            fStatus = kError; 
            PrintMsgAndAbort();
            return;
        }
    }

    
    //__________________________________________________________________________________________________________________________________
    //set whether an abort should be called when an error is encountered. 
    inline void SetAbortOnError(bool _val) { fAbortOnError=_val; }

    inline std::vector<std::string> GetOutputBranches() const { return fOutputBranches; }

    //__________________________________________________________________________________________________________________________________
    //assignment operator 
    inline ROOT::RDF::RNode& operator=(const ROOT::RDF::RNode& node) {
        fNode = node; 
        return fNode; 
    }
};

#endif 
