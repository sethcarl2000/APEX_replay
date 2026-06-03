#ifndef MetadataFetcher_h
#define MetadataFetcher_h

//ROOT headers
#include <ROOT/RDataFrame.hxx>
#include <TString.h>
//stdlib headers
#include <vector>
#include <string> 
#include <stdexcept>
#include <limits> 

//this struct will be used to associate meta-data with actual branch variables
struct metadata_branch_t {
  //three integers which identify a meta-data variable
  int run{-1}, segment{-1}, rawfile{-1};
  //the value to access
  double dt_center, dt_sigma;
};

//comparison operator
inline bool operator==(const metadata_branch_t& lhs, const metadata_branch_t& rhs) {
  return
    (lhs.run     == rhs.run) &&
    (lhs.segment == rhs.segment) &&
    (lhs.rawfile == rhs.rawfile);  
}

struct MetadataFetcher {
  
  //name of branch to produce 
  std::string fBranch;

  int fMinRun, fMaxRun; 
  
  std::vector<metadata_branch_t> fData;
  std::vector<int> fRunIndex; //each run's position in the 'data' vector
  
  MetadataFetcher(std::string branch, int min_run, int max_run)
    : fBranch{branch},
      fMinRun{min_run}, fMaxRun{max_run},
      fData{} {
    //initialize this list with all elements equal to -1 
    fRunIndex = std::vector<int>(fMaxRun-fMinRun+1, -1); 
  }; 
  
  //fetch meta-data value associated with the given run-segment-rawfile 
  metadata_branch_t Get(int run, int segment, int rawfile) const; 
};  

#endif

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ifndef MetadataFetcher_cxx
#define MetadataFetcher_cxx

//_________________________________________________________________________________________
metadata_branch_t MetadataFetcher::Get(int run, int segment, int rawfile) const
{
  static const double kNaN = std::numeric_limits<double>::quiet_NaN(); 
  
  metadata_branch_t target_data{ run, segment, rawfile, kNaN, kNaN }; 

  if (run > fMaxRun || run < fMinRun) {
    throw std::logic_error(Form("in <MetadataFetcher::Get> invalid run number: %i."
				" valid range is [%i,%i]",
				run, fMinRun, fMaxRun));
    return target_data;
  }

  //find our starting place to search
  int ind = fRunIndex[run-fMinRun];

  std::printf("<MetadataFetcher::Get>: looking for data:\n"
	      "  run=%i, segment=%i, rawfile=%i\n",
	      run, segment, rawfile
	      ); 
  
  std::printf("run %i has index: %i\n", run, ind); 
  
  //no meta-data for this run
  if (ind < 0) { return target_data; }

  auto data = fData[ind]; 
  while (ind < (int)fData.size() && data.run == run) {
    
    std::printf("data: run=%i, segment=%i, rawfile=%i, match? %s\n",
		data.run, data.segment, data.rawfile,
		(data==target_data ? "yes" : "no")); 
    
    //search all the meta-data for this run until we find our match 
    if (data == target_data) { return data; }
    ++ind;
    data = fData[ind]; 
  }

  return target_data; 
}
//_________________________________________________________________________________________
//_________________________________________________________________________________________
//_________________________________________________________________________________________
//_________________________________________________________________________________________
//_________________________________________________________________________________________
//_________________________________________________________________________________________
//_________________________________________________________________________________________
//_________________________________________________________________________________________


#endif
