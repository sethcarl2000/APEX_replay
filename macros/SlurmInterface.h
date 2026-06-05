#ifndef SlurmInterface_h
#define SlurmInterface_h

#include <cstdlib> 
#include <string> 
#include <thread> 
#include <stdexcept> 

//simple singleton class that reports basic info about slurm job environment
class SlurmInterface {
 private:

  int job_id{-1}; 
  
  SlurmInterface() {
    //check to see if we're in a slurm job
    const char* job_id_str = std::getenv("SLURM_JOB_ID");
    if (job_id_str) {
      //we ARE in a slurm job
      job_id = std::stoi( std::string(job_id_str) ); 
    } else {
      //we are NOT in a slurm job
      job_id = -1; 
    }
  };

  ~SlurmInterface() {}; 
  
  
 public: 

  SlurmInterface(const SlurmInterface&) = delete;
  SlurmInterface& operator=(const SlurmInterface&) = delete;

  //static access fcn
  static SlurmInterface& Instance() {
    static SlurmInterface inst;
    return inst;
  }; 

  bool inside_slurm_job() const { return job_id >= 0; }
  
  int get_n_threads() const {
    if (inside_slurm_job()) {
      
      char* slurm_n_cpus = std::getenv("SLURM_CPUS_PER_TASK");
      
      if (!slurm_n_cpus) {
	throw std::logic_error("in <SlurmInterface::n_threads>: "
			       " we are supposedly in a slurm-job, but the env. 'SLURM_CPUS_PER_TASK' is undefined.");
	return 0; 
      }
      
      return std::stoi( std::string(slurm_n_cpus) ); 
    } 
    
    //if we aren't in a slurm-job, just return the hardware concurrency 
    return std::thread::hardware_concurrency();    
  }
  
  int get_job_id() const { return job_id; }
  
}; 


#endif
