# APEX_replay

## Bash scripts
Unless otherwise noted, these files are in the 'old-bash-scripts' directory

### cache_rawFile.sh
Cache raw (CODA) APEX runs from tape onto the cache directory. 
Usage example:
```tsch
./cache_rawFile.sh 4000 4050
```
This will request all (extant) runs from 4000 4050 (and each run 'part file' for runs that are split over multiple files). They will be cached in ```/cache/halla/apex/raw```. An email will be sent to my email ```sethhall@jlab.org``` upon completion of the task (All files are requested as one ```jcache``` job). 

### run_decode.sh 
This automates the process of sumbitting slurm jobs to decode raw APEX runs. I will copy-paste the help message printed by this script, invoked by ```./run_decode.sh -h```: 
```
run_decode.sh
  Runs decode script on desired cached run.
  Will use default options in decode_run.C, and submit them as jobs.

 Usage: run_decode.sh -r {1st-run} [-l {last-run}] [options]

 Options:
  -r  -  Run number to decode (required)

  -l  -  Last run number to decode (all runs in-between will be checked)

  -h  -  displays this help message, then quits

  -n  -  print the number of raw files this run currently has in the cache,
         then exit.

  -s  -  Print the number of raw files, and their total size in the cache (in MB)

  -p  -  Print the hard-coded parameters you want to run, then exit.

  -t  -  Submit jobs for test only, so that the queue can be tested.

  -z  -  Submit jobs in the priority queue (default is production)

  -m  -  Request memory-per-cpu in MB (use '-p' to see default)
```
FIXME: make this interactive, so that we don't have to modify the script itself on each unique invocation
