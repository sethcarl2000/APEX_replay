# APEX_replay

## Bash scripts
Executable bash scripts in the `scripts` directory

### run_full_replay.sh
Usage example: 
```bash
./scripts/run_full_replay.sh 3000 3050
```
This would submit runs 3000 to 3050 for a full replay, granted they are cached (see `cache_rawFile.sh`). 

Each APEX run has one or more raw files. Regardless of how many runs are submitted for a replay (be it 1 or all ~1000), This scripts scans the cache for rawfiles corresponding to a runs in the range `run0` to `run1`, and submits a **separate slurm-job for each rawfile**. 

For each run that succeeds, the final output file (with reconstructed VDC tracks, along with all other misc. information) will be under the directory `${PATH_APEX_VOLATILE}/production` with the following format: 

```
"replay.[run number].raw-num-[rawfile number].seg-[segment number].root"
```
Here, the rawfile number is the index of the rawfile (starting at 0 for the first rawfile encountered for a given run, and counting up from there). The segment number is the index of the segment (large rawfiles may be split into mutliple _segments_ during the decoding step if they are above a given size). After a replay is finished, the segment and rawfile number have no practical importance; they are just used to organize input and output files during the replay process. All replay files for a given run (or multiple runs) may be combined using `hadd` once a replay is done without any loss of information.


Under the hood, the `run_full_replay.sh` script creates a 'slurm array file', on the volatile disk under `${PATH_APEX_VOLATILE}/slurm/array_${run0}_${run1}.list`. This file contains information that the next script, `run-full-replay-arrray` will make use of. This `.list` file is `|`-delimited and has the following format: 

```
[slurm array id]|[run number]|[rawfile number]|[rawfile path]
```
Where the slurm array id starts at `0` for the first raw file in this run range, and the `rawfile number` index starts at 0 for the first rawfile belonging to each run, and counts upward for each addtional rawfile in a given run. 

The reason to bother with making such a meta-data csv is so that we can take advantage of slurm arrays. Arrays are a nice way to submit a large batch of similar jobs at once. My understanding is that the job-queueing system likes arrays over manually submitting hundreds or thousands of jobs independently, which gives you better throughput for large arrays of hundredes of jobs. 

### cache_rawFile.sh
Cache raw (CODA) APEX runs from tape onto the cache directory. 
Usage example:
```bash
./scripts/cache_rawFile.sh 4000 4050
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


