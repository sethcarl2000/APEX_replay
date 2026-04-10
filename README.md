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
